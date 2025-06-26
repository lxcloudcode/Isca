module ml_interface_mod

#ifdef INTERNAL_FILE_NML
    use mpp_mod, only: input_nml_file
#else
    use fms_mod, only: open_namelist_file, close_file
#endif

use fms_mod, only: write_version_number, file_exist, close_file, stdlog, error_mesg, NOTE, FATAL, read_data, field_size, uppercase, mpp_pe, check_nml_error, mpp_root_pe

use time_manager_mod, only: time_type

use interpolator_mod, only: interpolate_type,interpolator_init&
     &,CONSTANT,interpolator

use ennuf_example_mod, only: example_ml_model
use ENNUF_RH_mod, only: ENNUF_RH_model
use ENNUF_Th_mod, only: ENNUF_Th_model
use constants_mod, only: cp_air, hlv, stefan, rdgas, rvgas, grav, vonkarm, dens_h2o, PSTD_MKS
use sat_vapor_pres_mod, only: lookup_es ! needed for relative humdity to be converted to q
use   monin_obukhov_mod, only: mo_drag, mo_profile
use  sat_vapor_pres_mod, only: escomp, descomp
use        diag_manager_mod, only: register_diag_field, send_data
use   spectral_dynamics_mod, only: get_axis_id
use        time_manager_mod, only: time_type

implicit none
private
!=================================================================================================================================

character(len=128) :: version= &
'$Id: ml_interface.F90,v 1.0'

character(len=128) :: tagname= &
'$Name:  $'
character(len=15), parameter :: mod_name='ml_interface'

!=================================================================================================================================

public :: ml_interface_init, read_2d_ml_generated_file, ENNUF_2d_test_prediction, ENNUF_2d_T_RH_prediction

character(len=256) :: conv_input_file  = 'ml_input'
character(len=256) :: tstd_field_name = 'tstd' 
character(len=256) :: qstd_field_name = 'qstd' 
logical :: alt_gustiness         = .false.
logical :: bucket = .true.
logical :: raoult_sat_vap        = .false.
logical :: do_simple             = .false.
real    :: gust_const            =  1.0
real    :: gust_min              =  0.0
logical :: no_neg_q              = .false.  ! for backwards compatibility
logical :: use_mixing_ratio      = .false.

namelist / ml_interface_nml / conv_input_file, tstd_field_name, qstd_field_name, &
                              alt_gustiness, bucket, do_simple, gust_const, gust_min, &
                              no_neg_q, raoult_sat_vap, use_mixing_ratio


logical :: module_is_initialized =.false.
type(interpolate_type),save :: conv_input_file_interp
real :: d622, d378, hlars, gcp, kappa, d608

integer ::           &
     id_theta_std,  &   ! predicted theta std
     id_t_std,  &   ! predicted t std     
     id_RH_std,  &   ! predicted theta std
     id_q_std   ! predicted t std     

integer, dimension(4) :: axes

contains


subroutine ml_interface_init(is, ie, js, je, rad_lonb_2d, rad_latb_2d, perturb_ml_using_input_file, Time)
    
    real, intent(in), dimension(:,:) :: rad_lonb_2d, rad_latb_2d
    integer, intent(in) :: is, ie, js, je
    logical, intent(in) :: perturb_ml_using_input_file
    integer :: nml_unit, ierr, io, nseed
    type(time_type) :: Time

    if(module_is_initialized) return

#ifdef INTERNAL_FILE_NML
        read (input_nml_file, nml=ml_interface_nml, iostat=io)
        ierr = check_nml_error (io,'ml_interface_nml')
#else
        if ( file_exist('input.nml') ) then
            nml_unit = open_namelist_file()
            read (nml_unit, ml_interface_nml, iostat=io)
            ierr = check_nml_error (io,'ml_interface_nml')
            call close_file(nml_unit)
        endif
#endif

    if ( mpp_pe() == mpp_root_pe() ) write (stdlog(), nml=ml_interface_nml)

    if (perturb_ml_using_input_file) then
        call interpolator_init( conv_input_file_interp, trim(conv_input_file)//'.nc', rad_lonb_2d, rad_latb_2d, data_out_of_bounds=(/CONSTANT/) )
    endif

    call random_seed(size=nseed)

    ! initialise variables for rh_calc
    d622 = RDGAS/RVGAS
    d378 = 1.-d622    
    hlars  = hlv/rvgas
    gcp    = grav/cp_air
    kappa  = rdgas/cp_air
    d608   = d378/d622

  axes = get_axis_id()

  id_theta_std = register_diag_field(mod_name, 'theta_std_ml', &
                    axes(1:2), Time, 'predicted std of Theta from ML', 'K')
  id_t_std = register_diag_field(mod_name, 't_std_ml', &
                    axes(1:2), Time, 'predicted std of T from ML', 'K')
  id_RH_std = register_diag_field(mod_name, 'RH_std_ml', &
                    axes(1:2), Time, 'predicted std of RH from ML', 'percent')
  id_q_std = register_diag_field(mod_name, 'q_std_ml', &
                    axes(1:2), Time, 'predicted std of sphum from ML', 'kg/kg')                                                            

    module_is_initialized=.true.

    return
end subroutine ml_interface_init


! subroutine read_ml_generated_file(p_half, num_levels, tstd, qstd)

!     real, dimension(:,:,:), intent(in)  :: p_half
!     integer, intent(in):: num_levels    
!     real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)), intent(out)                   :: tstd, qstd
!     real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)) :: sigma_half    

!     if(.not.module_is_initialized) then
!         call error_mesg('ml_interface','ml_interface module is not initialized',FATAL)
!       endif

!     do i in range(size(p_half,1))
!       do j in range(size(p_half,2))
!         sigma_half(i,j,:) = p_half(i,j,:) /p_half(i,j,num_levels+1) 
!       enddo
!     enddo

!     call interpolator( conv_input_file_interp, sigma_half, tstd, tstd_field_name)
!     call interpolator( conv_input_file_interp, sigma_half, qstd, qstd_field_name)

! end subroutine read_ml_generated_file

subroutine read_2d_ml_generated_file(tstd)

    real, dimension(:,:), intent(out)                   :: tstd

    if(.not.module_is_initialized) then
        call error_mesg('ml_interface','ml_interface module is not initialized',FATAL)
    endif

    call interpolator( conv_input_file_interp, tstd, tstd_field_name)

end subroutine read_2d_ml_generated_file

subroutine ENNUF_2d_test_prediction(temp_in, q_in, tstd)

    real, dimension(:,:), intent(in)   :: temp_in, q_in
    real, dimension(:,:), intent(out)  :: tstd

    integer :: i, j
    real, dimension(size(temp_in,1), size(temp_in, 2), 4) :: four_predictors
    real(kind=4), dimension(size(temp_in,1), size(temp_in, 2), 2) :: two_outputs

    if(.not.module_is_initialized) then
        call error_mesg('ml_interface','ml_interface module is not initialized',FATAL)
    endif

    do i = 1, size(temp_in,1)
        do j = 1, size(temp_in,2)    

            four_predictors(i,j,1) = 1.0 !temp_in(i,j)
            four_predictors(i,j,2) = 2.0 !temp_in(i,j)
            four_predictors(i,j,3) = 3.0 !q_in(i,j)
            four_predictors(i,j,4) = 4.0 !q_in(i,j)

            call example_ml_model(real(four_predictors(i,j,:),4), two_outputs(i,j,:))

            tstd(i,j) = two_outputs(i,j,1)
            write(6,*) real(four_predictors(i,j,:),4), two_outputs(i,j,1), two_outputs(i,j,2)
        enddo
    enddo

end subroutine ENNUF_2d_test_prediction

function rnorm() result( fn_val )
    !   Adapted from code released into the public domain by Alan Miller
    !   https://jblevins.org/mirror/amiller/rnorm.f90
    !   This version doubles the computations required for many calls
    !   but is thread-safe.
    !   Generate a random normal deviate using the polar method.
    !   Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
    !              normal variables', Siam Rev., vol.6, 260-264, 1964.

    implicit none
    real(kind=4)  :: fn_val

    ! Local variables
    real(kind=4)            :: u, v, sum, sln
    real(kind=4), parameter :: one = 1.0, vsmall = tiny( one )

    ! Generate a pair of random normals
    do
    call random_number( u )
    call random_number( v )
    u = scale( u, 1 ) - one
    v = scale( v, 1 ) - one
    sum = u*u + v*v + vsmall         ! vsmall added to prevent LOG(zero) / zero
    if (sum < one) exit
    end do
    sln = sqrt(- scale( log(sum), 1 ) / sum)
    fn_val = u*sln
return
end function rnorm

subroutine ENNUF_2d_T_RH_prediction(temp_in, q_in, u_in, v_in, t_surf_in, q_surf_in, u_surf_in, v_surf_in, rough_mom, rough_heat, rough_moist, rough_scale, gust, bucket_depth_in, land, seawater, avail_in, Time, ptemp_2m, rh_2m, sdor, surf_geopot, temp_2m, u_10m, v_10m, num_levels, p_full, p_half, z_full, z_surf, pert_t, pert_q)

    real, dimension(:,:,:), intent(in)     :: temp_in, q_in, u_in, v_in
    real, dimension(:,:),   intent(in)     :: ptemp_2m, rh_2m, sdor, surf_geopot, temp_2m, u_10m, v_10m, z_surf, t_surf_in, q_surf_in, u_surf_in, v_surf_in
    real, dimension(:,:),   intent(in)     :: rough_mom, rough_heat, rough_moist, rough_scale, gust, bucket_depth_in    
    real, dimension(:,:,:), intent(in)     :: p_full, p_half, z_full
    integer, intent(in)                    :: num_levels    
    logical, dimension(:,:), intent(in)    :: land, seawater, avail_in
    type(time_type),            intent(in) :: Time
    real, dimension(:,:,:), intent(out)    :: pert_t, pert_q


    integer :: i, j, z_tick
    real, dimension(size(temp_in,1), size(temp_in, 2))    :: tdewpoint_2m, q_sat, es_2m
    real, dimension(size(temp_in,1), size(temp_in, 2))    :: temp_2m_sv, u_10m_sv, v_10m_sv, q_2m_sv, rh_2m_sv, ptemp_2m_sv

    real, dimension(size(temp_in,1), size(temp_in, 2), 6) :: Th_predictors
    real, dimension(size(temp_in,1), size(temp_in, 2), 7) :: RH_predictors
    real(kind=4), dimension(size(temp_in,1), size(temp_in, 2)) :: Th_outputs, RH_outputs, qstd, Tstd, random_num_theta, random_num_rh

    real :: ptemp_2m_mean, ptemp_2m_std, rh2_mean, rh2_std, sdor_mean, sdor_std, surf_geopot_mean, surf_geopot_std, wind_10m_mean, wind_10m_std, &
            tdewpoint_2m_mean, tdewpoint_2m_std, temp_2m_mean, temp_2m_std, surf_p_mean, surf_p_std, test_input_value

    logical :: used

    if(.not.module_is_initialized) then
        call error_mesg('ml_interface','ml_interface module is not initialized',FATAL)
    endif

    ptemp_2m_mean    = 281.92596
    ptemp_2m_std     = 17.551208
    rh2_mean         = 74.51201 / 100.
    rh2_std          = 14.990538 / 100.
    sdor_mean        = 20.685488
    sdor_std         = 45.847626
    surf_geopot_mean = 3613.3135
    surf_geopot_std  = 7890.1587
    wind_10m_mean    = 5.9862623
    wind_10m_std     = 3.5270183
    tdewpoint_2m_mean= 274.4919
    tdewpoint_2m_std = 20.493927
    temp_2m_mean     = 279.14786
    temp_2m_std      = 21.044209
    surf_p_mean      = 96789.875
    surf_p_std       = 9142.519

    call calc_surface_variables(temp_in(:,:,num_levels),     q_in(:,:,num_levels),   u_in(:,:,num_levels),     v_in(:,:,num_levels),     p_full(:,:,num_levels),  &   
                                  z_full(:,:,num_levels)-z_surf(:,:),    &
                                  p_half(:,:,num_levels+1),    t_surf_in(:,:),     t_surf_in(:,:),      q_surf_in(:,:), u_surf_in,    v_surf_in,      &
                                  rough_mom, rough_heat, rough_moist, rough_scale, gust, bucket_depth_in, land, seawater, avail_in,        &                          
                                  temp_2m_sv, u_10m_sv, v_10m_sv,                                            &
                                  q_2m_sv, rh_2m_sv, ptemp_2m_sv)


    tdewpoint_2m = temp_2m_sv / (1. - ((RVGAS/HLV)*temp_2m_sv*log(rh_2m_sv)))

    ! write(6,*) MAXVAL(temp_2m), MINVAL(temp_2m), MAXVAL(temp_2m_sv), MINVAL(temp_2m_sv)
    ! write(6,*) MAXVAL(ptemp_2m), MINVAL(ptemp_2m), MAXVAL(ptemp_2m_sv), MINVAL(ptemp_2m_sv)
    ! write(6,*) MAXVAL(rh_2m), MINVAL(rh_2m), MAXVAL(rh_2m_sv), MINVAL(rh_2m_sv)
    ! ! write(6,*) MAXVAL(q_2m), MINVAL(q_2m), MAXVAL(q_2m_sv), MINVAL(q_2m_sv)
    ! write(6,*) MAXVAL(u_10m), MINVAL(u_10m), MAXVAL(u_10m_sv), MINVAL(u_10m_sv)
    ! write(6,*) MAXVAL(v_10m), MINVAL(v_10m), MAXVAL(v_10m_sv), MINVAL(v_10m_sv)

    CALL LOOKUP_ES(temp_2m_sv, es_2m)  !find saturation vapor pressure
    q_sat = d622 * es_2m / (p_half(i,j,num_levels+1) - d378*es_2m)

    Th_outputs = 0.
    RH_outputs = 0.

    test_input_value = 0.2

    do i = 1, size(temp_in,1)
        do j = 1, size(temp_in,2)    

            RH_predictors(i,j,1) = 0.15!test_input_value!(ptemp_2m_sv(i,j)-ptemp_2m_mean)/ptemp_2m_std !ptemp_2m 
            RH_predictors(i,j,2) = -0.03!test_input_value!(rh_2m_sv(i,j)-rh2_mean)/rh2_std !RH2
            RH_predictors(i,j,3) = 0.07!test_input_value!(sdor(i,j)-sdor_mean)/sdor_std !sdor
            RH_predictors(i,j,4) = -0.11!test_input_value!(surf_geopot(i,j)-surf_geopot_mean)/(surf_geopot_std) !surface geopotential
            RH_predictors(i,j,5) = 0.01!test_input_value!(tdewpoint_2m(i,j)-tdewpoint_2m_mean)/tdewpoint_2m_std !2m dew-point temperature
            RH_predictors(i,j,6) = -0.09!test_input_value!(temp_2m_sv(i,j)-temp_2m_mean)/temp_2m_std !2m temperature
            RH_predictors(i,j,7) = 0.17!test_input_value!(p_half(i,j,num_levels+1)-surf_p_mean)/surf_p_std !surface pressure

            call ENNUF_RH_model(real(RH_predictors(i,j,:),4), RH_outputs(i,j))

            qstd(i,j) = RH_outputs(i,j)*q_sat(i,j)


            Th_predictors(i,j,1) = 0.15!test_input_value!(ptemp_2m_sv(i,j)-ptemp_2m_mean)/ptemp_2m_std !ptemp_2m 
            Th_predictors(i,j,2) = -0.03!test_input_value!(rh_2m_sv(i,j)-rh2_mean)/rh2_std !RH2
            Th_predictors(i,j,3) = 0.07!test_input_value!(sdor(i,j)-sdor_mean)/sdor_std !sdor
            Th_predictors(i,j,4) = -0.11!test_input_value!(surf_geopot(i,j)-surf_geopot_mean)/(surf_geopot_std) !surface geopotential
            Th_predictors(i,j,5) = 0.01!test_input_value!(p_half(i,j,num_levels+1)-surf_p_mean)/surf_p_std !surface pressure
            Th_predictors(i,j,6) = -0.09!test_input_value!(sqrt(u_10m_sv(i,j)**2 + v_10m_sv(i,j)**2)-wind_10m_mean)/wind_10m_std !10m wind speed

            call ENNUF_Th_model(real(Th_predictors(i,j,:),4), Th_outputs(i,j))

            Tstd(i,j) = Th_outputs(i,j) * (p_half(i,j,num_levels+1)/PSTD_MKS)**KAPPA

            ! make sure random number is drawn from normal distribution with mean of 0 and std of 1. And now clip it between -3 and 3 standard deviations.
            random_num_theta(i,j) = MIN(rnorm(), 3.0)
            random_num_rh(i,j)    = MIN(rnorm(), 3.0)
            random_num_theta(i,j) = MAX(random_num_theta(i,j), -3.0)
            random_num_rh(i,j)    = MAX(random_num_rh(i,j), -3.0)


            write(6,*) RH_outputs(i,j), TH_outputs(i,j), qstd(i,j), Tstd(i,j)

        enddo
    enddo

  write(6,*) 'maxvals = ', maxval(RH_outputs), minval(RH_outputs), maxval(Th_outputs), minval(Th_outputs)

  if(id_theta_std > 0) used = send_data(id_theta_std, real(Th_outputs,8), Time)
  if(id_t_std > 0)     used = send_data(id_t_std, real(Tstd,8), Time)
  if(id_RH_std > 0)    used = send_data(id_theta_std, real(RH_outputs*100.,8), Time)
  if(id_q_std > 0)     used = send_data(id_t_std, real(qstd,8), Time)

  pert_t(:,:,:) = temp_in(:,:,:)
  pert_q(:,:,:) = q_in(:,:,:)

  ! do z_tick=1, num_levels
  !   pert_t(:,:,z_tick) = temp_in(:,:,z_tick) + 0.01*random_num_theta* Tstd*(p_full(:,:,z_tick)/p_half(:,:,num_levels+1))
  !   pert_q(:,:,z_tick) = q_in(:,:,z_tick)    + random_num_rh*    qstd*(p_full(:,:,z_tick)/p_half(:,:,num_levels+1))
  ! enddo

    pert_t(:,:,num_levels) = pert_t(:,:,num_levels) + 0.001*real(random_num_theta* Tstd,8)
    pert_q(:,:,num_levels) = pert_q(:,:,num_levels) + 0.001*real(random_num_theta* qstd,8)

    !Need to clip T and q perturbations so that they're not crazy big


end subroutine ENNUF_2d_T_RH_prediction

subroutine calc_surface_variables(t_atm,     q_atm_in,   u_atm,     v_atm,     p_atm,     z_atm,    &
                                  p_surf,    t_surf,     t_ca,      q_surf, u_surf,    v_surf,      &
                                  rough_mom, rough_heat, rough_moist, rough_scale, gust,            &      
                                  bucket_depth, land, seawater, avail,                                                            &                    
                                  temp_2m, u_10m, v_10m,                                            &
                                  q_2m, rh_2m, ptemp_2m &                                     
                                  )

  real, intent(in),  dimension(:,:) :: &
       t_atm,     q_atm_in,   u_atm,     v_atm,              &
       p_atm,     z_atm,      t_ca,                          &
       p_surf,    t_surf,     u_surf,    v_surf,             &
       q_surf, &
       rough_mom, rough_heat, rough_moist, rough_scale, gust, &
       bucket_depth    
  logical, intent(in), dimension(:,:) :: land, seawater, avail
  real, intent(out),  dimension(:,:) :: &       
       temp_2m, u_10m, v_10m,                                &
       q_2m, rh_2m, ptemp_2m 

 integer :: j

  do j = 1, size(t_atm,2)
    call calc_surface_variables_1d(t_atm(:,j),     q_atm_in(:,j),   u_atm(:,j),     v_atm(:,j),     p_atm(:,j),     z_atm(:,j),    &
                                  p_surf(:,j),    t_surf(:,j),     t_ca(:,j),      q_surf(:,j), u_surf(:,j),    v_surf(:,j),      &
                                  rough_mom(:,j), rough_heat(:,j), rough_moist(:,j), rough_scale(:,j), gust(:,j),            &      
                                  bucket_depth(:,j), land(:,j), seawater(:,j), avail(:,j),                                                            &                    
                                  temp_2m(:,j), u_10m(:,j), v_10m(:,j),                                            &
                                  q_2m(:,j), rh_2m(:,j), ptemp_2m(:,j) &                                     
                                  )
  enddo

end subroutine calc_surface_variables


subroutine calc_surface_variables_1d(t_atm,     q_atm_in,   u_atm,     v_atm,     p_atm,     z_atm,    &
                                  p_surf,    t_surf,     t_ca,      q_surf, u_surf,    v_surf,      &
                                  rough_mom, rough_heat, rough_moist, rough_scale, gust,            &      
                                  bucket_depth, land, seawater, avail,                                                            &                    
                                  temp_2m, u_10m, v_10m,                                            &
                                  q_2m, rh_2m, ptemp_2m &                                     
                                  )

  real, intent(in),  dimension(:) :: &
       t_atm,     q_atm_in,   u_atm,     v_atm,              &
       p_atm,     z_atm,      t_ca,                          &
       p_surf,    t_surf,     u_surf,    v_surf,             &
       q_surf, &
       rough_mom, rough_heat, rough_moist, rough_scale, gust, &
       bucket_depth    
  logical, intent(in), dimension(:) :: land, seawater, avail
  real, intent(out),  dimension(:) :: &       
       temp_2m, u_10m, v_10m,                                &
       q_2m, rh_2m, ptemp_2m 

  real, dimension(size(t_atm,1)) :: &
      w_atm,     u_star,     b_star,    q_star,             &
      cd_m,      cd_t,       cd_q, &        
      thv_atm,  th_atm,   tv_atm,    thv_surf,            &
      e_sat,    e_sat1,   q_sat,     q_sat1,    p_ratio,  &
      t_surf0,  t_surf1,  u_dif,     v_dif,               &
      rho_drag, drag_t,    drag_m,   drag_q,    rho,      &
      q_atm,    q_surf0,  dw_atmdu,  dw_atmdv,  w_gust,   &
      e_sat_2m, q_sat_2m, ex_del_h, ex_del_m, ex_del_q

  real, parameter:: del_temp=0.1, del_temp_inv=1.0/del_temp
  real :: zrefm, zrefh      
  integer :: i

!---- use local value of surf temp ----

  t_surf0 = 200.   !  avoids out-of-bounds in es lookup
  where (avail)
     where (land)
        t_surf0 = t_ca
     elsewhere
        t_surf0 = t_surf
     endwhere
  endwhere

  t_surf1 = t_surf0 + del_temp

  call escomp ( t_surf0, e_sat  )  ! saturation vapor pressure
  call escomp ( t_surf1, e_sat1 )  ! perturbed  vapor pressure

  if(use_mixing_ratio) then
    ! surface mixing ratio at saturation
    q_sat   = d622*e_sat /(p_surf-e_sat )
    q_sat1  = d622*e_sat1/(p_surf-e_sat1)
  elseif(do_simple) then                  !rif:(09/02/09)
    q_sat   = d622*e_sat / p_surf
    q_sat1  = d622*e_sat1/ p_surf
  else
    ! surface specific humidity at saturation
    q_sat   = d622*e_sat /(p_surf-d378*e_sat )
    q_sat1  = d622*e_sat1/(p_surf-d378*e_sat1)
  endif

  ! initilaize surface air humidity according to surface type
  where (land)
!     q_surf0 = q_surf ! land calculates it
     q_surf0 = q_sat ! our simplified land evaporation model does not calculate q_surf, so we specify it as q_sat.
  elsewhere
     q_surf0 = q_sat  ! everything else assumes saturated sfc humidity
  endwhere

  if (raoult_sat_vap) where (seawater) q_surf0 = 0.98 * q_surf0

  ! check for negative atmospheric humidities
  where(avail) q_atm = q_atm_in
  if(no_neg_q) then
     where(avail .and. q_atm_in < 0.0) q_atm = 0.0
  endif

  ! initialize surface air humidity depending on whether surface is dry or wet (bucket empty or not)
  if (bucket) then
  where (bucket_depth <= 0.0)
      q_surf0 = q_atm
  elsewhere
      q_surf0 = q_sat    ! everything else assumes saturated sfc humidity
  end where
  endif

  ! generate information needed by monin_obukhov
  where (avail)
     p_ratio = (p_surf/p_atm)**kappa

     tv_atm  = t_atm  * (1.0 + d608*q_atm)     ! virtual temperature
     th_atm  = t_atm  * p_ratio                ! potential T, using p_surf as reference
     thv_atm = tv_atm * p_ratio                ! virt. potential T, using p_surf as reference
     thv_surf= t_surf0 * (1.0 + d608*q_surf0 ) ! surface virtual (potential) T
!     thv_surf= t_surf0                        ! surface virtual (potential) T -- just for testing tun off the q_surf

     u_dif = u_surf - u_atm                    ! velocity components relative to surface
     v_dif = v_surf - v_atm
  endwhere

  if(alt_gustiness) then
     do i = 1, size(avail)
        if (.not.avail(i)) cycle
        w_atm(i) = max(sqrt(u_dif(i)**2 + v_dif(i)**2), gust_const)
        ! derivatives of surface wind w.r.t. atm. wind components
        if(w_atm(i) > gust_const) then
           dw_atmdu(i) = u_dif(i)/w_atm(i)
           dw_atmdv(i) = v_dif(i)/w_atm(i)
        else
           dw_atmdu(i) = 0.0
           dw_atmdv(i) = 0.0
        endif
     enddo
  else
     if (gust_min > 0.0) then
       where(avail)
         w_gust = max(gust,gust_min) ! minimum gustiness
       end where
     else
       where(avail)
         w_gust = gust
       end where
     endif

     where(avail)
        w_atm = sqrt(u_dif*u_dif + v_dif*v_dif + w_gust*w_gust)
        ! derivatives of surface wind w.r.t. atm. wind components
        dw_atmdu = u_dif/w_atm
        dw_atmdv = v_dif/w_atm
     endwhere
  endif

  !  monin-obukhov similarity theory
  call mo_drag (thv_atm, thv_surf, z_atm,                  &
       rough_mom, rough_heat, rough_moist, w_atm,          &
       cd_m, cd_t, cd_q, u_star, b_star, avail             )

! added for 10m winds and 2m temperature add mo_profile()

  zrefm = 10. !want winds at 10m
  zrefh = 2.  !want temp and q at 2m

  call mo_profile( zrefm, zrefh, z_atm,   rough_mom, &
       rough_heat, rough_moist,          &
       u_star, b_star, q_star,        &
       ex_del_m, ex_del_h, ex_del_q, avail  )

! adapted from https://github.com/mom-ocean/MOM5/blob/3702ad86f9653f4e315b98613eb824a47d89cf00/src/coupler/flux_exchange.F90#L1932

     !    ------- reference temp -----------
        where (avail) &
           temp_2m = t_surf + (t_atm - t_surf) * ex_del_h !t_ca = canopy temperature, assuming that there is no canopy (no difference between land and ocean), t_ca = t_surf

     !    ------- reference ptemp -----------
        where (avail) &
           ptemp_2m = (t_surf*(PSTD_MKS/p_surf)**kappa) + ((t_atm*(PSTD_MKS/p_atm)**kappa) - (t_surf*(PSTD_MKS/p_surf)**kappa)) * ex_del_h 

     !    ------- reference u comp -----------
        where (avail) &
           u_10m = u_atm * ex_del_m ! setting u at surface to 0.

     !    ------- reference v comp -----------
       where (avail) &
           v_10m = v_atm * ex_del_m ! setting v at surface to 0.

! end of low level wind additions


      ! Add 2m q and RH

      ! ------- reference specific humidity -----------
      where (avail) &
           q_2m = q_surf + (q_atm - q_surf) * ex_del_q

      call escomp ( temp_2m, e_sat_2m )

      if(use_mixing_ratio) then
         ! surface mixing ratio at saturation
         q_sat_2m   = d622 * e_sat_2m / (p_surf - e_sat_2m)
      elseif(do_simple) then
         q_sat_2m   = d622 * e_sat_2m / p_surf
      else
         ! surface specific humidity at saturation
         q_sat_2m   = d622 * e_sat_2m / (p_surf - d378*e_sat)
      endif

      ! ------- reference relative humidity -----------
      where (avail) &
         rh_2m = q_2m / q_sat_2m
end subroutine calc_surface_variables_1d

end module ml_interface_mod
