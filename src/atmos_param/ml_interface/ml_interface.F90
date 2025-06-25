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
use constants_mod, only: RVGAS, HLV

implicit none
private
!=================================================================================================================================

character(len=128) :: version= &
'$Id: ml_interface.F90,v 1.0'

character(len=128) :: tagname= &
'$Name:  $'
character(len=10), parameter :: mod_name='ml_interface'

!=================================================================================================================================

public :: ml_interface_init, read_2d_ml_generated_file, ENNUF_2d_test_prediction, ENNUF_2d_T_RH_prediction

character(len=256) :: conv_input_file  = 'ml_input'
character(len=256) :: tstd_field_name = 'tstd' 
character(len=256) :: qstd_field_name = 'qstd' 

namelist / ml_interface_nml / conv_input_file, tstd_field_name, qstd_field_name 


logical :: module_is_initialized =.false.
type(interpolate_type),save :: conv_input_file_interp


contains


subroutine ml_interface_init(is, ie, js, je, rad_lonb_2d, rad_latb_2d, perturb_ml_using_input_file)
    
    real, intent(in), dimension(:,:) :: rad_lonb_2d, rad_latb_2d
    integer, intent(in) :: is, ie, js, je
    logical, intent(in) :: perturb_ml_using_input_file
    integer :: nml_unit, ierr, io, nseed

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

subroutine ENNUF_2d_T_RH_prediction(temp_in, q_in, ptemp_2m, rh_2m, sdor, surf_geopot, temp_2m, u_10m, v_10m, num_levels, p_full, p_half, pert_t, pert_q)

    real, dimension(:,:,:), intent(in)     :: temp_in, q_in
    real, dimension(:,:),   intent(in)     :: ptemp_2m, rh_2m, sdor, surf_geopot, temp_2m, u_10m, v_10m
    real, dimension(:,:,:), intent(in)     :: p_full, p_half
    integer, intent(in)                    :: num_levels    
    real, dimension(:,:,:), intent(out)    :: pert_t, pert_q


    integer :: i, j, z_tick
    real, dimension(size(temp_in,1), size(temp_in, 2))    :: tdewpoint_2m
    real, dimension(size(temp_in,1), size(temp_in, 2), 6) :: Th_predictors
    real, dimension(size(temp_in,1), size(temp_in, 2), 7) :: RH_predictors
    real(kind=4), dimension(size(temp_in,1), size(temp_in, 2)) :: Th_outputs, RH_outputs, qstd, Tstd, random_num_theta, random_num_rh

    real :: ptemp_2m_mean, ptemp_2m_std, rh2_mean, rh2_std, sdor_mean, sdor_std, surf_geopot_mean, surf_geopot_std, wind_10m_mean, wind_10m_std, &
            tdewpoint_2m_mean, tdewpoint_2m_std, temp_2m_mean, temp_2m_std, surf_p_mean, surf_p_std

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

    tdewpoint_2m = temp_2m / (1. - ((RVGAS/HLV)*temp_2m*log(rh_2m)))

    do i = 1, size(temp_in,1)
        do j = 1, size(temp_in,2)    

            RH_predictors(i,j,1) = (ptemp_2m(i,j)-ptemp_2m_mean)/ptemp_2m_std !ptemp_2m 
            RH_predictors(i,j,2) = (rh_2m(i,j)-rh2_mean)/rh2_std !RH2
            RH_predictors(i,j,3) = (sdor(i,j)-sdor_mean)/sdor_std !sdor
            RH_predictors(i,j,4) = (surf_geopot(i,j)-surf_geopot_mean)/(surf_geopot_std) !surface geopotential
            RH_predictors(i,j,5) = (tdewpoint_2m(i,j)-tdewpoint_2m_mean)/tdewpoint_2m_std !2m dew-point temperature
            RH_predictors(i,j,6) = (temp_2m(i,j)-temp_2m_mean)/temp_2m_std !2m temperature
            RH_predictors(i,j,7) = (p_half(i,j,num_levels+1)-surf_p_mean)/surf_p_std !surface pressure

            call ENNUF_RH_model(real(RH_predictors(i,j,:),4), RH_outputs(i,j))

            qstd(i,j) = RH_outputs(i,j)


            Th_predictors(i,j,1) = (ptemp_2m(i,j)-ptemp_2m_mean)/ptemp_2m_std !ptemp_2m 
            Th_predictors(i,j,2) = (rh_2m(i,j)-rh2_mean)/rh2_std !RH2
            Th_predictors(i,j,3) = (sdor(i,j)-sdor_mean)/sdor_std !sdor
            Th_predictors(i,j,4) = (surf_geopot(i,j)-surf_geopot_mean)/(surf_geopot_std) !surface geopotential
            Th_predictors(i,j,5) = (p_half(i,j,num_levels+1)-surf_p_mean)/surf_p_std !surface pressure
            Th_predictors(i,j,6) = (sqrt(u_10m(i,j)**2 + v_10m(i,j)**2)-wind_10m_mean)/wind_10m_std !10m wind speed

            call ENNUF_Th_model(real(Th_predictors(i,j,:),4), Th_outputs(i,j))

            Tstd(i,j) = Th_outputs(i,j)

            ! make sure random number is drawn from normal distribution with mean of 0 and std of 1. And now clip it between -3 and 3 standard deviations.
            random_num_theta(i,j) = MIN(rnorm(), 3.0)
            random_num_rh(i,j)    = MIN(rnorm(), 3.0)
            random_num_theta(i,j) = MAX(random_num_theta(i,j), -3.0)
            random_num_rh(i,j)    = MAX(random_num_rh(i,j), -3.0)


            write(6,*) qstd(i,j), Tstd(i,j), random_num_theta(i,j), random_num_rh(i,j)

        enddo
    enddo

    !TODO:
    !Need to think about how to go from predictions of theta and RH to perturbing T and q

  do z_tick=1, num_levels
    pert_t(:,:,z_tick) = temp_in(:,:,z_tick) + random_num_theta* Tstd*(p_full(:,:,z_tick)/p_half(:,:,num_levels+1))
    pert_q(:,:,z_tick) = q_in(:,:,z_tick)    + random_num_rh*    qstd*(p_full(:,:,z_tick)/p_half(:,:,num_levels+1))
  enddo

    !Need to clip T and q perturbations so that they're not crazy big


end subroutine ENNUF_2d_T_RH_prediction


end module ml_interface_mod
