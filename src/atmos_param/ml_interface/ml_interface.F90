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
    integer :: nml_unit, ierr, io

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

subroutine ENNUF_2d_T_RH_prediction(temp_in, q_in, num_levels, p_full, p_half, pert_t, pert_q)

    real, dimension(:,:,:), intent(in)     :: temp_in, q_in
    real, dimension(:,:,:), intent(in) :: p_full, p_half
    integer, intent(in)                  :: num_levels    
    real, dimension(:,:,:), intent(out)    :: pert_t, pert_q


    integer :: i, j, z_tick
    real, dimension(size(temp_in,1), size(temp_in, 2), 6) :: Th_predictors
    real, dimension(size(temp_in,1), size(temp_in, 2), 7) :: RH_predictors
    real(kind=4), dimension(size(temp_in,1), size(temp_in, 2)) :: Th_outputs, RH_outputs, RHstd, Thstd

    if(.not.module_is_initialized) then
        call error_mesg('ml_interface','ml_interface module is not initialized',FATAL)
    endif

    do i = 1, size(temp_in,1)
        do j = 1, size(temp_in,2)    

            RH_predictors(i,j,1) = 1.0 !temp_in(i,j,num_levels)
            RH_predictors(i,j,2) = 2.0 !temp_in(i,j,num_levels)
            RH_predictors(i,j,3) = 3.0 !q_in(i,j,num_levels)
            RH_predictors(i,j,4) = 4.0 !q_in(i,j,num_levels)
            RH_predictors(i,j,5) = 4.0 !q_in(i,j,num_levels)
            RH_predictors(i,j,6) = 4.0 !q_in(i,j,num_levels)
            RH_predictors(i,j,7) = 4.0 !q_in(i,j,num_levels)                                    

            call ENNUF_RH_model(real(RH_predictors(i,j,:),4), RH_outputs(i,j))

            RHstd(i,j) = RH_outputs(i,j)


            Th_predictors(i,j,1) = 1.0 !temp_in(i,j,num_levels)
            Th_predictors(i,j,2) = 2.0 !temp_in(i,j,num_levels)
            Th_predictors(i,j,3) = 3.0 !q_in(i,j,num_levels)
            Th_predictors(i,j,4) = 4.0 !q_in(i,j,num_levels)
            Th_predictors(i,j,5) = 4.0 !q_in(i,j,num_levels)
            Th_predictors(i,j,6) = 4.0 !q_in(i,j,num_levels)

            call ENNUF_Th_model(real(Th_predictors(i,j,:),4), Th_outputs(i,j))

            Thstd(i,j) = Th_outputs(i,j)

            write(6,*) RHstd(i,j), Thstd(i,j)
        enddo
    enddo

  do z_tick=1, num_levels
    pert_t(:,:,z_tick) = temp_in(:,:,z_tick) + Thstd*(p_full(:,:,z_tick)/p_half(:,:,num_levels+1))
    pert_q(:,:,z_tick) = q_in(:,:,z_tick)    + RHstd*(p_full(:,:,z_tick)/p_half(:,:,num_levels+1))
  enddo

end subroutine ENNUF_2d_T_RH_prediction


end module ml_interface_mod
