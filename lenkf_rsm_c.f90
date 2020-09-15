module lenkf_rsm_c
  use lenkf_rsm, ONLY: lenkf_rsm_fortran => lenkf_analysis_rsm
  use mod_base_assimilation_manager, ONLY: base_assimilation_manager
  use exceptions, ONLY: error_status

  use iso_c_binding

  implicit none

  type, extends(base_assimilation_manager):: c_function_container
    type(c_ptr)::info_ptr
    type(c_funptr)::localize_fptr
    type(c_funptr)::add_obs_err_fptr
  contains
    procedure::localize => localize_wrapper
    procedure::add_obs_err => add_obs_err_wrapper
  end type c_function_container

contains

  subroutine localize_wrapper(this, &
                              istep, ibatch, dim_p, dim_obs, HP_p, HPH, status)

    abstract INTERFACE
      SUBROUTINE localize( &
        istep, ind_p, dim_p, dim_obs, HP_p, HPH, info_ptr) BIND(C)
        ! Apply localization to HP and HPH^T
        USE iso_c_binding
        INTEGER(c_int32_t), INTENT(in), value :: istep, ind_p, dim_p, dim_obs
        REAL(c_double), INTENT(inout) :: HP_p(dim_obs, dim_p)
        REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)
        type(c_ptr), value::info_ptr
      END SUBROUTINE localize
    end INTERFACE

    class(c_function_container), target :: this
    INTEGER(c_int32_t), INTENT(in), value :: istep, ibatch, dim_p, dim_obs
    REAL(c_double), INTENT(inout) :: HP_p(dim_obs, dim_p), HPH(dim_obs, dim_obs)
    class(error_status), intent(out), allocatable, optional :: status
    procedure(localize), pointer::U_localize

    call c_f_procpointer(this%localize_fptr, U_localize)
    call U_localize(istep, ibatch, dim_p, dim_obs, HP_p, HPH, this%info_ptr)

  end subroutine localize_wrapper

  subroutine add_obs_err_wrapper(this, istep, ibatch, dim_obs, HPH, status)

    abstract INTERFACE
      SUBROUTINE add_obs_err(step, ind_p, dim_obs, HPH, info_ptr) BIND(C)
        ! Add observation error covariance matrix
        USE iso_c_binding
        INTEGER(c_int32_t), INTENT(in), value :: step, ind_p, dim_obs
        REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)
        type(c_ptr), value::info_ptr
      END SUBROUTINE add_obs_err
    end INTERFACE

    class(c_function_container), target::this
    INTEGER(c_int32_t), INTENT(in), value :: istep, ibatch, dim_obs
    REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)
    class(error_status), intent(out), allocatable, optional :: status
    procedure(add_obs_err), pointer::U_add_obs_err

    call c_f_procpointer(this%add_obs_err_fptr, U_add_obs_err)
    call U_add_obs_err(istep, ibatch, dim_obs, HPH, this%info_ptr)
  end subroutine add_obs_err_wrapper

  subroutine lenkf_analysis_rsm_c( &
    step, ind_p, dim_p, dim_obs_p, dim_obs, &
    dim_ens, rank_ana, state_p, ens_p, predictions, innovations, &
    U_add_obs_err, U_localize, forget, flag, info_ptr) bind(c)

    abstract INTERFACE
      SUBROUTINE add_obs_err(step, ind_p, dim_obs, HPH) BIND(C)
        ! Add observation error covariance matrix
        USE iso_c_binding
        INTEGER(c_int32_t), INTENT(in), value :: step, ind_p, dim_obs
        REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)
      END SUBROUTINE add_obs_err
      SUBROUTINE localize(step, ind_p, dim_p, dim_obs, HP_p, HPH) BIND(C)
        ! Apply localization to HP and HPH^T
        USE iso_c_binding
        INTEGER(c_int32_t), INTENT(in), value :: step, ind_p, dim_p, dim_obs
        REAL(c_double), INTENT(inout) :: HP_p(dim_obs, dim_p)
        REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)
      END SUBROUTINE localize
    end INTERFACE

    ! !ARGUMENTS:
    INTEGER(c_int32_t), INTENT(in), value :: step      ! Iteration number
    INTEGER(c_int32_t), INTENT(in), value :: ind_p      ! Integer index of PE
    INTEGER(c_int32_t), INTENT(in), value  :: dim_p     ! PE-local dimension of model state
    INTEGER(c_int32_t), INTENT(in), value :: dim_obs_p  ! PE-local dimension of observation vector
    INTEGER(c_int32_t), INTENT(in), value :: dim_obs    ! Global dimension of observation vector
    INTEGER(c_int32_t), INTENT(in), value :: dim_ens   ! Size of state ensemble
    INTEGER(c_int32_t), INTENT(in), value :: rank_ana  ! Rank to be considered for inversion of HPH
    REAL(c_double), INTENT(inout)  :: state_p(dim_p)        ! PE-local ensemble mean state
    REAL(c_double), INTENT(inout)  :: ens_p(dim_p, dim_ens) ! PE-local state ensemble
    REAL(c_double), INTENT(inout)  :: predictions(dim_obs, dim_ens) ! PE-local state ensemble
    REAL(c_double), INTENT(in)     :: innovations(dim_obs, dim_ens) ! Global array of innovations
    REAL(c_double), INTENT(in), value     :: forget    ! Forgetting factor
    INTEGER(c_int32_t), INTENT(out) :: flag    ! Status flag
    type(c_ptr), value::info_ptr

    type(c_funptr), value :: U_add_obs_err
    type(c_funptr), value :: U_localize

    type(c_function_container)::c_functions_container

    c_functions_container%info_ptr = info_ptr
    c_functions_container%localize_fptr = U_localize
    c_functions_container%add_obs_err_fptr = U_add_obs_err

    call lenkf_rsm_fortran(step, ind_p, dim_p, dim_obs_p, dim_obs, dim_ens, &
                           rank_ana, state_p, ens_p, predictions, innovations, &
                           forget, flag, c_functions_container)

  end subroutine lenkf_analysis_rsm_c

end module lenkf_rsm_c
