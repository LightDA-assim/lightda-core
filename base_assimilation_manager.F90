module mod_base_assimilation_manager
  type, abstract::base_assimilation_manager
  contains
    procedure(localize), deferred :: localize
    procedure(add_obs_err), deferred :: add_obs_err
  end type base_assimilation_manager
  abstract INTERFACE
    SUBROUTINE add_obs_err(this, istep, ibatch, dim_obs, HPH, status)
      ! Add observation error covariance matrix
      USE iso_c_binding
      use exceptions, ONLY: error_status
      import base_assimilation_manager
      class(base_assimilation_manager), target::this
      INTEGER(c_int32_t), INTENT(in), value :: istep, ibatch, dim_obs
      REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)
      class(error_status), intent(out), allocatable, optional :: status
    END SUBROUTINE add_obs_err
    SUBROUTINE localize(this, istep, ibatch, dim_p, dim_obs, HP_p, HPH, status)
      ! Apply localization to HP and HPH^T
      USE iso_c_binding
      use exceptions, ONLY: error_status
      import base_assimilation_manager
      class(base_assimilation_manager), target::this
      INTEGER(c_int32_t), INTENT(in), value :: istep, ibatch, dim_p, dim_obs
      REAL(c_double), INTENT(inout) :: HP_p(dim_obs, dim_p)
      REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)
      class(error_status), intent(out), allocatable, optional :: status
    END SUBROUTINE localize

  END INTERFACE
end module mod_base_assimilation_manager
