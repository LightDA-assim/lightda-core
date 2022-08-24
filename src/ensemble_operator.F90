module mod_ensemble_operator

  use, intrinsic :: iso_c_binding, only: c_double

  implicit none

  type::ensemble_operator
  contains
    procedure::execute
  end type ensemble_operator

  type, extends(ensemble_operator)::inflate_ensemble
    real(c_double)::factor = 1.0
  contains
    procedure::execute => inflate_ensemble_execute
  end type inflate_ensemble

contains

  subroutine execute( &
    this, ibatch, dim_p, dim_ens, ens_p, mgr, status)

    !! Base ensemble operator. Leaves ensemble state unchanged. Used in tests.

    use exceptions, ONLY: error_container
    use mod_base_assimilation_manager, ONLY: base_assimilation_manager

    class(ensemble_operator) :: this
    integer, intent(in)::ibatch
    integer, intent(in)::dim_p
    integer, intent(in)::dim_ens
    real(c_double), intent(inout)::ens_p(dim_p, dim_ens)

    class(base_assimilation_manager)::mgr
    type(error_container), intent(out), optional::status

  end subroutine execute

  subroutine inflate_ensemble_execute( &
    this, ibatch, dim_p, dim_ens, ens_p, mgr, status)

    !! Inflate the ensemble by a factor given by this%factor

    use exceptions, ONLY: error_container
    use mod_base_assimilation_manager, ONLY: base_assimilation_manager

    class(inflate_ensemble) :: this
    integer, intent(in)::ibatch
    integer, intent(in)::dim_p
    integer, intent(in)::dim_ens
    real(c_double), intent(inout)::ens_p(dim_p, dim_ens)

    class(base_assimilation_manager)::mgr
    type(error_container), intent(out), optional::status

    real(c_double)::row_mean, resid
    integer::irow, imember

    do irow = 1, dim_p
      row_mean = sum(ens_p(irow, :))/real(dim_ens, c_double)

      do imember = 1, dim_ens
        resid = ens_p(irow, imember) - row_mean
        ens_p(irow, imember) = row_mean + resid*this%factor
      end do
    end do

  end subroutine inflate_ensemble_execute

end module mod_ensemble_operator
