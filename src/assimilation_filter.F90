module mod_assimilation_filter

  implicit none

  type::assimilation_filter
  contains
    procedure::assimilate
  end type assimilation_filter

contains

  subroutine assimilate( &
    this, ibatch, dim_p, dim_obs, dim_ens, &
    ens_p, predictions, observations, obs_errors, &
    mgr, status)

    !! Base assimilation filter. Leaves ensemble state unchanged. Used in tests.

    use exceptions, ONLY: error_container
    use mod_base_assimilation_manager, ONLY: base_assimilation_manager

    class(assimilation_filter) :: this
    integer, intent(in)::ibatch
    integer, intent(in)::dim_p
    integer, intent(in)::dim_obs
    integer, intent(in)::dim_ens
    real(kind=8), intent(inout)::ens_p(dim_p, dim_ens)
    real(kind=8), intent(in)::predictions(dim_obs, dim_ens)
    real(kind=8), intent(in)::observations(dim_obs)
    real(kind=8), intent(in)::obs_errors(dim_obs)

    class(base_assimilation_manager)::mgr
    class(error_container), intent(out), allocatable, optional::status

  end subroutine assimilate

end module mod_assimilation_filter
