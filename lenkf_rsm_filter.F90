module mod_lenkf_rsm_filter

  use exceptions, ONLY: error_status, throw, new_exception
  use mod_assimilation_filter, ONLY: assimilation_filter
  use lenkf_rsm, ONLY: lenkf_analysis_rsm
  use mod_base_assimilation_manager, ONLY: base_assimilation_manager

  implicit none

  type, extends(assimilation_filter)::lenkf_rsm_filter
    real(kind=8)::forget = 0.6
  contains
    procedure::assimilate
  end type lenkf_rsm_filter

contains
  subroutine get_innovations( &
    observations, predictions, obs_errors, innovations, status)

    !! Compute innovations for a given subset of the model domain.

    use random

    implicit none

    ! Arguments
    real(kind=8), intent(in)::observations(:)
        !! Observation values for the subset
    real(kind=8), intent(in)::obs_errors(:)
        !! Observation errors for the subset
    real(kind=8), intent(in)::predictions(:, :)
        !! Predictions for the subset
    real(kind=8), intent(out)::innovations(:, :)
        !! Innovations for the subset
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    integer::imember, iobs, obs_count, n_ensemble

    obs_count = size(observations)
    n_ensemble = size(predictions, 2)

    if (size(observations) /= obs_count) then
      print '(A,I0,A,I0)', 'Observations array has wrong length. Expected ', &
        obs_count, ', got ', size(observations)
      stop
    end if

    if (size(predictions, 1) /= obs_count .or. &
        size(predictions, 2) /= n_ensemble) then
      print '(A,I0,A,I0,A,I0,A,I0,A)', &
        'Predictions array has wrong shape. Expected (' &
        , obs_count, ',', n_ensemble, '), got (', &
        size(predictions, 1), ',', size(predictions, 2), ')'
      stop
    end if

    if (size(innovations, 1) /= obs_count .or. &
        size(innovations, 2) /= n_ensemble) then
      print '(A,I0,A,I0,A,I0,A,I0,A)', &
        'Innovations array has wrong shape. Expected (', &
        obs_count, ',', n_ensemble, '), got (' &
        , size(innovations, 1), ',', size(innovations, 2), ')'
      stop
    end if

    do imember = 1, n_ensemble
      do iobs = 1, obs_count
        innovations(iobs, imember) = observations(iobs) - &
                                     predictions(iobs, imember) + &
                                     random_normal()*obs_errors(iobs)
      end do
    end do

  end subroutine get_innovations

  subroutine assimilate( &
    this, istep, ibatch, dim_p, dim_obs_p, dim_obs, dim_ens, &
    ens_p, predictions, observations, obs_errors, &
    mgr, status)

    class(lenkf_rsm_filter) :: this
    integer, intent(in)::istep
    integer, intent(in)::ibatch
    integer, intent(in)::dim_p
    integer, intent(in)::dim_obs_p
    integer, intent(in)::dim_obs
    integer, intent(in)::dim_ens
    real(kind=8), intent(inout)::ens_p(dim_p, dim_ens)
    real(kind=8), intent(in)::predictions(dim_obs, dim_ens)
    real(kind=8), intent(in)::observations(dim_obs)
    real(kind=8), intent(in)::obs_errors(dim_obs)
    real(kind=8)::innovations(dim_obs, dim_ens)

    class(base_assimilation_manager)::mgr
    class(error_status), intent(out), allocatable, optional :: status

    real(kind=8)::state_p(dim_p)
    integer::flag

    call lenkf_analysis_rsm( &
      istep, ibatch, dim_p, dim_obs_p, dim_obs, &
      dim_ens, int(0), state_p, ens_p, predictions, innovations, &
      this%forget, flag, mgr)

    if (flag == 2) then
      call throw(status, new_exception( &
           'lenkf_analysis_rsm encountered a problem in solve for &
           &Kalman gain', 'assimilate'))
      return
    else if (flag /= 0) then
      call throw(status, new_exception( &
                 'Error in lenkf_analysis_rsm', 'assimilate'))
      return
    end if

  end subroutine assimilate

end module mod_lenkf_rsm_filter
