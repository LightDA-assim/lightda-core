#include "mpi_types.h"

module assimilation_model_interface

  use exceptions, ONLY: throw, new_exception, error_status
  use system_mpi

  implicit none

  type, abstract::base_model_interface

    !! Base class for model interfaces

    integer::n_ensemble
    MPI_COMM_TYPE::comm
  contains
    procedure(I_get_subset_predictions), deferred::get_subset_predictions
    procedure(I_get_subset_observations), deferred::get_subset_observations
    procedure(I_get_subset_obs_count), deferred::get_subset_obs_count
    procedure(I_get_subset_obs_err), deferred::get_subset_obs_err
    procedure::get_innovations
    procedure(I_get_state_size), deferred::get_state_size
    procedure(I_set_state_subset), deferred::set_state_subset
    procedure(I_get_state_darray), deferred::get_state_darray
    procedure::get_weight_obs_obs
    procedure::get_weight_model_obs
    procedure::read_state
    procedure::write_state
  end type base_model_interface

  abstract interface

    function I_get_state_darray(this, istep, imember) result(state_darray)

       !! Get the requested ensemble member state as a darray

      use distributed_array, ONLY: darray
      import base_model_interface

      ! Arguments
      class(base_model_interface)::this
           !! Model interface
      integer, intent(in)::istep
           !! Iteration number
      integer, intent(in)::imember
           !! Ensemble member index

      ! Returns
      type(darray)::state_darray
           !! State array represented as a darray object

    end function I_get_state_darray

    function I_get_state_size(this, istep, status) result(size)

       !! Returns the number elements in the model state array

      use exceptions, ONLY: error_status
      import base_model_interface

      implicit none

      ! Arguments
      class(base_model_interface)::this
           !! Model interface
      integer, intent(in)::istep
           !! Iteration number
      class(error_status), intent(out), allocatable, optional::status
           !! Error status

      ! Returns
      integer::size
           !! State size

    end function I_get_state_size

    subroutine I_get_subset_predictions( &
      this, istep, subset_offset, subset_size, predictions, status)
       !! Get prediction values for a subset of the model state

      use exceptions, ONLY: error_status
      import base_model_interface

      implicit none

      ! Arguments
      class(base_model_interface)::this
           !! Model interface
      integer, intent(in)::istep
           !! Iteration number
      integer, intent(in)::subset_offset
           !! Offset of subset from start of state array
      integer, intent(in)::subset_size
           !! Size of subset
      real(kind=8), intent(out)::predictions(:, :)
           !! Predicted values. Will have shape (n_observations,n_ensemble).
      class(error_status), intent(out), allocatable, optional::status
           !! Error status

    end subroutine I_get_subset_predictions

    function I_get_subset_obs_count( &
      this, istep, subset_offset, subset_size, status) result(obs_count)
       !! Get the number of observations affecting a given subset of the
       !! model domain. This will be the length of the array returned by I_get_subset_observations.

      use exceptions, ONLY: error_status
      import base_model_interface

      implicit none

      ! Arguments
      class(base_model_interface)::this
           !! Model interface
      integer, intent(in)::istep
           !! Iteration number
      integer, intent(in)::subset_offset
           !! Offset of subset from start of state array
      integer, intent(in)::subset_size
           !! Size of subset
      class(error_status), intent(out), allocatable, optional::status
           !! Error status

      integer::obs_count
           !! Number of observations

    end function I_get_subset_obs_count

    subroutine I_get_subset_observations( &
      this, istep, subset_offset, subset_size, observations, status)
       !! Get the values of observations affecting a given subset of the
       !! model domain.

      use exceptions, ONLY: error_status
      import base_model_interface

      implicit none

      ! Arguments
      class(base_model_interface)::this
           !! Model interface
      integer, intent(in)::istep
           !! Iteration number
      integer, intent(in)::subset_offset
           !! Offset of subset from start of state array
      integer, intent(in)::subset_size
           !! Size of subset
      real(kind=8), intent(out)::observations(:)
           !! Values of observations. Length must equal the value returned
           !! by get_subset_obs_count
      class(error_status), intent(out), allocatable, optional::status
           !! Error status

    end subroutine I_get_subset_observations

    subroutine I_get_subset_obs_err( &
      this, istep, subset_offset, subset_size, obs_err, status)
       !! Get the errors (uncertainties) associated with the observations
       !! affecting a given subset of the model domain.

      USE iso_c_binding
      use exceptions, ONLY: error_status
      import base_model_interface

      implicit none

      ! Arguments
      class(base_model_interface)::this
           !! Model interface
      integer, intent(in)::istep
           !! Iteration number
      integer, intent(in)::subset_offset
           !! Offset of subset from start of state array
      integer, intent(in)::subset_size
           !! Size of subset
      REAL(c_double), INTENT(out) :: obs_err(:)
           !! Values of observation errors. Length must equal the value
           !! returned by get_subset_obs_count.
      class(error_status), intent(out), allocatable, optional::status
           !! Error status

    end subroutine I_get_subset_obs_err

    subroutine I_set_state_subset( &
      this, istep, imember, subset_offset, subset_size, subset_state, status)

       !! Returns model state values for a given subset of the model state
       !! Note that ability to write data may depend on the rank of the calling
       !! processor, so this subroutine should only be called after a call to
       !! get_io_ranks, and the values of subset_offset and subset_size should
       !! be chosen to fall within a continguous segment of data stored on the
       !! calling processor rank as determined by the output from get_io_ranks.
       !! Failure to satisfy this requirement may result in an error.

      use exceptions, ONLY: error_status
      import base_model_interface

      implicit none

      ! Arguments
      class(base_model_interface)::this
           !! Model interface
      integer, intent(in)::istep
           !! Iteration number
      integer, intent(in)::imember
           !! Ensemble member index
      integer, intent(in)::subset_offset
           !! Offset of subset from start of state array
      integer, intent(in)::subset_size
           !! Size of subset
      real(kind=8), intent(in)::subset_state(subset_size)
           !! Values of the model state in the requested subset
      class(error_status), intent(out), allocatable, optional::status
           !! Error status
    end subroutine I_set_state_subset

  end interface

contains

  subroutine read_state(this, istep, status)

    !! Load the model state from disk for a given iteration `istep`. This
    !! subroutine is to be called before any calls to get_state_subset for the
    !! specified iteration `istep`. The default implementation is a no-op,
    !! enabling a model interface to load the model state incrementally as
    !! ensemble state is requested.

    ! Arguments
    class(base_model_interface)::this
        !! Model interface
    integer, intent(in)::istep
        !! Iteration number
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

  end subroutine read_state

  subroutine write_state(this, istep, status)

    !! Record the model state after assimilation of `istep`. This
    !! subroutine is to be called after all calls to set_state_subset have
    !! completed for the specified iteration `istep`. The default
    !! implementation is a no-op, enabling a model interface to load the model
    !! state incrementally as ensemble state is received from the assimilation
    !! workers.

    ! Arguments
    class(base_model_interface)::this
        !! Model interface
    integer, intent(in)::istep
        !! Iteration number
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

  end subroutine write_state

  subroutine get_innovations( &
    this, istep, subset_offset, subset_size, observations, &
    predictions, obs_errors, innovations, status)

    !! Compute innovations for a given subset of the model domain.

    use random

    implicit none

    ! Arguments
    class(base_model_interface)::this
        !! Model interface
    integer, intent(in)::istep
        !! Iteration number
    integer, intent(in)::subset_offset
        !! Offset of subset from start of state array
    integer, intent(in)::subset_size
        !! Size of subset
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

    integer::imember, iobs, obs_count

    obs_count = this%get_subset_obs_count(istep, subset_offset, subset_size)

    if (size(observations) /= obs_count) then
      print '(A,I0,A,I0)', 'Observations array has wrong length. Expected ', &
        obs_count, ', got ', size(observations)
      stop
    end if

    if (size(predictions, 1) /= obs_count .or. &
        size(predictions, 2) /= this%n_ensemble) then
      print '(A,I0,A,I0,A,I0,A,I0,A)', &
        'Predictions array has wrong shape. Expected (' &
        , obs_count, ',', this%n_ensemble, '), got (', &
        size(predictions, 1), ',', size(predictions, 2), ')'
      stop
    end if

    if (size(innovations, 1) /= obs_count .or. &
        size(innovations, 2) /= this%n_ensemble) then
      print '(A,I0,A,I0,A,I0,A,I0,A)', &
        'Innovations array has wrong shape. Expected (', &
        obs_count, ',', this%n_ensemble, '), got (' &
        , size(innovations, 1), ',', size(innovations, 2), ')'
      stop
    end if

    do imember = 1, this%n_ensemble
      do iobs = 1, obs_count
        innovations(iobs, imember) = observations(iobs) - &
                                     predictions(iobs, imember) + &
                                     random_normal()*obs_errors(iobs)
      end do
    end do

  end subroutine get_innovations

  function get_weight_obs_obs(this, istep, iobs1, iobs2, status) result(weight)

    !! Get localization weight for a given pair of observations. Default
    !! implementation returns 1 for any input (i.e., no localization).

    ! Arguments
    class(base_model_interface)::this
        !! Model interface
    integer, intent(in)::istep
        !! Iteration number
    integer, intent(in)::iobs1
        !! Index of the first observation
    integer, intent(in)::iobs2
        !! Index of the second observation
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Returns
    real(kind=8)::weight
        !! Localization weight

    weight = 1

  end function get_weight_obs_obs

  function get_weight_model_obs(this, istep, imodel, iobs, status) result(weight)

    !! Get localization weight for a given observation at a given index in the
    !! model state. Default implementation returns 1 for any input (i.e., no
    !! localization).

    ! Arguments
    class(base_model_interface)::this
        !! Model interface
    integer, intent(in)::istep
        !! Iteration number
    integer, intent(in)::imodel
        !! Index in the model state array
    integer, intent(in)::iobs
        !! Index in the observations array
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Returns
    real(kind=8)::weight
        !! Localization weight

    weight = 1

  end function get_weight_model_obs

end module assimilation_model_interface
