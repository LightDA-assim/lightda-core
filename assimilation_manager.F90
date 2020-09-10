#include "mpi_types.h"

module mod_assimilation_manager

  use system_mpi
  use assimilation_batch_manager, ONLY: assim_batch_manager, new_batch_manager
  use mod_observation_manager, ONLY: observation_manager, &
                                     new_observation_manager
  use observations, ONLY: observation_set
  use forward_operator, ONLY: base_forward_operator
  use localization, ONLY: base_localizer
  use assimilation_model_interface
  use mod_base_assimilation_manager, ONLY: base_assimilation_manager
  use iso_c_binding
  use distributed_array, ONLY: darray, darray_segment
  use mod_assimilation_filter, ONLY: assimilation_filter

  implicit none

  type, extends(base_assimilation_manager)::assimilation_manager
    class(assim_batch_manager), allocatable::batch_manager
    class(base_model_interface), pointer::model_interface
    class(assimilation_filter), pointer :: filter
    class(observation_manager), allocatable::obs_manager
    integer::n_ensemble
  contains
    procedure::assimilate
    procedure::localize
    procedure::add_obs_err
  end type assimilation_manager

contains

  function new_assimilation_manager( &
    model_interface, istep, n_ensemble, &
    forward_operator, observation_sets, max_batch_size, &
    localizer, filter, comm)

    !! Create a new assimilation_manager

    ! Arguments
    class(base_model_interface), target::model_interface
        !! Model interface
    integer::istep
        !! Assimilation step
    integer::n_ensemble
        !! Number of ensemble members
    class(assimilation_filter), target :: filter
        !! Assimilation filter
    class(observation_set)::observation_sets(:)
        !! Observation sets
    class(base_forward_operator), target::forward_operator
    class(base_localizer), optional::localizer
        !! Localizer
    integer, optional::max_batch_size
        !! Maximum batch size
    MPI_COMM_TYPE :: comm
        !! MPI communicator

    ! Result
    type(assimilation_manager)::new_assimilation_manager
        !! New assimilation manager

    integer :: actual_max_batch_size

    if (present(max_batch_size)) then
      actual_max_batch_size = max_batch_size
    else
      actual_max_batch_size = 500
    end if

    new_assimilation_manager%model_interface => model_interface
    new_assimilation_manager%filter => filter
    new_assimilation_manager%batch_manager = &
      new_batch_manager( &
      model_interface, n_ensemble, &
      model_interface%get_state_size(istep), &
      actual_max_batch_size, model_interface%comm)

    new_assimilation_manager%obs_manager = &
      new_observation_manager( &
      model_interface, forward_operator, observation_sets, localizer)

  end function new_assimilation_manager

  subroutine assimilate( &
    this, istep)

    class(assimilation_manager), intent(inout)::this
        !! Assimilation manager
    integer, intent(in) :: istep
        !! Iteration number
    real(kind=8), allocatable::local_batches(:, :, :)
    integer, allocatable::local_batch_inds(:)
    type(darray_segment), allocatable :: observations(:)
    type(darray_segment), allocatable :: obs_errors(:)
    type(darray_segment), allocatable :: predictions(:)
    real(kind=8), allocatable :: batch_observations(:)
    real(kind=8), allocatable :: batch_obs_err(:)
    real(kind=8), allocatable :: batch_predictions(:, :)
    real(kind=8), allocatable::batch_mean_state(:), batch_states(:, :)
    real(kind=8)::forget
    integer::rank, ierr, comm_size, n_batches, n_local_batches, ibatch, &
              n_obs_batch, ibatch_local, batch_offset, batch_length, &
              batch_size, state_size, n_ensemble
    type(darray)::batches
    MPI_COMM_TYPE::comm

    comm = this%batch_manager%get_comm()
    batch_size = this%batch_manager%get_batch_size()
    state_size = this%batch_manager%get_state_size()
    n_ensemble = this%batch_manager%get_n_ensemble()

    call mpi_comm_rank(comm, rank, ierr)
    call mpi_comm_size(comm, comm_size, ierr)

    forget = 0.8

    ! Get the number of batches and allocate batch arrays
    n_batches = this%batch_manager%get_n_batches()

    allocate (batch_mean_state(batch_size))
    allocate (batch_states(batch_size, n_ensemble))

    ! Get number of local batches
    n_local_batches = this%batch_manager%get_rank_batch_count(rank)

    ! Allocate array to hold local ensemble state
    allocate (local_batches(batch_size, n_local_batches, n_ensemble))

    ! Get local batch indices
    allocate (local_batch_inds(n_local_batches))
    call this%batch_manager%get_rank_batches(rank, local_batch_inds)

    batches = this%batch_manager%get_batches_darray()

    ! Load the ensemble state
    call this%batch_manager%load_ensemble_state(istep, local_batches)

    obs_manager = &
      new_observation_manager( &
      istep, this%model_interface, batches, this%forward_operator, this%observation_sets, this%localizer)

    observations = obs_manager%get_batches_obs_values()

    obs_errors = obs_manager%get_batches_obs_errors()

    predictions = obs_manager%get_batches_predictions()

    ! Assimilate local batches
    do ibatch_local = 1, n_local_batches
      ibatch = local_batch_inds(ibatch_local)

      batch_offset = this%batch_manager%get_batch_offset(ibatch)
      batch_length = this%batch_manager%get_batch_length(ibatch)

      batch_observations = observations(ibatch)%data
      batch_obs_err = obs_errors(ibatch)%data
      batch_predictions = reshape(predictions(ibatch)%data, &
                                  (/size(batch_observations), n_ensemble/))

      batch_states = local_batches(:, ibatch_local, :)

      call this%filter%assimilate( &
        istep, ibatch, batch_size, n_obs_batch, &
        n_obs_batch, n_ensemble, &
        batch_states, batch_predictions, &
        batch_observations, batch_obs_err, this)

      local_batches(:, ibatch_local, :) = batch_states

    end do

    ! Write the ensemble state
    call this%batch_manager%store_results(istep, local_batches)

  end subroutine assimilate

  SUBROUTINE add_obs_err(this, istep, ibatch, dim_obs, HPH)
    ! Add observation error covariance matrix
    USE iso_c_binding
    class(assimilation_manager)::this
    INTEGER(c_int32_t), INTENT(in), value :: istep, ibatch, dim_obs
    REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)
  END SUBROUTINE add_obs_err

  SUBROUTINE localize(this, istep, ibatch, dim_p, dim_obs, HP_p, HPH)
    ! Apply localization to HP and HPH^T
    USE iso_c_binding
    class(assimilation_manager)::this
    INTEGER(c_int32_t), INTENT(in), value :: istep, ibatch, dim_p, dim_obs
    REAL(c_double), INTENT(inout) :: HP_p(dim_obs, dim_p)
    REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)
  END SUBROUTINE localize

end module mod_assimilation_manager
