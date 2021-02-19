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
  use exceptions, ONLY: throw, new_exception, error_status
  use util, ONLY: str

  implicit none

  type, extends(base_assimilation_manager)::assimilation_manager
    class(assim_batch_manager), allocatable::batch_manager
    class(base_model_interface), pointer::model_interface
    class(assimilation_filter), pointer :: filter
    class(base_localizer), pointer :: localizer
        !! Localizer
    class(base_forward_operator), pointer::forward_operator
        !! Forward operator
    class(observation_set), pointer :: observation_sets(:)
        !! Observation sets

    integer::n_ensemble
    type(darray_segment), allocatable, private :: observations(:)
    type(darray_segment), allocatable, private :: obs_errors(:)
    type(darray_segment), allocatable, private :: predictions(:)
    type(observation_manager), private::obs_manager
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
    class(observation_set), target :: observation_sets(:)
        !! Observation sets
    class(base_forward_operator), target::forward_operator
        !! Forward operator
    class(base_localizer), target, optional::localizer
        !! Localizer
    integer, optional::max_batch_size
        !! Maximum batch size
    MPI_COMM_TYPE :: comm
        !! MPI communicator

    ! Result
    type(assimilation_manager)::new_assimilation_manager
        !! New assimilation manager

    integer :: actual_max_batch_size

    type(base_localizer), target::default_localizer

    if (present(max_batch_size)) then
      actual_max_batch_size = max_batch_size
    else
      actual_max_batch_size = 500
    end if

    new_assimilation_manager%model_interface => model_interface
    new_assimilation_manager%filter => filter
    new_assimilation_manager%forward_operator => forward_operator
    new_assimilation_manager%observation_sets => observation_sets
    new_assimilation_manager%batch_manager = &
      new_batch_manager( &
      model_interface, n_ensemble, &
      model_interface%get_state_size(), &
      actual_max_batch_size, model_interface%comm)

    if (present(localizer)) then
      new_assimilation_manager%localizer => localizer
    else
      new_assimilation_manager%localizer => default_localizer
    end if

  end function new_assimilation_manager

  subroutine assimilate( &
    this)

    class(assimilation_manager), intent(inout)::this
        !! Assimilation manager
    real(kind=8), allocatable::local_batches(:, :, :)
    integer, allocatable::local_batch_inds(:)
    real(kind=8), allocatable :: batch_observations(:)
    real(kind=8), allocatable :: batch_obs_err(:)
    real(kind=8), allocatable :: batch_predictions(:, :)
    real(kind=8), allocatable::batch_mean_state(:), batch_states(:, :)
    integer::rank, ierr, comm_size, n_batches, n_local_batches, ibatch, &
              ibatch_local, batch_size, state_size, n_ensemble
    type(darray)::batches
    MPI_COMM_TYPE::comm

    comm = this%batch_manager%get_comm()
    batch_size = this%batch_manager%get_batch_size()
    state_size = this%batch_manager%get_state_size()
    n_ensemble = this%batch_manager%get_n_ensemble()

    call mpi_comm_rank(comm, rank, ierr)
    call mpi_comm_size(comm, comm_size, ierr)

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
    call this%batch_manager%load_ensemble_state(local_batches)

    this%obs_manager = &
      new_observation_manager( &
      this%model_interface, batches, this%forward_operator, &
      this%observation_sets, this%localizer)

    this%observations = this%obs_manager%get_batches_obs_values()

    this%obs_errors = this%obs_manager%get_batches_obs_errors()

    this%predictions = this%obs_manager%get_batches_predictions(n_ensemble)

    ! Assimilate local batches
    do ibatch_local = 1, n_local_batches

      ibatch = local_batch_inds(ibatch_local)

      batch_observations = this%observations(ibatch)%data

      if (size(batch_observations) == 0) cycle

      batch_obs_err = this%obs_errors(ibatch)%data
      batch_predictions = reshape(this%predictions(ibatch)%data, &
                                  (/size(batch_observations), n_ensemble/))

      batch_states = local_batches(:, ibatch_local, :)

      call this%filter%assimilate( &
        ibatch, batch_size, size(batch_observations), n_ensemble, &
        batch_states, batch_predictions, &
        batch_observations, batch_obs_err, this)

      local_batches(:, ibatch_local, :) = batch_states

    end do

    ! Write the ensemble state
    call this%batch_manager%store_results(local_batches)

  end subroutine assimilate

  SUBROUTINE add_obs_err(this, ibatch, dim_obs, HPH, status)
    ! Add observation error covariance matrix
    USE iso_c_binding
    class(assimilation_manager), target::this
    INTEGER(c_int32_t), INTENT(in), value :: ibatch, dim_obs
    REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)
    class(error_status), intent(out), allocatable, optional :: status

    real(kind=8), pointer::batch_obs_err(:)
    ! Pointer to observation errors for this batch

    integer::iobs ! Loop counter

    batch_obs_err => this%obs_errors(ibatch)%data

    if (size(batch_obs_err) /= dim_obs) then
      call throw(status, new_exception( &
                 'Wrong number of observations for batch '//str(ibatch)// &
                 '. Expected '//str(size(batch_obs_err))// &
                 ', got '//str(dim_obs), &
                 'add_obs_err'))
      return
    end if

    do iobs = 1, dim_obs
      HPH(iobs, iobs) = HPH(iobs, iobs) + batch_obs_err(iobs)**2
    end do

  END SUBROUTINE add_obs_err

  SUBROUTINE localize(this, ibatch, dim_p, dim_obs, HP_p, HPH, status)
    ! Apply localization to HP and HPH^T
    USE iso_c_binding
    class(assimilation_manager), target::this
    INTEGER(c_int32_t), INTENT(in), value :: ibatch, dim_p, dim_obs
    REAL(c_double), INTENT(inout) :: HP_p(dim_obs, dim_p)
    REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)
    class(error_status), intent(out), allocatable, optional :: status

    integer, allocatable::batch_obs_set_inds(:)
       !! Indices of the originating observation set for each observation

    integer, allocatable::batch_obs_inds(:)
       !! Index of each observation in its observation set

    integer::batch_offset, batch_length
    ! Location of batch in the model state array

    class(observation_set), pointer :: obs_set1, obs_set2
    ! Observation set pointers

    integer::iobs_batch1, iobs_batch2, iobs_set1, iobs_set2, iobs_inset1, &
              iobs_inset2, ipos
    ! Loop counters

    real(kind=8)::w  ! Localization weight

    ! Locate batch in the state array
    batch_offset = this%batch_manager%get_batch_offset(ibatch)
    batch_length = this%batch_manager%get_batch_length(ibatch)

    ! Get an index into this%observation_sets for each observation
    batch_obs_set_inds = this%obs_manager%get_batch_obs_set_inds(ibatch)

    ! Get indices of the observations in their respective observation sets
    batch_obs_inds = this%obs_manager%get_batch_obs_inds(ibatch)

    if (size(batch_obs_set_inds) /= dim_obs) then
      call throw(status, new_exception( &
                 'Inconsistent observation count. Expected dim_obs='// &
                 str(size(batch_obs_set_inds))//', got dim_obs='// &
                 str(dim_obs)//'.', &
                 'localize'))
      return
    end if

    if (size(batch_obs_inds) /= dim_obs) then
      call throw(status, new_exception( &
                 'Inconsistent observation count. Expected dim_obs='// &
                 str(size(batch_obs_inds))// &
                 ', got dim_obs='//str(dim_obs)//'.', &
                 'localize'))
      return
    end if

    if (dim_p /= batch_length) then
      call throw(status, new_exception( &
                 'Inconsistent batch size for batch '//str(ibatch)// &
                 '. Expected '//str(dim_p)// &
                 ', got '//str(batch_length)//'.', &
                 'localize'))
      return
    end if

    do iobs_batch1 = 1, dim_obs

      iobs_set1 = batch_obs_set_inds(iobs_batch1)
      iobs_inset1 = batch_obs_inds(iobs_batch1)

      obs_set1 => this%observation_sets(iobs_set1)

      do iobs_batch2 = 1, dim_obs

        iobs_set2 = batch_obs_set_inds(iobs_batch2)
        iobs_inset2 = batch_obs_inds(iobs_batch2)

        obs_set2 => this%observation_sets(iobs_set2)

        ! Get localization weights
        w = this%localizer%get_weight_obs_obs( &
            obs_set1, iobs_inset1, obs_set2, iobs_inset2)

        ! Multiply HPH by the localization weights
        HPH(iobs_batch1, iobs_batch2) = HPH(iobs_batch1, iobs_batch2)*w

      end do

      do ipos = 1, dim_p

        ! Get localization weights
        w = this%localizer%get_weight_model_obs( &
            obs_set1, iobs_inset1, this%model_interface, &
            ipos + batch_offset)

        ! Multiply HP_p by the localization weights
        HP_p(iobs_batch1, ipos) = HP_p(iobs_batch1, ipos)*w

      end do

    end do

  END SUBROUTINE localize

end module mod_assimilation_manager