module mod_observation_manager

  use observations, ONLY: observation_set
  use localization, ONLY: base_localizer
  use assimilation_model_interface, ONLY: base_model_interface
  use forward_operator, ONLY: base_forward_operator
  use distributed_array, ONLY: darray, darray_segment
  use exceptions, ONLY: error_container, throw, new_exception
  use util, ONLY: str
  use fhash, ONLY: fhash_tbl_t, key => fhash_key

  use system_mpi

  implicit none

  type::segment_mask
    logical, allocatable::mask(:)
  end type segment_mask

  type::observation_manager

     !! Interface between the observation sets, forward operator, localizer,
     !! and assimilator. Obtains observation and prediction values for each
     !! batch.

    class(base_model_interface), pointer, public :: model_interface
        !! Interface to the model

    class(base_forward_operator), pointer :: forward_operator
        !! Forward operator

    class(observation_set), pointer :: observation_sets(:)
        !! Observation sets

    class(base_localizer), pointer :: localizer
        !! Localizer

    class(darray), pointer :: batches
        !! Batches

    type(fhash_tbl_t), allocatable :: batches_weight_masks(:)
        !! Mask arrays of observations with nonzero weights for each batch

    type(segment_mask), allocatable :: prediction_masks(:)
        !! Mask arrays of observations that can be predicted for each batch

    real(kind=8) :: min_weight = 1e-10

  contains

    procedure::get_batches_obs_values
    procedure::get_batches_obs_errors
    procedure::get_batches_predictions
    procedure::get_batch_obs_inds
    procedure::get_batch_obs_set_inds
    procedure, private::get_batches_weight_masks
    procedure, private::get_prediction_mask
    procedure, private::get_batch_weight_mask
    procedure, private::get_batch_obs_count

  end type observation_manager

contains
  !> Define custom getter for logical arrays
  subroutine fhash_get_segment_mask_ptr(tbl, k, ptr)

    use fhash, only: fhash_key_t
    type(fhash_tbl_t), intent(in) :: tbl
    class(fhash_key_t), intent(in) :: k
    type(segment_mask), intent(out), pointer :: ptr

    integer :: stat
    class(*), pointer :: data

    call tbl%get_raw_ptr(k, data, stat)

    if (stat /= 0) print *, 'error ', stat! Error handling: key not found

    select type (d => data)
    type is (segment_mask)
      ptr => d
    class default
      ! Error handling: found wrong type
      print *, 'fhash_get_logical_ptr called on a hash table value ' &
        //'that is not a logical array.'
      error stop
    end select

  end subroutine fhash_get_segment_mask_ptr

  function new_observation_manager( &
    model_interface, batches, forward_operator, observation_sets, &
    localizer, status)

    class(base_model_interface), intent(in), target::model_interface
        !! Model interface

    class(darray), target::batches
        !! Assimilation batches

    class(base_forward_operator), intent(in), target::forward_operator
        !! Forward operator

    class(observation_set), intent(in), target:: observation_sets(:)
        !! Observation sets

    class(base_localizer), intent(in), target, optional :: localizer
        !! Localizer

    type(error_container), intent(out), optional::status
        !! Error status

    ! Result
    type(observation_manager)::new_observation_manager
        !! New observation manager

    type(base_localizer), target::default_localizer

    new_observation_manager%model_interface => model_interface
    new_observation_manager%batches => batches
    new_observation_manager%forward_operator => forward_operator
    new_observation_manager%observation_sets => observation_sets

    if (present(localizer)) then
      new_observation_manager%localizer => localizer
    else
      new_observation_manager%localizer => default_localizer
    end if

    allocate (new_observation_manager%batches_weight_masks( &
              size(observation_sets)))
    allocate (new_observation_manager%prediction_masks(size(observation_sets)))

  end function new_observation_manager

  function get_batch_weight_mask(this, ibatch, iobs_set, status) &
    result(mask)

    !! Get a mask array indicating whether each point in `obs_set` has a
    !! positive localization weight for at least one point in the portion of
    !! the model domain contained in `batch`. Used to determine which
    !! observations should be provided to the batch manager when processing
    !! `batch`.

    class(observation_manager), intent(inout)::this
        !! Observation manager
    integer::ibatch
        !! Batch segment index
    integer::iobs_set
        !! Observation set index
    type(error_container), intent(out), optional::status
        !! Error status

    ! Result
    logical, pointer::mask(:)
        !! Array indicating whether each observation has a weight greater than
        !! this%min_weight for at least one point in `batch`

    integer :: iobs, imodel ! Loop counters

    type(segment_mask), pointer::mask_container

    class(observation_set), pointer::obs_set
        !! Observation set
    class(darray_segment), pointer::batch
        !! Batch segment

    real(kind=8)::w
    ! Localization weight

    integer::stat ! Status of key in table
    integer::nobs ! Number of observations

    call this%batches_weight_masks(iobs_set)%check_key(key(ibatch), stat)

    if (stat == 0) then
      ! Assign pointer and return
      call fhash_get_segment_mask_ptr( &
        this%batches_weight_masks(iobs_set), key(ibatch), mask_container)
      mask => mask_container%mask
      return
    end if

    ! Assign batch pointer
    batch => this%batches%segments(ibatch)

    ! Assign observation set pointer
    obs_set => this%observation_sets(iobs_set)

    nobs = obs_set%get_size()

    ! Allocate the array
    allocate (mask_container)
    allocate (mask_container%mask(nobs))
    mask => mask_container%mask

    ! Set mask to false initially
    mask = .false.

    do iobs = 1, nobs

      do imodel = batch%offset + 1, batch%offset + batch%length

        w = this%localizer%get_weight_model_obs( &
            obs_set, iobs, this%model_interface, imodel, status)

        if (w > this%min_weight) then

          mask(iobs) = .true.

          cycle

        end if

      end do

    end do

    call this%batches_weight_masks(iobs_set)%set( &
      key(ibatch), value=mask_container)

  end function get_batch_weight_mask

  function get_batches_weight_masks(this, status) &
    result(batches_weight_masks)

    !! Get the weight masks for this%batches

    class(observation_manager), intent(inout)::this
        !! Observation manager
    type(error_container), intent(out), optional::status
        !! Error status

    ! Result
    type(segment_mask), allocatable :: batches_weight_masks(:, :)
        !! Mask arrays of observations with nonzero weights for each batch

    integer::ibatch, iobs_set, iobs, iobs_batch, imodel
    ! Loop counters

    class(observation_set), pointer::obs_set
    ! Pointer to current observation set

    class(darray_segment), pointer::batch
    ! Pointer to current batch segment

    integer::rank ! MPI rank
    integer::ierr ! MPI status code

    call mpi_comm_rank(this%batches%comm, rank, ierr)

    allocate (batches_weight_masks( &
              size(this%batches%segments), size(this%observation_sets)))

    do ibatch = 1, size(this%batches%segments)

      batch => this%batches%segments(ibatch)

      if (batch%rank == rank) then

        do iobs_set = 1, size(this%observation_sets)

          obs_set => this%observation_sets(iobs_set)

          batches_weight_masks(ibatch, iobs_set)%mask = &
            this%get_batch_weight_mask(ibatch, iobs_set)

        end do
      end if
    end do

  end function get_batches_weight_masks

  function get_prediction_mask(this, iobs_set, status) &
    result(mask)

    !! Get the weight masks for this%batches

    class(observation_manager), target, intent(inout)::this
        !! Observation manager
    type(error_container), intent(out), optional::status
        !! Error status

    ! Result
    type(segment_mask), allocatable :: batches_prediction_masks(:, :)
        !! Mask arrays of observations with nonzero weights for each batch

    integer::ibatch, iobs_set, iobs, iobs_batch, imodel
    ! Loop counters

    class(observation_set), pointer::obs_set
    ! Pointer to current observation set

    class(darray_segment), pointer::batch
    ! Pointer to current batch segment

    ! Result
    logical, pointer::mask(:)

    integer::nobs ! Number of observations

    if (allocated(this%prediction_masks(iobs_set)%mask)) then
      ! Assign pointer and return
      mask => this%prediction_masks(iobs_set)%mask
      return
    end if

    obs_set => this%observation_sets(iobs_set)

    nobs = obs_set%get_size()

    ! Allocate the array
    allocate (this%prediction_masks(iobs_set)%mask(nobs))
    mask => this%prediction_masks(iobs_set)%mask

    ! Get prediction mask from forward operator
    mask = this%forward_operator%get_predictions_mask(obs_set)

  end function get_prediction_mask

  function get_batch_obs_count(this, ibatch) result(obs_count)

    !! Get the number of observations to be assimilated by the batch at index
    !! `ibatch`

    ! Arguments
    class(observation_manager), intent(inout)::this
        !! Observation manager
    integer, intent(in) :: ibatch
        !! Batch index

    ! Result
    integer::obs_count
        !! Number of observations associated with `ibatch`

    integer:: iobs_set ! Loop counter

    logical, pointer::weight_mask(:), prediction_mask(:)

    obs_count = 0

    do iobs_set = 1, size(this%observation_sets)
      obs_count = obs_count + count( &
                  this%get_batch_weight_mask(ibatch, iobs_set) .and. &
                  this%get_prediction_mask(iobs_set))
    end do

  end function get_batch_obs_count

  function get_batch_obs_inds(this, ibatch) result(obs_inds)

    !! Get the indices in the respective observation set for each observation
    !! to be assimilated into the batch `ibatch`.
    !!
    !! This function returns the indices to locate observations in the
    !! observation sets; use get_batch_obs_set_inds to determine which
    !! observation set each observation comes from

    ! Arguments
    class(observation_manager), intent(inout)::this
        !! Observation manager
    integer, intent(in) :: ibatch
        !! Batch index

    ! Result
    integer, allocatable::obs_inds(:)
        !! Number of observations associated with `ibatch`

    integer:: iobs_set, iobs_batch, iobs_inset
    ! Loop counters

    logical, allocatable::mask(:)

    iobs_batch = 1

    allocate (obs_inds(this%get_batch_obs_count(ibatch)))

    do iobs_set = 1, size(this%observation_sets)

      mask = (this%get_batch_weight_mask(ibatch, iobs_set) &
              .and. this%get_prediction_mask(iobs_set))

      do iobs_inset = 1, this%observation_sets(iobs_set)%get_size()

        if (mask(iobs_inset)) then

          obs_inds(iobs_batch) = iobs_inset
          iobs_batch = iobs_batch + 1

        end if

      end do
    end do

  end function get_batch_obs_inds

  function get_batch_obs_set_inds(this, ibatch) result(set_inds)

    !! Get the indices of the respective observation set for each observation
    !! to be assimilated into the batch `ibatch`
    !!
    !! This function returns the index of the observation set for each
    !! observation; use the indices returned from get_batch_obs_inds to locate
    !! the observations in their respective observation sets.

    ! Arguments
    class(observation_manager), intent(inout)::this
        !! Observation manager
    integer, intent(in) :: ibatch
        !! Batch index

    ! Result
    integer, allocatable::set_inds(:)
        !! Number of observations associated with `ibatch`

    integer:: iobs_set, iobs_batch
    ! Loop counters

    integer::set_batch_overlap
    ! Number of observations in a set that influence a given batch

    iobs_batch = 1

    allocate (set_inds(this%get_batch_obs_count(ibatch)))

    do iobs_set = 1, size(this%observation_sets)

      ! Determine how many observations from the set are used in this batch
      set_batch_overlap = &
        count( &
        this%get_batch_weight_mask(ibatch, iobs_set) .and. &
        this%get_prediction_mask(iobs_set))

      ! Store the set index in set_inds
      set_inds(iobs_batch:iobs_batch + set_batch_overlap - 1) = iobs_set

      ! Advance iobs_batch
      iobs_batch = iobs_batch + set_batch_overlap

    end do

  end function get_batch_obs_set_inds

  function get_batches_obs_values(this, status) &
    result(obs_values)

    !! Get observation values for each batch

    class(observation_manager), intent(inout)::this
        !! Observation manager
    type(error_container), intent(out), optional::status
        !! Error status

    type(darray_segment) :: obs_values(size(this%batches%segments))
        !! Observation values required for assimilation of each batch

    integer::ibatch, iobs_set, iobs, iobs_batch
    ! Loop counters

    class(observation_set), pointer::obs_set
    ! Pointer to current observation set

    class(darray_segment), pointer::batch
    ! Pointer to current batch segment

    integer::batch_obs_count ! Number of observations for a given batch
    real(kind=8), allocatable::values(:) ! Values for an observation set

    integer::rank ! MPI rank
    integer::ierr ! MPI status code
    integer::n_obs ! Number of observations in set
    integer::n_sets ! Number of observation sets

    logical, allocatable :: mask(:)

    call mpi_comm_rank(this%batches%comm, rank, ierr)

    n_sets = size(this%observation_sets)

    do iobs_set = 1, size(this%observation_sets)

      if (rank == 0) then
        print *, 'Loading values from observation set '//str(iobs_set)// &
          ' of '//str(n_sets)
      end if

      obs_set => this%observation_sets(iobs_set)
      values = obs_set%get_values()

      if (rank == 0) then
        n_obs = size(values)
        print *, 'Read '//str(n_obs)//' observations from observation set '// &
          str(iobs_set)
      end if

      do ibatch = 1, size(this%batches%segments)

        batch => this%batches%segments(ibatch)

        if (batch%rank == rank) then

          iobs_batch = 1

          batch_obs_count = this%get_batch_obs_count(ibatch)
          allocate (obs_values(ibatch)%data(batch_obs_count))

          do iobs = 1, this%observation_sets(iobs_set)%get_size()

            mask = this%get_batch_weight_mask(ibatch, iobs_set) .and. &
                   this%get_prediction_mask(iobs_set)

            if (any(mask)) then
              obs_values(ibatch)%data(iobs_batch) = values(iobs)
              iobs_batch = iobs_batch + 1
            end if

          end do
        end if
      end do
    end do

  end function get_batches_obs_values

  function get_batches_obs_errors(this, status) &
    result(obs_errors)

    !! Get observation errors for each batch

    ! Arguments
    class(observation_manager), intent(inout)::this
        !! Observation manager
    type(error_container), intent(out), optional::status
        !! Error status

    ! Result
    type(darray_segment) :: obs_errors(size(this%batches%segments))
        !! Observation errors required for assimilation of each batch

    integer::ibatch, iobs_set, iobs, iobs_batch
    ! Loop counters

    class(observation_set), pointer::obs_set
    ! Pointer to current observation set

    class(darray_segment), pointer::batch
    ! Pointer to current batch segment

    integer::batch_obs_count ! Number of observations for a given batch
    real(kind=8), allocatable::errors(:) ! Values for an observation set

    integer::rank ! MPI rank
    integer::ierr ! MPI status code

    logical, allocatable :: mask(:)

    call mpi_comm_rank(this%batches%comm, rank, ierr)

    do ibatch = 1, size(this%batches%segments)

      batch => this%batches%segments(ibatch)

      if (batch%rank == rank) then

        iobs_batch = 1

        batch_obs_count = this%get_batch_obs_count(ibatch)
        allocate (obs_errors(ibatch)%data(batch_obs_count))

        do iobs_set = 1, size(this%observation_sets)

          obs_set => this%observation_sets(iobs_set)
          errors = obs_set%get_errors()

          do iobs = 1, this%observation_sets(iobs_set)%get_size()

            mask = this%get_batch_weight_mask(ibatch, iobs_set) .and. &
                   this%get_prediction_mask(iobs_set)

            if (any(mask)) then
              obs_errors(ibatch)%data(iobs_batch) = errors(iobs)
              iobs_batch = iobs_batch + 1
            end if

          end do
        end do
      end if
    end do

  end function get_batches_obs_errors

  function get_batches_predictions(this, n_ensemble, status) &
    result(predictions)

    !! Get predictions for each batch

    class(observation_manager), intent(inout)::this
        !! Observation manager
    integer, intent(in) :: n_ensemble
    type(error_container), intent(out), optional::status
        !! Error status

    type(darray_segment) :: predictions(size(this%batches%segments)*n_ensemble)
        !! Observation errors required for assimilation of each batch

    integer::ibatch, iobs_set, iobs, iobs_batch
    ! Loop counters

    class(observation_set), pointer::obs_set
    ! Pointer to current observation set

    class(darray_segment), pointer::batch
    ! Pointer to current batch segment

    integer::batch_obs_count ! Number of observations for a given batch

    real(kind=8), allocatable::set_predictions(:, :)
    ! Predicted values for an observation set

    real(kind=8), allocatable::batch_predictions(:, :)
    ! Predicted values for a batch

    integer::rank ! MPI rank
    integer::ierr ! MPI status code

    logical, allocatable:: mask(:)

    call mpi_comm_rank(this%batches%comm, rank, ierr)

    do iobs_set = 1, size(this%observation_sets)

      obs_set => this%observation_sets(iobs_set)
      set_predictions = this%forward_operator%get_predictions(obs_set)

      do ibatch = 1, size(this%batches%segments)

        batch => this%batches%segments(ibatch)

        if (batch%rank == rank) then

          iobs_batch = 1

          batch_obs_count = this%get_batch_obs_count(ibatch)

          allocate (batch_predictions(batch_obs_count, n_ensemble))

          do iobs = 1, this%observation_sets(iobs_set)%get_size()

            mask = this%get_batch_weight_mask(ibatch, iobs_set) .and. &
                   this%get_prediction_mask(iobs_set)

            if (any(mask)) then
              batch_predictions(iobs_batch, :) = set_predictions(iobs, :)
              iobs_batch = iobs_batch + 1
            end if

          end do

          allocate (predictions(ibatch)%data, &
                    source=pack(batch_predictions, .true.))

          deallocate (batch_predictions)

        end if

      end do
    end do

  end function get_batches_predictions

end module mod_observation_manager
