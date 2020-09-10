module mod_observation_manager

  use observations, ONLY: observation_set
  use localization, ONLY: base_localizer
  use assimilation_model_interface, ONLY: base_model_interface
  use forward_operator, ONLY: base_forward_operator
  use distributed_array, ONLY: darray, darray_segment
  use exceptions, ONLY: error_status, throw, new_exception

  use system_mpi

  implicit none

  type::segment_mask
     logical, allocatable::mask(:)
  end type segment_mask

  type::observation_manager

     !! Interface between the observation sets, forward operator, localizer,
     !! and assimilator. Obtains observation and prediction values for each
     !! batch.

    integer :: istep
        !! Assimilation step

    class(base_model_interface), pointer, public :: model_interface
         !! Interface to the model

    class(base_forward_operator), pointer :: forward_operator

    class(observation_set), pointer :: observation_sets(:)

    class(base_localizer), pointer :: localizer
        !! Localizer

    class(darray), pointer :: batches
        !! Batches

    type(segment_mask), allocatable :: batches_weight_masks(:,:)
        !! Mask arrays of observations with nonzero weights for each batch

    type(segment_mask), allocatable :: batches_prediction_masks(:,:)
        !! Mask arrays of observations that can be predicted for each batch

    real(kind=8) :: min_weight = 1e-10

  contains

    procedure::get_batches_obs_values
    procedure::get_batches_obs_errors
    procedure::get_batches_predictions
    procedure, private::get_batches_weight_masks
    procedure, private::get_batches_prediction_masks
    procedure, private::get_batch_weight_mask
    procedure, private::get_batch_obs_count

  end type observation_manager

contains

  function new_observation_manager( &
       istep, model_interface, batches, forward_operator, observation_sets, &
       localizer, status)

    integer::istep
        !! Assimilation step

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

    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Result
    type(observation_manager)::new_observation_manager
        !! New observation manager

    type(base_localizer), target::default_localizer

    new_observation_manager%istep = istep
    new_observation_manager%model_interface => model_interface
    new_observation_manager%batches => batches
    new_observation_manager%forward_operator => forward_operator
    new_observation_manager%observation_sets => observation_sets

    if (present(localizer)) then
      new_observation_manager%localizer => localizer
    else
      new_observation_manager%localizer => default_localizer
    end if

    new_observation_manager%batches_weight_masks = &
         new_observation_manager%get_batches_weight_masks(status)

    new_observation_manager%batches_prediction_masks = &
         new_observation_manager%get_batches_prediction_masks(status)

  end function new_observation_manager

  function get_batch_weight_mask(this, batch, obs_set, status) &
       result(mask)

    !! Get a mask array indicating whether each point in `obs_set` has a
    !! positive localization weight for at least one point in the portion of
    !! the model domain contained in `batch`. Used to determine which
    !! observations should be provided to the batch manager when processing
    !! `batch`.

    class(observation_manager), intent(inout)::this
        !! Observation manager
    class(observation_set)::obs_set
        !! Observation set
    class(darray_segment)::batch
        !! Batch segment
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Result
    logical,allocatable::mask(:)
        !! Array indicating whether each observation has a weight greater than
        !! this%min_weight for at least one point in `batch`

    integer :: iobs, imodel ! Loop counters

    real(kind=8)::w
        ! Localization weight

    allocate(mask(obs_set%get_size()))

    do iobs=1,obs_set%get_size()

       do imodel = batch%offset + 1, batch%offset + batch%length

          w=this%localizer%get_weight_model_obs( &
               this%istep, obs_set, iobs, this%model_interface, imodel, status)

          if(w>this%min_weight) then

             mask(iobs)=.true.

             cycle

          end if

       end do
    end do

  end function get_batch_weight_mask

  function get_batches_weight_masks(this, status) &
       result(batches_weight_masks)

    !! Get the weight masks for this%batches

    class(observation_manager), intent(inout)::this
        !! Observation manager
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Result
    type(segment_mask), allocatable :: batches_weight_masks(:,:)
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

    allocate(batches_weight_masks( &
         size(this%batches%segments), size(this%observation_sets)))

    do ibatch=1,size(this%batches%segments)

       batch=>this%batches%segments(ibatch)

       if(batch%rank==rank) then

          do iobs_set=1,size(this%observation_sets)

             obs_set=>this%observation_sets(iobs_set)

             batches_weight_masks(ibatch, iobs_set)%mask = &
                  this%get_batch_weight_mask(batch, obs_set)

          end do
       end if
    end do

  end function get_batches_weight_masks

  function get_batches_prediction_masks(this, status) &
       result(batches_prediction_masks)

    !! Get the weight masks for this%batches

    class(observation_manager), intent(inout)::this
        !! Observation manager
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Result
    type(segment_mask), allocatable :: batches_prediction_masks(:,:)
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

    allocate(batches_prediction_masks( &
         size(this%batches%segments), size(this%observation_sets)))

    do ibatch=1,size(this%batches%segments)

       batch=>this%batches%segments(ibatch)

       if(batch%rank==rank) then

          do iobs_set=1,size(this%observation_sets)

             obs_set=>this%observation_sets(iobs_set)

             batches_prediction_masks(ibatch, iobs_set)%mask = &
                  this%forward_operator%get_predictions_mask( &
                  this%istep, obs_set)

          end do
       end if
    end do

  end function get_batches_prediction_masks

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

    integer:: iobs_set
        ! Loop counters

    obs_count=0

    do iobs_set=1, size(this%observation_sets)
       obs_count=obs_count + count( &
            this%batches_weight_masks(ibatch,iobs_set)%mask .and. &
            this%batches_prediction_masks(ibatch,iobs_set)%mask)
    end do

  end function get_batch_obs_count

  function get_batches_obs_values(this, status) &
    result(obs_values)

    !! Get observation values for each batch

    class(observation_manager), intent(inout)::this
        !! Observation manager
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    type(darray_segment) :: obs_values(size(this%batches%segments))
        !! Observation values required for assimilation of each batch

    integer::ibatch, iobs_set, iobs, iobs_batch, imodel
        ! Loop counters

    class(observation_set), pointer::obs_set
        ! Pointer to current observation set

    class(darray_segment), pointer::batch
        ! Pointer to current batch segment

    integer::batch_obs_count ! Number of observations for a given batch
    real(kind=8), allocatable::values(:) ! Values for an observation set

    integer::rank ! MPI rank
    integer::ierr ! MPI status code

    logical :: weight_mask, prediction_mask

    call mpi_comm_rank(this%batches%comm, rank, ierr)

    do ibatch=1,size(this%batches%segments)

       batch=>this%batches%segments(ibatch)

       if(batch%rank==rank) then

          iobs_batch=1

          batch_obs_count=this%get_batch_obs_count(ibatch)
          allocate(obs_values(ibatch)%data(batch_obs_count))

          do iobs_set=1,size(this%observation_sets)

             obs_set=>this%observation_sets(iobs_set)
             values = obs_set%get_values()


             do iobs=1,this%observation_sets(iobs_set)%get_size()

                weight_mask = this%batches_weight_masks( &
                     ibatch,iobs_set)%mask(iobs)

                prediction_mask = this%batches_prediction_masks( &
                     ibatch,iobs_set)%mask(iobs)

                if(weight_mask .and. prediction_mask) then
                   obs_values(ibatch)%data(iobs_batch) = values(iobs)
                   iobs_batch = iobs_batch + 1
                end if

             end do
          end do
       end if
    end do

  end function get_batches_obs_values

  function get_batches_obs_errors(this, status) &
    result(obs_errors)

    !! Get observation errors for each batch

    class(observation_manager), intent(inout)::this
        !! Observation manager
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    type(darray_segment) :: obs_errors(size(this%batches%segments))
        !! Observation errors required for assimilation of each batch

    call throw(status, &
               new_exception('Not yet implemented', 'get_batches_obs_errors'))

  end function get_batches_obs_errors

  function get_batches_predictions(this, status) &
    result(predictions)

    !! Get predictions for each batch

    class(observation_manager), intent(inout)::this
        !! Observation manager
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    type(darray_segment) :: predictions(size(this%batches%segments))
        !! Observation errors required for assimilation of each batch

    call throw(status, &
               new_exception('Not yet implemented', 'get_batches_predictions'))

  end function get_batches_predictions

end module mod_observation_manager
