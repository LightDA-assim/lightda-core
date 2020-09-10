module mod_observation_manager

  use observations, ONLY: observation_set
  use localization, ONLY: base_localizer
  use assimilation_model_interface, ONLY: base_model_interface
  use forward_operator, ONLY: base_forward_operator
  use distributed_array, ONLY: darray, darray_segment
  use exceptions, ONLY: error_status, throw, new_exception

  implicit none

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

  end function new_observation_manager

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
