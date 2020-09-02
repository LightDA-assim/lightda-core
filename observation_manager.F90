module mod_observation_manager

  use observations, ONLY: observation_set
  use localization, ONLY: base_localizer
  use assimilation_model_interface, ONLY: base_model_interface
  use forward_operator, ONLY: base_forward_operator
  use distributed_array, ONLY: darray
  use exceptions, ONLY: error_status, throw, new_exception

  implicit none

  type::observation_manager

     !! Interface between the observation sets, forward operator, localizer,
     !! and assimilator. Obtains observation and prediction values for each
     !! batch.

    class(base_model_interface), pointer, public :: model_interface
         !! Interface to the model

    class(base_forward_operator), pointer :: forward_operator

    class(observation_set), pointer :: observation_sets(:)

    class(base_localizer), pointer :: localizer

  contains

    procedure::get_batches_obs_values
    procedure::get_batches_obs_errors
    procedure::get_batches_predictions

  end type observation_manager

contains

  function new_observation_manager( &
    model_interface, forward_operator, observation_sets, localizer, status)

    class(base_model_interface), intent(in), target::model_interface
        !! Model interface

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

    new_observation_manager%model_interface => model_interface
    new_observation_manager%forward_operator => forward_operator
    new_observation_manager%observation_sets => observation_sets

    if (present(localizer)) then
      new_observation_manager%localizer => localizer
    else
      new_observation_manager%localizer => default_localizer
    end if

  end function new_observation_manager

  function get_batches_obs_values(this, istep, batches, status) &
    result(obs_values)

    !! Get observation values for each batch

    class(observation_manager), intent(inout)::this
        !! Observation manager
    integer, intent(in) :: istep
        !! Iteration number
    class(darray), intent(in) :: batches
        !! Model state array
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    type(darray) :: obs_values(size(batches%segments))
        !! Observation values required for assimilation of each batch

    call throw(status, &
               new_exception('Not yet implemented', 'get_batches_obs_values'))

  end function get_batches_obs_values

  function get_batches_obs_errors(this, istep, batches, status) &
    result(obs_errors)

    !! Get observation errors for each batch

    class(observation_manager), intent(inout)::this
        !! Observation manager
    integer, intent(in) :: istep
        !! Iteration number
    class(darray), intent(in) :: batches
        !! Model state array
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    type(darray) :: obs_errors(size(batches%segments))
        !! Observation errors required for assimilation of each batch

    call throw(status, &
               new_exception('Not yet implemented', 'get_batches_obs_errors'))

  end function get_batches_obs_errors

  function get_batches_predictions(this, istep, batches, status) &
    result(predictions)

    !! Get predictions for each batch

    class(observation_manager), intent(inout)::this
        !! Observation manager
    integer, intent(in) :: istep
        !! Iteration number
    class(darray), intent(in) :: batches
        !! Model state array
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    type(darray) :: predictions(size(batches%segments))
        !! Observation errors required for assimilation of each batch

    call throw(status, &
               new_exception('Not yet implemented', 'get_batches_predictions'))

  end function get_batches_predictions

end module mod_observation_manager
