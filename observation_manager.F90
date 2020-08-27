module mod_observation_manager

  use observations, ONLY: observation_set
  use localization, ONLY: base_localizer
  use assimilation_model_interface, ONLY: base_model_interface
  use forward_operator, ONLY: base_forward_operator

  type::observation_manager

     !! Interface between the observation sets, forward operator, localizer,
     !! and assimilator. Obtains observation and prediction values for each
     !! batch.

    class(base_model_interface), pointer, public :: model_interface
         !! Interface to the model

    class(base_forward_operator), pointer :: forward_operator

    class(observation_set), pointer :: observation_sets(:)

    class(base_localizer), pointer :: localizer

  end type observation_manager

contains

  function new_observation_manager( &
    model_interface, forward_operator, observation_sets, localizer)

    class(base_model_interface), intent(in), target::model_interface
        !! Model interface

    class(base_forward_operator), intent(in), target::forward_operator
        !! Forward operator

    class(observation_set), intent(in), target:: observation_sets(:)
        !! Observation sets

    class(base_localizer), intent(in), target, optional :: localizer

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

end module mod_observation_manager
