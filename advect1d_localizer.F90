module mod_advect1d_localization

  use observations, ONLY: observation_set
  use assimilation_model_interface, ONLY: base_model_interface
  use advect1d_observations, ONLY: advected_quantity_observation_set
  use advect1d_assimilate_interfaces, ONLY: advect1d_interface
  use exceptions, ONLY: error_status, throw, new_exception
  use localization, ONLY: base_localizer, localize_gaspari_cohn
  use util, ONLY: str

  implicit none

  type, extends(base_localizer)::advect1d_localizer
    real(kind=8)::cutoff = 0.2
    real(kind=8)::cutoff_u_a = 0.2
    integer::domain_size = 100
  contains
    procedure::get_weight_obs_obs
    procedure::get_weight_model_obs
  end type advect1d_localizer

contains

  function get_weight_obs_obs(this, obs_set1, iobs1, obs_set2, iobs2, &
                              status) result(weight)

    !! Get localization weight for a given pair of observations. Default
    !! implementation returns 1 for any input (i.e., no localization).

    ! Arguments
    class(advect1d_localizer)::this
        !! Localizer
    class(observation_set), pointer::obs_set1
        !! Observation set 1
    integer, intent(in)::iobs1
        !! Index of an observation in observation set 1
    class(observation_set), pointer::obs_set2
        !! Observation set 2
    integer, intent(in)::iobs2
        !! Index of an observation in observation set 2
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Returns
    real(kind=8)::weight
        !! Localization weight

    real(kind=8)::pos1, pos2, delta, distance
    integer::domain_size

    character(:), allocatable::errstr

    weight = 1

    select type (obs_set1)
    class is (advected_quantity_observation_set)
      pos1 = real(obs_set1%get_position(iobs1))/this%domain_size
    class default
      call throw(status, new_exception('Unknown observation type', &
                                       'get_weight_obs_obs'))
      return
    end select

    select type (obs_set2)
    class is (advected_quantity_observation_set)
      pos2 = real(obs_set2%get_position(iobs2))/this%domain_size
    class default
      call throw(status, new_exception('Unknown observation type', &
                                       'get_weight_obs_obs'))
      return
    end select

    delta = abs(pos1 - pos2)
    distance = min(delta, 1 - delta)

    if (distance < -1e-8) then
      errstr = 'Invalid distance'//str(distance)// &
               'computed for obs. positions'//str(pos1)//'and '//str(pos2)
      call throw(status, new_exception(errstr, 'get_weight_obs_obs'))
      return
    end if

    distance = max(distance, 0.0)

    weight = localize_gaspari_cohn(distance, this%cutoff)

  end function get_weight_obs_obs

  function get_weight_model_obs(this, obs_set, iobs, model_interface, &
                                imodel, status) result(weight)

    !! Get localization weight for a given observation at a given index in the
    !! model state. Default implementation returns 1 for any input (i.e., no
    !! localization).

    ! Arguments
    class(advect1d_localizer)::this
        !! Model interface
    class(base_model_interface)::model_interface
        !! Model interface
    integer, intent(in)::imodel
        !! Index in the model state array
    class(observation_set), pointer::obs_set
        !! Observation set
    integer, intent(in)::iobs
        !! Index in the observation set
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Returns
    real(kind=8)::weight
        !! Localization weight

    real(kind=8)::pos_obs, pos_model, delta, distance, cutoff
    integer::domain_size

    character(:), allocatable::errstr

    select type (obs_set)
    class is (advected_quantity_observation_set)
      pos_obs = real(obs_set%get_position(iobs))/this%domain_size
    class default
      call throw(status, new_exception('Unknown observation type', &
                                       'get_weight_model_obs'))
      return
    end select

    select type (model_interface)
    class is (advect1d_interface)
      domain_size = model_interface%get_state_size()/2

      pos_model = real(mod(imodel - 1, domain_size))/domain_size
    class default
      call throw(status, new_exception('Unknown model type', &
                                       'get_weight_model_obs'))
      return
    end select

    delta = abs(pos_obs - pos_model)
    distance = min(delta, 1 - delta)

    if (distance < -1e-8) then
      errstr = 'Invalid distance'//str(distance)// &
               'computed for obs. position='//str(pos_obs)// &
               'and model position='//str(pos_model)
      call throw(status, new_exception(errstr, 'get_weight_model_obs'))
      return
    end if

    distance = max(distance, 0.0)

    if (imodel < domain_size) then
      cutoff = this%cutoff
    else
      cutoff = this%cutoff_u_a
    end if

    weight = localize_gaspari_cohn(distance, cutoff)

    if (weight > 1 .or. weight < 0) then
      errstr = 'Invalid weight='//str(weight)// &
               'returned from localize_gaspari_cohn for distance='// &
               str(distance)//'and cutoff='//str(cutoff)
      call throw(status, new_exception(errstr, 'get_weight_model_obs'))
      return
    end if

  end function get_weight_model_obs

end module mod_advect1d_localization
