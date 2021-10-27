module localization

  use observations, ONLY: observation_set
  use assimilation_model_interface, ONLY: base_model_interface
  use exceptions, ONLY: error_container, throw, new_exception

  implicit none

  type::base_localizer
    real(kind=8) :: min_weight = 1e-10
  contains
    procedure::get_weight_obs_obs
    procedure::get_weight_model_obs
    procedure::get_weight_mask
  end type base_localizer

contains

  function get_weight_obs_obs(this, obs_set1, iobs1, obs_set2, iobs2, &
                              status) result(weight)

    !! Get localization weight for a given pair of observations. Default
    !! implementation returns 1 for any input (i.e., no localization).

    ! Arguments
    class(base_localizer)::this
        !! Model interface
    class(observation_set), pointer::obs_set1
        !! Observation set 1
    integer, intent(in)::iobs1
        !! Index of an observation in observation set 1
    class(observation_set), pointer::obs_set2
        !! Observation set 2
    integer, intent(in)::iobs2
        !! Index of an observation in observation set 2
    type(error_container), intent(out), optional::status
        !! Error status

    ! Returns
    real(kind=8)::weight
        !! Localization weight

    weight = 1

  end function get_weight_obs_obs

  function get_weight_model_obs(this, obs_set, iobs, model_interface, &
                                imodel, status) result(weight)

    !! Get localization weight for a given observation at a given index in the
    !! model state. Default implementation returns 1 for any input (i.e., no
    !! localization).

    ! Arguments
    class(base_localizer)::this
        !! Model interface
    class(base_model_interface)::model_interface
        !! Model interface
    integer, intent(in)::imodel
        !! Index in the model state array
    class(observation_set)::obs_set
        !! Observation set
    integer, intent(in)::iobs
        !! Index in the observation set
    type(error_container), intent(out), optional::status
        !! Error status

    ! Returns
    real(kind=8)::weight
        !! Localization weight

    weight = 1

  end function get_weight_model_obs

  subroutine get_weight_mask(this, mask, obs_set, model_interface, &
                             imodel_start, imodel_end, status)

    !! Get a mask array whose values indicate which observations have nonzero
    !! localization weights for the given subset of the model state array

    ! Arguments
    class(base_localizer)::this
        !! Model interface
    logical::mask(:)
    class(base_model_interface)::model_interface
        !! Model interface
    integer, intent(in)::imodel_start
        !! Subset start index in the model state array
    integer, intent(in)::imodel_end
        !! Subset end index in the model state array
    class(observation_set)::obs_set
        !! Observation set
    type(error_container), intent(out), optional::status
        !! Error status

    integer::iobs
        !! Index in the observation set
    integer::imodel
        !! Index in the model state array
    real(kind=8)::w
        !! Localization weight

    if (size(mask) /= obs_set%get_size()) then
      call throw(status, new_exception('Wrong size array passed to mask ' &
                                       //'argument of get_weight_mask'))
    end if

    do iobs = 1, obs_set%get_size()

      do imodel = imodel_start, imodel_end

        w = this%get_weight_model_obs( &
            obs_set, iobs, model_interface, imodel, status)

        if (w > this%min_weight) then

          mask(iobs) = .true.

          cycle

        end if

      end do

    end do

  end subroutine get_weight_mask

  function gaspari_cohn_mid(z, c) result(f)
    real(kind=8), intent(in)::z, c
    real(kind=8)::f
    f = 1./12*(z/c)**5 - 0.5*(z/c)**4 + 5./8*(z/c)**3 &
        + 5./3*(z/c)**2 - 5*z/c - 2./3*c/z + 4
  end function gaspari_cohn_mid

  function gaspari_cohn_close(z, c) result(f)
    real(kind=8), intent(in)::z, c
    real(kind=8)::f
    f = -0.25*(z/c)**5 + 0.5*(z/c)**4 + 5./8*(z/c)**3 - 5./3*(z/c)**2 + 1
  end function gaspari_cohn_close

  function localize_gaspari_cohn(z, c, status) result(f)
    real(kind=8), intent(in)::z, c
    real(kind=8)::f
    type(error_container), intent(out), optional::status
        !! Error status

    if (z == 0) then
      f = 1
    else if (z <= c) then
      f = gaspari_cohn_close(z, c)
    else if (z <= 2*c) then
      f = max(gaspari_cohn_mid(z, c), 0.0_8)
    else if (z < 0) then
      call throw(status, new_exception( &
                 'Negative distance passed to localize_gaspari_cohn', &
                 'localize_gaspari_cohn'))
      return
    else
      f = 0
    end if

  end function localize_gaspari_cohn

end module localization
