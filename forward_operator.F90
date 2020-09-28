module forward_operator

  implicit none

  type, abstract::base_forward_operator

  contains
    procedure(get_predictions_mask), deferred::get_predictions_mask
    procedure(get_predictions), deferred::get_predictions

  end type base_forward_operator

  abstract interface

    function get_predictions_mask(this, obs_set, status) result(mask)

       !! Returns a mask array indicating which observations in `obs_set`
       !! can be predicted by the model

      use observations, ONLY: observation_set
      use exceptions, ONLY: error_status

      import base_forward_operator

      ! Arguments
      class(base_forward_operator)::this
           !! Forward operator
      class(observation_set)::obs_set
           !! Observation set
      class(error_status), intent(out), allocatable, optional::status
           !! Error status

      ! Result
      logical, allocatable::mask(:)
           !! Mask array

    end function get_predictions_mask

    function get_predictions(this, obs_set, status) result(predictions)

       !! Get the predictions for the observations in `obs_set`

      use observations, ONLY: observation_set
      use exceptions, ONLY: error_status

      import base_forward_operator

      ! Arguments
      class(base_forward_operator)::this
           !! Forward operator
      class(observation_set)::obs_set
           !! Observation set
      class(error_status), intent(out), allocatable, optional::status
           !! Error status

      ! Result
      real(kind=8), allocatable::predictions(:, :)
           !! Mask array

    end function get_predictions

  end interface
end module forward_operator
