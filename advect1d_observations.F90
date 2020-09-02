module advect1d_observations

  use observations, ONLY: observation_set
  use exceptions, ONLY: error_status, throw, new_exception
  use util, ONLY: str

  implicit none

  type, extends(observation_set) :: advected_quantity_observation_set
    integer::n_observations = 0
  contains
    procedure::get_position
  end type advected_quantity_observation_set

contains

  function get_position(this, istep, iobs, status) result(position)

    !! Get the position of the observation at index `iobs`.

    ! Arguments
    class(advected_quantity_observation_set)::this
        !! Observation set
    integer::istep
        !! Assimilation step
    integer::iobs
        !! Index into the observation set
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Result
    real(kind=8)::position
        !! Coordinate in the model domain

    character(:), allocatable::errstr

    if (iobs < 1 .or. iobs > this%n_observations) then
      errstr = 'Invalid observation index'//str(iobs)
      call throw(status, new_exception(errstr, 'get_position'))
      return
    end if

    call throw(status, new_exception('Not yet implemented', 'get_position'))

  end function get_position

end module advect1d_observations
