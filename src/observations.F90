module observations

  use exceptions, only: throw, new_exception, error_container

  implicit none

  type, abstract::observation_set
    logical, allocatable, private::mask(:)
  contains
    procedure(get_size), deferred::get_size
    procedure(get_values), deferred::get_values
    procedure(get_values), deferred::get_errors
    procedure::get_mask
    procedure::set_mask
  end type observation_set

  abstract interface

    function get_size(this, status) result(size)

      use exceptions, ONLY: error_container

      import observation_set

      ! Arguments
      class(observation_set)::this
           !! Observation set
      type(error_container), intent(out), optional::status
           !! Error status

      ! Result
      integer::size
           !! Number of observation values in set

    end function get_size

    function get_values(this, status) result(values)

       !! Get the values of the observations

      use exceptions, ONLY: error_container

      import observation_set

      ! Arguments
      class(observation_set)::this
           !! Observation set
      type(error_container), intent(out), optional::status
           !! Error status

      ! Result
      real(kind=8), allocatable::values(:)
           !! Observation values

    end function get_values

    function get_errors(this, status) result(errors)

       !! Get errors (uncertainties) associated with the observations

      use exceptions, ONLY: error_container

      import observation_set

      ! Arguments
      class(observation_set)::this
           !! Observation set
      type(error_container), intent(out), optional::status
           !! Error status

      ! Result
      real(kind=8), allocatable::errors(:)
           !! Observation errors

    end function get_errors

  end interface

contains

  function get_mask(this, status) result(mask)

    ! Arguments
    class(observation_set), target::this
         !! Observation set
    type(error_container), intent(out), optional::status
         !! Error status

    ! Result
    logical, pointer::mask(:)
         !! Number of observation values in set

    if (allocated(this%mask)) then
      ! Mask already exists, just return the pointer
      mask => this%mask
      return
    end if

    ! Allocate the mask array
    allocate (this%mask(this%get_size(status)))

    ! Set all elements to .true. by default
    this%mask = .true.

    mask => this%mask

  end function get_mask

  subroutine set_mask(this, mask, status)

    ! Arguments
    class(observation_set)::this
         !! Observation set
    logical, intent(in)::mask(:)
    type(error_container), intent(out), optional::status
         !! Error status

    if (size(mask) /= this%get_size()) then

      call throw(status, new_exception( &
                 'Size of mask differs from size of observation set.'))

    end if

    this%mask = mask

  end subroutine set_mask

end module observations
