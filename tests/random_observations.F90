#include "mpi_types.h"
module random_observations

  use observations, ONLY: observation_set
  use exceptions, ONLY: error_status, throw, new_exception
  use hdf5_exceptions, ONLY: new_hdf5_exception
  use util, ONLY: str
  use system_mpi

  implicit none

  type, extends(observation_set) :: random_observation_set
    private
    integer::size = 0
    real(kind=8), allocatable::observations(:)
    real(kind=8), allocatable::positions(:)
    real(kind=8), allocatable::errors(:)
    MPI_COMM_TYPE::comm
  contains
    procedure::get_position
    procedure::get_size
    procedure::get_values
    procedure::get_errors
    procedure, private::generate
  end type random_observation_set

contains

  function new_random_observation_set(size, comm)

    ! Arguments
    integer, intent(in)::size
        !! Number of observations in set
    MPI_COMM_TYPE, intent(in)::comm
        !! MPI communicator

    ! Result
    type(random_observation_set):: &
      new_random_observation_set
        !! New observation set

    new_random_observation_set%size = size
    new_random_observation_set%comm = comm
    call new_random_observation_set%generate()

  end function new_random_observation_set

  subroutine generate(this)
    class(random_observation_set)::this

    allocate (this%observations(this%size))
    allocate (this%positions(this%size))
    allocate (this%errors(this%size))

    call random_number(this%observations)
    call random_number(this%positions)
    call random_number(this%errors)

  end subroutine generate

  function get_size(this, status) result(size)

    ! Arguments
    class(random_observation_set)::this
        !! Observation set
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Result
    integer::size
        !! Number of observation values in set

    size = this%size

  end function get_size

  function get_values(this, status) result(values)

    ! Arguments
    class(random_observation_set)::this
        !! Observation set
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Result
    real(kind=8), allocatable::values(:)
        !! Observation values

    values = this%observations

  end function get_values

  function get_errors(this, status) result(errors)

    ! Arguments
    class(random_observation_set)::this
        !! Observation set
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Result
    real(kind=8), allocatable::errors(:)
        !! Observation errors

    errors = this%errors

  end function get_errors

  function get_position(this, iobs, status) result(pos)

    ! Arguments
    class(random_observation_set)::this
        !! Observation set
    integer, intent(in)::iobs
        !! Observation index
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Result
    real(kind=8)::pos
        !! Observation errors

    pos = this%positions(iobs)

  end function get_position

end module random_observations
