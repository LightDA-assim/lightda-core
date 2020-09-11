#include "mpi_types.h"
module advect1d_observations

  use observations, ONLY: observation_set
  use exceptions, ONLY: error_status, throw, new_exception
  use util, ONLY: str
  use system_mpi

  implicit none

  type, extends(observation_set) :: advected_quantity_observation_set
     integer::istep
     integer::n_observations=0
     real(kind=8), allocatable::observations(:)
     real(kind=8), allocatable::errors(:)
     integer, allocatable::positions(:)
     MPI_COMM_TYPE::comm
  contains
    procedure::get_position
    procedure::get_size
    procedure::get_values
    procedure,private::read_observations
    procedure,private::load_observations_parallel
  end type advected_quantity_observation_set

contains

  function new_advected_quantity_observation_set(istep, comm)

    ! Arguments
    integer, intent(in)::istep
        !! Assimilation step
    MPI_COMM_TYPE, intent(in)::comm
        !! MPI communicator

    ! Result
    type(advected_quantity_observation_set):: &
         new_advected_quantity_observation_set
        !! New observation set

    new_advected_quantity_observation_set%istep = istep
    new_advected_quantity_observation_set%comm = comm
    call new_advected_quantity_observation_set%load_observations_parallel(istep)

  end function new_advected_quantity_observation_set

  subroutine read_observations(this, istep, status)

    use hdf5

    !! Read the observations from disk

    ! Arguments
    class(advected_quantity_observation_set)::this
        !! Model interface
    integer, intent(in)::istep
        !! Iteration number
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    character(len=50)::obs_filename
    integer(HID_T)::h5file_h, dset_h, dataspace
    integer(HSIZE_T)::dims(1), maxdims(1)
    integer::ierr, rank

    ! Set the HDF5 filename
    write (obs_filename, "(A,I0,A)") &
      'ensembles/', istep, '/observations.h5'

    ! Open the file
    call h5fopen_f(obs_filename, h5F_ACC_RDONLY_F, h5file_h, ierr)

    ! Open the observations dataset
    call h5dopen_f(h5file_h, 'observations', dset_h, ierr)

    ! Get the dataspace handle
    call h5dget_space_f(dset_h, dataspace, ierr)

    ! Get the dataset size
    call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, ierr)

    this%n_observations = dims(1)

    ! Close the dataspace
    call h5sclose_f(dataspace, ierr)

    allocate (this%observations(this%n_observations))

    ! Read the data
    call h5dread_f(dset_h, H5T_NATIVE_DOUBLE, this%observations, dims, ierr)

    ! Close the dataset
    call h5dclose_f(dset_h, ierr)

    ! Open the obs_positions dataset
    call h5dopen_f(h5file_h, 'obs_positions', dset_h, ierr)

    allocate (this%positions(this%n_observations))

    ! Read the data
    call h5dread_f(dset_h, H5T_NATIVE_INTEGER, this%positions, dims, ierr)

    ! Close the dataset
    call h5dclose_f(dset_h, ierr)

    ! Open the obs_errors dataset
    call h5dopen_f(h5file_h, 'obs_errors', dset_h, ierr)

    allocate (this%errors(this%n_observations))

    ! Read the data
    call h5dread_f(dset_h, H5T_NATIVE_DOUBLE, this%errors, dims, ierr)

    ! Close the dataset
    call h5dclose_f(dset_h, ierr)

    ! Close the file
    call h5fclose_f(h5file_h, ierr)

  end subroutine read_observations

  subroutine load_observations_parallel(this, istep)

    !! Load the observations from disk and broadcast them to all processors.

    ! Arguments
    class(advected_quantity_observation_set)::this
        !! Model interface
    integer, intent(in)::istep
        !! Iteration number

    integer::ierr, rank

    call mpi_comm_rank(this%comm, rank, ierr)

    if(rank==0) then
       call this%read_observations(istep)
    end if

    call mpi_bcast(this%n_observations, 1, MPI_INTEGER, 0, this%comm, ierr)

    if(rank>0) then
      deallocate (this%observations)
      allocate (this%observations(this%n_observations))
      deallocate (this%positions)
      allocate (this%positions(this%n_observations))
      deallocate (this%errors)
      allocate (this%errors(this%n_observations))
    end if

    call mpi_bcast( &
      this%observations, this%n_observations, MPI_DOUBLE_PRECISION, 0, &
      this%comm, ierr)
    call mpi_bcast( &
      this%positions, this%n_observations, MPI_INTEGER, 0, &
      this%comm, ierr)
    call mpi_bcast( &
      this%errors, this%n_observations, MPI_DOUBLE_PRECISION, 0, &
      this%comm, ierr)

  end subroutine load_observations_parallel

  function get_size(this,status) result(size)

    ! Arguments
    class(advected_quantity_observation_set)::this
        !! Observation set
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Result
    integer::size
        !! Number of observation values in set

    size = this%n_observations

  end function get_size

  function get_values(this,status) result(values)

    ! Arguments
    class(advected_quantity_observation_set)::this
        !! Observation set
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Result
    real(kind=8), allocatable::values(:)
        !! Observation values

    values = this%observations

  end function get_values

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
    integer::position
        !! Coordinate in the model domain

    character(:), allocatable::errstr

    if (iobs < 1 .or. iobs > this%n_observations) then
      errstr = 'Invalid observation index '//str(iobs)
      call throw(status, new_exception(errstr, 'get_position'))
      return
    end if

    position = this%positions(iobs)

  end function get_position

end module advect1d_observations
