#include "mpi_types.h"
module advect1d_observations

  use observations, ONLY: observation_set
  use exceptions, ONLY: error_status, throw, new_exception
  use hdf5_exceptions, ONLY: new_hdf5_exception
  use util, ONLY: str
  use system_mpi

  implicit none

  type, extends(observation_set) :: advected_quantity_observation_set
    integer::istep
    integer::n_observations = 0
    real(kind=8), allocatable::observations(:)
    real(kind=8), allocatable::errors(:)
    integer, allocatable::positions(:)
    MPI_COMM_TYPE::comm
  contains
    procedure::get_position
    procedure::get_size
    procedure::get_values
    procedure::get_errors
    procedure, private::read_observations
    procedure, private::load_observations_parallel
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

    character(:), allocatable::obs_filename
    integer(HID_T)::h5file_h, dset_h, dataspace
    integer(HSIZE_T)::dims(1), maxdims(1)
    integer::ierr, rank
    logical::file_exists

    call h5open_f(ierr)

#ifdef OVERRIDABLE_FINALIZERS
    call h5eset_auto_f(0, ierr)
#endif

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'Error initializing HDF5', &
                 procedure='read_observations'))
      return
    end if

    ! Set the HDF5 filename
    obs_filename = 'ensembles/'//str(istep)//'/observations.h5'

    inquire (FILE=obs_filename, exist=file_exists)

    if (.not. file_exists) then
      call throw(status, new_exception( &
                 'File '//obs_filename//' does not exist', 'read_observations'))
      return
    end if

    ! Open the file
    call h5fopen_f('ensembles/4/observations.h5', &
                   H5F_ACC_RDONLY_F, h5file_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error opening observations file', &
                 filename=obs_filename, &
                 procedure='read_observations'))
      return
    end if

    ! Open the observations dataset
    call h5dopen_f(h5file_h, 'observations', dset_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error opening observations dataset', &
                 filename=obs_filename, procedure='read_observations'))
      return
    end if

    ! Get the dataspace handle
    call h5dget_space_f(dset_h, dataspace, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error getting dataspace handle', &
                 filename=obs_filename, procedure='read_observations'))
      return
    end if

    ! Get the dataset size
    call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error getting dataset size', &
                 filename=obs_filename, procedure='read_observations'))
      return
    end if

    this%n_observations = dims(1)

    ! Close the dataspace
    call h5sclose_f(dataspace, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error closing dataspace', &
                 filename=obs_filename, procedure='read_observations'))
      return
    end if

    allocate (this%observations(this%n_observations))

    ! Read the data
    call h5dread_f(dset_h, H5T_NATIVE_DOUBLE, this%observations, dims, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error reading observations', &
                 filename=obs_filename, procedure='read_observations'))
      return
    end if

    ! Close the dataset
    call h5dclose_f(dset_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error', &
                 filename=obs_filename, procedure='read_observations'))
      return
    end if

    ! Open the obs_positions dataset
    call h5dopen_f(h5file_h, 'obs_positions', dset_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error opening the obs_positions dataset', &
                 filename=obs_filename, procedure='read_observations'))
      return
    end if

    allocate (this%positions(this%n_observations))

    ! Read the data
    call h5dread_f(dset_h, H5T_NATIVE_INTEGER, this%positions, dims, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error reading obs positions', &
                 filename=obs_filename, procedure='read_observations'))
      return
    end if

    ! Close the dataset
    call h5dclose_f(dset_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error closing dataset', &
                 filename=obs_filename, procedure='read_observations'))
      return
    end if

    ! Open the obs_errors dataset
    call h5dopen_f(h5file_h, 'obs_errors', dset_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error', &
                 filename=obs_filename, procedure='read_observations'))
      return
    end if

    allocate (this%errors(this%n_observations))

    ! Read the data
    call h5dread_f(dset_h, H5T_NATIVE_DOUBLE, this%errors, dims, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error reading observation errors', &
                 filename=obs_filename, procedure='read_observations'))
      return
    end if

    ! Close the dataset
    call h5dclose_f(dset_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error', &
                 filename=obs_filename, procedure='read_observations'))
      return
    end if

    ! Close the file
    call h5fclose_f(h5file_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error', &
                 filename=obs_filename, procedure='read_observations'))
      return
    end if

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

    if (rank == 0) then
      call this%read_observations(istep)
    end if

    call mpi_bcast(this%n_observations, 1, MPI_INTEGER, 0, this%comm, ierr)

    if (rank > 0) then
      allocate (this%observations(this%n_observations))
      allocate (this%positions(this%n_observations))
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

  function get_size(this, status) result(size)

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

  function get_values(this, status) result(values)

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

  function get_errors(this, status) result(errors)

    ! Arguments
    class(advected_quantity_observation_set)::this
        !! Observation set
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Result
    real(kind=8), allocatable::errors(:)
        !! Observation errors

    errors = this%errors

  end function get_errors

  function get_position(this, iobs, status) result(position)

    !! Get the position of the observation at index `iobs`.

    ! Arguments
    class(advected_quantity_observation_set)::this
        !! Observation set
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
