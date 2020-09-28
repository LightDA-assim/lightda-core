#include "mpi_types.h"

module advect1d_assimilate_interfaces
  use assimilation_model_interface, ONLY: base_model_interface
  use system_mpi
  use iso_c_binding
  use exceptions, ONLY: throw, new_exception, error_status
  use hdf5
  use hdf5_exceptions, ONLY: new_hdf5_exception
  use util, ONLY: str
  use distributed_array, ONLY: darray, darray_segment, new_darray

  implicit none

  type, extends(base_model_interface)::advect1d_interface
     !! I/O interface for the advect1d model

    private
    real(kind=8), allocatable::local_io_data(:, :)
    real(kind=8)::cutoff, cutoff_u_a
    integer, allocatable::io_ranks(:)
    integer::state_size, local_io_size
    logical::observations_read = .false., predictions_computed = .false., &
              state_loaded = .false.
    integer::istep
  contains
    procedure::get_state_size
    procedure::set_state_subset
    procedure::read_state
    procedure::write_state
    procedure::get_state_darray
  end type advect1d_interface

contains

  function new_advect1d_interface( &
    istep, n_ensemble, state_size, comm, status) result(this)
    !! Create a new advect1d_interface instance

    integer(c_int), intent(in)::istep, n_ensemble, state_size
    MPI_COMM_TYPE, intent(in)::comm
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    type(advect1d_interface)::this
    integer::ierr, rank, comm_size

    this%istep = istep
    this%n_ensemble = n_ensemble
    this%state_size = state_size
    this%comm = comm

    this%cutoff = 0.1
    this%cutoff_u_a = 0.2

    this%observations_read = .false.
    this%predictions_computed = .false.
    this%state_loaded = .false.

    call mpi_comm_rank(comm, rank, ierr)
    call mpi_comm_size(comm, comm_size, ierr)

    allocate (this%io_ranks(n_ensemble))

    ! Assign ensemble members to processors for i/o purposes
    this%io_ranks = get_member_ranks(comm_size, n_ensemble)

    ! Get the number of ensemble members read and written locally
    this%local_io_size = get_rank_io_size(n_ensemble, this%io_ranks, rank)

    ! Allocate array for local i/o data
    allocate (this%local_io_data(state_size, this%local_io_size))

    call h5open_f(ierr)

#ifdef OVERRIDABLE_FINALIZERS
    call h5eset_auto_f(0, ierr)
#endif

    if (ierr < 0) then
      call throw(status, new_hdf5_exception(ierr, &
                                            'Error initializing HDF5', &
                                            procedure='new_advect1d_interface'))
      return
    end if

  end function new_advect1d_interface

  function get_member_ranks(comm_size, n_ensemble) result(io_ranks)
    !! Returns an array of rank assignments for ensemble members

    ! Arguments
    integer, intent(in)::comm_size
         !! MPI communicator size
    integer, intent(in)::n_ensemble
         !! Number of ensemble members

    integer::io_ranks(n_ensemble)
         !! Array containing the MPI rank that handles I/O for each ensemble
         !! member

    integer::i, stride

    stride = max(comm_size/n_ensemble, 1)

    do i = 1, n_ensemble
      io_ranks(i) = mod(i*stride, comm_size)
    end do

  end function get_member_ranks

  function get_rank_io_size(n_ensemble, io_ranks, rank) result(size)
    !! Returns the number of ensemble members assigned to the given
    !! rank for i/o purposes

    ! Arguments
    integer, intent(in)::n_ensemble
         !! Number of ensemble members
    integer, intent(in)::io_ranks(n_ensemble)
         !! Array of rank assignments
    integer, intent(in)::rank
         !! MPI rank

    integer::i, size

    size = 0

    do i = 1, n_ensemble
      if (io_ranks(i) == rank) size = size + 1
    end do

  end function get_rank_io_size

  function get_local_io_index(n_ensemble, io_ranks, rank, imember) result(index)
    !! Returns the index of the specified ensemble member in the local i/o
    !! data array for the given MPI rank

    ! Araguments
    integer, intent(in)::n_ensemble
         !! Number of ensemble members
    integer, intent(in)::io_ranks(n_ensemble)
         !! Array of rank assignments
    integer, intent(in)::rank
         !! MPI rank
    integer, intent(in)::imember
         !! Global index of requested ensemble member

    integer::index
         !! Local i/o index of requested ensemble member

    integer::i, local_io_counter

    index = -1

    local_io_counter = 0

    do i = 1, n_ensemble
      if (io_ranks(i) == rank) then
        local_io_counter = local_io_counter + 1
        if (i == imember) then
          index = local_io_counter
          exit
        end if
      end if
    end do

  end function get_local_io_index

  subroutine read_member_state(istep, imember, member_state, state_size, status)

    !! Read the model state from disk for a given ensemble member

    ! Arguments
    integer, intent(in)::istep
        !! Iteration number
    integer, intent(in)::imember
        !! Ensemble member index
    real(kind=8), intent(inout)::member_state(state_size)
        !! On exit, array holding the model state values
    integer, intent(in)::state_size
        !! Size of the model state
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    character(len=50)::preassim_filename
    integer(HID_T)::h5file_h, dset_h, dataspace, memspace
    integer(HSIZE_T)::offset(2), count(2), stride(2), blocksize(2)
    integer::ierr

    character(:), allocatable::errstr_suffix

    errstr_suffix = ' for member '//str(imember)

    stride = (/1, 1/)
    offset = (/imember - 1, 0/)
    count = (/1, 1/)
    blocksize = (/1, state_size/)

    ! Set the HDF5 filename
    write (preassim_filename, "(A,I0,A)") &
      'ensembles/', istep, '/preassim.h5'

    ! Open the file
    call h5fopen_f(preassim_filename, h5F_ACC_RDONLY_F, h5file_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error opening ensemble state file'//errstr_suffix, &
                 filename=preassim_filename, &
                 procedure='read_member_state'))
      return
    end if

    ! Open the dataset
    call h5dopen_f(h5file_h, 'ensemble_state', dset_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error opening dataset'//errstr_suffix, &
                 filename=preassim_filename, &
                 procedure='read_member_state'))
      return
    end if

    ! Define a dataspace within the dataset so we only read the
    ! specified ensemble member
    call h5dget_space_f(dset_h, dataspace, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error configuring dataspace'//errstr_suffix, &
                 filename=preassim_filename, &
                 procedure='read_member_state'))
      return
    end if

    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                               offset, count, ierr, stride, blocksize)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error configuring hyperslab'//errstr_suffix, &
                 filename=preassim_filename, &
                 procedure='read_member_state'))
      return
    end if

    ! Memory dataspace (needed since the local array shape differs
    ! from the dataspace in the file)
    call h5screate_simple_f(2, blocksize, memspace, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error creating memory dataspace'//errstr_suffix, &
                 filename=preassim_filename, &
                 procedure='read_member_state'))
      return
    end if

    ! Read the data
    call h5dread_f(dset_h, H5T_IEEE_F64LE, member_state, blocksize, ierr, &
                   file_space_id=dataspace, mem_space_id=memspace)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error reading member state data'//errstr_suffix, &
                 filename=preassim_filename, &
                 procedure='read_member_state'))
      return
    end if

    ! Close the dataspaces
    call h5sclose_f(dataspace, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error closing dataspace'//errstr_suffix, &
                 filename=preassim_filename, &
                 procedure='read_member_state'))
      return
    end if

    call h5sclose_f(memspace, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error closing memory dataspace'//errstr_suffix, &
                 filename=preassim_filename, &
                 procedure='read_member_state'))
      return
    end if

    ! Close the dataset
    call h5dclose_f(dset_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error closing dataset'//errstr_suffix, &
                 filename=preassim_filename, &
                 procedure='read_member_state'))
      return
    end if

    ! Close the file
    call h5fclose_f(h5file_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception( &
                 ierr, &
                 'HDF5 error closing file'//errstr_suffix, &
                 filename=preassim_filename, &
                 procedure='read_member_state'))
      return
    end if

  end subroutine read_member_state

  subroutine read_state(this, status)

    !! Read the model state and observations from disk, compute predictions

    ! Arguments
    class(advect1d_interface)::this
        !! Model interface
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    integer::imember, local_io_counter, rank, ierr

    call mpi_comm_rank(this%comm, rank, ierr)

    local_io_counter = 1

    do imember = 1, this%n_ensemble
      if (this%io_ranks(imember) == rank) then
        call read_member_state(this%istep, imember, &
                               this%local_io_data(:, local_io_counter), &
                               this%state_size)
        local_io_counter = local_io_counter + 1
      end if
    end do

    this%state_loaded = .true.

  end subroutine read_state

  function get_state_darray(this, imember, status) result(state_darray)

    !! Get the requested ensemble member state as a darray

    ! Arguments
    class(advect1d_interface)::this
        !! Model interface
    integer, intent(in)::imember
        !! Ensemble member index
    class(error_status), intent(out), allocatable, optional :: status

    type(darray)::state_darray
        !! State array represented as a darray object

    type(darray_segment)::segments(1) ! Array of segments comprising the darray

    integer::rank !! MPI rank
    integer::ierr !! MPI status code

    if (.not. this%state_loaded) then
      call this%read_state(status)
    end if

    call mpi_comm_rank(this%comm, rank, ierr)

    ! Populate the segment indicating that it covers the entire model state
    ! and is stored on the processor rank found in this%io_ranks(imember)
    segments(1)%rank = this%io_ranks(imember)
    segments(1)%comm = this%comm
    segments(1)%offset = 0
    segments(1)%length = this%state_size

    if (this%io_ranks(imember) == rank) then
      ! Copy the member state data to the darray segment
      segments(1)%data = &
        this%local_io_data(1:this%state_size, &
                           get_local_io_index(this%n_ensemble, &
                                              this%io_ranks, rank, imember))
    end if

    state_darray = new_darray(segments, this%comm)

  end function get_state_darray

  subroutine set_state_subset( &
    this, imember, subset_offset, subset_size, subset_state, status)

    !! Update the state for a portion of the model domain.

    ! Arguments
    class(advect1d_interface)::this
        !! Model interface
    integer, intent(in)::imember
        !! Ensemble member index
    integer, intent(in)::subset_offset
        !! Offset of subset from start of state array
    integer, intent(in)::subset_size
        !! Size of subset
    real(kind=8), intent(in)::subset_state(subset_size)
        !! Values to set
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    integer::rank, ierr

    call mpi_comm_rank(this%comm, rank, ierr)

    if (get_local_io_index(this%n_ensemble, this%io_ranks, rank, imember) &
        > 0) then

      ! Write to this%local_io_data
      this%local_io_data( &
        subset_offset + 1:subset_offset + subset_size, &
        get_local_io_index(this%n_ensemble, this%io_ranks, rank, imember) &
        ) = subset_state

    else
      call throw(status, new_exception( &
             'Indices for non-local ensemble state &
             &passed to set_state_subset', &
             'set_state_subset'))
      return
    end if

  end subroutine set_state_subset

  function get_state_size(this, status) result(size)

    !! Get the size of the model state array

    ! Arguments
    class(advect1d_interface)::this
        !! Model interface
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    integer::size

    size = this%state_size

  end function get_state_size

  subroutine write_state(this, status)

    !! Write the new state to disk

    ! Arguments
    class(advect1d_interface)::this
        !! Model interface
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    integer::imember, rank, ierr, imember_local

    call MPI_Comm_rank(this%comm, rank, ierr)

    imember_local = 1

    do imember = 1, this%n_ensemble

      if (this%io_ranks(imember) == rank) then

        call write_ensemble_member( &
          this, imember, this%n_ensemble, &
          this%comm, this%local_io_data(:, imember_local), &
          this%state_size)

        imember_local = imember_local + 1

      end if

    end do

  end subroutine write_state

  subroutine write_ensemble_member( &
    this, imember, n_ensemble, comm, member_state, state_size, status)

    !! Write the new state to disk for one ensemble member

    ! Arguments
    class(advect1d_interface)::this
        !! Model interface
    integer, intent(in)::imember
        !! Ensemble member index
    integer, intent(in)::n_ensemble
        !! Number of ensemble members
    MPI_COMM_TYPE::comm
        !! MPI communicator
    real(kind=8), intent(in)::member_state(state_size)
        !! Model state values to write
    integer, intent(in)::state_size
        !! Size of model state
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    character(len=80)::postassim_filename
    integer(HID_T)::h5file_h, dset_h, dataspace, memspace
    integer(HSIZE_T)::offset(2), count(2), stride(2), data_block(2)
    logical::exists
    integer::ierr

    stride = (/2, 1/)
    offset = (/imember - 1, 0/)
    count = (/1, 1/)
    data_block = (/1, state_size/)

    ! Set the HDF5 filename
    write (postassim_filename, "(A,I0,A,I0,A)") &
      'ensembles/', this%istep, '/postassim_', imember - 1, '.h5'

    ! Create the file
    call h5fcreate_f(postassim_filename, H5F_ACC_TRUNC_F, h5file_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception(ierr, &
                                            'HDF5 error creating output file', &
                                            filename=postassim_filename, &
                                            procedure='write_member_state'))
      return
    end if

    ! Define a dataspace
    call h5screate_simple_f(1, (/int(state_size, 8)/), dataspace, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception(ierr, &
                                            'HDF5 error creating dataspace', &
                                            filename=postassim_filename, &
                                            procedure='write_member_state'))
      return
    end if

    ! Create the dataset
    call h5dcreate_f(h5file_h, 'ensemble_state', H5T_IEEE_F64LE, &
                     dataspace, dset_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception(ierr, &
                                            'HDF5 error creating dataset', &
                                            filename=postassim_filename, &
                                            procedure='write_member_state'))
      return
    end if

    ! Write the data
    call h5dwrite_f( &
      dset_h, H5T_IEEE_F64LE, member_state, data_block, ierr, &
      file_space_id=dataspace)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception(ierr, &
                                            'HDF5 error writing data', &
                                            filename=postassim_filename, &
                                            procedure='write_member_state'))
      return
    end if

    ! Close the dataspace
    call h5sclose_f(dataspace, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception(ierr, &
                                            'HDF5 error closing dataspace', &
                                            filename=postassim_filename, &
                                            procedure='write_member_state'))
      return
    end if

    ! Close the dataset
    call h5dclose_f(dset_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception(ierr, &
                                            'HDF5 error closing dataset', &
                                            filename=postassim_filename, &
                                            procedure='write_member_state'))
      return
    end if

    ! Close the file
    call h5fclose_f(h5file_h, ierr)

    if (ierr < 0) then
      call throw(status, new_hdf5_exception(ierr, &
                                            'HDF5 error closing file', &
                                            filename=postassim_filename, &
                                            procedure='write_member_state'))
      return
    end if

  end subroutine write_ensemble_member

end module advect1d_assimilate_interfaces
