#include "mpi_types.h"

module per_member_model_interfaces
  use assimilation_model_interface, ONLY: base_model_interface
  use system_mpi
  use iso_c_binding
  use exceptions, ONLY: throw, new_exception, error_container
  use util, ONLY: str
  use distributed_array, ONLY: darray, darray_segment, new_darray

  implicit none

  type, extends(base_model_interface), abstract::per_member_model_interface
     !! I/O interface for the advect1d model

    private
    real(kind=8), allocatable::local_io_data(:, :)
    integer, allocatable::io_ranks(:)
    integer::state_size, local_io_size
    logical::state_loaded = .false.
  contains
    procedure::initialize_per_member_model_interface
    procedure::initialize => initialize_per_member_model_interface
    procedure::get_state_size
    procedure::set_state_subset
    procedure::read_state
    procedure::write_state
    procedure::get_state_darray
    procedure::is_local_member
    procedure::get_local_member_state
    procedure(I_store_member_state), deferred::store_member_state
    procedure(I_load_member_state), deferred::load_member_state
  end type per_member_model_interface

  abstract interface

    subroutine I_store_member_state( &
      this, imember, n_ensemble, member_state, state_size, status)

      !! Store the new state for one ensemble member

      use exceptions, ONLY: error_container
      import per_member_model_interface

      ! Arguments
      class(per_member_model_interface)::this
          !! Model interface
      integer, intent(in)::imember
          !! Ensemble member index
      integer, intent(in)::n_ensemble
          !! Number of ensemble members
      integer, intent(in)::state_size
          !! Size of model state
      real(kind=8), intent(in)::member_state(state_size)
          !! Model state values to write
      type(error_container), intent(out), optional::status
          !! Error status
    end subroutine I_store_member_state

    subroutine I_load_member_state( &
      this, imember, member_state, state_size, status)

      !! Read the model state from disk for a given ensemble member

      use exceptions, ONLY: error_container
      import per_member_model_interface

      ! Arguments
      class(per_member_model_interface)::this
          !! Model interface
      integer, intent(in)::imember
          !! Ensemble member index
      integer, intent(in)::state_size
          !! Size of the model state
      real(kind=8), intent(inout)::member_state(state_size)
          !! On exit, array holding the model state values
      type(error_container), intent(out), optional::status
          !! Error status

    end subroutine I_load_member_state
  end interface

contains

  subroutine initialize_per_member_model_interface( &
    this, n_ensemble, state_size, comm)

    class(per_member_model_interface), intent(inout)::this
    integer(c_int), intent(in)::n_ensemble, state_size
    MPI_COMM_TYPE, intent(in)::comm
    integer::ierr, rank, comm_size

    this%n_ensemble = n_ensemble
    this%state_size = state_size
    this%comm = comm
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

  end subroutine initialize_per_member_model_interface

  function get_state_size(this, status) result(size)
    class(per_member_model_interface)::this
    integer::size
    type(error_container), intent(out), optional::status
        !! Error status

    size = this%state_size

  end function get_state_size

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

  subroutine read_state(this, status)

    !! Read the model state and observations from disk, compute predictions

    ! Arguments
    class(per_member_model_interface)::this
        !! Model interface
    type(error_container), intent(out), optional::status
        !! Error status

    integer::imember, local_io_counter, rank, ierr

    call mpi_comm_rank(this%comm, rank, ierr)

    local_io_counter = 1

    do imember = 1, this%n_ensemble
      if (this%io_ranks(imember) == rank) then
        call this%load_member_state(imember, &
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
    class(per_member_model_interface)::this
        !! Model interface
    integer, intent(in)::imember
        !! Ensemble member index
    type(error_container), intent(out), optional :: status

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
    class(per_member_model_interface)::this
        !! Model interface
    integer, intent(in)::imember
        !! Ensemble member index
    integer, intent(in)::subset_offset
        !! Offset of subset from start of state array
    integer, intent(in)::subset_size
        !! Size of subset
    real(kind=8), intent(in)::subset_state(subset_size)
        !! Values to set
    type(error_container), intent(out), optional::status
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

  subroutine write_state(this, status)

    !! Write the new state to disk

    ! Arguments
    class(per_member_model_interface)::this
        !! Model interface
    type(error_container), intent(out), optional::status
        !! Error status

    integer::imember, rank, ierr, imember_local

    call MPI_Comm_rank(this%comm, rank, ierr)

    imember_local = 1

    do imember = 1, this%n_ensemble

      if (this%io_ranks(imember) == rank) then

        call this%store_member_state( &
          imember, this%n_ensemble, &
          this%local_io_data(:, imember_local), &
          this%state_size)

        imember_local = imember_local + 1

      end if

    end do

  end subroutine write_state

  function is_local_member(this, imember, status) result(is_local)

    ! Arguments
    class(per_member_model_interface)::this
        !! Model interface
    integer, intent(in)::imember
        !! Ensemble member index
    type(error_container), intent(out), optional::status
        !! Error status

    logical::is_local
        !! Whether imember is stored locally

    integer::rank, ierr

    call MPI_Comm_rank(this%comm, rank, ierr)

    is_local = .false.

    if (this%io_ranks(imember) == rank) is_local = .true.

  end function is_local_member

  function get_local_member_state(this, imember, status) result(member_state)

    ! Arguments
    class(per_member_model_interface), target::this
        !! Model interface
    integer, intent(in)::imember
        !! Ensemble member index
    type(error_container), intent(out), optional::status
        !! Error status

    real(kind=8), pointer::member_state(:)
        !! Whether imember is stored locally

    integer:: rank, ierr

    if (.not. this%is_local_member(imember)) then
      call throw( &
        status, new_exception( &
        'get_local_member_state called on a non-local ensemble member. ' &
        //'Call is_local_member first to check whether a member is local.'))
    end if

    call mpi_comm_rank(this%comm, rank, ierr)

    member_state => this%local_io_data( &
                    :, &
                    get_local_io_index( &
                    this%n_ensemble, this%io_ranks, rank, imember))

  end function get_local_member_state

end module per_member_model_interfaces
