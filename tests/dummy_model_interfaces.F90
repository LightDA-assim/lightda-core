#include "mpi_types.h"

module dummy_model_interfaces
  use assimilation_model_interface, ONLY: base_model_interface
  use system_mpi
  use iso_c_binding
  use random, ONLY: random_normal
  use exceptions, ONLY: throw, new_exception, error_container
  use distributed_array, ONLY: darray, darray_segment, new_darray

  implicit none

  type, extends(base_model_interface)::dummy_model_interface
    private
    real(kind=8), allocatable::observations(:), obs_errors(:), predictions(:, :)
    real(kind=8), pointer::local_io_data(:, :)
    real(kind=8)::cutoff, cutoff_u_a
    integer, allocatable::obs_positions(:), io_ranks(:)
    integer::n_observations, state_size, local_io_size
    logical::observations_read, predictions_computed, state_loaded
  contains
    procedure::get_state_size
    procedure::set_state_subset
    procedure::read_state
    procedure::get_ensemble_state
    procedure::get_state_darray
  end type dummy_model_interface

contains

  function new_dummy_model( &
    n_ensemble, n_observations, state_size, comm) result(this)

    integer(c_int), intent(in)::n_ensemble, n_observations, state_size
    MPI_COMM_TYPE::comm
    type(dummy_model_interface)::this
    integer::ierr, rank, comm_size

    this%n_ensemble = n_ensemble
    this%n_observations = n_observations
    this%state_size = state_size
    this%comm = comm

    this%observations_read = .false.
    this%predictions_computed = .false.
    this%state_loaded = .false.

    call mpi_comm_rank(comm, rank, ierr)
    call mpi_comm_size(comm, comm_size, ierr)

    allocate (this%observations(this%n_observations))
    allocate (this%obs_errors(this%n_observations))
    allocate (this%obs_positions(this%n_observations))
    allocate (this%predictions(this%n_observations, this%n_ensemble))
    allocate (this%io_ranks(n_ensemble))

    this%io_ranks = get_member_ranks(comm_size, n_ensemble)

    ! Get the number of ensemble members read and written locally
    this%local_io_size = get_rank_io_size(n_ensemble, this%io_ranks, rank)

    ! Allocate array for local i/o data
    allocate (this%local_io_data(state_size, this%n_ensemble))

  end function new_dummy_model

  function get_state_size(this, status) result(size)
    class(dummy_model_interface)::this
    integer::size
    type(error_container), intent(out), optional::status
        !! Error status

    size = this%state_size

  end function get_state_size

  function get_member_ranks(comm_size, n_ensemble) result(io_ranks)
    integer, intent(in)::comm_size, n_ensemble
    integer::io_ranks(n_ensemble)
    integer::i, stride
    stride = max(comm_size/n_ensemble, 1)

    do i = 1, n_ensemble
      io_ranks(i) = mod(i*stride, comm_size)
    end do

  end function get_member_ranks

  function get_rank_io_size(n_ensemble, io_ranks, rank) result(size)
    integer, intent(in)::n_ensemble
    integer, intent(in)::io_ranks(n_ensemble)
    integer, intent(in)::rank
    integer::i, size

    size = 0

    do i = 1, n_ensemble
      if (io_ranks(i) == rank) size = size + 1
    end do

  end function get_rank_io_size

  subroutine read_state(this, status)
    class(dummy_model_interface)::this
    type(error_container), intent(out), optional::status
        !! Error status

    integer::imember, rank, ierr, i
    real(kind=8)::r

    call mpi_comm_rank(this%comm, rank, ierr)

    ! Populate local i/o data
    do imember = 1, this%n_ensemble
      do i = 1, this%state_size
        this%local_io_data(i, imember) = imember*this%state_size + i
      end do
    end do

    ! Broadcast so all processes have the same ensemble state
    call MPI_Bcast( &
      this%local_io_data, this%state_size*this%n_ensemble, &
      MPI_DOUBLE_PRECISION, 0, this%comm, ierr)

    this%state_loaded = .true.

  end subroutine read_state

  function get_state_darray(this, imember, status) result(state_darray)

    !! Get the requested ensemble member state as a darray

    ! Arguments
    class(dummy_model_interface)::this
        !! Model interface
    integer, intent(in)::imember
        !! Ensemble member index
    type(error_container), intent(out), optional::status
        !! Error status

    type(darray)::state_darray
        !! State array represented as a darray object

    type(darray_segment)::segments(1)

    integer::rank !! MPI rank
    integer::ierr !! MPI status code

    if (.not. this%state_loaded) call this%read_state()

    call mpi_comm_rank(this%comm, rank, ierr)

    ! Populate the segment indicating that it covers the entire model state
    ! and is stored on the processor rank found in this%io_ranks(imember)
    segments(1)%rank = this%io_ranks(imember)
    segments(1)%comm = this%comm
    segments(1)%offset = 0
    segments(1)%length = this%state_size

    if (this%io_ranks(imember) == rank) then
      ! Copy the member state data to the darray segment
      segments(1)%data = this%local_io_data(1:this%state_size, imember)
    end if

    state_darray = new_darray(segments, this%comm)

  end function get_state_darray

  subroutine set_state_subset( &
    this, imember, subset_offset, subset_size, subset_state, status)

    class(dummy_model_interface)::this
    integer, intent(in)::imember, subset_offset, subset_size
    real(kind=8), intent(in)::subset_state(subset_size)
    type(error_container), intent(out), optional::status
        !! Error status

    integer::rank, ierr

    call mpi_comm_rank(this%comm, rank, ierr)

    if (this%io_ranks(imember) == rank) then

      this%local_io_data( &
        subset_offset + 1:subset_offset + subset_size, &
        imember) = subset_state

    else
      call throw(status, new_exception( &
              'Indices for non-local ensemble &
              &state passed to set_state_subset'))
      return
    end if

  end subroutine set_state_subset

  subroutine after_ensemble_results_received(this, istep)
    class(dummy_model_interface)::this
    integer, intent(in)::istep
    integer::imember, rank, ierr, imember_local

    call MPI_Comm_rank(this%comm, rank, ierr)

    imember_local = 1

    do imember = 1, this%n_ensemble

      if (this%io_ranks(imember) == rank) then

        imember_local = imember_local + 1

      end if

    end do

  end subroutine after_ensemble_results_received

  function get_ensemble_state(this) result(local_io_data)
    class(dummy_model_interface)::this
    real(kind=8)::local_io_data(this%state_size, this%n_ensemble)
    integer::imember, ierr

    do imember = 1, this%n_ensemble

      call MPI_Bcast(this%local_io_data(:, imember), this%state_size, &
                     MPI_DOUBLE_PRECISION, this%io_ranks(imember), &
                     this%comm, ierr)

    end do

    local_io_data = this%local_io_data

  end function get_ensemble_state

end module dummy_model_interfaces
