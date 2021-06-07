#include "mpi_types.h"

module per_member_dummy_model_interfaces
  use per_member_model_interfaces, ONLY: per_member_model_interface
  use system_mpi
  use iso_c_binding
  use random, ONLY: random_normal
  use exceptions, ONLY: throw, new_exception, error_container
  use distributed_array, ONLY: darray, darray_segment, new_darray

  implicit none

  type, extends(per_member_model_interface)::per_member_dummy_model_interface
    private
    real(kind=8), allocatable::observations(:), obs_errors(:), predictions(:, :)
    integer, allocatable::obs_positions(:)
    integer::n_observations
    logical::observations_read, predictions_computed
  contains
    procedure::store_member_state
    procedure::load_member_state
  end type per_member_dummy_model_interface

contains

  function new_per_member_dummy_model( &
    n_ensemble, n_observations, state_size, comm) result(this)

    integer(c_int), intent(in)::n_ensemble, n_observations, state_size
    MPI_COMM_TYPE::comm
    type(per_member_dummy_model_interface)::this
    integer::ierr, rank, comm_size

    call this%initialize_per_member_model_interface( &
      n_ensemble, state_size, comm)

    this%observations_read = .false.
    this%predictions_computed = .false.

    call mpi_comm_rank(comm, rank, ierr)
    call mpi_comm_size(comm, comm_size, ierr)

    allocate (this%observations(this%n_observations))
    allocate (this%obs_errors(this%n_observations))
    allocate (this%obs_positions(this%n_observations))
    allocate (this%predictions(this%n_observations, this%n_ensemble))

  end function new_per_member_dummy_model

  subroutine load_member_state( &
    this, istep, imember, member_state, state_size, status)

    !! Read the model state from disk for a given ensemble member

    ! Arguments
    class(per_member_dummy_model_interface)::this
        !! Model interface
    integer, intent(in)::istep
        !! Iteration number
    integer, intent(in)::imember
        !! Ensemble member index
    integer, intent(in)::state_size
        !! Size of the model state
    real(kind=8), intent(inout)::member_state(state_size)
        !! On exit, array holding the model state values
    type(error_container), intent(out), optional::status
        !! Error status

    integer::i
      !! Loop counter

    do i = 1, state_size
       member_state(i) = imember*this%get_state_size()+i
    end do

  end subroutine load_member_state

  subroutine store_member_state( &
    this, imember, n_ensemble, member_state, state_size, status)

    !! Store the new state for one ensemble member

    use exceptions, ONLY: error_container

    ! Arguments
    class(per_member_dummy_model_interface)::this
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
  end subroutine store_member_state

end module per_member_dummy_model_interfaces
