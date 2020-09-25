#include "mpi_types.h"

module batch_manager_tests
  use system_mpi
  use assimilation_batch_manager, ONLY: assim_batch_manager, new_batch_manager
  use dummy_model_interfaces, ONLY: dummy_model_interface, new_dummy_model
  use exceptions, ONLY: throw, error_status, new_exception
  use dummy_assimilator
  use iso_c_binding
  use util, ONLY: str

  implicit none
contains

  subroutine test_batch_math()

    type(assim_batch_manager)::batch_manager
    type(dummy_model_interface)::model_interface
    integer, parameter::n_observations = 5
    integer, parameter::state_size = 100
    integer, parameter::batch_size = 7
    integer, parameter::n_ensemble = 15
    real(kind=8), allocatable::local_batches(:, :, :)
    integer, allocatable::local_batch_inds(:)
    integer::rank, comm_size, ierr, batch_length, batch_offset, ibatch, istep, &
              n_batches, sum_batch_lengths, offset, last_offset
    MPI_COMM_TYPE::comm
    real(kind=8), parameter::forget = 0.6

    comm = mpi_comm_world

    comm_size = 10

    model_interface = new_dummy_model( &
                      n_ensemble, n_observations, state_size, comm)

    batch_manager = new_batch_manager( &
                    model_interface, n_ensemble, state_size, &
                    batch_size, mpi_comm_world)

    n_batches = batch_manager%get_n_batches()
    sum_batch_lengths = 0
    last_offset = 0

    do ibatch = 1, n_batches
      batch_length = batch_manager%get_batch_length(ibatch)
      sum_batch_lengths = sum_batch_lengths + batch_length
      batch_offset = batch_manager%get_batch_length(ibatch)
      if (ibatch > 1 .and. batch_offset /= last_offset + batch_length) then
        print *, 'Offset incorrect for batch', ibatch, &
          '. Expected', last_offset + batch_length, 'got', batch_offset
        error stop
      end if
    end do

    if (sum_batch_lengths /= state_size) then
      print *, 'Sum of batch lengths was', sum_batch_lengths, &
        ' expected', state_size
      error stop
    end if

  end subroutine test_batch_math

end module batch_manager_tests

program test_batch_manager

  use system_mpi
  use batch_manager_tests
  implicit none

  integer ierr

  call mpi_init(ierr)
  call test_batch_math()
  call mpi_finalize(ierr)

end program test_batch_manager
