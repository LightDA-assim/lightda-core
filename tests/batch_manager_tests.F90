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

  subroutine localize(istep, ibatch, dim_p, dim_obs, HP_p, HPH, mgr)
    class(*), intent(in)::mgr
    INTEGER(c_int32_t), INTENT(in), value :: istep, ibatch, dim_p, dim_obs
    REAL(c_double), INTENT(inout) :: HP_p(dim_obs, dim_p), HPH(dim_obs, dim_obs)

    select type (mgr)
    class is (assim_batch_manager)
      call mgr%localize(istep, ibatch, dim_p, dim_obs, HP_p, HPH)
    class default
      print *, 'Could not determine argument type. &
           &Should be class assim_batch_manager'
    end select
  end subroutine localize

  subroutine add_obs_err(istep, ibatch, dim_obs, HPH, mgr)
    class(*), intent(in)::mgr
    INTEGER(c_int32_t), INTENT(in), value :: istep, ibatch, dim_obs
    REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)

    select type (mgr)
    class is (assim_batch_manager)
      call mgr%add_obs_err(istep, ibatch, dim_obs, HPH)
    class default
      print *, 'Could not determine argument type. &
           &Should be class assim_batch_manager'
    end select
  end subroutine add_obs_err

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

  subroutine test_empty_assimilator(status)

    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    type(assim_batch_manager)::batch_manager
    type(dummy_model_interface), target::model_interface
    integer, parameter::n_observations = 5
    integer, parameter::state_size = 100
    integer, parameter::batch_size = 5
    integer, parameter::n_ensemble = 15
    real(kind=8), allocatable::local_batches(:, :, :)
    integer, allocatable::local_batch_inds(:)
    real(kind=8)::innovations(n_observations, n_ensemble)
    real(kind=8)::predictions(n_observations, n_ensemble)
    real(kind=8)::observations(n_observations)
    real(kind=8)::obs_errors(n_observations)
    real(kind=8)::batch_mean_state(batch_size)
    real(kind=8)::batch_states(batch_size, n_ensemble)
    real(kind=8)::ensemble_before_assimilation(state_size, n_ensemble)
    real(kind=8)::ensemble_after_assimilation(state_size, n_ensemble)
    integer::rank, comm_size, ierr, batch_length, batch_offset, &
              ibatch, ibatch_local, istep, n_local_batches, n_obs_batch
    integer::imember, i
    real(kind=8)::delta
    MPI_COMM_TYPE::comm
    real(kind=8), parameter::forget = 0.6
    character(:), allocatable::errstr ! Error string

    istep = 1

    comm = mpi_comm_world

    call mpi_comm_rank(comm, rank, ierr)
    call mpi_comm_size(comm, comm_size, ierr)

    if (comm_size <= 1) then
      errstr = 'Not enough processes to run test. Have '//str(comm_size)// &
               ' process on MPI communicator, need at least 2.'
      call throw(status, new_exception(errstr, 'test_darray_transfer'))
      return
    end if

    model_interface = new_dummy_model( &
                      n_ensemble, n_observations, state_size, comm)

    batch_manager = new_batch_manager( &
                    model_interface, n_ensemble, state_size, &
                    batch_size, mpi_comm_world)

    ! Get number of local batches
    n_local_batches = batch_manager%get_rank_batch_count(rank)

    ! Allocate array to hold ensemble state
    allocate (local_batches(batch_size, n_local_batches, n_ensemble))

    ! Get local batch indices
    allocate (local_batch_inds(n_local_batches))
    call batch_manager%get_rank_batches(rank, local_batch_inds)

    ! Load the ensemble state
    call batch_manager%load_ensemble_state(istep, local_batches)

    ensemble_before_assimilation = model_interface%get_ensemble_state()

    ! Assimilate local batches
    do ibatch_local = 1, n_local_batches
      ibatch = local_batch_inds(ibatch_local)

      batch_offset = batch_manager%get_batch_offset(ibatch)
      batch_length = batch_manager%get_batch_length(ibatch)

      ! Get number of predictions for this batch
      n_obs_batch = model_interface%get_subset_obs_count( &
                    istep, batch_offset, batch_length)

      call model_interface%get_subset_predictions( &
        istep, batch_offset, batch_size, predictions)
      call model_interface%get_subset_observations( &
        istep, batch_offset, batch_size, observations)
      call model_interface%get_subset_obs_err( &
        istep, batch_offset, batch_size, obs_errors)
      call model_interface%get_innovations( &
        istep, batch_offset, batch_length, observations, predictions, &
        obs_errors, innovations)

      batch_states = local_batches(ibatch_local, :, :)

      call dummy_assimilator_assimilate( &
        istep, ibatch, batch_size, n_obs_batch, &
        n_obs_batch, n_ensemble, int(0), batch_mean_state, &
        batch_states, predictions, innovations, add_obs_err, &
        localize, forget, ierr, batch_manager)

      local_batches(ibatch_local, :, :) = batch_states

    end do

    call batch_manager%store_results(istep, local_batches)

    ensemble_after_assimilation = model_interface%get_ensemble_state()

    do imember = 1, n_ensemble
      do i = 1, state_size
        delta = ensemble_after_assimilation(i, imember) - &
                ensemble_before_assimilation(i, imember)
        if (abs(delta) > 1e-8) then
          print *, imember, i, ensemble_before_assimilation(i, imember), &
            ensemble_after_assimilation(i, imember), delta
        end if
      end do
    end do

    if (.not. all(abs(ensemble_before_assimilation - &
                      ensemble_after_assimilation) < 1e-8)) then
      print *, 'Ensemble state changed during assimilation or i/o on rank', rank
      error stop
    end if

  end subroutine test_empty_assimilator

end module batch_manager_tests

program test_batch_manager

  use system_mpi
  use batch_manager_tests
  implicit none

  integer ierr

  call mpi_init(ierr)
  call test_batch_math()
  call test_empty_assimilator()
  call mpi_finalize(ierr)

end program test_batch_manager
