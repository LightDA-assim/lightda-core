#include "mpi_types.h"

module assimilation_batch_manager
  use system_mpi
  use iso_c_binding
  use random, ONLY: random_normal
  use random_integer, ONLY: randint
  use assimilation_model_interface
  use distributed_array, ONLY: darray, darray_segment, new_darray
  use exceptions, ONLY: error_status, throw, new_exception
  use util, ONLY: str

  implicit none

  type :: assim_batch_manager

     !! Batch manager for data assimilation. Divides the model state into
     !! batches, assigns batches to MPI processor ranks, and handles
     !! transmission of model state data, observations, and predictions
     !! to the processor ranks assigned to assimilate each batch.

    private
    class(base_model_interface), pointer, public::model_interface
         !! Interface to the model
    integer, allocatable::batch_ranks(:)
         !! Array of MPI ranks assigned to compute assimilation on each batch
    integer::n_ensemble
         !! Number of ensemble members
    integer::state_size
         !! Size of model state
    integer::n_observations
         !! Number of observations
    integer::n_batches
         !! Number of batches
    integer::batch_size
         !! Maximum batch size
    integer::n_local_batches
         !! Number of batches assigned to the local processor rank
    MPI_COMM_TYPE::comm
         !! MPI communicator
    logical, allocatable::batch_results_received(:, :)
         !! Whether assimilation results have been received for each batch
  contains
    procedure::load_ensemble_state
    procedure::receive_results
    procedure::get_comm
    procedure::get_batch_size
    procedure::get_state_size
    procedure::get_n_batches
    procedure::get_n_ensemble
    procedure::add_obs_err
    procedure::localize
    procedure::store_results
    procedure::get_rank_batch_count
    procedure::get_rank_batches
    procedure::get_batch_offset
    procedure::get_batch_length
    procedure::get_batches_darray
  end type assim_batch_manager

contains

  function new_batch_manager(model_interface, n_ensemble, state_size, &
                             batch_size, comm)
    !! Create a new batch manager

    ! Arguments
    integer, intent(in)::n_ensemble
        !! Number of ensemble members
    integer, intent(in)::state_size
        !! Number of elements in model state
    integer, intent(in)::batch_size
        !! Maximum batch size
    MPI_COMM_TYPE::comm
        !! MPI communicator
    class(base_model_interface), intent(inout), target::model_interface
        !! Model interface

    type(assim_batch_manager)::new_batch_manager
        !! Newly created batch manager

    integer::ierr        ! MPI status code
    integer::rank        ! MPI rank
    integer::comm_size   ! Size of MPI communicator

    ! Initialize state info
    new_batch_manager%comm = comm
    new_batch_manager%model_interface => model_interface
    new_batch_manager%n_ensemble = n_ensemble
    new_batch_manager%state_size = state_size
    new_batch_manager%n_observations = 0
    new_batch_manager%batch_size = batch_size

    call mpi_comm_size(comm, comm_size, ierr)
    call mpi_comm_rank(comm, rank, ierr)

    new_batch_manager%n_batches = get_batch_count(state_size, batch_size)

    allocate (new_batch_manager%batch_ranks(new_batch_manager%n_batches))
    allocate (new_batch_manager%batch_results_received( &
              new_batch_manager%n_ensemble, new_batch_manager%n_batches))

    new_batch_manager%batch_results_received = .false.

    ! Assign batches to process ranks
    call get_batch_ranks(comm_size, new_batch_manager%batch_ranks)

    ! Get number of local batches
    new_batch_manager%n_local_batches = &
      new_batch_manager%get_rank_batch_count(rank)

  end function new_batch_manager

  function get_comm(this) result(comm)

    !! Get batch manager's MPI communicator

    ! Arguments
    class(assim_batch_manager)::this
        !! Batch manager

    ! Result
    MPI_COMM_TYPE::comm
        !! MPI communicator

    comm = this%comm

  end function get_comm

  function get_batch_size(this) result(batch_size)

    !! Get batch manager's maximum batch size

    ! Arguments
    class(assim_batch_manager)::this
        !! Batch manager

    ! Result
    integer::batch_size
        !! Batch size

    batch_size = this%batch_size

  end function get_batch_size

  function get_state_size(this) result(state_size)

    !! Get batch model state size

    ! Arguments
    class(assim_batch_manager)::this
        !! Batch manager

    ! Result
    integer::state_size
        !! State size

    state_size = this%state_size

  end function get_state_size

  function get_n_batches(this) result(n_batches)

    !! Get number of batches

    ! Arguments
    class(assim_batch_manager)::this
        !! Batch manager

    ! Result
    integer::n_batches
        !! Number of batches

    n_batches = this%n_batches

  end function get_n_batches

  function get_n_ensemble(this) result(n_ensemble)

    !! Get ensemble size

    ! Arguments
    class(assim_batch_manager)::this
        !! Batch manager

    ! Result
    integer::n_ensemble
        !! Ensemble size

    n_ensemble = this%n_ensemble

  end function get_n_ensemble

  subroutine shuffle(a)

    !! Randomize the order of elements in an array in-place

    ! Arguments
    integer, intent(inout) :: a(:)
        !! Array to reorder

    integer :: i, randpos, temp

    do i = size(a), 2, -1
      randpos = randint(i)
      temp = a(randpos)
      a(randpos) = a(i)
      a(i) = temp
    end do

  end subroutine shuffle

  subroutine get_batch_ranks(comm_size, batch_ranks)

    !! Assign MPI ranks to a set of assimilation batches

    ! Arguments
    integer, intent(in)::comm_size
        !! MPI communicator size

    integer, intent(inout)::batch_ranks(:)
        !! Array of batch ranks

    integer::random_size
        !! Size of random number seed
    integer, dimension(:), allocatable::seed
        !! Random number seed

    integer::i, ibatch !! Loop counters
    integer::rank !! MPI rank

    ! Start at highest rank (so if distribution is uneven, rank 0 gets
    ! fewer batches)
    rank = comm_size

    do ibatch = 1, size(batch_ranks)
      ! Decrement rank
      rank = rank - 1

      ! raise to comm_size-1 if we reach zero
      if (rank < 0) rank = comm_size - 1

      batch_ranks(ibatch) = rank
    end do

    ! Initialize random number generator
    call random_seed(size=random_size)
    allocate (seed(random_size))
    seed = 576834934 + 37*(/(i - 1, i=1, random_size)/)
    call random_seed(put=seed)

    ! Shuffle the batch ranks (to distribute load in case difficult segments
    ! of the domain are clustered near each other in the state array)
    call shuffle(batch_ranks)

  end subroutine get_batch_ranks

  function get_batch_count(state_size, batch_size) result(count)

    !! Get the number of batches required to cover the given model state size
    !! with the given batch size

    ! Arguments
    integer, intent(in)::state_size
        !! Model state size
    integer, intent(in)::batch_size
        !! Batch size

    ! Result
    integer::count
        !! Number of batches

    count = state_size/batch_size + min(mod(state_size, batch_size), 1)

  end function get_batch_count

  function get_batch_offset(this, ibatch) result(offset)

    !! Locate a batch in the model state array

    ! Arguments
    class(assim_batch_manager)::this
        !! Batch manager
    integer, intent(in)::ibatch
        !! Batch index

    ! Result
    integer::offset
        !! Offset from the start of the model state array

    offset = this%batch_size*(ibatch - 1)

  end function get_batch_offset

  function get_batch_length(this, ibatch) result(length)

    !! Get the length of a given batch

    ! Arguments
    class(assim_batch_manager)::this
        !! Batch manager
    integer, intent(in)::ibatch
        !! Batch index

    ! Result
    integer::length
        !! Length of batch

    integer::offset
        !! Offset of batch from start of the model state array

    offset = this%get_batch_offset(ibatch)
    length = min(this%state_size, offset + this%batch_size) - offset

  end function get_batch_length

  function get_rank_batch_count(this, rank) result(count)

    !! Get the number of batches assigned to a given MPI processor rank

    ! Arguments
    class(assim_batch_manager)::this
        !! Batch manager
    integer, intent(in)::rank
        !! MPI processor rank

    ! Result
    integer::count
        !! Number of batches assigned to `rank`

    integer::ibatch ! Loop counter

    count = 0

    do ibatch = 1, this%n_batches
      if (this%batch_ranks(ibatch) == rank) count = count + 1
    end do

  end function get_rank_batch_count

  subroutine get_rank_batches(this, rank, batches, status)

    !! Get an array of batches assigned to the MPI processor `rank`.

    ! Arguments
    class(assim_batch_manager)::this
        !! Batch manager
    integer, intent(in)::rank
        !! MPI processor rank
    integer, intent(out)::batches(:)
        !! Array of batch indices assigned to `rank`
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Result
    integer::n_batches
        !! Number of batches assigned to `rank`

    integer::ibatch_rank, ibatch ! Loop counters
    integer::ierr ! MPI status code

    character(:), allocatable :: errstr

    n_batches = this%get_rank_batch_count(rank)

    if (size(batches, 1) /= n_batches) then
      errstr = 'Wrong array size passed to get_rank_batches. &
           &Expected '//str(n_batches)//' got '//str(size(batches, 1))
      call throw(status, new_exception(errstr, 'get_rank_batches'))
      return
    end if

    ibatch_rank = 1

    do ibatch = 1, this%n_batches
      if (this%batch_ranks(ibatch) == rank) then
        batches(ibatch_rank) = ibatch
        ibatch_rank = ibatch_rank + 1
      end if
    end do

  end subroutine get_rank_batches

  function get_batches_darray(this) result(batches_darray)

    !! Distributed array with the size of the model state and

    class(assim_batch_manager)::this
        !! Batch manager

    type(darray)::batches_darray
        !! darray of model state with segments aligning to batch ranks

    type(darray_segment), target::batch_segments(this%n_batches)
    ! darray segments that make up batches_darray

    type(darray_segment), pointer::batch_segment
    ! Pointer to a member of batch_segments

    real(kind=8)::empty_batch(this%batch_size)
    ! Array of zeros to populate the array

    integer::rank ! MPI rank
    integer::ierr ! MPI status code

    integer::ibatch ! Loop counter

    call mpi_comm_rank(this%comm, rank, ierr)

    empty_batch = 0

    ! Build darray segments for each batch
    do ibatch = 1, this%n_batches

      batch_segment => batch_segments(ibatch)
      batch_segment%offset = this%get_batch_offset(ibatch)
      batch_segment%length = this%get_batch_length(ibatch)
      batch_segment%rank = this%batch_ranks(ibatch)
      batch_segment%comm = this%comm

      if (this%batch_ranks(ibatch) == rank) then

        ! A straightforward allocate() call doesn't work on members of
        ! derived types, so instead we use assignment which causes the
        ! array to be allocated implictly. Copying an array of the
        ! correct size triggers the required allocation.

        batch_segment%data = empty_batch(1:batch_segment%length)

      end if
    end do

    ! Create the darray from the array of segments
    batches_darray = new_darray(batch_segments, this%comm)

  end function get_batches_darray

  subroutine load_ensemble_state(this, istep, local_batches)

    !! Get the ensemble state from the model interface, divide into
    !! assimilation batches, and transmit the batch data to each processor

    ! Arguments
    class(assim_batch_manager)::this
        !! Batch manager
    integer, intent(in)::istep
        !! Assimilation step
    real(kind=8), intent(out), target::local_batches( &
                                        this%batch_size, this%n_local_batches, &
                                        this%n_ensemble)
        !! Array of local batch data

    integer::rank ! MPI rank
    integer::ierr ! MPI status code
    integer::imember, ibatch, ibatch_local ! Loop counters

    type(darray)::state_darray   ! Model state darray from the model interface
    type(darray)::batches_darray ! darray with segments aligning to batch ranks

    call mpi_comm_rank(this%comm, rank, ierr)

    call this%model_interface%read_state(istep)

    ! Get the assimilation batches darray
    batches_darray = this%get_batches_darray()

    do imember = 1, this%n_ensemble

      ! Get the model state darray from the model interface
      state_darray = this%model_interface%get_state_darray(istep, imember)

      ! Transfer the data to the batches array (this sends all the data to the
      ! correct processor ranks for processing)
      call state_darray%transfer_to_darray(batches_darray)

      ! Reset the local batch counter
      ibatch_local = 1

      do ibatch = 1, this%n_batches
        if (this%batch_ranks(ibatch) == rank) then

          ! Copy batch data to the local_batches array
          local_batches(:, ibatch_local, imember) = &
            batches_darray%segments(ibatch)%data

          ! Increment the local batch counter
          ibatch_local = ibatch_local + 1

        end if
      end do

    end do

  end subroutine load_ensemble_state

  subroutine receive_results(this, istep, local_batches)

    !! Receive assimilation results and transmit them to the model interface

    ! Arguments
    class(assim_batch_manager)::this
        !! Batch manager
    integer, intent(in)::istep
        !! Assimilation step
    real(kind=8), intent(in)::local_batches( &
                               this%batch_size, this%n_local_batches, &
                               this%n_ensemble)
        !! Array of local batch data

    integer::imember, ibatch, ibatch_local, isegment ! Loop counters
    integer::ierr ! MPI status code
    integer::rank ! MPI processor rank

    type(darray), target::state_darray
    ! Model state darray from the model interface
    type(darray)::batches_darray
    ! darray with segments aligning to batch ranks

    type(darray_segment), pointer::state_segment

    call mpi_comm_rank(this%comm, rank, ierr)

    ! Get the assimilation batches darray
    batches_darray = this%get_batches_darray()

    do imember = 1, this%n_ensemble

      ! Reset local batch counter
      ibatch_local = 1

      do ibatch = 1, this%n_batches
        if (this%batch_ranks(ibatch) == rank) then

          ! Copy batch data to the local_batches array
          batches_darray%segments(ibatch)%data = &
            local_batches(:, ibatch_local, imember)

          ! Increment the local batch counter
          ibatch_local = ibatch_local + 1

        end if

      end do

      ! Get the model state darray from the model interface
      state_darray = this%model_interface%get_state_darray(istep, imember)

      ! Transfer the data to state_darray
      call batches_darray%transfer_to_darray(state_darray)

      do isegment = 1, size(state_darray%segments)

        state_segment => state_darray%segments(isegment)

        if (state_segment%rank == rank) then

          ! Store the data in the model interface
          call this%model_interface%set_state_subset( &
            istep, imember, state_segment%offset, &
            state_segment%length, state_segment%data)

        end if

      end do

    end do

    this%batch_results_received = .true.

  end subroutine receive_results

  subroutine store_results(this, istep, local_batches)

    !! Store the assimilation results

    ! Arguments
    class(assim_batch_manager)::this
        !! Batch manager
    integer, intent(in)::istep
        !! Assimilation step
    real(kind=8), intent(in)::local_batches( &
                               this%batch_size, this%n_local_batches, &
                               this%n_ensemble)
        !! Assimilation results

    integer::imember, rank, ierr

    do while (any(this%batch_results_received .eqv. .false.))
      ! Receive any pending results
      call this%receive_results(istep, local_batches)
    end do

    ! Tell the model interface that we've finished receiving results
    call this%model_interface%write_state(istep)

  end subroutine store_results

  SUBROUTINE add_obs_err(this, istep, ibatch, dim_obs, HPH)

    !! Add observation error covariance matrix to the array `HPH`

    ! Arguments
    class(assim_batch_manager)::this
        !! Batch manager
    integer, intent(in), value :: istep
        !! Assimilation index
    integer, intent(in)::ibatch
        !! Batch index
    integer, intent(in)::dim_obs
        !! Dimensions of the observation array
    real(kind=8), INTENT(inout) :: HPH(dim_obs, dim_obs)
        !! Input array
    class(error_status), allocatable::status

    real(kind=8)::obs_err(dim_obs)
        !! Observation errors

    integer::iobs           !! Loop counter
    integer::batch_offset   !! Location of batch in the model state array
    integer::batch_length   !! Batch size

    ! Locate batch in the state array
    batch_offset = this%get_batch_offset(ibatch)
    batch_length = this%get_batch_length(ibatch)

    call throw(status, new_exception('Not yet implemented','add_obs_err'))
    return

    ! Get the observation errors
    !call this%get_subset_obs_err(istep, batch_offset, &
    !                                             batch_length, obs_err)

    do iobs = 1, dim_obs
      HPH(iobs, iobs) = HPH(iobs, iobs) + obs_err(iobs)**2
    end do

  END SUBROUTINE add_obs_err

  SUBROUTINE localize(this, istep, ibatch, dim_p, dim_obs, HP_p, HPH)
    !! Apply localization to HP and HPH^T

    ! Arguments
    class(assim_batch_manager)::this
        !! Batch manager
    integer, intent(in)::istep
        !! Assimilation step
    integer, intent(in)::ibatch
        !! Batch index
    integer, intent(in)::dim_p
        !! Batch size
    integer, intent(in)::dim_obs
        !! Number of observations

    real(kind=8), INTENT(inout) :: HP_p(dim_obs, dim_p)
        !! HP array
    real(kind=8), INTENT(inout) :: HPH(dim_obs, dim_obs)
        !! HPH array

    integer::domain_size

    integer::batch_offset     ! Location of batch in the model state array
    integer::batch_length     ! Batch size

    real(kind=8)::w           ! Localization weight

    integer::iobs1, iobs2, ipos ! Loop counters

    class(error_status), allocatable::status

    ! Locate batch in the state array
    batch_offset = this%get_batch_offset(ibatch)
    batch_length = this%get_batch_length(ibatch)

    if (dim_p /= batch_length) then
      print *, 'Inconsistent batch size'
      stop
    end if

    call throw(status, new_exception('Not yet implemented','localize'))
    return

    do iobs1 = 1, dim_obs

      do iobs2 = 1, dim_obs

        ! Get localization weights
        !w = this%model_interface%get_weight_obs_obs(istep, iobs1, iobs2)

        ! Multiply HPH by the localization weights
        HPH(iobs1, iobs2) = HPH(iobs1, iobs2)*w

      end do

      do ipos = 1, dim_p

        ! Get localization weights
        !w = this%model_interface%get_weight_model_obs( &
        !    istep, ipos + batch_offset, iobs1)

        ! Multiply HP_p by the localization weights
        HP_p(iobs1, ipos) = HP_p(iobs1, ipos)*w

      end do
    end do

  END SUBROUTINE localize

end module assimilation_batch_manager
