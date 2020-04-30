module assimilate
  implicit none

  include 'mpif.h'

contains

  subroutine shuffle(a)
    integer, intent(inout) :: a(:)
    integer :: i, randpos, temp
    real :: r
 
    do i = size(a), 2, -1
       call random_number(r)
       randpos = int(r * i) + 1
       temp = a(randpos)
       a(randpos) = a(i)
       a(i) = temp
    end do
 
  end subroutine shuffle

  subroutine get_batch_ranks(comm_size,batch_ranks)
    integer,intent(in)::comm_size
    integer,intent(inout)::batch_ranks(:)
    integer::random_size
    integer,dimension(:),allocatable::seed
    integer::i,ibatch,rank

    ! Start at highest rank (so if distribution is uneven, rank 0 gets
    ! fewer batches)
    rank=comm_size+1

    ! 
    do ibatch=1,size(batch_ranks)
       ! Decrement rank
       rank=rank-1

       ! raise to comm_size if we reach zero
       if(rank<1) rank=comm_size

       batch_ranks(ibatch)=rank
    end do

    ! Initialize random number generator
    call random_seed(size=random_size)
    allocate(seed(random_size))
    seed = 576834934 + 37 * (/ (i-1, i=1, random_size) /)
    call random_seed(put=seed)

    ! Shuffle the batch ranks (to distribute load in case difficult segments
    ! of the domain are clustered near each other in the state array)
    call shuffle(batch_ranks)

  end subroutine get_batch_ranks

  function get_batch_count(state_size,batch_size) result(count)
    integer,intent(in)::state_size,batch_size
    integer::count
    count=state_size/batch_size+min(mod(state_size,batch_size),1)
  end function get_batch_count

  function get_batch_offset(batch_size,ibatch) result(offset)
    integer,intent(in)::batch_size,ibatch
    integer::offset
    offset=batch_size*ibatch
  end function get_batch_offset

  subroutine assimilate_parallel(istep,hrut,n_ensemble,batch_size,state_size,comm,U_load_ensemble_state)

    abstract interface

       subroutine load_ensemble_state()
         
       end subroutine load_ensemble_state

       subroutine transmit_results()
         
       end subroutine transmit_results

    end interface

    integer,intent(in) :: istep,n_ensemble,batch_size,state_size,comm
    real(kind=8) :: hrut
    integer,dimension(:),allocatable :: io_ranks, batch_ranks
    !real(kind=8),dimension(:,:),allocatable::

    integer::rank,ierr,comm_size,n_batches

    procedure(load_ensemble_state) :: U_load_ensemble_state

    call mpi_comm_rank(comm, rank, ierr)
    call mpi_comm_size(comm, comm_size, ierr)

    ! Get the number of batches and allocate batch arrays
    n_batches=get_batch_count(n_ensemble,batch_size)
    allocate(batch_ranks(n_batches))

    ! Assign batches to process ranks
    call get_batch_ranks(comm_size,batch_ranks)

    ! Load the ensemble state
    call U_load_ensemble_state()

  end subroutine assimilate_parallel

end module assimilate
