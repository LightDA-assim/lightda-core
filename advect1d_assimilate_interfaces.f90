module advect1d_assimilate_interfaces
  use mpi

contains

  subroutine get_io_ranks(comm_size,n_ensemble,io_ranks)
    integer,intent(in)::comm_size
    integer,intent(inout)::io_ranks(n_ensemble)
    integer::i,stride
    stride=max(comm_size/n_ensemble,1)

    do i=1,n_ensemble
       io_ranks(i)=mod(i,comm_size)*stride+1
    end do

  end subroutine get_io_ranks

  function get_rank_io_size(n_ensemble,io_ranks,rank) result(size)
    integer,intent(in)::n_ensemble
    integer,intent(in)::io_ranks(n_ensemble)
    integer,intent(in)::rank
    integer::size

    size=0

    do i=1,n_ensemble
       if(io_ranks(i)==rank) size=size+1
    end do

  end function get_rank_io_size

  subroutine load_ensemble_state(istep,rank,comm,state_size,batch_ranks, &
       local_batches,n_batches,n_local_batches,batch_size,n_ensemble)

    integer,intent(in)::istep,rank,comm,n_local_batches,state_size,batch_size,n_batches
    integer,intent(in)::batch_ranks(n_batches)
    real(kind=8),intent(out)::local_batches(n_local_batches,state_size,n_ensemble)

    
    
  end subroutine load_ensemble_state

  subroutine transmit_results()
  end subroutine transmit_results

end module advect1d_assimilate_interfaces
