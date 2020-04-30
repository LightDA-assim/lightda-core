module advect1d_assimilate_interfaces
  include 'mpif.h'

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

  subroutine load_ensemble_state(comm)

    
    
  end subroutine load_ensemble_state

  subroutine transmit_results()
  end subroutine transmit_results

end module advect1d_assimilate_interfaces
