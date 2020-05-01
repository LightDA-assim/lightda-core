module advect1d_assimilate_interfaces
  use mpi
  use iso_c_binding

  implicit none

  type :: io_info
     integer::comm
     integer,allocatable::io_ranks(:)
  end type io_info

contains

  subroutine init_interface(info_ptr,io_comm,n_ensemble)
    type(c_ptr),intent(out)::info_ptr
    integer(c_int),intent(in)::io_comm,n_ensemble
    type(io_info), pointer ::info
    integer::ierr,rank,comm_size

    ! Initialize state info
    allocate(info)
    allocate(info%io_ranks(n_ensemble))
    info_ptr=c_loc(info)

    info%comm=io_comm

    call mpi_comm_size(io_comm,comm_size,ierr)

    ! Assign ranks to processors for i/o purposes
    call get_io_ranks(comm_size,n_ensemble,info%io_ranks)

  end subroutine init_interface

  subroutine get_io_ranks(comm_size,n_ensemble,io_ranks)
    integer,intent(in)::comm_size,n_ensemble
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
    integer::i,size

    size=0

    do i=1,n_ensemble
       if(io_ranks(i)==rank) size=size+1
    end do

  end function get_rank_io_size

  subroutine load_ensemble_state(info_c_ptr,istep,rank,comm,state_size, &
    batch_ranks,local_batches,n_batches,n_local_batches,batch_size,n_ensemble)

    type(c_ptr),intent(inout)::info_c_ptr
    integer,intent(in)::istep,rank,comm,n_local_batches,state_size,batch_size,n_batches,n_ensemble
    integer,intent(in)::batch_ranks(n_batches)
    real(kind=8),intent(out)::local_batches(n_local_batches,state_size,n_ensemble)

  end subroutine load_ensemble_state

  subroutine transmit_results()
  end subroutine transmit_results

  subroutine cleanup_interface(info_ptr)
    type(c_ptr),intent(in)::info_ptr
    type(io_info), pointer ::info

    call c_f_pointer(info_ptr,info)

    deallocate(info%io_ranks)
    deallocate(info)

  end subroutine cleanup_interface

end module advect1d_assimilate_interfaces
