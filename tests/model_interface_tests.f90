module model_interface_tests
  use assimilation_model_interface, ONLY: base_model_interface
  use system_mpi
contains

  subroutine test_buffer_length(iface)
    class(base_model_interface)::iface
    real(kind=8),pointer::buf(:)
    integer::length,sum_buffer_lengths,ierr

    length=min(iface%get_state_size(),5)
    buf=>iface%get_state_subset_buffer(1,1,0,length)

    call MPI_Allreduce(size(buf),sum_buffer_lengths,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)

    if(sum_buffer_lengths/=length) then
       print *,'Buffer of length',length,'expected, got length',size(buf)
       error stop
    end if

  end subroutine test_buffer_length

  subroutine run_all(iface)
    class(base_model_interface)::iface

    call test_buffer_length(iface)

  end subroutine run_all
  
end module model_interface_tests
