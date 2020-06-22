module model_interface_tests
  use assimilation_model_interface, ONLY: base_model_interface
  use system_mpi
contains

  subroutine test_buffer_length(iface)
    class(base_model_interface)::iface
    real(kind=8),pointer::buf(:)
    integer::length,sum_buffer_lengths,ierr

    ! Set requested buffer length
    length=min(iface%get_state_size(),5)

    ! Get buffer pointer
    buf=>iface%get_state_subset_buffer(1,1,0,length)

    ! Add up buffer lengths across all processors
    call MPI_Allreduce(size(buf),sum_buffer_lengths,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)

    ! Check that buffer lengths add up to the requested length
    if(sum_buffer_lengths/=length) then
       print *,'Buffer of length',length,'expected, got length',size(buf)
       error stop
    end if

  end subroutine test_buffer_length

  subroutine test_buffer_readwrite(iface)
    class(base_model_interface)::iface
    real(kind=8),pointer::buf(:)
    integer::length,ierr,i
    integer,parameter::shift=1

    ! Set requested buffer length
    length=min(iface%get_state_size(),5)

    if(length>0) then

       ! Get buffer pointer
       buf=>iface%get_state_subset_buffer(1,1,0,length)

       ! Write to buffer
       do i=1,size(buf)
          buf(i)=i
       end do

       ! Get new buffer, shifted relative to original
       buf=>iface%get_state_subset_buffer(1,1,shift,length)

       ! Check values in new buffer
       do i=1,size(buf)-1
          if(buf(i)/=i+shift) then
             print *,'Wrote value of',i+shift,'read back value of',buf(i+shift)
             error stop
          end if
       end do

    end if

  end subroutine test_buffer_readwrite

  subroutine run_all(iface)
    class(base_model_interface)::iface

    call test_buffer_length(iface)
    call test_buffer_readwrite(iface)

  end subroutine run_all
  
end module model_interface_tests
