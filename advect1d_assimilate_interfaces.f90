module advect1d_assimilate_interfaces
  use mpi
  use iso_c_binding

  use hdf5

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

    call h5open_f(ierr)

  end subroutine init_interface

  subroutine get_io_ranks(comm_size,n_ensemble,io_ranks)
    integer,intent(in)::comm_size,n_ensemble
    integer,intent(inout)::io_ranks(n_ensemble)
    integer::i,stride
    stride=max(comm_size/n_ensemble,1)

    do i=1,n_ensemble
       io_ranks(i)=mod(i,comm_size)*stride
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
    type(io_info),pointer::info
    integer::i,ierr,io_rank
    character(len=50)::preassim_filename
    integer(HID_T)::h5file_h,dset_h,dataspace,memspace
    real(kind=8)::member_state(state_size)
    integer(HSIZE_T)::offset(2),count(2),stride(2),block(2)

    call c_f_pointer(info_c_ptr,info)

    call mpi_comm_rank(info%comm,io_rank,ierr)

    do i=1,n_ensemble
       if(info%io_ranks(i)==io_rank) then

          stride=(/2,1/)
          offset=(/i-1,0/)
          count=(/1,1/)
          block=(/1,state_size/)

          ! Set the HDF5 filename
          write(preassim_filename,"(A,I0,A)") &
               'ensembles/',istep,'/preassim.h5'

          ! Open the file
          call h5fopen_f(preassim_filename,h5F_ACC_RDONLY_F,h5file_h,ierr)

          ! Open the dataset
          call h5dopen_f(h5file_h,'ensemble_state',dset_h,ierr)
          call h5dget_space_f(dset_h,dataspace,ierr)
          call h5screate_simple_f(2,block,memspace,ierr)
          call h5sselect_hyperslab_f(dataspace,H5S_SELECT_SET_F, &
               offset,count,ierr,stride,block)
          call h5dread_f(dset_h,H5T_IEEE_F64LE,member_state,block,ierr,file_space_id=dataspace,mem_space_id=memspace)
          call h5dclose_f(dset_h,ierr)

       end if
    end do

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
