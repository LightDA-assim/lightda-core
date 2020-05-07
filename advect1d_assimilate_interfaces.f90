module advect1d_assimilate_interfaces
  use mpi
  use iso_c_binding
  use assimilate, ONLY: get_batch_offset, get_batch_length

  use hdf5

  implicit none

  type :: io_info
     integer,allocatable::io_ranks(:)
     integer::comm_size,n_ensemble,state_size,local_io_size
     real(kind=8),allocatable::local_io_data(:,:)
  end type io_info

contains

  subroutine init_interface(info_ptr,n_ensemble,state_size,comm_size,rank)
    type(c_ptr),intent(out)::info_ptr
    integer(c_int),intent(in)::n_ensemble,state_size,comm_size
    type(io_info), pointer ::info
    integer::ierr,rank

    ! Initialize state info
    allocate(info)
    info%n_ensemble=n_ensemble
    info%state_size=state_size

    allocate(info%io_ranks(n_ensemble))

    ! Assign ranks to processors for i/o purposes
    call get_io_ranks(comm_size,n_ensemble,info%io_ranks)

    info%local_io_size=get_rank_io_size(n_ensemble,info%io_ranks,rank)

    allocate(info%local_io_data(state_size,info%local_io_size))

    info_ptr=c_loc(info)

    call h5open_f(ierr)

  end subroutine init_interface

  subroutine get_io_ranks(comm_size,n_ensemble,io_ranks)
    integer,intent(in)::comm_size,n_ensemble
    integer,intent(inout)::io_ranks(n_ensemble)
    integer::i,stride
    stride=max(comm_size/n_ensemble,1)

    do i=1,n_ensemble
       io_ranks(i)=mod(i*stride,comm_size)
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

  subroutine read_state(istep,imember,member_state,state_size)
    
    integer,intent(in)::istep,imember,state_size
    real(kind=8),intent(inout)::member_state(state_size)
    character(len=50)::preassim_filename
    integer(HID_T)::h5file_h,dset_h,dataspace,memspace
    integer(HSIZE_T)::offset(2),count(2),stride(2),block(2)
    integer::ierr

    stride=(/2,1/)
    offset=(/imember-1,0/)
    count=(/1,1/)
    block=(/1,state_size/)

    ! Set the HDF5 filename
      write(preassim_filename,"(A,I0,A)") &
           'ensembles/',istep,'/preassim.h5'

      ! Open the file
      call h5fopen_f(preassim_filename,h5F_ACC_RDONLY_F,h5file_h,ierr)

      ! Open the dataset
      call h5dopen_f(h5file_h,'ensemble_state',dset_h,ierr)

      ! Define a dataspace within the dataset so we only read the
      ! specified ensemble member
      call h5dget_space_f(dset_h,dataspace,ierr)
      call h5sselect_hyperslab_f(dataspace,H5S_SELECT_SET_F, &
           offset,count,ierr,stride,block)

      ! Memory dataspace (needed since the local array shape differs
      ! from the dataspace in the file)
      call h5screate_simple_f(2,block,memspace,ierr)

      ! Read the data
      call h5dread_f(dset_h,H5T_IEEE_F64LE,member_state,block,ierr,file_space_id=dataspace,mem_space_id=memspace)

      ! Close the dataspaces
      call h5sclose_f(dataspace,ierr)
      call h5sclose_f(memspace,ierr)

      ! Close the dataset
      call h5dclose_f(dset_h,ierr)

      ! Close the file
      call h5fclose_f(h5file_h,ierr)
  end subroutine read_state

  subroutine load_ensemble_state(info_c_ptr,istep,rank,comm,state_size, &
    batch_ranks,local_batches,n_batches,n_local_batches,batch_size,n_ensemble)

    type(c_ptr),intent(inout)::info_c_ptr
    integer,intent(in)::istep,rank,comm,n_local_batches,state_size,batch_size,n_batches,n_ensemble
    integer,intent(in)::batch_ranks(n_batches)
    real(kind=8),intent(out)::local_batches(n_local_batches,batch_size,n_ensemble)
    real(kind=8)::batch_state(batch_size)
    type(io_info),pointer::info
    integer::imember,ierr,ibatch,ibatch_local,intercomm,comm_size,batch_length,batch_offset,local_io_counter
    real(kind=8)::member_state(state_size)
    integer::status(MPI_STATUS_SIZE)
    integer::request

    call c_f_pointer(info_c_ptr,info)

    call mpi_comm_size(comm,comm_size,ierr)

    local_io_counter=1

    do imember=1,n_ensemble
       if(info%io_ranks(imember)==rank) then

          ! Read batch data from disk
          call read_state(istep,imember,member_state,state_size)

          ! Copy data to local io array for later use
          info%local_io_data(:,local_io_counter)=member_state

          local_io_counter=local_io_counter+1

          ! Reset the local batch counter
          ibatch_local=1

          do ibatch=1,n_batches

             ! Locate batch in the state array
             batch_offset=get_batch_offset(batch_size,ibatch)
             batch_length=get_batch_length(batch_size,ibatch,state_size)

             ! Array containing the state for this batch
             batch_state=member_state(ibatch)

             if(batch_ranks(ibatch)==rank) then

                ! Copy batch state into the local_batches array
                local_batches(ibatch_local,:,imember)=batch_state

                ! Increment the local batch counter
                ibatch_local=ibatch_local+1

             else
                ! Send the batch to the appropriate process
                call MPI_Isend(batch_state(1:batch_length),batch_length, &
                     MPI_DOUBLE_PRECISION,batch_ranks(ibatch),ibatch,comm, &
                     request,ierr)
             end if
          end do

       end if
    end do

    do imember=1,n_ensemble
       if(info%io_ranks(imember)/=rank) then

          ! Reset the local batch counter
          ibatch_local=1

          do ibatch=1,n_batches
             if(batch_ranks(ibatch)==rank) then

                ! Determine the length of this batch
                batch_length=get_batch_length(batch_size,ibatch,state_size)

                ! Receive the batch data from the i/o process
                call MPI_Recv( &
                     local_batches( ibatch_local,1:batch_length,imember), &
                     batch_length,MPI_DOUBLE_PRECISION,info%io_ranks(imember), &
                     ibatch,comm,status,ierr)

                ! Increment the local batch counter
                ibatch_local=ibatch_local+1

             end if
          end do
       end if
    end do

  end subroutine load_ensemble_state

  subroutine transmit_results()
  end subroutine transmit_results

  subroutine write_state(istep,imember,n_ensemble,comm,member_state,state_size)

    integer,intent(in)::istep,imember,n_ensemble,state_size
    real(kind=8),intent(in)::member_state(state_size)
    character(len=80)::postassim_filename
    integer(HID_T)::h5file_h,dset_h,dataspace,memspace
    integer(HSIZE_T)::offset(2),count(2),stride(2),data_block(2)
    logical::exists
    integer::mpi_info
    integer::comm,ierr

    mpi_info=MPI_INFO_NULL

    stride=(/2,1/)
    offset=(/imember-1,0/)
    count=(/1,1/)
    data_block=(/1,state_size/)

    ! Set the HDF5 filename
    write(postassim_filename,"(A,I0,A,I0,A)") &
         'ensembles/',istep,'/postassim_',imember-1,'.h5'

    ! Create the file
    call h5fcreate_f(postassim_filename,H5F_ACC_TRUNC_F,h5file_h,ierr)

    ! Define a dataspace
    call h5screate_simple_f(1,(/int(state_size,8)/),dataspace,ierr)

    ! Create the dataset
    call h5dcreate_f(h5file_h,'ensemble_state',H5T_IEEE_F64LE, &
         dataspace,dset_h,ierr)

    ! Write the data
    call h5dwrite_f(dset_h,H5T_IEEE_F64LE,member_state,data_block,ierr,file_space_id=dataspace)

    ! Close the dataspace
    call h5sclose_f(dataspace,ierr)

    ! Close the dataset
    call h5dclose_f(dset_h,ierr)

    ! Close the file
    call h5fclose_f(h5file_h,ierr)

  end subroutine write_state

  subroutine store_results(info_c_ptr,istep,rank,comm,state_size)
    implicit none

    type(c_ptr),intent(inout)::info_c_ptr
    integer,intent(in)::istep,rank,comm,state_size
    type(io_info),pointer::info
    integer::imember,local_io_counter

    call c_f_pointer(info_c_ptr,info)

    local_io_counter=1

    do imember=1,info%n_ensemble
       if(info%io_ranks(imember)==rank) then

          call write_state(istep,imember,info%n_ensemble,comm, &
               info%local_io_data(:,local_io_counter),state_size)

          local_io_counter=local_io_counter+1

       end if
    end do

  end subroutine store_results

  subroutine cleanup_interface(info_ptr)
    type(c_ptr),intent(in)::info_ptr
    type(io_info), pointer ::info

    call c_f_pointer(info_ptr,info)

    deallocate(info%io_ranks)
    deallocate(info%local_io_data)
    deallocate(info)

  end subroutine cleanup_interface

end module advect1d_assimilate_interfaces
