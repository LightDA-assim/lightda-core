module advect1d_assimilate_interfaces
  use mpi
  use iso_c_binding
  use assimilate, ONLY: get_batch_offset, get_batch_length, get_batch_count
  use random, ONLY: random_normal

  use hdf5

  implicit none

  type :: io_info
     integer,allocatable::io_ranks(:)
     integer::comm_size,n_ensemble,state_size,local_io_size,n_observations,n_batches
     real(kind=8),allocatable::local_io_data(:,:)
     real(kind=8),allocatable::observations(:),obs_errors(:),predictions(:,:)
     integer,allocatable::obs_positions(:)
     logical,allocatable::batch_results_received(:,:)
     logical::observations_read
     logical::predictions_computed
  end type io_info

contains

  subroutine init_interface(info_ptr,n_ensemble,state_size,batch_size,comm_size,rank)
    type(c_ptr),intent(out)::info_ptr
    integer(c_int),intent(in)::n_ensemble,state_size,batch_size,comm_size
    type(io_info), pointer ::info
    integer::ierr,rank

    ! Initialize state info
    allocate(info)
    info%n_ensemble=n_ensemble
    info%state_size=state_size
    info%n_observations=0

    info%observations_read=.false.
    info%predictions_computed=.false.

    allocate(info%io_ranks(n_ensemble))

    ! Assign ranks to processors for i/o purposes
    call get_io_ranks(comm_size,n_ensemble,info%io_ranks)

    info%local_io_size=get_rank_io_size(n_ensemble,info%io_ranks,rank)

    info%n_batches=get_batch_count(state_size,batch_size)

    allocate(info%observations(info%n_observations))
    allocate(info%obs_errors(info%n_observations))
    allocate(info%obs_positions(info%n_observations))

    allocate(info%local_io_data(state_size,info%local_io_size))
    allocate(info%batch_results_received(info%local_io_size,info%n_batches))

    info%batch_results_received=.false.

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

  function get_local_io_index(n_ensemble,io_ranks,rank,imember) result(index)
    integer,intent(in)::n_ensemble,rank,imember
    integer,intent(in)::io_ranks(n_ensemble)
    integer::i,local_io_counter,index

    index=-1

    local_io_counter=0

    do i=1,n_ensemble
       if(io_ranks(i)==rank) then
          local_io_counter=local_io_counter+1
          if(i==imember) then
             index=local_io_counter
             exit
          end if
       end if
    end do

  end function get_local_io_index

  function get_batch_observation_count(info_ptr,istep,ibatch, &
       batch_size,comm,rank) result(n_obs_batch)

    use iso_c_binding

    implicit none

    type(c_ptr),intent(inout)::info_ptr
    integer,intent(in)::istep,ibatch,batch_size,comm,rank
    integer::n_obs_batch
    type(io_info),pointer::info

    call c_f_pointer(info_ptr,info)

    if(.not.info%observations_read) call read_observations(info,istep)

    n_obs_batch=info%n_observations

  end function get_batch_observation_count

  subroutine get_batch_predictions(info_ptr,istep,ibatch,batch_size,n_ensemble,n_obs_batch,rank,comm,predictions)

    use iso_c_binding
    use mpi

    implicit none

    type(c_ptr),intent(inout)::info_ptr
    integer,intent(in)::istep,ibatch,rank,comm,batch_size,n_ensemble, &
         n_obs_batch
    real(kind=8),intent(inout)::predictions(n_obs_batch,n_ensemble)
    integer::i,imember,ierr
    type(io_info),pointer::info

    call c_f_pointer(info_ptr,info)

    if(n_obs_batch /= info%n_observations) then
       print '(A,I0,A,I0,A,I0)', &
            'Wrong n_obs_batch passed to get_batch_predictions for batch ', &
            ibatch,'. Expected ',info%n_observations,', got ',n_obs_batch
       call mpi_abort(comm,1,ierr)
    end if

    call c_f_pointer(info_ptr,info)

    predictions=info%predictions

  end subroutine get_batch_predictions

  subroutine get_batch_innovations(info_c_ptr,istep,ibatch,batch_size,n_ensemble,n_obs_batch,rank,comm,innovations)

    use iso_c_binding

    implicit none

    type(c_ptr),intent(inout)::info_c_ptr
    integer,intent(in)::istep,ibatch,rank,comm,batch_size,n_ensemble, &
         n_obs_batch
    real(kind=8),intent(inout)::innovations(n_obs_batch,n_ensemble)
    type(io_info),pointer::info
    integer::imember,iobs

    call c_f_pointer(info_c_ptr,info)

    do imember=1,n_ensemble
       do iobs=1,n_obs_batch
          innovations(iobs,imember)=info%observations(iobs) - &
               info%predictions(iobs,imember) + &
               random_normal()*info%obs_errors(iobs)
       end do
    end do

  end subroutine get_batch_innovations

  subroutine get_batch_observations(info_c_ptr,istep,ibatch,batch_size,n_obs_batch,rank,comm,observations)

    use iso_c_binding
    use mpi

    implicit none

    type(c_ptr),intent(inout)::info_c_ptr
    integer,intent(in)::istep,ibatch,rank,comm,batch_size,n_obs_batch
    real(kind=8),intent(inout)::observations(n_obs_batch)
    type(io_info),pointer::info
    integer::ierr

    call c_f_pointer(info_c_ptr,info)

    if(n_obs_batch /= info%n_observations) then
       print '(A,I0,A,I0,A,I0)', &
            'Wrong n_obs_batch passed to get_batch_observations for batch ', &
            ibatch,'. Expected ',info%n_observations,', got ',n_obs_batch
       call mpi_abort(comm,1,ierr)
    end if

    observations=info%observations

  end subroutine get_batch_observations

  SUBROUTINE add_obs_err(step,ind_p,dim_obs,HPH,info_ptr)
    ! Add observation error covariance matrix
    USE iso_c_binding
    type(c_ptr),intent(inout)::info_ptr
    INTEGER(c_int32_t), INTENT(in), value :: step, ind_p, dim_obs
    REAL(c_double), INTENT(inout) :: HPH(dim_obs,dim_obs)
    type(io_info),pointer::info
    integer::iobs

    call c_f_pointer(info_ptr,info)

    do iobs=1,dim_obs
       HPH(iobs,iobs)=HPH(iobs,iobs)+info%obs_errors(iobs)
    end do

  END SUBROUTINE add_obs_err

  SUBROUTINE localize(step,ind_p,dim_p,dim_obs,HP_p,HPH,info_ptr)
    ! Apply localization to HP and HPH^T
    USE iso_c_binding
    INTEGER(c_int32_t), INTENT(in), value :: step, ind_p, dim_p, dim_obs
    REAL(c_double), INTENT(inout) :: HP_p(dim_obs,dim_p), HPH(dim_obs,dim_obs)
    type(c_ptr),intent(inout)::info_ptr
  END SUBROUTINE localize

  subroutine read_observations(info,istep)
    type(io_info),intent(inout), pointer :: info
    integer,intent(in)::istep
    character(len=50)::obs_filename
    integer(HID_T)::h5file_h,dset_h,dataspace
    integer(HSIZE_T)::dims(1),maxdims(1)
    integer::ierr,rank

    ! Set the HDF5 filename
    write(obs_filename,"(A,I0,A)") &
         'ensembles/',istep,'/observations.h5'

    ! Open the file
    call h5fopen_f(obs_filename,h5F_ACC_RDONLY_F,h5file_h,ierr)

    ! Open the observations dataset
    call h5dopen_f(h5file_h,'observations',dset_h,ierr)

    ! Get the dataspace handle
    call h5dget_space_f(dset_h,dataspace,ierr)

    ! Get the dataset size
    call h5sget_simple_extent_dims_f(dataspace,dims,maxdims,ierr)

    info%n_observations=dims(1)

    ! Close the dataspace
    call h5sclose_f(dataspace,ierr)

    deallocate(info%observations)
    allocate(info%observations(info%n_observations))

    ! Read the data
    call h5dread_f(dset_h,H5T_NATIVE_DOUBLE,info%observations,dims,ierr)

    ! Close the dataset
    call h5dclose_f(dset_h,ierr)

    ! Open the obs_positions dataset
    call h5dopen_f(h5file_h,'obs_positions',dset_h,ierr)

    deallocate(info%obs_positions)
    allocate(info%obs_positions(info%n_observations))

    ! Read the data
    call h5dread_f(dset_h,H5T_NATIVE_INTEGER,info%obs_positions,dims,ierr)

    ! Close the dataset
    call h5dclose_f(dset_h,ierr)

    ! Open the obs_errors dataset
    call h5dopen_f(h5file_h,'obs_errors',dset_h,ierr)

    deallocate(info%obs_errors)
    allocate(info%obs_errors(info%n_observations))

    ! Read the data
    call h5dread_f(dset_h,H5T_NATIVE_DOUBLE,info%obs_errors,dims,ierr)

    ! Close the dataset
    call h5dclose_f(dset_h,ierr)

    ! Close the file
    call h5fclose_f(h5file_h,ierr)

    info%observations_read=.true.

  end subroutine read_observations

  subroutine load_observations_parallel(info,istep,comm,rank)
    type(io_info),intent(inout), pointer :: info
    integer,intent(in)::istep,comm
    integer::ierr,rank

    if(info%observations_read) return

    if(rank==0) call read_observations(info,istep)

    call mpi_bcast(info%n_observations,1,MPI_INTEGER,0,comm,ierr)

    if(rank>0) then
       deallocate(info%observations)
       allocate(info%observations(info%n_observations))
       deallocate(info%obs_positions)
       allocate(info%obs_positions(info%n_observations))
       deallocate(info%obs_errors)
       allocate(info%obs_errors(info%n_observations))
    end if

    call mpi_bcast(info%observations,info%n_observations,MPI_INTEGER,0,comm, &
         ierr)
    call mpi_bcast(info%obs_positions,info%n_observations,MPI_INTEGER,0,comm, &
         ierr)
    call mpi_bcast(info%obs_errors,info%n_observations,MPI_INTEGER,0,comm, &
         ierr)

    info%observations_read=.true.

  end subroutine load_observations_parallel

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
    integer::imember,ierr,ibatch,ibatch_local,intercomm,comm_size,batch_length,batch_offset,local_io_counter,iobs
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
             batch_state=member_state(batch_offset+1:batch_offset+batch_length+1)

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

    ! Load observations
    call load_observations_parallel(info,istep,comm,rank)
    allocate(info%predictions(info%n_observations,n_ensemble))

    local_io_counter=1
    do imember=1,n_ensemble
       if(info%io_ranks(imember)==rank) then
          do iobs=1,info%n_observations
             info%predictions(iobs,imember)=info%local_io_data( &
                  info%obs_positions(iobs)+1,local_io_counter)
          end do
          local_io_counter=local_io_counter+1
       end if
       call mpi_bcast(info%predictions(:,imember),info%n_observations, &
            MPI_DOUBLE,info%io_ranks(imember),comm,ierr)
    end do

    info%predictions_computed=.true.

  end subroutine load_ensemble_state

  subroutine receive_results(info,istep,batch_size,batch_ranks,n_batches,n_ensemble,rank,comm)
    type(io_info),pointer::info
    integer,intent(in)::istep,batch_size,rank,comm,n_ensemble
    integer,intent(in)::batch_ranks(n_batches)
    integer::imember,ibatch,n_batches,batch_rank,batch_offset,batch_length, &
         member_rank,local_io_index,status(MPI_STATUS_SIZE),ierr
    logical::flag

    do imember=1,n_ensemble
       member_rank=info%io_ranks(imember)
       if(member_rank==rank) then
          local_io_index=get_local_io_index(n_ensemble,info%io_ranks, &
               rank,imember)
          do ibatch=1,n_batches

             batch_rank=batch_ranks(ibatch)

             ! Check whether batch has already been received
             if(info%batch_results_received(local_io_index,ibatch)) cycle

             ! Check whether data is ready to receive
             call mpi_iprobe(batch_rank,ibatch,comm,flag,status,ierr)
             if(flag.eqv..false.) cycle

             ! Locate batch in the state array
             batch_offset=get_batch_offset(batch_size,ibatch)
             batch_length=get_batch_length(batch_size,ibatch,info%state_size)

             ! Receive data
             call mpi_recv(info%local_io_data( &
                  batch_offset+1:batch_offset+batch_length,local_io_index), &
                  batch_length,MPI_DOUBLE_PRECISION,batch_rank,ibatch,comm, &
                  status,ierr)

             ! Record that batch was received
             info%batch_results_received(local_io_index,ibatch)=.true.

          end do
       end if
    end do

  end subroutine receive_results

  subroutine transmit_results(info_ptr,istep,ibatch,batch_state,batch_size,n_ensemble,batch_ranks,n_batches,comm,rank)
    type(c_ptr),intent(inout)::info_ptr
    integer,intent(in)::istep,ibatch,comm,rank,batch_size,n_ensemble,n_batches
    integer,intent(in)::batch_ranks(n_batches)
    real(kind=8),intent(in)::batch_state(batch_size,n_ensemble)
    type(io_info),pointer::info
    integer::ierr,imember,member_rank,offset,batch_length,batch_offset,local_io_index,ibatch_recv
    integer::req

    call c_f_pointer(info_ptr,info)

    ! Locate batch in the state array
    batch_offset=get_batch_offset(batch_size,ibatch)
    batch_length=get_batch_length(batch_size,ibatch,info%state_size)

    do imember=1,n_ensemble
       member_rank=info%io_ranks(imember)
       if(member_rank==rank) then
          local_io_index=get_local_io_index(n_ensemble,info%io_ranks, &
               rank,imember)

          info%local_io_data(batch_offset+1:batch_offset+batch_length, &
               local_io_index)=batch_state(1:batch_length,imember)
       else
          call mpi_isend(batch_state(1:batch_length,imember),batch_length,MPI_DOUBLE_PRECISION, &
               member_rank,ibatch,comm,req,ierr)
       end if
    end do

    call receive_results(info,istep,batch_size,batch_ranks,n_batches,n_ensemble,rank,comm)

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

  subroutine store_results(info_c_ptr,istep,batch_size,batch_ranks,n_batches,rank,comm,state_size)
    implicit none

    type(c_ptr),intent(inout)::info_c_ptr
    integer,intent(in)::istep,rank,comm,state_size,batch_size,n_batches
    integer,intent(in)::batch_ranks(n_batches)
    type(io_info),pointer::info
    integer::imember,local_io_counter

    call c_f_pointer(info_c_ptr,info)

    local_io_counter=1

    do while(count(info%batch_results_received)<info%local_io_size)
       call receive_results(info,istep,batch_size,batch_ranks,n_batches,info%n_ensemble,rank,comm)
    end do

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

    deallocate(info%observations)
    deallocate(info%obs_positions)
    deallocate(info%obs_errors)
    deallocate(info%batch_results_received)

    if(info%predictions_computed) deallocate(info%predictions)

    deallocate(info%io_ranks)
    deallocate(info%local_io_data)
    deallocate(info)

  end subroutine cleanup_interface

end module advect1d_assimilate_interfaces
