module advect1d_assimilate_interfaces
  use assimilation_model_interface, ONLY: base_model_interface
  use system_mpi
  use iso_c_binding
  use random, ONLY: random_normal

  use hdf5

  implicit none

  type, extends(base_model_interface)::advect1d_interface
     private
     real(kind=8),allocatable::observations(:),obs_errors(:),predictions(:,:)
     real(kind=8),pointer::local_io_data(:,:)
     real(kind=8)::cutoff,cutoff_u_a
     integer,allocatable::obs_positions(:),io_ranks(:)
     integer::n_observations,state_size,comm,local_io_size
     logical::observations_read,predictions_computed,state_loaded
   contains
     procedure::get_state_size
     procedure::get_subset_io_segment_data
     procedure::get_state_subset_buffer
     procedure::get_subset_obs_count
     procedure::get_subset_predictions
     procedure::get_subset_observations
     procedure::get_subset_obs_err
     procedure::get_weight_obs_obs
     procedure::get_weight_model_obs
     procedure::read_observations
     procedure::before_loading_ensemble_state
     procedure,private::compute_predictions
     procedure,private::load_ensemble_state
     procedure::after_ensemble_results_received
  end type advect1d_interface

contains

  function new_advect1d_interface(n_ensemble,n_observations,state_size,comm) result(this)

    integer(c_int),intent(in)::n_ensemble,n_observations,state_size,comm
    type(advect1d_interface)::this
    integer::ierr,rank,comm_size

    this%n_ensemble=n_ensemble
    this%n_observations=n_observations
    this%state_size=state_size
    this%comm=comm

    this%cutoff=0.1
    this%cutoff_u_a=0.2

    this%observations_read=.false.
    this%predictions_computed=.false.
    this%state_loaded=.false.

    call mpi_comm_rank(comm,rank,ierr)
    call mpi_comm_size(comm,comm_size,ierr)

    allocate(this%observations(this%n_observations))
    allocate(this%obs_errors(this%n_observations))
    allocate(this%obs_positions(this%n_observations))
    allocate(this%predictions(this%n_observations,this%n_ensemble))
    allocate(this%io_ranks(n_ensemble))

    ! Assign ensemble members to processors for i/o purposes
    call get_io_ranks(comm_size,n_ensemble,this%io_ranks)

    ! Get the number of ensemble members read and written locally
    this%local_io_size=get_rank_io_size(n_ensemble,this%io_ranks,rank)

    ! Allocate array for local i/o data
    allocate(this%local_io_data(state_size,this%local_io_size))

    call h5open_f(ierr)

  end function new_advect1d_interface

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

  function get_subset_obs_count(this,istep,subset_offset,subset_size) result(obs_count)

    implicit none

    class(advect1d_interface)::this
    integer,intent(in)::istep,subset_offset,subset_size
    integer::obs_count

    obs_count=this%n_observations

  end function get_subset_obs_count

  subroutine get_subset_predictions(this,istep,subset_offset,subset_size,predictions)

    use iso_c_binding

    class(advect1d_interface)::this
    integer,intent(in)::istep,subset_offset,subset_size
    real(kind=8),intent(inout)::predictions(:,:)
    integer::i,imember,ierr

    if(size(predictions,1) /= this%n_observations .or. &
       size(predictions,2) /= this%n_ensemble) then
       print '(A,I0,A,I0,A,I0,A,I0,A)', &
            'Wrong shape passed to predictions argument of get_subset_predictions. Expected (', &
            this%n_observations,',',this%n_ensemble, &
            '), got (',size(predictions,1),',',size(predictions,2),').'
       call mpi_abort(this%comm,1,ierr)
    end if

    if(.not. this%observations_read) then
       print *,'Error: Observations not yet read'
       error stop
    end if

    if(.not. this%predictions_computed) then
       print *,'Error: Predictions not yet computed'
       error stop
    end if

    predictions=this%predictions

  end subroutine get_subset_predictions

  subroutine get_subset_observations(this,istep,subset_offset,subset_size,observations)

    use iso_c_binding

    implicit none

    class(advect1d_interface)::this
    integer,intent(in)::istep,subset_offset,subset_size
    real(kind=8),intent(out)::observations(:)
    integer::ierr

    if(size(observations) /= this%n_observations) then
       print '(A,I0,A,I0,A,I0)', &
            'Wrong size array passed to observations argument of get_batch_observations. Expected size=', &
            this%n_observations,', got size=',size(observations)
       call mpi_abort(this%comm,1,ierr)
    end if

    if(.not. this%observations_read) then
       print *,'Error observations not yet read'
       error stop
    end if

    observations=this%observations

  end subroutine get_subset_observations

  SUBROUTINE get_subset_obs_err(this,istep,subset_offset,subset_size,obs_err)
    USE iso_c_binding
    class(advect1d_interface)::this
    integer,intent(in)::istep,subset_offset,subset_size
    REAL(c_double), INTENT(out) :: obs_err(:)

    if(.not. this%observations_read) then
       print *,'Error observations not yet read'
       error stop
    end if

    obs_err=this%obs_errors

  END SUBROUTINE get_subset_obs_err

  function get_weight_obs_obs(this,istep,subset_offset,subset_size,iobs1,iobs2) result(weight)

    use localization, ONLY: localize_gaspari_cohn

    class(advect1d_interface)::this
    integer,intent(in)::istep,subset_offset,subset_size,iobs1,iobs2
    real(kind=8)::weight
    real(kind=8)::pos1,pos2,delta,distance
    integer::domain_size

    if(.not. this%observations_read) then
       print *,'Error observations not yet read'
       error stop
    end if

    domain_size=this%state_size/2

    pos1=real(this%obs_positions(iobs1)-1)/domain_size
    pos2=real(this%obs_positions(iobs2)-1)/domain_size
    delta=abs(pos1-pos2)
    distance=min(delta,1-delta)

    weight=localize_gaspari_cohn(distance,this%cutoff)

  end function get_weight_obs_obs

  function get_weight_model_obs(this,istep,subset_offset,subset_size,imodel,iobs) result(weight)

    use localization, ONLY: localize_gaspari_cohn

    class(advect1d_interface)::this
    integer,intent(in)::istep,subset_offset,subset_size,imodel,iobs
    real(kind=8)::weight
    real(kind=8)::pos_obs,pos_model,delta,distance,cutoff
    integer::domain_size

    domain_size=this%state_size/2

    pos_obs=real(this%obs_positions(iobs)-1)/domain_size
    pos_model=real(mod(imodel+subset_offset-1,domain_size))/domain_size
    delta=abs(pos_obs-pos_model)
    distance=min(delta,1-delta)

    if(imodel<domain_size) then
       cutoff=this%cutoff
    else
       cutoff=this%cutoff_u_a
    end if

    weight=localize_gaspari_cohn(distance,cutoff)

  end function get_weight_model_obs

  subroutine read_observations(this,istep)
    class(advect1d_interface)::this
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

    this%n_observations=dims(1)

    ! Close the dataspace
    call h5sclose_f(dataspace,ierr)

    deallocate(this%observations)
    allocate(this%observations(this%n_observations))

    ! Read the data
    call h5dread_f(dset_h,H5T_NATIVE_DOUBLE,this%observations,dims,ierr)

    ! Close the dataset
    call h5dclose_f(dset_h,ierr)

    ! Open the obs_positions dataset
    call h5dopen_f(h5file_h,'obs_positions',dset_h,ierr)

    deallocate(this%obs_positions)
    allocate(this%obs_positions(this%n_observations))

    ! Read the data
    call h5dread_f(dset_h,H5T_NATIVE_INTEGER,this%obs_positions,dims,ierr)

    ! Close the dataset
    call h5dclose_f(dset_h,ierr)

    ! Open the obs_errors dataset
    call h5dopen_f(h5file_h,'obs_errors',dset_h,ierr)

    deallocate(this%obs_errors)
    allocate(this%obs_errors(this%n_observations))

    ! Read the data
    call h5dread_f(dset_h,H5T_NATIVE_DOUBLE,this%obs_errors,dims,ierr)

    ! Close the dataset
    call h5dclose_f(dset_h,ierr)

    ! Close the file
    call h5fclose_f(h5file_h,ierr)

    this%observations_read=.true.

  end subroutine read_observations

  subroutine load_observations_parallel(this,istep)
    type(advect1d_interface)::this
    integer,intent(in)::istep
    integer::ierr,rank

    call mpi_comm_rank(this%comm,rank,ierr)

    if(rank==0) call read_observations(this,istep)

    call mpi_bcast(this%n_observations,1,MPI_INTEGER,0,this%comm,ierr)

    if(rank>0) then
       deallocate(this%observations)
       allocate(this%observations(this%n_observations))
       deallocate(this%obs_positions)
       allocate(this%obs_positions(this%n_observations))
       deallocate(this%obs_errors)
       allocate(this%obs_errors(this%n_observations))
    end if

    call mpi_bcast(this%observations,this%n_observations,MPI_DOUBLE_PRECISION,0, &
         this%comm,ierr)
    call mpi_bcast(this%obs_positions,this%n_observations,MPI_INTEGER,0, &
         this%comm,ierr)
    call mpi_bcast(this%obs_errors,this%n_observations,MPI_DOUBLE_PRECISION,0, &
         this%comm,ierr)

    this%observations_read=.true.

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

  subroutine load_ensemble_state(this,istep)
    class(advect1d_interface)::this
    integer,intent(in)::istep
    integer::imember,local_io_counter,rank,ierr

    call mpi_comm_rank(this%comm,rank,ierr)

    call load_observations_parallel(this,istep)

    local_io_counter=1

    do imember=1,this%n_ensemble
       if(this%io_ranks(imember)==rank) then
          call read_state(istep,imember, &
               this%local_io_data(:,local_io_counter),this%state_size)
          local_io_counter=local_io_counter+1
       end if
    end do

    this%state_loaded=.true.

    call this%compute_predictions(istep)

  end subroutine load_ensemble_state

  subroutine compute_predictions(this,istep)
    class(advect1d_interface)::this
    integer,intent(in)::istep
    integer::imember,local_io_counter,rank,ierr,iobs
    real(kind=8)::member_predictions(this%n_observations)

    call mpi_comm_rank(this%comm,rank,ierr)

    if(.not. this%observations_read) then
       print *,'Error: Observations not yet read'
       error stop
    end if

    if(.not. this%state_loaded) then
       print *,'Error: Ensemble state not yet loaded'
       error stop
    end if

    local_io_counter=1

    do imember=1,this%n_ensemble
       if(this%io_ranks(imember)==rank) then

          ! Compute predictions for this ensemble member
          do iobs=1,this%n_observations
             member_predictions(iobs)=this%local_io_data( &
                  this%obs_positions(iobs)+1,local_io_counter)
          end do

          local_io_counter=local_io_counter+1
       end if

       ! Broadcast to all processors
       call mpi_bcast(member_predictions,this%n_observations, &
            MPI_DOUBLE_PRECISION,this%io_ranks(imember),this%comm,ierr)

       this%predictions(:,imember)=member_predictions

    end do

    this%predictions_computed=.true.
  end subroutine compute_predictions

  subroutine before_loading_ensemble_state(this,istep)
    class(advect1d_interface)::this
    integer,intent(in)::istep
    integer::imember,local_io_counter,rank,ierr,iobs

    call mpi_comm_rank(this%comm,rank,ierr)

    call this%read_observations(istep)

    call this%load_ensemble_state(istep)

    call this%compute_predictions(istep)

  end subroutine before_loading_ensemble_state

  subroutine get_subset_io_segment_data(this,istep,imember,subset_offset,subset_size,counts,offsets)
    class(advect1d_interface)::this
    integer,intent(in)::istep,imember,subset_offset,subset_size
    integer,intent(out)::counts(:),offsets(:)
    integer::comm_size,ierr

    call mpi_comm_size(this%comm,comm_size,ierr)

    if(size(counts)/=comm_size) then
       print '(A,I0,A,I0,A)','Wrong size passed to counts argument of get_subset_io_segment_data. Expected ', &
            comm_size,', got ',size(counts),'.'
       error stop
    end if

    if(size(offsets)/=comm_size) then
       print '(A,I0,A,I0,A)','Wrong size passed to counts argument of get_subset_io_segment_data. Expected ', &
            comm_size,', got ',size(offsets),'.'
       error stop
    end if

    ! Fill counts with zeros
    counts=0
    offsets=0

    ! Set counts for the rank assigned to do i/o for requested ensemble member
    counts(this%io_ranks(imember)+1)=subset_size
    offsets(this%io_ranks(imember)+1)=0

  end subroutine get_subset_io_segment_data

  function get_state_subset_buffer(this,istep,imember,subset_offset,subset_size) result(buffer)

    class(advect1d_interface)::this
    integer,intent(in)::istep,imember,subset_offset,subset_size
    real(kind=8),pointer::buffer(:)
    real(kind=8),target::empty(0)
    integer::rank,ierr

    call mpi_comm_rank(this%comm,rank,ierr)

    if(get_local_io_index(this%n_ensemble,this%io_ranks,rank,imember)>0) then

       ! Point to the appropriate position in local_io_data
       buffer=>this%local_io_data( &
            subset_offset+1:subset_offset+subset_size, &
            get_local_io_index(this%n_ensemble,this%io_ranks,rank,imember))

    else
       ! Return a null pointer
       buffer=>empty
    end if

  end function get_state_subset_buffer

  function get_state_size(this) result(size)
    class(advect1d_interface)::this
    integer::size
    size=this%state_size
  end function get_state_size

  subroutine after_ensemble_results_received(this,istep)
    class(advect1d_interface)::this
    integer,intent(in)::istep
    integer::imember,rank,ierr,imember_local

    call MPI_Comm_rank(this%comm,rank,ierr)

    imember_local=1

    do imember=1,this%n_ensemble

       if(this%io_ranks(imember)==rank) then

          call write_state(this,istep,imember,this%n_ensemble,this%comm,this%local_io_data(:,imember_local),this%state_size)

          imember_local=imember_local+1

       end if

    end do

  end subroutine after_ensemble_results_received

  subroutine write_state(this,istep,imember,n_ensemble,comm,member_state,state_size)

    type(advect1d_interface)::this
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

end module advect1d_assimilate_interfaces
