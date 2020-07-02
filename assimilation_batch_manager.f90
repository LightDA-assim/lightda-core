module assimilation_batch_manager
  use system_mpi
  use iso_c_binding
  use random, ONLY: random_normal
  use random_integer, ONLY: randint
  use assimilation_model_interface

  implicit none

  type :: assim_batch_manager
     private
     class(base_model_interface),pointer,public::model_interface
     integer,allocatable::batch_ranks(:)
     integer::comm,n_ensemble,state_size,n_observations,n_batches,batch_size,n_local_batches
     logical,allocatable::batch_results_received(:,:)
     integer,allocatable::batch_receive_reqs(:,:),batch_send_reqs(:,:)
   contains
     procedure::load_ensemble_state
     procedure::receive_results
     procedure::transmit_results
     procedure::get_comm
     procedure::get_batch_size
     procedure::get_state_size
     procedure::get_n_batches
     procedure::get_n_ensemble
     procedure::add_obs_err
     procedure::localize
     procedure::store_results
     procedure::get_rank_batch_count
     procedure::get_rank_batches
     procedure::get_batch_offset
     procedure::get_batch_length
     final::cleanup
  end type assim_batch_manager

contains

  function new_batch_manager(model_interface,n_ensemble,state_size,batch_size,comm)
    integer(c_int),intent(in)::n_ensemble,state_size,batch_size,comm
    class(base_model_interface),intent(inout),target::model_interface
    type(assim_batch_manager)::new_batch_manager
    integer::ierr,rank,comm_size

    ! Initialize state info
    new_batch_manager%comm=comm
    new_batch_manager%model_interface=>model_interface
    new_batch_manager%n_ensemble=n_ensemble
    new_batch_manager%state_size=state_size
    new_batch_manager%n_observations=0
    new_batch_manager%batch_size=batch_size

    call mpi_comm_size(comm,comm_size,ierr)
    call mpi_comm_rank(comm,rank,ierr)

    new_batch_manager%n_batches=get_batch_count(state_size,batch_size)

    allocate(new_batch_manager%batch_ranks(new_batch_manager%n_batches))
    allocate(new_batch_manager%batch_results_received(new_batch_manager%n_ensemble,new_batch_manager%n_batches))
    allocate(new_batch_manager%batch_receive_reqs(new_batch_manager%n_ensemble,new_batch_manager%n_batches))
    allocate(new_batch_manager%batch_send_reqs(new_batch_manager%n_ensemble,new_batch_manager%n_batches))

    new_batch_manager%batch_receive_reqs=MPI_REQUEST_NULL

    new_batch_manager%batch_results_received=.false.

    ! Assign batches to process ranks
    call get_batch_ranks(comm_size,new_batch_manager%batch_ranks)

    ! Get number of local batches
    new_batch_manager%n_local_batches=new_batch_manager%get_rank_batch_count(rank)

  end function new_batch_manager

  function get_comm(this) result(comm)
    class(assim_batch_manager)::this
    integer::comm

    comm=this%comm
    
  end function get_comm

  function get_batch_size(this) result(batch_size)
    class(assim_batch_manager)::this
    integer::batch_size

    batch_size=this%batch_size
    
  end function get_batch_size

  function get_state_size(this) result(state_size)
    class(assim_batch_manager)::this
    integer::state_size

    state_size=this%state_size
    
  end function get_state_size

  function get_n_batches(this) result(n_batches)
    class(assim_batch_manager)::this
    integer::n_batches

    n_batches=this%n_batches
    
  end function get_n_batches

  function get_n_ensemble(this) result(n_ensemble)
    class(assim_batch_manager)::this
    integer::n_ensemble

    n_ensemble=this%n_ensemble
    
  end function get_n_ensemble

  subroutine shuffle(a)
    integer, intent(inout) :: a(:)
    integer :: i, randpos, temp
 
    do i = size(a), 2, -1
       randpos=randint(i)
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
    rank=comm_size

    do ibatch=1,size(batch_ranks)
       ! Decrement rank
       rank=rank-1

       ! raise to comm_size-1 if we reach zero
       if(rank<0) rank=comm_size-1

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

  function get_batch_offset(this,ibatch) result(offset)

    class(assim_batch_manager)::this
    integer,intent(in)::ibatch
    integer::offset
    offset=this%batch_size*(ibatch-1)
  end function get_batch_offset

  function get_batch_length(this,ibatch) result(length)

    class(assim_batch_manager)::this
    integer,intent(in)::ibatch
    integer::offset,length
    offset=this%get_batch_offset(ibatch)
    length=min(this%state_size,offset+this%batch_size)-offset
  end function get_batch_length

  function get_rank_batch_count(this,rank) result(count)

    class(assim_batch_manager)::this
    integer,intent(in)::rank
    integer::ibatch,count
    count=0

    do ibatch=1,this%n_batches
       if(this%batch_ranks(ibatch)==rank) count=count+1
    end do

  end function get_rank_batch_count

  subroutine get_rank_batches(this,rank,batches)
    class(assim_batch_manager)::this
    integer,intent(in)::rank
    integer,intent(out)::batches(:)
    integer::n_batches
    integer::ibatch_rank,ibatch,ierr

    n_batches=this%get_rank_batch_count(rank)

    if(size(batches,1)/=n_batches) then
       print '(A,I0,A,I0)','Wrong array size passed to get_rank_batches. Expected ',n_batches,' got ',size(batches,1)
       error stop
       call mpi_abort(this%comm,ierr)
    end if

    ibatch_rank=1

    do ibatch=1,this%n_batches
       if(this%batch_ranks(ibatch)==rank) then
          batches(ibatch_rank)=ibatch
          ibatch_rank=ibatch_rank+1
       end if
    end do

  end subroutine get_rank_batches

  subroutine load_ensemble_state(this,istep,local_batches)

    class(assim_batch_manager)::this

    integer,intent(in)::istep
    real(kind=8),intent(out)::local_batches(this%batch_size,this%n_local_batches,this%n_ensemble)
    real(kind=8),pointer::sendbuf(:)
    real(kind=8)::received_batch_state(this%batch_size)
    integer::imember,ierr,ibatch,ibatch_local,intercomm,comm_size, &
         batch_length,batch_offset,local_io_counter,iobs,rank,local_batch_length
    integer::status(MPI_STATUS_SIZE)
    integer,allocatable::batch_io_counts(:),batch_io_offsets(:)

    call this%model_interface%before_loading_ensemble_state(istep)

    call mpi_comm_size(this%comm,comm_size,ierr)

    call mpi_comm_rank(this%comm,rank,ierr)

    allocate(batch_io_counts(comm_size))
    allocate(batch_io_offsets(comm_size))

    local_io_counter=1

    batch_io_offsets=0

    do imember=1,this%n_ensemble

       local_io_counter=local_io_counter+1

       ! Reset the local batch counter
       ibatch_local=1

       do ibatch=1,this%n_batches

          ! Locate batch in the state array
          batch_offset=this%get_batch_offset(ibatch)
          batch_length=this%get_batch_length(ibatch)

          ! Get number of state array members read by each process
          call this%model_interface%get_subset_io_segment_data(istep, &
               imember,batch_offset,batch_length,batch_io_counts,batch_io_offsets)

          ! Part of the state array read by local process
          local_batch_length=batch_io_counts(rank+1)

          ! Get the batch data that is to be read by the local processor
          sendbuf=>this%model_interface%get_state_subset_buffer( &
               istep,imember,batch_offset,batch_length)

          ! Gather batch data to the assigned processor from batch_ranks
          call MPI_Igatherv(sendbuf,local_batch_length, &
               MPI_DOUBLE_PRECISION,received_batch_state,batch_io_counts, &
               batch_io_offsets, MPI_DOUBLE_PRECISION, &
               this%batch_ranks(ibatch), &
               this%comm, this%batch_send_reqs(imember,ibatch), ierr)

          if(this%batch_ranks(ibatch)==rank) then

             ! Wait for MPI_Igatherv to complete
             call MPI_Wait(this%batch_send_reqs(imember,ibatch),status,ierr)

             ! Copy batch state into the local_batches array
             local_batches(:,ibatch_local,imember)=received_batch_state

             ! Increment the local batch counter
             ibatch_local=ibatch_local+1

          end if

       end do

       if(ibatch_local-1/=this%n_local_batches) then
          print '(A,I0,A,I0)', &
               'Inconsistency in the local batch count. Expected ', &
               this%n_local_batches,' found ',ibatch_local-1
          error stop
       end if

    end do

  end subroutine load_ensemble_state

  subroutine scatter_batch(this,istep,imember,ibatch,sendbuf)
    class(assim_batch_manager)::this
    integer,intent(in)::istep,ibatch,imember
    real(kind=8),pointer::recvbuf(:)
    real(kind=8),intent(in)::sendbuf(:)
    integer,allocatable::batch_io_counts(:),batch_io_offsets(:)
    integer::batch_offset,batch_length,comm_size,rank,ierr

    call mpi_comm_size(this%comm,comm_size,ierr)

    call mpi_comm_rank(this%comm,rank,ierr)

    ! Locate batch in the state array
    batch_offset=this%get_batch_offset(ibatch)
    batch_length=this%get_batch_length(ibatch)

    if(rank==this%batch_ranks(ibatch)) then
       if(size(sendbuf)/=batch_length ) then
          print '(A,I0,A,I0,A)','Wrong buffer size passed to scatter_batch. Expected ', &
               batch_length,', got ',size(sendbuf),'.'
          error stop
       end if
    end if

    allocate(batch_io_counts(comm_size))
    allocate(batch_io_offsets(comm_size))

    ! Get number of state array members to be written by each process
    call this%model_interface%get_subset_io_segment_data(istep, &
         imember,batch_offset,batch_length,batch_io_counts,batch_io_offsets)

    ! Get receive buffer
    recvbuf=>this%model_interface%get_state_subset_buffer(istep,imember, &
         batch_offset,batch_length)

    if(batch_io_counts(rank+1)>0) then
       if(.not. associated(recvbuf)) then
          print *,'Unassociated receive buffer where a buffer of size',batch_io_counts(rank+1),'was expected.'
          error stop
       else if(size(recvbuf)<batch_io_counts(rank+1)) then
          print *,'Expected receive buffer of size',batch_io_counts(rank+1),'got size',size(recvbuf)
          error stop
       end if
    end if

    call MPI_Iscatterv(sendbuf,batch_io_counts,batch_io_offsets, &
         MPI_DOUBLE_PRECISION,recvbuf,batch_io_counts(rank+1), &
         MPI_DOUBLE_PRECISION,this%batch_ranks(ibatch),this%comm, &
         this%batch_receive_reqs(imember,ibatch),ierr)

  end subroutine scatter_batch

  subroutine receive_results(this,istep,local_batches)

    class(assim_batch_manager)::this
    integer,intent(in)::istep
    real(kind=8),intent(in),target::local_batches(this%batch_size,this%n_local_batches,this%n_ensemble)
    integer::imember,ibatch,n_batches,batch_rank,batch_offset,batch_length, &
         member_rank,local_io_index,status(MPI_STATUS_SIZE),ierr
    integer::rank,comm_size,ireq,completed_req_count,req_ind,ibatch_local
    integer::completed_req_inds(this%n_ensemble*this%n_batches)
    integer,allocatable::batch_io_counts(:),batch_io_offsets(:)
    integer::statuses(MPI_STATUS_SIZE,this%n_ensemble*this%n_batches)
    real(kind=8),pointer::recvbuf(:)
    real(kind=8),pointer::sendbuf(:)
    real(kind=8),target::empty(0)
    logical::flag

    call mpi_comm_size(this%comm,comm_size,ierr)

    call mpi_comm_rank(this%comm,rank,ierr)

    allocate(batch_io_counts(comm_size))
    allocate(batch_io_offsets(comm_size))

    sendbuf=>empty

    do imember=1,this%n_ensemble
       ibatch_local=1
       do ibatch=1,this%n_batches


          ! Locate batch in the state array
          batch_offset=this%get_batch_offset(ibatch)
          batch_length=this%get_batch_length(ibatch)

          ! Get number of state array members to be written by each process
          call this%model_interface%get_subset_io_segment_data(istep, &
               imember,batch_offset,batch_length,batch_io_counts,batch_io_offsets)

          ! Get receive buffer
          recvbuf=>this%model_interface%get_state_subset_buffer(istep,imember, &
               batch_offset,batch_length)

          ! Set send buffer
          if(this%batch_ranks(ibatch)==rank) then
             sendbuf=>local_batches(:batch_length,ibatch_local,imember)
             ibatch_local=ibatch_local+1
          else
             sendbuf=>empty
          end if

          ! Get number of state array members to be written by each process
          call this%model_interface%get_subset_io_segment_data(istep, &
               imember,batch_offset,batch_length,batch_io_counts,batch_io_offsets)

          call scatter_batch(this,istep,imember,ibatch,sendbuf)

          if(batch_io_counts(rank+1)>0) then
             call MPI_Wait(this%batch_receive_reqs(imember,ibatch),status,ierr)

             ! Tell the model interface that the batch has been received
             call this%model_interface%after_member_state_received(istep,imember,batch_offset,batch_length)

          end if

          ! Record that batch was received
          this%batch_results_received(imember,ibatch)=.true.

       end do

    end do

  end subroutine receive_results

  subroutine print_remaining_batches(this)
    class(assim_batch_manager)::this
    integer::imember,ibatch,rank,ierr

    call mpi_comm_rank(this%comm,rank,ierr)

    print *,'Batches still pending:'

    do imember=1,this%n_ensemble
       do ibatch=1,this%n_batches
          if(.not.this%batch_results_received(imember,ibatch)) then
             print *,'member',imember,'batch',ibatch,&
                  'from rank',this%batch_ranks(ibatch),&
                  'request',this%batch_receive_reqs(imember,ibatch)
          end if
       end do
    end do
  end subroutine print_remaining_batches

  subroutine request_batch_index(this,req_ind,imember,ibatch)
    class(assim_batch_manager)::this
    integer,intent(in)::req_ind
    integer,intent(out)::imember,ibatch

    imember=mod(req_ind,this%n_ensemble)+1
    ibatch=(req_ind/this%n_ensemble)+1
  end subroutine request_batch_index

  subroutine transmit_results(this,istep,ibatch,batch_state)
    class(assim_batch_manager),intent(inout)::this
    integer,intent(in)::istep,ibatch
    real(kind=8),intent(in)::batch_state(this%batch_size,this%n_ensemble)
    integer::ierr,imember,member_rank,offset,batch_length,batch_offset, &
         comm_size,rank
    integer::req

    do imember=1,this%n_ensemble
       call scatter_batch(this,istep,imember,ibatch,batch_state(:,imember))
    end do

    call this%receive_results(istep,batch_state)

  end subroutine transmit_results

  subroutine store_results(this,istep,local_batches)

    implicit none

    class(assim_batch_manager)::this
    integer,intent(in)::istep
    real(kind=8),intent(in)::local_batches(this%batch_size,this%n_local_batches,this%n_ensemble)
    integer::imember,local_io_counter,rank,ierr
    integer,allocatable::status(:,:)

    call mpi_comm_rank(this%comm,rank,ierr)

    allocate(status(MPI_STATUS_SIZE,size(this%batch_send_reqs)))

    call MPI_Waitall(size(this%batch_send_reqs),this%batch_send_reqs,status,ierr)

    local_io_counter=1

    do while(any(this%batch_results_received.eqv..false.))
       call this%receive_results(istep,local_batches)
    end do

    call this%model_interface%after_ensemble_results_received(istep)

  end subroutine store_results

  SUBROUTINE add_obs_err(this,istep,ibatch,dim_obs,HPH)
    ! Add observation error covariance matrix
    class(assim_batch_manager)::this
    INTEGER(c_int32_t), INTENT(in), value :: istep, ibatch, dim_obs
    REAL(c_double), INTENT(inout) :: HPH(dim_obs,dim_obs)
    real(kind=8)::obs_err(dim_obs)
    integer::iobs,batch_offset,batch_length

    ! Locate batch in the state array
    batch_offset=this%get_batch_offset(ibatch)
    batch_length=this%get_batch_length(ibatch)

    call this%model_interface%get_subset_obs_err(istep,batch_offset, &
         batch_length,obs_err)

    do iobs=1,dim_obs
       HPH(iobs,iobs)=HPH(iobs,iobs)+obs_err(iobs)**2
    end do

  END SUBROUTINE add_obs_err

  SUBROUTINE localize(this,istep,ibatch,dim_p,dim_obs,HP_p,HPH)
    ! Apply localization to HP and HPH^T
    USE iso_c_binding
    class(assim_batch_manager)::this
    INTEGER(c_int32_t), INTENT(in), value :: istep, ibatch, dim_p, dim_obs
    REAL(c_double), INTENT(inout) :: HP_p(dim_obs,dim_p), HPH(dim_obs,dim_obs)
    real(kind=8)::cutoff,cutoff_u_a,pos,pos_obs1,pos_obs2,pos_obs,c,distance,delta,w
    integer::domain_size,iobs1,iobs2,ipos,batch_offset,batch_length

    cutoff=0.1
    cutoff_u_a=0.2

    ! Locate batch in the state array
    batch_offset=this%get_batch_offset(ibatch)
    batch_length=this%get_batch_length(ibatch)

    if(dim_p/=batch_length) then
       print *,'Inconsistent batch size'
       stop
    end if

    do iobs1=1,dim_obs

       do iobs2=1,dim_obs
          w=this%model_interface%get_weight_obs_obs(istep,iobs1,iobs2)

          HPH(iobs1,iobs2)=HPH(iobs1,iobs2)*w

       end do

       do ipos=1,dim_p

          w=this%model_interface%get_weight_model_obs(istep,ipos+batch_offset,iobs1)

          HP_p(iobs1,ipos)=HP_p(iobs1,ipos)*w

       end do
    end do

  END SUBROUTINE localize

  subroutine cleanup(this)
    type(assim_batch_manager)::this
    integer,allocatable::status(:,:)
    integer::ierr

    allocate(status(MPI_STATUS_SIZE,size(this%batch_receive_reqs)))

    call MPI_Waitall(size(this%batch_receive_reqs),this%batch_receive_reqs,status,ierr)

  end subroutine cleanup

end module assimilation_batch_manager
