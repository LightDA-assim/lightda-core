#include "mpi_types.h"

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
     integer::n_ensemble,state_size,n_observations,n_batches,batch_size,n_local_batches
     MPI_COMM_TYPE::comm
     logical,allocatable::batch_results_received(:,:)
   contains
     procedure::load_ensemble_state
     procedure::receive_results
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
    integer(c_int),intent(in)::n_ensemble,state_size,batch_size
    MPI_COMM_TYPE::comm
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

    new_batch_manager%batch_results_received=.false.

    ! Assign batches to process ranks
    call get_batch_ranks(comm_size,new_batch_manager%batch_ranks)

    ! Get number of local batches
    new_batch_manager%n_local_batches=new_batch_manager%get_rank_batch_count(rank)

  end function new_batch_manager

  function get_comm(this) result(comm)
    class(assim_batch_manager)::this
    MPI_COMM_TYPE::comm

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
       call mpi_abort(this%comm,1,ierr)
    end if

    ibatch_rank=1

    do ibatch=1,this%n_batches
       if(this%batch_ranks(ibatch)==rank) then
          batches(ibatch_rank)=ibatch
          ibatch_rank=ibatch_rank+1
       end if
    end do

  end subroutine get_rank_batches

  subroutine segment_range_overlap(segment1_start,segment1_end,segment2_start,segment2_end,overlap_start,overlap_end)
    integer,intent(in)::segment1_start,segment1_end,segment2_start,segment2_end
    integer,intent(out)::overlap_start,overlap_end

    overlap_start=max(segment1_start,segment2_start)
    overlap_end=min(segment1_end,segment2_end)
  end subroutine segment_range_overlap

  subroutine load_ensemble_state(this,istep,local_batches)

    use util, ONLY: append_array

    class(assim_batch_manager)::this

    integer,intent(in)::istep
    real(kind=8),intent(out),target::local_batches(this%batch_size,this%n_local_batches,this%n_ensemble)
    real(kind=8),allocatable::sendbuf(:)
    real(kind=8),pointer::recvbuf(:)
    real(kind=8)::received_batch_state(this%batch_size)
    integer::imember,ierr,ibatch,ibatch_local,intercomm,comm_size, &
         batch_length,batch_offset,iobs,rank, &
         irank
#ifdef HAVE_MPI_F08_MODULE
    MPI_STATUS_TYPE,allocatable::statuses(:)
#else
    MPI_STATUS_TYPE,allocatable::statuses(:,:)
#endif
    MPI_REQUEST_TYPE,allocatable::receive_reqs(:),receive_req
    integer,allocatable::io_ranks(:),io_counts(:),io_offsets(:)
    integer::batch_io_offsets,batch_segment_start,batch_segment_end, &
         overlap_start,overlap_end,local_batch_length,ireq

    call mpi_comm_size(this%comm,comm_size,ierr)

    call mpi_comm_rank(this%comm,rank,ierr)

    call this%model_interface%read_state(istep)

    batch_io_offsets=0

    ! Initialize receive_reqs
    allocate(receive_reqs(0))

    do imember=1,this%n_ensemble

       ! Reset the local batch counter
       ibatch_local=1

       ! Get locations of data (i/o segments) for this ensemble member
       call this%model_interface%get_io_ranks(istep, &
            imember,io_ranks,io_counts,io_offsets)

       do ibatch=1,this%n_batches

          ! Locate batch in the state array
          batch_offset=this%get_batch_offset(ibatch)
          batch_length=this%get_batch_length(ibatch)

          ! Loop over i/o segments
          do irank=1,size(io_ranks)

             call segment_range_overlap( &
                  batch_offset+1,batch_offset+batch_length, &
                  io_offsets(irank)+1,io_offsets(irank)+io_counts(irank), &
                  overlap_start,overlap_end)

             if(overlap_start>overlap_end) then
                ! I/O segment does not overlap with batch
                cycle
             end if

             if(io_ranks(irank)==rank) then
                ! This segment will be read locally, no receive required
                cycle

             end if

             if(this%batch_ranks(ibatch)==rank) then

                ! Location of data in the batch
                batch_segment_start=overlap_start-batch_offset
                batch_segment_end=overlap_end-batch_offset

                ! Size of data to be read
                local_batch_length=batch_segment_end-batch_segment_start+1

                receive_reqs=append_array(receive_reqs,MPI_REQUEST_NULL)

                recvbuf=>local_batches(batch_segment_start: &
                     batch_segment_end,ibatch_local,imember)

                ! Receive this segment from its i/o processor
                call MPI_Irecv(recvbuf, &
                     local_batch_length, MPI_DOUBLE_PRECISION, &
                     io_ranks(irank),1,this%comm,receive_reqs(size(receive_reqs)),ierr)

             end if

          end do

          if(this%batch_ranks(ibatch)==rank) then
             ! Increment the local batch counter
             ibatch_local=ibatch_local+1
          end if

       end do

    end do

       if(ibatch_local-1/=this%n_local_batches) then
          print '(A,I0,A,I0,A,I0,A)', &
               'Inconsistency in the local batch count on rank ',rank, &
               '. Expected ', this%n_local_batches, &
               ', got ',ibatch_local-1,'.'
          error stop
       end if

    do imember=1,this%n_ensemble

       ! Reset the local batch counter
       ibatch_local=1

       ! Get locations of data (i/o segments) for this ensemble member
       call this%model_interface%get_io_ranks(istep, &
            imember,io_ranks,io_counts,io_offsets)

       do ibatch=1,this%n_batches

          ! Locate batch in the state array
          batch_offset=this%get_batch_offset(ibatch)
          batch_length=this%get_batch_length(ibatch)

          ! Loop over i/o segments
          do irank=1,size(io_ranks)

             call segment_range_overlap( &
                  batch_offset+1,batch_offset+batch_length, &
                  io_offsets(irank)+1,io_offsets(irank)+io_counts(irank), &
                  overlap_start,overlap_end)

             if(overlap_start>overlap_end) then
                ! I/O segment does not overlap with batch
                cycle
             end if

             ! Location of data in the model state array
             batch_segment_start=overlap_start-batch_offset
             batch_segment_end=overlap_end-batch_offset

             ! Size of data to be read
             local_batch_length=batch_segment_end-batch_segment_start+1

             if(io_ranks(irank)==rank) then

                ! This segment will be read locally, we will request it from
                ! the model interface

                ! Allocate buffer for locally read data in this batch
                allocate(sendbuf(local_batch_length))

                ! Get the data from the model
                sendbuf=this%model_interface%get_state_subset(istep,imember, &
                     overlap_start-1,local_batch_length)

                if(this%batch_ranks(ibatch)==rank) then

                   ! This segment will be assimilated locally, store it in
                   ! the local_batches array

                   local_batches(batch_segment_start: &
                        batch_segment_end,ibatch_local,imember)=sendbuf

                else

                   ! Send data to the process where it will be assimilated
                   ! Note that we do a blocking send because the send needs
                   ! to complete before we deallocate sendbuf
                   call MPI_Send(sendbuf,local_batch_length, &
                        MPI_DOUBLE_PRECISION,this%batch_ranks(ibatch),1, &
                        this%comm)

                end if

                deallocate(sendbuf)

             end if

          end do

          if(this%batch_ranks(ibatch)==rank) then
             ! Increment the local batch counter
             ibatch_local=ibatch_local+1
          end if

       end do

       if(ibatch_local-1/=this%n_local_batches) then
          print '(A,I0,A,I0,A,I0,A)', &
               'Inconsistency in the local batch count on rank ',rank, &
               '. Expected ', this%n_local_batches, &
               ', got ',ibatch_local-1,'.'
          error stop
       end if

    end do

#ifdef HAVE_MPI_F08_MODULE
    allocate(statuses(size(receive_reqs)))
#else
    allocate(statuses(MPI_STATUS_SIZE,size(receive_reqs)))
#endif

    ! Wait for receives to complete
    call MPI_Waitall(size(receive_reqs),receive_reqs,statuses,ierr)

  end subroutine load_ensemble_state

  subroutine receive_results(this,istep,local_batches)

    use util, ONLY: append_array

    class(assim_batch_manager)::this
    integer,intent(in)::istep
    real(kind=8),intent(in),target::local_batches(this%batch_size,this%n_local_batches,this%n_ensemble)
    integer::imember,ibatch,n_batches,batch_rank,batch_offset,batch_length, &
         member_rank,local_io_index,ierr
#ifdef HAVE_MPI_F08_MODULE
    MPI_STATUS_TYPE,allocatable::statuses(:)
    MPI_STATUS_TYPE::status
#else
    MPI_STATUS_TYPE,allocatable::statuses(:,:)
    MPI_STATUS_TYPE::status
#endif
    integer::rank,comm_size,ireq,completed_req_count,req_ind,ibatch_local
    integer::completed_req_inds(this%n_ensemble*this%n_batches)
    integer::batch_segment_start,batch_segment_end
    integer::irank
    integer::local_batch_length,overlap_start,overlap_end
    integer,allocatable::io_ranks(:),io_counts(:),io_offsets(:)
    real(kind=8),pointer::sendbuf(:)
    real(kind=8),allocatable::recvbuf(:)
    real(kind=8),target::empty(0)
    MPI_REQUEST_TYPE,allocatable::req
    MPI_REQUEST_TYPE,allocatable::reqs(:)
    logical::flag

    call mpi_comm_size(this%comm,comm_size,ierr)

    call mpi_comm_rank(this%comm,rank,ierr)

    sendbuf=>empty

    allocate(reqs(0))

    do imember=1,this%n_ensemble

       ! Reset local batch counter
       ibatch_local=1

       ! Get locations of data (i/o segments) for this ensemble member
       call this%model_interface%get_io_ranks(istep, &
            imember,io_ranks,io_counts,io_offsets)

       do ibatch=1,this%n_batches

          if(this%batch_ranks(ibatch)==rank) then

             ! Locate batch in the state array
             batch_offset=this%get_batch_offset(ibatch)
             batch_length=this%get_batch_length(ibatch)

             ! Loop over i/o segments
             do irank=1,size(io_ranks)

                call segment_range_overlap( &
                     batch_offset+1,batch_offset+batch_length, &
                     io_offsets(irank)+1,io_offsets(irank)+io_counts(irank), &
                     overlap_start,overlap_end)

                if(overlap_start>overlap_end) then
                   ! I/O segment does not overlap with batch
                   cycle
                end if

                ! Location of data in the model state array
                batch_segment_start=overlap_start-batch_offset
                batch_segment_end=overlap_end-batch_offset

                ! Size of data to be read
                local_batch_length=batch_segment_end-batch_segment_start+1

                sendbuf=>local_batches(batch_segment_start:batch_segment_end, &
                     ibatch_local,imember)

                if(io_ranks(irank)==rank) then

                   call this%model_interface%set_state_subset( &
                        istep,imember,overlap_start-1, &
                        local_batch_length,sendbuf)

                else

                   reqs=append_array(reqs,MPI_REQUEST_NULL)

                   call MPI_Isend(sendbuf,local_batch_length, &
                        MPI_DOUBLE_PRECISION,io_ranks(irank),1,this%comm, &
                        reqs(size(reqs)),ierr)

                end if

                ibatch_local=ibatch_local+1

             end do

          end if

       end do

    end do

    do imember=1,this%n_ensemble

       ! Reset local batch counter
       ibatch_local=1

       ! Get locations of data (i/o segments) for this ensemble member
       call this%model_interface%get_io_ranks(istep, &
            imember,io_ranks,io_counts,io_offsets)

       do ibatch=1,this%n_batches

          if(this%batch_ranks(ibatch)/=rank) then

             ! Locate batch in the state array
             batch_offset=this%get_batch_offset(ibatch)
             batch_length=this%get_batch_length(ibatch)

             ! Loop over i/o segments
             do irank=1,size(io_ranks)

                call segment_range_overlap( &
                     batch_offset+1,batch_offset+batch_length, &
                     io_offsets(irank)+1,io_offsets(irank)+io_counts(irank), &
                     overlap_start,overlap_end)

                if(overlap_start>overlap_end) then
                   ! I/O segment does not overlap with batch
                   cycle
                end if

                ! Size of data to be read
                local_batch_length=overlap_end-overlap_start+1

                if(io_ranks(irank)==rank) then

                   allocate(recvbuf(local_batch_length))

                   call MPI_Recv(recvbuf,local_batch_length, &
                        MPI_DOUBLE_PRECISION,this%batch_ranks(ibatch),1,this%comm, &
                        status,ierr)

                   call this%model_interface%set_state_subset( &
                        istep,imember,overlap_start-1, &
                        local_batch_length,recvbuf)

                   deallocate(recvbuf)

                end if

             end do

          end if

       end do

    end do

#ifdef HAVE_MPI_F08_MODULE
    allocate(statuses(size(reqs)))
#else
    allocate(statuses(MPI_STATUS_SIZE,size(reqs)))
#endif

    call MPI_Waitall(size(reqs),reqs,statuses,ierr)

    this%batch_results_received=.true.

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
                  'from rank',this%batch_ranks(ibatch)
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

  subroutine store_results(this,istep,local_batches)

    implicit none

    class(assim_batch_manager)::this
    integer,intent(in)::istep
    real(kind=8),intent(in)::local_batches(this%batch_size,this%n_local_batches,this%n_ensemble)
    integer::imember,rank,ierr

    do while(any(this%batch_results_received.eqv..false.))
       call this%receive_results(istep,local_batches)
    end do

    call this%model_interface%write_state(istep)

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
#ifdef HAVE_MPI_F08_MODULE
    MPI_STATUS_TYPE,allocatable::status(:)
#else
    MPI_STATUS_TYPE,allocatable::status(:,:)
#endif
    integer::ierr

  end subroutine cleanup

end module assimilation_batch_manager
