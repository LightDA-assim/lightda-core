#include "mpi_types.h"

module assimilation_batch_manager
  use system_mpi
  use iso_c_binding
  use random, ONLY: random_normal
  use random_integer, ONLY: randint
  use assimilation_model_interface
  use distributed_array, ONLY: darray, darray_segment, new_darray

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
     procedure::get_batches_darray
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

  function get_batches_darray(this) result(batches_darray)

    !! Distributed array with the size of the model state and 

    class(assim_batch_manager)::this
        !! Batch manager

    type(darray)::batches_darray
        !! darray of model state with segments aligning to batch ranks

    type(darray_segment),target::batch_segments(this%n_batches)
        ! darray segments that make up batches_darray

    type(darray_segment),pointer::batch_segment
        ! Pointer to a member of batch_segments

    real(kind=8)::empty_batch(this%batch_size)
        ! Array of zeros to populate the array

    integer::rank ! MPI rank
    integer::ierr ! MPI status code

    integer::ibatch ! Loop counter

    call mpi_comm_rank(this%comm,rank,ierr)

    empty_batch=0

    ! Build darray segments for each batch
    do ibatch=1,this%n_batches
       batch_segment=>batch_segments(ibatch)
       batch_segment%offset=this%get_batch_offset(ibatch)
       batch_segment%length=this%get_batch_length(ibatch)
       batch_segment%rank=this%batch_ranks(ibatch)
       batch_segment%comm=this%comm

       if(this%batch_ranks(ibatch)==rank) then

          ! A straightforward allocate() call doesn't work on members of
          ! derived types, so instead we use assignment which causes the
          ! array to be allocated implictly. Copying an array of the
          ! correct size triggers the required allocation.

          batch_segment%data=empty_batch(1:batch_segment%length)

       end if
    end do

    ! Create the darray from the array of segments
    batches_darray=new_darray(batch_segments,this%comm)

  end function get_batches_darray

  subroutine load_ensemble_state(this,istep,local_batches)

    !! Get the ensemble state from the model interface, divide into
    !! assimilation batches, and transmit the batch data to each processor

    class(assim_batch_manager)::this
        !! Batch manager

    integer,intent(in)::istep
        !! Assimilation step
    real(kind=8),intent(out),target::local_batches(this%batch_size,this%n_local_batches,this%n_ensemble)
        !! Array of local batch data

    integer::rank ! MPI rank
    integer::ierr ! MPI status code
    integer::imember,ibatch,ibatch_local ! Loop counters

    type(darray)::state_darray   ! Model state darray from the model interface
    type(darray)::batches_darray ! darray with segments aligning to batch ranks

    call mpi_comm_rank(this%comm,rank,ierr)

    call this%model_interface%read_state(istep)

    ! Get the assimilation batches darray
    batches_darray=this%get_batches_darray()

    do imember=1,this%n_ensemble

       ! Get the model state darray from the model interface
       state_darray=this%model_interface%get_state_darray(istep,imember)

       ! Transfer the data to the batches array (this sends all the data to the
       ! correct processor ranks for processing)
       call state_darray%transfer_to_darray(batches_darray)

       ! Reset the local batch counter
       ibatch_local=1

       do ibatch=1,this%n_batches
          if(this%batch_ranks(ibatch)==rank) then

             ! Copy batch data to the local_batches array
             local_batches(:,ibatch_local,imember)=batches_darray%segments(ibatch)%data

             ! Increment the local batch counter
             ibatch_local=ibatch_local+1

          end if
       end do

    end do

  end subroutine load_ensemble_state

  subroutine receive_results(this,istep,local_batches)

    class(assim_batch_manager)::this
    integer,intent(in)::istep
    real(kind=8),intent(in)::local_batches(this%batch_size,this%n_local_batches,this%n_ensemble)
    integer::imember,ibatch, ierr
    integer::rank,ibatch_local,isegment

    type(darray),target::state_darray   ! Model state darray from the model interface
    type(darray)::batches_darray ! darray with segments aligning to batch ranks

    type(darray_segment),pointer::state_segment

    call mpi_comm_rank(this%comm,rank,ierr)

    ! Get the assimilation batches darray
    batches_darray=this%get_batches_darray()

    do imember=1,this%n_ensemble

       ! Reset local batch counter
       ibatch_local=1

       do ibatch=1,this%n_batches
          if(this%batch_ranks(ibatch)==rank) then

             ! Copy batch data to the local_batches array
             batches_darray%segments(ibatch)%data=local_batches(:,ibatch_local,imember)

             ! Increment the local batch counter
             ibatch_local=ibatch_local+1

          end if

       end do

       ! Get the model state darray from the model interface
       state_darray=this%model_interface%get_state_darray(istep,imember)

       ! Transfer the data to state_darray
       call batches_darray%transfer_to_darray(state_darray)

       do isegment=1,size(state_darray%segments)

          state_segment=>state_darray%segments(isegment)

          if(state_segment%rank==rank) then

             ! Store the data in the model interface
             call this%model_interface%set_state_subset( &
                  istep,imember,state_segment%offset, &
                  state_segment%length,state_segment%data)

          end if

       end do

    end do

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
