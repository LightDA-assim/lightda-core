module assimilate

  use mpi
  use iso_c_binding

  implicit none

contains

  subroutine shuffle(a)
    integer, intent(inout) :: a(:)
    integer :: i, randpos, temp
    real :: r
 
    do i = size(a), 2, -1
       call random_number(r)
       randpos = int(r * i) + 1
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

    ! 
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

  function get_batch_offset(batch_size,ibatch) result(offset)
    integer,intent(in)::batch_size,ibatch
    integer::offset
    offset=batch_size*(ibatch-1)
  end function get_batch_offset

  function get_batch_length(batch_size,ibatch,state_size) result(length)
    integer,intent(in)::batch_size,ibatch,state_size
    integer::offset,length
    offset=get_batch_offset(batch_size,ibatch)
    length=min(state_size,offset+batch_size)-offset
  end function get_batch_length

  function get_rank_batch_count(n_batches,batch_ranks,rank) result(count)
    integer,intent(in)::n_batches,batch_ranks(n_batches),rank
    integer::ibatch,count
    count=0

    do ibatch=1,n_batches
       if(batch_ranks(ibatch)==rank) count=count+1
    end do

  end function get_rank_batch_count

  subroutine assimilate_parallel(interface_info,istep,n_ensemble,batch_size, &
       state_size,n_observations,n_obs_batch_max,comm,U_load_ensemble_state,U_transmit_results, &
       U_store_results,U_get_batch_observation_count,U_get_batch_observations, &
       U_get_batch_predictions,U_get_batch_innovations,U_add_obs_err,U_localize)

    use lenkf_rsm, ONLY: lenkf_analysis_rsm

    abstract interface

       subroutine load_ensemble_state(interface_info,istep,rank,comm, &
            state_size,batch_ranks,local_batches,n_batches,n_local_batches, &
            batch_size,n_ensemble)

         use iso_c_binding
         implicit none

         type(c_ptr),intent(inout)::interface_info
         integer,intent(in)::istep,rank,comm,n_local_batches,state_size,batch_size,n_batches,n_ensemble
         integer,intent(in)::batch_ranks(n_batches)
         real(kind=8),intent(out)::local_batches(n_local_batches,state_size,n_ensemble)
         
       end subroutine load_ensemble_state

       subroutine transmit_results(info_ptr,istep,ibatch,batch_state,batch_size,n_ensemble,batch_ranks,n_batches,comm,rank)
         use iso_c_binding
         implicit none
         type(c_ptr),intent(inout)::info_ptr
         real(kind=8),intent(in)::batch_state(batch_size,n_ensemble)
         integer,intent(in)::istep,ibatch,batch_size,n_ensemble,comm,rank,n_batches
         integer,intent(in)::batch_ranks(n_batches)
         
       end subroutine transmit_results

       subroutine store_results(interface_info,istep,batch_size,batch_ranks,n_batches,rank,comm,state_size)

         use iso_c_binding

         implicit none

         type(c_ptr),intent(inout)::interface_info
         integer,intent(in)::istep,rank,comm,state_size,batch_size,n_batches
         integer,intent(in)::batch_ranks(n_batches)

       end subroutine store_results

       function get_batch_observation_count(interface_info,istep,ibatch, &
            batch_size,comm,rank) result(n_obs_batch)

         use iso_c_binding

         implicit none

         type(c_ptr),intent(inout)::interface_info
         integer,intent(in)::istep,ibatch,batch_size,comm,rank
         integer::n_obs_batch

       end function get_batch_observation_count

       subroutine get_batch_predictions(interface_info,istep,ibatch,batch_size,n_ensemble,n_obs_batch,rank,comm,predictions)

         use iso_c_binding

         implicit none

         type(c_ptr),intent(inout)::interface_info
         integer,intent(in)::istep,ibatch,rank,comm,batch_size,n_ensemble, &
              n_obs_batch
         real(kind=8),intent(inout)::predictions(n_obs_batch,n_ensemble)

       end subroutine get_batch_predictions

       subroutine get_batch_innovations(interface_info,istep,ibatch,batch_size,n_ensemble,n_obs_batch,rank,comm,innovations)

         use iso_c_binding

         implicit none

         type(c_ptr),intent(inout)::interface_info
         integer,intent(in)::istep,ibatch,rank,comm,batch_size,n_ensemble, &
              n_obs_batch
         real(kind=8),intent(inout)::innovations(n_obs_batch,n_ensemble)

       end subroutine get_batch_innovations

       subroutine get_batch_observations(interface_info,istep,ibatch,batch_size,n_obs_batch,rank,comm,observations)

         use iso_c_binding

         implicit none

         type(c_ptr),intent(inout)::interface_info
         integer,intent(in)::istep,ibatch,rank,comm,batch_size,n_obs_batch
         real(kind=8),intent(inout)::observations(n_obs_batch)

       end subroutine get_batch_observations

       SUBROUTINE add_obs_err(step,ind_p,dim_obs,HPH,info_ptr)
         ! Add observation error covariance matrix
         USE iso_c_binding
         type(c_ptr),intent(inout)::info_ptr
         INTEGER(c_int32_t), INTENT(in), value :: step, ind_p, dim_obs
         REAL(c_double), INTENT(inout) :: HPH(dim_obs,dim_obs)
       END SUBROUTINE add_obs_err

       SUBROUTINE localize(step,ind_p,dim_p,dim_obs,HP_p,HPH,info_ptr)
         ! Apply localization to HP and HPH^T
         USE iso_c_binding
         INTEGER(c_int32_t), INTENT(in), value :: step, ind_p, dim_p, dim_obs
         REAL(c_double), INTENT(inout) :: HP_p(dim_obs,dim_p), HPH(dim_obs,dim_obs)
         type(c_ptr),intent(inout)::info_ptr
       END SUBROUTINE localize

  end interface

    integer,intent(in) :: istep,n_ensemble,batch_size,state_size,comm,n_observations,n_obs_batch_max
    type(c_ptr)::interface_info
    integer,dimension(:),allocatable :: io_ranks, batch_ranks
    real(kind=8),allocatable::local_batches(:,:,:)
    real(kind=8),allocatable::innovations(:,:),predictions(:,:),observations(:)
    real(kind=8),allocatable::batch_mean_state(:),batch_states(:,:)
    real(kind=8)::forget

    integer::rank,ierr,comm_size,n_batches,n_local_batches,ibatch, &
         n_obs_batch,ibatch_local, batch_length

    procedure(load_ensemble_state) :: U_load_ensemble_state
    procedure(transmit_results) :: U_transmit_results
    procedure(store_results) :: U_store_results
    procedure(get_batch_observation_count) :: U_get_batch_observation_count
    procedure(get_batch_predictions) :: U_get_batch_predictions
    procedure(get_batch_innovations) :: U_get_batch_innovations
    procedure(get_batch_observations) :: U_get_batch_observations
    procedure(add_obs_err) :: U_add_obs_err
    procedure(localize) :: U_localize

    call mpi_comm_rank(comm, rank, ierr)
    call mpi_comm_size(comm, comm_size, ierr)

    forget=1

    ! Get the number of batches and allocate batch arrays
    n_batches=get_batch_count(state_size,batch_size)

    allocate(batch_ranks(n_batches))
    allocate(innovations(n_obs_batch_max,n_ensemble))
    allocate(predictions(n_obs_batch_max,n_ensemble))
    allocate(observations(n_obs_batch_max))
    allocate(batch_mean_state(batch_size))
    allocate(batch_states(batch_size,n_ensemble))

    ! Assign batches to process ranks
    call get_batch_ranks(comm_size,batch_ranks)

    ! Get number of local batches
    n_local_batches=get_rank_batch_count(n_batches,batch_ranks,rank)

    ! Allocate array to hold local ensemble state
    allocate(local_batches(n_local_batches,state_size,n_ensemble))

    ! Load the ensemble state
    call U_load_ensemble_state(interface_info,istep,rank,comm,state_size, &
         batch_ranks,local_batches,n_batches,n_local_batches,batch_size, &
         n_ensemble)

    ibatch_local=1

    ! Assimilate local batches
    do ibatch=1,n_batches

       if(batch_ranks(ibatch)/=rank) cycle

       ! Get number of predictions for this batch
       n_obs_batch=U_get_batch_observation_count(interface_info,istep,ibatch,batch_size,comm,rank)

       batch_length=get_batch_length(batch_size,ibatch,state_size)

       if(size(innovations,1)<=n_obs_batch) then
          deallocate(observations)
          allocate(observations(n_obs_batch))
          deallocate(innovations)
          allocate(innovations(n_obs_batch,n_ensemble))
          deallocate(predictions)
          allocate(predictions(n_obs_batch,n_ensemble))
       end if

       call U_get_batch_predictions(interface_info,istep,ibatch,batch_size,n_ensemble,n_obs_batch,rank,comm,predictions)
       call U_get_batch_innovations(interface_info,istep,ibatch,batch_size,n_ensemble,n_obs_batch,rank,comm,innovations)
       call U_get_batch_observations(interface_info,istep,ibatch,batch_size,n_obs_batch,rank,comm,observations)

       batch_states=local_batches(ibatch_local,:,:)

       call lenkf_analysis_rsm(istep,ibatch,batch_size,n_obs_batch, &
            n_obs_batch,n_ensemble,int(0),batch_mean_state,batch_states, &
            predictions,innovations,U_add_obs_err,U_localize,forget,ierr,interface_info)

       call U_transmit_results(interface_info,istep,ibatch, &
            batch_states,batch_size,n_ensemble,batch_ranks, &
            n_batches,comm,rank)

       ibatch_local=ibatch_local+1

    end do

    ! Write the ensemble state
    call U_store_results(interface_info,istep,batch_size,batch_ranks,n_batches,rank,comm,state_size)

    deallocate(local_batches)
    deallocate(batch_ranks)
    deallocate(observations)
    deallocate(innovations)
    deallocate(predictions)
    deallocate(batch_states)

  end subroutine assimilate_parallel

end module assimilate
