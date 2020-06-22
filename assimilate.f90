module assimilate

  use system_mpi
  use assimilation_batch_manager,ONLY:assim_batch_manager
  use assimilation_model_interface
  use iso_c_binding

  implicit none

contains

  subroutine localize(istep,ibatch,dim_p,dim_obs,HP_p,HPH,mgr)
    class(*),intent(in)::mgr
    INTEGER(c_int32_t), INTENT(in), value :: istep, ibatch, dim_p, dim_obs
    REAL(c_double), INTENT(inout) :: HP_p(dim_obs,dim_p), HPH(dim_obs,dim_obs)

    select type(mgr)
    class is (assim_batch_manager)
       call mgr%localize(istep,ibatch,dim_p,dim_obs,HP_p,HPH)
    class default
       print *,'Could not determine argument type. Should be class assim_batch_manager'
    end select
  end subroutine localize

  subroutine add_obs_err(istep,ibatch,dim_obs,HPH,mgr)
    class(*),intent(in)::mgr
    INTEGER(c_int32_t), INTENT(in), value :: istep, ibatch, dim_obs
    REAL(c_double), INTENT(inout) :: HPH(dim_obs,dim_obs)

    select type(mgr)
    class is (assim_batch_manager)
       call mgr%add_obs_err(istep,ibatch,dim_obs,HPH)
    class default
       print *,'Could not determine argument type. Should be class assim_batch_manager'
    end select
  end subroutine add_obs_err

  subroutine assimilate_parallel(batch_manager,assimilator,istep,n_observations,n_obs_batch_max)

    integer,intent(in) :: istep,n_observations,n_obs_batch_max
    real(kind=8),allocatable::local_batches(:,:,:)
    integer,allocatable::local_batch_inds(:)
    real(kind=8),allocatable::innovations(:,:),predictions(:,:),observations(:),obs_errors(:)
    real(kind=8),allocatable::batch_mean_state(:),batch_states(:,:)
    real(kind=8)::forget
    class(assim_batch_manager)::batch_manager
    class(base_model_interface),pointer::model_interface
    integer::comm,rank,ierr,comm_size,n_batches,n_local_batches,ibatch, &
         n_obs_batch, ibatch_local, batch_offset, batch_length, batch_size, &
         state_size,n_ensemble

    abstract interface
       subroutine U_assimilator(step,ind_p,dim_p, dim_obs_p, dim_obs, dim_ens, rank_ana, &
       state_p, ens_p, predictions, innovations, U_add_obs_err, U_localize, forget, flag,info)
         USE iso_c_binding

         IMPLICIT NONE

         abstract INTERFACE
            SUBROUTINE add_obs_err(step,ind_p,dim_obs,HPH,info)
              ! Add observation error covariance matrix
              USE iso_c_binding
              INTEGER(c_int32_t), INTENT(in), value :: step, ind_p, dim_obs
              REAL(c_double), INTENT(inout) :: HPH(dim_obs,dim_obs)
              class(*),intent(in)::info
            END SUBROUTINE add_obs_err
            SUBROUTINE localize(step,ind_p,dim_p,dim_obs,HP_p,HPH,info)
              ! Apply localization to HP and HPH^T
              USE iso_c_binding
              INTEGER(c_int32_t), INTENT(in), value :: step, ind_p, dim_p, dim_obs
              REAL(c_double), INTENT(inout) :: HP_p(dim_obs,dim_p), HPH(dim_obs,dim_obs)
              class(*),intent(in)::info
            END SUBROUTINE localize

         END INTERFACE

         ! !ARGUMENTS:
         INTEGER(c_int32_t), INTENT(in), value :: step      ! Iteration number
         INTEGER(c_int32_t), INTENT(in), value :: ind_p      ! Integer index of PE
         INTEGER(c_int32_t), INTENT(in), value  :: dim_p     ! PE-local dimension of model state
         INTEGER(c_int32_t), INTENT(in), value :: dim_obs_p  ! PE-local dimension of observation vector
         INTEGER(c_int32_t), INTENT(in), value :: dim_obs    ! Global dimension of observation vector
         INTEGER(c_int32_t), INTENT(in), value :: dim_ens   ! Size of state ensemble
         INTEGER(c_int32_t), INTENT(in), value :: rank_ana  ! Rank to be considered for inversion of HPH
         REAL(c_double), INTENT(inout)  :: state_p(dim_p)        ! PE-local ensemble mean state
         REAL(c_double), INTENT(inout)  :: ens_p(dim_p, dim_ens) ! PE-local state ensemble
         REAL(c_double), INTENT(inout)  :: predictions(dim_obs, dim_ens) ! PE-local state ensemble
         REAL(c_double), INTENT(in)     :: innovations(dim_obs, dim_ens) ! Global array of innovations
         REAL(c_double), INTENT(in), value     :: forget    ! Forgetting factor
         INTEGER(c_int32_t), INTENT(out) :: flag    ! Status flag
         class(*),intent(inout)::info

         procedure(add_obs_err) :: U_add_obs_err
         procedure(localize) :: U_localize

       end subroutine U_assimilator
    end interface

    procedure(U_assimilator)::assimilator

    model_interface=>batch_manager%model_interface
    comm=batch_manager%get_comm()
    batch_size=batch_manager%get_batch_size()
    state_size=batch_manager%get_state_size()
    n_ensemble=batch_manager%get_n_ensemble()

    call mpi_comm_rank(comm, rank, ierr)
    call mpi_comm_size(comm, comm_size, ierr)

    forget=0.6

    ! Get the number of batches and allocate batch arrays
    n_batches=batch_manager%get_n_batches()

    allocate(innovations(n_obs_batch_max,n_ensemble))
    allocate(predictions(n_obs_batch_max,n_ensemble))
    allocate(observations(n_obs_batch_max))
    allocate(obs_errors(n_obs_batch_max))
    allocate(batch_mean_state(batch_size))
    allocate(batch_states(batch_size,n_ensemble))

    ! Get number of local batches
    n_local_batches=batch_manager%get_rank_batch_count(rank)

    ! Allocate array to hold local ensemble state
    allocate(local_batches(batch_size,n_local_batches,n_ensemble))

    ! Get local batch indices
    allocate(local_batch_inds(n_local_batches))
    call batch_manager%get_rank_batches(rank,local_batch_inds)

    ! Load the ensemble state
    call batch_manager%load_ensemble_state(istep,local_batches)

    ! Assimilate local batches
    do ibatch_local=1,n_local_batches
       ibatch=local_batch_inds(ibatch_local)

       batch_offset=batch_manager%get_batch_offset(ibatch)
       batch_length=batch_manager%get_batch_length(ibatch)

       ! Get number of predictions for this batch
       n_obs_batch=model_interface%get_subset_obs_count( &
            istep,batch_offset,batch_length)

       if(size(innovations,1)<=n_obs_batch) then
          deallocate(observations)
          allocate(observations(n_obs_batch))
          deallocate(innovations)
          allocate(innovations(n_obs_batch,n_ensemble))
          deallocate(predictions)
          allocate(predictions(n_obs_batch,n_ensemble))
       end if

       call model_interface%get_subset_predictions(istep,batch_offset,batch_size,predictions)
       call model_interface%get_subset_observations(istep,batch_offset,batch_size,observations)
       call model_interface%get_subset_obs_err(istep,batch_offset,batch_size,obs_errors)
       call model_interface%get_innovations(istep,batch_offset,batch_length,observations,predictions,obs_errors,innovations)

       batch_states=local_batches(:,ibatch_local,:)

       call assimilator(istep,ibatch,batch_size,n_obs_batch, &
            n_obs_batch,n_ensemble,int(0),batch_mean_state,batch_states, &
            predictions,innovations,add_obs_err,localize,forget,ierr,batch_manager)

       local_batches(:,ibatch_local,:)=batch_states

    end do

    ! Write the ensemble state
    print *,'Rank',rank,'completed all batches, storing results'
    call batch_manager%store_results(istep,local_batches)

    deallocate(local_batches)
    deallocate(observations)
    deallocate(innovations)
    deallocate(predictions)
    deallocate(batch_states)

  end subroutine assimilate_parallel

end module assimilate
