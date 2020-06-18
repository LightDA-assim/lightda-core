module dummy_assimilator_tests
  use dummy_assimilator, ONLY:dummy_assimilator_assimilate
  use iso_c_binding
contains
  subroutine localize(istep,ibatch,dim_p,dim_obs,HP_p,HPH,mgr)
    class(*),intent(in)::mgr
    INTEGER(c_int32_t), INTENT(in), value :: istep, ibatch, dim_p, dim_obs
    REAL(c_double), INTENT(inout) :: HP_p(dim_obs,dim_p), HPH(dim_obs,dim_obs)

  end subroutine localize

  subroutine add_obs_err(istep,ibatch,dim_obs,HPH,mgr)
    class(*),intent(in)::mgr
    INTEGER(c_int32_t), INTENT(in), value :: istep, ibatch, dim_obs
    REAL(c_double), INTENT(inout) :: HPH(dim_obs,dim_obs)

  end subroutine add_obs_err

  subroutine test_unchanged_data()
    integer,parameter::istep=1,ibatch=1,batch_size=10,n_obs_batch=80,n_ensemble=30
    real(kind=8),parameter::forget=0.6
    real(kind=8)::batch_mean_state(n_obs_batch)
    real(kind=8)::batch_states_orig(batch_size,n_ensemble)
    real(kind=8)::batch_states(batch_size,n_ensemble)
    real(kind=8)::predictions(n_obs_batch,n_ensemble)
    real(kind=8)::innovations(n_obs_batch,n_ensemble)
    integer::ierr,batch_manager=0

    call random_number(predictions)
    call random_number(innovations)
    call random_number(batch_mean_state)
    call random_number(batch_states)

    batch_states_orig=batch_states

    call dummy_assimilator_assimilate(istep,ibatch,batch_size,n_obs_batch, &
         n_obs_batch,n_ensemble,int(0),batch_mean_state,batch_states, &
         predictions,innovations,add_obs_err,localize,forget,ierr,batch_manager)

    if(.not. all(batch_states==batch_states_orig)) then
       print *,'Batch states were changed by dummy assimilator'
       error stop
    end if
  end subroutine test_unchanged_data
end module dummy_assimilator_tests

program test_dummy_assimilator
  use dummy_assimilator_tests

  call test_unchanged_data()

end program test_dummy_assimilator
