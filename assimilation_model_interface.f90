module assimilation_model_interface
  implicit none

  type,abstract::base_model_interface
     integer::n_ensemble
   contains
     procedure(I_get_subset_predictions), deferred::get_subset_predictions
     procedure(I_get_subset_observations), deferred::get_subset_observations
     procedure(I_get_subset_obs_count), deferred::get_subset_obs_count
     procedure(I_get_subset_obs_err), deferred::get_subset_obs_err
     procedure::get_innovations=>get_innovations
     procedure::before_loading_ensemble_state
     procedure::after_ensemble_state_loaded
     procedure(I_get_subset_io_segment_data), deferred::get_subset_io_segment_data
     procedure(I_get_member_state), deferred::get_member_state
     procedure(I_get_receive_buffer), deferred::get_receive_buffer
     procedure(I_get_weight_obs_obs), deferred::get_weight_obs_obs
     procedure(I_get_weight_model_obs), deferred::get_weight_model_obs
     procedure::after_member_state_received
     procedure::after_ensemble_results_received
  end type base_model_interface

  abstract interface

     subroutine I_get_subset_predictions(this,istep,subset_offset,subset_size,predictions)
       import base_model_interface

       implicit none

       class(base_model_interface)::this
       integer,intent(in)::istep,subset_offset,subset_size
       real(kind=8),intent(inout)::predictions(:,:)
     end subroutine I_get_subset_predictions

     function I_get_subset_obs_count(this,istep,subset_offset,subset_size) result(obs_count)
       import base_model_interface

       implicit none

       class(base_model_interface)::this
       integer,intent(in)::istep,subset_offset,subset_size
       integer::obs_count
     end function I_get_subset_obs_count

     subroutine I_get_subset_observations(this,istep,subset_offset,subset_size,observations)
       import base_model_interface

       implicit none

       class(base_model_interface)::this
       integer,intent(in)::istep,subset_offset,subset_size
       real(kind=8),intent(out)::observations(:)
     end subroutine I_get_subset_observations

     subroutine I_get_subset_obs_err(this,istep,subset_offset,subset_size,obs_err)
       USE iso_c_binding
       import base_model_interface

       implicit none

       class(base_model_interface)::this
       integer,intent(in)::istep,subset_offset,subset_size
       REAL(c_double), INTENT(out) :: obs_err(:)
     end subroutine I_get_subset_obs_err

     subroutine I_get_subset_io_segment_data(this,istep,imember,subset_offset,subset_size,counts,offsets)
       import base_model_interface

       implicit none

       class(base_model_interface)::this
       integer,intent(in)::istep,imember,subset_offset,subset_size
       integer,intent(out)::counts(:),offsets(:)
     end subroutine I_get_subset_io_segment_data

     function I_get_receive_buffer(this,istep,imember,subset_offset,subset_size) result(buffer)
       import base_model_interface

       implicit none

       class(base_model_interface)::this
       integer,intent(in)::istep,imember,subset_offset,subset_size
       real(kind=8),pointer::buffer(:)

     end function I_get_receive_buffer

     subroutine I_get_member_state(this,istep,imember,subset_offset,subset_size,subset_state)
       import base_model_interface

       implicit none

       class(base_model_interface)::this
       integer,intent(in)::istep,imember,subset_offset,subset_size
       real(kind=8),intent(out)::subset_state(subset_size)
     end subroutine I_get_member_state

     function I_get_weight_obs_obs(this,istep,subset_offset,subset_size,iobs1,iobs2) result(weight)
       import base_model_interface

       implicit none

       class(base_model_interface)::this
       integer,intent(in)::istep,subset_offset,subset_size,iobs1,iobs2
       real(kind=8)::weight
     end function I_get_weight_obs_obs

     function I_get_weight_model_obs(this,istep,subset_offset,subset_size,imodel,iobs) result(weight)
       import base_model_interface

       implicit none

       class(base_model_interface)::this
       integer,intent(in)::istep,subset_offset,subset_size,imodel,iobs
       real(kind=8)::weight
     end function I_get_weight_model_obs

  end interface

contains
  
  subroutine before_loading_ensemble_state(this,istep)
    class(base_model_interface)::this
    integer,intent(in)::istep

    ! Empty procedure which can be overridden by child classes that want to
    ! run certain tasks before ensemble state batches are requested

    ! Will be called before the first call to get_member_state
  end subroutine before_loading_ensemble_state

  subroutine after_ensemble_state_loaded(this,istep)
    class(base_model_interface)::this
    integer,intent(in)::istep

    ! Empty procedure which can be overridden by child classes that want to
    ! run certain tasks (e.g. loading observations or computing forward
    ! operators) after ensemble state is loaded

    ! Will be called after the last call to get_member_state and
    ! before assimilation begins
  end subroutine after_ensemble_state_loaded

  subroutine after_member_state_received(this,istep,imember,subset_offset,subset_size)

    implicit none

    class(base_model_interface)::this
    integer,intent(in)::istep,imember,subset_offset,subset_size

    ! Empty procedure which can be overridden by child classes that want to
    ! run certain tasks (e.g. loading observations or computing forward
    ! operators) after a member state is received

    ! Will be called after the buffer returned by the corresponding call to
    ! get_receive_buffer is populated with new data
  end subroutine after_member_state_received

  subroutine after_ensemble_results_received(this,istep)
    class(base_model_interface)::this
    integer,intent(in)::istep

    ! Empty procedure which can be overridden by child classes that want to
    ! run certain tasks before ensemble state batches are requested

    ! Will be called after the complete assimilated ensemble state has been
    ! received from the workers
  end subroutine after_ensemble_results_received

  subroutine get_innovations(this,istep,batch_offset,batch_length,observations,predictions,obs_errors,innovations)

    use random

    implicit none

    class(base_model_interface)::this
    integer,intent(in)::istep,batch_offset,batch_length
    real(kind=8),intent(in)::observations(:), &
         obs_errors(:),predictions(:,:)
    real(kind=8),intent(out)::innovations(:,:)
    integer::imember,iobs,obs_count

    obs_count=this%get_subset_obs_count(istep,batch_offset,batch_length)

    if(size(observations)/=obs_count) then
       print '(A,I0,A,I0)','Observations array has wrong length. Expected ',obs_count,', got ',size(observations)
       stop
    end if

    if(size(predictions,1)/=obs_count .or. &
         size(predictions,2)/=this%n_ensemble) then
       print '(A,I0,A,I0,A,I0,A,I0,A)', &
            'Predictions array has wrong shape. Expected (' &
            ,obs_count,',',this%n_ensemble,'), got (', &
            size(predictions,1),',',size(predictions,2),')'
       stop
    end if

    if(size(innovations,1)/=obs_count .or. &
         size(innovations,2)/=this%n_ensemble) then
       print '(A,I0,A,I0,A,I0,A,I0,A)', &
            'Innovations array has wrong shape. Expected (', &
            obs_count,',',this%n_ensemble,'), got (' &
            ,size(innovations,1),',',size(innovations,2),')'
       stop
    end if

    do imember=1,this%n_ensemble
       do iobs=1,obs_count
          innovations(iobs,imember)=observations(iobs) - &
               predictions(iobs,imember) + &
               random_normal()*obs_errors(iobs)
       end do
    end do

  end subroutine get_innovations

end module assimilation_model_interface
