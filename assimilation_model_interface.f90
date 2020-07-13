module assimilation_model_interface
  implicit none

  type,abstract::base_model_interface
     !! Base class for model interfaces
     integer::n_ensemble
   contains
     procedure(I_get_subset_predictions), deferred::get_subset_predictions
     procedure(I_get_subset_observations), deferred::get_subset_observations
     procedure(I_get_subset_obs_count), deferred::get_subset_obs_count
     procedure(I_get_subset_obs_err), deferred::get_subset_obs_err
     procedure::get_innovations
     procedure(I_get_state_size), deferred::get_state_size
     procedure(I_get_subset_io_segment_data), deferred::get_subset_io_segment_data
     procedure(I_get_state_subset_buffer), deferred::get_state_subset_buffer
     procedure::get_weight_obs_obs
     procedure::get_weight_model_obs
     procedure::write_state
  end type base_model_interface

  abstract interface

     function I_get_state_size(this) result(size)
       !! Returns the number elements in the model state array
       import base_model_interface

       implicit none

       ! Arguments
       !! Model interface
       class(base_model_interface)::this
       !! State size
       integer::size

     end function I_get_state_size

     subroutine I_get_subset_predictions(this,istep,subset_offset,subset_size,predictions)
       !! Get prediction values for a subset of the model state
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

     function I_get_state_subset_buffer(this,istep,imember,subset_offset,subset_size) result(buffer)
       import base_model_interface

       implicit none

       class(base_model_interface)::this
       integer,intent(in)::istep,imember,subset_offset,subset_size
       real(kind=8),pointer::buffer(:)

     end function I_get_state_subset_buffer

     subroutine I_get_member_state(this,istep,imember,subset_offset,subset_size,subset_state)
       import base_model_interface

       implicit none

       class(base_model_interface)::this
       integer,intent(in)::istep,imember,subset_offset,subset_size
       real(kind=8),intent(out)::subset_state(subset_size)
     end subroutine I_get_member_state

  end interface

contains
  
  subroutine write_state(this,istep)
    class(base_model_interface)::this
    integer,intent(in)::istep

    ! Subroutine to be called after the complete assimilated ensemble state has been
    ! received from the assimilation workers

    ! Left empty since some models may do disk i/o incrementally as results are received

  end subroutine write_state

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

  function get_weight_obs_obs(this,istep,iobs1,iobs2) result(weight)

    class(base_model_interface)::this
    integer,intent(in)::istep,iobs1,iobs2
    real(kind=8)::weight
    real(kind=8)::pos1,pos2,delta,distance
    integer::domain_size

    weight=1

  end function get_weight_obs_obs

  function get_weight_model_obs(this,istep,imodel,iobs) result(weight)

    class(base_model_interface)::this
    integer,intent(in)::istep,imodel,iobs
    real(kind=8)::weight
    real(kind=8)::pos_obs,pos_model,delta,distance,cutoff
    integer::domain_size

    weight=1

  end function get_weight_model_obs

end module assimilation_model_interface
