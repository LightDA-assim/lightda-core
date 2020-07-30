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
     procedure(I_get_io_ranks), deferred::get_io_ranks
     procedure(I_get_state_subset), deferred::get_state_subset
     procedure(I_set_state_subset), deferred::set_state_subset
     procedure::get_weight_obs_obs
     procedure::get_weight_model_obs
     procedure::read_state
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
       real(kind=8),intent(out)::predictions(:,:)
     end subroutine I_get_subset_predictions

     function I_get_subset_obs_count(this,istep,subset_offset,subset_size) result(obs_count)
       !! Get the number of observations impacting a given subset
       import base_model_interface

       implicit none

       ! Arguments
       class(base_model_interface)::this
           !! Model interface
       integer,intent(in)::istep
           !! Iteration number
       integer,intent(in)::subset_offset
           !! Offset of subset from start of state array
       integer,intent(in)::subset_size
           !! Size of subset

       integer::obs_count
           !! Number of observations

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

     subroutine I_get_io_ranks(this,istep,imember,ranks,counts,offsets)
       import base_model_interface

       implicit none

       class(base_model_interface)::this
       integer,intent(in)::istep,imember
       integer,intent(out),allocatable::ranks(:),counts(:),offsets(:)
     end subroutine I_get_io_ranks

     function I_get_state_subset(this,istep,imember,subset_offset,subset_size) result(buffer)
       import base_model_interface

       implicit none

       class(base_model_interface)::this
       integer,intent(in)::istep,imember,subset_offset,subset_size
       real(kind=8)::buffer(subset_size)

     end function I_get_state_subset

     subroutine I_set_state_subset(this,istep,imember,subset_offset,subset_size,subset_state)
       import base_model_interface

       implicit none

       class(base_model_interface)::this
       integer,intent(in)::istep,imember,subset_offset,subset_size
       real(kind=8),intent(in)::subset_state(subset_size)
     end subroutine I_set_state_subset

  end interface

contains
  
  subroutine read_state(this,istep)
    class(base_model_interface)::this
    integer,intent(in)::istep

    ! Subroutine to be called before requesting any model state data for the given iteration istep
    ! received from the assimilation workers

    ! Left empty since some models may choose to do disk i/o incrementally as ensemble state is requested

  end subroutine read_state

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
