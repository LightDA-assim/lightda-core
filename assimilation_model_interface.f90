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

       ! Arguments
       class(base_model_interface)::this
           !! Model interface
       integer,intent(in)::istep
           !! Iteration number
       integer,intent(in)::subset_offset
           !! Offset of subset from start of state array
       integer,intent(in)::subset_size
           !! Size of subset
       real(kind=8),intent(out)::predictions(:,:)
           !! Predicted values. Will have shape (n_observations,n_ensemble).

     end subroutine I_get_subset_predictions

     function I_get_subset_obs_count(this,istep,subset_offset,subset_size) result(obs_count)
       !! Get the number of observations affecting a given subset of the
       !! model domain. This will be the length of the array returned by I_get_subset_observations.

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
       !! Get the values of observations affecting a given subset of the
       !! model domain.

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
       real(kind=8),intent(out)::observations(:)
           !! Values of observations. Length must equal the value returned
           !! by get_subset_obs_count

     end subroutine I_get_subset_observations

     subroutine I_get_subset_obs_err(this,istep,subset_offset,subset_size,obs_err)
       !! Get the errors (uncertainties) associated with the observations
       !! affecting a given subset of the model domain.

       USE iso_c_binding
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
       REAL(c_double), INTENT(out) :: obs_err(:)
           !! Values of observation errors. Length must equal the value
           !! returned by get_subset_obs_count.

     end subroutine I_get_subset_obs_err

     subroutine I_get_io_ranks(this,istep,imember,ranks,counts,offsets)
       import base_model_interface

       !! Provides the rank assignments for model state i/o, in the form of
       !! an array of processor ranks, and arrays of lengths and offsets
       !! indicating what portion of the state array is to be read/written
       !! by each processor.
       !!
       !! For a parallel model, the model state may be distributed across
       !! multiple processors, with no single processor holding the entire
       !! model state in memory. Similarly, different processors may be
       !! responsible for processing different ensemble members. This subroutine
       !! lets other components know what processor holds each part of model
       !! state.
       !!
       !! The three arrays ranks, counts, and offsets returned by get_io_ranks
       !! are assumed to have the same length. The ranks array is a sequence
       !! of MPI processor ranks, each of which holds a segment of the model
       !! state array. The counts array is a sequence of integers indicating the size of the segment of the model state array held by each processor. The
       !! offsets array indicates the locations of these segments in the model
       !! state array.
       !!
       !! As an example, if a model state of size 10 is divided evenly across
       !! two processors with ranks 0 and 1, the return values of get_io_ranks
       !! would be
       !!
       !!    :::fortran
       !!    ranks   = ( 0, 1 )
       !!    counts  = ( 5, 5 )
       !!    offsets = ( 0, 5 )
       !!
       !! In the case of a serial model, the entire model state will be held by
       !! a single processor, and the return values from get_io_ranks should
       !! indicate this. For example, for a serial model with a state array of
       !! length 10 held on processor 0, the return value from get_io_ranks
       !! would be
       !!
       !!    :::fortran
       !!    ranks   = (/ 0 /)
       !!    counts  = (/ 10 /)
       !!    offsets = (/ 0 /)

       implicit none

       ! Arguments
       class(base_model_interface)::this
           !! Model interface
       integer,intent(in)::istep
           !! Iteration number
       integer, intent(in)::imember
           !! Ensemble member index
       integer,intent(out),allocatable::ranks(:)
           !! Array of processor ranks which hold portions of the model state
           !! for the requested ensemble member
       integer,intent(out),allocatable::counts(:)
           !! Length of the segments of the model state array
       integer,intent(out),allocatable::offsets(:)
           !! Offsets indicating where each segment begins from the start of the
           !! model state array
     end subroutine I_get_io_ranks

     function I_get_state_subset(this,istep,imember,subset_offset,subset_size) result(buffer)

       !! Returns model state values for a given subset of the model state.
       !! Note that availability of data may depend on the rank of the calling
       !! processor, so this subroutine should only be called after a call to
       !! get_io_ranks, and the values of subset_offset and subset_size should
       !! be chosen to fall within a continguous segment of data stored on the
       !! calling processor rank as determined by the output from get_io_ranks.
       !! Failure to satisfy this requirement may result in an error.

       import base_model_interface

       implicit none

       ! Arguments
       class(base_model_interface)::this
           !! Model interface
       integer,intent(in)::istep
           !! Iteration number
       integer, intent(in)::imember
           !! Ensemble member index
       integer,intent(in)::subset_offset
           !! Offset of subset from start of state array
       integer,intent(in)::subset_size
           !! Size of subset
       real(kind=8)::buffer(subset_size)
           !! Values of the model state in the requested subset

     end function I_get_state_subset

     subroutine I_set_state_subset(this,istep,imember,subset_offset,subset_size,subset_state)

       !! Returns model state values for a given subset of the model state
       !! Note that ability to write data may depend on the rank of the calling
       !! processor, so this subroutine should only be called after a call to
       !! get_io_ranks, and the values of subset_offset and subset_size should
       !! be chosen to fall within a continguous segment of data stored on the
       !! calling processor rank as determined by the output from get_io_ranks.
       !! Failure to satisfy this requirement may result in an error.

       import base_model_interface

       implicit none

       ! Arguments
       class(base_model_interface)::this
           !! Model interface
       integer,intent(in)::istep
           !! Iteration number
       integer, intent(in)::imember
           !! Ensemble member index
       integer,intent(in)::subset_offset
           !! Offset of subset from start of state array
       integer,intent(in)::subset_size
           !! Size of subset
       real(kind=8),intent(in)::subset_state(subset_size)
           !! Values of the model state in the requested subset
     end subroutine I_set_state_subset

  end interface

contains
  
  subroutine read_state(this,istep)

    !! Load the model state from disk for a given iteration `istep`. This
    !! subroutine is to be called before any calls to get_state_subset for the
    !! specified iteration `istep`. The default implementation is a no-op,
    !! enabling a model interface to load the model state incrementally as
    !! ensemble state is requested.

    ! Arguments
    class(base_model_interface)::this
        !! Model interface
    integer,intent(in)::istep
        !! Iteration number

  end subroutine read_state

  subroutine write_state(this,istep)

    !! Record the model state after assimilation of `istep`. This
    !! subroutine is to be called after all calls to set_state_subset have
    !! completed for the specified iteration `istep`. The default
    !! implementation is a no-op, enabling a model interface to load the model
    !! state incrementally as ensemble state is received from the assimilation
    !! workers.

    ! Arguments
    class(base_model_interface)::this
        !! Model interface
    integer,intent(in)::istep
        !! Iteration number

  end subroutine write_state

  subroutine get_innovations(this,istep,subset_offset,subset_size,observations,predictions,obs_errors,innovations)

    !! Compute innovations for a given subset of the model domain.

    use random

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
    real(kind=8),intent(in)::observations(:)
        !! Observation values for the subset
    real(kind=8),intent(in)::obs_errors(:)
        !! Observation errors for the subset
    real(kind=8),intent(in)::predictions(:,:)
        !! Predictions for the subset
    real(kind=8),intent(out)::innovations(:,:)
        !! Innovations for the subset
    integer::imember,iobs,obs_count

    obs_count=this%get_subset_obs_count(istep,subset_offset,subset_size)

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

    !! Get localization weight for a given pair of observations. Default
    !! implementation returns 1 for any input (i.e., no localization).

    ! Arguments
    class(base_model_interface)::this
        !! Model interface
    integer,intent(in)::istep
        !! Iteration number
    integer,intent(in)::iobs1
        !! Index of the first observation
    integer,intent(in)::iobs2
        !! Index of the second observation

    ! Returns
    real(kind=8)::weight
        !! Localization weight

    weight=1

  end function get_weight_obs_obs

  function get_weight_model_obs(this,istep,imodel,iobs) result(weight)

    !! Get localization weight for a given observation at a given index in the
    !! model state. Default implementation returns 1 for any input (i.e., no
    !! localization).

    ! Arguments
    class(base_model_interface)::this
        !! Model interface
    integer,intent(in)::istep
        !! Iteration number
    integer,intent(in)::imodel
        !! Index in the model state array
    integer,intent(in)::iobs
        !! Index in the observations array

    ! Returns
    real(kind=8)::weight
        !! Localization weight

    weight=1

  end function get_weight_model_obs

end module assimilation_model_interface
