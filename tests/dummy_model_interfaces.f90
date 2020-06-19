module dummy_model_interfaces
  use assimilation_model_interface, ONLY: base_model_interface
  use system_mpi
  use iso_c_binding
  use random, ONLY: random_normal

  implicit none

  type, extends(base_model_interface)::dummy_model_interface
     private
     real(kind=8),allocatable::observations(:),obs_errors(:),predictions(:,:)
     real(kind=8),pointer::local_io_data(:,:)
     real(kind=8)::cutoff,cutoff_u_a
     integer,allocatable::obs_positions(:),io_ranks(:)
     integer::n_observations,state_size,comm,local_io_size
     logical::observations_read,predictions_computed,state_loaded
   contains
     procedure::get_member_state
     procedure::get_subset_io_segment_data
     procedure::get_receive_buffer
     procedure::get_subset_obs_count
     procedure::get_subset_predictions
     procedure::get_subset_observations
     procedure::get_subset_obs_err
     procedure::get_weight_obs_obs
     procedure::get_weight_model_obs
     procedure::before_loading_ensemble_state
     procedure,private::compute_predictions
     procedure,private::load_ensemble_state
     procedure,private::load_observations
     procedure::after_ensemble_results_received
     procedure::get_ensemble_state
  end type dummy_model_interface

contains

  function new_dummy_model(n_ensemble,n_observations,state_size,comm) result(this)

    integer(c_int),intent(in)::n_ensemble,n_observations,state_size,comm
    type(dummy_model_interface)::this
    integer::ierr,rank,comm_size

    this%n_ensemble=n_ensemble
    this%n_observations=n_observations
    this%state_size=state_size
    this%comm=comm

    this%observations_read=.false.
    this%predictions_computed=.false.
    this%state_loaded=.false.

    call mpi_comm_rank(comm,rank,ierr)
    call mpi_comm_size(comm,comm_size,ierr)

    allocate(this%observations(this%n_observations))
    allocate(this%obs_errors(this%n_observations))
    allocate(this%obs_positions(this%n_observations))
    allocate(this%predictions(this%n_observations,this%n_ensemble))
    allocate(this%io_ranks(n_ensemble))

    call get_io_ranks(comm_size,n_ensemble,this%io_ranks)

    ! Get the number of ensemble members read and written locally
    this%local_io_size=get_rank_io_size(n_ensemble,this%io_ranks,rank)

    ! Allocate array for local i/o data
    allocate(this%local_io_data(state_size,this%n_ensemble))

  end function new_dummy_model

  subroutine get_io_ranks(comm_size,n_ensemble,io_ranks)
    integer,intent(in)::comm_size,n_ensemble
    integer,intent(inout)::io_ranks(n_ensemble)
    integer::i,stride
    stride=max(comm_size/n_ensemble,1)

    do i=1,n_ensemble
       io_ranks(i)=mod(i*stride,comm_size)
    end do

  end subroutine get_io_ranks

  function get_rank_io_size(n_ensemble,io_ranks,rank) result(size)
    integer,intent(in)::n_ensemble
    integer,intent(in)::io_ranks(n_ensemble)
    integer,intent(in)::rank
    integer::i,size

    size=0

    do i=1,n_ensemble
       if(io_ranks(i)==rank) size=size+1
    end do

  end function get_rank_io_size

  function get_subset_obs_count(this,istep,subset_offset,subset_size) result(obs_count)

    implicit none

    class(dummy_model_interface)::this
    integer,intent(in)::istep,subset_offset,subset_size
    integer::obs_count

    obs_count=this%n_observations

  end function get_subset_obs_count

  subroutine get_subset_predictions(this,istep,subset_offset,subset_size,predictions)

    use iso_c_binding

    class(dummy_model_interface)::this
    integer,intent(in)::istep,subset_offset,subset_size
    real(kind=8),intent(inout)::predictions(:,:)
    integer::i,imember,ierr

    if(size(predictions,1) /= this%n_observations .or. &
       size(predictions,2) /= this%n_ensemble) then
       print '(A,I0,A,I0,A,I0,A,I0,A)', &
            'Wrong shape passed to predictions argument of get_subset_predictions. Expected (', &
            this%n_observations,',',this%n_ensemble, &
            '), got (',size(predictions,1),',',size(predictions,2),').'
       call mpi_abort(this%comm,1,ierr)
    end if

    predictions=this%predictions

  end subroutine get_subset_predictions

  subroutine get_subset_observations(this,istep,subset_offset,subset_size,observations)

    use iso_c_binding

    implicit none

    class(dummy_model_interface)::this
    integer,intent(in)::istep,subset_offset,subset_size
    real(kind=8),intent(out)::observations(:)
    integer::ierr

    if(size(observations) /= this%n_observations) then
       print '(A,I0,A,I0,A,I0)', &
            'Wrong size array passed to observations argument of get_batch_observations. Expected size=', &
            this%n_observations,', got size=',size(observations)
       call mpi_abort(this%comm,1,ierr)
    end if

    observations=this%observations

  end subroutine get_subset_observations

  SUBROUTINE get_subset_obs_err(this,istep,subset_offset,subset_size,obs_err)
    USE iso_c_binding
    class(dummy_model_interface)::this
    integer,intent(in)::istep,subset_offset,subset_size
    REAL(c_double), INTENT(out) :: obs_err(:)

    obs_err=this%obs_errors

  END SUBROUTINE get_subset_obs_err

  function get_weight_obs_obs(this,istep,subset_offset,subset_size,iobs1,iobs2) result(weight)

    class(dummy_model_interface)::this
    integer,intent(in)::istep,subset_offset,subset_size,iobs1,iobs2
    real(kind=8)::weight
    real(kind=8)::pos1,pos2,delta,distance
    integer::domain_size

    weight=1

  end function get_weight_obs_obs

  function get_weight_model_obs(this,istep,subset_offset,subset_size,imodel,iobs) result(weight)

    class(dummy_model_interface)::this
    integer,intent(in)::istep,subset_offset,subset_size,imodel,iobs
    real(kind=8)::weight
    real(kind=8)::pos_obs,pos_model,delta,distance,cutoff
    integer::domain_size

    weight=1

  end function get_weight_model_obs

  subroutine load_observations(this,istep)
    class(dummy_model_interface)::this
    integer,intent(in)::istep
    integer::ierr,rank,iobs
    real(kind=8)::r

    call mpi_comm_rank(this%comm,rank,ierr)

    if(this%observations_read) return

    call mpi_bcast(this%n_observations,1,MPI_INTEGER,0,this%comm,ierr)

    do iobs=1,this%n_observations
       call random_number(r)
       this%obs_positions(iobs)=1
       this%obs_positions(iobs)=floor(1+r*this%n_observations)
    end do

    call random_number(this%observations)

    call mpi_bcast(this%observations,this%n_observations,MPI_INTEGER,0, &
         this%comm,ierr)
    call mpi_bcast(this%obs_positions,this%n_observations,MPI_INTEGER,0, &
         this%comm,ierr)
    call mpi_bcast(this%obs_errors,this%n_observations,MPI_INTEGER,0, &
         this%comm,ierr)

    this%observations_read=.true.

  end subroutine load_observations

  subroutine load_ensemble_state(this,istep)
    class(dummy_model_interface)::this
    integer,intent(in)::istep
    integer::imember,rank,ierr,i
    real(kind=8)::r

    call mpi_comm_rank(this%comm,rank,ierr)

    ! Populate local i/o data
    do imember=1,this%n_ensemble
       do i=1,this%state_size
          call random_number(r)
          this%local_io_data(i,imember)=i+r/2
       end do
    end do

    ! Broadcast so all processes have the same ensemble state
    call MPI_Bcast(this%local_io_data,this%state_size*this%n_ensemble,MPI_DOUBLE_PRECISION,0,this%comm,ierr)

    call this%load_observations(istep)

    this%state_loaded=.true.

  end subroutine load_ensemble_state

  subroutine compute_predictions(this,istep)
    class(dummy_model_interface)::this
    integer,intent(in)::istep
    integer::imember,local_io_counter,rank,ierr,iobs
    real(kind=8)::member_predictions(this%n_observations)

    call mpi_comm_rank(this%comm,rank,ierr)

    if(.not. this%state_loaded) call this%load_ensemble_state(istep)

    if(.not. this%observations_read) call this%load_observations(istep)

    local_io_counter=1

    do imember=1,this%n_ensemble
       if(this%io_ranks(imember)==rank) then

          ! Compute predictions for this ensemble member
          do iobs=1,this%n_observations
             member_predictions(iobs)=this%local_io_data( &
                  this%obs_positions(iobs)+1,imember)
          end do

       end if

       ! Broadcast to all processors
       call mpi_bcast(member_predictions,this%n_observations, &
            MPI_DOUBLE_PRECISION,this%io_ranks(imember),this%comm,ierr)

       this%predictions(:,imember)=member_predictions

    end do
  end subroutine compute_predictions

  subroutine before_loading_ensemble_state(this,istep)
    class(dummy_model_interface)::this
    integer,intent(in)::istep
    integer::imember,local_io_counter,rank,ierr,iobs

    call mpi_comm_rank(this%comm,rank,ierr)

    call this%load_ensemble_state(istep)

    call this%compute_predictions(istep)

  end subroutine before_loading_ensemble_state

  subroutine get_subset_io_segment_data(this,istep,imember,subset_offset,subset_size,counts,offsets)
    class(dummy_model_interface)::this
    integer,intent(in)::istep,imember,subset_offset,subset_size
    integer,intent(out)::counts(:),offsets(:)
    integer::comm_size,ierr

    call mpi_comm_size(this%comm,comm_size,ierr)

    if(size(counts)/=comm_size) then
       print '(A,I0,A,I0,A)','Wrong size passed to counts argument of get_subset_io_segment_data. Expected ', &
            comm_size,', got ',size(counts),'.'
       error stop
    end if

    if(size(offsets)/=comm_size) then
       print '(A,I0,A,I0,A)','Wrong size passed to counts argument of get_subset_io_segment_data. Expected ', &
            comm_size,', got ',size(offsets),'.'
       error stop
    end if

    ! Fill counts with zeros
    counts=0
    offsets=0

    ! Set counts for the rank assigned to do i/o for requested ensemble member
    counts(this%io_ranks(imember)+1)=subset_size
    offsets(this%io_ranks(imember)+1)=0

  end subroutine get_subset_io_segment_data

  function get_receive_buffer(this,istep,imember,subset_offset,subset_size) result(buffer)

    class(dummy_model_interface)::this
    integer,intent(in)::istep,imember,subset_offset,subset_size
    real(kind=8),pointer::buffer(:)
    real(kind=8),target::empty(0)
    integer::rank,ierr

    call mpi_comm_rank(this%comm,rank,ierr)

    if(this%io_ranks(imember)==rank) then

       ! Point to the appropriate position in local_io_data
       buffer=>this%local_io_data( &
            subset_offset+1:subset_offset+subset_size+1, &
            imember)

    else
       ! Return a null pointer
       buffer=>empty
    end if

  end function get_receive_buffer

  subroutine get_member_state(this,istep,imember,subset_offset,subset_size,subset_state)

    class(dummy_model_interface)::this
    integer,intent(in)::istep,imember,subset_offset,subset_size
    real(kind=8),intent(out)::subset_state(subset_size)
    integer::local_io_index,rank,ierr

    call mpi_comm_rank(this%comm,rank,ierr)

    subset_state=this%local_io_data( &
         subset_offset+1:subset_offset+subset_size+1, &
         imember)

  end subroutine get_member_state

  subroutine after_ensemble_results_received(this,istep)
    class(dummy_model_interface)::this
    integer,intent(in)::istep
    integer::imember,rank,ierr,imember_local

    call MPI_Comm_rank(this%comm,rank,ierr)

    imember_local=1

    do imember=1,this%n_ensemble

       if(this%io_ranks(imember)==rank) then

          imember_local=imember_local+1

       end if

    end do

  end subroutine after_ensemble_results_received

  function get_ensemble_state(this) result(local_io_data)
    class(dummy_model_interface)::this
    real(kind=8)::local_io_data(this%state_size,this%n_ensemble)
    local_io_data=this%local_io_data
  end function get_ensemble_state

end module dummy_model_interfaces
