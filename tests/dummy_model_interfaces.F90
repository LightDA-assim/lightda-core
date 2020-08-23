#include "mpi_types.h"

module dummy_model_interfaces
  use assimilation_model_interface, ONLY: base_model_interface
  use system_mpi
  use iso_c_binding
  use random, ONLY: random_normal
  use exceptions, ONLY: throw, new_exception, error_status
  use distributed_array, ONLY: darray, darray_segment, new_darray

  implicit none

  type, extends(base_model_interface)::dummy_model_interface
     private
     real(kind=8),allocatable::observations(:),obs_errors(:),predictions(:,:)
     real(kind=8),pointer::local_io_data(:,:)
     real(kind=8)::cutoff,cutoff_u_a
     integer,allocatable::obs_positions(:),io_ranks(:)
     integer::n_observations,state_size,local_io_size
     MPI_COMM_TYPE::comm
     logical::observations_read,predictions_computed,state_loaded
   contains
     procedure::get_state_size
     procedure::get_io_ranks
     procedure::get_state_subset
     procedure::set_state_subset
     procedure::get_subset_obs_count
     procedure::get_subset_predictions
     procedure::get_subset_observations
     procedure::get_subset_obs_err
     procedure,private::compute_predictions
     procedure::read_state
     procedure,private::load_observations
     procedure::get_ensemble_state
     procedure::get_state_darray
  end type dummy_model_interface

contains

  function new_dummy_model(n_ensemble,n_observations,state_size,comm) result(this)

    integer(c_int),intent(in)::n_ensemble,n_observations,state_size
    MPI_COMM_TYPE::comm
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

    this%io_ranks=get_member_ranks(comm_size,n_ensemble)

    ! Get the number of ensemble members read and written locally
    this%local_io_size=get_rank_io_size(n_ensemble,this%io_ranks,rank)

    ! Allocate array for local i/o data
    allocate(this%local_io_data(state_size,this%n_ensemble))

  end function new_dummy_model

  function get_state_size(this,istep,status) result(size)
    class(dummy_model_interface)::this
    integer,intent(in)::istep
    integer::size
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    size=this%state_size

  end function get_state_size

  function get_member_ranks(comm_size,n_ensemble) result(io_ranks)
    integer,intent(in)::comm_size,n_ensemble
    integer::io_ranks(n_ensemble)
    integer::i,stride
    stride=max(comm_size/n_ensemble,1)

    do i=1,n_ensemble
       io_ranks(i)=mod(i*stride,comm_size)
    end do

  end function get_member_ranks

  subroutine get_io_ranks(this,istep,imember,ranks,counts,offsets,status)
    class(dummy_model_interface)::this
    integer,intent(in)::istep,imember
    integer,intent(out),allocatable::ranks(:),counts(:),offsets(:)
    integer::comm_size,ierr
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    call mpi_comm_size(this%comm,comm_size,ierr)

    allocate(ranks(1),counts(1),offsets(1))

    ranks(1)=this%io_ranks(imember)
    counts(1)=this%state_size
    offsets(1)=0

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

  function get_subset_obs_count(this,istep,subset_offset,subset_size,status) result(obs_count)

    implicit none

    class(dummy_model_interface)::this
    integer,intent(in)::istep,subset_offset,subset_size
    integer::obs_count
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    obs_count=this%n_observations

  end function get_subset_obs_count

  subroutine get_subset_predictions(this,istep,subset_offset,subset_size,predictions,status)

    use iso_c_binding

    class(dummy_model_interface)::this
    integer,intent(in)::istep,subset_offset,subset_size
    real(kind=8),intent(out)::predictions(:,:)
    integer::i,imember,ierr
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    character(:),allocatable::errstr

    if(size(predictions,1) /= this%n_observations .or. &
       size(predictions,2) /= this%n_ensemble) then
       write(errstr,'(A,I0,A,I0,A,I0,A,I0,A)') &
            'Wrong shape passed to predictions argument of get_subset_predictions. Expected (', &
            this%n_observations,',',this%n_ensemble, &
            '), got (',size(predictions,1),',',size(predictions,2),').'
       call throw(status,new_exception(errstr,'get_subset_predictions'))
       return
    end if

    predictions=this%predictions

  end subroutine get_subset_predictions

  subroutine get_subset_observations(this,istep,subset_offset,subset_size,observations,status)

    use iso_c_binding

    implicit none

    class(dummy_model_interface)::this
    integer,intent(in)::istep,subset_offset,subset_size
    real(kind=8),intent(out)::observations(:)
    integer::ierr
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    character(:),allocatable::errstr

    if(size(observations) /= this%n_observations) then
       write(errstr,'(A,I0,A,I0,A,I0)') &
            'Wrong size array passed to observations argument of get_batch_observations. Expected size=', &
            this%n_observations,', got size=',size(observations)
       call throw(status,new_exception(errstr,'get_subset_observations'))
       return
    end if

    observations=this%observations

  end subroutine get_subset_observations

  SUBROUTINE get_subset_obs_err(this,istep,subset_offset,subset_size,obs_err,status)
    USE iso_c_binding
    class(dummy_model_interface)::this
    integer,intent(in)::istep,subset_offset,subset_size
    REAL(c_double), INTENT(out) :: obs_err(:)
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    obs_err=this%obs_errors

  END SUBROUTINE get_subset_obs_err

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

  subroutine read_state(this,istep,status)
    class(dummy_model_interface)::this
    integer,intent(in)::istep
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    integer::imember,rank,ierr,i
    real(kind=8)::r

    call mpi_comm_rank(this%comm,rank,ierr)

    ! Populate local i/o data
    do imember=1,this%n_ensemble
       do i=1,this%state_size
          this%local_io_data(i,imember)=imember*this%state_size+i
       end do
    end do

    ! Broadcast so all processes have the same ensemble state
    call MPI_Bcast(this%local_io_data,this%state_size*this%n_ensemble,MPI_DOUBLE_PRECISION,0,this%comm,ierr)

    call this%load_observations(istep)

    this%state_loaded=.true.

  end subroutine read_state

  subroutine compute_predictions(this,istep)
    class(dummy_model_interface)::this
    integer,intent(in)::istep
    integer::imember,local_io_counter,rank,ierr,iobs
    real(kind=8)::member_predictions(this%n_observations)

    call mpi_comm_rank(this%comm,rank,ierr)

    if(.not. this%state_loaded) call this%read_state(istep)

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

  function get_state_darray(this,istep,imember) result(state_darray)

    !! Get the requested ensemble member state as a darray

    ! Arguments
    class(dummy_model_interface)::this
        !! Model interface
    integer,intent(in)::istep
        !! Iteration number
    integer,intent(in)::imember
        !! Ensemble member index

    type(darray)::state_darray
        !! State array represented as a darray object

    type(darray_segment)::segments(1)

    integer::rank !! MPI rank
    integer::ierr !! MPI status code

    call mpi_comm_rank(this%comm,rank,ierr)

    ! Populate the segment indicating that it covers the entire model state
    ! and is stored on the processor rank found in this%io_ranks(imember)
    segments(1)%rank=this%io_ranks(imember)
    segments(1)%comm=this%comm
    segments(1)%offset=0
    segments(1)%length=this%state_size

    if(this%io_ranks(imember)==rank) then
       ! Copy the member state data to the darray segment
       segments(1)%data=this%local_io_data(1:this%state_size,imember)
    end if

    state_darray=new_darray(segments,this%comm)

  end function get_state_darray

  function get_state_subset(this,istep,imember,subset_offset,subset_size,status) result(buffer)

    class(dummy_model_interface)::this
    integer,intent(in)::istep,imember,subset_offset,subset_size
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    real(kind=8)::buffer(subset_size)

    integer::rank,ierr
    character(:),allocatable::errstr

    call mpi_comm_rank(this%comm,rank,ierr)

    if(this%io_ranks(imember)==rank) then

       buffer=this%local_io_data( &
            subset_offset+1:subset_offset+subset_size, &
            imember)

    else
       write(errstr,*) 'Indices for non-local ensemble state passed to get_state_subset. Data for member', &
            imember,'requested on rank',rank
       call throw(status,new_exception(errstr,'get_state_subset'))
       return
    end if

  end function get_state_subset

  subroutine set_state_subset(this,istep,imember,subset_offset,subset_size,subset_state,status)

    class(dummy_model_interface)::this
    integer,intent(in)::istep,imember,subset_offset,subset_size
    real(kind=8),intent(in)::subset_state(subset_size)
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    integer::rank,ierr

    call mpi_comm_rank(this%comm,rank,ierr)

    if(this%io_ranks(imember)==rank) then

       this%local_io_data( &
            subset_offset+1:subset_offset+subset_size, &
            imember)=subset_state

    else
       call throw(status,new_exception('Indices for non-local ensemble state passed to set_state_subset'))
       return
    end if

  end subroutine set_state_subset

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
    integer::imember,ierr

    do imember=1,this%n_ensemble

       call MPI_Bcast(this%local_io_data(:,imember),this%state_size, &
            MPI_DOUBLE_PRECISION,this%io_ranks(imember),this%comm,ierr)

    end do

    local_io_data=this%local_io_data

  end function get_ensemble_state

end module dummy_model_interfaces
