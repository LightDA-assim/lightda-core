module mod_dummy_model_forward_operator

  use forward_operator, ONLY: base_forward_operator
  use observations, ONLY: observation_set
  use random_observations, ONLY: random_observation_set
  use dummy_model_interfaces, ONLY: dummy_model_interface
  use exceptions, ONLY: error_status, throw, new_exception
  use system_mpi

  implicit none

  type, extends(base_forward_operator) :: dummy_model_forward_operator

    class(dummy_model_interface), pointer::model_interface => null()

  contains
    procedure::get_predictions_mask
    procedure::get_predictions
    procedure, private::get_random_observation_predictions

  end type dummy_model_forward_operator

contains

  function new_dummy_model_forward_operator(model_interface)

    ! Arguments
    class(dummy_model_interface), target :: model_interface

    ! Result
    type(dummy_model_forward_operator)::new_dummy_model_forward_operator

    new_dummy_model_forward_operator%model_interface => model_interface

  end function new_dummy_model_forward_operator

  function get_predictions_mask(this, istep, obs_set, status) result(mask)

    !! Returns a mask array indicating which observations in `obs_set`
    !! can be predicted by the model

    ! Arguments
    class(dummy_model_forward_operator)::this
        !! Forward operator
    integer, intent(in)::istep
        !! Assimilation step
    class(observation_set)::obs_set
        !! Observation set
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Result
    logical, allocatable::mask(:)
        !! Mask array

    integer::iobs

    allocate (mask(obs_set%get_size()))

    if (.not. associated(this%model_interface)) then
      call throw(status, new_exception('Model interface is undefined', &
                                       'get_predictions_mask'))
      return
    end if

    select type (obs_set)
    class is (random_observation_set)
      mask = .true.
    class default

      ! Unknown observation type, set mask to false
      mask = .false.

    end select

  end function get_predictions_mask

  function get_predictions(this, istep, obs_set, status) result(predictions)

    !! Get the predictions for the observations in `obs_set`

    ! Arguments
    class(dummy_model_forward_operator)::this
        !! Forward operator
    integer, intent(in)::istep
        !! Assimilation step
    class(observation_set)::obs_set
        !! Observation set
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    ! Result
    real(kind=8), allocatable::predictions(:, :)
        !! Mask array

    integer::iobs

    if (.not. associated(this%model_interface)) then
      call throw(status, new_exception('Model interface is undefined', &
                                       'get_predictions_mask'))
      return
    end if

    allocate (predictions(obs_set%get_size(), this%model_interface%n_ensemble))

    select type (obs_set)
    class is (random_observation_set)

      predictions = this%get_random_observation_predictions( &
                    istep, obs_set, status)

    class default

      ! Unknown observation type, set all values to 0
      predictions = 0

    end select

  end function get_predictions

  function get_random_observation_predictions(this, istep, obs_set, status) &
    result(predictions)

    use distributed_array, ONLY: darray

    !! Compute predictions for all ensemble members and broadcast to all
    !! processors

    ! Arguments
    class(dummy_model_forward_operator)::this
        !! Model interface
    integer, intent(in)::istep
        !! Iteration number
    type(random_observation_set), intent(in) :: obs_set
    class(error_status), intent(out), allocatable, optional::status
        !! Error status

    integer::imember, rank, ierr, iobs
    real(kind=8), allocatable::member_predictions(:)
    real(kind=8), allocatable::predictions(:, :)

    real(kind=8), pointer::member_state(:)

    class(darray), allocatable, target::state

    integer::state_size

    call mpi_comm_rank(this%model_interface%comm, rank, ierr)

    if (.not. associated(this%model_interface)) then
      call throw(status, new_exception('Model interface is undefined', &
                                       'get_predictions_mask'))
      return
    end if

    allocate (member_predictions(obs_set%get_size()))
    allocate (predictions(obs_set%get_size(), this%model_interface%n_ensemble))

    state_size = this%model_interface%get_state_size(istep)

    do imember = 1, this%model_interface%n_ensemble
      state = this%model_interface%get_state_darray(istep, imember)

      if (state%segments(1)%rank == rank) then

        member_state => state%segments(1)%data

        ! Compute predictions for this ensemble member
        do iobs = 1, obs_set%get_size()
          member_predictions(iobs) = &
            member_state( &
            max(min( &
                int(obs_set%get_position(iobs)*state_size + 1), &
                state_size), &
                1))
        end do

      end if

      ! Broadcast to all processors
      call mpi_bcast(member_predictions, obs_set%get_size(), &
                     MPI_DOUBLE_PRECISION, state%segments(1)%rank, &
                     state%comm, ierr)

      predictions(:, imember) = member_predictions

    end do

  end function get_random_observation_predictions

end module mod_dummy_model_forward_operator
