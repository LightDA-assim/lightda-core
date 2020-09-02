program advect1d_assimmilate

  use advect1d_assimilate_interfaces, ONLY: advect1d_interface, &
                                            new_advect1d_interface
  use assimilation_batch_manager, ONLY: assim_batch_manager, new_batch_manager
  use system_mpi
  use iso_c_binding
  use lenkf_rsm, ONLY: lenkf_analysis_rsm
  use mod_assimilation_manager, ONLY: assimilation_manager, &
                                      new_assimilation_manager
  use mod_lenkf_rsm_filter, ONLY: lenkf_rsm_filter
  use mod_advect1d_forward_operator, ONLY: advect1d_forward_operator
  use advect1d_observations, ONLY: advected_quantity_observation_set
  use mod_advect1d_localization, ONLY: advect1d_localizer

  implicit none

  integer::ierr, istep, batch_size, n_ensemble, state_size, comm_size, rank, &
            n_observations

  type(advect1d_interface), target::model_interface
  type(assim_batch_manager)::batch_manager
  type(assimilation_manager)::assim_mgr
  type(lenkf_rsm_filter)::filter
  type(advect1d_forward_operator)::forward_operator
  type(advected_quantity_observation_set)::observation_sets(1)
  type(advect1d_localizer)::localizer

  call mpi_init(ierr)

  call mpi_comm_size(mpi_comm_world, comm_size, ierr)

  call mpi_comm_rank(mpi_comm_world, rank, ierr)

  ! Parse command line arguments
  call parse_arguments( &
    istep, n_ensemble, n_observations, batch_size, state_size)

  ! Initialize i/o interface for accessing data for assimilation
  model_interface = new_advect1d_interface( &
                    n_ensemble, n_observations, state_size, mpi_comm_world)

  assim_mgr = new_assimilation_manager( &
              model_interface, istep, n_ensemble, forward_operator, &
              observation_sets, batch_size, localizer, filter, mpi_comm_world)

  ! Run the assimilation
  call assim_mgr%assimilate(istep)

  call mpi_finalize(ierr)

end program advect1d_assimmilate

subroutine parse_arguments( &
  istep, n_ensemble, n_observations, batch_size, state_size)

  integer::i, istep, batch_size, n_ensemble, state_size, n_observations
  character(len=32) :: arg

  state_size = 0
  n_ensemble = 0
  n_observations = 0
  batch_size = 1
  istep = 0

  i = 1

  do while (i <= command_argument_count())

    call get_command_argument(i, arg)

    select case (arg)

    case ('-i', '--istep')

      i = i + 1
      call get_command_argument(i, arg)
      read (arg, *) istep

    case ('-n', '--n_ensemble')

      i = i + 1
      call get_command_argument(i, arg)
      read (arg, *) n_ensemble

    case ('-o', '--n_observations')

      i = i + 1
      call get_command_argument(i, arg)
      read (arg, *) n_observations

    case ('-s', '--state_size')

      i = i + 1
      call get_command_argument(i, arg)
      read (arg, *) state_size

    case ('-b', '--batch_size')

      i = i + 1
      call get_command_argument(i, arg)
      read (arg, *) batch_size

      if (batch_size < 1) then
        print *, 'Batch size cannot be less than 1'
        stop
      end if

    case ('-h', '--help')
      call print_help()
      stop

    case default

      print '(a,a,/)', 'Unrecognized command-line option: ', arg
      call print_help()
      stop

    end select

    i = i + 1

  end do

end subroutine parse_arguments

subroutine print_help

  print *, 'Usage: advect1d_assimilate [OPTIONS]'
  print *, ''
  print *, 'Options:'
  print *, '  -h, --help        Print usage information and exit'
  print *, '  -i, --istep       Time step or assimilation index'
  print *, '  -b, --batch_size  Size of assimilation batches'
  print *, '  -n, --n_ensemble  Number of ensemble members'
  print *, '  -s, --state_size  Size of model state array'

end subroutine print_help
