program advect1d_assimmilate

  use advect1d_assimilate_interfaces
  use assimilate
  use mpi
  use iso_c_binding

  implicit none

  integer::ierr,istep,batch_size,n_ensemble,state_size,comm_size,rank, &
       n_observations

  type(c_ptr)::interface_ptr

  call mpi_init(ierr)

  call mpi_comm_size(mpi_comm_world,comm_size,ierr)

  call mpi_comm_rank(mpi_comm_world,rank,ierr)

  ! Parse command line arguments
  call parse_arguments(istep,n_ensemble,batch_size,state_size)

  ! Initialize i/o interface for accessing data for assimilation
  call init_interface(interface_ptr,n_ensemble,state_size,batch_size,comm_size,rank)

  n_observations=get_batch_observation_count(interface_ptr,istep,1,batch_size,MPI_COMM_WORLD,rank)

  ! Run the assimilation
  call assimilate_parallel(interface_ptr,istep,n_ensemble,batch_size, &
       state_size,n_observations,n_observations,MPI_COMM_WORLD, &
       load_ensemble_state,u_transmit_results=transmit_results, &
       u_store_results=store_results, &
       u_get_batch_observation_count=get_batch_observation_count, &
       u_get_batch_observations=get_batch_observations, &
       u_get_batch_predictions=get_batch_predictions, &
       u_get_batch_innovations=get_batch_innovations, &
       u_add_obs_err=add_obs_err, &
       u_localize=localize)

  ! Cleanup i/o interface
  call cleanup_interface(interface_ptr)

  call mpi_finalize(ierr)

end program advect1d_assimmilate

subroutine parse_arguments(istep,n_ensemble,batch_size,state_size)

  integer::i,istep,batch_size,n_ensemble,state_size
  character(len=32) :: arg

  state_size=0
  n_ensemble=0
  batch_size=1
  istep=0

  i=1

  do while (i <= command_argument_count())

     call get_command_argument(i,arg)

     select case(arg)

     case('-i','--istep')

        i=i+1
        call get_command_argument(i,arg)
        read(arg,*) istep

     case('-n','--n_ensemble')

        i=i+1
        call get_command_argument(i,arg)
        read(arg,*) n_ensemble

     case('-s','--state_size')

        i=i+1
        call get_command_argument(i,arg)
        read(arg,*) state_size

     case('-b','--batch_size')

        i=i+1
        call get_command_argument(i,arg)
        read(arg,*) batch_size

        if(batch_size<1) then
           print *,'Batch size cannot be less than 1'
           stop
        end if

     case('-h','--help')
        call print_help()
        stop

     case default

        print '(a,a,/)', 'Unrecognized command-line option: ', arg
        call print_help()
        stop

     end select

     i=i+1

  end do

end subroutine parse_arguments

subroutine print_help

  print *,'Usage: advect1d_assimilate [OPTIONS]'
  print *,''
  print *,'Options:'
  print *,'  -h, --help        Print usage information and exit'
  print *,'  -i, --istep       Time step or assimilation index'
  print *,'  -b, --batch_size  Size of assimilation batches'
  print *,'  -n, --n_ensemble  Number of ensemble members'
  print *,'  -s, --state_size  Size of model state array'

end subroutine print_help
