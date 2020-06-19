program dummy_model_interface_tests
  use model_interface_tests, ONLY: run_all
  use dummy_model_interfaces,ONLY:dummy_model_interface,new_dummy_model
  use system_mpi
  implicit none
  type(dummy_model_interface)::model_interface
  integer,parameter::n_observations=5
  integer,parameter::state_size=100
  integer,parameter::n_ensemble=15
  integer::ierr

  call mpi_init(ierr)

  model_interface=new_dummy_model(n_ensemble,n_observations,state_size,mpi_comm_world)

  call run_all(model_interface)
  call mpi_finalize(ierr)
end program dummy_model_interface_tests
