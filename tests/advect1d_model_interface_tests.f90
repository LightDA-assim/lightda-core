program dummy_model_interface_tests
  use model_interface_tests, ONLY: run_all
  use advect1d_assimilate_interfaces,ONLY:advect1d_interface,new_advect1d_interface
  use system_mpi
  implicit none
  type(advect1d_interface)::model_interface
  integer,parameter::n_observations=5
  integer,parameter::state_size=100
  integer,parameter::n_ensemble=15
  integer::ierr

  call mpi_init(ierr)

  model_interface=new_advect1d_interface(n_ensemble,n_observations,state_size,mpi_comm_world)

  call run_all(model_interface)
  call mpi_finalize(ierr)
end program dummy_model_interface_tests
