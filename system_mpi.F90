module system_mpi
#ifdef HAVE_MPI_F08_MODULE
  use mpi_f08
#elif HAVE_MPI_F90_MODULE
  use mpi
#else
  include 'mpif.h'
#endif
end module system_mpi
