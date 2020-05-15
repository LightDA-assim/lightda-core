module system_mpi
#if HAVE_MPIF90_MODULE
  use mpi
#else
  include 'mpif.h'
#endif
end module system_mpi
