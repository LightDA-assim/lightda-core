#include "mpi_types.h"

program test_darray

  use darray_tests, ONLY: test_darray_transfer
  use system_mpi

  integer::ierr ! MPI error code

  call mpi_init(ierr)

  call test_darray_transfer()

  call mpi_finalize(ierr)

end program test_darray
