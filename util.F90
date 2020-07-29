#include "mpi_types.h"

module util
  use system_mpi
  implicit none
  
  interface append_array
     module procedure append_array_real8, append_array_mpi_request
  end interface append_array
contains
  function append_array_real8(a,x) result(b)
    real(kind=8),intent(in)::a(:),x
    real(kind=8),allocatable::b(:)

    allocate(b(size(a)+1))
    b=[a,[x]]
  end function append_array_real8

  function append_array_mpi_request(a,x) result(b)
    MPI_REQUEST_TYPE,intent(in)::a(:),x
    MPI_REQUEST_TYPE,allocatable::b(:)

    allocate(b(size(a)+1))
    b=[a,[x]]
  end function append_array_mpi_request

end module util
