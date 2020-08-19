#include "mpi_types.h"

module system_mpi
#ifdef HAVE_MPI_F08_MODULE
  use mpi_f08
#elif HAVE_MPI_F90_MODULE
  use mpi
#else
  include 'mpif.h'
#endif

  implicit none

contains

  subroutine system_mpi_waitall(reqs)
    !! Wrapper for MPI_Waitall. Automatically allocates an arry to hold the
    !! of the appropriate type for the MPI interface in use, then calls
    !! MPI_Waitall on the provided array of requests

    MPI_REQUEST_TYPE,intent(inout)::reqs(:)
        !! Array of MPI requests

#ifdef HAVE_MPI_F08_MODULE
    MPI_STATUS_TYPE,allocatable::statuses(:)
        !! Statuses returned from MPI requests
#else
    MPI_STATUS_TYPE,allocatable::statuses(:,:)
        !! Statuses returned from MPI requests
#endif

    integer::ierr ! MPI error code

    ! Allocate the statuses array
#ifdef HAVE_MPI_F08_MODULE
    allocate(statuses(size(reqs)))
#else
    allocate(statuses(MPI_STATUS_SIZE,size(reqs)))
#endif

    call MPI_Waitall(size(reqs),reqs,statuses,ierr)

  end subroutine system_mpi_waitall

end module system_mpi
