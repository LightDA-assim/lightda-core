#include "mpi_types.h"

module distributed_array

  use system_mpi
  use exceptions, ONLY: throw, error_status, new_exception

  implicit none

  type :: darray_segment
     integer::rank
     integer::offset
     MPI_COMM_TYPE::comm
     real(kind=8),allocatable::data(:)
   contains
     procedure::write_data=>write_segment_data
     procedure::read_data=>write_segment_data
  end type darray_segment

  type :: darray
     MPI_COMM_TYPE::comm
     type(darray_segment),allocatable::segments
   contains
  end type darray

contains

  subroutine transfer_data(source,dest)
    !! Transfer data between two distributed arrays

    class(darray)::source
        !! Source array
    class(darray)::dest
        !! Destination array

  end subroutine transfer_data

  subroutine write_segment_data(this,offset,data,status)
    !! Write to the segment's data array.
    !! Note that the offset must correspond to a location within the extent
    !! of the segment and the length of data must be small enough to stay
    !! within the segment.
    !! This routine must only be called on the processor whose rank corresponds
    !! to `this%rank`.

    use exceptions, ONLY: throw, new_exception, error_status

    class(darray_segment)::this
        !! Distributed array segment
    integer,intent(in)::offset
        !! Offset relative to the start of the global distributed array
    real(kind=8),intent(in)::data(:)
        !! Data to write
    class(error_status),allocatable,optional::status
        !! Error status

    integer::rank !! MPI rank
    integer::ierr !! Error code returned from MPI

    call mpi_comm_rank(this%comm,rank,ierr)

    if(rank/=this%rank) then
       call throw(status, new_exception( &
            'write_segment_data called from a foreign processor rank', &
            'write_segment_data'))
       return
    end if

    if(offset<this%offset .or. &
         offset+size(data)>this%offset+size(this%data)) then
       call throw(status, &
            new_exception('Attempt to write to data outside of segment', &
            'write_segment_data'))
       return
    end if

  end subroutine write_segment_data

  function read_segment_data(this,offset,size) result(data)
    !! Read data from the segment starting at `offset` and continuing for
    !! `size` array elements.
    !! Note that the `offset` and `size` must correspond to locations within the
    !! extent of the segment.
    !! This routine must only be called on the processor whose rank corresponds
    !! to `this%rank`.

    ! Arguments
    class(darray_segment)::this
        !! Distributed array segment
    integer,intent(in)::offset
        !! Offset relative to the start of the global distributed array
    integer,intent(in)::size
        !! Number of elements to read

    real(kind=8)::data(size)
        !! Data that was read
    

  end function read_segment_data
  
end module distributed_array
