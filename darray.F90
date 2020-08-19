#include "mpi_types.h"

module distributed_array

  use system_mpi
  use exceptions, ONLY: throw, error_status, new_exception

  implicit none

  type :: darray_segment
     !! Process-local contiguous segment of a distributed array
     integer::rank
         !! Processor rank that holds the segment data
     integer::offset
         !! Offset of the segment from the start of the global array
     integer::length
         !! Length of data
     MPI_COMM_TYPE::comm
         !! MPI communicator
     real(kind=8),allocatable::data(:)
         !! Segment data
   contains
     procedure::write_data=>write_segment_data
     procedure::read_data=>write_segment_data
  end type darray_segment

  type :: darray
     !! Distributed array type

     MPI_COMM_TYPE::comm
         !! MPI communicator
     type(darray_segment),allocatable::segments(:)
         !! Array of segments, each stored on one processor
   contains
  end type darray

contains

  subroutine transfer_data(source,dest,status)
    !! Transfer data between two distributed arrays

    ! Arguments
    class(darray),intent(in)::source
        !! Source array
    class(darray),intent(inout)::dest
        !! Destination array
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    integer::rank !! MPI rank
    integer::ierr !! MPI error code

    if(source%comm/=dest%comm) then
       call throw(status,new_exception( &
            'source and destination arrays must use the same MPI communicator', &
            'transfer_data'))
       return
    end if

    call mpi_comm_rank(source%comm,rank,ierr)

  end subroutine transfer_data

  function new_darray(segments,comm,status)

    !! Create a new darray

    ! Arguments
    type(darray_segment),allocatable::segments(:)
        !! Array of darray segments
    MPI_COMM_TYPE,intent(in)::comm
        !! MPI communicator
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    type(darray)::new_darray
        !! New darray to create

    integer::i ! Loop counter
    character(:),allocatable::errstr ! Error string
    integer::expected_offset ! Expected offset for segment

    ! Loop over segments and do sanity checks
    do i=1,size(segments)

       if(segments(i)%comm/=comm) then

          write (errstr,'(A,I0,A)') &
               'Segment ',i,' uses a different communicator than the one specified for the darray'

          call throw(status,new_exception(errstr, 'new_darray'))

          return

       end if

       if(i==1) then

          ! Check that first segment is at start of array
          if(segments(i)%offset/=0) then
             call throw(status,new_exception( &
                  'First segment offset must be 0', &
                  'new_darray'))
             return
          end if

       else

          ! Check that subsequent segments cover the array without gaps

          expected_offset=segments(i-1)%offset+segments(i-1)%length

          if(expected_offset/=segments(i)%offset) then
             write(errstr,'(A,I0,A,I0,A,I0,A,I0,A)') &
                  'Expected segment(',i,')%offset=', &
                  expected_offset,', got segment(',i,')%offset=', &
                  segments(i)%offset,'.'
             call throw(status,new_exception(errstr,'new_darray'))
             return
          end if

       end if

       if(segments(i)%length>0 .and. &
            segments(i)%length/=size(segments(i)%data)) then

          write (errstr,'(A,I0,A,I0,A,I0,A,I0,A)') &
               'segment(',i,')%length is ',segments(i)%length, &
               ' but segment(',i,')%data has a size of ', &
               size(segments(i)%data),'.'

          call throw(status,new_exception(errstr,'new_darray'))

          return

       end if

    end do

    new_darray%segments=segments
    new_darray%comm=comm

  end function new_darray

  subroutine write_segment_data(this,offset,data,status)
    !! Write to the segment's data array.
    !! Note that the offset must correspond to a location within the extent
    !! of the segment and the length of data must be small enough to stay
    !! within the segment.
    !! This routine must only be called on the processor whose rank corresponds
    !! to `this%rank`.

    class(darray_segment),intent(inout)::this
        !! Distributed array segment
    integer,intent(in)::offset
        !! Offset relative to the start of the global distributed array
    real(kind=8),intent(in)::data(:)
        !! Data to write
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    integer::rank !! MPI rank
    integer::ierr !! Error code returned from MPI
    integer::local_start, local_end !! Range on the local array

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

    ! Determine segment-local array indices
    local_start=offset-this%offset+1
    local_end=offset-this%offset+size(data)

    ! Copy data to output array
    this%data(local_start:local_end)=data

  end subroutine write_segment_data

  function read_segment_data(this,offset,length,status) result(data)
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
    integer,intent(in)::length
        !! Number of elements to read
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    real(kind=8)::data(length)
        !! Data that was read

    integer::rank !! MPI rank
    integer::ierr !! Error code returned from MPI
    integer::local_start, local_end !! Range on the local array

    character(*),parameter::nameproc='read_segment_data'

    call mpi_comm_rank(this%comm,rank,ierr)

    if(rank/=this%rank) then
       call throw(status, new_exception( &
            'write_segment_data called from a foreign processor rank', &
            'write_segment_data'))
       return
    end if

    if( offset<this%offset .or. &
         offset+length > this%offset+size(this%data)) then
       call throw(status, &
            new_exception('Attempt to write to data outside of segment', &
            'write_segment_data'))
       return
    end if

    ! Determine segment-local array indices
    local_start=offset-this%offset+1
    local_end=offset-this%offset+length

    ! Copy data to output array
    data=this%data(local_start:local_end)

  end function read_segment_data
  
end module distributed_array
