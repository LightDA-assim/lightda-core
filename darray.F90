#include "mpi_types.h"

module distributed_array

  use system_mpi
  use exceptions, ONLY: throw, error_status, new_exception
  use util, ONLY: str

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
     procedure::get_segments_for_range
     procedure::get_segment_index_for_offset
     procedure::transfer_to_segment
     procedure::transfer_to_darray=>transfer_data
  end type darray

contains

  subroutine transfer_data(source,dest,status)
    !! Transfer data between two distributed arrays

    ! Arguments
    class(darray),intent(in)::source
        !! Source array
    class(darray),target,intent(inout)::dest
        !! Destination array
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    integer::rank !! MPI rank
    integer::ierr !! MPI error code

    class(darray_segment),pointer::dest_segment
        ! Pointer to destination segment
    class(darray_segment),pointer::overlapping_source_segments(:)
        ! Pointer to an array of source segments that overlap with a given
        ! destination segment
    integer::i ! Loop counter

    MPI_REQUEST_TYPE,allocatable::segment_reqs(:) ! Segment MPI requests
    MPI_REQUEST_TYPE,allocatable::reqs(:)         ! Combined MPI requests

    ! Initialize reqs
    allocate(reqs(0))

    if(source%comm/=dest%comm) then
       call throw(status,new_exception( &
            'source and destination arrays must use the same MPI communicator', &
            'transfer_data'))
       return
    end if

    call mpi_comm_rank(source%comm,rank,ierr)

    do i=1,size(dest%segments)

       dest_segment=>dest%segments(i)

       call source%transfer_to_segment(dest_segment,reqs_out=segment_reqs, &
            status=status)

       reqs=[reqs,segment_reqs]

    end do

    call system_mpi_waitall(reqs)

  end subroutine transfer_data

  subroutine transfer_to_segment(this,dest_segment,reqs_out,status)

    !! Transfer data to darray segment `dest_segment`

    ! Arguments
    class(darray),target::this
        !! Distributed array segment
    class(error_status),intent(out),allocatable,optional::status
        !! Error status
    class(darray_segment),intent(inout),pointer::dest_segment
        !! Pointer to destination segment
    MPI_REQUEST_TYPE,intent(out),allocatable,target,optional::reqs_out(:)
        !! MPI requests to be returned to sender

    MPI_REQUEST_TYPE,pointer::reqs(:)
        !! MPI requests

    class(darray_segment),pointer::source_segments(:)
        !! List of segments

    real(kind=8),pointer::source_data(:)
    real(kind=8),pointer::dest_data(:)
    integer::rank !! MPI rank
    integer::ierr !! MPI error code

    integer::i ! Loop counter

    integer::overlap_start,overlap_end
        ! Start/end offsets of an overlapping array section

    call mpi_comm_rank(this%comm,rank,ierr)

    ! Get the segments in this array that overlap with the destination segment
    source_segments=>this%get_segments_for_range( &
         dest_segment%offset,dest_segment%length,status)

    ! Allocate requests array
    if(present(reqs_out)) then
       allocate(reqs_out(size(source_segments)*2))
       reqs=>reqs_out
    else
       allocate(reqs(size(source_segments)*2))
    end if

    reqs=MPI_REQUEST_NULL

    do i=1,size(source_segments)

       ! Locate the overlapping region between source-segments(i) and
       ! dest_segment
       overlap_start=max(source_segments(i)%offset,dest_segment%offset)
       overlap_end=min( &
            source_segments(i)%offset+source_segments(i)%length, &
            dest_segment%offset+dest_segment%length)

       if(dest_segment%rank==rank) then

          ! Configure destination buffer
          dest_data=>dest_segment%data( &
               overlap_start-dest_segment%offset+1: &
               overlap_end-dest_segment%offset)
       end if

       if(source_segments(i)%rank==rank) then

          ! Configure source data buffer
          source_data=>source_segments(i)%data( &
               overlap_start-source_segments(i)%offset+1: &
               overlap_end-source_segments(i)%offset)

          if(dest_segment%rank==rank) then
             ! Copy from the source buffer to the destination buffer
             dest_data=source_data
          else
             ! Send to destination segment rank
             call MPI_Isend(source_data,size(source_data), &
                  MPI_DOUBLE_PRECISION,dest_segment%rank,1,this%comm, &
                  reqs(i),ierr)
          end if
       else


          if(dest_segment%rank==rank) then
             ! Receive from source segment rank
             call MPI_Irecv(dest_data,size(dest_data),MPI_DOUBLE_PRECISION, &
                  source_segments(i)%rank,1,this%comm, &
                  reqs(i+size(source_segments)),ierr)
          end if
       end if

    end do

    if(.not.present(reqs_out)) call system_mpi_waitall(reqs)

  end subroutine transfer_to_segment

  function new_darray(segments,comm,status)

    !! Create a new darray

    ! Arguments
    class(darray_segment),intent(in)::segments(:)
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

    integer::rank !! MPI rank
    integer::ierr !! MPI error code

    call mpi_comm_rank(comm,rank,ierr)

    ! Loop over segments and do sanity checks
    do i=1,size(segments)

       if(segments(i)%comm/=comm) then

          errstr='Segment '//str(i,'(I0)')//' uses a different communicator than the one specified for the darray'

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
             errstr='Expected segment('//str(i)//')%offset='// &
                  str(expected_offset)// &
                  ', got segment('//str(i)//')%offset='// &
                  str(segments(i)%offset)//'.'
             call throw(status,new_exception(errstr,'new_darray'))
             return
          end if

       end if

       if(segments(i)%rank==rank .and. &
            segments(i)%length/=size(segments(i)%data)) then

          errstr='segments('//str(i)//')%length is '//str(segments(i)%length)// &
               ' but segments('//str(i)//')%data has a size of '// &
               str(size(segments(i)%data))//'.'

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

  function get_segments_for_range(this,offset,length,status) result(segments)
    !! Get array segments that cover the range of array elements
    !! `offset`+1:`offset`+`length`

    ! Arguments
    class(darray),target::this
        !! Distributed array segment
    integer,intent(in)::offset
        !! Offset relative to the start of the global distributed array
    integer,intent(in)::length
        !! Length of range
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    class(darray_segment),pointer::segments(:)
        !! List of segments

    integer::istart,iend ! Start and end indices in this%segments

    istart=this%get_segment_index_for_offset(offset,status)
    iend=this%get_segment_index_for_offset(offset+length,status)

    segments=>this%segments(istart:iend)

  end function get_segments_for_range

  function get_segment_index_for_offset(this,offset,status) result(iseg)
    !! Get index of in this%segments for the segment containing the element with offset `offset`

    ! Arguments
    class(darray)::this
        !! Distributed array segment
    integer,intent(in)::offset
        !! Offset relative to the start of the global distributed array
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    integer::iseg !! Segment index

    integer::ilower, iupper ! Bounds of search bracket for matching segment
    integer::imid           ! Midpoint of bracket
    integer::offset_min, offset_max ! Range of valid offsets for the array

    character(:),allocatable::errstr ! Error string

    ilower=1
    iupper=size(this%segments)

    offset_min=this%segments(1)%offset
    offset_max=this%segments(size(this%segments))%offset+ &
         this%segments(size(this%segments))%length

    if(offset<offset_min .or. offset>offset_max) then

       errstr='Request for offset '//str(offset)// &
            ' which is out of the range of the array ('//str(offset_min)// &
            '-'//str(offset_max)//')'
       call throw(status,new_exception(errstr,'get_segment_index_for_offset'))
       return
    end if

    ! Bisection search for the matching segment
    do while(this%segments(ilower)%offset>offset .or. &
         this%segments(ilower)%offset+this%segments(ilower)%length<offset)

       ! Midpoint of bracket
       imid=(ilower+iupper)/2

       ! Reduce bracket size by comparing offset to this%segments(imid)%offset
       if(this%segments(imid)%offset>=offset) iupper=imid
       if(this%segments(imid)%offset<=offset) ilower=imid
       if(this%segments(iupper)%offset<=offset) ilower=iupper

    end do

    iseg=ilower

  end function get_segment_index_for_offset

end module distributed_array
