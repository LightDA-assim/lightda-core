#include "mpi_types.h"

module distributed_array

  use system_mpi
  use exceptions, ONLY: throw, error_container, new_exception, exception
  use util, ONLY: str
  use iso_fortran_env, ONLY: real64

  implicit none

  type, extends(exception)::darray_transfer_error
    !! Error during darray transfer

    integer::gfortran_workaround_1 = 0
      !! Workaround for gfortran bug #110987.
      !! (prevents a segmentation fault triggered when finalizing
      !! temporary objects that inherit from a parent type but don't
      !! have any data members aside from those in the parent)
  end type darray_transfer_error

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
    real(kind=8), allocatable::data(:)
         !! Segment data
  contains
    procedure::write_data => write_segment_data
    procedure::read_data => write_segment_data
  end type darray_segment

  type :: darray_segment_set
     !! Set of darray segments, which may be noncontiguous or overlapping.

    MPI_COMM_TYPE::comm
         !! MPI communicator
    type(darray_segment), allocatable::segments(:)
         !! Array of segments, each stored on one processor
    real(real64), allocatable, private::buffer(:)
        !! Send/receive buffer
    class(darray_segment_set), pointer::foreign => null()
    logical::transfer_pending = .false.
    logical::is_transfer_dest = .false.
  contains
    procedure::get_local_data
  end type darray_segment_set

  type, extends(darray_segment_set) :: darray
     !! Distributed array comprised of contiguous and nonoverlapping darray
     !! segments

    MPI_REQUEST_TYPE::req
        !! Array of MPI requests
    integer, allocatable, private::sendcounts(:), recvcounts(:)
        !! Send/receive counts for MPI_Ialltoallv
    integer, allocatable, private::sdspls(:), rdspls(:)
        !! Displacements for MPI_Ialltoallv

  contains
    procedure::get_segments_for_range
    procedure::get_segment_index_for_offset
    procedure::transfer_to_darray => transfer_data
    procedure::get_length
    procedure::finish_transfer
  end type darray

contains

  function new_darray_transfer_error(message, procedure) result(exc)

    !! Create a new darray transfer error

    ! Arguments
    character(*), intent(in)::message
        !! Error message
    character(*), intent(in), optional::procedure
        !! Procedure where error occured

    type(darray_transfer_error)::exc
        !! New exception object

    exc%message = message

    if (present(procedure)) exc%procedure = procedure

  end function new_darray_transfer_error

  function get_length(this) result(length)

    ! Arguments
    class(darray) :: this
        !! Darray

    ! Result
    integer::length

    integer::i

    length = 0
    do i = 1, size(this%segments)
      length = length + this%segments(i)%length
    end do

  end function get_length

  function get_local_data(this) result(local_data)

    !! Get the portion of the segment data that is stored on the local MPI
    !! process as a flattened array

    ! Arguments
    class(darray_segment_set) :: this
        !! Darray segment set

    ! Result
    real(kind=8), allocatable :: local_data(:)

    integer::rank !! MPI rank
    integer::ierr !! MPI error code

    integer::i ! Loop counter

    integer::local_data_size ! Size of local data
    integer::local_offset ! Offset in the local data array

    call mpi_comm_rank(this%comm, rank, ierr)

    ! Determine size of local data
    local_data_size = 0
    do i = 1, size(this%segments)
      if (rank == this%segments(i)%rank) then
        local_data_size = local_data_size + this%segments(i)%length
      end if
    end do

    ! Allocate the result array
    allocate (local_data(local_data_size))

    ! Copy the data into the array
    local_offset = 0
    do i = 1, size(this%segments)
      if (rank == this%segments(i)%rank) then

        ! Copy the data
        local_data(local_offset + 1:local_offset + this%segments(i)%length) &
          = this%segments(i)%data

        ! Move the local offset to the end of this segment
        local_offset = local_offset + this%segments(i)%length

      end if
    end do

  end function get_local_data

  subroutine get_segment_overlap(source_segment, dest_segment, &
                                 overlap_start, overlap_end)

    !! Locate the overlapping region of two darray segments overlap within the
    !! global array
    class(darray_segment), intent(in)::source_segment
        !! Source segment
    class(darray_segment), intent(in)::dest_segment
        !! Destination segment
    integer, intent(out)::overlap_start
        !! Start offset of overlap range
    integer, intent(out)::overlap_end
        !! End offset of overlap range

    overlap_start = max(source_segment%offset, dest_segment%offset)
    overlap_end = min( &
                  source_segment%offset + source_segment%length, &
                  dest_segment%offset + dest_segment%length)

  end subroutine get_segment_overlap

  subroutine transfer_data(source, dest, finish_immediately, status)
    !! Transfer data between two distributed arrays

    ! Arguments
    class(darray), target, intent(inout)::source
        !! Source array
    class(darray_segment_set), target, intent(inout)::dest
        !! Destination segments
    logical, optional::finish_immediately
        !! Whether transfer should be finished immediately
    type(error_container), intent(out), optional::status
        !! Error status

    integer::rank !! MPI rank
    integer::ierr !! MPI error code

    class(darray_segment), pointer::dest_segment
    ! Pointer to destination segment
    class(darray_segment), pointer::overlapping_source_segments(:)
    ! Pointer to an array of source segments that overlap with a given
    ! destination segment
    integer::i, irank, idest, isource ! Loop counters
    integer, pointer::crsr ! Cursor into the send buffer
    integer, allocatable, target::crsrs(:)
        !! Array of cursors into all the send buffers
    class(darray_segment), pointer::source_segments(:)
        !! List of segments
    real(kind=8), pointer::source_data(:)
    real(kind=8), pointer::dest_data(:)
    integer::dest_rank, source_rank
        !! Ranks of source and destination processes for a message

    integer::comm_size
    integer::overlap_start, overlap_end
    ! Start/end offsets of an overlapping array section
    logical::call_finish
        !! Whether the finish_transfer procedure
        !! should be called

    if (source%transfer_pending .or. dest%transfer_pending) then
      call throw(status, new_darray_transfer_error( &
                 'transfer_data called on an array or segment set ' &
                 //'that already has a transfer pending.'))
      return
    end if

    if (source%comm /= dest%comm) then
      call throw(status, new_exception( &
           'source and destination arrays &
           &must use the same MPI communicator', &
           'transfer_data'))
      return
    end if

    if (present(finish_immediately)) then
      call_finish = finish_immediately
    else
      call_finish = .true.
    end if

    call mpi_comm_rank(source%comm, rank, ierr)

    call mpi_comm_size(source%comm, comm_size, ierr)

    allocate (source%sendcounts(comm_size), source%recvcounts(comm_size), &
              source%sdspls(comm_size), source%rdspls(comm_size), &
              crsrs(comm_size))

    source%sendcounts = 0
    source%recvcounts = 0

    ! Loop over the segments to determine the required sizes for the send and receive buffers
    do idest = 1, size(dest%segments)

      dest_segment => dest%segments(idest)

      ! Get the segments in this array that overlap with the destination segment
      source_segments => source%get_segments_for_range( &
                         dest_segment%offset, dest_segment%length, status)

      dest_rank = dest_segment%rank

      do isource = 1, size(source_segments)

        source_rank = source_segments(isource)%rank

        if (source_rank == rank .or. dest_rank == rank) then

          ! Locate the overlapping region between source-segments(i) and
          ! dest_segment

          call get_segment_overlap(source_segments(isource), dest_segment, &
                                   overlap_start, overlap_end)

        end if

        if (source_rank == rank .and. dest_rank /= rank) then

          ! Add the overlap size to the send count
          source%sendcounts(dest_rank + 1) = &
            source%sendcounts(dest_rank + 1) &
            + overlap_end - overlap_start

        else if (dest_rank == rank .and. source_rank /= rank) then

          ! Add the overlap size to the receive count
          source%recvcounts(source_rank + 1) = &
            source%recvcounts(source_rank + 1) &
            + overlap_end - overlap_start

        end if

      end do

    end do

    ! Allocate send/receive buffers
    allocate (source%buffer(sum(source%sendcounts)))
    allocate (dest%buffer(sum(source%recvcounts)))

    ! Set send buffer displacements
    source%sdspls(1) = 0
    do dest_rank = 2, comm_size
      source%sdspls(dest_rank) = &
        source%sdspls(dest_rank - 1) &
        + source%sendcounts(dest_rank - 1)
    end do

    ! Set receive buffer displacements
    source%rdspls(1) = 0
    do source_rank = 2, comm_size
      source%rdspls(source_rank) = &
        source%rdspls(source_rank - 1) &
        + source%recvcounts(source_rank - 1)
    end do

    ! Initialize cursors
    crsrs = 0

    ! Loop over the segments and populate the send buffer
    do irank = 1, comm_size
      do idest = 1, size(dest%segments)

        dest_segment => dest%segments(idest)
        dest_rank = dest_segment%rank

        if (dest_rank /= irank - 1) cycle

        ! Get the segments in source array that overlap with the destination segment
        source_segments => source%get_segments_for_range( &
                           dest_segment%offset, dest_segment%length, status)

        do isource = 1, size(source_segments)

          source_rank = source_segments(isource)%rank

          if (source_rank == rank) then

            ! Locate the overlapping region between source-segments(i) and
            ! dest_segment
            call get_segment_overlap(source_segments(isource), dest_segment, &
                                     overlap_start, overlap_end)

            source_data => source_segments(isource)%data( &
                           overlap_start - source_segments(isource)%offset + 1: &
                           overlap_end - source_segments(isource)%offset)

            if (dest_rank == rank) then

              ! Configure destination buffer
              dest_data => dest_segment%data( &
                           overlap_start - dest_segment%offset + 1: &
                           overlap_end - dest_segment%offset)

              ! Copy from the source buffer to the destination buffer
              dest_data = source_data
            else
              crsr => crsrs(source_rank + 1)

              ! Copy segment data to the send buffer
              source%buffer( &
                crsr + 1:crsr + overlap_end - overlap_start) = source_data

              crsr = crsr + overlap_end - overlap_start
            end if

          end if

        end do

      end do

    end do

    call MPI_Ialltoallv( &
      source%buffer, source%sendcounts, source%sdspls, MPI_DOUBLE_PRECISION, &
      dest%buffer, source%recvcounts, source%rdspls, MPI_DOUBLE_PRECISION, &
      source%comm, source%req, ierr)

    source%foreign => dest
    dest%foreign => source

    source%transfer_pending = .true.
    source%is_transfer_dest = .false.
    dest%transfer_pending = .true.
    dest%is_transfer_dest = .true.

    if (call_finish) call source%finish_transfer()

  end subroutine transfer_data

  subroutine finish_transfer(source, status)

    !! Finish any pending transfers. Wait for any pending MPI requests to
    !! complete, then transfer data from the receive buffers, free
    !! arrays allocated by transfer_data, and set transfer_pending to .false.
    !! on source and dest darrays.
    !!
    !! Calling this on a darray with no transfers pending is a no-op.

    class(darray), intent(inout), target::source
    type(error_container), intent(out), optional::status
        !! Error status
    class(darray_segment_set), pointer::dest
    integer::crsr ! Cursor into the current receive buffer
    integer::idest, isource
    class(darray_segment), pointer::dest_segment
    ! Pointer to destination segment
    class(darray_segment), pointer::source_segments(:)
    ! Pointer to an array of source segments
    class(darray_segment), pointer::overlapping_source_segments(:)
    ! Pointer to an array of source segments that overlap with a given
    ! destination segment
    real(kind=8), pointer::dest_data(:)
    integer::overlap_start, overlap_end
    integer::rank
    integer::irank
    integer::ierr
    integer::comm_size
    integer::dest_rank, source_rank
        !! Ranks of source and destination processes for a message

    dest => source%foreign

    if (.not. source%transfer_pending) return

    if (source%is_transfer_dest) then
      call throw(status, new_darray_transfer_error( &
                 'finish_transfer should be called on ' &
                 //'the transfer source darray, not the destination darray ' &
                 //'or segment set.'))
      return
    end if

    call mpi_comm_rank(source%comm, rank, ierr)
    call mpi_comm_size(source%comm, comm_size, ierr)

    call MPI_Wait(source%req, MPI_STATUS_IGNORE, ierr)

    crsr = 0

    do irank = 1, comm_size

      do idest = 1, size(dest%segments)

        dest_segment => dest%segments(idest)
        dest_rank = dest_segment%rank

        if (dest_rank == rank) then
          ! Get the segments in this array that overlap with the destination segment
          source_segments => source%get_segments_for_range( &
                             dest_segment%offset, dest_segment%length, status)

          do isource = 1, size(source_segments)

            source_rank = source_segments(isource)%rank

            if (source_rank /= rank .and. source_rank == irank - 1) then

              ! Locate the overlapping region between source-segments(i) and
              ! dest_segment
              call get_segment_overlap(source_segments(isource), dest_segment, &
                                       overlap_start, overlap_end)

              ! Configure destination buffer
              dest_data => dest_segment%data( &
                           overlap_start - dest_segment%offset + 1: &
                           overlap_end - dest_segment%offset)

              ! Copy from the receive buffer to the destination segment
              dest_data = dest%buffer( &
                          crsr + 1:crsr + overlap_end - overlap_start)

              ! Move cursor
              crsr = crsr + overlap_end - overlap_start

            end if

          end do

        end if

      end do

    end do

    deallocate (source%buffer, source%sendcounts, source%recvcounts, &
                source%sdspls, source%rdspls, dest%buffer)

    source%transfer_pending = .false.
    dest%transfer_pending = .false.

    nullify (source%foreign, dest%foreign)

  end subroutine finish_transfer

  function new_darray(segments, comm, status)

    !! Create a new darray

    ! Arguments
    type(darray_segment), intent(in)::segments(:)
        !! Array of darray segments
    MPI_COMM_TYPE, intent(in)::comm
        !! MPI communicator
    type(error_container), intent(out), optional::status
        !! Error status

    type(darray)::new_darray
        !! New darray to create

    integer::i ! Loop counter
    character(:), allocatable::errstr ! Error string
    integer::expected_offset ! Expected offset for segment

    integer::rank !! MPI rank
    integer::ierr !! MPI error code

    call mpi_comm_rank(comm, rank, ierr)

    ! Loop over segments and do sanity checks
    do i = 1, size(segments)

      if (segments(i)%comm /= comm) then

        errstr = 'Segment '//str(i, '(I0)')// &
             ' uses a different communicator than the one &
             &specified for the darray'

        call throw(status, new_exception(errstr, 'new_darray'))

        return

      end if

      if (i == 1) then

        ! Check that first segment is at start of array
        if (segments(i)%offset /= 0) then
          call throw(status, new_exception( &
                     'First segment offset must be 0', &
                     'new_darray'))
          return
        end if

      else

        ! Check that subsequent segments cover the array without gaps

        expected_offset = segments(i - 1)%offset + segments(i - 1)%length

        if (expected_offset /= segments(i)%offset) then
          errstr = 'Expected segment('//str(i)//')%offset='// &
                   str(expected_offset)// &
                   ', got segment('//str(i)//')%offset='// &
                   str(segments(i)%offset)//'.'
          call throw(status, new_exception(errstr, 'new_darray'))
          return
        end if

      end if

      if (segments(i)%rank == rank) then
        if (segments(i)%length /= size(segments(i)%data)) then

          errstr = 'segments('//str(i)//')%length is ' &
                   //str(segments(i)%length)// &
                   ' but segments('//str(i)//')%data has a size of '// &
                   str(size(segments(i)%data))//'.'

          call throw(status, new_exception(errstr, 'new_darray'))

          return
        end if

      end if

    end do

    new_darray%segments = segments
    new_darray%comm = comm

  end function new_darray

  subroutine write_segment_data(this, offset, data, status)
    !! Write to the segment's data array.
    !! Note that the offset must correspond to a location within the extent
    !! of the segment and the length of data must be small enough to stay
    !! within the segment.
    !! This routine must only be called on the processor whose rank corresponds
    !! to `this%rank`.

    class(darray_segment), intent(inout)::this
        !! Distributed array segment
    integer, intent(in)::offset
        !! Offset relative to the start of the global distributed array
    real(kind=8), intent(in)::data(:)
        !! Data to write
    type(error_container), intent(out), optional::status
        !! Error status

    integer::rank !! MPI rank
    integer::ierr !! Error code returned from MPI
    integer::local_start, local_end !! Range on the local array

    call mpi_comm_rank(this%comm, rank, ierr)

    if (rank /= this%rank) then
      call throw(status, new_exception( &
                 'write_segment_data called from a foreign processor rank', &
                 'write_segment_data'))
      return
    end if

    if (offset < this%offset .or. &
        offset + size(data) > this%offset + size(this%data)) then
      call throw(status, &
                 new_exception('Attempt to write to data outside of segment', &
                               'write_segment_data'))
      return
    end if

    ! Determine segment-local array indices
    local_start = offset - this%offset + 1
    local_end = offset - this%offset + size(data)

    ! Copy data to output array
    this%data(local_start:local_end) = data

  end subroutine write_segment_data

  function read_segment_data(this, offset, length, status) result(data)
    !! Read data from the segment starting at `offset` and continuing for
    !! `size` array elements.
    !! Note that the `offset` and `size` must correspond to locations within the
    !! extent of the segment.
    !! This routine must only be called on the processor whose rank corresponds
    !! to `this%rank`.

    ! Arguments
    class(darray_segment)::this
        !! Distributed array segment
    integer, intent(in)::offset
        !! Offset relative to the start of the global distributed array
    integer, intent(in)::length
        !! Number of elements to read
    type(error_container), intent(out), optional::status
        !! Error status

    real(kind=8)::data(length)
        !! Data that was read

    integer::rank !! MPI rank
    integer::ierr !! Error code returned from MPI
    integer::local_start, local_end !! Range on the local array

    character(*), parameter::nameproc = 'read_segment_data'

    call mpi_comm_rank(this%comm, rank, ierr)

    if (rank /= this%rank) then
      call throw(status, new_exception( &
                 'write_segment_data called from a foreign processor rank', &
                 'write_segment_data'))
      return
    end if

    if (offset < this%offset .or. &
        offset + length > this%offset + size(this%data)) then
      call throw(status, &
                 new_exception('Attempt to write to data outside of segment', &
                               'write_segment_data'))
      return
    end if

    ! Determine segment-local array indices
    local_start = offset - this%offset + 1
    local_end = offset - this%offset + length

    ! Copy data to output array
    data = this%data(local_start:local_end)

  end function read_segment_data

  function get_segments_for_range(this, offset, length, status) result(segments)
    !! Get array segments that cover the range of array elements
    !! `offset`+1:`offset`+`length`

    ! Arguments
    class(darray), target::this
        !! Distributed array segment
    integer, intent(in)::offset
        !! Offset relative to the start of the global distributed array
    integer, intent(in)::length
        !! Length of range
    type(error_container), intent(out), optional::status
        !! Error status

    class(darray_segment), pointer::segments(:)
        !! List of segments

    integer::istart, iend ! Start and end indices in this%segments

    istart = this%get_segment_index_for_offset(offset, status)
    iend = this%get_segment_index_for_offset(offset + length, status)

    segments => this%segments(istart:iend)

  end function get_segments_for_range

  function get_segment_index_for_offset(this, offset, status) result(iseg)
    !! Get index of in this%segments for the segment containing the element with offset `offset`

    ! Arguments
    class(darray)::this
        !! Distributed array segment
    integer, intent(in)::offset
        !! Offset relative to the start of the global distributed array
    type(error_container), intent(out), optional::status
        !! Error status

    integer::iseg !! Segment index

    integer::ilower, iupper ! Bounds of search bracket for matching segment
    integer::imid           ! Midpoint of bracket
    integer::offset_min, offset_max ! Range of valid offsets for the array

    character(:), allocatable::errstr ! Error string

    ilower = 1
    iupper = size(this%segments)

    offset_min = this%segments(1)%offset
    offset_max = this%segments(size(this%segments))%offset + &
                 this%segments(size(this%segments))%length

    if (offset < offset_min .or. offset > offset_max) then

      errstr = 'Request for offset '//str(offset)// &
               ' which is out of the range of the array ('//str(offset_min)// &
               '-'//str(offset_max)//')'
      call throw(status, new_exception(errstr, 'get_segment_index_for_offset'))
      return
    end if

    ! Bisection search for the matching segment
    do while (this%segments(ilower)%offset > offset .or. &
              this%segments(ilower)%offset + &
              this%segments(ilower)%length < offset)

      ! Midpoint of bracket
      imid = (ilower + iupper)/2

      ! Reduce bracket size by comparing offset to this%segments(imid)%offset
      if (this%segments(imid)%offset >= offset) iupper = imid
      if (this%segments(imid)%offset <= offset) ilower = imid
      if (this%segments(iupper)%offset <= offset) ilower = iupper

    end do

    iseg = ilower

  end function get_segment_index_for_offset

end module distributed_array
