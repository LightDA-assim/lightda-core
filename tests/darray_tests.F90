#include "mpi_types.h"

module darray_tests
  use system_mpi
  use distributed_array, ONLY: darray, darray_segment, darray_segment_set, &
       new_darray, transfer_data, darray_transfer_error
  use exceptions, ONLY: throw, error_container, new_exception, transfer_error
  use util, ONLY: str

  implicit none

contains

  function build_darray(arr, np, comm) result(built_darray)

    !! Convert an ordinary array into a distributed array

    real(kind=8), intent(in)::arr(:)
        !! Array to distribute
    integer, intent(in)::np
        !! Number of processors
    MPI_COMM_TYPE, intent(in)::comm
        !! MPI communicator

    type(darray)::built_darray
        !! Built darray

    integer::rank                ! MPI processor rank
    integer::ierr                ! MPI error code

    type(darray_segment)::segments(np)   ! Array segments
    integer::n                    ! Array size

    integer::step                 ! Maximum size of an array segment

    integer::i    ! Loop counter

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    n = size(arr)

    step = n/np + min(1, mod(n, np))
    do i = 1, np
      segments(i)%offset = (i - 1)*step
      segments(i)%length = min(step, n - segments(i)%offset)
      segments(i)%rank = i - 1
      segments(i)%comm = comm
      if (rank == segments(i)%rank) then
        segments(i)%data = arr( &
                           segments(i)%offset + 1: &
                           segments(i)%offset + segments(i)%length)
      end if
    end do

    built_darray = new_darray(segments, comm)

  end function build_darray

  subroutine check_darray_contents(check_darray, arr, status)

    !! Check the contents of a darray

    ! Arguments
    real(kind=8), intent(in), target::arr(:)
        !! Array to distribute
    class(darray_segment_set), intent(in), target::check_darray
        !! Built darray
    type(error_container), intent(out), optional::status
        !! Error status

    integer::rank                ! MPI processor rank
    integer::ierr                ! MPI error code

    real(kind=8), pointer::expected_segment_data(:)
    ! Pointer to segment's portion of darray
    type(darray_segment), pointer::segment ! Pointer to a darray segment

    character(:), allocatable::errstr ! Error string

    integer::i    ! Loop counter

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    do i = 1, size(check_darray%segments)

      segment => check_darray%segments(i)

      if (segment%rank == rank) then

        expected_segment_data => arr(segment%offset + 1: &
                                     segment%offset + segment%length)

        if (any(segment%data /= expected_segment_data)) then

          errstr = 'Wrong data in segment '//str(i, "(I0)")//'.'

          call throw(status, new_exception(errstr, 'check_darray_contents'))

        end if

      end if

    end do

  end subroutine check_darray_contents

  subroutine test_darray_transfer(status)

    !! Test transferring data between darrays with different processor
    !! layouts

    ! Arguments
    type(error_container), intent(out), optional::status
        !! Error status

    integer, parameter::n = 98     ! Number of array elements
    integer, parameter::np_src = 4  ! Number of source processes
    integer, parameter::np_dest = 3 ! Number of dest processes
    integer::step                      ! Maximum segment size
    real(kind=8)::src_arr(n)           ! Source array data
    real(kind=8)::dest_arr(n)          ! Destination array data
    type(darray)::src_darray           ! Source darray
    type(darray)::dest_darray          ! Destination darray
    type(darray_segment_set)::dest_set ! Destination darray segment set
    type(darray_segment), allocatable::segments(:)
                                       ! Segments to go in dest_set
    integer::rank                      ! MPI processor rank
    integer::ierr                      ! MPI error code
    integer::nproc                     ! Number of processes

    character(:), allocatable::errstr ! Error string

    integer::i    ! Loop counter

    type(error_container)::child_status
        !! Error status of called procedures

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    if (nproc < max(np_src, np_dest)) then
      errstr = 'Not enough processes to run test. Have '//str(nproc)// &
               ' process on MPI communicator, need at least '// &
               str(max(np_src, np_dest))//'.'
      call throw(status, new_exception(errstr, 'test_darray_transfer'))
      return
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Test transferring to dest_darray
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, n
      src_arr(i) = i
      dest_arr(i) = 0
    end do

    src_darray = build_darray(src_arr, np_src, MPI_COMM_WORLD)
    dest_darray = build_darray(dest_arr, np_dest, MPI_COMM_WORLD)

    call transfer_data(src_darray, dest_darray)

    ! Call finish_transfer (should be a no-op since finish_immediately defaults to .false.)
    call src_darray%finish_transfer()
    call dest_darray%finish_transfer()

    call check_darray_contents(dest_darray, src_arr)

    do i = 1, n
      src_arr(i) = i
      dest_arr(i) = 0
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Test deferred call to finish_transfer
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, n
      src_arr(i) = i
      dest_arr(i) = 0
    end do

    src_darray = build_darray(src_arr, np_src, MPI_COMM_WORLD)
    dest_darray = build_darray(dest_arr, np_dest, MPI_COMM_WORLD)

    call transfer_data(src_darray, dest_darray, finish_immediately = .false.)

    ! Call finish_transfer on the destination darray (should throw a darray_transfer_error)
    call dest_darray%finish_transfer(child_status)

    ! Check the exception thrown by finish_transfer
    select type(error => child_status%info)
    class is (darray_transfer_error)
       ! We got the expected exception, mark it as handled
       error%handled = .true.
    class default
       call throw(status, new_exception('Did not catch a darray_transfer_error as expected.'))
       return
    end select

    ! Call finish_transfer on the source darray
    call src_darray%finish_transfer(status)

    call check_darray_contents(dest_darray, src_arr)

    do i = 1, n
      src_arr(i) = i
      dest_arr(i) = 0
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Test transferring back to src_darray
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    src_darray = build_darray(dest_arr, np_src, MPI_COMM_WORLD)

    call transfer_data(dest_darray, src_darray)

    call check_darray_contents(src_darray, src_arr)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Test transferring to a segment set with overlapping segment ranges
    ! and non-increasing ranks
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dest_arr = 0

    allocate (segments(2))
    segments(1)%offset = 60
    segments(1)%length = 10
    segments(1)%rank = 1
    segments(1)%comm = MPI_COMM_WORLD
    segments(1)%data = dest_arr(1:segments(1)%length)
    segments(2)%offset = 40
    segments(2)%length = 30
    segments(2)%rank = 0
    segments(2)%comm = MPI_COMM_WORLD
    segments(2)%data = dest_arr(1:segments(2)%length)
    dest_set%segments = segments
    dest_set%comm = MPI_COMM_WORLD

    call transfer_data(src_darray, dest_set)

    call check_darray_contents(dest_set, src_arr)

  end subroutine test_darray_transfer

end module darray_tests
