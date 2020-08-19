#include "mpi_types.h"

module darray_tests
  use system_mpi
  use distributed_array, ONLY: darray, darray_segment, new_darray, transfer_data
  use exceptions, ONLY: throw, error_status, new_exception
  use util, ONLY: str

  implicit none

contains

  function build_darray(arr,np,comm) result(built_darray)

    !! Convert an ordinary array into a distributed array

    real(kind=8),intent(in)::arr(:)
        !! Array to distribute
    integer,intent(in)::np
        !! Number of processors
    MPI_COMM_TYPE,intent(in)::comm
        !! MPI communicator

    type(darray)::built_darray
        !! Built darray

    integer::rank                ! MPI processor rank
    integer::ierr                ! MPI error code

    type(darray_segment)::segments(np)   ! Array segments
    integer::n                    ! Array size

    integer::step                 ! Maximum size of an array segment

    integer::i    ! Loop counter

    call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

    n=size(arr)

    step=n/np+min(1,mod(n,np))
    do i=1,np
       segments(i)%offset=(i-1)*step
       segments(i)%length=min(step,n-segments(i)%offset)
       segments(i)%rank=i-1
       segments(i)%comm=comm
       if(rank==segments(i)%rank) then
          segments(i)%data=arr( &
               segments(i)%offset+1: &
               segments(i)%offset+segments(i)%length)
       end if
    end do

    built_darray%segments=segments
    built_darray%comm=comm

  end function build_darray

  subroutine check_darray_contents(check_darray,arr,status)

    !! Check the contents of a darray

    ! Arguments
    real(kind=8),intent(in),target::arr(:)
        !! Array to distribute
    type(darray),intent(in),target::check_darray
        !! Built darray
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    integer::rank                ! MPI processor rank
    integer::ierr                ! MPI error code

    real(kind=8),pointer::expected_segment_data(:)
        ! Pointer to segment's portion of darray
    type(darray_segment),pointer::segment ! Pointer to a darray segment

    character(:),allocatable::errstr ! Error string

    integer::i    ! Loop counter

    call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

    do i=1,size(check_darray%segments)

       segment=>check_darray%segments(i)

       if(segment%rank==rank) then

          expected_segment_data=>arr(segment%offset+1: &
               segment%offset+segment%length)

          if(any(segment%data /= expected_segment_data)) then

             errstr='Wrong data in segment '//str(i,"(I0)")//'.'

             call throw(status,new_exception(errstr,'check_darray_contents'))

          end if

       end if

    end do

  end subroutine check_darray_contents

  subroutine test_darray_transfer(status)

    !! Test transferring data between darrays with different processor
    !! layouts

    ! Arguments
    class(error_status),intent(out),allocatable,optional::status
        !! Error status

    integer,parameter::n=98     ! Number of array elements
    integer,parameter::np_src=4  ! Number of source processes
    integer,parameter::np_dest=3 ! Number of dest processes
    integer::step                ! Maximum segment size
    real(kind=8)::src_arr(n)     ! Source array data
    real(kind=8)::dest_arr(n)    ! Destination array data
    type(darray)::src_darray     ! Source darray
    type(darray)::dest_darray    ! Destination darray
    integer::rank                ! MPI processor rank
    integer::ierr                ! MPI error code
    integer::nproc               ! Number of processes

    character(:),allocatable::errstr ! Error string

    integer::i    ! Loop counter

    call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

    call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)

    if(nproc<max(np_src,np_dest)) then
       errstr='Not enough processes to run test. Have '//str(nproc)// &
            ' process on MPI communicator, need at least '// &
            str(max(np_src,np_dest))//'.'
       call throw(status,new_exception(errstr,'test_darray_transfer'))
       return
    end if

    do i=1,n
       src_arr(i)=i
       dest_arr(i)=0
    end do

    src_darray=build_darray(src_arr,np_src,MPI_COMM_WORLD)
    dest_darray=build_darray(dest_arr,np_dest,MPI_COMM_WORLD)

    call transfer_data(src_darray,dest_darray)

    call check_darray_contents(dest_darray,src_arr)

  end subroutine test_darray_transfer

end module darray_tests
