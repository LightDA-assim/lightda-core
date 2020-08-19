#include "mpi_types.h"

module util
  use system_mpi
  use exceptions, ONLY: error_status,throw,new_exception

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

  function str(x,fmt,status)
    class(*),intent(in)::x
    character(*),optional,intent(in)::fmt
    class(error_status),intent(out),allocatable,optional::status
        !! Error status
    character(30)::str_tmp
    character(:),allocatable::str

    select type(x)
    type is(integer)
       if(present(fmt)) then
          write(str_tmp,fmt) x
       else
          write(str_tmp,'(I0)') x
       end if
    class default
       call throw(status, &
            new_exception('Variable of unknown type passed to str()','str'))
       return
    end select

    str=trim(str_tmp)

  end function str

end module util
