module util_tests
  use util
  implicit none
contains
  subroutine test_append_array
    real(kind=8),allocatable::a(:),b(:)
    real(kind=8)::x
    integer,parameter::n=10
    integer::i

    allocate(a(n))

    do i=1,n
       a(i)=i
    end do

    x=n+1

    b=append_array(a,x)

    do i=1,n
       if(b(i)/=a(i)) then
          print *,'Error: b(',i,')=',b(i),', expected b(',i,')=',a(i)
          error stop
       end if
    end do

    if(b(n+1)/=x) then
       print *,'Error: b(',n+1,')=',b(n+1),', expected b(',n+1,')=',x
       error stop
    end if

  end subroutine test_append_array
end module util_tests

program test_util
  use util_tests
  implicit none

  call test_append_array
end program test_util
