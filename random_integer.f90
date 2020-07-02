module random_integer

  implicit none

contains

  function randint(a,b)

    integer,intent(in)::a
    integer,intent(in),optional::b
    integer::rangemin
    integer::rangemax
    real::r
    integer::randint

    if(present(b)) then
       rangemin=a
       rangemax=b
    else
       rangemin=1
       rangemax=a
    end if

    call random_number(r)

    randint = int(r*(rangemax-rangemin-1))+rangemin

  end function randint
end module random_integer
