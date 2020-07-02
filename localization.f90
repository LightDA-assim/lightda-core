module localization

  implicit none

contains

  function gaspari_cohn_mid(z,c) result(f)
    real(kind=8),intent(in)::z,c
    real(kind=8)::f
    f=1./12*(z/c)**5 - 0.5*(z/c)**4 + 5./8*(z/c)**3 &
        + 5./3*(z/c)**2 - 5*z/c - 2./3*c/z + 4
  end function gaspari_cohn_mid

  function gaspari_cohn_close(z,c) result(f)
    real(kind=8),intent(in)::z,c
    real(kind=8)::f
    f=-0.25*(z/c)**5 + 0.5*(z/c)**4 + 5./8*(z/c)**3 - 5./3*(z/c)**2 + 1
  end function gaspari_cohn_close


  function localize_gaspari_cohn(z,c) result(f)
    real(kind=8),intent(in)::z,c
    real(kind=8)::f

    if(z==0) then
       f=1
    else if(z<=c) then
       f=gaspari_cohn_close(z,c)
    else if(z<=2*c) then
       f=max(gaspari_cohn_mid(z,c),0.0)
    else if(z<0) then
       print *,'Error: Negative distance passed to localize_gaspari_cohn'
       error stop
    else
       f=0
    end if

  end function localize_gaspari_cohn

end module localization
