module mod_advect1d_forward_operator

  use forward_operator, ONLY: base_forward_operator

  implicit none

  type, extends(base_forward_operator) :: advect1d_forward_operator

  end type advect1d_forward_operator

end module mod_advect1d_forward_operator
