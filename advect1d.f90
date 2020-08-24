module advect1d

  use iso_c_binding

contains

  function limiter(sm, sp, limiter_type) result(s)
    implicit none
    real(kind=8), intent(in)::sm, sp
    real(kind=8) :: s
    real(kind=8) :: epsilon
    integer, intent(in)::limiter_type

    epsilon = 1d-12

    s = 0

    select case (limiter_type)

    case (0) ! First order upwind
      s = 0

    case (1) ! Lax Wendroff
      s = sp

    case (2) ! MinMod
      s = max(0., min(sm, sp))*sign(1.0_8, sm + sp)

    case (3) ! Harmonic
      s = (sm*abs(sp) + abs(sm)*sp)/(abs(sm + abs(sp) + epsilon))

    case (4) ! Geometric
      s = sqrt((sm*sp + abs(sm*sp))/2)*sign(1.0_8 + 0, sm + sp)

    case (5) !superbee
      s = max(0., max(min(2*abs(sm), abs(sp)), min(abs(sm), 2*abs(sp)))) &
          *sign(1.0_8, sm + sp)

    end select

  end function limiter

  subroutine advect1d_step(u, a, dx, dt, limiter_type, nx) bind(c)
    implicit none

    real(c_double), dimension(nx), intent(inout) :: u, a
    real(c_double), value :: dx, dt
    real(kind=8) :: fr, fl

    integer(c_int32_t), value :: limiter_type, nx
    integer::i

    call flux(2, fr, u, a, dx, dt, limiter_type)

    do i = 2, size(a) - 1

      fl = fr

      call flux(i + 1, fr, u, a, dx, dt, limiter_type)

      u(i) = u(i) + dt/dx*(fl - fr)

    end do

  end subroutine advect1d_step

  subroutine flux(i, f, u, a, dx, dt, limiter_type)
    implicit none

    real(kind=8), dimension(:) :: u, a
    real(kind=8) :: sm, sp, ap, am, fp, fm, dx, dt, f
    integer, intent(in) :: i, limiter_type

    sm = (u(i) - u(i - 1))/dx

    sp = (u(i + 1) - u(i))/dx

    ap = max(a(i - 1), 0.)
    am = min(a(i), 0.)

    fp = ap*(u(i - 1) &
             + dx/2.0*(1.0 - a(i - 1)*dt/dx)*limiter(sm, sp, limiter_type))

    fm = am*(u(i) &
             - dx/2.0*(1.0 - abs(a(i))*dt/dx)*limiter(sm, sp, limiter_type))

    f = fp + fm

  end subroutine flux

end module advect1d
