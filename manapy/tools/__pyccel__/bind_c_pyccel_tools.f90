module bind_c_pyccel_tools

  use pyccel_tools, only: time_step
  use pyccel_tools, only: initialisation_gaussian_3d
  use pyccel_tools, only: update_new_value
  use pyccel_tools, only: initialisation_gaussian_2d

  use, intrinsic :: ISO_C_Binding, only : i64 => C_INT64_T , f64 => &
        C_DOUBLE
  implicit none

  contains

  !........................................
  subroutine bind_c_initialisation_gaussian_2d(n0_ne, ne, n0_u, u, n0_v, &
        v, n0_w, w, n0_P, P, n0_center, n1_center, center, Pinit) bind( &
        c)

    implicit none

    integer(i64), value :: n0_ne
    real(f64), intent(inout) :: ne(0:n0_ne - 1_i64)
    integer(i64), value :: n0_u
    real(f64), intent(inout) :: u(0:n0_u - 1_i64)
    integer(i64), value :: n0_v
    real(f64), intent(inout) :: v(0:n0_v - 1_i64)
    integer(i64), value :: n0_w
    real(f64), intent(inout) :: w(0:n0_w - 1_i64)
    integer(i64), value :: n0_P
    real(f64), intent(inout) :: P(0:n0_P - 1_i64)
    integer(i64), value :: n0_center
    integer(i64), value :: n1_center
    real(f64), intent(in) :: center(0:n1_center - 1_i64,0:n0_center - &
          1_i64)
    real(f64), value :: Pinit

    call initialisation_gaussian_2d(ne, u, v, w, P, center, Pinit)

  end subroutine bind_c_initialisation_gaussian_2d
  !........................................

  !........................................
  subroutine bind_c_initialisation_gaussian_3d(n0_ne, ne, n0_u, u, n0_v, &
        v, n0_w, w, n0_P, P, n0_center, n1_center, center, Pinit) bind( &
        c)

    implicit none

    integer(i64), value :: n0_ne
    real(f64), intent(inout) :: ne(0:n0_ne - 1_i64)
    integer(i64), value :: n0_u
    real(f64), intent(inout) :: u(0:n0_u - 1_i64)
    integer(i64), value :: n0_v
    real(f64), intent(inout) :: v(0:n0_v - 1_i64)
    integer(i64), value :: n0_w
    real(f64), intent(inout) :: w(0:n0_w - 1_i64)
    integer(i64), value :: n0_P
    real(f64), intent(inout) :: P(0:n0_P - 1_i64)
    integer(i64), value :: n0_center
    integer(i64), value :: n1_center
    real(f64), intent(in) :: center(0:n1_center - 1_i64,0:n0_center - &
          1_i64)
    real(f64), value :: Pinit

    call initialisation_gaussian_3d(ne, u, v, w, P, center, Pinit)

  end subroutine bind_c_initialisation_gaussian_3d
  !........................................

  !........................................
  function bind_c_time_step(n0_u, u, n0_v, v, n0_w, w, cfl, n0_normal, &
        n1_normal, normal, n0_mesure, mesure, n0_volume, volume, &
        n0_faceid, n1_faceid, faceid, dim, Dxx, Dyy, Dzz) bind(c) &
        result(dt)

    implicit none

    integer(i64), value :: n0_u
    real(f64), intent(in) :: u(0:n0_u - 1_i64)
    integer(i64), value :: n0_v
    real(f64), intent(in) :: v(0:n0_v - 1_i64)
    integer(i64), value :: n0_w
    real(f64), intent(in) :: w(0:n0_w - 1_i64)
    real(f64), value :: cfl
    integer(i64), value :: n0_normal
    integer(i64), value :: n1_normal
    real(f64), intent(in) :: normal(0:n1_normal - 1_i64,0:n0_normal - &
          1_i64)
    integer(i64), value :: n0_mesure
    real(f64), intent(in) :: mesure(0:n0_mesure - 1_i64)
    integer(i64), value :: n0_volume
    real(f64), intent(in) :: volume(0:n0_volume - 1_i64)
    integer(i64), value :: n0_faceid
    integer(i64), value :: n1_faceid
    integer(i64), intent(in) :: faceid(0:n1_faceid - 1_i64,0:n0_faceid - &
          1_i64)
    integer(i64), value :: dim
    real(f64), value :: Dxx
    real(f64), value :: Dyy
    real(f64), value :: Dzz
    real(f64) :: dt

    dt = time_step(u, v, w, cfl, normal, mesure, volume, faceid, dim, &
          Dxx, Dyy, Dzz)

  end function bind_c_time_step
  !........................................

  !........................................
  subroutine bind_c_update_new_value(n0_ne_c, ne_c, n0_u_c, u_c, n0_v_c, &
        v_c, n0_P_c, P_c, n0_rez_ne, rez_ne, n0_dissip_ne, dissip_ne, &
        n0_src_ne, src_ne, dtime, n0_vol, vol) bind(c)

    implicit none

    integer(i64), value :: n0_ne_c
    real(f64), intent(inout) :: ne_c(0:n0_ne_c - 1_i64)
    integer(i64), value :: n0_u_c
    real(f64), intent(in) :: u_c(0:n0_u_c - 1_i64)
    integer(i64), value :: n0_v_c
    real(f64), intent(in) :: v_c(0:n0_v_c - 1_i64)
    integer(i64), value :: n0_P_c
    real(f64), intent(in) :: P_c(0:n0_P_c - 1_i64)
    integer(i64), value :: n0_rez_ne
    real(f64), intent(in) :: rez_ne(0:n0_rez_ne - 1_i64)
    integer(i64), value :: n0_dissip_ne
    real(f64), intent(in) :: dissip_ne(0:n0_dissip_ne - 1_i64)
    integer(i64), value :: n0_src_ne
    real(f64), intent(in) :: src_ne(0:n0_src_ne - 1_i64)
    real(f64), value :: dtime
    integer(i64), value :: n0_vol
    real(f64), intent(in) :: vol(0:n0_vol - 1_i64)

    call update_new_value(ne_c, u_c, v_c, P_c, rez_ne, dissip_ne, src_ne &
          , dtime, vol)

  end subroutine bind_c_update_new_value
  !........................................

end module bind_c_pyccel_tools
