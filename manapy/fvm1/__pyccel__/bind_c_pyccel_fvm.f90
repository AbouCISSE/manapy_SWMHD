module bind_c_pyccel_fvm

  use pyccel_fvm, only: explicitscheme_convective_2d
  use pyccel_fvm, only: explicitscheme_dissipative
  use pyccel_fvm, only: explicitscheme_convective_3d

  use, intrinsic :: ISO_C_Binding, only : i64 => C_INT64_T , f64 => &
        C_DOUBLE
  implicit none

  contains

  !........................................
  subroutine bind_c_explicitscheme_dissipative(n0_wx_face, wx_face, &
        n0_wy_face, wy_face, n0_wz_face, wz_face, n0_cellidf, &
        n1_cellidf, cellidf, n0_normalf, n1_normalf, normalf, n0_namef, &
        namef, n0_dissip_w, dissip_w, Dxx, Dyy, Dzz) bind(c)

    implicit none

    integer(i64), value :: n0_wx_face
    real(f64), intent(in) :: wx_face(0:n0_wx_face - 1_i64)
    integer(i64), value :: n0_wy_face
    real(f64), intent(in) :: wy_face(0:n0_wy_face - 1_i64)
    integer(i64), value :: n0_wz_face
    real(f64), intent(in) :: wz_face(0:n0_wz_face - 1_i64)
    integer(i64), value :: n0_cellidf
    integer(i64), value :: n1_cellidf
    integer(i64), intent(inout) :: cellidf(0:n1_cellidf - 1_i64,0: &
          n0_cellidf - 1_i64)
    integer(i64), value :: n0_normalf
    integer(i64), value :: n1_normalf
    real(f64), intent(in) :: normalf(0:n1_normalf - 1_i64,0:n0_normalf - &
          1_i64)
    integer(i64), value :: n0_namef
    integer(i64), intent(in) :: namef(0:n0_namef - 1_i64)
    integer(i64), value :: n0_dissip_w
    real(f64), intent(inout) :: dissip_w(0:n0_dissip_w - 1_i64)
    real(f64), value :: Dxx
    real(f64), value :: Dyy
    real(f64), value :: Dzz

    call explicitscheme_dissipative(wx_face, wy_face, wz_face, cellidf, &
          normalf, namef, dissip_w, Dxx, Dyy, Dzz)

  end subroutine bind_c_explicitscheme_dissipative
  !........................................

  !........................................
  subroutine bind_c_compute_upwind_flux(w_l, w_r, u_face, v_face, w_face &
        , n0_normal, normal, n0_flux_w, flux_w) bind(c)

    implicit none

    real(f64), value :: w_l
    real(f64), value :: w_r
    real(f64), value :: u_face
    real(f64), value :: v_face
    real(f64), value :: w_face
    integer(i64), value :: n0_normal
    real(f64), intent(in) :: normal(0:n0_normal - 1_i64)
    integer(i64), value :: n0_flux_w
    real(f64), intent(inout) :: flux_w(0:n0_flux_w - 1_i64)
    real(f64) :: sol_0009
    real(f64) :: sign_0009

    sol_0009 = 0.0_f64
    sign_0009 = u_face * normal(0_i64) + v_face * normal(1_i64) + w_face &
          * normal(2_i64)
    if (sign_0009 >= 0_i64) then
      sol_0009 = w_l
    else
      sol_0009 = w_r
    end if
    flux_w(0_i64) = sign_0009 * sol_0009

  end subroutine bind_c_compute_upwind_flux
  !........................................

  !........................................
  subroutine bind_c_explicitscheme_convective_2d(n0_rez_w, rez_w, n0_w_c &
        , w_c, n0_w_ghost, w_ghost, n0_w_halo, w_halo, n0_u_face, &
        u_face, n0_v_face, v_face, n0_w_face, w_face, n0_w_x, w_x, &
        n0_w_y, w_y, n0_w_z, w_z, n0_wx_halo, wx_halo, n0_wy_halo, &
        wy_halo, n0_wz_halo, wz_halo, n0_psi, psi, n0_psi_halo, &
        psi_halo, n0_centerc, n1_centerc, centerc, n0_centerf, &
        n1_centerf, centerf, n0_centerh, n1_centerh, centerh, &
        n0_centerg, n1_centerg, centerg, n0_cellidf, n1_cellidf, &
        cellidf, n0_mesuref, mesuref, n0_normalf, n1_normalf, normalf, &
        n0_halofid, halofid, n0_name, name, n0_innerfaces, innerfaces, &
        n0_halofaces, halofaces, n0_boundaryfaces, boundaryfaces, &
        n0_periodicboundaryfaces, periodicboundaryfaces, n0_shift, &
        n1_shift, shift, order) bind(c)

    implicit none

    integer(i64), value :: n0_rez_w
    real(f64), intent(inout) :: rez_w(0:n0_rez_w - 1_i64)
    integer(i64), value :: n0_w_c
    real(f64), intent(in) :: w_c(0:n0_w_c - 1_i64)
    integer(i64), value :: n0_w_ghost
    real(f64), intent(in) :: w_ghost(0:n0_w_ghost - 1_i64)
    integer(i64), value :: n0_w_halo
    real(f64), intent(in) :: w_halo(0:n0_w_halo - 1_i64)
    integer(i64), value :: n0_u_face
    real(f64), intent(in) :: u_face(0:n0_u_face - 1_i64)
    integer(i64), value :: n0_v_face
    real(f64), intent(in) :: v_face(0:n0_v_face - 1_i64)
    integer(i64), value :: n0_w_face
    real(f64), intent(in) :: w_face(0:n0_w_face - 1_i64)
    integer(i64), value :: n0_w_x
    real(f64), intent(in) :: w_x(0:n0_w_x - 1_i64)
    integer(i64), value :: n0_w_y
    real(f64), intent(in) :: w_y(0:n0_w_y - 1_i64)
    integer(i64), value :: n0_w_z
    real(f64), intent(in) :: w_z(0:n0_w_z - 1_i64)
    integer(i64), value :: n0_wx_halo
    real(f64), intent(in) :: wx_halo(0:n0_wx_halo - 1_i64)
    integer(i64), value :: n0_wy_halo
    real(f64), intent(in) :: wy_halo(0:n0_wy_halo - 1_i64)
    integer(i64), value :: n0_wz_halo
    real(f64), intent(in) :: wz_halo(0:n0_wz_halo - 1_i64)
    integer(i64), value :: n0_psi
    real(f64), intent(in) :: psi(0:n0_psi - 1_i64)
    integer(i64), value :: n0_psi_halo
    real(f64), intent(in) :: psi_halo(0:n0_psi_halo - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_centerf
    integer(i64), value :: n1_centerf
    real(f64), intent(in) :: centerf(0:n1_centerf - 1_i64,0:n0_centerf - &
          1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: n0_centerg
    integer(i64), value :: n1_centerg
    real(f64), intent(in) :: centerg(0:n1_centerg - 1_i64,0:n0_centerg - &
          1_i64)
    integer(i64), value :: n0_cellidf
    integer(i64), value :: n1_cellidf
    integer(i64), intent(inout) :: cellidf(0:n1_cellidf - 1_i64,0: &
          n0_cellidf - 1_i64)
    integer(i64), value :: n0_mesuref
    real(f64), intent(in) :: mesuref(0:n0_mesuref - 1_i64)
    integer(i64), value :: n0_normalf
    integer(i64), value :: n1_normalf
    real(f64), intent(in) :: normalf(0:n1_normalf - 1_i64,0:n0_normalf - &
          1_i64)
    integer(i64), value :: n0_halofid
    integer(i64), intent(in) :: halofid(0:n0_halofid - 1_i64)
    integer(i64), value :: n0_name
    integer(i64), intent(in) :: name(0:n0_name - 1_i64)
    integer(i64), value :: n0_innerfaces
    integer(i64), intent(in) :: innerfaces(0:n0_innerfaces - 1_i64)
    integer(i64), value :: n0_halofaces
    integer(i64), intent(in) :: halofaces(0:n0_halofaces - 1_i64)
    integer(i64), value :: n0_boundaryfaces
    integer(i64), intent(in) :: boundaryfaces(0:n0_boundaryfaces - 1_i64 &
          )
    integer(i64), value :: n0_periodicboundaryfaces
    integer(i64), intent(in) :: periodicboundaryfaces(0: &
          n0_periodicboundaryfaces - 1_i64)
    integer(i64), value :: n0_shift
    integer(i64), value :: n1_shift
    real(f64), intent(in) :: shift(0:n1_shift - 1_i64,0:n0_shift - 1_i64 &
          )
    integer(i64), value :: order

    call explicitscheme_convective_2d(rez_w, w_c, w_ghost, w_halo, &
          u_face, v_face, w_face, w_x, w_y, w_z, wx_halo, wy_halo, &
          wz_halo, psi, psi_halo, centerc, centerf, centerh, centerg, &
          cellidf, mesuref, normalf, halofid, name, innerfaces, &
          halofaces, boundaryfaces, periodicboundaryfaces, shift, order &
          )

  end subroutine bind_c_explicitscheme_convective_2d
  !........................................

  !........................................
  subroutine bind_c_explicitscheme_convective_3d(n0_rez_w, rez_w, n0_w_c &
        , w_c, n0_w_ghost, w_ghost, n0_w_halo, w_halo, n0_u_face, &
        u_face, n0_v_face, v_face, n0_w_face, w_face, n0_w_x, w_x, &
        n0_w_y, w_y, n0_w_z, w_z, n0_wx_halo, wx_halo, n0_wy_halo, &
        wy_halo, n0_wz_halo, wz_halo, n0_psi, psi, n0_psi_halo, &
        psi_halo, n0_centerc, n1_centerc, centerc, n0_centerf, &
        n1_centerf, centerf, n0_centerh, n1_centerh, centerh, &
        n0_centerg, n1_centerg, centerg, n0_cellidf, n1_cellidf, &
        cellidf, n0_mesuref, mesuref, n0_normalf, n1_normalf, normalf, &
        n0_halofid, halofid, n0_name, name, n0_innerfaces, innerfaces, &
        n0_halofaces, halofaces, n0_boundaryfaces, boundaryfaces, &
        n0_periodicboundaryfaces, periodicboundaryfaces, n0_shift, &
        n1_shift, shift, order) bind(c)

    implicit none

    integer(i64), value :: n0_rez_w
    real(f64), intent(inout) :: rez_w(0:n0_rez_w - 1_i64)
    integer(i64), value :: n0_w_c
    real(f64), intent(in) :: w_c(0:n0_w_c - 1_i64)
    integer(i64), value :: n0_w_ghost
    real(f64), intent(in) :: w_ghost(0:n0_w_ghost - 1_i64)
    integer(i64), value :: n0_w_halo
    real(f64), intent(in) :: w_halo(0:n0_w_halo - 1_i64)
    integer(i64), value :: n0_u_face
    real(f64), intent(in) :: u_face(0:n0_u_face - 1_i64)
    integer(i64), value :: n0_v_face
    real(f64), intent(in) :: v_face(0:n0_v_face - 1_i64)
    integer(i64), value :: n0_w_face
    real(f64), intent(in) :: w_face(0:n0_w_face - 1_i64)
    integer(i64), value :: n0_w_x
    real(f64), intent(in) :: w_x(0:n0_w_x - 1_i64)
    integer(i64), value :: n0_w_y
    real(f64), intent(in) :: w_y(0:n0_w_y - 1_i64)
    integer(i64), value :: n0_w_z
    real(f64), intent(in) :: w_z(0:n0_w_z - 1_i64)
    integer(i64), value :: n0_wx_halo
    real(f64), intent(in) :: wx_halo(0:n0_wx_halo - 1_i64)
    integer(i64), value :: n0_wy_halo
    real(f64), intent(in) :: wy_halo(0:n0_wy_halo - 1_i64)
    integer(i64), value :: n0_wz_halo
    real(f64), intent(in) :: wz_halo(0:n0_wz_halo - 1_i64)
    integer(i64), value :: n0_psi
    real(f64), intent(in) :: psi(0:n0_psi - 1_i64)
    integer(i64), value :: n0_psi_halo
    real(f64), intent(in) :: psi_halo(0:n0_psi_halo - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_centerf
    integer(i64), value :: n1_centerf
    real(f64), intent(in) :: centerf(0:n1_centerf - 1_i64,0:n0_centerf - &
          1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: n0_centerg
    integer(i64), value :: n1_centerg
    real(f64), intent(in) :: centerg(0:n1_centerg - 1_i64,0:n0_centerg - &
          1_i64)
    integer(i64), value :: n0_cellidf
    integer(i64), value :: n1_cellidf
    integer(i64), intent(inout) :: cellidf(0:n1_cellidf - 1_i64,0: &
          n0_cellidf - 1_i64)
    integer(i64), value :: n0_mesuref
    real(f64), intent(in) :: mesuref(0:n0_mesuref - 1_i64)
    integer(i64), value :: n0_normalf
    integer(i64), value :: n1_normalf
    real(f64), intent(in) :: normalf(0:n1_normalf - 1_i64,0:n0_normalf - &
          1_i64)
    integer(i64), value :: n0_halofid
    integer(i64), intent(in) :: halofid(0:n0_halofid - 1_i64)
    integer(i64), value :: n0_name
    integer(i64), intent(in) :: name(0:n0_name - 1_i64)
    integer(i64), value :: n0_innerfaces
    integer(i64), intent(in) :: innerfaces(0:n0_innerfaces - 1_i64)
    integer(i64), value :: n0_halofaces
    integer(i64), intent(in) :: halofaces(0:n0_halofaces - 1_i64)
    integer(i64), value :: n0_boundaryfaces
    integer(i64), intent(in) :: boundaryfaces(0:n0_boundaryfaces - 1_i64 &
          )
    integer(i64), value :: n0_periodicboundaryfaces
    integer(i64), intent(in) :: periodicboundaryfaces(0: &
          n0_periodicboundaryfaces - 1_i64)
    integer(i64), value :: n0_shift
    integer(i64), value :: n1_shift
    real(f64), intent(in) :: shift(0:n1_shift - 1_i64,0:n0_shift - 1_i64 &
          )
    integer(i64), value :: order

    call explicitscheme_convective_3d(rez_w, w_c, w_ghost, w_halo, &
          u_face, v_face, w_face, w_x, w_y, w_z, wx_halo, wy_halo, &
          wz_halo, psi, psi_halo, centerc, centerf, centerh, centerg, &
          cellidf, mesuref, normalf, halofid, name, innerfaces, &
          halofaces, boundaryfaces, periodicboundaryfaces, shift, order &
          )

  end subroutine bind_c_explicitscheme_convective_3d
  !........................................

end module bind_c_pyccel_fvm
