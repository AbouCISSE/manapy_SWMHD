module pyccel_fvm


  use, intrinsic :: ISO_C_Binding, only : i64 => C_INT64_T , f64 => &
        C_DOUBLE
  implicit none

  contains

  !........................................
  subroutine explicitscheme_dissipative(wx_face, wy_face, wz_face, &
        cellidf, normalf, namef, dissip_w, Dxx, Dyy, Dzz)

    implicit none

    real(f64), intent(in) :: wx_face(0:)
    real(f64), intent(in) :: wy_face(0:)
    real(f64), intent(in) :: wz_face(0:)
    integer(i64), intent(inout) :: cellidf(0:,0:)
    real(f64), intent(in) :: normalf(0:,0:)
    integer(i64), intent(in) :: namef(0:)
    real(f64), intent(inout) :: dissip_w(0:)
    real(f64), value :: Dxx
    real(f64), value :: Dyy
    real(f64), value :: Dzz
    integer(i64) :: nbface
    real(f64) :: norm(0:2_i64)
    integer(i64) :: i
    real(f64) :: q
    real(f64) :: flux_w

    nbface = size(cellidf,2,i64)
    norm = 0.0_f64
    dissip_w(:) = 0.0_f64
    do i = 0_i64, nbface - 1_i64, 1_i64
      norm(:) = normalf(:, i)
      q = Dxx * wx_face(i) * norm(0_i64) + Dyy * wy_face(i) * norm(1_i64 &
            ) + Dzz * wz_face(i) * norm(2_i64)
      flux_w = q
      if (namef(i) == 0_i64) then
        dissip_w(cellidf(0_i64, i)) = dissip_w(cellidf(0_i64, i)) + &
              flux_w
        dissip_w(cellidf(1_i64, i)) = dissip_w(cellidf(1_i64, i)) - &
              flux_w
      else
        dissip_w(cellidf(0_i64, i)) = dissip_w(cellidf(0_i64, i)) + &
              flux_w
      end if
    end do

  end subroutine explicitscheme_dissipative
  !........................................

  !........................................
  !........................................

  !........................................
  subroutine explicitscheme_convective_2d(rez_w, w_c, w_ghost, w_halo, &
        u_face, v_face, w_face, w_x, w_y, w_z, wx_halo, wy_halo, &
        wz_halo, psi, psi_halo, centerc, centerf, centerh, centerg, &
        cellidf, mesuref, normalf, halofid, name, innerfaces, halofaces &
        , boundaryfaces, periodicboundaryfaces, shift, order)

    implicit none

    real(f64), intent(inout) :: rez_w(0:)
    real(f64), intent(in) :: w_c(0:)
    real(f64), intent(in) :: w_ghost(0:)
    real(f64), intent(in) :: w_halo(0:)
    real(f64), intent(in) :: u_face(0:)
    real(f64), intent(in) :: v_face(0:)
    real(f64), intent(in) :: w_face(0:)
    real(f64), intent(in) :: w_x(0:)
    real(f64), intent(in) :: w_y(0:)
    real(f64), intent(in) :: w_z(0:)
    real(f64), intent(in) :: wx_halo(0:)
    real(f64), intent(in) :: wy_halo(0:)
    real(f64), intent(in) :: wz_halo(0:)
    real(f64), intent(in) :: psi(0:)
    real(f64), intent(in) :: psi_halo(0:)
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: centerf(0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    real(f64), intent(in) :: centerg(0:,0:)
    integer(i64), intent(inout) :: cellidf(0:,0:)
    real(f64), intent(in) :: mesuref(0:)
    real(f64), intent(in) :: normalf(0:,0:)
    integer(i64), intent(in) :: halofid(0:)
    integer(i64), intent(in) :: name(0:)
    integer(i64), intent(in) :: innerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: boundaryfaces(0:)
    integer(i64), intent(in) :: periodicboundaryfaces(0:)
    real(f64), intent(in) :: shift(0:,0:)
    integer(i64), value :: order
    real(f64) :: center_left(0:1_i64)
    real(f64) :: center_right(0:1_i64)
    real(f64) :: r_l(0:1_i64)
    real(f64) :: r_r(0:1_i64)
    real(f64) :: normal(0:2_i64)
    real(f64) :: flux_w(0:0_i64)
    integer(i64) :: i
    real(f64) :: w_l
    real(f64) :: w_r
    real(f64) :: w_x_left
    real(f64) :: w_x_right
    real(f64) :: w_y_left
    real(f64) :: w_y_right
    real(f64) :: psi_left
    real(f64) :: psi_right
    integer(i64) :: Dummy_0009
    real(f64) :: sol_0001
    real(f64) :: sign_0001
    integer(i64) :: Dummy_0010
    real(f64) :: sol_0002
    real(f64) :: sign_0002
    integer(i64) :: Dummy_0011
    real(f64) :: sol_0003
    real(f64) :: sign_0003
    integer(i64) :: Dummy_0012
    real(f64) :: sol_0004
    real(f64) :: sign_0004

    center_left = 0.0_f64
    center_right = 0.0_f64
    r_l = 0.0_f64
    r_r = 0.0_f64
    normal = 0.0_f64
    flux_w = 0.0_f64
    rez_w(:) = 0.0_f64
    do Dummy_0009 = 0_i64, size(innerfaces, kind=i64) - 1_i64, 1_i64
      i = innerfaces(Dummy_0009)
      w_l = w_c(cellidf(0_i64, i))
      normal(:) = normalf(:, i)
      w_r = w_c(cellidf(1_i64, i))
      center_left(:) = centerc(0_i64:1_i64, cellidf(0_i64, i))
      center_right(:) = centerc(0_i64:1_i64, cellidf(1_i64, i))
      w_x_left = w_x(cellidf(0_i64, i))
      w_x_right = w_x(cellidf(1_i64, i))
      w_y_left = w_y(cellidf(0_i64, i))
      w_y_right = w_y(cellidf(1_i64, i))
      psi_left = psi(cellidf(0_i64, i))
      psi_right = psi(cellidf(1_i64, i))
      r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
      r_r(0_i64) = centerf(0_i64, i) - center_right(0_i64)
      r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
      r_r(1_i64) = centerf(1_i64, i) - center_right(1_i64)
      w_l = w_l + (order - 1_i64) * psi_left * (w_x_left * r_l(0_i64) + &
            w_y_left * r_l(1_i64))
      w_r = w_r + (order - 1_i64) * psi_right * (w_x_right * r_r(0_i64) &
            + w_y_right * r_r(1_i64))
      sol_0001 = 0.0_f64
      sign_0001 = u_face(i) * normal(0_i64) + v_face(i) * normal(1_i64) &
            + w_face(i) * normal(2_i64)
      if (sign_0001 >= 0_i64) then
        sol_0001 = w_l
      else
        sol_0001 = w_r
      end if
      flux_w(0_i64) = sign_0001 * sol_0001
      rez_w(cellidf(0_i64, i)) = rez_w(cellidf(0_i64, i)) - flux_w(0_i64 &
            )
      rez_w(cellidf(1_i64, i)) = rez_w(cellidf(1_i64, i)) + flux_w(0_i64 &
            )
    end do
    do Dummy_0010 = 0_i64, size(periodicboundaryfaces, kind=i64) - 1_i64 &
          , 1_i64
      i = periodicboundaryfaces(Dummy_0010)
      w_l = w_c(cellidf(0_i64, i))
      normal(:) = normalf(:, i)
      w_r = w_c(cellidf(1_i64, i))
      center_left(:) = centerc(0_i64:1_i64, cellidf(0_i64, i))
      center_right(:) = centerc(0_i64:1_i64, cellidf(1_i64, i))
      w_x_left = w_x(cellidf(0_i64, i))
      w_x_right = w_x(cellidf(1_i64, i))
      w_y_left = w_y(cellidf(0_i64, i))
      w_y_right = w_y(cellidf(1_i64, i))
      psi_left = psi(cellidf(0_i64, i))
      psi_right = psi(cellidf(1_i64, i))
      if (name(i) == 11_i64 .or. name(i) == 22_i64) then
        r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
        r_r(0_i64) = centerf(0_i64, i) - center_right(0_i64) - shift( &
              0_i64, cellidf(1_i64, i))
        r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
        r_r(1_i64) = centerf(1_i64, i) - center_right(1_i64)
      end if
      if (name(i) == 33_i64 .or. name(i) == 44_i64) then
        r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
        r_r(0_i64) = centerf(0_i64, i) - center_right(0_i64)
        r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
        r_r(1_i64) = centerf(1_i64, i) - center_right(1_i64) - shift( &
              1_i64, cellidf(1_i64, i))
      end if
      w_l = w_l + (order - 1_i64) * psi_left * (w_x_left * r_l(0_i64) + &
            w_y_left * r_l(1_i64))
      w_r = w_r + (order - 1_i64) * psi_right * (w_x_right * r_r(0_i64) &
            + w_y_right * r_r(1_i64))
      sol_0002 = 0.0_f64
      sign_0002 = u_face(i) * normal(0_i64) + v_face(i) * normal(1_i64) &
            + w_face(i) * normal(2_i64)
      if (sign_0002 >= 0_i64) then
        sol_0002 = w_l
      else
        sol_0002 = w_r
      end if
      flux_w(0_i64) = sign_0002 * sol_0002
      rez_w(cellidf(0_i64, i)) = rez_w(cellidf(0_i64, i)) - flux_w(0_i64 &
            )
    end do
    do Dummy_0011 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0011)
      w_l = w_c(cellidf(0_i64, i))
      normal(:) = normalf(:, i)
      w_r = w_halo(halofid(i))
      center_left(:) = centerc(0_i64:1_i64, cellidf(0_i64, i))
      center_right(:) = centerh(0_i64:1_i64, halofid(i))
      w_x_left = w_x(cellidf(0_i64, i))
      w_x_right = wx_halo(halofid(i))
      w_y_left = w_y(cellidf(0_i64, i))
      w_y_right = wy_halo(halofid(i))
      psi_left = psi(cellidf(0_i64, i))
      psi_right = psi_halo(halofid(i))
      r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
      r_r(0_i64) = centerf(0_i64, i) - center_right(0_i64)
      r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
      r_r(1_i64) = centerf(1_i64, i) - center_right(1_i64)
      w_l = w_l + (order - 1_i64) * psi_left * (w_x_left * r_l(0_i64) + &
            w_y_left * r_l(1_i64))
      w_r = w_r + (order - 1_i64) * psi_right * (w_x_right * r_r(0_i64) &
            + w_y_right * r_r(1_i64))
      sol_0003 = 0.0_f64
      sign_0003 = u_face(i) * normal(0_i64) + v_face(i) * normal(1_i64) &
            + w_face(i) * normal(2_i64)
      if (sign_0003 >= 0_i64) then
        sol_0003 = w_l
      else
        sol_0003 = w_r
      end if
      flux_w(0_i64) = sign_0003 * sol_0003
      rez_w(cellidf(0_i64, i)) = rez_w(cellidf(0_i64, i)) - flux_w(0_i64 &
            )
    end do
    do Dummy_0012 = 0_i64, size(boundaryfaces, kind=i64) - 1_i64, 1_i64
      i = boundaryfaces(Dummy_0012)
      w_l = w_c(cellidf(0_i64, i))
      normal(:) = normalf(:, i)
      w_r = w_ghost(i)
      center_left(:) = centerc(0_i64:1_i64, cellidf(0_i64, i))
      w_x_left = w_x(cellidf(0_i64, i))
      w_y_left = w_y(cellidf(0_i64, i))
      psi_left = psi(cellidf(0_i64, i))
      r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
      r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
      w_l = w_l + (order - 1_i64) * psi_left * (w_x_left * r_l(0_i64) + &
            w_y_left * r_l(1_i64))
      w_r = w_r
      sol_0004 = 0.0_f64
      sign_0004 = u_face(i) * normal(0_i64) + v_face(i) * normal(1_i64) &
            + w_face(i) * normal(2_i64)
      if (sign_0004 >= 0_i64) then
        sol_0004 = w_l
      else
        sol_0004 = w_r
      end if
      flux_w(0_i64) = sign_0004 * sol_0004
      rez_w(cellidf(0_i64, i)) = rez_w(cellidf(0_i64, i)) - flux_w(0_i64 &
            )
    end do

  end subroutine explicitscheme_convective_2d
  !........................................

  !........................................
  subroutine explicitscheme_convective_3d(rez_w, w_c, w_ghost, w_halo, &
        u_face, v_face, w_face, w_x, w_y, w_z, wx_halo, wy_halo, &
        wz_halo, psi, psi_halo, centerc, centerf, centerh, centerg, &
        cellidf, mesuref, normalf, halofid, name, innerfaces, halofaces &
        , boundaryfaces, periodicboundaryfaces, shift, order)

    implicit none

    real(f64), intent(inout) :: rez_w(0:)
    real(f64), intent(in) :: w_c(0:)
    real(f64), intent(in) :: w_ghost(0:)
    real(f64), intent(in) :: w_halo(0:)
    real(f64), intent(in) :: u_face(0:)
    real(f64), intent(in) :: v_face(0:)
    real(f64), intent(in) :: w_face(0:)
    real(f64), intent(in) :: w_x(0:)
    real(f64), intent(in) :: w_y(0:)
    real(f64), intent(in) :: w_z(0:)
    real(f64), intent(in) :: wx_halo(0:)
    real(f64), intent(in) :: wy_halo(0:)
    real(f64), intent(in) :: wz_halo(0:)
    real(f64), intent(in) :: psi(0:)
    real(f64), intent(in) :: psi_halo(0:)
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: centerf(0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    real(f64), intent(in) :: centerg(0:,0:)
    integer(i64), intent(inout) :: cellidf(0:,0:)
    real(f64), intent(in) :: mesuref(0:)
    real(f64), intent(in) :: normalf(0:,0:)
    integer(i64), intent(in) :: halofid(0:)
    integer(i64), intent(in) :: name(0:)
    integer(i64), intent(in) :: innerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: boundaryfaces(0:)
    integer(i64), intent(in) :: periodicboundaryfaces(0:)
    real(f64), intent(in) :: shift(0:,0:)
    integer(i64), value :: order
    real(f64) :: center_left(0:2_i64)
    real(f64) :: center_right(0:2_i64)
    real(f64) :: r_l(0:2_i64)
    real(f64) :: r_r(0:2_i64)
    real(f64) :: normal(0:2_i64)
    real(f64) :: flux_w(0:0_i64)
    integer(i64) :: i
    real(f64) :: w_l
    real(f64) :: w_r
    real(f64) :: w_x_left
    real(f64) :: w_x_right
    real(f64) :: w_y_left
    real(f64) :: w_y_right
    real(f64) :: w_z_left
    real(f64) :: w_z_right
    real(f64) :: psi_left
    real(f64) :: psi_right
    integer(i64) :: Dummy_0013
    real(f64) :: sol_0005
    real(f64) :: sign_0005
    integer(i64) :: Dummy_0014
    real(f64) :: sol_0006
    real(f64) :: sign_0006
    integer(i64) :: Dummy_0015
    real(f64) :: sol_0007
    real(f64) :: sign_0007
    integer(i64) :: Dummy_0016
    real(f64) :: sol_0008
    real(f64) :: sign_0008

    center_left = 0.0_f64
    center_right = 0.0_f64
    r_l = 0.0_f64
    r_r = 0.0_f64
    normal = 0.0_f64
    flux_w = 0.0_f64
    rez_w(:) = 0.0_f64
    do Dummy_0013 = 0_i64, size(innerfaces, kind=i64) - 1_i64, 1_i64
      i = innerfaces(Dummy_0013)
      w_l = w_c(cellidf(0_i64, i))
      normal(:) = normalf(:, i)
      w_r = w_c(cellidf(1_i64, i))
      center_left(:) = centerc(:, cellidf(0_i64, i))
      center_right(:) = centerc(:, cellidf(1_i64, i))
      w_x_left = w_x(cellidf(0_i64, i))
      w_x_right = w_x(cellidf(1_i64, i))
      w_y_left = w_y(cellidf(0_i64, i))
      w_y_right = w_y(cellidf(1_i64, i))
      w_z_left = w_z(cellidf(0_i64, i))
      w_z_right = w_z(cellidf(1_i64, i))
      psi_left = psi(cellidf(0_i64, i))
      psi_right = psi(cellidf(1_i64, i))
      r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
      r_r(0_i64) = centerf(0_i64, i) - center_right(0_i64)
      r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
      r_r(1_i64) = centerf(1_i64, i) - center_right(1_i64)
      r_l(2_i64) = centerf(2_i64, i) - center_left(2_i64)
      r_r(2_i64) = centerf(2_i64, i) - center_right(2_i64)
      w_l = w_l + (order - 1_i64) * psi_left * (w_x_left * r_l(0_i64) + &
            w_y_left * r_l(1_i64) + w_z_left * r_l(2_i64))
      w_r = w_r + (order - 1_i64) * psi_right * (w_x_right * r_r(0_i64) &
            + w_y_right * r_r(1_i64) + w_z_right * r_r(2_i64))
      sol_0005 = 0.0_f64
      sign_0005 = u_face(i) * normal(0_i64) + v_face(i) * normal(1_i64) &
            + w_face(i) * normal(2_i64)
      if (sign_0005 >= 0_i64) then
        sol_0005 = w_l
      else
        sol_0005 = w_r
      end if
      flux_w(0_i64) = sign_0005 * sol_0005
      rez_w(cellidf(0_i64, i)) = rez_w(cellidf(0_i64, i)) - flux_w(0_i64 &
            )
      rez_w(cellidf(1_i64, i)) = rez_w(cellidf(1_i64, i)) + flux_w(0_i64 &
            )
    end do
    do Dummy_0014 = 0_i64, size(periodicboundaryfaces, kind=i64) - 1_i64 &
          , 1_i64
      i = periodicboundaryfaces(Dummy_0014)
      w_l = w_c(cellidf(0_i64, i))
      normal(:) = normalf(:, i)
      w_r = w_c(cellidf(1_i64, i))
      center_left(:) = centerc(:, cellidf(0_i64, i))
      center_right(:) = centerc(:, cellidf(1_i64, i))
      w_x_left = w_x(cellidf(0_i64, i))
      w_x_right = w_x(cellidf(1_i64, i))
      w_y_left = w_y(cellidf(0_i64, i))
      w_y_right = w_y(cellidf(1_i64, i))
      w_z_left = w_z(cellidf(0_i64, i))
      w_z_right = w_z(cellidf(1_i64, i))
      psi_left = psi(cellidf(0_i64, i))
      psi_right = psi(cellidf(1_i64, i))
      if (name(i) == 11_i64 .or. name(i) == 22_i64) then
        r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
        r_r(0_i64) = centerf(0_i64, i) - center_right(0_i64) - shift( &
              0_i64, cellidf(1_i64, i))
        r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
        r_r(1_i64) = centerf(1_i64, i) - center_right(1_i64)
        r_l(2_i64) = centerf(2_i64, i) - center_left(2_i64)
        r_r(2_i64) = centerf(2_i64, i) - center_right(2_i64)
      end if
      if (name(i) == 33_i64 .or. name(i) == 44_i64) then
        r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
        r_r(0_i64) = centerf(0_i64, i) - center_right(0_i64)
        r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
        r_r(1_i64) = centerf(1_i64, i) - center_right(1_i64) - shift( &
              1_i64, cellidf(1_i64, i))
        r_l(2_i64) = centerf(2_i64, i) - center_left(2_i64)
        r_r(2_i64) = centerf(2_i64, i) - center_right(2_i64)
      end if
      if (name(i) == 55_i64 .or. name(i) == 66_i64) then
        r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
        r_r(0_i64) = centerf(0_i64, i) - center_right(0_i64)
        r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
        r_r(1_i64) = centerf(1_i64, i) - center_right(1_i64)
        r_l(2_i64) = centerf(2_i64, i) - center_left(2_i64)
        r_r(2_i64) = centerf(2_i64, i) - center_right(2_i64) - shift( &
              2_i64, cellidf(1_i64, i))
      end if
      w_l = w_l + (order - 1_i64) * psi_left * (w_x_left * r_l(0_i64) + &
            w_y_left * r_l(1_i64) + w_z_left * r_l(2_i64))
      w_r = w_r + (order - 1_i64) * psi_right * (w_x_right * r_r(0_i64) &
            + w_y_right * r_r(1_i64) + w_z_right * r_r(2_i64))
      sol_0006 = 0.0_f64
      sign_0006 = u_face(i) * normal(0_i64) + v_face(i) * normal(1_i64) &
            + w_face(i) * normal(2_i64)
      if (sign_0006 >= 0_i64) then
        sol_0006 = w_l
      else
        sol_0006 = w_r
      end if
      flux_w(0_i64) = sign_0006 * sol_0006
      rez_w(cellidf(0_i64, i)) = rez_w(cellidf(0_i64, i)) - flux_w(0_i64 &
            )
    end do
    do Dummy_0015 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0015)
      w_l = w_c(cellidf(0_i64, i))
      normal(:) = normalf(:, i)
      w_r = w_halo(halofid(i))
      center_left(:) = centerc(:, cellidf(0_i64, i))
      center_right(:) = centerh(:, halofid(i))
      w_x_left = w_x(cellidf(0_i64, i))
      w_x_right = wx_halo(halofid(i))
      w_y_left = w_y(cellidf(0_i64, i))
      w_y_right = wy_halo(halofid(i))
      w_z_left = w_z(cellidf(0_i64, i))
      w_z_right = wz_halo(halofid(i))
      psi_left = psi(cellidf(0_i64, i))
      psi_right = psi_halo(halofid(i))
      r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
      r_r(0_i64) = centerf(0_i64, i) - center_right(0_i64)
      r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
      r_r(1_i64) = centerf(1_i64, i) - center_right(1_i64)
      r_l(2_i64) = centerf(2_i64, i) - center_left(2_i64)
      r_r(2_i64) = centerf(2_i64, i) - center_right(2_i64)
      w_l = w_l + (order - 1_i64) * psi_left * (w_x_left * r_l(0_i64) + &
            w_y_left * r_l(1_i64) + w_z_left * r_l(2_i64))
      w_r = w_r + (order - 1_i64) * psi_right * (w_x_right * r_r(0_i64) &
            + w_y_right * r_r(1_i64) + w_z_right * r_r(2_i64))
      sol_0007 = 0.0_f64
      sign_0007 = u_face(i) * normal(0_i64) + v_face(i) * normal(1_i64) &
            + w_face(i) * normal(2_i64)
      if (sign_0007 >= 0_i64) then
        sol_0007 = w_l
      else
        sol_0007 = w_r
      end if
      flux_w(0_i64) = sign_0007 * sol_0007
      rez_w(cellidf(0_i64, i)) = rez_w(cellidf(0_i64, i)) - flux_w(0_i64 &
            )
    end do
    do Dummy_0016 = 0_i64, size(boundaryfaces, kind=i64) - 1_i64, 1_i64
      i = boundaryfaces(Dummy_0016)
      w_l = w_c(cellidf(0_i64, i))
      normal(:) = normalf(:, i)
      w_r = w_ghost(i)
      center_left(:) = centerc(:, cellidf(0_i64, i))
      w_x_left = w_x(cellidf(0_i64, i))
      w_y_left = w_y(cellidf(0_i64, i))
      w_z_left = w_z(cellidf(0_i64, i))
      psi_left = psi(cellidf(0_i64, i))
      r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
      r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
      r_l(2_i64) = centerf(2_i64, i) - center_left(2_i64)
      w_l = w_l + (order - 1_i64) * psi_left * (w_x_left * r_l(0_i64) + &
            w_y_left * r_l(1_i64) + w_z_left * r_l(2_i64))
      w_r = w_r
      sol_0008 = 0.0_f64
      sign_0008 = u_face(i) * normal(0_i64) + v_face(i) * normal(1_i64) &
            + w_face(i) * normal(2_i64)
      if (sign_0008 >= 0_i64) then
        sol_0008 = w_l
      else
        sol_0008 = w_r
      end if
      flux_w(0_i64) = sign_0008 * sol_0008
      rez_w(cellidf(0_i64, i)) = rez_w(cellidf(0_i64, i)) - flux_w(0_i64 &
            )
    end do

  end subroutine explicitscheme_convective_3d
  !........................................

end module pyccel_fvm
