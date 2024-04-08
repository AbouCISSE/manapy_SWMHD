module tools


  use, intrinsic :: ISO_C_Binding, only : i64 => C_INT64_T , f64 => &
        C_DOUBLE
  implicit none

  contains

  !........................................
  subroutine update_SW(h_c, hu_c, hv_c, hc_c, Z_c, rez_h, rez_hu, rez_hv &
        , rez_hc, rez_Z, src_h, src_hu, src_hv, src_hc, src_Z, &
        dissip_hc, corio_hu, corio_hv, wind_hu, wind_hv, dtime, vol)

    implicit none

    real(f64), intent(inout) :: h_c(0:)
    real(f64), intent(inout) :: hu_c(0:)
    real(f64), intent(inout) :: hv_c(0:)
    real(f64), intent(inout) :: hc_c(0:)
    real(f64), intent(inout) :: Z_c(0:)
    real(f64), intent(in) :: rez_h(0:)
    real(f64), intent(in) :: rez_hu(0:)
    real(f64), intent(in) :: rez_hv(0:)
    real(f64), intent(in) :: rez_hc(0:)
    real(f64), intent(in) :: rez_Z(0:)
    real(f64), intent(in) :: src_h(0:)
    real(f64), intent(in) :: src_hu(0:)
    real(f64), intent(in) :: src_hv(0:)
    real(f64), intent(in) :: src_hc(0:)
    real(f64), intent(in) :: src_Z(0:)
    real(f64), intent(in) :: dissip_hc(0:)
    real(f64), intent(in) :: corio_hu(0:)
    real(f64), intent(in) :: corio_hv(0:)
    real(f64), value :: wind_hu
    real(f64), value :: wind_hv
    real(f64), value :: dtime
    real(f64), intent(in) :: vol(0:)
    integer(i64) :: i

    !$omp parallel do
    do i = 0_i64, size(h_c, kind=i64) - 1_i64, 1_i64
      h_c(i) = h_c(i) + dtime * (rez_h(i) + src_h(i)) / vol(i)
      hu_c(i) = hu_c(i) + dtime * ((rez_hu(i) + src_hu(i)) / vol(i) + &
            corio_hu(i) + wind_hu)
      hv_c(i) = hv_c(i) + dtime * ((rez_hv(i) + src_hv(i)) / vol(i) + &
            corio_hv(i) + wind_hv)
      hc_c(i) = hc_c(i) + dtime * (rez_hc(i) + src_hc(i) - dissip_hc(i &
            )) / vol(i)
      Z_c(i) = Z_c(i) + dtime * (rez_Z(i) + src_Z(i)) / vol(i)
    end do
    !$omp end parallel do

  end subroutine update_SW
  !........................................

  !........................................
  function time_step_SW(h_c, hu_c, hv_c, cfl, normal, mesure, volume, &
        faceid, Dxx, Dyy) result(dt)

    implicit none

    real(f64) :: dt
    real(f64), intent(in) :: h_c(0:)
    real(f64), intent(in) :: hu_c(0:)
    real(f64), intent(in) :: hv_c(0:)
    real(f64), value :: cfl
    real(f64), intent(in) :: normal(0:,0:)
    real(f64), intent(in) :: mesure(0:)
    real(f64), intent(in) :: volume(0:)
    integer(i64), intent(in) :: faceid(0:,0:)
    real(f64), value :: Dxx
    real(f64), value :: Dyy
    real(f64) :: grav
    integer(i64) :: nbelement
    real(f64) :: u_n
    real(f64) :: norm(0:2_i64)
    integer(i64) :: i
    real(f64) :: velson
    real(f64) :: lam
    integer(i64) :: j
    real(f64) :: lam_convect
    real(f64) :: mes
    real(f64) :: lam_diff

    grav = 9.81_f64
    nbelement = size(faceid,2,i64)
    u_n = 0.0_f64
    norm = 0.0_f64
    dt = 1000000.0_f64
    do i = 0_i64, nbelement - 1_i64, 1_i64
      velson = sqrt(grav * h_c(i))
      lam = 0.0_f64
      do j = 0_i64, 2_i64, 1_i64
        norm(:) = normal(:, faceid(j, i))
        !convective part
        u_n = abs(hu_c(i) / h_c(i) * norm(0_i64) + hv_c(i) / h_c(i) * &
              norm(1_i64))
        lam_convect = u_n / mesure(faceid(j, i)) + velson
        lam = lam + lam_convect * mesure(faceid(j, i))
        !diffusion part
        mes = sqrt(norm(0_i64) * norm(0_i64) + norm(1_i64) * norm(1_i64 &
              ))
        lam_diff = Dxx * mes ** 2_i64 + Dyy * mes ** 2_i64
        lam = lam + lam_diff / volume(i)
      end do
      dt = minval([dt, cfl * volume(i) / lam])
    end do
    return

  end function time_step_SW
  !........................................

  !........................................
  function time_step_SWMHD(h_c, hu_c, hv_c, hB1_c, hB2_c, cfl, normal, &
        mesure, volume, faceid) result(dt)

    implicit none

    real(f64) :: dt
    real(f64), intent(in) :: h_c(0:)
    real(f64), intent(in) :: hu_c(0:)
    real(f64), intent(in) :: hv_c(0:)
    real(f64), intent(in) :: hB1_c(0:)
    real(f64), intent(in) :: hB2_c(0:)
    real(f64), value :: cfl
    real(f64), intent(in) :: normal(0:,0:)
    real(f64), intent(in) :: mesure(0:)
    real(f64), intent(in) :: volume(0:)
    integer(i64), intent(in) :: faceid(0:,0:)
    real(f64) :: grav
    integer(i64) :: nbelement
    real(f64) :: u_n
    real(f64) :: B_n
    integer(i64) :: i
    real(f64) :: lam
    integer(i64) :: j
    real(f64) :: wb
    real(f64) :: lam1
    real(f64) :: lam2
    real(f64) :: lam3
    real(f64) :: lam4
    real(f64) :: lam_convect

    grav = 9.81_f64
    nbelement = size(faceid,2,i64)
    u_n = 0.0_f64
    B_n = 0.0_f64
    dt = 1000000.0_f64
    do i = 0_i64, nbelement - 1_i64, 1_i64
      lam = 0.0_f64
      do j = 0_i64, 2_i64, 1_i64
        u_n = abs(hu_c(i) / h_c(i) * normal(0_i64, faceid(j, i)) + hv_c( &
              i) / h_c(i) * normal(1_i64, faceid(j, i))) / mesure( &
              faceid(j, i))
        B_n = abs(hB1_c(i) / h_c(i) * normal(0_i64, faceid(j, i)) + &
              hB2_c(i) / h_c(i) * normal(1_i64, faceid(j, i))) / mesure &
              (faceid(j, i))
        wb = sqrt(B_n ** 2_i64 + grav * h_c(i))
        lam1 = abs(u_n - wb)
        lam2 = abs(u_n - B_n)
        lam3 = abs(u_n + B_n)
        lam4 = abs(u_n + wb)
        !lam_convect = u_n + velson #max(lam1, lam2, lam3, lam4)
        lam_convect = maxval([lam1, lam2, lam3, lam4])
        lam = lam + lam_convect * mesure(faceid(j, i))
        !lam += lam_convect * mesure[faceid[i][j]]
      end do
      !dt_c[i]  = cfl * volume[i]/lam
      dt = minval([dt, cfl * volume(i) / lam])
    end do
    return

  end function time_step_SWMHD
  !........................................

  !........................................
  subroutine update_SWMHD(h_c, hu_c, hv_c, hB1_c, hB2_c, PSI_c, Z_c, &
        rez_h, rez_hu, rez_hv, rez_hB1, rez_hB2, rez_PSI, rez_Z, src_h, &
        src_hu, src_hv, src_hB1, src_hB2, src_PSI, src_Z, corio_hu, &
        corio_hv, dtime, vol, GLM, cpsi)

    implicit none

    real(f64), intent(inout) :: h_c(0:)
    real(f64), intent(inout) :: hu_c(0:)
    real(f64), intent(inout) :: hv_c(0:)
    real(f64), intent(inout) :: hB1_c(0:)
    real(f64), intent(inout) :: hB2_c(0:)
    real(f64), intent(inout) :: PSI_c(0:)
    real(f64), intent(inout) :: Z_c(0:)
    real(f64), intent(in) :: rez_h(0:)
    real(f64), intent(in) :: rez_hu(0:)
    real(f64), intent(in) :: rez_hv(0:)
    real(f64), intent(in) :: rez_hB1(0:)
    real(f64), intent(in) :: rez_hB2(0:)
    real(f64), intent(in) :: rez_PSI(0:)
    real(f64), intent(in) :: rez_Z(0:)
    real(f64), intent(in) :: src_h(0:)
    real(f64), intent(in) :: src_hu(0:)
    real(f64), intent(in) :: src_hv(0:)
    real(f64), intent(in) :: src_hB1(0:)
    real(f64), intent(in) :: src_hB2(0:)
    real(f64), intent(in) :: src_PSI(0:)
    real(f64), intent(in) :: src_Z(0:)
    real(f64), intent(in) :: corio_hu(0:)
    real(f64), intent(in) :: corio_hv(0:)
    real(f64), value :: dtime
    real(f64), intent(in) :: vol(0:)
    integer(i64), value :: GLM
    real(f64), value :: cpsi
    real(f64) :: cr

    h_c(:) = h_c(:) + dtime * (rez_h(:) + src_h(:)) / vol(:)
    hu_c(:) = hu_c(:) + dtime * ((rez_hu(:) + src_hu(:)) / vol(:) + &
          corio_hu(:))
    hv_c(:) = hv_c(:) + dtime * ((rez_hv(:) + src_hv(:)) / vol(:) + &
          corio_hv(:))
    hB1_c(:) = hB1_c(:) + dtime * (rez_hB1(:) + src_hB1(:)) / vol(:)
    hB2_c(:) = hB2_c(:) + dtime * (rez_hB2(:) + src_hB2(:)) / vol(:)
    PSI_c(:) = PSI_c(:) + dtime * (rez_PSI(:) + src_PSI(:)) / vol(:)
    Z_c(:) = Z_c(:) + dtime * (rez_Z(:) + src_Z(:)) / vol(:)
    if (GLM == 10_i64) then
      cr = 0.01_f64
      PSI_c(:) = exp((-dtime) * (cpsi / cr)) * PSI_c(:)
    end if
    !renvoi max

  end subroutine update_SWMHD
  !........................................

  !........................................
  function cpsi_global(h_c, hu_c, hv_c, hB1_c, hB2_c, cfl, normal, &
        mesure, volume, faceid) result(cpsiglobal)

    implicit none

    real(f64) :: cpsiglobal
    real(f64), intent(in) :: h_c(0:)
    real(f64), intent(in) :: hu_c(0:)
    real(f64), intent(in) :: hv_c(0:)
    real(f64), intent(in) :: hB1_c(0:)
    real(f64), intent(in) :: hB2_c(0:)
    real(f64), value :: cfl
    real(f64), intent(in) :: normal(0:,0:)
    real(f64), intent(in) :: mesure(0:)
    real(f64), intent(in) :: volume(0:)
    integer(i64), intent(in) :: faceid(0:,0:)
    integer(i64) :: grav
    integer(i64) :: nbelement
    real(f64) :: u_n
    real(f64) :: B_n
    integer(i64) :: i
    real(f64) :: lam
    integer(i64) :: j
    real(f64) :: wb
    real(f64) :: lam1
    real(f64) :: lam2
    real(f64) :: lam3
    real(f64) :: lam4
    real(f64) :: lam_convect

    grav = 1_i64
    nbelement = size(faceid,2,i64)
    u_n = 0.0_f64
    B_n = 0.0_f64
    cpsiglobal = 0.0_f64
    do i = 0_i64, nbelement - 1_i64, 1_i64
      lam = 0.0_f64
      do j = 0_i64, 2_i64, 1_i64
        u_n = abs(hu_c(i) / h_c(i) * normal(0_i64, faceid(j, i)) + hv_c( &
              i) / h_c(i) * normal(1_i64, faceid(j, i))) / mesure( &
              faceid(j, i))
        B_n = abs(hB1_c(i) / h_c(i) * normal(0_i64, faceid(j, i)) + &
              hB2_c(i) / h_c(i) * normal(1_i64, faceid(j, i))) / mesure &
              (faceid(j, i))
        wb = sqrt(B_n ** 2_i64 + grav * h_c(i))
        lam1 = abs(u_n - wb)
        lam2 = abs(u_n - B_n)
        lam3 = abs(u_n + B_n)
        lam4 = abs(u_n + wb)
        lam_convect = maxval([lam1, lam2, lam3, lam4])
        lam = lam + lam_convect * mesure(faceid(j, i))
      end do
      cpsiglobal = maxval([cpsiglobal, lam_convect])
    end do
    return

  end function cpsi_global
  !........................................

  !........................................
  subroutine initialisation_SW(h, hu, hv, hc, Z, center) 

    implicit none

    real(f64), intent(inout) :: h(0:)
    real(f64), intent(inout) :: hu(0:)
    real(f64), intent(inout) :: hv(0:)
    real(f64), intent(inout) :: hc(0:)
    real(f64), intent(inout) :: Z(0:)
    real(f64), intent(in) :: center(0:,0:)
    integer(i64) :: nbelements
    integer(i64) :: i
    real(f64) :: xcent

    nbelements = size(center,2,i64)
    do i = 0_i64, nbelements - 1_i64, 1_i64
      xcent = center(0_i64, i)
      h(i) = 2_i64
      Z(i) = 0.0_f64
      if (xcent < 0.5_f64) then
        h(i) = 5.0_f64
      end if
      hu(i) = 0.0_f64
      hv(i) = 0.0_f64
      hc(i) = 0.0_f64
    end do

  end subroutine initialisation_SW
  !........................................

  !........................................
  subroutine term_source_srnh_SWf(src_h, src_hu, src_hv, src_hc, src_Z, &
        h_c, hu_c, hv_c, hc_c, Z_c, h_ghost, hu_ghost, hv_ghost, &
        hc_ghost, Z_ghost, h_halo, hu_halo, hv_halo, hc_halo, Z_halo, &
        h_x, h_y, psi, hx_halo, hy_halo, psi_halo, nodeidc, faceidc, &
        cellidc, cellidf, centerc, normalc, namef, centerf, centerh, &
        vertexn, halofid, order)

    implicit none

    real(f64), intent(inout) :: src_h(0:)
    real(f64), intent(inout) :: src_hu(0:)
    real(f64), intent(inout) :: src_hv(0:)
    real(f64), intent(inout) :: src_hc(0:)
    real(f64), intent(inout) :: src_Z(0:)
    real(f64), intent(in) :: h_c(0:)
    real(f64), intent(in) :: hu_c(0:)
    real(f64), intent(in) :: hv_c(0:)
    real(f64), intent(in) :: hc_c(0:)
    real(f64), intent(in) :: Z_c(0:)
    real(f64), intent(in) :: h_ghost(0:)
    real(f64), intent(in) :: hu_ghost(0:)
    real(f64), intent(in) :: hv_ghost(0:)
    real(f64), intent(in) :: hc_ghost(0:)
    real(f64), intent(in) :: Z_ghost(0:)
    real(f64), intent(in) :: h_halo(0:)
    real(f64), intent(in) :: hu_halo(0:)
    real(f64), intent(in) :: hv_halo(0:)
    real(f64), intent(in) :: hc_halo(0:)
    real(f64), intent(in) :: Z_halo(0:)
    real(f64), intent(in) :: h_x(0:)
    real(f64), intent(in) :: h_y(0:)
    real(f64), intent(in) :: psi(0:)
    real(f64), intent(in) :: hx_halo(0:)
    real(f64), intent(in) :: hy_halo(0:)
    real(f64), intent(in) :: psi_halo(0:)
    integer(i64), intent(in) :: nodeidc(0:,0:)
    integer(i64), intent(in) :: faceidc(0:,0:)
    integer(i64), intent(in) :: cellidc(0:,0:)
    integer(i64), intent(in) :: cellidf(0:,0:)
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: normalc(0:,0:,0:)
    integer(i64), intent(in) :: namef(0:)
    real(f64), intent(in) :: centerf(0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    real(f64), intent(in) :: vertexn(0:,0:)
    integer(i64), intent(in) :: halofid(0:)
    integer(i64), value :: order
    real(f64) :: grav
    integer(i64) :: nbelement
    real(f64), allocatable :: hi_p(:)
    real(f64), allocatable :: zi_p(:)
    real(f64), allocatable :: zv(:)
    real(f64), allocatable :: mata(:)
    real(f64), allocatable :: matb(:)
    real(f64), allocatable :: ns(:,:)
    real(f64), allocatable :: ss(:,:)
    real(f64), allocatable :: s_1(:)
    real(f64), allocatable :: s_2(:)
    real(f64), allocatable :: s_3(:)
    real(f64), allocatable :: b(:)
    real(f64), allocatable :: G(:)
    integer(i64) :: i
    real(f64) :: c_1
    real(f64) :: c_2
    real(f64) :: c_3
    real(f64) :: delta
    real(f64) :: deltax
    real(f64) :: deltay
    real(f64) :: deltaz
    real(f64) :: h_1
    real(f64) :: h_2
    real(f64) :: h_3
    real(f64) :: z_1
    real(f64) :: z_2
    real(f64) :: z_3
    integer(i64) :: j
    integer(i64) :: f
    real(f64) :: h_1p
    real(f64) :: z_1p
    real(f64) :: h_p1
    real(f64) :: z_p1

    grav = 9.81_f64
    nbelement = size(h_c, kind=i64)
    allocate(hi_p(0:2_i64))
    hi_p = 0.0_f64
    allocate(zi_p(0:2_i64))
    zi_p = 0.0_f64
    allocate(zv(0:2_i64))
    zv = 0.0_f64
    allocate(mata(0:2_i64))
    mata = 0.0_f64
    allocate(matb(0:2_i64))
    matb = 0.0_f64
    allocate(ns(0:2_i64, 0:2_i64))
    ns = 0.0_f64
    allocate(ss(0:2_i64, 0:2_i64))
    ss = 0.0_f64
    allocate(s_1(0:2_i64))
    s_1 = 0.0_f64
    allocate(s_2(0:2_i64))
    s_2 = 0.0_f64
    allocate(s_3(0:2_i64))
    s_3 = 0.0_f64
    allocate(b(0:2_i64))
    b = 0.0_f64
    allocate(G(0:2_i64))
    G = 0.0_f64
    do i = 0_i64, nbelement - 1_i64, 1_i64
      G(:) = centerc(:, i)
      c_1 = 0.0_f64
      c_2 = 0.0_f64
      do j = 0_i64, 2_i64, 1_i64
        f = faceidc(j, i)
        ss(:, j) = normalc(:, j, i)
        if (namef(f) == 10_i64) then
          h_1p = h_c(i)
          z_1p = Z_c(i)
          h_p1 = h_halo(halofid(f))
          z_p1 = Z_halo(halofid(f))
        else if (namef(f) == 0_i64) then
          h_1p = h_c(i)
          z_1p = Z_c(i)
          h_p1 = h_c(cellidc(j, i))
          z_p1 = Z_c(cellidc(j, i))
        else
          h_1p = h_c(i)
          z_1p = Z_c(i)
          h_p1 = h_ghost(f)
          z_p1 = Z_ghost(f)
        end if
        zv(j) = z_p1
        mata(j) = h_p1 * ss(0_i64, j)
        matb(j) = h_p1 * ss(1_i64, j)
        c_1 = c_1 + 0.5_f64 * (h_1p + h_p1) * 0.5_f64 * (h_1p + h_p1) * &
              ss(0_i64, j)
        c_2 = c_2 + 0.5_f64 * (h_1p + h_p1) * 0.5_f64 * (h_1p + h_p1) * &
              ss(1_i64, j)
        hi_p(j) = h_1p
        zi_p(j) = z_1p
      end do
      c_3 = 3.0_f64 * h_1p
      delta = mata(1_i64) * matb(2_i64) - mata(2_i64) * matb(1_i64) - ( &
            mata(0_i64) * matb(2_i64) - matb(0_i64) * mata(2_i64)) + ( &
            mata(0_i64) * matb(1_i64) - matb(0_i64) * mata(1_i64))
      deltax = c_3 * (mata(1_i64) * matb(2_i64) - mata(2_i64) * matb( &
            1_i64)) - (c_1 * matb(2_i64) - c_2 * mata(2_i64)) + (c_1 * &
            matb(1_i64) - c_2 * mata(1_i64))
      deltay = c_1 * matb(2_i64) - c_2 * mata(2_i64) - c_3 * (mata(0_i64 &
            ) * matb(2_i64) - matb(0_i64) * mata(2_i64)) + (mata(0_i64 &
            ) * c_2 - matb(0_i64) * c_1)
      deltaz = mata(1_i64) * c_2 - matb(1_i64) * c_1 - (mata(0_i64) * &
            c_2 - matb(0_i64) * c_1) + c_3 * (mata(0_i64) * matb(1_i64 &
            ) - matb(0_i64) * mata(1_i64))
      h_1 = deltax / delta
      h_2 = deltay / delta
      h_3 = deltaz / delta
      z_1 = zi_p(0_i64) + hi_p(0_i64) - h_1
      z_2 = zi_p(1_i64) + hi_p(1_i64) - h_2
      z_3 = zi_p(2_i64) + hi_p(2_i64) - h_3
      b(:) = vertexn(0_i64:2_i64, nodeidc(1_i64, i))
      ns(:, 0_i64) = [G(1_i64) - b(1_i64), -(G(0_i64) - b(0_i64)), &
            0.0_f64]
      ns(:, 1_i64) = ns(:, 0_i64) - ss(:, 1_i64)
      ns(:, 2_i64) = ns(:, 0_i64) + ss(:, 0_i64)
      s_1 = 0.5_f64 * h_1 * (zv(0_i64) * ss(:, 0_i64) + z_2 * ns(:, &
            0_i64) + z_3 * (-1_i64) * ns(:, 2_i64))
      s_2 = 0.5_f64 * h_2 * (zv(1_i64) * ss(:, 1_i64) + z_1 * (-1_i64) * &
            ns(:, 0_i64) + z_3 * ns(:, 1_i64))
      s_3 = 0.5_f64 * h_3 * (zv(2_i64) * ss(:, 2_i64) + z_1 * ns(:, &
            2_i64) + z_2 * (-1_i64) * ns(:, 1_i64))
      !TODO
      src_h(i) = 0_i64
      src_hu(i) = (-grav) * (s_1(0_i64) + s_2(0_i64) + s_3(0_i64))
      src_hv(i) = (-grav) * (s_1(1_i64) + s_2(1_i64) + s_3(1_i64))
      src_hc(i) = 0.0_f64
      src_Z(i) = 0.0_f64
    end do
    if (allocated(hi_p)) then
      deallocate(hi_p)
    end if
    if (allocated(zi_p)) then
      deallocate(zi_p)
    end if
    if (allocated(zv)) then
      deallocate(zv)
    end if
    if (allocated(mata)) then
      deallocate(mata)
    end if
    if (allocated(matb)) then
      deallocate(matb)
    end if
    if (allocated(ns)) then
      deallocate(ns)
    end if
    if (allocated(ss)) then
      deallocate(ss)
    end if
    if (allocated(s_1)) then
      deallocate(s_1)
    end if
    if (allocated(s_2)) then
      deallocate(s_2)
    end if
    if (allocated(s_3)) then
      deallocate(s_3)
    end if
    if (allocated(b)) then
      deallocate(b)
    end if
    if (allocated(G)) then
      deallocate(G)
    end if

  end subroutine term_source_srnh_SWf
  !........................................

  !........................................
  subroutine term_source_srnh_SWMHD(src_h, src_hu, src_hv, src_hB1, &
        src_hB2, src_PSI, src_Z, h_c, hu_c, hv_c, hc_c, Z_c, h_ghost, &
        hu_ghost, hv_ghost, hc_ghost, Z_ghost, h_halo, hu_halo, hv_halo &
        , hc_halo, Z_halo, h_x, h_y, psi, hx_halo, hy_halo, psi_halo, &
        nodeidc, faceidc, cellidc, cellidf, centerc, normalc, namef, &
        centerf, centerh, vertexn, halofid, order)

    implicit none

    real(f64), intent(inout) :: src_h(0:)
    real(f64), intent(inout) :: src_hu(0:)
    real(f64), intent(inout) :: src_hv(0:)
    real(f64), intent(inout) :: src_hB1(0:)
    real(f64), intent(inout) :: src_hB2(0:)
    real(f64), intent(inout) :: src_PSI(0:)
    real(f64), intent(inout) :: src_Z(0:)
    real(f64), intent(in) :: h_c(0:)
    real(f64), intent(in) :: hu_c(0:)
    real(f64), intent(in) :: hv_c(0:)
    real(f64), intent(in) :: hc_c(0:)
    real(f64), intent(in) :: Z_c(0:)
    real(f64), intent(in) :: h_ghost(0:)
    real(f64), intent(in) :: hu_ghost(0:)
    real(f64), intent(in) :: hv_ghost(0:)
    real(f64), intent(in) :: hc_ghost(0:)
    real(f64), intent(in) :: Z_ghost(0:)
    real(f64), intent(in) :: h_halo(0:)
    real(f64), intent(in) :: hu_halo(0:)
    real(f64), intent(in) :: hv_halo(0:)
    real(f64), intent(in) :: hc_halo(0:)
    real(f64), intent(in) :: Z_halo(0:)
    real(f64), intent(in) :: h_x(0:)
    real(f64), intent(in) :: h_y(0:)
    real(f64), intent(in) :: psi(0:)
    real(f64), intent(in) :: hx_halo(0:)
    real(f64), intent(in) :: hy_halo(0:)
    real(f64), intent(in) :: psi_halo(0:)
    integer(i64), intent(in) :: nodeidc(0:,0:)
    integer(i64), intent(in) :: faceidc(0:,0:)
    integer(i64), intent(in) :: cellidc(0:,0:)
    integer(i64), intent(in) :: cellidf(0:,0:)
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: normalc(0:,0:,0:)
    integer(i64), intent(in) :: namef(0:)
    real(f64), intent(in) :: centerf(0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    real(f64), intent(in) :: vertexn(0:,0:)
    integer(i64), intent(in) :: halofid(0:)
    integer(i64), value :: order
    real(f64) :: grav
    integer(i64) :: nbelement
    real(f64), allocatable :: hi_p(:)
    real(f64), allocatable :: zi_p(:)
    real(f64), allocatable :: zv(:)
    real(f64), allocatable :: mata(:)
    real(f64), allocatable :: matb(:)
    real(f64), allocatable :: ns(:,:)
    real(f64), allocatable :: ss(:,:)
    real(f64), allocatable :: s_1(:)
    real(f64), allocatable :: s_2(:)
    real(f64), allocatable :: s_3(:)
    real(f64), allocatable :: b(:)
    real(f64), allocatable :: G(:)
    integer(i64) :: i
    real(f64) :: c_1
    real(f64) :: c_2
    real(f64) :: c_3
    real(f64) :: delta
    real(f64) :: deltax
    real(f64) :: deltay
    real(f64) :: deltaz
    real(f64) :: h_1
    real(f64) :: h_2
    real(f64) :: h_3
    real(f64) :: z_1
    real(f64) :: z_2
    real(f64) :: z_3
    integer(i64) :: j
    integer(i64) :: f
    real(f64) :: h_1p
    real(f64) :: z_1p
    real(f64) :: h_p1
    real(f64) :: z_p1

    grav = 9.81_f64
    nbelement = size(h_c, kind=i64)
    allocate(hi_p(0:2_i64))
    hi_p = 0.0_f64
    allocate(zi_p(0:2_i64))
    zi_p = 0.0_f64
    allocate(zv(0:2_i64))
    zv = 0.0_f64
    allocate(mata(0:2_i64))
    mata = 0.0_f64
    allocate(matb(0:2_i64))
    matb = 0.0_f64
    allocate(ns(0:2_i64, 0:2_i64))
    ns = 0.0_f64
    allocate(ss(0:2_i64, 0:2_i64))
    ss = 0.0_f64
    allocate(s_1(0:2_i64))
    s_1 = 0.0_f64
    allocate(s_2(0:2_i64))
    s_2 = 0.0_f64
    allocate(s_3(0:2_i64))
    s_3 = 0.0_f64
    allocate(b(0:2_i64))
    b = 0.0_f64
    allocate(G(0:2_i64))
    G = 0.0_f64
    do i = 0_i64, nbelement - 1_i64, 1_i64
      G(:) = centerc(:, i)
      c_1 = 0.0_f64
      c_2 = 0.0_f64
      do j = 0_i64, 2_i64, 1_i64
        f = faceidc(j, i)
        ss(:, j) = normalc(:, j, i)
        if (namef(f) == 10_i64) then
          h_1p = h_c(i)
          z_1p = Z_c(i)
          h_p1 = h_halo(halofid(f))
          z_p1 = Z_halo(halofid(f))
        else if (namef(f) == 0_i64) then
          h_1p = h_c(i)
          z_1p = Z_c(i)
          h_p1 = h_c(cellidc(j, i))
          z_p1 = Z_c(cellidc(j, i))
        else
          h_1p = h_c(i)
          z_1p = Z_c(i)
          h_p1 = h_ghost(f)
          z_p1 = Z_ghost(f)
        end if
        zv(j) = z_p1
        mata(j) = h_p1 * ss(0_i64, j)
        matb(j) = h_p1 * ss(1_i64, j)
        c_1 = c_1 + 0.5_f64 * (h_1p + h_p1) * 0.5_f64 * (h_1p + h_p1) * &
              ss(0_i64, j)
        c_2 = c_2 + 0.5_f64 * (h_1p + h_p1) * 0.5_f64 * (h_1p + h_p1) * &
              ss(1_i64, j)
        hi_p(j) = h_1p
        zi_p(j) = z_1p
      end do
      c_3 = 3.0_f64 * h_1p
      delta = mata(1_i64) * matb(2_i64) - mata(2_i64) * matb(1_i64) - ( &
            mata(0_i64) * matb(2_i64) - matb(0_i64) * mata(2_i64)) + ( &
            mata(0_i64) * matb(1_i64) - matb(0_i64) * mata(1_i64))
      deltax = c_3 * (mata(1_i64) * matb(2_i64) - mata(2_i64) * matb( &
            1_i64)) - (c_1 * matb(2_i64) - c_2 * mata(2_i64)) + (c_1 * &
            matb(1_i64) - c_2 * mata(1_i64))
      deltay = c_1 * matb(2_i64) - c_2 * mata(2_i64) - c_3 * (mata(0_i64 &
            ) * matb(2_i64) - matb(0_i64) * mata(2_i64)) + (mata(0_i64 &
            ) * c_2 - matb(0_i64) * c_1)
      deltaz = mata(1_i64) * c_2 - matb(1_i64) * c_1 - (mata(0_i64) * &
            c_2 - matb(0_i64) * c_1) + c_3 * (mata(0_i64) * matb(1_i64 &
            ) - matb(0_i64) * mata(1_i64))
      h_1 = deltax / delta
      h_2 = deltay / delta
      h_3 = deltaz / delta
      z_1 = zi_p(0_i64) + hi_p(0_i64) - h_1
      z_2 = zi_p(1_i64) + hi_p(1_i64) - h_2
      z_3 = zi_p(2_i64) + hi_p(2_i64) - h_3
      b(:) = vertexn(0_i64:2_i64, nodeidc(1_i64, i))
      ns(:, 0_i64) = [G(1_i64) - b(1_i64), -(G(0_i64) - b(0_i64)), &
            0.0_f64]
      ns(:, 1_i64) = ns(:, 0_i64) - ss(:, 1_i64)
      ns(:, 2_i64) = ns(:, 0_i64) + ss(:, 0_i64)
      s_1 = 0.5_f64 * h_1 * (zv(0_i64) * ss(:, 0_i64) + z_2 * ns(:, &
            0_i64) + z_3 * (-1_i64) * ns(:, 2_i64))
      s_2 = 0.5_f64 * h_2 * (zv(1_i64) * ss(:, 1_i64) + z_1 * (-1_i64) * &
            ns(:, 0_i64) + z_3 * ns(:, 1_i64))
      s_3 = 0.5_f64 * h_3 * (zv(2_i64) * ss(:, 2_i64) + z_1 * ns(:, &
            2_i64) + z_2 * (-1_i64) * ns(:, 1_i64))
      !TODO
      src_h(i) = 0.0_f64
      src_hu(i) = (-grav) * (s_1(0_i64) + s_2(0_i64) + s_3(0_i64))
      src_hv(i) = (-grav) * (s_1(1_i64) + s_2(1_i64) + s_3(1_i64))
      src_hB1(i) = 0.0_f64
      src_hB2(i) = 0.0_f64
      src_PSI(i) = 0.0_f64
      src_Z(i) = 0.0_f64
    end do
    if (allocated(hi_p)) then
      deallocate(hi_p)
    end if
    if (allocated(zi_p)) then
      deallocate(zi_p)
    end if
    if (allocated(zv)) then
      deallocate(zv)
    end if
    if (allocated(mata)) then
      deallocate(mata)
    end if
    if (allocated(matb)) then
      deallocate(matb)
    end if
    if (allocated(ns)) then
      deallocate(ns)
    end if
    if (allocated(ss)) then
      deallocate(ss)
    end if
    if (allocated(s_1)) then
      deallocate(s_1)
    end if
    if (allocated(s_2)) then
      deallocate(s_2)
    end if
    if (allocated(s_3)) then
      deallocate(s_3)
    end if
    if (allocated(b)) then
      deallocate(b)
    end if
    if (allocated(G)) then
      deallocate(G)
    end if

  end subroutine term_source_srnh_SWMHD
  !........................................

  !........................................
  !........................................

  !........................................
  !........................................

  !........................................
  subroutine explicitscheme_convective_SW(rez_h, rez_hu, rez_hv, rez_hc, &
        rez_Z, h_c, hu_c, hv_c, hc_c, Z_c, h_ghost, hu_ghost, hv_ghost, &
        hc_ghost, Z_ghost, h_halo, hu_halo, hv_halo, hc_halo, Z_halo, &
        h_x, h_y, hx_halo, hy_halo, hc_x, hc_y, hcx_halo, hcy_halo, psi &
        , psi_halo, centerc, centerf, centerh, centerg, cellidf, &
        mesuref, normalf, halofid, innerfaces, halofaces, boundaryfaces &
        , order)

    implicit none

    real(f64), intent(inout) :: rez_h(0:)
    real(f64), intent(inout) :: rez_hu(0:)
    real(f64), intent(inout) :: rez_hv(0:)
    real(f64), intent(inout) :: rez_hc(0:)
    real(f64), intent(inout) :: rez_Z(0:)
    real(f64), intent(in) :: h_c(0:)
    real(f64), intent(in) :: hu_c(0:)
    real(f64), intent(in) :: hv_c(0:)
    real(f64), intent(in) :: hc_c(0:)
    real(f64), intent(in) :: Z_c(0:)
    real(f64), intent(in) :: h_ghost(0:)
    real(f64), intent(in) :: hu_ghost(0:)
    real(f64), intent(in) :: hv_ghost(0:)
    real(f64), intent(in) :: hc_ghost(0:)
    real(f64), intent(in) :: Z_ghost(0:)
    real(f64), intent(in) :: h_halo(0:)
    real(f64), intent(in) :: hu_halo(0:)
    real(f64), intent(in) :: hv_halo(0:)
    real(f64), intent(in) :: hc_halo(0:)
    real(f64), intent(in) :: Z_halo(0:)
    real(f64), intent(in) :: h_x(0:)
    real(f64), intent(in) :: h_y(0:)
    real(f64), intent(in) :: hx_halo(0:)
    real(f64), intent(in) :: hy_halo(0:)
    real(f64), intent(in) :: hc_x(0:)
    real(f64), intent(in) :: hc_y(0:)
    real(f64), intent(in) :: hcx_halo(0:)
    real(f64), intent(in) :: hcy_halo(0:)
    real(f64), intent(in) :: psi(0:)
    real(f64), intent(in) :: psi_halo(0:)
    real(f64), intent(in), target :: centerc(0:,0:)
    real(f64), intent(in) :: centerf(0:,0:)
    real(f64), intent(in), target :: centerh(0:,0:)
    real(f64), intent(in), target :: centerg(0:,0:)
    integer(i64), intent(inout) :: cellidf(0:,0:)
    real(f64), intent(in) :: mesuref(0:)
    real(f64), intent(in), target :: normalf(0:,0:)
    integer(i64), intent(in) :: halofid(0:)
    integer(i64), intent(in) :: innerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: boundaryfaces(0:)
    integer(i64), value :: order
    real(f64) :: grav
    real(f64), allocatable :: flux(:)
    real(f64), allocatable :: r_l(:)
    real(f64), allocatable :: r_r(:)
    integer(i64) :: i
    real(f64) :: h_l
    real(f64) :: hu_l
    real(f64) :: hv_l
    real(f64) :: hc_l
    real(f64) :: Z_l
    real(f64), pointer :: normal(:)
    real(f64) :: mesure
    real(f64) :: h_r
    real(f64) :: hu_r
    real(f64) :: hv_r
    real(f64) :: hc_r
    real(f64) :: Z_r
    real(f64), pointer :: center_left(:)
    real(f64), pointer :: center_right(:)
    real(f64) :: h_x_left
    real(f64) :: h_x_right
    real(f64) :: h_y_left
    real(f64) :: h_y_right
    real(f64) :: hc_x_left
    real(f64) :: hc_x_right
    real(f64) :: hc_y_left
    real(f64) :: hc_y_right
    real(f64) :: psi_left
    real(f64) :: psi_right
    integer(i64) :: Dummy_0007
    real(f64), allocatable :: ninv_0001(:)
    real(f64), allocatable :: w_dif_0001(:)
    real(f64), allocatable :: rmat_0001(:,:)
    integer(i64) :: As_0001
    real(f64) :: p_0001
    real(f64) :: xi_0001
    real(f64) :: u_h_0001
    real(f64) :: v_h_0001
    real(f64) :: c_h_0001
    real(f64) :: un_h_0001
    real(f64) :: vn_h_0001
    real(f64) :: hroe_0001
    real(f64) :: uroe_0001
    real(f64) :: vroe_0001
    real(f64) :: croe_0001
    real(f64) :: uleft_0001
    real(f64) :: vleft_0001
    real(f64) :: uright_0001
    real(f64) :: vright_0001
    real(f64) :: w_lrh_0001
    real(f64) :: w_lrhu_0001
    real(f64) :: w_lrhv_0001
    real(f64) :: w_lrhc_0001
    real(f64) :: w_lrz_0001
    real(f64) :: d_0001
    real(f64) :: sound_0001
    real(f64) :: Q_0001
    real(f64) :: R_0001
    real(f64) :: theta_0001
    real(f64) :: lambda1_0001
    real(f64) :: lambda2_0001
    real(f64) :: lambda3_0001
    real(f64) :: lambda4_0001
    real(f64) :: lambda5_0001
    real(f64) :: alpha1_0001
    real(f64) :: alpha2_0001
    real(f64) :: alpha3_0001
    real(f64) :: beta_0001
    real(f64) :: gamma1_0001
    real(f64) :: gamma2_0001
    real(f64) :: gamma3_0001
    real(f64) :: sigma1_0001
    real(f64) :: sigma2_0001
    real(f64) :: sigma3_0001
    real(f64) :: epsilon_0001
    real(f64) :: sign1_0001
    real(f64) :: sign2_0001
    real(f64) :: sign3_0001
    real(f64) :: sign4_0001
    real(f64) :: sign5_0001
    real(f64) :: hnew_0001
    real(f64) :: unew_0001
    real(f64) :: vnew_0001
    real(f64) :: cnew_0001
    real(f64) :: znew_0001
    real(f64) :: u_hu_0001
    real(f64) :: u_hv_0001
    real(f64) :: u_hc_0001
    real(f64) :: u_z_0001
    real(f64) :: q_s_0001
    integer(i64) :: Dummy_0008
    real(f64), allocatable :: ninv_0002(:)
    real(f64), allocatable :: w_dif_0002(:)
    real(f64), allocatable :: rmat_0002(:,:)
    integer(i64) :: As_0002
    real(f64) :: p_0002
    real(f64) :: xi_0002
    real(f64) :: u_h_0002
    real(f64) :: v_h_0002
    real(f64) :: c_h_0002
    real(f64) :: un_h_0002
    real(f64) :: vn_h_0002
    real(f64) :: hroe_0002
    real(f64) :: uroe_0002
    real(f64) :: vroe_0002
    real(f64) :: croe_0002
    real(f64) :: uleft_0002
    real(f64) :: vleft_0002
    real(f64) :: uright_0002
    real(f64) :: vright_0002
    real(f64) :: w_lrh_0002
    real(f64) :: w_lrhu_0002
    real(f64) :: w_lrhv_0002
    real(f64) :: w_lrhc_0002
    real(f64) :: w_lrz_0002
    real(f64) :: d_0002
    real(f64) :: sound_0002
    real(f64) :: Q_0002
    real(f64) :: R_0002
    real(f64) :: theta_0002
    real(f64) :: lambda1_0002
    real(f64) :: lambda2_0002
    real(f64) :: lambda3_0002
    real(f64) :: lambda4_0002
    real(f64) :: lambda5_0002
    real(f64) :: alpha1_0002
    real(f64) :: alpha2_0002
    real(f64) :: alpha3_0002
    real(f64) :: beta_0002
    real(f64) :: gamma1_0002
    real(f64) :: gamma2_0002
    real(f64) :: gamma3_0002
    real(f64) :: sigma1_0002
    real(f64) :: sigma2_0002
    real(f64) :: sigma3_0002
    real(f64) :: epsilon_0002
    real(f64) :: sign1_0002
    real(f64) :: sign2_0002
    real(f64) :: sign3_0002
    real(f64) :: sign4_0002
    real(f64) :: sign5_0002
    real(f64) :: hnew_0002
    real(f64) :: unew_0002
    real(f64) :: vnew_0002
    real(f64) :: cnew_0002
    real(f64) :: znew_0002
    real(f64) :: u_hu_0002
    real(f64) :: u_hv_0002
    real(f64) :: u_hc_0002
    real(f64) :: u_z_0002
    real(f64) :: q_s_0002
    integer(i64) :: Dummy_0009
    real(f64), allocatable :: ninv_0003(:)
    real(f64), allocatable :: w_dif_0003(:)
    real(f64), allocatable :: rmat_0003(:,:)
    integer(i64) :: As_0003
    real(f64) :: p_0003
    real(f64) :: xi_0003
    real(f64) :: u_h_0003
    real(f64) :: v_h_0003
    real(f64) :: c_h_0003
    real(f64) :: un_h_0003
    real(f64) :: vn_h_0003
    real(f64) :: hroe_0003
    real(f64) :: uroe_0003
    real(f64) :: vroe_0003
    real(f64) :: croe_0003
    real(f64) :: uleft_0003
    real(f64) :: vleft_0003
    real(f64) :: uright_0003
    real(f64) :: vright_0003
    real(f64) :: w_lrh_0003
    real(f64) :: w_lrhu_0003
    real(f64) :: w_lrhv_0003
    real(f64) :: w_lrhc_0003
    real(f64) :: w_lrz_0003
    real(f64) :: d_0003
    real(f64) :: sound_0003
    real(f64) :: Q_0003
    real(f64) :: R_0003
    real(f64) :: theta_0003
    real(f64) :: lambda1_0003
    real(f64) :: lambda2_0003
    real(f64) :: lambda3_0003
    real(f64) :: lambda4_0003
    real(f64) :: lambda5_0003
    real(f64) :: alpha1_0003
    real(f64) :: alpha2_0003
    real(f64) :: alpha3_0003
    real(f64) :: beta_0003
    real(f64) :: gamma1_0003
    real(f64) :: gamma2_0003
    real(f64) :: gamma3_0003
    real(f64) :: sigma1_0003
    real(f64) :: sigma2_0003
    real(f64) :: sigma3_0003
    real(f64) :: epsilon_0003
    real(f64) :: sign1_0003
    real(f64) :: sign2_0003
    real(f64) :: sign3_0003
    real(f64) :: sign4_0003
    real(f64) :: sign5_0003
    real(f64) :: hnew_0003
    real(f64) :: unew_0003
    real(f64) :: vnew_0003
    real(f64) :: cnew_0003
    real(f64) :: znew_0003
    real(f64) :: u_hu_0003
    real(f64) :: u_hv_0003
    real(f64) :: u_hc_0003
    real(f64) :: u_z_0003
    real(f64) :: q_s_0003

    rez_h(:) = 0.0_f64
    rez_hu(:) = 0.0_f64
    rez_hv(:) = 0.0_f64
    rez_hc(:) = 0.0_f64
    rez_Z(:) = 0.0_f64
    grav = 9.81_f64
    allocate(flux(0:4_i64))
    flux = 0.0_f64
    allocate(r_l(0:1_i64))
    r_l = 0.0_f64
    allocate(r_r(0:1_i64))
    r_r = 0.0_f64
    do Dummy_0007 = 0_i64, size(innerfaces, kind=i64) - 1_i64, 1_i64
      i = innerfaces(Dummy_0007)
      h_l = h_c(cellidf(0_i64, i))
      hu_l = hu_c(cellidf(0_i64, i))
      hv_l = hv_c(cellidf(0_i64, i))
      hc_l = hc_c(cellidf(0_i64, i))
      Z_l = Z_c(cellidf(0_i64, i))
      normal(0:) => normalf(:, i)
      mesure = mesuref(i)
      h_r = h_c(cellidf(1_i64, i))
      hu_r = hu_c(cellidf(1_i64, i))
      hv_r = hv_c(cellidf(1_i64, i))
      hc_r = hc_c(cellidf(1_i64, i))
      Z_r = Z_c(cellidf(1_i64, i))
      center_left(0:) => centerc(:, cellidf(0_i64, i))
      center_right(0:) => centerc(:, cellidf(1_i64, i))
      h_x_left = h_x(cellidf(0_i64, i))
      h_x_right = h_x(cellidf(1_i64, i))
      h_y_left = h_y(cellidf(0_i64, i))
      h_y_right = h_y(cellidf(1_i64, i))
      hc_x_left = hc_x(cellidf(0_i64, i))
      hc_x_right = hc_x(cellidf(1_i64, i))
      hc_y_left = hc_y(cellidf(0_i64, i))
      hc_y_right = hc_y(cellidf(1_i64, i))
      psi_left = psi(cellidf(0_i64, i))
      psi_right = psi(cellidf(1_i64, i))
      r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
      r_r(0_i64) = centerf(0_i64, i) - center_right(0_i64)
      r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
      r_r(1_i64) = centerf(1_i64, i) - center_right(1_i64)
      h_l = h_l + (order - 1_i64) * psi_left * (h_x_left * r_l(0_i64) + &
            h_y_left * r_l(1_i64))
      h_r = h_r + (order - 1_i64) * psi_right * (h_x_right * r_r(0_i64) &
            + h_y_right * r_r(1_i64))
      hc_l = hc_l + (order - 1_i64) * psi_left * (hc_x_left * r_l(0_i64 &
            ) + hc_y_left * r_l(1_i64))
      hc_r = hc_r + (order - 1_i64) * psi_right * (hc_x_right * r_r( &
            0_i64) + hc_y_right * r_r(1_i64))
      allocate(ninv_0001(0:1_i64))
      ninv_0001 = 0.0_f64
      allocate(w_dif_0001(0:4_i64))
      w_dif_0001 = 0.0_f64
      allocate(rmat_0001(0:4_i64, 0:4_i64))
      rmat_0001 = 0.0_f64
      As_0001 = 0_i64
      p_0001 = 0.4_f64
      xi_0001 = 1_i64 / (1_i64 - p_0001)
      ninv_0001 = 0.0_f64
      ninv_0001(0_i64) = (-1_i64) * normal(1_i64)
      ninv_0001(1_i64) = normal(0_i64)
      u_h_0001 = (hu_l / h_l * sqrt(h_l) + hu_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      v_h_0001 = (hv_l / h_l * sqrt(h_l) + hv_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      c_h_0001 = (hc_l / h_l * sqrt(h_l) + hc_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      !uvh =  array([uh, vh])
      un_h_0001 = u_h_0001 * normal(0_i64) + v_h_0001 * normal(1_i64)
      un_h_0001 = un_h_0001 / mesure
      vn_h_0001 = u_h_0001 * ninv_0001(0_i64) + v_h_0001 * ninv_0001( &
            1_i64)
      vn_h_0001 = vn_h_0001 / mesure
      hroe_0001 = (h_l + h_r) / 2_i64
      uroe_0001 = un_h_0001
      vroe_0001 = vn_h_0001
      croe_0001 = c_h_0001
      uleft_0001 = hu_l * normal(0_i64) + hv_l * normal(1_i64)
      uleft_0001 = uleft_0001 / mesure
      vleft_0001 = hu_l * ninv_0001(0_i64) + hv_l * ninv_0001(1_i64)
      vleft_0001 = vleft_0001 / mesure
      uright_0001 = hu_r * normal(0_i64) + hv_r * normal(1_i64)
      uright_0001 = uright_0001 / mesure
      vright_0001 = hu_r * ninv_0001(0_i64) + hv_r * ninv_0001(1_i64)
      vright_0001 = vright_0001 / mesure
      w_lrh_0001 = (h_l + h_r) / 2_i64
      w_lrhu_0001 = (uleft_0001 + uright_0001) / 2_i64
      w_lrhv_0001 = (vleft_0001 + vright_0001) / 2_i64
      w_lrhc_0001 = (hc_l + hc_r) / 2_i64
      w_lrz_0001 = (Z_l + Z_r) / 2_i64
      w_dif_0001(0_i64) = h_r - h_l
      w_dif_0001(1_i64) = uright_0001 - uleft_0001
      w_dif_0001(2_i64) = vright_0001 - vleft_0001
      w_dif_0001(3_i64) = hc_r - hc_l
      w_dif_0001(4_i64) = Z_r - Z_l
      d_0001 = As_0001 * xi_0001 * (3_i64 * uroe_0001 ** 2_i64 + &
            vroe_0001 ** 2_i64)
      sound_0001 = sqrt(grav * hroe_0001)
      Q_0001 = (-(uroe_0001 ** 2_i64 + 3_i64 * grav * (hroe_0001 + &
            d_0001))) / 9_i64
      R_0001 = uroe_0001 * (9_i64 * grav * (2_i64 * hroe_0001 - d_0001) &
            - 2_i64 * uroe_0001 ** 2_i64) / 54_i64
      theta_0001 = acos(R_0001 / sqrt(-Q_0001 ** 3_i64))
      !Les valeurs propres
      lambda1_0001 = 2_i64 * sqrt(-Q_0001) * cos(theta_0001 / 3_i64) + &
            2.0_f64 / 3.0_f64 * uroe_0001
      lambda2_0001 = 2_i64 * sqrt(-Q_0001) * cos((theta_0001 + 2_i64 * &
            3.141592653589793_f64) / 3_i64) + 2.0_f64 / 3.0_f64 * &
            uroe_0001
      lambda3_0001 = 2_i64 * sqrt(-Q_0001) * cos((theta_0001 + 4_i64 * &
            3.141592653589793_f64) / 3_i64) + 2.0_f64 / 3.0_f64 * &
            uroe_0001
      lambda4_0001 = uroe_0001
      lambda5_0001 = uroe_0001
      !définition de alpha
      alpha1_0001 = lambda1_0001 - uroe_0001
      alpha2_0001 = lambda2_0001 - uroe_0001
      alpha3_0001 = lambda3_0001 - uroe_0001
      !définition de beta
      beta_0001 = 2_i64 * As_0001 * xi_0001 * vroe_0001 / hroe_0001
      !définition de gamma
      gamma1_0001 = sound_0001 ** 2_i64 - uroe_0001 ** 2_i64 + &
            lambda2_0001 * lambda3_0001 - beta_0001 * alpha2_0001 * &
            alpha3_0001 * vroe_0001
      gamma2_0001 = sound_0001 ** 2_i64 - uroe_0001 ** 2_i64 + &
            lambda1_0001 * lambda3_0001 - beta_0001 * alpha1_0001 * &
            alpha3_0001 * vroe_0001
      gamma3_0001 = sound_0001 ** 2_i64 - uroe_0001 ** 2_i64 + &
            lambda1_0001 * lambda2_0001 - beta_0001 * alpha1_0001 * &
            alpha2_0001 * vroe_0001
      !définition de sigma
      sigma1_0001 = (-alpha1_0001) * alpha2_0001 + alpha2_0001 * &
            alpha3_0001 - alpha1_0001 * alpha3_0001 + alpha1_0001 ** &
            2_i64
      sigma2_0001 = alpha1_0001 * alpha2_0001 + alpha2_0001 * &
            alpha3_0001 - alpha1_0001 * alpha3_0001 - alpha2_0001 ** &
            2_i64
      sigma3_0001 = alpha1_0001 * alpha2_0001 - alpha2_0001 * &
            alpha3_0001 - alpha1_0001 * alpha3_0001 + alpha3_0001 ** &
            2_i64
      epsilon_0001 = 1e-10_f64
      if (abs(lambda1_0001) < epsilon_0001) then
        sign1_0001 = 0.0_f64
      else
        sign1_0001 = lambda1_0001 / abs(lambda1_0001)
      end if
      if (abs(lambda2_0001) < epsilon_0001) then
        sign2_0001 = 0.0_f64
      else
        sign2_0001 = lambda2_0001 / abs(lambda2_0001)
      end if
      if (abs(lambda3_0001) < epsilon_0001) then
        sign3_0001 = 0.0_f64
      else
        sign3_0001 = lambda3_0001 / abs(lambda3_0001)
      end if
      if (abs(lambda4_0001) < epsilon_0001) then
        sign4_0001 = 0.0_f64
      else
        sign4_0001 = lambda4_0001 / abs(lambda4_0001)
      end if
      if (abs(lambda5_0001) < epsilon_0001) then
        sign5_0001 = 0.0_f64
      else
        sign5_0001 = lambda5_0001 / abs(lambda5_0001)
      end if
      !1ère colonne
      rmat_0001(0_i64, 0_i64) = sign1_0001 * (gamma1_0001 / sigma1_0001 &
            ) - sign2_0001 * (gamma2_0001 / sigma2_0001) + sign3_0001 * &
            (gamma3_0001 / sigma3_0001) + sign5_0001 * (beta_0001 * &
            vroe_0001)
      rmat_0001(0_i64, 1_i64) = lambda1_0001 * sign1_0001 * (gamma1_0001 &
            / sigma1_0001) - lambda2_0001 * sign2_0001 * (gamma2_0001 / &
            sigma2_0001) + lambda3_0001 * sign3_0001 * (gamma3_0001 / &
            sigma3_0001) + sign5_0001 * (beta_0001 * uroe_0001 * &
            vroe_0001)
      rmat_0001(0_i64, 2_i64) = vroe_0001 * sign1_0001 * (gamma1_0001 / &
            sigma1_0001) - vroe_0001 * sign2_0001 * (gamma2_0001 / &
            sigma2_0001) + vroe_0001 * sign3_0001 * (gamma3_0001 / &
            sigma3_0001) - vroe_0001 * sign5_0001 * (1_i64 - beta_0001 &
            * vroe_0001)
      rmat_0001(0_i64, 3_i64) = croe_0001 * sign1_0001 * (gamma1_0001 / &
            sigma1_0001) - croe_0001 * sign2_0001 * (gamma2_0001 / &
            sigma2_0001) + croe_0001 * sign3_0001 * (gamma3_0001 / &
            sigma3_0001) - croe_0001 * sign4_0001 + croe_0001 * &
            sign5_0001 * beta_0001 * vroe_0001
      rmat_0001(0_i64, 4_i64) = (alpha1_0001 ** 2_i64 / sound_0001 ** &
            2_i64 - 1_i64) * sign1_0001 * (gamma1_0001 / sigma1_0001) - &
            (alpha2_0001 ** 2_i64 / sound_0001 ** 2_i64 - 1_i64) * &
            sign2_0001 * (gamma2_0001 / sigma2_0001) + (alpha3_0001 ** &
            2_i64 / sound_0001 ** 2_i64 - 1_i64) * sign3_0001 * ( &
            gamma3_0001 / sigma3_0001) - sign5_0001 * (beta_0001 * &
            vroe_0001)
      !2ème colonne
      rmat_0001(1_i64, 0_i64) = (-sign1_0001) * (alpha2_0001 + &
            alpha3_0001) / sigma1_0001 + sign2_0001 * (alpha1_0001 + &
            alpha3_0001) / sigma2_0001 - sign3_0001 * (alpha1_0001 + &
            alpha2_0001) / sigma3_0001
      rmat_0001(1_i64, 1_i64) = (-lambda1_0001) * sign1_0001 * ( &
            alpha2_0001 + alpha3_0001) / sigma1_0001 + lambda2_0001 * &
            sign2_0001 * (alpha1_0001 + alpha3_0001) / sigma2_0001 - &
            lambda3_0001 * sign3_0001 * (alpha1_0001 + alpha2_0001) / &
            sigma3_0001
      rmat_0001(1_i64, 2_i64) = (-vroe_0001) * sign1_0001 * (alpha2_0001 &
            + alpha3_0001) / sigma1_0001 + vroe_0001 * sign2_0001 * ( &
            alpha1_0001 + alpha3_0001) / sigma2_0001 - vroe_0001 * &
            sign3_0001 * (alpha1_0001 + alpha2_0001) / sigma3_0001
      rmat_0001(1_i64, 3_i64) = (-croe_0001) * sign1_0001 * (alpha2_0001 &
            + alpha3_0001) / sigma1_0001 + croe_0001 * sign2_0001 * ( &
            alpha1_0001 + alpha3_0001) / sigma2_0001 - croe_0001 * &
            sign3_0001 * (alpha1_0001 + alpha2_0001) / sigma3_0001
      rmat_0001(1_i64, 4_i64) = (-(alpha1_0001 ** 2_i64 / sound_0001 ** &
            2_i64 - 1_i64)) * sign1_0001 * (alpha2_0001 + alpha3_0001) &
            / sigma1_0001 + (alpha2_0001 ** 2_i64 / sound_0001 ** 2_i64 &
            - 1_i64) * sign2_0001 * (alpha1_0001 + alpha3_0001) / &
            sigma2_0001 - (alpha3_0001 ** 2_i64 / sound_0001 ** 2_i64 - &
            1_i64) * sign3_0001 * (alpha1_0001 + alpha2_0001) / &
            sigma3_0001
      !3ème colonne
      rmat_0001(2_i64, 0_i64) = sign1_0001 * beta_0001 * alpha2_0001 * &
            alpha3_0001 / sigma1_0001 - sign2_0001 * beta_0001 * &
            alpha1_0001 * alpha3_0001 / sigma2_0001 + sign3_0001 * &
            beta_0001 * alpha1_0001 * alpha2_0001 / sigma3_0001 - &
            sign5_0001 * beta_0001
      rmat_0001(2_i64, 1_i64) = lambda1_0001 * sign1_0001 * beta_0001 * &
            alpha2_0001 * alpha3_0001 / sigma1_0001 - lambda2_0001 * &
            sign2_0001 * beta_0001 * alpha1_0001 * alpha3_0001 / &
            sigma2_0001 + lambda3_0001 * sign3_0001 * beta_0001 * &
            alpha1_0001 * alpha2_0001 / sigma3_0001 - sign5_0001 * &
            beta_0001 * uroe_0001
      rmat_0001(2_i64, 2_i64) = vroe_0001 * sign1_0001 * beta_0001 * &
            alpha2_0001 * alpha3_0001 / sigma1_0001 - vroe_0001 * &
            sign2_0001 * beta_0001 * alpha1_0001 * alpha3_0001 / &
            sigma2_0001 + vroe_0001 * sign3_0001 * beta_0001 * &
            alpha1_0001 * alpha2_0001 / sigma3_0001 + sign5_0001 * ( &
            1_i64 - beta_0001 * vroe_0001)
      rmat_0001(2_i64, 3_i64) = croe_0001 * sign1_0001 * beta_0001 * &
            alpha2_0001 * alpha3_0001 / sigma1_0001 - croe_0001 * &
            sign2_0001 * beta_0001 * alpha1_0001 * alpha3_0001 / &
            sigma2_0001 + croe_0001 * sign3_0001 * beta_0001 * &
            alpha1_0001 * alpha2_0001 / sigma3_0001 - croe_0001 * &
            sign5_0001 * beta_0001
      rmat_0001(2_i64, 4_i64) = (alpha1_0001 ** 2_i64 / sound_0001 ** &
            2_i64 - 1_i64) * sign1_0001 * beta_0001 * alpha2_0001 * &
            alpha3_0001 / sigma1_0001 - (alpha2_0001 ** 2_i64 / &
            sound_0001 ** 2_i64 - 1_i64) * sign2_0001 * beta_0001 * &
            alpha1_0001 * alpha3_0001 / sigma2_0001 + (alpha3_0001 ** &
            2_i64 / sound_0001 ** 2_i64 - 1_i64) * sign3_0001 * &
            beta_0001 * alpha1_0001 * alpha2_0001 / sigma3_0001 + &
            sign5_0001 * beta_0001
      !4ème colonne
      rmat_0001(3_i64, 0_i64) = 0.0_f64
      rmat_0001(3_i64, 1_i64) = 0.0_f64
      rmat_0001(3_i64, 2_i64) = 0.0_f64
      rmat_0001(3_i64, 3_i64) = sign4_0001
      rmat_0001(3_i64, 4_i64) = 0.0_f64
      !5ème colone
      rmat_0001(4_i64, 0_i64) = sign1_0001 * sound_0001 ** 2_i64 / &
            sigma1_0001 - sign2_0001 * sound_0001 ** 2_i64 / &
            sigma2_0001 + sign3_0001 * sound_0001 ** 2_i64 / &
            sigma3_0001
      rmat_0001(4_i64, 1_i64) = lambda1_0001 * sign1_0001 * sound_0001 &
            ** 2_i64 / sigma1_0001 - lambda2_0001 * sign2_0001 * &
            sound_0001 ** 2_i64 / sigma2_0001 + lambda3_0001 * &
            sign3_0001 * sound_0001 ** 2_i64 / sigma3_0001
      rmat_0001(4_i64, 2_i64) = vroe_0001 * sign1_0001 * sound_0001 ** &
            2_i64 / sigma1_0001 - vroe_0001 * sign2_0001 * sound_0001 &
            ** 2_i64 / sigma2_0001 + vroe_0001 * sign3_0001 * &
            sound_0001 ** 2_i64 / sigma3_0001
      rmat_0001(4_i64, 3_i64) = croe_0001 * sign1_0001 * sound_0001 ** &
            2_i64 / sigma1_0001 - croe_0001 * sign2_0001 * sound_0001 &
            ** 2_i64 / sigma2_0001 + croe_0001 * sign3_0001 * &
            sound_0001 ** 2_i64 / sigma3_0001
      rmat_0001(4_i64, 4_i64) = (alpha1_0001 ** 2_i64 / sound_0001 ** &
            2_i64 - 1_i64) * sign1_0001 * sound_0001 ** 2_i64 / &
            sigma1_0001 - (alpha2_0001 ** 2_i64 / sound_0001 ** 2_i64 - &
            1_i64) * sign2_0001 * sound_0001 ** 2_i64 / sigma2_0001 + ( &
            alpha3_0001 ** 2_i64 / sound_0001 ** 2_i64 - 1_i64) * &
            sign3_0001 * sound_0001 ** 2_i64 / sigma3_0001
      hnew_0001 = sum(rmat_0001(:, 0_i64) * w_dif_0001(:))
      unew_0001 = sum(rmat_0001(:, 1_i64) * w_dif_0001(:))
      vnew_0001 = sum(rmat_0001(:, 2_i64) * w_dif_0001(:))
      cnew_0001 = sum(rmat_0001(:, 3_i64) * w_dif_0001(:))
      znew_0001 = sum(rmat_0001(:, 4_i64) * w_dif_0001(:))
      u_h_0001 = hnew_0001 / 2_i64
      u_hu_0001 = unew_0001 / 2_i64
      u_hv_0001 = vnew_0001 / 2_i64
      u_hc_0001 = cnew_0001 / 2_i64
      u_z_0001 = znew_0001 / 2_i64
      w_lrh_0001 = w_lrh_0001 - u_h_0001
      w_lrhu_0001 = w_lrhu_0001 - u_hu_0001
      w_lrhv_0001 = w_lrhv_0001 - u_hv_0001
      w_lrhc_0001 = w_lrhc_0001 - u_hc_0001
      w_lrz_0001 = w_lrz_0001 - u_z_0001
      unew_0001 = 0.0_f64
      vnew_0001 = 0.0_f64
      unew_0001 = w_lrhu_0001 * normal(0_i64) + w_lrhv_0001 * (-1_i64) * &
            normal(1_i64)
      unew_0001 = unew_0001 / mesure
      vnew_0001 = w_lrhu_0001 * (-1_i64) * ninv_0001(0_i64) + &
            w_lrhv_0001 * ninv_0001(1_i64)
      vnew_0001 = vnew_0001 / mesure
      w_lrhu_0001 = unew_0001
      w_lrhv_0001 = vnew_0001
      q_s_0001 = normal(0_i64) * unew_0001 + normal(1_i64) * vnew_0001
      flux(0_i64) = q_s_0001
      flux(1_i64) = q_s_0001 * w_lrhu_0001 / w_lrh_0001 + 0.5_f64 * grav &
            * w_lrh_0001 * w_lrh_0001 * normal(0_i64)
      flux(2_i64) = q_s_0001 * w_lrhv_0001 / w_lrh_0001 + 0.5_f64 * grav &
            * w_lrh_0001 * w_lrh_0001 * normal(1_i64)
      flux(3_i64) = q_s_0001 * w_lrhc_0001 / w_lrh_0001
      flux(4_i64) = As_0001 * xi_0001 * normal(0_i64) * unew_0001 * ( &
            unew_0001 ** 2_i64 + vnew_0001 ** 2_i64) / w_lrh_0001 ** &
            3_i64 + As_0001 * xi_0001 * normal(1_i64) * vnew_0001 * ( &
            unew_0001 ** 2_i64 + vnew_0001 ** 2_i64) / w_lrh_0001 ** &
            3_i64
      if (allocated(ninv_0001)) then
        deallocate(ninv_0001)
      end if
      if (allocated(w_dif_0001)) then
        deallocate(w_dif_0001)
      end if
      if (allocated(rmat_0001)) then
        deallocate(rmat_0001)
      end if
      rez_h(cellidf(0_i64, i)) = rez_h(cellidf(0_i64, i)) - flux(0_i64)
      rez_hu(cellidf(0_i64, i)) = rez_hu(cellidf(0_i64, i)) - flux(1_i64 &
            )
      rez_hv(cellidf(0_i64, i)) = rez_hv(cellidf(0_i64, i)) - flux(2_i64 &
            )
      rez_hc(cellidf(0_i64, i)) = rez_hc(cellidf(0_i64, i)) - flux(3_i64 &
            )
      rez_Z(cellidf(0_i64, i)) = rez_Z(cellidf(0_i64, i)) - flux(4_i64)
      rez_h(cellidf(1_i64, i)) = rez_h(cellidf(1_i64, i)) + flux(0_i64)
      rez_hu(cellidf(1_i64, i)) = rez_hu(cellidf(1_i64, i)) + flux(1_i64 &
            )
      rez_hv(cellidf(1_i64, i)) = rez_hv(cellidf(1_i64, i)) + flux(2_i64 &
            )
      rez_hc(cellidf(1_i64, i)) = rez_hc(cellidf(1_i64, i)) + flux(3_i64 &
            )
      rez_Z(cellidf(1_i64, i)) = rez_Z(cellidf(1_i64, i)) + flux(4_i64)
    end do
    do Dummy_0008 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0008)
      h_l = h_c(cellidf(0_i64, i))
      hu_l = hu_c(cellidf(0_i64, i))
      hv_l = hv_c(cellidf(0_i64, i))
      hc_l = hc_c(cellidf(0_i64, i))
      Z_l = Z_c(cellidf(0_i64, i))
      normal(0:) => normalf(:, i)
      mesure = mesuref(i)
      h_r = h_halo(halofid(i))
      hu_r = hu_halo(halofid(i))
      hv_r = hv_halo(halofid(i))
      hc_r = hc_halo(halofid(i))
      Z_r = Z_halo(halofid(i))
      center_left(0:) => centerc(:, cellidf(0_i64, i))
      center_right(0:) => centerh(:, halofid(i))
      h_x_left = h_x(cellidf(0_i64, i))
      h_x_right = hx_halo(halofid(i))
      h_y_left = h_y(cellidf(0_i64, i))
      h_y_right = hy_halo(halofid(i))
      hc_x_left = hc_x(cellidf(0_i64, i))
      hc_x_right = hcx_halo(halofid(i))
      hc_y_left = hc_y(cellidf(0_i64, i))
      hc_y_right = hcy_halo(halofid(i))
      psi_left = psi(cellidf(0_i64, i))
      psi_right = psi_halo(halofid(i))
      r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
      r_r(0_i64) = centerf(0_i64, i) - center_right(0_i64)
      r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
      r_r(1_i64) = centerf(1_i64, i) - center_right(1_i64)
      h_l = h_l + (order - 1_i64) * psi_left * (h_x_left * r_l(0_i64) + &
            h_y_left * r_l(1_i64))
      h_r = h_r + (order - 1_i64) * psi_right * (h_x_right * r_r(0_i64) &
            + h_y_right * r_r(1_i64))
      hc_l = hc_l + (order - 1_i64) * psi_left * (hc_x_left * r_l(0_i64 &
            ) + hc_y_left * r_l(1_i64))
      hc_r = hc_r + (order - 1_i64) * psi_right * (hc_x_right * r_r( &
            0_i64) + hc_y_right * r_r(1_i64))
      allocate(ninv_0002(0:1_i64))
      ninv_0002 = 0.0_f64
      allocate(w_dif_0002(0:4_i64))
      w_dif_0002 = 0.0_f64
      allocate(rmat_0002(0:4_i64, 0:4_i64))
      rmat_0002 = 0.0_f64
      As_0002 = 0_i64
      p_0002 = 0.4_f64
      xi_0002 = 1_i64 / (1_i64 - p_0002)
      ninv_0002 = 0.0_f64
      ninv_0002(0_i64) = (-1_i64) * normal(1_i64)
      ninv_0002(1_i64) = normal(0_i64)
      u_h_0002 = (hu_l / h_l * sqrt(h_l) + hu_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      v_h_0002 = (hv_l / h_l * sqrt(h_l) + hv_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      c_h_0002 = (hc_l / h_l * sqrt(h_l) + hc_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      !uvh =  array([uh, vh])
      un_h_0002 = u_h_0002 * normal(0_i64) + v_h_0002 * normal(1_i64)
      un_h_0002 = un_h_0002 / mesure
      vn_h_0002 = u_h_0002 * ninv_0002(0_i64) + v_h_0002 * ninv_0002( &
            1_i64)
      vn_h_0002 = vn_h_0002 / mesure
      hroe_0002 = (h_l + h_r) / 2_i64
      uroe_0002 = un_h_0002
      vroe_0002 = vn_h_0002
      croe_0002 = c_h_0002
      uleft_0002 = hu_l * normal(0_i64) + hv_l * normal(1_i64)
      uleft_0002 = uleft_0002 / mesure
      vleft_0002 = hu_l * ninv_0002(0_i64) + hv_l * ninv_0002(1_i64)
      vleft_0002 = vleft_0002 / mesure
      uright_0002 = hu_r * normal(0_i64) + hv_r * normal(1_i64)
      uright_0002 = uright_0002 / mesure
      vright_0002 = hu_r * ninv_0002(0_i64) + hv_r * ninv_0002(1_i64)
      vright_0002 = vright_0002 / mesure
      w_lrh_0002 = (h_l + h_r) / 2_i64
      w_lrhu_0002 = (uleft_0002 + uright_0002) / 2_i64
      w_lrhv_0002 = (vleft_0002 + vright_0002) / 2_i64
      w_lrhc_0002 = (hc_l + hc_r) / 2_i64
      w_lrz_0002 = (Z_l + Z_r) / 2_i64
      w_dif_0002(0_i64) = h_r - h_l
      w_dif_0002(1_i64) = uright_0002 - uleft_0002
      w_dif_0002(2_i64) = vright_0002 - vleft_0002
      w_dif_0002(3_i64) = hc_r - hc_l
      w_dif_0002(4_i64) = Z_r - Z_l
      d_0002 = As_0002 * xi_0002 * (3_i64 * uroe_0002 ** 2_i64 + &
            vroe_0002 ** 2_i64)
      sound_0002 = sqrt(grav * hroe_0002)
      Q_0002 = (-(uroe_0002 ** 2_i64 + 3_i64 * grav * (hroe_0002 + &
            d_0002))) / 9_i64
      R_0002 = uroe_0002 * (9_i64 * grav * (2_i64 * hroe_0002 - d_0002) &
            - 2_i64 * uroe_0002 ** 2_i64) / 54_i64
      theta_0002 = acos(R_0002 / sqrt(-Q_0002 ** 3_i64))
      !Les valeurs propres
      lambda1_0002 = 2_i64 * sqrt(-Q_0002) * cos(theta_0002 / 3_i64) + &
            2.0_f64 / 3.0_f64 * uroe_0002
      lambda2_0002 = 2_i64 * sqrt(-Q_0002) * cos((theta_0002 + 2_i64 * &
            3.141592653589793_f64) / 3_i64) + 2.0_f64 / 3.0_f64 * &
            uroe_0002
      lambda3_0002 = 2_i64 * sqrt(-Q_0002) * cos((theta_0002 + 4_i64 * &
            3.141592653589793_f64) / 3_i64) + 2.0_f64 / 3.0_f64 * &
            uroe_0002
      lambda4_0002 = uroe_0002
      lambda5_0002 = uroe_0002
      !définition de alpha
      alpha1_0002 = lambda1_0002 - uroe_0002
      alpha2_0002 = lambda2_0002 - uroe_0002
      alpha3_0002 = lambda3_0002 - uroe_0002
      !définition de beta
      beta_0002 = 2_i64 * As_0002 * xi_0002 * vroe_0002 / hroe_0002
      !définition de gamma
      gamma1_0002 = sound_0002 ** 2_i64 - uroe_0002 ** 2_i64 + &
            lambda2_0002 * lambda3_0002 - beta_0002 * alpha2_0002 * &
            alpha3_0002 * vroe_0002
      gamma2_0002 = sound_0002 ** 2_i64 - uroe_0002 ** 2_i64 + &
            lambda1_0002 * lambda3_0002 - beta_0002 * alpha1_0002 * &
            alpha3_0002 * vroe_0002
      gamma3_0002 = sound_0002 ** 2_i64 - uroe_0002 ** 2_i64 + &
            lambda1_0002 * lambda2_0002 - beta_0002 * alpha1_0002 * &
            alpha2_0002 * vroe_0002
      !définition de sigma
      sigma1_0002 = (-alpha1_0002) * alpha2_0002 + alpha2_0002 * &
            alpha3_0002 - alpha1_0002 * alpha3_0002 + alpha1_0002 ** &
            2_i64
      sigma2_0002 = alpha1_0002 * alpha2_0002 + alpha2_0002 * &
            alpha3_0002 - alpha1_0002 * alpha3_0002 - alpha2_0002 ** &
            2_i64
      sigma3_0002 = alpha1_0002 * alpha2_0002 - alpha2_0002 * &
            alpha3_0002 - alpha1_0002 * alpha3_0002 + alpha3_0002 ** &
            2_i64
      epsilon_0002 = 1e-10_f64
      if (abs(lambda1_0002) < epsilon_0002) then
        sign1_0002 = 0.0_f64
      else
        sign1_0002 = lambda1_0002 / abs(lambda1_0002)
      end if
      if (abs(lambda2_0002) < epsilon_0002) then
        sign2_0002 = 0.0_f64
      else
        sign2_0002 = lambda2_0002 / abs(lambda2_0002)
      end if
      if (abs(lambda3_0002) < epsilon_0002) then
        sign3_0002 = 0.0_f64
      else
        sign3_0002 = lambda3_0002 / abs(lambda3_0002)
      end if
      if (abs(lambda4_0002) < epsilon_0002) then
        sign4_0002 = 0.0_f64
      else
        sign4_0002 = lambda4_0002 / abs(lambda4_0002)
      end if
      if (abs(lambda5_0002) < epsilon_0002) then
        sign5_0002 = 0.0_f64
      else
        sign5_0002 = lambda5_0002 / abs(lambda5_0002)
      end if
      !1ère colonne
      rmat_0002(0_i64, 0_i64) = sign1_0002 * (gamma1_0002 / sigma1_0002 &
            ) - sign2_0002 * (gamma2_0002 / sigma2_0002) + sign3_0002 * &
            (gamma3_0002 / sigma3_0002) + sign5_0002 * (beta_0002 * &
            vroe_0002)
      rmat_0002(0_i64, 1_i64) = lambda1_0002 * sign1_0002 * (gamma1_0002 &
            / sigma1_0002) - lambda2_0002 * sign2_0002 * (gamma2_0002 / &
            sigma2_0002) + lambda3_0002 * sign3_0002 * (gamma3_0002 / &
            sigma3_0002) + sign5_0002 * (beta_0002 * uroe_0002 * &
            vroe_0002)
      rmat_0002(0_i64, 2_i64) = vroe_0002 * sign1_0002 * (gamma1_0002 / &
            sigma1_0002) - vroe_0002 * sign2_0002 * (gamma2_0002 / &
            sigma2_0002) + vroe_0002 * sign3_0002 * (gamma3_0002 / &
            sigma3_0002) - vroe_0002 * sign5_0002 * (1_i64 - beta_0002 &
            * vroe_0002)
      rmat_0002(0_i64, 3_i64) = croe_0002 * sign1_0002 * (gamma1_0002 / &
            sigma1_0002) - croe_0002 * sign2_0002 * (gamma2_0002 / &
            sigma2_0002) + croe_0002 * sign3_0002 * (gamma3_0002 / &
            sigma3_0002) - croe_0002 * sign4_0002 + croe_0002 * &
            sign5_0002 * beta_0002 * vroe_0002
      rmat_0002(0_i64, 4_i64) = (alpha1_0002 ** 2_i64 / sound_0002 ** &
            2_i64 - 1_i64) * sign1_0002 * (gamma1_0002 / sigma1_0002) - &
            (alpha2_0002 ** 2_i64 / sound_0002 ** 2_i64 - 1_i64) * &
            sign2_0002 * (gamma2_0002 / sigma2_0002) + (alpha3_0002 ** &
            2_i64 / sound_0002 ** 2_i64 - 1_i64) * sign3_0002 * ( &
            gamma3_0002 / sigma3_0002) - sign5_0002 * (beta_0002 * &
            vroe_0002)
      !2ème colonne
      rmat_0002(1_i64, 0_i64) = (-sign1_0002) * (alpha2_0002 + &
            alpha3_0002) / sigma1_0002 + sign2_0002 * (alpha1_0002 + &
            alpha3_0002) / sigma2_0002 - sign3_0002 * (alpha1_0002 + &
            alpha2_0002) / sigma3_0002
      rmat_0002(1_i64, 1_i64) = (-lambda1_0002) * sign1_0002 * ( &
            alpha2_0002 + alpha3_0002) / sigma1_0002 + lambda2_0002 * &
            sign2_0002 * (alpha1_0002 + alpha3_0002) / sigma2_0002 - &
            lambda3_0002 * sign3_0002 * (alpha1_0002 + alpha2_0002) / &
            sigma3_0002
      rmat_0002(1_i64, 2_i64) = (-vroe_0002) * sign1_0002 * (alpha2_0002 &
            + alpha3_0002) / sigma1_0002 + vroe_0002 * sign2_0002 * ( &
            alpha1_0002 + alpha3_0002) / sigma2_0002 - vroe_0002 * &
            sign3_0002 * (alpha1_0002 + alpha2_0002) / sigma3_0002
      rmat_0002(1_i64, 3_i64) = (-croe_0002) * sign1_0002 * (alpha2_0002 &
            + alpha3_0002) / sigma1_0002 + croe_0002 * sign2_0002 * ( &
            alpha1_0002 + alpha3_0002) / sigma2_0002 - croe_0002 * &
            sign3_0002 * (alpha1_0002 + alpha2_0002) / sigma3_0002
      rmat_0002(1_i64, 4_i64) = (-(alpha1_0002 ** 2_i64 / sound_0002 ** &
            2_i64 - 1_i64)) * sign1_0002 * (alpha2_0002 + alpha3_0002) &
            / sigma1_0002 + (alpha2_0002 ** 2_i64 / sound_0002 ** 2_i64 &
            - 1_i64) * sign2_0002 * (alpha1_0002 + alpha3_0002) / &
            sigma2_0002 - (alpha3_0002 ** 2_i64 / sound_0002 ** 2_i64 - &
            1_i64) * sign3_0002 * (alpha1_0002 + alpha2_0002) / &
            sigma3_0002
      !3ème colonne
      rmat_0002(2_i64, 0_i64) = sign1_0002 * beta_0002 * alpha2_0002 * &
            alpha3_0002 / sigma1_0002 - sign2_0002 * beta_0002 * &
            alpha1_0002 * alpha3_0002 / sigma2_0002 + sign3_0002 * &
            beta_0002 * alpha1_0002 * alpha2_0002 / sigma3_0002 - &
            sign5_0002 * beta_0002
      rmat_0002(2_i64, 1_i64) = lambda1_0002 * sign1_0002 * beta_0002 * &
            alpha2_0002 * alpha3_0002 / sigma1_0002 - lambda2_0002 * &
            sign2_0002 * beta_0002 * alpha1_0002 * alpha3_0002 / &
            sigma2_0002 + lambda3_0002 * sign3_0002 * beta_0002 * &
            alpha1_0002 * alpha2_0002 / sigma3_0002 - sign5_0002 * &
            beta_0002 * uroe_0002
      rmat_0002(2_i64, 2_i64) = vroe_0002 * sign1_0002 * beta_0002 * &
            alpha2_0002 * alpha3_0002 / sigma1_0002 - vroe_0002 * &
            sign2_0002 * beta_0002 * alpha1_0002 * alpha3_0002 / &
            sigma2_0002 + vroe_0002 * sign3_0002 * beta_0002 * &
            alpha1_0002 * alpha2_0002 / sigma3_0002 + sign5_0002 * ( &
            1_i64 - beta_0002 * vroe_0002)
      rmat_0002(2_i64, 3_i64) = croe_0002 * sign1_0002 * beta_0002 * &
            alpha2_0002 * alpha3_0002 / sigma1_0002 - croe_0002 * &
            sign2_0002 * beta_0002 * alpha1_0002 * alpha3_0002 / &
            sigma2_0002 + croe_0002 * sign3_0002 * beta_0002 * &
            alpha1_0002 * alpha2_0002 / sigma3_0002 - croe_0002 * &
            sign5_0002 * beta_0002
      rmat_0002(2_i64, 4_i64) = (alpha1_0002 ** 2_i64 / sound_0002 ** &
            2_i64 - 1_i64) * sign1_0002 * beta_0002 * alpha2_0002 * &
            alpha3_0002 / sigma1_0002 - (alpha2_0002 ** 2_i64 / &
            sound_0002 ** 2_i64 - 1_i64) * sign2_0002 * beta_0002 * &
            alpha1_0002 * alpha3_0002 / sigma2_0002 + (alpha3_0002 ** &
            2_i64 / sound_0002 ** 2_i64 - 1_i64) * sign3_0002 * &
            beta_0002 * alpha1_0002 * alpha2_0002 / sigma3_0002 + &
            sign5_0002 * beta_0002
      !4ème colonne
      rmat_0002(3_i64, 0_i64) = 0.0_f64
      rmat_0002(3_i64, 1_i64) = 0.0_f64
      rmat_0002(3_i64, 2_i64) = 0.0_f64
      rmat_0002(3_i64, 3_i64) = sign4_0002
      rmat_0002(3_i64, 4_i64) = 0.0_f64
      !5ème colone
      rmat_0002(4_i64, 0_i64) = sign1_0002 * sound_0002 ** 2_i64 / &
            sigma1_0002 - sign2_0002 * sound_0002 ** 2_i64 / &
            sigma2_0002 + sign3_0002 * sound_0002 ** 2_i64 / &
            sigma3_0002
      rmat_0002(4_i64, 1_i64) = lambda1_0002 * sign1_0002 * sound_0002 &
            ** 2_i64 / sigma1_0002 - lambda2_0002 * sign2_0002 * &
            sound_0002 ** 2_i64 / sigma2_0002 + lambda3_0002 * &
            sign3_0002 * sound_0002 ** 2_i64 / sigma3_0002
      rmat_0002(4_i64, 2_i64) = vroe_0002 * sign1_0002 * sound_0002 ** &
            2_i64 / sigma1_0002 - vroe_0002 * sign2_0002 * sound_0002 &
            ** 2_i64 / sigma2_0002 + vroe_0002 * sign3_0002 * &
            sound_0002 ** 2_i64 / sigma3_0002
      rmat_0002(4_i64, 3_i64) = croe_0002 * sign1_0002 * sound_0002 ** &
            2_i64 / sigma1_0002 - croe_0002 * sign2_0002 * sound_0002 &
            ** 2_i64 / sigma2_0002 + croe_0002 * sign3_0002 * &
            sound_0002 ** 2_i64 / sigma3_0002
      rmat_0002(4_i64, 4_i64) = (alpha1_0002 ** 2_i64 / sound_0002 ** &
            2_i64 - 1_i64) * sign1_0002 * sound_0002 ** 2_i64 / &
            sigma1_0002 - (alpha2_0002 ** 2_i64 / sound_0002 ** 2_i64 - &
            1_i64) * sign2_0002 * sound_0002 ** 2_i64 / sigma2_0002 + ( &
            alpha3_0002 ** 2_i64 / sound_0002 ** 2_i64 - 1_i64) * &
            sign3_0002 * sound_0002 ** 2_i64 / sigma3_0002
      hnew_0002 = sum(rmat_0002(:, 0_i64) * w_dif_0002(:))
      unew_0002 = sum(rmat_0002(:, 1_i64) * w_dif_0002(:))
      vnew_0002 = sum(rmat_0002(:, 2_i64) * w_dif_0002(:))
      cnew_0002 = sum(rmat_0002(:, 3_i64) * w_dif_0002(:))
      znew_0002 = sum(rmat_0002(:, 4_i64) * w_dif_0002(:))
      u_h_0002 = hnew_0002 / 2_i64
      u_hu_0002 = unew_0002 / 2_i64
      u_hv_0002 = vnew_0002 / 2_i64
      u_hc_0002 = cnew_0002 / 2_i64
      u_z_0002 = znew_0002 / 2_i64
      w_lrh_0002 = w_lrh_0002 - u_h_0002
      w_lrhu_0002 = w_lrhu_0002 - u_hu_0002
      w_lrhv_0002 = w_lrhv_0002 - u_hv_0002
      w_lrhc_0002 = w_lrhc_0002 - u_hc_0002
      w_lrz_0002 = w_lrz_0002 - u_z_0002
      unew_0002 = 0.0_f64
      vnew_0002 = 0.0_f64
      unew_0002 = w_lrhu_0002 * normal(0_i64) + w_lrhv_0002 * (-1_i64) * &
            normal(1_i64)
      unew_0002 = unew_0002 / mesure
      vnew_0002 = w_lrhu_0002 * (-1_i64) * ninv_0002(0_i64) + &
            w_lrhv_0002 * ninv_0002(1_i64)
      vnew_0002 = vnew_0002 / mesure
      w_lrhu_0002 = unew_0002
      w_lrhv_0002 = vnew_0002
      q_s_0002 = normal(0_i64) * unew_0002 + normal(1_i64) * vnew_0002
      flux(0_i64) = q_s_0002
      flux(1_i64) = q_s_0002 * w_lrhu_0002 / w_lrh_0002 + 0.5_f64 * grav &
            * w_lrh_0002 * w_lrh_0002 * normal(0_i64)
      flux(2_i64) = q_s_0002 * w_lrhv_0002 / w_lrh_0002 + 0.5_f64 * grav &
            * w_lrh_0002 * w_lrh_0002 * normal(1_i64)
      flux(3_i64) = q_s_0002 * w_lrhc_0002 / w_lrh_0002
      flux(4_i64) = As_0002 * xi_0002 * normal(0_i64) * unew_0002 * ( &
            unew_0002 ** 2_i64 + vnew_0002 ** 2_i64) / w_lrh_0002 ** &
            3_i64 + As_0002 * xi_0002 * normal(1_i64) * vnew_0002 * ( &
            unew_0002 ** 2_i64 + vnew_0002 ** 2_i64) / w_lrh_0002 ** &
            3_i64
      if (allocated(ninv_0002)) then
        deallocate(ninv_0002)
      end if
      if (allocated(w_dif_0002)) then
        deallocate(w_dif_0002)
      end if
      if (allocated(rmat_0002)) then
        deallocate(rmat_0002)
      end if
      rez_h(cellidf(0_i64, i)) = rez_h(cellidf(0_i64, i)) - flux(0_i64)
      rez_hu(cellidf(0_i64, i)) = rez_hu(cellidf(0_i64, i)) - flux(1_i64 &
            )
      rez_hv(cellidf(0_i64, i)) = rez_hv(cellidf(0_i64, i)) - flux(2_i64 &
            )
      rez_hc(cellidf(0_i64, i)) = rez_hc(cellidf(0_i64, i)) - flux(3_i64 &
            )
      rez_Z(cellidf(0_i64, i)) = rez_Z(cellidf(0_i64, i)) - flux(4_i64)
    end do
    do Dummy_0009 = 0_i64, size(boundaryfaces, kind=i64) - 1_i64, 1_i64
      i = boundaryfaces(Dummy_0009)
      h_l = h_c(cellidf(0_i64, i))
      hu_l = hu_c(cellidf(0_i64, i))
      hv_l = hv_c(cellidf(0_i64, i))
      hc_l = hc_c(cellidf(0_i64, i))
      Z_l = Z_c(cellidf(0_i64, i))
      normal(0:) => normalf(:, i)
      mesure = mesuref(i)
      h_r = h_ghost(i)
      hu_r = hu_ghost(i)
      hv_r = hv_ghost(i)
      hc_r = hc_ghost(i)
      Z_r = Z_ghost(i)
      center_left(0:) => centerc(:, cellidf(0_i64, i))
      center_right(0:) => centerg(:, i)
      h_x_left = h_x(cellidf(0_i64, i))
      h_y_left = h_y(cellidf(0_i64, i))
      hc_x_left = hc_x(cellidf(0_i64, i))
      hc_y_left = hc_y(cellidf(0_i64, i))
      psi_left = psi(cellidf(0_i64, i))
      r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
      r_r(0_i64) = centerf(0_i64, i) - center_right(0_i64)
      r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
      r_r(1_i64) = centerf(1_i64, i) - center_right(1_i64)
      h_l = h_l + (order - 1_i64) * psi_left * (h_x_left * r_l(0_i64) + &
            h_y_left * r_l(1_i64))
      h_r = h_r
      hc_l = hc_l + (order - 1_i64) * psi_left * (hc_x_left * r_l(0_i64 &
            ) + hc_y_left * r_l(1_i64))
      hc_r = hc_r
      allocate(ninv_0003(0:1_i64))
      ninv_0003 = 0.0_f64
      allocate(w_dif_0003(0:4_i64))
      w_dif_0003 = 0.0_f64
      allocate(rmat_0003(0:4_i64, 0:4_i64))
      rmat_0003 = 0.0_f64
      As_0003 = 0_i64
      p_0003 = 0.4_f64
      xi_0003 = 1_i64 / (1_i64 - p_0003)
      ninv_0003 = 0.0_f64
      ninv_0003(0_i64) = (-1_i64) * normal(1_i64)
      ninv_0003(1_i64) = normal(0_i64)
      u_h_0003 = (hu_l / h_l * sqrt(h_l) + hu_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      v_h_0003 = (hv_l / h_l * sqrt(h_l) + hv_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      c_h_0003 = (hc_l / h_l * sqrt(h_l) + hc_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      !uvh =  array([uh, vh])
      un_h_0003 = u_h_0003 * normal(0_i64) + v_h_0003 * normal(1_i64)
      un_h_0003 = un_h_0003 / mesure
      vn_h_0003 = u_h_0003 * ninv_0003(0_i64) + v_h_0003 * ninv_0003( &
            1_i64)
      vn_h_0003 = vn_h_0003 / mesure
      hroe_0003 = (h_l + h_r) / 2_i64
      uroe_0003 = un_h_0003
      vroe_0003 = vn_h_0003
      croe_0003 = c_h_0003
      uleft_0003 = hu_l * normal(0_i64) + hv_l * normal(1_i64)
      uleft_0003 = uleft_0003 / mesure
      vleft_0003 = hu_l * ninv_0003(0_i64) + hv_l * ninv_0003(1_i64)
      vleft_0003 = vleft_0003 / mesure
      uright_0003 = hu_r * normal(0_i64) + hv_r * normal(1_i64)
      uright_0003 = uright_0003 / mesure
      vright_0003 = hu_r * ninv_0003(0_i64) + hv_r * ninv_0003(1_i64)
      vright_0003 = vright_0003 / mesure
      w_lrh_0003 = (h_l + h_r) / 2_i64
      w_lrhu_0003 = (uleft_0003 + uright_0003) / 2_i64
      w_lrhv_0003 = (vleft_0003 + vright_0003) / 2_i64
      w_lrhc_0003 = (hc_l + hc_r) / 2_i64
      w_lrz_0003 = (Z_l + Z_r) / 2_i64
      w_dif_0003(0_i64) = h_r - h_l
      w_dif_0003(1_i64) = uright_0003 - uleft_0003
      w_dif_0003(2_i64) = vright_0003 - vleft_0003
      w_dif_0003(3_i64) = hc_r - hc_l
      w_dif_0003(4_i64) = Z_r - Z_l
      d_0003 = As_0003 * xi_0003 * (3_i64 * uroe_0003 ** 2_i64 + &
            vroe_0003 ** 2_i64)
      sound_0003 = sqrt(grav * hroe_0003)
      Q_0003 = (-(uroe_0003 ** 2_i64 + 3_i64 * grav * (hroe_0003 + &
            d_0003))) / 9_i64
      R_0003 = uroe_0003 * (9_i64 * grav * (2_i64 * hroe_0003 - d_0003) &
            - 2_i64 * uroe_0003 ** 2_i64) / 54_i64
      theta_0003 = acos(R_0003 / sqrt(-Q_0003 ** 3_i64))
      !Les valeurs propres
      lambda1_0003 = 2_i64 * sqrt(-Q_0003) * cos(theta_0003 / 3_i64) + &
            2.0_f64 / 3.0_f64 * uroe_0003
      lambda2_0003 = 2_i64 * sqrt(-Q_0003) * cos((theta_0003 + 2_i64 * &
            3.141592653589793_f64) / 3_i64) + 2.0_f64 / 3.0_f64 * &
            uroe_0003
      lambda3_0003 = 2_i64 * sqrt(-Q_0003) * cos((theta_0003 + 4_i64 * &
            3.141592653589793_f64) / 3_i64) + 2.0_f64 / 3.0_f64 * &
            uroe_0003
      lambda4_0003 = uroe_0003
      lambda5_0003 = uroe_0003
      !définition de alpha
      alpha1_0003 = lambda1_0003 - uroe_0003
      alpha2_0003 = lambda2_0003 - uroe_0003
      alpha3_0003 = lambda3_0003 - uroe_0003
      !définition de beta
      beta_0003 = 2_i64 * As_0003 * xi_0003 * vroe_0003 / hroe_0003
      !définition de gamma
      gamma1_0003 = sound_0003 ** 2_i64 - uroe_0003 ** 2_i64 + &
            lambda2_0003 * lambda3_0003 - beta_0003 * alpha2_0003 * &
            alpha3_0003 * vroe_0003
      gamma2_0003 = sound_0003 ** 2_i64 - uroe_0003 ** 2_i64 + &
            lambda1_0003 * lambda3_0003 - beta_0003 * alpha1_0003 * &
            alpha3_0003 * vroe_0003
      gamma3_0003 = sound_0003 ** 2_i64 - uroe_0003 ** 2_i64 + &
            lambda1_0003 * lambda2_0003 - beta_0003 * alpha1_0003 * &
            alpha2_0003 * vroe_0003
      !définition de sigma
      sigma1_0003 = (-alpha1_0003) * alpha2_0003 + alpha2_0003 * &
            alpha3_0003 - alpha1_0003 * alpha3_0003 + alpha1_0003 ** &
            2_i64
      sigma2_0003 = alpha1_0003 * alpha2_0003 + alpha2_0003 * &
            alpha3_0003 - alpha1_0003 * alpha3_0003 - alpha2_0003 ** &
            2_i64
      sigma3_0003 = alpha1_0003 * alpha2_0003 - alpha2_0003 * &
            alpha3_0003 - alpha1_0003 * alpha3_0003 + alpha3_0003 ** &
            2_i64
      epsilon_0003 = 1e-10_f64
      if (abs(lambda1_0003) < epsilon_0003) then
        sign1_0003 = 0.0_f64
      else
        sign1_0003 = lambda1_0003 / abs(lambda1_0003)
      end if
      if (abs(lambda2_0003) < epsilon_0003) then
        sign2_0003 = 0.0_f64
      else
        sign2_0003 = lambda2_0003 / abs(lambda2_0003)
      end if
      if (abs(lambda3_0003) < epsilon_0003) then
        sign3_0003 = 0.0_f64
      else
        sign3_0003 = lambda3_0003 / abs(lambda3_0003)
      end if
      if (abs(lambda4_0003) < epsilon_0003) then
        sign4_0003 = 0.0_f64
      else
        sign4_0003 = lambda4_0003 / abs(lambda4_0003)
      end if
      if (abs(lambda5_0003) < epsilon_0003) then
        sign5_0003 = 0.0_f64
      else
        sign5_0003 = lambda5_0003 / abs(lambda5_0003)
      end if
      !1ère colonne
      rmat_0003(0_i64, 0_i64) = sign1_0003 * (gamma1_0003 / sigma1_0003 &
            ) - sign2_0003 * (gamma2_0003 / sigma2_0003) + sign3_0003 * &
            (gamma3_0003 / sigma3_0003) + sign5_0003 * (beta_0003 * &
            vroe_0003)
      rmat_0003(0_i64, 1_i64) = lambda1_0003 * sign1_0003 * (gamma1_0003 &
            / sigma1_0003) - lambda2_0003 * sign2_0003 * (gamma2_0003 / &
            sigma2_0003) + lambda3_0003 * sign3_0003 * (gamma3_0003 / &
            sigma3_0003) + sign5_0003 * (beta_0003 * uroe_0003 * &
            vroe_0003)
      rmat_0003(0_i64, 2_i64) = vroe_0003 * sign1_0003 * (gamma1_0003 / &
            sigma1_0003) - vroe_0003 * sign2_0003 * (gamma2_0003 / &
            sigma2_0003) + vroe_0003 * sign3_0003 * (gamma3_0003 / &
            sigma3_0003) - vroe_0003 * sign5_0003 * (1_i64 - beta_0003 &
            * vroe_0003)
      rmat_0003(0_i64, 3_i64) = croe_0003 * sign1_0003 * (gamma1_0003 / &
            sigma1_0003) - croe_0003 * sign2_0003 * (gamma2_0003 / &
            sigma2_0003) + croe_0003 * sign3_0003 * (gamma3_0003 / &
            sigma3_0003) - croe_0003 * sign4_0003 + croe_0003 * &
            sign5_0003 * beta_0003 * vroe_0003
      rmat_0003(0_i64, 4_i64) = (alpha1_0003 ** 2_i64 / sound_0003 ** &
            2_i64 - 1_i64) * sign1_0003 * (gamma1_0003 / sigma1_0003) - &
            (alpha2_0003 ** 2_i64 / sound_0003 ** 2_i64 - 1_i64) * &
            sign2_0003 * (gamma2_0003 / sigma2_0003) + (alpha3_0003 ** &
            2_i64 / sound_0003 ** 2_i64 - 1_i64) * sign3_0003 * ( &
            gamma3_0003 / sigma3_0003) - sign5_0003 * (beta_0003 * &
            vroe_0003)
      !2ème colonne
      rmat_0003(1_i64, 0_i64) = (-sign1_0003) * (alpha2_0003 + &
            alpha3_0003) / sigma1_0003 + sign2_0003 * (alpha1_0003 + &
            alpha3_0003) / sigma2_0003 - sign3_0003 * (alpha1_0003 + &
            alpha2_0003) / sigma3_0003
      rmat_0003(1_i64, 1_i64) = (-lambda1_0003) * sign1_0003 * ( &
            alpha2_0003 + alpha3_0003) / sigma1_0003 + lambda2_0003 * &
            sign2_0003 * (alpha1_0003 + alpha3_0003) / sigma2_0003 - &
            lambda3_0003 * sign3_0003 * (alpha1_0003 + alpha2_0003) / &
            sigma3_0003
      rmat_0003(1_i64, 2_i64) = (-vroe_0003) * sign1_0003 * (alpha2_0003 &
            + alpha3_0003) / sigma1_0003 + vroe_0003 * sign2_0003 * ( &
            alpha1_0003 + alpha3_0003) / sigma2_0003 - vroe_0003 * &
            sign3_0003 * (alpha1_0003 + alpha2_0003) / sigma3_0003
      rmat_0003(1_i64, 3_i64) = (-croe_0003) * sign1_0003 * (alpha2_0003 &
            + alpha3_0003) / sigma1_0003 + croe_0003 * sign2_0003 * ( &
            alpha1_0003 + alpha3_0003) / sigma2_0003 - croe_0003 * &
            sign3_0003 * (alpha1_0003 + alpha2_0003) / sigma3_0003
      rmat_0003(1_i64, 4_i64) = (-(alpha1_0003 ** 2_i64 / sound_0003 ** &
            2_i64 - 1_i64)) * sign1_0003 * (alpha2_0003 + alpha3_0003) &
            / sigma1_0003 + (alpha2_0003 ** 2_i64 / sound_0003 ** 2_i64 &
            - 1_i64) * sign2_0003 * (alpha1_0003 + alpha3_0003) / &
            sigma2_0003 - (alpha3_0003 ** 2_i64 / sound_0003 ** 2_i64 - &
            1_i64) * sign3_0003 * (alpha1_0003 + alpha2_0003) / &
            sigma3_0003
      !3ème colonne
      rmat_0003(2_i64, 0_i64) = sign1_0003 * beta_0003 * alpha2_0003 * &
            alpha3_0003 / sigma1_0003 - sign2_0003 * beta_0003 * &
            alpha1_0003 * alpha3_0003 / sigma2_0003 + sign3_0003 * &
            beta_0003 * alpha1_0003 * alpha2_0003 / sigma3_0003 - &
            sign5_0003 * beta_0003
      rmat_0003(2_i64, 1_i64) = lambda1_0003 * sign1_0003 * beta_0003 * &
            alpha2_0003 * alpha3_0003 / sigma1_0003 - lambda2_0003 * &
            sign2_0003 * beta_0003 * alpha1_0003 * alpha3_0003 / &
            sigma2_0003 + lambda3_0003 * sign3_0003 * beta_0003 * &
            alpha1_0003 * alpha2_0003 / sigma3_0003 - sign5_0003 * &
            beta_0003 * uroe_0003
      rmat_0003(2_i64, 2_i64) = vroe_0003 * sign1_0003 * beta_0003 * &
            alpha2_0003 * alpha3_0003 / sigma1_0003 - vroe_0003 * &
            sign2_0003 * beta_0003 * alpha1_0003 * alpha3_0003 / &
            sigma2_0003 + vroe_0003 * sign3_0003 * beta_0003 * &
            alpha1_0003 * alpha2_0003 / sigma3_0003 + sign5_0003 * ( &
            1_i64 - beta_0003 * vroe_0003)
      rmat_0003(2_i64, 3_i64) = croe_0003 * sign1_0003 * beta_0003 * &
            alpha2_0003 * alpha3_0003 / sigma1_0003 - croe_0003 * &
            sign2_0003 * beta_0003 * alpha1_0003 * alpha3_0003 / &
            sigma2_0003 + croe_0003 * sign3_0003 * beta_0003 * &
            alpha1_0003 * alpha2_0003 / sigma3_0003 - croe_0003 * &
            sign5_0003 * beta_0003
      rmat_0003(2_i64, 4_i64) = (alpha1_0003 ** 2_i64 / sound_0003 ** &
            2_i64 - 1_i64) * sign1_0003 * beta_0003 * alpha2_0003 * &
            alpha3_0003 / sigma1_0003 - (alpha2_0003 ** 2_i64 / &
            sound_0003 ** 2_i64 - 1_i64) * sign2_0003 * beta_0003 * &
            alpha1_0003 * alpha3_0003 / sigma2_0003 + (alpha3_0003 ** &
            2_i64 / sound_0003 ** 2_i64 - 1_i64) * sign3_0003 * &
            beta_0003 * alpha1_0003 * alpha2_0003 / sigma3_0003 + &
            sign5_0003 * beta_0003
      !4ème colonne
      rmat_0003(3_i64, 0_i64) = 0.0_f64
      rmat_0003(3_i64, 1_i64) = 0.0_f64
      rmat_0003(3_i64, 2_i64) = 0.0_f64
      rmat_0003(3_i64, 3_i64) = sign4_0003
      rmat_0003(3_i64, 4_i64) = 0.0_f64
      !5ème colone
      rmat_0003(4_i64, 0_i64) = sign1_0003 * sound_0003 ** 2_i64 / &
            sigma1_0003 - sign2_0003 * sound_0003 ** 2_i64 / &
            sigma2_0003 + sign3_0003 * sound_0003 ** 2_i64 / &
            sigma3_0003
      rmat_0003(4_i64, 1_i64) = lambda1_0003 * sign1_0003 * sound_0003 &
            ** 2_i64 / sigma1_0003 - lambda2_0003 * sign2_0003 * &
            sound_0003 ** 2_i64 / sigma2_0003 + lambda3_0003 * &
            sign3_0003 * sound_0003 ** 2_i64 / sigma3_0003
      rmat_0003(4_i64, 2_i64) = vroe_0003 * sign1_0003 * sound_0003 ** &
            2_i64 / sigma1_0003 - vroe_0003 * sign2_0003 * sound_0003 &
            ** 2_i64 / sigma2_0003 + vroe_0003 * sign3_0003 * &
            sound_0003 ** 2_i64 / sigma3_0003
      rmat_0003(4_i64, 3_i64) = croe_0003 * sign1_0003 * sound_0003 ** &
            2_i64 / sigma1_0003 - croe_0003 * sign2_0003 * sound_0003 &
            ** 2_i64 / sigma2_0003 + croe_0003 * sign3_0003 * &
            sound_0003 ** 2_i64 / sigma3_0003
      rmat_0003(4_i64, 4_i64) = (alpha1_0003 ** 2_i64 / sound_0003 ** &
            2_i64 - 1_i64) * sign1_0003 * sound_0003 ** 2_i64 / &
            sigma1_0003 - (alpha2_0003 ** 2_i64 / sound_0003 ** 2_i64 - &
            1_i64) * sign2_0003 * sound_0003 ** 2_i64 / sigma2_0003 + ( &
            alpha3_0003 ** 2_i64 / sound_0003 ** 2_i64 - 1_i64) * &
            sign3_0003 * sound_0003 ** 2_i64 / sigma3_0003
      hnew_0003 = sum(rmat_0003(:, 0_i64) * w_dif_0003(:))
      unew_0003 = sum(rmat_0003(:, 1_i64) * w_dif_0003(:))
      vnew_0003 = sum(rmat_0003(:, 2_i64) * w_dif_0003(:))
      cnew_0003 = sum(rmat_0003(:, 3_i64) * w_dif_0003(:))
      znew_0003 = sum(rmat_0003(:, 4_i64) * w_dif_0003(:))
      u_h_0003 = hnew_0003 / 2_i64
      u_hu_0003 = unew_0003 / 2_i64
      u_hv_0003 = vnew_0003 / 2_i64
      u_hc_0003 = cnew_0003 / 2_i64
      u_z_0003 = znew_0003 / 2_i64
      w_lrh_0003 = w_lrh_0003 - u_h_0003
      w_lrhu_0003 = w_lrhu_0003 - u_hu_0003
      w_lrhv_0003 = w_lrhv_0003 - u_hv_0003
      w_lrhc_0003 = w_lrhc_0003 - u_hc_0003
      w_lrz_0003 = w_lrz_0003 - u_z_0003
      unew_0003 = 0.0_f64
      vnew_0003 = 0.0_f64
      unew_0003 = w_lrhu_0003 * normal(0_i64) + w_lrhv_0003 * (-1_i64) * &
            normal(1_i64)
      unew_0003 = unew_0003 / mesure
      vnew_0003 = w_lrhu_0003 * (-1_i64) * ninv_0003(0_i64) + &
            w_lrhv_0003 * ninv_0003(1_i64)
      vnew_0003 = vnew_0003 / mesure
      w_lrhu_0003 = unew_0003
      w_lrhv_0003 = vnew_0003
      q_s_0003 = normal(0_i64) * unew_0003 + normal(1_i64) * vnew_0003
      flux(0_i64) = q_s_0003
      flux(1_i64) = q_s_0003 * w_lrhu_0003 / w_lrh_0003 + 0.5_f64 * grav &
            * w_lrh_0003 * w_lrh_0003 * normal(0_i64)
      flux(2_i64) = q_s_0003 * w_lrhv_0003 / w_lrh_0003 + 0.5_f64 * grav &
            * w_lrh_0003 * w_lrh_0003 * normal(1_i64)
      flux(3_i64) = q_s_0003 * w_lrhc_0003 / w_lrh_0003
      flux(4_i64) = As_0003 * xi_0003 * normal(0_i64) * unew_0003 * ( &
            unew_0003 ** 2_i64 + vnew_0003 ** 2_i64) / w_lrh_0003 ** &
            3_i64 + As_0003 * xi_0003 * normal(1_i64) * vnew_0003 * ( &
            unew_0003 ** 2_i64 + vnew_0003 ** 2_i64) / w_lrh_0003 ** &
            3_i64
      if (allocated(ninv_0003)) then
        deallocate(ninv_0003)
      end if
      if (allocated(w_dif_0003)) then
        deallocate(w_dif_0003)
      end if
      if (allocated(rmat_0003)) then
        deallocate(rmat_0003)
      end if
      rez_h(cellidf(0_i64, i)) = rez_h(cellidf(0_i64, i)) - flux(0_i64)
      rez_hu(cellidf(0_i64, i)) = rez_hu(cellidf(0_i64, i)) - flux(1_i64 &
            )
      rez_hv(cellidf(0_i64, i)) = rez_hv(cellidf(0_i64, i)) - flux(2_i64 &
            )
      rez_hc(cellidf(0_i64, i)) = rez_hc(cellidf(0_i64, i)) - flux(3_i64 &
            )
      rez_Z(cellidf(0_i64, i)) = rez_Z(cellidf(0_i64, i)) - flux(4_i64)
    end do
    if (allocated(flux)) then
      deallocate(flux)
    end if
    if (allocated(r_l)) then
      deallocate(r_l)
    end if
    if (allocated(r_r)) then
      deallocate(r_r)
    end if

  end subroutine explicitscheme_convective_SW
  !........................................

  !........................................
  subroutine explicitscheme_convective_SWMHD(rez_h, rez_hu, rez_hv, &
        rez_hB1, rez_hB2, rez_PSI, rez_Z, h_c, hu_c, hv_c, hB1_c, hB2_c &
        , hPSIc, Z_c, h_ghost, hu_ghost, hv_ghost, hB1_ghost, hB2_ghost &
        , hPSIghost, Z_ghost, h_halo, hu_halo, hv_halo, hB1_halo, &
        hB2_halo, hPSIhalo, Z_halo, h_x, h_y, hx_halo, hy_halo, psi, &
        psi_halo, centerc, centerf, centerh, centerg, cellidf, mesuref, &
        normalf, halofid, innerfaces, halofaces, boundaryfaces, order, &
        cpsi, hB1_cst, hB2_cst)

    implicit none

    real(f64), intent(inout) :: rez_h(0:)
    real(f64), intent(inout) :: rez_hu(0:)
    real(f64), intent(inout) :: rez_hv(0:)
    real(f64), intent(inout) :: rez_hB1(0:)
    real(f64), intent(inout) :: rez_hB2(0:)
    real(f64), intent(inout) :: rez_PSI(0:)
    real(f64), intent(inout) :: rez_Z(0:)
    real(f64), intent(in) :: h_c(0:)
    real(f64), intent(in) :: hu_c(0:)
    real(f64), intent(in) :: hv_c(0:)
    real(f64), intent(in) :: hB1_c(0:)
    real(f64), intent(in) :: hB2_c(0:)
    real(f64), intent(in) :: hPSIc(0:)
    real(f64), intent(in) :: Z_c(0:)
    real(f64), intent(in) :: h_ghost(0:)
    real(f64), intent(in) :: hu_ghost(0:)
    real(f64), intent(in) :: hv_ghost(0:)
    real(f64), intent(in) :: hB1_ghost(0:)
    real(f64), intent(in) :: hB2_ghost(0:)
    real(f64), intent(in) :: hPSIghost(0:)
    real(f64), intent(in) :: Z_ghost(0:)
    real(f64), intent(in) :: h_halo(0:)
    real(f64), intent(in) :: hu_halo(0:)
    real(f64), intent(in) :: hv_halo(0:)
    real(f64), intent(in) :: hB1_halo(0:)
    real(f64), intent(in) :: hB2_halo(0:)
    real(f64), intent(in) :: hPSIhalo(0:)
    real(f64), intent(in) :: Z_halo(0:)
    real(f64), intent(in) :: h_x(0:)
    real(f64), intent(in) :: h_y(0:)
    real(f64), intent(in) :: hx_halo(0:)
    real(f64), intent(in) :: hy_halo(0:)
    real(f64), intent(in) :: psi(0:)
    real(f64), intent(in) :: psi_halo(0:)
    real(f64), intent(in), target :: centerc(0:,0:)
    real(f64), intent(in) :: centerf(0:,0:)
    real(f64), intent(in), target :: centerh(0:,0:)
    real(f64), intent(in), target :: centerg(0:,0:)
    integer(i64), intent(inout) :: cellidf(0:,0:)
    real(f64), intent(in) :: mesuref(0:)
    real(f64), intent(in), target :: normalf(0:,0:)
    integer(i64), intent(in) :: halofid(0:)
    integer(i64), intent(in) :: innerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: boundaryfaces(0:)
    integer(i64), value :: order
    real(f64), value :: cpsi
    real(f64), intent(in) :: hB1_cst(0:)
    real(f64), intent(in) :: hB2_cst(0:)
    real(f64) :: grav
    real(f64), allocatable :: flux(:)
    real(f64), allocatable :: r_l(:)
    real(f64), allocatable :: r_r(:)
    integer(i64) :: i
    real(f64) :: h_l
    real(f64) :: hu_l
    real(f64) :: hv_l
    real(f64) :: hB1_l
    real(f64) :: hB2_l
    real(f64) :: hPSI_l
    real(f64) :: Z_l
    real(f64) :: hB1c_l
    real(f64) :: hB2c_l
    real(f64), pointer :: normal(:)
    real(f64) :: mesure
    real(f64) :: h_r
    real(f64) :: hu_r
    real(f64) :: hv_r
    real(f64) :: hB1_r
    real(f64) :: hB2_r
    real(f64) :: hPSI_r
    real(f64) :: Z_r
    real(f64) :: hB1c_r
    real(f64) :: hB2c_r
    real(f64), pointer :: center_left(:)
    real(f64), pointer :: center_right(:)
    real(f64) :: h_x_left
    real(f64) :: h_x_right
    real(f64) :: h_y_left
    real(f64) :: h_y_right
    real(f64) :: psi_left
    real(f64) :: psi_right
    integer(i64) :: Dummy_0010
    real(f64), allocatable :: ninv_0004(:)
    real(f64), allocatable :: w_dif_0004(:)
    real(f64) :: u_h_0004
    real(f64) :: v_h_0004
    real(f64) :: B1_h_0001
    real(f64) :: B2_h_0001
    real(f64) :: un_h_0004
    real(f64) :: vn_h_0004
    real(f64) :: B1n_h_0001
    real(f64) :: B2n_h_0001
    real(f64) :: hroe_0004
    real(f64) :: uroe_0004
    real(f64) :: vroe_0004
    real(f64) :: B1roe_0001
    real(f64) :: B2roe_0001
    real(f64) :: uleft_0004
    real(f64) :: vleft_0004
    real(f64) :: B1left_0001
    real(f64) :: B2left_0001
    real(f64) :: uright_0004
    real(f64) :: vright_0004
    real(f64) :: B1right_0001
    real(f64) :: B2right_0001
    real(f64) :: B1i_0001
    real(f64) :: B1j_0001
    real(f64) :: absy_0001
    real(f64) :: w_lrh_0004
    real(f64) :: w_lrhu_0004
    real(f64) :: w_lrhv_0004
    real(f64) :: w_lrhB1_0001
    real(f64) :: w_lrhB2_0001
    real(f64) :: w_lrhPSI_0001
    real(f64) :: w_lrz_0004
    real(f64), allocatable, target :: signA_0001(:,:)
    real(f64) :: sound_0004
    real(f64) :: w_0001
    real(f64) :: lambda1_0004
    real(f64) :: lambda2_0004
    real(f64) :: lambda3_0004
    real(f64) :: lambda4_0004
    real(f64) :: epsilon_0004
    real(f64) :: s1_0001
    real(f64) :: pi1_0001
    real(f64) :: s2_0001
    real(f64) :: pi2_0001
    real(f64) :: s3_0001
    real(f64) :: pi3_0001
    real(f64) :: s4_0001
    real(f64) :: pi4_0001
    real(f64) :: gamma1_0004
    real(f64) :: gamma2_0004
    real(f64) :: sigma1_0004
    real(f64) :: sigma2_0004
    real(f64) :: mu1_0001
    real(f64) :: mu2_0001
    integer(i64) :: ann_0001
    real(f64), pointer :: smmat_0001(:,:)
    real(f64) :: hnew_0004
    real(f64) :: unew_0004
    real(f64) :: vnew_0004
    real(f64) :: B1new_0001
    real(f64) :: B2new_0001
    real(f64) :: znew_0004
    real(f64) :: Pnew_0001
    real(f64) :: u_hu_0004
    real(f64) :: u_hv_0004
    real(f64) :: u_hP_0001
    real(f64) :: u_hB1_0001
    real(f64) :: u_hB2_0001
    real(f64) :: u_z_0004
    real(f64) :: w_lrhP_0001
    real(f64) :: w_hP_0001
    real(f64) :: mw_hB1_0001
    real(f64) :: mhP_0001
    real(f64), allocatable :: norm_0001(:)
    real(f64) :: q_s_0004
    real(f64) :: p_s_0001
    real(f64) :: Flux_B1psi_0001
    real(f64) :: Flux_B2psi_0001
    real(f64) :: Flux_hPpsi_0001
    integer(i64) :: i_0001
    integer(i64) :: Dummy_0011
    real(f64), allocatable :: ninv_0005(:)
    real(f64), allocatable :: w_dif_0005(:)
    real(f64) :: u_h_0005
    real(f64) :: v_h_0005
    real(f64) :: B1_h_0002
    real(f64) :: B2_h_0002
    real(f64) :: un_h_0005
    real(f64) :: vn_h_0005
    real(f64) :: B1n_h_0002
    real(f64) :: B2n_h_0002
    real(f64) :: hroe_0005
    real(f64) :: uroe_0005
    real(f64) :: vroe_0005
    real(f64) :: B1roe_0002
    real(f64) :: B2roe_0002
    real(f64) :: uleft_0005
    real(f64) :: vleft_0005
    real(f64) :: B1left_0002
    real(f64) :: B2left_0002
    real(f64) :: uright_0005
    real(f64) :: vright_0005
    real(f64) :: B1right_0002
    real(f64) :: B2right_0002
    real(f64) :: B1i_0002
    real(f64) :: B1j_0002
    real(f64) :: absy_0002
    real(f64) :: w_lrh_0005
    real(f64) :: w_lrhu_0005
    real(f64) :: w_lrhv_0005
    real(f64) :: w_lrhB1_0002
    real(f64) :: w_lrhB2_0002
    real(f64) :: w_lrhPSI_0002
    real(f64) :: w_lrz_0005
    real(f64), allocatable, target :: signA_0002(:,:)
    real(f64) :: sound_0005
    real(f64) :: w_0002
    real(f64) :: lambda1_0005
    real(f64) :: lambda2_0005
    real(f64) :: lambda3_0005
    real(f64) :: lambda4_0005
    real(f64) :: epsilon_0005
    real(f64) :: s1_0002
    real(f64) :: pi1_0002
    real(f64) :: s2_0002
    real(f64) :: pi2_0002
    real(f64) :: s3_0002
    real(f64) :: pi3_0002
    real(f64) :: s4_0002
    real(f64) :: pi4_0002
    real(f64) :: gamma1_0005
    real(f64) :: gamma2_0005
    real(f64) :: sigma1_0005
    real(f64) :: sigma2_0005
    real(f64) :: mu1_0002
    real(f64) :: mu2_0002
    integer(i64) :: ann_0002
    real(f64), pointer :: smmat_0002(:,:)
    real(f64) :: hnew_0005
    real(f64) :: unew_0005
    real(f64) :: vnew_0005
    real(f64) :: B1new_0002
    real(f64) :: B2new_0002
    real(f64) :: znew_0005
    real(f64) :: Pnew_0002
    real(f64) :: u_hu_0005
    real(f64) :: u_hv_0005
    real(f64) :: u_hP_0002
    real(f64) :: u_hB1_0002
    real(f64) :: u_hB2_0002
    real(f64) :: u_z_0005
    real(f64) :: w_lrhP_0002
    real(f64) :: w_hP_0002
    real(f64) :: mw_hB1_0002
    real(f64) :: mhP_0002
    real(f64), allocatable :: norm_0002(:)
    real(f64) :: q_s_0005
    real(f64) :: p_s_0002
    real(f64) :: Flux_B1psi_0002
    real(f64) :: Flux_B2psi_0002
    real(f64) :: Flux_hPpsi_0002
    integer(i64) :: i_0002
    integer(i64) :: Dummy_0012
    real(f64), allocatable :: ninv_0006(:)
    real(f64), allocatable :: w_dif_0006(:)
    real(f64) :: u_h_0006
    real(f64) :: v_h_0006
    real(f64) :: B1_h_0003
    real(f64) :: B2_h_0003
    real(f64) :: un_h_0006
    real(f64) :: vn_h_0006
    real(f64) :: B1n_h_0003
    real(f64) :: B2n_h_0003
    real(f64) :: hroe_0006
    real(f64) :: uroe_0006
    real(f64) :: vroe_0006
    real(f64) :: B1roe_0003
    real(f64) :: B2roe_0003
    real(f64) :: uleft_0006
    real(f64) :: vleft_0006
    real(f64) :: B1left_0003
    real(f64) :: B2left_0003
    real(f64) :: uright_0006
    real(f64) :: vright_0006
    real(f64) :: B1right_0003
    real(f64) :: B2right_0003
    real(f64) :: B1i_0003
    real(f64) :: B1j_0003
    real(f64) :: absy_0003
    real(f64) :: w_lrh_0006
    real(f64) :: w_lrhu_0006
    real(f64) :: w_lrhv_0006
    real(f64) :: w_lrhB1_0003
    real(f64) :: w_lrhB2_0003
    real(f64) :: w_lrhPSI_0003
    real(f64) :: w_lrz_0006
    real(f64), allocatable, target :: signA_0003(:,:)
    real(f64) :: sound_0006
    real(f64) :: w_0003
    real(f64) :: lambda1_0006
    real(f64) :: lambda2_0006
    real(f64) :: lambda3_0006
    real(f64) :: lambda4_0006
    real(f64) :: epsilon_0006
    real(f64) :: s1_0003
    real(f64) :: pi1_0003
    real(f64) :: s2_0003
    real(f64) :: pi2_0003
    real(f64) :: s3_0003
    real(f64) :: pi3_0003
    real(f64) :: s4_0003
    real(f64) :: pi4_0003
    real(f64) :: gamma1_0006
    real(f64) :: gamma2_0006
    real(f64) :: sigma1_0006
    real(f64) :: sigma2_0006
    real(f64) :: mu1_0003
    real(f64) :: mu2_0003
    integer(i64) :: ann_0003
    real(f64), pointer :: smmat_0003(:,:)
    real(f64) :: hnew_0006
    real(f64) :: unew_0006
    real(f64) :: vnew_0006
    real(f64) :: B1new_0003
    real(f64) :: B2new_0003
    real(f64) :: znew_0006
    real(f64) :: Pnew_0003
    real(f64) :: u_hu_0006
    real(f64) :: u_hv_0006
    real(f64) :: u_hP_0003
    real(f64) :: u_hB1_0003
    real(f64) :: u_hB2_0003
    real(f64) :: u_z_0006
    real(f64) :: w_lrhP_0003
    real(f64) :: w_hP_0003
    real(f64) :: mw_hB1_0003
    real(f64) :: mhP_0003
    real(f64), allocatable :: norm_0003(:)
    real(f64) :: q_s_0006
    real(f64) :: p_s_0003
    real(f64) :: Flux_B1psi_0003
    real(f64) :: Flux_B2psi_0003
    real(f64) :: Flux_hPpsi_0003
    integer(i64) :: i_0003

    rez_h(:) = 0.0_f64
    rez_hu(:) = 0.0_f64
    rez_hv(:) = 0.0_f64
    rez_hB1(:) = 0.0_f64
    rez_hB2(:) = 0.0_f64
    rez_PSI(:) = 0.0_f64
    rez_Z(:) = 0.0_f64
    grav = 1.0_f64
    allocate(flux(0:6_i64))
    flux = 0.0_f64
    allocate(r_l(0:1_i64))
    r_l = 0.0_f64
    allocate(r_r(0:1_i64))
    r_r = 0.0_f64
    do Dummy_0010 = 0_i64, size(innerfaces, kind=i64) - 1_i64, 1_i64
      i = innerfaces(Dummy_0010)
      h_l = h_c(cellidf(0_i64, i))
      hu_l = hu_c(cellidf(0_i64, i))
      hv_l = hv_c(cellidf(0_i64, i))
      hB1_l = hB1_c(cellidf(0_i64, i))
      hB2_l = hB2_c(cellidf(0_i64, i))
      hPSI_l = hPSIc(cellidf(0_i64, i))
      Z_l = Z_c(cellidf(0_i64, i))
      hB1c_l = hB1_cst(cellidf(0_i64, i))
      hB2c_l = hB2_cst(cellidf(0_i64, i))
      normal(0:) => normalf(:, i)
      mesure = mesuref(i)
      h_r = h_c(cellidf(1_i64, i))
      hu_r = hu_c(cellidf(1_i64, i))
      hv_r = hv_c(cellidf(1_i64, i))
      hB1_r = hB1_c(cellidf(1_i64, i))
      hB2_r = hB2_c(cellidf(1_i64, i))
      hPSI_r = hPSIc(cellidf(1_i64, i))
      Z_r = Z_c(cellidf(1_i64, i))
      hB1c_r = hB1_cst(cellidf(1_i64, i))
      hB2c_r = hB2_cst(cellidf(1_i64, i))
      center_left(0:) => centerc(:, cellidf(0_i64, i))
      center_right(0:) => centerc(:, cellidf(1_i64, i))
      h_x_left = h_x(cellidf(0_i64, i))
      h_x_right = h_x(cellidf(1_i64, i))
      h_y_left = h_y(cellidf(0_i64, i))
      h_y_right = h_y(cellidf(1_i64, i))
      psi_left = psi(cellidf(0_i64, i))
      psi_right = psi(cellidf(1_i64, i))
      r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
      r_r(0_i64) = centerf(0_i64, i) - center_right(0_i64)
      r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
      r_r(1_i64) = centerf(1_i64, i) - center_right(1_i64)
      h_l = h_l + (order - 1_i64) * psi_left * (h_x_left * r_l(0_i64) + &
            h_y_left * r_l(1_i64))
      h_r = h_r + (order - 1_i64) * psi_right * (h_x_right * r_r(0_i64) &
            + h_y_right * r_r(1_i64))
      allocate(ninv_0004(0:1_i64))
      ninv_0004 = 0.0_f64
      allocate(w_dif_0004(0:5_i64))
      w_dif_0004 = 0.0_f64
      ninv_0004(0_i64) = (-1_i64) * normal(1_i64)
      ninv_0004(1_i64) = normal(0_i64)
      u_h_0004 = (hu_l / h_l * sqrt(h_l) + hu_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      v_h_0004 = (hv_l / h_l * sqrt(h_l) + hv_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      B1_h_0001 = (hB1_l / h_l * sqrt(h_l) + hB1_r / h_r * sqrt(h_r)) / &
            (sqrt(h_l) + sqrt(h_r))
      B2_h_0001 = (hB2_l / h_l * sqrt(h_l) + hB2_r / h_r * sqrt(h_r)) / &
            (sqrt(h_l) + sqrt(h_r))
      !uvh =  array([uh, vh])
      un_h_0004 = u_h_0004 * normal(0_i64) + v_h_0004 * normal(1_i64)
      un_h_0004 = un_h_0004 / mesure
      vn_h_0004 = u_h_0004 * ninv_0004(0_i64) + v_h_0004 * ninv_0004( &
            1_i64)
      vn_h_0004 = vn_h_0004 / mesure
      B1n_h_0001 = B1_h_0001 * normal(0_i64) + B2_h_0001 * normal(1_i64)
      B1n_h_0001 = B1n_h_0001 / mesure
      B2n_h_0001 = B1_h_0001 * ninv_0004(0_i64) + B2_h_0001 * ninv_0004( &
            1_i64)
      B2n_h_0001 = B2n_h_0001 / mesure
      hroe_0004 = (h_l + h_r) / 2_i64
      uroe_0004 = un_h_0004
      vroe_0004 = vn_h_0004
      B1roe_0001 = B1n_h_0001
      B2roe_0001 = B2n_h_0001
      uleft_0004 = hu_l * normal(0_i64) + hv_l * normal(1_i64)
      uleft_0004 = uleft_0004 / mesure
      vleft_0004 = hu_l * ninv_0004(0_i64) + hv_l * ninv_0004(1_i64)
      vleft_0004 = vleft_0004 / mesure
      B1left_0001 = hB1_l * normal(0_i64) + hB2_l * normal(1_i64)
      B1left_0001 = B1left_0001 / mesure
      B2left_0001 = hB1_l * ninv_0004(0_i64) + hB2_l * ninv_0004(1_i64)
      B2left_0001 = B2left_0001 / mesure
      uright_0004 = hu_r * normal(0_i64) + hv_r * normal(1_i64)
      uright_0004 = uright_0004 / mesure
      vright_0004 = hu_r * ninv_0004(0_i64) + hv_r * ninv_0004(1_i64)
      vright_0004 = vright_0004 / mesure
      B1right_0001 = hB1_r * normal(0_i64) + hB2_r * normal(1_i64)
      B1right_0001 = B1right_0001 / mesure
      B2right_0001 = hB1_r * ninv_0004(0_i64) + hB2_r * ninv_0004(1_i64)
      B2right_0001 = B2right_0001 / mesure
      B1i_0001 = (hB1c_l * normal(0_i64) + hB2c_l * normal(1_i64)) / &
            mesure
      B1j_0001 = (hB1c_r * normal(0_i64) + hB1c_r * normal(1_i64)) / &
            mesure
      absy_0001 = (B1i_0001 + B1j_0001) / 2_i64
      w_lrh_0004 = (h_l + h_r) / 2_i64
      w_lrhu_0004 = (uleft_0004 + uright_0004) / 2_i64
      w_lrhv_0004 = (vleft_0004 + vright_0004) / 2_i64
      w_lrhB1_0001 = (B1left_0001 + B1right_0001) / 2_i64
      w_lrhB2_0001 = (B2left_0001 + B2right_0001) / 2_i64
      w_lrhPSI_0001 = (hPSI_l + hPSI_r) / 2_i64
      w_lrz_0004 = (Z_l + Z_r) / 2_i64
      w_dif_0004(0_i64) = h_r - h_l
      w_dif_0004(1_i64) = uright_0004 - uleft_0004
      w_dif_0004(2_i64) = vright_0004 - vleft_0004
      w_dif_0004(3_i64) = B1right_0001 - B1left_0001
      w_dif_0004(4_i64) = B2right_0001 - B2left_0001
      w_dif_0004(5_i64) = Z_r - Z_l
      allocate(signA_0001(0:5_i64, 0:5_i64))
      signA_0001 = 0.0_f64
      sound_0004 = sqrt(grav * hroe_0004)
      B1roe_0001 = absy_0001 / hroe_0004
      w_0001 = sqrt(B1roe_0001 * B1roe_0001 + grav * hroe_0004)
      lambda1_0004 = uroe_0004 - w_0001
      lambda2_0004 = uroe_0004 - B1roe_0001
      lambda3_0004 = uroe_0004 + B1roe_0001
      lambda4_0004 = uroe_0004 + w_0001
      !cpsi = max(fabs(lambda1) , fabs(lambda4))
      epsilon_0004 = 1e-10_f64
      if (abs(lambda1_0004) < epsilon_0004) then
        s1_0001 = 0.0_f64
        pi1_0001 = 0.0_f64
      else
        s1_0001 = lambda1_0004 / abs(lambda1_0004)
        pi1_0001 = s1_0001 / lambda1_0004
      end if
      if (abs(lambda2_0004) < epsilon_0004) then
        s2_0001 = 0.0_f64
        pi2_0001 = 0.0_f64
      else
        s2_0001 = lambda2_0004 / abs(lambda2_0004)
        pi2_0001 = 1.0_f64 / abs(lambda2_0004)
      end if
      if (abs(lambda3_0004) < epsilon_0004) then
        s3_0001 = 0.0_f64
        pi3_0001 = 0.0_f64
      else
        s3_0001 = lambda3_0004 / abs(lambda3_0004)
        pi3_0001 = 1.0_f64 / abs(lambda3_0004)
      end if
      if (abs(lambda4_0004) < epsilon_0004) then
        s4_0001 = 0.0_f64
        pi4_0001 = 0.0_f64
      else
        s4_0001 = lambda4_0004 / abs(lambda4_0004)
        pi4_0001 = 1.0_f64 / abs(lambda4_0004)
      end if
      gamma1_0004 = vroe_0004 + B2roe_0001
      gamma2_0004 = vroe_0004 - B2roe_0001
      sigma1_0004 = vroe_0004 * (s1_0001 * lambda4_0004 - s4_0001 * &
            lambda1_0004) - w_0001 * (s2_0001 * gamma1_0004 + s3_0001 * &
            gamma2_0004)
      sigma2_0004 = B2roe_0001 * (s1_0001 * lambda4_0004 - s4_0001 * &
            lambda1_0004) - w_0001 * (s2_0001 * gamma1_0004 - s3_0001 * &
            gamma2_0004)
      if (abs(lambda2_0004) < epsilon_0004 .and. abs(lambda3_0004) < &
            epsilon_0004) then
        mu1_0001 = B1roe_0001 * vroe_0004 * pi1_0001 / w_0001 - &
              B1roe_0001 * vroe_0004 * pi4_0001 / w_0001
        mu2_0001 = B1roe_0001 * B2roe_0001 * pi1_0001 / w_0001 - &
              B1roe_0001 * B2roe_0001 * pi4_0001 / w_0001
        ann_0001 = 0_i64
      else
        mu1_0001 = B1roe_0001 * vroe_0004 * pi1_0001 / w_0001 - &
              B1roe_0001 * vroe_0004 * pi4_0001 / w_0001 - 0.5_f64 * ( &
              gamma1_0004 * pi2_0001 - gamma2_0004 * pi3_0001)
        mu2_0001 = B1roe_0001 * B2roe_0001 * pi1_0001 / w_0001 - &
              B1roe_0001 * B2roe_0001 * pi4_0001 / w_0001 - 0.5_f64 * ( &
              gamma1_0004 * pi2_0001 + gamma2_0004 * pi3_0001)
        ann_0001 = 0_i64
      end if
      !1ère colonne de la matrice A
      signA_0001(0_i64, 0_i64) = (s1_0001 * lambda4_0004 - s4_0001 * &
            lambda1_0004) / (2_i64 * w_0001)
      signA_0001(0_i64, 1_i64) = lambda1_0004 * lambda4_0004 * (s1_0001 &
            - s4_0001) / (2_i64 * w_0001)
      signA_0001(0_i64, 2_i64) = sigma1_0004 / (2_i64 * w_0001)
      signA_0001(0_i64, 3_i64) = 0.0_f64
      signA_0001(0_i64, 4_i64) = sigma2_0004 / (2_i64 * w_0001)
      signA_0001(0_i64, 5_i64) = 0.0_f64
      !2ème colonne de la matrice A
      signA_0001(1_i64, 0_i64) = (s4_0001 - s1_0001) / (2_i64 * w_0001)
      signA_0001(1_i64, 1_i64) = (s4_0001 * lambda4_0004 - s1_0001 * &
            lambda1_0004) / (2_i64 * w_0001)
      signA_0001(1_i64, 2_i64) = vroe_0004 * (s4_0001 - s1_0001) / ( &
            2_i64 * w_0001)
      signA_0001(1_i64, 3_i64) = 0.0_f64
      signA_0001(1_i64, 4_i64) = B2roe_0001 * (s4_0001 - s1_0001) / ( &
            2_i64 * w_0001)
      signA_0001(1_i64, 5_i64) = 0.0_f64
      !3ème colonne de la matrice A
      signA_0001(2_i64, 0_i64) = 0.0_f64
      signA_0001(2_i64, 1_i64) = 0.0_f64
      signA_0001(2_i64, 2_i64) = (s2_0001 + s3_0001) / 2_i64
      signA_0001(2_i64, 3_i64) = 0.0_f64
      signA_0001(2_i64, 4_i64) = (s2_0001 - s3_0001) / 2_i64
      signA_0001(2_i64, 5_i64) = 0.0_f64
      !4ème colonne de la matrice A
      signA_0001(3_i64, 0_i64) = ann_0001 * B1roe_0001 * (pi1_0001 - &
            pi4_0001) / w_0001
      signA_0001(3_i64, 1_i64) = ann_0001 * B1roe_0001 * (s1_0001 - &
            s4_0001) / w_0001
      signA_0001(3_i64, 2_i64) = ann_0001 * mu1_0001
      signA_0001(3_i64, 3_i64) = 0.0_f64
      signA_0001(3_i64, 4_i64) = ann_0001 * mu2_0001
      signA_0001(3_i64, 5_i64) = 0.0_f64
      !5ème colonne de la matrice A
      signA_0001(4_i64, 0_i64) = 0.0_f64
      signA_0001(4_i64, 1_i64) = 0.0_f64
      signA_0001(4_i64, 2_i64) = (s2_0001 - s3_0001) / 2_i64
      signA_0001(4_i64, 3_i64) = 0.0_f64
      signA_0001(4_i64, 4_i64) = (s2_0001 + s3_0001) / 2_i64
      signA_0001(4_i64, 5_i64) = 0.0_f64
      !6ème colonne de la matrice A
      signA_0001(5_i64, 0_i64) = sound_0004 ** 2_i64 * (pi4_0001 - &
            pi1_0001) / (2_i64 * w_0001)
      signA_0001(5_i64, 1_i64) = sound_0004 ** 2_i64 * (s4_0001 - &
            s1_0001) / (2_i64 * w_0001)
      signA_0001(5_i64, 2_i64) = sound_0004 ** 2_i64 * vroe_0004 * ( &
            pi4_0001 - pi1_0001) / (2_i64 * w_0001)
      signA_0001(5_i64, 3_i64) = 0.0_f64
      signA_0001(5_i64, 4_i64) = sound_0004 ** 2_i64 * B2roe_0001 * ( &
            pi4_0001 - pi1_0001) / (2_i64 * w_0001)
      signA_0001(5_i64, 5_i64) = 0.0_f64
      smmat_0001(0:, 0:) => signA_0001
      hnew_0004 = 0.0_f64
      unew_0004 = 0.0_f64
      vnew_0004 = 0.0_f64
      B1new_0001 = 0.0_f64
      B2new_0001 = 0.0_f64
      znew_0004 = 0.0_f64
      do i_0001 = 0_i64, 5_i64, 1_i64
        hnew_0004 = hnew_0004 + smmat_0001(i_0001, 0_i64) * w_dif_0004( &
              i_0001)
        unew_0004 = unew_0004 + smmat_0001(i_0001, 1_i64) * w_dif_0004( &
              i_0001)
        vnew_0004 = vnew_0004 + smmat_0001(i_0001, 2_i64) * w_dif_0004( &
              i_0001)
        B1new_0001 = B1new_0001 + smmat_0001(i_0001, 3_i64) * w_dif_0004 &
              (i_0001)
        B2new_0001 = B2new_0001 + smmat_0001(i_0001, 4_i64) * w_dif_0004 &
              (i_0001)
        znew_0004 = znew_0004 + smmat_0001(i_0001, 5_i64) * w_dif_0004( &
              i_0001)
      end do
      Pnew_0001 = cpsi * (B1right_0001 - B1left_0001)
      u_h_0004 = hnew_0004 / 2_i64
      u_hu_0004 = unew_0004 / 2_i64
      u_hv_0004 = vnew_0004 / 2_i64
      u_hP_0001 = Pnew_0001 / 2_i64
      u_hB1_0001 = B1new_0001 / 2_i64
      u_hB2_0001 = B2new_0001 / 2_i64
      u_z_0004 = znew_0004 / 2_i64
      w_lrh_0004 = w_lrh_0004 - u_h_0004
      w_lrhu_0004 = w_lrhu_0004 - u_hu_0004
      w_lrhv_0004 = w_lrhv_0004 - u_hv_0004
      w_lrhP_0001 = w_lrhPSI_0001 - u_hP_0001
      w_lrhB1_0001 = absy_0001
      w_lrhB2_0001 = w_lrhB2_0001 - u_hB2_0001
      w_lrz_0004 = w_lrz_0004 - u_z_0004
      w_hP_0001 = hPSI_r - hPSI_l
      mw_hB1_0001 = w_lrhB1_0001 / 2_i64 - 1_i64 / (2_i64 * cpsi) * &
            w_hP_0001
      mhP_0001 = (hPSI_r - hPSI_l) / 2_i64 - cpsi * w_dif_0004(3_i64) / &
            2_i64
      unew_0004 = 0.0_f64
      vnew_0004 = 0.0_f64
      B1new_0001 = 0.0_f64
      B2new_0001 = 0.0_f64
      unew_0004 = w_lrhu_0004 * normal(0_i64) - w_lrhv_0004 * normal( &
            1_i64)
      unew_0004 = unew_0004 / mesure
      vnew_0004 = w_lrhu_0004 * normal(1_i64) + w_lrhv_0004 * normal( &
            0_i64)
      vnew_0004 = vnew_0004 / mesure
      B1new_0001 = w_lrhB1_0001 * normal(0_i64) - w_lrhB2_0001 * normal( &
            1_i64)
      B1new_0001 = B1new_0001 / mesure
      B2new_0001 = w_lrhB1_0001 * normal(1_i64) + w_lrhB2_0001 * normal( &
            0_i64)
      B2new_0001 = B2new_0001 / mesure
      w_lrhu_0004 = unew_0004
      w_lrhv_0004 = vnew_0004
      w_lrhB1_0001 = B1new_0001
      w_lrhB2_0001 = B2new_0001
      allocate(norm_0001(0:size(normal, kind=i64) - 1_i64))
      norm_0001 = normal / mesure
      q_s_0004 = normal(0_i64) * unew_0004 + normal(1_i64) * vnew_0004
      p_s_0001 = normal(0_i64) * B1new_0001 + normal(1_i64) * B2new_0001
      Flux_B1psi_0001 = mhP_0001 * norm_0001(0_i64) * mesure
      Flux_B2psi_0001 = mhP_0001 * norm_0001(1_i64) * mesure
      Flux_hPpsi_0001 = cpsi * cpsi * mw_hB1_0001 * mesure
      flux(0_i64) = q_s_0004
      flux(1_i64) = q_s_0004 * w_lrhu_0004 / w_lrh_0004 + 0.5_f64 * grav &
            * w_lrh_0004 * w_lrh_0004 * normal(0_i64) - p_s_0001 * &
            w_lrhB1_0001 / w_lrh_0004
      flux(2_i64) = q_s_0004 * w_lrhv_0004 / w_lrh_0004 + 0.5_f64 * grav &
            * w_lrh_0004 * w_lrh_0004 * normal(1_i64) - p_s_0001 * &
            w_lrhB2_0001 / w_lrh_0004
      flux(3_i64) = (w_lrhv_0004 * w_lrhB1_0001 / w_lrh_0004 - &
            w_lrhu_0004 * w_lrhB2_0001 / w_lrh_0004) * normal(1_i64) + &
            Flux_B1psi_0001
      flux(4_i64) = (w_lrhu_0004 * w_lrhB2_0001 / w_lrh_0004 - &
            w_lrhv_0004 * w_lrhB1_0001 / w_lrh_0004) * normal(0_i64) + &
            Flux_B2psi_0001
      flux(5_i64) = Flux_hPpsi_0001
      flux(6_i64) = 0_i64
      if (allocated(ninv_0004)) then
        deallocate(ninv_0004)
      end if
      if (allocated(w_dif_0004)) then
        deallocate(w_dif_0004)
      end if
      if (allocated(signA_0001)) then
        deallocate(signA_0001)
      end if
      if (allocated(norm_0001)) then
        deallocate(norm_0001)
      end if
      rez_h(cellidf(0_i64, i)) = rez_h(cellidf(0_i64, i)) - flux(0_i64)
      rez_hu(cellidf(0_i64, i)) = rez_hu(cellidf(0_i64, i)) - flux(1_i64 &
            )
      rez_hv(cellidf(0_i64, i)) = rez_hv(cellidf(0_i64, i)) - flux(2_i64 &
            )
      rez_hB1(cellidf(0_i64, i)) = rez_hB1(cellidf(0_i64, i)) - flux( &
            3_i64)
      rez_hB2(cellidf(0_i64, i)) = rez_hB2(cellidf(0_i64, i)) - flux( &
            4_i64)
      rez_PSI(cellidf(0_i64, i)) = rez_PSI(cellidf(0_i64, i)) - flux( &
            5_i64)
      rez_Z(cellidf(0_i64, i)) = rez_Z(cellidf(0_i64, i)) - flux(6_i64)
      rez_h(cellidf(1_i64, i)) = rez_h(cellidf(1_i64, i)) + flux(0_i64)
      rez_hu(cellidf(1_i64, i)) = rez_hu(cellidf(1_i64, i)) + flux(1_i64 &
            )
      rez_hv(cellidf(1_i64, i)) = rez_hv(cellidf(1_i64, i)) + flux(2_i64 &
            )
      rez_hB1(cellidf(1_i64, i)) = rez_hB1(cellidf(1_i64, i)) + flux( &
            3_i64)
      rez_hB2(cellidf(1_i64, i)) = rez_hB2(cellidf(1_i64, i)) + flux( &
            4_i64)
      rez_PSI(cellidf(1_i64, i)) = rez_PSI(cellidf(1_i64, i)) + flux( &
            5_i64)
      rez_Z(cellidf(1_i64, i)) = rez_Z(cellidf(1_i64, i)) + flux(6_i64)
    end do
    do Dummy_0011 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0011)
      h_l = h_c(cellidf(0_i64, i))
      hu_l = hu_c(cellidf(0_i64, i))
      hv_l = hv_c(cellidf(0_i64, i))
      hB1_l = hB1_c(cellidf(0_i64, i))
      hB2_l = hB2_c(cellidf(0_i64, i))
      hPSI_l = hPSIc(cellidf(0_i64, i))
      Z_l = Z_c(cellidf(0_i64, i))
      hB1c_l = hB1_cst(cellidf(0_i64, i))
      hB2c_l = hB2_cst(cellidf(0_i64, i))
      normal(0:) => normalf(:, i)
      mesure = mesuref(i)
      h_r = h_halo(halofid(i))
      hu_r = hu_halo(halofid(i))
      hv_r = hv_halo(halofid(i))
      hB1_r = hB1_halo(halofid(i))
      hB2_r = hB2_halo(halofid(i))
      hPSI_r = hPSIhalo(halofid(i))
      Z_r = Z_halo(halofid(i))
      hB1c_r = hB1_cst(cellidf(1_i64, i))
      hB2c_r = hB2_cst(cellidf(1_i64, i))
      center_left(0:) => centerc(:, cellidf(0_i64, i))
      center_right(0:) => centerh(:, halofid(i))
      h_x_left = h_x(cellidf(0_i64, i))
      h_x_right = hx_halo(halofid(i))
      h_y_left = h_y(cellidf(0_i64, i))
      h_y_right = hy_halo(halofid(i))
      psi_left = psi(cellidf(0_i64, i))
      psi_right = psi_halo(halofid(i))
      r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
      r_r(0_i64) = centerf(0_i64, i) - center_right(0_i64)
      r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
      r_r(1_i64) = centerf(1_i64, i) - center_right(1_i64)
      h_l = h_l + (order - 1_i64) * psi_left * (h_x_left * r_l(0_i64) + &
            h_y_left * r_l(1_i64))
      h_r = h_r + (order - 1_i64) * psi_right * (h_x_right * r_r(0_i64) &
            + h_y_right * r_r(1_i64))
      allocate(ninv_0005(0:1_i64))
      ninv_0005 = 0.0_f64
      allocate(w_dif_0005(0:5_i64))
      w_dif_0005 = 0.0_f64
      ninv_0005(0_i64) = (-1_i64) * normal(1_i64)
      ninv_0005(1_i64) = normal(0_i64)
      u_h_0005 = (hu_l / h_l * sqrt(h_l) + hu_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      v_h_0005 = (hv_l / h_l * sqrt(h_l) + hv_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      B1_h_0002 = (hB1_l / h_l * sqrt(h_l) + hB1_r / h_r * sqrt(h_r)) / &
            (sqrt(h_l) + sqrt(h_r))
      B2_h_0002 = (hB2_l / h_l * sqrt(h_l) + hB2_r / h_r * sqrt(h_r)) / &
            (sqrt(h_l) + sqrt(h_r))
      !uvh =  array([uh, vh])
      un_h_0005 = u_h_0005 * normal(0_i64) + v_h_0005 * normal(1_i64)
      un_h_0005 = un_h_0005 / mesure
      vn_h_0005 = u_h_0005 * ninv_0005(0_i64) + v_h_0005 * ninv_0005( &
            1_i64)
      vn_h_0005 = vn_h_0005 / mesure
      B1n_h_0002 = B1_h_0002 * normal(0_i64) + B2_h_0002 * normal(1_i64)
      B1n_h_0002 = B1n_h_0002 / mesure
      B2n_h_0002 = B1_h_0002 * ninv_0005(0_i64) + B2_h_0002 * ninv_0005( &
            1_i64)
      B2n_h_0002 = B2n_h_0002 / mesure
      hroe_0005 = (h_l + h_r) / 2_i64
      uroe_0005 = un_h_0005
      vroe_0005 = vn_h_0005
      B1roe_0002 = B1n_h_0002
      B2roe_0002 = B2n_h_0002
      uleft_0005 = hu_l * normal(0_i64) + hv_l * normal(1_i64)
      uleft_0005 = uleft_0005 / mesure
      vleft_0005 = hu_l * ninv_0005(0_i64) + hv_l * ninv_0005(1_i64)
      vleft_0005 = vleft_0005 / mesure
      B1left_0002 = hB1_l * normal(0_i64) + hB2_l * normal(1_i64)
      B1left_0002 = B1left_0002 / mesure
      B2left_0002 = hB1_l * ninv_0005(0_i64) + hB2_l * ninv_0005(1_i64)
      B2left_0002 = B2left_0002 / mesure
      uright_0005 = hu_r * normal(0_i64) + hv_r * normal(1_i64)
      uright_0005 = uright_0005 / mesure
      vright_0005 = hu_r * ninv_0005(0_i64) + hv_r * ninv_0005(1_i64)
      vright_0005 = vright_0005 / mesure
      B1right_0002 = hB1_r * normal(0_i64) + hB2_r * normal(1_i64)
      B1right_0002 = B1right_0002 / mesure
      B2right_0002 = hB1_r * ninv_0005(0_i64) + hB2_r * ninv_0005(1_i64)
      B2right_0002 = B2right_0002 / mesure
      B1i_0002 = (hB1c_l * normal(0_i64) + hB2c_l * normal(1_i64)) / &
            mesure
      B1j_0002 = (hB1c_r * normal(0_i64) + hB1c_r * normal(1_i64)) / &
            mesure
      absy_0002 = (B1i_0002 + B1j_0002) / 2_i64
      w_lrh_0005 = (h_l + h_r) / 2_i64
      w_lrhu_0005 = (uleft_0005 + uright_0005) / 2_i64
      w_lrhv_0005 = (vleft_0005 + vright_0005) / 2_i64
      w_lrhB1_0002 = (B1left_0002 + B1right_0002) / 2_i64
      w_lrhB2_0002 = (B2left_0002 + B2right_0002) / 2_i64
      w_lrhPSI_0002 = (hPSI_l + hPSI_r) / 2_i64
      w_lrz_0005 = (Z_l + Z_r) / 2_i64
      w_dif_0005(0_i64) = h_r - h_l
      w_dif_0005(1_i64) = uright_0005 - uleft_0005
      w_dif_0005(2_i64) = vright_0005 - vleft_0005
      w_dif_0005(3_i64) = B1right_0002 - B1left_0002
      w_dif_0005(4_i64) = B2right_0002 - B2left_0002
      w_dif_0005(5_i64) = Z_r - Z_l
      allocate(signA_0002(0:5_i64, 0:5_i64))
      signA_0002 = 0.0_f64
      sound_0005 = sqrt(grav * hroe_0005)
      B1roe_0002 = absy_0002 / hroe_0005
      w_0002 = sqrt(B1roe_0002 * B1roe_0002 + grav * hroe_0005)
      lambda1_0005 = uroe_0005 - w_0002
      lambda2_0005 = uroe_0005 - B1roe_0002
      lambda3_0005 = uroe_0005 + B1roe_0002
      lambda4_0005 = uroe_0005 + w_0002
      !cpsi = max(fabs(lambda1) , fabs(lambda4))
      epsilon_0005 = 1e-10_f64
      if (abs(lambda1_0005) < epsilon_0005) then
        s1_0002 = 0.0_f64
        pi1_0002 = 0.0_f64
      else
        s1_0002 = lambda1_0005 / abs(lambda1_0005)
        pi1_0002 = s1_0002 / lambda1_0005
      end if
      if (abs(lambda2_0005) < epsilon_0005) then
        s2_0002 = 0.0_f64
        pi2_0002 = 0.0_f64
      else
        s2_0002 = lambda2_0005 / abs(lambda2_0005)
        pi2_0002 = 1.0_f64 / abs(lambda2_0005)
      end if
      if (abs(lambda3_0005) < epsilon_0005) then
        s3_0002 = 0.0_f64
        pi3_0002 = 0.0_f64
      else
        s3_0002 = lambda3_0005 / abs(lambda3_0005)
        pi3_0002 = 1.0_f64 / abs(lambda3_0005)
      end if
      if (abs(lambda4_0005) < epsilon_0005) then
        s4_0002 = 0.0_f64
        pi4_0002 = 0.0_f64
      else
        s4_0002 = lambda4_0005 / abs(lambda4_0005)
        pi4_0002 = 1.0_f64 / abs(lambda4_0005)
      end if
      gamma1_0005 = vroe_0005 + B2roe_0002
      gamma2_0005 = vroe_0005 - B2roe_0002
      sigma1_0005 = vroe_0005 * (s1_0002 * lambda4_0005 - s4_0002 * &
            lambda1_0005) - w_0002 * (s2_0002 * gamma1_0005 + s3_0002 * &
            gamma2_0005)
      sigma2_0005 = B2roe_0002 * (s1_0002 * lambda4_0005 - s4_0002 * &
            lambda1_0005) - w_0002 * (s2_0002 * gamma1_0005 - s3_0002 * &
            gamma2_0005)
      if (abs(lambda2_0005) < epsilon_0005 .and. abs(lambda3_0005) < &
            epsilon_0005) then
        mu1_0002 = B1roe_0002 * vroe_0005 * pi1_0002 / w_0002 - &
              B1roe_0002 * vroe_0005 * pi4_0002 / w_0002
        mu2_0002 = B1roe_0002 * B2roe_0002 * pi1_0002 / w_0002 - &
              B1roe_0002 * B2roe_0002 * pi4_0002 / w_0002
        ann_0002 = 0_i64
      else
        mu1_0002 = B1roe_0002 * vroe_0005 * pi1_0002 / w_0002 - &
              B1roe_0002 * vroe_0005 * pi4_0002 / w_0002 - 0.5_f64 * ( &
              gamma1_0005 * pi2_0002 - gamma2_0005 * pi3_0002)
        mu2_0002 = B1roe_0002 * B2roe_0002 * pi1_0002 / w_0002 - &
              B1roe_0002 * B2roe_0002 * pi4_0002 / w_0002 - 0.5_f64 * ( &
              gamma1_0005 * pi2_0002 + gamma2_0005 * pi3_0002)
        ann_0002 = 0_i64
      end if
      !1ère colonne de la matrice A
      signA_0002(0_i64, 0_i64) = (s1_0002 * lambda4_0005 - s4_0002 * &
            lambda1_0005) / (2_i64 * w_0002)
      signA_0002(0_i64, 1_i64) = lambda1_0005 * lambda4_0005 * (s1_0002 &
            - s4_0002) / (2_i64 * w_0002)
      signA_0002(0_i64, 2_i64) = sigma1_0005 / (2_i64 * w_0002)
      signA_0002(0_i64, 3_i64) = 0.0_f64
      signA_0002(0_i64, 4_i64) = sigma2_0005 / (2_i64 * w_0002)
      signA_0002(0_i64, 5_i64) = 0.0_f64
      !2ème colonne de la matrice A
      signA_0002(1_i64, 0_i64) = (s4_0002 - s1_0002) / (2_i64 * w_0002)
      signA_0002(1_i64, 1_i64) = (s4_0002 * lambda4_0005 - s1_0002 * &
            lambda1_0005) / (2_i64 * w_0002)
      signA_0002(1_i64, 2_i64) = vroe_0005 * (s4_0002 - s1_0002) / ( &
            2_i64 * w_0002)
      signA_0002(1_i64, 3_i64) = 0.0_f64
      signA_0002(1_i64, 4_i64) = B2roe_0002 * (s4_0002 - s1_0002) / ( &
            2_i64 * w_0002)
      signA_0002(1_i64, 5_i64) = 0.0_f64
      !3ème colonne de la matrice A
      signA_0002(2_i64, 0_i64) = 0.0_f64
      signA_0002(2_i64, 1_i64) = 0.0_f64
      signA_0002(2_i64, 2_i64) = (s2_0002 + s3_0002) / 2_i64
      signA_0002(2_i64, 3_i64) = 0.0_f64
      signA_0002(2_i64, 4_i64) = (s2_0002 - s3_0002) / 2_i64
      signA_0002(2_i64, 5_i64) = 0.0_f64
      !4ème colonne de la matrice A
      signA_0002(3_i64, 0_i64) = ann_0002 * B1roe_0002 * (pi1_0002 - &
            pi4_0002) / w_0002
      signA_0002(3_i64, 1_i64) = ann_0002 * B1roe_0002 * (s1_0002 - &
            s4_0002) / w_0002
      signA_0002(3_i64, 2_i64) = ann_0002 * mu1_0002
      signA_0002(3_i64, 3_i64) = 0.0_f64
      signA_0002(3_i64, 4_i64) = ann_0002 * mu2_0002
      signA_0002(3_i64, 5_i64) = 0.0_f64
      !5ème colonne de la matrice A
      signA_0002(4_i64, 0_i64) = 0.0_f64
      signA_0002(4_i64, 1_i64) = 0.0_f64
      signA_0002(4_i64, 2_i64) = (s2_0002 - s3_0002) / 2_i64
      signA_0002(4_i64, 3_i64) = 0.0_f64
      signA_0002(4_i64, 4_i64) = (s2_0002 + s3_0002) / 2_i64
      signA_0002(4_i64, 5_i64) = 0.0_f64
      !6ème colonne de la matrice A
      signA_0002(5_i64, 0_i64) = sound_0005 ** 2_i64 * (pi4_0002 - &
            pi1_0002) / (2_i64 * w_0002)
      signA_0002(5_i64, 1_i64) = sound_0005 ** 2_i64 * (s4_0002 - &
            s1_0002) / (2_i64 * w_0002)
      signA_0002(5_i64, 2_i64) = sound_0005 ** 2_i64 * vroe_0005 * ( &
            pi4_0002 - pi1_0002) / (2_i64 * w_0002)
      signA_0002(5_i64, 3_i64) = 0.0_f64
      signA_0002(5_i64, 4_i64) = sound_0005 ** 2_i64 * B2roe_0002 * ( &
            pi4_0002 - pi1_0002) / (2_i64 * w_0002)
      signA_0002(5_i64, 5_i64) = 0.0_f64
      smmat_0002(0:, 0:) => signA_0002
      hnew_0005 = 0.0_f64
      unew_0005 = 0.0_f64
      vnew_0005 = 0.0_f64
      B1new_0002 = 0.0_f64
      B2new_0002 = 0.0_f64
      znew_0005 = 0.0_f64
      do i_0002 = 0_i64, 5_i64, 1_i64
        hnew_0005 = hnew_0005 + smmat_0002(i_0002, 0_i64) * w_dif_0005( &
              i_0002)
        unew_0005 = unew_0005 + smmat_0002(i_0002, 1_i64) * w_dif_0005( &
              i_0002)
        vnew_0005 = vnew_0005 + smmat_0002(i_0002, 2_i64) * w_dif_0005( &
              i_0002)
        B1new_0002 = B1new_0002 + smmat_0002(i_0002, 3_i64) * w_dif_0005 &
              (i_0002)
        B2new_0002 = B2new_0002 + smmat_0002(i_0002, 4_i64) * w_dif_0005 &
              (i_0002)
        znew_0005 = znew_0005 + smmat_0002(i_0002, 5_i64) * w_dif_0005( &
              i_0002)
      end do
      Pnew_0002 = cpsi * (B1right_0002 - B1left_0002)
      u_h_0005 = hnew_0005 / 2_i64
      u_hu_0005 = unew_0005 / 2_i64
      u_hv_0005 = vnew_0005 / 2_i64
      u_hP_0002 = Pnew_0002 / 2_i64
      u_hB1_0002 = B1new_0002 / 2_i64
      u_hB2_0002 = B2new_0002 / 2_i64
      u_z_0005 = znew_0005 / 2_i64
      w_lrh_0005 = w_lrh_0005 - u_h_0005
      w_lrhu_0005 = w_lrhu_0005 - u_hu_0005
      w_lrhv_0005 = w_lrhv_0005 - u_hv_0005
      w_lrhP_0002 = w_lrhPSI_0002 - u_hP_0002
      w_lrhB1_0002 = absy_0002
      w_lrhB2_0002 = w_lrhB2_0002 - u_hB2_0002
      w_lrz_0005 = w_lrz_0005 - u_z_0005
      w_hP_0002 = hPSI_r - hPSI_l
      mw_hB1_0002 = w_lrhB1_0002 / 2_i64 - 1_i64 / (2_i64 * cpsi) * &
            w_hP_0002
      mhP_0002 = (hPSI_r - hPSI_l) / 2_i64 - cpsi * w_dif_0005(3_i64) / &
            2_i64
      unew_0005 = 0.0_f64
      vnew_0005 = 0.0_f64
      B1new_0002 = 0.0_f64
      B2new_0002 = 0.0_f64
      unew_0005 = w_lrhu_0005 * normal(0_i64) - w_lrhv_0005 * normal( &
            1_i64)
      unew_0005 = unew_0005 / mesure
      vnew_0005 = w_lrhu_0005 * normal(1_i64) + w_lrhv_0005 * normal( &
            0_i64)
      vnew_0005 = vnew_0005 / mesure
      B1new_0002 = w_lrhB1_0002 * normal(0_i64) - w_lrhB2_0002 * normal( &
            1_i64)
      B1new_0002 = B1new_0002 / mesure
      B2new_0002 = w_lrhB1_0002 * normal(1_i64) + w_lrhB2_0002 * normal( &
            0_i64)
      B2new_0002 = B2new_0002 / mesure
      w_lrhu_0005 = unew_0005
      w_lrhv_0005 = vnew_0005
      w_lrhB1_0002 = B1new_0002
      w_lrhB2_0002 = B2new_0002
      allocate(norm_0002(0:size(normal, kind=i64) - 1_i64))
      norm_0002 = normal / mesure
      q_s_0005 = normal(0_i64) * unew_0005 + normal(1_i64) * vnew_0005
      p_s_0002 = normal(0_i64) * B1new_0002 + normal(1_i64) * B2new_0002
      Flux_B1psi_0002 = mhP_0002 * norm_0002(0_i64) * mesure
      Flux_B2psi_0002 = mhP_0002 * norm_0002(1_i64) * mesure
      Flux_hPpsi_0002 = cpsi * cpsi * mw_hB1_0002 * mesure
      flux(0_i64) = q_s_0005
      flux(1_i64) = q_s_0005 * w_lrhu_0005 / w_lrh_0005 + 0.5_f64 * grav &
            * w_lrh_0005 * w_lrh_0005 * normal(0_i64) - p_s_0002 * &
            w_lrhB1_0002 / w_lrh_0005
      flux(2_i64) = q_s_0005 * w_lrhv_0005 / w_lrh_0005 + 0.5_f64 * grav &
            * w_lrh_0005 * w_lrh_0005 * normal(1_i64) - p_s_0002 * &
            w_lrhB2_0002 / w_lrh_0005
      flux(3_i64) = (w_lrhv_0005 * w_lrhB1_0002 / w_lrh_0005 - &
            w_lrhu_0005 * w_lrhB2_0002 / w_lrh_0005) * normal(1_i64) + &
            Flux_B1psi_0002
      flux(4_i64) = (w_lrhu_0005 * w_lrhB2_0002 / w_lrh_0005 - &
            w_lrhv_0005 * w_lrhB1_0002 / w_lrh_0005) * normal(0_i64) + &
            Flux_B2psi_0002
      flux(5_i64) = Flux_hPpsi_0002
      flux(6_i64) = 0_i64
      if (allocated(ninv_0005)) then
        deallocate(ninv_0005)
      end if
      if (allocated(w_dif_0005)) then
        deallocate(w_dif_0005)
      end if
      if (allocated(signA_0002)) then
        deallocate(signA_0002)
      end if
      if (allocated(norm_0002)) then
        deallocate(norm_0002)
      end if
      rez_h(cellidf(0_i64, i)) = rez_h(cellidf(0_i64, i)) - flux(0_i64)
      rez_hu(cellidf(0_i64, i)) = rez_hu(cellidf(0_i64, i)) - flux(1_i64 &
            )
      rez_hv(cellidf(0_i64, i)) = rez_hv(cellidf(0_i64, i)) - flux(2_i64 &
            )
      rez_hB1(cellidf(0_i64, i)) = rez_hB1(cellidf(0_i64, i)) - flux( &
            3_i64)
      rez_hB2(cellidf(0_i64, i)) = rez_hB2(cellidf(0_i64, i)) - flux( &
            4_i64)
      rez_PSI(cellidf(0_i64, i)) = rez_PSI(cellidf(0_i64, i)) - flux( &
            5_i64)
      rez_Z(cellidf(0_i64, i)) = rez_Z(cellidf(0_i64, i)) - flux(6_i64)
    end do
    do Dummy_0012 = 0_i64, size(boundaryfaces, kind=i64) - 1_i64, 1_i64
      i = boundaryfaces(Dummy_0012)
      h_l = h_c(cellidf(0_i64, i))
      hu_l = hu_c(cellidf(0_i64, i))
      hv_l = hv_c(cellidf(0_i64, i))
      hB1_l = hB1_c(cellidf(0_i64, i))
      hB2_l = hB2_c(cellidf(0_i64, i))
      hPSI_l = hPSIc(cellidf(0_i64, i))
      Z_l = Z_c(cellidf(0_i64, i))
      hB1c_l = hB1_cst(cellidf(0_i64, i))
      hB2c_l = hB2_cst(cellidf(0_i64, i))
      normal(0:) => normalf(:, i)
      mesure = mesuref(i)
      h_r = h_ghost(i)
      hu_r = hu_ghost(i)
      hv_r = hv_ghost(i)
      hB1_r = hB1_ghost(i)
      hB2_r = hB2_ghost(i)
      hPSI_r = hPSIghost(i)
      Z_r = Z_ghost(i)
      hB1c_r = hB1_cst(cellidf(1_i64, i))
      hB2c_r = hB2_cst(cellidf(1_i64, i))
      center_left(0:) => centerc(:, cellidf(0_i64, i))
      center_right(0:) => centerg(:, i)
      h_x_left = h_x(cellidf(0_i64, i))
      h_y_left = h_y(cellidf(0_i64, i))
      psi_left = psi(cellidf(0_i64, i))
      r_l(0_i64) = centerf(0_i64, i) - center_left(0_i64)
      r_r(0_i64) = centerf(0_i64, i) - center_right(0_i64)
      r_l(1_i64) = centerf(1_i64, i) - center_left(1_i64)
      r_r(1_i64) = centerf(1_i64, i) - center_right(1_i64)
      h_l = h_l + (order - 1_i64) * psi_left * (h_x_left * r_l(0_i64) + &
            h_y_left * r_l(1_i64))
      h_r = h_r
      allocate(ninv_0006(0:1_i64))
      ninv_0006 = 0.0_f64
      allocate(w_dif_0006(0:5_i64))
      w_dif_0006 = 0.0_f64
      ninv_0006(0_i64) = (-1_i64) * normal(1_i64)
      ninv_0006(1_i64) = normal(0_i64)
      u_h_0006 = (hu_l / h_l * sqrt(h_l) + hu_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      v_h_0006 = (hv_l / h_l * sqrt(h_l) + hv_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      B1_h_0003 = (hB1_l / h_l * sqrt(h_l) + hB1_r / h_r * sqrt(h_r)) / &
            (sqrt(h_l) + sqrt(h_r))
      B2_h_0003 = (hB2_l / h_l * sqrt(h_l) + hB2_r / h_r * sqrt(h_r)) / &
            (sqrt(h_l) + sqrt(h_r))
      !uvh =  array([uh, vh])
      un_h_0006 = u_h_0006 * normal(0_i64) + v_h_0006 * normal(1_i64)
      un_h_0006 = un_h_0006 / mesure
      vn_h_0006 = u_h_0006 * ninv_0006(0_i64) + v_h_0006 * ninv_0006( &
            1_i64)
      vn_h_0006 = vn_h_0006 / mesure
      B1n_h_0003 = B1_h_0003 * normal(0_i64) + B2_h_0003 * normal(1_i64)
      B1n_h_0003 = B1n_h_0003 / mesure
      B2n_h_0003 = B1_h_0003 * ninv_0006(0_i64) + B2_h_0003 * ninv_0006( &
            1_i64)
      B2n_h_0003 = B2n_h_0003 / mesure
      hroe_0006 = (h_l + h_r) / 2_i64
      uroe_0006 = un_h_0006
      vroe_0006 = vn_h_0006
      B1roe_0003 = B1n_h_0003
      B2roe_0003 = B2n_h_0003
      uleft_0006 = hu_l * normal(0_i64) + hv_l * normal(1_i64)
      uleft_0006 = uleft_0006 / mesure
      vleft_0006 = hu_l * ninv_0006(0_i64) + hv_l * ninv_0006(1_i64)
      vleft_0006 = vleft_0006 / mesure
      B1left_0003 = hB1_l * normal(0_i64) + hB2_l * normal(1_i64)
      B1left_0003 = B1left_0003 / mesure
      B2left_0003 = hB1_l * ninv_0006(0_i64) + hB2_l * ninv_0006(1_i64)
      B2left_0003 = B2left_0003 / mesure
      uright_0006 = hu_r * normal(0_i64) + hv_r * normal(1_i64)
      uright_0006 = uright_0006 / mesure
      vright_0006 = hu_r * ninv_0006(0_i64) + hv_r * ninv_0006(1_i64)
      vright_0006 = vright_0006 / mesure
      B1right_0003 = hB1_r * normal(0_i64) + hB2_r * normal(1_i64)
      B1right_0003 = B1right_0003 / mesure
      B2right_0003 = hB1_r * ninv_0006(0_i64) + hB2_r * ninv_0006(1_i64)
      B2right_0003 = B2right_0003 / mesure
      B1i_0003 = (hB1c_l * normal(0_i64) + hB2c_l * normal(1_i64)) / &
            mesure
      B1j_0003 = (hB1c_r * normal(0_i64) + hB1c_r * normal(1_i64)) / &
            mesure
      absy_0003 = (B1i_0003 + B1j_0003) / 2_i64
      w_lrh_0006 = (h_l + h_r) / 2_i64
      w_lrhu_0006 = (uleft_0006 + uright_0006) / 2_i64
      w_lrhv_0006 = (vleft_0006 + vright_0006) / 2_i64
      w_lrhB1_0003 = (B1left_0003 + B1right_0003) / 2_i64
      w_lrhB2_0003 = (B2left_0003 + B2right_0003) / 2_i64
      w_lrhPSI_0003 = (hPSI_l + hPSI_r) / 2_i64
      w_lrz_0006 = (Z_l + Z_r) / 2_i64
      w_dif_0006(0_i64) = h_r - h_l
      w_dif_0006(1_i64) = uright_0006 - uleft_0006
      w_dif_0006(2_i64) = vright_0006 - vleft_0006
      w_dif_0006(3_i64) = B1right_0003 - B1left_0003
      w_dif_0006(4_i64) = B2right_0003 - B2left_0003
      w_dif_0006(5_i64) = Z_r - Z_l
      allocate(signA_0003(0:5_i64, 0:5_i64))
      signA_0003 = 0.0_f64
      sound_0006 = sqrt(grav * hroe_0006)
      B1roe_0003 = absy_0003 / hroe_0006
      w_0003 = sqrt(B1roe_0003 * B1roe_0003 + grav * hroe_0006)
      lambda1_0006 = uroe_0006 - w_0003
      lambda2_0006 = uroe_0006 - B1roe_0003
      lambda3_0006 = uroe_0006 + B1roe_0003
      lambda4_0006 = uroe_0006 + w_0003
      !cpsi = max(fabs(lambda1) , fabs(lambda4))
      epsilon_0006 = 1e-10_f64
      if (abs(lambda1_0006) < epsilon_0006) then
        s1_0003 = 0.0_f64
        pi1_0003 = 0.0_f64
      else
        s1_0003 = lambda1_0006 / abs(lambda1_0006)
        pi1_0003 = s1_0003 / lambda1_0006
      end if
      if (abs(lambda2_0006) < epsilon_0006) then
        s2_0003 = 0.0_f64
        pi2_0003 = 0.0_f64
      else
        s2_0003 = lambda2_0006 / abs(lambda2_0006)
        pi2_0003 = 1.0_f64 / abs(lambda2_0006)
      end if
      if (abs(lambda3_0006) < epsilon_0006) then
        s3_0003 = 0.0_f64
        pi3_0003 = 0.0_f64
      else
        s3_0003 = lambda3_0006 / abs(lambda3_0006)
        pi3_0003 = 1.0_f64 / abs(lambda3_0006)
      end if
      if (abs(lambda4_0006) < epsilon_0006) then
        s4_0003 = 0.0_f64
        pi4_0003 = 0.0_f64
      else
        s4_0003 = lambda4_0006 / abs(lambda4_0006)
        pi4_0003 = 1.0_f64 / abs(lambda4_0006)
      end if
      gamma1_0006 = vroe_0006 + B2roe_0003
      gamma2_0006 = vroe_0006 - B2roe_0003
      sigma1_0006 = vroe_0006 * (s1_0003 * lambda4_0006 - s4_0003 * &
            lambda1_0006) - w_0003 * (s2_0003 * gamma1_0006 + s3_0003 * &
            gamma2_0006)
      sigma2_0006 = B2roe_0003 * (s1_0003 * lambda4_0006 - s4_0003 * &
            lambda1_0006) - w_0003 * (s2_0003 * gamma1_0006 - s3_0003 * &
            gamma2_0006)
      if (abs(lambda2_0006) < epsilon_0006 .and. abs(lambda3_0006) < &
            epsilon_0006) then
        mu1_0003 = B1roe_0003 * vroe_0006 * pi1_0003 / w_0003 - &
              B1roe_0003 * vroe_0006 * pi4_0003 / w_0003
        mu2_0003 = B1roe_0003 * B2roe_0003 * pi1_0003 / w_0003 - &
              B1roe_0003 * B2roe_0003 * pi4_0003 / w_0003
        ann_0003 = 0_i64
      else
        mu1_0003 = B1roe_0003 * vroe_0006 * pi1_0003 / w_0003 - &
              B1roe_0003 * vroe_0006 * pi4_0003 / w_0003 - 0.5_f64 * ( &
              gamma1_0006 * pi2_0003 - gamma2_0006 * pi3_0003)
        mu2_0003 = B1roe_0003 * B2roe_0003 * pi1_0003 / w_0003 - &
              B1roe_0003 * B2roe_0003 * pi4_0003 / w_0003 - 0.5_f64 * ( &
              gamma1_0006 * pi2_0003 + gamma2_0006 * pi3_0003)
        ann_0003 = 0_i64
      end if
      !1ère colonne de la matrice A
      signA_0003(0_i64, 0_i64) = (s1_0003 * lambda4_0006 - s4_0003 * &
            lambda1_0006) / (2_i64 * w_0003)
      signA_0003(0_i64, 1_i64) = lambda1_0006 * lambda4_0006 * (s1_0003 &
            - s4_0003) / (2_i64 * w_0003)
      signA_0003(0_i64, 2_i64) = sigma1_0006 / (2_i64 * w_0003)
      signA_0003(0_i64, 3_i64) = 0.0_f64
      signA_0003(0_i64, 4_i64) = sigma2_0006 / (2_i64 * w_0003)
      signA_0003(0_i64, 5_i64) = 0.0_f64
      !2ème colonne de la matrice A
      signA_0003(1_i64, 0_i64) = (s4_0003 - s1_0003) / (2_i64 * w_0003)
      signA_0003(1_i64, 1_i64) = (s4_0003 * lambda4_0006 - s1_0003 * &
            lambda1_0006) / (2_i64 * w_0003)
      signA_0003(1_i64, 2_i64) = vroe_0006 * (s4_0003 - s1_0003) / ( &
            2_i64 * w_0003)
      signA_0003(1_i64, 3_i64) = 0.0_f64
      signA_0003(1_i64, 4_i64) = B2roe_0003 * (s4_0003 - s1_0003) / ( &
            2_i64 * w_0003)
      signA_0003(1_i64, 5_i64) = 0.0_f64
      !3ème colonne de la matrice A
      signA_0003(2_i64, 0_i64) = 0.0_f64
      signA_0003(2_i64, 1_i64) = 0.0_f64
      signA_0003(2_i64, 2_i64) = (s2_0003 + s3_0003) / 2_i64
      signA_0003(2_i64, 3_i64) = 0.0_f64
      signA_0003(2_i64, 4_i64) = (s2_0003 - s3_0003) / 2_i64
      signA_0003(2_i64, 5_i64) = 0.0_f64
      !4ème colonne de la matrice A
      signA_0003(3_i64, 0_i64) = ann_0003 * B1roe_0003 * (pi1_0003 - &
            pi4_0003) / w_0003
      signA_0003(3_i64, 1_i64) = ann_0003 * B1roe_0003 * (s1_0003 - &
            s4_0003) / w_0003
      signA_0003(3_i64, 2_i64) = ann_0003 * mu1_0003
      signA_0003(3_i64, 3_i64) = 0.0_f64
      signA_0003(3_i64, 4_i64) = ann_0003 * mu2_0003
      signA_0003(3_i64, 5_i64) = 0.0_f64
      !5ème colonne de la matrice A
      signA_0003(4_i64, 0_i64) = 0.0_f64
      signA_0003(4_i64, 1_i64) = 0.0_f64
      signA_0003(4_i64, 2_i64) = (s2_0003 - s3_0003) / 2_i64
      signA_0003(4_i64, 3_i64) = 0.0_f64
      signA_0003(4_i64, 4_i64) = (s2_0003 + s3_0003) / 2_i64
      signA_0003(4_i64, 5_i64) = 0.0_f64
      !6ème colonne de la matrice A
      signA_0003(5_i64, 0_i64) = sound_0006 ** 2_i64 * (pi4_0003 - &
            pi1_0003) / (2_i64 * w_0003)
      signA_0003(5_i64, 1_i64) = sound_0006 ** 2_i64 * (s4_0003 - &
            s1_0003) / (2_i64 * w_0003)
      signA_0003(5_i64, 2_i64) = sound_0006 ** 2_i64 * vroe_0006 * ( &
            pi4_0003 - pi1_0003) / (2_i64 * w_0003)
      signA_0003(5_i64, 3_i64) = 0.0_f64
      signA_0003(5_i64, 4_i64) = sound_0006 ** 2_i64 * B2roe_0003 * ( &
            pi4_0003 - pi1_0003) / (2_i64 * w_0003)
      signA_0003(5_i64, 5_i64) = 0.0_f64
      smmat_0003(0:, 0:) => signA_0003
      hnew_0006 = 0.0_f64
      unew_0006 = 0.0_f64
      vnew_0006 = 0.0_f64
      B1new_0003 = 0.0_f64
      B2new_0003 = 0.0_f64
      znew_0006 = 0.0_f64
      do i_0003 = 0_i64, 5_i64, 1_i64
        hnew_0006 = hnew_0006 + smmat_0003(i_0003, 0_i64) * w_dif_0006( &
              i_0003)
        unew_0006 = unew_0006 + smmat_0003(i_0003, 1_i64) * w_dif_0006( &
              i_0003)
        vnew_0006 = vnew_0006 + smmat_0003(i_0003, 2_i64) * w_dif_0006( &
              i_0003)
        B1new_0003 = B1new_0003 + smmat_0003(i_0003, 3_i64) * w_dif_0006 &
              (i_0003)
        B2new_0003 = B2new_0003 + smmat_0003(i_0003, 4_i64) * w_dif_0006 &
              (i_0003)
        znew_0006 = znew_0006 + smmat_0003(i_0003, 5_i64) * w_dif_0006( &
              i_0003)
      end do
      Pnew_0003 = cpsi * (B1right_0003 - B1left_0003)
      u_h_0006 = hnew_0006 / 2_i64
      u_hu_0006 = unew_0006 / 2_i64
      u_hv_0006 = vnew_0006 / 2_i64
      u_hP_0003 = Pnew_0003 / 2_i64
      u_hB1_0003 = B1new_0003 / 2_i64
      u_hB2_0003 = B2new_0003 / 2_i64
      u_z_0006 = znew_0006 / 2_i64
      w_lrh_0006 = w_lrh_0006 - u_h_0006
      w_lrhu_0006 = w_lrhu_0006 - u_hu_0006
      w_lrhv_0006 = w_lrhv_0006 - u_hv_0006
      w_lrhP_0003 = w_lrhPSI_0003 - u_hP_0003
      w_lrhB1_0003 = absy_0003
      w_lrhB2_0003 = w_lrhB2_0003 - u_hB2_0003
      w_lrz_0006 = w_lrz_0006 - u_z_0006
      w_hP_0003 = hPSI_r - hPSI_l
      mw_hB1_0003 = w_lrhB1_0003 / 2_i64 - 1_i64 / (2_i64 * cpsi) * &
            w_hP_0003
      mhP_0003 = (hPSI_r - hPSI_l) / 2_i64 - cpsi * w_dif_0006(3_i64) / &
            2_i64
      unew_0006 = 0.0_f64
      vnew_0006 = 0.0_f64
      B1new_0003 = 0.0_f64
      B2new_0003 = 0.0_f64
      unew_0006 = w_lrhu_0006 * normal(0_i64) - w_lrhv_0006 * normal( &
            1_i64)
      unew_0006 = unew_0006 / mesure
      vnew_0006 = w_lrhu_0006 * normal(1_i64) + w_lrhv_0006 * normal( &
            0_i64)
      vnew_0006 = vnew_0006 / mesure
      B1new_0003 = w_lrhB1_0003 * normal(0_i64) - w_lrhB2_0003 * normal( &
            1_i64)
      B1new_0003 = B1new_0003 / mesure
      B2new_0003 = w_lrhB1_0003 * normal(1_i64) + w_lrhB2_0003 * normal( &
            0_i64)
      B2new_0003 = B2new_0003 / mesure
      w_lrhu_0006 = unew_0006
      w_lrhv_0006 = vnew_0006
      w_lrhB1_0003 = B1new_0003
      w_lrhB2_0003 = B2new_0003
      allocate(norm_0003(0:size(normal, kind=i64) - 1_i64))
      norm_0003 = normal / mesure
      q_s_0006 = normal(0_i64) * unew_0006 + normal(1_i64) * vnew_0006
      p_s_0003 = normal(0_i64) * B1new_0003 + normal(1_i64) * B2new_0003
      Flux_B1psi_0003 = mhP_0003 * norm_0003(0_i64) * mesure
      Flux_B2psi_0003 = mhP_0003 * norm_0003(1_i64) * mesure
      Flux_hPpsi_0003 = cpsi * cpsi * mw_hB1_0003 * mesure
      flux(0_i64) = q_s_0006
      flux(1_i64) = q_s_0006 * w_lrhu_0006 / w_lrh_0006 + 0.5_f64 * grav &
            * w_lrh_0006 * w_lrh_0006 * normal(0_i64) - p_s_0003 * &
            w_lrhB1_0003 / w_lrh_0006
      flux(2_i64) = q_s_0006 * w_lrhv_0006 / w_lrh_0006 + 0.5_f64 * grav &
            * w_lrh_0006 * w_lrh_0006 * normal(1_i64) - p_s_0003 * &
            w_lrhB2_0003 / w_lrh_0006
      flux(3_i64) = (w_lrhv_0006 * w_lrhB1_0003 / w_lrh_0006 - &
            w_lrhu_0006 * w_lrhB2_0003 / w_lrh_0006) * normal(1_i64) + &
            Flux_B1psi_0003
      flux(4_i64) = (w_lrhu_0006 * w_lrhB2_0003 / w_lrh_0006 - &
            w_lrhv_0006 * w_lrhB1_0003 / w_lrh_0006) * normal(0_i64) + &
            Flux_B2psi_0003
      flux(5_i64) = Flux_hPpsi_0003
      flux(6_i64) = 0_i64
      if (allocated(ninv_0006)) then
        deallocate(ninv_0006)
      end if
      if (allocated(w_dif_0006)) then
        deallocate(w_dif_0006)
      end if
      if (allocated(signA_0003)) then
        deallocate(signA_0003)
      end if
      if (allocated(norm_0003)) then
        deallocate(norm_0003)
      end if
      rez_h(cellidf(0_i64, i)) = rez_h(cellidf(0_i64, i)) - flux(0_i64)
      rez_hu(cellidf(0_i64, i)) = rez_hu(cellidf(0_i64, i)) - flux(1_i64 &
            )
      rez_hv(cellidf(0_i64, i)) = rez_hv(cellidf(0_i64, i)) - flux(2_i64 &
            )
      rez_hB1(cellidf(0_i64, i)) = rez_hB1(cellidf(0_i64, i)) - flux( &
            3_i64)
      rez_hB2(cellidf(0_i64, i)) = rez_hB2(cellidf(0_i64, i)) - flux( &
            4_i64)
      rez_PSI(cellidf(0_i64, i)) = rez_PSI(cellidf(0_i64, i)) - flux( &
            5_i64)
      rez_Z(cellidf(0_i64, i)) = rez_Z(cellidf(0_i64, i)) - flux(6_i64)
    end do
    if (allocated(flux)) then
      deallocate(flux)
    end if
    if (allocated(r_l)) then
      deallocate(r_l)
    end if
    if (allocated(r_r)) then
      deallocate(r_r)
    end if

  end subroutine explicitscheme_convective_SWMHD
  !........................................

  !........................................
  subroutine term_coriolis_SW(hu_c, hv_c, corio_hu, corio_hv, f_c) 

    implicit none

    real(f64), intent(in) :: hu_c(0:)
    real(f64), intent(in) :: hv_c(0:)
    real(f64), intent(inout) :: corio_hu(0:)
    real(f64), intent(inout) :: corio_hv(0:)
    real(f64), value :: f_c
    integer(i64) :: i

    !$omp parallel do
    do i = 0_i64, size(hu_c, kind=i64) - 1_i64, 1_i64
      corio_hu(i) = f_c * hu_c(i)
      corio_hv(i) = (-f_c) * hv_c(i)
    end do
    !$omp end parallel do

  end subroutine term_coriolis_SW
  !........................................

  !........................................
  subroutine term_friction_SW(h_c, hu_c, hv_c, grav, eta, time) 

    implicit none

    real(f64), intent(in) :: h_c(0:)
    real(f64), intent(inout) :: hu_c(0:)
    real(f64), intent(inout) :: hv_c(0:)
    real(f64), value :: grav
    real(f64), value :: eta
    real(f64), value :: time
    integer(i64) :: nbelement
    integer(i64) :: i
    real(f64) :: ufric
    real(f64) :: vfric
    real(f64) :: hfric
    real(f64) :: A
    real(f64) :: hutild
    real(f64) :: hvtild

    nbelement = size(h_c, kind=i64)
    !$omp parallel do
    do i = 0_i64, nbelement - 1_i64, 1_i64
      ufric = hu_c(i) / h_c(i)
      vfric = hv_c(i) / h_c(i)
      hfric = h_c(i)
      A = 1_i64 + grav * time * eta ** 2_i64 * sqrt(ufric ** 2_i64 + &
            vfric ** 2_i64) / hfric ** (4.0_f64 / 3.0_f64)
      hutild = hu_c(i) / A
      hvtild = hv_c(i) / A
      hu_c(i) = hutild
      hv_c(i) = hvtild
    end do
    !$omp end parallel do

  end subroutine term_friction_SW
  !........................................

  !........................................
  subroutine term_wind_SW(uwind, vwind, TAUXWX, TAUXWY) 

    implicit none

    real(f64), intent(out) :: TAUXWX
    real(f64), intent(out) :: TAUXWY
    real(f64), value :: uwind
    real(f64), value :: vwind
    real(f64) :: RHO_air
    real(f64) :: WNORME
    real(f64) :: CW

    RHO_air = 1.28_f64
    WNORME = sqrt(uwind ** 2_i64 + vwind ** 2_i64)
    CW = RHO_air * (0.75_f64 + 0.067_f64 * WNORME) * 0.001_f64
    TAUXWX = CW * uwind * WNORME
    TAUXWY = CW * vwind * WNORME
    return

  end subroutine term_wind_SW
  !........................................

end module tools
