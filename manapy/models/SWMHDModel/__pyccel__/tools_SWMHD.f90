module tools_SWMHD


  use, intrinsic :: ISO_C_Binding, only : f64 => C_DOUBLE , i64 => &
        C_INT64_T
  implicit none

  contains

  !........................................
  subroutine Total_Energy(h_c, hu_c, hv_c, hB1_c, hB2_c, Z_c, grav, &
        volumec, numt, numc, numm, nump)

    implicit none

    real(f64), intent(out) :: numt
    real(f64), intent(out) :: numc
    real(f64), intent(out) :: numm
    real(f64), intent(out) :: nump
    real(f64), intent(in) :: h_c(0:)
    real(f64), intent(in) :: hu_c(0:)
    real(f64), intent(in) :: hv_c(0:)
    real(f64), intent(in) :: hB1_c(0:)
    real(f64), intent(in) :: hB2_c(0:)
    real(f64), intent(in) :: Z_c(0:)
    real(f64), value :: grav
    real(f64), intent(in) :: volumec(0:)
    integer(i64) :: nbelement
    real(f64) :: num_t
    real(f64) :: num_c
    real(f64) :: num_m
    real(f64) :: num_p
    real(f64), allocatable :: Et(:)
    real(f64), allocatable :: Ec(:)
    real(f64), allocatable :: Ep(:)
    real(f64), allocatable :: Em(:)
    integer(i64) :: i
    real(f64) :: hc
    real(f64) :: uc
    real(f64) :: vc
    real(f64) :: B1c
    real(f64) :: B2c
    real(f64) :: bc

    nbelement = size(h_c, kind=i64)
    num_t = 0.0_f64
    num_c = 0.0_f64
    num_m = 0.0_f64
    num_p = 0.0_f64
    !from numpy import zeros
    allocate(Et(0:nbelement - 1_i64))
    Et = 0.0_f64
    allocate(Ec(0:nbelement - 1_i64))
    Ec = 0.0_f64
    allocate(Ep(0:nbelement - 1_i64))
    Ep = 0.0_f64
    allocate(Em(0:nbelement - 1_i64))
    Em = 0.0_f64
    do i = 0_i64, nbelement - 1_i64, 1_i64
      hc = h_c(i)
      uc = hu_c(i) / h_c(i)
      vc = hv_c(i) / h_c(i)
      B1c = hB1_c(i) / h_c(i)
      B2c = hB2_c(i) / h_c(i)
      bc = Z_c(i)
      Ec(i) = 0.5_f64 * hc * (uc ** 2_i64 + vc ** 2_i64)
      Em(i) = 0.5_f64 * hc * (B1c ** 2_i64 + B2c ** 2_i64)
      Ep(i) = 0.5_f64 * grav * hc ** 2_i64 + grav * bc * hc
      Et(i) = Ec(i) + Em(i) + Ep(i)
      num_t = num_t + volumec(i) * Et(i)
      num_c = num_c + volumec(i) * Ec(i)
      num_m = num_m + volumec(i) * Em(i)
      num_p = num_p + volumec(i) * Ep(i)
    end do
    numt = num_t
    numc = num_c
    numm = num_m
    nump = num_p
    if (allocated(Et)) then
      deallocate(Et)
    end if
    if (allocated(Ec)) then
      deallocate(Ec)
    end if
    if (allocated(Ep)) then
      deallocate(Ep)
    end if
    if (allocated(Em)) then
      deallocate(Em)
    end if
    return

  end subroutine Total_Energy
  !........................................

  !........................................
  subroutine explicitscheme_stabilizater(wx_face, wy_face, cellidf, &
        normalf, namef, vepsilon, dissip_w, h_c)

    implicit none

    real(f64), intent(in) :: wx_face(0:)
    real(f64), intent(in) :: wy_face(0:)
    integer(i64), intent(inout) :: cellidf(0:,0:)
    real(f64), intent(in) :: normalf(0:,0:)
    integer(i64), intent(in) :: namef(0:)
    real(f64), value :: vepsilon
    real(f64), intent(inout) :: dissip_w(0:)
    real(f64), intent(in) :: h_c(0:)
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
      q = vepsilon * h_c(i) * (wx_face(i) * norm(0_i64) + wy_face(i) * &
            norm(1_i64))
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

  end subroutine explicitscheme_stabilizater
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

    grav = 1.0_f64
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
  subroutine initialisation_SWMHD(h, hu, hv, hB1, hB2, PSI, Z, center, &
        choix, k1, k2, eps, tol)

    implicit none

    real(f64), intent(inout) :: h(0:)
    real(f64), intent(inout) :: hu(0:)
    real(f64), intent(inout) :: hv(0:)
    real(f64), intent(inout) :: hB1(0:)
    real(f64), intent(inout) :: hB2(0:)
    real(f64), intent(inout) :: PSI(0:)
    real(f64), intent(inout) :: Z(0:)
    real(f64), intent(in) :: center(0:,0:)
    integer(i64), value :: choix
    real(f64), value :: k1
    real(f64), value :: k2
    real(f64), value :: eps
    real(f64), value :: tol
    integer(i64) :: nbelements
    integer(i64) :: i
    real(f64) :: xcent
    real(f64) :: ycent
    real(f64) :: AA
    real(f64) :: uu
    real(f64) :: bb
    real(f64) :: h0
    real(f64) :: u0
    real(f64) :: v0
    real(f64) :: B01
    real(f64) :: B02
    real(f64) :: grav
    real(f64) :: k
    real(f64) :: c0
    real(f64) :: ss
    real(f64) :: a1
    real(f64) :: b1
    real(f64) :: h_hot
    real(f64) :: u_hot
    real(f64) :: v_hot
    real(f64) :: B1_hot
    real(f64) :: B2_hot
    real(f64) :: hint
    real(f64) :: uint
    real(f64) :: vint
    real(f64) :: B1int
    real(f64) :: B2int
    real(f64) :: gamma
    real(f64) :: g
    real(f64) :: umax
    real(f64) :: Bmax
    real(f64) :: hmax
    real(f64) :: rcent
    real(f64) :: ee
    real(f64) :: e1
    real(f64) :: hin
    real(f64) :: uin
    real(f64) :: vin
    real(f64) :: B1in
    real(f64) :: B2in
    real(f64) :: s1
    real(f64) :: s2
    real(f64) :: s3

    !fabs, pi, cos, sin, exp
    nbelements = size(center,2,i64)
    if (choix == -1_i64) then
      do i = 0_i64, nbelements - 1_i64, 1_i64
        xcent = center(0_i64, i)
        ycent = center(1_i64, i)
        AA = sin(2_i64 * 3.141592653589793_f64 * xcent - 2_i64 * &
              3.141592653589793_f64 * ycent)
        uu = 2.0_f64
        bb = 5.0_f64
        h0 = 100.0_f64
        u0 = (uu + AA) / h0
        v0 = (uu + AA) / h0
        B01 = (bb + 2_i64 * AA) / h0
        B02 = (bb + 2_i64 * AA) / h0
        PSI(i) = 0.0_f64
        h(i) = 100_i64
        hu(i) = uu + AA
        hv(i) = uu + AA
        hB1(i) = bb + 2_i64 * AA
        hB2(i) = bb + 2_i64 * AA
        PSI(i) = 0.0_f64
        grav = 1.0_f64
        k = sqrt(k1 ** 2_i64 + k2 ** 2_i64)
        !B0 = sqrt(B01**2 + B02**2)
        c0 = sqrt(h0 * grav)
        !Produit scalaire de (k1,k2) et (B01, B02)
        ss = k1 * B01 + k2 * B02
        a1 = sqrt(ss ** 2_i64 + (k * c0) ** 2_i64)
        b1 = (-0.5_f64) * eps * k * k
        h_hot = h0 * k * k * cos(k1 * xcent + k2 * ycent)
        u_hot = k1 * (a1 * cos(k1 * xcent + k2 * ycent) - b1 * sin(k1 * &
              xcent + k2 * ycent))
        v_hot = k2 * (a1 * cos(k1 * xcent + k2 * ycent) - b1 * sin(k1 * &
              xcent + k2 * ycent))
        B1_hot = (-k1) * ss * cos(k1 * xcent + k2 * ycent)
        B2_hot = (-k2) * ss * cos(k1 * xcent + k2 * ycent)
        hint = h0 + tol * h_hot
        uint = u0 + tol * u_hot
        vint = v0 + tol * v_hot
        B1int = B01 + tol * B1_hot
        B2int = B02 + tol * B2_hot
        h(i) = hint
        hu(i) = hint * uint
        hv(i) = hint * vint
        hB1(i) = hint * B1int
        hB2(i) = hint * B2int
      end do
    end if
    if (choix == -3_i64) then
      do i = 0_i64, nbelements - 1_i64, 1_i64
        xcent = center(0_i64, i)
        ycent = center(1_i64, i)
        AA = 0.0_f64 * sin(2_i64 * 3.141592653589793_f64 * xcent - 2_i64 &
              * 3.141592653589793_f64 * ycent)
        uu = 0.0_f64
        bb = 3.0_f64
        h0 = 2.0_f64
        u0 = (uu + AA) / h0
        v0 = (uu + AA) / h0
        B01 = (bb + 2_i64 * AA) / h0
        B02 = (bb + 2_i64 * AA) / h0
        PSI(i) = 0.0_f64
        h(i) = 100_i64
        hu(i) = uu + AA
        hv(i) = uu + AA
        hB1(i) = bb + 2_i64 * AA
        hB2(i) = bb + 2_i64 * AA
        PSI(i) = 0.0_f64
        grav = 1.0_f64
        k = sqrt(k1 ** 2_i64 + k2 ** 2_i64)
        !B0 = sqrt(B01**2 + B02**2)
        c0 = sqrt(h0 * grav)
        !Produit scalaire de (k1,k2) et (B01, B02)
        ss = k1 * B01 + k2 * B02
        a1 = sqrt(ss ** 2_i64 + (k * c0) ** 2_i64)
        b1 = (-0.5_f64) * eps * k * k
        h_hot = h0 * k * k * cos(k1 * xcent + k2 * ycent)
        u_hot = k1 * (a1 * cos(k1 * xcent + k2 * ycent) - b1 * sin(k1 * &
              xcent + k2 * ycent))
        v_hot = k2 * (a1 * cos(k1 * xcent + k2 * ycent) - b1 * sin(k1 * &
              xcent + k2 * ycent))
        B1_hot = (-k1) * ss * cos(k1 * xcent + k2 * ycent)
        B2_hot = (-k2) * ss * cos(k1 * xcent + k2 * ycent)
        hint = h0 + tol * h_hot
        uint = u0 + tol * u_hot
        vint = v0 + tol * v_hot
        B1int = B01 + tol * B1_hot
        B2int = B02 + tol * B2_hot
        h(i) = hint
        hu(i) = hint * uint
        hv(i) = hint * vint
        hB1(i) = hint * B1int
        hB2(i) = hint * B2int
      end do
    end if
    if (choix == -2_i64) then
      do i = 0_i64, nbelements - 1_i64, 1_i64
        xcent = center(0_i64, i)
        ycent = center(1_i64, i)
        AA = sin(2_i64 * 3.141592653589793_f64 * xcent - 2_i64 * &
              3.141592653589793_f64 * ycent)
        uu = 2.0_f64
        bb = 5.0_f64
        h0 = 100.0_f64
        u0 = (uu + AA) / h0
        v0 = (uu + AA) / h0
        B01 = (bb + 2_i64 * AA) / h0
        B02 = (bb + 2_i64 * AA) / h0
        PSI(i) = 0.0_f64
        h(i) = 100_i64
        hu(i) = uu + AA
        hv(i) = uu + AA
        hB1(i) = bb + 2_i64 * AA
        hB2(i) = bb + 2_i64 * AA
        PSI(i) = 0.0_f64
        grav = 1.0_f64
        k = sqrt(k1 ** 2_i64 + k2 ** 2_i64)
        !B0 = sqrt(B01**2 + B02**2)
        c0 = sqrt(h0 * grav)
        !Produit scalaire de (k1,k2) et (B01, B02)
        ss = k1 * B01 + k2 * B02
        b1 = (-0.5_f64) * eps * k * k
        h_hot = 0.0_f64
        u_hot = k2 * ss * cos(k1 * xcent + k2 * ycent)
        v_hot = (-k1) * ss * cos(k1 * xcent + k2 * ycent)
        B1_hot = (-k2) * (ss * cos(k1 * xcent + k2 * ycent) + b1 * sin( &
              k1 * xcent + k2 * ycent))
        B2_hot = (-k1) * (ss * cos(k1 * xcent + k2 * ycent) + b1 * sin( &
              k1 * xcent + k2 * ycent))
        hint = h0 + tol * h_hot
        uint = u0 + tol * u_hot
        vint = v0 + tol * v_hot
        B1int = B01 + tol * B1_hot
        B2int = B02 + tol * B2_hot
        h(i) = hint
        hu(i) = hint * uint
        hv(i) = hint * vint
        hB1(i) = hint * B1int
        hB2(i) = hint * B2int
      end do
    else if (choix == 1_i64) then
      do i = 0_i64, nbelements - 1_i64, 1_i64
        xcent = center(0_i64, i)
        if (xcent <= 0_i64) then
          hu(i) = 0.0_f64
          h(i) = 1.0_f64
          hv(i) = 0.0_f64
          Z(i) = 0.0_f64
          hB1(i) = 1.0_f64
          hB2(i) = 0.0_f64
          PSI(i) = 0.0_f64
        else
          Z(i) = 0.0_f64
          h(i) = 2.0_f64
          hu(i) = 0.0_f64
          hv(i) = 0.0_f64
          hB1(i) = 1.0_f64
          hB2(i) = 2.0_f64
          PSI(i) = 0.0_f64
        end if
      end do
    else if (choix == 2_i64) then
      do i = 0_i64, nbelements - 1_i64, 1_i64
        xcent = center(0_i64, i)
        ycent = center(1_i64, i)
        if (sqrt(xcent ** 2_i64 + ycent ** 2_i64) <= 0.3_f64) then
          hu(i) = 0.0_f64
          h(i) = 1.0_f64
          hv(i) = 0.0_f64
          Z(i) = 0.0_f64
          hB1(i) = 0.1_f64
          hB2(i) = 0.0_f64
          PSI(i) = 0.0_f64
        else
          Z(i) = 0.0_f64
          h(i) = 0.1_f64
          hu(i) = 0.0_f64
          hv(i) = 0.0_f64
          hB1(i) = 0.1_f64
          hB2(i) = 0.0_f64
          PSI(i) = 0.0_f64
        end if
      end do
    else if (choix == 12_i64) then
      do i = 0_i64, nbelements - 1_i64, 1_i64
        xcent = center(0_i64, i)
        ycent = center(1_i64, i)
        gamma = 1.667_f64
        hint = gamma ** 2_i64
        h(i) = hint
        hu(i) = (-hint) * sin(ycent)
        hv(i) = hint * sin(xcent)
        hB1(i) = (-hint) * sin(ycent)
        hB2(i) = (-hint) * sin(2_i64 * xcent)
        PSI(i) = 0.0_f64
      end do
    else if (choix == 40_i64) then
      do i = 0_i64, nbelements - 1_i64, 1_i64
        xcent = center(0_i64, i)
        ycent = center(1_i64, i)
        g = 1.0_f64
        umax = 0.2_f64
        Bmax = 0.1_f64
        hmax = 1.0_f64
        rcent = sqrt(xcent ** 2_i64 + ycent ** 2_i64)
        ee = exp(1_i64 - rcent ** 2_i64)
        e1 = exp(0.5_f64 * (1_i64 - rcent ** 2_i64))
        hin = hmax - 1.0_f64 / (2.0_f64 * g) * (umax ** 2_i64 - Bmax ** &
              2_i64) * ee
        uin = 1.0_f64 - umax * e1 * ycent
        vin = 1.0_f64 + umax * e1 * xcent
        B1in = (-Bmax) * e1 * ycent
        B2in = Bmax * e1 * xcent
        h(i) = hin
        hu(i) = hin * uin
        hv(i) = hin * vin
        hB1(i) = hin * B1in
        hB2(i) = hin * B2in
        PSI(i) = 0.0_f64
      end do
    else if (choix == 60_i64) then
      do i = 0_i64, nbelements - 1_i64, 1_i64
        xcent = center(0_i64, i)
        ycent = center(1_i64, i)
        AA = sin(2_i64 * 3.141592653589793_f64 * xcent - 2_i64 * &
              3.141592653589793_f64 * ycent)
        uu = 2.0_f64
        bb = 5.0_f64
        h(i) = 100_i64
        hu(i) = uu + AA
        hv(i) = uu + AA
        hB1(i) = bb + 2_i64 * AA
        hB2(i) = bb + 2_i64 * AA
        PSI(i) = 0.0_f64
      end do
    else if (choix == 50_i64) then
      do i = 0_i64, nbelements - 1_i64, 1_i64
        xcent = center(0_i64, i)
        ycent = center(1_i64, i)
        Z(i) = 0.8_f64 * exp((-5_i64) * (xcent - 1_i64) ** 2_i64 - &
              50_i64 * (ycent - 0.5_f64) ** 2_i64)
        if (0.05_f64 < xcent .and. xcent < 0.15_f64) then
          h(i) = 1.0_f64 - Z(i) + 0.01_f64
        else
          h(i) = 1_i64 - Z(i)
        end if
        hu(i) = 0.0_f64
        hv(i) = 0.0_f64
        hB1(i) = 1.002_f64 * h(i)
        hB2(i) = 0.1_f64 * h(i)
        PSI(i) = 0.0_f64
      end do
    else if (choix == 51_i64) then
      do i = 0_i64, nbelements - 1_i64, 1_i64
        xcent = center(0_i64, i)
        ycent = center(1_i64, i)
        Z(i) = 0.8_f64 * exp((-5_i64) * (xcent - 1_i64) ** 2_i64 - &
              50_i64 * (ycent - 0.5_f64) ** 2_i64)
        if (0.05_f64 < xcent .and. xcent < 0.15_f64) then
          h(i) = 1.0_f64 - Z(i)
        else
          h(i) = 1_i64 - Z(i)
        end if
        hu(i) = 0.0_f64
        hv(i) = 0.0_f64
        hB1(i) = 1.0_f64
        hB2(i) = 0.0_f64
        PSI(i) = 0.0_f64
      end do
    else if (choix == 52_i64) then
      do i = 0_i64, nbelements - 1_i64, 1_i64
        xcent = center(0_i64, i)
        ycent = center(1_i64, i)
        Z(i) = 0.5_f64 * cos(10_i64 * 3.141592653589793_f64 * (xcent - &
              1.5_f64) + 1.0_f64)
        if (0.05_f64 < xcent .and. xcent < 0.15_f64) then
          h(i) = 1.0_f64 - Z(i)
        else
          h(i) = 1_i64 - Z(i)
        end if
        hu(i) = 0.0_f64
        hv(i) = 0.0_f64
        hB1(i) = 1.0_f64
        hB2(i) = 0.0_f64
        PSI(i) = 0.0_f64
      end do
    else if (choix == 53_i64) then
      do i = 0_i64, nbelements - 1_i64, 1_i64
        xcent = center(0_i64, i)
        ycent = center(1_i64, i)
        Z(i) = 0.2_f64 * exp((-0.5_f64) * (xcent + 1_i64) ** 2_i64) + &
              0.3_f64 * exp(-(xcent - 1.5_f64) ** 2_i64)
        if (0.05_f64 < xcent .and. xcent < 0.15_f64) then
          h(i) = 1.0_f64
        else
          h(i) = 1.0_f64
        end if
        hu(i) = 1.0_f64 * h(i)
        hv(i) = 0.0_f64
        if (xcent <= 0_i64) then
          hB1(i) = 0.05_f64 * h(i) * 0.0_f64
          hB2(i) = 0.0_f64
          PSI(i) = 0.0_f64
        else
          hB1(i) = 0.1_f64 * h(i) * 0.0_f64
          hB2(i) = 0.1_f64 * h(i) * 0.0_f64
          PSI(i) = 0.0_f64
        end if
      end do
    else if (choix == 6_i64) then
      do i = 0_i64, nbelements - 1_i64, 1_i64
        xcent = center(0_i64, i)
        ycent = center(1_i64, i)
        g = 9.81_f64
        s1 = sin(2_i64 * 3.141592653589793_f64 * xcent)
        s2 = s1 + 2_i64
        s3 = 2_i64 * g
        Z(i) = s1 - 1_i64 / (s3 * s2 ** 2_i64) + 2_i64
        h(i) = 2_i64 + s1
        hu(i) = 2_i64 + s1
        hv(i) = 2_i64 + s1
        hB1(i) = 1_i64
        hB2(i) = 4_i64 + 2_i64 * s1
        PSI(i) = 0.0_f64
      end do
    else if (choix == 7_i64) then
      do i = 0_i64, nbelements - 1_i64, 1_i64
        xcent = center(0_i64, i)
        ycent = center(1_i64, i)
        if (xcent <= 0_i64) then
          hu(i) = 0.0_f64
          h(i) = 1.0_f64
          hv(i) = 1.0_f64
          Z(i) = 0.0_f64
          hB1(i) = 1.0_f64
          hB2(i) = 1.0_f64
          PSI(i) = 0.0_f64
        else
          Z(i) = 0.0_f64
          h(i) = 1.0_f64 - 0.0001_f64
          hu(i) = 0.0_f64
          hv(i) = 1.0_f64 - 0.0001_f64
          hB1(i) = (1.0_f64 - 0.0001_f64) * (1_i64 + 0.0001_f64 / ( &
                1.0_f64 - 0.0001_f64))
          hB2(i) = (1_i64 + 0.0002_f64) * (1_i64 - 0.0001_f64)
          PSI(i) = 0.0_f64
        end if
      end do
    end if

  end subroutine initialisation_SWMHD
  !........................................

  !........................................
  subroutine term_source_srnh_SWMHD(src_h, src_hu, src_hv, src_hB1, &
        src_hB2, src_PSI, src_Z, h_c, hu_c, hv_c, Z_c, h_ghost, &
        hu_ghost, hv_ghost, Z_ghost, h_halo, hu_halo, hv_halo, Z_halo, &
        h_x, h_y, psi, hx_halo, hy_halo, psi_halo, nodeidc, faceidc, &
        cellidc, cellidf, centerc, normalc, namef, centerf, centerh, &
        vertexn, halofid, order)

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
    real(f64), intent(in) :: Z_c(0:)
    real(f64), intent(in) :: h_ghost(0:)
    real(f64), intent(in) :: hu_ghost(0:)
    real(f64), intent(in) :: hv_ghost(0:)
    real(f64), intent(in) :: Z_ghost(0:)
    real(f64), intent(in) :: h_halo(0:)
    real(f64), intent(in) :: hu_halo(0:)
    real(f64), intent(in) :: hv_halo(0:)
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
  subroutine explicitscheme_convective_SWMHD(rez_h, rez_hu, rez_hv, &
        rez_hB1, rez_hB2, rez_PSI, rez_Z, h_c, hu_c, hv_c, hB1_c, hB2_c &
        , hPSIc, Z_c, h_ghost, hu_ghost, hv_ghost, hB1_ghost, hB2_ghost &
        , hPSIghost, Z_ghost, h_halo, hu_halo, hv_halo, hB1_halo, &
        hB2_halo, hPSIhalo, Z_halo, h_x, h_y, hx_halo, hy_halo, psi, &
        psi_halo, centerc, centerf, centerh, centerg, cellidf, mesuref, &
        normalf, halofid, innerfaces, halofaces, boundaryfaces, &
        periodicboundaryfaces, shift, order, cpsi, hB1_cst, hB2_cst)

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
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: centerf(0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    real(f64), intent(in) :: centerg(0:,0:)
    integer(i64), intent(inout) :: cellidf(0:,0:)
    real(f64), intent(in) :: mesuref(0:)
    real(f64), intent(in), target :: normalf(0:,0:)
    integer(i64), intent(in) :: halofid(0:)
    integer(i64), intent(in) :: innerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: boundaryfaces(0:)
    integer(i64), intent(in) :: periodicboundaryfaces(0:)
    real(f64), intent(in) :: shift(0:,0:)
    integer(i64), value :: order
    real(f64), value :: cpsi
    real(f64), intent(in) :: hB1_cst(0:)
    real(f64), intent(in) :: hB2_cst(0:)
    real(f64) :: grav
    real(f64), allocatable :: flux(:)
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
    integer(i64) :: Dummy_0005
    real(f64), allocatable :: ninv_0001(:)
    real(f64), allocatable :: w_dif_0001(:)
    real(f64) :: u_h_0001
    real(f64) :: v_h_0001
    real(f64) :: B1_h_0001
    real(f64) :: B2_h_0001
    real(f64) :: un_h_0001
    real(f64) :: vn_h_0001
    real(f64) :: B1n_h_0001
    real(f64) :: B2n_h_0001
    real(f64) :: hroe_0001
    real(f64) :: uroe_0001
    real(f64) :: vroe_0001
    real(f64) :: B1roe_0001
    real(f64) :: B2roe_0001
    real(f64) :: uleft_0001
    real(f64) :: vleft_0001
    real(f64) :: B1left_0001
    real(f64) :: B2left_0001
    real(f64) :: uright_0001
    real(f64) :: vright_0001
    real(f64) :: B1right_0001
    real(f64) :: B2right_0001
    real(f64) :: B1i_0001
    real(f64) :: B1j_0001
    real(f64) :: absy_0001
    real(f64) :: w_lrh_0001
    real(f64) :: w_lrhu_0001
    real(f64) :: w_lrhv_0001
    real(f64) :: w_lrhB1_0001
    real(f64) :: w_lrhB2_0001
    real(f64) :: w_lrhPSI_0001
    real(f64) :: w_lrz_0001
    real(f64), allocatable, target :: signA_0001(:,:)
    real(f64) :: sound_0001
    real(f64) :: w_0001
    real(f64) :: lambda1_0001
    real(f64) :: lambda2_0001
    real(f64) :: lambda3_0001
    real(f64) :: lambda4_0001
    real(f64) :: epsilon_0001
    real(f64) :: s1_0001
    real(f64) :: pi1_0001
    real(f64) :: s2_0001
    real(f64) :: pi2_0001
    real(f64) :: s3_0001
    real(f64) :: pi3_0001
    real(f64) :: s4_0001
    real(f64) :: pi4_0001
    real(f64) :: gamma1_0001
    real(f64) :: gamma2_0001
    real(f64) :: sigma1_0001
    real(f64) :: sigma2_0001
    real(f64) :: mu1_0001
    real(f64) :: mu2_0001
    integer(i64) :: ann_0001
    real(f64), pointer :: smmat_0001(:,:)
    real(f64) :: hnew_0001
    real(f64) :: unew_0001
    real(f64) :: vnew_0001
    real(f64) :: B1new_0001
    real(f64) :: B2new_0001
    real(f64) :: znew_0001
    real(f64) :: Pnew_0001
    real(f64) :: u_hu_0001
    real(f64) :: u_hv_0001
    real(f64) :: u_hP_0001
    real(f64) :: u_hB1_0001
    real(f64) :: u_hB2_0001
    real(f64) :: u_z_0001
    real(f64) :: w_lrhP_0001
    real(f64) :: w_hP_0001
    real(f64) :: mw_hB1_0001
    real(f64) :: mhP_0001
    real(f64), allocatable :: norm_0001(:)
    real(f64) :: q_s_0001
    real(f64) :: p_s_0001
    real(f64) :: Flux_B1psi_0001
    real(f64) :: Flux_B2psi_0001
    real(f64) :: Flux_hPpsi_0001
    integer(i64) :: i_0001
    integer(i64) :: Dummy_0006
    real(f64), allocatable :: ninv_0002(:)
    real(f64), allocatable :: w_dif_0002(:)
    real(f64) :: u_h_0002
    real(f64) :: v_h_0002
    real(f64) :: B1_h_0002
    real(f64) :: B2_h_0002
    real(f64) :: un_h_0002
    real(f64) :: vn_h_0002
    real(f64) :: B1n_h_0002
    real(f64) :: B2n_h_0002
    real(f64) :: hroe_0002
    real(f64) :: uroe_0002
    real(f64) :: vroe_0002
    real(f64) :: B1roe_0002
    real(f64) :: B2roe_0002
    real(f64) :: uleft_0002
    real(f64) :: vleft_0002
    real(f64) :: B1left_0002
    real(f64) :: B2left_0002
    real(f64) :: uright_0002
    real(f64) :: vright_0002
    real(f64) :: B1right_0002
    real(f64) :: B2right_0002
    real(f64) :: B1i_0002
    real(f64) :: B1j_0002
    real(f64) :: absy_0002
    real(f64) :: w_lrh_0002
    real(f64) :: w_lrhu_0002
    real(f64) :: w_lrhv_0002
    real(f64) :: w_lrhB1_0002
    real(f64) :: w_lrhB2_0002
    real(f64) :: w_lrhPSI_0002
    real(f64) :: w_lrz_0002
    real(f64), allocatable, target :: signA_0002(:,:)
    real(f64) :: sound_0002
    real(f64) :: w_0002
    real(f64) :: lambda1_0002
    real(f64) :: lambda2_0002
    real(f64) :: lambda3_0002
    real(f64) :: lambda4_0002
    real(f64) :: epsilon_0002
    real(f64) :: s1_0002
    real(f64) :: pi1_0002
    real(f64) :: s2_0002
    real(f64) :: pi2_0002
    real(f64) :: s3_0002
    real(f64) :: pi3_0002
    real(f64) :: s4_0002
    real(f64) :: pi4_0002
    real(f64) :: gamma1_0002
    real(f64) :: gamma2_0002
    real(f64) :: sigma1_0002
    real(f64) :: sigma2_0002
    real(f64) :: mu1_0002
    real(f64) :: mu2_0002
    integer(i64) :: ann_0002
    real(f64), pointer :: smmat_0002(:,:)
    real(f64) :: hnew_0002
    real(f64) :: unew_0002
    real(f64) :: vnew_0002
    real(f64) :: B1new_0002
    real(f64) :: B2new_0002
    real(f64) :: znew_0002
    real(f64) :: Pnew_0002
    real(f64) :: u_hu_0002
    real(f64) :: u_hv_0002
    real(f64) :: u_hP_0002
    real(f64) :: u_hB1_0002
    real(f64) :: u_hB2_0002
    real(f64) :: u_z_0002
    real(f64) :: w_lrhP_0002
    real(f64) :: w_hP_0002
    real(f64) :: mw_hB1_0002
    real(f64) :: mhP_0002
    real(f64), allocatable :: norm_0002(:)
    real(f64) :: q_s_0002
    real(f64) :: p_s_0002
    real(f64) :: Flux_B1psi_0002
    real(f64) :: Flux_B2psi_0002
    real(f64) :: Flux_hPpsi_0002
    integer(i64) :: i_0002
    integer(i64) :: Dummy_0007
    real(f64), allocatable :: ninv_0003(:)
    real(f64), allocatable :: w_dif_0003(:)
    real(f64) :: u_h_0003
    real(f64) :: v_h_0003
    real(f64) :: B1_h_0003
    real(f64) :: B2_h_0003
    real(f64) :: un_h_0003
    real(f64) :: vn_h_0003
    real(f64) :: B1n_h_0003
    real(f64) :: B2n_h_0003
    real(f64) :: hroe_0003
    real(f64) :: uroe_0003
    real(f64) :: vroe_0003
    real(f64) :: B1roe_0003
    real(f64) :: B2roe_0003
    real(f64) :: uleft_0003
    real(f64) :: vleft_0003
    real(f64) :: B1left_0003
    real(f64) :: B2left_0003
    real(f64) :: uright_0003
    real(f64) :: vright_0003
    real(f64) :: B1right_0003
    real(f64) :: B2right_0003
    real(f64) :: B1i_0003
    real(f64) :: B1j_0003
    real(f64) :: absy_0003
    real(f64) :: w_lrh_0003
    real(f64) :: w_lrhu_0003
    real(f64) :: w_lrhv_0003
    real(f64) :: w_lrhB1_0003
    real(f64) :: w_lrhB2_0003
    real(f64) :: w_lrhPSI_0003
    real(f64) :: w_lrz_0003
    real(f64), allocatable, target :: signA_0003(:,:)
    real(f64) :: sound_0003
    real(f64) :: w_0003
    real(f64) :: lambda1_0003
    real(f64) :: lambda2_0003
    real(f64) :: lambda3_0003
    real(f64) :: lambda4_0003
    real(f64) :: epsilon_0003
    real(f64) :: s1_0003
    real(f64) :: pi1_0003
    real(f64) :: s2_0003
    real(f64) :: pi2_0003
    real(f64) :: s3_0003
    real(f64) :: pi3_0003
    real(f64) :: s4_0003
    real(f64) :: pi4_0003
    real(f64) :: gamma1_0003
    real(f64) :: gamma2_0003
    real(f64) :: sigma1_0003
    real(f64) :: sigma2_0003
    real(f64) :: mu1_0003
    real(f64) :: mu2_0003
    integer(i64) :: ann_0003
    real(f64), pointer :: smmat_0003(:,:)
    real(f64) :: hnew_0003
    real(f64) :: unew_0003
    real(f64) :: vnew_0003
    real(f64) :: B1new_0003
    real(f64) :: B2new_0003
    real(f64) :: znew_0003
    real(f64) :: Pnew_0003
    real(f64) :: u_hu_0003
    real(f64) :: u_hv_0003
    real(f64) :: u_hP_0003
    real(f64) :: u_hB1_0003
    real(f64) :: u_hB2_0003
    real(f64) :: u_z_0003
    real(f64) :: w_lrhP_0003
    real(f64) :: w_hP_0003
    real(f64) :: mw_hB1_0003
    real(f64) :: mhP_0003
    real(f64), allocatable :: norm_0003(:)
    real(f64) :: q_s_0003
    real(f64) :: p_s_0003
    real(f64) :: Flux_B1psi_0003
    real(f64) :: Flux_B2psi_0003
    real(f64) :: Flux_hPpsi_0003
    integer(i64) :: i_0003
    integer(i64) :: Dummy_0008
    real(f64), allocatable :: ninv_0004(:)
    real(f64), allocatable :: w_dif_0004(:)
    real(f64) :: u_h_0004
    real(f64) :: v_h_0004
    real(f64) :: B1_h_0004
    real(f64) :: B2_h_0004
    real(f64) :: un_h_0004
    real(f64) :: vn_h_0004
    real(f64) :: B1n_h_0004
    real(f64) :: B2n_h_0004
    real(f64) :: hroe_0004
    real(f64) :: uroe_0004
    real(f64) :: vroe_0004
    real(f64) :: B1roe_0004
    real(f64) :: B2roe_0004
    real(f64) :: uleft_0004
    real(f64) :: vleft_0004
    real(f64) :: B1left_0004
    real(f64) :: B2left_0004
    real(f64) :: uright_0004
    real(f64) :: vright_0004
    real(f64) :: B1right_0004
    real(f64) :: B2right_0004
    real(f64) :: B1i_0004
    real(f64) :: B1j_0004
    real(f64) :: absy_0004
    real(f64) :: w_lrh_0004
    real(f64) :: w_lrhu_0004
    real(f64) :: w_lrhv_0004
    real(f64) :: w_lrhB1_0004
    real(f64) :: w_lrhB2_0004
    real(f64) :: w_lrhPSI_0004
    real(f64) :: w_lrz_0004
    real(f64), allocatable, target :: signA_0004(:,:)
    real(f64) :: sound_0004
    real(f64) :: w_0004
    real(f64) :: lambda1_0004
    real(f64) :: lambda2_0004
    real(f64) :: lambda3_0004
    real(f64) :: lambda4_0004
    real(f64) :: epsilon_0004
    real(f64) :: s1_0004
    real(f64) :: pi1_0004
    real(f64) :: s2_0004
    real(f64) :: pi2_0004
    real(f64) :: s3_0004
    real(f64) :: pi3_0004
    real(f64) :: s4_0004
    real(f64) :: pi4_0004
    real(f64) :: gamma1_0004
    real(f64) :: gamma2_0004
    real(f64) :: sigma1_0004
    real(f64) :: sigma2_0004
    real(f64) :: mu1_0004
    real(f64) :: mu2_0004
    integer(i64) :: ann_0004
    real(f64), pointer :: smmat_0004(:,:)
    real(f64) :: hnew_0004
    real(f64) :: unew_0004
    real(f64) :: vnew_0004
    real(f64) :: B1new_0004
    real(f64) :: B2new_0004
    real(f64) :: znew_0004
    real(f64) :: Pnew_0004
    real(f64) :: u_hu_0004
    real(f64) :: u_hv_0004
    real(f64) :: u_hP_0004
    real(f64) :: u_hB1_0004
    real(f64) :: u_hB2_0004
    real(f64) :: u_z_0004
    real(f64) :: w_lrhP_0004
    real(f64) :: w_hP_0004
    real(f64) :: mw_hB1_0004
    real(f64) :: mhP_0004
    real(f64), allocatable :: norm_0004(:)
    real(f64) :: q_s_0004
    real(f64) :: p_s_0004
    real(f64) :: Flux_B1psi_0004
    real(f64) :: Flux_B2psi_0004
    real(f64) :: Flux_hPpsi_0004
    integer(i64) :: i_0004

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
    !r_l = zeros(2)
    !r_r = zeros(2)
    !
    do Dummy_0005 = 0_i64, size(innerfaces, kind=i64) - 1_i64, 1_i64
      i = innerfaces(Dummy_0005)
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
      !center_left = centerc[cellidf[i][0]]
      !center_right = centerc[cellidf[i][1]]
      !
      !h_x_left = h_x[cellidf[i][0]];   h_x_right = h_x[cellidf[i][1]]
      !h_y_left = h_y[cellidf[i][0]];   h_y_right = h_y[cellidf[i][1]]
      !
      !psi_left = psi[cellidf[i][0]]; psi_right = psi[cellidf[i][1]]
      !
      !r_l[0] = centerf[i][0] - center_left[0]; r_r[0] = centerf[i][0] - center_right[0];
      !r_l[1] = centerf[i][1] - center_left[1]; r_r[1] = centerf[i][1] - center_right[1];
      !h_l  = h_l  + (order - 1) * psi_left  * (h_x_left * r_l[0]  + h_y_left * r_l[1] )
      !h_r  = h_r  + (order - 1) * psi_right * (h_x_right* r_r[0]  + h_y_right* r_r[1] )
      allocate(ninv_0001(0:1_i64))
      ninv_0001 = 0.0_f64
      allocate(w_dif_0001(0:5_i64))
      w_dif_0001 = 0.0_f64
      ninv_0001(0_i64) = (-1_i64) * normal(1_i64)
      ninv_0001(1_i64) = normal(0_i64)
      u_h_0001 = (hu_l / h_l * sqrt(h_l) + hu_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      v_h_0001 = (hv_l / h_l * sqrt(h_l) + hv_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      B1_h_0001 = (hB1_l / h_l * sqrt(h_l) + hB1_r / h_r * sqrt(h_r)) / &
            (sqrt(h_l) + sqrt(h_r))
      B2_h_0001 = (hB2_l / h_l * sqrt(h_l) + hB2_r / h_r * sqrt(h_r)) / &
            (sqrt(h_l) + sqrt(h_r))
      !uvh =  array([uh, vh])
      un_h_0001 = u_h_0001 * normal(0_i64) + v_h_0001 * normal(1_i64)
      un_h_0001 = un_h_0001 / mesure
      vn_h_0001 = u_h_0001 * ninv_0001(0_i64) + v_h_0001 * ninv_0001( &
            1_i64)
      vn_h_0001 = vn_h_0001 / mesure
      B1n_h_0001 = B1_h_0001 * normal(0_i64) + B2_h_0001 * normal(1_i64)
      B1n_h_0001 = B1n_h_0001 / mesure
      B2n_h_0001 = B1_h_0001 * ninv_0001(0_i64) + B2_h_0001 * ninv_0001( &
            1_i64)
      B2n_h_0001 = B2n_h_0001 / mesure
      hroe_0001 = (h_l + h_r) / 2_i64
      uroe_0001 = un_h_0001
      vroe_0001 = vn_h_0001
      B1roe_0001 = B1n_h_0001
      B2roe_0001 = B2n_h_0001
      uleft_0001 = hu_l * normal(0_i64) + hv_l * normal(1_i64)
      uleft_0001 = uleft_0001 / mesure
      vleft_0001 = hu_l * ninv_0001(0_i64) + hv_l * ninv_0001(1_i64)
      vleft_0001 = vleft_0001 / mesure
      B1left_0001 = hB1_l * normal(0_i64) + hB2_l * normal(1_i64)
      B1left_0001 = B1left_0001 / mesure
      B2left_0001 = hB1_l * ninv_0001(0_i64) + hB2_l * ninv_0001(1_i64)
      B2left_0001 = B2left_0001 / mesure
      uright_0001 = hu_r * normal(0_i64) + hv_r * normal(1_i64)
      uright_0001 = uright_0001 / mesure
      vright_0001 = hu_r * ninv_0001(0_i64) + hv_r * ninv_0001(1_i64)
      vright_0001 = vright_0001 / mesure
      B1right_0001 = hB1_r * normal(0_i64) + hB2_r * normal(1_i64)
      B1right_0001 = B1right_0001 / mesure
      B2right_0001 = hB1_r * ninv_0001(0_i64) + hB2_r * ninv_0001(1_i64)
      B2right_0001 = B2right_0001 / mesure
      B1i_0001 = (hB1c_l * normal(0_i64) + hB2c_l * normal(1_i64)) / &
            mesure
      B1j_0001 = (hB1c_r * normal(0_i64) + hB2c_r * normal(1_i64)) / &
            mesure
      absy_0001 = (B1i_0001 + B1j_0001) / 2_i64
      w_lrh_0001 = (h_l + h_r) / 2_i64
      w_lrhu_0001 = (uleft_0001 + uright_0001) / 2_i64
      w_lrhv_0001 = (vleft_0001 + vright_0001) / 2_i64
      w_lrhB1_0001 = (B1left_0001 + B1right_0001) / 2_i64
      w_lrhB2_0001 = (B2left_0001 + B2right_0001) / 2_i64
      w_lrhPSI_0001 = (hPSI_l + hPSI_r) / 2_i64
      w_lrz_0001 = (Z_l + Z_r) / 2_i64
      w_dif_0001(0_i64) = h_r - h_l
      w_dif_0001(1_i64) = uright_0001 - uleft_0001
      w_dif_0001(2_i64) = vright_0001 - vleft_0001
      w_dif_0001(3_i64) = B1right_0001 - B1left_0001
      w_dif_0001(4_i64) = B2right_0001 - B2left_0001
      w_dif_0001(5_i64) = Z_r - Z_l
      allocate(signA_0001(0:5_i64, 0:5_i64))
      signA_0001 = 0.0_f64
      sound_0001 = sqrt(grav * hroe_0001)
      !B1roe = absy/hroe
      w_0001 = sqrt(B1roe_0001 * B1roe_0001 + grav * hroe_0001)
      lambda1_0001 = uroe_0001 - w_0001
      lambda2_0001 = uroe_0001 - B1roe_0001
      lambda3_0001 = uroe_0001 + B1roe_0001
      lambda4_0001 = uroe_0001 + w_0001
      epsilon_0001 = 1e-15_f64
      if (abs(lambda1_0001) < epsilon_0001) then
        s1_0001 = 0.0_f64
        pi1_0001 = 0.0_f64
      else
        s1_0001 = lambda1_0001 / abs(lambda1_0001)
        pi1_0001 = s1_0001 / lambda1_0001
      end if
      if (abs(lambda2_0001) < epsilon_0001) then
        s2_0001 = 0.0_f64
        pi2_0001 = 0.0_f64
      else
        s2_0001 = lambda2_0001 / abs(lambda2_0001)
        pi2_0001 = 1.0_f64 / abs(lambda2_0001)
      end if
      if (abs(lambda3_0001) < epsilon_0001) then
        s3_0001 = 0.0_f64
        pi3_0001 = 0.0_f64
      else
        s3_0001 = lambda3_0001 / abs(lambda3_0001)
        pi3_0001 = 1.0_f64 / abs(lambda3_0001)
      end if
      if (abs(lambda4_0001) < epsilon_0001) then
        s4_0001 = 0.0_f64
        pi4_0001 = 0.0_f64
      else
        s4_0001 = lambda4_0001 / abs(lambda4_0001)
        pi4_0001 = 1.0_f64 / abs(lambda4_0001)
      end if
      gamma1_0001 = vroe_0001 + B2roe_0001
      gamma2_0001 = vroe_0001 - B2roe_0001
      sigma1_0001 = vroe_0001 * (s1_0001 * lambda4_0001 - s4_0001 * &
            lambda1_0001) - w_0001 * (s2_0001 * gamma1_0001 + s3_0001 * &
            gamma2_0001)
      sigma2_0001 = B2roe_0001 * (s1_0001 * lambda4_0001 - s4_0001 * &
            lambda1_0001) - w_0001 * (s2_0001 * gamma1_0001 - s3_0001 * &
            gamma2_0001)
      if (abs(lambda2_0001) < epsilon_0001 .and. abs(lambda3_0001) < &
            epsilon_0001) then
        mu1_0001 = B1roe_0001 * vroe_0001 * pi1_0001 / w_0001 - &
              B1roe_0001 * vroe_0001 * pi4_0001 / w_0001
        mu2_0001 = B1roe_0001 * B2roe_0001 * pi1_0001 / w_0001 - &
              B1roe_0001 * B2roe_0001 * pi4_0001 / w_0001
        ann_0001 = 1_i64
      else
        mu1_0001 = B1roe_0001 * vroe_0001 * pi1_0001 / w_0001 - &
              B1roe_0001 * vroe_0001 * pi4_0001 / w_0001 - 0.5_f64 * ( &
              gamma1_0001 * pi2_0001 - gamma2_0001 * pi3_0001)
        mu2_0001 = B1roe_0001 * B2roe_0001 * pi1_0001 / w_0001 - &
              B1roe_0001 * B2roe_0001 * pi4_0001 / w_0001 - 0.5_f64 * ( &
              gamma1_0001 * pi2_0001 + gamma2_0001 * pi3_0001)
        ann_0001 = 1_i64
      end if
      !1ère colonne de la matrice A
      signA_0001(0_i64, 0_i64) = (s1_0001 * lambda4_0001 - s4_0001 * &
            lambda1_0001) / (2_i64 * w_0001)
      signA_0001(0_i64, 1_i64) = lambda1_0001 * lambda4_0001 * (s1_0001 &
            - s4_0001) / (2_i64 * w_0001)
      signA_0001(0_i64, 2_i64) = sigma1_0001 / (2_i64 * w_0001)
      signA_0001(0_i64, 3_i64) = 0.0_f64
      signA_0001(0_i64, 4_i64) = sigma2_0001 / (2_i64 * w_0001)
      signA_0001(0_i64, 5_i64) = 0.0_f64
      !2ème colonne de la matrice A
      signA_0001(1_i64, 0_i64) = (s4_0001 - s1_0001) / (2_i64 * w_0001)
      signA_0001(1_i64, 1_i64) = (s4_0001 * lambda4_0001 - s1_0001 * &
            lambda1_0001) / (2_i64 * w_0001)
      signA_0001(1_i64, 2_i64) = vroe_0001 * (s4_0001 - s1_0001) / ( &
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
      signA_0001(5_i64, 0_i64) = sound_0001 ** 2_i64 * (pi4_0001 - &
            pi1_0001) / (2_i64 * w_0001)
      signA_0001(5_i64, 1_i64) = sound_0001 ** 2_i64 * (s4_0001 - &
            s1_0001) / (2_i64 * w_0001)
      signA_0001(5_i64, 2_i64) = sound_0001 ** 2_i64 * vroe_0001 * ( &
            pi4_0001 - pi1_0001) / (2_i64 * w_0001)
      signA_0001(5_i64, 3_i64) = 0.0_f64
      signA_0001(5_i64, 4_i64) = sound_0001 ** 2_i64 * B2roe_0001 * ( &
            pi4_0001 - pi1_0001) / (2_i64 * w_0001)
      signA_0001(5_i64, 5_i64) = 0.0_f64
      smmat_0001(0:, 0:) => signA_0001
      hnew_0001 = 0.0_f64
      unew_0001 = 0.0_f64
      vnew_0001 = 0.0_f64
      B1new_0001 = 0.0_f64
      B2new_0001 = 0.0_f64
      znew_0001 = 0.0_f64
      do i_0001 = 0_i64, 5_i64, 1_i64
        hnew_0001 = hnew_0001 + smmat_0001(i_0001, 0_i64) * w_dif_0001( &
              i_0001)
        unew_0001 = unew_0001 + smmat_0001(i_0001, 1_i64) * w_dif_0001( &
              i_0001)
        vnew_0001 = vnew_0001 + smmat_0001(i_0001, 2_i64) * w_dif_0001( &
              i_0001)
        B1new_0001 = B1new_0001 + smmat_0001(i_0001, 3_i64) * w_dif_0001 &
              (i_0001)
        B2new_0001 = B2new_0001 + smmat_0001(i_0001, 4_i64) * w_dif_0001 &
              (i_0001)
        znew_0001 = znew_0001 + smmat_0001(i_0001, 5_i64) * w_dif_0001( &
              i_0001)
      end do
      Pnew_0001 = cpsi * (B1right_0001 - B1left_0001)
      u_h_0001 = hnew_0001 / 2_i64
      u_hu_0001 = unew_0001 / 2_i64
      u_hv_0001 = vnew_0001 / 2_i64
      u_hP_0001 = Pnew_0001 / 2_i64
      u_hB1_0001 = B1new_0001 / 2_i64
      u_hB2_0001 = B2new_0001 / 2_i64
      u_z_0001 = znew_0001 / 2_i64
      w_lrh_0001 = w_lrh_0001 - u_h_0001
      w_lrhu_0001 = w_lrhu_0001 - u_hu_0001
      w_lrhv_0001 = w_lrhv_0001 - u_hv_0001
      w_lrhP_0001 = w_lrhPSI_0001 - u_hP_0001
      w_lrhB1_0001 = w_lrhB1_0001 - u_hB1_0001
      w_lrhB2_0001 = w_lrhB2_0001 - u_hB2_0001
      w_lrz_0001 = w_lrz_0001 - u_z_0001
      w_hP_0001 = hPSI_r - hPSI_l
      mw_hB1_0001 = w_lrhB1_0001 / 2_i64 - 1_i64 / (2_i64 * cpsi) * &
            w_hP_0001
      mhP_0001 = (hPSI_r - hPSI_l) / 2_i64 - cpsi * w_dif_0001(3_i64) / &
            2_i64
      unew_0001 = 0.0_f64
      vnew_0001 = 0.0_f64
      B1new_0001 = 0.0_f64
      B2new_0001 = 0.0_f64
      unew_0001 = w_lrhu_0001 * normal(0_i64) - w_lrhv_0001 * normal( &
            1_i64)
      unew_0001 = unew_0001 / mesure
      vnew_0001 = w_lrhu_0001 * normal(1_i64) + w_lrhv_0001 * normal( &
            0_i64)
      vnew_0001 = vnew_0001 / mesure
      B1new_0001 = w_lrhB1_0001 * normal(0_i64) - w_lrhB2_0001 * normal( &
            1_i64)
      B1new_0001 = B1new_0001 / mesure
      B2new_0001 = w_lrhB1_0001 * normal(1_i64) + w_lrhB2_0001 * normal( &
            0_i64)
      B2new_0001 = B2new_0001 / mesure
      w_lrhu_0001 = unew_0001
      w_lrhv_0001 = vnew_0001
      w_lrhB1_0001 = B1new_0001
      w_lrhB2_0001 = B2new_0001
      allocate(norm_0001(0:size(normal, kind=i64) - 1_i64))
      norm_0001 = normal / mesure
      q_s_0001 = normal(0_i64) * unew_0001 + normal(1_i64) * vnew_0001
      p_s_0001 = normal(0_i64) * B1new_0001 + normal(1_i64) * B2new_0001
      Flux_B1psi_0001 = mhP_0001 * norm_0001(0_i64) * mesure
      Flux_B2psi_0001 = mhP_0001 * norm_0001(1_i64) * mesure
      Flux_hPpsi_0001 = cpsi * cpsi * mw_hB1_0001 * mesure
      flux(0_i64) = q_s_0001
      flux(1_i64) = q_s_0001 * w_lrhu_0001 / w_lrh_0001 + 0.5_f64 * grav &
            * w_lrh_0001 * w_lrh_0001 * normal(0_i64) - p_s_0001 * &
            w_lrhB1_0001 / w_lrh_0001
      flux(2_i64) = q_s_0001 * w_lrhv_0001 / w_lrh_0001 + 0.5_f64 * grav &
            * w_lrh_0001 * w_lrh_0001 * normal(1_i64) - p_s_0001 * &
            w_lrhB2_0001 / w_lrh_0001
      flux(3_i64) = (w_lrhv_0001 * w_lrhB1_0001 / w_lrh_0001 - &
            w_lrhu_0001 * w_lrhB2_0001 / w_lrh_0001) * normal(1_i64) + &
            Flux_B1psi_0001
      flux(4_i64) = (w_lrhu_0001 * w_lrhB2_0001 / w_lrh_0001 - &
            w_lrhv_0001 * w_lrhB1_0001 / w_lrh_0001) * normal(0_i64) + &
            Flux_B2psi_0001
      flux(5_i64) = Flux_hPpsi_0001
      flux(6_i64) = 0.0_f64
      if (allocated(ninv_0001)) then
        deallocate(ninv_0001)
      end if
      if (allocated(w_dif_0001)) then
        deallocate(w_dif_0001)
      end if
      if (allocated(signA_0001)) then
        deallocate(signA_0001)
      end if
      if (allocated(norm_0001)) then
        deallocate(norm_0001)
      end if
      !LF_scheme_MHD(hu_l, hu_r, hv_l, hv_r, h_l, h_r, hB1_l, hB1_r, hB2_l, hB2_r, Z_l, Z_r, normal, mesure, grav, flux)
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
    do Dummy_0006 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0006)
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
      !
      !center_left = centerc[cellidf[i][0]]
      !center_right = centerh[halofid[i]]
      !
      !h_x_left = h_x[cellidf[i][0]];   h_x_right = hx_halo[halofid[i]]
      !h_y_left = h_y[cellidf[i][0]];   h_y_right = hy_halo[halofid[i]]
      !
      !psi_left = psi[cellidf[i][0]]; psi_right = psi_halo[halofid[i]]
      !
      !r_l[0] = centerf[i][0] - center_left[0]; r_r[0] = centerf[i][0] - center_right[0];
      !r_l[1] = centerf[i][1] - center_left[1]; r_r[1] = centerf[i][1] - center_right[1];
      !h_l  = h_l  + (order - 1) * psi_left  * (h_x_left   * r_l[0] + h_y_left   * r_l[1] )
      !h_r  = h_r  + (order - 1) * psi_right * (h_x_right  * r_r[0] + h_y_right  * r_r[1] )
      allocate(ninv_0002(0:1_i64))
      ninv_0002 = 0.0_f64
      allocate(w_dif_0002(0:5_i64))
      w_dif_0002 = 0.0_f64
      ninv_0002(0_i64) = (-1_i64) * normal(1_i64)
      ninv_0002(1_i64) = normal(0_i64)
      u_h_0002 = (hu_l / h_l * sqrt(h_l) + hu_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      v_h_0002 = (hv_l / h_l * sqrt(h_l) + hv_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      B1_h_0002 = (hB1_l / h_l * sqrt(h_l) + hB1_r / h_r * sqrt(h_r)) / &
            (sqrt(h_l) + sqrt(h_r))
      B2_h_0002 = (hB2_l / h_l * sqrt(h_l) + hB2_r / h_r * sqrt(h_r)) / &
            (sqrt(h_l) + sqrt(h_r))
      !uvh =  array([uh, vh])
      un_h_0002 = u_h_0002 * normal(0_i64) + v_h_0002 * normal(1_i64)
      un_h_0002 = un_h_0002 / mesure
      vn_h_0002 = u_h_0002 * ninv_0002(0_i64) + v_h_0002 * ninv_0002( &
            1_i64)
      vn_h_0002 = vn_h_0002 / mesure
      B1n_h_0002 = B1_h_0002 * normal(0_i64) + B2_h_0002 * normal(1_i64)
      B1n_h_0002 = B1n_h_0002 / mesure
      B2n_h_0002 = B1_h_0002 * ninv_0002(0_i64) + B2_h_0002 * ninv_0002( &
            1_i64)
      B2n_h_0002 = B2n_h_0002 / mesure
      hroe_0002 = (h_l + h_r) / 2_i64
      uroe_0002 = un_h_0002
      vroe_0002 = vn_h_0002
      B1roe_0002 = B1n_h_0002
      B2roe_0002 = B2n_h_0002
      uleft_0002 = hu_l * normal(0_i64) + hv_l * normal(1_i64)
      uleft_0002 = uleft_0002 / mesure
      vleft_0002 = hu_l * ninv_0002(0_i64) + hv_l * ninv_0002(1_i64)
      vleft_0002 = vleft_0002 / mesure
      B1left_0002 = hB1_l * normal(0_i64) + hB2_l * normal(1_i64)
      B1left_0002 = B1left_0002 / mesure
      B2left_0002 = hB1_l * ninv_0002(0_i64) + hB2_l * ninv_0002(1_i64)
      B2left_0002 = B2left_0002 / mesure
      uright_0002 = hu_r * normal(0_i64) + hv_r * normal(1_i64)
      uright_0002 = uright_0002 / mesure
      vright_0002 = hu_r * ninv_0002(0_i64) + hv_r * ninv_0002(1_i64)
      vright_0002 = vright_0002 / mesure
      B1right_0002 = hB1_r * normal(0_i64) + hB2_r * normal(1_i64)
      B1right_0002 = B1right_0002 / mesure
      B2right_0002 = hB1_r * ninv_0002(0_i64) + hB2_r * ninv_0002(1_i64)
      B2right_0002 = B2right_0002 / mesure
      B1i_0002 = (hB1c_l * normal(0_i64) + hB2c_l * normal(1_i64)) / &
            mesure
      B1j_0002 = (hB1c_r * normal(0_i64) + hB2c_r * normal(1_i64)) / &
            mesure
      absy_0002 = (B1i_0002 + B1j_0002) / 2_i64
      w_lrh_0002 = (h_l + h_r) / 2_i64
      w_lrhu_0002 = (uleft_0002 + uright_0002) / 2_i64
      w_lrhv_0002 = (vleft_0002 + vright_0002) / 2_i64
      w_lrhB1_0002 = (B1left_0002 + B1right_0002) / 2_i64
      w_lrhB2_0002 = (B2left_0002 + B2right_0002) / 2_i64
      w_lrhPSI_0002 = (hPSI_l + hPSI_r) / 2_i64
      w_lrz_0002 = (Z_l + Z_r) / 2_i64
      w_dif_0002(0_i64) = h_r - h_l
      w_dif_0002(1_i64) = uright_0002 - uleft_0002
      w_dif_0002(2_i64) = vright_0002 - vleft_0002
      w_dif_0002(3_i64) = B1right_0002 - B1left_0002
      w_dif_0002(4_i64) = B2right_0002 - B2left_0002
      w_dif_0002(5_i64) = Z_r - Z_l
      allocate(signA_0002(0:5_i64, 0:5_i64))
      signA_0002 = 0.0_f64
      sound_0002 = sqrt(grav * hroe_0002)
      !B1roe = absy/hroe
      w_0002 = sqrt(B1roe_0002 * B1roe_0002 + grav * hroe_0002)
      lambda1_0002 = uroe_0002 - w_0002
      lambda2_0002 = uroe_0002 - B1roe_0002
      lambda3_0002 = uroe_0002 + B1roe_0002
      lambda4_0002 = uroe_0002 + w_0002
      epsilon_0002 = 1e-15_f64
      if (abs(lambda1_0002) < epsilon_0002) then
        s1_0002 = 0.0_f64
        pi1_0002 = 0.0_f64
      else
        s1_0002 = lambda1_0002 / abs(lambda1_0002)
        pi1_0002 = s1_0002 / lambda1_0002
      end if
      if (abs(lambda2_0002) < epsilon_0002) then
        s2_0002 = 0.0_f64
        pi2_0002 = 0.0_f64
      else
        s2_0002 = lambda2_0002 / abs(lambda2_0002)
        pi2_0002 = 1.0_f64 / abs(lambda2_0002)
      end if
      if (abs(lambda3_0002) < epsilon_0002) then
        s3_0002 = 0.0_f64
        pi3_0002 = 0.0_f64
      else
        s3_0002 = lambda3_0002 / abs(lambda3_0002)
        pi3_0002 = 1.0_f64 / abs(lambda3_0002)
      end if
      if (abs(lambda4_0002) < epsilon_0002) then
        s4_0002 = 0.0_f64
        pi4_0002 = 0.0_f64
      else
        s4_0002 = lambda4_0002 / abs(lambda4_0002)
        pi4_0002 = 1.0_f64 / abs(lambda4_0002)
      end if
      gamma1_0002 = vroe_0002 + B2roe_0002
      gamma2_0002 = vroe_0002 - B2roe_0002
      sigma1_0002 = vroe_0002 * (s1_0002 * lambda4_0002 - s4_0002 * &
            lambda1_0002) - w_0002 * (s2_0002 * gamma1_0002 + s3_0002 * &
            gamma2_0002)
      sigma2_0002 = B2roe_0002 * (s1_0002 * lambda4_0002 - s4_0002 * &
            lambda1_0002) - w_0002 * (s2_0002 * gamma1_0002 - s3_0002 * &
            gamma2_0002)
      if (abs(lambda2_0002) < epsilon_0002 .and. abs(lambda3_0002) < &
            epsilon_0002) then
        mu1_0002 = B1roe_0002 * vroe_0002 * pi1_0002 / w_0002 - &
              B1roe_0002 * vroe_0002 * pi4_0002 / w_0002
        mu2_0002 = B1roe_0002 * B2roe_0002 * pi1_0002 / w_0002 - &
              B1roe_0002 * B2roe_0002 * pi4_0002 / w_0002
        ann_0002 = 1_i64
      else
        mu1_0002 = B1roe_0002 * vroe_0002 * pi1_0002 / w_0002 - &
              B1roe_0002 * vroe_0002 * pi4_0002 / w_0002 - 0.5_f64 * ( &
              gamma1_0002 * pi2_0002 - gamma2_0002 * pi3_0002)
        mu2_0002 = B1roe_0002 * B2roe_0002 * pi1_0002 / w_0002 - &
              B1roe_0002 * B2roe_0002 * pi4_0002 / w_0002 - 0.5_f64 * ( &
              gamma1_0002 * pi2_0002 + gamma2_0002 * pi3_0002)
        ann_0002 = 1_i64
      end if
      !1ère colonne de la matrice A
      signA_0002(0_i64, 0_i64) = (s1_0002 * lambda4_0002 - s4_0002 * &
            lambda1_0002) / (2_i64 * w_0002)
      signA_0002(0_i64, 1_i64) = lambda1_0002 * lambda4_0002 * (s1_0002 &
            - s4_0002) / (2_i64 * w_0002)
      signA_0002(0_i64, 2_i64) = sigma1_0002 / (2_i64 * w_0002)
      signA_0002(0_i64, 3_i64) = 0.0_f64
      signA_0002(0_i64, 4_i64) = sigma2_0002 / (2_i64 * w_0002)
      signA_0002(0_i64, 5_i64) = 0.0_f64
      !2ème colonne de la matrice A
      signA_0002(1_i64, 0_i64) = (s4_0002 - s1_0002) / (2_i64 * w_0002)
      signA_0002(1_i64, 1_i64) = (s4_0002 * lambda4_0002 - s1_0002 * &
            lambda1_0002) / (2_i64 * w_0002)
      signA_0002(1_i64, 2_i64) = vroe_0002 * (s4_0002 - s1_0002) / ( &
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
      signA_0002(5_i64, 0_i64) = sound_0002 ** 2_i64 * (pi4_0002 - &
            pi1_0002) / (2_i64 * w_0002)
      signA_0002(5_i64, 1_i64) = sound_0002 ** 2_i64 * (s4_0002 - &
            s1_0002) / (2_i64 * w_0002)
      signA_0002(5_i64, 2_i64) = sound_0002 ** 2_i64 * vroe_0002 * ( &
            pi4_0002 - pi1_0002) / (2_i64 * w_0002)
      signA_0002(5_i64, 3_i64) = 0.0_f64
      signA_0002(5_i64, 4_i64) = sound_0002 ** 2_i64 * B2roe_0002 * ( &
            pi4_0002 - pi1_0002) / (2_i64 * w_0002)
      signA_0002(5_i64, 5_i64) = 0.0_f64
      smmat_0002(0:, 0:) => signA_0002
      hnew_0002 = 0.0_f64
      unew_0002 = 0.0_f64
      vnew_0002 = 0.0_f64
      B1new_0002 = 0.0_f64
      B2new_0002 = 0.0_f64
      znew_0002 = 0.0_f64
      do i_0002 = 0_i64, 5_i64, 1_i64
        hnew_0002 = hnew_0002 + smmat_0002(i_0002, 0_i64) * w_dif_0002( &
              i_0002)
        unew_0002 = unew_0002 + smmat_0002(i_0002, 1_i64) * w_dif_0002( &
              i_0002)
        vnew_0002 = vnew_0002 + smmat_0002(i_0002, 2_i64) * w_dif_0002( &
              i_0002)
        B1new_0002 = B1new_0002 + smmat_0002(i_0002, 3_i64) * w_dif_0002 &
              (i_0002)
        B2new_0002 = B2new_0002 + smmat_0002(i_0002, 4_i64) * w_dif_0002 &
              (i_0002)
        znew_0002 = znew_0002 + smmat_0002(i_0002, 5_i64) * w_dif_0002( &
              i_0002)
      end do
      Pnew_0002 = cpsi * (B1right_0002 - B1left_0002)
      u_h_0002 = hnew_0002 / 2_i64
      u_hu_0002 = unew_0002 / 2_i64
      u_hv_0002 = vnew_0002 / 2_i64
      u_hP_0002 = Pnew_0002 / 2_i64
      u_hB1_0002 = B1new_0002 / 2_i64
      u_hB2_0002 = B2new_0002 / 2_i64
      u_z_0002 = znew_0002 / 2_i64
      w_lrh_0002 = w_lrh_0002 - u_h_0002
      w_lrhu_0002 = w_lrhu_0002 - u_hu_0002
      w_lrhv_0002 = w_lrhv_0002 - u_hv_0002
      w_lrhP_0002 = w_lrhPSI_0002 - u_hP_0002
      w_lrhB1_0002 = w_lrhB1_0002 - u_hB1_0002
      w_lrhB2_0002 = w_lrhB2_0002 - u_hB2_0002
      w_lrz_0002 = w_lrz_0002 - u_z_0002
      w_hP_0002 = hPSI_r - hPSI_l
      mw_hB1_0002 = w_lrhB1_0002 / 2_i64 - 1_i64 / (2_i64 * cpsi) * &
            w_hP_0002
      mhP_0002 = (hPSI_r - hPSI_l) / 2_i64 - cpsi * w_dif_0002(3_i64) / &
            2_i64
      unew_0002 = 0.0_f64
      vnew_0002 = 0.0_f64
      B1new_0002 = 0.0_f64
      B2new_0002 = 0.0_f64
      unew_0002 = w_lrhu_0002 * normal(0_i64) - w_lrhv_0002 * normal( &
            1_i64)
      unew_0002 = unew_0002 / mesure
      vnew_0002 = w_lrhu_0002 * normal(1_i64) + w_lrhv_0002 * normal( &
            0_i64)
      vnew_0002 = vnew_0002 / mesure
      B1new_0002 = w_lrhB1_0002 * normal(0_i64) - w_lrhB2_0002 * normal( &
            1_i64)
      B1new_0002 = B1new_0002 / mesure
      B2new_0002 = w_lrhB1_0002 * normal(1_i64) + w_lrhB2_0002 * normal( &
            0_i64)
      B2new_0002 = B2new_0002 / mesure
      w_lrhu_0002 = unew_0002
      w_lrhv_0002 = vnew_0002
      w_lrhB1_0002 = B1new_0002
      w_lrhB2_0002 = B2new_0002
      allocate(norm_0002(0:size(normal, kind=i64) - 1_i64))
      norm_0002 = normal / mesure
      q_s_0002 = normal(0_i64) * unew_0002 + normal(1_i64) * vnew_0002
      p_s_0002 = normal(0_i64) * B1new_0002 + normal(1_i64) * B2new_0002
      Flux_B1psi_0002 = mhP_0002 * norm_0002(0_i64) * mesure
      Flux_B2psi_0002 = mhP_0002 * norm_0002(1_i64) * mesure
      Flux_hPpsi_0002 = cpsi * cpsi * mw_hB1_0002 * mesure
      flux(0_i64) = q_s_0002
      flux(1_i64) = q_s_0002 * w_lrhu_0002 / w_lrh_0002 + 0.5_f64 * grav &
            * w_lrh_0002 * w_lrh_0002 * normal(0_i64) - p_s_0002 * &
            w_lrhB1_0002 / w_lrh_0002
      flux(2_i64) = q_s_0002 * w_lrhv_0002 / w_lrh_0002 + 0.5_f64 * grav &
            * w_lrh_0002 * w_lrh_0002 * normal(1_i64) - p_s_0002 * &
            w_lrhB2_0002 / w_lrh_0002
      flux(3_i64) = (w_lrhv_0002 * w_lrhB1_0002 / w_lrh_0002 - &
            w_lrhu_0002 * w_lrhB2_0002 / w_lrh_0002) * normal(1_i64) + &
            Flux_B1psi_0002
      flux(4_i64) = (w_lrhu_0002 * w_lrhB2_0002 / w_lrh_0002 - &
            w_lrhv_0002 * w_lrhB1_0002 / w_lrh_0002) * normal(0_i64) + &
            Flux_B2psi_0002
      flux(5_i64) = Flux_hPpsi_0002
      flux(6_i64) = 0.0_f64
      if (allocated(ninv_0002)) then
        deallocate(ninv_0002)
      end if
      if (allocated(w_dif_0002)) then
        deallocate(w_dif_0002)
      end if
      if (allocated(signA_0002)) then
        deallocate(signA_0002)
      end if
      if (allocated(norm_0002)) then
        deallocate(norm_0002)
      end if
      !LF_scheme_MHD(hu_l, hu_r, hv_l, hv_r, h_l, h_r, hB1_l, hB1_r, hB2_l, hB2_r, Z_l, Z_r, normal, mesure, grav, flux)
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
    do Dummy_0007 = 0_i64, size(periodicboundaryfaces, kind=i64) - 1_i64 &
          , 1_i64
      i = periodicboundaryfaces(Dummy_0007)
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
      allocate(ninv_0003(0:1_i64))
      ninv_0003 = 0.0_f64
      allocate(w_dif_0003(0:5_i64))
      w_dif_0003 = 0.0_f64
      ninv_0003(0_i64) = (-1_i64) * normal(1_i64)
      ninv_0003(1_i64) = normal(0_i64)
      u_h_0003 = (hu_l / h_l * sqrt(h_l) + hu_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      v_h_0003 = (hv_l / h_l * sqrt(h_l) + hv_r / h_r * sqrt(h_r)) / ( &
            sqrt(h_l) + sqrt(h_r))
      B1_h_0003 = (hB1_l / h_l * sqrt(h_l) + hB1_r / h_r * sqrt(h_r)) / &
            (sqrt(h_l) + sqrt(h_r))
      B2_h_0003 = (hB2_l / h_l * sqrt(h_l) + hB2_r / h_r * sqrt(h_r)) / &
            (sqrt(h_l) + sqrt(h_r))
      !uvh =  array([uh, vh])
      un_h_0003 = u_h_0003 * normal(0_i64) + v_h_0003 * normal(1_i64)
      un_h_0003 = un_h_0003 / mesure
      vn_h_0003 = u_h_0003 * ninv_0003(0_i64) + v_h_0003 * ninv_0003( &
            1_i64)
      vn_h_0003 = vn_h_0003 / mesure
      B1n_h_0003 = B1_h_0003 * normal(0_i64) + B2_h_0003 * normal(1_i64)
      B1n_h_0003 = B1n_h_0003 / mesure
      B2n_h_0003 = B1_h_0003 * ninv_0003(0_i64) + B2_h_0003 * ninv_0003( &
            1_i64)
      B2n_h_0003 = B2n_h_0003 / mesure
      hroe_0003 = (h_l + h_r) / 2_i64
      uroe_0003 = un_h_0003
      vroe_0003 = vn_h_0003
      B1roe_0003 = B1n_h_0003
      B2roe_0003 = B2n_h_0003
      uleft_0003 = hu_l * normal(0_i64) + hv_l * normal(1_i64)
      uleft_0003 = uleft_0003 / mesure
      vleft_0003 = hu_l * ninv_0003(0_i64) + hv_l * ninv_0003(1_i64)
      vleft_0003 = vleft_0003 / mesure
      B1left_0003 = hB1_l * normal(0_i64) + hB2_l * normal(1_i64)
      B1left_0003 = B1left_0003 / mesure
      B2left_0003 = hB1_l * ninv_0003(0_i64) + hB2_l * ninv_0003(1_i64)
      B2left_0003 = B2left_0003 / mesure
      uright_0003 = hu_r * normal(0_i64) + hv_r * normal(1_i64)
      uright_0003 = uright_0003 / mesure
      vright_0003 = hu_r * ninv_0003(0_i64) + hv_r * ninv_0003(1_i64)
      vright_0003 = vright_0003 / mesure
      B1right_0003 = hB1_r * normal(0_i64) + hB2_r * normal(1_i64)
      B1right_0003 = B1right_0003 / mesure
      B2right_0003 = hB1_r * ninv_0003(0_i64) + hB2_r * ninv_0003(1_i64)
      B2right_0003 = B2right_0003 / mesure
      B1i_0003 = (hB1c_l * normal(0_i64) + hB2c_l * normal(1_i64)) / &
            mesure
      B1j_0003 = (hB1c_r * normal(0_i64) + hB2c_r * normal(1_i64)) / &
            mesure
      absy_0003 = (B1i_0003 + B1j_0003) / 2_i64
      w_lrh_0003 = (h_l + h_r) / 2_i64
      w_lrhu_0003 = (uleft_0003 + uright_0003) / 2_i64
      w_lrhv_0003 = (vleft_0003 + vright_0003) / 2_i64
      w_lrhB1_0003 = (B1left_0003 + B1right_0003) / 2_i64
      w_lrhB2_0003 = (B2left_0003 + B2right_0003) / 2_i64
      w_lrhPSI_0003 = (hPSI_l + hPSI_r) / 2_i64
      w_lrz_0003 = (Z_l + Z_r) / 2_i64
      w_dif_0003(0_i64) = h_r - h_l
      w_dif_0003(1_i64) = uright_0003 - uleft_0003
      w_dif_0003(2_i64) = vright_0003 - vleft_0003
      w_dif_0003(3_i64) = B1right_0003 - B1left_0003
      w_dif_0003(4_i64) = B2right_0003 - B2left_0003
      w_dif_0003(5_i64) = Z_r - Z_l
      allocate(signA_0003(0:5_i64, 0:5_i64))
      signA_0003 = 0.0_f64
      sound_0003 = sqrt(grav * hroe_0003)
      !B1roe = absy/hroe
      w_0003 = sqrt(B1roe_0003 * B1roe_0003 + grav * hroe_0003)
      lambda1_0003 = uroe_0003 - w_0003
      lambda2_0003 = uroe_0003 - B1roe_0003
      lambda3_0003 = uroe_0003 + B1roe_0003
      lambda4_0003 = uroe_0003 + w_0003
      epsilon_0003 = 1e-15_f64
      if (abs(lambda1_0003) < epsilon_0003) then
        s1_0003 = 0.0_f64
        pi1_0003 = 0.0_f64
      else
        s1_0003 = lambda1_0003 / abs(lambda1_0003)
        pi1_0003 = s1_0003 / lambda1_0003
      end if
      if (abs(lambda2_0003) < epsilon_0003) then
        s2_0003 = 0.0_f64
        pi2_0003 = 0.0_f64
      else
        s2_0003 = lambda2_0003 / abs(lambda2_0003)
        pi2_0003 = 1.0_f64 / abs(lambda2_0003)
      end if
      if (abs(lambda3_0003) < epsilon_0003) then
        s3_0003 = 0.0_f64
        pi3_0003 = 0.0_f64
      else
        s3_0003 = lambda3_0003 / abs(lambda3_0003)
        pi3_0003 = 1.0_f64 / abs(lambda3_0003)
      end if
      if (abs(lambda4_0003) < epsilon_0003) then
        s4_0003 = 0.0_f64
        pi4_0003 = 0.0_f64
      else
        s4_0003 = lambda4_0003 / abs(lambda4_0003)
        pi4_0003 = 1.0_f64 / abs(lambda4_0003)
      end if
      gamma1_0003 = vroe_0003 + B2roe_0003
      gamma2_0003 = vroe_0003 - B2roe_0003
      sigma1_0003 = vroe_0003 * (s1_0003 * lambda4_0003 - s4_0003 * &
            lambda1_0003) - w_0003 * (s2_0003 * gamma1_0003 + s3_0003 * &
            gamma2_0003)
      sigma2_0003 = B2roe_0003 * (s1_0003 * lambda4_0003 - s4_0003 * &
            lambda1_0003) - w_0003 * (s2_0003 * gamma1_0003 - s3_0003 * &
            gamma2_0003)
      if (abs(lambda2_0003) < epsilon_0003 .and. abs(lambda3_0003) < &
            epsilon_0003) then
        mu1_0003 = B1roe_0003 * vroe_0003 * pi1_0003 / w_0003 - &
              B1roe_0003 * vroe_0003 * pi4_0003 / w_0003
        mu2_0003 = B1roe_0003 * B2roe_0003 * pi1_0003 / w_0003 - &
              B1roe_0003 * B2roe_0003 * pi4_0003 / w_0003
        ann_0003 = 1_i64
      else
        mu1_0003 = B1roe_0003 * vroe_0003 * pi1_0003 / w_0003 - &
              B1roe_0003 * vroe_0003 * pi4_0003 / w_0003 - 0.5_f64 * ( &
              gamma1_0003 * pi2_0003 - gamma2_0003 * pi3_0003)
        mu2_0003 = B1roe_0003 * B2roe_0003 * pi1_0003 / w_0003 - &
              B1roe_0003 * B2roe_0003 * pi4_0003 / w_0003 - 0.5_f64 * ( &
              gamma1_0003 * pi2_0003 + gamma2_0003 * pi3_0003)
        ann_0003 = 1_i64
      end if
      !1ère colonne de la matrice A
      signA_0003(0_i64, 0_i64) = (s1_0003 * lambda4_0003 - s4_0003 * &
            lambda1_0003) / (2_i64 * w_0003)
      signA_0003(0_i64, 1_i64) = lambda1_0003 * lambda4_0003 * (s1_0003 &
            - s4_0003) / (2_i64 * w_0003)
      signA_0003(0_i64, 2_i64) = sigma1_0003 / (2_i64 * w_0003)
      signA_0003(0_i64, 3_i64) = 0.0_f64
      signA_0003(0_i64, 4_i64) = sigma2_0003 / (2_i64 * w_0003)
      signA_0003(0_i64, 5_i64) = 0.0_f64
      !2ème colonne de la matrice A
      signA_0003(1_i64, 0_i64) = (s4_0003 - s1_0003) / (2_i64 * w_0003)
      signA_0003(1_i64, 1_i64) = (s4_0003 * lambda4_0003 - s1_0003 * &
            lambda1_0003) / (2_i64 * w_0003)
      signA_0003(1_i64, 2_i64) = vroe_0003 * (s4_0003 - s1_0003) / ( &
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
      signA_0003(5_i64, 0_i64) = sound_0003 ** 2_i64 * (pi4_0003 - &
            pi1_0003) / (2_i64 * w_0003)
      signA_0003(5_i64, 1_i64) = sound_0003 ** 2_i64 * (s4_0003 - &
            s1_0003) / (2_i64 * w_0003)
      signA_0003(5_i64, 2_i64) = sound_0003 ** 2_i64 * vroe_0003 * ( &
            pi4_0003 - pi1_0003) / (2_i64 * w_0003)
      signA_0003(5_i64, 3_i64) = 0.0_f64
      signA_0003(5_i64, 4_i64) = sound_0003 ** 2_i64 * B2roe_0003 * ( &
            pi4_0003 - pi1_0003) / (2_i64 * w_0003)
      signA_0003(5_i64, 5_i64) = 0.0_f64
      smmat_0003(0:, 0:) => signA_0003
      hnew_0003 = 0.0_f64
      unew_0003 = 0.0_f64
      vnew_0003 = 0.0_f64
      B1new_0003 = 0.0_f64
      B2new_0003 = 0.0_f64
      znew_0003 = 0.0_f64
      do i_0003 = 0_i64, 5_i64, 1_i64
        hnew_0003 = hnew_0003 + smmat_0003(i_0003, 0_i64) * w_dif_0003( &
              i_0003)
        unew_0003 = unew_0003 + smmat_0003(i_0003, 1_i64) * w_dif_0003( &
              i_0003)
        vnew_0003 = vnew_0003 + smmat_0003(i_0003, 2_i64) * w_dif_0003( &
              i_0003)
        B1new_0003 = B1new_0003 + smmat_0003(i_0003, 3_i64) * w_dif_0003 &
              (i_0003)
        B2new_0003 = B2new_0003 + smmat_0003(i_0003, 4_i64) * w_dif_0003 &
              (i_0003)
        znew_0003 = znew_0003 + smmat_0003(i_0003, 5_i64) * w_dif_0003( &
              i_0003)
      end do
      Pnew_0003 = cpsi * (B1right_0003 - B1left_0003)
      u_h_0003 = hnew_0003 / 2_i64
      u_hu_0003 = unew_0003 / 2_i64
      u_hv_0003 = vnew_0003 / 2_i64
      u_hP_0003 = Pnew_0003 / 2_i64
      u_hB1_0003 = B1new_0003 / 2_i64
      u_hB2_0003 = B2new_0003 / 2_i64
      u_z_0003 = znew_0003 / 2_i64
      w_lrh_0003 = w_lrh_0003 - u_h_0003
      w_lrhu_0003 = w_lrhu_0003 - u_hu_0003
      w_lrhv_0003 = w_lrhv_0003 - u_hv_0003
      w_lrhP_0003 = w_lrhPSI_0003 - u_hP_0003
      w_lrhB1_0003 = w_lrhB1_0003 - u_hB1_0003
      w_lrhB2_0003 = w_lrhB2_0003 - u_hB2_0003
      w_lrz_0003 = w_lrz_0003 - u_z_0003
      w_hP_0003 = hPSI_r - hPSI_l
      mw_hB1_0003 = w_lrhB1_0003 / 2_i64 - 1_i64 / (2_i64 * cpsi) * &
            w_hP_0003
      mhP_0003 = (hPSI_r - hPSI_l) / 2_i64 - cpsi * w_dif_0003(3_i64) / &
            2_i64
      unew_0003 = 0.0_f64
      vnew_0003 = 0.0_f64
      B1new_0003 = 0.0_f64
      B2new_0003 = 0.0_f64
      unew_0003 = w_lrhu_0003 * normal(0_i64) - w_lrhv_0003 * normal( &
            1_i64)
      unew_0003 = unew_0003 / mesure
      vnew_0003 = w_lrhu_0003 * normal(1_i64) + w_lrhv_0003 * normal( &
            0_i64)
      vnew_0003 = vnew_0003 / mesure
      B1new_0003 = w_lrhB1_0003 * normal(0_i64) - w_lrhB2_0003 * normal( &
            1_i64)
      B1new_0003 = B1new_0003 / mesure
      B2new_0003 = w_lrhB1_0003 * normal(1_i64) + w_lrhB2_0003 * normal( &
            0_i64)
      B2new_0003 = B2new_0003 / mesure
      w_lrhu_0003 = unew_0003
      w_lrhv_0003 = vnew_0003
      w_lrhB1_0003 = B1new_0003
      w_lrhB2_0003 = B2new_0003
      allocate(norm_0003(0:size(normal, kind=i64) - 1_i64))
      norm_0003 = normal / mesure
      q_s_0003 = normal(0_i64) * unew_0003 + normal(1_i64) * vnew_0003
      p_s_0003 = normal(0_i64) * B1new_0003 + normal(1_i64) * B2new_0003
      Flux_B1psi_0003 = mhP_0003 * norm_0003(0_i64) * mesure
      Flux_B2psi_0003 = mhP_0003 * norm_0003(1_i64) * mesure
      Flux_hPpsi_0003 = cpsi * cpsi * mw_hB1_0003 * mesure
      flux(0_i64) = q_s_0003
      flux(1_i64) = q_s_0003 * w_lrhu_0003 / w_lrh_0003 + 0.5_f64 * grav &
            * w_lrh_0003 * w_lrh_0003 * normal(0_i64) - p_s_0003 * &
            w_lrhB1_0003 / w_lrh_0003
      flux(2_i64) = q_s_0003 * w_lrhv_0003 / w_lrh_0003 + 0.5_f64 * grav &
            * w_lrh_0003 * w_lrh_0003 * normal(1_i64) - p_s_0003 * &
            w_lrhB2_0003 / w_lrh_0003
      flux(3_i64) = (w_lrhv_0003 * w_lrhB1_0003 / w_lrh_0003 - &
            w_lrhu_0003 * w_lrhB2_0003 / w_lrh_0003) * normal(1_i64) + &
            Flux_B1psi_0003
      flux(4_i64) = (w_lrhu_0003 * w_lrhB2_0003 / w_lrh_0003 - &
            w_lrhv_0003 * w_lrhB1_0003 / w_lrh_0003) * normal(0_i64) + &
            Flux_B2psi_0003
      flux(5_i64) = Flux_hPpsi_0003
      flux(6_i64) = 0.0_f64
      if (allocated(ninv_0003)) then
        deallocate(ninv_0003)
      end if
      if (allocated(w_dif_0003)) then
        deallocate(w_dif_0003)
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
    do Dummy_0008 = 0_i64, size(boundaryfaces, kind=i64) - 1_i64, 1_i64
      i = boundaryfaces(Dummy_0008)
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
      !center_left = centerc[cellidf[i][0]]
      !center_right = centerg[i]
      !
      !h_x_left = h_x[cellidf[i][0]];   h_y_left = h_y[cellidf[i][0]];
      !
      !
      !psi_left = psi[cellidf[i][0]];
      !
      !r_l[0] = centerf[i][0] - center_left[0]; r_r[0] = centerf[i][0] - center_right[0];
      !r_l[1] = centerf[i][1] - center_left[1]; r_r[1] = centerf[i][1] - center_right[1];
      !h_l  = h_l  + (order - 1) * psi_left  * (h_x_left * r_l[0] + h_y_left * r_l[1] )
      !h_r  = h_r
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
      B1_h_0004 = (hB1_l / h_l * sqrt(h_l) + hB1_r / h_r * sqrt(h_r)) / &
            (sqrt(h_l) + sqrt(h_r))
      B2_h_0004 = (hB2_l / h_l * sqrt(h_l) + hB2_r / h_r * sqrt(h_r)) / &
            (sqrt(h_l) + sqrt(h_r))
      !uvh =  array([uh, vh])
      un_h_0004 = u_h_0004 * normal(0_i64) + v_h_0004 * normal(1_i64)
      un_h_0004 = un_h_0004 / mesure
      vn_h_0004 = u_h_0004 * ninv_0004(0_i64) + v_h_0004 * ninv_0004( &
            1_i64)
      vn_h_0004 = vn_h_0004 / mesure
      B1n_h_0004 = B1_h_0004 * normal(0_i64) + B2_h_0004 * normal(1_i64)
      B1n_h_0004 = B1n_h_0004 / mesure
      B2n_h_0004 = B1_h_0004 * ninv_0004(0_i64) + B2_h_0004 * ninv_0004( &
            1_i64)
      B2n_h_0004 = B2n_h_0004 / mesure
      hroe_0004 = (h_l + h_r) / 2_i64
      uroe_0004 = un_h_0004
      vroe_0004 = vn_h_0004
      B1roe_0004 = B1n_h_0004
      B2roe_0004 = B2n_h_0004
      uleft_0004 = hu_l * normal(0_i64) + hv_l * normal(1_i64)
      uleft_0004 = uleft_0004 / mesure
      vleft_0004 = hu_l * ninv_0004(0_i64) + hv_l * ninv_0004(1_i64)
      vleft_0004 = vleft_0004 / mesure
      B1left_0004 = hB1_l * normal(0_i64) + hB2_l * normal(1_i64)
      B1left_0004 = B1left_0004 / mesure
      B2left_0004 = hB1_l * ninv_0004(0_i64) + hB2_l * ninv_0004(1_i64)
      B2left_0004 = B2left_0004 / mesure
      uright_0004 = hu_r * normal(0_i64) + hv_r * normal(1_i64)
      uright_0004 = uright_0004 / mesure
      vright_0004 = hu_r * ninv_0004(0_i64) + hv_r * ninv_0004(1_i64)
      vright_0004 = vright_0004 / mesure
      B1right_0004 = hB1_r * normal(0_i64) + hB2_r * normal(1_i64)
      B1right_0004 = B1right_0004 / mesure
      B2right_0004 = hB1_r * ninv_0004(0_i64) + hB2_r * ninv_0004(1_i64)
      B2right_0004 = B2right_0004 / mesure
      B1i_0004 = (hB1c_l * normal(0_i64) + hB2c_l * normal(1_i64)) / &
            mesure
      B1j_0004 = (hB1c_r * normal(0_i64) + hB2c_r * normal(1_i64)) / &
            mesure
      absy_0004 = (B1i_0004 + B1j_0004) / 2_i64
      w_lrh_0004 = (h_l + h_r) / 2_i64
      w_lrhu_0004 = (uleft_0004 + uright_0004) / 2_i64
      w_lrhv_0004 = (vleft_0004 + vright_0004) / 2_i64
      w_lrhB1_0004 = (B1left_0004 + B1right_0004) / 2_i64
      w_lrhB2_0004 = (B2left_0004 + B2right_0004) / 2_i64
      w_lrhPSI_0004 = (hPSI_l + hPSI_r) / 2_i64
      w_lrz_0004 = (Z_l + Z_r) / 2_i64
      w_dif_0004(0_i64) = h_r - h_l
      w_dif_0004(1_i64) = uright_0004 - uleft_0004
      w_dif_0004(2_i64) = vright_0004 - vleft_0004
      w_dif_0004(3_i64) = B1right_0004 - B1left_0004
      w_dif_0004(4_i64) = B2right_0004 - B2left_0004
      w_dif_0004(5_i64) = Z_r - Z_l
      allocate(signA_0004(0:5_i64, 0:5_i64))
      signA_0004 = 0.0_f64
      sound_0004 = sqrt(grav * hroe_0004)
      !B1roe = absy/hroe
      w_0004 = sqrt(B1roe_0004 * B1roe_0004 + grav * hroe_0004)
      lambda1_0004 = uroe_0004 - w_0004
      lambda2_0004 = uroe_0004 - B1roe_0004
      lambda3_0004 = uroe_0004 + B1roe_0004
      lambda4_0004 = uroe_0004 + w_0004
      epsilon_0004 = 1e-15_f64
      if (abs(lambda1_0004) < epsilon_0004) then
        s1_0004 = 0.0_f64
        pi1_0004 = 0.0_f64
      else
        s1_0004 = lambda1_0004 / abs(lambda1_0004)
        pi1_0004 = s1_0004 / lambda1_0004
      end if
      if (abs(lambda2_0004) < epsilon_0004) then
        s2_0004 = 0.0_f64
        pi2_0004 = 0.0_f64
      else
        s2_0004 = lambda2_0004 / abs(lambda2_0004)
        pi2_0004 = 1.0_f64 / abs(lambda2_0004)
      end if
      if (abs(lambda3_0004) < epsilon_0004) then
        s3_0004 = 0.0_f64
        pi3_0004 = 0.0_f64
      else
        s3_0004 = lambda3_0004 / abs(lambda3_0004)
        pi3_0004 = 1.0_f64 / abs(lambda3_0004)
      end if
      if (abs(lambda4_0004) < epsilon_0004) then
        s4_0004 = 0.0_f64
        pi4_0004 = 0.0_f64
      else
        s4_0004 = lambda4_0004 / abs(lambda4_0004)
        pi4_0004 = 1.0_f64 / abs(lambda4_0004)
      end if
      gamma1_0004 = vroe_0004 + B2roe_0004
      gamma2_0004 = vroe_0004 - B2roe_0004
      sigma1_0004 = vroe_0004 * (s1_0004 * lambda4_0004 - s4_0004 * &
            lambda1_0004) - w_0004 * (s2_0004 * gamma1_0004 + s3_0004 * &
            gamma2_0004)
      sigma2_0004 = B2roe_0004 * (s1_0004 * lambda4_0004 - s4_0004 * &
            lambda1_0004) - w_0004 * (s2_0004 * gamma1_0004 - s3_0004 * &
            gamma2_0004)
      if (abs(lambda2_0004) < epsilon_0004 .and. abs(lambda3_0004) < &
            epsilon_0004) then
        mu1_0004 = B1roe_0004 * vroe_0004 * pi1_0004 / w_0004 - &
              B1roe_0004 * vroe_0004 * pi4_0004 / w_0004
        mu2_0004 = B1roe_0004 * B2roe_0004 * pi1_0004 / w_0004 - &
              B1roe_0004 * B2roe_0004 * pi4_0004 / w_0004
        ann_0004 = 1_i64
      else
        mu1_0004 = B1roe_0004 * vroe_0004 * pi1_0004 / w_0004 - &
              B1roe_0004 * vroe_0004 * pi4_0004 / w_0004 - 0.5_f64 * ( &
              gamma1_0004 * pi2_0004 - gamma2_0004 * pi3_0004)
        mu2_0004 = B1roe_0004 * B2roe_0004 * pi1_0004 / w_0004 - &
              B1roe_0004 * B2roe_0004 * pi4_0004 / w_0004 - 0.5_f64 * ( &
              gamma1_0004 * pi2_0004 + gamma2_0004 * pi3_0004)
        ann_0004 = 1_i64
      end if
      !1ère colonne de la matrice A
      signA_0004(0_i64, 0_i64) = (s1_0004 * lambda4_0004 - s4_0004 * &
            lambda1_0004) / (2_i64 * w_0004)
      signA_0004(0_i64, 1_i64) = lambda1_0004 * lambda4_0004 * (s1_0004 &
            - s4_0004) / (2_i64 * w_0004)
      signA_0004(0_i64, 2_i64) = sigma1_0004 / (2_i64 * w_0004)
      signA_0004(0_i64, 3_i64) = 0.0_f64
      signA_0004(0_i64, 4_i64) = sigma2_0004 / (2_i64 * w_0004)
      signA_0004(0_i64, 5_i64) = 0.0_f64
      !2ème colonne de la matrice A
      signA_0004(1_i64, 0_i64) = (s4_0004 - s1_0004) / (2_i64 * w_0004)
      signA_0004(1_i64, 1_i64) = (s4_0004 * lambda4_0004 - s1_0004 * &
            lambda1_0004) / (2_i64 * w_0004)
      signA_0004(1_i64, 2_i64) = vroe_0004 * (s4_0004 - s1_0004) / ( &
            2_i64 * w_0004)
      signA_0004(1_i64, 3_i64) = 0.0_f64
      signA_0004(1_i64, 4_i64) = B2roe_0004 * (s4_0004 - s1_0004) / ( &
            2_i64 * w_0004)
      signA_0004(1_i64, 5_i64) = 0.0_f64
      !3ème colonne de la matrice A
      signA_0004(2_i64, 0_i64) = 0.0_f64
      signA_0004(2_i64, 1_i64) = 0.0_f64
      signA_0004(2_i64, 2_i64) = (s2_0004 + s3_0004) / 2_i64
      signA_0004(2_i64, 3_i64) = 0.0_f64
      signA_0004(2_i64, 4_i64) = (s2_0004 - s3_0004) / 2_i64
      signA_0004(2_i64, 5_i64) = 0.0_f64
      !4ème colonne de la matrice A
      signA_0004(3_i64, 0_i64) = ann_0004 * B1roe_0004 * (pi1_0004 - &
            pi4_0004) / w_0004
      signA_0004(3_i64, 1_i64) = ann_0004 * B1roe_0004 * (s1_0004 - &
            s4_0004) / w_0004
      signA_0004(3_i64, 2_i64) = ann_0004 * mu1_0004
      signA_0004(3_i64, 3_i64) = 0.0_f64
      signA_0004(3_i64, 4_i64) = ann_0004 * mu2_0004
      signA_0004(3_i64, 5_i64) = 0.0_f64
      !5ème colonne de la matrice A
      signA_0004(4_i64, 0_i64) = 0.0_f64
      signA_0004(4_i64, 1_i64) = 0.0_f64
      signA_0004(4_i64, 2_i64) = (s2_0004 - s3_0004) / 2_i64
      signA_0004(4_i64, 3_i64) = 0.0_f64
      signA_0004(4_i64, 4_i64) = (s2_0004 + s3_0004) / 2_i64
      signA_0004(4_i64, 5_i64) = 0.0_f64
      !6ème colonne de la matrice A
      signA_0004(5_i64, 0_i64) = sound_0004 ** 2_i64 * (pi4_0004 - &
            pi1_0004) / (2_i64 * w_0004)
      signA_0004(5_i64, 1_i64) = sound_0004 ** 2_i64 * (s4_0004 - &
            s1_0004) / (2_i64 * w_0004)
      signA_0004(5_i64, 2_i64) = sound_0004 ** 2_i64 * vroe_0004 * ( &
            pi4_0004 - pi1_0004) / (2_i64 * w_0004)
      signA_0004(5_i64, 3_i64) = 0.0_f64
      signA_0004(5_i64, 4_i64) = sound_0004 ** 2_i64 * B2roe_0004 * ( &
            pi4_0004 - pi1_0004) / (2_i64 * w_0004)
      signA_0004(5_i64, 5_i64) = 0.0_f64
      smmat_0004(0:, 0:) => signA_0004
      hnew_0004 = 0.0_f64
      unew_0004 = 0.0_f64
      vnew_0004 = 0.0_f64
      B1new_0004 = 0.0_f64
      B2new_0004 = 0.0_f64
      znew_0004 = 0.0_f64
      do i_0004 = 0_i64, 5_i64, 1_i64
        hnew_0004 = hnew_0004 + smmat_0004(i_0004, 0_i64) * w_dif_0004( &
              i_0004)
        unew_0004 = unew_0004 + smmat_0004(i_0004, 1_i64) * w_dif_0004( &
              i_0004)
        vnew_0004 = vnew_0004 + smmat_0004(i_0004, 2_i64) * w_dif_0004( &
              i_0004)
        B1new_0004 = B1new_0004 + smmat_0004(i_0004, 3_i64) * w_dif_0004 &
              (i_0004)
        B2new_0004 = B2new_0004 + smmat_0004(i_0004, 4_i64) * w_dif_0004 &
              (i_0004)
        znew_0004 = znew_0004 + smmat_0004(i_0004, 5_i64) * w_dif_0004( &
              i_0004)
      end do
      Pnew_0004 = cpsi * (B1right_0004 - B1left_0004)
      u_h_0004 = hnew_0004 / 2_i64
      u_hu_0004 = unew_0004 / 2_i64
      u_hv_0004 = vnew_0004 / 2_i64
      u_hP_0004 = Pnew_0004 / 2_i64
      u_hB1_0004 = B1new_0004 / 2_i64
      u_hB2_0004 = B2new_0004 / 2_i64
      u_z_0004 = znew_0004 / 2_i64
      w_lrh_0004 = w_lrh_0004 - u_h_0004
      w_lrhu_0004 = w_lrhu_0004 - u_hu_0004
      w_lrhv_0004 = w_lrhv_0004 - u_hv_0004
      w_lrhP_0004 = w_lrhPSI_0004 - u_hP_0004
      w_lrhB1_0004 = w_lrhB1_0004 - u_hB1_0004
      w_lrhB2_0004 = w_lrhB2_0004 - u_hB2_0004
      w_lrz_0004 = w_lrz_0004 - u_z_0004
      w_hP_0004 = hPSI_r - hPSI_l
      mw_hB1_0004 = w_lrhB1_0004 / 2_i64 - 1_i64 / (2_i64 * cpsi) * &
            w_hP_0004
      mhP_0004 = (hPSI_r - hPSI_l) / 2_i64 - cpsi * w_dif_0004(3_i64) / &
            2_i64
      unew_0004 = 0.0_f64
      vnew_0004 = 0.0_f64
      B1new_0004 = 0.0_f64
      B2new_0004 = 0.0_f64
      unew_0004 = w_lrhu_0004 * normal(0_i64) - w_lrhv_0004 * normal( &
            1_i64)
      unew_0004 = unew_0004 / mesure
      vnew_0004 = w_lrhu_0004 * normal(1_i64) + w_lrhv_0004 * normal( &
            0_i64)
      vnew_0004 = vnew_0004 / mesure
      B1new_0004 = w_lrhB1_0004 * normal(0_i64) - w_lrhB2_0004 * normal( &
            1_i64)
      B1new_0004 = B1new_0004 / mesure
      B2new_0004 = w_lrhB1_0004 * normal(1_i64) + w_lrhB2_0004 * normal( &
            0_i64)
      B2new_0004 = B2new_0004 / mesure
      w_lrhu_0004 = unew_0004
      w_lrhv_0004 = vnew_0004
      w_lrhB1_0004 = B1new_0004
      w_lrhB2_0004 = B2new_0004
      allocate(norm_0004(0:size(normal, kind=i64) - 1_i64))
      norm_0004 = normal / mesure
      q_s_0004 = normal(0_i64) * unew_0004 + normal(1_i64) * vnew_0004
      p_s_0004 = normal(0_i64) * B1new_0004 + normal(1_i64) * B2new_0004
      Flux_B1psi_0004 = mhP_0004 * norm_0004(0_i64) * mesure
      Flux_B2psi_0004 = mhP_0004 * norm_0004(1_i64) * mesure
      Flux_hPpsi_0004 = cpsi * cpsi * mw_hB1_0004 * mesure
      flux(0_i64) = q_s_0004
      flux(1_i64) = q_s_0004 * w_lrhu_0004 / w_lrh_0004 + 0.5_f64 * grav &
            * w_lrh_0004 * w_lrh_0004 * normal(0_i64) - p_s_0004 * &
            w_lrhB1_0004 / w_lrh_0004
      flux(2_i64) = q_s_0004 * w_lrhv_0004 / w_lrh_0004 + 0.5_f64 * grav &
            * w_lrh_0004 * w_lrh_0004 * normal(1_i64) - p_s_0004 * &
            w_lrhB2_0004 / w_lrh_0004
      flux(3_i64) = (w_lrhv_0004 * w_lrhB1_0004 / w_lrh_0004 - &
            w_lrhu_0004 * w_lrhB2_0004 / w_lrh_0004) * normal(1_i64) + &
            Flux_B1psi_0004
      flux(4_i64) = (w_lrhu_0004 * w_lrhB2_0004 / w_lrh_0004 - &
            w_lrhv_0004 * w_lrhB1_0004 / w_lrh_0004) * normal(0_i64) + &
            Flux_B2psi_0004
      flux(5_i64) = Flux_hPpsi_0004
      flux(6_i64) = 0.0_f64
      if (allocated(ninv_0004)) then
        deallocate(ninv_0004)
      end if
      if (allocated(w_dif_0004)) then
        deallocate(w_dif_0004)
      end if
      if (allocated(signA_0004)) then
        deallocate(signA_0004)
      end if
      if (allocated(norm_0004)) then
        deallocate(norm_0004)
      end if
      !LF_scheme_MHD(hu_l, hu_r, hv_l, hv_r, h_l, h_r, hB1_l, hB1_r, hB2_l, hB2_r, Z_l, Z_r, normal, mesure, grav, flux)
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

  end subroutine explicitscheme_convective_SWMHD
  !........................................

end module tools_SWMHD
