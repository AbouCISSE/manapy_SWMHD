module bind_c_tools_SWMHD

  use tools_SWMHD, only: term_source_srnh_SWMHD
  use tools_SWMHD, only: cpsi_global
  use tools_SWMHD, only: update_SWMHD
  use tools_SWMHD, only: explicitscheme_stabilizater
  use tools_SWMHD, only: time_step_SWMHD
  use tools_SWMHD, only: explicitscheme_convective_SWMHD
  use tools_SWMHD, only: initialisation_SWMHD
  use tools_SWMHD, only: Total_Energy

  use, intrinsic :: ISO_C_Binding, only : f64 => C_DOUBLE , i64 => &
        C_INT64_T
  implicit none

  contains

  !........................................
  subroutine bind_c_total_energy(n0_h_c, h_c, n0_hu_c, hu_c, n0_hv_c, &
        hv_c, n0_hB1_c, hB1_c, n0_hB2_c, hB2_c, n0_Z_c, Z_c, grav, &
        n0_volumec, volumec, numt, numc, numm, nump) bind(c)

    implicit none

    integer(i64), value :: n0_h_c
    real(f64), intent(in) :: h_c(0:n0_h_c - 1_i64)
    integer(i64), value :: n0_hu_c
    real(f64), intent(in) :: hu_c(0:n0_hu_c - 1_i64)
    integer(i64), value :: n0_hv_c
    real(f64), intent(in) :: hv_c(0:n0_hv_c - 1_i64)
    integer(i64), value :: n0_hB1_c
    real(f64), intent(in) :: hB1_c(0:n0_hB1_c - 1_i64)
    integer(i64), value :: n0_hB2_c
    real(f64), intent(in) :: hB2_c(0:n0_hB2_c - 1_i64)
    integer(i64), value :: n0_Z_c
    real(f64), intent(in) :: Z_c(0:n0_Z_c - 1_i64)
    real(f64), value :: grav
    integer(i64), value :: n0_volumec
    real(f64), intent(in) :: volumec(0:n0_volumec - 1_i64)
    real(f64), intent(out) :: numt
    real(f64), intent(out) :: numc
    real(f64), intent(out) :: numm
    real(f64), intent(out) :: nump

    call Total_Energy(h_c, hu_c, hv_c, hB1_c, hB2_c, Z_c, grav, volumec, &
          numt = numt, numc = numc, numm = numm, nump = nump)

  end subroutine bind_c_total_energy
  !........................................

  !........................................
  subroutine bind_c_explicitscheme_stabilizater(n0_wx_face, wx_face, &
        n0_wy_face, wy_face, n0_cellidf, n1_cellidf, cellidf, &
        n0_normalf, n1_normalf, normalf, n0_namef, namef, vepsilon, &
        n0_dissip_w, dissip_w, n0_h_c, h_c) bind(c)

    implicit none

    integer(i64), value :: n0_wx_face
    real(f64), intent(in) :: wx_face(0:n0_wx_face - 1_i64)
    integer(i64), value :: n0_wy_face
    real(f64), intent(in) :: wy_face(0:n0_wy_face - 1_i64)
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
    real(f64), value :: vepsilon
    integer(i64), value :: n0_dissip_w
    real(f64), intent(inout) :: dissip_w(0:n0_dissip_w - 1_i64)
    integer(i64), value :: n0_h_c
    real(f64), intent(in) :: h_c(0:n0_h_c - 1_i64)

    call explicitscheme_stabilizater(wx_face, wy_face, cellidf, normalf, &
          namef, vepsilon, dissip_w, h_c)

  end subroutine bind_c_explicitscheme_stabilizater
  !........................................

  !........................................
  function bind_c_time_step_swmhd(n0_h_c, h_c, n0_hu_c, hu_c, n0_hv_c, &
        hv_c, n0_hB1_c, hB1_c, n0_hB2_c, hB2_c, cfl, n0_normal, &
        n1_normal, normal, n0_mesure, mesure, n0_volume, volume, &
        n0_faceid, n1_faceid, faceid) bind(c) result(dt)

    implicit none

    integer(i64), value :: n0_h_c
    real(f64), intent(in) :: h_c(0:n0_h_c - 1_i64)
    integer(i64), value :: n0_hu_c
    real(f64), intent(in) :: hu_c(0:n0_hu_c - 1_i64)
    integer(i64), value :: n0_hv_c
    real(f64), intent(in) :: hv_c(0:n0_hv_c - 1_i64)
    integer(i64), value :: n0_hB1_c
    real(f64), intent(in) :: hB1_c(0:n0_hB1_c - 1_i64)
    integer(i64), value :: n0_hB2_c
    real(f64), intent(in) :: hB2_c(0:n0_hB2_c - 1_i64)
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
    real(f64) :: dt

    dt = time_step_SWMHD(h_c, hu_c, hv_c, hB1_c, hB2_c, cfl, normal, &
          mesure, volume, faceid)

  end function bind_c_time_step_swmhd
  !........................................

  !........................................
  subroutine bind_c_update_swmhd(n0_h_c, h_c, n0_hu_c, hu_c, n0_hv_c, &
        hv_c, n0_hB1_c, hB1_c, n0_hB2_c, hB2_c, n0_PSI_c, PSI_c, n0_Z_c &
        , Z_c, n0_rez_h, rez_h, n0_rez_hu, rez_hu, n0_rez_hv, rez_hv, &
        n0_rez_hB1, rez_hB1, n0_rez_hB2, rez_hB2, n0_rez_PSI, rez_PSI, &
        n0_rez_Z, rez_Z, n0_src_h, src_h, n0_src_hu, src_hu, n0_src_hv, &
        src_hv, n0_src_hB1, src_hB1, n0_src_hB2, src_hB2, n0_src_PSI, &
        src_PSI, n0_src_Z, src_Z, n0_corio_hu, corio_hu, n0_corio_hv, &
        corio_hv, dtime, n0_vol, vol, GLM, cpsi) bind(c)

    implicit none

    integer(i64), value :: n0_h_c
    real(f64), intent(inout) :: h_c(0:n0_h_c - 1_i64)
    integer(i64), value :: n0_hu_c
    real(f64), intent(inout) :: hu_c(0:n0_hu_c - 1_i64)
    integer(i64), value :: n0_hv_c
    real(f64), intent(inout) :: hv_c(0:n0_hv_c - 1_i64)
    integer(i64), value :: n0_hB1_c
    real(f64), intent(inout) :: hB1_c(0:n0_hB1_c - 1_i64)
    integer(i64), value :: n0_hB2_c
    real(f64), intent(inout) :: hB2_c(0:n0_hB2_c - 1_i64)
    integer(i64), value :: n0_PSI_c
    real(f64), intent(inout) :: PSI_c(0:n0_PSI_c - 1_i64)
    integer(i64), value :: n0_Z_c
    real(f64), intent(inout) :: Z_c(0:n0_Z_c - 1_i64)
    integer(i64), value :: n0_rez_h
    real(f64), intent(in) :: rez_h(0:n0_rez_h - 1_i64)
    integer(i64), value :: n0_rez_hu
    real(f64), intent(in) :: rez_hu(0:n0_rez_hu - 1_i64)
    integer(i64), value :: n0_rez_hv
    real(f64), intent(in) :: rez_hv(0:n0_rez_hv - 1_i64)
    integer(i64), value :: n0_rez_hB1
    real(f64), intent(in) :: rez_hB1(0:n0_rez_hB1 - 1_i64)
    integer(i64), value :: n0_rez_hB2
    real(f64), intent(in) :: rez_hB2(0:n0_rez_hB2 - 1_i64)
    integer(i64), value :: n0_rez_PSI
    real(f64), intent(in) :: rez_PSI(0:n0_rez_PSI - 1_i64)
    integer(i64), value :: n0_rez_Z
    real(f64), intent(in) :: rez_Z(0:n0_rez_Z - 1_i64)
    integer(i64), value :: n0_src_h
    real(f64), intent(in) :: src_h(0:n0_src_h - 1_i64)
    integer(i64), value :: n0_src_hu
    real(f64), intent(in) :: src_hu(0:n0_src_hu - 1_i64)
    integer(i64), value :: n0_src_hv
    real(f64), intent(in) :: src_hv(0:n0_src_hv - 1_i64)
    integer(i64), value :: n0_src_hB1
    real(f64), intent(in) :: src_hB1(0:n0_src_hB1 - 1_i64)
    integer(i64), value :: n0_src_hB2
    real(f64), intent(in) :: src_hB2(0:n0_src_hB2 - 1_i64)
    integer(i64), value :: n0_src_PSI
    real(f64), intent(in) :: src_PSI(0:n0_src_PSI - 1_i64)
    integer(i64), value :: n0_src_Z
    real(f64), intent(in) :: src_Z(0:n0_src_Z - 1_i64)
    integer(i64), value :: n0_corio_hu
    real(f64), intent(in) :: corio_hu(0:n0_corio_hu - 1_i64)
    integer(i64), value :: n0_corio_hv
    real(f64), intent(in) :: corio_hv(0:n0_corio_hv - 1_i64)
    real(f64), value :: dtime
    integer(i64), value :: n0_vol
    real(f64), intent(in) :: vol(0:n0_vol - 1_i64)
    integer(i64), value :: GLM
    real(f64), value :: cpsi

    call update_SWMHD(h_c, hu_c, hv_c, hB1_c, hB2_c, PSI_c, Z_c, rez_h, &
          rez_hu, rez_hv, rez_hB1, rez_hB2, rez_PSI, rez_Z, src_h, &
          src_hu, src_hv, src_hB1, src_hB2, src_PSI, src_Z, corio_hu, &
          corio_hv, dtime, vol, GLM, cpsi)

  end subroutine bind_c_update_swmhd
  !........................................

  !........................................
  function bind_c_cpsi_global(n0_h_c, h_c, n0_hu_c, hu_c, n0_hv_c, hv_c, &
        n0_hB1_c, hB1_c, n0_hB2_c, hB2_c, cfl, n0_normal, n1_normal, &
        normal, n0_mesure, mesure, n0_volume, volume, n0_faceid, &
        n1_faceid, faceid) bind(c) result(cpsiglobal)

    implicit none

    integer(i64), value :: n0_h_c
    real(f64), intent(in) :: h_c(0:n0_h_c - 1_i64)
    integer(i64), value :: n0_hu_c
    real(f64), intent(in) :: hu_c(0:n0_hu_c - 1_i64)
    integer(i64), value :: n0_hv_c
    real(f64), intent(in) :: hv_c(0:n0_hv_c - 1_i64)
    integer(i64), value :: n0_hB1_c
    real(f64), intent(in) :: hB1_c(0:n0_hB1_c - 1_i64)
    integer(i64), value :: n0_hB2_c
    real(f64), intent(in) :: hB2_c(0:n0_hB2_c - 1_i64)
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
    real(f64) :: cpsiglobal

    cpsiglobal = cpsi_global(h_c, hu_c, hv_c, hB1_c, hB2_c, cfl, normal, &
          mesure, volume, faceid)

  end function bind_c_cpsi_global
  !........................................

  !........................................
  subroutine bind_c_initialisation_swmhd(n0_h, h, n0_hu, hu, n0_hv, hv, &
        n0_hB1, hB1, n0_hB2, hB2, n0_PSI, PSI, n0_Z, Z, n0_center, &
        n1_center, center, choix, k1, k2, eps, tol) bind(c)

    implicit none

    integer(i64), value :: n0_h
    real(f64), intent(inout) :: h(0:n0_h - 1_i64)
    integer(i64), value :: n0_hu
    real(f64), intent(inout) :: hu(0:n0_hu - 1_i64)
    integer(i64), value :: n0_hv
    real(f64), intent(inout) :: hv(0:n0_hv - 1_i64)
    integer(i64), value :: n0_hB1
    real(f64), intent(inout) :: hB1(0:n0_hB1 - 1_i64)
    integer(i64), value :: n0_hB2
    real(f64), intent(inout) :: hB2(0:n0_hB2 - 1_i64)
    integer(i64), value :: n0_PSI
    real(f64), intent(inout) :: PSI(0:n0_PSI - 1_i64)
    integer(i64), value :: n0_Z
    real(f64), intent(inout) :: Z(0:n0_Z - 1_i64)
    integer(i64), value :: n0_center
    integer(i64), value :: n1_center
    real(f64), intent(in) :: center(0:n1_center - 1_i64,0:n0_center - &
          1_i64)
    integer(i64), value :: choix
    real(f64), value :: k1
    real(f64), value :: k2
    real(f64), value :: eps
    real(f64), value :: tol

    call initialisation_SWMHD(h, hu, hv, hB1, hB2, PSI, Z, center, choix &
          , k1, k2, eps, tol)

  end subroutine bind_c_initialisation_swmhd
  !........................................

  !........................................
  subroutine bind_c_term_source_srnh_swmhd(n0_src_h, src_h, n0_src_hu, &
        src_hu, n0_src_hv, src_hv, n0_src_hB1, src_hB1, n0_src_hB2, &
        src_hB2, n0_src_PSI, src_PSI, n0_src_Z, src_Z, n0_h_c, h_c, &
        n0_hu_c, hu_c, n0_hv_c, hv_c, n0_Z_c, Z_c, n0_h_ghost, h_ghost, &
        n0_hu_ghost, hu_ghost, n0_hv_ghost, hv_ghost, n0_Z_ghost, &
        Z_ghost, n0_h_halo, h_halo, n0_hu_halo, hu_halo, n0_hv_halo, &
        hv_halo, n0_Z_halo, Z_halo, n0_h_x, h_x, n0_h_y, h_y, n0_psi, &
        psi, n0_hx_halo, hx_halo, n0_hy_halo, hy_halo, n0_psi_halo, &
        psi_halo, n0_nodeidc, n1_nodeidc, nodeidc, n0_faceidc, &
        n1_faceidc, faceidc, n0_cellidc, n1_cellidc, cellidc, &
        n0_cellidf, n1_cellidf, cellidf, n0_centerc, n1_centerc, &
        centerc, n0_normalc, n1_normalc, n2_normalc, normalc, n0_namef, &
        namef, n0_centerf, n1_centerf, centerf, n0_centerh, n1_centerh, &
        centerh, n0_vertexn, n1_vertexn, vertexn, n0_halofid, halofid, &
        order) bind(c)

    implicit none

    integer(i64), value :: n0_src_h
    real(f64), intent(inout) :: src_h(0:n0_src_h - 1_i64)
    integer(i64), value :: n0_src_hu
    real(f64), intent(inout) :: src_hu(0:n0_src_hu - 1_i64)
    integer(i64), value :: n0_src_hv
    real(f64), intent(inout) :: src_hv(0:n0_src_hv - 1_i64)
    integer(i64), value :: n0_src_hB1
    real(f64), intent(inout) :: src_hB1(0:n0_src_hB1 - 1_i64)
    integer(i64), value :: n0_src_hB2
    real(f64), intent(inout) :: src_hB2(0:n0_src_hB2 - 1_i64)
    integer(i64), value :: n0_src_PSI
    real(f64), intent(inout) :: src_PSI(0:n0_src_PSI - 1_i64)
    integer(i64), value :: n0_src_Z
    real(f64), intent(inout) :: src_Z(0:n0_src_Z - 1_i64)
    integer(i64), value :: n0_h_c
    real(f64), intent(in) :: h_c(0:n0_h_c - 1_i64)
    integer(i64), value :: n0_hu_c
    real(f64), intent(in) :: hu_c(0:n0_hu_c - 1_i64)
    integer(i64), value :: n0_hv_c
    real(f64), intent(in) :: hv_c(0:n0_hv_c - 1_i64)
    integer(i64), value :: n0_Z_c
    real(f64), intent(in) :: Z_c(0:n0_Z_c - 1_i64)
    integer(i64), value :: n0_h_ghost
    real(f64), intent(in) :: h_ghost(0:n0_h_ghost - 1_i64)
    integer(i64), value :: n0_hu_ghost
    real(f64), intent(in) :: hu_ghost(0:n0_hu_ghost - 1_i64)
    integer(i64), value :: n0_hv_ghost
    real(f64), intent(in) :: hv_ghost(0:n0_hv_ghost - 1_i64)
    integer(i64), value :: n0_Z_ghost
    real(f64), intent(in) :: Z_ghost(0:n0_Z_ghost - 1_i64)
    integer(i64), value :: n0_h_halo
    real(f64), intent(in) :: h_halo(0:n0_h_halo - 1_i64)
    integer(i64), value :: n0_hu_halo
    real(f64), intent(in) :: hu_halo(0:n0_hu_halo - 1_i64)
    integer(i64), value :: n0_hv_halo
    real(f64), intent(in) :: hv_halo(0:n0_hv_halo - 1_i64)
    integer(i64), value :: n0_Z_halo
    real(f64), intent(in) :: Z_halo(0:n0_Z_halo - 1_i64)
    integer(i64), value :: n0_h_x
    real(f64), intent(in) :: h_x(0:n0_h_x - 1_i64)
    integer(i64), value :: n0_h_y
    real(f64), intent(in) :: h_y(0:n0_h_y - 1_i64)
    integer(i64), value :: n0_psi
    real(f64), intent(in) :: psi(0:n0_psi - 1_i64)
    integer(i64), value :: n0_hx_halo
    real(f64), intent(in) :: hx_halo(0:n0_hx_halo - 1_i64)
    integer(i64), value :: n0_hy_halo
    real(f64), intent(in) :: hy_halo(0:n0_hy_halo - 1_i64)
    integer(i64), value :: n0_psi_halo
    real(f64), intent(in) :: psi_halo(0:n0_psi_halo - 1_i64)
    integer(i64), value :: n0_nodeidc
    integer(i64), value :: n1_nodeidc
    integer(i64), intent(in) :: nodeidc(0:n1_nodeidc - 1_i64,0: &
          n0_nodeidc - 1_i64)
    integer(i64), value :: n0_faceidc
    integer(i64), value :: n1_faceidc
    integer(i64), intent(in) :: faceidc(0:n1_faceidc - 1_i64,0: &
          n0_faceidc - 1_i64)
    integer(i64), value :: n0_cellidc
    integer(i64), value :: n1_cellidc
    integer(i64), intent(in) :: cellidc(0:n1_cellidc - 1_i64,0: &
          n0_cellidc - 1_i64)
    integer(i64), value :: n0_cellidf
    integer(i64), value :: n1_cellidf
    integer(i64), intent(in) :: cellidf(0:n1_cellidf - 1_i64,0: &
          n0_cellidf - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_normalc
    integer(i64), value :: n1_normalc
    integer(i64), value :: n2_normalc
    real(f64), intent(in) :: normalc(0:n2_normalc - 1_i64,0:n1_normalc - &
          1_i64,0:n0_normalc - 1_i64)
    integer(i64), value :: n0_namef
    integer(i64), intent(in) :: namef(0:n0_namef - 1_i64)
    integer(i64), value :: n0_centerf
    integer(i64), value :: n1_centerf
    real(f64), intent(in) :: centerf(0:n1_centerf - 1_i64,0:n0_centerf - &
          1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: n0_vertexn
    integer(i64), value :: n1_vertexn
    real(f64), intent(in) :: vertexn(0:n1_vertexn - 1_i64,0:n0_vertexn - &
          1_i64)
    integer(i64), value :: n0_halofid
    integer(i64), intent(in) :: halofid(0:n0_halofid - 1_i64)
    integer(i64), value :: order

    call term_source_srnh_SWMHD(src_h, src_hu, src_hv, src_hB1, src_hB2, &
          src_PSI, src_Z, h_c, hu_c, hv_c, Z_c, h_ghost, hu_ghost, &
          hv_ghost, Z_ghost, h_halo, hu_halo, hv_halo, Z_halo, h_x, h_y &
          , psi, hx_halo, hy_halo, psi_halo, nodeidc, faceidc, cellidc, &
          cellidf, centerc, normalc, namef, centerf, centerh, vertexn, &
          halofid, order)

  end subroutine bind_c_term_source_srnh_swmhd
  !........................................

  !........................................
  subroutine bind_c_lf_scheme_mhd(hu_l, hu_r, hv_l, hv_r, h_l, h_r, &
        hB1_l, hB1_r, hB2_l, hB2_r, Z_l, Z_r, n0_normal, normal, mesure &
        , grav, n0_flux, flux) bind(c)

    implicit none

    real(f64), value :: hu_l
    real(f64), value :: hu_r
    real(f64), value :: hv_l
    real(f64), value :: hv_r
    real(f64), value :: h_l
    real(f64), value :: h_r
    real(f64), value :: hB1_l
    real(f64), value :: hB1_r
    real(f64), value :: hB2_l
    real(f64), value :: hB2_r
    real(f64), value :: Z_l
    real(f64), value :: Z_r
    integer(i64), value :: n0_normal
    real(f64), intent(in) :: normal(0:n0_normal - 1_i64)
    real(f64), value :: mesure
    real(f64), value :: grav
    integer(i64), value :: n0_flux
    real(f64), intent(inout) :: flux(0:n0_flux - 1_i64)
    real(f64), allocatable :: norm_0005(:)
    real(f64), allocatable :: ninv_0005(:)
    real(f64) :: hl_0001
    real(f64) :: ul_0001
    real(f64) :: vl_0001
    real(f64) :: B1l_0001
    real(f64) :: B2l_0001
    real(f64) :: hr_0001
    real(f64) :: ur_0001
    real(f64) :: vr_0001
    real(f64) :: B1r_0001
    real(f64) :: B2r_0001
    real(f64) :: u_h_0005
    real(f64) :: v_h_0005
    real(f64) :: B1_h_0005
    real(f64) :: B2_h_0005
    real(f64) :: un_h_0005
    real(f64) :: vn_h_0005
    real(f64) :: B1n_h_0005
    real(f64) :: B2n_h_0005
    real(f64) :: hroe_0005
    real(f64) :: uroe_0005
    real(f64) :: vroe_0005
    real(f64) :: B1roe_0005
    real(f64) :: B2roe_0005
    real(f64) :: n1_0001
    real(f64) :: n2_0001
    real(f64) :: Unl_0001
    real(f64) :: Bnl_0001
    real(f64) :: wl_0001
    real(f64) :: Unr_0001
    real(f64) :: Bnr_0001
    real(f64) :: wr_0001
    real(f64) :: lambda1l_0001
    real(f64) :: lambda2l_0001
    real(f64) :: lambda3l_0001
    real(f64) :: lambda4l_0001
    real(f64) :: lambda1r_0001
    real(f64) :: lambda2r_0001
    real(f64) :: lambda3r_0001
    real(f64) :: lambda4r_0001
    real(f64) :: ll_0001
    real(f64) :: lr_0001
    real(f64) :: lambda_star_0001
    real(f64) :: q_l_0001
    real(f64) :: m_l_0001
    real(f64) :: p_l_0001
    real(f64) :: q_r_0001
    real(f64) :: m_r_0001
    real(f64) :: p_r_0001
    real(f64) :: fleft_h_0001
    real(f64) :: fleft_hu_0001
    real(f64) :: fleft_hv_0001
    real(f64) :: fleft_hB1_0001
    real(f64) :: fleft_hB2_0001
    real(f64) :: fright_h_0001
    real(f64) :: fright_hu_0001
    real(f64) :: fright_hv_0001
    real(f64) :: fright_hB1_0001
    real(f64) :: fright_hB2_0001
    real(f64) :: f_h_0001
    real(f64) :: f_hu_0001
    real(f64) :: f_hv_0001
    real(f64) :: f_hB1_0001
    real(f64) :: f_hB2_0001

    !@inline
    !def compute_flux_swmhd_LF(flux, fleft, fright, w_l, w_r, normal, mesure, grav):
    !print("Im Roe scheme")
    allocate(norm_0005(0:size(normal, kind=i64) - 1_i64))
    norm_0005 = normal / mesure
    allocate(ninv_0005(0:1_i64))
    ninv_0005 = 0.0_f64
    ninv_0005(0_i64) = (-1_i64) * norm_0005(1_i64)
    ninv_0005(1_i64) = norm_0005(0_i64)
    hl_0001 = h_l
    ul_0001 = hu_l / h_l
    vl_0001 = hv_l / h_l
    B1l_0001 = hB1_l / h_l
    B2l_0001 = hB2_l / h_l
    hr_0001 = h_r
    ur_0001 = hu_r / h_r
    vr_0001 = hv_r / h_r
    B1r_0001 = hB1_r / h_r
    B2r_0001 = hB2_r / h_r
    u_h_0005 = (hu_l / h_l * sqrt(h_l) + hu_r / h_r * sqrt(h_r)) / (sqrt &
          (h_l) + sqrt(h_r))
    v_h_0005 = (hv_l / h_l * sqrt(h_l) + hv_r / h_r * sqrt(h_r)) / (sqrt &
          (h_l) + sqrt(h_r))
    B1_h_0005 = (hB1_l / h_l * sqrt(h_l) + hB1_r / h_r * sqrt(h_r)) / ( &
          sqrt(h_l) + sqrt(h_r))
    B2_h_0005 = (hB2_l / h_l * sqrt(h_l) + hB2_r / h_r * sqrt(h_r)) / ( &
          sqrt(h_l) + sqrt(h_r))
    un_h_0005 = u_h_0005 * normal(0_i64) + v_h_0005 * normal(1_i64)
    un_h_0005 = un_h_0005 / mesure
    vn_h_0005 = u_h_0005 * ninv_0005(0_i64) + v_h_0005 * ninv_0005(1_i64 &
          )
    vn_h_0005 = vn_h_0005 / mesure
    B1n_h_0005 = B1_h_0005 * normal(0_i64) + B2_h_0005 * normal(1_i64)
    B1n_h_0005 = B1n_h_0005 / mesure
    B2n_h_0005 = B1_h_0005 * ninv_0005(0_i64) + B2_h_0005 * ninv_0005( &
          1_i64)
    B2n_h_0005 = B2n_h_0005 / mesure
    hroe_0005 = (h_l + h_r) / 2_i64
    uroe_0005 = un_h_0005
    vroe_0005 = vn_h_0005
    B1roe_0005 = B1n_h_0005
    B2roe_0005 = B2n_h_0005
    n1_0001 = norm_0005(0_i64)
    n2_0001 = norm_0005(1_i64)
    Unl_0001 = ul_0001 * n1_0001 + vl_0001 * n2_0001
    Bnl_0001 = B1l_0001 * n1_0001 + B2l_0001 * n2_0001
    wl_0001 = sqrt(grav * hl_0001 + Bnl_0001 ** 2_i64)
    Unr_0001 = ur_0001 * n1_0001 + vr_0001 * n2_0001
    Bnr_0001 = B1r_0001 * n1_0001 + B2r_0001 * n2_0001
    wr_0001 = sqrt(grav * hr_0001 + Bnr_0001 ** 2_i64)
    !Utr     = vr*n1  - ur*n2
    !Btr     = B2r*n1 - B1r*n2
    !Utl     = vl*n1  - ul*n2
    !Btl     = B2l*n1 - B1l*n2
    lambda1l_0001 = abs(Unl_0001 - wl_0001)
    lambda2l_0001 = abs(Unl_0001 - Bnl_0001)
    lambda3l_0001 = abs(Unl_0001 + Bnl_0001)
    lambda4l_0001 = abs(Unl_0001 + wl_0001)
    lambda1r_0001 = abs(Unr_0001 - wr_0001)
    lambda2r_0001 = abs(Unr_0001 - Bnr_0001)
    lambda3r_0001 = abs(Unr_0001 + Bnr_0001)
    lambda4r_0001 = abs(Unr_0001 + wr_0001)
    ll_0001 = maxval([lambda1l_0001, lambda2l_0001, lambda3l_0001, &
          lambda4l_0001])
    lr_0001 = maxval([lambda1r_0001, lambda2r_0001, lambda3r_0001, &
          lambda4r_0001])
    lambda_star_0001 = maxval([ll_0001, lr_0001])
    q_l_0001 = hu_l * norm_0005(0_i64) + hv_l * norm_0005(1_i64)
    m_l_0001 = hB1_l * norm_0005(0_i64) + hB2_l * norm_0005(1_i64)
    p_l_0001 = 0.5_f64 * grav * h_l * h_l
    q_r_0001 = hu_r * norm_0005(0_i64) + hv_r * norm_0005(1_i64)
    m_r_0001 = hB1_r * norm_0005(0_i64) + hB2_r * norm_0005(1_i64)
    p_r_0001 = 0.5_f64 * grav * h_r * h_r
    fleft_h_0001 = q_l_0001
    fleft_hu_0001 = q_l_0001 * hu_l / h_l + p_l_0001 * norm_0005(0_i64) &
          - m_l_0001 * hB1_l / h_l
    fleft_hv_0001 = q_l_0001 * hv_l / h_l + p_l_0001 * norm_0005(1_i64) &
          - m_l_0001 * hB2_l / h_l
    fleft_hB1_0001 = (hv_l * hB1_l / h_l - hu_l * hB2_l / h_l) * &
          norm_0005(1_i64)
    fleft_hB2_0001 = (hu_l * hB2_l / h_l - hv_l * hB1_l / h_l) * &
          norm_0005(0_i64)
    fright_h_0001 = q_r_0001
    fright_hu_0001 = q_r_0001 * hu_r / h_r + p_r_0001 * norm_0005(0_i64 &
          ) - m_r_0001 * hB1_r / h_r
    fright_hv_0001 = q_r_0001 * hv_r / h_r + p_r_0001 * norm_0005(1_i64 &
          ) - m_r_0001 * hB2_r / h_r
    fright_hB1_0001 = (hv_r * hB1_r / h_r - hu_r * hB2_r / h_r) * &
          norm_0005(1_i64)
    fright_hB2_0001 = (hu_r * hB2_r / h_r - hv_r * hB1_r / h_r) * &
          norm_0005(0_i64)
    f_h_0001 = 0.5_f64 * (fleft_h_0001 + fright_h_0001) - 0.5_f64 * &
          lambda_star_0001 * (h_r - h_l)
    f_hu_0001 = 0.5_f64 * (fleft_hu_0001 + fright_hu_0001) - 0.5_f64 * &
          lambda_star_0001 * (hu_r - hu_l)
    f_hv_0001 = 0.5_f64 * (fleft_hv_0001 + fright_hv_0001) - 0.5_f64 * &
          lambda_star_0001 * (hv_r - hv_l)
    f_hB1_0001 = 0.5_f64 * (fleft_hB1_0001 + fright_hB1_0001) - 0.5_f64 &
          * lambda_star_0001 * (hB1_r - hB1_l)
    f_hB2_0001 = 0.5_f64 * (fleft_hB2_0001 + fright_hB2_0001) - 0.5_f64 &
          * lambda_star_0001 * (hB2_r - hB2_l)
    flux(0_i64) = f_h_0001 * mesure
    flux(1_i64) = f_hu_0001 * mesure
    flux(2_i64) = f_hv_0001 * mesure
    flux(3_i64) = f_hB1_0001 * mesure
    flux(4_i64) = f_hB2_0001 * mesure
    flux(5_i64) = 0.0_f64
    flux(6_i64) = 0.0_f64
    if (allocated(norm_0005)) then
      deallocate(norm_0005)
    end if
    if (allocated(ninv_0005)) then
      deallocate(ninv_0005)
    end if

  end subroutine bind_c_lf_scheme_mhd
  !........................................

  !........................................
  subroutine bind_c_srnh_scheme_mhd(hu_l, hu_r, hv_l, hv_r, h_l, h_r, &
        hB1_l, hB1_r, hB2_l, hB2_r, hPSI_l, hPSI_r, hB1c_l, hB1c_r, &
        hB2c_l, hB2c_r, Z_l, Z_r, n0_normal, normal, mesure, grav, &
        n0_flux, flux, cpsi) bind(c)

    implicit none

    real(f64), value :: hu_l
    real(f64), value :: hu_r
    real(f64), value :: hv_l
    real(f64), value :: hv_r
    real(f64), value :: h_l
    real(f64), value :: h_r
    real(f64), value :: hB1_l
    real(f64), value :: hB1_r
    real(f64), value :: hB2_l
    real(f64), value :: hB2_r
    real(f64), value :: hPSI_l
    real(f64), value :: hPSI_r
    real(f64), value :: hB1c_l
    real(f64), value :: hB1c_r
    real(f64), value :: hB2c_l
    real(f64), value :: hB2c_r
    real(f64), value :: Z_l
    real(f64), value :: Z_r
    integer(i64), value :: n0_normal
    real(f64), intent(in) :: normal(0:n0_normal - 1_i64)
    real(f64), value :: mesure
    real(f64), value :: grav
    integer(i64), value :: n0_flux
    real(f64), intent(inout) :: flux(0:n0_flux - 1_i64)
    real(f64), value :: cpsi
    real(f64), allocatable :: ninv_0006(:)
    real(f64), allocatable :: w_dif_0005(:)
    real(f64) :: u_h_0006
    real(f64) :: v_h_0006
    real(f64) :: B1_h_0006
    real(f64) :: B2_h_0006
    real(f64) :: un_h_0006
    real(f64) :: vn_h_0006
    real(f64) :: B1n_h_0006
    real(f64) :: B2n_h_0006
    real(f64) :: hroe_0006
    real(f64) :: uroe_0006
    real(f64) :: vroe_0006
    real(f64) :: B1roe_0006
    real(f64) :: B2roe_0006
    real(f64) :: uleft_0005
    real(f64) :: vleft_0005
    real(f64) :: B1left_0005
    real(f64) :: B2left_0005
    real(f64) :: uright_0005
    real(f64) :: vright_0005
    real(f64) :: B1right_0005
    real(f64) :: B2right_0005
    real(f64) :: B1i_0005
    real(f64) :: B1j_0005
    real(f64) :: absy_0005
    real(f64) :: w_lrh_0005
    real(f64) :: w_lrhu_0005
    real(f64) :: w_lrhv_0005
    real(f64) :: w_lrhB1_0005
    real(f64) :: w_lrhB2_0005
    real(f64) :: w_lrhPSI_0005
    real(f64) :: w_lrz_0005
    real(f64), allocatable, target :: signA_0005(:,:)
    real(f64) :: sound_0005
    real(f64) :: w_0005
    real(f64) :: lambda1_0005
    real(f64) :: lambda2_0005
    real(f64) :: lambda3_0005
    real(f64) :: lambda4_0005
    real(f64) :: epsilon_0005
    real(f64) :: s1_0005
    real(f64) :: pi1_0005
    real(f64) :: s2_0005
    real(f64) :: pi2_0005
    real(f64) :: s3_0005
    real(f64) :: pi3_0005
    real(f64) :: s4_0005
    real(f64) :: pi4_0005
    real(f64) :: gamma1_0005
    real(f64) :: gamma2_0005
    real(f64) :: sigma1_0005
    real(f64) :: sigma2_0005
    real(f64) :: mu1_0005
    real(f64) :: mu2_0005
    integer(i64) :: ann_0005
    real(f64), pointer :: smmat_0005(:,:)
    real(f64) :: hnew_0005
    real(f64) :: unew_0005
    real(f64) :: vnew_0005
    real(f64) :: B1new_0005
    real(f64) :: B2new_0005
    real(f64) :: znew_0005
    real(f64) :: Pnew_0005
    real(f64) :: u_hu_0005
    real(f64) :: u_hv_0005
    real(f64) :: u_hP_0005
    real(f64) :: u_hB1_0005
    real(f64) :: u_hB2_0005
    real(f64) :: u_z_0005
    real(f64) :: w_lrhP_0005
    real(f64) :: w_hP_0005
    real(f64) :: mw_hB1_0005
    real(f64) :: mhP_0005
    real(f64), allocatable :: norm_0006(:)
    real(f64) :: q_s_0005
    real(f64) :: p_s_0005
    real(f64) :: Flux_B1psi_0005
    real(f64) :: Flux_B2psi_0005
    real(f64) :: Flux_hPpsi_0005
    integer(i64) :: i_0005

    allocate(ninv_0006(0:1_i64))
    ninv_0006 = 0.0_f64
    allocate(w_dif_0005(0:5_i64))
    w_dif_0005 = 0.0_f64
    ninv_0006(0_i64) = (-1_i64) * normal(1_i64)
    ninv_0006(1_i64) = normal(0_i64)
    u_h_0006 = (hu_l / h_l * sqrt(h_l) + hu_r / h_r * sqrt(h_r)) / (sqrt &
          (h_l) + sqrt(h_r))
    v_h_0006 = (hv_l / h_l * sqrt(h_l) + hv_r / h_r * sqrt(h_r)) / (sqrt &
          (h_l) + sqrt(h_r))
    B1_h_0006 = (hB1_l / h_l * sqrt(h_l) + hB1_r / h_r * sqrt(h_r)) / ( &
          sqrt(h_l) + sqrt(h_r))
    B2_h_0006 = (hB2_l / h_l * sqrt(h_l) + hB2_r / h_r * sqrt(h_r)) / ( &
          sqrt(h_l) + sqrt(h_r))
    !uvh =  array([uh, vh])
    un_h_0006 = u_h_0006 * normal(0_i64) + v_h_0006 * normal(1_i64)
    un_h_0006 = un_h_0006 / mesure
    vn_h_0006 = u_h_0006 * ninv_0006(0_i64) + v_h_0006 * ninv_0006(1_i64 &
          )
    vn_h_0006 = vn_h_0006 / mesure
    B1n_h_0006 = B1_h_0006 * normal(0_i64) + B2_h_0006 * normal(1_i64)
    B1n_h_0006 = B1n_h_0006 / mesure
    B2n_h_0006 = B1_h_0006 * ninv_0006(0_i64) + B2_h_0006 * ninv_0006( &
          1_i64)
    B2n_h_0006 = B2n_h_0006 / mesure
    hroe_0006 = (h_l + h_r) / 2_i64
    uroe_0006 = un_h_0006
    vroe_0006 = vn_h_0006
    B1roe_0006 = B1n_h_0006
    B2roe_0006 = B2n_h_0006
    uleft_0005 = hu_l * normal(0_i64) + hv_l * normal(1_i64)
    uleft_0005 = uleft_0005 / mesure
    vleft_0005 = hu_l * ninv_0006(0_i64) + hv_l * ninv_0006(1_i64)
    vleft_0005 = vleft_0005 / mesure
    B1left_0005 = hB1_l * normal(0_i64) + hB2_l * normal(1_i64)
    B1left_0005 = B1left_0005 / mesure
    B2left_0005 = hB1_l * ninv_0006(0_i64) + hB2_l * ninv_0006(1_i64)
    B2left_0005 = B2left_0005 / mesure
    uright_0005 = hu_r * normal(0_i64) + hv_r * normal(1_i64)
    uright_0005 = uright_0005 / mesure
    vright_0005 = hu_r * ninv_0006(0_i64) + hv_r * ninv_0006(1_i64)
    vright_0005 = vright_0005 / mesure
    B1right_0005 = hB1_r * normal(0_i64) + hB2_r * normal(1_i64)
    B1right_0005 = B1right_0005 / mesure
    B2right_0005 = hB1_r * ninv_0006(0_i64) + hB2_r * ninv_0006(1_i64)
    B2right_0005 = B2right_0005 / mesure
    B1i_0005 = (hB1c_l * normal(0_i64) + hB2c_l * normal(1_i64)) / &
          mesure
    B1j_0005 = (hB1c_r * normal(0_i64) + hB2c_r * normal(1_i64)) / &
          mesure
    absy_0005 = (B1i_0005 + B1j_0005) / 2_i64
    w_lrh_0005 = (h_l + h_r) / 2_i64
    w_lrhu_0005 = (uleft_0005 + uright_0005) / 2_i64
    w_lrhv_0005 = (vleft_0005 + vright_0005) / 2_i64
    w_lrhB1_0005 = (B1left_0005 + B1right_0005) / 2_i64
    w_lrhB2_0005 = (B2left_0005 + B2right_0005) / 2_i64
    w_lrhPSI_0005 = (hPSI_l + hPSI_r) / 2_i64
    w_lrz_0005 = (Z_l + Z_r) / 2_i64
    w_dif_0005(0_i64) = h_r - h_l
    w_dif_0005(1_i64) = uright_0005 - uleft_0005
    w_dif_0005(2_i64) = vright_0005 - vleft_0005
    w_dif_0005(3_i64) = B1right_0005 - B1left_0005
    w_dif_0005(4_i64) = B2right_0005 - B2left_0005
    w_dif_0005(5_i64) = Z_r - Z_l
    allocate(signA_0005(0:5_i64, 0:5_i64))
    signA_0005 = 0.0_f64
    sound_0005 = sqrt(grav * hroe_0006)
    !B1roe = absy/hroe
    w_0005 = sqrt(B1roe_0006 * B1roe_0006 + grav * hroe_0006)
    lambda1_0005 = uroe_0006 - w_0005
    lambda2_0005 = uroe_0006 - B1roe_0006
    lambda3_0005 = uroe_0006 + B1roe_0006
    lambda4_0005 = uroe_0006 + w_0005
    epsilon_0005 = 1e-15_f64
    if (abs(lambda1_0005) < epsilon_0005) then
      s1_0005 = 0.0_f64
      pi1_0005 = 0.0_f64
    else
      s1_0005 = lambda1_0005 / abs(lambda1_0005)
      pi1_0005 = s1_0005 / lambda1_0005
    end if
    if (abs(lambda2_0005) < epsilon_0005) then
      s2_0005 = 0.0_f64
      pi2_0005 = 0.0_f64
    else
      s2_0005 = lambda2_0005 / abs(lambda2_0005)
      pi2_0005 = 1.0_f64 / abs(lambda2_0005)
    end if
    if (abs(lambda3_0005) < epsilon_0005) then
      s3_0005 = 0.0_f64
      pi3_0005 = 0.0_f64
    else
      s3_0005 = lambda3_0005 / abs(lambda3_0005)
      pi3_0005 = 1.0_f64 / abs(lambda3_0005)
    end if
    if (abs(lambda4_0005) < epsilon_0005) then
      s4_0005 = 0.0_f64
      pi4_0005 = 0.0_f64
    else
      s4_0005 = lambda4_0005 / abs(lambda4_0005)
      pi4_0005 = 1.0_f64 / abs(lambda4_0005)
    end if
    gamma1_0005 = vroe_0006 + B2roe_0006
    gamma2_0005 = vroe_0006 - B2roe_0006
    sigma1_0005 = vroe_0006 * (s1_0005 * lambda4_0005 - s4_0005 * &
          lambda1_0005) - w_0005 * (s2_0005 * gamma1_0005 + s3_0005 * &
          gamma2_0005)
    sigma2_0005 = B2roe_0006 * (s1_0005 * lambda4_0005 - s4_0005 * &
          lambda1_0005) - w_0005 * (s2_0005 * gamma1_0005 - s3_0005 * &
          gamma2_0005)
    if (abs(lambda2_0005) < epsilon_0005 .and. abs(lambda3_0005) < &
          epsilon_0005) then
      mu1_0005 = B1roe_0006 * vroe_0006 * pi1_0005 / w_0005 - B1roe_0006 &
            * vroe_0006 * pi4_0005 / w_0005
      mu2_0005 = B1roe_0006 * B2roe_0006 * pi1_0005 / w_0005 - &
            B1roe_0006 * B2roe_0006 * pi4_0005 / w_0005
      ann_0005 = 1_i64
    else
      mu1_0005 = B1roe_0006 * vroe_0006 * pi1_0005 / w_0005 - B1roe_0006 &
            * vroe_0006 * pi4_0005 / w_0005 - 0.5_f64 * (gamma1_0005 * &
            pi2_0005 - gamma2_0005 * pi3_0005)
      mu2_0005 = B1roe_0006 * B2roe_0006 * pi1_0005 / w_0005 - &
            B1roe_0006 * B2roe_0006 * pi4_0005 / w_0005 - 0.5_f64 * ( &
            gamma1_0005 * pi2_0005 + gamma2_0005 * pi3_0005)
      ann_0005 = 1_i64
    end if
    !1ère colonne de la matrice A
    signA_0005(0_i64, 0_i64) = (s1_0005 * lambda4_0005 - s4_0005 * &
          lambda1_0005) / (2_i64 * w_0005)
    signA_0005(0_i64, 1_i64) = lambda1_0005 * lambda4_0005 * (s1_0005 - &
          s4_0005) / (2_i64 * w_0005)
    signA_0005(0_i64, 2_i64) = sigma1_0005 / (2_i64 * w_0005)
    signA_0005(0_i64, 3_i64) = 0.0_f64
    signA_0005(0_i64, 4_i64) = sigma2_0005 / (2_i64 * w_0005)
    signA_0005(0_i64, 5_i64) = 0.0_f64
    !2ème colonne de la matrice A
    signA_0005(1_i64, 0_i64) = (s4_0005 - s1_0005) / (2_i64 * w_0005)
    signA_0005(1_i64, 1_i64) = (s4_0005 * lambda4_0005 - s1_0005 * &
          lambda1_0005) / (2_i64 * w_0005)
    signA_0005(1_i64, 2_i64) = vroe_0006 * (s4_0005 - s1_0005) / (2_i64 &
          * w_0005)
    signA_0005(1_i64, 3_i64) = 0.0_f64
    signA_0005(1_i64, 4_i64) = B2roe_0006 * (s4_0005 - s1_0005) / (2_i64 &
          * w_0005)
    signA_0005(1_i64, 5_i64) = 0.0_f64
    !3ème colonne de la matrice A
    signA_0005(2_i64, 0_i64) = 0.0_f64
    signA_0005(2_i64, 1_i64) = 0.0_f64
    signA_0005(2_i64, 2_i64) = (s2_0005 + s3_0005) / 2_i64
    signA_0005(2_i64, 3_i64) = 0.0_f64
    signA_0005(2_i64, 4_i64) = (s2_0005 - s3_0005) / 2_i64
    signA_0005(2_i64, 5_i64) = 0.0_f64
    !4ème colonne de la matrice A
    signA_0005(3_i64, 0_i64) = ann_0005 * B1roe_0006 * (pi1_0005 - &
          pi4_0005) / w_0005
    signA_0005(3_i64, 1_i64) = ann_0005 * B1roe_0006 * (s1_0005 - &
          s4_0005) / w_0005
    signA_0005(3_i64, 2_i64) = ann_0005 * mu1_0005
    signA_0005(3_i64, 3_i64) = 0.0_f64
    signA_0005(3_i64, 4_i64) = ann_0005 * mu2_0005
    signA_0005(3_i64, 5_i64) = 0.0_f64
    !5ème colonne de la matrice A
    signA_0005(4_i64, 0_i64) = 0.0_f64
    signA_0005(4_i64, 1_i64) = 0.0_f64
    signA_0005(4_i64, 2_i64) = (s2_0005 - s3_0005) / 2_i64
    signA_0005(4_i64, 3_i64) = 0.0_f64
    signA_0005(4_i64, 4_i64) = (s2_0005 + s3_0005) / 2_i64
    signA_0005(4_i64, 5_i64) = 0.0_f64
    !6ème colonne de la matrice A
    signA_0005(5_i64, 0_i64) = sound_0005 ** 2_i64 * (pi4_0005 - &
          pi1_0005) / (2_i64 * w_0005)
    signA_0005(5_i64, 1_i64) = sound_0005 ** 2_i64 * (s4_0005 - s1_0005 &
          ) / (2_i64 * w_0005)
    signA_0005(5_i64, 2_i64) = sound_0005 ** 2_i64 * vroe_0006 * ( &
          pi4_0005 - pi1_0005) / (2_i64 * w_0005)
    signA_0005(5_i64, 3_i64) = 0.0_f64
    signA_0005(5_i64, 4_i64) = sound_0005 ** 2_i64 * B2roe_0006 * ( &
          pi4_0005 - pi1_0005) / (2_i64 * w_0005)
    signA_0005(5_i64, 5_i64) = 0.0_f64
    smmat_0005(0:, 0:) => signA_0005
    hnew_0005 = 0.0_f64
    unew_0005 = 0.0_f64
    vnew_0005 = 0.0_f64
    B1new_0005 = 0.0_f64
    B2new_0005 = 0.0_f64
    znew_0005 = 0.0_f64
    do i_0005 = 0_i64, 5_i64, 1_i64
      hnew_0005 = hnew_0005 + smmat_0005(i_0005, 0_i64) * w_dif_0005( &
            i_0005)
      unew_0005 = unew_0005 + smmat_0005(i_0005, 1_i64) * w_dif_0005( &
            i_0005)
      vnew_0005 = vnew_0005 + smmat_0005(i_0005, 2_i64) * w_dif_0005( &
            i_0005)
      B1new_0005 = B1new_0005 + smmat_0005(i_0005, 3_i64) * w_dif_0005( &
            i_0005)
      B2new_0005 = B2new_0005 + smmat_0005(i_0005, 4_i64) * w_dif_0005( &
            i_0005)
      znew_0005 = znew_0005 + smmat_0005(i_0005, 5_i64) * w_dif_0005( &
            i_0005)
    end do
    Pnew_0005 = cpsi * (B1right_0005 - B1left_0005)
    u_h_0006 = hnew_0005 / 2_i64
    u_hu_0005 = unew_0005 / 2_i64
    u_hv_0005 = vnew_0005 / 2_i64
    u_hP_0005 = Pnew_0005 / 2_i64
    u_hB1_0005 = B1new_0005 / 2_i64
    u_hB2_0005 = B2new_0005 / 2_i64
    u_z_0005 = znew_0005 / 2_i64
    w_lrh_0005 = w_lrh_0005 - u_h_0006
    w_lrhu_0005 = w_lrhu_0005 - u_hu_0005
    w_lrhv_0005 = w_lrhv_0005 - u_hv_0005
    w_lrhP_0005 = w_lrhPSI_0005 - u_hP_0005
    w_lrhB1_0005 = w_lrhB1_0005 - u_hB1_0005
    w_lrhB2_0005 = w_lrhB2_0005 - u_hB2_0005
    w_lrz_0005 = w_lrz_0005 - u_z_0005
    w_hP_0005 = hPSI_r - hPSI_l
    mw_hB1_0005 = w_lrhB1_0005 / 2_i64 - 1_i64 / (2_i64 * cpsi) * &
          w_hP_0005
    mhP_0005 = (hPSI_r - hPSI_l) / 2_i64 - cpsi * w_dif_0005(3_i64) / &
          2_i64
    unew_0005 = 0.0_f64
    vnew_0005 = 0.0_f64
    B1new_0005 = 0.0_f64
    B2new_0005 = 0.0_f64
    unew_0005 = w_lrhu_0005 * normal(0_i64) - w_lrhv_0005 * normal(1_i64 &
          )
    unew_0005 = unew_0005 / mesure
    vnew_0005 = w_lrhu_0005 * normal(1_i64) + w_lrhv_0005 * normal(0_i64 &
          )
    vnew_0005 = vnew_0005 / mesure
    B1new_0005 = w_lrhB1_0005 * normal(0_i64) - w_lrhB2_0005 * normal( &
          1_i64)
    B1new_0005 = B1new_0005 / mesure
    B2new_0005 = w_lrhB1_0005 * normal(1_i64) + w_lrhB2_0005 * normal( &
          0_i64)
    B2new_0005 = B2new_0005 / mesure
    w_lrhu_0005 = unew_0005
    w_lrhv_0005 = vnew_0005
    w_lrhB1_0005 = B1new_0005
    w_lrhB2_0005 = B2new_0005
    allocate(norm_0006(0:size(normal, kind=i64) - 1_i64))
    norm_0006 = normal / mesure
    q_s_0005 = normal(0_i64) * unew_0005 + normal(1_i64) * vnew_0005
    p_s_0005 = normal(0_i64) * B1new_0005 + normal(1_i64) * B2new_0005
    Flux_B1psi_0005 = mhP_0005 * norm_0006(0_i64) * mesure
    Flux_B2psi_0005 = mhP_0005 * norm_0006(1_i64) * mesure
    Flux_hPpsi_0005 = cpsi * cpsi * mw_hB1_0005 * mesure
    flux(0_i64) = q_s_0005
    flux(1_i64) = q_s_0005 * w_lrhu_0005 / w_lrh_0005 + 0.5_f64 * grav * &
          w_lrh_0005 * w_lrh_0005 * normal(0_i64) - p_s_0005 * &
          w_lrhB1_0005 / w_lrh_0005
    flux(2_i64) = q_s_0005 * w_lrhv_0005 / w_lrh_0005 + 0.5_f64 * grav * &
          w_lrh_0005 * w_lrh_0005 * normal(1_i64) - p_s_0005 * &
          w_lrhB2_0005 / w_lrh_0005
    flux(3_i64) = (w_lrhv_0005 * w_lrhB1_0005 / w_lrh_0005 - w_lrhu_0005 &
          * w_lrhB2_0005 / w_lrh_0005) * normal(1_i64) + &
          Flux_B1psi_0005
    flux(4_i64) = (w_lrhu_0005 * w_lrhB2_0005 / w_lrh_0005 - w_lrhv_0005 &
          * w_lrhB1_0005 / w_lrh_0005) * normal(0_i64) + &
          Flux_B2psi_0005
    flux(5_i64) = Flux_hPpsi_0005
    flux(6_i64) = 0.0_f64
    if (allocated(ninv_0006)) then
      deallocate(ninv_0006)
    end if
    if (allocated(w_dif_0005)) then
      deallocate(w_dif_0005)
    end if
    if (allocated(signA_0005)) then
      deallocate(signA_0005)
    end if
    if (allocated(norm_0006)) then
      deallocate(norm_0006)
    end if

  end subroutine bind_c_srnh_scheme_mhd
  !........................................

  !........................................
  subroutine bind_c_explicitscheme_convective_swmhd(n0_rez_h, rez_h, &
        n0_rez_hu, rez_hu, n0_rez_hv, rez_hv, n0_rez_hB1, rez_hB1, &
        n0_rez_hB2, rez_hB2, n0_rez_PSI, rez_PSI, n0_rez_Z, rez_Z, &
        n0_h_c, h_c, n0_hu_c, hu_c, n0_hv_c, hv_c, n0_hB1_c, hB1_c, &
        n0_hB2_c, hB2_c, n0_hPSIc, hPSIc, n0_Z_c, Z_c, n0_h_ghost, &
        h_ghost, n0_hu_ghost, hu_ghost, n0_hv_ghost, hv_ghost, &
        n0_hB1_ghost, hB1_ghost, n0_hB2_ghost, hB2_ghost, n0_hPSIghost, &
        hPSIghost, n0_Z_ghost, Z_ghost, n0_h_halo, h_halo, n0_hu_halo, &
        hu_halo, n0_hv_halo, hv_halo, n0_hB1_halo, hB1_halo, &
        n0_hB2_halo, hB2_halo, n0_hPSIhalo, hPSIhalo, n0_Z_halo, Z_halo &
        , n0_h_x, h_x, n0_h_y, h_y, n0_hx_halo, hx_halo, n0_hy_halo, &
        hy_halo, n0_psi, psi, n0_psi_halo, psi_halo, n0_centerc, &
        n1_centerc, centerc, n0_centerf, n1_centerf, centerf, &
        n0_centerh, n1_centerh, centerh, n0_centerg, n1_centerg, &
        centerg, n0_cellidf, n1_cellidf, cellidf, n0_mesuref, mesuref, &
        n0_normalf, n1_normalf, normalf, n0_halofid, halofid, &
        n0_innerfaces, innerfaces, n0_halofaces, halofaces, &
        n0_boundaryfaces, boundaryfaces, n0_periodicboundaryfaces, &
        periodicboundaryfaces, n0_shift, n1_shift, shift, order, cpsi, &
        n0_hB1_cst, hB1_cst, n0_hB2_cst, hB2_cst) bind(c)

    implicit none

    integer(i64), value :: n0_rez_h
    real(f64), intent(inout) :: rez_h(0:n0_rez_h - 1_i64)
    integer(i64), value :: n0_rez_hu
    real(f64), intent(inout) :: rez_hu(0:n0_rez_hu - 1_i64)
    integer(i64), value :: n0_rez_hv
    real(f64), intent(inout) :: rez_hv(0:n0_rez_hv - 1_i64)
    integer(i64), value :: n0_rez_hB1
    real(f64), intent(inout) :: rez_hB1(0:n0_rez_hB1 - 1_i64)
    integer(i64), value :: n0_rez_hB2
    real(f64), intent(inout) :: rez_hB2(0:n0_rez_hB2 - 1_i64)
    integer(i64), value :: n0_rez_PSI
    real(f64), intent(inout) :: rez_PSI(0:n0_rez_PSI - 1_i64)
    integer(i64), value :: n0_rez_Z
    real(f64), intent(inout) :: rez_Z(0:n0_rez_Z - 1_i64)
    integer(i64), value :: n0_h_c
    real(f64), intent(in) :: h_c(0:n0_h_c - 1_i64)
    integer(i64), value :: n0_hu_c
    real(f64), intent(in) :: hu_c(0:n0_hu_c - 1_i64)
    integer(i64), value :: n0_hv_c
    real(f64), intent(in) :: hv_c(0:n0_hv_c - 1_i64)
    integer(i64), value :: n0_hB1_c
    real(f64), intent(in) :: hB1_c(0:n0_hB1_c - 1_i64)
    integer(i64), value :: n0_hB2_c
    real(f64), intent(in) :: hB2_c(0:n0_hB2_c - 1_i64)
    integer(i64), value :: n0_hPSIc
    real(f64), intent(in) :: hPSIc(0:n0_hPSIc - 1_i64)
    integer(i64), value :: n0_Z_c
    real(f64), intent(in) :: Z_c(0:n0_Z_c - 1_i64)
    integer(i64), value :: n0_h_ghost
    real(f64), intent(in) :: h_ghost(0:n0_h_ghost - 1_i64)
    integer(i64), value :: n0_hu_ghost
    real(f64), intent(in) :: hu_ghost(0:n0_hu_ghost - 1_i64)
    integer(i64), value :: n0_hv_ghost
    real(f64), intent(in) :: hv_ghost(0:n0_hv_ghost - 1_i64)
    integer(i64), value :: n0_hB1_ghost
    real(f64), intent(in) :: hB1_ghost(0:n0_hB1_ghost - 1_i64)
    integer(i64), value :: n0_hB2_ghost
    real(f64), intent(in) :: hB2_ghost(0:n0_hB2_ghost - 1_i64)
    integer(i64), value :: n0_hPSIghost
    real(f64), intent(in) :: hPSIghost(0:n0_hPSIghost - 1_i64)
    integer(i64), value :: n0_Z_ghost
    real(f64), intent(in) :: Z_ghost(0:n0_Z_ghost - 1_i64)
    integer(i64), value :: n0_h_halo
    real(f64), intent(in) :: h_halo(0:n0_h_halo - 1_i64)
    integer(i64), value :: n0_hu_halo
    real(f64), intent(in) :: hu_halo(0:n0_hu_halo - 1_i64)
    integer(i64), value :: n0_hv_halo
    real(f64), intent(in) :: hv_halo(0:n0_hv_halo - 1_i64)
    integer(i64), value :: n0_hB1_halo
    real(f64), intent(in) :: hB1_halo(0:n0_hB1_halo - 1_i64)
    integer(i64), value :: n0_hB2_halo
    real(f64), intent(in) :: hB2_halo(0:n0_hB2_halo - 1_i64)
    integer(i64), value :: n0_hPSIhalo
    real(f64), intent(in) :: hPSIhalo(0:n0_hPSIhalo - 1_i64)
    integer(i64), value :: n0_Z_halo
    real(f64), intent(in) :: Z_halo(0:n0_Z_halo - 1_i64)
    integer(i64), value :: n0_h_x
    real(f64), intent(in) :: h_x(0:n0_h_x - 1_i64)
    integer(i64), value :: n0_h_y
    real(f64), intent(in) :: h_y(0:n0_h_y - 1_i64)
    integer(i64), value :: n0_hx_halo
    real(f64), intent(in) :: hx_halo(0:n0_hx_halo - 1_i64)
    integer(i64), value :: n0_hy_halo
    real(f64), intent(in) :: hy_halo(0:n0_hy_halo - 1_i64)
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
    real(f64), value :: cpsi
    integer(i64), value :: n0_hB1_cst
    real(f64), intent(in) :: hB1_cst(0:n0_hB1_cst - 1_i64)
    integer(i64), value :: n0_hB2_cst
    real(f64), intent(in) :: hB2_cst(0:n0_hB2_cst - 1_i64)

    call explicitscheme_convective_SWMHD(rez_h, rez_hu, rez_hv, rez_hB1, &
          rez_hB2, rez_PSI, rez_Z, h_c, hu_c, hv_c, hB1_c, hB2_c, hPSIc &
          , Z_c, h_ghost, hu_ghost, hv_ghost, hB1_ghost, hB2_ghost, &
          hPSIghost, Z_ghost, h_halo, hu_halo, hv_halo, hB1_halo, &
          hB2_halo, hPSIhalo, Z_halo, h_x, h_y, hx_halo, hy_halo, psi, &
          psi_halo, centerc, centerf, centerh, centerg, cellidf, &
          mesuref, normalf, halofid, innerfaces, halofaces, &
          boundaryfaces, periodicboundaryfaces, shift, order, cpsi, &
          hB1_cst, hB2_cst)

  end subroutine bind_c_explicitscheme_convective_swmhd
  !........................................

end module bind_c_tools_SWMHD
