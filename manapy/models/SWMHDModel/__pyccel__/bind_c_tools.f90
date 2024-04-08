module bind_c_tools

  use tools, only: term_source_srnh_SWf
  use tools, only: term_coriolis_SW
  use tools, only: term_wind_SW
  use tools, only: update_SW
  use tools, only: explicitscheme_convective_SW
  use tools, only: term_friction_SW
  use tools, only: update_SWMHD
  use tools, only: initialisation_SW
  use tools, only: term_source_srnh_SWMHD
  use tools, only: time_step_SW
  use tools, only: time_step_SWMHD
  use tools, only: explicitscheme_convective_SWMHD
  use tools, only: cpsi_global

  use, intrinsic :: ISO_C_Binding, only : i64 => C_INT64_T , f64 => &
        C_DOUBLE
  implicit none

  contains

  !........................................
  subroutine bind_c_update_sw(n0_h_c, h_c, n0_hu_c, hu_c, n0_hv_c, hv_c, &
        n0_hc_c, hc_c, n0_Z_c, Z_c, n0_rez_h, rez_h, n0_rez_hu, rez_hu, &
        n0_rez_hv, rez_hv, n0_rez_hc, rez_hc, n0_rez_Z, rez_Z, n0_src_h &
        , src_h, n0_src_hu, src_hu, n0_src_hv, src_hv, n0_src_hc, &
        src_hc, n0_src_Z, src_Z, n0_dissip_hc, dissip_hc, n0_corio_hu, &
        corio_hu, n0_corio_hv, corio_hv, wind_hu, wind_hv, dtime, &
        n0_vol, vol) bind(c)

    implicit none

    integer(i64), value :: n0_h_c
    real(f64), intent(inout) :: h_c(0:n0_h_c - 1_i64)
    integer(i64), value :: n0_hu_c
    real(f64), intent(inout) :: hu_c(0:n0_hu_c - 1_i64)
    integer(i64), value :: n0_hv_c
    real(f64), intent(inout) :: hv_c(0:n0_hv_c - 1_i64)
    integer(i64), value :: n0_hc_c
    real(f64), intent(inout) :: hc_c(0:n0_hc_c - 1_i64)
    integer(i64), value :: n0_Z_c
    real(f64), intent(inout) :: Z_c(0:n0_Z_c - 1_i64)
    integer(i64), value :: n0_rez_h
    real(f64), intent(in) :: rez_h(0:n0_rez_h - 1_i64)
    integer(i64), value :: n0_rez_hu
    real(f64), intent(in) :: rez_hu(0:n0_rez_hu - 1_i64)
    integer(i64), value :: n0_rez_hv
    real(f64), intent(in) :: rez_hv(0:n0_rez_hv - 1_i64)
    integer(i64), value :: n0_rez_hc
    real(f64), intent(in) :: rez_hc(0:n0_rez_hc - 1_i64)
    integer(i64), value :: n0_rez_Z
    real(f64), intent(in) :: rez_Z(0:n0_rez_Z - 1_i64)
    integer(i64), value :: n0_src_h
    real(f64), intent(in) :: src_h(0:n0_src_h - 1_i64)
    integer(i64), value :: n0_src_hu
    real(f64), intent(in) :: src_hu(0:n0_src_hu - 1_i64)
    integer(i64), value :: n0_src_hv
    real(f64), intent(in) :: src_hv(0:n0_src_hv - 1_i64)
    integer(i64), value :: n0_src_hc
    real(f64), intent(in) :: src_hc(0:n0_src_hc - 1_i64)
    integer(i64), value :: n0_src_Z
    real(f64), intent(in) :: src_Z(0:n0_src_Z - 1_i64)
    integer(i64), value :: n0_dissip_hc
    real(f64), intent(in) :: dissip_hc(0:n0_dissip_hc - 1_i64)
    integer(i64), value :: n0_corio_hu
    real(f64), intent(in) :: corio_hu(0:n0_corio_hu - 1_i64)
    integer(i64), value :: n0_corio_hv
    real(f64), intent(in) :: corio_hv(0:n0_corio_hv - 1_i64)
    real(f64), value :: wind_hu
    real(f64), value :: wind_hv
    real(f64), value :: dtime
    integer(i64), value :: n0_vol
    real(f64), intent(in) :: vol(0:n0_vol - 1_i64)

    call update_SW(h_c, hu_c, hv_c, hc_c, Z_c, rez_h, rez_hu, rez_hv, &
          rez_hc, rez_Z, src_h, src_hu, src_hv, src_hc, src_Z, &
          dissip_hc, corio_hu, corio_hv, wind_hu, wind_hv, dtime, vol)

  end subroutine bind_c_update_sw
  !........................................

  !........................................
  function bind_c_time_step_sw(n0_h_c, h_c, n0_hu_c, hu_c, n0_hv_c, hv_c &
        , cfl, n0_normal, n1_normal, normal, n0_mesure, mesure, &
        n0_volume, volume, n0_faceid, n1_faceid, faceid, Dxx, Dyy) bind &
        (c) result(dt)

    implicit none

    integer(i64), value :: n0_h_c
    real(f64), intent(in) :: h_c(0:n0_h_c - 1_i64)
    integer(i64), value :: n0_hu_c
    real(f64), intent(in) :: hu_c(0:n0_hu_c - 1_i64)
    integer(i64), value :: n0_hv_c
    real(f64), intent(in) :: hv_c(0:n0_hv_c - 1_i64)
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
    real(f64), value :: Dxx
    real(f64), value :: Dyy
    real(f64) :: dt

    dt = time_step_SW(h_c, hu_c, hv_c, cfl, normal, mesure, volume, &
          faceid, Dxx, Dyy)

  end function bind_c_time_step_sw
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
  subroutine bind_c_initialisation_sw(n0_h, h, n0_hu, hu, n0_hv, hv, &
        n0_hc, hc, n0_Z, Z, n0_center, n1_center, center) bind(c)

    implicit none

    integer(i64), value :: n0_h
    real(f64), intent(inout) :: h(0:n0_h - 1_i64)
    integer(i64), value :: n0_hu
    real(f64), intent(inout) :: hu(0:n0_hu - 1_i64)
    integer(i64), value :: n0_hv
    real(f64), intent(inout) :: hv(0:n0_hv - 1_i64)
    integer(i64), value :: n0_hc
    real(f64), intent(inout) :: hc(0:n0_hc - 1_i64)
    integer(i64), value :: n0_Z
    real(f64), intent(inout) :: Z(0:n0_Z - 1_i64)
    integer(i64), value :: n0_center
    integer(i64), value :: n1_center
    real(f64), intent(in) :: center(0:n1_center - 1_i64,0:n0_center - &
          1_i64)

    call initialisation_SW(h, hu, hv, hc, Z, center)

  end subroutine bind_c_initialisation_sw
  !........................................

  !........................................
  subroutine bind_c_term_source_srnh_swf(n0_src_h, src_h, n0_src_hu, &
        src_hu, n0_src_hv, src_hv, n0_src_hc, src_hc, n0_src_Z, src_Z, &
        n0_h_c, h_c, n0_hu_c, hu_c, n0_hv_c, hv_c, n0_hc_c, hc_c, &
        n0_Z_c, Z_c, n0_h_ghost, h_ghost, n0_hu_ghost, hu_ghost, &
        n0_hv_ghost, hv_ghost, n0_hc_ghost, hc_ghost, n0_Z_ghost, &
        Z_ghost, n0_h_halo, h_halo, n0_hu_halo, hu_halo, n0_hv_halo, &
        hv_halo, n0_hc_halo, hc_halo, n0_Z_halo, Z_halo, n0_h_x, h_x, &
        n0_h_y, h_y, n0_psi, psi, n0_hx_halo, hx_halo, n0_hy_halo, &
        hy_halo, n0_psi_halo, psi_halo, n0_nodeidc, n1_nodeidc, nodeidc &
        , n0_faceidc, n1_faceidc, faceidc, n0_cellidc, n1_cellidc, &
        cellidc, n0_cellidf, n1_cellidf, cellidf, n0_centerc, &
        n1_centerc, centerc, n0_normalc, n1_normalc, n2_normalc, &
        normalc, n0_namef, namef, n0_centerf, n1_centerf, centerf, &
        n0_centerh, n1_centerh, centerh, n0_vertexn, n1_vertexn, &
        vertexn, n0_halofid, halofid, order) bind(c)

    implicit none

    integer(i64), value :: n0_src_h
    real(f64), intent(inout) :: src_h(0:n0_src_h - 1_i64)
    integer(i64), value :: n0_src_hu
    real(f64), intent(inout) :: src_hu(0:n0_src_hu - 1_i64)
    integer(i64), value :: n0_src_hv
    real(f64), intent(inout) :: src_hv(0:n0_src_hv - 1_i64)
    integer(i64), value :: n0_src_hc
    real(f64), intent(inout) :: src_hc(0:n0_src_hc - 1_i64)
    integer(i64), value :: n0_src_Z
    real(f64), intent(inout) :: src_Z(0:n0_src_Z - 1_i64)
    integer(i64), value :: n0_h_c
    real(f64), intent(in) :: h_c(0:n0_h_c - 1_i64)
    integer(i64), value :: n0_hu_c
    real(f64), intent(in) :: hu_c(0:n0_hu_c - 1_i64)
    integer(i64), value :: n0_hv_c
    real(f64), intent(in) :: hv_c(0:n0_hv_c - 1_i64)
    integer(i64), value :: n0_hc_c
    real(f64), intent(in) :: hc_c(0:n0_hc_c - 1_i64)
    integer(i64), value :: n0_Z_c
    real(f64), intent(in) :: Z_c(0:n0_Z_c - 1_i64)
    integer(i64), value :: n0_h_ghost
    real(f64), intent(in) :: h_ghost(0:n0_h_ghost - 1_i64)
    integer(i64), value :: n0_hu_ghost
    real(f64), intent(in) :: hu_ghost(0:n0_hu_ghost - 1_i64)
    integer(i64), value :: n0_hv_ghost
    real(f64), intent(in) :: hv_ghost(0:n0_hv_ghost - 1_i64)
    integer(i64), value :: n0_hc_ghost
    real(f64), intent(in) :: hc_ghost(0:n0_hc_ghost - 1_i64)
    integer(i64), value :: n0_Z_ghost
    real(f64), intent(in) :: Z_ghost(0:n0_Z_ghost - 1_i64)
    integer(i64), value :: n0_h_halo
    real(f64), intent(in) :: h_halo(0:n0_h_halo - 1_i64)
    integer(i64), value :: n0_hu_halo
    real(f64), intent(in) :: hu_halo(0:n0_hu_halo - 1_i64)
    integer(i64), value :: n0_hv_halo
    real(f64), intent(in) :: hv_halo(0:n0_hv_halo - 1_i64)
    integer(i64), value :: n0_hc_halo
    real(f64), intent(in) :: hc_halo(0:n0_hc_halo - 1_i64)
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

    call term_source_srnh_SWf(src_h, src_hu, src_hv, src_hc, src_Z, h_c, &
          hu_c, hv_c, hc_c, Z_c, h_ghost, hu_ghost, hv_ghost, hc_ghost, &
          Z_ghost, h_halo, hu_halo, hv_halo, hc_halo, Z_halo, h_x, h_y, &
          psi, hx_halo, hy_halo, psi_halo, nodeidc, faceidc, cellidc, &
          cellidf, centerc, normalc, namef, centerf, centerh, vertexn, &
          halofid, order)

  end subroutine bind_c_term_source_srnh_swf
  !........................................

  !........................................
  subroutine bind_c_term_source_srnh_swmhd(n0_src_h, src_h, n0_src_hu, &
        src_hu, n0_src_hv, src_hv, n0_src_hB1, src_hB1, n0_src_hB2, &
        src_hB2, n0_src_PSI, src_PSI, n0_src_Z, src_Z, n0_h_c, h_c, &
        n0_hu_c, hu_c, n0_hv_c, hv_c, n0_hc_c, hc_c, n0_Z_c, Z_c, &
        n0_h_ghost, h_ghost, n0_hu_ghost, hu_ghost, n0_hv_ghost, &
        hv_ghost, n0_hc_ghost, hc_ghost, n0_Z_ghost, Z_ghost, n0_h_halo &
        , h_halo, n0_hu_halo, hu_halo, n0_hv_halo, hv_halo, n0_hc_halo, &
        hc_halo, n0_Z_halo, Z_halo, n0_h_x, h_x, n0_h_y, h_y, n0_psi, &
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
    integer(i64), value :: n0_hc_c
    real(f64), intent(in) :: hc_c(0:n0_hc_c - 1_i64)
    integer(i64), value :: n0_Z_c
    real(f64), intent(in) :: Z_c(0:n0_Z_c - 1_i64)
    integer(i64), value :: n0_h_ghost
    real(f64), intent(in) :: h_ghost(0:n0_h_ghost - 1_i64)
    integer(i64), value :: n0_hu_ghost
    real(f64), intent(in) :: hu_ghost(0:n0_hu_ghost - 1_i64)
    integer(i64), value :: n0_hv_ghost
    real(f64), intent(in) :: hv_ghost(0:n0_hv_ghost - 1_i64)
    integer(i64), value :: n0_hc_ghost
    real(f64), intent(in) :: hc_ghost(0:n0_hc_ghost - 1_i64)
    integer(i64), value :: n0_Z_ghost
    real(f64), intent(in) :: Z_ghost(0:n0_Z_ghost - 1_i64)
    integer(i64), value :: n0_h_halo
    real(f64), intent(in) :: h_halo(0:n0_h_halo - 1_i64)
    integer(i64), value :: n0_hu_halo
    real(f64), intent(in) :: hu_halo(0:n0_hu_halo - 1_i64)
    integer(i64), value :: n0_hv_halo
    real(f64), intent(in) :: hv_halo(0:n0_hv_halo - 1_i64)
    integer(i64), value :: n0_hc_halo
    real(f64), intent(in) :: hc_halo(0:n0_hc_halo - 1_i64)
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
          src_PSI, src_Z, h_c, hu_c, hv_c, hc_c, Z_c, h_ghost, hu_ghost &
          , hv_ghost, hc_ghost, Z_ghost, h_halo, hu_halo, hv_halo, &
          hc_halo, Z_halo, h_x, h_y, psi, hx_halo, hy_halo, psi_halo, &
          nodeidc, faceidc, cellidc, cellidf, centerc, normalc, namef, &
          centerf, centerh, vertexn, halofid, order)

  end subroutine bind_c_term_source_srnh_swmhd
  !........................................

  !........................................
  subroutine bind_c_srnh_scheme(hu_l, hu_r, hv_l, hv_r, h_l, h_r, hc_l, &
        hc_r, Z_l, Z_r, n0_normal, normal, mesure, grav, n0_flux, flux &
        ) bind(c)

    implicit none

    real(f64), value :: hu_l
    real(f64), value :: hu_r
    real(f64), value :: hv_l
    real(f64), value :: hv_r
    real(f64), value :: h_l
    real(f64), value :: h_r
    real(f64), value :: hc_l
    real(f64), value :: hc_r
    real(f64), value :: Z_l
    real(f64), value :: Z_r
    integer(i64), value :: n0_normal
    real(f64), intent(in) :: normal(0:n0_normal - 1_i64)
    real(f64), value :: mesure
    real(f64), value :: grav
    integer(i64), value :: n0_flux
    real(f64), intent(inout) :: flux(0:n0_flux - 1_i64)
    real(f64), allocatable :: ninv_0007(:)
    real(f64), allocatable :: w_dif_0007(:)
    real(f64), allocatable :: rmat_0004(:,:)
    integer(i64) :: As_0004
    real(f64) :: p_0004
    real(f64) :: xi_0004
    real(f64) :: u_h_0007
    real(f64) :: v_h_0007
    real(f64) :: c_h_0004
    real(f64) :: un_h_0007
    real(f64) :: vn_h_0007
    real(f64) :: hroe_0007
    real(f64) :: uroe_0007
    real(f64) :: vroe_0007
    real(f64) :: croe_0004
    real(f64) :: uleft_0007
    real(f64) :: vleft_0007
    real(f64) :: uright_0007
    real(f64) :: vright_0007
    real(f64) :: w_lrh_0007
    real(f64) :: w_lrhu_0007
    real(f64) :: w_lrhv_0007
    real(f64) :: w_lrhc_0004
    real(f64) :: w_lrz_0007
    real(f64) :: d_0004
    real(f64) :: sound_0007
    real(f64) :: Q_0004
    real(f64) :: R_0004
    real(f64) :: theta_0004
    real(f64) :: lambda1_0007
    real(f64) :: lambda2_0007
    real(f64) :: lambda3_0007
    real(f64) :: lambda4_0007
    real(f64) :: lambda5_0004
    real(f64) :: alpha1_0004
    real(f64) :: alpha2_0004
    real(f64) :: alpha3_0004
    real(f64) :: beta_0004
    real(f64) :: gamma1_0007
    real(f64) :: gamma2_0007
    real(f64) :: gamma3_0004
    real(f64) :: sigma1_0007
    real(f64) :: sigma2_0007
    real(f64) :: sigma3_0004
    real(f64) :: epsilon_0007
    real(f64) :: sign1_0004
    real(f64) :: sign2_0004
    real(f64) :: sign3_0004
    real(f64) :: sign4_0004
    real(f64) :: sign5_0004
    real(f64) :: hnew_0007
    real(f64) :: unew_0007
    real(f64) :: vnew_0007
    real(f64) :: cnew_0004
    real(f64) :: znew_0007
    real(f64) :: u_hu_0007
    real(f64) :: u_hv_0007
    real(f64) :: u_hc_0004
    real(f64) :: u_z_0007
    real(f64) :: q_s_0007

    allocate(ninv_0007(0:1_i64))
    ninv_0007 = 0.0_f64
    allocate(w_dif_0007(0:4_i64))
    w_dif_0007 = 0.0_f64
    allocate(rmat_0004(0:4_i64, 0:4_i64))
    rmat_0004 = 0.0_f64
    As_0004 = 0_i64
    p_0004 = 0.4_f64
    xi_0004 = 1_i64 / (1_i64 - p_0004)
    ninv_0007 = 0.0_f64
    ninv_0007(0_i64) = (-1_i64) * normal(1_i64)
    ninv_0007(1_i64) = normal(0_i64)
    u_h_0007 = (hu_l / h_l * sqrt(h_l) + hu_r / h_r * sqrt(h_r)) / (sqrt &
          (h_l) + sqrt(h_r))
    v_h_0007 = (hv_l / h_l * sqrt(h_l) + hv_r / h_r * sqrt(h_r)) / (sqrt &
          (h_l) + sqrt(h_r))
    c_h_0004 = (hc_l / h_l * sqrt(h_l) + hc_r / h_r * sqrt(h_r)) / (sqrt &
          (h_l) + sqrt(h_r))
    !uvh =  array([uh, vh])
    un_h_0007 = u_h_0007 * normal(0_i64) + v_h_0007 * normal(1_i64)
    un_h_0007 = un_h_0007 / mesure
    vn_h_0007 = u_h_0007 * ninv_0007(0_i64) + v_h_0007 * ninv_0007(1_i64 &
          )
    vn_h_0007 = vn_h_0007 / mesure
    hroe_0007 = (h_l + h_r) / 2_i64
    uroe_0007 = un_h_0007
    vroe_0007 = vn_h_0007
    croe_0004 = c_h_0004
    uleft_0007 = hu_l * normal(0_i64) + hv_l * normal(1_i64)
    uleft_0007 = uleft_0007 / mesure
    vleft_0007 = hu_l * ninv_0007(0_i64) + hv_l * ninv_0007(1_i64)
    vleft_0007 = vleft_0007 / mesure
    uright_0007 = hu_r * normal(0_i64) + hv_r * normal(1_i64)
    uright_0007 = uright_0007 / mesure
    vright_0007 = hu_r * ninv_0007(0_i64) + hv_r * ninv_0007(1_i64)
    vright_0007 = vright_0007 / mesure
    w_lrh_0007 = (h_l + h_r) / 2_i64
    w_lrhu_0007 = (uleft_0007 + uright_0007) / 2_i64
    w_lrhv_0007 = (vleft_0007 + vright_0007) / 2_i64
    w_lrhc_0004 = (hc_l + hc_r) / 2_i64
    w_lrz_0007 = (Z_l + Z_r) / 2_i64
    w_dif_0007(0_i64) = h_r - h_l
    w_dif_0007(1_i64) = uright_0007 - uleft_0007
    w_dif_0007(2_i64) = vright_0007 - vleft_0007
    w_dif_0007(3_i64) = hc_r - hc_l
    w_dif_0007(4_i64) = Z_r - Z_l
    d_0004 = As_0004 * xi_0004 * (3_i64 * uroe_0007 ** 2_i64 + vroe_0007 &
          ** 2_i64)
    sound_0007 = sqrt(grav * hroe_0007)
    Q_0004 = (-(uroe_0007 ** 2_i64 + 3_i64 * grav * (hroe_0007 + d_0004 &
          ))) / 9_i64
    R_0004 = uroe_0007 * (9_i64 * grav * (2_i64 * hroe_0007 - d_0004) - &
          2_i64 * uroe_0007 ** 2_i64) / 54_i64
    theta_0004 = acos(R_0004 / sqrt(-Q_0004 ** 3_i64))
    !Les valeurs propres
    lambda1_0007 = 2_i64 * sqrt(-Q_0004) * cos(theta_0004 / 3_i64) + &
          2.0_f64 / 3.0_f64 * uroe_0007
    lambda2_0007 = 2_i64 * sqrt(-Q_0004) * cos((theta_0004 + 2_i64 * &
          3.141592653589793_f64) / 3_i64) + 2.0_f64 / 3.0_f64 * &
          uroe_0007
    lambda3_0007 = 2_i64 * sqrt(-Q_0004) * cos((theta_0004 + 4_i64 * &
          3.141592653589793_f64) / 3_i64) + 2.0_f64 / 3.0_f64 * &
          uroe_0007
    lambda4_0007 = uroe_0007
    lambda5_0004 = uroe_0007
    !définition de alpha
    alpha1_0004 = lambda1_0007 - uroe_0007
    alpha2_0004 = lambda2_0007 - uroe_0007
    alpha3_0004 = lambda3_0007 - uroe_0007
    !définition de beta
    beta_0004 = 2_i64 * As_0004 * xi_0004 * vroe_0007 / hroe_0007
    !définition de gamma
    gamma1_0007 = sound_0007 ** 2_i64 - uroe_0007 ** 2_i64 + &
          lambda2_0007 * lambda3_0007 - beta_0004 * alpha2_0004 * &
          alpha3_0004 * vroe_0007
    gamma2_0007 = sound_0007 ** 2_i64 - uroe_0007 ** 2_i64 + &
          lambda1_0007 * lambda3_0007 - beta_0004 * alpha1_0004 * &
          alpha3_0004 * vroe_0007
    gamma3_0004 = sound_0007 ** 2_i64 - uroe_0007 ** 2_i64 + &
          lambda1_0007 * lambda2_0007 - beta_0004 * alpha1_0004 * &
          alpha2_0004 * vroe_0007
    !définition de sigma
    sigma1_0007 = (-alpha1_0004) * alpha2_0004 + alpha2_0004 * &
          alpha3_0004 - alpha1_0004 * alpha3_0004 + alpha1_0004 ** &
          2_i64
    sigma2_0007 = alpha1_0004 * alpha2_0004 + alpha2_0004 * alpha3_0004 &
          - alpha1_0004 * alpha3_0004 - alpha2_0004 ** 2_i64
    sigma3_0004 = alpha1_0004 * alpha2_0004 - alpha2_0004 * alpha3_0004 &
          - alpha1_0004 * alpha3_0004 + alpha3_0004 ** 2_i64
    epsilon_0007 = 1e-10_f64
    if (abs(lambda1_0007) < epsilon_0007) then
      sign1_0004 = 0.0_f64
    else
      sign1_0004 = lambda1_0007 / abs(lambda1_0007)
    end if
    if (abs(lambda2_0007) < epsilon_0007) then
      sign2_0004 = 0.0_f64
    else
      sign2_0004 = lambda2_0007 / abs(lambda2_0007)
    end if
    if (abs(lambda3_0007) < epsilon_0007) then
      sign3_0004 = 0.0_f64
    else
      sign3_0004 = lambda3_0007 / abs(lambda3_0007)
    end if
    if (abs(lambda4_0007) < epsilon_0007) then
      sign4_0004 = 0.0_f64
    else
      sign4_0004 = lambda4_0007 / abs(lambda4_0007)
    end if
    if (abs(lambda5_0004) < epsilon_0007) then
      sign5_0004 = 0.0_f64
    else
      sign5_0004 = lambda5_0004 / abs(lambda5_0004)
    end if
    !1ère colonne
    rmat_0004(0_i64, 0_i64) = sign1_0004 * (gamma1_0007 / sigma1_0007) - &
          sign2_0004 * (gamma2_0007 / sigma2_0007) + sign3_0004 * ( &
          gamma3_0004 / sigma3_0004) + sign5_0004 * (beta_0004 * &
          vroe_0007)
    rmat_0004(0_i64, 1_i64) = lambda1_0007 * sign1_0004 * (gamma1_0007 / &
          sigma1_0007) - lambda2_0007 * sign2_0004 * (gamma2_0007 / &
          sigma2_0007) + lambda3_0007 * sign3_0004 * (gamma3_0004 / &
          sigma3_0004) + sign5_0004 * (beta_0004 * uroe_0007 * &
          vroe_0007)
    rmat_0004(0_i64, 2_i64) = vroe_0007 * sign1_0004 * (gamma1_0007 / &
          sigma1_0007) - vroe_0007 * sign2_0004 * (gamma2_0007 / &
          sigma2_0007) + vroe_0007 * sign3_0004 * (gamma3_0004 / &
          sigma3_0004) - vroe_0007 * sign5_0004 * (1_i64 - beta_0004 * &
          vroe_0007)
    rmat_0004(0_i64, 3_i64) = croe_0004 * sign1_0004 * (gamma1_0007 / &
          sigma1_0007) - croe_0004 * sign2_0004 * (gamma2_0007 / &
          sigma2_0007) + croe_0004 * sign3_0004 * (gamma3_0004 / &
          sigma3_0004) - croe_0004 * sign4_0004 + croe_0004 * &
          sign5_0004 * beta_0004 * vroe_0007
    rmat_0004(0_i64, 4_i64) = (alpha1_0004 ** 2_i64 / sound_0007 ** &
          2_i64 - 1_i64) * sign1_0004 * (gamma1_0007 / sigma1_0007) - ( &
          alpha2_0004 ** 2_i64 / sound_0007 ** 2_i64 - 1_i64) * &
          sign2_0004 * (gamma2_0007 / sigma2_0007) + (alpha3_0004 ** &
          2_i64 / sound_0007 ** 2_i64 - 1_i64) * sign3_0004 * ( &
          gamma3_0004 / sigma3_0004) - sign5_0004 * (beta_0004 * &
          vroe_0007)
    !2ème colonne
    rmat_0004(1_i64, 0_i64) = (-sign1_0004) * (alpha2_0004 + alpha3_0004 &
          ) / sigma1_0007 + sign2_0004 * (alpha1_0004 + alpha3_0004) / &
          sigma2_0007 - sign3_0004 * (alpha1_0004 + alpha2_0004) / &
          sigma3_0004
    rmat_0004(1_i64, 1_i64) = (-lambda1_0007) * sign1_0004 * ( &
          alpha2_0004 + alpha3_0004) / sigma1_0007 + lambda2_0007 * &
          sign2_0004 * (alpha1_0004 + alpha3_0004) / sigma2_0007 - &
          lambda3_0007 * sign3_0004 * (alpha1_0004 + alpha2_0004) / &
          sigma3_0004
    rmat_0004(1_i64, 2_i64) = (-vroe_0007) * sign1_0004 * (alpha2_0004 + &
          alpha3_0004) / sigma1_0007 + vroe_0007 * sign2_0004 * ( &
          alpha1_0004 + alpha3_0004) / sigma2_0007 - vroe_0007 * &
          sign3_0004 * (alpha1_0004 + alpha2_0004) / sigma3_0004
    rmat_0004(1_i64, 3_i64) = (-croe_0004) * sign1_0004 * (alpha2_0004 + &
          alpha3_0004) / sigma1_0007 + croe_0004 * sign2_0004 * ( &
          alpha1_0004 + alpha3_0004) / sigma2_0007 - croe_0004 * &
          sign3_0004 * (alpha1_0004 + alpha2_0004) / sigma3_0004
    rmat_0004(1_i64, 4_i64) = (-(alpha1_0004 ** 2_i64 / sound_0007 ** &
          2_i64 - 1_i64)) * sign1_0004 * (alpha2_0004 + alpha3_0004) / &
          sigma1_0007 + (alpha2_0004 ** 2_i64 / sound_0007 ** 2_i64 - &
          1_i64) * sign2_0004 * (alpha1_0004 + alpha3_0004) / &
          sigma2_0007 - (alpha3_0004 ** 2_i64 / sound_0007 ** 2_i64 - &
          1_i64) * sign3_0004 * (alpha1_0004 + alpha2_0004) / &
          sigma3_0004
    !3ème colonne
    rmat_0004(2_i64, 0_i64) = sign1_0004 * beta_0004 * alpha2_0004 * &
          alpha3_0004 / sigma1_0007 - sign2_0004 * beta_0004 * &
          alpha1_0004 * alpha3_0004 / sigma2_0007 + sign3_0004 * &
          beta_0004 * alpha1_0004 * alpha2_0004 / sigma3_0004 - &
          sign5_0004 * beta_0004
    rmat_0004(2_i64, 1_i64) = lambda1_0007 * sign1_0004 * beta_0004 * &
          alpha2_0004 * alpha3_0004 / sigma1_0007 - lambda2_0007 * &
          sign2_0004 * beta_0004 * alpha1_0004 * alpha3_0004 / &
          sigma2_0007 + lambda3_0007 * sign3_0004 * beta_0004 * &
          alpha1_0004 * alpha2_0004 / sigma3_0004 - sign5_0004 * &
          beta_0004 * uroe_0007
    rmat_0004(2_i64, 2_i64) = vroe_0007 * sign1_0004 * beta_0004 * &
          alpha2_0004 * alpha3_0004 / sigma1_0007 - vroe_0007 * &
          sign2_0004 * beta_0004 * alpha1_0004 * alpha3_0004 / &
          sigma2_0007 + vroe_0007 * sign3_0004 * beta_0004 * &
          alpha1_0004 * alpha2_0004 / sigma3_0004 + sign5_0004 * (1_i64 &
          - beta_0004 * vroe_0007)
    rmat_0004(2_i64, 3_i64) = croe_0004 * sign1_0004 * beta_0004 * &
          alpha2_0004 * alpha3_0004 / sigma1_0007 - croe_0004 * &
          sign2_0004 * beta_0004 * alpha1_0004 * alpha3_0004 / &
          sigma2_0007 + croe_0004 * sign3_0004 * beta_0004 * &
          alpha1_0004 * alpha2_0004 / sigma3_0004 - croe_0004 * &
          sign5_0004 * beta_0004
    rmat_0004(2_i64, 4_i64) = (alpha1_0004 ** 2_i64 / sound_0007 ** &
          2_i64 - 1_i64) * sign1_0004 * beta_0004 * alpha2_0004 * &
          alpha3_0004 / sigma1_0007 - (alpha2_0004 ** 2_i64 / &
          sound_0007 ** 2_i64 - 1_i64) * sign2_0004 * beta_0004 * &
          alpha1_0004 * alpha3_0004 / sigma2_0007 + (alpha3_0004 ** &
          2_i64 / sound_0007 ** 2_i64 - 1_i64) * sign3_0004 * beta_0004 &
          * alpha1_0004 * alpha2_0004 / sigma3_0004 + sign5_0004 * &
          beta_0004
    !4ème colonne
    rmat_0004(3_i64, 0_i64) = 0.0_f64
    rmat_0004(3_i64, 1_i64) = 0.0_f64
    rmat_0004(3_i64, 2_i64) = 0.0_f64
    rmat_0004(3_i64, 3_i64) = sign4_0004
    rmat_0004(3_i64, 4_i64) = 0.0_f64
    !5ème colone
    rmat_0004(4_i64, 0_i64) = sign1_0004 * sound_0007 ** 2_i64 / &
          sigma1_0007 - sign2_0004 * sound_0007 ** 2_i64 / sigma2_0007 &
          + sign3_0004 * sound_0007 ** 2_i64 / sigma3_0004
    rmat_0004(4_i64, 1_i64) = lambda1_0007 * sign1_0004 * sound_0007 ** &
          2_i64 / sigma1_0007 - lambda2_0007 * sign2_0004 * sound_0007 &
          ** 2_i64 / sigma2_0007 + lambda3_0007 * sign3_0004 * &
          sound_0007 ** 2_i64 / sigma3_0004
    rmat_0004(4_i64, 2_i64) = vroe_0007 * sign1_0004 * sound_0007 ** &
          2_i64 / sigma1_0007 - vroe_0007 * sign2_0004 * sound_0007 ** &
          2_i64 / sigma2_0007 + vroe_0007 * sign3_0004 * sound_0007 ** &
          2_i64 / sigma3_0004
    rmat_0004(4_i64, 3_i64) = croe_0004 * sign1_0004 * sound_0007 ** &
          2_i64 / sigma1_0007 - croe_0004 * sign2_0004 * sound_0007 ** &
          2_i64 / sigma2_0007 + croe_0004 * sign3_0004 * sound_0007 ** &
          2_i64 / sigma3_0004
    rmat_0004(4_i64, 4_i64) = (alpha1_0004 ** 2_i64 / sound_0007 ** &
          2_i64 - 1_i64) * sign1_0004 * sound_0007 ** 2_i64 / &
          sigma1_0007 - (alpha2_0004 ** 2_i64 / sound_0007 ** 2_i64 - &
          1_i64) * sign2_0004 * sound_0007 ** 2_i64 / sigma2_0007 + ( &
          alpha3_0004 ** 2_i64 / sound_0007 ** 2_i64 - 1_i64) * &
          sign3_0004 * sound_0007 ** 2_i64 / sigma3_0004
    hnew_0007 = sum(rmat_0004(:, 0_i64) * w_dif_0007(:))
    unew_0007 = sum(rmat_0004(:, 1_i64) * w_dif_0007(:))
    vnew_0007 = sum(rmat_0004(:, 2_i64) * w_dif_0007(:))
    cnew_0004 = sum(rmat_0004(:, 3_i64) * w_dif_0007(:))
    znew_0007 = sum(rmat_0004(:, 4_i64) * w_dif_0007(:))
    u_h_0007 = hnew_0007 / 2_i64
    u_hu_0007 = unew_0007 / 2_i64
    u_hv_0007 = vnew_0007 / 2_i64
    u_hc_0004 = cnew_0004 / 2_i64
    u_z_0007 = znew_0007 / 2_i64
    w_lrh_0007 = w_lrh_0007 - u_h_0007
    w_lrhu_0007 = w_lrhu_0007 - u_hu_0007
    w_lrhv_0007 = w_lrhv_0007 - u_hv_0007
    w_lrhc_0004 = w_lrhc_0004 - u_hc_0004
    w_lrz_0007 = w_lrz_0007 - u_z_0007
    unew_0007 = 0.0_f64
    vnew_0007 = 0.0_f64
    unew_0007 = w_lrhu_0007 * normal(0_i64) + w_lrhv_0007 * (-1_i64) * &
          normal(1_i64)
    unew_0007 = unew_0007 / mesure
    vnew_0007 = w_lrhu_0007 * (-1_i64) * ninv_0007(0_i64) + w_lrhv_0007 &
          * ninv_0007(1_i64)
    vnew_0007 = vnew_0007 / mesure
    w_lrhu_0007 = unew_0007
    w_lrhv_0007 = vnew_0007
    q_s_0007 = normal(0_i64) * unew_0007 + normal(1_i64) * vnew_0007
    flux(0_i64) = q_s_0007
    flux(1_i64) = q_s_0007 * w_lrhu_0007 / w_lrh_0007 + 0.5_f64 * grav * &
          w_lrh_0007 * w_lrh_0007 * normal(0_i64)
    flux(2_i64) = q_s_0007 * w_lrhv_0007 / w_lrh_0007 + 0.5_f64 * grav * &
          w_lrh_0007 * w_lrh_0007 * normal(1_i64)
    flux(3_i64) = q_s_0007 * w_lrhc_0004 / w_lrh_0007
    flux(4_i64) = As_0004 * xi_0004 * normal(0_i64) * unew_0007 * ( &
          unew_0007 ** 2_i64 + vnew_0007 ** 2_i64) / w_lrh_0007 ** &
          3_i64 + As_0004 * xi_0004 * normal(1_i64) * vnew_0007 * ( &
          unew_0007 ** 2_i64 + vnew_0007 ** 2_i64) / w_lrh_0007 ** &
          3_i64
    if (allocated(ninv_0007)) then
      deallocate(ninv_0007)
    end if
    if (allocated(w_dif_0007)) then
      deallocate(w_dif_0007)
    end if
    if (allocated(rmat_0004)) then
      deallocate(rmat_0004)
    end if

  end subroutine bind_c_srnh_scheme
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
    real(f64), allocatable :: ninv_0008(:)
    real(f64), allocatable :: w_dif_0008(:)
    real(f64) :: u_h_0008
    real(f64) :: v_h_0008
    real(f64) :: B1_h_0004
    real(f64) :: B2_h_0004
    real(f64) :: un_h_0008
    real(f64) :: vn_h_0008
    real(f64) :: B1n_h_0004
    real(f64) :: B2n_h_0004
    real(f64) :: hroe_0008
    real(f64) :: uroe_0008
    real(f64) :: vroe_0008
    real(f64) :: B1roe_0004
    real(f64) :: B2roe_0004
    real(f64) :: uleft_0008
    real(f64) :: vleft_0008
    real(f64) :: B1left_0004
    real(f64) :: B2left_0004
    real(f64) :: uright_0008
    real(f64) :: vright_0008
    real(f64) :: B1right_0004
    real(f64) :: B2right_0004
    real(f64) :: B1i_0004
    real(f64) :: B1j_0004
    real(f64) :: absy_0004
    real(f64) :: w_lrh_0008
    real(f64) :: w_lrhu_0008
    real(f64) :: w_lrhv_0008
    real(f64) :: w_lrhB1_0004
    real(f64) :: w_lrhB2_0004
    real(f64) :: w_lrhPSI_0004
    real(f64) :: w_lrz_0008
    real(f64), allocatable, target :: signA_0004(:,:)
    real(f64) :: sound_0008
    real(f64) :: w_0004
    real(f64) :: lambda1_0008
    real(f64) :: lambda2_0008
    real(f64) :: lambda3_0008
    real(f64) :: lambda4_0008
    real(f64) :: epsilon_0008
    real(f64) :: s1_0004
    real(f64) :: pi1_0004
    real(f64) :: s2_0004
    real(f64) :: pi2_0004
    real(f64) :: s3_0004
    real(f64) :: pi3_0004
    real(f64) :: s4_0004
    real(f64) :: pi4_0004
    real(f64) :: gamma1_0008
    real(f64) :: gamma2_0008
    real(f64) :: sigma1_0008
    real(f64) :: sigma2_0008
    real(f64) :: mu1_0004
    real(f64) :: mu2_0004
    integer(i64) :: ann_0004
    real(f64), pointer :: smmat_0004(:,:)
    real(f64) :: hnew_0008
    real(f64) :: unew_0008
    real(f64) :: vnew_0008
    real(f64) :: B1new_0004
    real(f64) :: B2new_0004
    real(f64) :: znew_0008
    real(f64) :: Pnew_0004
    real(f64) :: u_hu_0008
    real(f64) :: u_hv_0008
    real(f64) :: u_hP_0004
    real(f64) :: u_hB1_0004
    real(f64) :: u_hB2_0004
    real(f64) :: u_z_0008
    real(f64) :: w_lrhP_0004
    real(f64) :: w_hP_0004
    real(f64) :: mw_hB1_0004
    real(f64) :: mhP_0004
    real(f64), allocatable :: norm_0004(:)
    real(f64) :: q_s_0008
    real(f64) :: p_s_0004
    real(f64) :: Flux_B1psi_0004
    real(f64) :: Flux_B2psi_0004
    real(f64) :: Flux_hPpsi_0004
    integer(i64) :: i_0004

    allocate(ninv_0008(0:1_i64))
    ninv_0008 = 0.0_f64
    allocate(w_dif_0008(0:5_i64))
    w_dif_0008 = 0.0_f64
    ninv_0008(0_i64) = (-1_i64) * normal(1_i64)
    ninv_0008(1_i64) = normal(0_i64)
    u_h_0008 = (hu_l / h_l * sqrt(h_l) + hu_r / h_r * sqrt(h_r)) / (sqrt &
          (h_l) + sqrt(h_r))
    v_h_0008 = (hv_l / h_l * sqrt(h_l) + hv_r / h_r * sqrt(h_r)) / (sqrt &
          (h_l) + sqrt(h_r))
    B1_h_0004 = (hB1_l / h_l * sqrt(h_l) + hB1_r / h_r * sqrt(h_r)) / ( &
          sqrt(h_l) + sqrt(h_r))
    B2_h_0004 = (hB2_l / h_l * sqrt(h_l) + hB2_r / h_r * sqrt(h_r)) / ( &
          sqrt(h_l) + sqrt(h_r))
    !uvh =  array([uh, vh])
    un_h_0008 = u_h_0008 * normal(0_i64) + v_h_0008 * normal(1_i64)
    un_h_0008 = un_h_0008 / mesure
    vn_h_0008 = u_h_0008 * ninv_0008(0_i64) + v_h_0008 * ninv_0008(1_i64 &
          )
    vn_h_0008 = vn_h_0008 / mesure
    B1n_h_0004 = B1_h_0004 * normal(0_i64) + B2_h_0004 * normal(1_i64)
    B1n_h_0004 = B1n_h_0004 / mesure
    B2n_h_0004 = B1_h_0004 * ninv_0008(0_i64) + B2_h_0004 * ninv_0008( &
          1_i64)
    B2n_h_0004 = B2n_h_0004 / mesure
    hroe_0008 = (h_l + h_r) / 2_i64
    uroe_0008 = un_h_0008
    vroe_0008 = vn_h_0008
    B1roe_0004 = B1n_h_0004
    B2roe_0004 = B2n_h_0004
    uleft_0008 = hu_l * normal(0_i64) + hv_l * normal(1_i64)
    uleft_0008 = uleft_0008 / mesure
    vleft_0008 = hu_l * ninv_0008(0_i64) + hv_l * ninv_0008(1_i64)
    vleft_0008 = vleft_0008 / mesure
    B1left_0004 = hB1_l * normal(0_i64) + hB2_l * normal(1_i64)
    B1left_0004 = B1left_0004 / mesure
    B2left_0004 = hB1_l * ninv_0008(0_i64) + hB2_l * ninv_0008(1_i64)
    B2left_0004 = B2left_0004 / mesure
    uright_0008 = hu_r * normal(0_i64) + hv_r * normal(1_i64)
    uright_0008 = uright_0008 / mesure
    vright_0008 = hu_r * ninv_0008(0_i64) + hv_r * ninv_0008(1_i64)
    vright_0008 = vright_0008 / mesure
    B1right_0004 = hB1_r * normal(0_i64) + hB2_r * normal(1_i64)
    B1right_0004 = B1right_0004 / mesure
    B2right_0004 = hB1_r * ninv_0008(0_i64) + hB2_r * ninv_0008(1_i64)
    B2right_0004 = B2right_0004 / mesure
    B1i_0004 = (hB1c_l * normal(0_i64) + hB2c_l * normal(1_i64)) / &
          mesure
    B1j_0004 = (hB1c_r * normal(0_i64) + hB1c_r * normal(1_i64)) / &
          mesure
    absy_0004 = (B1i_0004 + B1j_0004) / 2_i64
    w_lrh_0008 = (h_l + h_r) / 2_i64
    w_lrhu_0008 = (uleft_0008 + uright_0008) / 2_i64
    w_lrhv_0008 = (vleft_0008 + vright_0008) / 2_i64
    w_lrhB1_0004 = (B1left_0004 + B1right_0004) / 2_i64
    w_lrhB2_0004 = (B2left_0004 + B2right_0004) / 2_i64
    w_lrhPSI_0004 = (hPSI_l + hPSI_r) / 2_i64
    w_lrz_0008 = (Z_l + Z_r) / 2_i64
    w_dif_0008(0_i64) = h_r - h_l
    w_dif_0008(1_i64) = uright_0008 - uleft_0008
    w_dif_0008(2_i64) = vright_0008 - vleft_0008
    w_dif_0008(3_i64) = B1right_0004 - B1left_0004
    w_dif_0008(4_i64) = B2right_0004 - B2left_0004
    w_dif_0008(5_i64) = Z_r - Z_l
    allocate(signA_0004(0:5_i64, 0:5_i64))
    signA_0004 = 0.0_f64
    sound_0008 = sqrt(grav * hroe_0008)
    B1roe_0004 = absy_0004 / hroe_0008
    w_0004 = sqrt(B1roe_0004 * B1roe_0004 + grav * hroe_0008)
    lambda1_0008 = uroe_0008 - w_0004
    lambda2_0008 = uroe_0008 - B1roe_0004
    lambda3_0008 = uroe_0008 + B1roe_0004
    lambda4_0008 = uroe_0008 + w_0004
    !cpsi = max(fabs(lambda1) , fabs(lambda4))
    epsilon_0008 = 1e-10_f64
    if (abs(lambda1_0008) < epsilon_0008) then
      s1_0004 = 0.0_f64
      pi1_0004 = 0.0_f64
    else
      s1_0004 = lambda1_0008 / abs(lambda1_0008)
      pi1_0004 = s1_0004 / lambda1_0008
    end if
    if (abs(lambda2_0008) < epsilon_0008) then
      s2_0004 = 0.0_f64
      pi2_0004 = 0.0_f64
    else
      s2_0004 = lambda2_0008 / abs(lambda2_0008)
      pi2_0004 = 1.0_f64 / abs(lambda2_0008)
    end if
    if (abs(lambda3_0008) < epsilon_0008) then
      s3_0004 = 0.0_f64
      pi3_0004 = 0.0_f64
    else
      s3_0004 = lambda3_0008 / abs(lambda3_0008)
      pi3_0004 = 1.0_f64 / abs(lambda3_0008)
    end if
    if (abs(lambda4_0008) < epsilon_0008) then
      s4_0004 = 0.0_f64
      pi4_0004 = 0.0_f64
    else
      s4_0004 = lambda4_0008 / abs(lambda4_0008)
      pi4_0004 = 1.0_f64 / abs(lambda4_0008)
    end if
    gamma1_0008 = vroe_0008 + B2roe_0004
    gamma2_0008 = vroe_0008 - B2roe_0004
    sigma1_0008 = vroe_0008 * (s1_0004 * lambda4_0008 - s4_0004 * &
          lambda1_0008) - w_0004 * (s2_0004 * gamma1_0008 + s3_0004 * &
          gamma2_0008)
    sigma2_0008 = B2roe_0004 * (s1_0004 * lambda4_0008 - s4_0004 * &
          lambda1_0008) - w_0004 * (s2_0004 * gamma1_0008 - s3_0004 * &
          gamma2_0008)
    if (abs(lambda2_0008) < epsilon_0008 .and. abs(lambda3_0008) < &
          epsilon_0008) then
      mu1_0004 = B1roe_0004 * vroe_0008 * pi1_0004 / w_0004 - B1roe_0004 &
            * vroe_0008 * pi4_0004 / w_0004
      mu2_0004 = B1roe_0004 * B2roe_0004 * pi1_0004 / w_0004 - &
            B1roe_0004 * B2roe_0004 * pi4_0004 / w_0004
      ann_0004 = 0_i64
    else
      mu1_0004 = B1roe_0004 * vroe_0008 * pi1_0004 / w_0004 - B1roe_0004 &
            * vroe_0008 * pi4_0004 / w_0004 - 0.5_f64 * (gamma1_0008 * &
            pi2_0004 - gamma2_0008 * pi3_0004)
      mu2_0004 = B1roe_0004 * B2roe_0004 * pi1_0004 / w_0004 - &
            B1roe_0004 * B2roe_0004 * pi4_0004 / w_0004 - 0.5_f64 * ( &
            gamma1_0008 * pi2_0004 + gamma2_0008 * pi3_0004)
      ann_0004 = 0_i64
    end if
    !1ère colonne de la matrice A
    signA_0004(0_i64, 0_i64) = (s1_0004 * lambda4_0008 - s4_0004 * &
          lambda1_0008) / (2_i64 * w_0004)
    signA_0004(0_i64, 1_i64) = lambda1_0008 * lambda4_0008 * (s1_0004 - &
          s4_0004) / (2_i64 * w_0004)
    signA_0004(0_i64, 2_i64) = sigma1_0008 / (2_i64 * w_0004)
    signA_0004(0_i64, 3_i64) = 0.0_f64
    signA_0004(0_i64, 4_i64) = sigma2_0008 / (2_i64 * w_0004)
    signA_0004(0_i64, 5_i64) = 0.0_f64
    !2ème colonne de la matrice A
    signA_0004(1_i64, 0_i64) = (s4_0004 - s1_0004) / (2_i64 * w_0004)
    signA_0004(1_i64, 1_i64) = (s4_0004 * lambda4_0008 - s1_0004 * &
          lambda1_0008) / (2_i64 * w_0004)
    signA_0004(1_i64, 2_i64) = vroe_0008 * (s4_0004 - s1_0004) / (2_i64 &
          * w_0004)
    signA_0004(1_i64, 3_i64) = 0.0_f64
    signA_0004(1_i64, 4_i64) = B2roe_0004 * (s4_0004 - s1_0004) / (2_i64 &
          * w_0004)
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
    signA_0004(5_i64, 0_i64) = sound_0008 ** 2_i64 * (pi4_0004 - &
          pi1_0004) / (2_i64 * w_0004)
    signA_0004(5_i64, 1_i64) = sound_0008 ** 2_i64 * (s4_0004 - s1_0004 &
          ) / (2_i64 * w_0004)
    signA_0004(5_i64, 2_i64) = sound_0008 ** 2_i64 * vroe_0008 * ( &
          pi4_0004 - pi1_0004) / (2_i64 * w_0004)
    signA_0004(5_i64, 3_i64) = 0.0_f64
    signA_0004(5_i64, 4_i64) = sound_0008 ** 2_i64 * B2roe_0004 * ( &
          pi4_0004 - pi1_0004) / (2_i64 * w_0004)
    signA_0004(5_i64, 5_i64) = 0.0_f64
    smmat_0004(0:, 0:) => signA_0004
    hnew_0008 = 0.0_f64
    unew_0008 = 0.0_f64
    vnew_0008 = 0.0_f64
    B1new_0004 = 0.0_f64
    B2new_0004 = 0.0_f64
    znew_0008 = 0.0_f64
    do i_0004 = 0_i64, 5_i64, 1_i64
      hnew_0008 = hnew_0008 + smmat_0004(i_0004, 0_i64) * w_dif_0008( &
            i_0004)
      unew_0008 = unew_0008 + smmat_0004(i_0004, 1_i64) * w_dif_0008( &
            i_0004)
      vnew_0008 = vnew_0008 + smmat_0004(i_0004, 2_i64) * w_dif_0008( &
            i_0004)
      B1new_0004 = B1new_0004 + smmat_0004(i_0004, 3_i64) * w_dif_0008( &
            i_0004)
      B2new_0004 = B2new_0004 + smmat_0004(i_0004, 4_i64) * w_dif_0008( &
            i_0004)
      znew_0008 = znew_0008 + smmat_0004(i_0004, 5_i64) * w_dif_0008( &
            i_0004)
    end do
    Pnew_0004 = cpsi * (B1right_0004 - B1left_0004)
    u_h_0008 = hnew_0008 / 2_i64
    u_hu_0008 = unew_0008 / 2_i64
    u_hv_0008 = vnew_0008 / 2_i64
    u_hP_0004 = Pnew_0004 / 2_i64
    u_hB1_0004 = B1new_0004 / 2_i64
    u_hB2_0004 = B2new_0004 / 2_i64
    u_z_0008 = znew_0008 / 2_i64
    w_lrh_0008 = w_lrh_0008 - u_h_0008
    w_lrhu_0008 = w_lrhu_0008 - u_hu_0008
    w_lrhv_0008 = w_lrhv_0008 - u_hv_0008
    w_lrhP_0004 = w_lrhPSI_0004 - u_hP_0004
    w_lrhB1_0004 = absy_0004
    w_lrhB2_0004 = w_lrhB2_0004 - u_hB2_0004
    w_lrz_0008 = w_lrz_0008 - u_z_0008
    w_hP_0004 = hPSI_r - hPSI_l
    mw_hB1_0004 = w_lrhB1_0004 / 2_i64 - 1_i64 / (2_i64 * cpsi) * &
          w_hP_0004
    mhP_0004 = (hPSI_r - hPSI_l) / 2_i64 - cpsi * w_dif_0008(3_i64) / &
          2_i64
    unew_0008 = 0.0_f64
    vnew_0008 = 0.0_f64
    B1new_0004 = 0.0_f64
    B2new_0004 = 0.0_f64
    unew_0008 = w_lrhu_0008 * normal(0_i64) - w_lrhv_0008 * normal(1_i64 &
          )
    unew_0008 = unew_0008 / mesure
    vnew_0008 = w_lrhu_0008 * normal(1_i64) + w_lrhv_0008 * normal(0_i64 &
          )
    vnew_0008 = vnew_0008 / mesure
    B1new_0004 = w_lrhB1_0004 * normal(0_i64) - w_lrhB2_0004 * normal( &
          1_i64)
    B1new_0004 = B1new_0004 / mesure
    B2new_0004 = w_lrhB1_0004 * normal(1_i64) + w_lrhB2_0004 * normal( &
          0_i64)
    B2new_0004 = B2new_0004 / mesure
    w_lrhu_0008 = unew_0008
    w_lrhv_0008 = vnew_0008
    w_lrhB1_0004 = B1new_0004
    w_lrhB2_0004 = B2new_0004
    allocate(norm_0004(0:size(normal, kind=i64) - 1_i64))
    norm_0004 = normal / mesure
    q_s_0008 = normal(0_i64) * unew_0008 + normal(1_i64) * vnew_0008
    p_s_0004 = normal(0_i64) * B1new_0004 + normal(1_i64) * B2new_0004
    Flux_B1psi_0004 = mhP_0004 * norm_0004(0_i64) * mesure
    Flux_B2psi_0004 = mhP_0004 * norm_0004(1_i64) * mesure
    Flux_hPpsi_0004 = cpsi * cpsi * mw_hB1_0004 * mesure
    flux(0_i64) = q_s_0008
    flux(1_i64) = q_s_0008 * w_lrhu_0008 / w_lrh_0008 + 0.5_f64 * grav * &
          w_lrh_0008 * w_lrh_0008 * normal(0_i64) - p_s_0004 * &
          w_lrhB1_0004 / w_lrh_0008
    flux(2_i64) = q_s_0008 * w_lrhv_0008 / w_lrh_0008 + 0.5_f64 * grav * &
          w_lrh_0008 * w_lrh_0008 * normal(1_i64) - p_s_0004 * &
          w_lrhB2_0004 / w_lrh_0008
    flux(3_i64) = (w_lrhv_0008 * w_lrhB1_0004 / w_lrh_0008 - w_lrhu_0008 &
          * w_lrhB2_0004 / w_lrh_0008) * normal(1_i64) + &
          Flux_B1psi_0004
    flux(4_i64) = (w_lrhu_0008 * w_lrhB2_0004 / w_lrh_0008 - w_lrhv_0008 &
          * w_lrhB1_0004 / w_lrh_0008) * normal(0_i64) + &
          Flux_B2psi_0004
    flux(5_i64) = Flux_hPpsi_0004
    flux(6_i64) = 0_i64
    if (allocated(ninv_0008)) then
      deallocate(ninv_0008)
    end if
    if (allocated(w_dif_0008)) then
      deallocate(w_dif_0008)
    end if
    if (allocated(signA_0004)) then
      deallocate(signA_0004)
    end if
    if (allocated(norm_0004)) then
      deallocate(norm_0004)
    end if

  end subroutine bind_c_srnh_scheme_mhd
  !........................................

  !........................................
  subroutine bind_c_explicitscheme_convective_sw(n0_rez_h, rez_h, &
        n0_rez_hu, rez_hu, n0_rez_hv, rez_hv, n0_rez_hc, rez_hc, &
        n0_rez_Z, rez_Z, n0_h_c, h_c, n0_hu_c, hu_c, n0_hv_c, hv_c, &
        n0_hc_c, hc_c, n0_Z_c, Z_c, n0_h_ghost, h_ghost, n0_hu_ghost, &
        hu_ghost, n0_hv_ghost, hv_ghost, n0_hc_ghost, hc_ghost, &
        n0_Z_ghost, Z_ghost, n0_h_halo, h_halo, n0_hu_halo, hu_halo, &
        n0_hv_halo, hv_halo, n0_hc_halo, hc_halo, n0_Z_halo, Z_halo, &
        n0_h_x, h_x, n0_h_y, h_y, n0_hx_halo, hx_halo, n0_hy_halo, &
        hy_halo, n0_hc_x, hc_x, n0_hc_y, hc_y, n0_hcx_halo, hcx_halo, &
        n0_hcy_halo, hcy_halo, n0_psi, psi, n0_psi_halo, psi_halo, &
        n0_centerc, n1_centerc, centerc, n0_centerf, n1_centerf, &
        centerf, n0_centerh, n1_centerh, centerh, n0_centerg, &
        n1_centerg, centerg, n0_cellidf, n1_cellidf, cellidf, &
        n0_mesuref, mesuref, n0_normalf, n1_normalf, normalf, &
        n0_halofid, halofid, n0_innerfaces, innerfaces, n0_halofaces, &
        halofaces, n0_boundaryfaces, boundaryfaces, order) bind(c)

    implicit none

    integer(i64), value :: n0_rez_h
    real(f64), intent(inout) :: rez_h(0:n0_rez_h - 1_i64)
    integer(i64), value :: n0_rez_hu
    real(f64), intent(inout) :: rez_hu(0:n0_rez_hu - 1_i64)
    integer(i64), value :: n0_rez_hv
    real(f64), intent(inout) :: rez_hv(0:n0_rez_hv - 1_i64)
    integer(i64), value :: n0_rez_hc
    real(f64), intent(inout) :: rez_hc(0:n0_rez_hc - 1_i64)
    integer(i64), value :: n0_rez_Z
    real(f64), intent(inout) :: rez_Z(0:n0_rez_Z - 1_i64)
    integer(i64), value :: n0_h_c
    real(f64), intent(in) :: h_c(0:n0_h_c - 1_i64)
    integer(i64), value :: n0_hu_c
    real(f64), intent(in) :: hu_c(0:n0_hu_c - 1_i64)
    integer(i64), value :: n0_hv_c
    real(f64), intent(in) :: hv_c(0:n0_hv_c - 1_i64)
    integer(i64), value :: n0_hc_c
    real(f64), intent(in) :: hc_c(0:n0_hc_c - 1_i64)
    integer(i64), value :: n0_Z_c
    real(f64), intent(in) :: Z_c(0:n0_Z_c - 1_i64)
    integer(i64), value :: n0_h_ghost
    real(f64), intent(in) :: h_ghost(0:n0_h_ghost - 1_i64)
    integer(i64), value :: n0_hu_ghost
    real(f64), intent(in) :: hu_ghost(0:n0_hu_ghost - 1_i64)
    integer(i64), value :: n0_hv_ghost
    real(f64), intent(in) :: hv_ghost(0:n0_hv_ghost - 1_i64)
    integer(i64), value :: n0_hc_ghost
    real(f64), intent(in) :: hc_ghost(0:n0_hc_ghost - 1_i64)
    integer(i64), value :: n0_Z_ghost
    real(f64), intent(in) :: Z_ghost(0:n0_Z_ghost - 1_i64)
    integer(i64), value :: n0_h_halo
    real(f64), intent(in) :: h_halo(0:n0_h_halo - 1_i64)
    integer(i64), value :: n0_hu_halo
    real(f64), intent(in) :: hu_halo(0:n0_hu_halo - 1_i64)
    integer(i64), value :: n0_hv_halo
    real(f64), intent(in) :: hv_halo(0:n0_hv_halo - 1_i64)
    integer(i64), value :: n0_hc_halo
    real(f64), intent(in) :: hc_halo(0:n0_hc_halo - 1_i64)
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
    integer(i64), value :: n0_hc_x
    real(f64), intent(in) :: hc_x(0:n0_hc_x - 1_i64)
    integer(i64), value :: n0_hc_y
    real(f64), intent(in) :: hc_y(0:n0_hc_y - 1_i64)
    integer(i64), value :: n0_hcx_halo
    real(f64), intent(in) :: hcx_halo(0:n0_hcx_halo - 1_i64)
    integer(i64), value :: n0_hcy_halo
    real(f64), intent(in) :: hcy_halo(0:n0_hcy_halo - 1_i64)
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
    integer(i64), value :: order

    call explicitscheme_convective_SW(rez_h, rez_hu, rez_hv, rez_hc, &
          rez_Z, h_c, hu_c, hv_c, hc_c, Z_c, h_ghost, hu_ghost, &
          hv_ghost, hc_ghost, Z_ghost, h_halo, hu_halo, hv_halo, &
          hc_halo, Z_halo, h_x, h_y, hx_halo, hy_halo, hc_x, hc_y, &
          hcx_halo, hcy_halo, psi, psi_halo, centerc, centerf, centerh, &
          centerg, cellidf, mesuref, normalf, halofid, innerfaces, &
          halofaces, boundaryfaces, order)

  end subroutine bind_c_explicitscheme_convective_sw
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
        n0_boundaryfaces, boundaryfaces, order, cpsi, n0_hB1_cst, &
        hB1_cst, n0_hB2_cst, hB2_cst) bind(c)

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
          boundaryfaces, order, cpsi, hB1_cst, hB2_cst)

  end subroutine bind_c_explicitscheme_convective_swmhd
  !........................................

  !........................................
  subroutine bind_c_term_coriolis_sw(n0_hu_c, hu_c, n0_hv_c, hv_c, &
        n0_corio_hu, corio_hu, n0_corio_hv, corio_hv, f_c) bind(c)

    implicit none

    integer(i64), value :: n0_hu_c
    real(f64), intent(in) :: hu_c(0:n0_hu_c - 1_i64)
    integer(i64), value :: n0_hv_c
    real(f64), intent(in) :: hv_c(0:n0_hv_c - 1_i64)
    integer(i64), value :: n0_corio_hu
    real(f64), intent(inout) :: corio_hu(0:n0_corio_hu - 1_i64)
    integer(i64), value :: n0_corio_hv
    real(f64), intent(inout) :: corio_hv(0:n0_corio_hv - 1_i64)
    real(f64), value :: f_c

    call term_coriolis_SW(hu_c, hv_c, corio_hu, corio_hv, f_c)

  end subroutine bind_c_term_coriolis_sw
  !........................................

  !........................................
  subroutine bind_c_term_friction_sw(n0_h_c, h_c, n0_hu_c, hu_c, n0_hv_c &
        , hv_c, grav, eta, time) bind(c)

    implicit none

    integer(i64), value :: n0_h_c
    real(f64), intent(in) :: h_c(0:n0_h_c - 1_i64)
    integer(i64), value :: n0_hu_c
    real(f64), intent(inout) :: hu_c(0:n0_hu_c - 1_i64)
    integer(i64), value :: n0_hv_c
    real(f64), intent(inout) :: hv_c(0:n0_hv_c - 1_i64)
    real(f64), value :: grav
    real(f64), value :: eta
    real(f64), value :: time

    call term_friction_SW(h_c, hu_c, hv_c, grav, eta, time)

  end subroutine bind_c_term_friction_sw
  !........................................

  !........................................
  subroutine bind_c_term_wind_sw(uwind, vwind, TAUXWX, TAUXWY) bind(c) 

    implicit none

    real(f64), value :: uwind
    real(f64), value :: vwind
    real(f64), intent(out) :: TAUXWX
    real(f64), intent(out) :: TAUXWY

    call term_wind_SW(uwind, vwind, TAUXWX = TAUXWX, TAUXWY = TAUXWY)

  end subroutine bind_c_term_wind_sw
  !........................................

end module bind_c_tools
