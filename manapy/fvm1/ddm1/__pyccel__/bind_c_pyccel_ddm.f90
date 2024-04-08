module bind_c_pyccel_ddm

  use pyccel_ddm, only: create_cell_faceid
  use pyccel_ddm, only: variables_3d
  use pyccel_ddm, only: compute_K
  use pyccel_ddm, only: create_info_2dfaces
  use pyccel_ddm, only: variables
  use pyccel_ddm, only: create_NormalFacesOfCell
  use pyccel_ddm, only: create_3dfaces
  use pyccel_ddm, only: create_info_3dfaces
  use pyccel_ddm, only: create_2dfaces
  use pyccel_ddm, only: create_cellsOfFace
  use pyccel_ddm, only: face_gradient_info_2d
  use pyccel_ddm, only: Compute_3dcentervolumeOfCell
  use pyccel_ddm, only: Compute_2dcentervolumeOfCell

  use, intrinsic :: ISO_C_Binding, only : i64 => C_INT64_T , f64 => &
        C_DOUBLE
  implicit none

  contains

  !........................................
  subroutine bind_c_compute_k(n0_cellid, n1_cellid, cellid, n0_namef, &
        namef, n0_ghostcenterf, n1_ghostcenterf, ghostcenterf, &
        n0_centerc, n1_centerc, centerc, n0_K, n1_K, K, dim) bind(c)

    implicit none

    integer(i64), value :: n0_cellid
    integer(i64), value :: n1_cellid
    integer(i64), intent(in) :: cellid(0:n1_cellid - 1_i64,0:n0_cellid - &
          1_i64)
    integer(i64), value :: n0_namef
    integer(i64), intent(in) :: namef(0:n0_namef - 1_i64)
    integer(i64), value :: n0_ghostcenterf
    integer(i64), value :: n1_ghostcenterf
    real(f64), intent(in) :: ghostcenterf(0:n1_ghostcenterf - 1_i64,0: &
          n0_ghostcenterf - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_K
    integer(i64), value :: n1_K
    real(f64), intent(inout) :: K(0:n1_K - 1_i64,0:n0_K - 1_i64)
    integer(i64), value :: dim

    call compute_K(cellid, namef, ghostcenterf, centerc, K, dim)

  end subroutine bind_c_compute_k
  !........................................

  !........................................
  subroutine bind_c_create_info_2dfaces(n0_cellid, n1_cellid, cellid, &
        n0_nodeid, n1_nodeid, nodeid, n0_namen, namen, n0_vertex, &
        n1_vertex, vertex, n0_centerc, n1_centerc, centerc, nbfaces, &
        n0_normalf, n1_normalf, normalf, n0_mesuref, mesuref, &
        n0_centerf, n1_centerf, centerf, n0_namef, namef) bind(c)

    implicit none

    integer(i64), value :: n0_cellid
    integer(i64), value :: n1_cellid
    integer(i64), intent(in) :: cellid(0:n1_cellid - 1_i64,0:n0_cellid - &
          1_i64)
    integer(i64), value :: n0_nodeid
    integer(i64), value :: n1_nodeid
    integer(i64), intent(in) :: nodeid(0:n1_nodeid - 1_i64,0:n0_nodeid - &
          1_i64)
    integer(i64), value :: n0_namen
    integer(i64), intent(in) :: namen(0:n0_namen - 1_i64)
    integer(i64), value :: n0_vertex
    integer(i64), value :: n1_vertex
    real(f64), intent(in) :: vertex(0:n1_vertex - 1_i64,0:n0_vertex - &
          1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: nbfaces
    integer(i64), value :: n0_normalf
    integer(i64), value :: n1_normalf
    real(f64), intent(inout) :: normalf(0:n1_normalf - 1_i64,0: &
          n0_normalf - 1_i64)
    integer(i64), value :: n0_mesuref
    real(f64), intent(inout) :: mesuref(0:n0_mesuref - 1_i64)
    integer(i64), value :: n0_centerf
    integer(i64), value :: n1_centerf
    real(f64), intent(inout) :: centerf(0:n1_centerf - 1_i64,0: &
          n0_centerf - 1_i64)
    integer(i64), value :: n0_namef
    integer(i64), intent(inout) :: namef(0:n0_namef - 1_i64)

    call create_info_2dfaces(cellid, nodeid, namen, vertex, centerc, &
          nbfaces, normalf, mesuref, centerf, namef)

  end subroutine bind_c_create_info_2dfaces
  !........................................

  !........................................
  subroutine bind_c_create_info_3dfaces(n0_cellid, n1_cellid, cellid, &
        n0_nodeid, n1_nodeid, nodeid, n0_namen, namen, n0_vertex, &
        n1_vertex, vertex, n0_centerc, n1_centerc, centerc, nbfaces, &
        n0_normalf, n1_normalf, normalf, n0_mesuref, mesuref, &
        n0_centerf, n1_centerf, centerf, n0_namef, namef) bind(c)

    implicit none

    integer(i64), value :: n0_cellid
    integer(i64), value :: n1_cellid
    integer(i64), intent(in) :: cellid(0:n1_cellid - 1_i64,0:n0_cellid - &
          1_i64)
    integer(i64), value :: n0_nodeid
    integer(i64), value :: n1_nodeid
    integer(i64), intent(in) :: nodeid(0:n1_nodeid - 1_i64,0:n0_nodeid - &
          1_i64)
    integer(i64), value :: n0_namen
    integer(i64), intent(in) :: namen(0:n0_namen - 1_i64)
    integer(i64), value :: n0_vertex
    integer(i64), value :: n1_vertex
    real(f64), intent(in) :: vertex(0:n1_vertex - 1_i64,0:n0_vertex - &
          1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: nbfaces
    integer(i64), value :: n0_normalf
    integer(i64), value :: n1_normalf
    real(f64), intent(inout) :: normalf(0:n1_normalf - 1_i64,0: &
          n0_normalf - 1_i64)
    integer(i64), value :: n0_mesuref
    real(f64), intent(inout) :: mesuref(0:n0_mesuref - 1_i64)
    integer(i64), value :: n0_centerf
    integer(i64), value :: n1_centerf
    real(f64), intent(inout) :: centerf(0:n1_centerf - 1_i64,0: &
          n0_centerf - 1_i64)
    integer(i64), value :: n0_namef
    integer(i64), intent(inout) :: namef(0:n0_namef - 1_i64)

    call create_info_3dfaces(cellid, nodeid, namen, vertex, centerc, &
          nbfaces, normalf, mesuref, centerf, namef)

  end subroutine bind_c_create_info_3dfaces
  !........................................

  !........................................
  subroutine bind_c_compute_2dcentervolumeofcell(n0_nodeid, n1_nodeid, &
        nodeid, n0_vertex, n1_vertex, vertex, nbelements, n0_center, &
        n1_center, center, n0_volume, volume) bind(c)

    implicit none

    integer(i64), value :: n0_nodeid
    integer(i64), value :: n1_nodeid
    integer(i64), intent(inout) :: nodeid(0:n1_nodeid - 1_i64,0: &
          n0_nodeid - 1_i64)
    integer(i64), value :: n0_vertex
    integer(i64), value :: n1_vertex
    real(f64), intent(in) :: vertex(0:n1_vertex - 1_i64,0:n0_vertex - &
          1_i64)
    integer(i64), value :: nbelements
    integer(i64), value :: n0_center
    integer(i64), value :: n1_center
    real(f64), intent(inout) :: center(0:n1_center - 1_i64,0:n0_center - &
          1_i64)
    integer(i64), value :: n0_volume
    real(f64), intent(inout) :: volume(0:n0_volume - 1_i64)

    call Compute_2dcentervolumeOfCell(nodeid, vertex, nbelements, center &
          , volume)

  end subroutine bind_c_compute_2dcentervolumeofcell
  !........................................

  !........................................
  subroutine bind_c_compute_3dcentervolumeofcell(n0_nodeid, n1_nodeid, &
        nodeid, n0_vertex, n1_vertex, vertex, nbelements, n0_center, &
        n1_center, center, n0_volume, volume) bind(c)

    implicit none

    integer(i64), value :: n0_nodeid
    integer(i64), value :: n1_nodeid
    integer(i64), intent(in) :: nodeid(0:n1_nodeid - 1_i64,0:n0_nodeid - &
          1_i64)
    integer(i64), value :: n0_vertex
    integer(i64), value :: n1_vertex
    real(f64), intent(in) :: vertex(0:n1_vertex - 1_i64,0:n0_vertex - &
          1_i64)
    integer(i64), value :: nbelements
    integer(i64), value :: n0_center
    integer(i64), value :: n1_center
    real(f64), intent(inout) :: center(0:n1_center - 1_i64,0:n0_center - &
          1_i64)
    integer(i64), value :: n0_volume
    real(f64), intent(inout) :: volume(0:n0_volume - 1_i64)

    call Compute_3dcentervolumeOfCell(nodeid, vertex, nbelements, center &
          , volume)

  end subroutine bind_c_compute_3dcentervolumeofcell
  !........................................

  !........................................
  subroutine bind_c_create_cellsofface(n0_faceid, n1_faceid, faceid, &
        nbelements, nbfaces, n0_cellid, n1_cellid, cellid, dim) bind(c &
        )

    implicit none

    integer(i64), value :: n0_faceid
    integer(i64), value :: n1_faceid
    integer(i64), intent(inout) :: faceid(0:n1_faceid - 1_i64,0: &
          n0_faceid - 1_i64)
    integer(i64), value :: nbelements
    integer(i64), value :: nbfaces
    integer(i64), value :: n0_cellid
    integer(i64), value :: n1_cellid
    integer(i64), intent(inout) :: cellid(0:n1_cellid - 1_i64,0: &
          n0_cellid - 1_i64)
    integer(i64), value :: dim

    call create_cellsOfFace(faceid, nbelements, nbfaces, cellid, dim)

  end subroutine bind_c_create_cellsofface
  !........................................

  !........................................
  subroutine bind_c_create_2dfaces(n0_nodeidc, n1_nodeidc, nodeidc, &
        nbelements, n0_faces, n1_faces, faces, n0_cellf, n1_cellf, &
        cellf) bind(c)

    implicit none

    integer(i64), value :: n0_nodeidc
    integer(i64), value :: n1_nodeidc
    integer(i64), intent(in) :: nodeidc(0:n1_nodeidc - 1_i64,0: &
          n0_nodeidc - 1_i64)
    integer(i64), value :: nbelements
    integer(i64), value :: n0_faces
    integer(i64), value :: n1_faces
    integer(i64), intent(inout) :: faces(0:n1_faces - 1_i64,0:n0_faces - &
          1_i64)
    integer(i64), value :: n0_cellf
    integer(i64), value :: n1_cellf
    integer(i64), intent(inout) :: cellf(0:n1_cellf - 1_i64,0:n0_cellf - &
          1_i64)

    call create_2dfaces(nodeidc, nbelements, faces, cellf)

  end subroutine bind_c_create_2dfaces
  !........................................

  !........................................
  subroutine bind_c_create_cell_faceid(nbelements, n0_oldTonewIndex, &
        oldTonewIndex, n0_cellf, n1_cellf, cellf, n0_faceid, n1_faceid, &
        faceid, dim) bind(c)

    implicit none

    integer(i64), value :: nbelements
    integer(i64), value :: n0_oldTonewIndex
    integer(i64), intent(in) :: oldTonewIndex(0:n0_oldTonewIndex - 1_i64 &
          )
    integer(i64), value :: n0_cellf
    integer(i64), value :: n1_cellf
    integer(i64), intent(in) :: cellf(0:n1_cellf - 1_i64,0:n0_cellf - &
          1_i64)
    integer(i64), value :: n0_faceid
    integer(i64), value :: n1_faceid
    integer(i64), intent(inout) :: faceid(0:n1_faceid - 1_i64,0: &
          n0_faceid - 1_i64)
    integer(i64), value :: dim

    call create_cell_faceid(nbelements, oldTonewIndex, cellf, faceid, &
          dim)

  end subroutine bind_c_create_cell_faceid
  !........................................

  !........................................
  subroutine bind_c_create_3dfaces(n0_nodeidc, n1_nodeidc, nodeidc, &
        nbelements, n0_faces, n1_faces, faces, n0_cellf, n1_cellf, &
        cellf) bind(c)

    implicit none

    integer(i64), value :: n0_nodeidc
    integer(i64), value :: n1_nodeidc
    integer(i64), intent(in) :: nodeidc(0:n1_nodeidc - 1_i64,0: &
          n0_nodeidc - 1_i64)
    integer(i64), value :: nbelements
    integer(i64), value :: n0_faces
    integer(i64), value :: n1_faces
    integer(i64), intent(inout) :: faces(0:n1_faces - 1_i64,0:n0_faces - &
          1_i64)
    integer(i64), value :: n0_cellf
    integer(i64), value :: n1_cellf
    integer(i64), intent(inout) :: cellf(0:n1_cellf - 1_i64,0:n0_cellf - &
          1_i64)

    call create_3dfaces(nodeidc, nbelements, faces, cellf)

  end subroutine bind_c_create_3dfaces
  !........................................

  !........................................
  subroutine bind_c_create_normalfacesofcell(n0_centerc, n1_centerc, &
        centerc, n0_centerf, n1_centerf, centerf, n0_faceid, n1_faceid, &
        faceid, n0_normal, n1_normal, normal, nbelements, n0_nf, n1_nf, &
        n2_nf, nf, dim) bind(c)

    implicit none

    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_centerf
    integer(i64), value :: n1_centerf
    real(f64), intent(in) :: centerf(0:n1_centerf - 1_i64,0:n0_centerf - &
          1_i64)
    integer(i64), value :: n0_faceid
    integer(i64), value :: n1_faceid
    integer(i64), intent(in) :: faceid(0:n1_faceid - 1_i64,0:n0_faceid - &
          1_i64)
    integer(i64), value :: n0_normal
    integer(i64), value :: n1_normal
    real(f64), intent(in) :: normal(0:n1_normal - 1_i64,0:n0_normal - &
          1_i64)
    integer(i64), value :: nbelements
    integer(i64), value :: n0_nf
    integer(i64), value :: n1_nf
    integer(i64), value :: n2_nf
    real(f64), intent(inout) :: nf(0:n2_nf - 1_i64,0:n1_nf - 1_i64,0: &
          n0_nf - 1_i64)
    integer(i64), value :: dim

    call create_NormalFacesOfCell(centerc, centerf, faceid, normal, &
          nbelements, nf, dim)

  end subroutine bind_c_create_normalfacesofcell
  !........................................

  !........................................
  subroutine bind_c_face_gradient_info_2d(n0_cellidf, n1_cellidf, &
        cellidf, n0_nodeidf, n1_nodeidf, nodeidf, n0_centergf, &
        n1_centergf, centergf, n0_namef, namef, n0_normalf, n1_normalf, &
        normalf, n0_centerc, n1_centerc, centerc, n0_centerh, &
        n1_centerh, centerh, n0_halofid, halofid, n0_vertexn, &
        n1_vertexn, vertexn, n0_airDiamond, airDiamond, n0_param1, &
        param1, n0_param2, param2, n0_param3, param3, n0_param4, param4 &
        , n0_f_1, n1_f_1, f_1, n0_f_2, n1_f_2, f_2, n0_f_3, n1_f_3, f_3 &
        , n0_f_4, n1_f_4, f_4, n0_shift, n1_shift, shift, dim) bind(c)

    implicit none

    integer(i64), value :: n0_cellidf
    integer(i64), value :: n1_cellidf
    integer(i64), intent(in) :: cellidf(0:n1_cellidf - 1_i64,0: &
          n0_cellidf - 1_i64)
    integer(i64), value :: n0_nodeidf
    integer(i64), value :: n1_nodeidf
    integer(i64), intent(in) :: nodeidf(0:n1_nodeidf - 1_i64,0: &
          n0_nodeidf - 1_i64)
    integer(i64), value :: n0_centergf
    integer(i64), value :: n1_centergf
    real(f64), intent(in) :: centergf(0:n1_centergf - 1_i64,0: &
          n0_centergf - 1_i64)
    integer(i64), value :: n0_namef
    integer(i64), intent(in) :: namef(0:n0_namef - 1_i64)
    integer(i64), value :: n0_normalf
    integer(i64), value :: n1_normalf
    real(f64), intent(in) :: normalf(0:n1_normalf - 1_i64,0:n0_normalf - &
          1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: n0_halofid
    integer(i64), intent(in) :: halofid(0:n0_halofid - 1_i64)
    integer(i64), value :: n0_vertexn
    integer(i64), value :: n1_vertexn
    real(f64), intent(in) :: vertexn(0:n1_vertexn - 1_i64,0:n0_vertexn - &
          1_i64)
    integer(i64), value :: n0_airDiamond
    real(f64), intent(inout) :: airDiamond(0:n0_airDiamond - 1_i64)
    integer(i64), value :: n0_param1
    real(f64), intent(inout) :: param1(0:n0_param1 - 1_i64)
    integer(i64), value :: n0_param2
    real(f64), intent(inout) :: param2(0:n0_param2 - 1_i64)
    integer(i64), value :: n0_param3
    real(f64), intent(inout) :: param3(0:n0_param3 - 1_i64)
    integer(i64), value :: n0_param4
    real(f64), intent(inout) :: param4(0:n0_param4 - 1_i64)
    integer(i64), value :: n0_f_1
    integer(i64), value :: n1_f_1
    real(f64), intent(inout) :: f_1(0:n1_f_1 - 1_i64,0:n0_f_1 - 1_i64)
    integer(i64), value :: n0_f_2
    integer(i64), value :: n1_f_2
    real(f64), intent(inout) :: f_2(0:n1_f_2 - 1_i64,0:n0_f_2 - 1_i64)
    integer(i64), value :: n0_f_3
    integer(i64), value :: n1_f_3
    real(f64), intent(inout) :: f_3(0:n1_f_3 - 1_i64,0:n0_f_3 - 1_i64)
    integer(i64), value :: n0_f_4
    integer(i64), value :: n1_f_4
    real(f64), intent(inout) :: f_4(0:n1_f_4 - 1_i64,0:n0_f_4 - 1_i64)
    integer(i64), value :: n0_shift
    integer(i64), value :: n1_shift
    real(f64), intent(in) :: shift(0:n1_shift - 1_i64,0:n0_shift - 1_i64 &
          )
    integer(i64), value :: dim

    call face_gradient_info_2d(cellidf, nodeidf, centergf, namef, &
          normalf, centerc, centerh, halofid, vertexn, airDiamond, &
          param1, param2, param3, param4, f_1, f_2, f_3, f_4, shift, &
          dim)

  end subroutine bind_c_face_gradient_info_2d
  !........................................

  !........................................
  subroutine bind_c_variables(n0_centerc, n1_centerc, centerc, &
        n0_cellidn, n1_cellidn, cellidn, n0_haloidn, n1_haloidn, &
        haloidn, n0_periodicn, n1_periodicn, periodicn, n0_vertexn, &
        n1_vertexn, vertexn, n0_namen, namen, n0_centergn, n1_centergn, &
        n2_centergn, centergn, n0_halocentergn, n1_halocentergn, &
        n2_halocentergn, halocentergn, n0_centerh, n1_centerh, centerh, &
        nbproc, n0_R_x, R_x, n0_R_y, R_y, n0_lambda_x, lambda_x, &
        n0_lambda_y, lambda_y, n0_number, number, n0_shift, n1_shift, &
        shift) bind(c)

    implicit none

    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_cellidn
    integer(i64), value :: n1_cellidn
    integer(i64), intent(in) :: cellidn(0:n1_cellidn - 1_i64,0: &
          n0_cellidn - 1_i64)
    integer(i64), value :: n0_haloidn
    integer(i64), value :: n1_haloidn
    integer(i64), intent(in) :: haloidn(0:n1_haloidn - 1_i64,0: &
          n0_haloidn - 1_i64)
    integer(i64), value :: n0_periodicn
    integer(i64), value :: n1_periodicn
    integer(i64), intent(in) :: periodicn(0:n1_periodicn - 1_i64,0: &
          n0_periodicn - 1_i64)
    integer(i64), value :: n0_vertexn
    integer(i64), value :: n1_vertexn
    real(f64), intent(in) :: vertexn(0:n1_vertexn - 1_i64,0:n0_vertexn - &
          1_i64)
    integer(i64), value :: n0_namen
    integer(i64), intent(in) :: namen(0:n0_namen - 1_i64)
    integer(i64), value :: n0_centergn
    integer(i64), value :: n1_centergn
    integer(i64), value :: n2_centergn
    real(f64), intent(in) :: centergn(0:n2_centergn - 1_i64,0: &
          n1_centergn - 1_i64,0:n0_centergn - 1_i64)
    integer(i64), value :: n0_halocentergn
    integer(i64), value :: n1_halocentergn
    integer(i64), value :: n2_halocentergn
    real(f64), intent(in) :: halocentergn(0:n2_halocentergn - 1_i64,0: &
          n1_halocentergn - 1_i64,0:n0_halocentergn - 1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: nbproc
    integer(i64), value :: n0_R_x
    real(f64), intent(inout) :: R_x(0:n0_R_x - 1_i64)
    integer(i64), value :: n0_R_y
    real(f64), intent(inout) :: R_y(0:n0_R_y - 1_i64)
    integer(i64), value :: n0_lambda_x
    real(f64), intent(inout) :: lambda_x(0:n0_lambda_x - 1_i64)
    integer(i64), value :: n0_lambda_y
    real(f64), intent(inout) :: lambda_y(0:n0_lambda_y - 1_i64)
    integer(i64), value :: n0_number
    integer(i64), intent(inout) :: number(0:n0_number - 1_i64)
    integer(i64), value :: n0_shift
    integer(i64), value :: n1_shift
    real(f64), intent(in) :: shift(0:n1_shift - 1_i64,0:n0_shift - 1_i64 &
          )

    call variables(centerc, cellidn, haloidn, periodicn, vertexn, namen, &
          centergn, halocentergn, centerh, nbproc, R_x, R_y, lambda_x, &
          lambda_y, number, shift)

  end subroutine bind_c_variables
  !........................................

  !........................................
  subroutine bind_c_variables_3d(n0_centerc, n1_centerc, centerc, &
        n0_cellidn, n1_cellidn, cellidn, n0_haloidn, n1_haloidn, &
        haloidn, n0_periodicn, n1_periodicn, periodicn, n0_vertexn, &
        n1_vertexn, vertexn, n0_namen, namen, n0_centergn, n1_centergn, &
        n2_centergn, centergn, n0_halocenterg, n1_halocenterg, &
        n2_halocenterg, halocenterg, n0_centerh, n1_centerh, centerh, &
        nbproc, n0_R_x, R_x, n0_R_y, R_y, n0_R_z, R_z, n0_lambda_x, &
        lambda_x, n0_lambda_y, lambda_y, n0_lambda_z, lambda_z, &
        n0_number, number, n0_shift, n1_shift, shift) bind(c)

    implicit none

    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_cellidn
    integer(i64), value :: n1_cellidn
    integer(i64), intent(in) :: cellidn(0:n1_cellidn - 1_i64,0: &
          n0_cellidn - 1_i64)
    integer(i64), value :: n0_haloidn
    integer(i64), value :: n1_haloidn
    integer(i64), intent(in) :: haloidn(0:n1_haloidn - 1_i64,0: &
          n0_haloidn - 1_i64)
    integer(i64), value :: n0_periodicn
    integer(i64), value :: n1_periodicn
    integer(i64), intent(in) :: periodicn(0:n1_periodicn - 1_i64,0: &
          n0_periodicn - 1_i64)
    integer(i64), value :: n0_vertexn
    integer(i64), value :: n1_vertexn
    real(f64), intent(in) :: vertexn(0:n1_vertexn - 1_i64,0:n0_vertexn - &
          1_i64)
    integer(i64), value :: n0_namen
    integer(i64), intent(in) :: namen(0:n0_namen - 1_i64)
    integer(i64), value :: n0_centergn
    integer(i64), value :: n1_centergn
    integer(i64), value :: n2_centergn
    real(f64), intent(in) :: centergn(0:n2_centergn - 1_i64,0: &
          n1_centergn - 1_i64,0:n0_centergn - 1_i64)
    integer(i64), value :: n0_halocenterg
    integer(i64), value :: n1_halocenterg
    integer(i64), value :: n2_halocenterg
    real(f64), intent(in) :: halocenterg(0:n2_halocenterg - 1_i64,0: &
          n1_halocenterg - 1_i64,0:n0_halocenterg - 1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: nbproc
    integer(i64), value :: n0_R_x
    real(f64), intent(inout) :: R_x(0:n0_R_x - 1_i64)
    integer(i64), value :: n0_R_y
    real(f64), intent(inout) :: R_y(0:n0_R_y - 1_i64)
    integer(i64), value :: n0_R_z
    real(f64), intent(inout) :: R_z(0:n0_R_z - 1_i64)
    integer(i64), value :: n0_lambda_x
    real(f64), intent(inout) :: lambda_x(0:n0_lambda_x - 1_i64)
    integer(i64), value :: n0_lambda_y
    real(f64), intent(inout) :: lambda_y(0:n0_lambda_y - 1_i64)
    integer(i64), value :: n0_lambda_z
    real(f64), intent(inout) :: lambda_z(0:n0_lambda_z - 1_i64)
    integer(i64), value :: n0_number
    integer(i64), intent(inout) :: number(0:n0_number - 1_i64)
    integer(i64), value :: n0_shift
    integer(i64), value :: n1_shift
    real(f64), intent(in) :: shift(0:n1_shift - 1_i64,0:n0_shift - 1_i64 &
          )

    call variables_3d(centerc, cellidn, haloidn, periodicn, vertexn, &
          namen, centergn, halocenterg, centerh, nbproc, R_x, R_y, R_z, &
          lambda_x, lambda_y, lambda_z, number, shift)

  end subroutine bind_c_variables_3d
  !........................................

end module bind_c_pyccel_ddm
