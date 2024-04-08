module bind_c_module_ddp

  use module_ddp, only: create_NormalFacesOfCell
  use module_ddp, only: create_2dfaces
  use module_ddp, only: create_info_2dfaces
  use module_ddp, only: create_cell_faceid
  use module_ddp, only: Compute_2dcentervolumeOfCell
  use module_ddp, only: create_cellsOfFace
  use module_ddp, only: create_3dfaces
  use module_ddp, only: create_info_3dfaces
  use module_ddp, only: Compute_3dcentervolumeOfCell

  use, intrinsic :: ISO_C_Binding, only : f64 => C_DOUBLE , i64 => &
      C_INT64_T
  implicit none

  contains

  !........................................
  subroutine bind_c_create_info_2dfaces(n0_cellid, n1_cellid, cellid, &
      n0_nodeid, n1_nodeid, nodeid, n0_namen, namen, n0_vertex, &
      n1_vertex, vertex, n0_centerc, n1_centerc, centerc, nbfaces, &
      n0_normalf, n1_normalf, normalf, n0_mesuref, mesuref, n0_centerf, &
      n1_centerf, centerf, n0_namef, namef) bind(c)

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
      n0_normalf, n1_normalf, normalf, n0_mesuref, mesuref, n0_centerf, &
      n1_centerf, centerf, n0_namef, namef) bind(c)

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
      nbelements, nbfaces, n0_cellid, n1_cellid, cellid, dim) bind(c)

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
      nbelements, n0_faces, n1_faces, faces, n0_cellf, n1_cellf, cellf &
      ) bind(c)

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
      nbelements, n0_faces, n1_faces, faces, n0_cellf, n1_cellf, cellf &
      ) bind(c)

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

end module bind_c_module_ddp
