module bind_c_pyccel_functions

  use pyccel_functions, only: compute_2dmatrix_size
  use pyccel_functions, only: facetocell
  use pyccel_functions, only: face_gradient_2d
  use pyccel_functions, only: haloghost_value_neumann
  use pyccel_functions, only: haloghost_value_dirichlet
  use pyccel_functions, only: haloghost_value_nonslip
  use pyccel_functions, only: haloghost_value_slip
  use pyccel_functions, only: get_rhs_loc_3d
  use pyccel_functions, only: compute_3dmatrix_size
  use pyccel_functions, only: search_element
  use pyccel_functions, only: barthlimiter_3d
  use pyccel_functions, only: compute_P_gradient_3d
  use pyccel_functions, only: get_rhs_glob_3d
  use pyccel_functions, only: get_triplet_3d
  use pyccel_functions, only: get_rhs_glob_2d
  use pyccel_functions, only: face_gradient_3d
  use pyccel_functions, only: barthlimiter_2d
  use pyccel_functions, only: centertovertex_2d
  use pyccel_functions, only: cell_gradient_2d
  use pyccel_functions, only: get_triplet_2d
  use pyccel_functions, only: cell_gradient_3d
  use pyccel_functions, only: face_gradient_2d_uv
  use pyccel_functions, only: compute_P_gradient_2d
  use pyccel_functions, only: centertovertex_3d
  use pyccel_functions, only: get_rhs_loc_2d

  use, intrinsic :: ISO_C_Binding, only : i32 => C_INT32_T , i64 => &
        C_INT64_T , f64 => C_DOUBLE
  implicit none

  contains

  !........................................
  subroutine bind_c_convert_solution(n0_x1, x1, n0_x1converted, &
        x1converted, n0_tc, tc, b0Size) bind(c)

    implicit none

    integer(i64), value :: n0_x1
    real(f64), intent(in) :: x1(0:n0_x1 - 1_i64)
    integer(i64), value :: n0_x1converted
    real(f64), intent(inout) :: x1converted(0:n0_x1converted - 1_i64)
    integer(i64), value :: n0_tc
    integer(i64), intent(in) :: tc(0:n0_tc - 1_i64)
    integer(i64), value :: b0Size
    integer(i64) :: i_0001

    !$omp parallel do
    do i_0001 = 0_i64, b0Size - 1_i64, 1_i64
      x1converted(i_0001) = x1(tc(i_0001))
    end do
    !$omp end parallel do

  end subroutine bind_c_convert_solution
  !........................................

  !........................................
  subroutine bind_c_rhs_value_dirichlet_node(n0_Pbordnode, Pbordnode, &
        n0_nodes, nodes, n0_value, value) bind(c)

    implicit none

    integer(i64), value :: n0_Pbordnode
    real(f64), intent(inout) :: Pbordnode(0:n0_Pbordnode - 1_i64)
    integer(i64), value :: n0_nodes
    integer(i64), intent(in) :: nodes(0:n0_nodes - 1_i64)
    integer(i64), value :: n0_value
    real(f64), intent(in) :: value(0:n0_value - 1_i64)
    integer(i64) :: i_0002
    integer(i64) :: Dummy_0001

    !$omp parallel do
    do Dummy_0001 = 0_i64, size(nodes, kind=i64) - 1_i64, 1_i64
      i_0002 = nodes(Dummy_0001)
      Pbordnode(i_0002) = value(i_0002)
    end do
    !$omp end parallel do

  end subroutine bind_c_rhs_value_dirichlet_node
  !........................................

  !........................................
  subroutine bind_c_rhs_value_dirichlet_face(n0_Pbordface, Pbordface, &
        n0_faces, faces, n0_value, value) bind(c)

    implicit none

    integer(i64), value :: n0_Pbordface
    real(f64), intent(inout) :: Pbordface(0:n0_Pbordface - 1_i64)
    integer(i64), value :: n0_faces
    integer(i64), intent(in) :: faces(0:n0_faces - 1_i64)
    integer(i64), value :: n0_value
    real(f64), intent(in) :: value(0:n0_value - 1_i64)
    integer(i64) :: i_0003
    integer(i64) :: Dummy_0002

    !$omp parallel do
    do Dummy_0002 = 0_i64, size(faces, kind=i64) - 1_i64, 1_i64
      i_0003 = faces(Dummy_0002)
      Pbordface(i_0003) = value(i_0003)
    end do
    !$omp end parallel do

  end subroutine bind_c_rhs_value_dirichlet_face
  !........................................

  !........................................
  subroutine bind_c_ghost_value_slip(n0_u_c, u_c, n0_v_c, v_c, &
        n0_w_ghost, w_ghost, n0_cellid, n1_cellid, cellid, n0_faces, &
        faces, n0_normal, n1_normal, normal, n0_mesure, mesure) bind(c &
        )

    implicit none

    integer(i64), value :: n0_u_c
    real(f64), intent(in) :: u_c(0:n0_u_c - 1_i64)
    integer(i64), value :: n0_v_c
    real(f64), intent(in) :: v_c(0:n0_v_c - 1_i64)
    integer(i64), value :: n0_w_ghost
    real(f64), intent(inout) :: w_ghost(0:n0_w_ghost - 1_i64)
    integer(i64), value :: n0_cellid
    integer(i64), value :: n1_cellid
    integer(i64), intent(in) :: cellid(0:n1_cellid - 1_i64,0:n0_cellid - &
          1_i64)
    integer(i64), value :: n0_faces
    integer(i64), intent(in) :: faces(0:n0_faces - 1_i64)
    integer(i64), value :: n0_normal
    integer(i64), value :: n1_normal
    real(f64), intent(in) :: normal(0:n1_normal - 1_i64,0:n0_normal - &
          1_i64)
    integer(i64), value :: n0_mesure
    real(f64), intent(in) :: mesure(0:n0_mesure - 1_i64)
    real(f64), allocatable :: s_n_0001(:)
    integer(i64) :: i_0004
    real(f64) :: u_i_0001
    real(f64) :: v_i_0001
    real(f64) :: u_g_0001
    integer(i64) :: Dummy_0003

    allocate(s_n_0001(0:2_i64))
    s_n_0001 = 0.0_f64
    !$omp parallel do
    do Dummy_0003 = 0_i64, size(faces, kind=i64) - 1_i64, 1_i64
      i_0004 = faces(Dummy_0003)
      u_i_0001 = u_c(cellid(0_i64, i_0004))
      v_i_0001 = v_c(cellid(0_i64, i_0004))
      s_n_0001(:) = normal(:, i_0004) / mesure(i_0004)
      u_g_0001 = u_i_0001 * (s_n_0001(1_i64) * s_n_0001(1_i64) - &
            s_n_0001(0_i64) * s_n_0001(0_i64)) - 2.0_f64 * v_i_0001 * &
            s_n_0001(0_i64) * s_n_0001(1_i64)
      w_ghost(i_0004) = u_c(cellid(0_i64, i_0004)) * u_g_0001
    end do
    !$omp end parallel do
    if (allocated(s_n_0001)) then
      deallocate(s_n_0001)
    end if

  end subroutine bind_c_ghost_value_slip
  !........................................

  !........................................
  subroutine bind_c_ghost_value_nonslip(n0_w_c, w_c, n0_w_ghost, w_ghost &
        , n0_cellid, n1_cellid, cellid, n0_faces, faces) bind(c)

    implicit none

    integer(i64), value :: n0_w_c
    real(f64), intent(in) :: w_c(0:n0_w_c - 1_i64)
    integer(i64), value :: n0_w_ghost
    real(f64), intent(inout) :: w_ghost(0:n0_w_ghost - 1_i64)
    integer(i64), value :: n0_cellid
    integer(i64), value :: n1_cellid
    integer(i64), intent(in) :: cellid(0:n1_cellid - 1_i64,0:n0_cellid - &
          1_i64)
    integer(i64), value :: n0_faces
    integer(i64), intent(in) :: faces(0:n0_faces - 1_i64)
    integer(i64) :: i_0005
    integer(i64) :: Dummy_0004

    !$omp parallel do
    do Dummy_0004 = 0_i64, size(faces, kind=i64) - 1_i64, 1_i64
      i_0005 = faces(Dummy_0004)
      w_ghost(i_0005) = (-1_i64) * w_c(cellid(0_i64, i_0005))
    end do
    !$omp end parallel do

  end subroutine bind_c_ghost_value_nonslip
  !........................................

  !........................................
  subroutine bind_c_ghost_value_neumann(n0_w_c, w_c, n0_w_ghost, w_ghost &
        , n0_cellid, n1_cellid, cellid, n0_faces, faces) bind(c)

    implicit none

    integer(i64), value :: n0_w_c
    real(f64), intent(in) :: w_c(0:n0_w_c - 1_i64)
    integer(i64), value :: n0_w_ghost
    real(f64), intent(inout) :: w_ghost(0:n0_w_ghost - 1_i64)
    integer(i64), value :: n0_cellid
    integer(i64), value :: n1_cellid
    integer(i64), intent(in) :: cellid(0:n1_cellid - 1_i64,0:n0_cellid - &
          1_i64)
    integer(i64), value :: n0_faces
    integer(i64), intent(in) :: faces(0:n0_faces - 1_i64)
    integer(i64) :: i_0006
    integer(i64) :: Dummy_0005

    !$omp parallel do
    do Dummy_0005 = 0_i64, size(faces, kind=i64) - 1_i64, 1_i64
      i_0006 = faces(Dummy_0005)
      w_ghost(i_0006) = w_c(cellid(0_i64, i_0006))
    end do
    !$omp end parallel do

  end subroutine bind_c_ghost_value_neumann
  !........................................

  !........................................
  subroutine bind_c_ghost_value_dirichlet(n0_value, value, n0_w_ghost, &
        w_ghost, n0_cellid, n1_cellid, cellid, n0_faces, faces) bind(c &
        )

    implicit none

    integer(i64), value :: n0_value
    real(f64), intent(in) :: value(0:n0_value - 1_i64)
    integer(i64), value :: n0_w_ghost
    real(f64), intent(inout) :: w_ghost(0:n0_w_ghost - 1_i64)
    integer(i64), value :: n0_cellid
    integer(i64), value :: n1_cellid
    integer(i64), intent(in) :: cellid(0:n1_cellid - 1_i64,0:n0_cellid - &
          1_i64)
    integer(i64), value :: n0_faces
    integer(i64), intent(in) :: faces(0:n0_faces - 1_i64)
    integer(i64) :: i_0007
    integer(i64) :: Dummy_0006

    !$omp parallel do
    do Dummy_0006 = 0_i64, size(faces, kind=i64) - 1_i64, 1_i64
      i_0007 = faces(Dummy_0006)
      w_ghost(i_0007) = value(i_0007)
    end do
    !$omp end parallel do

  end subroutine bind_c_ghost_value_dirichlet
  !........................................

  !........................................
  subroutine bind_c_haloghost_value_neumann(n0_w_halo, w_halo, &
        n0_w_haloghost, w_haloghost, n0_haloghostcenter, &
        n1_haloghostcenter, n2_haloghostcenter, haloghostcenter, &
        BCindex, n0_halonodes, halonodes) bind(c)

    implicit none

    integer(i64), value :: n0_w_halo
    real(f64), intent(in) :: w_halo(0:n0_w_halo - 1_i64)
    integer(i64), value :: n0_w_haloghost
    real(f64), intent(inout) :: w_haloghost(0:n0_w_haloghost - 1_i64)
    integer(i64), value :: n0_haloghostcenter
    integer(i64), value :: n1_haloghostcenter
    integer(i64), value :: n2_haloghostcenter
    real(f64), intent(in) :: haloghostcenter(0:n2_haloghostcenter - &
          1_i64,0:n1_haloghostcenter - 1_i64,0:n0_haloghostcenter - &
          1_i64)
    integer(i64), value :: BCindex
    integer(i64), value :: n0_halonodes
    integer(i64), intent(in) :: halonodes(0:n0_halonodes - 1_i64)

    call haloghost_value_neumann(w_halo, w_haloghost, haloghostcenter, &
          BCindex, halonodes)

  end subroutine bind_c_haloghost_value_neumann
  !........................................

  !........................................
  subroutine bind_c_haloghost_value_dirichlet(n0_value, value, &
        n0_w_haloghost, w_haloghost, n0_haloghostcenter, &
        n1_haloghostcenter, n2_haloghostcenter, haloghostcenter, &
        BCindex, n0_halonodes, halonodes) bind(c)

    implicit none

    integer(i64), value :: n0_value
    real(f64), intent(in) :: value(0:n0_value - 1_i64)
    integer(i64), value :: n0_w_haloghost
    real(f64), intent(inout) :: w_haloghost(0:n0_w_haloghost - 1_i64)
    integer(i64), value :: n0_haloghostcenter
    integer(i64), value :: n1_haloghostcenter
    integer(i64), value :: n2_haloghostcenter
    real(f64), intent(in) :: haloghostcenter(0:n2_haloghostcenter - &
          1_i64,0:n1_haloghostcenter - 1_i64,0:n0_haloghostcenter - &
          1_i64)
    integer(i64), value :: BCindex
    integer(i64), value :: n0_halonodes
    integer(i64), intent(in) :: halonodes(0:n0_halonodes - 1_i64)

    call haloghost_value_dirichlet(value, w_haloghost, haloghostcenter, &
          BCindex, halonodes)

  end subroutine bind_c_haloghost_value_dirichlet
  !........................................

  !........................................
  subroutine bind_c_haloghost_value_nonslip(n0_w_halo, w_halo, &
        n0_w_haloghost, w_haloghost, n0_haloghostcenter, &
        n1_haloghostcenter, n2_haloghostcenter, haloghostcenter, &
        BCindex, n0_halonodes, halonodes) bind(c)

    implicit none

    integer(i64), value :: n0_w_halo
    real(f64), intent(in) :: w_halo(0:n0_w_halo - 1_i64)
    integer(i64), value :: n0_w_haloghost
    real(f64), intent(inout) :: w_haloghost(0:n0_w_haloghost - 1_i64)
    integer(i64), value :: n0_haloghostcenter
    integer(i64), value :: n1_haloghostcenter
    integer(i64), value :: n2_haloghostcenter
    real(f64), intent(in) :: haloghostcenter(0:n2_haloghostcenter - &
          1_i64,0:n1_haloghostcenter - 1_i64,0:n0_haloghostcenter - &
          1_i64)
    integer(i64), value :: BCindex
    integer(i64), value :: n0_halonodes
    integer(i64), intent(in) :: halonodes(0:n0_halonodes - 1_i64)

    call haloghost_value_nonslip(w_halo, w_haloghost, haloghostcenter, &
          BCindex, halonodes)

  end subroutine bind_c_haloghost_value_nonslip
  !........................................

  !........................................
  subroutine bind_c_haloghost_value_slip(n0_u_halo, u_halo, n0_v_halo, &
        v_halo, n0_w_haloghost, w_haloghost, n0_haloghostcenter, &
        n1_haloghostcenter, n2_haloghostcenter, haloghostcenter, &
        BCindex, n0_halonodes, halonodes, n0_haloghostfaceinfo, &
        n1_haloghostfaceinfo, n2_haloghostfaceinfo, haloghostfaceinfo) &
        bind(c)

    implicit none

    integer(i64), value :: n0_u_halo
    real(f64), intent(in) :: u_halo(0:n0_u_halo - 1_i64)
    integer(i64), value :: n0_v_halo
    real(f64), intent(in) :: v_halo(0:n0_v_halo - 1_i64)
    integer(i64), value :: n0_w_haloghost
    real(f64), intent(inout) :: w_haloghost(0:n0_w_haloghost - 1_i64)
    integer(i64), value :: n0_haloghostcenter
    integer(i64), value :: n1_haloghostcenter
    integer(i64), value :: n2_haloghostcenter
    real(f64), intent(in) :: haloghostcenter(0:n2_haloghostcenter - &
          1_i64,0:n1_haloghostcenter - 1_i64,0:n0_haloghostcenter - &
          1_i64)
    integer(i64), value :: BCindex
    integer(i64), value :: n0_halonodes
    integer(i64), intent(in) :: halonodes(0:n0_halonodes - 1_i64)
    integer(i64), value :: n0_haloghostfaceinfo
    integer(i64), value :: n1_haloghostfaceinfo
    integer(i64), value :: n2_haloghostfaceinfo
    real(f64), intent(in) :: haloghostfaceinfo(0:n2_haloghostfaceinfo - &
          1_i64,0:n1_haloghostfaceinfo - 1_i64,0:n0_haloghostfaceinfo - &
          1_i64)

    call haloghost_value_slip(u_halo, v_halo, w_haloghost, &
          haloghostcenter, BCindex, halonodes, haloghostfaceinfo)

  end subroutine bind_c_haloghost_value_slip
  !........................................

  !........................................
  subroutine bind_c_cell_gradient_2d(n0_w_c, w_c, n0_w_ghost, w_ghost, &
        n0_w_halo, w_halo, n0_w_haloghost, w_haloghost, n0_centerc, &
        n1_centerc, centerc, n0_cellnid, n1_cellnid, cellnid, &
        n0_halonid, n1_halonid, halonid, n0_nodecid, n1_nodecid, &
        nodecid, n0_periodicn, n1_periodicn, periodicn, n0_periodic, &
        n1_periodic, periodic, n0_namen, namen, n0_centerg, n1_centerg, &
        n2_centerg, centerg, n0_halocenterg, n1_halocenterg, &
        n2_halocenterg, halocenterg, n0_vertexn, n1_vertexn, vertexn, &
        n0_centerh, n1_centerh, centerh, n0_shift, n1_shift, shift, &
        nbproc, n0_w_x, w_x, n0_w_y, w_y, n0_w_z, w_z) bind(c)

    implicit none

    integer(i64), value :: n0_w_c
    real(f64), intent(in) :: w_c(0:n0_w_c - 1_i64)
    integer(i64), value :: n0_w_ghost
    real(f64), intent(in) :: w_ghost(0:n0_w_ghost - 1_i64)
    integer(i64), value :: n0_w_halo
    real(f64), intent(in) :: w_halo(0:n0_w_halo - 1_i64)
    integer(i64), value :: n0_w_haloghost
    real(f64), intent(in) :: w_haloghost(0:n0_w_haloghost - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_cellnid
    integer(i64), value :: n1_cellnid
    integer(i64), intent(in) :: cellnid(0:n1_cellnid - 1_i64,0: &
          n0_cellnid - 1_i64)
    integer(i64), value :: n0_halonid
    integer(i64), value :: n1_halonid
    integer(i64), intent(in) :: halonid(0:n1_halonid - 1_i64,0: &
          n0_halonid - 1_i64)
    integer(i64), value :: n0_nodecid
    integer(i64), value :: n1_nodecid
    integer(i64), intent(in) :: nodecid(0:n1_nodecid - 1_i64,0: &
          n0_nodecid - 1_i64)
    integer(i64), value :: n0_periodicn
    integer(i64), value :: n1_periodicn
    integer(i64), intent(in) :: periodicn(0:n1_periodicn - 1_i64,0: &
          n0_periodicn - 1_i64)
    integer(i64), value :: n0_periodic
    integer(i64), value :: n1_periodic
    integer(i64), intent(in) :: periodic(0:n1_periodic - 1_i64,0: &
          n0_periodic - 1_i64)
    integer(i64), value :: n0_namen
    integer(i64), intent(in) :: namen(0:n0_namen - 1_i64)
    integer(i64), value :: n0_centerg
    integer(i64), value :: n1_centerg
    integer(i64), value :: n2_centerg
    real(f64), intent(in) :: centerg(0:n2_centerg - 1_i64,0:n1_centerg - &
          1_i64,0:n0_centerg - 1_i64)
    integer(i64), value :: n0_halocenterg
    integer(i64), value :: n1_halocenterg
    integer(i64), value :: n2_halocenterg
    real(f64), intent(in) :: halocenterg(0:n2_halocenterg - 1_i64,0: &
          n1_halocenterg - 1_i64,0:n0_halocenterg - 1_i64)
    integer(i64), value :: n0_vertexn
    integer(i64), value :: n1_vertexn
    real(f64), intent(in) :: vertexn(0:n1_vertexn - 1_i64,0:n0_vertexn - &
          1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: n0_shift
    integer(i64), value :: n1_shift
    real(f64), intent(in) :: shift(0:n1_shift - 1_i64,0:n0_shift - 1_i64 &
          )
    integer(i64), value :: nbproc
    integer(i64), value :: n0_w_x
    real(f64), intent(inout) :: w_x(0:n0_w_x - 1_i64)
    integer(i64), value :: n0_w_y
    real(f64), intent(inout) :: w_y(0:n0_w_y - 1_i64)
    integer(i64), value :: n0_w_z
    real(f64), intent(inout) :: w_z(0:n0_w_z - 1_i64)

    call cell_gradient_2d(w_c, w_ghost, w_halo, w_haloghost, centerc, &
          cellnid, halonid, nodecid, periodicn, periodic, namen, &
          centerg, halocenterg, vertexn, centerh, shift, nbproc, w_x, &
          w_y, w_z)

  end subroutine bind_c_cell_gradient_2d
  !........................................

  !........................................
  subroutine bind_c_cell_gradient_3d(n0_w_c, w_c, n0_w_ghost, w_ghost, &
        n0_w_halo, w_halo, n0_w_haloghost, w_haloghost, n0_centerc, &
        n1_centerc, centerc, n0_cellnid, n1_cellnid, cellnid, &
        n0_halonid, n1_halonid, halonid, n0_nodecid, n1_nodecid, &
        nodecid, n0_periodicn, n1_periodicn, periodicn, n0_periodic, &
        n1_periodic, periodic, n0_namen, namen, n0_centerg, n1_centerg, &
        n2_centerg, centerg, n0_halocenterg, n1_halocenterg, &
        n2_halocenterg, halocenterg, n0_vertexn, n1_vertexn, vertexn, &
        n0_centerh, n1_centerh, centerh, n0_shift, n1_shift, shift, &
        nbproc, n0_w_x, w_x, n0_w_y, w_y, n0_w_z, w_z) bind(c)

    implicit none

    integer(i64), value :: n0_w_c
    real(f64), intent(in) :: w_c(0:n0_w_c - 1_i64)
    integer(i64), value :: n0_w_ghost
    real(f64), intent(in) :: w_ghost(0:n0_w_ghost - 1_i64)
    integer(i64), value :: n0_w_halo
    real(f64), intent(in) :: w_halo(0:n0_w_halo - 1_i64)
    integer(i64), value :: n0_w_haloghost
    real(f64), intent(in) :: w_haloghost(0:n0_w_haloghost - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_cellnid
    integer(i64), value :: n1_cellnid
    integer(i64), intent(in) :: cellnid(0:n1_cellnid - 1_i64,0: &
          n0_cellnid - 1_i64)
    integer(i64), value :: n0_halonid
    integer(i64), value :: n1_halonid
    integer(i64), intent(in) :: halonid(0:n1_halonid - 1_i64,0: &
          n0_halonid - 1_i64)
    integer(i64), value :: n0_nodecid
    integer(i64), value :: n1_nodecid
    integer(i64), intent(in) :: nodecid(0:n1_nodecid - 1_i64,0: &
          n0_nodecid - 1_i64)
    integer(i64), value :: n0_periodicn
    integer(i64), value :: n1_periodicn
    integer(i64), intent(in) :: periodicn(0:n1_periodicn - 1_i64,0: &
          n0_periodicn - 1_i64)
    integer(i64), value :: n0_periodic
    integer(i64), value :: n1_periodic
    integer(i64), intent(in) :: periodic(0:n1_periodic - 1_i64,0: &
          n0_periodic - 1_i64)
    integer(i64), value :: n0_namen
    integer(i64), intent(in) :: namen(0:n0_namen - 1_i64)
    integer(i64), value :: n0_centerg
    integer(i64), value :: n1_centerg
    integer(i64), value :: n2_centerg
    real(f64), intent(in) :: centerg(0:n2_centerg - 1_i64,0:n1_centerg - &
          1_i64,0:n0_centerg - 1_i64)
    integer(i64), value :: n0_halocenterg
    integer(i64), value :: n1_halocenterg
    integer(i64), value :: n2_halocenterg
    real(f64), intent(in) :: halocenterg(0:n2_halocenterg - 1_i64,0: &
          n1_halocenterg - 1_i64,0:n0_halocenterg - 1_i64)
    integer(i64), value :: n0_vertexn
    integer(i64), value :: n1_vertexn
    real(f64), intent(in) :: vertexn(0:n1_vertexn - 1_i64,0:n0_vertexn - &
          1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: n0_shift
    integer(i64), value :: n1_shift
    real(f64), intent(in) :: shift(0:n1_shift - 1_i64,0:n0_shift - 1_i64 &
          )
    integer(i64), value :: nbproc
    integer(i64), value :: n0_w_x
    real(f64), intent(inout) :: w_x(0:n0_w_x - 1_i64)
    integer(i64), value :: n0_w_y
    real(f64), intent(inout) :: w_y(0:n0_w_y - 1_i64)
    integer(i64), value :: n0_w_z
    real(f64), intent(inout) :: w_z(0:n0_w_z - 1_i64)

    call cell_gradient_3d(w_c, w_ghost, w_halo, w_haloghost, centerc, &
          cellnid, halonid, nodecid, periodicn, periodic, namen, &
          centerg, halocenterg, vertexn, centerh, shift, nbproc, w_x, &
          w_y, w_z)

  end subroutine bind_c_cell_gradient_3d
  !........................................

  !........................................
  subroutine bind_c_face_gradient_2d(n0_w_c, w_c, n0_w_ghost, w_ghost, &
        n0_w_halo, w_halo, n0_w_node, w_node, n0_cellidf, n1_cellidf, &
        cellidf, n0_nodeidf, n1_nodeidf, nodeidf, n0_centergf, &
        n1_centergf, centergf, n0_namef, namef, n0_halofid, halofid, &
        n0_centerc, n1_centerc, centerc, n0_centerh, n1_centerh, &
        centerh, n0_vertexn, n1_vertexn, vertexn, n0_airDiamond, &
        airDiamond, n0_normalf, n1_normalf, normalf, n0_f_1, n1_f_1, &
        f_1, n0_f_2, n1_f_2, f_2, n0_f_3, n1_f_3, f_3, n0_f_4, n1_f_4, &
        f_4, n0_shift, n1_shift, shift, n0_wx_face, wx_face, n0_wy_face &
        , wy_face, n0_wz_face, wz_face, n0_innerfaces, innerfaces, &
        n0_halofaces, halofaces, n0_dirichletfaces, dirichletfaces, &
        n0_neumann, neumann, n0_periodicfaces, periodicfaces) bind(c)

    implicit none

    integer(i64), value :: n0_w_c
    real(f64), intent(in) :: w_c(0:n0_w_c - 1_i64)
    integer(i64), value :: n0_w_ghost
    real(f64), intent(in) :: w_ghost(0:n0_w_ghost - 1_i64)
    integer(i64), value :: n0_w_halo
    real(f64), intent(in) :: w_halo(0:n0_w_halo - 1_i64)
    integer(i64), value :: n0_w_node
    real(f64), intent(in) :: w_node(0:n0_w_node - 1_i64)
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
    integer(i64), value :: n0_halofid
    integer(i64), intent(in) :: halofid(0:n0_halofid - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: n0_vertexn
    integer(i64), value :: n1_vertexn
    real(f64), intent(in) :: vertexn(0:n1_vertexn - 1_i64,0:n0_vertexn - &
          1_i64)
    integer(i64), value :: n0_airDiamond
    real(f64), intent(in) :: airDiamond(0:n0_airDiamond - 1_i64)
    integer(i64), value :: n0_normalf
    integer(i64), value :: n1_normalf
    real(f64), intent(in) :: normalf(0:n1_normalf - 1_i64,0:n0_normalf - &
          1_i64)
    integer(i64), value :: n0_f_1
    integer(i64), value :: n1_f_1
    real(f64), intent(in) :: f_1(0:n1_f_1 - 1_i64,0:n0_f_1 - 1_i64)
    integer(i64), value :: n0_f_2
    integer(i64), value :: n1_f_2
    real(f64), intent(in) :: f_2(0:n1_f_2 - 1_i64,0:n0_f_2 - 1_i64)
    integer(i64), value :: n0_f_3
    integer(i64), value :: n1_f_3
    real(f64), intent(in) :: f_3(0:n1_f_3 - 1_i64,0:n0_f_3 - 1_i64)
    integer(i64), value :: n0_f_4
    integer(i64), value :: n1_f_4
    real(f64), intent(in) :: f_4(0:n1_f_4 - 1_i64,0:n0_f_4 - 1_i64)
    integer(i64), value :: n0_shift
    integer(i64), value :: n1_shift
    real(f64), intent(in) :: shift(0:n1_shift - 1_i64,0:n0_shift - 1_i64 &
          )
    integer(i64), value :: n0_wx_face
    real(f64), intent(inout) :: wx_face(0:n0_wx_face - 1_i64)
    integer(i64), value :: n0_wy_face
    real(f64), intent(inout) :: wy_face(0:n0_wy_face - 1_i64)
    integer(i64), value :: n0_wz_face
    real(f64), intent(in) :: wz_face(0:n0_wz_face - 1_i64)
    integer(i64), value :: n0_innerfaces
    integer(i64), intent(in) :: innerfaces(0:n0_innerfaces - 1_i64)
    integer(i64), value :: n0_halofaces
    integer(i64), intent(in) :: halofaces(0:n0_halofaces - 1_i64)
    integer(i64), value :: n0_dirichletfaces
    integer(i64), intent(in) :: dirichletfaces(0:n0_dirichletfaces - &
          1_i64)
    integer(i64), value :: n0_neumann
    integer(i64), intent(in) :: neumann(0:n0_neumann - 1_i64)
    integer(i64), value :: n0_periodicfaces
    integer(i64), intent(in) :: periodicfaces(0:n0_periodicfaces - 1_i64 &
          )

    call face_gradient_2d(w_c, w_ghost, w_halo, w_node, cellidf, nodeidf &
          , centergf, namef, halofid, centerc, centerh, vertexn, &
          airDiamond, normalf, f_1, f_2, f_3, f_4, shift, wx_face, &
          wy_face, wz_face, innerfaces, halofaces, dirichletfaces, &
          neumann, periodicfaces)

  end subroutine bind_c_face_gradient_2d
  !........................................

  !........................................
  subroutine bind_c_face_gradient_2d_uv(n0_h_c, h_c, n0_h_ghost, h_ghost &
        , n0_h_halo, h_halo, n0_h_node, h_node, n0_w_c, w_c, n0_w_ghost &
        , w_ghost, n0_w_halo, w_halo, n0_w_node, w_node, n0_cellidf, &
        n1_cellidf, cellidf, n0_nodeidf, n1_nodeidf, nodeidf, &
        n0_centergf, n1_centergf, centergf, n0_namef, namef, n0_halofid &
        , halofid, n0_centerc, n1_centerc, centerc, n0_centerh, &
        n1_centerh, centerh, n0_vertexn, n1_vertexn, vertexn, &
        n0_airDiamond, airDiamond, n0_normalf, n1_normalf, normalf, &
        n0_f_1, n1_f_1, f_1, n0_f_2, n1_f_2, f_2, n0_f_3, n1_f_3, f_3, &
        n0_f_4, n1_f_4, f_4, n0_shift, n1_shift, shift, n0_wx_face, &
        wx_face, n0_wy_face, wy_face, n0_wz_face, wz_face, &
        n0_innerfaces, innerfaces, n0_halofaces, halofaces, &
        n0_dirichletfaces, dirichletfaces, n0_neumann, neumann, &
        n0_periodicfaces, periodicfaces) bind(c)

    implicit none

    integer(i64), value :: n0_h_c
    real(f64), intent(in) :: h_c(0:n0_h_c - 1_i64)
    integer(i64), value :: n0_h_ghost
    real(f64), intent(in) :: h_ghost(0:n0_h_ghost - 1_i64)
    integer(i64), value :: n0_h_halo
    real(f64), intent(in) :: h_halo(0:n0_h_halo - 1_i64)
    integer(i64), value :: n0_h_node
    real(f64), intent(in) :: h_node(0:n0_h_node - 1_i64)
    integer(i64), value :: n0_w_c
    real(f64), intent(in) :: w_c(0:n0_w_c - 1_i64)
    integer(i64), value :: n0_w_ghost
    real(f64), intent(in) :: w_ghost(0:n0_w_ghost - 1_i64)
    integer(i64), value :: n0_w_halo
    real(f64), intent(in) :: w_halo(0:n0_w_halo - 1_i64)
    integer(i64), value :: n0_w_node
    real(f64), intent(in) :: w_node(0:n0_w_node - 1_i64)
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
    integer(i64), value :: n0_halofid
    integer(i64), intent(in) :: halofid(0:n0_halofid - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: n0_vertexn
    integer(i64), value :: n1_vertexn
    real(f64), intent(in) :: vertexn(0:n1_vertexn - 1_i64,0:n0_vertexn - &
          1_i64)
    integer(i64), value :: n0_airDiamond
    real(f64), intent(in) :: airDiamond(0:n0_airDiamond - 1_i64)
    integer(i64), value :: n0_normalf
    integer(i64), value :: n1_normalf
    real(f64), intent(in) :: normalf(0:n1_normalf - 1_i64,0:n0_normalf - &
          1_i64)
    integer(i64), value :: n0_f_1
    integer(i64), value :: n1_f_1
    real(f64), intent(in) :: f_1(0:n1_f_1 - 1_i64,0:n0_f_1 - 1_i64)
    integer(i64), value :: n0_f_2
    integer(i64), value :: n1_f_2
    real(f64), intent(in) :: f_2(0:n1_f_2 - 1_i64,0:n0_f_2 - 1_i64)
    integer(i64), value :: n0_f_3
    integer(i64), value :: n1_f_3
    real(f64), intent(in) :: f_3(0:n1_f_3 - 1_i64,0:n0_f_3 - 1_i64)
    integer(i64), value :: n0_f_4
    integer(i64), value :: n1_f_4
    real(f64), intent(in) :: f_4(0:n1_f_4 - 1_i64,0:n0_f_4 - 1_i64)
    integer(i64), value :: n0_shift
    integer(i64), value :: n1_shift
    real(f64), intent(in) :: shift(0:n1_shift - 1_i64,0:n0_shift - 1_i64 &
          )
    integer(i64), value :: n0_wx_face
    real(f64), intent(inout) :: wx_face(0:n0_wx_face - 1_i64)
    integer(i64), value :: n0_wy_face
    real(f64), intent(inout) :: wy_face(0:n0_wy_face - 1_i64)
    integer(i64), value :: n0_wz_face
    real(f64), intent(in) :: wz_face(0:n0_wz_face - 1_i64)
    integer(i64), value :: n0_innerfaces
    integer(i64), intent(in) :: innerfaces(0:n0_innerfaces - 1_i64)
    integer(i64), value :: n0_halofaces
    integer(i64), intent(in) :: halofaces(0:n0_halofaces - 1_i64)
    integer(i64), value :: n0_dirichletfaces
    integer(i64), intent(in) :: dirichletfaces(0:n0_dirichletfaces - &
          1_i64)
    integer(i64), value :: n0_neumann
    integer(i64), intent(in) :: neumann(0:n0_neumann - 1_i64)
    integer(i64), value :: n0_periodicfaces
    integer(i64), intent(in) :: periodicfaces(0:n0_periodicfaces - 1_i64 &
          )

    call face_gradient_2d_uv(h_c, h_ghost, h_halo, h_node, w_c, w_ghost, &
          w_halo, w_node, cellidf, nodeidf, centergf, namef, halofid, &
          centerc, centerh, vertexn, airDiamond, normalf, f_1, f_2, f_3 &
          , f_4, shift, wx_face, wy_face, wz_face, innerfaces, &
          halofaces, dirichletfaces, neumann, periodicfaces)

  end subroutine bind_c_face_gradient_2d_uv
  !........................................

  !........................................
  subroutine bind_c_face_gradient_3d(n0_w_c, w_c, n0_w_ghost, w_ghost, &
        n0_w_halo, w_halo, n0_w_node, w_node, n0_cellidf, n1_cellidf, &
        cellidf, n0_nodeidf, n1_nodeidf, nodeidf, n0_centergf, &
        n1_centergf, centergf, n0_namef, namef, n0_halofid, halofid, &
        n0_centerc, n1_centerc, centerc, n0_centerh, n1_centerh, &
        centerh, n0_vertexn, n1_vertexn, vertexn, n0_airDiamond, &
        airDiamond, n0_normalf, n1_normalf, normalf, n0_f_1, n1_f_1, &
        f_1, n0_f_2, n1_f_2, f_2, n0_f_3, n1_f_3, f_3, n0_f_4, n1_f_4, &
        f_4, n0_shift, n1_shift, shift, n0_wx_face, wx_face, n0_wy_face &
        , wy_face, n0_wz_face, wz_face, n0_innerfaces, innerfaces, &
        n0_halofaces, halofaces, n0_dirichletfaces, dirichletfaces, &
        n0_neumann, neumann, n0_periodicfaces, periodicfaces) bind(c)

    implicit none

    integer(i64), value :: n0_w_c
    real(f64), intent(in) :: w_c(0:n0_w_c - 1_i64)
    integer(i64), value :: n0_w_ghost
    real(f64), intent(in) :: w_ghost(0:n0_w_ghost - 1_i64)
    integer(i64), value :: n0_w_halo
    real(f64), intent(in) :: w_halo(0:n0_w_halo - 1_i64)
    integer(i64), value :: n0_w_node
    real(f64), intent(in) :: w_node(0:n0_w_node - 1_i64)
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
    integer(i64), value :: n0_halofid
    integer(i64), intent(in) :: halofid(0:n0_halofid - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: n0_vertexn
    integer(i64), value :: n1_vertexn
    real(f64), intent(in) :: vertexn(0:n1_vertexn - 1_i64,0:n0_vertexn - &
          1_i64)
    integer(i64), value :: n0_airDiamond
    real(f64), intent(in) :: airDiamond(0:n0_airDiamond - 1_i64)
    integer(i64), value :: n0_normalf
    integer(i64), value :: n1_normalf
    real(f64), intent(in) :: normalf(0:n1_normalf - 1_i64,0:n0_normalf - &
          1_i64)
    integer(i64), value :: n0_f_1
    integer(i64), value :: n1_f_1
    real(f64), intent(in) :: f_1(0:n1_f_1 - 1_i64,0:n0_f_1 - 1_i64)
    integer(i64), value :: n0_f_2
    integer(i64), value :: n1_f_2
    real(f64), intent(in) :: f_2(0:n1_f_2 - 1_i64,0:n0_f_2 - 1_i64)
    integer(i64), value :: n0_f_3
    integer(i64), value :: n1_f_3
    real(f64), intent(in) :: f_3(0:n1_f_3 - 1_i64,0:n0_f_3 - 1_i64)
    integer(i64), value :: n0_f_4
    integer(i64), value :: n1_f_4
    real(f64), intent(in) :: f_4(0:n1_f_4 - 1_i64,0:n0_f_4 - 1_i64)
    integer(i64), value :: n0_shift
    integer(i64), value :: n1_shift
    real(f64), intent(in) :: shift(0:n1_shift - 1_i64,0:n0_shift - 1_i64 &
          )
    integer(i64), value :: n0_wx_face
    real(f64), intent(inout) :: wx_face(0:n0_wx_face - 1_i64)
    integer(i64), value :: n0_wy_face
    real(f64), intent(inout) :: wy_face(0:n0_wy_face - 1_i64)
    integer(i64), value :: n0_wz_face
    real(f64), intent(inout) :: wz_face(0:n0_wz_face - 1_i64)
    integer(i64), value :: n0_innerfaces
    integer(i64), intent(in) :: innerfaces(0:n0_innerfaces - 1_i64)
    integer(i64), value :: n0_halofaces
    integer(i64), intent(in) :: halofaces(0:n0_halofaces - 1_i64)
    integer(i64), value :: n0_dirichletfaces
    integer(i64), intent(in) :: dirichletfaces(0:n0_dirichletfaces - &
          1_i64)
    integer(i64), value :: n0_neumann
    integer(i64), intent(in) :: neumann(0:n0_neumann - 1_i64)
    integer(i64), value :: n0_periodicfaces
    integer(i64), intent(in) :: periodicfaces(0:n0_periodicfaces - 1_i64 &
          )

    call face_gradient_3d(w_c, w_ghost, w_halo, w_node, cellidf, nodeidf &
          , centergf, namef, halofid, centerc, centerh, vertexn, &
          airDiamond, normalf, f_1, f_2, f_3, f_4, shift, wx_face, &
          wy_face, wz_face, innerfaces, halofaces, dirichletfaces, &
          neumann, periodicfaces)

  end subroutine bind_c_face_gradient_3d
  !........................................

  !........................................
  subroutine bind_c_centertovertex_2d(n0_w_c, w_c, n0_w_ghost, w_ghost, &
        n0_w_halo, w_halo, n0_w_haloghost, w_haloghost, n0_centerc, &
        n1_centerc, centerc, n0_centerh, n1_centerh, centerh, &
        n0_cellidn, n1_cellidn, cellidn, n0_periodicn, n1_periodicn, &
        periodicn, n0_haloidn, n1_haloidn, haloidn, n0_vertexn, &
        n1_vertexn, vertexn, n0_namen, namen, n0_centergn, n1_centergn, &
        n2_centergn, centergn, n0_halocentergn, n1_halocentergn, &
        n2_halocentergn, halocentergn, n0_R_x, R_x, n0_R_y, R_y, &
        n0_lambda_x, lambda_x, n0_lambda_y, lambda_y, n0_number, number &
        , n0_shift, n1_shift, shift, nbproc, n0_w_n, w_n) bind(c)

    implicit none

    integer(i64), value :: n0_w_c
    real(f64), intent(in) :: w_c(0:n0_w_c - 1_i64)
    integer(i64), value :: n0_w_ghost
    real(f64), intent(in) :: w_ghost(0:n0_w_ghost - 1_i64)
    integer(i64), value :: n0_w_halo
    real(f64), intent(in) :: w_halo(0:n0_w_halo - 1_i64)
    integer(i64), value :: n0_w_haloghost
    real(f64), intent(in) :: w_haloghost(0:n0_w_haloghost - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: n0_cellidn
    integer(i64), value :: n1_cellidn
    integer(i64), intent(in) :: cellidn(0:n1_cellidn - 1_i64,0: &
          n0_cellidn - 1_i64)
    integer(i64), value :: n0_periodicn
    integer(i64), value :: n1_periodicn
    integer(i64), intent(in) :: periodicn(0:n1_periodicn - 1_i64,0: &
          n0_periodicn - 1_i64)
    integer(i64), value :: n0_haloidn
    integer(i64), value :: n1_haloidn
    integer(i64), intent(in) :: haloidn(0:n1_haloidn - 1_i64,0: &
          n0_haloidn - 1_i64)
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
    integer(i64), value :: n0_R_x
    real(f64), intent(in) :: R_x(0:n0_R_x - 1_i64)
    integer(i64), value :: n0_R_y
    real(f64), intent(in) :: R_y(0:n0_R_y - 1_i64)
    integer(i64), value :: n0_lambda_x
    real(f64), intent(in) :: lambda_x(0:n0_lambda_x - 1_i64)
    integer(i64), value :: n0_lambda_y
    real(f64), intent(in) :: lambda_y(0:n0_lambda_y - 1_i64)
    integer(i64), value :: n0_number
    integer(i64), intent(in) :: number(0:n0_number - 1_i64)
    integer(i64), value :: n0_shift
    integer(i64), value :: n1_shift
    real(f64), intent(in) :: shift(0:n1_shift - 1_i64,0:n0_shift - 1_i64 &
          )
    integer(i64), value :: nbproc
    integer(i64), value :: n0_w_n
    real(f64), intent(inout) :: w_n(0:n0_w_n - 1_i64)

    call centertovertex_2d(w_c, w_ghost, w_halo, w_haloghost, centerc, &
          centerh, cellidn, periodicn, haloidn, vertexn, namen, &
          centergn, halocentergn, R_x, R_y, lambda_x, lambda_y, number, &
          shift, nbproc, w_n)

  end subroutine bind_c_centertovertex_2d
  !........................................

  !........................................
  subroutine bind_c_centertovertex_3d(n0_w_c, w_c, n0_w_ghost, w_ghost, &
        n0_w_halo, w_halo, n0_w_haloghost, w_haloghost, n0_centerc, &
        n1_centerc, centerc, n0_centerh, n1_centerh, centerh, &
        n0_cellidn, n1_cellidn, cellidn, n0_periodicn, n1_periodicn, &
        periodicn, n0_haloidn, n1_haloidn, haloidn, n0_vertexn, &
        n1_vertexn, vertexn, n0_namen, namen, n0_centergn, n1_centergn, &
        n2_centergn, centergn, n0_halocentergn, n1_halocentergn, &
        n2_halocentergn, halocentergn, n0_R_x, R_x, n0_R_y, R_y, n0_R_z &
        , R_z, n0_lambda_x, lambda_x, n0_lambda_y, lambda_y, &
        n0_lambda_z, lambda_z, n0_number, number, n0_shift, n1_shift, &
        shift, nbproc, n0_w_n, w_n) bind(c)

    implicit none

    integer(i64), value :: n0_w_c
    real(f64), intent(in) :: w_c(0:n0_w_c - 1_i64)
    integer(i64), value :: n0_w_ghost
    real(f64), intent(in) :: w_ghost(0:n0_w_ghost - 1_i64)
    integer(i64), value :: n0_w_halo
    real(f64), intent(in) :: w_halo(0:n0_w_halo - 1_i64)
    integer(i64), value :: n0_w_haloghost
    real(f64), intent(in) :: w_haloghost(0:n0_w_haloghost - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: n0_cellidn
    integer(i64), value :: n1_cellidn
    integer(i64), intent(in) :: cellidn(0:n1_cellidn - 1_i64,0: &
          n0_cellidn - 1_i64)
    integer(i64), value :: n0_periodicn
    integer(i64), value :: n1_periodicn
    integer(i64), intent(in) :: periodicn(0:n1_periodicn - 1_i64,0: &
          n0_periodicn - 1_i64)
    integer(i64), value :: n0_haloidn
    integer(i64), value :: n1_haloidn
    integer(i64), intent(in) :: haloidn(0:n1_haloidn - 1_i64,0: &
          n0_haloidn - 1_i64)
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
    integer(i64), value :: n0_R_x
    real(f64), intent(in) :: R_x(0:n0_R_x - 1_i64)
    integer(i64), value :: n0_R_y
    real(f64), intent(in) :: R_y(0:n0_R_y - 1_i64)
    integer(i64), value :: n0_R_z
    real(f64), intent(in) :: R_z(0:n0_R_z - 1_i64)
    integer(i64), value :: n0_lambda_x
    real(f64), intent(in) :: lambda_x(0:n0_lambda_x - 1_i64)
    integer(i64), value :: n0_lambda_y
    real(f64), intent(in) :: lambda_y(0:n0_lambda_y - 1_i64)
    integer(i64), value :: n0_lambda_z
    real(f64), intent(in) :: lambda_z(0:n0_lambda_z - 1_i64)
    integer(i64), value :: n0_number
    integer(i64), intent(in) :: number(0:n0_number - 1_i64)
    integer(i64), value :: n0_shift
    integer(i64), value :: n1_shift
    real(f64), intent(in) :: shift(0:n1_shift - 1_i64,0:n0_shift - 1_i64 &
          )
    integer(i64), value :: nbproc
    integer(i64), value :: n0_w_n
    real(f64), intent(inout) :: w_n(0:n0_w_n - 1_i64)

    call centertovertex_3d(w_c, w_ghost, w_halo, w_haloghost, centerc, &
          centerh, cellidn, periodicn, haloidn, vertexn, namen, &
          centergn, halocentergn, R_x, R_y, R_z, lambda_x, lambda_y, &
          lambda_z, number, shift, nbproc, w_n)

  end subroutine bind_c_centertovertex_3d
  !........................................

  !........................................
  subroutine bind_c_barthlimiter_2d(n0_w_c, w_c, n0_w_ghost, w_ghost, &
        n0_w_halo, w_halo, n0_w_x, w_x, n0_w_y, w_y, n0_w_z, w_z, &
        n0_psi, psi, n0_cellid, n1_cellid, cellid, n0_faceid, n1_faceid &
        , faceid, n0_namef, namef, n0_halofid, halofid, n0_centerc, &
        n1_centerc, centerc, n0_centerf, n1_centerf, centerf) bind(c)

    implicit none

    integer(i64), value :: n0_w_c
    real(f64), intent(in) :: w_c(0:n0_w_c - 1_i64)
    integer(i64), value :: n0_w_ghost
    real(f64), intent(in) :: w_ghost(0:n0_w_ghost - 1_i64)
    integer(i64), value :: n0_w_halo
    real(f64), intent(in) :: w_halo(0:n0_w_halo - 1_i64)
    integer(i64), value :: n0_w_x
    real(f64), intent(in) :: w_x(0:n0_w_x - 1_i64)
    integer(i64), value :: n0_w_y
    real(f64), intent(in) :: w_y(0:n0_w_y - 1_i64)
    integer(i64), value :: n0_w_z
    real(f64), intent(in) :: w_z(0:n0_w_z - 1_i64)
    integer(i64), value :: n0_psi
    real(f64), intent(inout) :: psi(0:n0_psi - 1_i64)
    integer(i64), value :: n0_cellid
    integer(i64), value :: n1_cellid
    integer(i64), intent(in) :: cellid(0:n1_cellid - 1_i64,0:n0_cellid - &
          1_i64)
    integer(i64), value :: n0_faceid
    integer(i64), value :: n1_faceid
    integer(i64), intent(in) :: faceid(0:n1_faceid - 1_i64,0:n0_faceid - &
          1_i64)
    integer(i64), value :: n0_namef
    integer(i64), intent(in) :: namef(0:n0_namef - 1_i64)
    integer(i64), value :: n0_halofid
    integer(i64), intent(in) :: halofid(0:n0_halofid - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_centerf
    integer(i64), value :: n1_centerf
    real(f64), intent(in) :: centerf(0:n1_centerf - 1_i64,0:n0_centerf - &
          1_i64)

    call barthlimiter_2d(w_c, w_ghost, w_halo, w_x, w_y, w_z, psi, &
          cellid, faceid, namef, halofid, centerc, centerf)

  end subroutine bind_c_barthlimiter_2d
  !........................................

  !........................................
  subroutine bind_c_barthlimiter_3d(n0_h_c, h_c, n0_h_ghost, h_ghost, &
        n0_h_halo, h_halo, n0_h_x, h_x, n0_h_y, h_y, n0_h_z, h_z, &
        n0_psi, psi, n0_cellid, n1_cellid, cellid, n0_faceid, n1_faceid &
        , faceid, n0_namef, namef, n0_halofid, halofid, n0_centerc, &
        n1_centerc, centerc, n0_centerf, n1_centerf, centerf) bind(c)

    implicit none

    integer(i64), value :: n0_h_c
    real(f64), intent(in) :: h_c(0:n0_h_c - 1_i64)
    integer(i64), value :: n0_h_ghost
    real(f64), intent(in) :: h_ghost(0:n0_h_ghost - 1_i64)
    integer(i64), value :: n0_h_halo
    real(f64), intent(in) :: h_halo(0:n0_h_halo - 1_i64)
    integer(i64), value :: n0_h_x
    real(f64), intent(in) :: h_x(0:n0_h_x - 1_i64)
    integer(i64), value :: n0_h_y
    real(f64), intent(in) :: h_y(0:n0_h_y - 1_i64)
    integer(i64), value :: n0_h_z
    real(f64), intent(in) :: h_z(0:n0_h_z - 1_i64)
    integer(i64), value :: n0_psi
    real(f64), intent(inout) :: psi(0:n0_psi - 1_i64)
    integer(i64), value :: n0_cellid
    integer(i64), value :: n1_cellid
    integer(i64), intent(in) :: cellid(0:n1_cellid - 1_i64,0:n0_cellid - &
          1_i64)
    integer(i64), value :: n0_faceid
    integer(i64), value :: n1_faceid
    integer(i64), intent(in) :: faceid(0:n1_faceid - 1_i64,0:n0_faceid - &
          1_i64)
    integer(i64), value :: n0_namef
    integer(i64), intent(in) :: namef(0:n0_namef - 1_i64)
    integer(i64), value :: n0_halofid
    integer(i64), intent(in) :: halofid(0:n0_halofid - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_centerf
    integer(i64), value :: n1_centerf
    real(f64), intent(in) :: centerf(0:n1_centerf - 1_i64,0:n0_centerf - &
          1_i64)

    call barthlimiter_3d(h_c, h_ghost, h_halo, h_x, h_y, h_z, psi, &
          cellid, faceid, namef, halofid, centerc, centerf)

  end subroutine bind_c_barthlimiter_3d
  !........................................

  !........................................
  function bind_c_search_element(n0_a, a, target_value) bind(c) result( &
        find)

    implicit none

    integer(i64), value :: n0_a
    integer(i64), intent(in) :: a(0:n0_a - 1_i64)
    integer(i64), value :: target_value
    integer(i64) :: find

    find = search_element(a, target_value)

  end function bind_c_search_element
  !........................................

  !........................................
  subroutine bind_c_get_triplet_2d(n0_cellfid, n1_cellfid, cellfid, &
        n0_nodeidf, n1_nodeidf, nodeidf, n0_vertexn, n1_vertexn, &
        vertexn, n0_halofid, halofid, n0_haloext, n1_haloext, haloext, &
        n0_namen, namen, n0_oldnamen, oldnamen, n0_volume, volume, &
        n0_cellnid, n1_cellnid, cellnid, n0_centerc, n1_centerc, &
        centerc, n0_centerh, n1_centerh, centerh, n0_halonid, &
        n1_halonid, halonid, n0_periodicnid, n1_periodicnid, &
        periodicnid, n0_centergn, n1_centergn, n2_centergn, centergn, &
        n0_halocentergn, n1_halocentergn, n2_halocentergn, halocentergn &
        , n0_airDiamond, airDiamond, n0_lambda_x, lambda_x, n0_lambda_y &
        , lambda_y, n0_number, number, n0_R_x, R_x, n0_R_y, R_y, &
        n0_param1, param1, n0_param2, param2, n0_param3, param3, &
        n0_param4, param4, n0_shift, n1_shift, shift, nbelements, &
        n0_loctoglob, loctoglob, n0_BCdirichlet, BCdirichlet, n0_a_loc, &
        a_loc, n0_irn_loc, irn_loc, n0_jcn_loc, jcn_loc, &
        n0_matrixinnerfaces, matrixinnerfaces, n0_halofaces, halofaces, &
        n0_dirichletfaces, dirichletfaces) bind(c)

    implicit none

    integer(i64), value :: n0_cellfid
    integer(i64), value :: n1_cellfid
    integer(i64), intent(in) :: cellfid(0:n1_cellfid - 1_i64,0: &
          n0_cellfid - 1_i64)
    integer(i64), value :: n0_nodeidf
    integer(i64), value :: n1_nodeidf
    integer(i64), intent(in) :: nodeidf(0:n1_nodeidf - 1_i64,0: &
          n0_nodeidf - 1_i64)
    integer(i64), value :: n0_vertexn
    integer(i64), value :: n1_vertexn
    real(f64), intent(in) :: vertexn(0:n1_vertexn - 1_i64,0:n0_vertexn - &
          1_i64)
    integer(i64), value :: n0_halofid
    integer(i64), intent(in) :: halofid(0:n0_halofid - 1_i64)
    integer(i64), value :: n0_haloext
    integer(i64), value :: n1_haloext
    integer(i64), intent(in) :: haloext(0:n1_haloext - 1_i64,0: &
          n0_haloext - 1_i64)
    integer(i64), value :: n0_namen
    integer(i64), intent(in) :: namen(0:n0_namen - 1_i64)
    integer(i64), value :: n0_oldnamen
    integer(i64), intent(in) :: oldnamen(0:n0_oldnamen - 1_i64)
    integer(i64), value :: n0_volume
    real(f64), intent(in) :: volume(0:n0_volume - 1_i64)
    integer(i64), value :: n0_cellnid
    integer(i64), value :: n1_cellnid
    integer(i64), intent(in) :: cellnid(0:n1_cellnid - 1_i64,0: &
          n0_cellnid - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: n0_halonid
    integer(i64), value :: n1_halonid
    integer(i64), intent(in) :: halonid(0:n1_halonid - 1_i64,0: &
          n0_halonid - 1_i64)
    integer(i64), value :: n0_periodicnid
    integer(i64), value :: n1_periodicnid
    integer(i64), intent(in) :: periodicnid(0:n1_periodicnid - 1_i64,0: &
          n0_periodicnid - 1_i64)
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
    integer(i64), value :: n0_airDiamond
    real(f64), intent(in) :: airDiamond(0:n0_airDiamond - 1_i64)
    integer(i64), value :: n0_lambda_x
    real(f64), intent(in) :: lambda_x(0:n0_lambda_x - 1_i64)
    integer(i64), value :: n0_lambda_y
    real(f64), intent(in) :: lambda_y(0:n0_lambda_y - 1_i64)
    integer(i64), value :: n0_number
    integer(i64), intent(in) :: number(0:n0_number - 1_i64)
    integer(i64), value :: n0_R_x
    real(f64), intent(in) :: R_x(0:n0_R_x - 1_i64)
    integer(i64), value :: n0_R_y
    real(f64), intent(in) :: R_y(0:n0_R_y - 1_i64)
    integer(i64), value :: n0_param1
    real(f64), intent(in) :: param1(0:n0_param1 - 1_i64)
    integer(i64), value :: n0_param2
    real(f64), intent(in) :: param2(0:n0_param2 - 1_i64)
    integer(i64), value :: n0_param3
    real(f64), intent(in) :: param3(0:n0_param3 - 1_i64)
    integer(i64), value :: n0_param4
    real(f64), intent(in) :: param4(0:n0_param4 - 1_i64)
    integer(i64), value :: n0_shift
    integer(i64), value :: n1_shift
    real(f64), intent(in) :: shift(0:n1_shift - 1_i64,0:n0_shift - 1_i64 &
          )
    integer(i64), value :: nbelements
    integer(i64), value :: n0_loctoglob
    integer(i64), intent(in) :: loctoglob(0:n0_loctoglob - 1_i64)
    integer(i64), value :: n0_BCdirichlet
    integer(i64), intent(in) :: BCdirichlet(0:n0_BCdirichlet - 1_i64)
    integer(i64), value :: n0_a_loc
    real(f64), intent(inout) :: a_loc(0:n0_a_loc - 1_i64)
    integer(i64), value :: n0_irn_loc
    integer(i32), intent(inout) :: irn_loc(0:n0_irn_loc - 1_i64)
    integer(i64), value :: n0_jcn_loc
    integer(i32), intent(inout) :: jcn_loc(0:n0_jcn_loc - 1_i64)
    integer(i64), value :: n0_matrixinnerfaces
    integer(i64), intent(in) :: matrixinnerfaces(0:n0_matrixinnerfaces - &
          1_i64)
    integer(i64), value :: n0_halofaces
    integer(i64), intent(in) :: halofaces(0:n0_halofaces - 1_i64)
    integer(i64), value :: n0_dirichletfaces
    integer(i64), intent(in) :: dirichletfaces(0:n0_dirichletfaces - &
          1_i64)

    call get_triplet_2d(cellfid, nodeidf, vertexn, halofid, haloext, &
          namen, oldnamen, volume, cellnid, centerc, centerh, halonid, &
          periodicnid, centergn, halocentergn, airDiamond, lambda_x, &
          lambda_y, number, R_x, R_y, param1, param2, param3, param4, &
          shift, nbelements, loctoglob, BCdirichlet, a_loc, irn_loc, &
          jcn_loc, matrixinnerfaces, halofaces, dirichletfaces)

  end subroutine bind_c_get_triplet_2d
  !........................................

  !........................................
  function bind_c_compute_2dmatrix_size(n0_nodeidf, n1_nodeidf, nodeidf, &
        n0_halofid, halofid, n0_cellnid, n1_cellnid, cellnid, &
        n0_halonid, n1_halonid, halonid, n0_periodicnid, n1_periodicnid &
        , periodicnid, n0_centergn, n1_centergn, n2_centergn, centergn, &
        n0_halocentergn, n1_halocentergn, n2_halocentergn, halocentergn &
        , n0_oldnamen, oldnamen, n0_BCdirichlet, BCdirichlet, &
        n0_matrixinnerfaces, matrixinnerfaces, n0_halofaces, halofaces, &
        n0_dirichletfaces, dirichletfaces) bind(c) result(cmpt)

    implicit none

    integer(i64), value :: n0_nodeidf
    integer(i64), value :: n1_nodeidf
    integer(i64), intent(in) :: nodeidf(0:n1_nodeidf - 1_i64,0: &
          n0_nodeidf - 1_i64)
    integer(i64), value :: n0_halofid
    integer(i64), intent(in) :: halofid(0:n0_halofid - 1_i64)
    integer(i64), value :: n0_cellnid
    integer(i64), value :: n1_cellnid
    integer(i64), intent(in) :: cellnid(0:n1_cellnid - 1_i64,0: &
          n0_cellnid - 1_i64)
    integer(i64), value :: n0_halonid
    integer(i64), value :: n1_halonid
    integer(i64), intent(in) :: halonid(0:n1_halonid - 1_i64,0: &
          n0_halonid - 1_i64)
    integer(i64), value :: n0_periodicnid
    integer(i64), value :: n1_periodicnid
    integer(i64), intent(in) :: periodicnid(0:n1_periodicnid - 1_i64,0: &
          n0_periodicnid - 1_i64)
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
    integer(i64), value :: n0_oldnamen
    integer(i64), intent(in) :: oldnamen(0:n0_oldnamen - 1_i64)
    integer(i64), value :: n0_BCdirichlet
    integer(i64), intent(in) :: BCdirichlet(0:n0_BCdirichlet - 1_i64)
    integer(i64), value :: n0_matrixinnerfaces
    integer(i64), intent(in) :: matrixinnerfaces(0:n0_matrixinnerfaces - &
          1_i64)
    integer(i64), value :: n0_halofaces
    integer(i64), intent(in) :: halofaces(0:n0_halofaces - 1_i64)
    integer(i64), value :: n0_dirichletfaces
    integer(i64), intent(in) :: dirichletfaces(0:n0_dirichletfaces - &
          1_i64)
    integer(i64) :: cmpt

    cmpt = compute_2dmatrix_size(nodeidf, halofid, cellnid, halonid, &
          periodicnid, centergn, halocentergn, oldnamen, BCdirichlet, &
          matrixinnerfaces, halofaces, dirichletfaces)

  end function bind_c_compute_2dmatrix_size
  !........................................

  !........................................
  function bind_c_compute_3dmatrix_size(n0_nodeidf, n1_nodeidf, nodeidf, &
        n0_halofid, halofid, n0_cellnid, n1_cellnid, cellnid, &
        n0_halonid, n1_halonid, halonid, n0_periodicnid, n1_periodicnid &
        , periodicnid, n0_centergn, n1_centergn, n2_centergn, centergn, &
        n0_halocentergn, n1_halocentergn, n2_halocentergn, halocentergn &
        , n0_oldnamen, oldnamen, n0_BCdirichlet, BCdirichlet, &
        n0_matrixinnerfaces, matrixinnerfaces, n0_halofaces, halofaces, &
        n0_dirichletfaces, dirichletfaces) bind(c) result(cmpt)

    implicit none

    integer(i64), value :: n0_nodeidf
    integer(i64), value :: n1_nodeidf
    integer(i64), intent(in) :: nodeidf(0:n1_nodeidf - 1_i64,0: &
          n0_nodeidf - 1_i64)
    integer(i64), value :: n0_halofid
    integer(i64), intent(in) :: halofid(0:n0_halofid - 1_i64)
    integer(i64), value :: n0_cellnid
    integer(i64), value :: n1_cellnid
    integer(i64), intent(in) :: cellnid(0:n1_cellnid - 1_i64,0: &
          n0_cellnid - 1_i64)
    integer(i64), value :: n0_halonid
    integer(i64), value :: n1_halonid
    integer(i64), intent(in) :: halonid(0:n1_halonid - 1_i64,0: &
          n0_halonid - 1_i64)
    integer(i64), value :: n0_periodicnid
    integer(i64), value :: n1_periodicnid
    integer(i64), intent(in) :: periodicnid(0:n1_periodicnid - 1_i64,0: &
          n0_periodicnid - 1_i64)
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
    integer(i64), value :: n0_oldnamen
    integer(i64), intent(in) :: oldnamen(0:n0_oldnamen - 1_i64)
    integer(i64), value :: n0_BCdirichlet
    integer(i64), intent(in) :: BCdirichlet(0:n0_BCdirichlet - 1_i64)
    integer(i64), value :: n0_matrixinnerfaces
    integer(i64), intent(in) :: matrixinnerfaces(0:n0_matrixinnerfaces - &
          1_i64)
    integer(i64), value :: n0_halofaces
    integer(i64), intent(in) :: halofaces(0:n0_halofaces - 1_i64)
    integer(i64), value :: n0_dirichletfaces
    integer(i64), intent(in) :: dirichletfaces(0:n0_dirichletfaces - &
          1_i64)
    integer(i64) :: cmpt

    cmpt = compute_3dmatrix_size(nodeidf, halofid, cellnid, halonid, &
          periodicnid, centergn, halocentergn, oldnamen, BCdirichlet, &
          matrixinnerfaces, halofaces, dirichletfaces)

  end function bind_c_compute_3dmatrix_size
  !........................................

  !........................................
  subroutine bind_c_get_triplet_3d(n0_cellfid, n1_cellfid, cellfid, &
        n0_nodeidf, n1_nodeidf, nodeidf, n0_vertexn, n1_vertexn, &
        vertexn, n0_halofid, halofid, n0_haloext, n1_haloext, haloext, &
        n0_namen, namen, n0_oldnamen, oldnamen, n0_volume, volume, &
        n0_centergn, n1_centergn, n2_centergn, centergn, &
        n0_halocentergn, n1_halocentergn, n2_halocentergn, halocentergn &
        , n0_periodicnid, n1_periodicnid, periodicnid, n0_cellnid, &
        n1_cellnid, cellnid, n0_centerc, n1_centerc, centerc, &
        n0_centerh, n1_centerh, centerh, n0_halonid, n1_halonid, &
        halonid, n0_airDiamond, airDiamond, n0_lambda_x, lambda_x, &
        n0_lambda_y, lambda_y, n0_lambda_z, lambda_z, n0_number, number &
        , n0_R_x, R_x, n0_R_y, R_y, n0_R_z, R_z, n0_param1, param1, &
        n0_param2, param2, n0_param3, param3, n0_shift, n1_shift, shift &
        , n0_loctoglob, loctoglob, n0_BCdirichlet, BCdirichlet, &
        n0_a_loc, a_loc, n0_irn_loc, irn_loc, n0_jcn_loc, jcn_loc, &
        n0_matrixinnerfaces, matrixinnerfaces, n0_halofaces, halofaces, &
        n0_dirichletfaces, dirichletfaces) bind(c)

    implicit none

    integer(i64), value :: n0_cellfid
    integer(i64), value :: n1_cellfid
    integer(i64), intent(in) :: cellfid(0:n1_cellfid - 1_i64,0: &
          n0_cellfid - 1_i64)
    integer(i64), value :: n0_nodeidf
    integer(i64), value :: n1_nodeidf
    integer(i64), intent(in) :: nodeidf(0:n1_nodeidf - 1_i64,0: &
          n0_nodeidf - 1_i64)
    integer(i64), value :: n0_vertexn
    integer(i64), value :: n1_vertexn
    real(f64), intent(in) :: vertexn(0:n1_vertexn - 1_i64,0:n0_vertexn - &
          1_i64)
    integer(i64), value :: n0_halofid
    integer(i64), intent(in) :: halofid(0:n0_halofid - 1_i64)
    integer(i64), value :: n0_haloext
    integer(i64), value :: n1_haloext
    integer(i64), intent(in) :: haloext(0:n1_haloext - 1_i64,0: &
          n0_haloext - 1_i64)
    integer(i64), value :: n0_namen
    integer(i64), intent(in) :: namen(0:n0_namen - 1_i64)
    integer(i64), value :: n0_oldnamen
    integer(i64), intent(in) :: oldnamen(0:n0_oldnamen - 1_i64)
    integer(i64), value :: n0_volume
    real(f64), intent(in) :: volume(0:n0_volume - 1_i64)
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
    integer(i64), value :: n0_periodicnid
    integer(i64), value :: n1_periodicnid
    integer(i64), intent(in) :: periodicnid(0:n1_periodicnid - 1_i64,0: &
          n0_periodicnid - 1_i64)
    integer(i64), value :: n0_cellnid
    integer(i64), value :: n1_cellnid
    integer(i64), intent(in) :: cellnid(0:n1_cellnid - 1_i64,0: &
          n0_cellnid - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: n0_halonid
    integer(i64), value :: n1_halonid
    integer(i64), intent(in) :: halonid(0:n1_halonid - 1_i64,0: &
          n0_halonid - 1_i64)
    integer(i64), value :: n0_airDiamond
    real(f64), intent(in) :: airDiamond(0:n0_airDiamond - 1_i64)
    integer(i64), value :: n0_lambda_x
    real(f64), intent(in) :: lambda_x(0:n0_lambda_x - 1_i64)
    integer(i64), value :: n0_lambda_y
    real(f64), intent(in) :: lambda_y(0:n0_lambda_y - 1_i64)
    integer(i64), value :: n0_lambda_z
    real(f64), intent(in) :: lambda_z(0:n0_lambda_z - 1_i64)
    integer(i64), value :: n0_number
    integer(i64), intent(in) :: number(0:n0_number - 1_i64)
    integer(i64), value :: n0_R_x
    real(f64), intent(in) :: R_x(0:n0_R_x - 1_i64)
    integer(i64), value :: n0_R_y
    real(f64), intent(in) :: R_y(0:n0_R_y - 1_i64)
    integer(i64), value :: n0_R_z
    real(f64), intent(in) :: R_z(0:n0_R_z - 1_i64)
    integer(i64), value :: n0_param1
    real(f64), intent(in) :: param1(0:n0_param1 - 1_i64)
    integer(i64), value :: n0_param2
    real(f64), intent(in) :: param2(0:n0_param2 - 1_i64)
    integer(i64), value :: n0_param3
    real(f64), intent(in) :: param3(0:n0_param3 - 1_i64)
    integer(i64), value :: n0_shift
    integer(i64), value :: n1_shift
    real(f64), intent(in) :: shift(0:n1_shift - 1_i64,0:n0_shift - 1_i64 &
          )
    integer(i64), value :: n0_loctoglob
    integer(i64), intent(in) :: loctoglob(0:n0_loctoglob - 1_i64)
    integer(i64), value :: n0_BCdirichlet
    integer(i64), intent(in) :: BCdirichlet(0:n0_BCdirichlet - 1_i64)
    integer(i64), value :: n0_a_loc
    real(f64), intent(inout) :: a_loc(0:n0_a_loc - 1_i64)
    integer(i64), value :: n0_irn_loc
    integer(i32), intent(inout) :: irn_loc(0:n0_irn_loc - 1_i64)
    integer(i64), value :: n0_jcn_loc
    integer(i32), intent(inout) :: jcn_loc(0:n0_jcn_loc - 1_i64)
    integer(i64), value :: n0_matrixinnerfaces
    integer(i64), intent(in) :: matrixinnerfaces(0:n0_matrixinnerfaces - &
          1_i64)
    integer(i64), value :: n0_halofaces
    integer(i64), intent(in) :: halofaces(0:n0_halofaces - 1_i64)
    integer(i64), value :: n0_dirichletfaces
    integer(i64), intent(in) :: dirichletfaces(0:n0_dirichletfaces - &
          1_i64)

    call get_triplet_3d(cellfid, nodeidf, vertexn, halofid, haloext, &
          namen, oldnamen, volume, centergn, halocentergn, periodicnid, &
          cellnid, centerc, centerh, halonid, airDiamond, lambda_x, &
          lambda_y, lambda_z, number, R_x, R_y, R_z, param1, param2, &
          param3, shift, loctoglob, BCdirichlet, a_loc, irn_loc, &
          jcn_loc, matrixinnerfaces, halofaces, dirichletfaces)

  end subroutine bind_c_get_triplet_3d
  !........................................

  !........................................
  subroutine bind_c_get_rhs_loc_2d(n0_cellfid, n1_cellfid, cellfid, &
        n0_nodeidf, n1_nodeidf, nodeidf, n0_oldname, oldname, n0_volume &
        , volume, n0_centergn, n1_centergn, n2_centergn, centergn, &
        n0_param1, param1, n0_param2, param2, n0_param3, param3, &
        n0_param4, param4, n0_Pbordnode, Pbordnode, n0_Pbordface, &
        Pbordface, n0_rhs_loc, rhs_loc, n0_BCdirichlet, BCdirichlet, &
        n0_centergf, n1_centergf, centergf, n0_matrixinnerfaces, &
        matrixinnerfaces, n0_halofaces, halofaces, n0_dirichletfaces, &
        dirichletfaces) bind(c)

    implicit none

    integer(i64), value :: n0_cellfid
    integer(i64), value :: n1_cellfid
    integer(i64), intent(in) :: cellfid(0:n1_cellfid - 1_i64,0: &
          n0_cellfid - 1_i64)
    integer(i64), value :: n0_nodeidf
    integer(i64), value :: n1_nodeidf
    integer(i64), intent(in) :: nodeidf(0:n1_nodeidf - 1_i64,0: &
          n0_nodeidf - 1_i64)
    integer(i64), value :: n0_oldname
    integer(i64), intent(in) :: oldname(0:n0_oldname - 1_i64)
    integer(i64), value :: n0_volume
    real(f64), intent(in) :: volume(0:n0_volume - 1_i64)
    integer(i64), value :: n0_centergn
    integer(i64), value :: n1_centergn
    integer(i64), value :: n2_centergn
    real(f64), intent(in) :: centergn(0:n2_centergn - 1_i64,0: &
          n1_centergn - 1_i64,0:n0_centergn - 1_i64)
    integer(i64), value :: n0_param1
    real(f64), intent(in) :: param1(0:n0_param1 - 1_i64)
    integer(i64), value :: n0_param2
    real(f64), intent(in) :: param2(0:n0_param2 - 1_i64)
    integer(i64), value :: n0_param3
    real(f64), intent(in) :: param3(0:n0_param3 - 1_i64)
    integer(i64), value :: n0_param4
    real(f64), intent(in) :: param4(0:n0_param4 - 1_i64)
    integer(i64), value :: n0_Pbordnode
    real(f64), intent(in) :: Pbordnode(0:n0_Pbordnode - 1_i64)
    integer(i64), value :: n0_Pbordface
    real(f64), intent(in) :: Pbordface(0:n0_Pbordface - 1_i64)
    integer(i64), value :: n0_rhs_loc
    real(f64), intent(inout) :: rhs_loc(0:n0_rhs_loc - 1_i64)
    integer(i64), value :: n0_BCdirichlet
    integer(i64), intent(in) :: BCdirichlet(0:n0_BCdirichlet - 1_i64)
    integer(i64), value :: n0_centergf
    integer(i64), value :: n1_centergf
    real(f64), intent(in) :: centergf(0:n1_centergf - 1_i64,0: &
          n0_centergf - 1_i64)
    integer(i64), value :: n0_matrixinnerfaces
    integer(i64), intent(in) :: matrixinnerfaces(0:n0_matrixinnerfaces - &
          1_i64)
    integer(i64), value :: n0_halofaces
    integer(i64), intent(in) :: halofaces(0:n0_halofaces - 1_i64)
    integer(i64), value :: n0_dirichletfaces
    integer(i64), intent(in) :: dirichletfaces(0:n0_dirichletfaces - &
          1_i64)

    call get_rhs_loc_2d(cellfid, nodeidf, oldname, volume, centergn, &
          param1, param2, param3, param4, Pbordnode, Pbordface, rhs_loc &
          , BCdirichlet, centergf, matrixinnerfaces, halofaces, &
          dirichletfaces)

  end subroutine bind_c_get_rhs_loc_2d
  !........................................

  !........................................
  subroutine bind_c_get_rhs_glob_2d(n0_cellfid, n1_cellfid, cellfid, &
        n0_nodeidf, n1_nodeidf, nodeidf, n0_oldname, oldname, n0_volume &
        , volume, n0_centergn, n1_centergn, n2_centergn, centergn, &
        n0_loctoglob, loctoglob, n0_param1, param1, n0_param2, param2, &
        n0_param3, param3, n0_param4, param4, n0_Pbordnode, Pbordnode, &
        n0_Pbordface, Pbordface, n0_rhs, rhs, n0_BCdirichlet, &
        BCdirichlet, n0_centergf, n1_centergf, centergf, &
        n0_matrixinnerfaces, matrixinnerfaces, n0_halofaces, halofaces, &
        n0_dirichletfaces, dirichletfaces) bind(c)

    implicit none

    integer(i64), value :: n0_cellfid
    integer(i64), value :: n1_cellfid
    integer(i64), intent(in) :: cellfid(0:n1_cellfid - 1_i64,0: &
          n0_cellfid - 1_i64)
    integer(i64), value :: n0_nodeidf
    integer(i64), value :: n1_nodeidf
    integer(i64), intent(in) :: nodeidf(0:n1_nodeidf - 1_i64,0: &
          n0_nodeidf - 1_i64)
    integer(i64), value :: n0_oldname
    integer(i64), intent(in) :: oldname(0:n0_oldname - 1_i64)
    integer(i64), value :: n0_volume
    real(f64), intent(in) :: volume(0:n0_volume - 1_i64)
    integer(i64), value :: n0_centergn
    integer(i64), value :: n1_centergn
    integer(i64), value :: n2_centergn
    real(f64), intent(in) :: centergn(0:n2_centergn - 1_i64,0: &
          n1_centergn - 1_i64,0:n0_centergn - 1_i64)
    integer(i64), value :: n0_loctoglob
    integer(i64), intent(in) :: loctoglob(0:n0_loctoglob - 1_i64)
    integer(i64), value :: n0_param1
    real(f64), intent(in) :: param1(0:n0_param1 - 1_i64)
    integer(i64), value :: n0_param2
    real(f64), intent(in) :: param2(0:n0_param2 - 1_i64)
    integer(i64), value :: n0_param3
    real(f64), intent(in) :: param3(0:n0_param3 - 1_i64)
    integer(i64), value :: n0_param4
    real(f64), intent(in) :: param4(0:n0_param4 - 1_i64)
    integer(i64), value :: n0_Pbordnode
    real(f64), intent(in) :: Pbordnode(0:n0_Pbordnode - 1_i64)
    integer(i64), value :: n0_Pbordface
    real(f64), intent(in) :: Pbordface(0:n0_Pbordface - 1_i64)
    integer(i64), value :: n0_rhs
    real(f64), intent(inout) :: rhs(0:n0_rhs - 1_i64)
    integer(i64), value :: n0_BCdirichlet
    integer(i64), intent(in) :: BCdirichlet(0:n0_BCdirichlet - 1_i64)
    integer(i64), value :: n0_centergf
    integer(i64), value :: n1_centergf
    real(f64), intent(in) :: centergf(0:n1_centergf - 1_i64,0: &
          n0_centergf - 1_i64)
    integer(i64), value :: n0_matrixinnerfaces
    integer(i64), intent(in) :: matrixinnerfaces(0:n0_matrixinnerfaces - &
          1_i64)
    integer(i64), value :: n0_halofaces
    integer(i64), intent(in) :: halofaces(0:n0_halofaces - 1_i64)
    integer(i64), value :: n0_dirichletfaces
    integer(i64), intent(in) :: dirichletfaces(0:n0_dirichletfaces - &
          1_i64)

    call get_rhs_glob_2d(cellfid, nodeidf, oldname, volume, centergn, &
          loctoglob, param1, param2, param3, param4, Pbordnode, &
          Pbordface, rhs, BCdirichlet, centergf, matrixinnerfaces, &
          halofaces, dirichletfaces)

  end subroutine bind_c_get_rhs_glob_2d
  !........................................

  !........................................
  subroutine bind_c_get_rhs_loc_3d(n0_cellfid, n1_cellfid, cellfid, &
        n0_nodeidf, n1_nodeidf, nodeidf, n0_oldname, oldname, n0_volume &
        , volume, n0_centergn, n1_centergn, n2_centergn, centergn, &
        n0_param1, param1, n0_param2, param2, n0_param3, param3, &
        n0_Pbordnode, Pbordnode, n0_Pbordface, Pbordface, n0_rhs_loc, &
        rhs_loc, n0_BCdirichlet, BCdirichlet, n0_matrixinnerfaces, &
        matrixinnerfaces, n0_halofaces, halofaces, n0_dirichletfaces, &
        dirichletfaces) bind(c)

    implicit none

    integer(i64), value :: n0_cellfid
    integer(i64), value :: n1_cellfid
    integer(i64), intent(in) :: cellfid(0:n1_cellfid - 1_i64,0: &
          n0_cellfid - 1_i64)
    integer(i64), value :: n0_nodeidf
    integer(i64), value :: n1_nodeidf
    integer(i64), intent(in) :: nodeidf(0:n1_nodeidf - 1_i64,0: &
          n0_nodeidf - 1_i64)
    integer(i64), value :: n0_oldname
    integer(i64), intent(in) :: oldname(0:n0_oldname - 1_i64)
    integer(i64), value :: n0_volume
    real(f64), intent(in) :: volume(0:n0_volume - 1_i64)
    integer(i64), value :: n0_centergn
    integer(i64), value :: n1_centergn
    integer(i64), value :: n2_centergn
    real(f64), intent(in) :: centergn(0:n2_centergn - 1_i64,0: &
          n1_centergn - 1_i64,0:n0_centergn - 1_i64)
    integer(i64), value :: n0_param1
    real(f64), intent(in) :: param1(0:n0_param1 - 1_i64)
    integer(i64), value :: n0_param2
    real(f64), intent(in) :: param2(0:n0_param2 - 1_i64)
    integer(i64), value :: n0_param3
    real(f64), intent(in) :: param3(0:n0_param3 - 1_i64)
    integer(i64), value :: n0_Pbordnode
    real(f64), intent(in) :: Pbordnode(0:n0_Pbordnode - 1_i64)
    integer(i64), value :: n0_Pbordface
    real(f64), intent(in) :: Pbordface(0:n0_Pbordface - 1_i64)
    integer(i64), value :: n0_rhs_loc
    real(f64), intent(inout) :: rhs_loc(0:n0_rhs_loc - 1_i64)
    integer(i64), value :: n0_BCdirichlet
    integer(i64), intent(in) :: BCdirichlet(0:n0_BCdirichlet - 1_i64)
    integer(i64), value :: n0_matrixinnerfaces
    integer(i64), intent(in) :: matrixinnerfaces(0:n0_matrixinnerfaces - &
          1_i64)
    integer(i64), value :: n0_halofaces
    integer(i64), intent(in) :: halofaces(0:n0_halofaces - 1_i64)
    integer(i64), value :: n0_dirichletfaces
    integer(i64), intent(in) :: dirichletfaces(0:n0_dirichletfaces - &
          1_i64)

    call get_rhs_loc_3d(cellfid, nodeidf, oldname, volume, centergn, &
          param1, param2, param3, Pbordnode, Pbordface, rhs_loc, &
          BCdirichlet, matrixinnerfaces, halofaces, dirichletfaces)

  end subroutine bind_c_get_rhs_loc_3d
  !........................................

  !........................................
  subroutine bind_c_get_rhs_glob_3d(n0_cellfid, n1_cellfid, cellfid, &
        n0_nodeidf, n1_nodeidf, nodeidf, n0_oldname, oldname, n0_volume &
        , volume, n0_centergn, n1_centergn, n2_centergn, centergn, &
        n0_loctoglob, loctoglob, n0_param1, param1, n0_param2, param2, &
        n0_param3, param3, n0_Pbordnode, Pbordnode, n0_Pbordface, &
        Pbordface, n0_rhs, rhs, n0_BCdirichlet, BCdirichlet, &
        n0_matrixinnerfaces, matrixinnerfaces, n0_halofaces, halofaces, &
        n0_dirichletfaces, dirichletfaces) bind(c)

    implicit none

    integer(i64), value :: n0_cellfid
    integer(i64), value :: n1_cellfid
    integer(i64), intent(in) :: cellfid(0:n1_cellfid - 1_i64,0: &
          n0_cellfid - 1_i64)
    integer(i64), value :: n0_nodeidf
    integer(i64), value :: n1_nodeidf
    integer(i64), intent(in) :: nodeidf(0:n1_nodeidf - 1_i64,0: &
          n0_nodeidf - 1_i64)
    integer(i64), value :: n0_oldname
    integer(i64), intent(in) :: oldname(0:n0_oldname - 1_i64)
    integer(i64), value :: n0_volume
    real(f64), intent(in) :: volume(0:n0_volume - 1_i64)
    integer(i64), value :: n0_centergn
    integer(i64), value :: n1_centergn
    integer(i64), value :: n2_centergn
    real(f64), intent(in) :: centergn(0:n2_centergn - 1_i64,0: &
          n1_centergn - 1_i64,0:n0_centergn - 1_i64)
    integer(i64), value :: n0_loctoglob
    integer(i64), intent(in) :: loctoglob(0:n0_loctoglob - 1_i64)
    integer(i64), value :: n0_param1
    real(f64), intent(in) :: param1(0:n0_param1 - 1_i64)
    integer(i64), value :: n0_param2
    real(f64), intent(in) :: param2(0:n0_param2 - 1_i64)
    integer(i64), value :: n0_param3
    real(f64), intent(in) :: param3(0:n0_param3 - 1_i64)
    integer(i64), value :: n0_Pbordnode
    real(f64), intent(in) :: Pbordnode(0:n0_Pbordnode - 1_i64)
    integer(i64), value :: n0_Pbordface
    real(f64), intent(in) :: Pbordface(0:n0_Pbordface - 1_i64)
    integer(i64), value :: n0_rhs
    real(f64), intent(inout) :: rhs(0:n0_rhs - 1_i64)
    integer(i64), value :: n0_BCdirichlet
    integer(i64), intent(in) :: BCdirichlet(0:n0_BCdirichlet - 1_i64)
    integer(i64), value :: n0_matrixinnerfaces
    integer(i64), intent(in) :: matrixinnerfaces(0:n0_matrixinnerfaces - &
          1_i64)
    integer(i64), value :: n0_halofaces
    integer(i64), intent(in) :: halofaces(0:n0_halofaces - 1_i64)
    integer(i64), value :: n0_dirichletfaces
    integer(i64), intent(in) :: dirichletfaces(0:n0_dirichletfaces - &
          1_i64)

    call get_rhs_glob_3d(cellfid, nodeidf, oldname, volume, centergn, &
          loctoglob, param1, param2, param3, Pbordnode, Pbordface, rhs, &
          BCdirichlet, matrixinnerfaces, halofaces, dirichletfaces)

  end subroutine bind_c_get_rhs_glob_3d
  !........................................

  !........................................
  subroutine bind_c_compute_p_gradient_2d(n0_P_c, P_c, n0_P_ghost, &
        P_ghost, n0_P_halo, P_halo, n0_P_node, P_node, n0_cellidf, &
        n1_cellidf, cellidf, n0_nodeidf, n1_nodeidf, nodeidf, &
        n0_centergf, n1_centergf, centergf, n0_namef, namef, n0_halofid &
        , halofid, n0_centerc, n1_centerc, centerc, n0_centerh, &
        n1_centerh, centerh, n0_oldname, oldname, n0_airDiamond, &
        airDiamond, n0_f_1, n1_f_1, f_1, n0_f_2, n1_f_2, f_2, n0_f_3, &
        n1_f_3, f_3, n0_f_4, n1_f_4, f_4, n0_normalf, n1_normalf, &
        normalf, n0_shift, n1_shift, shift, n0_Pbordnode, Pbordnode, &
        n0_Pbordface, Pbordface, n0_Px_face, Px_face, n0_Py_face, &
        Py_face, n0_Pz_face, Pz_face, n0_BCdirichlet, BCdirichlet, &
        n0_innerfaces, innerfaces, n0_halofaces, halofaces, &
        n0_neumannfaces, neumannfaces, n0_dirichletfaces, &
        dirichletfaces, n0_periodicfaces, periodicfaces) bind(c)

    implicit none

    integer(i64), value :: n0_P_c
    real(f64), intent(in) :: P_c(0:n0_P_c - 1_i64)
    integer(i64), value :: n0_P_ghost
    real(f64), intent(in) :: P_ghost(0:n0_P_ghost - 1_i64)
    integer(i64), value :: n0_P_halo
    real(f64), intent(in) :: P_halo(0:n0_P_halo - 1_i64)
    integer(i64), value :: n0_P_node
    real(f64), intent(in) :: P_node(0:n0_P_node - 1_i64)
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
    integer(i64), value :: n0_halofid
    integer(i64), intent(in) :: halofid(0:n0_halofid - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: n0_oldname
    integer(i64), intent(in) :: oldname(0:n0_oldname - 1_i64)
    integer(i64), value :: n0_airDiamond
    real(f64), intent(in) :: airDiamond(0:n0_airDiamond - 1_i64)
    integer(i64), value :: n0_f_1
    integer(i64), value :: n1_f_1
    real(f64), intent(in) :: f_1(0:n1_f_1 - 1_i64,0:n0_f_1 - 1_i64)
    integer(i64), value :: n0_f_2
    integer(i64), value :: n1_f_2
    real(f64), intent(in) :: f_2(0:n1_f_2 - 1_i64,0:n0_f_2 - 1_i64)
    integer(i64), value :: n0_f_3
    integer(i64), value :: n1_f_3
    real(f64), intent(in) :: f_3(0:n1_f_3 - 1_i64,0:n0_f_3 - 1_i64)
    integer(i64), value :: n0_f_4
    integer(i64), value :: n1_f_4
    real(f64), intent(in) :: f_4(0:n1_f_4 - 1_i64,0:n0_f_4 - 1_i64)
    integer(i64), value :: n0_normalf
    integer(i64), value :: n1_normalf
    real(f64), intent(in) :: normalf(0:n1_normalf - 1_i64,0:n0_normalf - &
          1_i64)
    integer(i64), value :: n0_shift
    integer(i64), value :: n1_shift
    real(f64), intent(in) :: shift(0:n1_shift - 1_i64,0:n0_shift - 1_i64 &
          )
    integer(i64), value :: n0_Pbordnode
    real(f64), intent(in) :: Pbordnode(0:n0_Pbordnode - 1_i64)
    integer(i64), value :: n0_Pbordface
    real(f64), intent(in) :: Pbordface(0:n0_Pbordface - 1_i64)
    integer(i64), value :: n0_Px_face
    real(f64), intent(inout) :: Px_face(0:n0_Px_face - 1_i64)
    integer(i64), value :: n0_Py_face
    real(f64), intent(inout) :: Py_face(0:n0_Py_face - 1_i64)
    integer(i64), value :: n0_Pz_face
    real(f64), intent(in) :: Pz_face(0:n0_Pz_face - 1_i64)
    integer(i64), value :: n0_BCdirichlet
    integer(i64), intent(in) :: BCdirichlet(0:n0_BCdirichlet - 1_i64)
    integer(i64), value :: n0_innerfaces
    integer(i64), intent(in) :: innerfaces(0:n0_innerfaces - 1_i64)
    integer(i64), value :: n0_halofaces
    integer(i64), intent(in) :: halofaces(0:n0_halofaces - 1_i64)
    integer(i64), value :: n0_neumannfaces
    integer(i64), intent(in) :: neumannfaces(0:n0_neumannfaces - 1_i64)
    integer(i64), value :: n0_dirichletfaces
    integer(i64), intent(in) :: dirichletfaces(0:n0_dirichletfaces - &
          1_i64)
    integer(i64), value :: n0_periodicfaces
    integer(i64), intent(in) :: periodicfaces(0:n0_periodicfaces - 1_i64 &
          )

    call compute_P_gradient_2d(P_c, P_ghost, P_halo, P_node, cellidf, &
          nodeidf, centergf, namef, halofid, centerc, centerh, oldname, &
          airDiamond, f_1, f_2, f_3, f_4, normalf, shift, Pbordnode, &
          Pbordface, Px_face, Py_face, Pz_face, BCdirichlet, innerfaces &
          , halofaces, neumannfaces, dirichletfaces, periodicfaces)

  end subroutine bind_c_compute_p_gradient_2d
  !........................................

  !........................................
  subroutine bind_c_compute_p_gradient_3d(n0_val_c, val_c, n0_v_ghost, &
        v_ghost, n0_v_halo, v_halo, n0_v_node, v_node, n0_cellidf, &
        n1_cellidf, cellidf, n0_nodeidf, n1_nodeidf, nodeidf, &
        n0_centergf, n1_centergf, centergf, n0_namef, namef, n0_halofid &
        , halofid, n0_centerc, n1_centerc, centerc, n0_centerh, &
        n1_centerh, centerh, n0_oldname, oldname, n0_airDiamond, &
        airDiamond, n0_n1, n1_n1, n1, n0_n2, n1_n2, n2, n0_n3, n1_n3, &
        n3, n0_n4, n1_n4, n4, n0_normalf, n1_normalf, normalf, n0_shift &
        , n1_shift, shift, n0_Pbordnode, Pbordnode, n0_Pbordface, &
        Pbordface, n0_Px_face, Px_face, n0_Py_face, Py_face, n0_Pz_face &
        , Pz_face, n0_BCdirichlet, BCdirichlet, n0_innerfaces, &
        innerfaces, n0_halofaces, halofaces, n0_neumannfaces, &
        neumannfaces, n0_dirichletfaces, dirichletfaces, &
        n0_periodicfaces, periodicfaces) bind(c)

    implicit none

    integer(i64), value :: n0_val_c
    real(f64), intent(in) :: val_c(0:n0_val_c - 1_i64)
    integer(i64), value :: n0_v_ghost
    real(f64), intent(in) :: v_ghost(0:n0_v_ghost - 1_i64)
    integer(i64), value :: n0_v_halo
    real(f64), intent(in) :: v_halo(0:n0_v_halo - 1_i64)
    integer(i64), value :: n0_v_node
    real(f64), intent(in) :: v_node(0:n0_v_node - 1_i64)
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
    integer(i64), value :: n0_halofid
    integer(i64), intent(in) :: halofid(0:n0_halofid - 1_i64)
    integer(i64), value :: n0_centerc
    integer(i64), value :: n1_centerc
    real(f64), intent(in) :: centerc(0:n1_centerc - 1_i64,0:n0_centerc - &
          1_i64)
    integer(i64), value :: n0_centerh
    integer(i64), value :: n1_centerh
    real(f64), intent(in) :: centerh(0:n1_centerh - 1_i64,0:n0_centerh - &
          1_i64)
    integer(i64), value :: n0_oldname
    integer(i64), intent(in) :: oldname(0:n0_oldname - 1_i64)
    integer(i64), value :: n0_airDiamond
    real(f64), intent(in) :: airDiamond(0:n0_airDiamond - 1_i64)
    integer(i64), value :: n0_n1
    integer(i64), value :: n1_n1
    real(f64), intent(in) :: n1(0:n1_n1 - 1_i64,0:n0_n1 - 1_i64)
    integer(i64), value :: n0_n2
    integer(i64), value :: n1_n2
    real(f64), intent(in) :: n2(0:n1_n2 - 1_i64,0:n0_n2 - 1_i64)
    integer(i64), value :: n0_n3
    integer(i64), value :: n1_n3
    real(f64), intent(in) :: n3(0:n1_n3 - 1_i64,0:n0_n3 - 1_i64)
    integer(i64), value :: n0_n4
    integer(i64), value :: n1_n4
    real(f64), intent(in) :: n4(0:n1_n4 - 1_i64,0:n0_n4 - 1_i64)
    integer(i64), value :: n0_normalf
    integer(i64), value :: n1_normalf
    real(f64), intent(in) :: normalf(0:n1_normalf - 1_i64,0:n0_normalf - &
          1_i64)
    integer(i64), value :: n0_shift
    integer(i64), value :: n1_shift
    real(f64), intent(in) :: shift(0:n1_shift - 1_i64,0:n0_shift - 1_i64 &
          )
    integer(i64), value :: n0_Pbordnode
    real(f64), intent(in) :: Pbordnode(0:n0_Pbordnode - 1_i64)
    integer(i64), value :: n0_Pbordface
    real(f64), intent(in) :: Pbordface(0:n0_Pbordface - 1_i64)
    integer(i64), value :: n0_Px_face
    real(f64), intent(inout) :: Px_face(0:n0_Px_face - 1_i64)
    integer(i64), value :: n0_Py_face
    real(f64), intent(inout) :: Py_face(0:n0_Py_face - 1_i64)
    integer(i64), value :: n0_Pz_face
    real(f64), intent(inout) :: Pz_face(0:n0_Pz_face - 1_i64)
    integer(i64), value :: n0_BCdirichlet
    integer(i64), intent(in) :: BCdirichlet(0:n0_BCdirichlet - 1_i64)
    integer(i64), value :: n0_innerfaces
    integer(i64), intent(in) :: innerfaces(0:n0_innerfaces - 1_i64)
    integer(i64), value :: n0_halofaces
    integer(i64), intent(in) :: halofaces(0:n0_halofaces - 1_i64)
    integer(i64), value :: n0_neumannfaces
    integer(i64), intent(in) :: neumannfaces(0:n0_neumannfaces - 1_i64)
    integer(i64), value :: n0_dirichletfaces
    integer(i64), intent(in) :: dirichletfaces(0:n0_dirichletfaces - &
          1_i64)
    integer(i64), value :: n0_periodicfaces
    integer(i64), intent(in) :: periodicfaces(0:n0_periodicfaces - 1_i64 &
          )

    call compute_P_gradient_3d(val_c, v_ghost, v_halo, v_node, cellidf, &
          nodeidf, centergf, namef, halofid, centerc, centerh, oldname, &
          airDiamond, n1, n2, n3, n4, normalf, shift, Pbordnode, &
          Pbordface, Px_face, Py_face, Pz_face, BCdirichlet, innerfaces &
          , halofaces, neumannfaces, dirichletfaces, periodicfaces)

  end subroutine bind_c_compute_p_gradient_3d
  !........................................

  !........................................
  subroutine bind_c_facetocell(n0_u_face, u_face, n0_u_c, u_c, &
        n0_faceidc, n1_faceidc, faceidc, dim) bind(c)

    implicit none

    integer(i64), value :: n0_u_face
    real(f64), intent(in) :: u_face(0:n0_u_face - 1_i64)
    integer(i64), value :: n0_u_c
    real(f64), intent(inout) :: u_c(0:n0_u_c - 1_i64)
    integer(i64), value :: n0_faceidc
    integer(i64), value :: n1_faceidc
    integer(i64), intent(in) :: faceidc(0:n1_faceidc - 1_i64,0: &
          n0_faceidc - 1_i64)
    integer(i64), value :: dim

    call facetocell(u_face, u_c, faceidc, dim)

  end subroutine bind_c_facetocell
  !........................................

end module bind_c_pyccel_functions
