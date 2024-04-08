module module_ddp


  use, intrinsic :: ISO_C_Binding, only : f64 => C_DOUBLE , i64 => &
      C_INT64_T
  implicit none

  contains

  !........................................
  pure subroutine create_info_2dfaces(cellid, nodeid, namen, vertex, &
      centerc, nbfaces, normalf, mesuref, centerf, namef)

    implicit none

    integer(i64), intent(in) :: cellid(0:,0:)
    integer(i64), intent(in) :: nodeid(0:,0:)
    integer(i64), intent(in) :: namen(0:)
    real(f64), intent(in) :: vertex(0:,0:)
    real(f64), intent(in) :: centerc(0:,0:)
    integer(i64), value :: nbfaces
    real(f64), intent(inout) :: normalf(0:,0:)
    real(f64), intent(inout) :: mesuref(0:)
    real(f64), intent(inout) :: centerf(0:,0:)
    integer(i64), intent(inout) :: namef(0:)
    real(f64) :: norm(0:2_i64)
    real(f64), allocatable :: snorm(:)
    integer(i64) :: i
    integer(i64) :: i_0001
    integer(i64) :: i_0002
    integer(i64) :: i_0003

    norm = 0.0_f64
    allocate(snorm(0:2_i64))
    snorm = 0.0_f64
    !Faces aux bords (1,2,3,4), Faces à l'interieur 0    A VOIR !!!!!
    do i = 0_i64, nbfaces - 1_i64, 1_i64
      if (cellid(1_i64, i) == -1_i64 .and. cellid(1_i64, i) /= -10_i64) &
      then
        if (namen(nodeid(0_i64, i)) == namen(nodeid(1_i64, i))) then
          namef(i) = namen(nodeid(0_i64, i))
        else if ((namen(nodeid(0_i64, i)) == 3_i64 .and. namen(nodeid( &
      1_i64, i)) /= 0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. &
      namen(nodeid(1_i64, i)) == 3_i64)) then
          namef(i) = 3_i64
        else if ((namen(nodeid(0_i64, i)) == 4_i64 .and. namen(nodeid( &
      1_i64, i)) /= 0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. &
      namen(nodeid(1_i64, i)) == 4_i64)) then
          namef(i) = 4_i64
        else if ((namen(nodeid(0_i64, i)) == 7_i64 .and. namen(nodeid( &
      1_i64, i)) /= 0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. &
      namen(nodeid(1_i64, i)) == 7_i64)) then
          namef(i) = 7_i64
        else if ((namen(nodeid(0_i64, i)) == 8_i64 .and. namen(nodeid( &
      1_i64, i)) /= 0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. &
      namen(nodeid(1_i64, i)) == 8_i64)) then
          namef(i) = 8_i64
        else
          namef(i) = 100_i64
        end if
      end if
      norm(0_i64) = vertex(1_i64, nodeid(0_i64, i)) - vertex(1_i64, &
      nodeid(1_i64, i))
      norm(1_i64) = vertex(0_i64, nodeid(1_i64, i)) - vertex(0_i64, &
      nodeid(0_i64, i))
      do i_0001 = 0_i64, size(centerf, 1_i64, i64) - 1_i64, 1_i64
        centerf(i_0001, i) = 0.5_f64 * (vertex(i_0001, nodeid(0_i64, i &
      )) + vertex(i_0001, nodeid(1_i64, i)))
      end do
      do i_0001 = 0_i64, 2_i64, 1_i64
        snorm(i_0001) = centerc(i_0001, cellid(0_i64, i)) - centerf( &
      i_0001, i)
      end do
      if (snorm(0_i64) * norm(0_i64) + snorm(1_i64) * norm(1_i64) > &
      0_i64) then
        do i_0002 = 0_i64, size(normalf, 1_i64, i64) - 1_i64, 1_i64
          normalf(i_0002, i) = (-1_i64) * norm(i_0002)
        end do
      else
        do i_0003 = 0_i64, size(normalf, 1_i64, i64) - 1_i64, 1_i64
          normalf(i_0003, i) = norm(i_0003)
        end do
      end if
      mesuref(i) = sqrt(normalf(0_i64, i) ** 2_i64 + normalf(1_i64, i) &
      ** 2_i64)
    end do

  end subroutine create_info_2dfaces
  !........................................

  !........................................
  pure subroutine create_info_3dfaces(cellid, nodeid, namen, vertex, &
      centerc, nbfaces, normalf, mesuref, centerf, namef)

    implicit none

    integer(i64), intent(in) :: cellid(0:,0:)
    integer(i64), intent(in) :: nodeid(0:,0:)
    integer(i64), intent(in) :: namen(0:)
    real(f64), intent(in) :: vertex(0:,0:)
    real(f64), intent(in) :: centerc(0:,0:)
    integer(i64), value :: nbfaces
    real(f64), intent(inout) :: normalf(0:,0:)
    real(f64), intent(inout) :: mesuref(0:)
    real(f64), intent(inout) :: centerf(0:,0:)
    integer(i64), intent(inout) :: namef(0:)
    real(f64) :: norm(0:2_i64)
    real(f64), allocatable :: snorm(:)
    real(f64) :: u(0:2_i64)
    real(f64) :: v(0:2_i64)
    integer(i64) :: i
    integer(i64) :: i_0004
    integer(i64) :: i_0005
    integer(i64) :: i_0006

    norm = 0.0_f64
    allocate(snorm(0:2_i64))
    snorm = 0.0_f64
    u = 0.0_f64
    v = 0.0_f64
    do i = 0_i64, nbfaces - 1_i64, 1_i64
      if (cellid(1_i64, i) == -1_i64) then
        if (namen(nodeid(0_i64, i)) == namen(nodeid(1_i64, i)) .and. &
      namen(nodeid(0_i64, i)) == namen(nodeid(2_i64, i))) then
          namef(i) = namen(nodeid(0_i64, i))
        else if ((namen(nodeid(0_i64, i)) == 3_i64 .and. namen(nodeid( &
      1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64, i)) /= 0_i64) .or. &
      (namen(nodeid(0_i64, i)) /= 0_i64 .and. namen(nodeid(1_i64, i)) &
      == 3_i64 .and. namen(nodeid(2_i64, i)) /= 0_i64) .or. (namen( &
      nodeid(0_i64, i)) /= 0_i64 .and. namen(nodeid(1_i64, i)) /= 0_i64 &
      .and. namen(nodeid(2_i64, i)) /= 3_i64)) then
          namef(i) = 3_i64
        else if ((namen(nodeid(0_i64, i)) == 4_i64 .and. namen(nodeid( &
      1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64, i)) /= 0_i64) .or. &
      (namen(nodeid(0_i64, i)) /= 0_i64 .and. namen(nodeid(1_i64, i)) &
      == 4_i64 .and. namen(nodeid(2_i64, i)) /= 0_i64) .or. (namen( &
      nodeid(0_i64, i)) /= 0_i64 .and. namen(nodeid(1_i64, i)) /= 0_i64 &
      .and. namen(nodeid(2_i64, i)) /= 4_i64)) then
          namef(i) = 4_i64
        else
          namef(i) = 100_i64
        end if
      end if
      u(:) = vertex(0_i64:2_i64, nodeid(1_i64, i)) - vertex(0_i64:2_i64, &
      nodeid(0_i64, i))
      v(:) = vertex(0_i64:2_i64, nodeid(2_i64, i)) - vertex(0_i64:2_i64, &
      nodeid(0_i64, i))
      norm(0_i64) = 0.5_f64 * (u(1_i64) * v(2_i64) - u(2_i64) * v(1_i64 &
      ))
      norm(1_i64) = 0.5_f64 * (u(2_i64) * v(0_i64) - u(0_i64) * v(2_i64 &
      ))
      norm(2_i64) = 0.5_f64 * (u(0_i64) * v(1_i64) - u(1_i64) * v(0_i64 &
      ))
      do i_0004 = 0_i64, size(centerf, 1_i64, i64) - 1_i64, 1_i64
        centerf(i_0004, i) = 1.0_f64 / 3_i64 * (vertex(i_0004, nodeid( &
      0_i64, i)) + vertex(i_0004, nodeid(1_i64, i)) + vertex(i_0004, &
      nodeid(2_i64, i)))
      end do
      do i_0004 = 0_i64, 2_i64, 1_i64
        snorm(i_0004) = centerc(i_0004, cellid(0_i64, i)) - centerf( &
      i_0004, i)
      end do
      if (snorm(0_i64) * norm(0_i64) + snorm(1_i64) * norm(1_i64) + &
      snorm(2_i64) * norm(2_i64) > 0_i64) then
        do i_0005 = 0_i64, size(normalf, 1_i64, i64) - 1_i64, 1_i64
          normalf(i_0005, i) = (-1_i64) * norm(i_0005)
        end do
      else
        do i_0006 = 0_i64, size(normalf, 1_i64, i64) - 1_i64, 1_i64
          normalf(i_0006, i) = norm(i_0006)
        end do
      end if
      mesuref(i) = sqrt(normalf(0_i64, i) ** 2_i64 + normalf(1_i64, i) &
      ** 2_i64 + normalf(2_i64, i) ** 2_i64)
    end do

  end subroutine create_info_3dfaces
  !........................................

  !........................................
  pure subroutine Compute_2dcentervolumeOfCell(nodeid, vertex, &
      nbelements, center, volume)

    implicit none

    integer(i64), intent(inout) :: nodeid(0:,0:)
    real(f64), intent(in) :: vertex(0:,0:)
    integer(i64), value :: nbelements
    real(f64), intent(inout) :: center(0:,0:)
    real(f64), intent(inout) :: volume(0:)
    integer(i64) :: i
    integer(i64) :: s_1
    integer(i64) :: s_2
    integer(i64) :: s_3
    real(f64) :: x_1
    real(f64) :: y_1
    real(f64) :: x_2
    real(f64) :: y_2
    real(f64) :: x_3
    real(f64) :: y_3
    real(f64) :: var1

    !calcul du barycentre et volume
    do i = 0_i64, nbelements - 1_i64, 1_i64
      s_1 = nodeid(0_i64, i)
      s_2 = nodeid(1_i64, i)
      s_3 = nodeid(2_i64, i)
      x_1 = vertex(0_i64, s_1)
      y_1 = vertex(1_i64, s_1)
      x_2 = vertex(0_i64, s_2)
      y_2 = vertex(1_i64, s_2)
      x_3 = vertex(0_i64, s_3)
      y_3 = vertex(1_i64, s_3)
      center(0_i64, i) = 1.0_f64 / 3_i64 * (x_1 + x_2 + x_3)
      center(1_i64, i) = 1.0_f64 / 3_i64 * (y_1 + y_2 + y_3)
      center(2_i64, i) = 0.0_f64
      volume(i) = 1.0_f64 / 2_i64 * abs((x_1 - x_2) * (y_1 - y_3) - (x_1 &
      - x_3) * (y_1 - y_2))
      var1 = (x_2 - x_1) * (y_3 - y_1) - (y_2 - y_1) * (x_3 - x_1)
      if (var1 < 0_i64) then
        nodeid(0_i64, i) = s_1
        nodeid(1_i64, i) = s_3
        nodeid(2_i64, i) = s_2
      end if
    end do

  end subroutine Compute_2dcentervolumeOfCell
  !........................................

  !........................................
  pure subroutine Compute_3dcentervolumeOfCell(nodeid, vertex, &
      nbelements, center, volume)

    implicit none

    integer(i64), intent(in) :: nodeid(0:,0:)
    real(f64), intent(in) :: vertex(0:,0:)
    integer(i64), value :: nbelements
    real(f64), intent(inout) :: center(0:,0:)
    real(f64), intent(inout) :: volume(0:)
    real(f64), allocatable :: wedge(:)
    real(f64) :: u(0:2_i64)
    real(f64) :: v(0:2_i64)
    real(f64) :: w(0:2_i64)
    integer(i64) :: i
    integer(i64) :: s_1
    integer(i64) :: s_2
    integer(i64) :: s_3
    integer(i64) :: s_4
    real(f64) :: x_1
    real(f64) :: y_1
    real(f64) :: z_1
    real(f64) :: x_2
    real(f64) :: y_2
    real(f64) :: z_2
    real(f64) :: x_3
    real(f64) :: y_3
    real(f64) :: z_3
    real(f64) :: x_4
    real(f64) :: y_4
    real(f64) :: z_4
    integer(i64) :: i_0007

    allocate(wedge(0:2_i64))
    wedge = 0.0_f64
    u = 0.0_f64
    v = 0.0_f64
    w = 0.0_f64
    !calcul du barycentre et volume
    do i = 0_i64, nbelements - 1_i64, 1_i64
      s_1 = nodeid(0_i64, i)
      s_2 = nodeid(1_i64, i)
      s_3 = nodeid(2_i64, i)
      s_4 = nodeid(3_i64, i)
      x_1 = vertex(0_i64, s_1)
      y_1 = vertex(1_i64, s_1)
      z_1 = vertex(2_i64, s_1)
      x_2 = vertex(0_i64, s_2)
      y_2 = vertex(1_i64, s_2)
      z_2 = vertex(2_i64, s_2)
      x_3 = vertex(0_i64, s_3)
      y_3 = vertex(1_i64, s_3)
      z_3 = vertex(2_i64, s_3)
      x_4 = vertex(0_i64, s_4)
      y_4 = vertex(1_i64, s_4)
      z_4 = vertex(2_i64, s_4)
      center(0_i64, i) = 1.0_f64 / 4_i64 * (x_1 + x_2 + x_3 + x_4)
      center(1_i64, i) = 1.0_f64 / 4_i64 * (y_1 + y_2 + y_3 + y_4)
      center(2_i64, i) = 1.0_f64 / 4_i64 * (z_1 + z_2 + z_3 + z_4)
      do i_0007 = 0_i64, 2_i64, 1_i64
        u(i_0007) = vertex(i_0007, s_2) - vertex(i_0007, s_1)
        v(i_0007) = vertex(i_0007, s_3) - vertex(i_0007, s_1)
        w(i_0007) = vertex(i_0007, s_4) - vertex(i_0007, s_1)
      end do
      wedge(0_i64) = v(1_i64) * w(2_i64) - v(2_i64) * w(1_i64)
      wedge(1_i64) = v(2_i64) * w(0_i64) - v(0_i64) * w(2_i64)
      wedge(2_i64) = v(0_i64) * w(1_i64) - v(1_i64) * w(0_i64)
      volume(i) = 1.0_f64 / 6_i64 * abs(u(0_i64) * wedge(0_i64) + u( &
      1_i64) * wedge(1_i64) + u(2_i64) * wedge(2_i64))
    end do

  end subroutine Compute_3dcentervolumeOfCell
  !........................................

  !........................................
  pure subroutine create_cellsOfFace(faceid, nbelements, nbfaces, cellid &
      , dim)

    implicit none

    integer(i64), intent(inout) :: faceid(0:,0:)
    integer(i64), value :: nbelements
    integer(i64), value :: nbfaces
    integer(i64), intent(inout) :: cellid(0:,0:)
    integer(i64), value :: dim
    integer(i64) :: i
    integer(i64) :: j

    do i = 0_i64, nbelements - 1_i64, 1_i64
      do j = 0_i64, dim + 1_i64 - 1_i64, 1_i64
        if (cellid(0_i64, faceid(j, i)) == -1_i64) then
          cellid(0_i64, faceid(j, i)) = i
        end if
        if (cellid(0_i64, faceid(j, i)) /= i) then
          cellid(0_i64, faceid(j, i)) = cellid(0_i64, faceid(j, i))
          cellid(1_i64, faceid(j, i)) = i
        end if
      end do
    end do

  end subroutine create_cellsOfFace
  !........................................

  !........................................
  pure subroutine create_2dfaces(nodeidc, nbelements, faces, cellf) 

    implicit none

    integer(i64), intent(in) :: nodeidc(0:,0:)
    integer(i64), value :: nbelements
    integer(i64), intent(inout) :: faces(0:,0:)
    integer(i64), intent(inout) :: cellf(0:,0:)
    integer(i64) :: k
    integer(i64) :: i

    !Create 2d faces
    k = 0_i64
    do i = 0_i64, nbelements - 1_i64, 1_i64
      faces(0_i64, k) = nodeidc(0_i64, i)
      faces(1_i64, k) = nodeidc(1_i64, i)
      faces(0_i64, k + 1_i64) = nodeidc(1_i64, i)
      faces(1_i64, k + 1_i64) = nodeidc(2_i64, i)
      faces(0_i64, k + 2_i64) = nodeidc(2_i64, i)
      faces(1_i64, k + 2_i64) = nodeidc(0_i64, i)
      cellf(0_i64, i) = k
      cellf(1_i64, i) = k + 1_i64
      cellf(2_i64, i) = k + 2_i64
      k = k + 3_i64
    end do

  end subroutine create_2dfaces
  !........................................

  !........................................
  pure subroutine create_cell_faceid(nbelements, oldTonewIndex, cellf, &
      faceid, dim)

    implicit none

    integer(i64), value :: nbelements
    integer(i64), intent(in) :: oldTonewIndex(0:)
    integer(i64), intent(in) :: cellf(0:,0:)
    integer(i64), intent(inout) :: faceid(0:,0:)
    integer(i64), value :: dim
    integer(i64) :: i
    integer(i64) :: j

    do i = 0_i64, nbelements - 1_i64, 1_i64
      do j = 0_i64, dim + 1_i64 - 1_i64, 1_i64
        faceid(j, i) = oldTonewIndex(cellf(j, i))
      end do
    end do

  end subroutine create_cell_faceid
  !........................................

  !........................................
  pure subroutine create_3dfaces(nodeidc, nbelements, faces, cellf) 

    implicit none

    integer(i64), intent(in) :: nodeidc(0:,0:)
    integer(i64), value :: nbelements
    integer(i64), intent(inout) :: faces(0:,0:)
    integer(i64), intent(inout) :: cellf(0:,0:)
    integer(i64) :: k
    integer(i64) :: i

    !Create 3d faces
    k = 0_i64
    do i = 0_i64, nbelements - 1_i64, 1_i64
      faces(0_i64, k) = nodeidc(0_i64, i)
      faces(1_i64, k) = nodeidc(1_i64, i)
      faces(2_i64, k) = nodeidc(2_i64, i)
      faces(0_i64, k + 1_i64) = nodeidc(2_i64, i)
      faces(1_i64, k + 1_i64) = nodeidc(3_i64, i)
      faces(2_i64, k + 1_i64) = nodeidc(0_i64, i)
      faces(0_i64, k + 2_i64) = nodeidc(0_i64, i)
      faces(1_i64, k + 2_i64) = nodeidc(1_i64, i)
      faces(2_i64, k + 2_i64) = nodeidc(3_i64, i)
      faces(0_i64, k + 3_i64) = nodeidc(3_i64, i)
      faces(1_i64, k + 3_i64) = nodeidc(1_i64, i)
      faces(2_i64, k + 3_i64) = nodeidc(2_i64, i)
      cellf(0_i64, i) = k
      cellf(1_i64, i) = k + 1_i64
      cellf(2_i64, i) = k + 2_i64
      cellf(3_i64, i) = k + 3_i64
      k = k + 4_i64
    end do

  end subroutine create_3dfaces
  !........................................

  !........................................
  pure subroutine create_NormalFacesOfCell(centerc, centerf, faceid, &
      normal, nbelements, nf, dim)

    implicit none

    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: centerf(0:,0:)
    integer(i64), intent(in) :: faceid(0:,0:)
    real(f64), intent(in) :: normal(0:,0:)
    integer(i64), value :: nbelements
    real(f64), intent(inout) :: nf(0:,0:,0:)
    integer(i64), value :: dim
    real(f64) :: ss(0:2_i64)
    real(f64) :: G(0:2_i64)
    real(f64) :: c(0:2_i64)
    integer(i64) :: i
    integer(i64) :: j
    integer(i64) :: f
    integer(i64) :: i_0008
    integer(i64) :: i_0009
    integer(i64) :: i_0010
    integer(i64) :: i_0011

    ss = 0.0_f64
    G = 0.0_f64
    c = 0.0_f64
    !compute the outgoing normal faces for each cell
    do i = 0_i64, nbelements - 1_i64, 1_i64
      do i_0008 = 0_i64, 2_i64, 1_i64
        G(i_0008) = centerc(i_0008, i)
      end do
      do j = 0_i64, dim + 1_i64 - 1_i64, 1_i64
        f = faceid(j, i)
        do i_0009 = 0_i64, 2_i64, 1_i64
          c(i_0009) = centerf(i_0009, f)
        end do
        if ((G(0_i64) - c(0_i64)) * normal(0_i64, f) + (G(1_i64) - c( &
      1_i64)) * normal(1_i64, f) + (G(2_i64) - c(2_i64)) * normal(2_i64 &
      , f) < 0.0_f64) then
          do i_0010 = 0_i64, 2_i64, 1_i64
            ss(i_0010) = normal(i_0010, f)
          end do
        else
          do i_0011 = 0_i64, 2_i64, 1_i64
            ss(i_0011) = (-1.0_f64) * normal(i_0011, f)
          end do
        end if
        do i_0009 = 0_i64, size(nf, 1_i64, i64) - 1_i64, 1_i64
          nf(i_0009, j, i) = ss(i_0009)
        end do
      end do
    end do

  end subroutine create_NormalFacesOfCell
  !........................................

end module module_ddp
