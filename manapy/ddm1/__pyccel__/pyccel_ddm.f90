module pyccel_ddm


  use, intrinsic :: ISO_C_Binding, only : i64 => C_INT64_T , f64 => &
        C_DOUBLE
  implicit none

  contains

  !........................................
  pure subroutine compute_K(cellid, namef, ghostcenterf, centerc, K, dim &
        )

    implicit none

    integer(i64), intent(in) :: cellid(0:,0:)
    integer(i64), intent(in) :: namef(0:)
    real(f64), intent(in) :: ghostcenterf(0:,0:)
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(inout) :: K(0:,0:)
    integer(i64), value :: dim
    integer(i64) :: nbfaces
    integer(i64) :: i
    integer(i64) :: c_left

    nbfaces = size(cellid,2,i64)
    do i = 0_i64, nbfaces - 1_i64, 1_i64
      if (namef(i) <= 4_i64 .and. namef(i) /= 0_i64) then
        c_left = cellid(0_i64, i)
        K(0_i64:dim - 1_i64, i) = 0.5_f64 * (centerc(0_i64:dim - 1_i64, &
              c_left) + ghostcenterf(0_i64:dim - 1_i64, i))
      end if
    end do

  end subroutine compute_K
  !........................................

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

    !from numpy import double, zeros, sqrt
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
              1_i64, i)) /= 0_i64) .or. (namen(nodeid(0_i64, i)) /= &
              0_i64 .and. namen(nodeid(1_i64, i)) == 3_i64)) then
          namef(i) = 3_i64
        else if ((namen(nodeid(0_i64, i)) == 4_i64 .and. namen(nodeid( &
              1_i64, i)) /= 0_i64) .or. (namen(nodeid(0_i64, i)) /= &
              0_i64 .and. namen(nodeid(1_i64, i)) == 4_i64)) then
          namef(i) = 4_i64
        else if ((namen(nodeid(0_i64, i)) == 33_i64 .and. namen(nodeid( &
              1_i64, i)) /= 0_i64) .or. (namen(nodeid(0_i64, i)) /= &
              0_i64 .and. namen(nodeid(1_i64, i)) == 33_i64)) then
          namef(i) = 33_i64
        else if ((namen(nodeid(0_i64, i)) == 44_i64 .and. namen(nodeid( &
              1_i64, i)) /= 0_i64) .or. (namen(nodeid(0_i64, i)) /= &
              0_i64 .and. namen(nodeid(1_i64, i)) == 44_i64)) then
          namef(i) = 44_i64
        else
          namef(i) = 100_i64
        end if
      end if
      norm(0_i64) = vertex(1_i64, nodeid(0_i64, i)) - vertex(1_i64, &
            nodeid(1_i64, i))
      norm(1_i64) = vertex(0_i64, nodeid(1_i64, i)) - vertex(0_i64, &
            nodeid(0_i64, i))
      centerf(:, i) = 0.5_f64 * (vertex(0_i64:2_i64, nodeid(0_i64, i)) + &
            vertex(0_i64:2_i64, nodeid(1_i64, i)))
      snorm(:) = centerc(:, cellid(0_i64, i)) - centerf(:, i)
      if (snorm(0_i64) * norm(0_i64) + snorm(1_i64) * norm(1_i64) > &
            0_i64) then
        normalf(:, i) = (-1_i64) * norm(:)
      else
        normalf(:, i) = norm(:)
      end if
      mesuref(i) = sqrt(normalf(0_i64, i) ** 2_i64 + normalf(1_i64, i) &
            ** 2_i64)
    end do
    if (allocated(snorm)) then
      deallocate(snorm)
    end if

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

    !from numpy import double, zeros, sqrt
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
        else if ((namen(nodeid(0_i64, i)) == 5_i64 .and. namen(nodeid( &
              1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64, i)) /= &
              0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. namen &
              (nodeid(1_i64, i)) == 5_i64 .and. namen(nodeid(2_i64, i &
              )) /= 0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. &
              namen(nodeid(1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64 &
              , i)) == 5_i64)) then
          namef(i) = 5_i64
        else if ((namen(nodeid(0_i64, i)) == 6_i64 .and. namen(nodeid( &
              1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64, i)) /= &
              0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. namen &
              (nodeid(1_i64, i)) == 6_i64 .and. namen(nodeid(2_i64, i &
              )) /= 0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. &
              namen(nodeid(1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64 &
              , i)) == 6_i64)) then
          namef(i) = 6_i64
        else if ((namen(nodeid(0_i64, i)) == 3_i64 .and. namen(nodeid( &
              1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64, i)) /= &
              0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. namen &
              (nodeid(1_i64, i)) == 3_i64 .and. namen(nodeid(2_i64, i &
              )) /= 0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. &
              namen(nodeid(1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64 &
              , i)) == 3_i64)) then
          namef(i) = 3_i64
        else if ((namen(nodeid(0_i64, i)) == 4_i64 .and. namen(nodeid( &
              1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64, i)) /= &
              0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. namen &
              (nodeid(1_i64, i)) == 4_i64 .and. namen(nodeid(2_i64, i &
              )) /= 0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. &
              namen(nodeid(1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64 &
              , i)) == 4_i64)) then
          namef(i) = 4_i64
        else if ((namen(nodeid(0_i64, i)) == 55_i64 .and. namen(nodeid( &
              1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64, i)) /= &
              0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. namen &
              (nodeid(1_i64, i)) == 55_i64 .and. namen(nodeid(2_i64, i &
              )) /= 0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. &
              namen(nodeid(1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64 &
              , i)) == 55_i64)) then
          namef(i) = 55_i64
        else if ((namen(nodeid(0_i64, i)) == 66_i64 .and. namen(nodeid( &
              1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64, i)) /= &
              0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. namen &
              (nodeid(1_i64, i)) == 66_i64 .and. namen(nodeid(2_i64, i &
              )) /= 0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. &
              namen(nodeid(1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64 &
              , i)) == 66_i64)) then
          namef(i) = 66_i64
        else if ((namen(nodeid(0_i64, i)) == 33_i64 .and. namen(nodeid( &
              1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64, i)) /= &
              0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. namen &
              (nodeid(1_i64, i)) == 33_i64 .and. namen(nodeid(2_i64, i &
              )) /= 0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. &
              namen(nodeid(1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64 &
              , i)) == 33_i64)) then
          namef(i) = 33_i64
        else if ((namen(nodeid(0_i64, i)) == 44_i64 .and. namen(nodeid( &
              1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64, i)) /= &
              0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. namen &
              (nodeid(1_i64, i)) == 44_i64 .and. namen(nodeid(2_i64, i &
              )) /= 0_i64) .or. (namen(nodeid(0_i64, i)) /= 0_i64 .and. &
              namen(nodeid(1_i64, i)) /= 0_i64 .and. namen(nodeid(2_i64 &
              , i)) == 44_i64)) then
          namef(i) = 44_i64
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
      centerf(:, i) = 1.0_f64 / 3_i64 * (vertex(:2_i64, nodeid(0_i64, i &
            )) + vertex(:2_i64, nodeid(1_i64, i)) + vertex(:2_i64, &
            nodeid(2_i64, i)))
      snorm(:) = centerc(:, cellid(0_i64, i)) - centerf(:, i)
      if (snorm(0_i64) * norm(0_i64) + snorm(1_i64) * norm(1_i64) + &
            snorm(2_i64) * norm(2_i64) > 0_i64) then
        normalf(:, i) = (-1_i64) * norm(:)
      else
        normalf(:, i) = norm(:)
      end if
      mesuref(i) = sqrt(normalf(0_i64, i) ** 2_i64 + normalf(1_i64, i) &
            ** 2_i64 + normalf(2_i64, i) ** 2_i64)
    end do
    if (allocated(snorm)) then
      deallocate(snorm)
    end if

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

    !from numpy import zeros, fabs
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
      u(:) = vertex(0_i64:2_i64, s_2) - vertex(0_i64:2_i64, s_1)
      v(:) = vertex(0_i64:2_i64, s_3) - vertex(0_i64:2_i64, s_1)
      w(:) = vertex(0_i64:2_i64, s_4) - vertex(0_i64:2_i64, s_1)
      wedge(0_i64) = v(1_i64) * w(2_i64) - v(2_i64) * w(1_i64)
      wedge(1_i64) = v(2_i64) * w(0_i64) - v(0_i64) * w(2_i64)
      wedge(2_i64) = v(0_i64) * w(1_i64) - v(1_i64) * w(0_i64)
      volume(i) = 1.0_f64 / 6_i64 * abs(u(0_i64) * wedge(0_i64) + u( &
            1_i64) * wedge(1_i64) + u(2_i64) * wedge(2_i64))
    end do
    if (allocated(wedge)) then
      deallocate(wedge)
    end if

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

    !from numpy import zeros
    ss = 0.0_f64
    G = 0.0_f64
    c = 0.0_f64
    !compute the outgoing normal faces for each cell
    do i = 0_i64, nbelements - 1_i64, 1_i64
      G(:) = centerc(:, i)
      do j = 0_i64, dim + 1_i64 - 1_i64, 1_i64
        f = faceid(j, i)
        c(:) = centerf(:, f)
        if ((G(0_i64) - c(0_i64)) * normal(0_i64, f) + (G(1_i64) - c( &
              1_i64)) * normal(1_i64, f) + (G(2_i64) - c(2_i64)) * &
              normal(2_i64, f) < 0.0_f64) then
          ss(:) = normal(:, f)
        else
          ss(:) = (-1.0_f64) * normal(:, f)
        end if
        nf(:, j, i) = ss(:)
      end do
    end do

  end subroutine create_NormalFacesOfCell
  !........................................

  !........................................
  subroutine face_gradient_info_2d(cellidf, nodeidf, centergf, namef, &
        normalf, centerc, centerh, halofid, vertexn, airDiamond, param1 &
        , param2, param3, param4, f_1, f_2, f_3, f_4, shift, dim)

    implicit none

    integer(i64), intent(in) :: cellidf(0:,0:)
    integer(i64), intent(in) :: nodeidf(0:,0:)
    real(f64), intent(in) :: centergf(0:,0:)
    integer(i64), intent(in) :: namef(0:)
    real(f64), intent(in) :: normalf(0:,0:)
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    integer(i64), intent(in) :: halofid(0:)
    real(f64), intent(in) :: vertexn(0:,0:)
    real(f64), intent(inout) :: airDiamond(0:)
    real(f64), intent(inout) :: param1(0:)
    real(f64), intent(inout) :: param2(0:)
    real(f64), intent(inout) :: param3(0:)
    real(f64), intent(inout) :: param4(0:)
    real(f64), intent(inout) :: f_1(0:,0:)
    real(f64), intent(inout) :: f_2(0:,0:)
    real(f64), intent(inout) :: f_3(0:,0:)
    real(f64), intent(inout) :: f_4(0:,0:)
    real(f64), intent(in) :: shift(0:,0:)
    integer(i64), value :: dim
    integer(i64) :: nbface
    real(f64) :: xy_1(0:dim - 1_i64)
    real(f64) :: xy_2(0:dim - 1_i64)
    real(f64) :: v_1(0:dim - 1_i64)
    real(f64) :: v_2(0:dim - 1_i64)
    integer(i64) :: i
    integer(i64) :: c_left
    integer(i64) :: c_right
    integer(i64) :: i_1
    integer(i64) :: i_2
    real(f64) :: n1
    real(f64) :: n2

    nbface = size(cellidf,2,i64)
    xy_1 = 0.0_f64
    xy_2 = 0.0_f64
    v_1 = 0.0_f64
    v_2 = 0.0_f64
    do i = 0_i64, nbface - 1_i64, 1_i64
      c_left = cellidf(0_i64, i)
      c_right = cellidf(1_i64, i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      xy_1(:) = vertexn(0_i64:dim - 1_i64, i_1)
      xy_2(:) = vertexn(0_i64:dim - 1_i64, i_2)
      v_1(:) = centerc(0_i64:dim - 1_i64, c_left)
      if (namef(i) == 0_i64) then
        v_2(:) = centerc(0_i64:dim - 1_i64, c_right)
      else if (namef(i) == 11_i64 .or. namef(i) == 22_i64) then
        v_2(0_i64) = centerc(0_i64, c_right) + shift(0_i64, c_right)
        v_2(1_i64) = centerc(1_i64, c_right)
      else if (namef(i) == 33_i64 .or. namef(i) == 44_i64) then
        v_2(0_i64) = centerc(0_i64, c_right)
        v_2(1_i64) = centerc(1_i64, c_right) + shift(1_i64, c_right)
      else if (namef(i) == 10_i64) then
        v_2(:) = centerh(0_i64:dim - 1_i64, halofid(i))
      else
        v_2(:) = centergf(0_i64:dim - 1_i64, i)
      end if
      f_1(:, i) = v_1(:) - xy_1(:)
      f_2(:, i) = xy_2(:) - v_1(:)
      f_3(:, i) = v_2(:) - xy_2(:)
      f_4(:, i) = xy_1(:) - v_2(:)
      n1 = normalf(0_i64, i)
      n2 = normalf(1_i64, i)
      airDiamond(i) = 0.5_f64 * ((xy_2(0_i64) - xy_1(0_i64)) * (v_2( &
            1_i64) - v_1(1_i64)) + (v_1(0_i64) - v_2(0_i64)) * (xy_2( &
            1_i64) - xy_1(1_i64)))
      param1(i) = 1.0_f64 / (2.0_f64 * airDiamond(i)) * ((f_1(1_i64, i) &
            + f_2(1_i64, i)) * n1 - (f_1(0_i64, i) + f_2(0_i64, i)) * &
            n2)
      param2(i) = 1.0_f64 / (2.0_f64 * airDiamond(i)) * ((f_2(1_i64, i) &
            + f_3(1_i64, i)) * n1 - (f_2(0_i64, i) + f_3(0_i64, i)) * &
            n2)
      param3(i) = 1.0_f64 / (2.0_f64 * airDiamond(i)) * ((f_3(1_i64, i) &
            + f_4(1_i64, i)) * n1 - (f_3(0_i64, i) + f_4(0_i64, i)) * &
            n2)
      param4(i) = 1.0_f64 / (2.0_f64 * airDiamond(i)) * ((f_4(1_i64, i) &
            + f_1(1_i64, i)) * n1 - (f_4(0_i64, i) + f_1(0_i64, i)) * &
            n2)
    end do

  end subroutine face_gradient_info_2d
  !........................................

  !........................................
  subroutine variables(centerc, cellidn, haloidn, periodicn, vertexn, &
        namen, centergn, halocentergn, centerh, nbproc, R_x, R_y, &
        lambda_x, lambda_y, number, shift)

    implicit none

    real(f64), intent(in) :: centerc(0:,0:)
    integer(i64), intent(in) :: cellidn(0:,0:)
    integer(i64), intent(in) :: haloidn(0:,0:)
    integer(i64), intent(in) :: periodicn(0:,0:)
    real(f64), intent(in) :: vertexn(0:,0:)
    integer(i64), intent(in) :: namen(0:)
    real(f64), intent(in) :: centergn(0:,0:,0:)
    real(f64), intent(in) :: halocentergn(0:,0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    integer(i64), value :: nbproc
    real(f64), intent(inout) :: R_x(0:)
    real(f64), intent(inout) :: R_y(0:)
    real(f64), intent(inout) :: lambda_x(0:)
    real(f64), intent(inout) :: lambda_y(0:)
    integer(i64), intent(inout) :: number(0:)
    real(f64), intent(in) :: shift(0:,0:)
    integer(i64) :: nbnode
    real(f64), allocatable :: I_xx(:)
    real(f64), allocatable :: I_yy(:)
    real(f64), allocatable :: I_xy(:)
    real(f64) :: center(0:2_i64)
    integer(i64) :: i
    real(f64) :: D
    integer(i64) :: j
    real(f64) :: Rx
    real(f64) :: Ry
    integer(i64) :: cell

    nbnode = size(R_x, kind=i64)
    allocate(I_xx(0:nbnode - 1_i64))
    I_xx = 0.0_f64
    allocate(I_yy(0:nbnode - 1_i64))
    I_yy = 0.0_f64
    allocate(I_xy(0:nbnode - 1_i64))
    I_xy = 0.0_f64
    center = 0.0_f64
    do i = 0_i64, nbnode - 1_i64, 1_i64
      do j = 0_i64, cellidn(size(cellidn, 1_i64, i64) - 1_i64, i) - &
            1_i64, 1_i64
        center(:) = centerc(0_i64:2_i64, cellidn(j, i))
        Rx = center(0_i64) - vertexn(0_i64, i)
        Ry = center(1_i64) - vertexn(1_i64, i)
        I_xx(i) = I_xx(i) + Rx * Rx
        I_yy(i) = I_yy(i) + Ry * Ry
        I_xy(i) = I_xy(i) + Rx * Ry
        R_x(i) = R_x(i) + Rx
        R_y(i) = R_y(i) + Ry
        number(i) = number(i) + 1_i64
      end do
      !ghost boundary (old vertex names)
      !if vertexn[i][3] == 1 or vertexn[i][3] == 2 or vertexn[i][3] == 3 or vertexn[i][3] == 4:
      if (centergn(2_i64, 0_i64, i) /= -1_i64) then
        do j = 0_i64, size(centergn(:, :, i),2,i64) - 1_i64, 1_i64
          cell = Int(centergn(size(centergn, 1_i64, i64) - 1_i64, j, i), &
                i64)
          if (cell /= -1_i64) then
            center(:) = centergn(0_i64:2_i64, j, i)
            Rx = center(0_i64) - vertexn(0_i64, i)
            Ry = center(1_i64) - vertexn(1_i64, i)
            I_xx(i) = I_xx(i) + Rx * Rx
            I_yy(i) = I_yy(i) + Ry * Ry
            I_xy(i) = I_xy(i) + Rx * Ry
            R_x(i) = R_x(i) + Rx
            R_y(i) = R_y(i) + Ry
            number(i) = number(i) + 1_i64
          end if
        end do
      end if
      !periodic boundary old vertex names)
      if (vertexn(3_i64, i) == 11_i64 .or. vertexn(3_i64, i) == 22_i64) &
            then
        do j = 0_i64, periodicn(size(periodicn, 1_i64, i64) - 1_i64, i) &
              - 1_i64, 1_i64
          cell = periodicn(j, i)
          center(0_i64) = centerc(0_i64, cell) + shift(0_i64, cell)
          center(1_i64) = centerc(1_i64, cell)
          Rx = center(0_i64) - vertexn(0_i64, i)
          Ry = center(1_i64) - vertexn(1_i64, i)
          I_xx(i) = I_xx(i) + Rx * Rx
          I_yy(i) = I_yy(i) + Ry * Ry
          I_xy(i) = I_xy(i) + Rx * Ry
          R_x(i) = R_x(i) + Rx
          R_y(i) = R_y(i) + Ry
          number(i) = number(i) + 1_i64
        end do
      else if (vertexn(3_i64, i) == 33_i64 .or. vertexn(3_i64, i) == &
            44_i64) then
        do j = 0_i64, periodicn(size(periodicn, 1_i64, i64) - 1_i64, i) &
              - 1_i64, 1_i64
          cell = periodicn(j, i)
          center(0_i64) = centerc(0_i64, cell)
          center(1_i64) = centerc(1_i64, cell) + shift(1_i64, cell)
          Rx = center(0_i64) - vertexn(0_i64, i)
          Ry = center(1_i64) - vertexn(1_i64, i)
          I_xx(i) = I_xx(i) + Rx * Rx
          I_yy(i) = I_yy(i) + Ry * Ry
          I_xy(i) = I_xy(i) + Rx * Ry
          R_x(i) = R_x(i) + Rx
          R_y(i) = R_y(i) + Ry
          number(i) = number(i) + 1_i64
        end do
      end if
      if (namen(i) == 10_i64) then
        do j = 0_i64, size(halocentergn(:, :, i),2,i64) - 1_i64, 1_i64
          cell = Int(halocentergn(size(halocentergn, 1_i64, i64) - 1_i64 &
                , j, i), i64)
          if (cell /= -1_i64) then
            center(:) = halocentergn(0_i64:2_i64, j, i)
            Rx = center(0_i64) - vertexn(0_i64, i)
            Ry = center(1_i64) - vertexn(1_i64, i)
            I_xx(i) = I_xx(i) + Rx * Rx
            I_yy(i) = I_yy(i) + Ry * Ry
            I_xy(i) = I_xy(i) + Rx * Ry
            R_x(i) = R_x(i) + Rx
            R_y(i) = R_y(i) + Ry
            number(i) = number(i) + 1_i64
          end if
        end do
        !if haloidn[i][-1] > 0:
        do j = 0_i64, haloidn(size(haloidn, 1_i64, i64) - 1_i64, i) - &
              1_i64, 1_i64
          cell = haloidn(j, i)
          center(:) = centerh(0_i64:2_i64, cell)
          Rx = center(0_i64) - vertexn(0_i64, i)
          Ry = center(1_i64) - vertexn(1_i64, i)
          I_xx(i) = I_xx(i) + Rx * Rx
          I_yy(i) = I_yy(i) + Ry * Ry
          I_xy(i) = I_xy(i) + Rx * Ry
          R_x(i) = R_x(i) + Rx
          R_y(i) = R_y(i) + Ry
          number(i) = number(i) + 1_i64
        end do
      end if
      D = I_xx(i) * I_yy(i) - I_xy(i) * I_xy(i)
      lambda_x(i) = (I_xy(i) * R_y(i) - I_yy(i) * R_x(i)) / D
      lambda_y(i) = (I_xy(i) * R_x(i) - I_xx(i) * R_y(i)) / D
    end do
    if (allocated(I_xx)) then
      deallocate(I_xx)
    end if
    if (allocated(I_yy)) then
      deallocate(I_yy)
    end if
    if (allocated(I_xy)) then
      deallocate(I_xy)
    end if

  end subroutine variables
  !........................................

  !........................................
  subroutine variables_3d(centerc, cellidn, haloidn, periodicn, vertexn, &
        namen, centergn, halocenterg, centerh, nbproc, R_x, R_y, R_z, &
        lambda_x, lambda_y, lambda_z, number, shift)

    implicit none

    real(f64), intent(in) :: centerc(0:,0:)
    integer(i64), intent(in) :: cellidn(0:,0:)
    integer(i64), intent(in) :: haloidn(0:,0:)
    integer(i64), intent(in) :: periodicn(0:,0:)
    real(f64), intent(in) :: vertexn(0:,0:)
    integer(i64), intent(in) :: namen(0:)
    real(f64), intent(in) :: centergn(0:,0:,0:)
    real(f64), intent(in) :: halocenterg(0:,0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    integer(i64), value :: nbproc
    real(f64), intent(inout) :: R_x(0:)
    real(f64), intent(inout) :: R_y(0:)
    real(f64), intent(inout) :: R_z(0:)
    real(f64), intent(inout) :: lambda_x(0:)
    real(f64), intent(inout) :: lambda_y(0:)
    real(f64), intent(inout) :: lambda_z(0:)
    integer(i64), intent(inout) :: number(0:)
    real(f64), intent(in) :: shift(0:,0:)
    integer(i64) :: nbnode
    real(f64), allocatable :: I_xx(:)
    real(f64), allocatable :: I_yy(:)
    real(f64), allocatable :: I_zz(:)
    real(f64), allocatable :: I_xy(:)
    real(f64), allocatable :: I_xz(:)
    real(f64), allocatable :: I_yz(:)
    real(f64) :: center(0:2_i64)
    integer(i64) :: i
    real(f64) :: D
    integer(i64) :: j
    real(f64) :: Rx
    real(f64) :: Ry
    real(f64) :: Rz
    integer(i64) :: cell

    nbnode = size(R_x, kind=i64)
    allocate(I_xx(0:nbnode - 1_i64))
    I_xx = 0.0_f64
    allocate(I_yy(0:nbnode - 1_i64))
    I_yy = 0.0_f64
    allocate(I_zz(0:nbnode - 1_i64))
    I_zz = 0.0_f64
    allocate(I_xy(0:nbnode - 1_i64))
    I_xy = 0.0_f64
    allocate(I_xz(0:nbnode - 1_i64))
    I_xz = 0.0_f64
    allocate(I_yz(0:nbnode - 1_i64))
    I_yz = 0.0_f64
    center = 0.0_f64
    do i = 0_i64, nbnode - 1_i64, 1_i64
      do j = 0_i64, cellidn(size(cellidn, 1_i64, i64) - 1_i64, i) - &
            1_i64, 1_i64
        center(:) = centerc(0_i64:2_i64, cellidn(j, i))
        Rx = center(0_i64) - vertexn(0_i64, i)
        Ry = center(1_i64) - vertexn(1_i64, i)
        Rz = center(2_i64) - vertexn(2_i64, i)
        I_xx(i) = I_xx(i) + Rx * Rx
        I_yy(i) = I_yy(i) + Ry * Ry
        I_zz(i) = I_zz(i) + Rz * Rz
        I_xy(i) = I_xy(i) + Rx * Ry
        I_xz(i) = I_xz(i) + Rx * Rz
        I_yz(i) = I_yz(i) + Ry * Rz
        R_x(i) = R_x(i) + Rx
        R_y(i) = R_y(i) + Ry
        R_z(i) = R_z(i) + Rz
        number(i) = number(i) + 1_i64
      end do
      if (centergn(3_i64, 0_i64, i) /= -1_i64) then
        do j = 0_i64, size(centergn(:, :, i),2,i64) - 1_i64, 1_i64
          cell = Int(centergn(size(centergn, 1_i64, i64) - 1_i64, j, i), &
                i64)
          if (cell /= -1_i64) then
            center(:) = centergn(0_i64:2_i64, j, i)
            Rx = center(0_i64) - vertexn(0_i64, i)
            Ry = center(1_i64) - vertexn(1_i64, i)
            Rz = center(2_i64) - vertexn(2_i64, i)
            I_xx(i) = I_xx(i) + Rx * Rx
            I_yy(i) = I_yy(i) + Ry * Ry
            I_zz(i) = I_zz(i) + Rz * Rz
            I_xy(i) = I_xy(i) + Rx * Ry
            I_xz(i) = I_xz(i) + Rx * Rz
            I_yz(i) = I_yz(i) + Ry * Rz
            R_x(i) = R_x(i) + Rx
            R_y(i) = R_y(i) + Ry
            R_z(i) = R_z(i) + Rz
            number(i) = number(i) + 1_i64
          end if
        end do
      end if
      !periodic boundary old vertex names)
      if (vertexn(3_i64, i) == 11_i64 .or. vertexn(3_i64, i) == 22_i64) &
            then
        do j = 0_i64, periodicn(size(periodicn, 1_i64, i64) - 1_i64, i) &
              - 1_i64, 1_i64
          cell = periodicn(j, i)
          center(0_i64) = centerc(0_i64, cell) + shift(0_i64, cell)
          center(1_i64) = centerc(1_i64, cell)
          center(2_i64) = centerc(2_i64, cell)
          Rx = center(0_i64) - vertexn(0_i64, i)
          Ry = center(1_i64) - vertexn(1_i64, i)
          Rz = center(2_i64) - vertexn(2_i64, i)
          I_xx(i) = I_xx(i) + Rx * Rx
          I_yy(i) = I_yy(i) + Ry * Ry
          I_zz(i) = I_zz(i) + Rz * Rz
          I_xy(i) = I_xy(i) + Rx * Ry
          I_xz(i) = I_xz(i) + Rx * Rz
          I_yz(i) = I_yz(i) + Ry * Rz
          R_x(i) = R_x(i) + Rx
          R_y(i) = R_y(i) + Ry
          R_z(i) = R_z(i) + Rz
          number(i) = number(i) + 1_i64
        end do
      else if (vertexn(3_i64, i) == 33_i64 .or. vertexn(3_i64, i) == &
            44_i64) then
        do j = 0_i64, periodicn(size(periodicn, 1_i64, i64) - 1_i64, i) &
              - 1_i64, 1_i64
          cell = periodicn(j, i)
          center(0_i64) = centerc(0_i64, cell)
          center(1_i64) = centerc(1_i64, cell) + shift(1_i64, cell)
          center(2_i64) = centerc(2_i64, cell)
          Rx = center(0_i64) - vertexn(0_i64, i)
          Ry = center(1_i64) - vertexn(1_i64, i)
          Rz = center(2_i64) - vertexn(2_i64, i)
          I_xx(i) = I_xx(i) + Rx * Rx
          I_yy(i) = I_yy(i) + Ry * Ry
          I_zz(i) = I_zz(i) + Rz * Rz
          I_xy(i) = I_xy(i) + Rx * Ry
          I_xz(i) = I_xz(i) + Rx * Rz
          I_yz(i) = I_yz(i) + Ry * Rz
          R_x(i) = R_x(i) + Rx
          R_y(i) = R_y(i) + Ry
          R_z(i) = R_z(i) + Rz
          number(i) = number(i) + 1_i64
        end do
      else if (vertexn(3_i64, i) == 55_i64 .or. vertexn(3_i64, i) == &
            66_i64) then
        do j = 0_i64, periodicn(size(periodicn, 1_i64, i64) - 1_i64, i) &
              - 1_i64, 1_i64
          cell = periodicn(j, i)
          center(0_i64) = centerc(0_i64, cell)
          center(1_i64) = centerc(1_i64, cell)
          center(2_i64) = centerc(2_i64, cell) + shift(2_i64, cell)
          Rx = center(0_i64) - vertexn(0_i64, i)
          Ry = center(1_i64) - vertexn(1_i64, i)
          Rz = center(2_i64) - vertexn(2_i64, i)
          I_xx(i) = I_xx(i) + Rx * Rx
          I_yy(i) = I_yy(i) + Ry * Ry
          I_zz(i) = I_zz(i) + Rz * Rz
          I_xy(i) = I_xy(i) + Rx * Ry
          I_xz(i) = I_xz(i) + Rx * Rz
          I_yz(i) = I_yz(i) + Ry * Rz
          R_x(i) = R_x(i) + Rx
          R_y(i) = R_y(i) + Ry
          R_z(i) = R_z(i) + Rz
          number(i) = number(i) + 1_i64
        end do
      end if
      if (namen(i) == 10_i64) then
        do j = 0_i64, haloidn(size(haloidn, 1_i64, i64) - 1_i64, i) - &
              1_i64, 1_i64
          cell = haloidn(j, i)
          center(:) = centerh(0_i64:2_i64, cell)
          Rx = center(0_i64) - vertexn(0_i64, i)
          Ry = center(1_i64) - vertexn(1_i64, i)
          Rz = center(2_i64) - vertexn(2_i64, i)
          I_xx(i) = I_xx(i) + Rx * Rx
          I_yy(i) = I_yy(i) + Ry * Ry
          I_zz(i) = I_zz(i) + Rz * Rz
          I_xy(i) = I_xy(i) + Rx * Ry
          I_xz(i) = I_xz(i) + Rx * Rz
          I_yz(i) = I_yz(i) + Ry * Rz
          R_x(i) = R_x(i) + Rx
          R_y(i) = R_y(i) + Ry
          R_z(i) = R_z(i) + Rz
          number(i) = number(i) + 1_i64
        end do
        !if namen[i] == 10:
        do j = 0_i64, size(halocenterg(:, :, i),2,i64) - 1_i64, 1_i64
          cell = Int(halocenterg(size(halocenterg, 1_i64, i64) - 1_i64, &
                j, i), i64)
          if (cell /= -1_i64) then
            center(:) = halocenterg(0_i64:2_i64, j, i)
            Rx = center(0_i64) - vertexn(0_i64, i)
            Ry = center(1_i64) - vertexn(1_i64, i)
            Rz = center(2_i64) - vertexn(2_i64, i)
            I_xx(i) = I_xx(i) + Rx * Rx
            I_yy(i) = I_yy(i) + Ry * Ry
            I_zz(i) = I_zz(i) + Rz * Rz
            I_xy(i) = I_xy(i) + Rx * Ry
            I_xz(i) = I_xz(i) + Rx * Rz
            I_yz(i) = I_yz(i) + Ry * Rz
            R_x(i) = R_x(i) + Rx
            R_y(i) = R_y(i) + Ry
            R_z(i) = R_z(i) + Rz
            number(i) = number(i) + 1_i64
          end if
        end do
      end if
      D = I_xx(i) * I_yy(i) * I_zz(i) + 2_i64 * I_xy(i) * I_xz(i) * I_yz &
            (i) - I_xx(i) * I_yz(i) * I_yz(i) - I_yy(i) * I_xz(i) * &
            I_xz(i) - I_zz(i) * I_xy(i) * I_xy(i)
      lambda_x(i) = ((I_yz(i) * I_yz(i) - I_yy(i) * I_zz(i)) * R_x(i) + &
            (I_xy(i) * I_zz(i) - I_xz(i) * I_yz(i)) * R_y(i) + (I_xz(i &
            ) * I_yy(i) - I_xy(i) * I_yz(i)) * R_z(i)) / D
      lambda_y(i) = ((I_xy(i) * I_zz(i) - I_xz(i) * I_yz(i)) * R_x(i) + &
            (I_xz(i) * I_xz(i) - I_xx(i) * I_zz(i)) * R_y(i) + (I_yz(i &
            ) * I_xx(i) - I_xz(i) * I_xy(i)) * R_z(i)) / D
      lambda_z(i) = ((I_xz(i) * I_yy(i) - I_xy(i) * I_yz(i)) * R_x(i) + &
            (I_yz(i) * I_xx(i) - I_xz(i) * I_xy(i)) * R_y(i) + (I_xy(i &
            ) * I_xy(i) - I_xx(i) * I_yy(i)) * R_z(i)) / D
    end do
    if (allocated(I_xx)) then
      deallocate(I_xx)
    end if
    if (allocated(I_yy)) then
      deallocate(I_yy)
    end if
    if (allocated(I_zz)) then
      deallocate(I_zz)
    end if
    if (allocated(I_xy)) then
      deallocate(I_xy)
    end if
    if (allocated(I_xz)) then
      deallocate(I_xz)
    end if
    if (allocated(I_yz)) then
      deallocate(I_yz)
    end if

  end subroutine variables_3d
  !........................................

end module pyccel_ddm
