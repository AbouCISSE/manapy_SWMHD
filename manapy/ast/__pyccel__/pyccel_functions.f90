module pyccel_functions


  use, intrinsic :: ISO_C_Binding, only : i32 => C_INT32_T , i64 => &
        C_INT64_T , f64 => C_DOUBLE
  implicit none

  contains

  !........................................
  !........................................

  !........................................
  !........................................

  !........................................
  !........................................

  !........................................
  !........................................

  !........................................
  !........................................

  !........................................
  !........................................

  !........................................
  !........................................

  !........................................
  subroutine haloghost_value_neumann(w_halo, w_haloghost, &
        haloghostcenter, BCindex, halonodes)

    implicit none

    real(f64), intent(in) :: w_halo(0:)
    real(f64), intent(inout) :: w_haloghost(0:)
    real(f64), intent(in) :: haloghostcenter(0:,0:,0:)
    integer(i64), value :: BCindex
    integer(i64), intent(in) :: halonodes(0:)
    integer(i64) :: i
    integer(i64) :: j
    integer(i64) :: cellhalo
    integer(i64) :: cellghost
    integer(i64) :: Dummy_0007

    do Dummy_0007 = 0_i64, size(halonodes, kind=i64) - 1_i64, 1_i64
      i = halonodes(Dummy_0007)
      do j = 0_i64, size(haloghostcenter(:, :, i),2,i64) - 1_i64, 1_i64
        if (haloghostcenter(size(haloghostcenter, 1_i64, i64) - 1_i64, j &
              , i) /= -1_i64) then
          if (haloghostcenter(size(haloghostcenter, 1_i64, i64) - 2_i64, &
                j, i) == BCindex) then
            cellhalo = Int(haloghostcenter(size(haloghostcenter, 1_i64, &
                  i64) - 3_i64, j, i), i64)
            cellghost = Int(haloghostcenter(size(haloghostcenter, 1_i64, &
                  i64) - 1_i64, j, i), i64)
            w_haloghost(cellghost) = w_halo(cellhalo)
          end if
        end if
      end do
    end do

  end subroutine haloghost_value_neumann
  !........................................

  !........................................
  subroutine haloghost_value_dirichlet(value, w_haloghost, &
        haloghostcenter, BCindex, halonodes)

    implicit none

    real(f64), intent(in) :: value(0:)
    real(f64), intent(inout) :: w_haloghost(0:)
    real(f64), intent(in) :: haloghostcenter(0:,0:,0:)
    integer(i64), value :: BCindex
    integer(i64), intent(in) :: halonodes(0:)
    integer(i64) :: i
    integer(i64) :: j
    integer(i64) :: cellghost
    integer(i64) :: Dummy_0008

    do Dummy_0008 = 0_i64, size(halonodes, kind=i64) - 1_i64, 1_i64
      i = halonodes(Dummy_0008)
      do j = 0_i64, size(haloghostcenter(:, :, i),2,i64) - 1_i64, 1_i64
        if (haloghostcenter(size(haloghostcenter, 1_i64, i64) - 1_i64, j &
              , i) /= -1_i64) then
          if (haloghostcenter(size(haloghostcenter, 1_i64, i64) - 2_i64, &
                j, i) == BCindex) then
            cellghost = Int(haloghostcenter(size(haloghostcenter, 1_i64, &
                  i64) - 1_i64, j, i), i64)
            w_haloghost(cellghost) = value(cellghost)
          end if
        end if
      end do
    end do

  end subroutine haloghost_value_dirichlet
  !........................................

  !........................................
  subroutine haloghost_value_nonslip(w_halo, w_haloghost, &
        haloghostcenter, BCindex, halonodes)

    implicit none

    real(f64), intent(in) :: w_halo(0:)
    real(f64), intent(inout) :: w_haloghost(0:)
    real(f64), intent(in) :: haloghostcenter(0:,0:,0:)
    integer(i64), value :: BCindex
    integer(i64), intent(in) :: halonodes(0:)
    integer(i64) :: i
    integer(i64) :: j
    integer(i64) :: cellghost
    integer(i64) :: Dummy_0009

    do Dummy_0009 = 0_i64, size(halonodes, kind=i64) - 1_i64, 1_i64
      i = halonodes(Dummy_0009)
      do j = 0_i64, size(haloghostcenter(:, :, i),2,i64) - 1_i64, 1_i64
        if (haloghostcenter(size(haloghostcenter, 1_i64, i64) - 1_i64, j &
              , i) /= -1_i64) then
          if (haloghostcenter(size(haloghostcenter, 1_i64, i64) - 2_i64, &
                j, i) == BCindex) then
            cellghost = Int(haloghostcenter(size(haloghostcenter, 1_i64, &
                  i64) - 1_i64, j, i), i64)
            w_haloghost(cellghost) = (-1_i64) * w_halo(cellghost)
          end if
        end if
      end do
    end do

  end subroutine haloghost_value_nonslip
  !........................................

  !........................................
  subroutine haloghost_value_slip(u_halo, v_halo, w_haloghost, &
        haloghostcenter, BCindex, halonodes, haloghostfaceinfo)

    implicit none

    real(f64), intent(in) :: u_halo(0:)
    real(f64), intent(in) :: v_halo(0:)
    real(f64), intent(inout) :: w_haloghost(0:)
    real(f64), intent(in) :: haloghostcenter(0:,0:,0:)
    integer(i64), value :: BCindex
    integer(i64), intent(in) :: halonodes(0:)
    real(f64), intent(in) :: haloghostfaceinfo(0:,0:,0:)
    real(f64), allocatable :: s_n(:)
    integer(i64) :: i
    integer(i64) :: j
    integer(i64) :: cellghost
    real(f64) :: u_i
    real(f64) :: v_i
    real(f64) :: mesure
    real(f64) :: u_g
    integer(i64) :: Dummy_0010

    allocate(s_n(0:1_i64))
    s_n = 0.0_f64
    do Dummy_0010 = 0_i64, size(halonodes, kind=i64) - 1_i64, 1_i64
      i = halonodes(Dummy_0010)
      do j = 0_i64, size(haloghostcenter(:, :, i),2,i64) - 1_i64, 1_i64
        if (haloghostcenter(size(haloghostcenter, 1_i64, i64) - 1_i64, j &
              , i) /= -1_i64) then
          if (haloghostcenter(size(haloghostcenter, 1_i64, i64) - 2_i64, &
                j, i) == BCindex) then
            cellghost = Int(haloghostcenter(size(haloghostcenter, 1_i64, &
                  i64) - 1_i64, j, i), i64)
            u_i = u_halo(cellghost)
            v_i = v_halo(cellghost)
            mesure = sqrt(haloghostfaceinfo(2_i64, j, i) ** 2_i64 + &
                  haloghostfaceinfo(3_i64, j, i) ** 2_i64)
            s_n(0_i64) = haloghostfaceinfo(2_i64, j, i) / mesure
            s_n(1_i64) = haloghostfaceinfo(3_i64, j, i) / mesure
            u_g = u_i * (s_n(1_i64) * s_n(1_i64) - s_n(0_i64) * s_n( &
                  0_i64)) - 2.0_f64 * v_i * s_n(0_i64) * s_n(1_i64)
            w_haloghost(i) = u_halo(cellghost) * u_g
          end if
        end if
      end do
    end do
    if (allocated(s_n)) then
      deallocate(s_n)
    end if

  end subroutine haloghost_value_slip
  !........................................

  !........................................
  subroutine cell_gradient_2d(w_c, w_ghost, w_halo, w_haloghost, centerc &
        , cellnid, halonid, nodecid, periodicn, periodic, namen, &
        centerg, halocenterg, vertexn, centerh, shift, nbproc, w_x, w_y &
        , w_z)

    implicit none

    real(f64), intent(in) :: w_c(0:)
    real(f64), intent(in) :: w_ghost(0:)
    real(f64), intent(in) :: w_halo(0:)
    real(f64), intent(in) :: w_haloghost(0:)
    real(f64), intent(in) :: centerc(0:,0:)
    integer(i64), intent(in) :: cellnid(0:,0:)
    integer(i64), intent(in) :: halonid(0:,0:)
    integer(i64), intent(in) :: nodecid(0:,0:)
    integer(i64), intent(in) :: periodicn(0:,0:)
    integer(i64), intent(in) :: periodic(0:,0:)
    integer(i64), intent(in) :: namen(0:)
    real(f64), intent(in) :: centerg(0:,0:,0:)
    real(f64), intent(in) :: halocenterg(0:,0:,0:)
    real(f64), intent(in) :: vertexn(0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    real(f64), intent(in) :: shift(0:,0:)
    integer(i64), value :: nbproc
    real(f64), intent(inout) :: w_x(0:)
    real(f64), intent(inout) :: w_y(0:)
    real(f64), intent(inout) :: w_z(0:)
    real(f64) :: center(0:2_i64)
    integer(i64) :: nbelement
    integer(i64) :: i
    real(f64) :: i_xx
    real(f64) :: i_yy
    real(f64) :: i_xy
    real(f64) :: j_xw
    real(f64) :: j_yw
    real(f64) :: dia
    integer(i64) :: j
    integer(i64) :: cell
    real(f64) :: j_x
    real(f64) :: j_y
    integer(i64) :: k
    integer(i64) :: nod

    center = 0.0_f64
    nbelement = size(w_c, kind=i64)
    do i = 0_i64, nbelement - 1_i64, 1_i64
      i_xx = 0.0_f64
      i_yy = 0.0_f64
      i_xy = 0.0_f64
      j_xw = 0.0_f64
      j_yw = 0.0_f64
      do j = 0_i64, cellnid(size(cellnid, 1_i64, i64) - 1_i64, i) - &
            1_i64, 1_i64
        cell = cellnid(j, i)
        j_x = centerc(0_i64, cell) - centerc(0_i64, i)
        j_y = centerc(1_i64, cell) - centerc(1_i64, i)
        i_xx = i_xx + j_x * j_x
        i_yy = i_yy + j_y * j_y
        i_xy = i_xy + j_x * j_y
        j_xw = j_xw + j_x * (w_c(cell) - w_c(i))
        j_yw = j_yw + j_y * (w_c(cell) - w_c(i))
      end do
      do k = 0_i64, 2_i64, 1_i64
        nod = nodecid(k, i)
        if (vertexn(3_i64, nod) == 11_i64 .or. vertexn(3_i64, nod) == &
              22_i64) then
          do j = 0_i64, periodic(size(periodic, 1_i64, i64) - 1_i64, nod &
                ) - 1_i64, 1_i64
            cell = Int(periodic(j, nod), i64)
            center(:) = centerc(0_i64:2_i64, cell)
            j_x = center(0_i64) + shift(0_i64, cell) - centerc(0_i64, i)
            j_y = center(1_i64) - centerc(1_i64, i)
            i_xx = i_xx + j_x * j_x
            i_yy = i_yy + j_y * j_y
            i_xy = i_xy + j_x * j_y
            j_xw = j_xw + j_x * (w_c(cell) - w_c(i))
            j_yw = j_yw + j_y * (w_c(cell) - w_c(i))
          end do
        end if
        if (vertexn(3_i64, nod) == 33_i64 .or. vertexn(3_i64, nod) == &
              44_i64) then
          do j = 0_i64, periodic(size(periodic, 1_i64, i64) - 1_i64, nod &
                ) - 1_i64, 1_i64
            cell = Int(periodic(j, nod), i64)
            center(:) = centerc(0_i64:2_i64, cell)
            j_x = center(0_i64) - centerc(0_i64, i)
            j_y = center(1_i64) + shift(1_i64, cell) - centerc(1_i64, i)
            i_xx = i_xx + j_x * j_x
            i_yy = i_yy + j_y * j_y
            i_xy = i_xy + j_x * j_y
            j_xw = j_xw + j_x * (w_c(cell) - w_c(i))
            j_yw = j_yw + j_y * (w_c(cell) - w_c(i))
          end do
        end if
      end do
      do j = 0_i64, halonid(size(halonid, 1_i64, i64) - 1_i64, i) - &
            1_i64, 1_i64
        cell = halonid(j, i)
        j_x = centerh(0_i64, cell) - centerc(0_i64, i)
        j_y = centerh(1_i64, cell) - centerc(1_i64, i)
        i_xx = i_xx + j_x * j_x
        i_yy = i_yy + j_y * j_y
        i_xy = i_xy + j_x * j_y
        j_xw = j_xw + j_x * (w_halo(cell) - w_c(i))
        j_yw = j_yw + j_y * (w_halo(cell) - w_c(i))
      end do
      do k = 0_i64, 2_i64, 1_i64
        nod = nodecid(k, i)
        if (vertexn(3_i64, nod) <= 4_i64) then
          do j = 0_i64, size(centerg(:, :, nod),2,i64) - 1_i64, 1_i64
            cell = Int(centerg(size(centerg, 1_i64, i64) - 1_i64, j, nod &
                  ), i64)
            if (cell /= -1_i64) then
              center(:) = centerg(0_i64:2_i64, j, nod)
              j_x = center(0_i64) - centerc(0_i64, i)
              j_y = center(1_i64) - centerc(1_i64, i)
              i_xx = i_xx + j_x * j_x
              i_yy = i_yy + j_y * j_y
              i_xy = i_xy + j_x * j_y
              j_xw = j_xw + j_x * (w_ghost(cell) - w_c(i))
              j_yw = j_yw + j_y * (w_ghost(cell) - w_c(i))
            end if
          end do
        end if
        do j = 0_i64, size(halocenterg(:, :, nod),2,i64) - 1_i64, 1_i64
          !-3 the index of global face
          cell = Int(halocenterg(size(halocenterg, 1_i64, i64) - 1_i64, &
                j, nod), i64)
          if (cell /= -1_i64) then
            center(:) = halocenterg(0_i64:2_i64, j, nod)
            j_x = center(0_i64) - centerc(0_i64, i)
            j_y = center(1_i64) - centerc(1_i64, i)
            i_xx = i_xx + j_x * j_x
            i_yy = i_yy + j_y * j_y
            i_xy = i_xy + j_x * j_y
            j_xw = j_xw + j_x * (w_haloghost(cell) - w_c(i))
            j_yw = j_yw + j_y * (w_haloghost(cell) - w_c(i))
          end if
        end do
      end do
      dia = i_xx * i_yy - i_xy * i_xy
      w_x(i) = (i_yy * j_xw - i_xy * j_yw) / dia
      w_y(i) = (i_xx * j_yw - i_xy * j_xw) / dia
      w_z(i) = 0.0_f64
    end do

  end subroutine cell_gradient_2d
  !........................................

  !........................................
  subroutine cell_gradient_3d(w_c, w_ghost, w_halo, w_haloghost, centerc &
        , cellnid, halonid, nodecid, periodicn, periodic, namen, &
        centerg, halocenterg, vertexn, centerh, shift, nbproc, w_x, w_y &
        , w_z)

    implicit none

    real(f64), intent(in) :: w_c(0:)
    real(f64), intent(in) :: w_ghost(0:)
    real(f64), intent(in) :: w_halo(0:)
    real(f64), intent(in) :: w_haloghost(0:)
    real(f64), intent(in) :: centerc(0:,0:)
    integer(i64), intent(in) :: cellnid(0:,0:)
    integer(i64), intent(in) :: halonid(0:,0:)
    integer(i64), intent(in) :: nodecid(0:,0:)
    integer(i64), intent(in) :: periodicn(0:,0:)
    integer(i64), intent(in) :: periodic(0:,0:)
    integer(i64), intent(in) :: namen(0:)
    real(f64), intent(in) :: centerg(0:,0:,0:)
    real(f64), intent(in) :: halocenterg(0:,0:,0:)
    real(f64), intent(in) :: vertexn(0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    real(f64), intent(in) :: shift(0:,0:)
    integer(i64), value :: nbproc
    real(f64), intent(inout) :: w_x(0:)
    real(f64), intent(inout) :: w_y(0:)
    real(f64), intent(inout) :: w_z(0:)
    integer(i64) :: nbelement
    real(f64) :: center(0:2_i64)
    integer(i64) :: i
    real(f64) :: i_xx
    real(f64) :: i_yy
    real(f64) :: i_zz
    real(f64) :: i_xy
    real(f64) :: i_xz
    real(f64) :: i_yz
    real(f64) :: j_x
    real(f64) :: j_y
    real(f64) :: j_z
    real(f64) :: dia
    integer(i64) :: j
    integer(i64) :: cell
    real(f64) :: jx
    real(f64) :: jy
    real(f64) :: jz
    integer(i64) :: k
    integer(i64) :: nod

    nbelement = size(w_c, kind=i64)
    center = 0.0_f64
    do i = 0_i64, nbelement - 1_i64, 1_i64
      i_xx = 0.0_f64
      i_yy = 0.0_f64
      i_zz = 0.0_f64
      i_xy = 0.0_f64
      i_xz = 0.0_f64
      i_yz = 0.0_f64
      j_x = 0.0_f64
      j_y = 0.0_f64
      j_z = 0.0_f64
      do j = 0_i64, cellnid(size(cellnid, 1_i64, i64) - 1_i64, i) - &
            1_i64, 1_i64
        cell = cellnid(j, i)
        jx = centerc(0_i64, cell) - centerc(0_i64, i)
        jy = centerc(1_i64, cell) - centerc(1_i64, i)
        jz = centerc(2_i64, cell) - centerc(2_i64, i)
        i_xx = i_xx + jx * jx
        i_yy = i_yy + jy * jy
        i_zz = i_zz + jz * jz
        i_xy = i_xy + jx * jy
        i_xz = i_xz + jx * jz
        i_yz = i_yz + jy * jz
        j_x = j_x + jx * (w_c(cell) - w_c(i))
        j_y = j_y + jy * (w_c(cell) - w_c(i))
        j_z = j_z + jz * (w_c(cell) - w_c(i))
      end do
      do j = 0_i64, periodicn(size(periodicn, 1_i64, i64) - 1_i64, i) - &
            1_i64, 1_i64
        cell = periodicn(j, i)
        center(:) = centerc(0_i64:2_i64, cell)
        j_x = center(0_i64) + shift(0_i64, cell) - centerc(0_i64, i)
        j_y = center(1_i64) + shift(1_i64, cell) - centerc(1_i64, i)
        j_z = center(2_i64) + shift(2_i64, cell) - centerc(2_i64, i)
        i_xx = i_xx + jx * jx
        i_yy = i_yy + jy * jy
        i_zz = i_zz + jz * jz
        i_xy = i_xy + jx * jy
        i_xz = i_xz + jx * jz
        i_yz = i_yz + jy * jz
        j_x = j_x + jx * (w_c(cell) - w_c(i))
        j_y = j_y + jy * (w_c(cell) - w_c(i))
        j_z = j_z + jz * (w_c(cell) - w_c(i))
      end do
      !if nbproc > 1:
      do j = 0_i64, halonid(size(halonid, 1_i64, i64) - 1_i64, i) - &
            1_i64, 1_i64
        cell = halonid(j, i)
        jx = centerh(0_i64, cell) - centerc(0_i64, i)
        jy = centerh(1_i64, cell) - centerc(1_i64, i)
        jz = centerh(2_i64, cell) - centerc(2_i64, i)
        i_xx = i_xx + jx * jx
        i_yy = i_yy + jy * jy
        i_zz = i_zz + jz * jz
        i_xy = i_xy + jx * jy
        i_xz = i_xz + jx * jz
        i_yz = i_yz + jy * jz
        j_x = j_x + jx * (w_halo(cell) - w_c(i))
        j_y = j_y + jy * (w_halo(cell) - w_c(i))
        j_z = j_z + jz * (w_halo(cell) - w_c(i))
      end do
      !TODO verify ghost center
      do k = 0_i64, 3_i64, 1_i64
        nod = nodecid(k, i)
        if (vertexn(3_i64, nod) <= 6_i64) then
          do j = 0_i64, size(centerg(:, :, nod),2,i64) - 1_i64, 1_i64
            cell = Int(centerg(size(centerg, 1_i64, i64) - 1_i64, j, nod &
                  ), i64)
            if (cell /= -1_i64) then
              center(:) = centerg(0_i64:2_i64, j, nod)
              jx = center(0_i64) - centerc(0_i64, i)
              jy = center(1_i64) - centerc(1_i64, i)
              jz = center(2_i64) - centerc(2_i64, i)
              i_xx = i_xx + jx * jx
              i_yy = i_yy + jy * jy
              i_zz = i_zz + jz * jz
              i_xy = i_xy + jx * jy
              i_xz = i_xz + jx * jz
              i_yz = i_yz + jy * jz
              j_x = j_x + jx * (w_ghost(cell) - w_c(i))
              j_y = j_y + jy * (w_ghost(cell) - w_c(i))
              j_z = j_z + jz * (w_ghost(cell) - w_c(i))
            end if
          end do
        end if
        !if namen[nod] == 10 :
        do j = 0_i64, size(halocenterg(:, :, nod),2,i64) - 1_i64, 1_i64
          !-3 the index of global face
          cell = Int(halocenterg(size(halocenterg, 1_i64, i64) - 1_i64, &
                j, nod), i64)
          if (cell /= -1_i64) then
            center(:) = halocenterg(0_i64:2_i64, j, nod)
            jx = center(0_i64) - centerc(0_i64, i)
            jy = center(1_i64) - centerc(1_i64, i)
            jz = center(2_i64) - centerc(2_i64, i)
            i_xx = i_xx + jx * jx
            i_yy = i_yy + jy * jy
            i_zz = i_zz + jz * jz
            i_xy = i_xy + jx * jy
            i_xz = i_xz + jx * jz
            i_yz = i_yz + jy * jz
            j_x = j_x + jx * (w_haloghost(cell) - w_c(i))
            j_y = j_y + jy * (w_haloghost(cell) - w_c(i))
            j_z = j_z + jz * (w_haloghost(cell) - w_c(i))
          end if
        end do
      end do
      dia = i_xx * i_yy * i_zz + 2.0_f64 * i_xy * i_xz * i_yz - i_xx * &
            i_yz ** 2_i64 - i_yy * i_xz ** 2_i64 - i_zz * i_xy ** 2_i64
      w_x(i) = ((i_yy * i_zz - i_yz ** 2_i64) * j_x + (i_xz * i_yz - &
            i_xy * i_zz) * j_y + (i_xy * i_yz - i_xz * i_yy) * j_z) / &
            dia
      w_y(i) = ((i_xz * i_yz - i_xy * i_zz) * j_x + (i_xx * i_zz - i_xz &
            ** 2_i64) * j_y + (i_xy * i_xz - i_yz * i_xx) * j_z) / dia
      w_z(i) = ((i_xy * i_yz - i_xz * i_yy) * j_x + (i_xy * i_xz - i_yz &
            * i_xx) * j_y + (i_xx * i_yy - i_xy ** 2_i64) * j_z) / dia
    end do

  end subroutine cell_gradient_3d
  !........................................

  !........................................
  subroutine face_gradient_2d(w_c, w_ghost, w_halo, w_node, cellidf, &
        nodeidf, centergf, namef, halofid, centerc, centerh, vertexn, &
        airDiamond, normalf, f_1, f_2, f_3, f_4, shift, wx_face, &
        wy_face, wz_face, innerfaces, halofaces, dirichletfaces, &
        neumann, periodicfaces)

    implicit none

    real(f64), intent(in) :: w_c(0:)
    real(f64), intent(in) :: w_ghost(0:)
    real(f64), intent(in) :: w_halo(0:)
    real(f64), intent(in) :: w_node(0:)
    integer(i64), intent(in) :: cellidf(0:,0:)
    integer(i64), intent(in) :: nodeidf(0:,0:)
    real(f64), intent(in) :: centergf(0:,0:)
    integer(i64), intent(in) :: namef(0:)
    integer(i64), intent(in) :: halofid(0:)
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    real(f64), intent(in) :: vertexn(0:,0:)
    real(f64), intent(in) :: airDiamond(0:)
    real(f64), intent(in) :: normalf(0:,0:)
    real(f64), intent(in) :: f_1(0:,0:)
    real(f64), intent(in) :: f_2(0:,0:)
    real(f64), intent(in) :: f_3(0:,0:)
    real(f64), intent(in) :: f_4(0:,0:)
    real(f64), intent(in) :: shift(0:,0:)
    real(f64), intent(inout) :: wx_face(0:)
    real(f64), intent(inout) :: wy_face(0:)
    real(f64), intent(in) :: wz_face(0:)
    integer(i64), intent(in) :: innerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: dirichletfaces(0:)
    integer(i64), intent(in) :: neumann(0:)
    integer(i64), intent(in) :: periodicfaces(0:)
    integer(i64) :: i
    integer(i64) :: c_left
    integer(i64) :: c_right
    integer(i64) :: i_1
    integer(i64) :: i_2
    real(f64) :: vi1
    real(f64) :: vi2
    real(f64) :: vv1
    real(f64) :: vv2
    integer(i64) :: Dummy_0011
    integer(i64) :: Dummy_0012
    integer(i64) :: Dummy_0013
    integer(i64) :: Dummy_0014
    integer(i64) :: Dummy_0015

    do Dummy_0011 = 0_i64, size(innerfaces, kind=i64) - 1_i64, 1_i64
      i = innerfaces(Dummy_0011)
      c_left = cellidf(0_i64, i)
      c_right = cellidf(1_i64, i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      vi1 = w_node(i_1)
      vi2 = w_node(i_2)
      vv1 = w_c(c_left)
      vv2 = w_c(c_right)
      wx_face(i) = 1_i64 / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * f_1( &
            1_i64, i) + (vv1 + vi2) * f_2(1_i64, i) + (vi2 + vv2) * f_3 &
            (1_i64, i) + (vv2 + vi1) * f_4(1_i64, i))
      wy_face(i) = (-1_i64) / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * &
            f_1(0_i64, i) + (vv1 + vi2) * f_2(0_i64, i) + (vi2 + vv2) * &
            f_3(0_i64, i) + (vv2 + vi1) * f_4(0_i64, i))
    end do
    do Dummy_0012 = 0_i64, size(periodicfaces, kind=i64) - 1_i64, 1_i64
      i = periodicfaces(Dummy_0012)
      c_left = cellidf(0_i64, i)
      c_right = cellidf(1_i64, i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      vi1 = w_node(i_1)
      vi2 = w_node(i_2)
      vv1 = w_c(c_left)
      vv2 = w_c(c_right)
      wx_face(i) = 1_i64 / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * f_1( &
            1_i64, i) + (vv1 + vi2) * f_2(1_i64, i) + (vi2 + vv2) * f_3 &
            (1_i64, i) + (vv2 + vi1) * f_4(1_i64, i))
      wy_face(i) = (-1_i64) / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * &
            f_1(0_i64, i) + (vv1 + vi2) * f_2(0_i64, i) + (vi2 + vv2) * &
            f_3(0_i64, i) + (vv2 + vi1) * f_4(0_i64, i))
    end do
    do Dummy_0013 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0013)
      c_left = cellidf(0_i64, i)
      c_right = halofid(i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      vi1 = w_node(i_1)
      vi2 = w_node(i_2)
      vv1 = w_c(c_left)
      vv2 = w_halo(c_right)
      wx_face(i) = 1_i64 / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * f_1( &
            1_i64, i) + (vv1 + vi2) * f_2(1_i64, i) + (vi2 + vv2) * f_3 &
            (1_i64, i) + (vv2 + vi1) * f_4(1_i64, i))
      wy_face(i) = (-1_i64) / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * &
            f_1(0_i64, i) + (vv1 + vi2) * f_2(0_i64, i) + (vi2 + vv2) * &
            f_3(0_i64, i) + (vv2 + vi1) * f_4(0_i64, i))
    end do
    do Dummy_0014 = 0_i64, size(dirichletfaces, kind=i64) - 1_i64, 1_i64
      i = dirichletfaces(Dummy_0014)
      c_left = cellidf(0_i64, i)
      c_right = i
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      vi1 = w_node(i_1)
      vi2 = w_node(i_2)
      vv1 = w_c(c_left)
      vv2 = w_ghost(c_right)
      wx_face(i) = 1_i64 / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * f_1( &
            1_i64, i) + (vv1 + vi2) * f_2(1_i64, i) + (vi2 + vv2) * f_3 &
            (1_i64, i) + (vv2 + vi1) * f_4(1_i64, i))
      wy_face(i) = (-1_i64) / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * &
            f_1(0_i64, i) + (vv1 + vi2) * f_2(0_i64, i) + (vi2 + vv2) * &
            f_3(0_i64, i) + (vv2 + vi1) * f_4(0_i64, i))
    end do
    do Dummy_0015 = 0_i64, size(neumann, kind=i64) - 1_i64, 1_i64
      i = neumann(Dummy_0015)
      c_left = cellidf(0_i64, i)
      c_right = i
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      vi1 = w_node(i_1)
      vi2 = w_node(i_2)
      vv1 = w_c(c_left)
      vv2 = w_ghost(c_right)
      wx_face(i) = 1_i64 / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * f_1( &
            1_i64, i) + (vv1 + vi2) * f_2(1_i64, i) + (vi2 + vv2) * f_3 &
            (1_i64, i) + (vv2 + vi1) * f_4(1_i64, i))
      wy_face(i) = (-1_i64) / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * &
            f_1(0_i64, i) + (vv1 + vi2) * f_2(0_i64, i) + (vi2 + vv2) * &
            f_3(0_i64, i) + (vv2 + vi1) * f_4(0_i64, i))
    end do

  end subroutine face_gradient_2d
  !........................................

  !........................................
  subroutine face_gradient_2d_uv(h_c, h_ghost, h_halo, h_node, w_c, &
        w_ghost, w_halo, w_node, cellidf, nodeidf, centergf, namef, &
        halofid, centerc, centerh, vertexn, airDiamond, normalf, f_1, &
        f_2, f_3, f_4, shift, wx_face, wy_face, wz_face, innerfaces, &
        halofaces, dirichletfaces, neumann, periodicfaces)

    implicit none

    real(f64), intent(in) :: h_c(0:)
    real(f64), intent(in) :: h_ghost(0:)
    real(f64), intent(in) :: h_halo(0:)
    real(f64), intent(in) :: h_node(0:)
    real(f64), intent(in) :: w_c(0:)
    real(f64), intent(in) :: w_ghost(0:)
    real(f64), intent(in) :: w_halo(0:)
    real(f64), intent(in) :: w_node(0:)
    integer(i64), intent(in) :: cellidf(0:,0:)
    integer(i64), intent(in) :: nodeidf(0:,0:)
    real(f64), intent(in) :: centergf(0:,0:)
    integer(i64), intent(in) :: namef(0:)
    integer(i64), intent(in) :: halofid(0:)
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    real(f64), intent(in) :: vertexn(0:,0:)
    real(f64), intent(in) :: airDiamond(0:)
    real(f64), intent(in) :: normalf(0:,0:)
    real(f64), intent(in) :: f_1(0:,0:)
    real(f64), intent(in) :: f_2(0:,0:)
    real(f64), intent(in) :: f_3(0:,0:)
    real(f64), intent(in) :: f_4(0:,0:)
    real(f64), intent(in) :: shift(0:,0:)
    real(f64), intent(inout) :: wx_face(0:)
    real(f64), intent(inout) :: wy_face(0:)
    real(f64), intent(in) :: wz_face(0:)
    integer(i64), intent(in) :: innerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: dirichletfaces(0:)
    integer(i64), intent(in) :: neumann(0:)
    integer(i64), intent(in) :: periodicfaces(0:)
    integer(i64) :: i
    integer(i64) :: c_left
    integer(i64) :: c_right
    integer(i64) :: i_1
    integer(i64) :: i_2
    real(f64) :: vi1
    real(f64) :: vi2
    real(f64) :: vv1
    real(f64) :: vv2
    integer(i64) :: Dummy_0016
    integer(i64) :: Dummy_0017
    integer(i64) :: Dummy_0018
    integer(i64) :: Dummy_0019
    integer(i64) :: Dummy_0020

    do Dummy_0016 = 0_i64, size(innerfaces, kind=i64) - 1_i64, 1_i64
      i = innerfaces(Dummy_0016)
      c_left = cellidf(0_i64, i)
      c_right = cellidf(1_i64, i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      vi1 = w_node(i_1) / h_node(i_1)
      vi2 = w_node(i_2) / h_node(i_2)
      vv1 = w_c(c_left) / h_c(c_left)
      vv2 = w_c(c_right) / h_c(c_right)
      wx_face(i) = 1_i64 / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * f_1( &
            1_i64, i) + (vv1 + vi2) * f_2(1_i64, i) + (vi2 + vv2) * f_3 &
            (1_i64, i) + (vv2 + vi1) * f_4(1_i64, i))
      wy_face(i) = (-1_i64) / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * &
            f_1(0_i64, i) + (vv1 + vi2) * f_2(0_i64, i) + (vi2 + vv2) * &
            f_3(0_i64, i) + (vv2 + vi1) * f_4(0_i64, i))
    end do
    do Dummy_0017 = 0_i64, size(periodicfaces, kind=i64) - 1_i64, 1_i64
      i = periodicfaces(Dummy_0017)
      c_left = cellidf(0_i64, i)
      c_right = cellidf(1_i64, i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      vi1 = w_node(i_1) / h_node(i_1)
      vi2 = w_node(i_2) / h_node(i_2)
      vv1 = w_c(c_left) / h_c(c_left)
      vv2 = w_c(c_right) / h_c(c_right)
      wx_face(i) = 1_i64 / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * f_1( &
            1_i64, i) + (vv1 + vi2) * f_2(1_i64, i) + (vi2 + vv2) * f_3 &
            (1_i64, i) + (vv2 + vi1) * f_4(1_i64, i))
      wy_face(i) = (-1_i64) / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * &
            f_1(0_i64, i) + (vv1 + vi2) * f_2(0_i64, i) + (vi2 + vv2) * &
            f_3(0_i64, i) + (vv2 + vi1) * f_4(0_i64, i))
    end do
    do Dummy_0018 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0018)
      c_left = cellidf(0_i64, i)
      c_right = halofid(i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      vi1 = w_node(i_1) / h_node(i_1)
      vi2 = w_node(i_2) / h_node(i_2)
      vv1 = w_c(c_left) / h_c(c_left)
      vv2 = w_halo(c_right) / h_halo(c_right)
      wx_face(i) = 1_i64 / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * f_1( &
            1_i64, i) + (vv1 + vi2) * f_2(1_i64, i) + (vi2 + vv2) * f_3 &
            (1_i64, i) + (vv2 + vi1) * f_4(1_i64, i))
      wy_face(i) = (-1_i64) / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * &
            f_1(0_i64, i) + (vv1 + vi2) * f_2(0_i64, i) + (vi2 + vv2) * &
            f_3(0_i64, i) + (vv2 + vi1) * f_4(0_i64, i))
    end do
    do Dummy_0019 = 0_i64, size(dirichletfaces, kind=i64) - 1_i64, 1_i64
      i = dirichletfaces(Dummy_0019)
      c_left = cellidf(0_i64, i)
      c_right = i
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      vi1 = w_node(i_1) / h_node(i_1)
      vi2 = w_node(i_2) / h_node(i_2)
      vv1 = w_c(c_left) / h_c(c_left)
      vv2 = w_ghost(c_right) / h_ghost(c_right)
      wx_face(i) = 1_i64 / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * f_1( &
            1_i64, i) + (vv1 + vi2) * f_2(1_i64, i) + (vi2 + vv2) * f_3 &
            (1_i64, i) + (vv2 + vi1) * f_4(1_i64, i))
      wy_face(i) = (-1_i64) / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * &
            f_1(0_i64, i) + (vv1 + vi2) * f_2(0_i64, i) + (vi2 + vv2) * &
            f_3(0_i64, i) + (vv2 + vi1) * f_4(0_i64, i))
    end do
    do Dummy_0020 = 0_i64, size(neumann, kind=i64) - 1_i64, 1_i64
      i = neumann(Dummy_0020)
      c_left = cellidf(0_i64, i)
      c_right = i
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      vi1 = w_node(i_1) / h_node(i_1)
      vi2 = w_node(i_2) / h_node(i_2)
      vv1 = w_c(c_left) / h_c(c_left)
      vv2 = w_ghost(c_right) / h_ghost(c_right)
      wx_face(i) = 1_i64 / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * f_1( &
            1_i64, i) + (vv1 + vi2) * f_2(1_i64, i) + (vi2 + vv2) * f_3 &
            (1_i64, i) + (vv2 + vi1) * f_4(1_i64, i))
      wy_face(i) = (-1_i64) / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * &
            f_1(0_i64, i) + (vv1 + vi2) * f_2(0_i64, i) + (vi2 + vv2) * &
            f_3(0_i64, i) + (vv2 + vi1) * f_4(0_i64, i))
    end do

  end subroutine face_gradient_2d_uv
  !........................................

  !........................................
  subroutine face_gradient_3d(w_c, w_ghost, w_halo, w_node, cellidf, &
        nodeidf, centergf, namef, halofid, centerc, centerh, vertexn, &
        airDiamond, normalf, f_1, f_2, f_3, f_4, shift, wx_face, &
        wy_face, wz_face, innerfaces, halofaces, dirichletfaces, &
        neumann, periodicfaces)

    implicit none

    real(f64), intent(in) :: w_c(0:)
    real(f64), intent(in) :: w_ghost(0:)
    real(f64), intent(in) :: w_halo(0:)
    real(f64), intent(in) :: w_node(0:)
    integer(i64), intent(in) :: cellidf(0:,0:)
    integer(i64), intent(in) :: nodeidf(0:,0:)
    real(f64), intent(in) :: centergf(0:,0:)
    integer(i64), intent(in) :: namef(0:)
    integer(i64), intent(in) :: halofid(0:)
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    real(f64), intent(in) :: vertexn(0:,0:)
    real(f64), intent(in) :: airDiamond(0:)
    real(f64), intent(in) :: normalf(0:,0:)
    real(f64), intent(in) :: f_1(0:,0:)
    real(f64), intent(in) :: f_2(0:,0:)
    real(f64), intent(in) :: f_3(0:,0:)
    real(f64), intent(in) :: f_4(0:,0:)
    real(f64), intent(in) :: shift(0:,0:)
    real(f64), intent(inout) :: wx_face(0:)
    real(f64), intent(inout) :: wy_face(0:)
    real(f64), intent(inout) :: wz_face(0:)
    integer(i64), intent(in) :: innerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: dirichletfaces(0:)
    integer(i64), intent(in) :: neumann(0:)
    integer(i64), intent(in) :: periodicfaces(0:)
    integer(i64) :: i
    integer(i64) :: c_left
    integer(i64) :: c_right
    integer(i64) :: i_1
    integer(i64) :: i_2
    integer(i64) :: i_3
    integer(i64) :: i_4
    real(f64) :: V_A
    real(f64) :: V_B
    real(f64) :: V_C
    real(f64) :: V_D
    real(f64) :: V_L
    real(f64) :: V_R
    integer(i64) :: Dummy_0021
    integer(i64) :: Dummy_0022
    integer(i64) :: Dummy_0023
    integer(i64) :: Dummy_0024
    integer(i64) :: Dummy_0025

    do Dummy_0021 = 0_i64, size(innerfaces, kind=i64) - 1_i64, 1_i64
      i = innerfaces(Dummy_0021)
      c_left = cellidf(0_i64, i)
      c_right = cellidf(1_i64, i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      i_3 = nodeidf(2_i64, i)
      i_4 = i_3
      V_A = w_node(i_1)
      V_B = w_node(i_2)
      V_C = w_node(i_3)
      V_D = w_node(i_4)
      V_L = w_c(c_left)
      V_R = w_c(c_right)
      wx_face(i) = (f_1(0_i64, i) * (V_A - V_C) + f_2(0_i64, i) * (V_B - &
            V_D) + normalf(0_i64, i) * (V_R - V_L)) / airDiamond(i)
      wy_face(i) = (f_1(1_i64, i) * (V_A - V_C) + f_2(1_i64, i) * (V_B - &
            V_D) + normalf(1_i64, i) * (V_R - V_L)) / airDiamond(i)
      wz_face(i) = (f_1(2_i64, i) * (V_A - V_C) + f_2(2_i64, i) * (V_B - &
            V_D) + normalf(2_i64, i) * (V_R - V_L)) / airDiamond(i)
    end do
    do Dummy_0022 = 0_i64, size(periodicfaces, kind=i64) - 1_i64, 1_i64
      i = periodicfaces(Dummy_0022)
      c_left = cellidf(0_i64, i)
      c_right = cellidf(1_i64, i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      i_3 = nodeidf(2_i64, i)
      i_4 = i_3
      V_A = w_node(i_1)
      V_B = w_node(i_2)
      V_C = w_node(i_3)
      V_D = w_node(i_4)
      V_L = w_c(c_left)
      V_R = w_c(c_right)
      wx_face(i) = (f_1(0_i64, i) * (V_A - V_C) + f_2(0_i64, i) * (V_B - &
            V_D) + normalf(0_i64, i) * (V_R - V_L)) / airDiamond(i)
      wy_face(i) = (f_1(1_i64, i) * (V_A - V_C) + f_2(1_i64, i) * (V_B - &
            V_D) + normalf(1_i64, i) * (V_R - V_L)) / airDiamond(i)
      wz_face(i) = (f_1(2_i64, i) * (V_A - V_C) + f_2(2_i64, i) * (V_B - &
            V_D) + normalf(2_i64, i) * (V_R - V_L)) / airDiamond(i)
    end do
    do Dummy_0023 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0023)
      c_left = cellidf(0_i64, i)
      c_right = halofid(i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      i_3 = nodeidf(2_i64, i)
      i_4 = i_3
      V_A = w_node(i_1)
      V_B = w_node(i_2)
      V_C = w_node(i_3)
      V_D = w_node(i_4)
      V_L = w_c(c_left)
      V_R = w_halo(c_right)
      wx_face(i) = (f_1(0_i64, i) * (V_A - V_C) + f_2(0_i64, i) * (V_B - &
            V_D) + normalf(0_i64, i) * (V_R - V_L)) / airDiamond(i)
      wy_face(i) = (f_1(1_i64, i) * (V_A - V_C) + f_2(1_i64, i) * (V_B - &
            V_D) + normalf(1_i64, i) * (V_R - V_L)) / airDiamond(i)
      wz_face(i) = (f_1(2_i64, i) * (V_A - V_C) + f_2(2_i64, i) * (V_B - &
            V_D) + normalf(2_i64, i) * (V_R - V_L)) / airDiamond(i)
    end do
    do Dummy_0024 = 0_i64, size(dirichletfaces, kind=i64) - 1_i64, 1_i64
      i = dirichletfaces(Dummy_0024)
      c_left = cellidf(0_i64, i)
      c_right = i
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      i_3 = nodeidf(2_i64, i)
      i_4 = i_3
      V_A = w_node(i_1)
      V_B = w_node(i_2)
      V_C = w_node(i_3)
      V_D = w_node(i_4)
      V_L = w_c(c_left)
      V_R = w_ghost(c_right)
      wx_face(i) = (f_1(0_i64, i) * (V_A - V_C) + f_2(0_i64, i) * (V_B - &
            V_D) + normalf(0_i64, i) * (V_R - V_L)) / airDiamond(i)
      wy_face(i) = (f_1(1_i64, i) * (V_A - V_C) + f_2(1_i64, i) * (V_B - &
            V_D) + normalf(1_i64, i) * (V_R - V_L)) / airDiamond(i)
      wz_face(i) = (f_1(2_i64, i) * (V_A - V_C) + f_2(2_i64, i) * (V_B - &
            V_D) + normalf(2_i64, i) * (V_R - V_L)) / airDiamond(i)
    end do
    do Dummy_0025 = 0_i64, size(neumann, kind=i64) - 1_i64, 1_i64
      i = neumann(Dummy_0025)
      c_left = cellidf(0_i64, i)
      c_right = i
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      i_3 = nodeidf(2_i64, i)
      i_4 = i_3
      V_A = w_node(i_1)
      V_B = w_node(i_2)
      V_C = w_node(i_3)
      V_D = w_node(i_4)
      V_L = w_c(c_left)
      V_R = w_ghost(c_right)
      wx_face(i) = (f_1(0_i64, i) * (V_A - V_C) + f_2(0_i64, i) * (V_B - &
            V_D) + normalf(0_i64, i) * (V_R - V_L)) / airDiamond(i)
      wy_face(i) = (f_1(1_i64, i) * (V_A - V_C) + f_2(1_i64, i) * (V_B - &
            V_D) + normalf(1_i64, i) * (V_R - V_L)) / airDiamond(i)
      wz_face(i) = (f_1(2_i64, i) * (V_A - V_C) + f_2(2_i64, i) * (V_B - &
            V_D) + normalf(2_i64, i) * (V_R - V_L)) / airDiamond(i)
    end do

  end subroutine face_gradient_3d
  !........................................

  !........................................
  subroutine centertovertex_2d(w_c, w_ghost, w_halo, w_haloghost, &
        centerc, centerh, cellidn, periodicn, haloidn, vertexn, namen, &
        centergn, halocentergn, R_x, R_y, lambda_x, lambda_y, number, &
        shift, nbproc, w_n)

    implicit none

    real(f64), intent(in) :: w_c(0:)
    real(f64), intent(in) :: w_ghost(0:)
    real(f64), intent(in) :: w_halo(0:)
    real(f64), intent(in) :: w_haloghost(0:)
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    integer(i64), intent(in) :: cellidn(0:,0:)
    integer(i64), intent(in) :: periodicn(0:,0:)
    integer(i64), intent(in) :: haloidn(0:,0:)
    real(f64), intent(in) :: vertexn(0:,0:)
    integer(i64), intent(in) :: namen(0:)
    real(f64), intent(in) :: centergn(0:,0:,0:)
    real(f64), intent(in) :: halocentergn(0:,0:,0:)
    real(f64), intent(in) :: R_x(0:)
    real(f64), intent(in) :: R_y(0:)
    real(f64), intent(in) :: lambda_x(0:)
    real(f64), intent(in) :: lambda_y(0:)
    integer(i64), intent(in) :: number(0:)
    real(f64), intent(in) :: shift(0:,0:)
    integer(i64), value :: nbproc
    real(f64), intent(inout) :: w_n(0:)
    integer(i64) :: nbnode
    real(f64) :: center(0:2_i64)
    integer(i64) :: i
    integer(i64) :: j
    integer(i64) :: cell
    real(f64) :: xdiff
    real(f64) :: ydiff
    real(f64) :: alpha

    w_n(:) = 0.0_f64
    nbnode = size(vertexn,2,i64)
    center = 0.0_f64
    do i = 0_i64, nbnode - 1_i64, 1_i64
      do j = 0_i64, cellidn(size(cellidn, 1_i64, i64) - 1_i64, i) - &
            1_i64, 1_i64
        cell = cellidn(j, i)
        center(:) = centerc(:, cell)
        xdiff = center(0_i64) - vertexn(0_i64, i)
        ydiff = center(1_i64) - vertexn(1_i64, i)
        alpha = (1.0_f64 + lambda_x(i) * xdiff + lambda_y(i) * ydiff) / &
              (number(i) + lambda_x(i) * R_x(i) + lambda_y(i) * R_y(i))
        w_n(i) = w_n(i) + alpha * w_c(cell)
      end do
      if (centergn(2_i64, 0_i64, i) /= -1_i64) then
        do j = 0_i64, size(centergn(:, :, i),2,i64) - 1_i64, 1_i64
          cell = Int(centergn(size(centergn, 1_i64, i64) - 1_i64, j, i), &
                i64)
          if (cell /= -1_i64) then
            center(:) = centergn(0_i64:2_i64, j, i)
            xdiff = center(0_i64) - vertexn(0_i64, i)
            ydiff = center(1_i64) - vertexn(1_i64, i)
            alpha = (1.0_f64 + lambda_x(i) * xdiff + lambda_y(i) * ydiff &
                  ) / (number(i) + lambda_x(i) * R_x(i) + lambda_y(i) * &
                  R_y(i))
            w_n(i) = w_n(i) + alpha * w_ghost(cell)
          end if
        end do
      end if
      !TODO Must be keeped like that checked ok ;)
      if (vertexn(3_i64, i) == 11_i64 .or. vertexn(3_i64, i) == 22_i64) &
            then
        do j = 0_i64, periodicn(size(periodicn, 1_i64, i64) - 1_i64, i) &
              - 1_i64, 1_i64
          cell = periodicn(j, i)
          center(:) = centerc(0_i64:2_i64, cell)
          xdiff = center(0_i64) + shift(0_i64, cell) - vertexn(0_i64, i)
          ydiff = center(1_i64) - vertexn(1_i64, i)
          alpha = (1.0_f64 + lambda_x(i) * xdiff + lambda_y(i) * ydiff) &
                / (number(i) + lambda_x(i) * R_x(i) + lambda_y(i) * R_y &
                (i))
          w_n(i) = w_n(i) + alpha * w_c(cell)
        end do
      else if (vertexn(3_i64, i) == 33_i64 .or. vertexn(3_i64, i) == &
            44_i64) then
        do j = 0_i64, periodicn(size(periodicn, 1_i64, i64) - 1_i64, i) &
              - 1_i64, 1_i64
          cell = periodicn(j, i)
          center(:) = centerc(0_i64:2_i64, cell)
          xdiff = center(0_i64) - vertexn(0_i64, i)
          ydiff = center(1_i64) + shift(1_i64, cell) - vertexn(1_i64, i)
          alpha = (1.0_f64 + lambda_x(i) * xdiff + lambda_y(i) * ydiff) &
                / (number(i) + lambda_x(i) * R_x(i) + lambda_y(i) * R_y &
                (i))
          w_n(i) = w_n(i) + alpha * w_c(cell)
        end do
      end if
      if (namen(i) == 10_i64) then
        do j = 0_i64, size(halocentergn(:, :, i),2,i64) - 1_i64, 1_i64
          cell = Int(halocentergn(size(halocentergn, 1_i64, i64) - 1_i64 &
                , j, i), i64)
          if (cell /= -1_i64) then
            center(:) = halocentergn(0_i64:2_i64, j, i)
            xdiff = center(0_i64) - vertexn(0_i64, i)
            ydiff = center(1_i64) - vertexn(1_i64, i)
            alpha = (1.0_f64 + lambda_x(i) * xdiff + lambda_y(i) * ydiff &
                  ) / (number(i) + lambda_x(i) * R_x(i) + lambda_y(i) * &
                  R_y(i))
            w_n(i) = w_n(i) + alpha * w_haloghost(cell)
          end if
        end do
        !if haloidn[i][-1] > 0 :
        do j = 0_i64, haloidn(size(haloidn, 1_i64, i64) - 1_i64, i) - &
              1_i64, 1_i64
          cell = haloidn(j, i)
          center(:) = centerh(0_i64:2_i64, cell)
          xdiff = center(0_i64) - vertexn(0_i64, i)
          ydiff = center(1_i64) - vertexn(1_i64, i)
          alpha = (1.0_f64 + lambda_x(i) * xdiff + lambda_y(i) * ydiff) &
                / (number(i) + lambda_x(i) * R_x(i) + lambda_y(i) * R_y &
                (i))
          w_n(i) = w_n(i) + alpha * w_halo(cell)
        end do
      end if
    end do

  end subroutine centertovertex_2d
  !........................................

  !........................................
  subroutine centertovertex_3d(w_c, w_ghost, w_halo, w_haloghost, &
        centerc, centerh, cellidn, periodicn, haloidn, vertexn, namen, &
        centergn, halocentergn, R_x, R_y, R_z, lambda_x, lambda_y, &
        lambda_z, number, shift, nbproc, w_n)

    implicit none

    real(f64), intent(in) :: w_c(0:)
    real(f64), intent(in) :: w_ghost(0:)
    real(f64), intent(in) :: w_halo(0:)
    real(f64), intent(in) :: w_haloghost(0:)
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    integer(i64), intent(in) :: cellidn(0:,0:)
    integer(i64), intent(in) :: periodicn(0:,0:)
    integer(i64), intent(in) :: haloidn(0:,0:)
    real(f64), intent(in) :: vertexn(0:,0:)
    integer(i64), intent(in) :: namen(0:)
    real(f64), intent(in) :: centergn(0:,0:,0:)
    real(f64), intent(in) :: halocentergn(0:,0:,0:)
    real(f64), intent(in) :: R_x(0:)
    real(f64), intent(in) :: R_y(0:)
    real(f64), intent(in) :: R_z(0:)
    real(f64), intent(in) :: lambda_x(0:)
    real(f64), intent(in) :: lambda_y(0:)
    real(f64), intent(in) :: lambda_z(0:)
    integer(i64), intent(in) :: number(0:)
    real(f64), intent(in) :: shift(0:,0:)
    integer(i64), value :: nbproc
    real(f64), intent(inout) :: w_n(0:)
    integer(i64) :: nbnode
    real(f64) :: center(0:2_i64)
    integer(i64) :: i
    integer(i64) :: j
    integer(i64) :: cell
    real(f64) :: xdiff
    real(f64) :: ydiff
    real(f64) :: zdiff
    real(f64) :: alpha

    w_n(:) = 0.0_f64
    nbnode = size(vertexn,2,i64)
    center = 0.0_f64
    do i = 0_i64, nbnode - 1_i64, 1_i64
      do j = 0_i64, cellidn(size(cellidn, 1_i64, i64) - 1_i64, i) - &
            1_i64, 1_i64
        cell = cellidn(j, i)
        center(:) = centerc(:, cell)
        xdiff = center(0_i64) - vertexn(0_i64, i)
        ydiff = center(1_i64) - vertexn(1_i64, i)
        zdiff = center(2_i64) - vertexn(2_i64, i)
        alpha = (1.0_f64 + lambda_x(i) * xdiff + lambda_y(i) * ydiff + &
              lambda_z(i) * zdiff) / (number(i) + lambda_x(i) * R_x(i) &
              + lambda_y(i) * R_y(i) + lambda_z(i) * R_z(i))
        w_n(i) = w_n(i) + alpha * w_c(cell)
      end do
      if (centergn(3_i64, 0_i64, i) /= -1_i64) then
        do j = 0_i64, size(centergn(:, :, i),2,i64) - 1_i64, 1_i64
          cell = Int(centergn(size(centergn, 1_i64, i64) - 1_i64, j, i), &
                i64)
          if (cell /= -1_i64) then
            center(:) = centergn(0_i64:2_i64, j, i)
            xdiff = center(0_i64) - vertexn(0_i64, i)
            ydiff = center(1_i64) - vertexn(1_i64, i)
            zdiff = center(2_i64) - vertexn(2_i64, i)
            alpha = (1.0_f64 + lambda_x(i) * xdiff + lambda_y(i) * ydiff &
                  + lambda_z(i) * zdiff) / (number(i) + lambda_x(i) * &
                  R_x(i) + lambda_y(i) * R_y(i) + lambda_z(i) * R_z(i))
            w_n(i) = w_n(i) + alpha * w_ghost(cell)
          end if
        end do
      end if
      if (vertexn(3_i64, i) == 11_i64 .or. vertexn(3_i64, i) == 22_i64) &
            then
        do j = 0_i64, periodicn(size(periodicn, 1_i64, i64) - 1_i64, i) &
              - 1_i64, 1_i64
          cell = periodicn(j, i)
          center(:) = centerc(0_i64:2_i64, cell)
          xdiff = center(0_i64) + shift(0_i64, cell) - vertexn(0_i64, i)
          ydiff = center(1_i64) - vertexn(1_i64, i)
          zdiff = center(2_i64) - vertexn(2_i64, i)
          alpha = (1.0_f64 + lambda_x(i) * xdiff + lambda_y(i) * ydiff + &
                lambda_z(i) * zdiff) / (number(i) + lambda_x(i) * R_x(i &
                ) + lambda_y(i) * R_y(i) + lambda_z(i) * R_z(i))
          w_n(i) = w_n(i) + alpha * w_c(cell)
        end do
      else if (vertexn(3_i64, i) == 33_i64 .or. vertexn(3_i64, i) == &
            44_i64) then
        do j = 0_i64, periodicn(size(periodicn, 1_i64, i64) - 1_i64, i) &
              - 1_i64, 1_i64
          cell = periodicn(j, i)
          center(:) = centerc(0_i64:2_i64, cell)
          xdiff = center(0_i64) - vertexn(0_i64, i)
          ydiff = center(1_i64) + shift(1_i64, cell) - vertexn(1_i64, i)
          zdiff = center(2_i64) - vertexn(2_i64, i)
          alpha = (1.0_f64 + lambda_x(i) * xdiff + lambda_y(i) * ydiff + &
                lambda_z(i) * zdiff) / (number(i) + lambda_x(i) * R_x(i &
                ) + lambda_y(i) * R_y(i) + lambda_z(i) * R_z(i))
          w_n(i) = w_n(i) + alpha * w_c(cell)
        end do
      else if (vertexn(3_i64, i) == 55_i64 .or. vertexn(3_i64, i) == &
            66_i64) then
        do j = 0_i64, periodicn(size(periodicn, 1_i64, i64) - 1_i64, i) &
              - 1_i64, 1_i64
          cell = periodicn(j, i)
          center(:) = centerc(0_i64:2_i64, cell)
          xdiff = center(0_i64) - vertexn(0_i64, i)
          ydiff = center(1_i64) - vertexn(1_i64, i)
          zdiff = center(2_i64) + shift(2_i64, cell) - vertexn(2_i64, i)
          alpha = (1.0_f64 + lambda_x(i) * xdiff + lambda_y(i) * ydiff + &
                lambda_z(i) * zdiff) / (number(i) + lambda_x(i) * R_x(i &
                ) + lambda_y(i) * R_y(i) + lambda_z(i) * R_z(i))
          w_n(i) = w_n(i) + alpha * w_c(cell)
        end do
      end if
      if (namen(i) == 10_i64) then
        do j = 0_i64, size(halocentergn(:, :, i),2,i64) - 1_i64, 1_i64
          cell = Int(halocentergn(size(halocentergn, 1_i64, i64) - 1_i64 &
                , j, i), i64)
          if (cell /= -1_i64) then
            center(:) = halocentergn(0_i64:2_i64, j, i)
            xdiff = center(0_i64) - vertexn(0_i64, i)
            ydiff = center(1_i64) - vertexn(1_i64, i)
            zdiff = center(2_i64) - vertexn(2_i64, i)
            alpha = (1.0_f64 + lambda_x(i) * xdiff + lambda_y(i) * ydiff &
                  + lambda_z(i) * zdiff) / (number(i) + lambda_x(i) * &
                  R_x(i) + lambda_y(i) * R_y(i) + lambda_z(i) * R_z(i))
            w_n(i) = w_n(i) + alpha * w_haloghost(cell)
          end if
        end do
        !if haloidn[i][-1] > 0 :
        do j = 0_i64, haloidn(size(haloidn, 1_i64, i64) - 1_i64, i) - &
              1_i64, 1_i64
          cell = haloidn(j, i)
          center(:) = centerh(0_i64:2_i64, cell)
          xdiff = center(0_i64) - vertexn(0_i64, i)
          ydiff = center(1_i64) - vertexn(1_i64, i)
          zdiff = center(2_i64) - vertexn(2_i64, i)
          alpha = (1.0_f64 + lambda_x(i) * xdiff + lambda_y(i) * ydiff + &
                lambda_z(i) * zdiff) / (number(i) + lambda_x(i) * R_x(i &
                ) + lambda_y(i) * R_y(i) + lambda_z(i) * R_z(i))
          w_n(i) = w_n(i) + alpha * w_halo(cell)
        end do
      end if
    end do

  end subroutine centertovertex_3d
  !........................................

  !........................................
  subroutine barthlimiter_2d(w_c, w_ghost, w_halo, w_x, w_y, w_z, psi, &
        cellid, faceid, namef, halofid, centerc, centerf)

    implicit none

    real(f64), intent(in) :: w_c(0:)
    real(f64), intent(in) :: w_ghost(0:)
    real(f64), intent(in) :: w_halo(0:)
    real(f64), intent(in) :: w_x(0:)
    real(f64), intent(in) :: w_y(0:)
    real(f64), intent(in) :: w_z(0:)
    real(f64), intent(inout) :: psi(0:)
    integer(i64), intent(in) :: cellid(0:,0:)
    integer(i64), intent(in) :: faceid(0:,0:)
    integer(i64), intent(in) :: namef(0:)
    integer(i64), intent(in) :: halofid(0:)
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: centerf(0:,0:)
    integer(i64) :: nbelement
    real(f64) :: val
    integer(i64) :: i
    real(f64) :: w_max
    real(f64) :: w_min
    integer(i64) :: j
    integer(i64) :: face
    real(f64) :: r_xyz1
    real(f64) :: r_xyz2
    real(f64) :: delta2
    real(f64) :: psi_ij
    real(f64) :: value

    nbelement = size(w_c, kind=i64)
    val = 1.0_f64
    psi(:) = val
    do i = 0_i64, nbelement - 1_i64, 1_i64
      w_max = w_c(i)
      w_min = w_c(i)
      do j = 0_i64, 2_i64, 1_i64
        face = faceid(j, i)
        if (namef(face) == 0_i64 .or. namef(face) > 10_i64) then
          !11 or namef[face] == 22 or namef[face] == 33 or namef[face] == 44:
          w_max = maxval([w_max, w_c(cellid(0_i64, face)), w_c(cellid( &
                1_i64, face))])
          w_min = minval([w_min, w_c(cellid(0_i64, face)), w_c(cellid( &
                1_i64, face))])
        else if (namef(face) == 1_i64 .or. namef(face) == 2_i64 .or. &
              namef(face) == 3_i64 .or. namef(face) == 4_i64) then
          w_max = maxval([w_max, w_c(cellid(0_i64, face)), w_ghost(face) &
                ])
          w_min = minval([w_min, w_c(cellid(0_i64, face)), w_ghost(face) &
                ])
        else
          w_max = maxval([w_max, w_c(cellid(0_i64, face)), w_halo( &
                halofid(face))])
          w_min = minval([w_min, w_c(cellid(0_i64, face)), w_halo( &
                halofid(face))])
        end if
      end do
      do j = 0_i64, 2_i64, 1_i64
        face = faceid(j, i)
        r_xyz1 = centerf(0_i64, face) - centerc(0_i64, i)
        r_xyz2 = centerf(1_i64, face) - centerc(1_i64, i)
        delta2 = w_x(i) * r_xyz1 + w_y(i) * r_xyz2
        !TODO choice of epsilon
        if (abs(delta2) < 1e-08_f64) then
          psi_ij = 1.0_f64
        else
          if (delta2 > 0.0_f64) then
            value = (w_max - w_c(i)) / delta2
            psi_ij = minval([val, value])
          end if
          if (delta2 < 0.0_f64) then
            value = (w_min - w_c(i)) / delta2
            psi_ij = minval([val, value])
          end if
        end if
        psi(i) = minval([psi(i), psi_ij])
      end do
    end do

  end subroutine barthlimiter_2d
  !........................................

  !........................................
  subroutine barthlimiter_3d(h_c, h_ghost, h_halo, h_x, h_y, h_z, psi, &
        cellid, faceid, namef, halofid, centerc, centerf)

    implicit none

    real(f64), intent(in) :: h_c(0:)
    real(f64), intent(in) :: h_ghost(0:)
    real(f64), intent(in) :: h_halo(0:)
    real(f64), intent(in) :: h_x(0:)
    real(f64), intent(in) :: h_y(0:)
    real(f64), intent(in) :: h_z(0:)
    real(f64), intent(inout) :: psi(0:)
    integer(i64), intent(in) :: cellid(0:,0:)
    integer(i64), intent(in) :: faceid(0:,0:)
    integer(i64), intent(in) :: namef(0:)
    integer(i64), intent(in) :: halofid(0:)
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: centerf(0:,0:)
    integer(i64) :: nbelement
    integer(i64) :: i
    real(f64) :: w_max
    real(f64) :: w_min
    integer(i64) :: j
    integer(i64) :: face
    real(f64) :: r_xyz1
    real(f64) :: r_xyz2
    real(f64) :: r_xyz3
    real(f64) :: delta2
    real(f64) :: psi_ij
    real(f64) :: value

    nbelement = size(h_c, kind=i64)
    psi(:) = 1.0_f64
    do i = 0_i64, nbelement - 1_i64, 1_i64
      w_max = h_c(i)
      w_min = h_c(i)
      do j = 0_i64, 3_i64, 1_i64
        face = faceid(j, i)
        if (namef(face) == 0_i64 .or. namef(face) > 10_i64) then
          w_max = maxval([w_max, h_c(cellid(0_i64, face)), h_c(cellid( &
                1_i64, face))])
          w_min = minval([w_min, h_c(cellid(0_i64, face)), h_c(cellid( &
                1_i64, face))])
        else if (namef(face) == 10_i64) then
          w_max = maxval([w_max, h_c(cellid(0_i64, face)), h_halo( &
                halofid(face))])
          w_min = minval([w_min, h_c(cellid(0_i64, face)), h_halo( &
                halofid(face))])
        else
          w_max = maxval([w_max, h_c(cellid(0_i64, face)), h_ghost(face) &
                ])
          w_min = minval([w_min, h_c(cellid(0_i64, face)), h_ghost(face) &
                ])
        end if
      end do
      do j = 0_i64, 3_i64, 1_i64
        face = faceid(j, i)
        r_xyz1 = centerf(0_i64, face) - centerc(0_i64, i)
        r_xyz2 = centerf(1_i64, face) - centerc(1_i64, i)
        r_xyz3 = centerf(2_i64, face) - centerc(2_i64, i)
        delta2 = h_x(i) * r_xyz1 + h_y(i) * r_xyz2 + h_z(i) * r_xyz3
        !TODO choice of epsilon
        if (abs(delta2) < 1e-10_f64) then
          psi_ij = 1.0_f64
        else
          if (delta2 > 0.0_f64) then
            value = (w_max - h_c(i)) / delta2
            psi_ij = minval([1.0_f64, value])
          end if
          if (delta2 < 0.0_f64) then
            value = (w_min - h_c(i)) / delta2
            psi_ij = minval([1.0_f64, value])
          end if
        end if
        psi(i) = minval([psi(i), psi_ij])
      end do
    end do

  end subroutine barthlimiter_3d
  !........................................

  !........................................
  function search_element(a, target_value) result(find)

    implicit none

    integer(i64) :: find
    integer(i64), intent(in) :: a(0:)
    integer(i64), value :: target_value
    integer(i64) :: val
    integer(i64) :: Dummy_0026

    find = 0_i64
    do Dummy_0026 = 0_i64, size(a, kind=i64) - 1_i64, 1_i64
      val = a(Dummy_0026)
      if (val == target_value) then
        find = 1_i64
        exit
      end if
    end do
    return

  end function search_element
  !........................................

  !........................................
  subroutine get_triplet_2d(cellfid, nodeidf, vertexn, halofid, haloext, &
        namen, oldnamen, volume, cellnid, centerc, centerh, halonid, &
        periodicnid, centergn, halocentergn, airDiamond, lambda_x, &
        lambda_y, number, R_x, R_y, param1, param2, param3, param4, &
        shift, nbelements, loctoglob, BCdirichlet, a_loc, irn_loc, &
        jcn_loc, matrixinnerfaces, halofaces, dirichletfaces)

    implicit none

    integer(i64), intent(in) :: cellfid(0:,0:)
    integer(i64), intent(in) :: nodeidf(0:,0:)
    real(f64), intent(in) :: vertexn(0:,0:)
    integer(i64), intent(in) :: halofid(0:)
    integer(i64), intent(in) :: haloext(0:,0:)
    integer(i64), intent(in) :: namen(0:)
    integer(i64), intent(in) :: oldnamen(0:)
    real(f64), intent(in) :: volume(0:)
    integer(i64), intent(in) :: cellnid(0:,0:)
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    integer(i64), intent(in) :: halonid(0:,0:)
    integer(i64), intent(in) :: periodicnid(0:,0:)
    real(f64), intent(in) :: centergn(0:,0:,0:)
    real(f64), intent(in) :: halocentergn(0:,0:,0:)
    real(f64), intent(in) :: airDiamond(0:)
    real(f64), intent(in) :: lambda_x(0:)
    real(f64), intent(in) :: lambda_y(0:)
    integer(i64), intent(in) :: number(0:)
    real(f64), intent(in) :: R_x(0:)
    real(f64), intent(in) :: R_y(0:)
    real(f64), intent(in) :: param1(0:)
    real(f64), intent(in) :: param2(0:)
    real(f64), intent(in) :: param3(0:)
    real(f64), intent(in) :: param4(0:)
    real(f64), intent(in) :: shift(0:,0:)
    integer(i64), value :: nbelements
    integer(i64), intent(in) :: loctoglob(0:)
    integer(i64), intent(in) :: BCdirichlet(0:)
    real(f64), intent(inout) :: a_loc(0:)
    integer(i32), intent(inout) :: irn_loc(0:)
    integer(i32), intent(inout) :: jcn_loc(0:)
    integer(i64), intent(in) :: matrixinnerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: dirichletfaces(0:)
    real(f64), allocatable :: center(:)
    real(f64), allocatable :: parameters(:)
    integer(i64) :: cmpt
    integer(i64) :: i
    integer(i64) :: c_left
    integer(i64) :: c_leftglob
    integer(i64) :: c_right
    integer(i64) :: c_rightglob
    real(f64) :: value
    integer(i64) :: cmptparam
    integer(i64) :: nod
    integer(i64) :: j
    real(f64) :: xdiff
    real(f64) :: ydiff
    real(f64) :: alpha
    integer(i64) :: index
    integer(i64) :: Dummy_0027
    integer(i64) :: Dummy_0028
    integer(i64) :: Dummy_0029
    integer(i64) :: Dummy_0030
    integer(i64) :: Dummy_0031

    allocate(center(0:1_i64))
    center = 0.0_f64
    allocate(parameters(0:1_i64))
    parameters = 0.0_f64
    cmpt = 0_i64
    do Dummy_0027 = 0_i64, size(matrixinnerfaces, kind=i64) - 1_i64, &
          1_i64
      i = matrixinnerfaces(Dummy_0027)
      c_left = cellfid(0_i64, i)
      c_leftglob = loctoglob(c_left)
      parameters(0_i64) = param4(i)
      parameters(1_i64) = param2(i)
      c_right = cellfid(1_i64, i)
      c_rightglob = loctoglob(c_right)
      irn_loc(cmpt) = c_leftglob
      jcn_loc(cmpt) = c_leftglob
      value = param1(i) / volume(c_left)
      a_loc(cmpt) = value
      cmpt = cmpt + 1_i64
      cmptparam = 0_i64
      do Dummy_0028 = 0_i64, size(nodeidf(:, i), kind=i64) - 1_i64, &
            1_i64
        nod = nodeidf(Dummy_0028, i)
        if (search_element(BCdirichlet, oldnamen(nod)) == 0_i64) then
          do j = 0_i64, cellnid(size(cellnid, 1_i64, i64) - 1_i64, nod) &
                - 1_i64, 1_i64
            center(:) = centerc(0_i64:1_i64, cellnid(j, nod))
            xdiff = center(0_i64) - vertexn(0_i64, nod)
            ydiff = center(1_i64) - vertexn(1_i64, nod)
            alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                  ydiff) / (number(nod) + lambda_x(nod) * R_x(nod) + &
                  lambda_y(nod) * R_y(nod))
            value = alpha / volume(c_left) * parameters(cmptparam)
            irn_loc(cmpt) = c_leftglob
            jcn_loc(cmpt) = loctoglob(cellnid(j, nod))
            a_loc(cmpt) = value
            cmpt = cmpt + 1_i64
            !right cell-----------------------------------
            value = (-1.0_f64) * alpha / volume(c_right) * parameters( &
                  cmptparam)
            irn_loc(cmpt) = c_rightglob
            jcn_loc(cmpt) = loctoglob(cellnid(j, nod))
            a_loc(cmpt) = value
            cmpt = cmpt + 1_i64
          end do
          do j = 0_i64, size(centergn(:, :, nod),2,i64) - 1_i64, 1_i64
            if (centergn(size(centergn, 1_i64, i64) - 1_i64, j, nod) /= &
                  -1_i64) then
              center(:) = centergn(0_i64:1_i64, j, nod)
              xdiff = center(0_i64) - vertexn(0_i64, nod)
              ydiff = center(1_i64) - vertexn(1_i64, nod)
              alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                    ydiff) / (number(nod) + lambda_x(nod) * R_x(nod) + &
                    lambda_y(nod) * R_y(nod))
              index = Int(centergn(2_i64, j, nod), i64)
              value = alpha / volume(c_left) * parameters(cmptparam)
              irn_loc(cmpt) = c_leftglob
              jcn_loc(cmpt) = loctoglob(index)
              a_loc(cmpt) = value
              cmpt = cmpt + 1_i64
              !right cell-----------------------------------
              value = (-1.0_f64) * alpha / volume(c_right) * parameters( &
                    cmptparam)
              irn_loc(cmpt) = c_rightglob
              jcn_loc(cmpt) = loctoglob(index)
              a_loc(cmpt) = value
              cmpt = cmpt + 1_i64
            end if
          end do
          do j = 0_i64, size(halocentergn(:, :, nod),2,i64) - 1_i64, &
                1_i64
            if (halocentergn(size(halocentergn, 1_i64, i64) - 1_i64, j, &
                  nod) /= -1_i64) then
              center(:) = halocentergn(0_i64:1_i64, j, nod)
              xdiff = center(0_i64) - vertexn(0_i64, nod)
              ydiff = center(1_i64) - vertexn(1_i64, nod)
              alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                    ydiff) / (number(nod) + lambda_x(nod) * R_x(nod) + &
                    lambda_y(nod) * R_y(nod))
              index = Int(halocentergn(2_i64, j, nod), i64)
              value = alpha / volume(c_left) * parameters(cmptparam)
              irn_loc(cmpt) = c_leftglob
              jcn_loc(cmpt) = haloext(0_i64, index)
              a_loc(cmpt) = value
              cmpt = cmpt + 1_i64
              !right cell-----------------------------------
              value = (-1.0_f64) * alpha / volume(c_right) * parameters( &
                    cmptparam)
              irn_loc(cmpt) = c_rightglob
              jcn_loc(cmpt) = haloext(0_i64, index)
              a_loc(cmpt) = value
              cmpt = cmpt + 1_i64
            end if
          end do
          do j = 0_i64, periodicnid(size(periodicnid, 1_i64, i64) - &
                1_i64, nod) - 1_i64, 1_i64
            if (vertexn(3_i64, nod) == 11_i64 .or. vertexn(3_i64, nod) &
                  == 22_i64) then
              center(0_i64) = centerc(0_i64, periodicnid(j, nod)) + &
                    shift(0_i64, periodicnid(j, nod))
              center(1_i64) = centerc(1_i64, periodicnid(j, nod))
            end if
            if (vertexn(3_i64, nod) == 33_i64 .or. vertexn(3_i64, nod) &
                  == 44_i64) then
              center(0_i64) = centerc(0_i64, periodicnid(j, nod))
              center(1_i64) = centerc(1_i64, periodicnid(j, nod)) + &
                    shift(1_i64, periodicnid(j, nod))
            end if
            xdiff = center(0_i64) - vertexn(0_i64, nod)
            ydiff = center(1_i64) - vertexn(1_i64, nod)
            alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                  ydiff) / (number(nod) + lambda_x(nod) * R_x(nod) + &
                  lambda_y(nod) * R_y(nod))
            value = alpha / volume(c_left) * parameters(cmptparam)
            irn_loc(cmpt) = c_leftglob
            jcn_loc(cmpt) = loctoglob(periodicnid(j, nod))
            a_loc(cmpt) = value
            cmpt = cmpt + 1_i64
            !right cell-----------------------------------
            value = (-1.0_f64) * alpha / volume(c_right) * parameters( &
                  cmptparam)
            irn_loc(cmpt) = c_rightglob
            jcn_loc(cmpt) = loctoglob(periodicnid(j, nod))
            a_loc(cmpt) = value
            cmpt = cmpt + 1_i64
          end do
          do j = 0_i64, halonid(size(halonid, 1_i64, i64) - 1_i64, nod) &
                - 1_i64, 1_i64
            center(:) = centerh(0_i64:1_i64, halonid(j, nod))
            xdiff = center(0_i64) - vertexn(0_i64, nod)
            ydiff = center(1_i64) - vertexn(1_i64, nod)
            alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                  ydiff) / (number(nod) + lambda_x(nod) * R_x(nod) + &
                  lambda_y(nod) * R_y(nod))
            value = alpha / volume(c_left) * parameters(cmptparam)
            irn_loc(cmpt) = c_leftglob
            jcn_loc(cmpt) = haloext(0_i64, halonid(j, nod))
            a_loc(cmpt) = value
            cmpt = cmpt + 1_i64
            !right cell-----------------------------------
            value = (-1.0_f64) * alpha / volume(c_right) * parameters( &
                  cmptparam)
            irn_loc(cmpt) = c_rightglob
            jcn_loc(cmpt) = haloext(0_i64, halonid(j, nod))
            a_loc(cmpt) = value
            cmpt = cmpt + 1_i64
          end do
        end if
        cmptparam = +1_i64
      end do
      irn_loc(cmpt) = c_leftglob
      jcn_loc(cmpt) = c_rightglob
      value = param3(i) / volume(c_left)
      a_loc(cmpt) = value
      cmpt = cmpt + 1_i64
      !right cell------------------------------------------------------
      irn_loc(cmpt) = c_rightglob
      jcn_loc(cmpt) = c_leftglob
      value = (-1.0_f64) * param1(i) / volume(c_right)
      a_loc(cmpt) = value
      cmpt = cmpt + 1_i64
      irn_loc(cmpt) = c_rightglob
      jcn_loc(cmpt) = c_rightglob
      value = (-1.0_f64) * param3(i) / volume(c_right)
      a_loc(cmpt) = value
      cmpt = cmpt + 1_i64
    end do
    do Dummy_0029 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0029)
      c_left = cellfid(0_i64, i)
      c_leftglob = loctoglob(c_left)
      parameters(0_i64) = param4(i)
      parameters(1_i64) = param2(i)
      c_rightglob = haloext(0_i64, halofid(i))
      c_right = halofid(i)
      irn_loc(cmpt) = c_leftglob
      jcn_loc(cmpt) = c_leftglob
      value = param1(i) / volume(c_left)
      a_loc(cmpt) = value
      cmpt = cmpt + 1_i64
      irn_loc(cmpt) = c_leftglob
      jcn_loc(cmpt) = c_rightglob
      value = param3(i) / volume(c_left)
      a_loc(cmpt) = value
      cmpt = cmpt + 1_i64
      cmptparam = 0_i64
      do Dummy_0030 = 0_i64, size(nodeidf(:, i), kind=i64) - 1_i64, &
            1_i64
        nod = nodeidf(Dummy_0030, i)
        if (search_element(BCdirichlet, oldnamen(nod)) == 0_i64) then
          do j = 0_i64, cellnid(size(cellnid, 1_i64, i64) - 1_i64, nod) &
                - 1_i64, 1_i64
            center(:) = centerc(0_i64:1_i64, cellnid(j, nod))
            xdiff = center(0_i64) - vertexn(0_i64, nod)
            ydiff = center(1_i64) - vertexn(1_i64, nod)
            alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                  ydiff) / (number(nod) + lambda_x(nod) * R_x(nod) + &
                  lambda_y(nod) * R_y(nod))
            value = alpha / volume(c_left) * parameters(cmptparam)
            irn_loc(cmpt) = c_leftglob
            jcn_loc(cmpt) = loctoglob(cellnid(j, nod))
            a_loc(cmpt) = value
            cmpt = cmpt + 1_i64
          end do
          do j = 0_i64, size(centergn(:, :, nod),2,i64) - 1_i64, 1_i64
            if (centergn(size(centergn, 1_i64, i64) - 1_i64, j, nod) /= &
                  -1_i64) then
              center(:) = centergn(0_i64:1_i64, j, nod)
              xdiff = center(0_i64) - vertexn(0_i64, nod)
              ydiff = center(1_i64) - vertexn(1_i64, nod)
              alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                    ydiff) / (number(nod) + lambda_x(nod) * R_x(nod) + &
                    lambda_y(nod) * R_y(nod))
              index = Int(centergn(2_i64, j, nod), i64)
              value = alpha / volume(c_left) * parameters(cmptparam)
              irn_loc(cmpt) = c_leftglob
              jcn_loc(cmpt) = loctoglob(index)
              a_loc(cmpt) = value
              cmpt = cmpt + 1_i64
            end if
          end do
          do j = 0_i64, size(halocentergn(:, :, nod),2,i64) - 1_i64, &
                1_i64
            if (halocentergn(size(halocentergn, 1_i64, i64) - 1_i64, j, &
                  nod) /= -1_i64) then
              center(:) = halocentergn(0_i64:1_i64, j, nod)
              xdiff = center(0_i64) - vertexn(0_i64, nod)
              ydiff = center(1_i64) - vertexn(1_i64, nod)
              alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                    ydiff) / (number(nod) + lambda_x(nod) * R_x(nod) + &
                    lambda_y(nod) * R_y(nod))
              index = Int(halocentergn(2_i64, j, nod), i64)
              value = alpha / volume(c_left) * parameters(cmptparam)
              irn_loc(cmpt) = c_leftglob
              jcn_loc(cmpt) = haloext(0_i64, index)
              a_loc(cmpt) = value
              cmpt = cmpt + 1_i64
            end if
          end do
          do j = 0_i64, halonid(size(halonid, 1_i64, i64) - 1_i64, nod) &
                - 1_i64, 1_i64
            center(:) = centerh(0_i64:1_i64, halonid(j, nod))
            xdiff = center(0_i64) - vertexn(0_i64, nod)
            ydiff = center(1_i64) - vertexn(1_i64, nod)
            alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                  ydiff) / (number(nod) + lambda_x(nod) * R_x(nod) + &
                  lambda_y(nod) * R_y(nod))
            value = alpha / volume(c_left) * parameters(cmptparam)
            irn_loc(cmpt) = c_leftglob
            jcn_loc(cmpt) = haloext(0_i64, halonid(j, nod))
            a_loc(cmpt) = value
            cmpt = cmpt + 1_i64
          end do
        end if
        cmptparam = cmptparam + 1_i64
      end do
    end do
    do Dummy_0031 = 0_i64, size(dirichletfaces, kind=i64) - 1_i64, 1_i64
      i = dirichletfaces(Dummy_0031)
      c_left = cellfid(0_i64, i)
      c_leftglob = loctoglob(c_left)
      irn_loc(cmpt) = c_leftglob
      jcn_loc(cmpt) = c_leftglob
      value = param1(i) / volume(c_left)
      a_loc(cmpt) = value
      cmpt = cmpt + 1_i64
      irn_loc(cmpt) = c_leftglob
      jcn_loc(cmpt) = c_leftglob
      value = (-1.0_f64) * param3(i) / volume(c_left)
      a_loc(cmpt) = value
      cmpt = cmpt + 1_i64
    end do
    if (allocated(center)) then
      deallocate(center)
    end if
    if (allocated(parameters)) then
      deallocate(parameters)
    end if

  end subroutine get_triplet_2d
  !........................................

  !........................................
  function compute_2dmatrix_size(nodeidf, halofid, cellnid, halonid, &
        periodicnid, centergn, halocentergn, oldnamen, BCdirichlet, &
        matrixinnerfaces, halofaces, dirichletfaces) result(cmpt)

    implicit none

    integer(i64) :: cmpt
    integer(i64), intent(in) :: nodeidf(0:,0:)
    integer(i64), intent(in) :: halofid(0:)
    integer(i64), intent(in) :: cellnid(0:,0:)
    integer(i64), intent(in) :: halonid(0:,0:)
    integer(i64), intent(in) :: periodicnid(0:,0:)
    real(f64), intent(in) :: centergn(0:,0:,0:)
    real(f64), intent(in) :: halocentergn(0:,0:,0:)
    integer(i64), intent(in) :: oldnamen(0:)
    integer(i64), intent(in) :: BCdirichlet(0:)
    integer(i64), intent(in) :: matrixinnerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: dirichletfaces(0:)
    integer(i64) :: i
    integer(i64) :: nod
    integer(i64) :: j
    integer(i64) :: Dummy_0032
    integer(i64) :: Dummy_0033
    integer(i64) :: Dummy_0034
    integer(i64) :: Dummy_0035
    integer(i64) :: Dummy_0036

    cmpt = 0_i64
    do Dummy_0032 = 0_i64, size(matrixinnerfaces, kind=i64) - 1_i64, &
          1_i64
      i = matrixinnerfaces(Dummy_0032)
      cmpt = cmpt + 1_i64
      do Dummy_0033 = 0_i64, size(nodeidf(:, i), kind=i64) - 1_i64, &
            1_i64
        nod = nodeidf(Dummy_0033, i)
        if (search_element(BCdirichlet, oldnamen(nod)) == 0_i64) then
          !if vertexn[nod][3] not in BCdirichlet:
          do j = 0_i64, cellnid(size(cellnid, 1_i64, i64) - 1_i64, nod) &
                - 1_i64, 1_i64
            cmpt = cmpt + 1_i64
            !right cell-----------------------------------
            cmpt = cmpt + 1_i64
          end do
          do j = 0_i64, size(centergn(:, :, nod),2,i64) - 1_i64, 1_i64
            if (centergn(size(centergn, 1_i64, i64) - 1_i64, j, nod) /= &
                  -1_i64) then
              cmpt = cmpt + 1_i64
              !right cell-----------------------------------
              cmpt = cmpt + 1_i64
            end if
          end do
          do j = 0_i64, size(halocentergn(:, :, nod),2,i64) - 1_i64, &
                1_i64
            if (halocentergn(size(halocentergn, 1_i64, i64) - 1_i64, j, &
                  nod) /= -1_i64) then
              cmpt = cmpt + 1_i64
              !right cell-----------------------------------
              cmpt = cmpt + 1_i64
            end if
          end do
          do j = 0_i64, periodicnid(size(periodicnid, 1_i64, i64) - &
                1_i64, nod) - 1_i64, 1_i64
            cmpt = cmpt + 1_i64
            !right cell-----------------------------------
            cmpt = cmpt + 1_i64
          end do
          do j = 0_i64, halonid(size(halonid, 1_i64, i64) - 1_i64, nod) &
                - 1_i64, 1_i64
            cmpt = cmpt + 1_i64
            !right cell-----------------------------------
            cmpt = cmpt + 1_i64
          end do
        end if
      end do
      cmpt = cmpt + 1_i64
      !right cell------------------------------------------------------
      cmpt = cmpt + 1_i64
      cmpt = cmpt + 1_i64
    end do
    !elif namef[i] == 10:
    do Dummy_0034 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0034)
      cmpt = cmpt + 1_i64
      cmpt = cmpt + 1_i64
      cmpt = cmpt + 1_i64
      do Dummy_0035 = 0_i64, size(nodeidf(:, i), kind=i64) - 1_i64, &
            1_i64
        nod = nodeidf(Dummy_0035, i)
        if (search_element(BCdirichlet, oldnamen(nod)) == 0_i64) then
          do j = 0_i64, cellnid(size(cellnid, 1_i64, i64) - 1_i64, nod) &
                - 1_i64, 1_i64
            cmpt = cmpt + 1_i64
          end do
          do j = 0_i64, size(centergn(:, :, nod),2,i64) - 1_i64, 1_i64
            if (centergn(size(centergn, 1_i64, i64) - 1_i64, j, nod) /= &
                  -1_i64) then
              cmpt = cmpt + 1_i64
            end if
          end do
          do j = 0_i64, size(halocentergn(:, :, nod),2,i64) - 1_i64, &
                1_i64
            if (halocentergn(size(halocentergn, 1_i64, i64) - 1_i64, j, &
                  nod) /= -1_i64) then
              cmpt = cmpt + 1_i64
            end if
          end do
          do j = 0_i64, halonid(size(halonid, 1_i64, i64) - 1_i64, nod) &
                - 1_i64, 1_i64
            cmpt = cmpt + 1_i64
          end do
        end if
      end do
    end do
    do Dummy_0036 = 0_i64, size(dirichletfaces, kind=i64) - 1_i64, 1_i64
      i = dirichletfaces(Dummy_0036)
      cmpt = cmpt + 1_i64
      cmpt = cmpt + 1_i64
    end do
    return

  end function compute_2dmatrix_size
  !........................................

  !........................................
  function compute_3dmatrix_size(nodeidf, halofid, cellnid, halonid, &
        periodicnid, centergn, halocentergn, oldnamen, BCdirichlet, &
        matrixinnerfaces, halofaces, dirichletfaces) result(cmpt)

    implicit none

    integer(i64) :: cmpt
    integer(i64), intent(in) :: nodeidf(0:,0:)
    integer(i64), intent(in) :: halofid(0:)
    integer(i64), intent(in) :: cellnid(0:,0:)
    integer(i64), intent(in) :: halonid(0:,0:)
    integer(i64), intent(in) :: periodicnid(0:,0:)
    real(f64), intent(in) :: centergn(0:,0:,0:)
    real(f64), intent(in) :: halocentergn(0:,0:,0:)
    integer(i64), intent(in) :: oldnamen(0:)
    integer(i64), intent(in) :: BCdirichlet(0:)
    integer(i64), intent(in) :: matrixinnerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: dirichletfaces(0:)
    integer(i64), allocatable :: nodes(:)
    integer(i64) :: i
    integer(i64) :: nod
    integer(i64) :: j
    integer(i64) :: Dummy_0037
    integer(i64) :: Dummy_0038
    integer(i64) :: Dummy_0039
    integer(i64) :: Dummy_0040
    integer(i64) :: Dummy_0041

    cmpt = 0_i64
    allocate(nodes(0:3_i64))
    nodes = 0_i64
    do Dummy_0037 = 0_i64, size(matrixinnerfaces, kind=i64) - 1_i64, &
          1_i64
      i = matrixinnerfaces(Dummy_0037)
      nodes(0_i64:2_i64) = nodeidf(:, i)
      nodes(3_i64) = nodeidf(2_i64, i)
      cmpt = cmpt + 1_i64
      do Dummy_0038 = 0_i64, size(nodes, kind=i64) - 1_i64, 1_i64
        nod = nodes(Dummy_0038)
        if (search_element(BCdirichlet, oldnamen(nod)) == 0_i64) then
          do j = 0_i64, cellnid(size(cellnid, 1_i64, i64) - 1_i64, nod) &
                - 1_i64, 1_i64
            cmpt = cmpt + 1_i64
            !right cell-----------------------------------
            cmpt = cmpt + 1_i64
          end do
          do j = 0_i64, size(centergn(:, :, nod),2,i64) - 1_i64, 1_i64
            if (centergn(size(centergn, 1_i64, i64) - 1_i64, j, nod) /= &
                  -1_i64) then
              cmpt = cmpt + 1_i64
              !right cell-----------------------------------
              cmpt = cmpt + 1_i64
            end if
          end do
          do j = 0_i64, size(halocentergn(:, :, nod),2,i64) - 1_i64, &
                1_i64
            if (halocentergn(size(halocentergn, 1_i64, i64) - 1_i64, j, &
                  nod) /= -1_i64) then
              cmpt = cmpt + 1_i64
              !right cell-----------------------------------
              cmpt = cmpt + 1_i64
            end if
          end do
          do j = 0_i64, periodicnid(size(periodicnid, 1_i64, i64) - &
                1_i64, nod) - 1_i64, 1_i64
            cmpt = cmpt + 1_i64
            !right cell-----------------------------------
            cmpt = cmpt + 1_i64
          end do
          do j = 0_i64, halonid(size(halonid, 1_i64, i64) - 1_i64, nod) &
                - 1_i64, 1_i64
            cmpt = cmpt + 1_i64
            !right cell-----------------------------------
            cmpt = cmpt + 1_i64
          end do
        end if
      end do
      cmpt = cmpt + 1_i64
      !right cell------------------------------------------------------
      cmpt = cmpt + 1_i64
      cmpt = cmpt + 1_i64
    end do
    do Dummy_0039 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0039)
      nodes(0_i64:2_i64) = nodeidf(:, i)
      nodes(3_i64) = nodeidf(2_i64, i)
      cmpt = cmpt + 1_i64
      cmpt = cmpt + 1_i64
      do Dummy_0040 = 0_i64, size(nodes, kind=i64) - 1_i64, 1_i64
        nod = nodes(Dummy_0040)
        if (search_element(BCdirichlet, oldnamen(nod)) == 0_i64) then
          do j = 0_i64, cellnid(size(cellnid, 1_i64, i64) - 1_i64, nod) &
                - 1_i64, 1_i64
            cmpt = cmpt + 1_i64
          end do
          do j = 0_i64, size(centergn(:, :, nod),2,i64) - 1_i64, 1_i64
            if (centergn(size(centergn, 1_i64, i64) - 1_i64, j, nod) /= &
                  -1_i64) then
              cmpt = cmpt + 1_i64
            end if
          end do
          do j = 0_i64, size(halocentergn(:, :, nod),2,i64) - 1_i64, &
                1_i64
            if (halocentergn(size(halocentergn, 1_i64, i64) - 1_i64, j, &
                  nod) /= -1_i64) then
              cmpt = cmpt + 1_i64
            end if
          end do
          do j = 0_i64, halonid(size(halonid, 1_i64, i64) - 1_i64, nod) &
                - 1_i64, 1_i64
            cmpt = cmpt + 1_i64
          end do
        end if
      end do
    end do
    do Dummy_0041 = 0_i64, size(dirichletfaces, kind=i64) - 1_i64, 1_i64
      i = dirichletfaces(Dummy_0041)
      cmpt = cmpt + 1_i64
      cmpt = cmpt + 1_i64
    end do
    if (allocated(nodes)) then
      deallocate(nodes)
    end if
    return

  end function compute_3dmatrix_size
  !........................................

  !........................................
  subroutine get_triplet_3d(cellfid, nodeidf, vertexn, halofid, haloext, &
        namen, oldnamen, volume, centergn, halocentergn, periodicnid, &
        cellnid, centerc, centerh, halonid, airDiamond, lambda_x, &
        lambda_y, lambda_z, number, R_x, R_y, R_z, param1, param2, &
        param3, shift, loctoglob, BCdirichlet, a_loc, irn_loc, jcn_loc, &
        matrixinnerfaces, halofaces, dirichletfaces)

    implicit none

    integer(i64), intent(in) :: cellfid(0:,0:)
    integer(i64), intent(in) :: nodeidf(0:,0:)
    real(f64), intent(in) :: vertexn(0:,0:)
    integer(i64), intent(in) :: halofid(0:)
    integer(i64), intent(in) :: haloext(0:,0:)
    integer(i64), intent(in) :: namen(0:)
    integer(i64), intent(in) :: oldnamen(0:)
    real(f64), intent(in) :: volume(0:)
    real(f64), intent(in), target :: centergn(0:,0:,0:)
    real(f64), intent(in), target :: halocentergn(0:,0:,0:)
    integer(i64), intent(in) :: periodicnid(0:,0:)
    integer(i64), intent(in) :: cellnid(0:,0:)
    real(f64), intent(in), target :: centerc(0:,0:)
    real(f64), intent(in), target :: centerh(0:,0:)
    integer(i64), intent(in) :: halonid(0:,0:)
    real(f64), intent(in) :: airDiamond(0:)
    real(f64), intent(in) :: lambda_x(0:)
    real(f64), intent(in) :: lambda_y(0:)
    real(f64), intent(in) :: lambda_z(0:)
    integer(i64), intent(in) :: number(0:)
    real(f64), intent(in) :: R_x(0:)
    real(f64), intent(in) :: R_y(0:)
    real(f64), intent(in) :: R_z(0:)
    real(f64), intent(in) :: param1(0:)
    real(f64), intent(in) :: param2(0:)
    real(f64), intent(in) :: param3(0:)
    real(f64), intent(in) :: shift(0:,0:)
    integer(i64), intent(in) :: loctoglob(0:)
    integer(i64), intent(in) :: BCdirichlet(0:)
    real(f64), intent(inout) :: a_loc(0:)
    integer(i32), intent(inout) :: irn_loc(0:)
    integer(i32), intent(inout) :: jcn_loc(0:)
    integer(i64), intent(in) :: matrixinnerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: dirichletfaces(0:)
    real(f64), allocatable :: parameters(:)
    integer(i64), allocatable :: nodes(:)
    integer(i64) :: cmpt
    integer(i64) :: i
    integer(i64) :: c_left
    integer(i64) :: c_leftglob
    integer(i64) :: c_right
    integer(i64) :: c_rightglob
    real(f64) :: value
    integer(i64) :: cmptparam
    integer(i64) :: nod
    integer(i64) :: j
    real(f64), pointer :: center(:)
    real(f64) :: xdiff
    real(f64) :: ydiff
    real(f64) :: zdiff
    real(f64) :: alpha
    integer(i64) :: index
    integer(i64) :: Dummy_0042
    integer(i64) :: Dummy_0043
    integer(i64) :: Dummy_0044
    integer(i64) :: Dummy_0045
    integer(i64) :: Dummy_0046

    allocate(parameters(0:3_i64))
    parameters = 0.0_f64
    allocate(nodes(0:3_i64))
    nodes = 0_i64
    cmpt = 0_i64
    do Dummy_0042 = 0_i64, size(matrixinnerfaces, kind=i64) - 1_i64, &
          1_i64
      i = matrixinnerfaces(Dummy_0042)
      c_left = cellfid(0_i64, i)
      c_leftglob = loctoglob(c_left)
      nodes(0_i64:2_i64) = nodeidf(:, i)
      nodes(3_i64) = nodeidf(2_i64, i)
      parameters(0_i64) = param1(i)
      parameters(1_i64) = param2(i)
      parameters(2_i64) = (-1.0_f64) * param1(i)
      parameters(3_i64) = (-1.0_f64) * param2(i)
      c_right = cellfid(1_i64, i)
      c_rightglob = loctoglob(c_right)
      irn_loc(cmpt) = c_leftglob
      jcn_loc(cmpt) = c_leftglob
      value = (-1_i64) * param3(i) / volume(c_left)
      a_loc(cmpt) = value
      cmpt = cmpt + 1_i64
      cmptparam = 0_i64
      do Dummy_0043 = 0_i64, size(nodes, kind=i64) - 1_i64, 1_i64
        nod = nodes(Dummy_0043)
        if (search_element(BCdirichlet, oldnamen(nod)) == 0_i64) then
          do j = 0_i64, cellnid(size(cellnid, 1_i64, i64) - 1_i64, nod) &
                - 1_i64, 1_i64
            center(0:) => centerc(:, cellnid(j, nod))
            xdiff = center(0_i64) - vertexn(0_i64, nod)
            ydiff = center(1_i64) - vertexn(1_i64, nod)
            zdiff = center(2_i64) - vertexn(2_i64, nod)
            alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                  ydiff + lambda_z(nod) * zdiff) / (number(nod) + &
                  lambda_x(nod) * R_x(nod) + lambda_y(nod) * R_y(nod) + &
                  lambda_z(nod) * R_z(nod))
            value = alpha / volume(c_left) * parameters(cmptparam)
            irn_loc(cmpt) = c_leftglob
            jcn_loc(cmpt) = loctoglob(cellnid(j, nod))
            a_loc(cmpt) = value
            cmpt = cmpt + 1_i64
            !right cell-----------------------------------
            value = (-1.0_f64) * alpha / volume(c_right) * parameters( &
                  cmptparam)
            irn_loc(cmpt) = c_rightglob
            jcn_loc(cmpt) = loctoglob(cellnid(j, nod))
            a_loc(cmpt) = value
            cmpt = cmpt + 1_i64
          end do
          do j = 0_i64, halonid(size(halonid, 1_i64, i64) - 1_i64, nod) &
                - 1_i64, 1_i64
            center(0:) => centerh(:, halonid(j, nod))
            xdiff = center(0_i64) - vertexn(0_i64, nod)
            ydiff = center(1_i64) - vertexn(1_i64, nod)
            zdiff = center(2_i64) - vertexn(2_i64, nod)
            alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                  ydiff + lambda_z(nod) * zdiff) / (number(nod) + &
                  lambda_x(nod) * R_x(nod) + lambda_y(nod) * R_y(nod) + &
                  lambda_z(nod) * R_z(nod))
            value = alpha / volume(c_left) * parameters(cmptparam)
            irn_loc(cmpt) = c_leftglob
            jcn_loc(cmpt) = haloext(0_i64, halonid(j, nod))
            a_loc(cmpt) = value
            cmpt = cmpt + 1_i64
            !right cell-----------------------------------
            value = (-1.0_f64) * alpha / volume(c_right) * parameters( &
                  cmptparam)
            irn_loc(cmpt) = c_rightglob
            jcn_loc(cmpt) = haloext(0_i64, halonid(j, nod))
            a_loc(cmpt) = value
            cmpt = cmpt + 1_i64
          end do
          do j = 0_i64, periodicnid(size(periodicnid, 1_i64, i64) - &
                1_i64, nod) - 1_i64, 1_i64
            if (vertexn(3_i64, nod) == 11_i64 .or. vertexn(3_i64, nod) &
                  == 22_i64) then
              center(0_i64) = centerc(0_i64, periodicnid(j, nod)) + &
                    shift(0_i64, periodicnid(j, nod))
              center(1_i64) = centerc(1_i64, periodicnid(j, nod))
              center(2_i64) = centerc(2_i64, periodicnid(j, nod))
            end if
            if (vertexn(3_i64, nod) == 33_i64 .or. vertexn(3_i64, nod) &
                  == 44_i64) then
              center(0_i64) = centerc(0_i64, periodicnid(j, nod))
              center(1_i64) = centerc(1_i64, periodicnid(j, nod)) + &
                    shift(1_i64, periodicnid(j, nod))
              center(2_i64) = centerc(2_i64, periodicnid(j, nod))
            end if
            if (vertexn(3_i64, nod) == 55_i64 .or. vertexn(3_i64, nod) &
                  == 66_i64) then
              center(0_i64) = centerc(0_i64, periodicnid(j, nod))
              center(1_i64) = centerc(1_i64, periodicnid(j, nod))
              center(2_i64) = centerc(2_i64, periodicnid(j, nod)) + &
                    shift(2_i64, periodicnid(j, nod))
            end if
            xdiff = center(0_i64) - vertexn(0_i64, nod)
            ydiff = center(1_i64) - vertexn(1_i64, nod)
            zdiff = center(2_i64) - vertexn(2_i64, nod)
            alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                  ydiff + lambda_z(nod) * zdiff) / (number(nod) + &
                  lambda_x(nod) * R_x(nod) + lambda_y(nod) * R_y(nod) + &
                  lambda_z(nod) * R_z(nod))
            value = alpha / volume(c_left) * parameters(cmptparam)
            irn_loc(cmpt) = c_leftglob
            jcn_loc(cmpt) = loctoglob(periodicnid(j, nod))
            a_loc(cmpt) = value
            cmpt = cmpt + 1_i64
            !right cell-----------------------------------
            value = (-1.0_f64) * alpha / volume(c_right) * parameters( &
                  cmptparam)
            irn_loc(cmpt) = c_rightglob
            jcn_loc(cmpt) = loctoglob(periodicnid(j, nod))
            a_loc(cmpt) = value
            cmpt = cmpt + 1_i64
          end do
          do j = 0_i64, size(centergn(:, :, nod),2,i64) - 1_i64, 1_i64
            if (centergn(size(centergn, 1_i64, i64) - 1_i64, j, nod) /= &
                  -1_i64) then
              center(0:) => centergn(0_i64:2_i64, j, nod)
              xdiff = center(0_i64) - vertexn(0_i64, nod)
              ydiff = center(1_i64) - vertexn(1_i64, nod)
              zdiff = center(2_i64) - vertexn(2_i64, nod)
              alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                    ydiff + lambda_z(nod) * zdiff) / (number(nod) + &
                    lambda_x(nod) * R_x(nod) + lambda_y(nod) * R_y(nod &
                    ) + lambda_z(nod) * R_z(nod))
              value = alpha / volume(c_left) * parameters(cmptparam)
              index = Int(centergn(3_i64, j, nod), i64)
              irn_loc(cmpt) = c_leftglob
              jcn_loc(cmpt) = loctoglob(index)
              a_loc(cmpt) = value
              cmpt = cmpt + 1_i64
              !right cell-----------------------------------
              value = (-1.0_f64) * alpha / volume(c_right) * parameters( &
                    cmptparam)
              irn_loc(cmpt) = c_rightglob
              jcn_loc(cmpt) = loctoglob(index)
              a_loc(cmpt) = value
              cmpt = cmpt + 1_i64
            end if
          end do
          do j = 0_i64, size(halocentergn(:, :, nod),2,i64) - 1_i64, &
                1_i64
            if (halocentergn(size(halocentergn, 1_i64, i64) - 1_i64, j, &
                  nod) /= -1_i64) then
              center(0:) => halocentergn(0_i64:2_i64, j, nod)
              xdiff = center(0_i64) - vertexn(0_i64, nod)
              ydiff = center(1_i64) - vertexn(1_i64, nod)
              zdiff = center(2_i64) - vertexn(2_i64, nod)
              alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                    ydiff + lambda_z(nod) * zdiff) / (number(nod) + &
                    lambda_x(nod) * R_x(nod) + lambda_y(nod) * R_y(nod &
                    ) + lambda_z(nod) * R_z(nod))
              value = alpha / volume(c_left) * parameters(cmptparam)
              index = Int(halocentergn(3_i64, j, nod), i64)
              irn_loc(cmpt) = c_leftglob
              jcn_loc(cmpt) = haloext(0_i64, index)
              a_loc(cmpt) = value
              cmpt = cmpt + 1_i64
              !right cell-----------------------------------
              value = (-1.0_f64) * alpha / volume(c_right) * parameters( &
                    cmptparam)
              irn_loc(cmpt) = c_rightglob
              jcn_loc(cmpt) = haloext(0_i64, index)
              a_loc(cmpt) = value
              cmpt = cmpt + 1_i64
            end if
          end do
        end if
        cmptparam = cmptparam + 1_i64
      end do
      irn_loc(cmpt) = c_leftglob
      jcn_loc(cmpt) = c_rightglob
      value = param3(i) / volume(c_left)
      a_loc(cmpt) = value
      cmpt = cmpt + 1_i64
      !right cell------------------------------------------------------
      irn_loc(cmpt) = c_rightglob
      jcn_loc(cmpt) = c_leftglob
      value = param3(i) / volume(c_right)
      a_loc(cmpt) = value
      cmpt = cmpt + 1_i64
      irn_loc(cmpt) = c_rightglob
      jcn_loc(cmpt) = c_rightglob
      value = (-1.0_f64) * param3(i) / volume(c_right)
      a_loc(cmpt) = value
      cmpt = cmpt + 1_i64
    end do
    do Dummy_0044 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0044)
      c_left = cellfid(0_i64, i)
      c_leftglob = loctoglob(c_left)
      nodes(0_i64:2_i64) = nodeidf(:, i)
      nodes(3_i64) = nodeidf(2_i64, i)
      parameters(0_i64) = param1(i)
      parameters(1_i64) = param2(i)
      parameters(2_i64) = (-1.0_f64) * param1(i)
      parameters(3_i64) = (-1.0_f64) * param2(i)
      c_rightglob = haloext(0_i64, halofid(i))
      c_right = halofid(i)
      cmptparam = 0_i64
      do Dummy_0045 = 0_i64, size(nodes, kind=i64) - 1_i64, 1_i64
        nod = nodes(Dummy_0045)
        if (search_element(BCdirichlet, oldnamen(nod)) == 0_i64) then
          do j = 0_i64, cellnid(size(cellnid, 1_i64, i64) - 1_i64, nod) &
                - 1_i64, 1_i64
            center(0:) => centerc(:, cellnid(j, nod))
            xdiff = center(0_i64) - vertexn(0_i64, nod)
            ydiff = center(1_i64) - vertexn(1_i64, nod)
            zdiff = center(2_i64) - vertexn(2_i64, nod)
            alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                  ydiff + lambda_z(nod) * zdiff) / (number(nod) + &
                  lambda_x(nod) * R_x(nod) + lambda_y(nod) * R_y(nod) + &
                  lambda_z(nod) * R_z(nod))
            value = alpha / volume(c_left) * parameters(cmptparam)
            irn_loc(cmpt) = c_leftglob
            jcn_loc(cmpt) = loctoglob(cellnid(j, nod))
            a_loc(cmpt) = value
            cmpt = cmpt + 1_i64
          end do
          do j = 0_i64, halonid(size(halonid, 1_i64, i64) - 1_i64, nod) &
                - 1_i64, 1_i64
            center(0:) => centerh(:, halonid(j, nod))
            xdiff = center(0_i64) - vertexn(0_i64, nod)
            ydiff = center(1_i64) - vertexn(1_i64, nod)
            zdiff = center(2_i64) - vertexn(2_i64, nod)
            alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                  ydiff + lambda_z(nod) * zdiff) / (number(nod) + &
                  lambda_x(nod) * R_x(nod) + lambda_y(nod) * R_y(nod) + &
                  lambda_z(nod) * R_z(nod))
            value = alpha / volume(c_left) * parameters(cmptparam)
            irn_loc(cmpt) = c_leftglob
            jcn_loc(cmpt) = haloext(0_i64, halonid(j, nod))
            a_loc(cmpt) = value
            cmpt = cmpt + 1_i64
          end do
          do j = 0_i64, size(centergn(:, :, nod),2,i64) - 1_i64, 1_i64
            if (centergn(size(centergn, 1_i64, i64) - 1_i64, j, nod) /= &
                  -1_i64) then
              center(0:) => centergn(0_i64:2_i64, j, nod)
              xdiff = center(0_i64) - vertexn(0_i64, nod)
              ydiff = center(1_i64) - vertexn(1_i64, nod)
              zdiff = center(2_i64) - vertexn(2_i64, nod)
              alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                    ydiff + lambda_z(nod) * zdiff) / (number(nod) + &
                    lambda_x(nod) * R_x(nod) + lambda_y(nod) * R_y(nod &
                    ) + lambda_z(nod) * R_z(nod))
              value = alpha / volume(c_left) * parameters(cmptparam)
              index = Int(centergn(3_i64, j, nod), i64)
              irn_loc(cmpt) = c_leftglob
              jcn_loc(cmpt) = loctoglob(index)
              a_loc(cmpt) = value
              cmpt = cmpt + 1_i64
            end if
          end do
          do j = 0_i64, size(halocentergn(:, :, nod),2,i64) - 1_i64, &
                1_i64
            if (halocentergn(size(halocentergn, 1_i64, i64) - 1_i64, j, &
                  nod) /= -1_i64) then
              center(0:) => halocentergn(0_i64:2_i64, j, nod)
              xdiff = center(0_i64) - vertexn(0_i64, nod)
              ydiff = center(1_i64) - vertexn(1_i64, nod)
              zdiff = center(2_i64) - vertexn(2_i64, nod)
              alpha = (1.0_f64 + lambda_x(nod) * xdiff + lambda_y(nod) * &
                    ydiff + lambda_z(nod) * zdiff) / (number(nod) + &
                    lambda_x(nod) * R_x(nod) + lambda_y(nod) * R_y(nod &
                    ) + lambda_z(nod) * R_z(nod))
              value = alpha / volume(c_left) * parameters(cmptparam)
              index = Int(halocentergn(3_i64, j, nod), i64)
              irn_loc(cmpt) = c_leftglob
              jcn_loc(cmpt) = haloext(0_i64, index)
              a_loc(cmpt) = value
              cmpt = cmpt + 1_i64
            end if
          end do
        end if
        cmptparam = cmptparam + 1_i64
      end do
      irn_loc(cmpt) = c_leftglob
      jcn_loc(cmpt) = c_leftglob
      value = (-1_i64) * param3(i) / volume(c_left)
      a_loc(cmpt) = value
      cmpt = cmpt + 1_i64
      irn_loc(cmpt) = c_leftglob
      jcn_loc(cmpt) = c_rightglob
      value = param3(i) / volume(c_left)
      a_loc(cmpt) = value
      cmpt = cmpt + 1_i64
    end do
    do Dummy_0046 = 0_i64, size(dirichletfaces, kind=i64) - 1_i64, 1_i64
      i = dirichletfaces(Dummy_0046)
      c_left = cellfid(0_i64, i)
      c_leftglob = loctoglob(c_left)
      parameters(0_i64) = param1(i)
      parameters(1_i64) = param2(i)
      parameters(2_i64) = (-1.0_f64) * param1(i)
      parameters(3_i64) = (-1.0_f64) * param2(i)
      irn_loc(cmpt) = c_leftglob
      jcn_loc(cmpt) = c_leftglob
      value = (-1_i64) * param3(i) / volume(c_left)
      a_loc(cmpt) = value
      cmpt = cmpt + 1_i64
      irn_loc(cmpt) = c_leftglob
      jcn_loc(cmpt) = c_leftglob
      value = (-1.0_f64) * param3(i) / volume(c_left)
      a_loc(cmpt) = value
      cmpt = cmpt + 1_i64
    end do
    if (allocated(parameters)) then
      deallocate(parameters)
    end if
    if (allocated(nodes)) then
      deallocate(nodes)
    end if

  end subroutine get_triplet_3d
  !........................................

  !........................................
  subroutine get_rhs_loc_2d(cellfid, nodeidf, oldname, volume, centergn, &
        param1, param2, param3, param4, Pbordnode, Pbordface, rhs_loc, &
        BCdirichlet, centergf, matrixinnerfaces, halofaces, &
        dirichletfaces)

    implicit none

    integer(i64), intent(in) :: cellfid(0:,0:)
    integer(i64), intent(in) :: nodeidf(0:,0:)
    integer(i64), intent(in) :: oldname(0:)
    real(f64), intent(in) :: volume(0:)
    real(f64), intent(in) :: centergn(0:,0:,0:)
    real(f64), intent(in) :: param1(0:)
    real(f64), intent(in) :: param2(0:)
    real(f64), intent(in) :: param3(0:)
    real(f64), intent(in) :: param4(0:)
    real(f64), intent(in) :: Pbordnode(0:)
    real(f64), intent(in) :: Pbordface(0:)
    real(f64), intent(inout) :: rhs_loc(0:)
    integer(i64), intent(in) :: BCdirichlet(0:)
    real(f64), intent(in) :: centergf(0:,0:)
    integer(i64), intent(in) :: matrixinnerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: dirichletfaces(0:)
    integer(i64) :: i
    integer(i64) :: c_right
    integer(i64) :: c_left
    integer(i64) :: i_1
    integer(i64) :: i_2
    real(f64) :: V
    real(f64) :: value_left
    real(f64) :: value_right
    real(f64) :: V_K
    real(f64) :: value
    integer(i64) :: Dummy_0047
    integer(i64) :: Dummy_0048
    integer(i64) :: Dummy_0049

    do Dummy_0047 = 0_i64, size(matrixinnerfaces, kind=i64) - 1_i64, &
          1_i64
      i = matrixinnerfaces(Dummy_0047)
      c_right = cellfid(1_i64, i)
      c_left = cellfid(0_i64, i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      if (search_element(BCdirichlet, oldname(i_1)) == 1_i64) then
        V = Pbordnode(i_1)
        value_left = (-1.0_f64) * V * param4(i) / volume(c_left)
        rhs_loc(c_left) = rhs_loc(c_left) + value_left
        value_right = V * param4(i) / volume(c_right)
        rhs_loc(c_right) = rhs_loc(c_right) + value_right
      end if
      if (search_element(BCdirichlet, oldname(i_2)) == 1_i64) then
        V = Pbordnode(i_2)
        value_left = (-1.0_f64) * V * param2(i) / volume(c_left)
        rhs_loc(c_left) = rhs_loc(c_left) + value_left
        value_right = V * param2(i) / volume(c_right)
        rhs_loc(c_right) = rhs_loc(c_right) + value_right
      end if
    end do
    do Dummy_0048 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0048)
      c_left = cellfid(0_i64, i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      if (search_element(BCdirichlet, oldname(i_1)) == 1_i64) then
        V = Pbordnode(i_1)
        value_left = (-1.0_f64) * V * param4(i) / volume(c_left)
        rhs_loc(c_left) = rhs_loc(c_left) + value_left
      end if
      if (search_element(BCdirichlet, oldname(i_2)) == 1_i64) then
        V = Pbordnode(i_2)
        value_left = (-1.0_f64) * V * param2(i) / volume(c_left)
        rhs_loc(c_left) = rhs_loc(c_left) + value_left
      end if
    end do
    !TODO verify
    do Dummy_0049 = 0_i64, size(dirichletfaces, kind=i64) - 1_i64, 1_i64
      i = dirichletfaces(Dummy_0049)
      c_left = cellfid(0_i64, i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      if (centergn(2_i64, 0_i64, i_1) /= -1_i64) then
        V = Pbordnode(i_1)
        value_left = (-1.0_f64) * V * param4(i) / volume(c_left)
        rhs_loc(c_left) = rhs_loc(c_left) + value_left
      end if
      if (centergn(2_i64, 0_i64, i_2) /= -1_i64) then
        V = Pbordnode(i_2)
        value_left = (-1.0_f64) * V * param2(i) / volume(c_left)
        rhs_loc(c_left) = rhs_loc(c_left) + value_left
      end if
      V_K = Pbordface(i)
      value = (-2.0_f64) * param3(i) / volume(c_left) * V_K
      rhs_loc(c_left) = rhs_loc(c_left) + value
    end do

  end subroutine get_rhs_loc_2d
  !........................................

  !........................................
  subroutine get_rhs_glob_2d(cellfid, nodeidf, oldname, volume, centergn &
        , loctoglob, param1, param2, param3, param4, Pbordnode, &
        Pbordface, rhs, BCdirichlet, centergf, matrixinnerfaces, &
        halofaces, dirichletfaces)

    implicit none

    integer(i64), intent(in) :: cellfid(0:,0:)
    integer(i64), intent(in) :: nodeidf(0:,0:)
    integer(i64), intent(in) :: oldname(0:)
    real(f64), intent(in) :: volume(0:)
    real(f64), intent(in) :: centergn(0:,0:,0:)
    integer(i64), intent(in) :: loctoglob(0:)
    real(f64), intent(in) :: param1(0:)
    real(f64), intent(in) :: param2(0:)
    real(f64), intent(in) :: param3(0:)
    real(f64), intent(in) :: param4(0:)
    real(f64), intent(in) :: Pbordnode(0:)
    real(f64), intent(in) :: Pbordface(0:)
    real(f64), intent(inout) :: rhs(0:)
    integer(i64), intent(in) :: BCdirichlet(0:)
    real(f64), intent(in) :: centergf(0:,0:)
    integer(i64), intent(in) :: matrixinnerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: dirichletfaces(0:)
    integer(i64) :: i
    integer(i64) :: c_left
    integer(i64) :: c_leftglob
    integer(i64) :: i_1
    integer(i64) :: i_2
    integer(i64) :: c_right
    integer(i64) :: c_rightglob
    real(f64) :: V
    real(f64) :: value_left
    real(f64) :: value_right
    real(f64) :: V_K
    real(f64) :: value
    integer(i64) :: Dummy_0050
    integer(i64) :: Dummy_0051
    integer(i64) :: Dummy_0052

    do Dummy_0050 = 0_i64, size(matrixinnerfaces, kind=i64) - 1_i64, &
          1_i64
      i = matrixinnerfaces(Dummy_0050)
      c_left = cellfid(0_i64, i)
      c_leftglob = loctoglob(c_left)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      c_right = cellfid(1_i64, i)
      c_rightglob = loctoglob(c_right)
      if (search_element(BCdirichlet, oldname(i_1)) == 1_i64) then
        V = Pbordnode(i_1)
        value_left = (-1.0_f64) * V * param4(i) / volume(c_left)
        rhs(c_leftglob) = rhs(c_leftglob) + value_left
        value_right = V * param4(i) / volume(c_right)
        rhs(c_rightglob) = rhs(c_rightglob) + value_right
      end if
      if (search_element(BCdirichlet, oldname(i_2)) == 1_i64) then
        V = Pbordnode(i_2)
        value_left = (-1.0_f64) * V * param2(i) / volume(c_left)
        rhs(c_leftglob) = rhs(c_leftglob) + value_left
        value_right = V * param2(i) / volume(c_right)
        rhs(c_rightglob) = rhs(c_rightglob) + value_right
      end if
    end do
    do Dummy_0051 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0051)
      c_left = cellfid(0_i64, i)
      c_leftglob = loctoglob(c_left)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      if (search_element(BCdirichlet, oldname(i_1)) == 1_i64) then
        V = Pbordnode(i_1)
        value_left = (-1.0_f64) * V * param4(i) / volume(c_left)
        rhs(c_leftglob) = rhs(c_leftglob) + value_left
      end if
      if (search_element(BCdirichlet, oldname(i_2)) == 1_i64) then
        V = Pbordnode(i_2)
        value_left = (-1.0_f64) * V * param2(i) / volume(c_left)
        rhs(c_leftglob) = rhs(c_leftglob) + value_left
      end if
    end do
    do Dummy_0052 = 0_i64, size(dirichletfaces, kind=i64) - 1_i64, 1_i64
      i = dirichletfaces(Dummy_0052)
      c_left = cellfid(0_i64, i)
      c_leftglob = loctoglob(c_left)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      if (centergn(2_i64, 0_i64, i_1) /= -1_i64) then
        V = Pbordnode(i_1)
        value_left = (-1.0_f64) * V * param4(i) / volume(c_left)
        rhs(c_leftglob) = rhs(c_leftglob) + value_left
      end if
      if (centergn(2_i64, 0_i64, i_2) /= -1_i64) then
        V = Pbordnode(i_2)
        value_left = (-1.0_f64) * V * param2(i) / volume(c_left)
        rhs(c_leftglob) = rhs(c_leftglob) + value_left
      end if
      V_K = Pbordface(i)
      value = (-2.0_f64) * param3(i) / volume(c_left) * V_K
      rhs(c_leftglob) = rhs(c_leftglob) + value
    end do

  end subroutine get_rhs_glob_2d
  !........................................

  !........................................
  subroutine get_rhs_loc_3d(cellfid, nodeidf, oldname, volume, centergn, &
        param1, param2, param3, Pbordnode, Pbordface, rhs_loc, &
        BCdirichlet, matrixinnerfaces, halofaces, dirichletfaces)

    implicit none

    integer(i64), intent(in) :: cellfid(0:,0:)
    integer(i64), intent(in) :: nodeidf(0:,0:)
    integer(i64), intent(in) :: oldname(0:)
    real(f64), intent(in) :: volume(0:)
    real(f64), intent(in) :: centergn(0:,0:,0:)
    real(f64), intent(in) :: param1(0:)
    real(f64), intent(in) :: param2(0:)
    real(f64), intent(in) :: param3(0:)
    real(f64), intent(in) :: Pbordnode(0:)
    real(f64), intent(in) :: Pbordface(0:)
    real(f64), intent(inout) :: rhs_loc(0:)
    integer(i64), intent(in) :: BCdirichlet(0:)
    integer(i64), intent(in) :: matrixinnerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: dirichletfaces(0:)
    real(f64), allocatable :: parameters(:)
    integer(i64), allocatable :: nodes(:)
    integer(i64) :: i
    integer(i64) :: c_left
    integer(i64) :: c_right
    integer(i64) :: cmpt
    integer(i64) :: nod
    real(f64) :: V
    real(f64) :: value_left
    real(f64) :: value_right
    real(f64) :: V_K
    real(f64) :: value
    integer(i64) :: Dummy_0053
    integer(i64) :: Dummy_0054
    integer(i64) :: Dummy_0055
    integer(i64) :: Dummy_0056
    integer(i64) :: Dummy_0057
    integer(i64) :: Dummy_0058

    allocate(parameters(0:3_i64))
    parameters = 0.0_f64
    allocate(nodes(0:3_i64))
    nodes = 0_i64
    do Dummy_0053 = 0_i64, size(matrixinnerfaces, kind=i64) - 1_i64, &
          1_i64
      i = matrixinnerfaces(Dummy_0053)
      c_left = cellfid(0_i64, i)
      c_right = cellfid(1_i64, i)
      nodes(0_i64:2_i64) = nodeidf(0_i64:2_i64, i)
      nodes(3_i64) = nodeidf(2_i64, i)
      parameters(0_i64) = param1(i)
      parameters(1_i64) = param2(i)
      parameters(2_i64) = (-1.0_f64) * param1(i)
      parameters(3_i64) = (-1.0_f64) * param2(i)
      cmpt = 0_i64
      do Dummy_0054 = 0_i64, size(nodes, kind=i64) - 1_i64, 1_i64
        nod = nodes(Dummy_0054)
        if (search_element(BCdirichlet, oldname(nod)) == 1_i64) then
          V = Pbordnode(nod)
          value_left = (-1.0_f64) * V * parameters(cmpt) / volume(c_left &
                )
          rhs_loc(c_left) = rhs_loc(c_left) + value_left
          value_right = V * parameters(cmpt) / volume(c_right)
          rhs_loc(c_right) = rhs_loc(c_right) + value_right
        end if
        cmpt = cmpt + 1_i64
      end do
    end do
    do Dummy_0055 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0055)
      c_left = cellfid(0_i64, i)
      nodes(0_i64:2_i64) = nodeidf(0_i64:2_i64, i)
      nodes(3_i64) = nodeidf(2_i64, i)
      parameters(0_i64) = param1(i)
      parameters(1_i64) = param2(i)
      parameters(2_i64) = (-1.0_f64) * param1(i)
      parameters(3_i64) = (-1.0_f64) * param2(i)
      cmpt = 0_i64
      do Dummy_0056 = 0_i64, size(nodes, kind=i64) - 1_i64, 1_i64
        nod = nodes(Dummy_0056)
        if (search_element(BCdirichlet, oldname(nod)) == 1_i64) then
          V = Pbordnode(nod)
          value_left = (-1.0_f64) * V * parameters(cmpt) / volume(c_left &
                )
          rhs_loc(c_left) = rhs_loc(c_left) + value_left
        end if
        cmpt = cmpt + 1_i64
      end do
    end do
    do Dummy_0057 = 0_i64, size(dirichletfaces, kind=i64) - 1_i64, 1_i64
      i = dirichletfaces(Dummy_0057)
      c_left = cellfid(0_i64, i)
      nodes(0_i64:2_i64) = nodeidf(0_i64:2_i64, i)
      nodes(3_i64) = nodeidf(2_i64, i)
      parameters(0_i64) = param1(i)
      parameters(1_i64) = param2(i)
      parameters(2_i64) = (-1.0_f64) * param1(i)
      parameters(3_i64) = (-1.0_f64) * param2(i)
      cmpt = 0_i64
      do Dummy_0058 = 0_i64, size(nodes, kind=i64) - 1_i64, 1_i64
        nod = nodes(Dummy_0058)
        if (centergn(3_i64, 0_i64, nod) /= -1_i64) then
          V = Pbordnode(nod)
          value_left = (-1.0_f64) * V * parameters(cmpt) / volume(c_left &
                )
          rhs_loc(c_left) = rhs_loc(c_left) + value_left
        end if
        cmpt = cmpt + 1_i64
      end do
      V_K = Pbordface(i)
      value = (-2.0_f64) * param3(i) / volume(c_left) * V_K
      rhs_loc(c_left) = rhs_loc(c_left) + value
    end do
    if (allocated(parameters)) then
      deallocate(parameters)
    end if
    if (allocated(nodes)) then
      deallocate(nodes)
    end if

  end subroutine get_rhs_loc_3d
  !........................................

  !........................................
  subroutine get_rhs_glob_3d(cellfid, nodeidf, oldname, volume, centergn &
        , loctoglob, param1, param2, param3, Pbordnode, Pbordface, rhs, &
        BCdirichlet, matrixinnerfaces, halofaces, dirichletfaces)

    implicit none

    integer(i64), intent(in) :: cellfid(0:,0:)
    integer(i64), intent(in) :: nodeidf(0:,0:)
    integer(i64), intent(in) :: oldname(0:)
    real(f64), intent(in) :: volume(0:)
    real(f64), intent(in) :: centergn(0:,0:,0:)
    integer(i64), intent(in) :: loctoglob(0:)
    real(f64), intent(in) :: param1(0:)
    real(f64), intent(in) :: param2(0:)
    real(f64), intent(in) :: param3(0:)
    real(f64), intent(in) :: Pbordnode(0:)
    real(f64), intent(in) :: Pbordface(0:)
    real(f64), intent(inout) :: rhs(0:)
    integer(i64), intent(in) :: BCdirichlet(0:)
    integer(i64), intent(in) :: matrixinnerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: dirichletfaces(0:)
    real(f64), allocatable :: parameters(:)
    integer(i64), allocatable :: nodes(:)
    integer(i64) :: i
    integer(i64) :: c_left
    integer(i64) :: c_leftglob
    integer(i64) :: c_right
    integer(i64) :: c_rightglob
    integer(i64) :: cmpt
    integer(i64) :: nod
    real(f64) :: V
    real(f64) :: value_left
    real(f64) :: value_right
    real(f64) :: V_K
    real(f64) :: value
    integer(i64) :: Dummy_0059
    integer(i64) :: Dummy_0060
    integer(i64) :: Dummy_0061
    integer(i64) :: Dummy_0062
    integer(i64) :: Dummy_0063
    integer(i64) :: Dummy_0064

    allocate(parameters(0:3_i64))
    parameters = 0.0_f64
    allocate(nodes(0:3_i64))
    nodes = 0_i64
    do Dummy_0059 = 0_i64, size(matrixinnerfaces, kind=i64) - 1_i64, &
          1_i64
      i = matrixinnerfaces(Dummy_0059)
      c_left = cellfid(0_i64, i)
      c_leftglob = loctoglob(c_left)
      nodes(0_i64:2_i64) = nodeidf(0_i64:2_i64, i)
      nodes(3_i64) = nodeidf(2_i64, i)
      parameters(0_i64) = param1(i)
      parameters(1_i64) = param2(i)
      parameters(2_i64) = (-1.0_f64) * param1(i)
      parameters(3_i64) = (-1.0_f64) * param2(i)
      c_right = cellfid(1_i64, i)
      c_rightglob = loctoglob(c_right)
      cmpt = 0_i64
      do Dummy_0060 = 0_i64, size(nodes, kind=i64) - 1_i64, 1_i64
        nod = nodes(Dummy_0060)
        if (search_element(BCdirichlet, oldname(nod)) == 1_i64) then
          V = Pbordnode(nod)
          value_left = (-1.0_f64) * V * parameters(cmpt) / volume(c_left &
                )
          rhs(c_leftglob) = rhs(c_leftglob) + value_left
          value_right = V * parameters(cmpt) / volume(c_right)
          rhs(c_rightglob) = rhs(c_rightglob) + value_right
        end if
        cmpt = cmpt + 1_i64
      end do
    end do
    do Dummy_0061 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0061)
      c_left = cellfid(0_i64, i)
      c_leftglob = loctoglob(c_left)
      nodes(0_i64:2_i64) = nodeidf(0_i64:2_i64, i)
      nodes(3_i64) = nodeidf(2_i64, i)
      parameters(0_i64) = param1(i)
      parameters(1_i64) = param2(i)
      parameters(2_i64) = (-1.0_f64) * param1(i)
      parameters(3_i64) = (-1.0_f64) * param2(i)
      cmpt = 0_i64
      do Dummy_0062 = 0_i64, size(nodes, kind=i64) - 1_i64, 1_i64
        nod = nodes(Dummy_0062)
        if (search_element(BCdirichlet, oldname(nod)) == 1_i64) then
          V = Pbordnode(nod)
          value_left = (-1.0_f64) * V * parameters(cmpt) / volume(c_left &
                )
          rhs(c_leftglob) = rhs(c_leftglob) + value_left
        end if
        cmpt = cmpt + 1_i64
      end do
    end do
    do Dummy_0063 = 0_i64, size(dirichletfaces, kind=i64) - 1_i64, 1_i64
      i = dirichletfaces(Dummy_0063)
      c_left = cellfid(0_i64, i)
      c_leftglob = loctoglob(c_left)
      nodes(0_i64:2_i64) = nodeidf(0_i64:2_i64, i)
      nodes(3_i64) = nodeidf(2_i64, i)
      parameters(0_i64) = param1(i)
      parameters(1_i64) = param2(i)
      parameters(2_i64) = (-1.0_f64) * param1(i)
      parameters(3_i64) = (-1.0_f64) * param2(i)
      cmpt = 0_i64
      do Dummy_0064 = 0_i64, size(nodes, kind=i64) - 1_i64, 1_i64
        nod = nodes(Dummy_0064)
        if (centergn(3_i64, 0_i64, nod) /= -1_i64) then
          V = Pbordnode(nod)
          value_left = (-1.0_f64) * V * parameters(cmpt) / volume(c_left &
                )
          rhs(c_leftglob) = rhs(c_leftglob) + value_left
        end if
        cmpt = cmpt + 1_i64
      end do
      V_K = Pbordface(i)
      value = (-2.0_f64) * param3(i) / volume(c_left) * V_K
      rhs(c_leftglob) = rhs(c_leftglob) + value
    end do
    if (allocated(parameters)) then
      deallocate(parameters)
    end if
    if (allocated(nodes)) then
      deallocate(nodes)
    end if

  end subroutine get_rhs_glob_3d
  !........................................

  !........................................
  subroutine compute_P_gradient_2d(P_c, P_ghost, P_halo, P_node, cellidf &
        , nodeidf, centergf, namef, halofid, centerc, centerh, oldname, &
        airDiamond, f_1, f_2, f_3, f_4, normalf, shift, Pbordnode, &
        Pbordface, Px_face, Py_face, Pz_face, BCdirichlet, innerfaces, &
        halofaces, neumannfaces, dirichletfaces, periodicfaces)

    implicit none

    real(f64), intent(in) :: P_c(0:)
    real(f64), intent(in) :: P_ghost(0:)
    real(f64), intent(in) :: P_halo(0:)
    real(f64), intent(in) :: P_node(0:)
    integer(i64), intent(in) :: cellidf(0:,0:)
    integer(i64), intent(in) :: nodeidf(0:,0:)
    real(f64), intent(in) :: centergf(0:,0:)
    integer(i64), intent(in) :: namef(0:)
    integer(i64), intent(in) :: halofid(0:)
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    integer(i64), intent(in) :: oldname(0:)
    real(f64), intent(in) :: airDiamond(0:)
    real(f64), intent(in) :: f_1(0:,0:)
    real(f64), intent(in) :: f_2(0:,0:)
    real(f64), intent(in) :: f_3(0:,0:)
    real(f64), intent(in) :: f_4(0:,0:)
    real(f64), intent(in) :: normalf(0:,0:)
    real(f64), intent(in) :: shift(0:,0:)
    real(f64), intent(in) :: Pbordnode(0:)
    real(f64), intent(in) :: Pbordface(0:)
    real(f64), intent(inout) :: Px_face(0:)
    real(f64), intent(inout) :: Py_face(0:)
    real(f64), intent(in) :: Pz_face(0:)
    integer(i64), intent(in) :: BCdirichlet(0:)
    integer(i64), intent(in) :: innerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: neumannfaces(0:)
    integer(i64), intent(in) :: dirichletfaces(0:)
    integer(i64), intent(in) :: periodicfaces(0:)
    integer(i64) :: i
    integer(i64) :: c_left
    integer(i64) :: c_right
    integer(i64) :: i_1
    integer(i64) :: i_2
    real(f64) :: vi1
    real(f64) :: vi2
    real(f64) :: vv1
    real(f64) :: vv2
    real(f64) :: VK
    integer(i64) :: Dummy_0065
    integer(i64) :: Dummy_0066
    integer(i64) :: Dummy_0067
    integer(i64) :: Dummy_0068
    integer(i64) :: Dummy_0069

    do Dummy_0065 = 0_i64, size(innerfaces, kind=i64) - 1_i64, 1_i64
      i = innerfaces(Dummy_0065)
      c_left = cellidf(0_i64, i)
      c_right = cellidf(1_i64, i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      vi1 = P_node(i_1)
      if (search_element(BCdirichlet, oldname(i_1)) == 1_i64) then
        vi1 = Pbordnode(i_1)
      end if
      vi2 = P_node(i_2)
      if (search_element(BCdirichlet, oldname(i_2)) == 1_i64) then
        vi2 = Pbordnode(i_2)
      end if
      vv1 = P_c(c_left)
      vv2 = P_c(c_right)
      Px_face(i) = (-1_i64) / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * &
            f_1(1_i64, i) + (vv1 + vi2) * f_2(1_i64, i) + (vi2 + vv2) * &
            f_3(1_i64, i) + (vv2 + vi1) * f_4(1_i64, i))
      Py_face(i) = 1_i64 / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * f_1( &
            0_i64, i) + (vv1 + vi2) * f_2(0_i64, i) + (vi2 + vv2) * f_3 &
            (0_i64, i) + (vv2 + vi1) * f_4(0_i64, i))
    end do
    do Dummy_0066 = 0_i64, size(periodicfaces, kind=i64) - 1_i64, 1_i64
      i = periodicfaces(Dummy_0066)
      c_left = cellidf(0_i64, i)
      c_right = cellidf(1_i64, i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      vi1 = P_node(i_1)
      if (search_element(BCdirichlet, oldname(i_1)) == 1_i64) then
        vi1 = Pbordnode(i_1)
      end if
      vi2 = P_node(i_2)
      if (search_element(BCdirichlet, oldname(i_2)) == 1_i64) then
        vi2 = Pbordnode(i_2)
      end if
      vv1 = P_c(c_left)
      vv2 = P_c(c_right)
      Px_face(i) = (-1_i64) / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * &
            f_1(1_i64, i) + (vv1 + vi2) * f_2(1_i64, i) + (vi2 + vv2) * &
            f_3(1_i64, i) + (vv2 + vi1) * f_4(1_i64, i))
      Py_face(i) = 1_i64 / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * f_1( &
            0_i64, i) + (vv1 + vi2) * f_2(0_i64, i) + (vi2 + vv2) * f_3 &
            (0_i64, i) + (vv2 + vi1) * f_4(0_i64, i))
    end do
    do Dummy_0067 = 0_i64, size(neumannfaces, kind=i64) - 1_i64, 1_i64
      i = neumannfaces(Dummy_0067)
      c_left = cellidf(0_i64, i)
      c_right = i
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      vi1 = P_node(i_1)
      if (search_element(BCdirichlet, oldname(i_1)) == 1_i64) then
        vi1 = Pbordnode(i_1)
      end if
      vi2 = P_node(i_2)
      if (search_element(BCdirichlet, oldname(i_2)) == 1_i64) then
        vi2 = Pbordnode(i_2)
      end if
      vv1 = P_c(c_left)
      vv2 = P_ghost(c_right)
      Px_face(i) = (-1_i64) / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * &
            f_1(1_i64, i) + (vv1 + vi2) * f_2(1_i64, i) + (vi2 + vv2) * &
            f_3(1_i64, i) + (vv2 + vi1) * f_4(1_i64, i))
      Py_face(i) = 1_i64 / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * f_1( &
            0_i64, i) + (vv1 + vi2) * f_2(0_i64, i) + (vi2 + vv2) * f_3 &
            (0_i64, i) + (vv2 + vi1) * f_4(0_i64, i))
    end do
    do Dummy_0068 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0068)
      c_left = cellidf(0_i64, i)
      c_right = halofid(i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      vi1 = P_node(i_1)
      if (search_element(BCdirichlet, oldname(i_1)) == 1_i64) then
        vi1 = Pbordnode(i_1)
      end if
      vi2 = P_node(i_2)
      if (search_element(BCdirichlet, oldname(i_2)) == 1_i64) then
        vi2 = Pbordnode(i_2)
      end if
      vv1 = P_c(c_left)
      vv2 = P_halo(c_right)
      Px_face(i) = (-1_i64) / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * &
            f_1(1_i64, i) + (vv1 + vi2) * f_2(1_i64, i) + (vi2 + vv2) * &
            f_3(1_i64, i) + (vv2 + vi1) * f_4(1_i64, i))
      Py_face(i) = 1_i64 / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * f_1( &
            0_i64, i) + (vv1 + vi2) * f_2(0_i64, i) + (vi2 + vv2) * f_3 &
            (0_i64, i) + (vv2 + vi1) * f_4(0_i64, i))
    end do
    do Dummy_0069 = 0_i64, size(dirichletfaces, kind=i64) - 1_i64, 1_i64
      i = dirichletfaces(Dummy_0069)
      c_left = cellidf(0_i64, i)
      c_right = i
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      vi1 = Pbordnode(i_1)
      vi2 = Pbordnode(i_2)
      vv1 = P_c(c_left)
      VK = Pbordface(i)
      vv2 = 2.0_f64 * VK - vv1
      Px_face(i) = (-1_i64) / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * &
            f_1(1_i64, i) + (vv1 + vi2) * f_2(1_i64, i) + (vi2 + vv2) * &
            f_3(1_i64, i) + (vv2 + vi1) * f_4(1_i64, i))
      Py_face(i) = 1_i64 / (2_i64 * airDiamond(i)) * ((vi1 + vv1) * f_1( &
            0_i64, i) + (vv1 + vi2) * f_2(0_i64, i) + (vi2 + vv2) * f_3 &
            (0_i64, i) + (vv2 + vi1) * f_4(0_i64, i))
    end do

  end subroutine compute_P_gradient_2d
  !........................................

  !........................................
  subroutine compute_P_gradient_3d(val_c, v_ghost, v_halo, v_node, &
        cellidf, nodeidf, centergf, namef, halofid, centerc, centerh, &
        oldname, airDiamond, n1, n2, n3, n4, normalf, shift, Pbordnode, &
        Pbordface, Px_face, Py_face, Pz_face, BCdirichlet, innerfaces, &
        halofaces, neumannfaces, dirichletfaces, periodicfaces)

    implicit none

    real(f64), intent(in) :: val_c(0:)
    real(f64), intent(in) :: v_ghost(0:)
    real(f64), intent(in) :: v_halo(0:)
    real(f64), intent(in) :: v_node(0:)
    integer(i64), intent(in) :: cellidf(0:,0:)
    integer(i64), intent(in) :: nodeidf(0:,0:)
    real(f64), intent(in) :: centergf(0:,0:)
    integer(i64), intent(in) :: namef(0:)
    integer(i64), intent(in) :: halofid(0:)
    real(f64), intent(in) :: centerc(0:,0:)
    real(f64), intent(in) :: centerh(0:,0:)
    integer(i64), intent(in) :: oldname(0:)
    real(f64), intent(in) :: airDiamond(0:)
    real(f64), intent(in) :: n1(0:,0:)
    real(f64), intent(in) :: n2(0:,0:)
    real(f64), intent(in) :: n3(0:,0:)
    real(f64), intent(in) :: n4(0:,0:)
    real(f64), intent(in) :: normalf(0:,0:)
    real(f64), intent(in) :: shift(0:,0:)
    real(f64), intent(in) :: Pbordnode(0:)
    real(f64), intent(in) :: Pbordface(0:)
    real(f64), intent(inout) :: Px_face(0:)
    real(f64), intent(inout) :: Py_face(0:)
    real(f64), intent(inout) :: Pz_face(0:)
    integer(i64), intent(in) :: BCdirichlet(0:)
    integer(i64), intent(in) :: innerfaces(0:)
    integer(i64), intent(in) :: halofaces(0:)
    integer(i64), intent(in) :: neumannfaces(0:)
    integer(i64), intent(in) :: dirichletfaces(0:)
    integer(i64), intent(in) :: periodicfaces(0:)
    integer(i64) :: i
    integer(i64) :: c_left
    integer(i64) :: c_right
    integer(i64) :: i_1
    integer(i64) :: i_2
    integer(i64) :: i_3
    integer(i64) :: i_4
    real(f64) :: V_A
    real(f64) :: V_B
    real(f64) :: V_C
    real(f64) :: V_D
    real(f64) :: V_L
    real(f64) :: V_R
    real(f64) :: V_K
    integer(i64) :: Dummy_0070
    integer(i64) :: Dummy_0071
    integer(i64) :: Dummy_0072
    integer(i64) :: Dummy_0073
    integer(i64) :: Dummy_0074

    do Dummy_0070 = 0_i64, size(innerfaces, kind=i64) - 1_i64, 1_i64
      i = innerfaces(Dummy_0070)
      c_left = cellidf(0_i64, i)
      c_right = cellidf(1_i64, i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      i_3 = nodeidf(2_i64, i)
      i_4 = i_3
      V_A = v_node(i_1)
      if (search_element(BCdirichlet, oldname(i_1)) == 1_i64) then
        V_A = Pbordnode(i_1)
      end if
      V_B = v_node(i_2)
      if (search_element(BCdirichlet, oldname(i_2)) == 1_i64) then
        V_B = Pbordnode(i_2)
      end if
      V_C = v_node(i_3)
      if (search_element(BCdirichlet, oldname(i_3)) == 1_i64) then
        V_C = Pbordnode(i_3)
      end if
      V_D = v_node(i_4)
      if (search_element(BCdirichlet, oldname(i_4)) == 1_i64) then
        V_D = Pbordnode(i_4)
      end if
      V_L = val_c(c_left)
      V_R = val_c(c_right)
      Px_face(i) = (-1.0_f64) * (n1(0_i64, i) * (V_A - V_C) + n2(0_i64, &
            i) * (V_B - V_D) + normalf(0_i64, i) * (V_R - V_L)) / &
            airDiamond(i)
      Py_face(i) = (-1.0_f64) * (n1(1_i64, i) * (V_A - V_C) + n2(1_i64, &
            i) * (V_B - V_D) + normalf(1_i64, i) * (V_R - V_L)) / &
            airDiamond(i)
      Pz_face(i) = (-1.0_f64) * (n1(2_i64, i) * (V_A - V_C) + n2(2_i64, &
            i) * (V_B - V_D) + normalf(2_i64, i) * (V_R - V_L)) / &
            airDiamond(i)
    end do
    do Dummy_0071 = 0_i64, size(periodicfaces, kind=i64) - 1_i64, 1_i64
      i = periodicfaces(Dummy_0071)
      c_left = cellidf(0_i64, i)
      c_right = cellidf(1_i64, i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      i_3 = nodeidf(2_i64, i)
      i_4 = i_3
      V_A = v_node(i_1)
      if (search_element(BCdirichlet, oldname(i_1)) == 1_i64) then
        V_A = Pbordnode(i_1)
      end if
      V_B = v_node(i_2)
      if (search_element(BCdirichlet, oldname(i_2)) == 1_i64) then
        V_B = Pbordnode(i_2)
      end if
      V_C = v_node(i_3)
      if (search_element(BCdirichlet, oldname(i_3)) == 1_i64) then
        V_C = Pbordnode(i_3)
      end if
      V_D = v_node(i_4)
      if (search_element(BCdirichlet, oldname(i_4)) == 1_i64) then
        V_D = Pbordnode(i_4)
      end if
      V_L = val_c(c_left)
      V_R = val_c(c_right)
      Px_face(i) = (-1.0_f64) * (n1(0_i64, i) * (V_A - V_C) + n2(0_i64, &
            i) * (V_B - V_D) + normalf(0_i64, i) * (V_R - V_L)) / &
            airDiamond(i)
      Py_face(i) = (-1.0_f64) * (n1(1_i64, i) * (V_A - V_C) + n2(1_i64, &
            i) * (V_B - V_D) + normalf(1_i64, i) * (V_R - V_L)) / &
            airDiamond(i)
      Pz_face(i) = (-1.0_f64) * (n1(2_i64, i) * (V_A - V_C) + n2(2_i64, &
            i) * (V_B - V_D) + normalf(2_i64, i) * (V_R - V_L)) / &
            airDiamond(i)
    end do
    do Dummy_0072 = 0_i64, size(neumannfaces, kind=i64) - 1_i64, 1_i64
      i = neumannfaces(Dummy_0072)
      c_left = cellidf(0_i64, i)
      c_right = i
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      i_3 = nodeidf(2_i64, i)
      i_4 = i_3
      V_A = v_node(i_1)
      if (search_element(BCdirichlet, oldname(i_1)) == 1_i64) then
        V_A = Pbordnode(i_1)
      end if
      V_B = v_node(i_2)
      if (search_element(BCdirichlet, oldname(i_2)) == 1_i64) then
        V_B = Pbordnode(i_2)
      end if
      V_C = v_node(i_3)
      if (search_element(BCdirichlet, oldname(i_3)) == 1_i64) then
        V_C = Pbordnode(i_3)
      end if
      V_D = v_node(i_4)
      if (search_element(BCdirichlet, oldname(i_4)) == 1_i64) then
        V_D = Pbordnode(i_4)
      end if
      V_L = val_c(c_left)
      V_R = v_ghost(c_right)
      Px_face(i) = (-1.0_f64) * (n1(0_i64, i) * (V_A - V_C) + n2(0_i64, &
            i) * (V_B - V_D) + normalf(0_i64, i) * (V_R - V_L)) / &
            airDiamond(i)
      Py_face(i) = (-1.0_f64) * (n1(1_i64, i) * (V_A - V_C) + n2(1_i64, &
            i) * (V_B - V_D) + normalf(1_i64, i) * (V_R - V_L)) / &
            airDiamond(i)
      Pz_face(i) = (-1.0_f64) * (n1(2_i64, i) * (V_A - V_C) + n2(2_i64, &
            i) * (V_B - V_D) + normalf(2_i64, i) * (V_R - V_L)) / &
            airDiamond(i)
    end do
    do Dummy_0073 = 0_i64, size(halofaces, kind=i64) - 1_i64, 1_i64
      i = halofaces(Dummy_0073)
      c_left = cellidf(0_i64, i)
      c_right = halofid(i)
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      i_3 = nodeidf(2_i64, i)
      i_4 = i_3
      V_A = v_node(i_1)
      if (search_element(BCdirichlet, oldname(i_1)) == 1_i64) then
        V_A = Pbordnode(i_1)
      end if
      V_B = v_node(i_2)
      if (search_element(BCdirichlet, oldname(i_2)) == 1_i64) then
        V_B = Pbordnode(i_2)
      end if
      V_C = v_node(i_3)
      if (search_element(BCdirichlet, oldname(i_3)) == 1_i64) then
        V_C = Pbordnode(i_3)
      end if
      V_D = v_node(i_4)
      if (search_element(BCdirichlet, oldname(i_4)) == 1_i64) then
        V_D = Pbordnode(i_4)
      end if
      V_L = val_c(c_left)
      V_R = v_halo(c_right)
      Px_face(i) = (-1.0_f64) * (n1(0_i64, i) * (V_A - V_C) + n2(0_i64, &
            i) * (V_B - V_D) + normalf(0_i64, i) * (V_R - V_L)) / &
            airDiamond(i)
      Py_face(i) = (-1.0_f64) * (n1(1_i64, i) * (V_A - V_C) + n2(1_i64, &
            i) * (V_B - V_D) + normalf(1_i64, i) * (V_R - V_L)) / &
            airDiamond(i)
      Pz_face(i) = (-1.0_f64) * (n1(2_i64, i) * (V_A - V_C) + n2(2_i64, &
            i) * (V_B - V_D) + normalf(2_i64, i) * (V_R - V_L)) / &
            airDiamond(i)
    end do
    do Dummy_0074 = 0_i64, size(dirichletfaces, kind=i64) - 1_i64, 1_i64
      i = dirichletfaces(Dummy_0074)
      c_left = cellidf(0_i64, i)
      c_right = i
      i_1 = nodeidf(0_i64, i)
      i_2 = nodeidf(1_i64, i)
      i_3 = nodeidf(2_i64, i)
      i_4 = i_3
      V_A = Pbordnode(i_1)
      V_B = Pbordnode(i_2)
      V_C = Pbordnode(i_3)
      V_D = Pbordnode(i_4)
      V_L = val_c(c_left)
      V_K = Pbordface(i)
      V_R = 2.0_f64 * V_K - V_L
      Px_face(i) = (-1.0_f64) * (n1(0_i64, i) * (V_A - V_C) + n2(0_i64, &
            i) * (V_B - V_D) + normalf(0_i64, i) * (V_R - V_L)) / &
            airDiamond(i)
      Py_face(i) = (-1.0_f64) * (n1(1_i64, i) * (V_A - V_C) + n2(1_i64, &
            i) * (V_B - V_D) + normalf(1_i64, i) * (V_R - V_L)) / &
            airDiamond(i)
      Pz_face(i) = (-1.0_f64) * (n1(2_i64, i) * (V_A - V_C) + n2(2_i64, &
            i) * (V_B - V_D) + normalf(2_i64, i) * (V_R - V_L)) / &
            airDiamond(i)
    end do

  end subroutine compute_P_gradient_3d
  !........................................

  !........................................
  subroutine facetocell(u_face, u_c, faceidc, dim) 

    implicit none

    real(f64), intent(in) :: u_face(0:)
    real(f64), intent(inout) :: u_c(0:)
    integer(i64), intent(in) :: faceidc(0:,0:)
    integer(i64), value :: dim
    integer(i64) :: nbelements
    integer(i64) :: i
    integer(i64) :: j

    nbelements = size(u_c, kind=i64)
    u_c(:) = 0.0_f64
    do i = 0_i64, nbelements - 1_i64, 1_i64
      do j = 0_i64, dim + 1_i64 - 1_i64, 1_i64
        u_c(i) = u_c(i) + u_face(faceidc(j, i))
      end do
    end do
    u_c(:) = u_c(:) / (dim + 1_i64)

  end subroutine facetocell
  !........................................

end module pyccel_functions
