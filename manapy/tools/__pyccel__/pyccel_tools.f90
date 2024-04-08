module pyccel_tools


  use, intrinsic :: ISO_C_Binding, only : i64 => C_INT64_T , f64 => &
        C_DOUBLE
  implicit none

  contains

  !........................................
  subroutine initialisation_gaussian_2d(ne, u, v, w, P, center, Pinit) 

    implicit none

    real(f64), intent(inout) :: ne(0:)
    real(f64), intent(inout) :: u(0:)
    real(f64), intent(inout) :: v(0:)
    real(f64), intent(inout) :: w(0:)
    real(f64), intent(inout) :: P(0:)
    real(f64), intent(in) :: center(0:,0:)
    real(f64), value :: Pinit
    integer(i64) :: nbelements
    real(f64) :: sigma
    integer(i64) :: i
    real(f64) :: xcent
    real(f64) :: ycent

    nbelements = size(center,2,i64)
    sigma = 0.05_f64
    do i = 0_i64, nbelements - 1_i64, 1_i64
      xcent = center(0_i64, i)
      ycent = center(1_i64, i)
      ne(i) = 5_i64 * exp((-1.0_f64) * ((xcent - 0.2_f64) ** 2_i64 + ( &
            ycent - 0.25_f64) ** 2_i64) / sigma ** 2_i64) + 1_i64
      u(i) = 0.0_f64
      v(i) = 0.0_f64
      w(i) = 0.0_f64
      P(i) = Pinit * (0.5_f64 - xcent)
    end do

  end subroutine initialisation_gaussian_2d
  !........................................

  !........................................
  subroutine initialisation_gaussian_3d(ne, u, v, w, P, center, Pinit) 

    implicit none

    real(f64), intent(inout) :: ne(0:)
    real(f64), intent(inout) :: u(0:)
    real(f64), intent(inout) :: v(0:)
    real(f64), intent(inout) :: w(0:)
    real(f64), intent(inout) :: P(0:)
    real(f64), intent(in) :: center(0:,0:)
    real(f64), value :: Pinit
    integer(i64) :: nbelements
    real(f64) :: sigma
    integer(i64) :: i
    real(f64) :: xcent
    real(f64) :: ycent
    real(f64) :: zcent

    nbelements = size(center,2,i64)
    sigma = 0.05_f64
    do i = 0_i64, nbelements - 1_i64, 1_i64
      xcent = center(0_i64, i)
      ycent = center(1_i64, i)
      zcent = center(2_i64, i)
      ne(i) = 5_i64 * exp((-1.0_f64) * ((xcent - 0.2_f64) ** 2_i64 + ( &
            ycent - 0.25_f64) ** 2_i64 + (zcent - 0.45_f64) ** 2_i64) / &
            sigma ** 2_i64) + 1_i64
      u(i) = 0.0_f64
      v(i) = 0.0_f64
      w(i) = 0.0_f64
      P(i) = Pinit * (0.5_f64 - xcent)
    end do

  end subroutine initialisation_gaussian_3d
  !........................................

  !........................................
  function time_step(u, v, w, cfl, normal, mesure, volume, faceid, dim, &
        Dxx, Dyy, Dzz) result(dt)

    implicit none

    real(f64) :: dt
    real(f64), intent(in) :: u(0:)
    real(f64), intent(in) :: v(0:)
    real(f64), intent(in) :: w(0:)
    real(f64), value :: cfl
    real(f64), intent(in) :: normal(0:,0:)
    real(f64), intent(in) :: mesure(0:)
    real(f64), intent(in) :: volume(0:)
    integer(i64), intent(in) :: faceid(0:,0:)
    integer(i64), value :: dim
    real(f64), value :: Dxx
    real(f64), value :: Dyy
    real(f64), value :: Dzz
    integer(i64) :: nbelement
    real(f64) :: u_n
    real(f64) :: norm(0:2_i64)
    integer(i64) :: i
    real(f64) :: lam
    integer(i64) :: j
    real(f64) :: lam_convect
    real(f64) :: mes
    real(f64) :: lam_diff

    nbelement = size(faceid,2,i64)
    u_n = 0.0_f64
    norm = 0.0_f64
    dt = 1000000.0_f64
    do i = 0_i64, nbelement - 1_i64, 1_i64
      lam = 0.0_f64
      do j = 0_i64, dim + 1_i64 - 1_i64, 1_i64
        norm(:) = normal(:, faceid(j, i))
        u_n = abs(u(i) * norm(0_i64) + v(i) * norm(1_i64) + w(i) * norm( &
              2_i64))
        lam_convect = u_n / mesure(faceid(j, i))
        lam = lam + lam_convect * mesure(faceid(j, i))
        mes = sqrt(norm(0_i64) * norm(0_i64) + norm(1_i64) * norm(1_i64 &
              ) + norm(2_i64) * norm(2_i64))
        lam_diff = Dxx * mes ** 2_i64 + Dyy * mes ** 2_i64 + Dzz * mes &
              ** 2_i64
        lam = lam + lam_diff / volume(i)
      end do
      dt = minval([dt, cfl * volume(i) / lam])
    end do
    return

  end function time_step
  !........................................

  !........................................
  subroutine update_new_value(ne_c, u_c, v_c, P_c, rez_ne, dissip_ne, &
        src_ne, dtime, vol)

    implicit none

    real(f64), intent(inout) :: ne_c(0:)
    real(f64), intent(in) :: u_c(0:)
    real(f64), intent(in) :: v_c(0:)
    real(f64), intent(in) :: P_c(0:)
    real(f64), intent(in) :: rez_ne(0:)
    real(f64), intent(in) :: dissip_ne(0:)
    real(f64), intent(in) :: src_ne(0:)
    real(f64), value :: dtime
    real(f64), intent(in) :: vol(0:)

    ne_c(:) = ne_c(:) + dtime * ((rez_ne(:) + dissip_ne(:)) / vol(:) + &
          src_ne(:))

  end subroutine update_new_value
  !........................................

end module pyccel_tools
