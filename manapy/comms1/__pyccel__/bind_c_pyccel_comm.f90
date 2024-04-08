module bind_c_pyccel_comm

  use pyccel_comm, only: define_halosend

  use, intrinsic :: ISO_C_Binding, only : i64 => C_INT64_T , f64 => &
        C_DOUBLE
  implicit none

  contains

  !........................................
  subroutine bind_c_define_halosend(n0_w_c, w_c, n0_w_halosend, &
        w_halosend, n0_indsend, indsend) bind(c)

    implicit none

    integer(i64), value :: n0_w_c
    real(f64), intent(in) :: w_c(0:n0_w_c - 1_i64)
    integer(i64), value :: n0_w_halosend
    real(f64), intent(inout) :: w_halosend(0:n0_w_halosend - 1_i64)
    integer(i64), value :: n0_indsend
    integer(i64), intent(in) :: indsend(0:n0_indsend - 1_i64)

    call define_halosend(w_c, w_halosend, indsend)

  end subroutine bind_c_define_halosend
  !........................................

end module bind_c_pyccel_comm
