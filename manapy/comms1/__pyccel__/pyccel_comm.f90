module pyccel_comm


  use, intrinsic :: ISO_C_Binding, only : i64 => C_INT64_T , f64 => &
        C_DOUBLE
  implicit none

  contains

  !........................................
  subroutine define_halosend(w_c, w_halosend, indsend) 

    implicit none

    real(f64), intent(in) :: w_c(0:)
    real(f64), intent(inout) :: w_halosend(0:)
    integer(i64), intent(in) :: indsend(0:)

    w_halosend(:) = w_c(indsend(:))

  end subroutine define_halosend
  !........................................

end module pyccel_comm
