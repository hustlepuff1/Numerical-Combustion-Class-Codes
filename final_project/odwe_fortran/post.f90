module post_mod
  use state_mod,    only: NVAR, NSPEC
  use chemistry_mod, only: Yspec
  implicit none
contains

  subroutine write_snapshot(U, ni, nj, x, y, isnap)
    implicit none
    integer, intent(in) :: ni, nj, isnap
    real(8), intent(in) :: U(NVAR, ni, nj)
    real(8), intent(in) :: x(ni), y(nj)
    integer :: i, j
    character(len=32) :: fname, fnameY

      !----------------------------------------
      ! Prevent "unused variable" remarks (#7712)
      !----------------------------------------
      if (.false.) then
        print *, x(1), y(1)
      end if

    ! ---------- conservative variables ----------
    write(fname, '("snap_", I5.5, ".txt")') isnap
    open(20, file=fname, status='replace', action='write')
    do j = 1, nj
      do i = 1, ni
        write(20, '(2I6, 4(1X, ES23.15E3))') i, j, &
             U(1,i,j), U(2,i,j), U(3,i,j), U(4,i,j)
      end do
    end do
    close(20)

    ! ---------- species mass fractions ----------
    ! Format: i, j, Y(1:NSPEC)
    write(fnameY, '("snapY_", I5.5, ".txt")') isnap
    open(21, file=fnameY, status='replace', action='write')
    do j = 1, nj
      do i = 1, ni
        write(21, '(2I6, 7(1X, ES23.15E3))') i, j, &
             Yspec(1,i,j), Yspec(2,i,j), Yspec(3,i,j), &
             Yspec(4,i,j), Yspec(5,i,j), Yspec(6,i,j), Yspec(7,i,j)
      end do
    end do
    close(21)

    write(*,'(" Wrote snapshot file:",A)') trim(fname)
  end subroutine write_snapshot

end module post_mod
