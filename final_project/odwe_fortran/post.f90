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
    write(fnameY, '("snapY_", I5.5, ".txt")') isnap
    open(21, file=fnameY, status='replace', action='write')
    do j = 1, nj
      do i = 1, ni
        ! Fixed to write all 8 species
        write(21, '(2I6, 8(1X, ES23.15E3))') i, j, &
             Yspec(1,i,j), Yspec(2,i,j), Yspec(3,i,j), &
             Yspec(4,i,j), Yspec(5,i,j), Yspec(6,i,j), &
             Yspec(7,i,j), Yspec(8,i,j)
      end do
    end do
    close(21)

    write(*,'(" Wrote snapshot file:",A)') trim(fname)
  end subroutine write_snapshot

end module post_mod