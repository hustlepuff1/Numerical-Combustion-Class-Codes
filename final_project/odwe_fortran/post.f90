module post_mod
  implicit none
contains

  subroutine write_snapshot(U, ni, nj, x, y, isnap)
    implicit none
    integer, intent(in) :: ni, nj, isnap
    real(8), intent(in) :: U(4, ni, nj)
    real(8), intent(in) :: x(ni), y(nj)
    integer :: i, j
    character(len=32) :: fname

    ! file name: snap_00001.txt
    write(fname, '("snap_", I5.5, ".txt")') isnap

    open(20, file=fname, status='replace', action='write')

    ! Format:
    !   2 integers + 4 floating values (scientific notation, fixed width)
    do j = 1, nj
      do i = 1, ni
        write(20, '(2I6, 4(1X, ES23.15E3))') i, j, U(1,i,j), U(2,i,j), U(3,i,j), U(4,i,j)
      end do
    end do

    close(20)
    write(*,'(" Wrote snapshot file:",A)') trim(fname)
  end subroutine write_snapshot

end module post_mod
