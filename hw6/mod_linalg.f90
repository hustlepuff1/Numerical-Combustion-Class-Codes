module mod_linalg
    use mod_global
    implicit none

contains

    ! Simple Gaussian Elimination solver for Ax = b
    ! N is 5 for this problem
    subroutine solve_linear_system(N, A, b, x)
        integer, intent(in) :: N
        real(dp), dimension(N,N), intent(inout) :: A ! Will be modified
        real(dp), dimension(N), intent(inout) :: b   ! Will be modified
        real(dp), dimension(N), intent(out) :: x
        
        integer :: i, j, k, max_row
        real(dp) :: factor, temp, max_val

        ! Forward Elimination
        do k = 1, N-1
            ! Partial Pivoting (Find largest pivot to avoid division by zero)
            max_val = abs(A(k,k))
            max_row = k
            do i = k+1, N
                if (abs(A(i,k)) > max_val) then
                    max_val = abs(A(i,k))
                    max_row = i
                end if
            end do

            ! Swap rows if needed
            if (max_row /= k) then
                do j = k, N
                    temp = A(k,j)
                    A(k,j) = A(max_row,j)
                    A(max_row,j) = temp
                end do
                temp = b(k)
                b(k) = b(max_row)
                b(max_row) = temp
            end if

            ! Eliminate
            do i = k+1, N
                factor = A(i,k) / A(k,k)
                do j = k, N
                    A(i,j) = A(i,j) - factor * A(k,j)
                end do
                b(i) = b(i) - factor * b(k)
            end do
        end do

        ! Backward Substitution
        x(N) = b(N) / A(N,N)
        do i = N-1, 1, -1
            temp = b(i)
            do j = i+1, N
                temp = temp - A(i,j) * x(j)
            end do
            x(i) = temp / A(i,i)
        end do

    end subroutine solve_linear_system

end module mod_linalg