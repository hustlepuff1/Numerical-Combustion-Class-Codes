! --------------------------------------------------------------------
! Program: p2 solver from book equation 5.6
! Purpose: Solves for pressure p2 using the secant method.
! --------------------------------------------------------------------
program hw3_solveForP2
    implicit none
    integer, parameter :: dp = kind(1.0d0)

    ! Secant Method variables
    real(dp) :: p_old, p_current, p_next
    real(dp) :: f_old, f_current
    real(dp) :: tolerance = 1.0d-6
    integer  :: max_iter = 30
    integer  :: i

    p_old = 10000.0_dp
    p_current = 80000.0_dp

    write(*,*) "--- Starting Secant Method (Corrected Eq.) ---"
    write(*,'(A, F12.4)') "Initial guess 1 (p_old)     = ", p_old
    write(*,'(A, F12.4)') "Initial guess 2 (p_current) = ", p_current
    write(*,*) "----------------------------------"

    do i = 1, max_iter
        f_old = F(p_old)
        f_current = F(p_current)

        if (abs(f_current - f_old) < 1.0d-12) then
            write(*,*) "Error: Stalled. Denominator is too small."
            exit
        end if

        p_next = p_current - f_current * (p_current - p_old) / (f_current - f_old)

        if (abs(p_next - p_current) < tolerance) then
            p_current = p_next
            write(*,'(A, I3, A)') "Convergence reached in ", i, " iterations."
            exit
        end if

        p_old = p_current
        p_current = p_next

        write(*,'(A, I3, A, F12.4)') "Iter ", i, ": p2 = ", p_current
    end do

    write(*,*) "----------------------------------"
    write(*,'(A, F12.4, A)') "Final p2 = ", p_current, " Pa"
    write(*,*) "----------------------------------"
    write(*,'(A, F12.4, A)') "Calculated u2 = ", calc_u2_shock(p_current), " m/s"
    write(*,*) "----------------------------------"

contains
    function F(p2)
        implicit none
        real(dp), intent(in) :: p2
        real(dp)             :: F
        
        real(dp), parameter :: p1 = 1000.0_dp
        real(dp), parameter :: p4_over_p1 = 100.0_dp
        real(dp), parameter :: gam = 1.4_dp
        
        real(dp), parameter :: c1 = (gam - 1.0_dp) / (2.0_dp * gam)
        real(dp), parameter :: c2 = (gam + 1.0_dp) / (2.0_dp * gam)
        real(dp) :: x, term_in_sqrt, bracket_term, pressure_term

        x = p2 / p1

        pressure_term = (x / p4_over_p1)**c1

        term_in_sqrt = c2 * (x - 1.0_dp) + 1.0_dp
        if (term_in_sqrt < 0) then
            write(*,*) "Error: sqrt of negative value in F(p2)."
            stop
        end if

        bracket_term = (x - 1.0_dp) / sqrt(term_in_sqrt)

        F = pressure_term - 1.0_dp + c1 * bracket_term
    end function F

    function calc_u2_shock(p2)
        implicit none
        real(dp), intent(in) :: p2
        real(dp)             :: calc_u2_shock
        
        real(dp), parameter :: p1 = 1000.0_dp
        real(dp), parameter :: rho1 = 0.01_dp
        real(dp), parameter :: u1 = 0.0_dp
        real(dp), parameter :: gam = 1.4_dp

        calc_u2_shock = u1 + (p2 - p1) / sqrt(rho1 * (((gam + 1.0_dp) / 2.0_dp) * p2 + ((gam - 1.0_dp) / 2.0_dp) * p1))
    end function calc_u2_shock

end program hw3_solveForP2