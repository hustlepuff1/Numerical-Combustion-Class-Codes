! Program to solve for p2 for Example 5.1 using the user's provided equation.
program hw3_example_5_1
    implicit none
    integer, parameter :: dp = kind(1.0d0)

    ! --- Secant Method variables ---
    real(dp) :: p_old, p_current, p_next
    real(dp) :: f_old, f_current
    real(dp) :: tolerance = 1.0d-6
    integer  :: max_iter = 30, i
    real(dp) :: final_u2

    ! --- Main Program ---
    ! Initial guesses for p2 (must be between p1 and p4)
    p_old = 20000.0_dp
    p_current = 80000.0_dp

    write(*,*) "--- Solving Example 5.1 ---"
    do i = 1, max_iter
        f_old = F(p_old)
        f_current = F(p_current)
        if (abs(f_current - f_old) < 1.0d-12) then
            write(*,*) "Error: Stalled."; exit
        end if
        p_next = p_current - f_current * (p_current - p_old) / (f_current - f_old)
        if (abs(p_next - p_current) < tolerance) then
            p_current = p_next
            write(*,'(A, I3, A)') "Convergence reached in ", i, " iterations."
            exit
        end if
        p_old = p_current
        p_current = p_next
    end do

    write(*,'(A, F12.4, A)') "Final p2 = ", p_current, " Pa"
    
    ! --- Verification Step ---
    ! Calculate u2 using the solved p2 and the shock relation (Eq. 5.4)
    final_u2 = calc_u2_shock(p_current)
    write(*,'(A, F12.4, A)') "Calculated u2 = ", final_u2, " m/s"
    write(*,'(A, F12.4, A)') "Book's answer for u2 = 329.5 m/s"

contains

    ! The function F(p2)=0 based on the user's typed equation.
    function F(p2)
        implicit none
        real(dp), intent(in) :: p2
        real(dp)             :: F
        real(dp), parameter :: p4=1d5, p1=1d4, u4=100d0, u1=-50d0, gam=1.4d0
        real(dp), parameter :: a4=374.1657d0, a1=334.6640d0
        real(dp) :: x, term_in_sqrt, bracket_term
        x = p2 / p1
        term_in_sqrt = ((gam+1d0)/(2d0*gam)) * (x - 1d0) + 1d0
        if (term_in_sqrt < 0) then; write(*,*) "Error: sqrt of negative"; stop; end if
        bracket_term = u4 - u1 - (a1/gam) * (x-1d0)/sqrt(term_in_sqrt)
        F = x * (1d0 + ((gam-1d0)/(2d0*a4)) * bracket_term)**(-2d0*gam/(gam-1d0)) - (p4/p1)
    end function F

    ! Function to calculate u2 from p2 using the shock relation (Eq. 5.4)
    function calc_u2_shock(p2)
        implicit none
        real(dp), intent(in) :: p2
        real(dp)             :: calc_u2_shock
        real(dp), parameter :: p1=1d4, rho1=0.125d0, u1=-50d0, gam=1.4d0
        calc_u2_shock = u1 + (p2-p1) / sqrt(rho1 * (((gam+1d0)/2d0)*p2 + ((gam-1d0)/2d0)*p1))
    end function calc_u2_shock

end program hw3_example_5_1