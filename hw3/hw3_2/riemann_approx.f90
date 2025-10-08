! ====================================================================
! Program: Approximate Riemann Solver for Homework #3
! Purpose: Uses Roe's method to construct a snapshot of the solution.
! Note:    Follows book pages 85 of Computational Gasdynamics
! ====================================================================
PROGRAM riemann_approx_solver_hw3
    IMPLICIT NONE
    integer, parameter :: dp = kind(1.0d0)

    ! --- Variable Declarations ---
    real(dp) :: gam, R
    real(dp) :: tmax
    real(dp) :: rho, u, p, a, Ma, s
    real(dp) :: pL, rhoL, uL, htL, etL
    real(dp) :: pR, rhoR, uR, htR, etR
    real(dp) :: rhoRL, uRL, hRL, aRL
    real(dp) :: ld1, ld2, ld3
    real(dp) :: drho, dp_jump, du
    real(dp) :: dv1, dv2, dv3
    real(dp) :: r1(3), r2(3), r3(3)
    real(dp) :: u1, u2, u3
    real(dp) :: x_pos, x_diaphragm, self_similar_var
    integer :: i
    
    ! --- Set Initial Conditions for Homework #3 ---
    gam = 1.4_dp
    R = 287.0_dp
    tmax = 3.9d-3
    pL = 1.0d5
    rhoL = 1.0_dp
    uL = 0.0_dp
    pR = 1.0d3
    rhoR = 0.01_dp
    uR = 0.0_dp
    x_diaphragm = 5.0_dp
    
    ! --- Step 1: Compute derived quantities ---
    etL = pL/(rhoL*(gam-1.0_dp)) + 0.5_dp*uL**2
    htL = etL + pL/rhoL
    etR = pR/(rhoR*(gam-1.0_dp)) + 0.5_dp*uR**2
    htR = etR + pR/rhoR
    
    ! --- Step 2: Compute Roe-averaged quantities ---
    rhoRL = sqrt(rhoR*rhoL)
    uRL = (sqrt(rhoR)*uR + sqrt(rhoL)*uL) / (sqrt(rhoR)+sqrt(rhoL))
    hRL = (sqrt(rhoR)*htR + sqrt(rhoL)*htL) / (sqrt(rhoR)+sqrt(rhoL))
    aRL = sqrt((gam-1.0_dp)*(hRL - 0.5_dp*uRL**2))
    
    ! --- Step 3: Compute the Roe-average wave speeds ---
    ld1 = uRL
    ld2 = uRL + aRL
    ld3 = uRL - aRL
    
    ! --- Step 4: Compute wave strengths ---
    drho = rhoR - rhoL
    dp_jump = pR - pL
    du = uR - uL
    dv1 = drho - dp_jump/aRL**2
    dv2 = du + dp_jump/(rhoRL*aRL)
    dv3 = du - dp_jump/(rhoRL*aRL)

    ! --- Step 5: Compute right characteristic vectors ---
    r1(1) = 1.0_dp; r1(2) = uRL; r1(3) = 0.5_dp*uRL**2
    r2(1) = 1.0_dp; r2(2) = uRL + aRL; r2(3) = hRL + aRL*uRL
    r2(:) = 0.5_dp*rhoRL/aRL * r2(:)
    r3(1) = 1.0_dp; r3(2) = uRL - aRL; r3(3) = hRL - aRL*uRL
    r3(:) = -0.5_dp*rhoRL/aRL * r3(:)
    
    ! --- Step 6: Compute the solution and write to file ---
    OPEN(11, File='Riemann_Approx_HW3.dat', status='replace')
    WRITE(11,*) 'X, density, velocity, pressure, speed_of_sound, Mach_number, entropy'
    
    do i = 0, 400
        x_pos = real(i, dp) * (10.0_dp / 400.0_dp)
        self_similar_var = (x_pos - x_diaphragm) / tmax

        if (self_similar_var <= ld3) then
            u1 = rhoL; u2 = rhoL*uL; u3 = rhoL*etL
        else if (self_similar_var <= ld1) then
            u1 = rhoL + r3(1)*dv3
            u2 = rhoL*uL + r3(2)*dv3
            u3 = rhoL*etL + r3(3)*dv3
        else if (self_similar_var <= ld2) then
            u1 = rhoR - r2(1)*dv2
            u2 = rhoR*uR - r2(2)*dv2
            u3 = rhoR*etR - r2(3)*dv2
        else
            u1 = rhoR; u2 = rhoR*uR; u3 = rhoR*etR
        end if
        
        rho = u1
        u = u2 / u1
        p = (gam-1.0_dp)*(u3 - 0.5_dp*(u2**2/u1))
        a = sqrt(gam*p/rho)
        s = (R/(gam-1.0_dp))*log(p/(rho**gam))
        if (a > 1.0d-12) then
          Ma = u/a
        else
          Ma = 0.0_dp
        end if
        write(11, '(7E22.14)') x_pos, rho, u, p, a, Ma, s
    end do
    CLOSE(11)

    write(*,*) "Approximate solution written to Riemann_Approx_HW3.dat"
    
END PROGRAM riemann_approx_solver_hw3
