! ====================================================================
! Program: main_fvm_solver
! Purpose: Solves the 1D Euler equations for the Riemann problem
!          using a Finite Volume Method with Roe's solver.
! ====================================================================
program main_fvm_solver
    use roe_solver_module
    implicit none

    ! --- Grid and Simulation Parameters ---
    integer, parameter :: nx = 400         ! Number of cells
    real(dp), parameter :: x_min = 0.0_dp, x_max = 10.0_dp
    real(dp), parameter :: t_final = 3.9d-3 ! Final time (3.9 ms)
    real(dp), parameter :: cfl = 0.8        ! Courant number
    real(dp) :: x_diaphragm = 5.0_dp      ! Diaphragm position

    ! --- Arrays ---
    real(dp) :: x(nx)                   ! Cell centers
    real(dp) :: rho(nx), vel(nx), p(nx) ! Primitive variables
    real(dp) :: U_cons(3, nx)           ! Conservative variables (rho, rho*u, rho*E)
    real(dp) :: F(3, nx+1)              ! Fluxes at interfaces
    
    ! --- Variables ---
    real(dp) :: dx, dt, t
    integer :: i
    integer :: file_unit
    real(dp) :: max_wave_speed, a

    ! =================================================================
    ! --- 1. SETUP ---
    ! =================================================================
    dx = (x_max - x_min) / real(nx, dp)
    t = 0.0_dp

    ! Create cell-centered grid and initial conditions
    do i = 1, nx
        x(i) = x_min + (real(i, dp) - 0.5_dp) * dx
        if (x(i) < x_diaphragm) then ! Left state
            p(i)   = 1.0d5
            rho(i) = 1.0d0
            vel(i) = 0.0d0
        else ! Right state
            p(i)   = 1.0d3
            rho(i) = 0.01d0
            vel(i) = 0.0d0
        end if
    end do
    
    ! Convert initial primitive to conservative variables
    U_cons(1,:) = rho
    U_cons(2,:) = rho * vel
    U_cons(3,:) = p / (gamma - 1.0_dp) + 0.5_dp * rho * vel**2

    ! =================================================================
    ! --- 2. MAIN TIME-STEPPING LOOP ---
    ! =================================================================
    do while (t < t_final)
        ! --- Calculate time step from CFL condition ---
        max_wave_speed = 0.0_dp
        do i = 1, nx
            a = sqrt(gamma * p(i) / rho(i))
            max_wave_speed = max(max_wave_speed, abs(vel(i)) + a)
        end do
        dt = cfl * dx / max_wave_speed
        ! Ensure we don't step past the final time
        if (t + dt > t_final) dt = t_final - t
        
        ! --- Calculate fluxes at interfaces ---
        ! Boundary conditions: copy interior state to ghost cells
        ! Left boundary
        call roe_flux(rho(1), vel(1), p(1), rho(1), vel(1), p(1), F(:,1))
        ! Interior interfaces
        do i = 2, nx
            call roe_flux(rho(i-1), vel(i-1), p(i-1), rho(i), vel(i), p(i), F(:,i))
        end do
        ! Right boundary
        call roe_flux(rho(nx), vel(nx), p(nx), rho(nx), vel(nx), p(nx), F(:,nx+1))

        ! --- Update conservative variables in each cell ---
        do i = 1, nx
            U_cons(:,i) = U_cons(:,i) - (dt / dx) * (F(:,i+1) - F(:,i))
        end do

        ! --- Update primitive variables from new conservative variables ---
        rho = U_cons(1,:)
        vel = U_cons(2,:) / rho
        p   = (gamma - 1.0_dp) * (U_cons(3,:) - 0.5_dp * rho * vel**2)
        
        t = t + dt
        write(*,*) "Time = ", t*1000.0_dp, " ms"
    end do

    ! =================================================================
    ! --- 3. WRITE FINAL RESULTS ---
    ! =================================================================
    open(newunit=file_unit, file='results_roe.dat', status='replace')
    write(file_unit, '(A8, A15, A15, A15)') "x", "pressure", "density", "velocity"
    do i = 1, nx
        write(file_unit, '(4E15.7)') x(i), p(i), rho(i), vel(i)
    end do
    close(file_unit)
    write(*,*) "Final results written to results_roe.dat"

end program main_fvm_solver
