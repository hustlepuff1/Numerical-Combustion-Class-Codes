program combustion_solver
    use mod_global
    use mod_chem
    use mod_linalg
    implicit none

    ! --- Variables ---
    real(dp) :: t, t_end, dt
    real(dp), dimension(4) :: Y      ! Species [A, B, C, D]
    real(dp) :: T_curr               ! Temperature
    real(dp), dimension(5) :: U      ! State Vector [YA, YB, YC, YD, T]
    real(dp), dimension(5) :: S      ! Source Vector
    real(dp), dimension(5) :: dU     ! Delta U (Implicit update)
    real(dp), dimension(5,5) :: J    ! Jacobian Matrix
    real(dp), dimension(5,5) :: LHS  ! Left Hand Side Matrix (I - dt*J)
    
    integer :: step, max_steps, method, i
    
    ! --- Setup Simulation ---
    ! method = 1: Explicit Euler
    ! method = 2: Linearized Implicit Euler
    
    print *, "Enter Solver Method (1=Explicit, 2=Implicit):"
    read *, method
    
    if (method == 1) then
        dt = 1.0e-11_dp   ! Tiny step for Explicit (Stability limited)
        print *, "Using Explicit Scheme with dt =", dt
    else
        dt = 1.0e-9_dp   ! Larger step for Implicit (Stable)
        print *, "Using Implicit Scheme with dt =", dt
    end if

    t = 0.0_dp
    t_end = 0.00003_dp      ! Simulate 0.1 seconds (adjust if needed)
    
    ! --- Initial Conditions ---
    ! Y = [1.0, 0.0, 0.0, 0.0] (Pure Fuel)
    U(1) = 1.0_dp
    U(2) = 0.0_dp
    U(3) = 0.0_dp
    U(4) = 0.0_dp
    ! Initial Temperature (Must be high enough to ignite!)
    U(5) = 1500.0_dp     
    
    ! --- Open Output File ---
    open(unit=10, file='results.csv', status='replace')
    write(10,*) 'Time,Y_A,Y_B,Y_C,Y_D,Temp' ! Header

    ! --- Time Loop ---
    step = 0
    do while (t < t_end)
        
        ! Write output every 10 steps (to save file size)
        if (mod(step, 10) == 0) then
             write(10, "(F12.6, 5(','F16.8))") t, U(1), U(2), U(3), U(4), U(5)
        end if
        
        ! Extract current state for helper functions
        Y = U(1:4)
        T_curr = U(5)
        
        ! --- SOLVER LOGIC ---
        call get_source_term(T_curr, Y, S)
        
        if (method == 1) then
            ! --- Explicit Euler: U_new = U_old + dt * S ---
            U = U + dt * S
            
        else
            ! --- Implicit Euler: (I - dt*J) * dU = dt * S ---
            
            ! 1. Get Jacobian J evaluated at current step
            call get_jacobian(T_curr, Y, J)
            
            ! 2. Form LHS Matrix A = (Identity - dt * J)
            LHS = 0.0_dp
            do i = 1, 5
                LHS(i,i) = 1.0_dp
            end do
            LHS = LHS - (dt * J)
            
            ! 3. Form RHS Vector b = dt * S
            dU = dt * S  ! Re-using dU variable as 'b' momentarily
            
            ! 4. Solve System A * x = b
            ! Result 'dU' will hold the change vector
            call solve_linear_system(5, LHS, dU, dU)
            
            ! 5. Update U_new = U_old + dU
            U = U + dU
            
        end if
        
        ! --- Check for Sanity ---
        if (U(5) > 5000.0_dp .or. U(5) < 300.0_dp) then
            print *, "Simulation Unstable! T =", U(5)
            exit
        end if
        
        ! Update Time
        t = t + dt
        step = step + 1
        
        ! Progress bar (simple)
        if (mod(step, 1000) == 0) print *, "Time: ", t, " Temp: ", U(5)

    end do
    
    close(10)
    print *, "Simulation Complete. Results written to results.csv"

end program combustion_solver