program combustion_solver
    use mod_global
    use mod_chem
    use mod_linalg
    implicit none

    ! --- Variables ---
    real(dp) :: t, t_end, dt
    ! NEW: Variables for output control
    real(dp) :: t_next_write, dt_write 
    
    real(dp), dimension(4) :: Y      ! Species [A, B, C, D]
    real(dp) :: T_curr               ! Temperature
    real(dp), dimension(5) :: U      ! State Vector [YA, YB, YC, YD, T]
    real(dp), dimension(5) :: S      ! Source Vector
    real(dp), dimension(5) :: dU     ! Delta U (Implicit update)
    real(dp), dimension(5,5) :: J    ! Jacobian Matrix
    real(dp), dimension(5,5) :: LHS  ! Left Hand Side Matrix (I - dt*J)
    
    integer :: step, max_steps, method, i
    
    ! --- Setup Simulation ---
    print *, "Enter Solver Method (1=Explicit, 2=Implicit):"
    read *, method
    
    t = 0.0_dp
    t_end = 0.3e-9_dp     ! 0.3 nanoseconds total time

    if (method == 1) then
        dt = 1.0e-19_dp   ! 1 femtosecond (Explicit)
        print *, "Using Explicit Scheme with dt =", dt
    else
        dt = 1.0e-13_dp   ! 0.1 picoseconds (Implicit)
        print *, "Using Implicit Scheme with dt =", dt
    end if

    ! --- INTELLIGENT OUTPUT SETUP ---
    ! We want exactly 1000 points in our CSV, regardless of the solver dt
    dt_write = t_end / 1000.0_dp
    t_next_write = 0.0_dp
    
    ! --- Initial Conditions ---
    U(1) = 1.0_dp
    U(2) = 0.0_dp
    U(3) = 0.0_dp
    U(4) = 0.0_dp
    U(5) = 1500.0_dp     
    
    open(unit=10, file='results.csv', status='replace')
    write(10,*) 'Time,Y_A,Y_B,Y_C,Y_D,Temp' 

    ! --- Time Loop ---
    step = 0
    do while (t < t_end)
        
        ! --- NEW WRITE LOGIC ---
        ! Write if we passed the target time OR it's the very first step
        if (t >= t_next_write .or. step == 0) then
             write(10, "(ES14.6, 5(','ES16.8))") t, U(1), U(2), U(3), U(4), U(5)
             t_next_write = t_next_write + dt_write
        end if
        
        ! Extract current state
        Y = U(1:4)
        T_curr = U(5)
        
        if (method == 1) then
            ! --- Explicit Euler ---
            call get_source_term(T_curr, Y, S)
            U = U + dt * S
            
        else
            ! --- Implicit Euler ---
            call get_source_term(T_curr, Y, S)
            call get_jacobian(T_curr, Y, J)
            
            ! LHS = I - dt*J
            LHS = 0.0_dp
            do i = 1, 5
                LHS(i,i) = 1.0_dp
            end do
            LHS = LHS - (dt * J)
            
            ! RHS = dt * S
            dU = dt * S
            
            call solve_linear_system(5, LHS, dU, dU)
            U = U + dU
        end if
        
        ! --- Check for Sanity ---
        if (U(5) > 6000.0_dp .or. U(5) < 0.0_dp) then
            print *, "Simulation Unstable! T =", U(5)
            exit
        end if
        
        t = t + dt
        step = step + 1
        
    end do
    
    close(10)
    print *, "Simulation Complete. Results written to results.csv"

end program combustion_solver