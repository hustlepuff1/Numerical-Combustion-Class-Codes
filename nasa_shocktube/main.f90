program nasa_solver
    use constants_mod
    use flux_mod
    use solver_mod
    implicit none
    
    integer :: step, frame_count
    real(dp) :: time, dt
    
    print *, "==================================="
    print *, "    NASA Channel: PRODUCTION       "
    print *, "==================================="
    
    call read_input()
    call allocate_memory()
    call initialize_flow()
    
    time = 0.0_dp
    frame_count = 0
    
    ! Write Initial State
    call write_output(frame_count)
    
    print *, "Starting Time Integration..."
    
    !$OMP PARALLEL
    !$OMP MASTER
        print *, "OpenMP Threads Active."
    !$OMP END MASTER
    !$OMP END PARALLEL

    do step = 1, max_steps
        call compute_time_step(dt)
        if (time + dt > t_final) dt = t_final - time
        
        call compute_residuals()
        call update_solution(dt, step)
        call apply_boundary_conditions()
        
        time = time + dt
        
        if (mod(step, 100) == 0) then
            print '(A, I6, A, ES10.3, A, ES10.3)', "Step:", step, "  Time:", time, "  dt:", dt
            frame_count = frame_count + 1
            call write_output(frame_count)
        end if
        
        if (time >= t_final) exit
    end do
    
    print *, "Done."
end program nasa_solver