module solver_mod
    use constants_mod, only: dp, gamma, R_univ, MW_H2, MW_O2, MW_H2O
    use flux_mod, only: get_flux_vectors
    use thermo_mod, only: get_temperature_from_energy
    use reaction_mod, only: reaction_subcycle
    implicit none

    integer :: ni, nj
    real(dp) :: L, H, dx, dy
    real(dp) :: CFL_target 
    
    real(dp) :: t_final      
    integer  :: max_steps
    integer  :: chem_switch
    
    real(dp) :: P_init, T_init, M_init
    real(dp) :: Y_H2_init, Y_O2_init
    
    ! Global Arrays
    real(dp) :: Q_inlet(10) 
    real(dp), allocatable :: x(:,:), y(:,:)
    real(dp), allocatable :: Q(:,:,:)    
    real(dp), allocatable :: Res(:,:,:)  

contains

    subroutine read_input()
        integer :: ios
        character(len=200) :: line_buf 
        
        open(10, file='input.txt', status='old', action='read', iostat=ios)
        read(10, '(A)') line_buf; read(line_buf, *) ni
        read(10, '(A)') line_buf; read(line_buf, *) nj
        read(10, '(A)') line_buf; read(line_buf, *) L
        read(10, '(A)') line_buf; read(line_buf, *) H
        read(10, '(A)') line_buf; read(line_buf, *) P_init
        read(10, '(A)') line_buf; read(line_buf, *) T_init
        read(10, '(A)') line_buf; read(line_buf, *) M_init
        read(10, '(A)') line_buf; read(line_buf, *) Y_H2_init
        read(10, '(A)') line_buf; read(line_buf, *) Y_O2_init
        read(10, '(A)') line_buf; read(line_buf, *) t_final
        read(10, '(A)') line_buf; read(line_buf, *) max_steps
        read(10, '(A)') line_buf; read(line_buf, *) chem_switch
        read(10, '(A)') line_buf; read(line_buf, *) CFL_target
        close(10)
        
        dx = L / real(ni - 1, dp); dy = H / real(nj - 1, dp)
        
        print *, "--- Input Read ---"
        print *, "CFL Target:", CFL_target
        if (chem_switch == 1) then
            print *, "Chemistry: ON"
        else
            print *, "Chemistry: OFF"
        end if
        
        ! Create History File (Overwrite if exists)
        open(30, file='history.dat', status='replace')
        write(30, '(A)') "Step, Time, L2_Norm_Rho"
        close(30)
    end subroutine read_input

    subroutine allocate_memory()
        allocate(Q(ni,nj,10)); allocate(Res(ni,nj,10))
        allocate(x(ni,nj)); allocate(y(ni,nj))
        Q=0.0_dp; Res=0.0_dp
    end subroutine allocate_memory

    subroutine initialize_flow()
        integer :: i, j
        real(dp) :: R_mix, rho, u, Et
        R_mix = R_univ * ( (Y_H2_init / MW_H2) + (Y_O2_init / MW_O2) )
        rho = P_init / (R_mix * T_init)
        u = M_init * sqrt(gamma * P_init / rho)
        Et = (P_init / ((gamma - 1.0_dp) * rho)) + 0.5_dp * u**2

        Q_inlet = 0.0_dp
        Q_inlet(1) = rho; Q_inlet(2) = rho * u; Q_inlet(4) = rho * Et
        Q_inlet(5) = rho * Y_H2_init; Q_inlet(6) = rho * Y_O2_init

        do j = 1, nj
            do i = 1, ni
                x(i,j) = real(i-1, dp) * dx
                y(i,j) = real(j-1, dp) * dy
                Q(i,j,:) = Q_inlet(:)
            end do
        end do
    end subroutine initialize_flow

    subroutine compute_time_step(dt_out)
        real(dp), intent(out) :: dt_out
        integer :: i, j
        real(dp) :: rho, u, v, P, a, max_eigen, global_max, e_int, Y_H2, Y_O2, T_local
        global_max = 0.0_dp
        do j = 1, nj
            do i = 1, ni
                rho = max(Q(i,j,1), 1.0e-10_dp)
                u = Q(i,j,2)/rho; v = Q(i,j,3)/rho
                e_int = Q(i,j,4)/rho - 0.5_dp*(u**2 + v**2)
                Y_H2 = Q(i,j,5)/rho; Y_O2 = Q(i,j,6)/rho
                
                call get_temperature_from_energy(e_int, Y_H2, Y_O2, 1500.0_dp, T_local)
                P = rho * (R_univ * (Y_H2/MW_H2 + Y_O2/MW_O2)) * T_local
                a = sqrt(1.4_dp * max(P, 1.0e-10_dp) / rho)
                max_eigen = sqrt(u**2 + v**2) + a
                if (max_eigen > global_max) global_max = max_eigen
            end do
        end do
        dt_out = CFL_target * min(dx, dy) / max(global_max, 1.0e-10_dp)
    end subroutine compute_time_step

    subroutine compute_residuals()
        integer :: i, j, k
        real(dp) :: Q_L(4), Q_R(4), E_L(4), E_R(4), F_L(4), F_R(4), Flux(4)
        real(dp) :: c_L, c_R, wave_speed, P_L, P_R, rL, rR, Y_L, Y_R, Flux_Species
        real(dp) :: e_int_L, e_int_R, T_L, T_R, Y_H2, Y_O2, R_mix_L, R_mix_R
        
        Res = 0.0_dp

        ! X-Fluxes
        do j = 2, nj-1
            do i = 1, ni-1
                Q_L = Q(i,j,1:4); Q_R = Q(i+1,j,1:4)
                rL = max(Q_L(1), 1.0e-10_dp); rR = max(Q_R(1), 1.0e-10_dp)
                
                ! Real Gas P (Left)
                e_int_L = Q_L(4)/rL - 0.5_dp*(Q_L(2)**2+Q_L(3)**2)/rL**2
                Y_H2 = Q(i,j,5)/rL; Y_O2 = Q(i,j,6)/rL
                call get_temperature_from_energy(e_int_L, Y_H2, Y_O2, 1500.0_dp, T_L)
                R_mix_L = R_univ * (Y_H2/MW_H2 + Y_O2/MW_O2) 
                P_L = rL * R_mix_L * T_L

                ! Real Gas P (Right)
                e_int_R = Q_R(4)/rR - 0.5_dp*(Q_R(2)**2+Q_R(3)**2)/rR**2
                Y_H2 = Q(i+1,j,5)/rR; Y_O2 = Q(i+1,j,6)/rR
                call get_temperature_from_energy(e_int_R, Y_H2, Y_O2, 1500.0_dp, T_R)
                R_mix_R = R_univ * (Y_H2/MW_H2 + Y_O2/MW_O2)
                P_R = rR * R_mix_R * T_R

                ! Manual Flux
                E_L(1)=Q_L(2); E_L(2)=Q_L(2)**2/rL+P_L; E_L(3)=Q_L(2)*Q_L(3)/rL; E_L(4)=(Q_L(4)+P_L)*Q_L(2)/rL
                E_R(1)=Q_R(2); E_R(2)=Q_R(2)**2/rR+P_R; E_R(3)=Q_R(2)*Q_R(3)/rR; E_R(4)=(Q_R(4)+P_R)*Q_R(2)/rR

                c_L = sqrt(1.4_dp*P_L/rL) + abs(Q_L(2)/rL)
                c_R = sqrt(1.4_dp*P_R/rR) + abs(Q_R(2)/rR)
                wave_speed = max(c_L, c_R)
                
                do k = 1, 4
                    Flux(k) = 0.5_dp*(E_L(k)+E_R(k)) - 0.5_dp*wave_speed*(Q_R(k)-Q_L(k))
                end do
                Res(i,j,1:4) = Res(i,j,1:4) - Flux(:)/dx
                Res(i+1,j,1:4) = Res(i+1,j,1:4) + Flux(:)/dx

                ! Species
                do k = 5, 10
                    Y_L = Q(i,j,k); Y_R = Q(i+1,j,k)
                    Flux_Species = 0.5_dp * ( (Y_L * Q_L(2)/rL) + (Y_R * Q_R(2)/rR) ) &
                                 - 0.5_dp * wave_speed * (Y_R - Y_L)
                    Res(i,j,k)   = Res(i,j,k)   - Flux_Species/dx
                    Res(i+1,j,k) = Res(i+1,j,k) + Flux_Species/dx
                end do
            end do
        end do

        ! Y-Fluxes 
        do i = 2, ni-1
            do j = 1, nj-1
                Q_L = Q(i,j,1:4); Q_R = Q(i,j+1,1:4)
                rL = max(Q_L(1), 1.0e-10_dp); rR = max(Q_R(1), 1.0e-10_dp)
                
                ! Use approximate pressure for Y-sweep speed (Optimization)
                P_L = 0.4_dp*(Q_L(4) - 0.5_dp*(Q_L(2)**2+Q_L(3)**2)/rL) 
                P_R = 0.4_dp*(Q_R(4) - 0.5_dp*(Q_R(2)**2+Q_R(3)**2)/rR)
                
                call get_flux_vectors(Q_L, E_L, F_L)
                call get_flux_vectors(Q_R, E_R, F_R)
                
                c_L = sqrt(1.4_dp*P_L/rL) + abs(Q_L(3)/rL)
                c_R = sqrt(1.4_dp*P_R/rR) + abs(Q_R(3)/rR)
                wave_speed = max(c_L, c_R)
                
                do k = 1, 4
                    Flux(k) = 0.5_dp*(F_L(k)+F_R(k)) - 0.5_dp*wave_speed*(Q_R(k)-Q_L(k))
                end do
                Res(i,j,1:4) = Res(i,j,1:4) - Flux(:)/dy
                Res(i,j+1,1:4) = Res(i,j+1,1:4) + Flux(:)/dy

                do k = 5, 10
                    Y_L = Q(i,j,k); Y_R = Q(i,j+1,k)
                    Flux_Species = 0.5_dp * ( (Y_L * Q_L(3)/rL) + (Y_R * Q_R(3)/rR) ) &
                                 - 0.5_dp * wave_speed * (Y_R - Y_L)
                    Res(i,j,k)   = Res(i,j,k)   - Flux_Species/dy
                    Res(i,j+1,k) = Res(i,j+1,k) + Flux_Species/dy
                end do
            end do
        end do
    end subroutine compute_residuals

    subroutine update_solution(dt_in, step_num)
        real(dp), intent(in) :: dt_in
        integer, intent(in) :: step_num
        integer :: i, j, k
        real(dp) :: rho, e_int, u, v, T, Y(6), L2_Norm
        
        Q = Q + dt_in * Res
        
        if (chem_switch == 1) then
            !$OMP PARALLEL DO PRIVATE(i, j, k, rho, u, v, e_int, T, Y) COLLAPSE(2)
            do j = 1, nj
                do i = 1, ni
                    rho = max(Q(i,j,1), 1.0e-10_dp)
                    do k = 1, 6; Y(k) = max(0.0_dp, Q(i,j,4+k) / rho); end do
                    u = Q(i,j,2)/rho; v = Q(i,j,3)/rho
                    e_int = Q(i,j,4)/rho - 0.5_dp*(u**2 + v**2)
                    call get_temperature_from_energy(e_int, Y(1), Y(2), 1500.0_dp, T)
                    
                    call reaction_subcycle(rho, T, Y, dt_in)
                    
                    do k = 1, 6; Q(i,j,4+k) = rho * Y(k); end do
                end do
            end do
            !$OMP END PARALLEL DO
        end if
        
        if (mod(step_num, 10) == 0) then
            L2_Norm = sqrt(sum(Res(:,:,1)**2) / real(ni*nj, dp))
            ! --- FIX: Safe Open/Append logic ---
            open(30, file='history.dat', status='unknown', position='append')
            write(30, '(I8, 1X, ES14.6, 1X, ES14.6)') step_num, dt_in*real(step_num,dp), L2_Norm
            close(30)
        end if
    end subroutine update_solution

    subroutine apply_boundary_conditions()
        integer :: i, j
        real(dp) :: rho_jet, u_jet, v_jet, P_jet, Et_jet, R_jet
        
        do j = 1, nj; Q(1, j, 1:10) = Q_inlet(1:10); Q(ni, j, 1:10) = Q(ni-1, j, 1:10); end do
        do i = 1, ni
            Q(i, nj, 1:10) = Q(i, nj-1, 1:10); Q(i, nj, 3) = 0.0_dp
            if (Y_H2_init < 1.0e-5_dp .and. x(i,1) >= 0.02_dp .and. x(i,1) <= 0.03_dp) then
                P_jet = 2.0_dp * P_init; R_jet = R_univ / MW_H2; rho_jet = P_jet / (R_jet * 300.0_dp)
                u_jet = 0.0_dp; v_jet = 500.0_dp          
                Et_jet = P_jet/((1.41_dp-1.0_dp)*rho_jet) + 0.5_dp*(u_jet**2 + v_jet**2)
                Q(i, 1, 1) = rho_jet; Q(i, 1, 2) = rho_jet * u_jet; Q(i, 1, 3) = rho_jet * v_jet
                Q(i, 1, 4) = rho_jet * Et_jet; Q(i, 1, 5) = rho_jet * 1.0_dp; Q(i, 1, 6:10) = 0.0_dp         
            else
                Q(i, 1, 1:10) = Q(i, 2, 1:10); Q(i, 1, 3) = 0.0_dp
            end if
        end do
    end subroutine apply_boundary_conditions

    subroutine write_output(frame)
        integer, intent(in) :: frame
        integer :: i, j
        real(dp) :: r, ru, rv, rE, uu, vv, P, T, R_mix
        real(dp) :: Y_H2, Y_O2, Y_H2O, Y_OH, Y_O, Y_H, e_int
        character(len=32) :: filename
        
        write(filename, '("flow_", I0.3, ".dat")') frame
        open(20, file=trim(filename), status='unknown')
        write(20, '(A)') "X Y Rho U V P T Y_H2 Y_O2 Y_H2O Y_OH Y_O Y_H Mach"
        
        do j = 1, nj
            do i = 1, ni
                r  = max(Q(i,j,1), 1.0e-10_dp)
                ru = Q(i,j,2); rv = Q(i,j,3); rE = Q(i,j,4)
                Y_H2 = Q(i,j,5)/r; Y_O2 = Q(i,j,6)/r; Y_H2O = Q(i,j,7)/r
                Y_OH = Q(i,j,8)/r; Y_O  = Q(i,j,9)/r; Y_H   = Q(i,j,10)/r
                uu = ru/r; vv = rv/r
                e_int = rE/r - 0.5_dp*(uu**2 + vv**2)
                call get_temperature_from_energy(e_int, Y_H2, Y_O2, 1500.0_dp, T)
                R_mix = R_univ * (Y_H2/MW_H2 + Y_O2/MW_O2 + Y_H2O/MW_H2O + (1.0-Y_H2-Y_O2-Y_H2O)/28.0_dp)
                P = r * R_mix * T
                write(20, '(14(ES14.6, 1X))') x(i,j), y(i,j), r, uu, vv, P, T, &
                                              Y_H2, Y_O2, Y_H2O, Y_OH, Y_O, Y_H, &
                                              sqrt(uu**2+vv**2)/sqrt(1.4_dp*R_mix*T)
            end do
        end do
        close(20)
    end subroutine write_output

end module solver_mod