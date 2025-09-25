! --------------------------------------------------------------------
! Program: 1d Advection solver
! Purpose: Solves u_t + u_x = 0 using various finite difference schemes.
! --------------------------------------------------------------------
program hw2_advectionSolver
    implicit none

    real, parameter :: x_min = -1.0
    real, parameter :: x_max = 3.0
    real, parameter :: t_final = 2.4
    real, parameter :: dx = 1.0 / 40.0
    real, parameter :: cfl = 0.8
    real, parameter :: pi = acos(-1.0)

    real    :: dt
    integer :: nx
    integer :: nt

    real, allocatable :: x(:)
    real, allocatable :: u_old(:)
    real, allocatable :: u_new(:)
    real, allocatable :: u_exact(:)
    real, allocatable :: u_temp(:) ! for leapfrog scheme, u n+1

    ! Loops
    integer :: i
    integer :: n

    ! File handling
    integer :: file_unit

    dt = cfl * dx
    nx = int((x_max - x_min) / dx) + 1
    nt = int(t_final / dt)

    write(*,*) "--- Simulation parameters ---"
    write(*,*) "Spatial points (nx): ", nx
    write(*,*) "Time steps   (nt):   ", nt
    write(*,*) "dx = ", dx
    write(*,*) "dt = ", dt
    write(*,*) "---------------------------"

    allocate(x(nx), u_old(nx), u_new(nx), u_exact(nx), u_temp(nx))

    do i = 1, nx
        x(i) = x_min + real(i-1) * dx

        if (abs(x(i)) <= 0.5) then
            u_old(i) = cos(pi * x(i))**2
        else
            u_old(i) = 0.0
        end if

    end do
    
    write(*,*) "Spatial grid and initial condition created."

    ! ! Main time-stepping loop
    ! do n = 1, nt
    !     ! call ftbsScheme(u_old, u_new, nx, cfl)
    !     ! call ftcsScheme(u_old, u_new, nx, cfl)
    !     ! call ftfsScheme(u_old, u_new, nx, cfl)
    !     call laxFriedrichsScheme(u_old, u_new, nx, cfl)
    !     u_old = u_new
    ! end do

    ! ## main time-stepping loop for leapfrog
    call ftcsScheme(u_old, u_new, nx, cfl)

    do n = 2, nt
        call leapfrogScheme(u_old, u_new, u_temp, nx, cfl)
        u_old = u_new
        u_new = u_temp
    end do
    ! ## end of main time-stepping loop for leapfrog


    write(*,*) "time-stepping complete. final time t = ", real(nt)*dt

    ! Exact solution
    do i = 1, nx
        if (abs(x(i) - t_final) <= 0.5) then
            u_exact(i) = cos(pi * (x(i) - t_final))**2
        else
            u_exact(i) = 0.0
        end if
    end do




    ! open(newunit=file_unit, file='results_ftbs.dat', status='replace')
    ! open(newunit=file_unit, file='results_ftcs.dat', status='replace')
    ! open(newunit=file_unit, file='results_ftfs.dat', status='replace')
    ! open(newunit=file_unit, file='results_laxfriedrichs.dat', status='replace')
    open(newunit=file_unit, file='results_leapfrog.dat', status='replace')

    write(file_unit, '(A8, A15, A15, A15)') "x", "u_numerical", "u_exact", "error"
    do i = 1, nx
        write(file_unit, '(4F15.7)') x(i), u_new(i), u_exact(i), u_new(i)-u_exact(i)
    end do
    close(file_unit)
    ! write(*,*) "results written to results_ftbs.dat"
    ! write(*,*) "results written to results_ftcs.dat"
    ! write(*,*) "results written to results_ftfs.dat"
    ! write(*,*) "results written to results_laxfriedrichs.dat"
    write(*,*) "results written to results_leapfrog.dat"


    deallocate(x, u_old, u_new, u_exact)

contains

subroutine ftbsScheme(u_old, u_new, nx, cfl)
    
    integer, intent(in)           :: nx
    real, intent(in)              :: cfl
    real, dimension(nx), intent(in)  :: u_old
    real, dimension(nx), intent(out) :: u_new

    integer :: i

    u_new(1) = 0.0

    do i = 2, nx - 1 
        u_new(i) = u_old(i) - cfl * (u_old(i) - u_old(i-1))
    end do

    u_new(nx) = u_new(nx - 1)

end subroutine ftbsScheme

subroutine ftcsScheme(u_old, u_new, nx, cfl)

    integer, intent(in)           :: nx
    real, intent(in)              :: cfl
    real, dimension(nx), intent(in)  :: u_old
    real, dimension(nx), intent(out) :: u_new

    integer :: i

    u_new(1) = 0.0

    do i = 2, nx - 1 
        u_new(i) = u_old(i) - (cfl / 2.0) * (u_old(i+1) - u_old(i-1))
    end do

    u_new(nx) = u_new(nx - 1)

end subroutine ftcsScheme

subroutine ftfsScheme(u_old, u_new, nx, cfl)

    integer, intent(in)           :: nx
    real, intent(in)              :: cfl
    real, dimension(nx), intent(in)  :: u_old
    real, dimension(nx), intent(out) :: u_new

    integer :: i

    u_new(1) = 0.0

    do i = 2, nx - 1 
        u_new(i) = u_old(i) - cfl * (u_old(i+1) - u_old(i))
    end do

    u_new(nx) = u_new(nx - 1)

end subroutine ftfsScheme

subroutine laxFriedrichsScheme(u_old, u_new, nx, cfl)

    integer, intent(in)           :: nx
    real, intent(in)              :: cfl
    real, dimension(nx), intent(in)  :: u_old
    real, dimension(nx), intent(out) :: u_new

    integer :: i

    u_new(1) = 0.0

    do i = 2, nx - 1 
        u_new(i) = 0.5 * (u_old(i+1) + u_old(i-1)) - &
                   (cfl / 2.0) * (u_old(i+1) - u_old(i-1))
    end do

    u_new(nx) = u_new(nx - 1)

end subroutine laxFriedrichsScheme

subroutine leapfrogScheme(u_oldest, u_current, u_new, nx, cfl)

    integer, intent(in)           :: nx
    real, intent(in)              :: cfl
    real, dimension(nx), intent(in)  :: u_oldest, u_current
    real, dimension(nx), intent(out) :: u_new

    integer :: i

    u_new(1) = 0.0

    do i = 2, nx - 1 
        u_new(i) = u_oldest(i) - cfl * (u_current(i+1) - u_current(i-1))
    end do

    u_new(nx) = u_new(nx - 1)

end subroutine leapfrogScheme



end program hw2_advectionSolver