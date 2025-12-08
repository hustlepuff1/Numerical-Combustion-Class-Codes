module state_mod
  implicit none
  private
  public :: NVAR, NSPEC
  public :: prim_to_cons, cons_to_prim, sound_speed
  public :: R_univ, Wspec, species_name
  public :: temperature_from_primitive
  public :: mixture_cp, mixture_cv
  public :: Yspec, allocate_state, deallocate_state
  public :: get_species, set_species
  ! --- NEW: NASA Poly Publics ---
  public :: init_nasa_coeffs, get_species_H_RT, get_species_S_R, get_species_cp_R

  ! ---------------------------
  ! Core flow + Species
  ! ---------------------------
  ! 8 species, order matches senior's C++ code:
  ! 1: H2, 2: O2, 3: O, 4: H, 5: H2O, 6: HO2, 7: OH, 8: N2
  integer, parameter :: NSPEC = 8

  ! NVAR is now 12 (4 flow vars + 8 species)
  integer, parameter :: NVAR  = 4 + NSPEC

  character(len=*), parameter :: species_name(NSPEC) = &
       [ character(len=4) :: 'H2  ', 'O2  ', 'O   ', 'H   ', &
                             'H2O ', 'HO2 ', 'OH  ', 'N2  ' ]

  ! Molar masses [kg/mol]
  real(8), parameter :: Wspec(NSPEC) = (/ &
       2.0d-3,  & ! H2
      32.0d-3,  & ! O2
      16.0d-3,  & ! O
       1.0d-3,  & ! H
      18.0d-3,  & ! H2O
      33.0d-3,  & ! HO2
      17.0d-3,  & ! OH
      28.0d-3   & ! N2
  /)

  real(8), parameter :: R_univ = 8.31446261815324d0  ! [J/(mol K)]

  ! Legacy pointer for chemistry module
  real(8), allocatable :: Yspec(:,:,:)

  ! --- NASA 7-Coefficient Polynomials Storage ---
  ! 7 coeffs for Low Temp (T < 1000K), 7 for High Temp (T >= 1000K)
  real(8), save :: nasa_lo(7, NSPEC)
  real(8), save :: nasa_hi(7, NSPEC)
  logical, save :: nasa_initialized = .false.

contains

  ! ==================================================================
  !  NASA POLYNOMIAL INITIALIZATION & GETTERS
  ! ==================================================================
  subroutine init_nasa_coeffs()
    if (nasa_initialized) return

    ! -----------------------------------------------------
    ! ORDER: 1:H2, 2:O2, 3:O, 4:H, 5:H2O, 6:HO2, 7:OH, 8:N2
    ! Data from GRI-Mech 3.0 / NASA 7-coeff format
    ! -----------------------------------------------------

    ! --- 1. H2 ---
    nasa_lo(:,1) = (/ 2.34433112E+00, 7.98052075E-03, -1.94781510E-05, &
                      2.01572094E-08, -7.37611761E-12, -9.17935173E+02, 6.83010238E-01 /)
    nasa_hi(:,1) = (/ 3.33727920E+00, -4.94024731E-05, 4.99456778E-07, &
                     -1.79566394E-10, 2.00255376E-14, -9.50158922E+02, -3.20502331E+00 /)

    ! --- 2. O2 ---
    nasa_lo(:,2) = (/ 3.78245636E+00, -2.99673416E-03, 9.84730201E-06, &
                     -9.68129509E-09, 3.24372837E-12, -1.06394356E+03, 3.65767573E+00 /)
    nasa_hi(:,2) = (/ 3.28253784E+00, 1.48308754E-03, -7.57966669E-07, &
                      2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129E+00 /)

    ! --- 3. O ---
    nasa_lo(:,3) = (/ 3.16826710E+00, -3.27931884E-03, 6.64306396E-06, &
                     -6.12806624E-09, 2.11265971E-12, 2.91222592E+04, 2.05193346E+00 /)
    nasa_hi(:,3) = (/ 2.56942078E+00, -8.59741137E-05, 4.19484589E-08, &
                     -1.00177799E-11, 1.22833691E-15, 2.92175791E+04, 4.78433864E+00 /)

    ! --- 4. H ---
    nasa_lo(:,4) = (/ 2.50000000E+00, 7.05332819E-13, -1.99591964E-15, &
                      2.30081632E-18, -9.27732332E-22, 2.54736599E+04, -4.46682853E-01 /)
    nasa_hi(:,4) = (/ 2.50000001E+00, -2.30842973E-11, 1.61561948E-14, &
                     -4.73515235E-18, 4.98197357E-22, 2.54736599E+04, -4.46682914E-01 /)

    ! --- 5. H2O ---
    nasa_lo(:,5) = (/ 4.19864056E+00, -2.03643410E-03, 6.52040211E-06, &
                     -5.48797062E-09, 1.77197817E-12, -3.02937267E+04, -8.49032208E-01 /)
    nasa_hi(:,5) = (/ 3.03399249E+00, 2.17691804E-03, -1.64072518E-07, &
                     -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00 /)

    ! --- 6. HO2 ---
    nasa_lo(:,6) = (/ 4.30179801E+00, -4.74912051E-03, 2.11582891E-05, &
                     -2.42763894E-08, 9.29225124E-12, 2.94808040E+02, 3.71666245E+00 /)
    nasa_hi(:,6) = (/ 4.01721090E+00, 2.23982013E-03, -6.33658150E-07, &
                      1.14246370E-10, -1.07908535E-14, 1.11856713E+02, 3.78510215E+00 /)

    ! --- 7. OH ---
    nasa_lo(:,7) = (/ 3.99201543E+00, -2.40131752E-03, 4.61793841E-06, &
                     -3.88113333E-09, 1.36411470E-12, 3.61508056E+03, -1.03925458E-01 /)
    nasa_hi(:,7) = (/ 3.09288767E+00, 5.48429716E-04, 1.26505228E-07, &
                     -8.79461556E-11, 1.17412376E-14, 3.85865700E+03, 4.47669610E+00 /)

    ! --- 8. N2 ---
    nasa_lo(:,8) = (/ 0.32986770E+01, 0.14082404E-02, -0.39632220E-05, &
                      0.56415150E-08, -0.24448540E-11, -0.10208999E+04, 0.39503720E+01 /)
    nasa_hi(:,8) = (/ 0.29266400E+01, 0.14879768E-02, -0.56847600E-06, &
                      0.10097038E-09, -0.67533510E-14, -0.92279770E+03, 0.59805280E+01 /)

    nasa_initialized = .true.
  end subroutine init_nasa_coeffs

  ! --- Specific Heat Cp/R ---
  function get_species_cp_R(k, T) result(cp_R)
    integer, intent(in) :: k
    real(8), intent(in) :: T
    real(8) :: cp_R, T2, T3, T4, a(7)

    if (.not. nasa_initialized) call init_nasa_coeffs()

    if (T < 1000.d0) then
       a = nasa_lo(:, k)
    else
       a = nasa_hi(:, k)
    end if

    T2 = T*T; T3 = T2*T; T4 = T3*T
    cp_R = a(1) + a(2)*T + a(3)*T2 + a(4)*T3 + a(5)*T4
  end function get_species_cp_R

  ! --- Enthalpy H/RT ---
  function get_species_H_RT(k, T) result(h_RT)
    integer, intent(in) :: k
    real(8), intent(in) :: T
    real(8) :: h_RT, T2, T3, T4, a(7)

    if (.not. nasa_initialized) call init_nasa_coeffs()

    if (T < 1000.d0) then
       a = nasa_lo(:, k)
    else
       a = nasa_hi(:, k)
    end if

    T2 = T*T; T3 = T2*T; T4 = T3*T
    h_RT = a(1) + (a(2)*T)/2.d0 + (a(3)*T2)/3.d0 + &
           (a(4)*T3)/4.d0 + (a(5)*T4)/5.d0 + a(6)/T
  end function get_species_H_RT

  ! --- Entropy S/R ---
  function get_species_S_R(k, T) result(s_R)
    integer, intent(in) :: k
    real(8), intent(in) :: T
    real(8) :: s_R, T2, T3, T4, a(7)

    if (.not. nasa_initialized) call init_nasa_coeffs()

    if (T < 1000.d0) then
       a = nasa_lo(:, k)
    else
       a = nasa_hi(:, k)
    end if

    T2 = T*T; T3 = T2*T; T4 = T3*T
    s_R = a(1)*log(T) + a(2)*T + (a(3)*T2)/2.d0 + &
          (a(4)*T3)/3.d0 + (a(5)*T4)/4.d0 + a(7)
  end function get_species_S_R

  ! ==================================================================
  !  EXISTING UTILITIES
  ! ==================================================================
  subroutine allocate_state(U, ni, nj)
    real(8), allocatable, intent(out) :: U(:,:,:)
    integer, intent(in) :: ni, nj
    allocate(U(NVAR, ni, nj))
    allocate(Yspec(NSPEC, ni, nj))
  end subroutine allocate_state

  subroutine deallocate_state(U)
    real(8), allocatable, intent(inout) :: U(:,:,:)
    if (allocated(U))     deallocate(U)
    if (allocated(Yspec)) deallocate(Yspec)
  end subroutine deallocate_state

  subroutine get_species(i, j, Yout)
    integer, intent(in)  :: i, j
    real(8), intent(out) :: Yout(NSPEC)
    Yout = Yspec(:, i, j)
  end subroutine get_species

  subroutine set_species(i, j, Yin)
    integer, intent(in)  :: i, j
    real(8), intent(in)  :: Yin(NSPEC)
    Yspec(:, i, j) = Yin
  end subroutine set_species

  subroutine prim_to_cons(rho, u, v, p, gamma, Ucons)
    real(8), intent(in)  :: rho, u, v, p, gamma
    real(8), intent(out) :: Ucons(NVAR)
    real(8) :: E
    E = p / ((gamma - 1.d0) * rho) + 0.5d0 * (u*u + v*v)
    Ucons(1) = rho
    Ucons(2) = rho * u
    Ucons(3) = rho * v
    Ucons(4) = rho * E
  end subroutine prim_to_cons

  subroutine cons_to_prim(Ucons, gamma, rho, u, v, p)
    real(8), intent(in)  :: Ucons(NVAR), gamma
    real(8), intent(out) :: rho, u, v, p
    real(8) :: E, q2
    rho = Ucons(1)
    if (rho < 1.0d-30) rho = 1.0d-30
    u   = Ucons(2) / rho
    v   = Ucons(3) / rho
    E   = Ucons(4) / rho
    q2  = 0.5d0 * (u*u + v*v)
    p   = (gamma - 1.d0) * rho * (E - q2)
  end subroutine cons_to_prim

  function sound_speed(Ucons, gamma) result(a)
    real(8), intent(in) :: Ucons(NVAR), gamma
    real(8) :: a
    real(8) :: rho, u, v, p
    call cons_to_prim(Ucons, gamma, rho, u, v, p)
    if (p < 1.0d-9) p = 1.0d-9
    a = sqrt(gamma * p / rho)
  end function sound_speed

  function temperature_from_primitive(rho, p, Rgas) result(T)
    real(8), intent(in) :: rho, p, Rgas
    real(8) :: T
    T = p / (rho * Rgas)
  end function temperature_from_primitive

  subroutine mixture_cp_cv(gamma, Rgas, cp, cv)
    real(8), intent(in)  :: gamma, Rgas
    real(8), intent(out) :: cp, cv
    cv = Rgas / (gamma - 1.d0)
    cp = gamma * cv
  end subroutine mixture_cp_cv

  function mixture_cp(gamma, Rgas) result(cp)
    real(8), intent(in) :: gamma, Rgas
    real(8) :: cp, cv_local
    call mixture_cp_cv(gamma, Rgas, cp, cv_local)
  end function mixture_cp

  function mixture_cv(gamma, Rgas) result(cv)
    real(8), intent(in) :: gamma, Rgas
    real(8) :: cv, cp_local
    call mixture_cp_cv(gamma, Rgas, cp_local, cv)
  end function mixture_cv

end module state_mod