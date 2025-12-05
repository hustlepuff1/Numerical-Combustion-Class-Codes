module state_mod
  implicit none
  private
  public :: NVAR, NSPEC
  public :: prim_to_cons, cons_to_prim, sound_speed
  public :: R_univ, Wspec, species_name
  public :: temperature_from_primitive
  public :: mixture_cp, mixture_cv

  ! ---------------------------
  ! Core flow: still 4 variables
  ! ---------------------------
  integer, parameter :: NVAR  = 4      ! [rho, rho*u, rho*v, rho*E]

  ! ---------------------------
  ! Chemistry / species data
  ! ---------------------------
  integer, parameter :: NSPEC = 7
  ! Index convention (you can change names/order if you like)
  ! 1: H2, 2: O2, 3: H, 4: O, 5: OH, 6: H2O, 7: N2
  character(len=*), parameter :: species_name(NSPEC) = &
       [ character(len=3) :: 'H2 ', 'O2 ', 'H  ', 'O  ', 'OH ', 'H2O', 'N2 ' ]

  ! Molar masses [kg/mol]
  real(8), parameter :: Wspec(NSPEC) = &
       (/ 2.0d-3, 32.0d-3, 1.0d-3, 16.0d-3, 17.0d-3, 18.0d-3, 28.0d-3 /)

  ! Universal gas constant
  real(8), parameter :: R_univ = 8.31446261815324d0  ! [J/(mol K)]

contains

  !======================================================
  ! 4-eq primitive -> conservative (same as before)
  !======================================================
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

  !======================================================
  ! 4-eq conservative -> primitive (same as before)
  !======================================================
  subroutine cons_to_prim(Ucons, gamma, rho, u, v, p)
    real(8), intent(in)  :: Ucons(NVAR), gamma
    real(8), intent(out) :: rho, u, v, p
    real(8) :: E, q2

    rho = Ucons(1)
    u   = Ucons(2) / rho
    v   = Ucons(3) / rho
    E   = Ucons(4) / rho
    q2  = 0.5d0 * (u*u + v*v)
    p   = (gamma - 1.d0) * rho * (E - q2)
  end subroutine cons_to_prim

  !======================================================
  ! sound speed a = sqrt(gamma * p / rho)
  !======================================================
  function sound_speed(Ucons, gamma) result(a)
    real(8), intent(in) :: Ucons(NVAR), gamma
    real(8) :: a
    real(8) :: rho, u, v, p

    call cons_to_prim(Ucons, gamma, rho, u, v, p)
    a = sqrt(gamma * p / rho)
  end function sound_speed

  !======================================================
  ! Helper: T from primitive (rho, p) and gas constant Rgas
  !======================================================
  function temperature_from_primitive(rho, p, Rgas) result(T)
    real(8), intent(in) :: rho, p, Rgas
    real(8) :: T
    T = p / (rho * Rgas)
  end function temperature_from_primitive

  !======================================================
  ! Simple mixture cp, cv (for now: constant gamma mix)
  !  -> Useful if later you want variable gamma
  !======================================================
  subroutine mixture_cp_cv(gamma, Rgas, cp, cv)
    real(8), intent(in)  :: gamma, Rgas
    real(8), intent(out) :: cp, cv

    cv = Rgas / (gamma - 1.d0)
    cp = gamma * cv
  end subroutine mixture_cp_cv

  function mixture_cp(gamma, Rgas) result(cp)
    real(8), intent(in) :: gamma, Rgas
    real(8) :: cp, cv
    call mixture_cp_cv(gamma, Rgas, cp, cv)
  end function mixture_cp

  function mixture_cv(gamma, Rgas) result(cv)
    real(8), intent(in) :: gamma, Rgas
    real(8) :: cp, cv_local
    real(8) :: cv
    call mixture_cp_cv(gamma, Rgas, cp, cv_local)
    cv = cv_local
  end function mixture_cv

end module state_mod
