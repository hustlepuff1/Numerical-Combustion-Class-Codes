module state_mod
  implicit none
  private
  public :: NVAR, NSPEC
  public :: prim_to_cons, cons_to_prim, sound_speed
  public :: R_univ, Wspec, species_name
  public :: temperature_from_primitive
  public :: mixture_cp, mixture_cv
  ! NEW:
  public :: Yspec, allocate_state, deallocate_state
  public :: get_species, set_species

  ! ---------------------------
  ! Core flow: 4 variables
  ! ---------------------------
  integer, parameter :: NVAR  = 4      ! [rho, rho*u, rho*v, rho*E]

  ! ---------------------------
  ! Chemistry / species data
  ! ---------------------------
  ! 8 species, order matches senior's C++ code:
  ! 1: H2, 2: O2, 3: O, 4: H, 5: H2O, 6: HO2, 7: OH, 8: N2
  integer, parameter :: NSPEC = 8

  character(len=*), parameter :: species_name(NSPEC) = &
       [ character(len=4) :: 'H2  ', 'O2  ', 'O   ', 'H   ', &
                             'H2O ', 'HO2 ', 'OH  ', 'N2  ' ]

  ! Molar masses [kg/mol]
  ! (kg/mol, so 2e-3 = 2 kg/kmol, etc.)
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

  ! Universal gas constant
  real(8), parameter :: R_univ = 8.31446261815324d0  ! [J/(mol K)]

  ! Species mass fractions on the CFD grid:
  !   Yspec(k, i, j), k = 1..NSPEC
  real(8), allocatable :: Yspec(:,:,:)

contains

  !======================================================
  ! Allocate / deallocate state + species fields
  !======================================================
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

  !======================================================
  ! Species getters / setters for BC & chemistry
  !======================================================
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

  !======================================================
  ! 4-eq primitive -> conservative
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
  ! 4-eq conservative -> primitive
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
  ! T from primitive (rho, p) and gas constant Rgas
  !======================================================
  function temperature_from_primitive(rho, p, Rgas) result(T)
    real(8), intent(in) :: rho, p, Rgas
    real(8) :: T
    T = p / (rho * Rgas)
  end function temperature_from_primitive

  !======================================================
  ! Simple mixture cp, cv (constant gamma mix for now)
  !======================================================
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
