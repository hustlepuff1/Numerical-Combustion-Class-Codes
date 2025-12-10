module constants_mod
    implicit none
    public
    
    ! --- Precision ---
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    ! --- Physical Constants ---
    real(dp), parameter :: gamma   = 1.4_dp
    real(dp), parameter :: R_univ  = 8314.462618_dp
    real(dp), parameter :: PI      = 3.14159265358979323846_dp
    
    ! --- Molecular Weights (kg/kmol) ---
    ! Consolidating ALL species here for safety
    real(dp), parameter :: MW_H2   = 2.016_dp   
    real(dp), parameter :: MW_O2   = 31.999_dp
    real(dp), parameter :: MW_H2O  = 18.015_dp
    real(dp), parameter :: MW_OH   = 17.007_dp
    real(dp), parameter :: MW_O    = 15.999_dp
    real(dp), parameter :: MW_H    = 1.008_dp
    real(dp), parameter :: MW_N2   = 28.013_dp

end module constants_mod