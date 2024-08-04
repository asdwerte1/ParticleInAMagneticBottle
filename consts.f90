module consts
    implicit none
    integer, parameter :: dp=selected_real_kind(15,300)
    integer, parameter :: i64=selected_int_kind(18)
    real(kind=dp), parameter :: pi=3.14159265358979323846, mu_nought=1.26e-6, epsilon_nought = 8.85e-12, SOL = 299792458
    real(kind=dp), parameter :: boltzmann=1.3806e-23
    complex(kind=dp) :: imaginary = cmplx(0.0_dp, 1.0_dp,dp)
end module consts
