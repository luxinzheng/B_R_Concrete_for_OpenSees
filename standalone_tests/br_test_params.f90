
module br_test_params
  implicit none
  integer, parameter :: dp = selected_real_kind(15,307)
contains

  subroutine set_props_default(props)
    real(dp), intent(out) :: props(37)
    ! Units: MPa, strains dimensionless
    real(dp) :: E0, nu, ft, fc, epsc, sig_u, epsu, stifac, shefac
    real(dp) :: beta_s, gama_s, rkapa_s, alfa_s
    real(dp) :: sp1(6), sp31(6), sp32(6), sp33(6)

    E0  = 30000.0_dp
    nu  = 0.20_dp
    ft  = 2.5_dp
    fc  = -30.0_dp
    epsc = -0.0020_dp
    sig_u = 0.20_dp*fc      ! ultimate compressive stress magnitude (negative)
    epsu  = -0.0060_dp

    stifac = 1.0e-2_dp      ! forumat.f90 enforces >=1e-2
    shefac = 0.50_dp

    beta_s = 0.90_dp
    gama_s = 1.00_dp
    rkapa_s = 0.70_dp
    alfa_s  = 0.12_dp

    ! Generic normalized tables (monotone in compression)
    sp1  = (/0.00_dp, 0.20_dp, 0.40_dp, 0.60_dp, 0.80_dp, 1.00_dp/)
    sp31 = (/1.00_dp, 1.05_dp, 1.10_dp, 1.12_dp, 1.10_dp, 1.00_dp/)
    sp32 = (/1.00_dp, 0.95_dp, 0.90_dp, 0.85_dp, 0.80_dp, 0.75_dp/)
    sp33 = (/1.00_dp, 0.90_dp, 0.82_dp, 0.75_dp, 0.70_dp, 0.65_dp/)

    props = 0.0_dp
    props(1)=E0; props(2)=nu; props(3)=ft; props(4)=fc; props(5)=epsc
    props(6)=sig_u; props(7)=epsu; props(8)=stifac; props(9)=shefac
    props(10)=beta_s; props(11)=gama_s; props(12)=rkapa_s; props(13)=alfa_s
    props(14:19)=sp1
    props(20:25)=sp31
    props(26:31)=sp32
    props(32:37)=sp33
  end subroutine set_props_default

  subroutine set_props_high_strength(props)
    real(dp), intent(out) :: props(37)
    call set_props_default(props)
    props(1)=42000.0_dp   ! E0
    props(3)=3.5_dp       ! ft
    props(4)=-60.0_dp     ! fc
    props(5)=-0.0022_dp   ! epsc
    props(6)=0.15_dp*props(4) ! sig_u
    props(7)=-0.0075_dp   ! epsu
    props(9)=0.35_dp      ! shefac (more degradation)
  end subroutine set_props_high_strength

  subroutine set_props_bathe1979(props)
    real(dp), intent(out) :: props(37)
    ! Bathe & Ramaswamy (1979) Figs. 8-9 material
    ! Units: ksi, in/in
    real(dp) :: sp1(6), sp31(6), sp32(6), sp33(6)

    props = 0.0_dp
    props(1) = 4600.0_dp        ! E0 (ksi)
    props(2) = 0.0_dp           ! nu (zero — pure uniaxial)
    props(3) = 0.45_dp          ! ft (ksi)
    props(4) = -4.45_dp         ! fc (ksi)
    props(5) = -0.00218_dp      ! epsc (in/in)
    props(6) = -3.84_dp         ! sigma_u (ksi)
    props(7) = -0.00515_dp      ! epsu (in/in)
    props(8) = 1.0e-6_dp        ! STIFAC (paper eta_n=1e-6, clamped to 0.01)
    props(9) = 0.50_dp          ! SHEFAC (eta_s)
    props(10) = 0.90_dp         ! beta
    props(11) = 1.00_dp         ! gama
    props(12) = 1.00_dp         ! rkapa (K in paper)
    props(13) = 0.12_dp         ! alfa

    sp1  = (/0.00_dp, 0.20_dp, 0.40_dp, 0.60_dp, 0.80_dp, 1.00_dp/)
    sp31 = (/1.00_dp, 1.05_dp, 1.10_dp, 1.12_dp, 1.10_dp, 1.00_dp/)
    sp32 = (/1.00_dp, 0.95_dp, 0.90_dp, 0.85_dp, 0.80_dp, 0.75_dp/)
    sp33 = (/1.00_dp, 0.90_dp, 0.82_dp, 0.75_dp, 0.70_dp, 0.65_dp/)
    props(14:19) = sp1
    props(20:25) = sp31
    props(26:31) = sp32
    props(32:37) = sp33
  end subroutine set_props_bathe1979

end module br_test_params
