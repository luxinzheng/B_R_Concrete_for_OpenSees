!===============================================================================
!  ADINA Concrete Plane-Stress Constitutive Model
!  Ported to OpenSees PlaneStress UserMaterial Interface
!
!  Based on: Bathe & Ramaswamy (1979), ADINA MODEL=5
!  Faithfully follows the original all_adina.for flow with one addition:
!    - Linear tension softening in crack direction (not in original ADINA)
!
!  Features:
!    - Saenz-type uniaxial compressive stress-strain curve
!    - Biaxial failure envelope (Kupfer-type)
!    - Fixed smeared crack model (0/1/2/3-dir cracks + crushing)
!    - Shear retention after cracking
!    - Hypoelastic stress update (incremental)
!    - Tension softening: sigma_t = ft*(1 - eps_open/eps_tu), linear to zero
!
!  Props layout (nprops=37):
!    1: E0        Young's modulus
!    2: VNU       Poisson's ratio
!    3: SIGMAT    Uniaxial tensile strength (> 0)
!    4: SIGMAC    Uniaxial compressive strength (< 0)
!    5: EPSC      Strain at peak compressive stress (< 0)
!    6: SIGMAU    Ultimate compressive stress (<= 0)
!    7: EPSU      Ultimate compressive strain (< 0)
!    8: STIFAC    Post-crack stiffness factor (e.g. 1e-4)
!    9: SHEFAC    Shear retention factor (0~1, e.g. 0.5)
!   10: BETA_S    Biaxial envelope transition parameter (0~1)
!   11: GAMA_S    Strain ratio parameter (default 1.0)
!   12: RKAPA_S   Anisotropic transition ratio (default 0.7)
!   13: ALFA_S    Drucker-Prager alpha (default 0.12)
!   14-19: SP1(6)    Biaxial envelope sigma1/fc table
!   20-25: SP31(6)   Biaxial envelope failure ratio table 1
!   26-31: SP32(6)   Biaxial envelope failure ratio table 2
!   32-37: SP33(6)   Biaxial envelope failure ratio table 3
!
!  State variables (nstatevs >= 8):
!    1:   EVMAX   Maximum equivalent stress measure
!    2:   ANGLE   Crack angle encoding (1000 = no crack)
!    3-5: CRKSTR(3) Crack strains / crush parameters
!    6:   PGRAV   Failure flag (0=normal, 100=crushed)
!    7:   init_flag (0=not initialized, 1=initialized)
!    8:   (spare)
!===============================================================================
module adina_concrete_mod
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: PI_VAL  = 3.141592653589793238d0
    real(dp), parameter :: DEG2RAD = 3.141592653589793238d0 / 180.0d0
    real(dp), parameter :: RAD2DEG = 180.0d0 / 3.141592653589793238d0
    real(dp), parameter :: ANG_NOCRACK = 1.0d3
    real(dp), parameter :: PGRAV_CRUSH = 1.0d2
    
    ! Module-level working variables (replaces COMMON blocks)
    real(dp) :: m_STIFAC, m_SHEFAC, m_SIGMAT, m_SIGMAC, m_EPSC
    real(dp) :: m_SIGMAU, m_EPSU, m_RAM5, m_RBM5, m_RCM5
    real(dp) :: m_BETA, m_GAMA, m_RKAPA, m_ALFA
    real(dp) :: m_SIGP(4), m_TEP(4), m_EP(4), m_YP(3)
    real(dp) :: m_E, m_VNU, m_RK, m_G
    real(dp) :: m_E12, m_E14, m_E24
    real(dp) :: m_EPSCP, m_SIGCP, m_FALSTR
    integer  :: m_MOD45, m_ILFSET
    
contains

!-----------------------------------------------------------------------
! compute_principal_2d: principal stresses and angle from 2D stress
!   STR = [sig_xx, sig_yy, tau_xy, sig_zz]
!   Sets m_SIGP(1:4), returns ANG in [0, 180)
!-----------------------------------------------------------------------
subroutine compute_principal_2d(STR, ANG)
    real(dp), intent(in)  :: STR(4)
    real(dp), intent(out) :: ANG
    real(dp) :: AA, BB, CC, DUM
    
    AA = (STR(1) + STR(2)) * 0.5d0
    BB = (STR(1) - STR(2)) * 0.5d0
    CC = sqrt(BB*BB + STR(3)*STR(3))
    m_SIGP(1) = AA + CC
    m_SIGP(2) = AA - CC
    m_SIGP(3) = 0.0d0
    m_SIGP(4) = STR(4)
    
    ANG = 45.0d0
    if (STR(3) == 0.0d0) ANG = 1.0d-4
    if (abs(BB) < 1.0d-7) return
    DUM = abs(STR(3) / BB)
    ANG = RAD2DEG * atan(DUM)
    
    if (BB > 0.0d0 .and. STR(3) > 0.0d0) then
        ! quadrant 1
    else if (BB < 0.0d0 .and. STR(3) > 0.0d0) then
        ANG = 180.0d0 - ANG
    else if (BB < 0.0d0 .and. STR(3) <= 0.0d0) then
        ANG = 180.0d0 + ANG
    else if (BB > 0.0d0 .and. STR(3) <= 0.0d0) then
        ANG = 360.0d0 - ANG
    end if
    ANG = ANG / 2.0d0
    if (ANG >= 180.0d0) ANG = ANG - 180.0d0
end subroutine

!-----------------------------------------------------------------------
! rotate_strain_to_principal: EPS → principal coords using simple angle
!-----------------------------------------------------------------------
subroutine rotate_strain_to_principal(EPS, ANG, EPSL)
    real(dp), intent(in)  :: EPS(4), ANG
    real(dp), intent(out) :: EPSL(4)
    real(dp) :: GAM, SG, CG
    
    GAM = ANG * DEG2RAD
    SG = sin(GAM); CG = cos(GAM)
    EPSL(1) = EPS(1)*CG*CG + EPS(2)*SG*SG + EPS(3)*SG*CG
    EPSL(2) = EPS(1)*SG*SG + EPS(2)*CG*CG - EPS(3)*SG*CG
    EPSL(3) = 0.0d0
    EPSL(4) = EPS(4)
end subroutine

!-----------------------------------------------------------------------
! rotate_strain_crack: strain → crack coords using ADINA angle encoding
!   Uses double-angle formulas matching CRAKID in all_adina.for
!-----------------------------------------------------------------------
subroutine rotate_strain_crack(EPS, ANGLE, EPSL)
    real(dp), intent(in)  :: EPS(4), ANGLE
    real(dp), intent(out) :: EPSL(4)
    real(dp) :: TANG, GAM, SG, CG, D11, D12
    
    TANG = ANGLE
    if (TANG < -541.0d0) TANG = TANG + 722.0d0
    if (TANG < -180.0d0) TANG = TANG + 361.0d0
    if (TANG >= 180.0d0) TANG = TANG - 180.0d0
    GAM = 2.0d0 * abs(TANG) * DEG2RAD
    SG = sin(GAM); CG = cos(GAM)
    
    D11 = (EPS(1) + EPS(2)) * 0.5d0
    D12 = (EPS(1) - EPS(2)) * 0.5d0
    EPSL(1) = D11 + D12*CG + EPS(3)*SG*0.5d0
    EPSL(2) = D11 - D12*CG - EPS(3)*SG*0.5d0
    EPSL(3) = -D12*SG + EPS(3)*CG
    EPSL(4) = EPS(4)
end subroutine

!-----------------------------------------------------------------------
! rotate_stress_crack: stress → crack coords using ADINA angle encoding
!-----------------------------------------------------------------------
subroutine rotate_stress_crack(STR, ANGLE, STRP)
    real(dp), intent(in)  :: STR(4), ANGLE
    real(dp), intent(out) :: STRP(4)
    real(dp) :: TANG, GAM, SG, CG, R11, R12
    
    TANG = ANGLE
    if (TANG < -541.0d0) TANG = TANG + 722.0d0
    if (TANG < -180.0d0) TANG = TANG + 361.0d0
    if (TANG >= 180.0d0) TANG = TANG - 180.0d0
    GAM = 2.0d0 * abs(TANG) * DEG2RAD
    SG = sin(GAM); CG = cos(GAM)
    
    R11 = (STR(1) + STR(2)) * 0.5d0
    R12 = (STR(1) - STR(2)) * 0.5d0
    STRP(1) = R11 + R12*CG + STR(3)*SG
    STRP(2) = R11 - R12*CG - STR(3)*SG
    STRP(3) = -R12*SG + STR(3)*CG
    STRP(4) = STR(4)
end subroutine

!-----------------------------------------------------------------------
! biaxial_envelope:
!   Compute biaxial failure envelope: m_SIGCP, m_EPSCP, m_FALSTR
!   P1 >= P2 >= P3 are sorted principal stresses
!-----------------------------------------------------------------------
subroutine biaxial_envelope(P1, P2, P3, SP1_T, SP31_T, SP32_T, SP33_T)
    real(dp), intent(in) :: P1, P2, P3
    real(dp), intent(in) :: SP1_T(6), SP31_T(6), SP32_T(6), SP33_T(6)
    real(dp) :: TEMP, RATIO, DSP, DSPI, FRAC, SLOPE
    real(dp) :: SP31I, SP32I, SP33I
    integer  :: I, J
    
    m_SIGCP = m_SIGMAC
    m_EPSCP = m_EPSC
    m_FALSTR = m_SIGMAT
    
    if (P3 >= 0.0d0) return
    
    if (P1 >= 0.0d0 .and. P2 >= 0.0d0) then
        m_FALSTR = m_SIGMAT * (1.0d0 - P3 / m_SIGMAC)
        if (m_FALSTR < 0.01d0 * m_SIGMAT) m_FALSTR = 0.01d0 * m_SIGMAT
        return
    end if
    
    if (P1 < 0.0d0) then
        TEMP = P1 / m_SIGMAC
        do I = 2, 6
            J = I - 1
            if (TEMP < SP1_T(I)) exit
        end do
        if (I > 6) I = 6
        J = I - 1
        DSP  = SP1_T(I) - SP1_T(J)
        DSPI = TEMP - SP1_T(J)
        if (abs(DSP) < 1.0d-15) then
            FRAC = 0.0d0
        else
            FRAC = DSPI / DSP
        end if
        SP31I = SP31_T(J) + FRAC * (SP31_T(I) - SP31_T(J))
        SP32I = SP32_T(J) + FRAC * (SP32_T(I) - SP32_T(J))
        SP33I = SP33_T(J) + FRAC * (SP33_T(I) - SP33_T(J))
    else
        SP31I = SP31_T(1)
        SP32I = SP32_T(1)
        SP33I = SP33_T(1)
        TEMP  = 0.0d0
    end if
    
    RATIO = P2 / m_SIGMAC
    if (RATIO > m_BETA * SP32I) then
        if (abs(SP33I - m_BETA * SP32I) < 1.0d-15) then
            SLOPE = 0.0d0
        else
            SLOPE = (SP33I - SP32I) / (SP33I - m_BETA * SP32I)
        end if
        m_SIGCP = SP32I * m_SIGMAC + SLOPE * (P2 - m_BETA * SP32I * m_SIGMAC)
    else
        if (abs(m_BETA * SP32I - TEMP) < 1.0d-15) then
            SLOPE = 0.0d0
        else
            SLOPE = (SP32I - SP31I) / (m_BETA * SP32I - TEMP)
        end if
        m_SIGCP = SP31I * m_SIGMAC + SLOPE * (P2 - P1)
        if (P1 > 0.0d0) m_SIGCP = SP31I * m_SIGMAC + SLOPE * P2
    end if
    
    if (m_SIGCP < m_SIGMAC) then
        m_EPSCP = m_GAMA * m_EPSC * m_SIGCP / m_SIGMAC
    end if
    
    if (P1 >= 0.0d0) then
        m_FALSTR = m_SIGMAT * (1.0d0 - P2/m_SIGCP) * (1.0d0 - P3/m_SIGCP)
        if (m_FALSTR < 0.001d0 * m_SIGMAT) m_FALSTR = m_SIGMAT
    else
        m_FALSTR = m_SIGCP
    end if
end subroutine

!-----------------------------------------------------------------------
! compute_saenz_params: Saenz curve shape parameters RAM5, RBM5, RCM5
!-----------------------------------------------------------------------
subroutine compute_saenz_params()
    real(dp) :: RP, ES, EU
    
    RP = m_EPSU / m_EPSC
    ES = m_SIGCP / m_EPSCP
    if (abs(m_EPSU * m_EPSCP / m_EPSC) > 1.0d-30) then
        EU = (m_SIGMAU * m_SIGCP / m_SIGMAC) / (m_EPSU * m_EPSCP / m_EPSC)
    else
        EU = m_E
    end if
    
    if (abs(RP * (RP - 1.0d0)**2) < 1.0d-30) then
        m_RAM5 = 0.0d0
    else
        m_RAM5 = (m_E/EU + (RP-2.0d0)*RP*RP*m_E/ES - (2.0d0*RP+1.0d0)*(RP-1.0d0)**2) &
                 / (RP * (RP - 1.0d0)**2)
    end if
    m_RBM5 = 2.0d0 * m_E / ES - 3.0d0 - 2.0d0 * m_RAM5
    m_RCM5 = 2.0d0 - m_E / ES + m_RAM5
end subroutine

!-----------------------------------------------------------------------
! saenz_tangent: tangent modulus at normalized strain DE = eps/epscp
!-----------------------------------------------------------------------
function saenz_tangent(DE) result(TY)
    real(dp), intent(in) :: DE
    real(dp) :: TY, denom
    
    if (DE <= 0.0d0) then
        TY = m_E
        return
    end if
    denom = 1.0d0 + m_RAM5*DE + m_RBM5*DE*DE + m_RCM5*DE**3
    if (abs(denom) < 1.0d-30) then
        TY = 0.0d0
        return
    end if
    TY = m_E * (1.0d0 - m_RBM5*DE*DE - 2.0d0*m_RCM5*DE**3) / (denom*denom)
end function

!-----------------------------------------------------------------------
! saenz_stress: stress on Saenz curve at strain eps_val
!-----------------------------------------------------------------------
function saenz_stress(eps_val) result(sigma)
    real(dp), intent(in) :: eps_val
    real(dp) :: sigma, DE, denom
    
    DE = eps_val / m_EPSCP
    if (DE <= 0.0d0) then
        sigma = m_E * eps_val
        return
    end if
    denom = 1.0d0 + m_RAM5*DE + m_RBM5*DE*DE + m_RCM5*DE**3
    if (abs(denom) < 1.0d-30) then
        sigma = 0.0d0
        return
    end if
    sigma = m_E * eps_val / denom
end function

!-----------------------------------------------------------------------
! compute_tangent_moduli:
!   Gauss quadrature over [TEP, EP] for effective tangent modulus
!   NUMINT=3 for stress update, NUMINT=1 for stiffness (point tangent)
!   Matches all_adina.for labels 7-8 exactly
!-----------------------------------------------------------------------
subroutine compute_tangent_moduli(NUMINT, is_crushed)
    integer, intent(in)  :: NUMINT
    logical, intent(in)  :: is_crushed
    
    real(dp), parameter :: XG3(3) = (/ -0.774596669241483d0, 0.0d0, 0.774596669241483d0 /)
    real(dp), parameter :: WG3(3) = (/  0.555555555555556d0, 0.888888888888889d0, 0.555555555555556d0 /)
    real(dp), parameter :: XG1(1) = (/ 0.0d0 /)
    real(dp), parameter :: WG1(1) = (/ 2.0d0 /)
    
    real(dp) :: DENM, DE, TY, E1
    integer  :: J, I, L, N1, N2
    
    call compute_saenz_params()
    
    N1 = 1; N2 = 4
    if (is_crushed) then
        N1 = 4
        if (m_SIGP(2) < m_SIGP(4)) N1 = 2
        N2 = N1
    end if
    
    DENM = abs(m_SIGP(1)) + abs(m_SIGP(2)) + abs(m_SIGP(4))
    if (DENM <= 0.00001d0 * m_SIGMAT) then
        m_MOD45 = 1
        return
    end if
    
    do J = N1, N2
        if (J == 3) cycle
        I = J
        if (J == 4) I = 3
        m_YP(I) = m_E
        
        if (m_SIGP(J) >= 0.001d0 * m_SIGCP) cycle
        
        m_YP(I) = 0.0d0
        do L = 1, NUMINT
            if (NUMINT == 3) then
                E1 = XG3(L)
            else
                E1 = XG1(1)
            end if
            DE = m_TEP(J) + (1.0d0 + E1) * (m_EP(J) - m_TEP(J)) / 2.0d0
            DE = DE / m_EPSCP
            TY = m_E
            if (DE > 0.0d0) TY = saenz_tangent(DE)
            if (NUMINT == 3) then
                m_YP(I) = m_YP(I) + 0.5d0 * WG3(L) * TY
            else
                m_YP(I) = m_YP(I) + 0.5d0 * WG1(1) * TY
            end if
        end do
    end do
    
    if (is_crushed) then
        m_E = m_YP(I)
        if (m_E > 0.0d0) m_E = 0.0d0
        m_MOD45 = 1
        return
    end if
    
    m_MOD45 = 1
    DE = m_RKAPA * m_SIGCP
    if (m_SIGP(1) < DE .or. m_SIGP(2) < DE .or. m_SIGP(4) < DE) m_MOD45 = 2
    
    m_E = (abs(m_SIGP(1))*m_YP(1) + abs(m_SIGP(2))*m_YP(2) + abs(m_SIGP(4))*m_YP(3)) / DENM
    
    if (m_MOD45 == 1) return
    
    m_E12 = m_E
    DENM = abs(m_SIGP(1)) + abs(m_SIGP(2))
    if (DENM > 0.0d0) m_E12 = (abs(m_SIGP(1))*m_YP(1) + abs(m_SIGP(2))*m_YP(2)) / DENM
    m_E14 = m_E
    DENM = abs(m_SIGP(1)) + abs(m_SIGP(4))
    if (DENM > 0.0d0) m_E14 = (abs(m_SIGP(1))*m_YP(1) + abs(m_SIGP(4))*m_YP(3)) / DENM
    m_E24 = m_E
    DENM = abs(m_SIGP(2)) + abs(m_SIGP(4))
    if (DENM > 0.0d0) m_E24 = (abs(m_SIGP(2))*m_YP(2) + abs(m_SIGP(4))*m_YP(3)) / DENM
end subroutine

!-----------------------------------------------------------------------
! build_constitutive_4x4:
!   Build 4x4 constitutive matrix in crack/principal coords.
!   NUMCRK encoding matches all_adina.for DCRACK exactly:
!     0 = uncracked, 1 = crack dir1, 2 = crack dir1+dir2,
!     3 = crushed,   4 = crack dir4(z), 5 = crack dir1+dir4,
!     6 = crack dir1+dir2+dir4
!-----------------------------------------------------------------------
subroutine build_constitutive_4x4(C, NUMCRK)
    real(dp), intent(out) :: C(4,4)
    integer,  intent(in)  :: NUMCRK
    real(dp) :: DUM, A1, B1, C1, D1
    real(dp) :: A2, B2, C2, D2, R2, S2, T2, W2, Z2
    integer  :: I, J
    
    C = 0.0d0
    
    if (m_MOD45 == 2) then
        A1 = (1.0d0 + m_VNU) * (1.0d0 - 2.0d0*m_VNU)
        B1 = (1.0d0 - m_VNU) / A1
        D1 = m_VNU / A1
        C1 = 1.0d0 / (2.0d0 * (1.0d0 + m_VNU))
        C(1,1) = B1 * m_YP(1)
        C(1,2) = m_E12 * D1
        C(1,4) = m_E14 * D1
        C(2,2) = B1 * m_YP(2)
        C(2,4) = m_E24 * D1
        C(3,3) = C1 * m_E12
        C(4,4) = B1 * m_YP(3)
    else
        DUM = 2.0d0 * m_G / 3.0d0
        A1 = m_RK + 2.0d0 * DUM
        B1 = m_RK - DUM
        C(1,1) = A1
        C(1,2) = B1
        C(1,4) = B1
        C(2,2) = A1
        C(2,4) = B1
        C(3,3) = m_G
        C(4,4) = A1
    end if
    
    C(2,1) = C(1,2)
    C(3,1) = 0.0d0; C(3,2) = 0.0d0; C(3,4) = 0.0d0
    C(4,1) = C(1,4); C(4,2) = C(2,4); C(4,3) = 0.0d0
    
    if (NUMCRK == 0 .or. NUMCRK == 3) return
    
    C(3,3) = C(3,3) * m_SHEFAC
    
    ! Compute condensed 2D stiffnesses (for remaining intact directions)
    if (m_MOD45 == 2) then
        A1 = (1.0d0 + m_VNU) * (1.0d0 - 2.0d0*m_VNU)
        D1 = m_VNU / A1
        B1 = (1.0d0 - m_VNU) / A1
        C1 = 1.0d0 / (2.0d0 * (1.0d0 + m_VNU))
        A2 = m_YP(2) / (1.0d0 - m_VNU*m_VNU)
        B2 = m_VNU * m_E24 / (1.0d0 - m_VNU*m_VNU)
        C2 = m_YP(3) / (1.0d0 - m_VNU*m_VNU)
        D2 = m_YP(3)
        R2 = m_YP(1) / (1.0d0 - m_VNU*m_VNU)
        S2 = m_VNU * m_E12 / (1.0d0 - m_VNU*m_VNU)
        T2 = A2
        W2 = C1 * m_E12
        Z2 = m_E12
    else
        A2 = 4.0d0*m_G*(m_RK + m_G/3.0d0) / (m_RK + 4.0d0*m_G/3.0d0)
        B2 = 2.0d0*m_G*(m_RK - 2.0d0*m_G/3.0d0) / (m_RK + 4.0d0*m_G/3.0d0)
        C2 = A2
        D2 = (9.0d0*m_RK*m_G) / (3.0d0*m_RK + m_G)
        R2 = A2
        S2 = B2
        T2 = A2
        W2 = m_G
        Z2 = D2
    end if
    
    ! Apply crack reduction — matches all_adina.for DCRACK labels 15-33
    select case (NUMCRK)
    case (1)
        ! One crack in direction 1: reduce row/col 1
        do I = 1, 4
            C(1,I) = C(1,I) * m_STIFAC
        end do
        C(2,2) = A2;  C(2,4) = B2;  C(4,4) = C2
        
    case (2)
        ! Two cracks in directions 1 and 2: reduce rows/cols 1,2
        do I = 1, 2
            do J = I, 4
                C(I,J) = C(I,J) * m_STIFAC
            end do
        end do
        C(4,4) = D2
        
    case (4)
        ! One crack in direction 4(z): reduce row/col 4
        do I = 1, 4
            C(I,4) = C(I,4) * m_STIFAC
        end do
        C(1,1) = R2;  C(1,2) = S2;  C(2,2) = T2;  C(3,3) = W2
        
    case (5)
        ! Cracks in directions 1 and 4: reduce ALL, restore dir 2
        do I = 1, 4
            do J = I, 4
                C(I,J) = C(I,J) * m_STIFAC
            end do
        end do
        C(2,2) = Z2
        C(3,3) = W2 * m_SHEFAC
        
    case (6)
        ! Cracks in directions 1, 2, and 4: reduce ALL, restore shear
        do I = 1, 4
            do J = I, 4
                C(I,J) = C(I,J) * m_STIFAC
            end do
        end do
        C(3,3) = W2 * m_SHEFAC
    end select
    
    ! Enforce symmetry
    do I = 1, 3
        do J = I+1, 4
            C(J,I) = C(I,J)
        end do
    end do
end subroutine

!-----------------------------------------------------------------------
! transform_to_global: C_global = T^T * C_local * T
!-----------------------------------------------------------------------
subroutine transform_to_global(C, ANGLE)
    real(dp), intent(inout) :: C(4,4)
    real(dp), intent(in)    :: ANGLE
    real(dp) :: T(4,4), D(4,4), TANGLE, GAM, SG, CG
    integer  :: IR, IC, INN
    
    TANGLE = ANGLE
    if (TANGLE < -541.0d0) TANGLE = TANGLE + 722.0d0
    if (TANGLE < -180.0d0) TANGLE = TANGLE + 361.0d0
    if (TANGLE >  180.0d0) TANGLE = TANGLE - 180.0d0
    GAM = abs(TANGLE) * DEG2RAD
    SG = sin(GAM); CG = cos(GAM)
    
    T = 0.0d0
    T(1,1) =  CG*CG;       T(1,2) =  SG*SG;       T(1,3) =  CG*SG
    T(2,1) =  SG*SG;       T(2,2) =  CG*CG;       T(2,3) = -CG*SG
    T(3,1) = -2.0d0*CG*SG; T(3,2) =  2.0d0*CG*SG; T(3,3) =  CG*CG - SG*SG
    T(4,4) =  1.0d0
    
    D = 0.0d0
    do IR = 1, 4
        do IC = 1, 4
            do INN = 1, 4
                D(IR,IC) = D(IR,IC) + T(INN,IR) * C(INN,IC)
            end do
        end do
    end do
    
    C = 0.0d0
    do IR = 1, 4
        do IC = IR, 4
            do INN = 1, 4
                C(IR,IC) = C(IR,IC) + D(IR,INN) * T(INN,IC)
            end do
            C(IC,IR) = C(IR,IC)
        end do
    end do
end subroutine

!-----------------------------------------------------------------------
! transform_stress_to_global: SIG_out = T^T * SIGP_in
!-----------------------------------------------------------------------
subroutine transform_stress_to_global(SIGP_in, ANGLE, SIG_out)
    real(dp), intent(in)  :: SIGP_in(4), ANGLE
    real(dp), intent(out) :: SIG_out(4)
    real(dp) :: T(4,4), TANGLE, GAM, SG, CG
    integer  :: IR, IC
    
    TANGLE = ANGLE
    if (TANGLE < -541.0d0) TANGLE = TANGLE + 722.0d0
    if (TANGLE < -180.0d0) TANGLE = TANGLE + 361.0d0
    if (TANGLE >  180.0d0) TANGLE = TANGLE - 180.0d0
    GAM = abs(TANGLE) * DEG2RAD
    SG = sin(GAM); CG = cos(GAM)
    
    T = 0.0d0
    T(1,1) =  CG*CG;  T(1,2) = SG*SG;  T(1,3) =  CG*SG
    T(2,1) =  SG*SG;  T(2,2) = CG*CG;  T(2,3) = -CG*SG
    T(3,1) = -2.0d0*CG*SG; T(3,2) = 2.0d0*CG*SG; T(3,3) = CG*CG - SG*SG
    T(4,4) =  1.0d0
    
    SIG_out = 0.0d0
    do IR = 1, 4
        do IC = 1, 4
            SIG_out(IR) = SIG_out(IR) + T(IC,IR) * SIGP_in(IC)
        end do
    end do
end subroutine

!-----------------------------------------------------------------------
! condense_plane_stress: 4x4 → 3x3 static condensation (eliminate zz)
!-----------------------------------------------------------------------
subroutine condense_plane_stress(C4, C3)
    real(dp), intent(in)  :: C4(4,4)
    real(dp), intent(out) :: C3(3,3)
    real(dp) :: A_fac
    integer  :: I, J
    
    C3 = 0.0d0
    if (abs(C4(4,4)) < 1.0d-30) then
        C3(1:3,1:3) = C4(1:3,1:3)
        return
    end if
    
    do I = 1, 3
        A_fac = C4(I,4) / C4(4,4)
        do J = I, 3
            C3(I,J) = C4(I,J) - C4(4,J) * A_fac
            C3(J,I) = C3(I,J)
        end do
    end do
end subroutine

!-----------------------------------------------------------------------
! sort_principal: sort 3 values descending → P1 >= P2 >= P3
!-----------------------------------------------------------------------
subroutine sort_principal(SIG3, P1, P2, P3)
    real(dp), intent(in)  :: SIG3(3)
    real(dp), intent(out) :: P1, P2, P3
    real(dp) :: S(3), TEMP
    integer  :: I, IS
    
    S(1) = SIG3(1); S(2) = SIG3(2); S(3) = SIG3(3)
    do
        IS = 0
        do I = 1, 2
            if (S(I+1) > S(I)) then
                TEMP = S(I+1); S(I+1) = S(I); S(I) = TEMP
                IS = IS + 1
            end if
        end do
        if (IS == 0) exit
    end do
    P1 = S(1); P2 = S(2); P3 = S(3)
end subroutine

end module adina_concrete_mod


!===============================================================================
!  PSUMAT: OpenSees Plane Stress User Material Interface
!
!  Flow follows all_adina.for for MODEL=5 (concrete) with ITYP2D=2 (plane stress)
!  Addition: linear tension softening in crack directions
!===============================================================================
subroutine PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
    use adina_concrete_mod
    implicit none
    integer,  intent(in)    :: nstatevs, nprops
    real(dp), intent(in)    :: props(nprops)
    real(dp), intent(in)    :: strain0(3), strain1(3), dstrain(3)
    real(dp), intent(inout) :: stress(3), statev(nstatevs)
    real(dp), intent(out)   :: tangent(3,3)

    real(dp) :: E0
    real(dp) :: SP1_T(6), SP31_T(6), SP32_T(6), SP33_T(6)
    real(dp) :: SIG(4), EPS(4), STRAIN(4), STRESS_NEW(4)
    real(dp) :: C4(4,4)
    real(dp) :: EVMAX, ANGLE, CRKSTR(3), PGRAV
    real(dp) :: TMM, T11, T22, T33, T12
    real(dp) :: EVV, EMM, E11, E22, E33, E12_old
    real(dp) :: SBAR, EVI
    real(dp) :: ANG, ANGPRI
    real(dp) :: P1, P2, P3, SIG3(3)
    real(dp) :: PEE, S11, S22, S33, S12
    real(dp) :: EPSMM, EPS11, EPS22, EPS33, EPS12
    real(dp) :: eps_tu, eps_open1, eps_open2, ETATAU
    real(dp) :: eps_comp_max, eps_min_principal, eps_center, eps_radius, m_E_save
    real(dp) :: E_sec_b, fac_b, sig_t_b(3), alpha_b
    real(dp) :: C_11_sec, C_12_sec, G_sec
    real(dp) :: t_floor
    integer  :: IKAS, NUMCRK, ii_tf
    logical  :: is_cracked, is_crushed

    ! ==========================================================
    ! 1. EXTRACT MATERIAL PROPERTIES
    ! ==========================================================
    E0       = props(1)
    m_VNU    = props(2)
    m_SIGMAT = props(3)
    m_SIGMAC = props(4)
    m_EPSC   = props(5)
    m_SIGMAU = props(6)
    m_EPSU   = props(7)
    m_STIFAC = props(8)
    if (m_STIFAC < 1.0d-2) m_STIFAC = 1.0d-2
    m_SHEFAC = props(9)
    m_BETA   = props(10)
    m_GAMA   = props(11)
    m_RKAPA  = props(12)
    m_ALFA   = props(13)
    SP1_T(1:6)  = props(14:19)
    SP31_T(1:6) = props(20:25)
    SP32_T(1:6) = props(26:31)
    SP33_T(1:6) = props(32:37)
    
    m_E = E0
    m_RK = m_E / (3.0d0 * (1.0d0 - 2.0d0*m_VNU))
    m_G  = m_E / (2.0d0 * (1.0d0 + m_VNU))

    ! ==========================================================
    ! 2. INITIALIZE STATE VARIABLES (first call)
    ! ==========================================================
    if (nstatevs >= 8 .and. statev(7) < 0.5d0) then
        statev(1)   = 0.0d0
        statev(2)   = ANG_NOCRACK
        statev(3:5) = 0.0d0
        statev(6)   = 0.0d0
        statev(7)   = 1.0d0
        statev(8)   = 0.0d0
    end if
    
    EVMAX     = statev(1)
    ANGLE     = statev(2)
    CRKSTR(1) = statev(3)
    CRKSTR(2) = statev(4)
    CRKSTR(3) = statev(5)
    PGRAV     = statev(6)
    eps_comp_max = 0.0d0
    if (nstatevs >= 8) eps_comp_max = statev(8)

    ! ==========================================================
    ! 3. BUILD 4-COMPONENT VECTORS
    ! ==========================================================
    SIG(1) = stress(1);  SIG(2) = stress(2);  SIG(3) = stress(3)
    SIG(4) = 0.0d0
    
    EPS(1) = strain0(1);  EPS(2) = strain0(2);  EPS(3) = strain0(3)
    EPS(4) = (strain0(1) + strain0(2)) * m_VNU / (m_VNU - 1.0d0)
    
    STRAIN(1) = strain1(1);  STRAIN(2) = strain1(2);  STRAIN(3) = strain1(3)
    STRAIN(4) = (strain1(1) + strain1(2)) * m_VNU / (m_VNU - 1.0d0)

    ! ==========================================================
    ! 4. DECOMPOSE OLD STRESS/STRAIN (matches all_adina.for)
    ! ==========================================================
    TMM = (SIG(1) + SIG(2) + SIG(4)) / 3.0d0
    T11 = SIG(1) - TMM;  T22 = SIG(2) - TMM
    T33 = SIG(4) - TMM;  T12 = SIG(3)
    
    EVV    = -(EPS(1) + EPS(2) + EPS(4))
    EMM    = -EVV / 3.0d0
    E11    = EPS(1) - EMM;  E22 = EPS(2) - EMM
    E33    = EPS(4) - EMM;  E12_old = EPS(3)
    
    SBAR = ((SIG(1)-SIG(2))**2 + (SIG(1)-SIG(4))**2 + (SIG(2)-SIG(4))**2) / 6.0d0
    EVV  = sqrt(SBAR + SIG(3)**2) + 3.0d0 * m_ALFA * TMM
    
    IKAS = 1
    if (abs(EVV - EVMAX) > abs(EVMAX) * 1.0d-8) IKAS = -1

    ! ==========================================================
    ! 5. INITIALIZE WORKING VARIABLES
    ! ==========================================================
    m_E      = E0
    m_G      = E0 / (2.0d0*(1.0d0 + m_VNU))
    m_RK     = E0 / (3.0d0*(1.0d0 - 2.0d0*m_VNU))
    m_MOD45  = 1
    NUMCRK   = 0
    ANGPRI   = ANG_NOCRACK
    m_ILFSET = 0
    m_SIGCP  = m_SIGMAC
    m_EPSCP  = m_EPSC
    m_FALSTR = m_SIGMAT
    
    is_cracked = (ANGLE < 361.0d0)
    is_crushed = (PGRAV == PGRAV_CRUSH)
    
    eps_tu = 10.0d0 * m_SIGMAT / E0

    ! ==========================================================
    ! BRANCH A: CRUSHED STATE (PGRAV = 100)
    !   Saenz-curve secant modulus for gradual compressive softening.
    !   Total-strain approach: sigma = C(E_sec) * eps
    ! ==========================================================
    if (is_crushed) then
        ! Saenz-curve secant approach for post-peak softening
        ! (replaces near-zero stiffness that caused abrupt stress drop)
        m_SIGCP = m_SIGMAC; m_EPSCP = m_EPSC
        call compute_saenz_params()
        
        eps_center = 0.5d0 * (STRAIN(1) + STRAIN(2))
        eps_radius = sqrt((0.5d0*(STRAIN(1)-STRAIN(2)))**2 + (0.5d0*STRAIN(3))**2)
        eps_min_principal = eps_center - eps_radius
        if (eps_min_principal < eps_comp_max) eps_comp_max = eps_min_principal
        
        if (abs(eps_comp_max) > 1.0d-15) then
            SBAR = saenz_stress(eps_comp_max)
            m_E = SBAR / eps_comp_max
        else
            m_E = E0
        end if
        if (m_E < E0 * m_STIFAC) m_E = E0 * m_STIFAC
        
        m_RK = m_E / (3.0d0*(1.0d0 - 2.0d0*m_VNU))
        m_G  = m_E / (2.0d0*(1.0d0 + m_VNU))
        m_MOD45 = 1
        
        call build_constitutive_4x4(C4, 0)
        call condense_plane_stress(C4, tangent)
        
        STRESS_NEW(1) = tangent(1,1)*STRAIN(1) + tangent(1,2)*STRAIN(2) + tangent(1,3)*STRAIN(3)
        STRESS_NEW(2) = tangent(2,1)*STRAIN(1) + tangent(2,2)*STRAIN(2) + tangent(2,3)*STRAIN(3)
        STRESS_NEW(3) = tangent(3,1)*STRAIN(1) + tangent(3,2)*STRAIN(2) + tangent(3,3)*STRAIN(3)
        STRESS_NEW(4) = 0.0d0
        
        stress(1) = STRESS_NEW(1)
        stress(2) = STRESS_NEW(2)
        stress(3) = STRESS_NEW(3)
        if (nstatevs >= 8) statev(8) = eps_comp_max
        return
    end if

    ! ==========================================================
    ! BRANCH B: CRACKED STATE (ANGLE < 361)
    !   matches all_adina.for labels 46-47 → DCRACK path
    !   Hypoelastic update in crack coordinates + tension softening
    ! ==========================================================
    if (is_cracked) then
        alpha_b = 0.0d0
        E_sec_b = E0
        
        ! --- B1: Determine NUMCRK from ANGLE encoding ---
        ! Matches CRAKID KKK=3 label 107
        NUMCRK = 1
        if (ANGLE < 0.0d0 .and. ANGLE > -181.0d0) NUMCRK = 2
        if (ANGLE >= 180.0d0 .and. ANGLE < 361.0d0) NUMCRK = 0
        if (ANGLE < -180.0d0 .and. ANGLE > -362.0d0) NUMCRK = 4
        if (ANGLE < -361.0d0 .and. ANGLE > -542.0d0) NUMCRK = 6
        if (ANGLE < -541.0d0) NUMCRK = 5
        
        ! --- B2: Rotate old stress and strains to crack coordinates ---
        call rotate_stress_crack(SIG, ANGLE, m_SIGP)
        call rotate_strain_crack(EPS, ANGLE, m_TEP)
        call rotate_strain_crack(STRAIN, ANGLE, m_EP)
        
        ! --- B3: Check crack closing ---
        ! Matches CRAKID labels 111, 180
        if (NUMCRK == 1) then
            if (m_EP(1) < 0.0d0 .and. m_EP(1) < CRKSTR(1)) then
                NUMCRK = 0
                ANGLE = ANGLE + 180.0d0
            end if
        else if (NUMCRK == 2) then
            if (m_EP(1) < 0.0d0 .and. m_EP(1) < CRKSTR(1)) NUMCRK = 1
            if (m_EP(2) < 0.0d0 .and. m_EP(2) < CRKSTR(2)) NUMCRK = max(NUMCRK-1, 0)
            if (NUMCRK < 2) ANGLE = abs(ANGLE)
            if (NUMCRK == 0) ANGLE = ANGLE + 180.0d0
        else if (NUMCRK == 5 .or. NUMCRK == 6) then
            if (NUMCRK == 6 .and. m_EP(2) < 0.0d0 .and. m_EP(2) < CRKSTR(2)) NUMCRK = 5
            if (NUMCRK == 5 .and. m_EP(1) < 0.0d0 .and. m_EP(1) < CRKSTR(1)) NUMCRK = 4
            if (m_EP(4) < 0.0d0 .and. m_EP(4) < CRKSTR(3)) NUMCRK = NUMCRK - 4
        end if
        
        ! --- B4: Compute E for stress update ---
        SIG3(1) = m_SIGP(1); SIG3(2) = m_SIGP(2); SIG3(3) = m_SIGP(4)
        call sort_principal(SIG3, P1, P2, P3)
        call biaxial_envelope(P1, P2, P3, SP1_T, SP31_T, SP32_T, SP33_T)
        
        if (IKAS > 0) then
            call compute_tangent_moduli(3, .false.)
        end if
        
        ! B4a: Update eps_comp_max with current strain (needed for secant)
        eps_center = 0.5d0 * (STRAIN(1) + STRAIN(2))
        eps_radius = sqrt((0.5d0*(STRAIN(1)-STRAIN(2)))**2 + (0.5d0*STRAIN(3))**2)
        eps_min_principal = eps_center - eps_radius
        if (eps_min_principal < eps_comp_max) eps_comp_max = eps_min_principal
        
        ! Compression damage memory: degrade modulus during unloading
        ! only (IKAS <= 0).  During loading (IKAS > 0) the Saenz
        ! tangent from compute_tangent_moduli is already correct.
        E_sec_b = E0
        if (abs(eps_comp_max) > 1.0d-15) then
            m_E_save = m_E
            m_E = E0
            m_SIGCP = m_SIGMAC; m_EPSCP = m_EPSC
            call compute_saenz_params()
            SBAR = saenz_stress(eps_comp_max)
            m_E = m_E_save
            E_sec_b = SBAR / eps_comp_max
            if (E_sec_b < E0 * m_STIFAC) E_sec_b = E0 * m_STIFAC
            if (IKAS <= 0 .and. m_E > E_sec_b) m_E = E_sec_b
        end if
        
        if (NUMCRK == 0) m_MOD45 = 1
        if (m_E <= 0.0d0) m_E = E0 * m_STIFAC
        
        m_RK = m_E / (3.0d0*(1.0d0 - 2.0d0*m_VNU))
        m_G  = m_E / (2.0d0*(1.0d0 + m_VNU))
        
        ! --- B5: Build C4 in crack coords and hypoelastic update ---
        call build_constitutive_4x4(C4, NUMCRK)
        
        ! Plane stress constraint: dsig_4 = 0 → compute deps_4
        if (abs(C4(4,4)) > 1.0d-30) then
            PEE = -(C4(4,1)*(m_EP(1)-m_TEP(1)) + C4(4,2)*(m_EP(2)-m_TEP(2))) / C4(4,4)
        else
            PEE = m_EP(4) - m_TEP(4)
        end if
        
        m_SIGP(3) = m_SIGP(3) + C4(3,3) * (m_EP(3) - m_TEP(3))
        m_SIGP(1) = m_SIGP(1) + C4(1,1)*(m_EP(1)-m_TEP(1)) + C4(1,2)*(m_EP(2)-m_TEP(2)) &
                   + C4(1,4)*PEE
        m_SIGP(2) = m_SIGP(2) + C4(2,1)*(m_EP(1)-m_TEP(1)) + C4(2,2)*(m_EP(2)-m_TEP(2)) &
                   + C4(2,4)*PEE
        
        ! --- B5b: Secant stress correction for intact directions ---
        ! For non-cracked directions, replace the incrementally accumulated
        ! stress with the total-strain secant value: sigma = C(E_sec) * eps.
        ! E_sec = saenz(eps_comp_max)/eps_comp_max captures damage history.
        ! This eliminates hypoelastic drift while preserving tension softening.
        C_11_sec = E_sec_b / (1.0d0 - m_VNU * m_VNU)
        C_12_sec = m_VNU * C_11_sec
        G_sec    = E_sec_b / (2.0d0 * (1.0d0 + m_VNU))
        
        select case (NUMCRK)
        case (0, 4)
            ! All in-plane directions intact (closed crack / z-crack only)
            m_SIGP(1) = C_11_sec * m_EP(1) + C_12_sec * m_EP(2)
            m_SIGP(2) = C_12_sec * m_EP(1) + C_11_sec * m_EP(2)
            m_SIGP(3) = G_sec * m_EP(3)
            m_SIGP(4) = 0.0d0
        case (1, 5)
            ! Dir 1 cracked (handled by softening below); dir 2 intact
            m_SIGP(2) = C_11_sec * m_EP(2)
            m_SIGP(3) = G_sec * m_SHEFAC * m_EP(3)
        case (2, 6)
            ! Dirs 1,2 cracked; only shear remains as intact
            m_SIGP(3) = G_sec * m_SHEFAC * m_EP(3)
        end select
        
        ! --- B6: Tension softening + crack stress zeroing ---
        ETATAU = 0.0d0
        if (m_SHEFAC > 1.0d-3) ETATAU = 1.0d0
        
        if (NUMCRK == 1 .or. NUMCRK == 5) then
            eps_open1 = m_EP(1) - CRKSTR(1)
            if (eps_open1 >= 0.0d0 .and. eps_open1 < eps_tu) then
                m_SIGP(1) = m_SIGMAT * (1.0d0 - eps_open1 / eps_tu)
            else
                m_SIGP(1) = 0.0d0
            end if
            m_SIGP(3) = m_SIGP(3) * ETATAU
        end if
        
        if (NUMCRK == 2 .or. NUMCRK == 6) then
            eps_open1 = m_EP(1) - CRKSTR(1)
            if (eps_open1 >= 0.0d0 .and. eps_open1 < eps_tu) then
                m_SIGP(1) = m_SIGMAT * (1.0d0 - eps_open1 / eps_tu)
            else if (eps_open1 >= eps_tu) then
                m_SIGP(1) = 0.0d0
            end if
            eps_open2 = m_EP(2) - CRKSTR(2)
            if (eps_open2 >= 0.0d0 .and. eps_open2 < eps_tu) then
                m_SIGP(2) = m_SIGMAT * (1.0d0 - eps_open2 / eps_tu)
            else if (eps_open2 >= eps_tu) then
                m_SIGP(2) = 0.0d0
            end if
            m_SIGP(3) = m_SIGP(3) * ETATAU
        end if
        
        if (NUMCRK >= 4) m_SIGP(4) = 0.0d0
        
        ! --- B6b: Re-cracking check for closed cracks (NUMCRK=0) ---
        if (NUMCRK == 0 .and. m_SIGP(1) > m_SIGMAT) then
            NUMCRK = 1
            ANGLE = ANGLE - 180.0d0
            eps_open1 = m_EP(1) - CRKSTR(1)
            if (eps_open1 >= 0.0d0 .and. eps_open1 < eps_tu) then
                m_SIGP(1) = m_SIGMAT * (1.0d0 - eps_open1 / eps_tu)
            else if (eps_open1 >= eps_tu) then
                m_SIGP(1) = 0.0d0
            else
                m_SIGP(1) = m_SIGMAT
            end if
            m_SIGP(3) = 0.0d0
        end if
        
        ! --- B6c: Cap direction 2 at ft for intact directions ---
        if (NUMCRK <= 1 .and. m_SIGP(2) > m_SIGMAT) then
            m_SIGP(2) = m_SIGMAT
        end if
        
        ! --- B7: Transform stress to global ---
        call transform_stress_to_global(m_SIGP, ANGLE, STRESS_NEW)
        STRESS_NEW(4) = 0.0d0
        
        ! --- B7a: Shear-only isotropic blend for fixed-crack model ---
        !   Blends only the shear stress with an isotropic (damage) estimate.
        !   Using a uniform alpha_b (independent of NUMCRK) ensures no
        !   discontinuous jumps when cracks open/close during unloading.
        !   Normal stresses are NOT blended to preserve tension softening.
        if (NUMCRK /= 3) then
            alpha_b = 0.15d0
        else
            alpha_b = 0.0d0
        end if
        
        if (alpha_b > 0.0d0) then
            eps_center = 0.5d0 * (STRAIN(1) + STRAIN(2))
            eps_radius = sqrt((0.5d0*(STRAIN(1)-STRAIN(2)))**2 + (0.5d0*STRAIN(3))**2)
            eps_min_principal = eps_center - eps_radius
            
            if (eps_min_principal < -1.0d-10) then
                m_E_save = m_E
                m_E = E0
                m_SIGCP = m_SIGMAC; m_EPSCP = m_EPSC
                call compute_saenz_params()
                if (eps_min_principal <= eps_comp_max) eps_comp_max = eps_min_principal
                SBAR = saenz_stress(eps_comp_max)
                m_E = m_E_save
                if (abs(eps_comp_max) > 1.0d-15) then
                    E_sec_b = SBAR / eps_comp_max
                else
                    E_sec_b = E0
                end if
                if (E_sec_b < E0 * m_STIFAC) E_sec_b = E0 * m_STIFAC
            else
                E_sec_b = E0
            end if
            
            sig_t_b(3) = E_sec_b / (2.0d0*(1.0d0 + m_VNU)) * STRAIN(3)
            STRESS_NEW(3) = (1.0d0-alpha_b)*STRESS_NEW(3) + alpha_b*sig_t_b(3)
        end if
        
        ! --- B8: Check for crushing (strain-based) ---
        ! Only crush when compressive strain exceeds ultimate strain EPSU,
        ! allowing the Saenz softening branch to develop naturally.
        eps_center = 0.5d0 * (STRAIN(1) + STRAIN(2))
        eps_radius = sqrt((0.5d0*(STRAIN(1)-STRAIN(2)))**2 + (0.5d0*STRAIN(3))**2)
        eps_min_principal = eps_center - eps_radius
        if (eps_min_principal < eps_comp_max) eps_comp_max = eps_min_principal
        
        if (eps_comp_max < m_EPSU) then
            NUMCRK = 3
            PGRAV  = PGRAV_CRUSH
            CRKSTR(1) = m_SIGCP
            CRKSTR(2) = m_EPSCP
        end if
        
        ! --- B9: Update EVMAX ---
        SBAR = ((STRESS_NEW(1)-STRESS_NEW(2))**2 + (STRESS_NEW(1)-STRESS_NEW(4))**2 &
              + (STRESS_NEW(2)-STRESS_NEW(4))**2) / 6.0d0
        TMM  = (STRESS_NEW(1) + STRESS_NEW(2) + STRESS_NEW(4)) / 3.0d0
        EVI  = sqrt(SBAR + STRESS_NEW(3)**2) + 3.0d0*m_ALFA*TMM
        if (EVI > EVMAX) EVMAX = EVI
        
        ! --- B10: Build tangent for OpenSees ---
        if (IKAS > 0) then
            call compute_tangent_moduli(1, .false.)
        end if
        
        ! Use secant modulus for tangent: consistent with B5b secant stress.
        ! E_tangent can vanish at Saenz peak → singular K → Newton diverges.
        ! E_sec = saenz(eps_max)/eps_max > 0 always → stable convergence.
        if (abs(eps_comp_max) > 1.0d-15) then
            m_E = E_sec_b
        end if
        if (m_E < E0 * m_STIFAC) m_E = E0 * m_STIFAC
        
        m_RK = m_E / (3.0d0*(1.0d0 - 2.0d0*m_VNU))
        m_G  = m_E / (2.0d0*(1.0d0 + m_VNU))
        
        call build_constitutive_4x4(C4, NUMCRK)
        
        if (NUMCRK /= 3) then
            call transform_to_global(C4, ANGLE)
        end if
        call condense_plane_stress(C4, tangent)
        
        if (alpha_b > 0.0d0 .and. NUMCRK /= 3) then
            tangent(3,3) = (1.0d0-alpha_b)*tangent(3,3) &
                         + alpha_b*E_sec_b/(2.0d0*(1.0d0+m_VNU))
        end if
        
        ! --- B11: Output ---
        ! Enforce minimum tangent diagonal to prevent singular global K
        t_floor = E0 * m_STIFAC * m_STIFAC
        do ii_tf = 1, 3
            if (tangent(ii_tf,ii_tf) < t_floor) tangent(ii_tf,ii_tf) = t_floor
        end do
        
        stress(1) = STRESS_NEW(1)
        stress(2) = STRESS_NEW(2)
        stress(3) = STRESS_NEW(3)
        
        statev(1) = EVMAX
        statev(2) = ANGLE
        statev(3) = CRKSTR(1)
        statev(4) = CRKSTR(2)
        statev(5) = CRKSTR(3)
        statev(6) = PGRAV
        if (nstatevs >= 8) statev(8) = eps_comp_max
        return
    end if

    ! ==========================================================
    ! BRANCH C: UNCRACKED STATE
    !   For compression: total-strain approach (σ = C(E_sec) * ε)
    !     Gives correct secant unloading and avoids hypoelastic drift.
    !   For pure tension: hypoelastic with E0 (follows ADINA)
    ! ==========================================================
    
    ! --- C1: Determine if compression is present ---
    eps_center = 0.5d0 * (STRAIN(1) + STRAIN(2))
    eps_radius = sqrt((0.5d0*(STRAIN(1) - STRAIN(2)))**2 + (0.5d0*STRAIN(3))**2)
    eps_min_principal = eps_center - eps_radius
    
    if (eps_min_principal < -1.0d-10) then
        ! ============================================================
        ! C-COMPRESSION: Total-strain approach with secant modulus
        !   σ = C(E_sec) * ε, where E_sec = saenz(eps_max)/eps_max
        !   Uses UNIAXIAL Saenz curve (no biaxial envelope feedback).
        !   Biaxial effects emerge naturally from plane stress coupling.
        ! ============================================================
        m_SIGCP = m_SIGMAC
        m_EPSCP = m_EPSC
        call compute_saenz_params()
        
        ! Compute secant modulus
        if (eps_min_principal <= eps_comp_max) then
            eps_comp_max = eps_min_principal
        end if
        
        SBAR = saenz_stress(eps_comp_max)
        if (abs(eps_comp_max) > 1.0d-15) then
            m_E = SBAR / eps_comp_max
        else
            m_E = E0
        end if
        if (m_E < E0 * m_STIFAC) m_E = E0 * m_STIFAC
        
        ! Build isotropic constitutive matrix and compute stress
        m_RK = m_E / (3.0d0*(1.0d0 - 2.0d0*m_VNU))
        m_G  = m_E / (2.0d0*(1.0d0 + m_VNU))
        m_MOD45 = 1
        call build_constitutive_4x4(C4, 0)
        call condense_plane_stress(C4, tangent)
        
        STRESS_NEW(1) = tangent(1,1)*STRAIN(1) + tangent(1,2)*STRAIN(2) + tangent(1,3)*STRAIN(3)
        STRESS_NEW(2) = tangent(2,1)*STRAIN(1) + tangent(2,2)*STRAIN(2) + tangent(2,3)*STRAIN(3)
        STRESS_NEW(3) = tangent(3,1)*STRAIN(1) + tangent(3,2)*STRAIN(2) + tangent(3,3)*STRAIN(3)
        STRESS_NEW(4) = 0.0d0
        
        ! Check for cracking / crushing using biaxial envelope
        call compute_principal_2d(STRESS_NEW, ANG)
        SIG3(1) = m_SIGP(1); SIG3(2) = m_SIGP(2); SIG3(3) = m_SIGP(4)
        call sort_principal(SIG3, P1, P2, P3)
        call biaxial_envelope(P1, P2, P3, SP1_T, SP31_T, SP32_T, SP33_T)
        
        NUMCRK = 0
        if (.not. (P1 < 0.0d0 .and. P3 < 0.0d0)) then
            if (P1 > m_FALSTR) NUMCRK = 1
            if (P2 > m_FALSTR) NUMCRK = 2
        end if
        if (eps_comp_max < m_EPSU) NUMCRK = 3
        
        if (NUMCRK > 0 .and. NUMCRK /= 3) then
            ANGLE = ANG
            if (NUMCRK == 2) ANGLE = -ANG
            call rotate_strain_to_principal(STRAIN, ANG, m_EP)
            CRKSTR(1) = m_EP(1); CRKSTR(2) = m_EP(2); CRKSTR(3) = m_EP(4)
            m_ILFSET = 1
            
            ! Rebuild stress with cracked C matrix (crack-direction ≈ 0)
            m_RK = m_E / (3.0d0*(1.0d0 - 2.0d0*m_VNU))
            m_G  = m_E / (2.0d0*(1.0d0 + m_VNU))
            call build_constitutive_4x4(C4, NUMCRK)
            call transform_to_global(C4, ANGLE)
            call condense_plane_stress(C4, tangent)
            STRESS_NEW(1) = tangent(1,1)*STRAIN(1) + tangent(1,2)*STRAIN(2) + tangent(1,3)*STRAIN(3)
            STRESS_NEW(2) = tangent(2,1)*STRAIN(1) + tangent(2,2)*STRAIN(2) + tangent(2,3)*STRAIN(3)
            STRESS_NEW(3) = tangent(3,1)*STRAIN(1) + tangent(3,2)*STRAIN(2) + tangent(3,3)*STRAIN(3)
            STRESS_NEW(4) = 0.0d0
        else if (NUMCRK == 3) then
            PGRAV = PGRAV_CRUSH
            CRKSTR(1) = m_SIGCP; CRKSTR(2) = m_EPSCP
            m_ILFSET = 1
        end if
        
        ! EVMAX update
        SBAR = ((STRESS_NEW(1)-STRESS_NEW(2))**2 + (STRESS_NEW(1)-STRESS_NEW(4))**2 &
              + (STRESS_NEW(2)-STRESS_NEW(4))**2) / 6.0d0
        TMM  = (STRESS_NEW(1) + STRESS_NEW(2) + STRESS_NEW(4)) / 3.0d0
        EVI  = sqrt(SBAR + STRESS_NEW(3)**2) + 3.0d0*m_ALFA*TMM
        if (m_ILFSET == 1) EVMAX = EVI
        if (EVI > EVMAX) EVMAX = EVI
        
        ! Output
        stress(1) = STRESS_NEW(1)
        stress(2) = STRESS_NEW(2)
        stress(3) = STRESS_NEW(3)
        
        statev(1) = EVMAX; statev(2) = ANGLE
        statev(3) = CRKSTR(1); statev(4) = CRKSTR(2); statev(5) = CRKSTR(3)
        statev(6) = PGRAV
        if (nstatevs >= 8) statev(8) = eps_comp_max
        return
    end if
    
    ! ============================================================
    ! C-TENSION: Pure tension (eps_min_principal >= 0)
    !   Hypoelastic with E0, following ADINA label 50 path
    ! ============================================================
    m_E = E0
    m_RK = m_E / (3.0d0*(1.0d0 - 2.0d0*m_VNU))
    m_G  = m_E / (2.0d0*(1.0d0 + m_VNU))
    
    ! Total-strain approach with E0 (eliminates potential drift
    ! at compression→tension transitions where old stress used E_sec < E0)
    C_11_sec = E0 / (1.0d0 - m_VNU * m_VNU)
    C_12_sec = m_VNU * C_11_sec
    
    STRESS_NEW(1) = C_11_sec * STRAIN(1) + C_12_sec * STRAIN(2)
    STRESS_NEW(2) = C_12_sec * STRAIN(1) + C_11_sec * STRAIN(2)
    STRESS_NEW(3) = m_G * STRAIN(3)
    STRESS_NEW(4) = 0.0d0

    ! Check for cracking
    call compute_principal_2d(STRESS_NEW, ANG)
    call rotate_strain_to_principal(STRAIN, ANG, m_EP)
    SIG3(1) = m_SIGP(1); SIG3(2) = m_SIGP(2); SIG3(3) = m_SIGP(4)
    call sort_principal(SIG3, P1, P2, P3)
    call biaxial_envelope(P1, P2, P3, SP1_T, SP31_T, SP32_T, SP33_T)
    
    NUMCRK = 0
    if (P1 > m_FALSTR) NUMCRK = 1
    if (P2 > m_FALSTR) NUMCRK = 2
    
    if (NUMCRK > 0) then
        ANGLE = ANG
        if (NUMCRK == 2) ANGLE = -ANG
        CRKSTR(1) = m_EP(1); CRKSTR(2) = m_EP(2); CRKSTR(3) = m_EP(4)
        m_ILFSET = 1
        
        ! Reduce stress in crack direction to ft at onset
        if (NUMCRK >= 1 .and. m_SIGP(1) > m_SIGMAT) m_SIGP(1) = m_SIGMAT
        if (NUMCRK == 2 .and. m_SIGP(2) > m_SIGMAT) m_SIGP(2) = m_SIGMAT
        m_SIGP(3) = 0.0d0
        call transform_stress_to_global(m_SIGP, ANG, STRESS_NEW)
        STRESS_NEW(4) = 0.0d0
    end if
    
    SBAR = ((STRESS_NEW(1)-STRESS_NEW(2))**2 + (STRESS_NEW(1)-STRESS_NEW(4))**2 &
          + (STRESS_NEW(2)-STRESS_NEW(4))**2) / 6.0d0
    TMM  = (STRESS_NEW(1) + STRESS_NEW(2) + STRESS_NEW(4)) / 3.0d0
    EVI  = sqrt(SBAR + STRESS_NEW(3)**2) + 3.0d0*m_ALFA*TMM
    if (m_ILFSET == 1) EVMAX = EVI
    if (EVI > EVMAX) EVMAX = EVI
    
    ! Tangent: use cracked C if cracked, else elastic
    if (NUMCRK > 0) then
        call build_constitutive_4x4(C4, NUMCRK)
        call transform_to_global(C4, ANGLE)
    else
        call build_constitutive_4x4(C4, 0)
    end if
    call condense_plane_stress(C4, tangent)

    ! NaN/Inf guard
    if (any(tangent /= tangent) .or. any(abs(tangent) > 1.0d20)) then
        tangent = 0.0d0
        tangent(1,1) = E0 / (1.0d0 - m_VNU*m_VNU)
        tangent(2,2) = tangent(1,1)
        tangent(1,2) = m_VNU * tangent(1,1); tangent(2,1) = tangent(1,2)
        tangent(3,3) = E0 / (2.0d0*(1.0d0 + m_VNU))
    end if
    if (any(STRESS_NEW(1:3) /= STRESS_NEW(1:3)) .or. any(abs(STRESS_NEW(1:3)) > 1.0d20)) then
        STRESS_NEW(1:4) = 0.0d0
    end if

    ! Enforce minimum tangent diagonal to prevent singular global K
    t_floor = E0 * m_STIFAC * m_STIFAC
    do ii_tf = 1, 3
        if (tangent(ii_tf,ii_tf) < t_floor) tangent(ii_tf,ii_tf) = t_floor
    end do

    stress(1) = STRESS_NEW(1)
    stress(2) = STRESS_NEW(2)
    stress(3) = STRESS_NEW(3)
    
    statev(1) = EVMAX; statev(2) = ANGLE
    statev(3) = CRKSTR(1); statev(4) = CRKSTR(2); statev(5) = CRKSTR(3)
    statev(6) = PGRAV
    if (nstatevs >= 8) statev(8) = eps_comp_max
    
end subroutine PSUMAT
