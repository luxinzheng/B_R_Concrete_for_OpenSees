!===============================================================================
!
!  OpenSees User-Defined Material Subroutine (Fortran 90 Format)
!  Integrated with ADINA Plane Stress Concrete Constitutive Model
!
!  Author: Integrated from ADINA concrete model
!  Date: 2026-02-13
!  Language: Fortran 90/95
!
!  Description:
!    This module implements ADINA's concrete cracking model
!    for plane stress conditions in OpenSees framework
!
!===============================================================================

module concrete_constants
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: PI = 4.0_dp * atan(1.0_dp)
    real(dp), parameter :: SMALL = 1.0e-20_dp
    real(dp), parameter :: TOLERANCE = 1.0e-8_dp
end module concrete_constants

!===============================================================================

module concrete_data_types
    use concrete_constants
    implicit none
    
    ! Material properties structure
    type :: MaterialProps
        real(dp) :: E          ! Young's modulus
        real(dp) :: VNU        ! Poisson's ratio
        real(dp) :: SIGMAT     ! Tensile strength
        real(dp) :: SIGMAC     ! Compressive strength
        real(dp) :: EPSC       ! Peak compressive strain
        real(dp) :: SIGMAU     ! Ultimate compressive stress
        real(dp) :: EPSU       ! Ultimate compressive strain
        real(dp) :: BETA       ! Constitutive parameter
        real(dp) :: GAMA       ! Constitutive parameter
        real(dp) :: RKAPA      ! Constitutive parameter
        real(dp) :: ALFA       ! Constitutive parameter
        real(dp) :: STIFAC     ! Normal stiffness reduction
        real(dp) :: SHEFAC     ! Shear stiffness reduction
    end type MaterialProps
    
    ! State variables structure
    type :: StateVars
        real(dp) :: EVMAX      ! Maximum equivalent strain
        real(dp) :: ANGLE      ! Crack angle
        real(dp) :: PGRAV      ! Crack indicator
        real(dp) :: EPSCP      ! Current compressive strain
        real(dp) :: SIGCP      ! Current compressive stress
        real(dp) :: FALSTR     ! Failure strength
        real(dp) :: SIGP(4)    ! Principal stresses
        real(dp) :: TEP(4)     ! Previous principal strains
        real(dp) :: EP(4)      ! Current principal strains
        integer  :: NUMCRK     ! Number of cracks
        integer  :: MOD45      ! Modified stiffness flag
        integer  :: ILFSET     ! Lifecycle flag
    end type StateVars
    
end module concrete_data_types

!===============================================================================

module concrete_utilities
    use concrete_constants
    implicit none
    
contains

    !---------------------------------------------------------------------------
    ! Gauss integration points and weights (3-point rule)
    !---------------------------------------------------------------------------
    subroutine get_gauss_points(xg, wgt)
        real(dp), intent(out) :: xg(3), wgt(3)
        
        xg(1)  = -0.7745966692_dp
        xg(2)  =  0.0_dp
        xg(3)  =  0.7745966692_dp
        
        wgt(1) =  0.5555555556_dp
        wgt(2) =  0.8888888889_dp
        wgt(3) =  0.5555555556_dp
    end subroutine get_gauss_points

    !---------------------------------------------------------------------------
    ! Calculate bulk and shear moduli
    !---------------------------------------------------------------------------
    subroutine calc_elastic_moduli(E, VNU, RK, G)
        real(dp), intent(in)  :: E, VNU
        real(dp), intent(out) :: RK, G
        
        RK = E / (3.0_dp * (1.0_dp - 2.0_dp * VNU))
        G  = E / (2.0_dp * (1.0_dp + VNU))
    end subroutine calc_elastic_moduli

    !---------------------------------------------------------------------------
    ! Calculate stress deviator and volumetric components
    !---------------------------------------------------------------------------
    subroutine calc_stress_components(stress, TMM, T11, T22, T33, T12)
        real(dp), intent(in)  :: stress(4)
        real(dp), intent(out) :: TMM, T11, T22, T33, T12
        
        TMM = (stress(1) + stress(2) + stress(4)) / 3.0_dp
        T11 = stress(1) - TMM
        T22 = stress(2) - TMM
        T33 = stress(4) - TMM
        T12 = stress(3)
    end subroutine calc_stress_components

    !---------------------------------------------------------------------------
    ! Calculate strain deviator and volumetric components
    !---------------------------------------------------------------------------
    subroutine calc_strain_components(strain, EVV, EMM, E11, E22, E33, E12)
        real(dp), intent(in)  :: strain(4)
        real(dp), intent(out) :: EVV, EMM, E11, E22, E33, E12
        
        EVV = -(strain(1) + strain(2) + strain(4))
        EMM = -EVV / 3.0_dp
        E11 = strain(1) - EMM
        E22 = strain(2) - EMM
        E33 = strain(4) - EMM
        E12 = strain(3)
    end subroutine calc_strain_components

end module concrete_utilities

!===============================================================================

module principal_stress_strain
    use concrete_constants
    implicit none
    
contains

    !---------------------------------------------------------------------------
    ! Calculate principal stresses and strains
    !---------------------------------------------------------------------------
    subroutine prncpl_concrete(sig, eps, ang, sigp, epsp)
        real(dp), intent(in)  :: sig(4), eps(4)
        real(dp), intent(out) :: ang, sigp(4), epsp(4)
        
        real(dp) :: A, B, S1, S2, S3, E1, E2, E3
        real(dp) :: SXX, SYY, SZZ, TXY
        real(dp) :: EXX, EYY, EZZ, GXY
        
        ! Stress components
        SXX = sig(1)
        SYY = sig(2)
        TXY = sig(3)
        SZZ = sig(4)
        
        ! Strain components
        EXX = eps(1)
        EYY = eps(2)
        GXY = eps(3)
        EZZ = eps(4)
        
        ! 2D principal stress calculation
        A = (SXX + SYY) / 2.0_dp
        B = sqrt(((SXX - SYY) / 2.0_dp)**2 + TXY**2)
        
        S1 = A + B
        S2 = A - B
        S3 = SZZ
        
        ! Store principal stresses (algebraically largest first)
        sigp(1) = S1
        sigp(2) = S2
        sigp(4) = S3
        sigp(3) = 0.0_dp
        
        ! Principal angle (in degrees)
        if (abs(SXX - SYY) > 1.0e-10_dp) then
            ang = 0.5_dp * atan2(2.0_dp * TXY, SXX - SYY) * 180.0_dp / PI
        else
            if (abs(TXY) > 1.0e-10_dp) then
                ang = 45.0_dp
            else
                ang = 0.0_dp
            end if
        end if
        
        ! 2D principal strain calculation
        A = (EXX + EYY) / 2.0_dp
        B = sqrt(((EXX - EYY) / 2.0_dp)**2 + (GXY / 2.0_dp)**2)
        
        E1 = A + B
        E2 = A - B
        E3 = EZZ
        
        epsp(1) = E1
        epsp(2) = E2
        epsp(4) = E3
        epsp(3) = 0.0_dp
        
    end subroutine prncpl_concrete

end module principal_stress_strain

!===============================================================================

module crack_identification
    use concrete_constants
    use principal_stress_strain
    implicit none
    
contains

    !---------------------------------------------------------------------------
    ! Identify and track crack state
    !---------------------------------------------------------------------------
    subroutine crakid_concrete(sig, eps, pgrav, crkstr, angle, tep, ep, numcrk)
        real(dp), intent(in)    :: sig(4), eps(4), pgrav
        real(dp), intent(inout) :: crkstr(3), angle
        real(dp), intent(inout) :: tep(4), ep(4)
        integer,  intent(inout) :: numcrk
        
        real(dp) :: sigp(4), epsp(4), ang
        integer  :: i
        
        ! Get principal stresses and strains
        call prncpl_concrete(sig, eps, ang, sigp, epsp)
        
        ! Update principal strains
        do i = 1, 4
            ep(i) = epsp(i)
        end do
        
        ! Track crack parameters
        if (numcrk > 0) then
            crkstr(1) = epsp(1)
            crkstr(2) = epsp(2)
            crkstr(3) = epsp(4)
        end if
        
    end subroutine crakid_concrete

end module crack_identification

!===============================================================================

module cracked_stiffness
    use concrete_constants
    implicit none
    
contains

    !---------------------------------------------------------------------------
    ! Calculate cracked stiffness matrix
    !---------------------------------------------------------------------------
    subroutine dcrack_concrete(C, sig, angle, numcrk, mod45, &
                               E, VNU, RK, G, YP, E12, E14, E24, &
                               STIFAC, SHEFAC)
        real(dp), intent(out) :: C(4,4)
        real(dp), intent(in)  :: sig(4), angle
        integer,  intent(in)  :: numcrk, mod45
        real(dp), intent(in)  :: E, VNU, RK, G, YP(3)
        real(dp), intent(in)  :: E12, E14, E24, STIFAC, SHEFAC
        
        real(dp) :: DUM, A1, B1, C1, D1
        real(dp) :: A2, B2, C2, D2, R2, S2, T2, W2, Z2
        real(dp) :: TANGLE, GAM, SG, CG
        real(dp) :: T_MATRIX(4,4), D_MATRIX(4,4)
        integer  :: i, j, ir, ic, in
        
        !-----------------------------------------------------------------------
        ! Build elastic stiffness matrix
        !-----------------------------------------------------------------------
        if (mod45 == 1) then
            ! Isotropic case
            DUM = 2.0_dp * G / 3.0_dp
            A1 = RK + 2.0_dp * DUM
            B1 = RK - DUM
            
            C(1,1) = A1
            C(1,2) = B1
            C(1,3) = 0.0_dp
            C(1,4) = B1
            C(2,2) = A1
            C(2,3) = 0.0_dp
            C(2,4) = B1
            C(3,3) = G
            C(3,4) = 0.0_dp
            C(4,4) = A1
        else
            ! Orthotropic case
            A1 = (1.0_dp + VNU) * (1.0_dp - 2.0_dp * VNU)
            B1 = (1.0_dp - VNU) / A1
            D1 = VNU / A1
            C1 = 1.0_dp / (2.0_dp * (1.0_dp + VNU))
            
            C(1,1) = B1 * YP(1)
            C(1,2) = E12 * D1
            C(1,3) = 0.0_dp
            C(1,4) = E14 * D1
            C(2,2) = B1 * YP(2)
            C(2,3) = 0.0_dp
            C(2,4) = E24 * D1
            C(3,3) = C1 * E12
            C(3,4) = 0.0_dp
            C(4,4) = B1 * YP(3)
        end if

        !-----------------------------------------------------------------------
        ! Apply cracking reduction factors
        !-----------------------------------------------------------------------
        if (numcrk == 0 .or. numcrk == 3) goto 35
        
        ! Reduce shear stiffness
        C(3,3) = C(3,3) * SHEFAC
        
        if (mod45 == 2) then
            ! Orthotropic reduction
            A1 = 1.0_dp / (1.0_dp - VNU * VNU)
            B1 = VNU * A1
            A2 = A1 * YP(2)
            B2 = B1 * E24
            C2 = A1 * YP(3)
            D2 = YP(3)
            R2 = A1 * YP(1)
            S2 = B1 * E12
            T2 = A2
            W2 = C1 * E12
            Z2 = E12
        else
            ! Isotropic reduction
            A2 = 4.0_dp*G*(RK+G/3.0_dp)/(RK+4.0_dp*G/3.0_dp)
            B2 = 2.0_dp*G*(RK-2.0_dp*G/3.0_dp)/(RK+4.0_dp*G/3.0_dp)
            C2 = A2
            D2 = (9.0_dp*RK*G)/(3.0_dp*RK + G)
            R2 = A2
            S2 = B2
            T2 = A2
            W2 = G
            Z2 = D2
        end if
        
        ! Apply reduction based on crack pattern
        select case (numcrk)
        case (1)
            ! Single crack in direction 1
            do i = 1, 4
                C(1,i) = C(1,i) * STIFAC
            end do
            C(2,2) = A2
            C(2,4) = B2
            C(4,4) = C2
            
        case (2)
            ! Cracks in directions 1 and 2
            do i = 1, 2
                do j = i, 4
                    C(i,j) = C(i,j) * STIFAC
                end do
            end do
            C(4,4) = D2
            
        case (4)
            ! Crack in direction 3
            do i = 1, 4
                C(i,4) = C(i,4) * STIFAC
            end do
            C(1,1) = R2
            C(1,2) = S2
            C(2,2) = T2
            C(3,3) = W2
            
        case (5)
            ! Cracks in directions 1 and 3
            do i = 1, 4
                do j = i, 4
                    C(i,j) = C(i,j) * STIFAC
                end do
            end do
            C(2,2) = Z2
            C(3,3) = W2 * SHEFAC
            
        case (6)
            ! Complete cracking
            do i = 1, 4
                do j = i, 4
                    C(i,j) = C(i,j) * STIFAC
                end do
            end do
            C(3,3) = W2 * SHEFAC
        end select

35      continue

        !-----------------------------------------------------------------------
        ! Symmetrize matrix
        !-----------------------------------------------------------------------
        do i = 1, 3
            do j = i+1, 4
                C(j,i) = C(i,j)
            end do
        end do

        !-----------------------------------------------------------------------
        ! Coordinate transformation (if cracked)
        !-----------------------------------------------------------------------
        if (numcrk > 0 .and. numcrk /= 3 .and. angle < 3.61e2_dp) then
            
            TANGLE = angle
            
            ! Normalize angle
            if (TANGLE < -5.41e2_dp) TANGLE = TANGLE + 722.0_dp
            if (TANGLE < -1.8e2_dp)  TANGLE = TANGLE + 361.0_dp
            if (TANGLE >  1.8e2_dp)  TANGLE = TANGLE - 180.0_dp
            
            GAM = abs(TANGLE) * PI / 180.0_dp
            SG = sin(GAM)
            CG = cos(GAM)
            
            ! Build transformation matrix
            T_MATRIX(1,1) =  CG**2
            T_MATRIX(1,2) =  SG**2
            T_MATRIX(1,3) =  CG * SG
            T_MATRIX(1,4) =  0.0_dp
            T_MATRIX(2,1) =  T_MATRIX(1,2)
            T_MATRIX(2,2) =  T_MATRIX(1,1)
            T_MATRIX(2,3) = -T_MATRIX(1,3)
            T_MATRIX(2,4) =  0.0_dp
            T_MATRIX(3,1) =  T_MATRIX(2,3) * 2.0_dp
            T_MATRIX(3,2) = -T_MATRIX(3,1)
            T_MATRIX(3,3) =  T_MATRIX(1,1) - T_MATRIX(1,2)
            T_MATRIX(3,4) =  0.0_dp
            T_MATRIX(4,1) =  0.0_dp
            T_MATRIX(4,2) =  0.0_dp
            T_MATRIX(4,3) =  0.0_dp
            T_MATRIX(4,4) =  1.0_dp
            
            ! First transformation: D = T^T * C
            do ir = 1, 4
                do ic = 1, 4
                    D_MATRIX(ir,ic) = 0.0_dp
                    do in = 1, 4
                        D_MATRIX(ir,ic) = D_MATRIX(ir,ic) + &
                                         T_MATRIX(in,ir) * C(in,ic)
                    end do
                end do
            end do
            
            ! Second transformation: C = D * T
            do ir = 1, 4
                do ic = ir, 4
                    C(ir,ic) = 0.0_dp
                    do in = 1, 4
                        C(ir,ic) = C(ir,ic) + &
                                  D_MATRIX(ir,in) * T_MATRIX(in,ic)
                    end do
                    C(ic,ir) = C(ir,ic)
                end do
            end do
        end if

        !-----------------------------------------------------------------------
        ! Plane stress condensation (eliminate z-direction)
        !-----------------------------------------------------------------------
        do i = 1, 3
            A1 = C(i,4) / C(4,4)
            do j = i, 3
                C(i,j) = C(i,j) - C(4,j) * A1
                C(j,i) = C(i,j)
            end do
        end do

    end subroutine dcrack_concrete

end module cracked_stiffness

!===============================================================================
! MAIN MATERIAL SUBROUTINE
!===============================================================================

subroutine PSUMAT(nstatevs, nprops, props, &
                  stress, strain0, strain1, dstrain, &
                  statev, tangent)
    
    use concrete_constants
    use concrete_data_types
    use concrete_utilities
    use principal_stress_strain
    use crack_identification
    use cracked_stiffness
    
    implicit none
    
    ! Arguments
    integer,  intent(in)    :: nstatevs, nprops
    real(dp), intent(in)    :: props(nprops)
    real(dp), intent(inout) :: stress(4), strain0(4), strain1(4), dstrain(4)
    real(dp), intent(inout) :: statev(nstatevs)
    real(dp), intent(out)   :: tangent(4,4)
    
    ! Local variables - Material properties
    type(MaterialProps) :: mat
    type(StateVars)     :: state
    
    ! Local variables - Computation
    real(dp) :: RK, G, ANGLE, EVMAX, EVGRAV, PGRAV
    real(dp) :: SIGP(4), TEP(4), EP(4), YP(3)
    real(dp) :: E12, E14, E24, EPSCP, SIGCP, FALSTR
    real(dp) :: TMM, T11, T22, T33, T12
    real(dp) :: EVV, EMM, E11, E22, E33, E12_DEV
    real(dp) :: SBAR, EVI, EPSMM, EPS11, EPS22, EPS33, EPS12
    real(dp) :: PEE, S11, S22, S33, S12
    real(dp) :: RP, ES, EU, RAM5, RBM5, RCM5
    real(dp) :: DENM, DE, TY, XG(3), WGT(3)
    real(dp) :: ANG, CRKSTR(3)
    integer  :: IKAS, NUMCRK, MOD45, ILFSET
    integer  :: I, J, K, L, NUMINT, N1, N2
    
    !---------------------------------------------------------------------------
    ! Read material properties from props array
    !---------------------------------------------------------------------------
    mat%E      = props(1)
    mat%VNU    = props(2)
    mat%SIGMAT = props(3)
    mat%SIGMAC = props(4)
    mat%EPSC   = props(5)
    mat%SIGMAU = props(6)
    mat%EPSU   = props(7)
    mat%BETA   = props(8)
    mat%GAMA   = props(9)
    mat%RKAPA  = props(10)
    mat%ALFA   = props(11)
    mat%STIFAC = props(12)
    mat%SHEFAC = props(13)

    !---------------------------------------------------------------------------
    ! Read state variables
    !---------------------------------------------------------------------------
    EVMAX   = statev(1)
    ANGLE   = statev(2)
    PGRAV   = statev(3)
    EPSCP   = statev(4)
    SIGCP   = statev(5)
    FALSTR  = statev(6)
    
    do I = 1, 4
        SIGP(I) = statev(6+I)
        TEP(I)  = statev(10+I)
        EP(I)   = statev(14+I)
    end do
    
    NUMCRK  = int(statev(19))
    MOD45   = int(statev(20))
    ILFSET  = int(statev(21))

    !---------------------------------------------------------------------------
    ! Initialize for first call
    !---------------------------------------------------------------------------
    if (EVMAX < SMALL) then
        EVMAX  = 0.0_dp
        ANGLE  = 1.0e3_dp
        PGRAV  = 0.0_dp
        EPSCP  = mat%EPSC
        SIGCP  = mat%SIGMAC
        FALSTR = 0.0_dp
        NUMCRK = 0
        MOD45  = 1
        ILFSET = 0
        SIGP   = 0.0_dp
        TEP    = 0.0_dp
        EP     = 0.0_dp
    end if

    !---------------------------------------------------------------------------
    ! Calculate elastic moduli
    !---------------------------------------------------------------------------
    call calc_elastic_moduli(mat%E, mat%VNU, RK, G)

    !---------------------------------------------------------------------------
    ! Calculate stress and strain components
    !---------------------------------------------------------------------------
    call calc_stress_components(stress, TMM, T11, T22, T33, T12)
    call calc_strain_components(strain1, EVV, EMM, E11, E22, E33, E12_DEV)

    !---------------------------------------------------------------------------
    ! Calculate equivalent strain (Concrete Model)
    !---------------------------------------------------------------------------
    SBAR = ((stress(1)-stress(2))**2 + (stress(1)-stress(4))**2 + &
            (stress(2)-stress(4))**2) / 6.0_dp
    EVI = sqrt(SBAR + stress(3)**2) + 3.0_dp * mat%ALFA * TMM

    !---------------------------------------------------------------------------
    ! Determine loading/unloading
    !---------------------------------------------------------------------------
    IKAS = 1
    if (abs(EVI - EVMAX) > abs(EVMAX) * TOLERANCE) IKAS = -1

    !---------------------------------------------------------------------------
    ! Calculate principal stresses and strains
    !---------------------------------------------------------------------------
    if (IKAS > 0 .and. (ANGLE > 3.61e2_dp .or. PGRAV == 1.0e2_dp)) then
        call prncpl_concrete(stress, strain1, ANG, SIGP, EP)
    end if
    
    if (ANGLE < 3.61e2_dp .and. PGRAV /= 1.0e2_dp) then
        call crakid_concrete(stress, strain1, PGRAV, CRKSTR, ANGLE, &
                           TEP, EP, NUMCRK)
    end if

    !---------------------------------------------------------------------------
    ! Calculate tangent moduli
    !---------------------------------------------------------------------------
    call get_gauss_points(XG, WGT)
    
    if (PGRAV /= 1.0e2_dp) then
        NUMINT = 3
        DENM = abs(SIGP(1)) + abs(SIGP(2)) + abs(SIGP(4))
        
        if (DENM > 0.00001_dp * mat%SIGMAT) then
            ! Calculate softening parameters
            RP = mat%EPSU / mat%EPSC
            ES = SIGCP / EPSCP
            EU = (mat%SIGMAU * SIGCP / mat%SIGMAC) / (mat%EPSU * EPSCP / mat%EPSC)
            
            RAM5 = mat%E/EU + (RP-2.0_dp)*RP*RP*mat%E/ES - &
                   (2.0_dp*RP+1.0_dp)*(RP-1.0_dp)**2
            RAM5 = RAM5 / (RP * (RP - 1.0_dp)**2)
            RBM5 = 2.0_dp * mat%E / ES - 3.0_dp - 2.0_dp * RAM5
            RCM5 = 2.0_dp - mat%E / ES + RAM5
            
            ! Calculate secant moduli for each principal direction
            N1 = 1
            N2 = 4
            if (PGRAV == 1.0e2_dp) then
                N1 = 4
                if (SIGP(2) < SIGP(4)) N1 = 2
                N2 = N1
            end if
            
            do J = N1, N2
                if (J == 3) cycle
                I = J
                if (J == 4) I = 3
                YP(I) = mat%E
                
                if (SIGP(J) < 0.001_dp * SIGCP) then
                    YP(I) = 0.0_dp
                    do L = 1, NUMINT
                        DE = TEP(J) + (1.0_dp + XG(L)) * (EP(J) - TEP(J)) / 2.0_dp
                        DE = DE / EPSCP
                        TY = mat%E
                        if (DE > 0.0_dp) then
                            TY = mat%E * (1.0_dp - RBM5*DE*DE - 2.0_dp*RCM5*DE**3)
                            TY = TY / (1.0_dp + RAM5*DE + RBM5*DE*DE + RCM5*DE**3)**2
                        end if
                        YP(I) = YP(I) + 0.5_dp * WGT(L) * TY
                    end do
                end if
            end do
            
            ! Check for different cracking states
            MOD45 = 1
            DE = mat%RKAPA * SIGCP
            if (SIGP(1) < DE .or. SIGP(2) < DE .or. SIGP(4) < DE) MOD45 = 2
            
            if (MOD45 == 1) then
                mat%E = (abs(SIGP(1))*YP(1) + abs(SIGP(2))*YP(2) + &
                        abs(SIGP(4))*YP(3)) / DENM
            else
                ! Calculate directional moduli
                E12 = mat%E
                DENM = abs(SIGP(1)) + abs(SIGP(2))
                if (DENM /= 0.0_dp) &
                    E12 = (abs(SIGP(1))*YP(1) + abs(SIGP(2))*YP(2)) / DENM
                
                E14 = mat%E
                DENM = abs(SIGP(1)) + abs(SIGP(4))
                if (DENM /= 0.0_dp) &
                    E14 = (abs(SIGP(1))*YP(1) + abs(SIGP(4))*YP(3)) / DENM
                
                E24 = mat%E
                DENM = abs(SIGP(2)) + abs(SIGP(4))
                if (DENM /= 0.0_dp) &
                    E24 = (abs(SIGP(2))*YP(2) + abs(SIGP(4))*YP(3)) / DENM
            end if
        else
            mat%E = props(1)
        end if
    else
        mat%E = 0.0_dp
    end if

    !---------------------------------------------------------------------------
    ! Update stresses (Plane Stress)
    !---------------------------------------------------------------------------
    EPSMM = -EVV / 3.0_dp
    EPS11 = strain1(1) - EPSMM
    EPS22 = strain1(2) - EPSMM
    EPS33 = strain1(4) - EPSMM
    EPS12 = strain1(3)
    
    PEE = 3.0_dp * RK * (EPSMM - EMM) + TMM
    
    S11 = 2.0_dp * G * (EPS11 - E11) + T11
    S22 = 2.0_dp * G * (EPS22 - E22) + T22
    S33 = 2.0_dp * G * (EPS33 - E33) + T33
    S12 = G * (EPS12 - E12_DEV) + T12
    
    stress(1) = S11 + PEE
    stress(2) = S22 + PEE
    stress(3) = S12
    stress(4) = 0.0_dp  ! Plane stress condition

    !---------------------------------------------------------------------------
    ! Build tangent stiffness matrix
    !---------------------------------------------------------------------------
    call dcrack_concrete(tangent, stress, ANGLE, NUMCRK, MOD45, &
                        mat%E, mat%VNU, RK, G, YP, E12, E14, E24, &
                        mat%STIFAC, mat%SHEFAC)

    !---------------------------------------------------------------------------
    ! Update principal stresses
    !---------------------------------------------------------------------------
    call prncpl_concrete(stress, strain1, ANG, SIGP, EP)

    !---------------------------------------------------------------------------
    ! Check for cracking
    !---------------------------------------------------------------------------
    if (ANGLE > 3.61e2_dp) then
        if (SIGP(1) < 0.0_dp .and. SIGP(4) < 0.0_dp) then
            NUMCRK = 0
        else
            if (SIGP(1) > FALSTR) NUMCRK = 1
            if (SIGP(2) > FALSTR) NUMCRK = 2
            if (SIGP(4) > FALSTR) NUMCRK = NUMCRK + 4
            if (SIGP(2) < SIGCP .or. SIGP(4) < SIGCP) NUMCRK = 3
            
            if (NUMCRK > 0) ANGLE = ANG
            
            if (NUMCRK == 3) then
                PGRAV = 1.0e2_dp
                FALSTR = SIGCP
                CRKSTR(1) = SIGCP
                CRKSTR(2) = EPSCP
            else if (NUMCRK > 0) then
                ANGLE = ANG
                if (NUMCRK == 2) ANGLE = -ANGLE
                if (NUMCRK >= 4) ANGLE = ANGLE - 361.0_dp
                CRKSTR(1) = EP(1)
                CRKSTR(2) = EP(2)
                CRKSTR(3) = EP(4)
            end if
        end if
    else
        call crakid_concrete(stress, strain1, PGRAV, CRKSTR, ANGLE, &
                           TEP, EP, NUMCRK)
    end if

    !---------------------------------------------------------------------------
    ! Update state at end of step
    !---------------------------------------------------------------------------
    TEP = EP
    
    if (EVI > EVMAX) EVMAX = EVI

    !---------------------------------------------------------------------------
    ! Write state variables
    !---------------------------------------------------------------------------
    statev(1) = EVMAX
    statev(2) = ANGLE
    statev(3) = PGRAV
    statev(4) = EPSCP
    statev(5) = SIGCP
    statev(6) = FALSTR
    
    do I = 1, 4
        statev(6+I)  = SIGP(I)
        statev(10+I) = TEP(I)
        statev(14+I) = EP(I)
    end do
    
    statev(19) = real(NUMCRK, dp)
    statev(20) = real(MOD45, dp)
    statev(21) = real(ILFSET, dp)

end subroutine PSUMAT
