!===============================================================================
!
!  ADINA Plane Stress Concrete Constitutive Model
!  Pure PSUMAT Interface (Fortran 90)
!
!  Interface: psumat (standard user material interface)
!  Author: Based on ADINA concrete model
!  Date: 2026-02-13
!
!===============================================================================

module concrete_constants_psumat
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: PI = 4.0_dp * atan(1.0_dp)
    real(dp), parameter :: SMALL = 1.0e-20_dp
    real(dp), parameter :: TOLERANCE = 1.0e-8_dp
end module concrete_constants_psumat

!===============================================================================

module concrete_storage_psumat
    use concrete_constants_psumat
    implicit none
    
    ! 状态变量存储 (每个积分点)
    integer, parameter :: MAX_ELEMENTS = 10000
    integer, parameter :: MAX_IPT = 10
    
    type :: StateVariables
        real(dp) :: EVMAX      ! 最大等效应变
        real(dp) :: ANGLE      ! 裂缝角度
        real(dp) :: PGRAV      ! 裂缝指示符
        real(dp) :: EPSCP      ! 当前压应变
        real(dp) :: SIGCP      ! 当前压应力
        real(dp) :: FALSTR     ! 破坏强度
        real(dp) :: SIGP(4)    ! 主应力
        real(dp) :: TEP(4)     ! 上一步主应变
        real(dp) :: EP(4)      ! 当前主应变
        real(dp) :: YP(3)      ! 切线模量
        integer  :: NUMCRK     ! 裂缝数量
        integer  :: MOD45      ! 刚度模式
        integer  :: ILFSET     ! 生命周期标志
        logical  :: initialized
    end type StateVariables
    
    ! 全局状态数组
    type(StateVariables), save :: state_array(MAX_ELEMENTS, MAX_IPT)
    logical, save :: storage_initialized = .false.
    
contains

    subroutine init_storage()
        integer :: i, j
        if (.not. storage_initialized) then
            do i = 1, MAX_ELEMENTS
                do j = 1, MAX_IPT
                    state_array(i,j)%EVMAX = 0.0_dp
                    state_array(i,j)%ANGLE = 1.0e3_dp
                    state_array(i,j)%PGRAV = 0.0_dp
                    state_array(i,j)%EPSCP = 0.002_dp
                    state_array(i,j)%SIGCP = 30.0e6_dp
                    state_array(i,j)%FALSTR = 0.0_dp
                    state_array(i,j)%SIGP = 0.0_dp
                    state_array(i,j)%TEP = 0.0_dp
                    state_array(i,j)%EP = 0.0_dp
                    state_array(i,j)%YP = 0.0_dp
                    state_array(i,j)%NUMCRK = 0
                    state_array(i,j)%MOD45 = 1
                    state_array(i,j)%ILFSET = 0
                    state_array(i,j)%initialized = .false.
                end do
            end do
            storage_initialized = .true.
        end if
    end subroutine init_storage

end module concrete_storage_psumat

!===============================================================================

module concrete_utilities_psumat
    use concrete_constants_psumat
    implicit none
    
contains

    subroutine calc_elastic_moduli(E, VNU, RK, G)
        real(dp), intent(in)  :: E, VNU
        real(dp), intent(out) :: RK, G
        RK = E / (3.0_dp * (1.0_dp - 2.0_dp * VNU))
        G  = E / (2.0_dp * (1.0_dp + VNU))
    end subroutine calc_elastic_moduli

    subroutine get_gauss_points(xg, wgt)
        real(dp), intent(out) :: xg(3), wgt(3)
        xg(1)  = -0.7745966692_dp
        xg(2)  =  0.0_dp
        xg(3)  =  0.7745966692_dp
        wgt(1) =  0.5555555556_dp
        wgt(2) =  0.8888888889_dp
        wgt(3) =  0.5555555556_dp
    end subroutine get_gauss_points

    subroutine prncpl_stress_strain(sig, eps, ang, sigp, epsp)
        real(dp), intent(in)  :: sig(4), eps(4)
        real(dp), intent(out) :: ang, sigp(4), epsp(4)
        real(dp) :: A, B, S1, S2, S3, E1, E2, E3
        
        ! 2D主应力
        A = (sig(1) + sig(2)) / 2.0_dp
        B = sqrt(((sig(1) - sig(2)) / 2.0_dp)**2 + sig(3)**2)
        S1 = A + B
        S2 = A - B
        S3 = sig(4)
        
        sigp(1) = S1
        sigp(2) = S2
        sigp(4) = S3
        sigp(3) = 0.0_dp
        
        ! 主应力角度
        if (abs(sig(1) - sig(2)) > 1.0e-10_dp) then
            ang = 0.5_dp * atan2(2.0_dp * sig(3), sig(1) - sig(2)) * 180.0_dp / PI
        else
            ang = 0.0_dp
            if (abs(sig(3)) > 1.0e-10_dp) ang = 45.0_dp
        end if
        
        ! 2D主应变
        A = (eps(1) + eps(2)) / 2.0_dp
        B = sqrt(((eps(1) - eps(2)) / 2.0_dp)**2 + (eps(3) / 2.0_dp)**2)
        E1 = A + B
        E2 = A - B
        E3 = eps(4)
        
        epsp(1) = E1
        epsp(2) = E2
        epsp(4) = E3
        epsp(3) = 0.0_dp
        
    end subroutine prncpl_stress_strain

end module concrete_utilities_psumat

!===============================================================================

module concrete_stiffness_psumat
    use concrete_constants_psumat
    implicit none
    
contains

    subroutine build_tangent_matrix(C, ANGLE, NUMCRK, MOD45, &
                                    E, VNU, RK, G, YP, E12, E14, E24, &
                                    STIFAC, SHEFAC, ITYP2D)
        real(dp), intent(out) :: C(4,4)
        real(dp), intent(in)  :: ANGLE, E, VNU, RK, G, YP(3)
        real(dp), intent(in)  :: E12, E14, E24, STIFAC, SHEFAC
        integer,  intent(in)  :: NUMCRK, MOD45, ITYP2D
        
        real(dp) :: DUM, A1, B1, C1, D1
        real(dp) :: A2, B2, C2, D2, R2, S2, T2, W2, Z2
        real(dp) :: TANGLE, GAM, SG, CG
        real(dp) :: T_MATRIX(4,4), D_MATRIX(4,4)
        integer  :: i, j, ir, ic, in
        
        ! 初始化
        C = 0.0_dp
        
        ! 构建弹性刚度矩阵
        if (MOD45 == 1) then
            ! 各向同性
            DUM = 2.0_dp * G / 3.0_dp
            A1 = RK + 2.0_dp * DUM
            B1 = RK - DUM
            
            C(1,1) = A1; C(1,2) = B1; C(1,3) = 0.0_dp; C(1,4) = B1
            C(2,2) = A1; C(2,3) = 0.0_dp; C(2,4) = B1
            C(3,3) = G;  C(3,4) = 0.0_dp
            C(4,4) = A1
        else
            ! 正交各向异性
            A1 = (1.0_dp + VNU) * (1.0_dp - 2.0_dp * VNU)
            B1 = (1.0_dp - VNU) / A1
            D1 = VNU / A1
            C1 = 1.0_dp / (2.0_dp * (1.0_dp + VNU))
            
            C(1,1) = B1 * YP(1); C(1,2) = E12 * D1
            C(1,3) = 0.0_dp;     C(1,4) = E14 * D1
            C(2,2) = B1 * YP(2); C(2,3) = 0.0_dp
            C(2,4) = E24 * D1
            C(3,3) = C1 * E12;   C(3,4) = 0.0_dp
            C(4,4) = B1 * YP(3)
        end if

        ! 开裂刚度折减
        if (NUMCRK /= 0 .and. NUMCRK /= 3) then
            C(3,3) = C(3,3) * SHEFAC
            
            if (MOD45 == 2) then
                A1 = 1.0_dp / (1.0_dp - VNU * VNU)
                B1 = VNU * A1
                A2 = A1 * YP(2); B2 = B1 * E24; C2 = A1 * YP(3)
                D2 = YP(3); R2 = A1 * YP(1); S2 = B1 * E12
                T2 = A2; W2 = C1 * E12; Z2 = E12
            else
                A2 = 4.0_dp*G*(RK+G/3.0_dp)/(RK+4.0_dp*G/3.0_dp)
                B2 = 2.0_dp*G*(RK-2.0_dp*G/3.0_dp)/(RK+4.0_dp*G/3.0_dp)
                C2 = A2; D2 = (9.0_dp*RK*G)/(3.0_dp*RK + G)
                R2 = A2; S2 = B2; T2 = A2; W2 = G; Z2 = D2
            end if
            
            select case (NUMCRK)
            case (1)
                do i = 1, 4
                    C(1,i) = C(1,i) * STIFAC
                end do
                C(2,2) = A2; C(2,4) = B2; C(4,4) = C2
            case (2)
                do i = 1, 2
                    do j = i, 4
                        C(i,j) = C(i,j) * STIFAC
                    end do
                end do
                C(4,4) = D2
            case (4)
                do i = 1, 4
                    C(i,4) = C(i,4) * STIFAC
                end do
                C(1,1) = R2; C(1,2) = S2; C(2,2) = T2; C(3,3) = W2
            case (5, 6)
                do i = 1, 4
                    do j = i, 4
                        C(i,j) = C(i,j) * STIFAC
                    end do
                end do
                if (NUMCRK == 5) C(2,2) = Z2
                C(3,3) = W2 * SHEFAC
            end select
        end if

        ! 对称化
        do i = 1, 3
            do j = i+1, 4
                C(j,i) = C(i,j)
            end do
        end do

        ! 坐标转换
        if (NUMCRK > 0 .and. NUMCRK /= 3 .and. ANGLE < 3.61e2_dp) then
            TANGLE = ANGLE
            if (TANGLE < -5.41e2_dp) TANGLE = TANGLE + 722.0_dp
            if (TANGLE < -1.8e2_dp)  TANGLE = TANGLE + 361.0_dp
            if (TANGLE >  1.8e2_dp)  TANGLE = TANGLE - 180.0_dp
            
            GAM = abs(TANGLE) * PI / 180.0_dp
            SG = sin(GAM); CG = cos(GAM)
            
            T_MATRIX(1,1) =  CG**2;     T_MATRIX(1,2) =  SG**2
            T_MATRIX(1,3) =  CG * SG;   T_MATRIX(1,4) =  0.0_dp
            T_MATRIX(2,1) =  SG**2;     T_MATRIX(2,2) =  CG**2
            T_MATRIX(2,3) = -CG * SG;   T_MATRIX(2,4) =  0.0_dp
            T_MATRIX(3,1) = -2.0_dp*CG*SG; T_MATRIX(3,2) = 2.0_dp*CG*SG
            T_MATRIX(3,3) =  CG**2 - SG**2; T_MATRIX(3,4) = 0.0_dp
            T_MATRIX(4,:) = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
            
            ! D = T^T * C
            do ir = 1, 4
                do ic = 1, 4
                    D_MATRIX(ir,ic) = 0.0_dp
                    do in = 1, 4
                        D_MATRIX(ir,ic) = D_MATRIX(ir,ic) + T_MATRIX(in,ir) * C(in,ic)
                    end do
                end do
            end do
            
            ! C = D * T
            do ir = 1, 4
                do ic = ir, 4
                    C(ir,ic) = 0.0_dp
                    do in = 1, 4
                        C(ir,ic) = C(ir,ic) + D_MATRIX(ir,in) * T_MATRIX(in,ic)
                    end do
                    C(ic,ir) = C(ir,ic)
                end do
            end do
        end if

        ! 平面应力缩减
        if (ITYP2D >= 2) then
            do i = 1, 3
                A1 = C(i,4) / C(4,4)
                do j = i, 3
                    C(i,j) = C(i,j) - C(4,j) * A1
                    C(j,i) = C(i,j)
                end do
            end do
        end if

    end subroutine build_tangent_matrix

end module concrete_stiffness_psumat

!===============================================================================
! PSUMAT主接口
!===============================================================================

subroutine PSUMAT(nstatevs, nprops, props, &
                  stress, strain0, strain1, dstrain, &
                  statev, tangent)
    
    use concrete_constants_psumat
    use concrete_storage_psumat
    use concrete_utilities_psumat
    use concrete_stiffness_psumat
    
    implicit none
    
    ! 参数
    integer,  intent(in)    :: nstatevs, nprops
    real(dp), intent(in)    :: props(nprops)
    real(dp), intent(inout) :: stress(4), strain0(4), strain1(4), dstrain(4)
    real(dp), intent(inout) :: statev(nstatevs)
    real(dp), intent(out)   :: tangent(4,4)
    
    ! 局部变量 - 材料参数
    real(dp) :: E, VNU, ft, fc, epsc, fu, epsu
    real(dp) :: beta, gama, rkapa, alfa, stifac, shefac
    
    ! 局部变量 - 计算
    real(dp) :: RK, G, EVMAX, ANGLE, PGRAV, EPSCP, SIGCP, FALSTR
    real(dp) :: SIGP(4), TEP(4), EP(4), YP(3)
    real(dp) :: E12, E14, E24
    real(dp) :: TMM, T11, T22, T33, T12
    real(dp) :: EVV, EMM, E11, E22, E33, E12_dev
    real(dp) :: SBAR, EVI, DENM, DE, TY
    real(dp) :: EPSMM, EPS11, EPS22, EPS33, EPS12
    real(dp) :: PEE, S11, S22, S33, S12
    real(dp) :: RP, ES, EU, RAM5, RBM5, RCM5
    real(dp) :: XG(3), WGT(3), ANG, CRKSTR(3)
    integer  :: NUMCRK, MOD45, ILFSET, IKAS
    integer  :: I, J, L, N1, N2, NUMINT, ITYP2D
    
    ! 初始化存储
    call init_storage()
    
    ! 读取材料参数
    E      = props(1)
    VNU    = props(2)
    ft     = props(3)
    fc     = props(4)
    epsc   = props(5)
    fu     = props(6)
    epsu   = props(7)
    beta   = props(8)
    gama   = props(9)
    rkapa  = props(10)
    alfa   = props(11)
    stifac = props(12)
    shefac = props(13)
    
    ! 平面应力类型
    ITYP2D = 2  ! 平面应力
    
    ! 初始化或读取状态变量
    if (nstatevs >= 21) then
        ! 使用外部状态变量
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
        
        ! 初次调用初始化
        if (EVMAX < SMALL) then
            EVMAX = 0.0_dp
            ANGLE = 1.0e3_dp
            PGRAV = 0.0_dp
            EPSCP = epsc
            SIGCP = fc
            FALSTR = 0.0_dp
            SIGP = 0.0_dp
            TEP = 0.0_dp
            EP = 0.0_dp
            NUMCRK = 0
            MOD45 = 1
            ILFSET = 0
        end if
    else
        ! 内部状态管理 (简化版)
        EVMAX = 0.0_dp
        ANGLE = 1.0e3_dp
        PGRAV = 0.0_dp
        EPSCP = epsc
        SIGCP = fc
        FALSTR = 0.0_dp
        SIGP = 0.0_dp
        TEP = 0.0_dp
        EP = 0.0_dp
        NUMCRK = 0
        MOD45 = 1
        ILFSET = 0
    end if
    
    ! 计算弹性模量
    call calc_elastic_moduli(E, VNU, RK, G)
    
    ! 应力应变分解
    TMM = (stress(1) + stress(2) + stress(4)) / 3.0_dp
    T11 = stress(1) - TMM
    T22 = stress(2) - TMM
    T33 = stress(4) - TMM
    T12 = stress(3)
    
    EVV = -(strain1(1) + strain1(2) + strain1(4))
    EMM = -EVV / 3.0_dp
    E11 = strain1(1) - EMM
    E22 = strain1(2) - EMM
    E33 = strain1(4) - EMM
    E12_dev = strain1(3)
    
    ! 等效应变
    SBAR = ((stress(1)-stress(2))**2 + (stress(1)-stress(4))**2 + &
            (stress(2)-stress(4))**2) / 6.0_dp
    EVI = sqrt(SBAR + stress(3)**2) + 3.0_dp * alfa * TMM
    
    ! 加载/卸载判断
    IKAS = 1
    if (abs(EVI - EVMAX) > abs(EVMAX) * TOLERANCE) IKAS = -1
    
    ! 主应力主应变
    if (IKAS > 0 .and. (ANGLE > 3.61e2_dp .or. PGRAV == 1.0e2_dp)) then
        call prncpl_stress_strain(stress, strain1, ANG, SIGP, EP)
    end if
    
    ! 计算切线模量
    call get_gauss_points(XG, WGT)
    
    if (PGRAV /= 1.0e2_dp) then
        NUMINT = 3
        DENM = abs(SIGP(1)) + abs(SIGP(2)) + abs(SIGP(4))
        
        if (DENM > 0.00001_dp * ft) then
            RP = epsu / epsc
            ES = SIGCP / EPSCP
            EU = (fu * SIGCP / fc) / (epsu * EPSCP / epsc)
            
            RAM5 = E/EU + (RP-2.0_dp)*RP*RP*E/ES - (2.0_dp*RP+1.0_dp)*(RP-1.0_dp)**2
            RAM5 = RAM5 / (RP * (RP - 1.0_dp)**2)
            RBM5 = 2.0_dp * E / ES - 3.0_dp - 2.0_dp * RAM5
            RCM5 = 2.0_dp - E / ES + RAM5
            
            N1 = 1; N2 = 4
            
            do J = N1, N2
                if (J == 3) cycle
                I = J
                if (J == 4) I = 3
                YP(I) = E
                
                if (SIGP(J) < 0.001_dp * SIGCP) then
                    YP(I) = 0.0_dp
                    do L = 1, NUMINT
                        DE = TEP(J) + (1.0_dp + XG(L)) * (EP(J) - TEP(J)) / 2.0_dp
                        DE = DE / EPSCP
                        TY = E
                        if (DE > 0.0_dp) then
                            TY = E * (1.0_dp - RBM5*DE*DE - 2.0_dp*RCM5*DE**3)
                            TY = TY / (1.0_dp + RAM5*DE + RBM5*DE*DE + RCM5*DE**3)**2
                        end if
                        YP(I) = YP(I) + 0.5_dp * WGT(L) * TY
                    end do
                end if
            end do
            
            MOD45 = 1
            DE = rkapa * SIGCP
            if (SIGP(1) < DE .or. SIGP(2) < DE .or. SIGP(4) < DE) MOD45 = 2
            
            if (MOD45 == 1) then
                E = (abs(SIGP(1))*YP(1) + abs(SIGP(2))*YP(2) + &
                     abs(SIGP(4))*YP(3)) / DENM
            else
                E12 = E; E14 = E; E24 = E
                DENM = abs(SIGP(1)) + abs(SIGP(2))
                if (DENM /= 0.0_dp) E12 = (abs(SIGP(1))*YP(1) + abs(SIGP(2))*YP(2))/DENM
                DENM = abs(SIGP(1)) + abs(SIGP(4))
                if (DENM /= 0.0_dp) E14 = (abs(SIGP(1))*YP(1) + abs(SIGP(4))*YP(3))/DENM
                DENM = abs(SIGP(2)) + abs(SIGP(4))
                if (DENM /= 0.0_dp) E24 = (abs(SIGP(2))*YP(2) + abs(SIGP(4))*YP(3))/DENM
            end if
        else
            E = props(1)
        end if
    else
        E = 0.0_dp
    end if
    
    ! 更新应力
    EPSMM = -EVV / 3.0_dp
    EPS11 = strain1(1) - EPSMM
    EPS22 = strain1(2) - EPSMM
    EPS33 = strain1(4) - EPSMM
    EPS12 = strain1(3)
    
    PEE = 3.0_dp * RK * (EPSMM - EMM) + TMM
    S11 = 2.0_dp * G * (EPS11 - E11) + T11
    S22 = 2.0_dp * G * (EPS22 - E22) + T22
    S33 = 2.0_dp * G * (EPS33 - E33) + T33
    S12 = G * (EPS12 - E12_dev) + T12
    
    stress(1) = S11 + PEE
    stress(2) = S22 + PEE
    stress(3) = S12
    stress(4) = 0.0_dp  ! 平面应力
    
    ! 构建切线刚度矩阵
    call build_tangent_matrix(tangent, ANGLE, NUMCRK, MOD45, &
                             E, VNU, RK, G, YP, E12, E14, E24, &
                             stifac, shefac, ITYP2D)
    
    ! 更新主应力
    call prncpl_stress_strain(stress, strain1, ANG, SIGP, EP)
    
    ! 检查开裂
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
            end if
        end if
    end if
    
    ! 更新状态
    TEP = EP
    if (EVI > EVMAX) EVMAX = EVI
    
    ! 保存状态变量
    if (nstatevs >= 21) then
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
    end if
    
end subroutine PSUMAT
