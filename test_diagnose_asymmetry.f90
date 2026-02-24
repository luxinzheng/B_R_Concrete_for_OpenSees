!===============================================================================
!  Diagnose asymmetry: detailed trace of crack state and stress during
!  cyclic shear with axial compression.
!  Focus: does direction-2 stress exceed ft after crack closing?
!===============================================================================
program test_diagnose_asymmetry
    use adina_concrete_mod
    implicit none
    integer, parameter :: NP=37, NS=13
    real(dp) :: props(NP), stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    real(dp) :: ANGLE, CRKSTR1, CRKSTR2
    real(dp) :: peak_pos, peak_neg, ft
    integer :: i, NUMCRK_decoded
    
    call setup_props(props)
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ft = props(3)
    peak_pos = 0.0d0; peak_neg = 0.0d0
    
    write(*,'(A)') '================================================================'
    write(*,'(A)') 'ASYMMETRY DIAGNOSIS: cyclic shear + axial compression'
    write(*,'(A)') '================================================================'
    write(*,'(A)') 'Phase  Step  gamma     sig_xx     sig_yy     tau_xy     ANGLE   NUMCRK  PGRAV  C11'
    
    ! Phase 1: Axial compression (eps_yy = -3e-4)
    do i = 1, 3
        strain1 = (/0.0d0, -1.0d-4*i, 0.0d0/)
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
    end do
    call print_state('AX', 3, strain1(3), stress, statev, tangent)
    
    ! Phase 2: Forward shear to gamma=1e-3 (small, avoid crushing)
    do i = 1, 10
        strain1 = (/0.0d0, -3.0d-4, 1.0d-4*i/)
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        call print_state('F+', i, strain1(3), stress, statev, tangent)
        if (stress(3) > peak_pos) peak_pos = stress(3)
    end do
    
    ! Phase 3: Return to zero
    do i = 1, 10
        strain1 = (/0.0d0, -3.0d-4, 1.0d-3 - 1.0d-4*i/)
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (mod(i,2)==0) call print_state('RT', i, strain1(3), stress, statev, tangent)
    end do
    
    ! Phase 4: Reverse shear to gamma=-1e-3
    write(*,'(A)') '--- REVERSE LOADING (key phase) ---'
    do i = 1, 10
        strain1 = (/0.0d0, -3.0d-4, -1.0d-4*i/)
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        call print_state('F-', i, strain1(3), stress, statev, tangent)
        if (stress(3) < peak_neg) peak_neg = stress(3)
    end do
    
    ! Phase 5: Return to zero (from negative)
    do i = 1, 10
        strain1 = (/0.0d0, -3.0d-4, -1.0d-3 + 1.0d-4*i/)
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (mod(i,2)==0) call print_state('R2', i, strain1(3), stress, statev, tangent)
    end do
    
    ! Phase 6: Second forward to gamma=1e-3
    write(*,'(A)') '--- SECOND FORWARD (compare with first) ---'
    do i = 1, 10
        strain1 = (/0.0d0, -3.0d-4, 1.0d-4*i/)
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        call print_state('2+', i, strain1(3), stress, statev, tangent)
    end do
    
    write(*,*)
    write(*,'(A,ES12.4)') '  Peak +tau (first):  ', peak_pos
    write(*,'(A,ES12.4)') '  Peak -tau (reverse): ', peak_neg
    write(*,'(A,F6.1,A)') '  Asymmetry ratio: ', abs(peak_neg/peak_pos)*100.0d0, '%'

contains
subroutine setup_props(props)
    real(dp), intent(out) :: props(NP)
    props = 0.0d0
    props(1)=21.4d9; props(2)=0.2d0; props(3)=2.07d6; props(4)=-20.7d6
    props(5)=-0.002d0; props(6)=-4.14d6; props(7)=-0.006d0
    props(8)=1.0d-4; props(9)=0.5d0; props(10)=0.75d0; props(11)=1.0d0
    props(12)=0.7d0; props(13)=0.12d0
    props(14:19)=(/0d0,0.25d0,0.5d0,0.75d0,1d0,1.2d0/)
    props(20:25)=(/1d0,1.1d0,1.2d0,1.15d0,1d0,0.9d0/)
    props(26:31)=(/1.3d0,1.5d0,2d0,2.3d0,2.7d0,3.2d0/)
    props(32:37)=(/1.25d0,1.45d0,1.95d0,2.25d0,2.65d0,3.15d0/)
end subroutine

subroutine print_state(phase, step, gam, sig, sv, tang)
    character(len=*), intent(in) :: phase
    integer, intent(in) :: step
    real(dp), intent(in) :: gam, sig(3), sv(NS), tang(3,3)
    real(dp) :: ang
    integer :: nc
    ang = sv(2)
    nc = -1
    if (ang >= 999.0d0) nc = -1
    if (ang >= 0.0d0 .and. ang < 180.0d0) nc = 1
    if (ang < 0.0d0 .and. ang > -181.0d0) nc = 2
    if (ang >= 180.0d0 .and. ang < 361.0d0) nc = 0
    write(*,'(A2,I5,F10.6,3ES11.3,F8.1,I5,F7.1,ES10.2)') &
        phase, step, gam, sig(1), sig(2), sig(3), sv(2), nc, sv(6), tang(3,3)
end subroutine
end program
