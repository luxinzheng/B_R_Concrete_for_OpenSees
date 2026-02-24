!===============================================================================
! Deep diagnosis: why is reverse direction force weaker?
! Track crack state evolution, direction-2 stress overshoot, and crushing
! under realistic wall element strain paths.
!===============================================================================
program test_diagnose_reverse
    use adina_concrete_mod
    implicit none
    integer, parameter :: NP=37, NS=13
    real(dp) :: props(NP), stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    real(dp) :: peak_fwd, peak_rev
    integer :: i, nsteps
    real(dp) :: t, gam_max
    real(dp) :: sig_p1, sig_p2, ang_p
    
    call setup_props(props)
    nsteps = 50
    gam_max = 1.5d-3
    
    write(*,'(A)') '================================================================'
    write(*,'(A)') '  REVERSE DIRECTION DIAGNOSIS'
    write(*,'(A)') '================================================================'
    
    ! ---- Test A: Boundary element path (flexure-dominated) ----
    write(*,'(A)') ''
    write(*,'(A)') '=== TEST A: Boundary Element (tension -> compression -> tension) ==='
    write(*,'(A)') 'Step  eps_xx     sig_xx     sig_yy     tau_xy   ANGLE   NC  PGRAV  dir2_overshoot'
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    
    ! A1: Axial precompression
    strain1 = (/0.0d0, -2.0d-4, 0.0d0/)
    dstrain = strain1 - strain0
    call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
    strain0 = strain1
    
    ! A2: Tension (boundary stretches) — crack forms
    do i = 1, nsteps
        t = dble(i)/dble(nsteps)
        strain1(1) = 5.0d-4 * t
        strain1(2) = -2.0d-4
        strain1(3) = 0.0d0
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (mod(i,10)==0) call pstate('A2', i, stress, statev)
    end do
    
    ! A3: Compression (boundary compresses) — crack closes
    do i = 1, nsteps
        t = dble(i)/dble(nsteps)
        strain1(1) = 5.0d-4 * (1.0d0-t) + (-8.0d-4) * t
        strain1(2) = -2.0d-4
        strain1(3) = 0.0d0
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (mod(i,10)==0) call pstate('A3', i, stress, statev)
    end do
    
    ! A4: Re-tension (reverse cycle) 
    do i = 1, nsteps
        t = dble(i)/dble(nsteps)
        strain1(1) = -8.0d-4 * (1.0d0-t) + 5.0d-4 * t
        strain1(2) = -2.0d-4
        strain1(3) = 0.0d0
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (mod(i,10)==0) call pstate('A4', i, stress, statev)
    end do
    
    ! ---- Test B: Web element path (shear-dominated) ----
    write(*,'(A)') ''
    write(*,'(A)') '=== TEST B: Web Element (axial + forward shear -> reverse shear) ==='
    write(*,'(A)') 'Step  gamma      sig_xx     sig_yy     tau_xy   ANGLE   NC  PGRAV  dir2_overshoot'
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    peak_fwd = 0.0d0; peak_rev = 0.0d0
    
    ! B1: Axial precompression
    strain1 = (/0.0d0, -3.0d-4, 0.0d0/)
    dstrain = strain1 - strain0
    call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
    strain0 = strain1
    
    ! B2: Forward shear
    do i = 1, nsteps
        t = dble(i)/dble(nsteps)
        strain1(1) = 0.0d0
        strain1(2) = -3.0d-4
        strain1(3) = gam_max * t
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (stress(3) > peak_fwd) peak_fwd = stress(3)
        if (mod(i,10)==0) call pstate('B2', i, stress, statev)
    end do
    
    ! B3: Return to zero
    do i = 1, nsteps
        t = dble(i)/dble(nsteps)
        strain1(3) = gam_max * (1.0d0-t)
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (mod(i,10)==0) call pstate('B3', i, stress, statev)
    end do
    
    ! B4: Reverse shear (KEY PHASE)
    write(*,'(A)') '--- B4: REVERSE SHEAR ---'
    do i = 1, nsteps
        t = dble(i)/dble(nsteps)
        strain1(3) = -gam_max * t
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (stress(3) < peak_rev) peak_rev = stress(3)
        if (mod(i,5)==0) call pstate('B4', i, stress, statev)
    end do
    
    write(*,*)
    write(*,'(A,ES12.4)') '  Peak +tau (forward):  ', peak_fwd
    write(*,'(A,ES12.4)') '  Peak -tau (reverse):  ', peak_rev
    if (abs(peak_fwd) > 1.0d-10) then
        write(*,'(A,F6.1,A)') '  Ratio |rev/fwd|: ', abs(peak_rev/peak_fwd)*100.0d0, '%'
    end if
    
    ! ---- Test C: Track principal stress and direction-2 overshoot ----
    write(*,'(A)') ''
    write(*,'(A)') '=== TEST C: Direction-2 Stress Tracking (NUMCRK=1 during reverse) ==='
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    
    ! C1: Crack via tension
    strain1 = (/2.0d-4, -1.0d-4, 3.0d-4/)
    dstrain = strain1 - strain0
    call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
    strain0 = strain1
    write(*,'(A,3ES11.3,A,F7.1,I3)') '  After crack:  ', stress, '  ANG=', statev(2), nint(statev(6))
    
    ! C2: Partial unload (crack stays open, NUMCRK=1)
    strain1 = (/1.0d-4, -1.0d-4, 1.0d-4/)
    dstrain = strain1 - strain0
    call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
    strain0 = strain1
    write(*,'(A,3ES11.3,A,F7.1,I3)') '  Partial unload:', stress, '  ANG=', statev(2), nint(statev(6))
    
    ! C3: Reverse shear (NUMCRK=1, direction 2 tension develops)
    do i = 1, 10
        strain1 = (/1.0d-4, -1.0d-4, 1.0d-4 - 1.5d-4*i/)
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        
        ! Check direction-2 stress in crack coords
        call check_dir2(stress, statev, props(3))
    end do

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

subroutine pstate(phase, step, sig, sv)
    character(len=*), intent(in) :: phase
    integer, intent(in) :: step
    real(dp), intent(in) :: sig(3), sv(NS)
    real(dp) :: ang, sigp(4), overshoot
    integer :: nc
    ang = sv(2)
    nc = -1
    if (ang >= 999.0d0) nc = -1
    if (ang >= 0.0d0 .and. ang < 180.0d0) nc = 1
    if (ang < 0.0d0 .and. ang > -181.0d0) nc = 2
    if (ang >= 180.0d0 .and. ang < 361.0d0) nc = 0
    
    overshoot = 0.0d0
    if (nc >= 0 .and. nc <= 1 .and. ang < 999.0d0) then
        call rotate_stress_crack((/sig(1),sig(2),sig(3),0.0d0/), ang, sigp)
        if (sigp(2) > 2.07d6) overshoot = sigp(2) - 2.07d6
    end if
    
    write(*,'(A2,I4,F10.6,3ES11.3,F7.1,I4,F6.1,ES10.2)') &
        phase, step, strain0(3), sig(1), sig(2), sig(3), sv(2), nc, sv(6), overshoot
end subroutine

subroutine check_dir2(sig, sv, ft)
    real(dp), intent(in) :: sig(3), sv(NS), ft
    real(dp) :: sigp(4), ang
    integer :: nc
    ang = sv(2)
    nc = -1
    if (ang >= 0.0d0 .and. ang < 180.0d0) nc = 1
    if (ang < 0.0d0 .and. ang > -181.0d0) nc = 2
    if (ang >= 180.0d0 .and. ang < 361.0d0) nc = 0
    
    if (nc >= 0 .and. ang < 999.0d0) then
        call rotate_stress_crack((/sig(1),sig(2),sig(3),0.0d0/), ang, sigp)
        write(*,'(A,3ES10.2,A,F7.1,I3,A,ES10.2,L3)') &
            '  crack-coord sig: ', sigp(1), sigp(2), sigp(3), &
            '  ANG=', ang, nc, '  dir2_over_ft=', sigp(2)-ft, sigp(2)>ft
    end if
end subroutine
end program
