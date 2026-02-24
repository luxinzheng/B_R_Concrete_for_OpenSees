!===============================================================================
!  Wall symmetry test: verify force recovery after crack closing
!  This test mimics what happens at a single material point in a shear wall:
!    1. Apply shear (causes diagonal cracking)
!    2. Reverse shear (cracks close, new cracks should form symmetrically)
!  The force in both directions should be similar for a symmetric loading.
!===============================================================================
program test_wall_symmetry
    use adina_concrete_mod
    implicit none
    integer, parameter :: NP=37, NS=13
    real(dp) :: props(NP)
    integer :: npass, nfail
    
    npass = 0; nfail = 0
    call setup_props(props)
    
    write(*,'(A)') '========================================================'
    write(*,'(A)') '  Wall Symmetry & Crack Recovery Tests'
    write(*,'(A)') '========================================================'
    write(*,*)
    
    call test_cyclic_shear_symmetry(props, npass, nfail)
    call test_crack_close_stiffness_recovery(props, npass, nfail)
    call test_axial_plus_cyclic_shear(props, npass, nfail)
    call test_tangent_after_crack_close(props, npass, nfail)
    call test_multiple_reversal_forces(props, npass, nfail)
    
    write(*,*)
    write(*,'(A)') '========================================================'
    write(*,'(A,I4,A,I4,A,I4)') '  TOTAL: ', npass+nfail, ', PASS=', npass, ', FAIL=', nfail
    write(*,'(A)') '========================================================'
    if (nfail > 0) then
        write(*,'(A)') '  *** SOME TESTS FAILED ***'
    else
        write(*,'(A)') '  ALL TESTS PASSED'
    end if

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

subroutine check(name, cond, npass, nfail)
    character(len=*), intent(in) :: name
    logical, intent(in) :: cond
    integer, intent(inout) :: npass, nfail
    if (cond) then; npass=npass+1; else; nfail=nfail+1; write(*,'(A,A)') '  FAIL: ', name; end if
end subroutine

! ===================================================================
! Test 1: Cyclic shear symmetry
!   Apply +gamma, then -gamma. Peak shear stress should be similar.
! ===================================================================
subroutine test_cyclic_shear_symmetry(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    real(dp) :: peak_pos, peak_neg
    integer :: i
    logical :: ok
    
    write(*,'(A)') 'TEST S1: Cyclic shear symmetry'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ok = .true.
    peak_pos = 0.0d0; peak_neg = 0.0d0
    
    ! Forward shear to gamma = 2e-3
    do i = 1, 20
        strain1 = 0.0d0; strain1(3) = 1.0d-4 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (any(stress /= stress)) ok = .false.
        if (abs(stress(3)) > abs(peak_pos)) peak_pos = stress(3)
    end do
    
    ! Return to zero
    do i = 1, 20
        strain1 = 0.0d0; strain1(3) = 2.0d-3 - 1.0d-4 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (any(stress /= stress)) ok = .false.
    end do
    
    ! Reverse shear to gamma = -2e-3
    do i = 1, 20
        strain1 = 0.0d0; strain1(3) = -1.0d-4 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (any(stress /= stress)) ok = .false.
        if (abs(stress(3)) > abs(peak_neg)) peak_neg = stress(3)
    end do
    
    write(*,'(A,ES12.4,A,ES12.4)') '  Peak +shear: ', peak_pos, '  Peak -shear: ', peak_neg
    write(*,'(A,F6.1,A)') '  Ratio |neg/pos|: ', abs(peak_neg/peak_pos)*100.0d0, '%'
    
    call check('ShearSym: no NaN', ok, npass, nfail)
    call check('ShearSym: ratio > 50%', abs(peak_neg) > 0.5d0 * abs(peak_pos), npass, nfail)
    call check('ShearSym: ratio > 30%', abs(peak_neg) > 0.3d0 * abs(peak_pos), npass, nfail)
    write(*,*)
end subroutine

! ===================================================================
! Test 2: Stiffness recovery after crack closing
!   Crack in tension, close, then compress.
!   Compare compressive stiffness with virgin material.
! ===================================================================
subroutine test_crack_close_stiffness_recovery(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress_v(3), stress_c(3)
    real(dp) :: strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev_v(NS), statev_c(NS), tangent_v(3,3), tangent_c(3,3)
    integer :: i
    
    write(*,'(A)') 'TEST S2: Stiffness recovery after crack closing'
    
    ! Virgin material: compress to eps = -1e-3
    stress_v = 0.0d0; strain0 = 0.0d0; statev_v = 0.0d0
    do i = 1, 10
        strain1 = 0.0d0; strain1(1) = -1.0d-4 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress_v,strain0,strain1,dstrain,statev_v,tangent_v)
        strain0 = strain1
    end do
    
    ! Cracked material: crack, close, then compress to same eps
    stress_c = 0.0d0; strain0 = 0.0d0; statev_c = 0.0d0
    ! Crack
    strain1 = 0.0d0; strain1(1) = 3.0d-4
    dstrain = strain1
    call PSUMAT(NS,NP,props,stress_c,strain0,strain1,dstrain,statev_c,tangent_c)
    strain0 = strain1
    ! Close and compress
    do i = 1, 13
        strain1 = 0.0d0; strain1(1) = 3.0d-4 - 1.0d-4 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress_c,strain0,strain1,dstrain,statev_c,tangent_c)
        strain0 = strain1
    end do
    
    write(*,'(A,ES12.4,A,ES12.4)') '  Virgin sig_xx:  ', stress_v(1), '  C11: ', tangent_v(1,1)
    write(*,'(A,ES12.4,A,ES12.4)') '  Cracked sig_xx: ', stress_c(1), '  C11: ', tangent_c(1,1)
    write(*,'(A,F6.1,A)') '  Stress ratio: ', abs(stress_c(1)/stress_v(1))*100.0d0, '%'
    write(*,'(A,F6.1,A)') '  Tangent ratio: ', tangent_c(1,1)/tangent_v(1,1)*100.0d0, '%'
    
    call check('Recovery: compressive stress', stress_c(1) < -1.0d6, npass, nfail)
    call check('Recovery: tangent > 50% virgin', tangent_c(1,1) > 0.5d0*tangent_v(1,1), npass, nfail)
    write(*,*)
end subroutine

! ===================================================================
! Test 3: Axial compression + cyclic shear (wall-like loading)
!   Check that forces are reasonably symmetric
! ===================================================================
subroutine test_axial_plus_cyclic_shear(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    real(dp) :: peak_tau_pos, peak_tau_neg
    integer :: i
    logical :: ok
    
    write(*,'(A)') 'TEST S3: Axial compression + cyclic shear (wall path)'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ok = .true.
    peak_tau_pos = 0.0d0; peak_tau_neg = 0.0d0
    
    ! Apply light axial compression first (avoid crushing)
    do i = 1, 3
        strain1(1) = 0.0d0
        strain1(2) = -1.0d-4 * i
        strain1(3) = 0.0d0
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
    end do
    
    ! Cyclic shear with constant axial compression (smaller amplitude)
    ! Forward
    do i = 1, 15
        strain1(1) = 0.0d0
        strain1(2) = -3.0d-4
        strain1(3) = 1.0d-4 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (any(stress /= stress)) ok = .false.
        if (stress(3) > peak_tau_pos) peak_tau_pos = stress(3)
    end do
    ! Return
    do i = 1, 15
        strain1(1) = 0.0d0
        strain1(2) = -3.0d-4
        strain1(3) = 1.5d-3 - 1.0d-4 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (any(stress /= stress)) ok = .false.
    end do
    ! Reverse
    do i = 1, 15
        strain1(1) = 0.0d0
        strain1(2) = -3.0d-4
        strain1(3) = -1.0d-4 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (any(stress /= stress)) ok = .false.
        if (stress(3) < peak_tau_neg) peak_tau_neg = stress(3)
    end do
    
    write(*,'(A,ES12.4,A,ES12.4)') '  Peak +tau: ', peak_tau_pos, '  Peak -tau: ', peak_tau_neg
    write(*,'(A,F6.1,A)') '  Symmetry ratio: ', abs(peak_tau_neg/peak_tau_pos)*100.0d0, '%'
    write(*,'(A,F8.1,A,F6.1)') '  ANGLE: ', statev(2), '  PGRAV: ', statev(6)
    
    call check('WallPath: no NaN', ok, npass, nfail)
    call check('WallPath: not crushed', statev(6) < 50.0d0, npass, nfail)
    call check('WallPath: symmetry > 40%', abs(peak_tau_neg) > 0.4d0*abs(peak_tau_pos), npass, nfail)
    write(*,*)
end subroutine

! ===================================================================
! Test 4: Tangent after crack close should be close to elastic
! ===================================================================
subroutine test_tangent_after_crack_close(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    real(dp) :: E_elastic
    integer :: i
    
    write(*,'(A)') 'TEST S4: Tangent after crack closing'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    E_elastic = props(1) / (1.0d0 - props(2)**2)
    
    ! Crack
    strain1 = 0.0d0; strain1(1) = 3.0d-4
    dstrain = strain1
    call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
    strain0 = strain1
    
    ! Close crack
    do i = 1, 6
        strain1 = 0.0d0; strain1(1) = 3.0d-4 - 1.0d-4 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
    end do
    
    write(*,'(A,ES12.4,A,ES12.4)') '  C11 after close: ', tangent(1,1), '  Elastic C11: ', E_elastic
    write(*,'(A,F6.1,A)') '  Ratio: ', tangent(1,1)/E_elastic*100.0d0, '%'
    write(*,'(A,F8.1)') '  ANGLE: ', statev(2)
    
    call check('TangClose: C11 > 50% elastic', tangent(1,1) > 0.5d0*E_elastic, npass, nfail)
    call check('TangClose: C11 positive', tangent(1,1) > 0.0d0, npass, nfail)
    write(*,*)
end subroutine

! ===================================================================
! Test 5: Multiple reversal — check force levels at each peak
! ===================================================================
subroutine test_multiple_reversal_forces(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    real(dp) :: forces(6)
    integer :: i, cyc
    logical :: ok
    
    write(*,'(A)') 'TEST S5: Multiple reversal force levels'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ok = .true.
    
    do cyc = 1, 3
        ! Forward to +5e-4
        do i = 1, 10
            strain1 = 0.0d0; strain1(1) = 5.0d-5 * i
            dstrain = strain1 - strain0
            call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
            strain0 = strain1
            if (any(stress /= stress)) ok = .false.
        end do
        forces(2*cyc-1) = stress(1)
        
        ! Reverse to -1e-3
        do i = 1, 15
            strain1 = 0.0d0; strain1(1) = 5.0d-4 - 1.0d-4 * i
            dstrain = strain1 - strain0
            call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
            strain0 = strain1
            if (any(stress /= stress)) ok = .false.
        end do
        forces(2*cyc) = stress(1)
        
        ! Return to start
        do i = 1, 10
            strain1 = 0.0d0; strain1(1) = -1.0d-3 + 1.5d-4 * i
            dstrain = strain1 - strain0
            call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
            strain0 = strain1
        end do
    end do
    
    write(*,'(A)') '  Peak forces at each reversal:'
    do cyc = 1, 3
        write(*,'(A,I1,A,ES12.4,A,ES12.4)') '    Cycle ', cyc, &
            ': +peak=', forces(2*cyc-1), '  -peak=', forces(2*cyc)
    end do
    
    call check('MultiRev: no NaN', ok, npass, nfail)
    call check('MultiRev: compressive force > 5 MPa', abs(forces(2)) > 5.0d6, npass, nfail)
    call check('MultiRev: stress bounded', all(abs(forces) < 1.0d9), npass, nfail)
    write(*,*)
end subroutine

end program test_wall_symmetry
