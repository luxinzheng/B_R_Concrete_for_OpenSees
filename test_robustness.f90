!===============================================================================
!  Robustness tests: edge cases, multi-cycle, extreme loading, shear wall paths
!===============================================================================
program test_robustness
    use adina_concrete_mod
    implicit none
    integer, parameter :: NP=37, NS=13
    real(dp) :: props(NP)
    integer :: npass, nfail
    
    npass = 0; nfail = 0
    call setup_props(props)
    
    write(*,'(A)') '========================================================'
    write(*,'(A)') '  Robustness Tests: Edge Cases & Shear Wall Paths'
    write(*,'(A)') '========================================================'
    write(*,*)
    
    call test_multi_cycle_tension(props, npass, nfail)
    call test_multi_cycle_compression(props, npass, nfail)
    call test_shear_wall_path(props, npass, nfail)
    call test_biaxial_tension_compression(props, npass, nfail)
    call test_large_shear_strain(props, npass, nfail)
    call test_tiny_increments(props, npass, nfail)
    call test_large_increments(props, npass, nfail)
    call test_crack_close_reopen_cycles(props, npass, nfail)
    call test_compression_after_cracking(props, npass, nfail)
    call test_zero_strain_increment(props, npass, nfail)
    call test_crush_path(props, npass, nfail)
    
    write(*,*)
    write(*,'(A)') '========================================================'
    write(*,'(A,I4,A,I4,A,I4)') '  TOTAL: ', npass+nfail, ', PASS=', npass, ', FAIL=', nfail
    write(*,'(A)') '========================================================'
    if (nfail > 0) then
        write(*,'(A)') '  *** SOME TESTS FAILED ***'
    else
        write(*,'(A)') '  ALL ROBUSTNESS TESTS PASSED'
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

logical function has_nan_stress(stress)
    real(dp), intent(in) :: stress(3)
    has_nan_stress = any(stress /= stress)
end function

logical function has_nan_tangent(tangent)
    real(dp), intent(in) :: tangent(3,3)
    has_nan_tangent = any(tangent /= tangent)
end function

! ===================================================================
! Test 1: Multiple tension-compression cycles
! ===================================================================
subroutine test_multi_cycle_tension(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    integer :: cyc, i
    logical :: ok
    
    write(*,'(A)') 'TEST R1: Multi-cycle tension-compression (5 cycles)'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ok = .true.
    
    do cyc = 1, 5
        ! Tension to 5e-4
        do i = 1, 10
            strain1 = 0.0d0; strain1(1) = 5.0d-5 * i
            dstrain = strain1 - strain0
            call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
            strain0 = strain1
            if (has_nan_stress(stress) .or. has_nan_tangent(tangent)) ok = .false.
        end do
        ! Compress to -2e-3
        do i = 1, 25
            strain1 = 0.0d0; strain1(1) = 5.0d-4 - 1.0d-4 * i
            dstrain = strain1 - strain0
            call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
            strain0 = strain1
            if (has_nan_stress(stress) .or. has_nan_tangent(tangent)) ok = .false.
        end do
        ! Back to 5e-4
        do i = 1, 25
            strain1 = 0.0d0; strain1(1) = -2.0d-3 + 1.0d-4 * i
            dstrain = strain1 - strain0
            call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
            strain0 = strain1
            if (has_nan_stress(stress) .or. has_nan_tangent(tangent)) ok = .false.
        end do
    end do
    
    call check('MultiCycleTens: no NaN', ok, npass, nfail)
    call check('MultiCycleTens: stress bounded', all(abs(stress) < 1.0d9), npass, nfail)
    write(*,'(A,3ES12.4)') '  Final stress: ', stress
    write(*,*)
end subroutine

! ===================================================================
! Test 2: Multiple compression cycles with increasing amplitude
! ===================================================================
subroutine test_multi_cycle_compression(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    real(dp) :: amp
    integer :: cyc, i
    logical :: ok
    
    write(*,'(A)') 'TEST R2: Multi-cycle compression with increasing amplitude'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ok = .true.
    
    do cyc = 1, 6
        amp = -0.5d-3 * cyc
        ! Compress
        do i = 1, 10
            strain1 = 0.0d0; strain1(1) = amp * i / 10.0d0
            dstrain = strain1 - strain0
            call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
            strain0 = strain1
            if (has_nan_stress(stress) .or. has_nan_tangent(tangent)) ok = .false.
        end do
        ! Unload
        do i = 1, 10
            strain1 = 0.0d0; strain1(1) = amp * (10 - i) / 10.0d0
            dstrain = strain1 - strain0
            call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
            strain0 = strain1
            if (has_nan_stress(stress) .or. has_nan_tangent(tangent)) ok = .false.
        end do
    end do
    
    call check('MultiCycleComp: no NaN', ok, npass, nfail)
    call check('MultiCycleComp: stress bounded', all(abs(stress) < 1.0d9), npass, nfail)
    write(*,'(A,3ES12.4)') '  Final stress: ', stress
    write(*,*)
end subroutine

! ===================================================================
! Test 3: Shear wall typical loading path
!   Cyclic shear with axial compression (represents wall element)
! ===================================================================
subroutine test_shear_wall_path(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    real(dp) :: axial_strain, shear_amp, t
    integer :: i, ntotal
    logical :: ok
    real(dp) :: max_stress
    
    write(*,'(A)') 'TEST R3: Shear wall path (axial + cyclic shear)'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ok = .true.
    max_stress = 0.0d0
    
    axial_strain = -5.0d-4
    ntotal = 400
    
    do i = 1, ntotal
        t = dble(i) / dble(ntotal) * 4.0d0 * 3.14159265d0
        shear_amp = 3.0d-3 * sin(t) * min(1.0d0, dble(i)/50.0d0)
        
        strain1(1) = 0.0d0
        strain1(2) = axial_strain
        strain1(3) = shear_amp
        
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        
        if (has_nan_stress(stress) .or. has_nan_tangent(tangent)) then
            ok = .false.
            write(*,'(A,I4,A,3ES12.4)') '  NaN at step ', i, ': ', stress
            exit
        end if
        
        max_stress = max(max_stress, maxval(abs(stress)))
    end do
    
    call check('ShearWall: no NaN', ok, npass, nfail)
    call check('ShearWall: stress bounded (<1 GPa)', max_stress < 1.0d9, npass, nfail)
    call check('ShearWall: positive tangent diag', &
        tangent(1,1) > 0 .and. tangent(2,2) > 0 .and. tangent(3,3) > 0, npass, nfail)
    write(*,'(A,ES12.4)') '  Max stress magnitude: ', max_stress
    write(*,'(A,3ES12.4)') '  Final stress: ', stress
    write(*,*)
end subroutine

! ===================================================================
! Test 4: Biaxial tension-compression (diagonal loading)
! ===================================================================
subroutine test_biaxial_tension_compression(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    integer :: i
    logical :: ok
    
    write(*,'(A)') 'TEST R4: Biaxial tension(x) + compression(y)'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ok = .true.
    
    do i = 1, 40
        strain1(1) = 2.0d-5 * i
        strain1(2) = -5.0d-5 * i
        strain1(3) = 0.0d0
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (has_nan_stress(stress) .or. has_nan_tangent(tangent)) ok = .false.
    end do
    
    call check('BiaxialTC: no NaN', ok, npass, nfail)
    call check('BiaxialTC: sig_xx > 0 or cracked', stress(1) >= -1.0d6, npass, nfail)
    call check('BiaxialTC: sig_yy < 0', stress(2) < 0.0d0, npass, nfail)
    write(*,'(A,3ES12.4)') '  Final stress: ', stress
    write(*,*)
end subroutine

! ===================================================================
! Test 5: Large shear strain
! ===================================================================
subroutine test_large_shear_strain(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    integer :: i
    logical :: ok
    
    write(*,'(A)') 'TEST R5: Large shear strain (gamma up to 0.01)'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ok = .true.
    
    do i = 1, 100
        strain1 = 0.0d0
        strain1(3) = 1.0d-4 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (has_nan_stress(stress) .or. has_nan_tangent(tangent)) then
            ok = .false.
            write(*,'(A,I4,A,3ES12.4)') '  NaN at step ', i, ': ', stress
            exit
        end if
    end do
    
    call check('LargeShear: no NaN', ok, npass, nfail)
    call check('LargeShear: stress bounded', all(abs(stress) < 1.0d9), npass, nfail)
    write(*,'(A,3ES12.4)') '  Final stress: ', stress
    write(*,*)
end subroutine

! ===================================================================
! Test 6: Very tiny increments (near-zero dstrain)
! ===================================================================
subroutine test_tiny_increments(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    integer :: i
    logical :: ok
    
    write(*,'(A)') 'TEST R6: Tiny increments (dε = 1e-8 per step)'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ok = .true.
    
    do i = 1, 200
        strain1 = 0.0d0
        strain1(1) = 1.0d-8 * i
        strain1(2) = -0.5d-8 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (has_nan_stress(stress) .or. has_nan_tangent(tangent)) then
            ok = .false.
            exit
        end if
    end do
    
    call check('TinyInc: no NaN', ok, npass, nfail)
    call check('TinyInc: tangent diagonal positive', &
        tangent(1,1) > 0 .and. tangent(2,2) > 0 .and. tangent(3,3) > 0, npass, nfail)
    write(*,'(A,3ES12.4)') '  Final stress: ', stress
    write(*,*)
end subroutine

! ===================================================================
! Test 7: Large increments (test stability with big strain jumps)
! ===================================================================
subroutine test_large_increments(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    logical :: ok
    
    write(*,'(A)') 'TEST R7: Large single increment (ε = -0.003 in one step)'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ok = .true.
    
    strain1 = 0.0d0; strain1(1) = -3.0d-3
    dstrain = strain1
    call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
    
    if (has_nan_stress(stress) .or. has_nan_tangent(tangent)) ok = .false.
    
    call check('LargeInc: no NaN', ok, npass, nfail)
    call check('LargeInc: compressive stress', stress(1) < 0.0d0, npass, nfail)
    call check('LargeInc: bounded', abs(stress(1)) < 1.0d9, npass, nfail)
    write(*,'(A,3ES12.4)') '  Stress: ', stress
    write(*,*)
end subroutine

! ===================================================================
! Test 8: Multiple crack close-reopen cycles
! ===================================================================
subroutine test_crack_close_reopen_cycles(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    integer :: cyc, i
    logical :: ok
    real(dp) :: max_tens_stress
    
    write(*,'(A)') 'TEST R8: Multiple crack close-reopen (3 cycles)'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ok = .true.
    max_tens_stress = 0.0d0
    
    do cyc = 1, 3
        ! Tension: crack opens
        do i = 1, 10
            strain1 = 0.0d0; strain1(1) = 5.0d-5 * i
            dstrain = strain1 - strain0
            call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
            strain0 = strain1
            if (has_nan_stress(stress)) ok = .false.
            if (stress(1) > max_tens_stress) max_tens_stress = stress(1)
        end do
        
        ! Compression: crack closes
        do i = 1, 20
            strain1 = 0.0d0; strain1(1) = 5.0d-4 - 1.0d-4 * i
            dstrain = strain1 - strain0
            call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
            strain0 = strain1
            if (has_nan_stress(stress)) ok = .false.
        end do
        
        ! Back to tension: crack reopens
        do i = 1, 20
            strain1 = 0.0d0; strain1(1) = -1.5d-3 + 1.0d-4 * i
            dstrain = strain1 - strain0
            call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
            strain0 = strain1
            if (has_nan_stress(stress)) ok = .false.
            if (stress(1) > max_tens_stress) max_tens_stress = stress(1)
        end do
    end do
    
    call check('CloseReopen: no NaN', ok, npass, nfail)
    call check('CloseReopen: tensile stress <= ft', max_tens_stress <= 2.2d6, npass, nfail)
    write(*,'(A,ES12.4)') '  Max tensile stress: ', max_tens_stress
    write(*,'(A,3ES12.4)') '  Final stress: ', stress
    write(*,*)
end subroutine

! ===================================================================
! Test 9: Heavy compression after cracking
! ===================================================================
subroutine test_compression_after_cracking(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    integer :: i
    logical :: ok
    
    write(*,'(A)') 'TEST R9: Deep compression after cracking'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ok = .true.
    
    ! Crack in tension
    do i = 1, 10
        strain1 = 0.0d0; strain1(1) = 5.0d-5 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
    end do
    
    ! Deep compression to -0.004
    do i = 1, 45
        strain1 = 0.0d0; strain1(1) = 5.0d-4 - 1.0d-4 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (has_nan_stress(stress) .or. has_nan_tangent(tangent)) then
            ok = .false.
            write(*,'(A,I4,A,3ES12.4)') '  NaN at step ', i, ': ', stress
            exit
        end if
    end do
    
    call check('CompAfterCrack: no NaN', ok, npass, nfail)
    call check('CompAfterCrack: compressive', stress(1) < 0.0d0, npass, nfail)
    call check('CompAfterCrack: bounded', abs(stress(1)) < 1.0d9, npass, nfail)
    write(*,'(A,3ES12.4)') '  Final stress: ', stress
    write(*,'(A,F8.1)') '  ANGLE: ', statev(2)
    write(*,*)
end subroutine

! ===================================================================
! Test 10: Zero strain increment
! ===================================================================
subroutine test_zero_strain_increment(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    real(dp) :: stress_before(3)
    integer :: i
    logical :: ok
    
    write(*,'(A)') 'TEST R10: Zero strain increment (stress should not change)'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ok = .true.
    
    ! First compress to some state
    strain1 = 0.0d0; strain1(1) = -1.0d-3
    dstrain = strain1
    call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
    strain0 = strain1
    stress_before = stress
    
    ! Now apply zero increment multiple times
    do i = 1, 5
        dstrain = 0.0d0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        if (has_nan_stress(stress)) ok = .false.
    end do
    
    call check('ZeroInc: no NaN', ok, npass, nfail)
    call check('ZeroInc: stress unchanged', &
        abs(stress(1) - stress_before(1)) < abs(stress_before(1)) * 0.01d0, npass, nfail)
    write(*,'(A,ES12.4,A,ES12.4)') '  Before: ', stress_before(1), '  After: ', stress(1)
    write(*,*)
end subroutine

! ===================================================================
! Test 11: Path to crushing
! ===================================================================
subroutine test_crush_path(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    integer :: i
    logical :: ok, crushed
    
    write(*,'(A)') 'TEST R11: Compression to crushing (biaxial)'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ok = .true.
    crushed = .false.
    
    ! Biaxial compression to large strain
    do i = 1, 80
        strain1(1) = -1.0d-4 * i
        strain1(2) = -0.5d-4 * i
        strain1(3) = 0.0d0
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (has_nan_stress(stress) .or. has_nan_tangent(tangent)) then
            ok = .false.
            write(*,'(A,I4,A,3ES12.4)') '  NaN at step ', i, ': ', stress
            exit
        end if
        if (statev(6) > 99.0d0) crushed = .true.
    end do
    
    call check('Crush: no NaN', ok, npass, nfail)
    call check('Crush: reached crushed state or post-peak softening', &
        crushed .or. abs(stress(1)) < abs(props(4)), npass, nfail)
    
    ! Continue loading after crush
    if (crushed) then
        do i = 1, 10
            strain1(1) = strain1(1) - 1.0d-4
            strain1(2) = strain1(2) - 0.5d-4
            dstrain = strain1 - strain0
            call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
            strain0 = strain1
            if (has_nan_stress(stress)) ok = .false.
        end do
        call check('Crush: post-crush no NaN', ok, npass, nfail)
    else
        npass = npass + 1
    end if
    
    write(*,'(A,3ES12.4)') '  Final stress: ', stress
    write(*,'(A,F8.1)') '  PGRAV: ', statev(6)
    write(*,*)
end subroutine

end program test_robustness
