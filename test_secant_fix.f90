!===============================================================================
! Comprehensive standalone test for secant stress correction in forumat.f90
! Tests: uniaxial, cyclic, crack open/close, ML scenario, biaxial
!===============================================================================
program test_secant_fix
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    integer :: n_pass, n_fail, n_total
    n_pass = 0; n_fail = 0; n_total = 0
    
    write(*,'(A)') '======================================================'
    write(*,'(A)') ' Secant Stress Correction: Comprehensive Test Suite'
    write(*,'(A)') '======================================================'
    
    call test_uniaxial_compression(n_pass, n_fail, n_total)
    call test_uniaxial_tension(n_pass, n_fail, n_total)
    call test_uniaxial_cyclic(n_pass, n_fail, n_total)
    call test_crack_close_reopen(n_pass, n_fail, n_total)
    call test_biaxial_compression(n_pass, n_fail, n_total)
    call test_ml_scenario_drift(n_pass, n_fail, n_total)
    call test_shear_cyclic(n_pass, n_fail, n_total)
    call test_endurance_500cycles(n_pass, n_fail, n_total)
    
    write(*,'(/,A)') '======================================================'
    write(*,'(A,I3,A,I3,A,I3)') ' RESULTS: ', n_pass, ' PASS / ', n_fail, ' FAIL / ', n_total, ' TOTAL'
    if (n_fail == 0) then
        write(*,'(A)') ' >>> ALL TESTS PASSED <<<'
    else
        write(*,'(A)') ' >>> SOME TESTS FAILED <<<'
    end if
    write(*,'(A)') '======================================================'
    
contains

subroutine init_props_ml(props)
    real(dp), intent(out) :: props(37)
    props = 0.0d0
    props(1)  = 2.40254d10   ! E0
    props(2)  = 0.2d0        ! VNU
    props(3)  = 2.0d6        ! SIGMAT (ft)
    props(4)  = -2.58d7      ! SIGMAC (fc)
    props(5)  = -0.003d0     ! EPSC
    props(6)  = -5.16d6      ! SIGMAU
    props(7)  = -0.021d0     ! EPSU
    props(8)  = 0.0001d0     ! STIFAC
    props(9)  = 0.5d0        ! SHEFAC
    props(10) = 0.75d0       ! BETA
    props(11) = 1.0d0        ! GAMA
    props(12) = 0.7d0        ! RKAPA
    props(13) = 0.12d0       ! ALFA
    props(14:19) = (/ 0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.2d0 /)
    props(20:25) = (/ 1.0d0, 1.4d0, 1.7d0, 2.2d0, 2.5d0, 2.8d0 /)
    props(26:31) = (/ 1.3d0, 1.5d0, 2.0d0, 2.3d0, 2.7d0, 3.2d0 /)
    props(32:37) = (/ 1.25d0, 1.45d0, 1.95d0, 2.25d0, 2.65d0, 3.15d0 /)
end subroutine

subroutine call_psumat(props, stress, strain0, strain1, statev, tangent)
    real(dp), intent(in)    :: props(37)
    real(dp), intent(inout) :: stress(3), statev(40)
    real(dp), intent(in)    :: strain0(3), strain1(3)
    real(dp), intent(out)   :: tangent(3,3)
    real(dp) :: dstrain(3)
    dstrain = strain1 - strain0
    call PSUMAT(40, 37, props, stress, strain0, strain1, dstrain, statev, tangent)
end subroutine

subroutine check(name, condition, n_pass, n_fail, n_total)
    character(len=*), intent(in) :: name
    logical, intent(in) :: condition
    integer, intent(inout) :: n_pass, n_fail, n_total
    n_total = n_total + 1
    if (condition) then
        n_pass = n_pass + 1
        write(*,'(A,A,A)') '  [PASS] ', name, ''
    else
        n_fail = n_fail + 1
        write(*,'(A,A,A)') '  [FAIL] ', name, ''
    end if
end subroutine

! =================================================================
! Test 1: Uniaxial compression — stress should follow Saenz curve
! =================================================================
subroutine test_uniaxial_compression(np, nf, nt)
    integer, intent(inout) :: np, nf, nt
    real(dp) :: props(37), stress(3), statev(40), tangent(3,3)
    real(dp) :: strain0(3), strain1(3)
    real(dp) :: eps, sig_xx, E0, fc, epsc, sig_saenz
    integer :: i
    
    write(*,'(/,A)') '--- Test 1: Uniaxial Compression (Saenz curve) ---'
    call init_props_ml(props)
    stress = 0.0d0; statev = 0.0d0; strain0 = 0.0d0
    E0 = props(1); fc = props(4); epsc = props(5)
    
    do i = 1, 30
        eps = -1.0d-4 * dble(i)
        strain1 = (/ eps, 0.0d0, 0.0d0 /)
        call call_psumat(props, stress, strain0, strain1, statev, tangent)
        strain0 = strain1
    end do
    
    sig_xx = stress(1)
    call check('Comp stress negative', sig_xx < 0.0d0, np, nf, nt)
    call check('Comp stress bounded |sig|<|fc|*1.5', abs(sig_xx) < abs(fc)*1.5d0, np, nf, nt)
    call check('No NaN in stress', stress(1)==stress(1) .and. stress(2)==stress(2), np, nf, nt)
    write(*,'(A,ES12.4,A,ES12.4)') '    sig_xx = ', sig_xx, '  fc = ', fc
end subroutine

! =================================================================
! Test 2: Uniaxial tension — should crack at ft
! =================================================================
subroutine test_uniaxial_tension(np, nf, nt)
    integer, intent(inout) :: np, nf, nt
    real(dp) :: props(37), stress(3), statev(40), tangent(3,3)
    real(dp) :: strain0(3), strain1(3)
    real(dp) :: eps, ft
    integer :: i
    logical :: cracked
    
    write(*,'(/,A)') '--- Test 2: Uniaxial Tension (cracking) ---'
    call init_props_ml(props)
    stress = 0.0d0; statev = 0.0d0; strain0 = 0.0d0
    ft = props(3)
    cracked = .false.
    
    do i = 1, 20
        eps = 5.0d-6 * dble(i)
        strain1 = (/ eps, 0.0d0, 0.0d0 /)
        call call_psumat(props, stress, strain0, strain1, statev, tangent)
        strain0 = strain1
        if (statev(2) < 361.0d0) cracked = .true.
    end do
    
    call check('Tension cracking occurred', cracked, np, nf, nt)
    call check('Post-crack stress <= ft', stress(1) <= ft * 1.01d0, np, nf, nt)
    call check('No NaN', stress(1)==stress(1), np, nf, nt)
    write(*,'(A,ES12.4,A,L1)') '    sig_xx = ', stress(1), '  cracked = ', cracked
end subroutine

! =================================================================
! Test 3: Uniaxial cyclic (compress → unload → tension → compress)
! This is the KEY test for stress drift prevention
! =================================================================
subroutine test_uniaxial_cyclic(np, nf, nt)
    integer, intent(inout) :: np, nf, nt
    real(dp) :: props(37), stress(3), statev(40), tangent(3,3)
    real(dp) :: strain0(3), strain1(3)
    real(dp) :: eps, fc, ft, E0, max_sig, min_sig
    integer :: i, cycle
    
    write(*,'(/,A)') '--- Test 3: Uniaxial Cyclic (drift prevention) ---'
    call init_props_ml(props)
    stress = 0.0d0; statev = 0.0d0; strain0 = 0.0d0
    fc = props(4); ft = props(3); E0 = props(1)
    max_sig = -1.0d30; min_sig = 1.0d30
    
    do cycle = 1, 20
        ! Compress to -0.002
        do i = 1, 20
            eps = -1.0d-4 * dble(i)
            strain1 = (/ eps, 0.0d0, 0.0d0 /)
            call call_psumat(props, stress, strain0, strain1, statev, tangent)
            strain0 = strain1
            if (stress(1) > max_sig) max_sig = stress(1)
            if (stress(1) < min_sig) min_sig = stress(1)
        end do
        ! Unload to +0.0001 (tension, may crack)
        do i = 19, -1, -1
            eps = -1.0d-4 * dble(i)
            strain1 = (/ eps, 0.0d0, 0.0d0 /)
            call call_psumat(props, stress, strain0, strain1, statev, tangent)
            strain0 = strain1
            if (stress(1) > max_sig) max_sig = stress(1)
            if (stress(1) < min_sig) min_sig = stress(1)
        end do
    end do
    
    call check('Max comp stress bounded', min_sig > fc * 2.0d0, np, nf, nt)
    call check('Max tens stress bounded', max_sig < ft * 2.0d0, np, nf, nt)
    call check('No NaN', stress(1)==stress(1), np, nf, nt)
    write(*,'(A,ES12.4,A,ES12.4)') '    min_sig = ', min_sig, '  max_sig = ', max_sig
    write(*,'(A,ES12.4,A,ES12.4)') '    fc = ', fc, '  ft = ', ft
end subroutine

! =================================================================
! Test 4: Crack open → close → reopen cycle (Multi-layer Shell pattern)
! =================================================================
subroutine test_crack_close_reopen(np, nf, nt)
    integer, intent(inout) :: np, nf, nt
    real(dp) :: props(37), stress(3), statev(40), tangent(3,3)
    real(dp) :: strain0(3), strain1(3)
    real(dp) :: eps_x, eps_y, fc, ft, max_abs_sig
    integer :: i, cycle
    
    write(*,'(/,A)') '--- Test 4: Crack Open/Close/Reopen Cycles ---'
    call init_props_ml(props)
    stress = 0.0d0; statev = 0.0d0; strain0 = 0.0d0
    fc = props(4); ft = props(3)
    max_abs_sig = 0.0d0
    
    do cycle = 1, 30
        ! Phase 1: Tension (crack opens)
        do i = 1, 10
            eps_x = 2.0d-5 * dble(i)
            eps_y = -0.2d0 * eps_x
            strain1 = (/ eps_x, eps_y, 0.0d0 /)
            call call_psumat(props, stress, strain0, strain1, statev, tangent)
            strain0 = strain1
        end do
        ! Phase 2: Compression (crack closes)
        do i = 10, -10, -1
            eps_x = 2.0d-5 * dble(i)
            eps_y = -0.2d0 * eps_x
            strain1 = (/ eps_x, eps_y, 1.0d-5 /)
            call call_psumat(props, stress, strain0, strain1, statev, tangent)
            strain0 = strain1
            if (abs(stress(1)) > max_abs_sig) max_abs_sig = abs(stress(1))
            if (abs(stress(2)) > max_abs_sig) max_abs_sig = abs(stress(2))
        end do
        ! Phase 3: Back to tension (crack reopens)
        do i = -10, 0
            eps_x = 2.0d-5 * dble(i)
            eps_y = -0.2d0 * eps_x
            strain1 = (/ eps_x, eps_y, 0.0d0 /)
            call call_psumat(props, stress, strain0, strain1, statev, tangent)
            strain0 = strain1
        end do
    end do
    
    call check('Max |stress| bounded (< 5*|fc|)', max_abs_sig < abs(fc)*5.0d0, np, nf, nt)
    call check('No NaN', all(stress == stress), np, nf, nt)
    write(*,'(A,ES12.4,A,ES12.4)') '    max_abs_sig = ', max_abs_sig, '  |fc| = ', abs(fc)
end subroutine

! =================================================================
! Test 5: Biaxial compression — enhanced strength
! =================================================================
subroutine test_biaxial_compression(np, nf, nt)
    integer, intent(inout) :: np, nf, nt
    real(dp) :: props(37), stress(3), statev(40), tangent(3,3)
    real(dp) :: strain0(3), strain1(3)
    real(dp) :: eps, fc
    integer :: i
    
    write(*,'(/,A)') '--- Test 5: Biaxial Compression ---'
    call init_props_ml(props)
    stress = 0.0d0; statev = 0.0d0; strain0 = 0.0d0
    fc = props(4)
    
    do i = 1, 30
        eps = -1.0d-4 * dble(i)
        strain1 = (/ eps, eps, 0.0d0 /)
        call call_psumat(props, stress, strain0, strain1, statev, tangent)
        strain0 = strain1
    end do
    
    call check('Biaxial stress < fc (enhanced)', stress(1) < fc, np, nf, nt)
    call check('Biaxial sig_xx ≈ sig_yy', abs(stress(1)-stress(2)) < abs(fc)*0.1d0, np, nf, nt)
    call check('No NaN', all(stress == stress), np, nf, nt)
    write(*,'(A,ES12.4,A,ES12.4)') '    sig_xx = ', stress(1), '  sig_yy = ', stress(2)
end subroutine

! =================================================================
! Test 6: Multi-layer_Shell scenario — large cyclic with shear
! This simulates the DisplacementControl-induced strain patterns
! =================================================================
subroutine test_ml_scenario_drift(np, nf, nt)
    integer, intent(inout) :: np, nf, nt
    real(dp) :: props(37), stress(3), statev(40), tangent(3,3)
    real(dp) :: strain0(3), strain1(3)
    real(dp) :: fc, ft, amp, eps_x, eps_y, gam_xy
    real(dp) :: max_sig, min_sig
    integer :: block, i, nstep
    
    write(*,'(/,A)') '--- Test 6: ML Scenario (cyclic drift test) ---'
    call init_props_ml(props)
    stress = 0.0d0; statev = 0.0d0; strain0 = 0.0d0
    fc = props(4); ft = props(3)
    max_sig = -1.0d30; min_sig = 1.0d30
    
    do block = 1, 11
        amp = 1.0d-3 * dble(block)
        nstep = 10
        
        ! Forward (tension dominant)
        do i = 1, nstep
            eps_x = amp * dble(i) / dble(nstep)
            eps_y = -0.15d0 * eps_x
            gam_xy = 0.3d0 * eps_x
            strain1 = (/ eps_x, eps_y, gam_xy /)
            call call_psumat(props, stress, strain0, strain1, statev, tangent)
            strain0 = strain1
            if (stress(1) > max_sig) max_sig = stress(1)
            if (stress(1) < min_sig) min_sig = stress(1)
        end do
        ! Reverse (compression dominant)
        do i = nstep, -nstep, -1
            eps_x = amp * dble(i) / dble(nstep)
            eps_y = -0.15d0 * eps_x
            gam_xy = 0.3d0 * eps_x
            strain1 = (/ eps_x, eps_y, gam_xy /)
            call call_psumat(props, stress, strain0, strain1, statev, tangent)
            strain0 = strain1
            if (stress(1) > max_sig) max_sig = stress(1)
            if (stress(1) < min_sig) min_sig = stress(1)
        end do
        ! Return to zero
        do i = -nstep, 0
            eps_x = amp * dble(i) / dble(nstep)
            eps_y = -0.15d0 * eps_x
            gam_xy = 0.3d0 * eps_x
            strain1 = (/ eps_x, eps_y, gam_xy /)
            call call_psumat(props, stress, strain0, strain1, statev, tangent)
            strain0 = strain1
            if (stress(1) > max_sig) max_sig = stress(1)
            if (stress(1) < min_sig) min_sig = stress(1)
        end do
    end do
    
    call check('ML: min_sig bounded (>10*fc?)', min_sig > fc * 10.0d0, np, nf, nt)
    call check('ML: max_sig bounded (<10*ft?)', max_sig < ft * 10.0d0, np, nf, nt)
    call check('ML: no NaN', all(stress == stress), np, nf, nt)
    write(*,'(A,ES12.4,A,ES12.4)') '    min_sig = ', min_sig, '  max_sig = ', max_sig
    write(*,'(A,ES12.4,A,ES12.4)') '    10*fc = ', fc*10.0d0, '  10*ft = ', ft*10.0d0
end subroutine

! =================================================================
! Test 7: Pure shear cyclic
! =================================================================
subroutine test_shear_cyclic(np, nf, nt)
    integer, intent(inout) :: np, nf, nt
    real(dp) :: props(37), stress(3), statev(40), tangent(3,3)
    real(dp) :: strain0(3), strain1(3)
    real(dp) :: gam, fc, max_tau
    integer :: cycle, i
    
    write(*,'(/,A)') '--- Test 7: Pure Shear Cyclic ---'
    call init_props_ml(props)
    stress = 0.0d0; statev = 0.0d0; strain0 = 0.0d0
    fc = props(4); max_tau = 0.0d0
    
    do cycle = 1, 50
        do i = 0, 20
            gam = 5.0d-4 * dble(i)
            strain1 = (/ 0.0d0, 0.0d0, gam /)
            call call_psumat(props, stress, strain0, strain1, statev, tangent)
            strain0 = strain1
            if (abs(stress(3)) > max_tau) max_tau = abs(stress(3))
        end do
        do i = 20, -20, -1
            gam = 5.0d-4 * dble(i)
            strain1 = (/ 0.0d0, 0.0d0, gam /)
            call call_psumat(props, stress, strain0, strain1, statev, tangent)
            strain0 = strain1
            if (abs(stress(3)) > max_tau) max_tau = abs(stress(3))
        end do
        do i = -20, 0
            gam = 5.0d-4 * dble(i)
            strain1 = (/ 0.0d0, 0.0d0, gam /)
            call call_psumat(props, stress, strain0, strain1, statev, tangent)
            strain0 = strain1
            if (abs(stress(3)) > max_tau) max_tau = abs(stress(3))
        end do
    end do
    
    call check('Shear: max |tau| bounded (< |fc|)', max_tau < abs(fc), np, nf, nt)
    call check('Shear: no NaN', all(stress == stress), np, nf, nt)
    write(*,'(A,ES12.4)') '    max_tau = ', max_tau
end subroutine

! =================================================================
! Test 8: 500-cycle endurance at ML amplitude
! =================================================================
subroutine test_endurance_500cycles(np, nf, nt)
    integer, intent(inout) :: np, nf, nt
    real(dp) :: props(37), stress(3), statev(40), tangent(3,3)
    real(dp) :: strain0(3), strain1(3)
    real(dp) :: fc, ft, amp, eps_x, gam
    real(dp) :: max_abs
    integer :: cycle, i
    
    write(*,'(/,A)') '--- Test 8: 500-Cycle Endurance ---'
    call init_props_ml(props)
    stress = 0.0d0; statev = 0.0d0; strain0 = 0.0d0
    fc = props(4); ft = props(3)
    max_abs = 0.0d0
    amp = 5.0d-3
    
    do cycle = 1, 500
        ! Forward
        do i = 1, 10
            eps_x = amp * dble(i) / 10.0d0
            gam = 0.2d0 * eps_x
            strain1 = (/ eps_x, -0.1d0*eps_x, gam /)
            call call_psumat(props, stress, strain0, strain1, statev, tangent)
            strain0 = strain1
            if (abs(stress(1)) > max_abs) max_abs = abs(stress(1))
            if (abs(stress(2)) > max_abs) max_abs = abs(stress(2))
            if (abs(stress(3)) > max_abs) max_abs = abs(stress(3))
        end do
        ! Reverse
        do i = 10, -10, -1
            eps_x = amp * dble(i) / 10.0d0
            gam = 0.2d0 * eps_x
            strain1 = (/ eps_x, -0.1d0*eps_x, gam /)
            call call_psumat(props, stress, strain0, strain1, statev, tangent)
            strain0 = strain1
            if (abs(stress(1)) > max_abs) max_abs = abs(stress(1))
            if (abs(stress(2)) > max_abs) max_abs = abs(stress(2))
            if (abs(stress(3)) > max_abs) max_abs = abs(stress(3))
        end do
        ! Return
        do i = -10, 0
            eps_x = amp * dble(i) / 10.0d0
            gam = 0.2d0 * eps_x
            strain1 = (/ eps_x, -0.1d0*eps_x, gam /)
            call call_psumat(props, stress, strain0, strain1, statev, tangent)
            strain0 = strain1
            if (abs(stress(1)) > max_abs) max_abs = abs(stress(1))
            if (abs(stress(2)) > max_abs) max_abs = abs(stress(2))
            if (abs(stress(3)) > max_abs) max_abs = abs(stress(3))
        end do
        
        if (any(stress /= stress)) then
            write(*,'(A,I5)') '  NaN at cycle ', cycle
            exit
        end if
    end do
    
    call check('Endurance: 500 cycles no NaN', all(stress == stress), np, nf, nt)
    call check('Endurance: max |stress| < 5*|fc|', max_abs < abs(fc)*5.0d0, np, nf, nt)
    write(*,'(A,ES12.4,A,ES12.4)') '    max_abs = ', max_abs, '  5*|fc| = ', abs(fc)*5.0d0
end subroutine

end program test_secant_fix
