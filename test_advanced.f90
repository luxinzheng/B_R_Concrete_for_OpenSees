!===============================================================================
!  Advanced tests: NR convergence, cracked state tangent, edge cases
!===============================================================================
program test_advanced
    use adina_concrete_mod
    implicit none
    integer, parameter :: NP=37, NS=13
    real(dp) :: props(NP)
    integer :: npass, nfail
    
    npass = 0; nfail = 0
    call setup_props(props)
    
    write(*,'(A)') '========================================================'
    write(*,'(A)') '  Advanced Tests: NR Convergence & Edge Cases'
    write(*,'(A)') '========================================================'
    write(*,*)
    
    call test_nr_convergence_compression(props, npass, nfail)
    call test_nr_convergence_tension(props, npass, nfail)
    call test_nr_convergence_cracked_shear(props, npass, nfail)
    call test_cracked_cyclic(props, npass, nfail)
    call test_tangent_positive_definite(props, npass, nfail)
    call test_stress_strain_csv(props)
    
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
! Newton-Raphson convergence test: compression
!   Apply target strain, iterate using tangent to find equilibrium
! ===================================================================
subroutine test_nr_convergence_compression(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    real(dp) :: sig_target(3), residual(3), dstrain_corr(3)
    real(dp) :: rnorm, rnorm0, SBAR
    integer :: iter, max_iter
    
    write(*,'(A)') 'TEST A1: NR convergence — compression (target σ_xx=-15 MPa)'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    sig_target = 0.0d0
    sig_target(1) = -15.0d6
    
    strain1 = 0.0d0
    strain1(1) = sig_target(1) * (1.0d0 - props(2)**2) / props(1)
    
    max_iter = 20
    do iter = 1, max_iter
        dstrain = strain1 - strain0
        statev = 0.0d0
        stress = 0.0d0
        call PSUMAT(NS, NP, props, stress, strain0, strain1, dstrain, statev, tangent)
        
        residual = stress - sig_target
        rnorm = sqrt(sum(residual**2))
        if (iter == 1) rnorm0 = rnorm
        
        if (rnorm < 1.0d-3) exit
        
        dstrain_corr(1) =  (tangent(2,2)*tangent(3,3)-tangent(2,3)*tangent(3,2))*residual(1) &
                         - (tangent(1,2)*tangent(3,3)-tangent(1,3)*tangent(3,2))*residual(2) &
                         + (tangent(1,2)*tangent(2,3)-tangent(1,3)*tangent(2,2))*residual(3)
        dstrain_corr(2) = -(tangent(2,1)*tangent(3,3)-tangent(2,3)*tangent(3,1))*residual(1) &
                         + (tangent(1,1)*tangent(3,3)-tangent(1,3)*tangent(3,1))*residual(2) &
                         - (tangent(1,1)*tangent(2,3)-tangent(1,3)*tangent(2,1))*residual(3)
        dstrain_corr(3) =  (tangent(2,1)*tangent(3,2)-tangent(2,2)*tangent(3,1))*residual(1) &
                         - (tangent(1,1)*tangent(3,2)-tangent(1,2)*tangent(3,1))*residual(2) &
                         + (tangent(1,1)*tangent(2,2)-tangent(1,2)*tangent(2,1))*residual(3)
        
        SBAR = tangent(1,1)*(tangent(2,2)*tangent(3,3)-tangent(2,3)*tangent(3,2)) &
              -tangent(1,2)*(tangent(2,1)*tangent(3,3)-tangent(2,3)*tangent(3,1)) &
              +tangent(1,3)*(tangent(2,1)*tangent(3,2)-tangent(2,2)*tangent(3,1))
        if (abs(SBAR) > 1.0d-30) dstrain_corr = dstrain_corr / SBAR
        
        strain1 = strain1 - dstrain_corr
    end do
    
    write(*,'(A,I3,A,ES10.3)') '  Converged in ', iter, ' iters, |R|=', rnorm
    call check('NR-Comp: converged', iter <= max_iter, npass, nfail)
    call check('NR-Comp: stress match', abs(stress(1) - sig_target(1)) / abs(sig_target(1)) < 0.01d0, npass, nfail)
    write(*,*)
end subroutine

! ===================================================================
! NR convergence test: tension (cracking)
! ===================================================================
subroutine test_nr_convergence_tension(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    real(dp) :: sig_target(3), residual(3), dstrain_corr(3)
    real(dp) :: rnorm, SBAR
    integer :: iter, max_iter
    
    write(*,'(A)') 'TEST A2: NR convergence — tension (target σ_xx=1.0 MPa post-crack)'
    
    ! First drive to cracking
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    strain1 = 0.0d0; strain1(1) = 2.0d-4
    dstrain = strain1 - strain0
    call PSUMAT(NS, NP, props, stress, strain0, strain1, dstrain, statev, tangent)
    strain0 = strain1
    
    write(*,'(A,ES12.4,A,F8.1)') '  Post-crack stress=', stress(1), ', ANGLE=', statev(2)
    
    ! Now try to find strain for σ_xx = 1.0 MPa (in softening zone)
    sig_target = 0.0d0
    sig_target(1) = 1.0d6
    
    strain1(1) = 2.5d-4
    
    max_iter = 20
    do iter = 1, max_iter
        dstrain = strain1 - strain0
        stress = 0.0d0  ! reset for fresh call
        call PSUMAT(NS, NP, props, stress, strain0, strain1, dstrain, statev, tangent)
        
        residual = stress - sig_target
        rnorm = sqrt(sum(residual**2))
        
        if (rnorm < 1.0d2) exit
        
        ! Solve dε = C^{-1} * R (using Cramer's rule)
        SBAR = tangent(1,1)*(tangent(2,2)*tangent(3,3)-tangent(2,3)*tangent(3,2)) &
              -tangent(1,2)*(tangent(2,1)*tangent(3,3)-tangent(2,3)*tangent(3,1)) &
              +tangent(1,3)*(tangent(2,1)*tangent(3,2)-tangent(2,2)*tangent(3,1))
        if (abs(SBAR) < 1.0d-30) exit
        
        dstrain_corr(1) =  (tangent(2,2)*tangent(3,3)-tangent(2,3)*tangent(3,2))*residual(1) &
                         - (tangent(1,2)*tangent(3,3)-tangent(1,3)*tangent(3,2))*residual(2) &
                         + (tangent(1,2)*tangent(2,3)-tangent(1,3)*tangent(2,2))*residual(3)
        dstrain_corr(2) = -(tangent(2,1)*tangent(3,3)-tangent(2,3)*tangent(3,1))*residual(1) &
                         + (tangent(1,1)*tangent(3,3)-tangent(1,3)*tangent(3,1))*residual(2) &
                         - (tangent(1,1)*tangent(2,3)-tangent(1,3)*tangent(2,1))*residual(3)
        dstrain_corr(3) =  (tangent(2,1)*tangent(3,2)-tangent(2,2)*tangent(3,1))*residual(1) &
                         - (tangent(1,1)*tangent(3,2)-tangent(1,2)*tangent(3,1))*residual(2) &
                         + (tangent(1,1)*tangent(2,2)-tangent(1,2)*tangent(2,1))*residual(3)
        dstrain_corr = dstrain_corr / SBAR
        
        strain1 = strain1 - 0.5d0 * dstrain_corr  ! damped update
    end do
    
    write(*,'(A,I3,A,ES10.3)') '  Iterations: ', iter, ', |R|=', rnorm
    call check('NR-Tens: converged or bounded', rnorm < 1.0d5, npass, nfail)
    write(*,*)
end subroutine

! ===================================================================
! NR convergence: cracked element under shear
! ===================================================================
subroutine test_nr_convergence_cracked_shear(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    integer :: i
    logical :: has_nan
    
    write(*,'(A)') 'TEST A3: Cracked element under shear loading'
    
    ! First crack the element in tension
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    strain1 = 0.0d0; strain1(1) = 5.0d-4
    dstrain = strain1 - strain0
    call PSUMAT(NS, NP, props, stress, strain0, strain1, dstrain, statev, tangent)
    strain0 = strain1
    write(*,'(A,F8.1)') '  After cracking: ANGLE=', statev(2)
    
    ! Now apply shear to the cracked element
    has_nan = .false.
    do i = 1, 20
        strain1(1) = 5.0d-4
        strain1(2) = 0.0d0
        strain1(3) = 1.0d-5 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS, NP, props, stress, strain0, strain1, dstrain, statev, tangent)
        strain0 = strain1
        if (any(stress /= stress) .or. any(tangent /= tangent)) has_nan = .true.
    end do
    
    call check('CrackShear: no NaN', .not. has_nan, npass, nfail)
    call check('CrackShear: positive tangent diag', &
        tangent(1,1) > 0 .and. tangent(2,2) > 0 .and. tangent(3,3) > 0, npass, nfail)
    
    write(*,'(A,3ES12.4)') '  Final stress: ', stress
    write(*,'(A,3ES12.4)') '  Tangent diag: ', tangent(1,1), tangent(2,2), tangent(3,3)
    write(*,*)
end subroutine

! ===================================================================
! Cracked cyclic: crack → close → reopen
! ===================================================================
subroutine test_cracked_cyclic(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    integer :: i
    logical :: has_nan
    
    write(*,'(A)') 'TEST A4: Cracked cyclic (crack → close → reopen)'
    
    ! Crack in tension
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    strain1 = 0.0d0; strain1(1) = 3.0d-4
    dstrain = strain1; call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
    strain0 = strain1
    write(*,'(A,ES12.4,A,F8.1)') '  After crack: sig=', stress(1), ', ANGLE=', statev(2)
    
    ! Close crack by compressing
    has_nan = .false.
    do i = 1, 10
        strain1 = 0.0d0
        strain1(1) = 3.0d-4 - 1.0d-4 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (any(stress /= stress) .or. any(tangent /= tangent)) has_nan = .true.
    end do
    write(*,'(A,ES12.4,A,F8.1)') '  After closing: sig=', stress(1), ', ANGLE=', statev(2)
    
    ! Reopen crack
    do i = 1, 10
        strain1 = 0.0d0
        strain1(1) = -7.0d-4 + 2.0d-4 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        if (any(stress /= stress) .or. any(tangent /= tangent)) has_nan = .true.
    end do
    write(*,'(A,ES12.4,A,F8.1)') '  After reopen: sig=', stress(1), ', ANGLE=', statev(2)
    
    call check('CrackCycl: no NaN', .not. has_nan, npass, nfail)
    call check('CrackCycl: reopened stress ≈ 0', abs(stress(1)) < 2.0d6, npass, nfail)
    write(*,*)
end subroutine

! ===================================================================
! Check tangent positive definiteness across loading states
! ===================================================================
subroutine test_tangent_positive_definite(props, npass, nfail)
    real(dp), intent(in) :: props(NP)
    integer, intent(inout) :: npass, nfail
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    real(dp) :: det_1, det_2, det_3
    integer :: i, n_neg
    
    write(*,'(A)') 'TEST A5: Tangent positive definiteness (compression path)'
    
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    n_neg = 0
    
    do i = 1, 40
        strain1 = 0.0d0; strain1(1) = -0.0001d0 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        
        det_1 = tangent(1,1)
        det_2 = tangent(1,1)*tangent(2,2) - tangent(1,2)*tangent(2,1)
        det_3 = tangent(1,1)*(tangent(2,2)*tangent(3,3)-tangent(2,3)*tangent(3,2)) &
               -tangent(1,2)*(tangent(2,1)*tangent(3,3)-tangent(2,3)*tangent(3,1)) &
               +tangent(1,3)*(tangent(2,1)*tangent(3,2)-tangent(2,2)*tangent(3,1))
        
        if (det_1 < 0 .or. det_2 < 0 .or. det_3 < 0) then
            n_neg = n_neg + 1
            if (n_neg <= 3) write(*,'(A,I3,A,3ES12.4)') '  Non-PD at step ', i, ': minors=', det_1, det_2, det_3
        end if
    end do
    
    call check('PD-Comp: tangent PD for compression path', n_neg == 0, npass, nfail)
    write(*,'(A,I3,A,I3)') '  Non-PD count: ', n_neg, ' / 40'
    
    ! Tension path (softening tangent may be non-PD, that's expected)
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    n_neg = 0
    do i = 1, 30
        strain1 = 0.0d0; strain1(1) = 1.0d-5 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        
        det_1 = tangent(1,1)
        if (det_1 < 0) n_neg = n_neg + 1
    end do
    write(*,'(A,I3,A,I3)') '  Tension: neg C11 count: ', n_neg, ' / 30'
    write(*,*)
end subroutine

! ===================================================================
! Write stress-strain CSV for external plotting
! ===================================================================
subroutine test_stress_strain_csv(props)
    real(dp), intent(in) :: props(NP)
    real(dp) :: stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: statev(NS), tangent(3,3)
    real(dp) :: eps_val
    integer :: i, unit_id
    
    write(*,'(A)') 'Writing stress-strain curves to CSV files...'
    
    ! Uniaxial compression
    unit_id = 20
    open(unit=unit_id, file='curve_compression.csv', status='replace')
    write(unit_id,'(A)') 'eps_xx,sig_xx,sig_yy,C11'
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    do i = 1, 60
        strain1 = 0.0d0; strain1(1) = -0.0001d0 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        write(unit_id,'(ES15.7,A,ES15.7,A,ES15.7,A,ES15.7)') &
            strain1(1), ',', stress(1), ',', stress(2), ',', tangent(1,1)
    end do
    close(unit_id)
    
    ! Uniaxial tension
    unit_id = 21
    open(unit=unit_id, file='curve_tension.csv', status='replace')
    write(unit_id,'(A)') 'eps_xx,sig_xx,sig_yy,C11'
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    do i = 1, 150
        strain1 = 0.0d0; strain1(1) = 1.0d-5 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        write(unit_id,'(ES15.7,A,ES15.7,A,ES15.7,A,ES15.7)') &
            strain1(1), ',', stress(1), ',', stress(2), ',', tangent(1,1)
    end do
    close(unit_id)
    
    ! Cyclic compression-tension
    unit_id = 22
    open(unit=unit_id, file='curve_cyclic.csv', status='replace')
    write(unit_id,'(A)') 'eps_xx,sig_xx,sig_yy,C11'
    stress = 0.0d0; strain0 = 0.0d0; statev = 0.0d0
    ! Compress to -0.003
    do i = 1, 30
        strain1 = 0.0d0; strain1(1) = -0.0001d0 * i
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        write(unit_id,'(ES15.7,A,ES15.7,A,ES15.7,A,ES15.7)') &
            strain1(1), ',', stress(1), ',', stress(2), ',', tangent(1,1)
    end do
    ! Unload to +0.0005
    do i = 1, 35
        eps_val = -0.003d0 + 0.0001d0 * i
        strain1 = 0.0d0; strain1(1) = eps_val
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        write(unit_id,'(ES15.7,A,ES15.7,A,ES15.7,A,ES15.7)') &
            strain1(1), ',', stress(1), ',', stress(2), ',', tangent(1,1)
    end do
    ! Reload to -0.004
    do i = 1, 45
        eps_val = 0.0005d0 - 0.0001d0 * i
        strain1 = 0.0d0; strain1(1) = eps_val
        dstrain = strain1 - strain0
        call PSUMAT(NS,NP,props,stress,strain0,strain1,dstrain,statev,tangent)
        strain0 = strain1
        write(unit_id,'(ES15.7,A,ES15.7,A,ES15.7,A,ES15.7)') &
            strain1(1), ',', stress(1), ',', stress(2), ',', tangent(1,1)
    end do
    close(unit_id)
    
    write(*,'(A)') '  Done: curve_compression.csv, curve_tension.csv, curve_cyclic.csv'
    write(*,*)
end subroutine

end program test_advanced
