!===========================================================================
! Wall-Like Simulation Test
!
! Simulates what OpenSees PlateFromPlaneStress does:
! - Plane stress condition (sig_yy = 0, solve for eps_yy via NR)
! - Cyclic eps_xx loading (±0.05% to ±0.5%)
! - Combined compression + shear (like wall bottom element)
!
! This is the most realistic standalone test for the wall analysis.
!===========================================================================
program test_wall_sim
    use adina_concrete_mod, only: dp
    implicit none

    integer, parameter :: NP = 37, NS = 13
    integer, parameter :: MAX_ITER = 80
    real(dp), parameter :: TOL = 1.0d-10
    real(dp), parameter :: deps = 5.0d-5  ! 0.005% per step

    real(dp) :: props(NP)
    real(dp) :: statev(NS), stress(3), strain0(3)
    real(dp) :: tangent(3,3)
    real(dp) :: eps_xx, eps_yy, sig_xx
    real(dp) :: stress_trial(3), statev_trial(NS)
    real(dp) :: strain1(3), dstrain(3)
    real(dp) :: d_eps_yy, residual, max_step, abs_res_old
    real(dp) :: amp_list(10), amp, target
    real(dp) :: sig_xx_prev, dsig_max
    real(dp) :: peak_pos_sig, peak_neg_sig
    integer  :: iamp, phase, istep, iter, total_steps
    integer  :: max_iters, nr_fail_count, total_iters
    integer  :: nsteps

    call setup_material(props)

    statev = 0.0d0
    stress = 0.0d0
    strain0 = 0.0d0
    eps_xx = 0.0d0
    eps_yy = 0.0d0
    sig_xx_prev = 0.0d0
    dsig_max = 0.0d0
    peak_pos_sig = 0.0d0
    peak_neg_sig = 0.0d0
    total_steps = 0
    max_iters = 0
    nr_fail_count = 0
    total_iters = 0

    ! Amplitude protocol (similar to reference model)
    amp_list = (/0.0005d0, 0.001d0, 0.0015d0, 0.002d0, 0.0025d0, &
                 0.003d0, 0.004d0, 0.005d0, 0.007d0, 0.01d0/)

    open(20, file='test_wall_sim.csv', status='replace')
    write(20,'(A)') 'step,eps_xx,sig_xx,eps_yy,nr_iters'

    write(*,'(A)') '========================================================'
    write(*,'(A)') '  Wall-Like Simulation: Cyclic Plane Stress'
    write(*,'(A)') '  Protocol: ±0.05% to ±1.0%'
    write(*,'(A)') '========================================================'

    do iamp = 1, 10
        amp = amp_list(iamp)

        do phase = 1, 2
            if (phase == 1) then
                target = -amp
            else
                target = amp
            end if

            nsteps = nint(abs(target - eps_xx) / deps)
            if (nsteps < 1) nsteps = 1

            do istep = 1, nsteps
                total_steps = total_steps + 1
                if (target > eps_xx) then
                    eps_xx = eps_xx + deps
                else
                    eps_xx = eps_xx - deps
                end if

                ! NR iteration for plane stress (sig_yy = 0)
                abs_res_old = 1.0d30
                max_step = max(abs(eps_xx) * 0.5d0, 1.0d-5)

                do iter = 1, MAX_ITER
                    strain1(1) = eps_xx
                    strain1(2) = eps_yy
                    strain1(3) = 0.0d0
                    dstrain = strain1 - strain0

                    stress_trial = stress
                    statev_trial = statev

                    call PSUMAT(NS, NP, props, stress_trial, strain0, &
                                strain1, dstrain, statev_trial, tangent)

                    residual = stress_trial(2)
                    if (abs(residual) < TOL) exit
                    if (abs(tangent(2,2)) < 1.0d-30) exit

                    d_eps_yy = -residual / tangent(2,2)
                    if (abs(d_eps_yy) > max_step) d_eps_yy = sign(max_step, d_eps_yy)
                    if (abs(residual) > abs_res_old * 1.1d0 .and. iter > 1) then
                        d_eps_yy = d_eps_yy * 0.5d0
                    end if
                    abs_res_old = abs(residual)
                    eps_yy = eps_yy + d_eps_yy
                    if (abs(d_eps_yy) < 1.0d-15) exit
                end do

                ! Commit converged step
                strain1(1) = eps_xx
                strain1(2) = eps_yy
                strain1(3) = 0.0d0
                dstrain = strain1 - strain0
                call PSUMAT(NS, NP, props, stress, strain0, &
                            strain1, dstrain, statev, tangent)
                strain0 = strain1
                sig_xx = stress(1)

                total_iters = total_iters + iter
                if (iter > max_iters) max_iters = iter
                if (iter >= MAX_ITER) nr_fail_count = nr_fail_count + 1

                dsig_max = max(dsig_max, abs(sig_xx - sig_xx_prev))
                sig_xx_prev = sig_xx
                if (sig_xx < peak_neg_sig) peak_neg_sig = sig_xx
                if (sig_xx > peak_pos_sig) peak_pos_sig = sig_xx

                write(20,'(I6,A,ES14.6,A,ES14.6,A,ES14.6,A,I3)') &
                    total_steps, ',', eps_xx, ',', sig_xx, ',', eps_yy, ',', iter
            end do

            write(*,'(A,I2,A,F7.4,A,F7.4,A,F8.3,A)') &
                '  Amp ', iamp, ': eps=', target*100, &
                '% → sig=', sig_xx, ' MPa (', eps_yy*100, '%)'
        end do
    end do
    close(20)

    write(*,*)
    write(*,'(A,I6)')    '  Total steps: ', total_steps
    write(*,'(A,I3)')     '  Max NR iterations/step: ', max_iters
    write(*,'(A,F6.1)')   '  Avg NR iterations/step: ', dble(total_iters)/dble(total_steps)
    write(*,'(A,I3)')     '  NR failures: ', nr_fail_count
    write(*,'(A,F10.3,A,F10.3,A)') '  Stress peaks: [', peak_neg_sig, ', ', peak_pos_sig, '] MPa'
    write(*,'(A,F8.4,A)') '  Max stress jump: ', dsig_max, ' MPa'

    if (peak_neg_sig < -0.5d0 .and. peak_pos_sig > 0.5d0) then
        write(*,'(A,F6.3)') '  Symmetry ratio: ', abs(peak_pos_sig / peak_neg_sig)
    end if

    if (nr_fail_count == 0) then
        write(*,'(A)') '  RESULT: ALL STEPS CONVERGED'
    else
        write(*,'(A,I3,A)') '  RESULT: ', nr_fail_count, ' STEPS FAILED TO CONVERGE'
    end if

contains

subroutine setup_material(props)
    real(dp), intent(out) :: props(NP)
    props = 0.0d0
    props(1)  = 30000.0d0     ! E0
    props(2)  = 0.2d0         ! nu
    props(3)  = 3.0d0         ! ft
    props(4)  = -30.0d0       ! fc
    props(5)  = -0.002d0      ! eps_c
    props(6)  = -6.0d0        ! fu
    props(7)  = -0.005d0      ! eps_u
    props(8)  = 1.0d-4        ! STIFAC
    props(9)  = 0.5d0         ! SHEFAC
    props(10) = 0.75d0; props(11) = 1.0d0; props(12) = 0.7d0; props(13) = 0.12d0
    props(14) = 0.0d0;  props(15) = 0.25d0; props(16) = 0.5d0
    props(17) = 0.75d0; props(18) = 1.0d0;  props(19) = 1.2d0
    props(20) = 1.0d0;  props(21) = 1.4d0;  props(22) = 1.7d0
    props(23) = 2.2d0;  props(24) = 2.5d0;  props(25) = 2.8d0
    props(26) = 1.3d0;  props(27) = 1.5d0;  props(28) = 2.0d0
    props(29) = 2.3d0;  props(30) = 2.7d0;  props(31) = 3.2d0
    props(32) = 1.25d0; props(33) = 1.45d0; props(34) = 1.95d0
    props(35) = 2.25d0; props(36) = 2.65d0; props(37) = 3.15d0
end subroutine

end program test_wall_sim
