!===============================================================================
! Generate stress-strain data files for plotting (before/after comparison)
! Outputs CSV files for each test scenario
!===============================================================================
program test_secant_plot
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    call gen_uniaxial_cyclic()
    call gen_crack_cycle()
    call gen_ml_scenario()
    
contains

subroutine init_props_ml(props)
    real(dp), intent(out) :: props(37)
    props = 0.0d0
    props(1)  = 2.40254d10; props(2)  = 0.2d0; props(3)  = 2.0d6
    props(4)  = -2.58d7;    props(5)  = -0.003d0; props(6)  = -5.16d6
    props(7)  = -0.021d0;   props(8)  = 0.0001d0; props(9)  = 0.5d0
    props(10) = 0.75d0;     props(11) = 1.0d0; props(12) = 0.7d0
    props(13) = 0.12d0
    props(14:19) = (/ 0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.2d0 /)
    props(20:25) = (/ 1.0d0, 1.4d0, 1.7d0, 2.2d0, 2.5d0, 2.8d0 /)
    props(26:31) = (/ 1.3d0, 1.5d0, 2.0d0, 2.3d0, 2.7d0, 3.2d0 /)
    props(32:37) = (/ 1.25d0, 1.45d0, 1.95d0, 2.25d0, 2.65d0, 3.15d0 /)
end subroutine

subroutine cpsumat(props, stress, s0, s1, sv, tg)
    real(dp), intent(in) :: props(37)
    real(dp), intent(inout) :: stress(3), sv(40)
    real(dp), intent(in) :: s0(3), s1(3)
    real(dp), intent(out) :: tg(3,3)
    real(dp) :: ds(3)
    ds = s1 - s0
    call PSUMAT(40, 37, props, stress, s0, s1, ds, sv, tg)
end subroutine

! Uniaxial cyclic: 20 cycles of compress → unload → tension
subroutine gen_uniaxial_cyclic()
    real(dp) :: props(37), stress(3), statev(40), tg(3,3)
    real(dp) :: s0(3), s1(3), eps
    integer :: cycle, i, step
    
    call init_props_ml(props)
    stress = 0.0d0; statev = 0.0d0; s0 = 0.0d0; step = 0
    
    open(unit=10, file='secant_uniaxial_cyclic.csv', status='replace')
    write(10,'(A)') 'step,eps_xx,sig_xx,sig_yy,tau_xy'
    
    do cycle = 1, 20
        do i = 1, 20
            eps = -1.0d-4 * dble(i)
            s1 = (/ eps, 0.0d0, 0.0d0 /)
            call cpsumat(props, stress, s0, s1, statev, tg)
            s0 = s1; step = step + 1
            write(10,'(I6,",",4(ES15.6,","))') step, eps, stress(1), stress(2), stress(3)
        end do
        do i = 19, -1, -1
            eps = -1.0d-4 * dble(i)
            s1 = (/ eps, 0.0d0, 0.0d0 /)
            call cpsumat(props, stress, s0, s1, statev, tg)
            s0 = s1; step = step + 1
            write(10,'(I6,",",4(ES15.6,","))') step, eps, stress(1), stress(2), stress(3)
        end do
    end do
    close(10)
    write(*,'(A)') 'Generated: secant_uniaxial_cyclic.csv'
end subroutine

! Crack open/close/reopen: 50 cycles
subroutine gen_crack_cycle()
    real(dp) :: props(37), stress(3), statev(40), tg(3,3)
    real(dp) :: s0(3), s1(3), ex, ey
    integer :: cycle, i, step
    
    call init_props_ml(props)
    stress = 0.0d0; statev = 0.0d0; s0 = 0.0d0; step = 0
    
    open(unit=11, file='secant_crack_cycle.csv', status='replace')
    write(11,'(A)') 'step,eps_xx,sig_xx,sig_yy,tau_xy'
    
    do cycle = 1, 50
        do i = 1, 10
            ex = 2.0d-5 * dble(i); ey = -0.2d0 * ex
            s1 = (/ ex, ey, 0.0d0 /)
            call cpsumat(props, stress, s0, s1, statev, tg)
            s0 = s1; step = step + 1
            write(11,'(I6,",",4(ES15.6,","))') step, ex, stress(1), stress(2), stress(3)
        end do
        do i = 10, -10, -1
            ex = 2.0d-5 * dble(i); ey = -0.2d0 * ex
            s1 = (/ ex, ey, 1.0d-5 /)
            call cpsumat(props, stress, s0, s1, statev, tg)
            s0 = s1; step = step + 1
            write(11,'(I6,",",4(ES15.6,","))') step, ex, stress(1), stress(2), stress(3)
        end do
        do i = -10, 0
            ex = 2.0d-5 * dble(i); ey = -0.2d0 * ex
            s1 = (/ ex, ey, 0.0d0 /)
            call cpsumat(props, stress, s0, s1, statev, tg)
            s0 = s1; step = step + 1
            write(11,'(I6,",",4(ES15.6,","))') step, ex, stress(1), stress(2), stress(3)
        end do
    end do
    close(11)
    write(*,'(A)') 'Generated: secant_crack_cycle.csv'
end subroutine

! ML scenario: 11 blocks with increasing amplitude + shear
subroutine gen_ml_scenario()
    real(dp) :: props(37), stress(3), statev(40), tg(3,3)
    real(dp) :: s0(3), s1(3), amp, ex, ey, gxy
    integer :: block, i, nstep, step
    
    call init_props_ml(props)
    stress = 0.0d0; statev = 0.0d0; s0 = 0.0d0; step = 0; nstep = 10
    
    open(unit=12, file='secant_ml_scenario.csv', status='replace')
    write(12,'(A)') 'step,eps_xx,sig_xx,sig_yy,tau_xy'
    
    do block = 1, 11
        amp = 1.0d-3 * dble(block)
        do i = 1, nstep
            ex = amp * dble(i) / dble(nstep)
            ey = -0.15d0 * ex; gxy = 0.3d0 * ex
            s1 = (/ ex, ey, gxy /)
            call cpsumat(props, stress, s0, s1, statev, tg)
            s0 = s1; step = step + 1
            write(12,'(I6,",",4(ES15.6,","))') step, ex, stress(1), stress(2), stress(3)
        end do
        do i = nstep, -nstep, -1
            ex = amp * dble(i) / dble(nstep)
            ey = -0.15d0 * ex; gxy = 0.3d0 * ex
            s1 = (/ ex, ey, gxy /)
            call cpsumat(props, stress, s0, s1, statev, tg)
            s0 = s1; step = step + 1
            write(12,'(I6,",",4(ES15.6,","))') step, ex, stress(1), stress(2), stress(3)
        end do
        do i = -nstep, 0
            ex = amp * dble(i) / dble(nstep)
            ey = -0.15d0 * ex; gxy = 0.3d0 * ex
            s1 = (/ ex, ey, gxy /)
            call cpsumat(props, stress, s0, s1, statev, tg)
            s0 = s1; step = step + 1
            write(12,'(I6,",",4(ES15.6,","))') step, ex, stress(1), stress(2), stress(3)
        end do
    end do
    close(12)
    write(*,'(A)') 'Generated: secant_ml_scenario.csv'
end subroutine

end program test_secant_plot
