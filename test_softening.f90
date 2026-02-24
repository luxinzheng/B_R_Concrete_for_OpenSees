program test_softening
    use adina_concrete_mod, only: dp
    implicit none
    
    integer, parameter :: NPROPS = 37, NSTATEVS = 13, MAX_NR = 100
    real(dp) :: props(NPROPS), statev(NSTATEVS), stress(3), strain0(3)
    real(dp) :: strain1(3), dstrain(3), tangent(3,3)
    real(dp) :: statev_t(NSTATEVS), stress_t(3)
    real(dp) :: eps_xx, eps_yy, deps, residual, d_eps
    real(dp) :: peak_sig, peak_eps
    integer :: i, nsteps, iter, iunit
    
    props = 0.0d0
    props(1)  = 21.4d9;  props(2) = 0.2d0
    props(3)  = 2.07d6;  props(4) = -20.7d6
    props(5)  = -0.002d0; props(6) = -4.14d6; props(7) = -0.006d0
    props(8)  = 1.0d-4;  props(9)  = 0.5d0
    props(10) = 0.75d0;  props(11) = 1.0d0;   props(12) = 0.7d0; props(13) = 0.12d0
    props(14:19) = (/ 0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.2d0 /)
    props(20:25) = (/ 1.0d0, 1.4d0,  1.7d0, 2.2d0,  2.5d0, 2.8d0 /)
    props(26:31) = (/ 1.3d0, 1.5d0,  2.0d0, 2.3d0,  2.7d0, 3.2d0 /)
    props(32:37) = (/ 1.25d0,1.45d0, 1.95d0,2.25d0, 2.65d0,3.15d0 /)
    
    statev = 0.0d0; stress = 0.0d0; strain0 = 0.0d0
    nsteps = 200
    deps = -0.008d0 / dble(nsteps)
    
    iunit = 20
    open(unit=iunit, file='test_softening.csv', status='replace')
    write(iunit,'(A)') 'eps_xx,sig_xx,eps_yy,E_sec_approx'
    write(iunit,'(4ES18.8)') 0.0d0, 0.0d0, 0.0d0, props(1)
    
    eps_yy = 0.0d0
    peak_sig = 0.0d0; peak_eps = 0.0d0
    
    do i = 1, nsteps
        eps_xx = deps * dble(i)
        
        do iter = 1, MAX_NR
            strain1(1) = eps_xx; strain1(2) = eps_yy; strain1(3) = 0.0d0
            dstrain = strain1 - strain0
            stress_t = stress; statev_t = statev
            call PSUMAT(NSTATEVS, NPROPS, props, stress_t, strain0, &
                        strain1, dstrain, statev_t, tangent)
            residual = stress_t(2)
            if (abs(residual) < 1.0d-6) exit
            if (abs(tangent(2,2)) < 1.0d-20) exit
            d_eps = -residual / tangent(2,2)
            if (abs(d_eps) > 5.0d-4) d_eps = sign(5.0d-4, d_eps)
            eps_yy = eps_yy + d_eps
        end do
        
        strain1(1) = eps_xx; strain1(2) = eps_yy; strain1(3) = 0.0d0
        dstrain = strain1 - strain0
        call PSUMAT(NSTATEVS, NPROPS, props, stress, strain0, &
                    strain1, dstrain, statev, tangent)
        strain0 = strain1
        
        if (stress(1) < peak_sig) then
            peak_sig = stress(1); peak_eps = eps_xx
        end if
        
        write(iunit,'(4ES18.8)') eps_xx, stress(1), eps_yy, tangent(1,1)
        
        if (mod(i, 20) == 0) then
            write(*,'(A,I4,A,ES11.3,A,ES11.3,A,ES11.3,A,I3)') &
                '  step', i, '  eps=', eps_xx, '  sig=', stress(1), &
                '  eps_yy=', eps_yy, '  NR=', iter
        end if
    end do
    
    close(iunit)
    write(*,*)
    write(*,'(A,ES12.4,A,ES12.4)') '  Peak stress = ', peak_sig, ' at eps = ', peak_eps
    write(*,'(A,ES12.4)') '  Final stress = ', stress(1)
    write(*,'(A,ES12.4)') '  fc (input)   = ', props(4)
    write(*,'(A,ES12.4)') '  fu (input)   = ', props(6)
    write(*,*) 'Output: test_softening.csv'
    
end program test_softening
