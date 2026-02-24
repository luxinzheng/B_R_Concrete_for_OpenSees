!===============================================================================
!  Pure Shear Test
!
!  Loading: eps_xx = 0, eps_yy = 0, gamma_xy from 0 to 0.5% in 100 steps
!
!  Theoretical response:
!    Phase 1 (elastic):  tau = G * gamma,  G = E/(2(1+nu)) = 12500 MPa
!    Phase 2 (cracked):  cracking at ~45 deg when sigma_1 > ft
!                         reduced shear stiffness: beta * G
!
!  Output: test_shear.csv
!===============================================================================
program test_shear
    use adina_concrete_mod, only: dp
    implicit none
    
    integer, parameter :: NPROPS   = 37
    integer, parameter :: NSTATEVS = 13
    integer, parameter :: NSTEPS   = 100
    
    real(dp) :: props(NPROPS)
    real(dp) :: statev(NSTATEVS), stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: tangent(3,3)
    real(dp) :: gamma, dgamma
    integer  :: istep, iunit
    
    call setup_material(props)
    
    statev  = 0.0d0
    stress  = 0.0d0
    strain0 = 0.0d0
    
    dgamma = 0.005d0 / dble(NSTEPS)   ! 0.5% / 100 steps = 0.005% per step
    
    write(*,*) '============================================'
    write(*,*) ' Pure Shear Test'
    write(*,*) ' gamma_xy: 0 to 0.5% in 100 steps'
    write(*,*) ' G = E/(2(1+nu)) = 12500 MPa'
    write(*,*) ' gamma_cr = ft/G = 0.024%'
    write(*,*) '============================================'
    
    open(newunit=iunit, file='test_shear.csv', status='replace')
    write(iunit, '(A)') 'step, gamma_xy, tau_xy, sig_xx, sig_yy'
    
    ! Initial state
    write(iunit, '(I6,A,4(ES20.12,A))') 0, ',', &
        0.0d0, ',', 0.0d0, ',', 0.0d0, ',', 0.0d0, ''
    
    do istep = 1, NSTEPS
        gamma = dgamma * dble(istep)
        
        strain1(1) = 0.0d0    ! eps_xx = 0
        strain1(2) = 0.0d0    ! eps_yy = 0
        strain1(3) = gamma     ! gamma_xy
        dstrain = strain1 - strain0
        
        call PSUMAT(NSTATEVS, NPROPS, props, stress, strain0, &
                    strain1, dstrain, statev, tangent)
        
        strain0 = strain1
        
        write(iunit, '(I6,A,4(ES20.12,A))') istep, ',', &
            gamma, ',', stress(3), ',', stress(1), ',', stress(2), ''
        
        if (mod(istep, 10) == 0 .or. istep <= 10) then
            write(*,'(A,I4,A,F8.4,A,F10.4,A,F10.4,A,F10.4)') &
                '  Step', istep, ': gamma%=', gamma*100, &
                ' tau=', stress(3), ' sx=', stress(1), ' sy=', stress(2)
        end if
    end do
    
    close(iunit)
    write(*,*) 'Output: test_shear.csv'
    write(*,*) 'Done.'
    
contains

subroutine setup_material(props)
    real(dp), intent(out) :: props(NPROPS)
    props = 0.0d0
    props(1)  = 30000.0d0     ! E0
    props(2)  = 0.2d0         ! nu
    props(3)  = 3.0d0         ! ft
    props(4)  = -30.0d0       ! fc
    props(5)  = -0.002d0      ! eps_c
    props(6)  = -6.0d0        ! fu
    props(7)  = -0.005d0      ! eps_u
    props(8)  = 1.0d-4        ! STIFAC
    props(9)  = 0.5d0         ! SHEFAC (shear retention factor)
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

end program test_shear
