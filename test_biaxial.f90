!===============================================================================
!  Biaxial Monotonic Loading Test
!
!  5 loading cases with strain ratios (eps_xx : eps_yy):
!    Case 1: -1 : -1   (equal biaxial compression)
!    Case 2:  1 :  1   (equal biaxial tension)
!    Case 3: -1 : -0.5 (unequal biaxial compression)
!    Case 4: -1 :  0.1 (compression + tension)
!    Case 5: -1 :  0.2 (compression + tension, sigma_yy ~ 0)
!
!  At least 40 steps, max |strain| = 1%
!  Output: test_biaxial_case{1..5}.csv
!===============================================================================
program test_biaxial
    use adina_concrete_mod, only: dp
    implicit none
    
    integer, parameter :: NPROPS   = 37
    integer, parameter :: NSTATEVS = 13
    integer, parameter :: NSTEPS   = 40
    integer, parameter :: NCASES   = 5
    
    real(dp) :: props(NPROPS)
    real(dp) :: statev(NSTATEVS), stress(3), strain0(3), strain1(3), dstrain(3)
    real(dp) :: tangent(3,3)
    real(dp) :: ratio_x(NCASES), ratio_y(NCASES)
    real(dp) :: max_strain, deps_scale, eps_xx, eps_yy
    integer  :: icase, istep, iunit
    character(len=80) :: filename
    character(len=80) :: case_names(NCASES)
    
    call setup_material(props)
    
    ! Loading ratios
    ratio_x(1) = -1.0d0;  ratio_y(1) = -1.0d0
    ratio_x(2) =  1.0d0;  ratio_y(2) =  1.0d0
    ratio_x(3) = -1.0d0;  ratio_y(3) = -0.5d0
    ratio_x(4) = -1.0d0;  ratio_y(4) =  0.1d0
    ratio_x(5) = -1.0d0;  ratio_y(5) =  0.2d0
    
    case_names(1) = '-1:-1   (biaxial compression)'
    case_names(2) = ' 1:1    (biaxial tension)'
    case_names(3) = '-1:-0.5 (unequal biaxial comp)'
    case_names(4) = '-1:0.1  (comp + tension)'
    case_names(5) = '-1:0.2  (comp + tension, sig_yy~0)'
    
    max_strain = 0.01d0   ! 1%
    
    write(*,*) '============================================'
    write(*,*) ' Biaxial Monotonic Loading Test'
    write(*,*) ' 5 cases, 40 steps each, max |strain|=1%'
    write(*,*) '============================================'
    
    do icase = 1, NCASES
        ! Strain increment: the larger component reaches 1% after NSTEPS steps
        deps_scale = max_strain / dble(NSTEPS) / max(abs(ratio_x(icase)), abs(ratio_y(icase)))
        
        ! Reset state
        statev  = 0.0d0
        stress  = 0.0d0
        strain0 = 0.0d0
        
        ! Open output CSV
        write(filename, '(A,I1,A)') 'test_biaxial_case', icase, '.csv'
        open(newunit=iunit, file=trim(filename), status='replace')
        write(iunit, '(A)') 'step, eps_xx, eps_yy, sig_xx, sig_yy, tau_xy'
        
        ! Initial state
        write(iunit, '(I6,A,5(ES20.12,A))') 0, ',', &
            0.0d0, ',', 0.0d0, ',', 0.0d0, ',', 0.0d0, ',', 0.0d0, ''
        
        do istep = 1, NSTEPS
            eps_xx = ratio_x(icase) * deps_scale * dble(istep)
            eps_yy = ratio_y(icase) * deps_scale * dble(istep)
            
            strain1(1) = eps_xx
            strain1(2) = eps_yy
            strain1(3) = 0.0d0
            dstrain = strain1 - strain0
            
            call PSUMAT(NSTATEVS, NPROPS, props, stress, strain0, &
                        strain1, dstrain, statev, tangent)
            
            strain0 = strain1
            
            write(iunit, '(I6,A,5(ES20.12,A))') istep, ',', &
                eps_xx, ',', eps_yy, ',', stress(1), ',', stress(2), ',', stress(3), ''
        end do
        
        close(iunit)
        write(*,'(A,I1,A,A,A,A)') '  Case ', icase, ': ', &
            trim(case_names(icase)), ' -> ', trim(filename)
        write(*,'(A,2ES12.4)') '    Final stress (xx,yy): ', stress(1), stress(2)
    end do
    
    write(*,*) 'Done.'
    
contains

subroutine setup_material(props)
    real(dp), intent(out) :: props(NPROPS)
    props = 0.0d0
    props(1)  = 30000.0d0     ! E0
    props(2)  = 0.2d0         ! nu
    props(3)  = 3.0d0         ! ft  (> 0)
    props(4)  = -30.0d0       ! fc  (< 0)
    props(5)  = -0.002d0      ! eps_c
    props(6)  = -6.0d0        ! fu
    props(7)  = -0.005d0      ! eps_u
    props(8)  = 1.0d-4        ! STIFAC
    props(9)  = 0.5d0         ! SHEFAC
    props(10) = 0.75d0; props(11) = 1.0d0; props(12) = 0.7d0; props(13) = 0.12d0
    ! Biaxial envelope spline points
    props(14) = 0.0d0;  props(15) = 0.25d0; props(16) = 0.5d0
    props(17) = 0.75d0; props(18) = 1.0d0;  props(19) = 1.2d0
    props(20) = 1.0d0;  props(21) = 1.4d0;  props(22) = 1.7d0
    props(23) = 2.2d0;  props(24) = 2.5d0;  props(25) = 2.8d0
    props(26) = 1.3d0;  props(27) = 1.5d0;  props(28) = 2.0d0
    props(29) = 2.3d0;  props(30) = 2.7d0;  props(31) = 3.2d0
    props(32) = 1.25d0; props(33) = 1.45d0; props(34) = 1.95d0
    props(35) = 2.25d0; props(36) = 2.65d0; props(37) = 3.15d0
end subroutine

end program test_biaxial
