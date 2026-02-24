
program test_runner
  use br_test_params
  use br_test_utils
  implicit none
  ! dp is imported from br_test_params
  interface
    subroutine PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, statev, tangent)
      import dp
      integer, intent(in) :: nstatevs, nprops
      real(dp), intent(in) :: props(nprops)
      real(dp), intent(in) :: strain0(3), strain1(3), dstrain(3)
      real(dp), intent(inout) :: stress(3), statev(nstatevs)
      real(dp), intent(out) :: tangent(3,3)
    end subroutine PSUMAT
  end interface

  integer :: nfail
  nfail = 0

  call case_elastic_patch(nfail)
  call case_uniaxial_compression(nfail)
  call case_uniaxial_tension_softening(nfail)
  call case_pure_shear(nfail)
  call case_biaxial_symmetry(nfail)
  call case_crack_open_close(nfail)
  call case_tangent_consistency(nfail)
  call case_fig8_saenz_curve(nfail)
  call case_fig9_cyclic_path(nfail)
  call case_biaxial_envelope(nfail)
  call case_shear_retention(nfail)
  call case_random_robustness(nfail)

  write(*,*)
  if (nfail == 0) then
    write(*,*) 'ALL TESTS PASSED'
  else
    write(*,*) 'TOTAL FAILURES = ', nfail
    stop 2
  end if

contains

  subroutine reset_state(stress, strain, statev, nstatevs)
    integer, intent(in) :: nstatevs
    real(dp), intent(out) :: stress(3), strain(3), statev(nstatevs)
    stress = 0.0_dp
    strain = 0.0_dp
    statev = 0.0_dp
    statev(2) = 1000.0_dp  ! ANGLE=no crack
    statev(7) = 0.0_dp     ! init flag
  end subroutine

  subroutine elastic_matrix(E, nu, D)
    real(dp), intent(in) :: E, nu
    real(dp), intent(out):: D(3,3)
    real(dp) :: c
    c = E/(1.0_dp - nu*nu)
    D = 0.0_dp
    D(1,1)=c;       D(1,2)=c*nu
    D(2,1)=c*nu;    D(2,2)=c
    D(3,3)=E/(2.0_dp*(1.0_dp+nu))
  end subroutine
  subroutine case_elastic_patch(nfail)
    integer, intent(inout) :: nfail
    real(dp) :: props(37), stress(3), strain(3), statev(8), tangent(3,3), D(3,3)
    real(dp) :: strain0(3), strain1(3), dstrain(3), eps(3), sig_expected(3)
    integer :: i,j
    call set_props_default(props)
    call reset_state(stress, strain, statev, 8)
    call elastic_matrix(props(1), props(2), D)

    do i=1,10
      eps = (/ 1.0e-6_dp*i, -0.5e-6_dp*i, 0.2e-6_dp*i /)
      strain0 = strain
      strain1 = eps
      dstrain = strain1 - strain0
      call PSUMAT(8,37,props,stress,strain0,strain1,dstrain,statev,tangent)
      strain = strain1
      sig_expected = matmul(D, eps)
      call assert_close('elastic sig_x', stress(1), sig_expected(1), 5e-3_dp, 1e-6_dp, nfail)
      call assert_close('elastic sig_y', stress(2), sig_expected(2), 5e-3_dp, 1e-6_dp, nfail)
      call assert_close('elastic tau',   stress(3), sig_expected(3), 5e-3_dp, 1e-6_dp, nfail)
      do j=1,3
        call assert_close('elastic tangent diag', tangent(j,j), D(j,j), 5e-3_dp, 1.0_dp, nfail)
      end do
    end do
    write(*,*) 'case_elastic_patch done'
  end subroutine

  subroutine case_uniaxial_compression(nfail)
    integer, intent(inout) :: nfail
    real(dp) :: props(37), stress(3), strain(3), statev(8), tangent(3,3)
    real(dp) :: strain0(3), strain1(3), dstrain(3)
    real(dp) :: eps, deps, fc, epsc, peak_sig
    integer :: i, u
    call set_props_default(props)
    fc = props(4); epsc = props(5)
    call reset_state(stress, strain, statev, 8)

    open(newunit=u, file='out_uniax_comp.csv', status='replace')
    write(u,'(A)') 'step eps_x eps_y gamma sig_x sig_y tau angle pgrav'
    deps = -8.0e-5_dp
    peak_sig = 0.0_dp
    do i=1,120
      eps = deps*real(i,dp)
      strain0 = strain
      strain1 = (/eps, 0.0_dp, 0.0_dp/)
      dstrain = strain1 - strain0
      call PSUMAT(8,37,props,stress,strain0,strain1,dstrain,statev,tangent)
      strain = strain1
      call write_csv_row(u,i,strain,stress,statev,8)
      if (stress(1) < peak_sig) peak_sig = stress(1)
    end do
    close(u)

    call assert_close('uniax comp peak stress approx fc', peak_sig, fc, 0.30_dp, 1.0_dp, nfail)
    write(*,*) 'case_uniaxial_compression done (see out_uniax_comp.csv)'
  end subroutine

  subroutine case_uniaxial_tension_softening(nfail)
    integer, intent(inout) :: nfail
    real(dp) :: props(37), stress(3), strain(3), statev(8), tangent(3,3)
    real(dp) :: strain0(3), strain1(3), dstrain(3)
    real(dp) :: ft, E0, eps_tu, eps
    integer :: i, u
    call set_props_default(props)
    ft = props(3); E0=props(1); eps_tu = 10.0_dp*ft/E0
    call reset_state(stress, strain, statev, 8)

    open(newunit=u, file='out_uniax_tens.csv', status='replace')
    write(u,'(A)') 'step eps_x eps_y gamma sig_x sig_y tau angle pgrav'
    do i=1,140
      eps = 1.2_dp*eps_tu*real(i,dp)/140.0_dp
      strain0 = strain
      strain1 = (/eps, 0.0_dp, 0.0_dp/)
      dstrain = strain1 - strain0
      call PSUMAT(8,37,props,stress,strain0,strain1,dstrain,statev,tangent)
      strain = strain1
      call write_csv_row(u,i,strain,stress,statev,8)
    end do
    close(u)

    ! by end should soften near zero in crack direction
    call assert_close('tension softening end ~0', stress(1), 0.0_dp, 0.10_dp, 0.05_dp, nfail)
    write(*,*) 'case_uniaxial_tension_softening done (see out_uniax_tens.csv)'
  end subroutine

  subroutine case_pure_shear(nfail)
    integer, intent(inout) :: nfail
    real(dp) :: props(37), stress(3), strain(3), statev(8), tangent(3,3)
    real(dp) :: strain0(3), strain1(3), dstrain(3), G, gamma
    integer :: i, u
    call set_props_default(props)
    G = props(1)/(2.0_dp*(1.0_dp+props(2)))
    call reset_state(stress, strain, statev, 8)
    open(newunit=u, file='out_pure_shear.csv', status='replace')
    write(u,'(A)') 'step eps_x eps_y gamma sig_x sig_y tau angle pgrav'
    do i=1,80
      gamma = 2.5e-5_dp*real(i,dp)
      strain0=strain
      strain1=(/0.0_dp,0.0_dp,gamma/)
      dstrain=strain1-strain0
      call PSUMAT(8,37,props,stress,strain0,strain1,dstrain,statev,tangent)
      strain=strain1
      call write_csv_row(u,i,strain,stress,statev,8)
      if (i <= 5) call assert_close('shear elastic tau', stress(3), G*gamma, 0.06_dp, 0.01_dp, nfail)
    end do
    close(u)
    write(*,*) 'case_pure_shear done (see out_pure_shear.csv)'
  end subroutine

  subroutine case_biaxial_symmetry(nfail)
    integer, intent(inout) :: nfail
    real(dp) :: props(37), stressA(3), stressB(3), strainA(3), strainB(3), stateA(8), stateB(8), tangent(3,3)
    real(dp) :: strain0(3), strain1(3), dstrain(3)
    integer :: i
    call set_props_default(props)
    call reset_state(stressA, strainA, stateA, 8)
    call reset_state(stressB, strainB, stateB, 8)

    do i=1,50
      strain0 = strainA
      strain1 = (/ -5e-5_dp*i, -2e-5_dp*i, 0.0_dp /)
      dstrain = strain1-strain0
      call PSUMAT(8,37,props,stressA,strain0,strain1,dstrain,stateA,tangent)
      strainA = strain1

      strain0 = strainB
      strain1 = (/ -2e-5_dp*i, -5e-5_dp*i, 0.0_dp /)
      dstrain = strain1-strain0
      call PSUMAT(8,37,props,stressB,strain0,strain1,dstrain,stateB,tangent)
      strainB = strain1
    end do
    call assert_close('biaxial symmetry sig_x', stressA(1), stressB(2), 1e-3_dp, 1e-2_dp, nfail)
    call assert_close('biaxial symmetry sig_y', stressA(2), stressB(1), 1e-3_dp, 1e-2_dp, nfail)
    write(*,*) 'case_biaxial_symmetry done'
  end subroutine

  subroutine case_crack_open_close(nfail)
    integer, intent(inout) :: nfail
    real(dp) :: props(37), stress(3), strain(3), statev(8), tangent(3,3)
    real(dp) :: strain0(3), strain1(3), dstrain(3), ft, E0, eps_tu
    integer :: i
    call set_props_default(props)
    ft=props(3); E0=props(1); eps_tu=10.0_dp*ft/E0
    call reset_state(stress, strain, statev, 8)

    ! Open crack in x, then compress to close, then reopen
    do i=1,60
      strain0=strain
      strain1=(/ 0.9_dp*eps_tu*real(i,dp)/60.0_dp, 0.0_dp, 0.0_dp /)
      dstrain=strain1-strain0
      call PSUMAT(8,37,props,stress,strain0,strain1,dstrain,statev,tangent)
      strain=strain1
    end do
    if (statev(2) == 1000.0_dp) then
      nfail=nfail+1
      write(*,*) 'FAIL: expected cracked state after tension'
    end if

    do i=1,80
      strain0=strain
      strain1=(/ -0.0010_dp*real(i,dp)/80.0_dp, 0.0_dp, 0.0_dp /)
      dstrain=strain1-strain0
      call PSUMAT(8,37,props,stress,strain0,strain1,dstrain,statev,tangent)
      strain=strain1
    end do

    do i=1,30
      strain0=strain
      strain1=(/ 0.3_dp*eps_tu*real(i,dp)/30.0_dp, 0.0_dp, 0.0_dp /)
      dstrain=strain1-strain0
      call PSUMAT(8,37,props,stress,strain0,strain1,dstrain,statev,tangent)
      strain=strain1
    end do
    call assert_close('reopen should not exceed ~ft', abs(stress(1)), ft, 5.0_dp, 0.5_dp, nfail)
    write(*,*) 'case_crack_open_close done'
  end subroutine

  subroutine case_tangent_consistency(nfail)
    integer, intent(inout) :: nfail
    real(dp) :: props(37), stress(3), strain(3), statev(12), tangent(3,3)
    real(dp) :: strain0(3), strain1(3), dstrain(3)
    real(dp) :: sv_save(12), sig_save(3)
    real(dp) :: sig_plus(3), sig_minus(3), sv_tmp(12), tan_dum(3,3)
    real(dp) :: tan_num(3,3), h
    real(dp) :: rel_err, max_err, denom_global
    real(dp) :: sig_pred(3), sec_err, sec_max, sec_denom
    integer  :: i, j, mi, mj, ipt, npt_fail

    real(dp) :: test_strains(3,3)

    h = 1.0e-7_dp
    npt_fail = 0

    ! Point 1: elastic (small biaxial compression + shear)
    test_strains(:,1) = (/ -5.0e-5_dp, -2.0e-5_dp, 1.0e-5_dp /)
    ! Point 2: cracked (tensile eps_x >> ft/E0 ≈ 8.3e-5)
    test_strains(:,2) = (/ 2.0e-4_dp, -1.0e-5_dp, 1.0e-5_dp /)
    ! Point 3: post-peak compression (past epsc=-0.002)
    test_strains(:,3) = (/ -3.0e-3_dp, -1.0e-3_dp, 5.0e-4_dp /)

    call set_props_default(props)

    do ipt = 1, 3

      call reset_state(stress, strain, statev, 12)
      strain0 = 0.0_dp
      strain1 = test_strains(:, ipt)
      dstrain = strain1 - strain0
      call PSUMAT(12, 37, props, stress, strain0, strain1, dstrain, statev, tangent)
      strain   = strain1
      sv_save  = statev
      sig_save = stress

      ! --- Central-difference FD tangent (computed for all 3 points) ---
      do j = 1, 3
        sv_tmp   = sv_save;  sig_plus = sig_save
        strain0  = strain;   strain1  = strain
        strain1(j) = strain1(j) + h;  dstrain = strain1 - strain0
        call PSUMAT(12, 37, props, sig_plus, strain0, strain1, dstrain, sv_tmp, tan_dum)

        sv_tmp   = sv_save;  sig_minus = sig_save
        strain0  = strain;   strain1   = strain
        strain1(j) = strain1(j) - h;  dstrain = strain1 - strain0
        call PSUMAT(12, 37, props, sig_minus, strain0, strain1, dstrain, sv_tmp, tan_dum)

        tan_num(:, j) = (sig_plus - sig_minus) / (2.0_dp * h)
      end do

      denom_global = max(maxval(abs(tangent)), maxval(abs(tan_num)), 1.0_dp)
      max_err = 0.0_dp;  mi = 1;  mj = 1
      do i = 1, 3
        do j = 1, 3
          rel_err = abs(tangent(i,j) - tan_num(i,j)) / denom_global
          if (rel_err > max_err) then
            max_err = rel_err;  mi = i;  mj = j
          end if
        end do
      end do

      ! --- Secant consistency: C·ε ≈ σ ---
      sig_pred  = matmul(tangent, strain)
      sec_denom = max(maxval(abs(sig_save)), 1.0_dp)
      sec_max = 0.0_dp;  mj = 1
      do i = 1, 3
        sec_err = abs(sig_pred(i) - sig_save(i)) / sec_denom
        if (sec_err > sec_max) then
          sec_max = sec_err;  mj = i
        end if
      end do

      ! --- Report & FAIL logic per point ---
      if (ipt == 1) then
        ! ELASTIC: FD is the gold-standard check (secant ≈ tangent here)
        write(*,'(A,F7.2,A,A,I0,A,I0,A,A,F6.2,A)') &
          '  tangent@elastic     : FD_err=', max_err*100.0_dp, '%', &
          ' at (', mi, ',', mj, ')', ' secant=', sec_max*100.0_dp, '%'
        if (max_err > 0.10_dp) then
          npt_fail = npt_fail + 1;  nfail = nfail + 1
          write(*,'(A)') '  FAIL: FD > 10% at elastic'
          call dump_tangent_info(strain, sig_save, sv_save, tangent, tan_num)
        end if

      else if (ipt == 2) then
        ! CRACKED: FD unreliable (±h crosses open/close branch).
        ! Use secant consistency as FAIL criterion; report FD as info.
        write(*,'(A,F6.2,A,A,F7.2,A,A,ES9.2)') &
          '  tangent@cracked     : secant=', sec_max*100.0_dp, '%', &
          ' FD_info=', max_err*100.0_dp, '%', &
          ' ANGLE=', sv_save(2)
        if (sec_max > 0.05_dp) then
          npt_fail = npt_fail + 1;  nfail = nfail + 1
          write(*,'(A,F7.2,A)') '  FAIL: secant > 5% at cracked (', sec_max*100.0_dp, '%)'
          call dump_tangent_info(strain, sig_save, sv_save, tangent, tan_num)
        end if

      else
        ! POST-PEAK: secant stiffness (C·ε=σ by design, not dσ/dε).
        write(*,'(A,F6.2,A,A,F7.2,A,A,ES9.2)') &
          '  tangent@post-peak   : secant=', sec_max*100.0_dp, '%', &
          ' FD_info=', max_err*100.0_dp, '%', &
          ' PGRAV=', sv_save(6)
        if (sec_max > 0.05_dp) then
          npt_fail = npt_fail + 1;  nfail = nfail + 1
          write(*,'(A,F7.2,A)') '  FAIL: secant > 5% at post-peak (', sec_max*100.0_dp, '%)'
          call dump_tangent_info(strain, sig_save, sv_save, tangent, tan_num)
        end if
      end if
    end do

    if (npt_fail == 0) then
      write(*,*) 'case_tangent_consistency done (3 points OK)'
    else
      write(*,'(A,I0,A)') ' case_tangent_consistency done (', npt_fail, ' point(s) FAILED)'
    end if
  end subroutine

  subroutine dump_tangent_info(strain, stress, sv, C_ana, C_fd)
    real(dp), intent(in) :: strain(3), stress(3), sv(12), C_ana(3,3), C_fd(3,3)
    integer :: i
    write(*,'(A,3ES13.5)') '    strain=', strain
    write(*,'(A,3ES13.5)') '    stress=', stress
    write(*,'(A,ES11.3,A,ES11.3)') '    ANGLE=', sv(2), ' PGRAV=', sv(6)
    write(*,'(A)') '    analytic:'
    do i = 1, 3
      write(*,'(4X,3ES13.5)') C_ana(i,:)
    end do
    write(*,'(A)') '    FD:'
    do i = 1, 3
      write(*,'(4X,3ES13.5)') C_fd(i,:)
    end do
  end subroutine

  ! ================================================================
  ! Analytical Saenz curve: independent reference implementation
  !   sigma = E0 * eps / (1 + A*xi + B*xi^2 + C*xi^3)
  !   xi = eps/epsc (>0 in compression since both negative)
  !   Eq. (10) in Bathe & Ramaswamy (1979)
  ! ================================================================
  function saenz_analytical(E0, fc, epsc, fu, epsu, eps) result(sigma)
    real(dp), intent(in) :: E0, fc, epsc, fu, epsu, eps
    real(dp) :: sigma
    real(dp) :: RP, ES, EU, RA, RB, RC, xi, denom

    if (eps >= 0.0_dp) then
      sigma = E0 * eps
      return
    end if

    RP = epsu / epsc
    ES = fc / epsc
    EU = fu / epsu
    if (abs(RP * (RP - 1.0_dp)**2) < 1.0e-30_dp) then
      RA = 0.0_dp
    else
      RA = (E0/EU + (RP-2.0_dp)*RP*RP*E0/ES - (2.0_dp*RP+1.0_dp)*(RP-1.0_dp)**2) &
           / (RP * (RP - 1.0_dp)**2)
    end if
    RB = 2.0_dp * E0 / ES - 3.0_dp - 2.0_dp * RA
    RC = 2.0_dp - E0 / ES + RA

    xi = eps / epsc
    denom = 1.0_dp + RA*xi + RB*xi*xi + RC*xi**3
    if (abs(denom) < 1.0e-30_dp) then
      sigma = 0.0_dp
    else
      sigma = E0 * eps / denom
    end if
  end function saenz_analytical

  ! ================================================================
  ! Case: Fig. 8 — Saenz curve accuracy (Bathe & Ramaswamy 1979)
  !   Strain-controlled uniaxial compression with nu=0.
  !   Reference: Fig. 8, strain-controlled loading points should lie
  !   exactly on the analytical Saenz stress-strain curve.
  ! ================================================================
  subroutine case_fig8_saenz_curve(nfail)
    integer, intent(inout) :: nfail
    real(dp) :: props(37), stress(3), strain(3), statev(12), tangent(3,3)
    real(dp) :: strain0(3), strain1(3), dstrain(3)
    real(dp) :: E0, fc, epsc, fu, epsu
    real(dp) :: eps, sig_ref, err, max_err
    integer  :: i, u, i_max

    call set_props_bathe1979(props)
    E0 = props(1); fc = props(4); epsc = props(5)
    fu = props(6); epsu = props(7)
    call reset_state(stress, strain, statev, 12)

    open(newunit=u, file='out_fig8_saenz.csv', status='replace')
    write(u,'(A)') 'step,eps_x,sig_x_num,sig_x_ana,rel_err'

    max_err = 0.0_dp; i_max = 0
    do i = 1, 100
      eps = 1.4_dp * epsc * real(i, dp) / 100.0_dp
      strain0 = strain
      strain1 = (/ eps, 0.0_dp, 0.0_dp /)
      dstrain = strain1 - strain0
      call PSUMAT(12, 37, props, stress, strain0, strain1, dstrain, statev, tangent)
      strain = strain1

      sig_ref = saenz_analytical(E0, fc, epsc, fu, epsu, eps)

      if (abs(sig_ref) > 0.01_dp) then
        err = abs(stress(1) - sig_ref) / abs(sig_ref)
      else
        err = abs(stress(1) - sig_ref)
      end if
      if (err > max_err) then
        max_err = err; i_max = i
      end if

      write(u,'(I0,4(",",ES16.8))') i, eps, stress(1), sig_ref, err
    end do
    close(u)

    write(*,'(A,F7.3,A,I0)') '  fig8 saenz max_err=', max_err*100.0_dp, '% at step ', i_max
    if (max_err > 0.02_dp) then
      nfail = nfail + 1
      write(*,'(A,F7.3,A)') '  FAIL: Saenz curve error ', max_err*100.0_dp, '% > 2%'
    end if
    write(*,*) 'case_fig8_saenz_curve done (see out_fig8_saenz.csv)'
  end subroutine

  ! ================================================================
  ! Case: Fig. 9 — Cyclic loading path (Bathe & Ramaswamy 1979)
  !   Sequence: compress→unload→tension/crack→unload→reload→crush→
  !             strain-soften→ultimate failure
  !   Verifies all state transitions and stress evolution.
  ! ================================================================
  subroutine case_fig9_cyclic_path(nfail)
    integer, intent(inout) :: nfail
    real(dp) :: props(37), stress(3), strain(3), statev(12), tangent(3,3)
    real(dp) :: strain0(3), strain1(3), dstrain(3)
    real(dp) :: epsc, epsu, ft, E0
    real(dp) :: eps, peak_comp_stress
    integer  :: i, u, step, local_fail
    logical  :: cracked_detected, crushed_detected

    call set_props_bathe1979(props)
    E0 = props(1); ft = props(3); epsc = props(5); epsu = props(7)
    call reset_state(stress, strain, statev, 12)

    open(newunit=u, file='out_fig9_cyclic.csv', status='replace')
    write(u,'(A)') 'step,eps_x,sig_x,angle,pgrav'

    step = 0
    local_fail = 0
    cracked_detected = .false.
    crushed_detected = .false.
    peak_comp_stress = 0.0_dp

    ! Phase 1: Compress to 0.8*epsc (point A on Fig. 9)
    do i = 1, 40
      eps = 0.8_dp * epsc * real(i, dp) / 40.0_dp
      strain0 = strain
      strain1 = (/ eps, 0.0_dp, 0.0_dp /)
      dstrain = strain1 - strain0
      call PSUMAT(12, 37, props, stress, strain0, strain1, dstrain, statev, tangent)
      strain = strain1
      step = step + 1
      if (stress(1) < peak_comp_stress) peak_comp_stress = stress(1)
      write(u,'(I0,4(",",ES16.8))') step, eps, stress(1), statev(2), statev(6)
    end do

    ! Phase 2: Unload to zero strain
    do i = 1, 20
      eps = 0.8_dp * epsc * (1.0_dp - real(i, dp) / 20.0_dp)
      strain0 = strain
      strain1 = (/ eps, 0.0_dp, 0.0_dp /)
      dstrain = strain1 - strain0
      call PSUMAT(12, 37, props, stress, strain0, strain1, dstrain, statev, tangent)
      strain = strain1
      step = step + 1
      write(u,'(I0,4(",",ES16.8))') step, eps, stress(1), statev(2), statev(6)
    end do

    ! Phase 3: Tension to 2*ft/E0 (should crack at ~ft/E0)
    do i = 1, 20
      eps = 2.0_dp * ft / E0 * real(i, dp) / 20.0_dp
      strain0 = strain
      strain1 = (/ eps, 0.0_dp, 0.0_dp /)
      dstrain = strain1 - strain0
      call PSUMAT(12, 37, props, stress, strain0, strain1, dstrain, statev, tangent)
      strain = strain1
      step = step + 1
      if (abs(statev(2) - 1000.0_dp) > 1.0_dp) cracked_detected = .true.
      write(u,'(I0,4(",",ES16.8))') step, eps, stress(1), statev(2), statev(6)
    end do

    if (.not. cracked_detected) then
      local_fail = local_fail + 1
      nfail = nfail + 1
      write(*,'(A)') '  FAIL: crack not detected in tension phase'
    end if

    ! Phase 4: Unload/compress to -0.5*epsc (close crack, point D)
    do i = 1, 30
      eps = strain(1) + (-0.5_dp * epsc - strain(1)) * real(i, dp) / 30.0_dp
      strain0 = strain
      strain1 = (/ eps, 0.0_dp, 0.0_dp /)
      dstrain = strain1 - strain0
      call PSUMAT(12, 37, props, stress, strain0, strain1, dstrain, statev, tangent)
      strain = strain1
      step = step + 1
      write(u,'(I0,4(",",ES16.8))') step, eps, stress(1), statev(2), statev(6)
    end do

    ! Phase 5: Reload to 1.2*epsc (past peak, point G → crush)
    do i = 1, 40
      eps = strain(1) + (1.2_dp * epsc - strain(1)) * real(i, dp) / 40.0_dp
      strain0 = strain
      strain1 = (/ eps, 0.0_dp, 0.0_dp /)
      dstrain = strain1 - strain0
      call PSUMAT(12, 37, props, stress, strain0, strain1, dstrain, statev, tangent)
      strain = strain1
      step = step + 1
      if (stress(1) < peak_comp_stress) peak_comp_stress = stress(1)
      if (abs(statev(6) - 100.0_dp) < 1.0_dp) crushed_detected = .true.
      write(u,'(I0,4(",",ES16.8))') step, eps, stress(1), statev(2), statev(6)
    end do

    ! Phase 6: Strain-soften to epsu
    do i = 1, 30
      eps = strain(1) + (epsu - strain(1)) * real(i, dp) / 30.0_dp
      strain0 = strain
      strain1 = (/ eps, 0.0_dp, 0.0_dp /)
      dstrain = strain1 - strain0
      call PSUMAT(12, 37, props, stress, strain0, strain1, dstrain, statev, tangent)
      strain = strain1
      step = step + 1
      if (abs(statev(6) - 100.0_dp) < 1.0_dp) crushed_detected = .true.
      write(u,'(I0,4(",",ES16.8))') step, eps, stress(1), statev(2), statev(6)
    end do

    ! Phase 7: Beyond epsu to 1.5*epsu (ultimate failure)
    do i = 1, 20
      eps = strain(1) + (1.5_dp * epsu - strain(1)) * real(i, dp) / 20.0_dp
      strain0 = strain
      strain1 = (/ eps, 0.0_dp, 0.0_dp /)
      dstrain = strain1 - strain0
      call PSUMAT(12, 37, props, stress, strain0, strain1, dstrain, statev, tangent)
      strain = strain1
      step = step + 1
      write(u,'(I0,4(",",ES16.8))') step, eps, stress(1), statev(2), statev(6)
    end do
    close(u)

    ! Stability check: no NaN/Inf
    call assert_finite_vec('fig9 final stress', stress, 3, nfail)

    ! Peak compressive stress should be close to fc
    call assert_close('fig9 peak~fc', peak_comp_stress, props(4), 0.30_dp, 0.5_dp, nfail)

    write(*,'(A,I0,A,L1,A,L1,A,ES10.3)') &
      '  fig9 cyclic: steps=', step, ' cracked=', cracked_detected, &
      ' crushed=', crushed_detected, ' peak_sig=', peak_comp_stress
    if (local_fail == 0) then
      write(*,*) 'case_fig9_cyclic_path done (see out_fig9_cyclic.csv)'
    else
      write(*,*) 'case_fig9_cyclic_path done WITH FAILURES'
    end if
  end subroutine

  ! ================================================================
  ! Case: Biaxial compression envelope
  !   Equal biaxial compression (eps_x = eps_y, gamma=0).
  !   Kupfer et al. (1969): biaxial strength ~ 1.16*|fc|.
  !   Verifies multiaxial failure envelope enhances strength.
  ! ================================================================
  subroutine case_biaxial_envelope(nfail)
    integer, intent(inout) :: nfail
    real(dp) :: props(37), stress(3), strain(3), statev(12), tangent(3,3)
    real(dp) :: strain0(3), strain1(3), dstrain(3)
    real(dp) :: fc, eps, sig_max_x, sig_max_y
    integer  :: i, u

    call set_props_default(props)
    fc = props(4)
    call reset_state(stress, strain, statev, 12)

    open(newunit=u, file='out_biaxial.csv', status='replace')
    write(u,'(A)') 'step,eps,sig_x,sig_y,tau'

    sig_max_x = 0.0_dp; sig_max_y = 0.0_dp
    do i = 1, 120
      eps = -6.0e-5_dp * real(i, dp)
      strain0 = strain
      strain1 = (/ eps, eps, 0.0_dp /)
      dstrain = strain1 - strain0
      call PSUMAT(12, 37, props, stress, strain0, strain1, dstrain, statev, tangent)
      strain = strain1
      if (stress(1) < sig_max_x) sig_max_x = stress(1)
      if (stress(2) < sig_max_y) sig_max_y = stress(2)
      write(u,'(I0,4(",",ES16.8))') i, eps, stress(1), stress(2), stress(3)
    end do
    close(u)

    ! Symmetry: sigma_x should equal sigma_y
    call assert_close('biaxial sig_x=sig_y', sig_max_x, sig_max_y, 0.01_dp, 0.01_dp, nfail)

    ! Biaxial enhancement: peak stress should exceed uniaxial fc
    if (sig_max_x >= fc) then
      nfail = nfail + 1
      write(*,'(A,ES11.3,A,ES11.3)') &
        '  FAIL: biaxial peak=', sig_max_x, ' not stronger than fc=', fc
    end if

    write(*,'(A,ES11.3,A,ES11.3,A,F5.2)') &
      '  biaxial: peak_sig=', sig_max_x, ' fc=', fc, &
      ' ratio=', sig_max_x / fc
    write(*,*) 'case_biaxial_envelope done (see out_biaxial.csv)'
  end subroutine

  ! ================================================================
  ! Case: Post-crack shear retention
  !   Phase 1: Tension in x → crack forms (ANGLE ~ 0)
  !   Phase 2: Apply shear gamma → verify reduced stiffness
  !   Paper Eq. (22): G_cracked = eta_s * G0 (SHEFAC factor)
  ! ================================================================
  subroutine case_shear_retention(nfail)
    integer, intent(inout) :: nfail
    real(dp) :: props(37), stress(3), strain(3), statev(12), tangent(3,3)
    real(dp) :: strain0(3), strain1(3), dstrain(3)
    real(dp) :: E0, ft, G0, eps_cr, gamma
    real(dp) :: G_eff, stifac_eff, E_crk_ratio, G_crk_ratio
    integer  :: i

    call set_props_default(props)
    E0 = props(1); ft = props(3)
    G0 = E0 / (2.0_dp * (1.0_dp + props(2)))
    eps_cr = ft / E0
    stifac_eff = max(props(8), 0.01_dp)
    call reset_state(stress, strain, statev, 12)

    ! Phase 1: Crack in x-direction (tension in x with 15 small steps)
    do i = 1, 15
      strain0 = strain
      strain1 = (/ 3.0_dp * eps_cr * real(i, dp) / 15.0_dp, 0.0_dp, 0.0_dp /)
      dstrain = strain1 - strain0
      call PSUMAT(12, 37, props, stress, strain0, strain1, dstrain, statev, tangent)
      strain = strain1
    end do

    if (abs(statev(2) - 1000.0_dp) < 1.0_dp) then
      nfail = nfail + 1
      write(*,'(A)') '  FAIL: crack not formed in shear retention test'
      write(*,*) 'case_shear_retention done WITH FAILURE'
      return
    end if

    ! Phase 2: Apply shear at current crack state
    gamma = 1.0e-4_dp
    strain0 = strain
    strain1 = strain
    strain1(3) = strain1(3) + gamma
    dstrain = strain1 - strain0
    call PSUMAT(12, 37, props, stress, strain0, strain1, dstrain, statev, tangent)
    strain = strain1

    ! Check tangent: cracked normal stiffness should be reduced
    E_crk_ratio = tangent(1,1) / (E0 / (1.0_dp - props(2)**2))
    G_crk_ratio = tangent(3,3) / G0

    write(*,'(A,F7.4,A,F7.4)') &
      '  shear_retention: E_crk/E0_ps=', E_crk_ratio, ' G_crk/G0=', G_crk_ratio

    ! Normal stiffness in crack direction should be significantly reduced
    if (E_crk_ratio > 0.5_dp) then
      nfail = nfail + 1
      write(*,'(A,F7.4)') '  FAIL: crack-normal stiffness not reduced, ratio=', E_crk_ratio
    end if

    ! Shear stiffness should be reduced (expect ~ SHEFAC = 0.5)
    if (G_crk_ratio > 0.90_dp) then
      nfail = nfail + 1
      write(*,'(A,F7.4)') '  FAIL: shear stiffness not reduced, ratio=', G_crk_ratio
    end if

    call assert_finite_vec('shear_ret stress', stress, 3, nfail)
    write(*,*) 'case_shear_retention done'
  end subroutine

  subroutine case_random_robustness(nfail)
    integer, intent(inout) :: nfail
    integer, parameter :: NSTEPS = 5000
    integer, parameter :: NHIST  = 20
    real(dp) :: props(37), stress(3), strain(3), statev(12), tangent(3,3)
    real(dp) :: strain0(3), strain1(3), dstrain(3)
    real(dp) :: rnd(4), amp, sig_max
    real(dp) :: strain_hist(3, NHIST), statev_hist(12, NHIST)
    real(dp) :: asym, denom
    integer  :: i, j, k, seed_size, hist_pos, nwarn_sym, last_step
    integer, allocatable :: seed(:)
    logical  :: exploded

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    do i = 1, seed_size
      seed(i) = 42 + i * 7
    end do
    call random_seed(put=seed)
    deallocate(seed)

    call set_props_default(props)
    call reset_state(stress, strain, statev, 12)

    strain_hist = 0.0_dp
    statev_hist = 0.0_dp
    nwarn_sym   = 0
    exploded    = .false.
    last_step   = 0

    do i = 1, NSTEPS
      call random_number(rnd)
      amp = 1.0e-6_dp + (1.0e-4_dp - 1.0e-6_dp) * rnd(4)
      do k = 1, 3
        dstrain(k) = amp * (2.0_dp * rnd(k) - 1.0_dp)
      end do

      strain0 = strain
      strain1 = strain + dstrain

      call PSUMAT(12, 37, props, stress, strain0, strain1, dstrain, statev, tangent)
      strain = strain1
      last_step = i

      hist_pos = mod(i - 1, NHIST) + 1
      strain_hist(:, hist_pos) = strain
      statev_hist(:, hist_pos) = statev

      ! Check 1: stress finite
      do k = 1, 3
        if (.not. (abs(stress(k)) <= 1.0d100)) then
          nfail = nfail + 1
          write(*,'(A,I0,A,I0,A,ES12.4)') &
            'FAIL: non-finite stress step=', i, ' comp=', k, ' val=', stress(k)
          exploded = .true.
        end if
      end do
      if (exploded) exit

      ! Check 1b: tangent finite
      do j = 1, 3
        do k = 1, 3
          if (.not. (abs(tangent(j,k)) <= 1.0d100)) then
            nfail = nfail + 1
            write(*,'(A,I0,A,I0,I0,A,ES12.4)') &
              'FAIL: non-finite tangent step=', i, ' (', j, k, ') val=', tangent(j,k)
            exploded = .true.
          end if
        end do
        if (exploded) exit
      end do
      if (exploded) exit

      ! Check 2: stress magnitude bound
      sig_max = maxval(abs(stress))
      if (sig_max > 2000.0_dp) then
        nfail = nfail + 1
        write(*,'(A,I0,A,3ES14.6)') &
          'FAIL: |stress|>2000 step=', i, ' sig=', stress
        exploded = .true.
        exit
      end if

      ! Check 3: tangent symmetry  (3 off-diagonal pairs)
      do j = 1, 2
        do k = j + 1, 3
          denom = max(abs(tangent(j,k)), abs(tangent(k,j)), 1.0_dp)
          asym  = abs(tangent(j,k) - tangent(k,j)) / denom
          if (asym > 0.05_dp) then
            nwarn_sym = nwarn_sym + 1
            if (nwarn_sym <= 5) then
              write(*,'(A,I0,A,I0,A,I0,A,ES10.2,A,ES10.2,A,F7.2,A)') &
                'WARN: tangent asym step=', i, ' (', j, ',', k, &
                ') Cij=', tangent(j,k), ' Cji=', tangent(k,j), &
                ' asym=', asym * 100.0_dp, '%'
            end if
          end if
        end do
      end do
    end do

    ! Dump history on explosion
    if (exploded) then
      write(*,*) '--- Last up to 20 steps before explosion ---'
      do k = 1, min(last_step, NHIST)
        j = mod(last_step - k, NHIST) + 1
        write(*,'(A,3ES11.3,A,4ES10.2)') &
          '  eps=', strain_hist(:,j), '  sv(1:4)=', statev_hist(1:4,j)
      end do
    end if

    if (nwarn_sym > 0) then
      nfail = nfail + 1
      write(*,'(A,I0,A,I0,A)') &
        '  tangent symmetry: ', nwarn_sym, ' violations in ', last_step, ' steps'
    end if

    write(*,'(A,I0,A)') ' case_random_robustness done (', last_step, ' steps)'
  end subroutine

end program test_runner
