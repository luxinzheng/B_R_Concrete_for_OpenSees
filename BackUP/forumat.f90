!===============================================================================
! Module: MHEOM2D
! Description: HEOM2D 2D Concrete Constitutive Model
!              (Hypo-Elastic Orthotropic Model)
! Based on: Darwin-Pecknold equivalent uniaxial strain approach
!           Kupfer biaxial failure surface
!           Noguchi rotating crack model
! Note: All empirical formulas are unit-independent via auto-detection.
!
! [FIX] Change log:
!   1. Fixed parameter mapping: compatible with standard PlaneStressUserMaterial
!      7-parameter interface (fc, ft, fcu, epsc0, epscu, epstu, stc)
!   2. Fixed state variable collision: initialization flag no longer overwrites
!      Cc(3,3). Data uses slots 1..38, flag at statev(nstatevs).
!   3. Fixed array out-of-bounds: props(8..16) accessed only if nprops allows.
!   4. Fixed negative Eo: auto-computed as 2*fcp/|ecu|, always positive.
!   5. Fixed tangent modulus floor: uses DABS(Eo) to guarantee positive minimum.
!   6. Removed redundant eud storage (eu=eud at save time), saving 2 slots.
!===============================================================================
MODULE MHEOM2D
  IMPLICIT NONE

  ! Number of material properties (max)
  INTEGER, PARAMETER :: NPROPS_HEOM2D = 80
  ! Number of state variables (max)
  INTEGER, PARAMETER :: NSTATV_HEOM2D = 60
  ! [FIX] Actual data state variables (excluding initialization flag)
  INTEGER, PARAMETER :: NSTV_DATA = 38
  ! Constants
  REAL(8), PARAMETER :: PI = 3.141592653589793D0
  REAL(8), PARAMETER :: RAD2DEG = 57.29577951308232D0
  REAL(8), PARAMETER :: DEG2RAD = 0.01745329252D0

  ! HEOM2D material data type
  TYPE :: typ_HEOM2D
    ! ========== Concrete material parameters ==========
    REAL(8) :: ft       ! Tensile strength (+)
    REAL(8) :: fcp      ! Compressive strength (+), note: fc = -fcp
    REAL(8) :: Eo       ! Initial elastic modulus
    REAL(8) :: ecu      ! Peak compressive strain (-)
    REAL(8) :: epcu     ! Ultimate compressive strain (-)
    REAL(8) :: epcus    ! Residual stress at epcu
    REAL(8) :: ebu      ! Ultimate tensile strain (+)

    ! ========== Model options ==========
    INTEGER :: iconc    ! Compression curve: 0-Saenz, 1-Popovics, 2-Shah-Fafitis, 3-Muguruma
    INTEGER :: icrack   ! Tension softening: 0-Shirai, 1-Vecchio-Collins
    INTEGER :: ired     ! Compression reduction: 0-Noguchi, 1-Vecchio-Collins
    INTEGER :: imod     ! Model mode: 0-Full model(with shear slip), 1-Basic model
    INTEGER :: ndivs    ! Tension softening subdivisions

    ! ========== Rebar data (up to 10 layers) ==========
    INTEGER :: ns       ! Number of rebar layers
    REAL(8) :: Eso(10)  ! Rebar elastic modulus
    REAL(8) :: eh(10)   ! Hardening onset strain
    REAL(8) :: Esh(10)  ! Hardening modulus
    REAL(8) :: fys(10)  ! Yield strength
    REAL(8) :: alfas(10)! Rebar angle (degrees)
    REAL(8) :: ro(10)   ! Reinforcement ratio

    ! ========== Crack spacing parameters ==========
    REAL(8) :: sx       ! Crack spacing in x-direction
    REAL(8) :: sy       ! Crack spacing in y-direction
    REAL(8) :: Ql       ! Shear transfer parameter

    ! ========== State indicators ==========
    INTEGER :: ntst     ! Phase indicator
    INTEGER :: iul(2)   ! Stress-strain curve position indicator
    INTEGER :: irek     ! Crack direction rotation indicator
    REAL(8) :: sni      ! Current Poisson's ratio
    REAL(8) :: oan      ! Initial principal angle
    REAL(8) :: angl     ! Previous principal angle
    REAL(8) :: eu(2)    ! Equivalent uniaxial strains
    REAL(8) :: eud(2)   ! Updated equivalent uniaxial strains
    REAL(8) :: etc(2)   ! Tangent moduli

    ! ========== Stresses and strains ==========
    REAL(8) :: G(3)     ! Current stress (X-Y coordinates)
    REAL(8) :: tep(3)   ! Current total strain (X-Y coordinates)
    REAL(8) :: ted(3)   ! Updated total strain (X-Y coordinates)
    REAL(8) :: tedc(3)  ! Total strain without slip
    REAL(8) :: tedcd(3) ! Updated total strain without slip
    REAL(8) :: de(3)    ! Strain increment
    REAL(8) :: tepc(3)  ! Principal total strain
    REAL(8) :: tepd(3)  ! Updated principal total strain
    REAL(8) :: Gd(3)    ! Updated concrete stress (X-Y coordinates)
    REAL(8) :: Gds(3)   ! Rebar stress contribution
    REAL(8) :: Gsl(3)   ! Shear slip stress contribution
    REAL(8) :: eslip(3) ! Slip strain

    ! ========== Stiffness matrices ==========
    REAL(8) :: C(3,3)   ! Total tangent stiffness matrix
    REAL(8) :: Cc(3,3)  ! Concrete tangent stiffness matrix

  END TYPE typ_HEOM2D

CONTAINS

  !-----------------------------------------------------------------------------
  ! Helper: get unit conversion factor (stress -> MPa)
  ! Detects unit system from elastic modulus magnitude
  !-----------------------------------------------------------------------------
  FUNCTION get_sfac(Eval) RESULT(sfac)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: Eval
    REAL(8) :: sfac
    IF (Eval .GT. 1.D8) THEN
      sfac = 1.D-6       ! Pa to MPa
    ELSE IF (Eval .GT. 1.D5) THEN
      sfac = 1.D-3       ! kPa to MPa
    ELSE
      sfac = 1.D0        ! already MPa
    END IF
  END FUNCTION get_sfac

  !-----------------------------------------------------------------------------
  ! [FIX] Initialize material properties
  ! Now compatible with standard PlaneStressUserMaterial 7-parameter interface:
  !   props(1) = fc    : compressive strength (positive)
  !   props(2) = ft    : tensile strength (positive)
  !   props(3) = fcu   : crushing/residual stress (negative)
  !   props(4) = epsc0 : peak compressive strain (negative)
  !   props(5) = epscu : ultimate compressive strain (negative)
  !   props(6) = epstu : ultimate tensile strain (positive)
  !   props(7) = stc   : shear transfer coefficient
  ! Extended parameters (optional, props 8+):
  !   props(8)  = iconc, props(9)  = icrack, props(10) = ired
  !   props(11) = imod,  props(12) = ndivs
  !   props(13) = sx,    props(14) = sy,     props(15) = Ql
  !   props(16) = ns,    props(17+)= rebar data
  !-----------------------------------------------------------------------------
  SUBROUTINE HEOM2D_InitProps(mat, props, nprops)
    TYPE(typ_HEOM2D), INTENT(INOUT) :: mat
    INTEGER, INTENT(IN) :: nprops
    REAL(8), INTENT(IN) :: props(nprops)
    INTEGER :: i, ich

    ! Initialize all rebar arrays to zero
    mat%Eso   = 0.D0
    mat%eh    = 0.D0
    mat%Esh   = 0.D0
    mat%fys   = 0.D0
    mat%alfas = 0.D0
    mat%ro    = 0.D0

    ! [FIX] Standard 7-parameter interface mapping
    mat%fcp   = DABS(props(1))             ! Compressive strength (positive)
    mat%ft    = DABS(props(2))             ! Tensile strength (positive)
    mat%epcus = props(3)                   ! Residual/crushing stress (negative)
    mat%ecu   = props(4)                   ! Peak compressive strain (negative)
    mat%epcu  = props(5)                   ! Ultimate compressive strain (negative)
    IF (nprops .GE. 6) THEN
      mat%ebu = DABS(props(6))             ! Ultimate tensile strain (positive)
    ELSE
      mat%ebu = 0.D0
    END IF
    IF (nprops .GE. 7) THEN
      mat%Ql  = props(7)                   ! Shear transfer coefficient
    ELSE
      mat%Ql  = 0.D0
    END IF

    ! [FIX] Ensure correct signs for strains and stress
    IF (mat%ecu  .GT. 0.D0) mat%ecu  = -mat%ecu
    IF (mat%epcu .GT. 0.D0) mat%epcu = -mat%epcu
    IF (mat%epcus .GT. 0.D0) mat%epcus = -mat%epcus

    ! [FIX] Auto-compute elastic modulus: Eo = 2*fcp / |ecu|
    IF (DABS(mat%ecu) .GT. 1.D-20) THEN
      mat%Eo = 2.D0 * mat%fcp / DABS(mat%ecu)
    ELSE
      mat%Eo = 0.D0
    END IF

    ! Fill remaining defaults via conc_param if needed
    IF (mat%Eo .LE. 0.D0 .OR. DABS(mat%ecu) .LT. 1.D-20) THEN
      CALL conc_param(mat%fcp, mat%Eo, mat%ecu, mat%epcu, mat%epcus, &
                      mat%ft/mat%fcp, mat%ebu, ich)
    END IF

    ! [FIX] Final safety: ensure Eo > 0
    IF (mat%Eo .LE. 0.D0) THEN
      mat%Eo = 5000.D0 * DSQRT(mat%fcp * get_sfac(mat%fcp)) / get_sfac(mat%fcp)
    END IF

    ! Ensure ebu is reasonable
    IF (mat%ebu .LE. mat%ft / mat%Eo) THEN
      mat%ebu = 20.D0 * mat%ft / mat%Eo
    END IF

    ! Ensure |epcus| <= fcp
    IF (DABS(mat%epcus) .GT. mat%fcp) THEN
      mat%epcus = -0.2D0 * mat%fcp
    END IF

    ! [FIX] Safe default model options (no out-of-bounds access)
    mat%iconc  = 0      ! Saenz compression curve
    mat%icrack = 0      ! Shirai tension softening
    mat%ired   = 0      ! Noguchi compression reduction
    mat%imod   = 1      ! Basic model (no shear slip)
    mat%ndivs  = 10     ! Tension softening subdivisions
    mat%sx     = 0.D0
    mat%sy     = 0.D0
    mat%ns     = 0

    ! [FIX] Override only if extended parameters are actually provided
    IF (nprops .GE. 8)  mat%iconc  = NINT(props(8))
    IF (nprops .GE. 9)  mat%icrack = NINT(props(9))
    IF (nprops .GE. 10) mat%ired   = NINT(props(10))
    IF (nprops .GE. 11) mat%imod   = NINT(props(11))
    IF (nprops .GE. 12) mat%ndivs  = NINT(props(12))
    IF (nprops .GE. 13) mat%sx     = props(13)
    IF (nprops .GE. 14) mat%sy     = props(14)
    IF (nprops .GE. 15) mat%Ql     = props(15)

    ! Rebar data (if provided)
    IF (nprops .GE. 16) THEN
      mat%ns = NINT(props(16))
      IF (mat%ns .GT. 0 .AND. mat%ns .LE. 10) THEN
        DO i = 1, mat%ns
          IF (16 + i*6 .LE. nprops) THEN
            mat%Eso(i)   = props(16 + (i-1)*6 + 1)
            mat%eh(i)    = props(16 + (i-1)*6 + 2)
            mat%Esh(i)   = props(16 + (i-1)*6 + 3)
            mat%fys(i)   = props(16 + (i-1)*6 + 4)
            mat%alfas(i) = props(16 + (i-1)*6 + 5)
            mat%ro(i)    = props(16 + (i-1)*6 + 6)
          END IF
        END DO
      ELSE
        mat%ns = 0
      END IF
    END IF

  END SUBROUTINE HEOM2D_InitProps

  !-----------------------------------------------------------------------------
  ! [FIX] Restore state from state variable array
  ! Layout: data in slots 1..NSTV_DATA (=38), flag at statev(nstatv)
  ! Removed eud storage (eu = eud at save time, saving 2 slots)
  !-----------------------------------------------------------------------------
  SUBROUTINE HEOM2D_GetState(mat, statev, nstatv)
    TYPE(typ_HEOM2D), INTENT(INOUT) :: mat
    INTEGER, INTENT(IN) :: nstatv
    REAL(8), INTENT(IN) :: statev(nstatv)
    INTEGER :: idx

    idx = 1

    ! Phase indicator
    mat%ntst = NINT(statev(idx)); idx = idx + 1       ! 1

    ! Position indicators
    mat%iul(1) = NINT(statev(idx)); idx = idx + 1     ! 2
    mat%iul(2) = NINT(statev(idx)); idx = idx + 1     ! 3
    mat%irek   = NINT(statev(idx)); idx = idx + 1     ! 4

    ! Poisson's ratio and angles
    mat%sni  = statev(idx); idx = idx + 1              ! 5
    mat%oan  = statev(idx); idx = idx + 1              ! 6
    mat%angl = statev(idx); idx = idx + 1              ! 7

    ! [FIX] Equivalent uniaxial strains (eu only; eud initialized = eu)
    mat%eu(1)  = statev(idx); idx = idx + 1            ! 8
    mat%eu(2)  = statev(idx); idx = idx + 1            ! 9
    mat%eud    = mat%eu   ! eud will be recomputed by concr

    ! Tangent moduli
    mat%etc(1) = statev(idx); idx = idx + 1            ! 10
    mat%etc(2) = statev(idx); idx = idx + 1            ! 11

    ! Stress (3 components)
    mat%G(1:3) = statev(idx:idx+2); idx = idx + 3     ! 12-14

    ! Strain (3 components)
    mat%tep(1:3) = statev(idx:idx+2); idx = idx + 3   ! 15-17

    ! Total strain without slip
    mat%tedc(1:3) = statev(idx:idx+2); idx = idx + 3  ! 18-20

    ! Principal total strain
    mat%tepc(1:3) = statev(idx:idx+2); idx = idx + 3  ! 21-23

    ! Slip strain
    mat%eslip(1:3) = statev(idx:idx+2); idx = idx + 3 ! 24-26

    ! Shear slip stress
    mat%Gsl(1:3) = statev(idx:idx+2); idx = idx + 3   ! 27-29

    ! Concrete tangent stiffness matrix (9 components)
    mat%Cc(1,1:3) = statev(idx:idx+2); idx = idx + 3  ! 30-32
    mat%Cc(2,1:3) = statev(idx:idx+2); idx = idx + 3  ! 33-35
    mat%Cc(3,1:3) = statev(idx:idx+2); idx = idx + 3  ! 36-38
    ! idx = NSTV_DATA + 1 = 39 (flag is at nstatv, not touched here)

  END SUBROUTINE HEOM2D_GetState

  !-----------------------------------------------------------------------------
  ! [FIX] Save current state to state variable array
  ! Writes only slots 1..NSTV_DATA (=38), preserves flag at nstatv
  !-----------------------------------------------------------------------------
  SUBROUTINE HEOM2D_SetState(mat, statev, nstatv)
    TYPE(typ_HEOM2D), INTENT(IN) :: mat
    INTEGER, INTENT(IN) :: nstatv
    REAL(8), INTENT(INOUT) :: statev(nstatv)
    INTEGER :: idx

    ! [FIX] Only zero data slots, preserve flag slot at nstatv
    statev(1:NSTV_DATA) = 0.D0

    idx = 1

    ! Phase indicator
    statev(idx) = DBLE(mat%ntst); idx = idx + 1

    ! Position indicators
    statev(idx) = DBLE(mat%iul(1)); idx = idx + 1
    statev(idx) = DBLE(mat%iul(2)); idx = idx + 1
    statev(idx) = DBLE(mat%irek);   idx = idx + 1

    ! Poisson's ratio and angles
    statev(idx) = mat%sni;  idx = idx + 1
    statev(idx) = mat%oan;  idx = idx + 1
    statev(idx) = mat%angl; idx = idx + 1

    ! [FIX] Equivalent uniaxial strains (eu only, already = eud)
    statev(idx) = mat%eu(1);  idx = idx + 1
    statev(idx) = mat%eu(2);  idx = idx + 1

    ! Tangent moduli
    statev(idx) = mat%etc(1); idx = idx + 1
    statev(idx) = mat%etc(2); idx = idx + 1

    ! Stress (updated)
    statev(idx:idx+2) = mat%Gd(1:3); idx = idx + 3

    ! Strain (updated)
    statev(idx:idx+2) = mat%ted(1:3); idx = idx + 3

    ! Total strain without slip (updated)
    statev(idx:idx+2) = mat%tedcd(1:3); idx = idx + 3

    ! Principal total strain (updated)
    statev(idx:idx+2) = mat%tepd(1:3); idx = idx + 3

    ! Slip strain
    statev(idx:idx+2) = mat%eslip(1:3); idx = idx + 3

    ! Shear slip stress
    statev(idx:idx+2) = mat%Gsl(1:3); idx = idx + 3

    ! Concrete tangent stiffness matrix
    statev(idx:idx+2) = mat%Cc(1,1:3); idx = idx + 3
    statev(idx:idx+2) = mat%Cc(2,1:3); idx = idx + 3
    statev(idx:idx+2) = mat%Cc(3,1:3); idx = idx + 3
    ! idx = NSTV_DATA + 1 = 39 (flag at nstatv is NOT touched)

  END SUBROUTINE HEOM2D_SetState

  !-----------------------------------------------------------------------------
  ! Initialize state (first call)
  !-----------------------------------------------------------------------------
  SUBROUTINE HEOM2D_InitState(mat)
    TYPE(typ_HEOM2D), INTENT(INOUT) :: mat

    mat%ntst  = 1
    mat%iul   = 1
    mat%irek  = 0
    mat%sni   = 0.2D0
    mat%oan   = 0.D0
    mat%angl  = 0.D0
    mat%eu    = 0.D0
    mat%eud   = 0.D0
    mat%etc   = mat%Eo
    mat%G     = 0.D0
    mat%tep   = 0.D0
    mat%ted   = 0.D0
    mat%tedc  = 0.D0
    mat%tedcd = 0.D0
    mat%tepc  = 0.D0
    mat%tepd  = 0.D0
    mat%Gd    = 0.D0
    mat%Gds   = 0.D0
    mat%Gsl   = 0.D0
    mat%eslip = 0.D0

    ! Initial elastic stiffness matrix
    mat%Cc = 0.D0
    mat%Cc(1,1) = mat%Eo / (1.D0 - mat%sni**2)
    mat%Cc(2,2) = mat%Cc(1,1)
    mat%Cc(1,2) = mat%sni * mat%Cc(1,1)
    mat%Cc(2,1) = mat%Cc(1,2)
    mat%Cc(3,3) = mat%Eo / (2.D0 * (1.D0 + mat%sni))
    mat%C = mat%Cc

  END SUBROUTINE HEOM2D_InitState

  !=============================================================================
  ! Core HEOM2D subroutine
  !=============================================================================
  SUBROUTINE HEOM2D_Core(ft, fcp, ebu, ecu, epcu, epcus, Eo, eu, eud, iul, etc, &
                    ntst, oan, angl, sni, G, tep, ted, tedc, tedcd, de, tepc, tepd, &
                    Gd, Gds, Gsl, C, Cc, ndivs, iconc, icrack, ired, imod, ns, &
                    Eso, eh, Esh, fys, alfas, ro, sx, sy, Ql, eslip, irek)
    IMPLICIT NONE
    INTEGER(KIND=4), PARAMETER    :: maxit = 100
    INTEGER(KIND=4)               :: ntst, ndivs, iconc, icrack, ired, imod, i, j, ij, &
                                     ns, is, iter, irek, mfail, msw
    INTEGER(KIND=4), DIMENSION(2) :: iul
    REAL(KIND=8), PARAMETER       :: eps0 = 0.01D0
    REAL(KIND=8)                  :: ft, fcp, ebu, ecu, epcu, epcus, Eo, sni, ang1, a, b, G2, &
                                     cs, ss, angl, oan, Qs, ec1, ec2, fc1, fc2, de1cr, vci, &
                                     sx, sy, Ql, s, w, ds, dsa, gsa, gsb, gs, &
                                     de1crp, diff, Em1, Em2, ecr, raz, &
                                     ros, fyh, Fp, Fx, denom_s, sfac_w
    REAL(KIND=8), DIMENSION(3)    :: G, tep, ted, tedc, tedcd, tepc, tepl, tepd, dG, de, tsi, &
                                     Gd, Gds, Gsl, eslip, deslip, eslipp
    REAL(KIND=8), DIMENSION(4)    :: psi
    REAL(KIND=8), DIMENSION(3,3)  :: C, Cc, Cr, tra, ua, ub, uc, uti, Cc1
    REAL(KIND=8), DIMENSION(2)    :: etc, eu, eud
    REAL(KIND=8), DIMENSION(10)   :: Eso, eh, Esh, fys, alfas, ro, Est, es, fs, esy, Qn, escr, &
                                     fscr, escrp, fscrp

    ! Unit conversion factor for Walraven equation
    sfac_w = get_sfac(Eo)

    ! I) Contribution from concrete elements
    DO i=1, 3
      dG(i) = 0.D0
      DO j=1, 3
        dG(i) = dG(i) + Cc(i, j)*de(j)
      END DO
    END DO
    DO i=1, 3
      tsi(i)    = G(i) + dG(i)
      ted(i)    = tep(i) + de(i)     ! Total strains with slip
      tedcd(i)  = tedc(i) + de(i)    ! Total strains without slip
      eslipp(i) = eslip(i)
    END DO

    ! Initialize tepd(3) to avoid undefined value
    tepd(3) = 0.D0
    tepl(3) = 0.D0

    ! 2: In principal stress direction
    CALL prin(tsi(1), tsi(2), tsi(3), psi(1), psi(2), psi(3), psi(4))
    ang1 = psi(4)
    IF(ntst.EQ.1) oan = ang1

    ! Calculate principal strains and their increments
    CALL chgcrd(tra, ang1, 1)
    DO i=1, 2
      tepl(i) = tepc(i)
      tepc(i) = 0.D0
      DO j=1, 3
        tepc(i) = tepc(i) + tra(i, j)*tedcd(j)
      END DO
    END DO

    ! Compute ros and fyh from rebar data before calling concr
    ros = 0.D0
    fyh = 0.D0
    IF (ns > 0) THEN
      DO is = 1, ns
        ros = ros + ro(is)
      END DO
      IF (ros > 0.D0) THEN
        DO is = 1, ns
          fyh = fyh + ro(is) * fys(is)
        END DO
        fyh = fyh / ros
      END IF
    END IF

    ! Calculate tangent modulus of concrete
    CALL concr(ft, fcp, Eo, ebu, ecu, epcu, epcus, ang1, oan, angl, tepl, tepc, eu, eud, &
               psi, iul, sni, etc, ndivs, iconc, icrack, ired, ns, ros, fyh, irek)

    ! Form D-matrix of concrete
    DO i=1, 3
      DO j=1, 3
        ua(i, j)  = 0.D0
        ub(i, j)  = 0.D0
        uc(i, j)  = 0.D0
        uti(i, j) = 0.D0
        Cc(i, j)  = 0.D0
      END DO
    END DO

    ! Equivalent uniaxial strain model by Darwin & Pecknold
    b  = 1.D0/(1.D0-sni*sni)
    IF (DABS(tepc(1)-tepc(2)) .GT. 1.D-20) THEN
        G2 = 0.5D0 * (psi(1)-psi(2)) / (tepc(1)-tepc(2))
    ELSE
        G2 = Eo / (2.D0*(1.D0+sni))
    END IF

    ! [FIX] Ensure non-negative tangent moduli using DABS(Eo)
    IF (etc(1) .LT. 0.D0) etc(1) = DABS(Eo) * 0.001D0
    IF (etc(2) .LT. 0.D0) etc(2) = DABS(Eo) * 0.001D0

    uti(1, 1) = etc(1)*b
    uti(1, 2) = sni*(etc(1)*etc(2))**0.5D0*b
    uti(2, 1) = uti(1,2)
    uti(2, 2) = etc(2)*b
    uti(3, 3) = (1.D0-sni**2) * G2

    CALL chgcrd(ua, -ang1, 0)
    CALL chgcrd(ub,  ang1, 1)
    DO i=1, 3
      DO j=1, 3
        DO ij=1, 3
          uc(i, j) = uc(i, j) + ua(i, ij)*uti(ij, j)
        END DO
      END DO
    END DO
    DO i=1, 3
      DO j=1, 3
        DO ij=1, 3
          Cc(i, j) = Cc(i, j) + uc(i, ij)*ub(ij, j)
        END DO
      END DO
    END DO

    ! Calculate stresses of concrete in global coordinates
    cs = DCOS(ang1*DEG2RAD)
    ss = DSIN(ang1*DEG2RAD)

    psi(3) = (psi(1)-psi(2))*0.5D0
    Gd(1) = psi(1)*cs**2 + psi(2)*ss**2
    Gd(2) = psi(1)*ss**2 + psi(2)*cs**2
    Gd(3) = (psi(1)-psi(2))*ss*cs

    ! Save updated values
    angl = psi(4)
    DO i=1, 2
      tepd(i) = tepc(i)
      tepc(i) = tepl(i)
    END DO

    ! Initialize rebar contributions
    DO i=1, 3
      Gds(i) = 0.D0
      DO j=1, 3
        Cr(i, j)  = 0.D0
      END DO
    END DO

    ! Initialize Cc1 and Gsl
    DO i=1, 3
      Gsl(i) = 0.D0
      DO j=1, 3
        Cc1(i, j) = 0.D0
      END DO
    END DO

    ! Initialize slip strain increments
    DO i=1, 3
      deslip(i) = 0.D0
    END DO

    IF(ns.EQ.0) THEN
      DO i=1, 3
        DO j=1, 3
          C(i, j) = Cc(i, j)
        END DO
      END DO
      ! Keep eslip unchanged
      DO i=1, 3
        eslip(i) = eslipp(i)
      END DO
      RETURN
    END IF

    ! II) Contribution from Re-bar elements
    raz = tedcd(1)-tedcd(2)
    IF(raz.GT.0.D0.AND.raz.LT.1.D-20) raz = 1.D-20
    IF(raz.LT.0.D0.AND.DABS(raz).LT.1.D-20) raz = -1.D-20
    Qs = 0.5D0*DATAN2(tedcd(3), raz)   ! [rad]

    DO is=1, ns
      es(is)  = 0.5D0*(ted(1)+ted(2)) + &
                0.5D0*(ted(1)-ted(2))*DCOS(2.D0*alfas(is)*DEG2RAD) + &
                0.5D0*ted(3)*DSIN(2.D0*alfas(is)*DEG2RAD)
      Qn(is)  = Qs - alfas(is)*DEG2RAD   ! [rad]
      esy(is)  = fys(is) / Eso(is)

      IF(DABS(es(is)).LT.esy(is)) THEN
        Est(is) = Eso(is)
        fs(is)  = Eso(is)*es(is)
      ELSE
        fs(is)  = fys(is)*DSIGN(1.D0,es(is)) + &
                  Esh(is)*(DABS(es(is))-esy(is))*DSIGN(1.D0,es(is))
        Est(is) = Esh(is)
      END IF

      DO i=1, 3
        DO j=1, 3
          ua(i, j)  = 0.D0
          ub(i, j)  = 0.D0
          uc(i, j)  = 0.D0
          uti(i, j) = 0.D0
        END DO
      END DO

      uti(1, 1) = ro(is)*Est(is)

      CALL chgcrd(ua, -alfas(is), 0)
      CALL chgcrd(ub,  alfas(is), 1)
      DO i=1, 3
        DO j=1, 3
          DO ij=1, 3
            uc(i, j) = uc(i, j) + ua(i, ij)*uti(ij, j)
          END DO
        END DO
      END DO
      DO i=1, 3
        DO j=1, 3
          DO ij=1, 3
            Cr(i, j) = Cr(i, j) + uc(i, ij)*ub(ij, j)
          END DO
        END DO
      END DO

      cs = DCOS(alfas(is)*DEG2RAD)
      ss = DSIN(alfas(is)*DEG2RAD)

      Gds(1) = Gds(1) + ro(is)*fs(is)*cs**2
      Gds(2) = Gds(2) + ro(is)*fs(is)*ss**2
      Gds(3) = Gds(3) + ro(is)*fs(is)*ss*cs

    END DO

    ! Total [D] matrix
    DO i=1, 3
      DO j=1, 3
        C(i, j) = Cc(i, j) + Cr(i, j)
      END DO
    END DO

    IF(imod.EQ.1) THEN
      ! Keep eslip unchanged for basic model
      DO i=1, 3
        eslip(i) = eslipp(i)
      END DO
      RETURN
    END IF

    ! III) Influence of shear-transfer across the crack
    IF(psi(1).GT.psi(2)) THEN
      fc1 = psi(1)
      fc2 = psi(2)
      ec1 = tepd(1)
      ec2 = tepd(2)
      Em1 = etc(1)
      IF(iul(1).EQ.4) Em1 = 0.D0
      Em2 = etc(2)
      IF(iul(2).EQ.4) Em2 = 0.D0
      msw = 1
    ELSE
      fc1 = psi(2)
      fc2 = psi(1)
      ec1 = tepd(2)
      ec2 = tepd(1)
      Em1 = etc(2)
      IF(iul(2).EQ.4) Em1 = 0.D0
      Em2 = etc(1)
      IF(iul(1).EQ.4) Em2 = 0.D0
      msw = 2
    END IF

    ! Influence of the shear slip on the crack surface
    ecr = ft/Eo
    IF(ec1.GT.ecr.AND.fc1.GT.0.D0) THEN

      ! Calculate "de1cr" via Newton-Raphson iteration
      mfail = 0
      iter = 0
      de1cr = 0.D0
      de1crp = de1cr

      88 iter = iter + 1
      IF(iter.GT.maxit) THEN
        mfail = 1
        GO TO 99
      END IF

      Fx = -fc1
      Fp = 0.D0
      DO i=1, ns
        escr(i) = es(i) + de1cr*DCOS(Qn(i))**2
        escrp(i) = DCOS(Qn(i))**2
        IF(DABS(escr(i)).LT.esy(i)) THEN
          fscr(i) = Eso(i)*escr(i)
          fscrp(i) = Eso(i)*escrp(i)
        ELSE
          fscr(i)  = fys(i)*DSIGN(1.D0,escr(i)) + &
                     Esh(i)*(DABS(escr(i))-esy(i))*DSIGN(1.D0,escr(i))
          fscrp(i) = Esh(i)*escrp(i)
        END IF
        Fx = Fx + ro(i)*(fscr(i)-fs(i))*DCOS(Qn(i))**2
        Fp = Fp + ro(i)*fscrp(i)*DCOS(Qn(i))**2
      END DO
      IF (DABS(Fp) .GT. 1.D-30) THEN
         de1cr = de1cr - Fx/Fp
      END IF

      IF (DABS(de1cr) .GT. 1.D-20) THEN
        diff = DABS(de1cr-de1crp)/DABS(de1cr)
      ELSE
        diff = 0.D0
      END IF
      de1crp = de1cr
      IF(diff.GT.eps0) GO TO 88

      ! Calculate induced shear force on the crack surface "vci"
      99 IF(mfail.EQ.0) THEN
        vci = 0.D0
        DO i=1, ns
          vci = vci + ro(i)*(fscr(i)-fs(i))*DCOS(Qn(i))*DSIN(Qn(i))
        END DO
      ELSE
        vci = 0.D0
      END IF

      ! Vecchio's method for calculating "s" and "w"
      denom_s = sy*DABS(DSIN(Qs)) + sx*DABS(DCOS(Qs))
      IF (denom_s .LT. 1.D-20) denom_s = 1.D-20
      s = sx*sy / denom_s
      w = ec1*DABS(s)*1.D3  ! Convert to [mm]

      IF(w.GT.1.D-6) THEN
        ! Walraven's equation (all stress quantities must be in MPa)
        ds = vci*sfac_w / (1.8D0*w**(-.8D0) + &
             (0.234D0*w**(-.707D0) - 0.2D0)*(1.17D0*fcp*sfac_w) )
        dsa = ds*1.D-3   ! mm to m
        gsa = dsa/s
      ELSE
        gsa = 0.D0
      END IF

      gsb = ted(3)*DCOS(2.D0*Qs) + (ted(2)-ted(1))*DSIN(2.D0*Qs)
      gs = gsb
      IF(DABS(gsa).GT.DABS(gs)) gs = gsa

      ! Slip strain vector
      eslip(1) = -0.5D0*gs*DSIN(2.D0*Qs)
      eslip(2) =  0.5D0*gs*DSIN(2.D0*Qs)
      eslip(3) =  gs*DCOS(2.D0*Qs)

      DO i=1, 3
        deslip(i) = eslip(i) - eslipp(i)
      END DO

      ! Build Cc1() matrix for shear slip contribution
      DO i=1, 3
        DO j=1, 3
          ua(i, j)  = 0.D0
          ub(i, j)  = 0.D0
          uc(i, j)  = 0.D0
          uti(i, j) = 0.D0
          Cc1(i, j) = 0.D0
        END DO
      END DO

      b  = 1.D0/(1.D0-sni*sni)

      IF (DABS(ec1-ec2) .GT. 1.D-20) THEN
          G2 = 0.5D0 * (fc1-fc2) / (ec1-ec2)
      ELSE
          G2 = 0.D0
      END IF

      uti(1, 1) = Em1*b
      uti(1, 2) = sni*(DABS(Em1*Em2))**0.5D0*b
      uti(2, 1) = uti(1,2)
      uti(2, 2) = Em2*b
      uti(3, 3) = G2 * (1.D0-sni**2)

      CALL chgcrd(ua, -ang1, 0)
      CALL chgcrd(ub,  ang1, 1)
      DO i=1, 3
        DO j=1, 3
          DO ij=1, 3
            uc(i, j) = uc(i, j) + ua(i, ij)*uti(ij, j)
          END DO
        END DO
      END DO
      DO i=1, 3
        DO j=1, 3
          DO ij=1, 3
            Cc1(i, j) = Cc1(i, j) + uc(i, ij)*ub(ij, j)
          END DO
        END DO
      END DO

    ELSE

      DO i=1, 3
        eslip(i)  = eslipp(i)
        deslip(i) = 0.D0
      END DO

    END IF

    ! Shear slip stress contribution
    DO i=1, 3
      DO j=1, 3
        Gsl(i) = Gsl(i) - Cc1(i, j)*deslip(j)
      END DO
    END DO

    RETURN
  END SUBROUTINE HEOM2D_Core


  !-----------------------------------------------------------------------------
  ! Coordinate transformation matrix
  ! N = 0 : STRESS , N = 1 : STRAIN
  !-----------------------------------------------------------------------------
  SUBROUTINE chgcrd(a, tha, n)
    IMPLICIT NONE
    INTEGER(KIND=4)              :: n
    REAL(KIND=8), DIMENSION(3,3) :: A
    REAL(KIND=8)                 :: tha, cs, ss

    cs = DCOS(tha*DEG2RAD)
    ss = DSIN(tha*DEG2RAD)

    IF(n.EQ.0) THEN
      A(1, 1) =  cs**2
      A(1, 2) =  ss**2
      A(1, 3) =  2.D0*cs*ss
      A(2, 1) =  A(1, 2)
      A(2, 2) =  A(1, 1)
      A(2, 3) = -A(1, 3)
      A(3, 1) = -cs*ss
      A(3, 2) =  cs*ss
      A(3, 3) =  cs**2-ss**2
    ELSE IF(n.EQ.1) THEN
      A(1, 1) =  cs**2
      A(1, 2) =  ss**2
      A(1, 3) =  cs*ss
      A(2, 1) =  A(1, 2)
      A(2, 2) =  A(1, 1)
      A(2, 3) = -A(1, 3)
      A(3, 1) = -cs*ss*2.D0
      A(3, 2) =  cs*ss*2.D0
      A(3, 3) =  cs**2-ss**2
    END IF

    RETURN
  END SUBROUTINE chgcrd


  !-----------------------------------------------------------------------------
  ! Swap two integers
  !-----------------------------------------------------------------------------
  SUBROUTINE rvers1(idata1, idata2)
    IMPLICIT NONE
    INTEGER(KIND=4) :: idata1, idata2, idummy
    idummy = idata1
    idata1 = idata2
    idata2 = idummy
    RETURN
  END SUBROUTINE rvers1


  !-----------------------------------------------------------------------------
  ! Swap two reals
  !-----------------------------------------------------------------------------
  SUBROUTINE rvers2(data1, data2)
    IMPLICIT NONE
    REAL(KIND=8) :: data1, data2, dummy
    dummy = data1
    data1 = data2
    data2 = dummy
    RETURN
  END SUBROUTINE rvers2


  !-----------------------------------------------------------------------------
  ! Constitutive model for concrete
  ! Based on equivalent uniaxial stress-strain relations
  ! Kupfer biaxial failure surface
  !-----------------------------------------------------------------------------
  SUBROUTINE concr(ft, fcp, Eo, ebu, ecu, epcu, epcus, ang1, oan, angl, tepl, tepc, &
                   eu, eud, psi, iul, sni, etc, ndivs, iconc, icrack, ired, ns, ros, fyh, irek)
    IMPLICIT NONE
    INTEGER(KIND=4), DIMENSION(2) :: iul
    INTEGER(KIND=4)               :: iul1, iul2, ndivs, iconc, icrack, ired, ns, irek
    REAL(KIND=8)                  :: ft, fcp, Eo, ebu, ecu, ang1, oan, angl, sni, abo, ab1, &
                                     tep1, tep2, eul1, eul2, dep1, dep2, tepl1, tepl2, &
                                     deu1, deu2, eu1, eu2, sn1, sn2, sc1, sc2, ec1, ec2, &
                                     fc, beta, beta1, alp, r1, r2, et1, et2, epcu, epcus, &
                                     Css, Cd, beta2, ros, fyh, dum
    REAL(KIND=8), DIMENSION(2)    :: eu, eud, etc
    REAL(KIND=8), DIMENSION(3)    :: tepl, tepc
    REAL(KIND=8), DIMENSION(4)    :: psi

    fc = -fcp

    ! Check rotation of principal direction (Noguchi's Rotating Crack Model)
    abo = DABS(ang1-oan)
    ab1 = DABS(ang1-angl)
    irek = 0
    IF(abo.GT.45.D0.AND.abo.LT.135.D0) THEN
      oan = DABS(oan+90.D0)
      IF(oan.GT.180.D0) oan = DABS(oan-180.D0)
      ab1 = DABS(ang1-angl-90.D0)
      ang1 = DABS(ang1+90.D0)
      IF(ang1.GT.180.D0) ang1 = DABS(ang1-180.D0)
      angl = DABS(angl+90.D0)
      IF(angl.GT.180.D0) angl = DABS(angl-180.D0)
      CALL rvers1(iul(1), iul(2))
      CALL rvers2(psi(1), psi(2))
      CALL rvers2(eu(1), eu(2))
      CALL rvers2(etc(1), etc(2))
      CALL rvers2(tepl(1), tepl(2))
      CALL rvers2(tepc(1), tepc(2))
      irek = 1
    ELSE
      IF(ab1.LT.45.D0.OR.ab1.GT.135.D0) THEN
        CONTINUE
      ELSE
        oan = DABS(oan+90.D0)
        IF(oan.GT.180.D0) oan = DABS(oan-180.D0)
        ab1 = DABS(ang1-angl-90.D0)
        ang1 = DABS(ang1+90.D0)
        IF(ang1.GT.180.D0) ang1 = DABS(ang1-180.D0)
        angl = DABS(angl+90.D0)
        IF(angl.GT.180.D0) angl = DABS(angl-180.D0)
        CALL rvers1(iul(1), iul(2))
        CALL rvers2(psi(1), psi(2))
        CALL rvers2(eu(1), eu(2))
        CALL rvers2(etc(1), etc(2))
        CALL rvers2(tepl(1), tepl(2))
        CALL rvers2(tepc(1), tepc(2))
        irek = 1
      END IF
    END IF

    sn1   = psi(1)
    sn2   = psi(2)
    eu1   = eu(1)
    eu2   = eu(2)
    iul1  = iul(1)
    iul2  = iul(2)
    et1   = etc(1)
    et2   = etc(2)
    tepl1 = tepl(1)
    tepl2 = tepl(2)
    tep1  = tepc(1)
    tep2  = tepc(2)

    ! Calculate equivalent uniaxial strains (Darwin-Pecknold/Chen)
    eul1 = eu1
    eul2 = eu2
    dep1 = tep1-tepl1
    dep2 = tep2-tepl2

    IF (et1 .GT. 0.D0 .AND. et2 .GT. 0.D0) THEN
        deu1 = (dep1+sni*(et2/et1)**0.5D0*dep2)/(1.D0-sni**2)
        deu2 = (dep2+sni*(et1/et2)**0.5D0*dep1)/(1.D0-sni**2)
    ELSE
        deu1 = dep1
        deu2 = dep2
    END IF

    ! Total uniaxial strains
    eu1 = eul1+deu1
    eu2 = eul2+deu2

    beta1 = 1.D0
    beta2 = 1.D0

    ! Initialize sc1, sc2
    sc1 = fc
    sc2 = fc

    ! Calculate sc1, sc2 using Kupfer's surface
    IF(iul1.EQ.8.OR.iul2.EQ.8.OR.iul1.EQ.9.OR.iul2.EQ.9) THEN
      GO TO 1
    END IF

    IF(iul1.EQ.4.OR.iul2.EQ.4) THEN
      IF(iul1.EQ.4) THEN
        IF(eu1.LT.0.D0) THEN
          sc2 = fc
        ELSE
          SELECT CASE (ired)
            CASE (0)
              ! Noguchi's equation
              beta = -eu1/ecu
              beta1 = 0.27D0+0.96D0*beta**0.167D0
              IF(beta1.LE.1.D0) beta1 = 1.D0
              beta2 = 1.D0
            CASE (1)
              ! Vecchio & Collins (1986)
              Css = 1.D0
              IF (DABS(eu2) .GT. 1.D-20) THEN
                 dum = -eu1/eu2 - 0.28D0
                 IF (dum .GT. 0.D0) THEN
                   Cd = 0.35D0 * dum**0.8D0
                 ELSE
                   Cd = 0.D0
                 END IF
              ELSE
                 Cd = 0.D0
              END IF
              beta1 = 1.D0+Css*Cd
              IF(beta1.LE.1.D0) beta1 = 1.D0
              beta2 = beta1
            CASE DEFAULT
              STOP 'Input error: ired!'
          END SELECT
          sc2 = fc/beta1
        END IF
        sc1 = sc2
      END IF

      IF(iul2.EQ.4) THEN
        IF(eu2.LT.0.D0) THEN
          sc1 = fc
        ELSE
          SELECT CASE (ired)
            CASE (0)
              ! Noguchi's equation
              beta = -eu2/ecu
              beta1 = 0.27D0+0.96D0*beta**0.167D0
              IF(beta1.LE.1.D0) beta1 = 1.D0
              beta2 = 1.D0
            CASE (1)
              ! Vecchio & Collins (1986)
              Css = 1.D0
              IF (DABS(eu1) .GT. 1.D-20) THEN
                 dum = -eu2/eu1 - 0.28D0
                 IF (dum .GT. 0.D0) THEN
                   Cd = 0.35D0 * dum**0.8D0
                 ELSE
                   Cd = 0.D0
                 END IF
              ELSE
                 Cd = 0.D0
              END IF
              beta1 = 1.D0+Css*Cd
              IF(beta1.LE.1.D0) beta1 = 1.D0
              beta2 = beta1
            CASE DEFAULT
              STOP 'Input error: ired!'
          END SELECT

          IF(iul1.EQ.4.AND.iul2.EQ.4) THEN
            IF(fc/beta1.LT.sc1) sc1 = fc/beta1
          ELSE
            sc1 = fc/beta1
          END IF
        END IF
        sc2 = sc1
      END IF
      GO TO 1
    END IF

    IF(sn1.LT.0.D0.AND.sn2.LT.0.D0) THEN
      IF(sn1.GT.sn2) THEN
        alp = sn1/sn2
        sc2 = (1.D0+3.65D0*alp)*fc/(1.D0+alp)**2
        sc1 = sc2
        IF(sc1.GT.0.65D0*fc) sc1 = 0.65D0*fc
        IF(sc2.GT.0.65D0*fc) sc2 = 0.65D0*fc
      ELSE
        alp = sn2/sn1
        sc1 = (1.D0+3.65D0*alp)*fc/(1.D0+alp)**2
        sc2 = sc1
        IF(sc1.GT.0.65D0*fc) sc1 = 0.65D0*fc
        IF(sc2.GT.0.65D0*fc) sc2 = 0.65D0*fc
      END IF
      GO TO 1
    END IF

    IF(sn1.GE.0.D0.AND.sn2.LT.0.D0) THEN
      alp = sn1/sn2
      IF(sn1.GT.ft) alp = ft/sn2
      IF(alp.NE.-1.D0) THEN
        sc2 = (1.D0+3.28D0*alp)*fc/(1.D0+alp)**2
        sc1 = sc2
        IF(sc1.GT.0.65D0*fc) sc1 = 0.65D0*fc
        IF(sc2.GT.0.65D0*fc) sc2 = 0.65D0*fc
      ELSE
        sc1 = fc
        sc2 = fc
      END IF
      GO TO 1
    END IF

    IF(sn2.GE.0.D0.AND.sn1.LT.0.D0) THEN
      alp = sn2/sn1
      IF(sn2.GT.ft) alp = ft/sn1
      IF(alp.NE.-1.D0) THEN
        sc1 = (1.D0+3.28D0*alp)*fc/(1.D0+alp)**2
        sc2 = sc1
        IF(sc1.GT.0.65D0*fc) sc1 = 0.65D0*fc
        IF(sc2.GT.0.65D0*fc) sc2 = 0.65D0*fc
      ELSE
        sc1 = fc
        sc2 = fc
      END IF
      GO TO 1
    END IF

    ! Tension - tension (no crack)
    sc1 = fc
    sc2 = fc

    1 CONTINUE

    IF(DABS(sc1).LT.-fc) THEN
      ec1 = ecu*(-1.6D0*(sc1/fc)**3+2.25D0*(sc1/fc)**2+0.35D0*sc1/fc)
    ELSE
      ec1 = ecu*(3.15D0*sc1/fc-2.15D0)
    END IF
    IF(DABS(sc2).LT.-fc) THEN
      ec2 = ecu*(-1.6D0*(sc2/fc)**3+2.25D0*(sc2/fc)**2+0.35D0*sc2/fc)
    ELSE
      ec2 = ecu*(3.15D0*sc2/fc-2.15D0)
    END IF
    ec1 = ec1/beta2
    ec2 = ec2/beta2

    ! Calculate tangent modulus of elasticity
    CALL secon(iul1, eu1, sn1, et1, ec1, sc1, fc, ft, Eo, epcu, epcus, ebu, ndivs, &
               iconc, icrack, ns)
    CALL secon(iul2, eu2, sn2, et2, ec2, sc2, fc, ft, Eo, epcu, epcus, ebu, ndivs, &
               iconc, icrack, ns)

    ! Update Poisson's ratio
    IF(iul1.NE.1.AND.iul2.NE.1) THEN
      sni = 0.2D0
    ELSE IF(sn1.LT.0.D0) THEN
      sni = 0.2D0
    ELSE IF(sn2.GE.0.D0) THEN
      sni = 0.2D0
    ELSE
      r1 = sn2/fc
      r2 = sn1/ft
      IF(r1.GT.1.D0) r1 = 1.D0
      IF(r2.GT.1.D0) r2 = 1.D0
      sni = 0.2D0+0.6D0*r1**4+0.4D0*r2**4
    END IF
    IF(iul1.EQ.4.OR.iul2.EQ.4.OR.iul1.EQ.8.OR.iul2.EQ.8.OR.iul1.EQ.9.OR.iul2.EQ.9) sni = 0.D0
    IF(sni.GT.0.5D0) sni = 0.5D0

    ! Return updated values
    eud(1)  = eu1
    eud(2)  = eu2
    psi(1)  = sn1
    psi(2)  = sn2
    iul(1)  = iul1
    iul(2)  = iul2
    etc(1)  = et1
    etc(2)  = et2

    RETURN
  END SUBROUTINE concr


  !-----------------------------------------------------------------------------
  ! Secant/tangent modulus computation for uniaxial stress-strain curves
  ! Unit-independent via auto-detection
  ! [FIX] Tangent modulus floor uses DABS(eo) to guarantee positive minimum
  !-----------------------------------------------------------------------------
  SUBROUTINE secon(iul, eu, sn, et, ec, sc, fc, ft, eo, epcu, epcus, ebu, ndivs, &
                   iconc, icrack, ns)
    IMPLICIT NONE
    INTEGER(KIND=4)               :: iul, ndivs, iconc, icrack, ns, iii
    REAL(KIND=8)                  :: eu, sn, et, ec, sc, fc, ft, eo, epcu, epcus, ebu, &
                                     ecr, dum, dum1, eex1, ssx1, es, Cebu, n, k, em, qc, eun
    REAL(KIND=8), DIMENSION(1001) :: gg, ee, Ecc
    REAL(KIND=8)                  :: sfac, sc_mpa, exp_val

    ! Unit conversion factor
    sfac = get_sfac(eo)
    sc_mpa = DABS(sc) * sfac   ! Peak stress in MPa (positive)

    ! Crack strain
    ecr = ft/eo

    ! Limit strain
    IF(ebu.LT.ecr) ebu = ecr

    IF(eu.GE.ec.AND.eu.LE.ecr) THEN
      iul = 1
      IF(eu.GE.0.D0) THEN
        et = eo
        sn = et*eu
      ELSE
        SELECT CASE (iconc)
          CASE (0)
            ! Saenz's equation (unit-independent)
            dum = 1.D0+(eo*ec/sc-2.D0)*eu/ec+(eu/ec)**2
            IF (DABS(dum) .LT. 1.D-20) dum = 1.D-20
            et = eo*(1.D0-(eu/ec)**2)/dum**2
            sn = eo*eu/dum
          CASE (1)
            ! Popovics's equation (ascending branch)
            n = 0.8D0 + sc_mpa / 17.D0
            IF (n .LT. 1.01D0) n = 1.01D0
            IF (n .GT. 20.D0) n = 20.D0
            k = 1.D0
            IF (DABS(eu/ec) .GT. 1.D-20) THEN
              exp_val = n*k*DLOG(DABS(eu/ec))
              IF (exp_val .GT. 500.D0) THEN
                dum = DEXP(500.D0)
              ELSE IF (exp_val .LT. -500.D0) THEN
                dum = (n-1.D0)
              ELSE
                dum = (n-1.D0) + DEXP(exp_val)
              END IF
            ELSE
              dum = (n-1.D0)
            END IF
            IF (DABS(dum) .LT. 1.D-20) dum = 1.D-20
            sn = sc * n * eu/ec / dum
            et = sc*n/(ec*dum**2) * ( dum - eu/ec*n*k*DEXP((n*k-1.D0)*DLOG(DABS(eu/ec)+1.D-30)) )
          CASE (2)
            ! Fafitis-Shah's equation (unit-independent)
            dum1 = eu/ec
            dum = eo*ec/sc
            IF (DABS(1.D0-dum1) .GT. 1.D-20) THEN
              et = eo*(1.D0-dum1)**(dum-1.D0)
              sn = sc*(1.D0-(1.D0-dum1)**dum)
            ELSE
              et = eo*0.01D0
              sn = sc
            END IF
          CASE (3)
            ! Modified Muguruma & Watanabe's equation (unit-independent)
            em = DABS(ec)
            qc = DABS(sc)
            eun = DABS(eu)
            sn = eo*eun + (qc-eo*em)/em**2 *eun**2
            sn = -sn
            et = eo + 2.D0*eun*(qc-eo*em)/em**2
          CASE DEFAULT
              STOP 'Input error: iconc!'
        END SELECT
        ! [FIX] Floor uses DABS(eo)
        IF(et.LT.DABS(eo)/100.D0) et = DABS(eo)/100.D0
      END IF

    ELSE IF(eu.GT.ecr) THEN
      iul = 4
      IF(eu.LE.ebu) THEN
        SELECT CASE (icrack)
          CASE (0)
            IF(ndivs.EQ.0) THEN
                et = DABS(eo)
                iii = 0
            ELSE IF(ndivs.GT.1000) THEN
              et = DABS(eo)*0.01D0
              iii = 0
            ELSE IF(ndivs.GT.0.AND.ndivs.LE.1000) THEN
              Cebu = ebu/ecr
              CALL tstif(ecr, Cebu, ft, ndivs, gg, ee, Ecc)
              CALL findGE(ndivs, gg, ee, Ecc, eu, sn, et, iii)
            END IF

            ! Shirai's equation for stress
            eex1 = (eu-ecr)/(ebu-ecr)
            ssx1 = 1.D0-2.748D0*eex1+2.654D0*eex1**2-0.906D0*eex1**3
            sn = ssx1*ft
            ! [FIX] Floor uses DABS(eo)
            IF(et.LT.DABS(eo)*0.01D0) et = DABS(eo)*0.01D0
            IF(iii.EQ.1) sn = 0.0D0

          CASE (1)
            ! Vecchio & Collins (MCFT)
            sn = ft / ( 1.D0+(200.D0*eu)**0.5D0 )
            et = DABS(eo)

          CASE DEFAULT
            STOP ' Input error: icrack!'
        END SELECT
        IF(sn.LT.0.D0) sn = 0.D0
      ELSE
        IF(ndivs.EQ.0) THEN
            et = DABS(eo)
        ELSE
          et = DABS(eo)*0.01D0
        END IF
        sn = 0.D0
      END IF

    ELSE IF(eu.GE.epcu.AND.eu.LT.ec) THEN
      iul = 8
      SELECT CASE (iconc)
        CASE (1)
          ! Popovics descending branch
          n = 0.8D0 + sc_mpa / 17.D0
          k = 0.67D0 + sc_mpa / 62.D0
          IF (n .LT. 1.01D0) n = 1.01D0
          IF (n .GT. 20.D0) n = 20.D0
          IF (k .GT. 5.D0) k = 5.D0
          IF (DABS(eu/ec) .GT. 1.D-20) THEN
            exp_val = n*k*DLOG(DABS(eu/ec))
            IF (exp_val .GT. 500.D0) THEN
              sn = epcus
            ELSE IF (exp_val .LT. -500.D0) THEN
              dum = (n-1.D0)
              sn = sc * n * eu/ec / dum
            ELSE
              dum = (n-1.D0) + DEXP(exp_val)
              IF (DABS(dum) .LT. 1.D-20) dum = 1.D-20
              sn = sc * n * eu/ec / dum
            END IF
          ELSE
            sn = 0.D0
          END IF
          IF(ndivs.EQ.0) THEN
              et = DABS(eo)
          ELSE
            et = DABS(eo)*0.01D0
          END IF
          IF(sn.GT.epcus) sn = epcus
        CASE DEFAULT
          IF (DABS(epcu-ec) .GT. 1.D-20) THEN
            es = (epcus-sc) / (epcu-ec)
          ELSE
            es = 0.D0
          END IF
          sn = sc+es*(eu-ec)
          IF(ndivs.EQ.0) THEN
              et = DABS(eo)
          ELSE
            et = DABS(eo)*0.01D0
          END IF
          IF(sn.GT.epcus) sn = epcus
      END SELECT

    ELSE IF(eu.LT.epcu) THEN
      iul = 9
      IF(ndivs.EQ.0) THEN
          et = DABS(eo)
      ELSE
        et = DABS(eo)*0.01D0
      END IF
      sn = epcus
    END IF

    RETURN
  END SUBROUTINE secon


  !-----------------------------------------------------------------------------
  ! Compute tension softening stiffness array
  !-----------------------------------------------------------------------------
  SUBROUTINE tstif(ecr, Cebu, ft, ndivs, GG, ee, Ec)
    IMPLICIT NONE
    INTEGER(KIND=4)            :: ndivs, i
    REAL(KIND=8)               :: ecr, Cebu, ft, ebu, dxx, xx
    REAL(KIND=8), DIMENSION(*) :: ee, GG, Ec

    ebu = Cebu*ecr
    dxx = 1.0D0/(ndivs*1.0D0)
    DO i=1, ndivs+1
      xx = (i-1)*dxx
      GG(i) = ft*(1.D0-2.748D0*xx+2.654D0*xx**2-0.906D0*xx**3)
      ee(i) = xx*(ebu-ecr)+ecr
    END DO
    DO i=1, ndivs-1
      IF (DABS(ee(i+1)-ee(i)) .GT. 1.D-30) THEN
        Ec(i) = (GG(i+1)-GG(i+2))/(ee(i+1)-ee(i))
      ELSE
        Ec(i) = 0.D0
      END IF
    END DO
    IF (ndivs .GE. 2) THEN
      Ec(ndivs) = Ec(ndivs-1)
    ELSE IF (ndivs .EQ. 1) THEN
      Ec(ndivs) = 0.D0
    END IF

    RETURN
  END SUBROUTINE tstif


  !-----------------------------------------------------------------------------
  ! Find tangent stiffness for given strain in tension softening curve
  !-----------------------------------------------------------------------------
  SUBROUTINE findGE(ndivs, GG, ee, Ec, eu, Gu, E, iii)
    IMPLICIT NONE
    INTEGER(KIND=4)            :: ndivs, iii, i
    REAL(KIND=8)               :: eu, Gu, E
    REAL(KIND=8), DIMENSION(*) :: ee, GG, Ec

    iii = 0
    DO i=1, ndivs
      IF(eu.GE.ee(i).AND.eu.LE.ee(i+1)) THEN
        E = Ec(i)
        IF(i.NE.ndivs) THEN
          Gu = GG(i+2)+E*(eu-ee(i))
        ELSE
          Gu = E*(eu-ee(i))
        END IF
        RETURN
      END IF
    END DO
    iii = 1

    RETURN
  END SUBROUTINE findGE


  !-----------------------------------------------------------------------------
  ! Auto-compute concrete parameters from fcp (unit-independent)
  !-----------------------------------------------------------------------------
  SUBROUTINE conc_param(fcp, Eo, ecu, epcu, epcus, ftp, ebu, ich)
    IMPLICIT NONE
    INTEGER(KIND=4) :: ich
    REAL(KIND=8)    :: fcp, Eo, ecu, epcu, epcus, ftp, ebu, ecr, ft
    REAL(KIND=8)    :: sfac, fcp_mpa

    ich = 0
    IF(fcp.EQ.0.D0) THEN
      ich = 999
      RETURN
    END IF

    ! Auto-detect unit from fcp magnitude
    IF(fcp .GT. 1.D6) THEN
      sfac = 1.D-6       ! Pa to MPa
    ELSE IF(fcp .GT. 1.D3) THEN
      sfac = 1.D-3       ! kPa to MPa
    ELSE
      sfac = 1.D0        ! already MPa
    END IF
    fcp_mpa = fcp * sfac

    IF(Eo.EQ.0.D0) THEN
      Eo = 5000.D0 * DSQRT(fcp_mpa) / sfac
    END IF

    IF(ecu.EQ.0.D0) THEN
      ecu = -2.D0*fcp/Eo
    END IF

    IF(epcu.EQ.0.D0) THEN
      epcu = 4.D0*ecu
    END IF

    IF(epcus.EQ.0.D0) THEN
      epcus = -0.2D0*fcp
    END IF

    IF(ftp.EQ.0.D0) THEN
      ft = 0.33D0 * DSQRT(fcp_mpa) / sfac
      ftp = ft/fcp
    END IF

    IF(ebu.EQ.0.D0) THEN
      ft = ftp*fcp
      ecr = ft/Eo
      ebu = 20.D0*ecr
    END IF

    RETURN
  END SUBROUTINE conc_param


  !-----------------------------------------------------------------------------
  ! Confinement effect on concrete parameters (unit-independent)
  !-----------------------------------------------------------------------------
  SUBROUTINE confine_sub(ros, fyh, hp, Sh, ecu, fcp, epcu, epcus, Eo, ftp, ebu)
    IMPLICIT NONE
    REAL(KIND=8) :: ros, fyh, hp, Sh, ecu, epcu, epcus, a1, a2, k, Z, &
                    ecc, ef, fcc, fc02, ecp, fcp, ft, ftp, ebu, ecr, Eo
    REAL(KIND=8) :: sfac, fcc_mpa

    IF(ros.NE.0.D0) THEN
      ! Unit conversion
      sfac = get_sfac(Eo)

      ecp = DABS(ecu)
      k   = 1.D0 + ros*fyh/fcp
      ecc = ecp*k
      fcc = fcp*k

      ! Empirical formula uses MPa
      fcc_mpa = fcc * sfac
      a1 = (3.D0 + 0.29D0*fcc_mpa) / (145.D0*fcc_mpa - 1000.D0)
      a2 = 0.75D0*ros*(hp/Sh)**0.5D0
      Z  = 0.5D0 / (a1+a2-ecc)
      IF (DABS(Z) .LT. 1.D-20) Z = 1.D-20
      ef = (0.8D0+Z*ecc) / Z
      IF(ef.LT.4.D0*ecc) ef = 4.D0*ecc
      fc02 = 0.2D0*fcc

      ecu   = -ecc
      fcp   = fcc
      epcu  = -ef
      epcus = -fc02

      ! Correction of tension concrete strength (unit-independent)
      ft  = 0.33D0 * DSQRT(fcp*sfac) / sfac
      ftp = ft/fcp
      ecr = ft/Eo
      ebu = 20.D0*ecr
    END IF

    RETURN
  END SUBROUTINE confine_sub

  !=============================================================================
  ! Principal stress computation (Mohr's circle)
  !=============================================================================
  SUBROUTINE prin(s11, s22, s12, p1, p2, p3, p4)
    IMPLICIT NONE
    REAL(8), INTENT(IN)  :: s11, s22, s12
    REAL(8), INTENT(OUT) :: p1, p2, p3, p4
    REAL(8) :: savg, sdif, R

    savg = 0.5D0 * (s11 + s22)
    sdif = 0.5D0 * (s11 - s22)
    R = DSQRT(sdif**2 + s12**2)

    p1 = savg + R
    p2 = savg - R
    p3 = R

    IF (DABS(sdif) .LT. 1.D-20 .AND. DABS(s12) .LT. 1.D-20) THEN
      p4 = 0.D0
    ELSE
      p4 = 0.5D0 * DATAN2(s12, sdif) * RAD2DEG
    END IF

  END SUBROUTINE prin

END MODULE MHEOM2D

!===============================================================================
! [FIX] OpenSees User Material Interface - HEOM2D
! Fixed: state variable flag no longer collides with Cc(3,3)
! Data uses statev(1..38), flag at statev(nstatevs)
! Minimum nstatevs = 39 (recommended: 40)
!===============================================================================
SUBROUTINE PSUMAT(nstatevs, nprops, props, stress, strain0, strain1, dstrain, &
                  statev, tangent)
  !DEC$ ATTRIBUTES DLLEXPORT :: PSUMAT
  !DEC$ ATTRIBUTES ALIAS:'PSUMAT' :: PSUMAT
  USE MHEOM2D
  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN) :: nstatevs, nprops
  REAL(8), INTENT(IN) :: props(nprops)
  REAL(8), INTENT(IN) :: strain0(3), strain1(3), dstrain(3)
  REAL(8), INTENT(INOUT) :: stress(3)
  REAL(8), INTENT(INOUT) :: statev(nstatevs)
  REAL(8), INTENT(OUT) :: tangent(3,3)

  ! Local variables
  TYPE(typ_HEOM2D) :: mat
  INTEGER :: i, is_initialized

  ! 1. Initialize material properties
  CALL HEOM2D_InitProps(mat, props, nprops)

  ! [FIX] 2. Check initialization flag at last slot (separate from data)
  is_initialized = NINT(statev(nstatevs))

  IF (is_initialized .NE. 12345) THEN
    ! First call: initialize state
    CALL HEOM2D_InitState(mat)
    statev(nstatevs) = 12345.D0
  ELSE
    ! Restore previous state (reads slots 1..NSTV_DATA only)
    CALL HEOM2D_GetState(mat, statev, nstatevs)
  END IF

  ! 3. Set current strain and stress
  mat%tep = strain0
  mat%de  = dstrain
  mat%G   = stress

  ! 4. Call HEOM2D core computation
  CALL HEOM2D_Core(mat%ft, mat%fcp, mat%ebu, mat%ecu, mat%epcu, &
              mat%epcus, mat%Eo, mat%eu, mat%eud, mat%iul, &
              mat%etc, mat%ntst, mat%oan, mat%angl, mat%sni, &
              mat%G, mat%tep, mat%ted, mat%tedc, mat%tedcd, &
              mat%de, mat%tepc, mat%tepd, mat%Gd, mat%Gds, &
              mat%Gsl, mat%C, mat%Cc, mat%ndivs, mat%iconc, &
              mat%icrack, mat%ired, mat%imod, mat%ns, mat%Eso, &
              mat%eh, mat%Esh, mat%fys, mat%alfas, mat%ro, &
              mat%sx, mat%sy, mat%Ql, mat%eslip, mat%irek)

  ! 5. Compute total stress = concrete + rebar + shear slip
  DO i = 1, 3
    stress(i) = mat%Gd(i) + mat%Gds(i) + mat%Gsl(i)
  END DO

  ! Return tangent stiffness matrix
  tangent = mat%C

  ! Update phase indicator (first step done)
  IF (mat%ntst .EQ. 1) mat%ntst = 0

  ! Update equivalent uniaxial strains
  mat%eu = mat%eud

  ! [FIX] 6. Save state variables (writes slots 1..NSTV_DATA only)
  CALL HEOM2D_SetState(mat, statev, nstatevs)
  ! Re-write flag at last slot (SetState does NOT touch nstatevs)
  statev(nstatevs) = 12345.D0

  RETURN
END SUBROUTINE PSUMAT