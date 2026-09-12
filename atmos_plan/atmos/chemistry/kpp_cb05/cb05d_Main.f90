! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  cb05d_Main.f90 - SCENARIO driver, DUMMY-REACTION forcing (replaces the
!  operator-splitting driver cb05s_Main.f90).
!
!  Emissions, ventilation and background inflow are NOT applied outside the
!  integrator any more: they are ordinary reactions in the mechanism
!    FSRC = X      : SRC_X(TEMP)       emission        [molec/cm3/s]
!    X    = FSNK   : KDILK(TEMP)       ventilation     [1/s]
!    FBCK = X      : BGIN_X(TEMP)      background in   [molec/cm3/s]
!  (see mk_dummy.py).  The driver only refreshes the forcing globals
!  (E_*, K_DIL, BG_*) once per outer step, before Update_RCONST.
!
!  Why this matters for the Moebius enumeration: with operator splitting the
!  inter-step coupling introduced by the numerical splitting would be mixed
!  into the interaction spectrum along with genuine chemical non-linearity,
!  and the source/sink terms would stay outside the stoichiometric matrix.
!
!  Switches (namelist RUNCFG):
!    NOXW/VOCW/COW = per-window (4 x 6h) multipliers on the NOx / VOC / CO
!                    anthropogenic emission: 1 = on, 0 = switched off.
!                    CO is the third intervention sector (12-player game).
!    DUMMYF  = .FALSE.  -> every forcing rate is 0: closed box
!                          (used to re-run the G3 conservation gate)
!    CONSTJ  = .TRUE.   -> freeze J/T/C_M at their midday values: autonomous
!                          closed box (used for the structural regression)
!
!  Solver-failure probe uses KPP semantics: IERR < 0 is a failure.
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM cb05_Driver

      USE cb05_Model
      USE cb05_Initialize, ONLY: Initialize

      IMPLICIT NONE

!~~~> run configuration (namelist)
      REAL(kind=dp)    :: TEND_S, DT_S
      CHARACTER(len=64):: STATE_IN
      REAL(kind=dp)    :: EMIS_NOX, EMIS_VOC, EMIS_ISO, EMIS_CO
      REAL(kind=dp)    :: NOX_FACTOR, VOC_FACTOR
      REAL(kind=dp)    :: NOXW(4), VOCW(4), COW(4)
      LOGICAL          :: WKEND, DUMMYF, CONSTJ
      REAL(kind=dp)    :: C_H2O
      REAL(kind=dp)    :: KDIL_DAY, KDIL_NIGHT
      REAL(kind=dp)    :: BG_O3, BG_NO2, BG_NO, BG_CO, BG_VOC, BG_ISO
      REAL(kind=dp)    :: TMIN_K, TMAX_K
      NAMELIST /RUNCFG/ TEND_S, DT_S, STATE_IN, &
           EMIS_NOX, EMIS_VOC, EMIS_ISO, EMIS_CO, &
           NOX_FACTOR, VOC_FACTOR, NOXW, VOCW, COW, WKEND, DUMMYF, CONSTJ, C_H2O, &
           KDIL_DAY, KDIL_NIGHT, &
           BG_O3, BG_NO2, BG_NO, BG_CO, BG_VOC, BG_ISO, &
           TMIN_K, TMAX_K

!~~~> locals
      LOGICAL :: lex, lstate
      INTEGER :: iu, i, IERR, ISTAT(20), NF, nfail
      INTEGER :: IW
      CHARACTER(len=64) :: SNAME
      REAL(kind=dp) :: SV, T, T0, HH, SDAY, PPB2N
      REAL(kind=dp) :: ENOX, EVOCE, EISOE, tmpd
      REAL(kind=dp) :: JNOON(17), DVAL(NSPEC), RSTATE(20)
      INTEGER :: ICNTRL(20)
      REAL(kind=dp) :: RCNTRL(20)
      REAL(kind=dp) :: GA

!~~~> defaults (urban summer weekday)
      TEND_S = 259200.0d0
      DT_S   = 300.0d0
      STATE_IN = 'none'
      EMIS_NOX = 16.0d0
      EMIS_VOC = 12.0d0
      EMIS_ISO = 1.0d0
      EMIS_CO  = 40.0d0
      NOX_FACTOR = 1.0d0
      VOC_FACTOR = 1.0d0
      NOXW = 1.0d0
      VOCW = 1.0d0
      COW  = 1.0d0
      WKEND  = .FALSE.
      DUMMYF = .TRUE.
      CONSTJ = .FALSE.
      C_H2O  = 5.0d+18
      KDIL_DAY   = 1500.0d0
      KDIL_NIGHT = 24.0d0*3600.0d0
      BG_O3  = 30.0d0
      BG_NO2 = 5.0d0
      BG_NO  = 1.0d0
      BG_CO  = 120.0d0
      BG_VOC = 20.0d0
      BG_ISO = 0.2d0
      TMIN_K = 288.0d0
      TMAX_K = 301.0d0

!~~~> read namelist if present
      INQUIRE(FILE='cb05.nml', EXIST=lex)
      IF (lex) THEN
         OPEN(NEWUNIT=iu, FILE='cb05.nml', STATUS='OLD', ACTION='READ')
         READ(iu, RUNCFG)
         CLOSE(iu)
      END IF

!~~~> Initialization (sets C_M, j = noon values, C defaults)
      CALL Initialize()
      JNOON = j                      ! keep noon photolysis values

!~~~> optional state import (order = dump order, C(1..NSPEC))
      IF (TRIM(STATE_IN) /= 'none') THEN
         INQUIRE(FILE=TRIM(STATE_IN), EXIST=lstate)
         IF (lstate) THEN
            OPEN(NEWUNIT=iu, FILE=TRIM(STATE_IN), STATUS='OLD')
            READ(iu,*) NF
            DO i=1,NF
               READ(iu,*) SNAME, SV
               C(i) = SV
            END DO
            CLOSE(iu)
         END IF
      END IF

!~~~> overwrite time config from namelist
      TSTART = 0.0d0
      TEND   = TEND_S
      DT     = DT_S
      STEPMIN = 0.0d0
      STEPMAX = 0.0d0

      DO i=1,NVAR
        RTOL(i) = 1.0d-4
        ATOL(i) = 1.0d+0
      END DO

      ICNTRL  = 0
      RCNTRL  = 0.d0
      nfail   = 0

!~~~> fixed species (they live at the tail of C; INTEGRATE maps
!~~~> FIX => C(NVAR+1:NSPEC), and the driver has no FIX pointer of its own)
      C(ind_H2O)  = C_H2O
      C(ind_M)    = C_M
      C(ind_FSRC) = 1.0d0            ! forcing pools: 1 molec/cm3 each
      C(ind_FBCK) = 1.0d0
      C(ind_FSNK) = 1.0d0

!~~~> forcing globals: zero them once so a closed run is well defined
      E_NO = 0.0d0
      E_NO2 = 0.0d0
      E_CO = 0.0d0
      E_PAR = 0.0d0
      E_TOL = 0.0d0
      E_XYL = 0.0d0
      E_PNA = 0.0d0
      E_ISOP = 0.0d0
      CBG_O3 = 0.0d0
      CBG_NO2 = 0.0d0
      CBG_NO = 0.0d0
      CBG_CO = 0.0d0
      CBG_PAR = 0.0d0
      CBG_TOL = 0.0d0
      CBG_XYL = 0.0d0
      CBG_ISOP = 0.0d0
      K_DIL = 0.0d0

!~~~> open outputs
      OPEN(UNIT=20, FILE='traj.dat', STATUS='REPLACE')
      WRITE(20,'(A)') '# t_h  o3_ppb  nox_ppb  temp_K  Sday  kdil_per_h  emis_nox_ppbh'

!~~~> time loop
      T = TSTART
kron: DO WHILE (T < TEND)

        T0   = T
        HH   = MOD(T, 86400.0d0)/3600.0d0
        IW   = IWIN(HH)             ! 6-hour intervention window 1..4

!~~~~> analytic day cycle (sunrise ~6h, sunset ~18h)
        IF (CONSTJ) THEN
           SDAY = 1.0d0
           TEMP = 298.0d0
           C_M  = 2.45d+19
           j    = JNOON
        ELSE
           SDAY = SIN(3.14159265358979d0*(HH-6.0d0)/12.0d0)
           IF (SDAY < 0.0d0) SDAY = 0.0d0
           SDAY = SDAY**1.7d0
           TEMP = TMIN_K + (TMAX_K-TMIN_K)*SQRT(SDAY)
           C_M  = 2.45d+19*(298.0d0/TEMP)
           j    = JNOON*SDAY
        END IF
        PPB2N = C_M/1.0d9

!~~~~> emission time profile
        IF (.NOT. WKEND) THEN
           GA = GPOSS(8.0d0,1.6d0,HH) + GPOSS(18.5d0,2.0d0,HH)
           tmpd = 0.45d0 + 1.30d0*GA
        ELSE
           GA = GPOSS(11.0d0,3.0d0,HH) + GPOSS(17.0d0,3.0d0,HH)
           tmpd = 0.55d0 + 0.70d0*GA
        END IF
        ENOX  = EMIS_NOX*NOX_FACTOR*NOXW(IW)*tmpd
        EVOCE = EMIS_VOC*VOC_FACTOR*VOCW(IW)*tmpd
        EISOE = EMIS_ISO

!~~~~> forcing globals: emissions [molec/cm3/s], ventilation [1/s],
!~~~~> background inflow = K_DIL * C_bg [molec/cm3/s]
        IF (DUMMYF) THEN
           IF (KDIL_DAY > 0.0d0 .AND. SDAY > 0.0d0) THEN
              K_DIL = 1.0d0/(KDIL_DAY + (KDIL_NIGHT-KDIL_DAY)*(1.0d0-SDAY))
           ELSE
              K_DIL = 1.0d0/KDIL_NIGHT
           END IF

           E_NO   = 0.75d0*ENOX*PPB2N/3600.0d0
           E_NO2  = 0.25d0*ENOX*PPB2N/3600.0d0
           E_CO   = EMIS_CO*COW(IW)*PPB2N/3600.0d0
           E_PAR  = 0.50d0*EVOCE*PPB2N/3600.0d0
           E_TOL  = 0.25d0*EVOCE*PPB2N/3600.0d0
           E_XYL  = 0.20d0*EVOCE*PPB2N/3600.0d0
           E_PNA  = 0.05d0*EVOCE*PPB2N/3600.0d0
           E_ISOP = EISOE*PPB2N/3600.0d0*(SDAY**1.5d0)

           CBG_O3   = BG_O3*PPB2N
           CBG_NO2  = BG_NO2*PPB2N
           CBG_NO   = BG_NO*PPB2N
           CBG_CO   = BG_CO*PPB2N
           CBG_PAR  = 0.50d0*BG_VOC*PPB2N
           CBG_TOL  = 0.25d0*BG_VOC*PPB2N
           CBG_XYL  = 0.25d0*BG_VOC*PPB2N
           CBG_ISOP = BG_ISO*PPB2N
        ELSE
           K_DIL = 0.0d0
           E_NO = 0.0d0
           E_NO2 = 0.0d0
           E_CO = 0.0d0
           E_PAR = 0.0d0
           E_TOL = 0.0d0
           E_XYL = 0.0d0
           E_PNA = 0.0d0
           E_ISOP = 0.0d0
           CBG_O3 = 0.0d0
           CBG_NO2 = 0.0d0
           CBG_NO = 0.0d0
           CBG_CO = 0.0d0
           CBG_PAR = 0.0d0
           CBG_TOL = 0.0d0
           CBG_XYL = 0.0d0
           CBG_ISOP = 0.0d0
        END IF

!~~~~> keep the fixed-species block consistent with the forcing
        C(ind_H2O)  = C_H2O
        C(ind_M)    = C_M
        C(ind_FSRC) = 1.0d0
        C(ind_FBCK) = 1.0d0
        C(ind_FSNK) = 1.0d0

!~~~~> update rates (must happen AFTER the forcing globals are refreshed)
        TIME = T
        CALL Update_SUN()
        CALL Update_RCONST()

!~~~~> integrate one outer step: chemistry AND forcing together
        CALL INTEGRATE( TIN       = T,        &
                        TOUT      = T+DT,     &
                        ICNTRL_U  = ICNTRL,   &
                        RCNTRL_U  = RCNTRL,   &
                        ISTATUS_U = ISTAT,    &
                        RSTATUS_U = RSTATE,   &
                        IERR_U    = IERR      )
        IF ( IERR .LT. 0 ) THEN
           nfail = nfail + 1
           WRITE(6,*) 'WARN solver failure at T=', T, ' IERR=', IERR
        END IF
        T = RSTATE(1)

!~~~~> trajectory line
        WRITE(20,'(7ES13.5)') T/3600.0d0, &
              C(ind_o3)/PPB2N, &
              (C(ind_no)+C(ind_no2))/PPB2N, &
              TEMP, SDAY, K_DIL*3600.0d0, ENOX

      END DO kron

!~~~> final state dump (same format as the pure-chemistry driver)
      CLOSE(20)
      OPEN(UNIT=97, FILE='state_out.dat', STATUS='REPLACE')
      WRITE(97,'(I4)') NSPEC
      DO i=1,NSPEC
         WRITE(97,'(A,1X,E25.15)') TRIM(SPC_NAMES(i)), C(i)/CFACTOR
      END DO
      CLOSE(97)

      WRITE(6,*) 'NFAIL=', nfail
      WRITE(6,*) 'DONE'

CONTAINS

      INTEGER FUNCTION IWIN(HH)
!       6-hour intervention windows: [0,6) [6,12) [12,18) [18,24)
        REAL(kind=dp) :: HH
        IWIN = MIN(4, INT(HH/6.0d0) + 1)
        IF (IWIN < 1) IWIN = 1
      END FUNCTION IWIN

      REAL(kind=dp) FUNCTION GPOSS(MU,SD,HH)
        REAL(kind=dp) :: MU,SD,HH,X
        X = HH-MU
        IF (X >  12.0d0) X = X - 24.0d0
        IF (X < -12.0d0) X = X + 24.0d0
        GPOSS = EXP(-(X/SD)**2)
      END FUNCTION GPOSS

END PROGRAM cb05_Driver
