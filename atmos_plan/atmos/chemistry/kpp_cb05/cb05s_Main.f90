! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  cb05s_Main.f90 - SCENARIO driver for the CB05 standalone box (G4).
!
!  Adds, on top of the pure-chemistry driver:
!    * namelist cb05.nml    : run length, emissions, dilution, weekend flags
!    * diurnal forcing      : J(t), T(t), k_dil(t) from an analytic day cycle
!    * sources/sinks        : operator splitting after each outer step
!                             C += dt*( E_i(t) - k(t)*(C_i - Cbg_i) )
!      (splitting is fine for G4's qualitative tests at dt=300 s;
!       the n=8 enumeration will move forcing into dummy equations)
!    * state input/output   : state_in.dat (optional) / state_out.dat
!    * traj.dat             : one line per outer step
!  Solver-failure probe uses KPP semantics: IERR < 0 is a failure
!  (positive IERR means success).
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
      REAL(kind=dp)    :: NOXW(4), VOCW(4)
      LOGICAL          :: WKEND
      REAL(kind=dp)    :: KDIL_DAY, KDIL_NIGHT
      REAL(kind=dp)    :: BG_O3, BG_NO2, BG_NO, BG_CO, BG_VOC, BG_ISO
      REAL(kind=dp)    :: TMIN_K, TMAX_K
      NAMELIST /RUNCFG/ TEND_S, DT_S, STATE_IN, &
           EMIS_NOX, EMIS_VOC, EMIS_ISO, EMIS_CO, &
           NOX_FACTOR, VOC_FACTOR, NOXW, VOCW, WKEND, &
           KDIL_DAY, KDIL_NIGHT, &
           BG_O3, BG_NO2, BG_NO, BG_CO, BG_VOC, BG_ISO, &
           TMIN_K, TMAX_K

!~~~> locals
      LOGICAL :: lex, lstate
      INTEGER :: iu, i, IERR, ISTAT(20), NF, nfail
      INTEGER :: IW
      CHARACTER(len=64) :: SNAME
      REAL(kind=dp) :: SV, T, T0, DTT, HH, SDAY, KDIL
      REAL(kind=dp) :: ENOX, EVOCE, EISOE, PPB2N, tmpd
      REAL(kind=dp) :: JNOON(17), E(NVAR), CBG(NVAR)
      REAL(kind=dp) :: DVAL(NSPEC), RSTATE(20)
      INTEGER :: ICNTRL(20)
      REAL(kind=dp) :: RCNTRL(20)
      REAL(kind=dp) :: GA, GB

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
      WKEND  = .FALSE.
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

!~~~> background vector (molecules/cm3), zero elsewhere
      CBG = 0.0d0
      PPB2N = C_M/1.0d9
      CBG(ind_o3)   = BG_O3*PPB2N
      CBG(ind_no2)  = BG_NO2*PPB2N
      CBG(ind_no)   = BG_NO*PPB2N
      CBG(ind_co)   = BG_CO*PPB2N
      CBG(ind_par)  = 0.50d0*BG_VOC*PPB2N
      CBG(ind_tol)  = 0.25d0*BG_VOC*PPB2N
      CBG(ind_xyl)  = 0.25d0*BG_VOC*PPB2N
      CBG(ind_isop) = BG_ISO*PPB2N

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
        SDAY = SIN(3.14159265358979d0*(HH-6.0d0)/12.0d0)
        IF (SDAY < 0.0d0) SDAY = 0.0d0
        SDAY = SDAY**1.7d0

!~~~~> temperature and air density
        TEMP = TMIN_K + (TMAX_K-TMIN_K)*SQRT(SDAY)
        C_M  = 2.45d+19*(298.0d0/TEMP)
        PPB2N = C_M/1.0d9

!~~~~> photolysis frequencies follow sunlight
        j = JNOON*SDAY

!~~~~> dilution rate (day/night boundary layer)
        IF (KDIL_DAY > 0.0d0 .AND. SDAY > 0.0d0) THEN
           KDIL = 1.0d0/(KDIL_DAY + (KDIL_NIGHT-KDIL_DAY)*(1.0d0-SDAY))
        ELSE
           KDIL = 1.0d0/KDIL_NIGHT
        END IF

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

!~~~~> source vector (molecules/cm3/s)
        E = 0.0d0
        E(ind_no)   = 0.75d0*ENOX*PPB2N/3600.0d0
        E(ind_no2)  = 0.25d0*ENOX*PPB2N/3600.0d0
        E(ind_co)   = EMIS_CO*PPB2N/3600.0d0
        E(ind_par)  = 0.50d0*EVOCE*PPB2N/3600.0d0
        E(ind_tol)  = 0.25d0*EVOCE*PPB2N/3600.0d0
        E(ind_xyl)  = 0.20d0*EVOCE*PPB2N/3600.0d0
        E(ind_pna)  = 0.05d0*EVOCE*PPB2N/3600.0d0
!        biogenic isoprene: light & temperature driven
        E(ind_isop) = EISOE*PPB2N/3600.0d0*(SDAY**1.5d0)

!~~~~> update rates, integrate one outer step
        TIME = T
        CALL Update_SUN()
        CALL Update_RCONST()
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

!~~~~> operator splitting: emissions in, dilution toward background out
        DTT = T - T0
        DO i=1,NVAR
           C(i) = C(i) + DTT*( E(i) - KDIL*( C(i) - CBG(i) ) )
           IF (C(i) < 0.0d0) C(i) = 0.0d0
        END DO
!       fixed species live at the tail of C (INTEGRATE maps FIX=>C(NVAR+1:NSPEC));
!       the driver-level FIX pointer is not associated here -> write C directly
        C(ind_H2O) = 5.0d+18        ! water vapor ~2% of air (OH source via O1D+H2O)
        C(ind_M)   = C_M

!~~~~> trajectory line
        WRITE(20,'(7ES13.5)') T/3600.0d0, &
              C(ind_o3)/PPB2N, &
              (C(ind_no)+C(ind_no2))/PPB2N, &
              TEMP, SDAY, KDIL*3600.0d0, ENOX

      END DO kron

!~~~> final state dump (same format as pure driver)
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
