! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! MAIN - Custom minimal driver for the CB05 standalone box.
!      - constant J / T / C_M (no forcing, closed box)
!      - dumps ALL species at t0 and tEND for the G3 conservation check
!      - tolerance lines kept in the KPP "RTOL(i) = .../ATOL(i) = ..." form
!        so ladder-style sed edits still work
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM cb05_Driver

!~~~> Global variables and routines
      USE cb05_Model
      USE cb05_Initialize, ONLY: Initialize

!~~~> Local variables
      REAL(kind=dp) :: T, DVAL(NSPEC)
      REAL(kind=dp) :: RSTATE(20)
      INTEGER :: i
      INTEGER :: IERR, ISTAT(20)

!~~~> Control (in) arguments for the integration
      INTEGER :: ICNTRL(20)
      REAL(kind=dp) :: RCNTRL(20)

      ICNTRL  = 0
      RCNTRL  = 0.d0

!~~~> Initialization
      STEPMIN = 0.0d0
      STEPMAX = 0.0d0

      DO i=1,NVAR
        RTOL(i) = 1.0d-4
        ATOL(i) = 1.0d+0
      END DO

      CALL Initialize()

!~~~> Open log file for write
      CALL InitSaveData()

!~~~> Dump initial state (all species, molecules/cm3)
      OPEN(UNIT=98, FILE='cb05_t0.dat', STATUS='REPLACE')
      WRITE(98,'(I4)') NSPEC
      DO i=1,NSPEC
        WRITE(98,'(A,1X,E25.15)') TRIM(SPC_NAMES(i)), C(i)/CFACTOR
      END DO
      CLOSE(98)

!~~~> Time loop
      T = TSTART
kron: DO WHILE (T < TEND)

        TIME = T

!~~~> Write values of monitored species at each iteration
        CALL GetMass( C, DVAL )
        CALL SaveData()

!~~~> Update sunlight intensity and reaction rates
        CALL Update_SUN()
        CALL Update_RCONST()

!~~~> Call the integrator
        CALL INTEGRATE( TIN       = T,        &
                        TOUT      = T+DT,     &
                        ICNTRL_U  = ICNTRL,   &
                        RCNTRL_U  = RCNTRL,   &
                        ISTATUS_U = ISTAT,    &
                        RSTATUS_U = RSTATE,   &
                        IERR_U    = IERR      )
        IF ( IERR .LT. 0 ) THEN
          WRITE(6,*) 'WARN solver failure at T=', T, ' IERR=', IERR
        END IF
        T = RSTATE(1)

      END DO kron
!~~~> End Time loop

!~~~> Write final values and close log file
      CALL GetMass( C, DVAL )
      TIME = T
      CALL SaveData()
      CALL CloseSaveData()

!~~~> Dump final state (all species)
      OPEN(UNIT=97, FILE='cb05_tend.dat', STATUS='REPLACE')
      WRITE(97,'(I4)') NSPEC
      DO i=1,NSPEC
        WRITE(97,'(A,1X,E25.15)') TRIM(SPC_NAMES(i)), C(i)/CFACTOR
      END DO
      CLOSE(97)

      WRITE(6,*) 'DONE'

END PROGRAM cb05_Driver
