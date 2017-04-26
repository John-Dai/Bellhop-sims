MODULE sspmod

  ! holds SSP input by user and associated variables

  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: MaxSSP = 2001, MaxMedia = 501
  REAL      (KIND=8) :: alphaR = 1500, betaR = 0, alphaI = 0, betaI = 0, rhoR = 1
  INTEGER               iz, ILoc, Lay

  TYPE SSPStructure
     INTEGER           :: Loc( MaxMedia ), NPts( MaxMedia ), NMedia
     REAL     (KIND=8) :: z( MaxSSP ), rho( MaxSSP )
     REAL     (KIND=8) :: Depth( MaxMedia ), sigma( MaxMedia )
     COMPLEX  (KIND=8) :: cp( MaxSSP ), cs( MaxSSP ), n2( MaxSSP ),  &
                          cpSpline( 4, MaxSSP ), csSpline( 4, MaxSSP ), rhoSpline( 4, MaxSSP )
     CHARACTER (LEN=1) :: Type
     CHARACTER (LEN=2) :: AttenUnit
     CHARACTER (LEN=8) :: Material( MaxMedia )
  END TYPE SSPStructure

  TYPE( SSPStructure ) :: SSP

CONTAINS

  SUBROUTINE EvaluateSSP( cP, cS, rho, Medium, N1, Freq, Task, ENVFile, PRTFile )

    ! Call the particular SSP routine specified by SSPType
    ! Performs two Tasks:
    !    Task = 'TAB'  then tabulate cP, cS, rho
    !    Task = 'INIT' then initialize
    ! Note that Freq is only needed if Task = 'INIT'

    INTEGER,           INTENT(IN)    :: ENVFile, PRTFile, Medium
    INTEGER,           INTENT(INOUT) :: N1
    REAL     (KIND=8), INTENT(OUT)   :: rho( * )
    REAL     (KIND=8), INTENT(IN)    :: Freq
    COMPLEX  (KIND=8), INTENT(OUT)   :: cP( * ), cS( * )
    CHARACTER (LEN=8), INTENT(IN)    :: Task
    REAL     (KIND=8)                :: h, z
    COMPLEX  (KIND=8)                :: cPT, cST

    SELECT CASE ( SSP%Type )
    CASE ( 'A' )  !  Analytic profile option 
       IF ( Task( 1 : 4 ) == 'INIT' ) THEN
          N1 = 21
          CALL ANALYT( cP, cS, rho, Medium, N1, Freq, Task )
          h = ( SSP%Depth( Medium+1 ) - SSP%Depth( Medium ) ) / ( N1 - 1 )

          DO iz = 1, N1
             z   = SSP%Depth( Medium ) + ( iz - 1 ) * h
             cPT =  cP( iz )
             cST =  cS( iz )
             WRITE( PRTFile, FMT="( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )" ) &
                  z,  REAL( cPT ),  REAL( cST ), rho( iz ), AIMAG( cPT ), AIMAG( cST )
          END DO
       ELSE
          CALL ANALYT( cP, cS, rho, Medium, N1, Freq, Task )
       ENDIF
    CASE ( 'N' )  !  N2-linear profile option
       CALL N2Linear( cP, cS, rho, Medium, N1, Freq, Task, ENVFile, PRTFile )
    CASE ( 'C' )  !  C-linear profile option 
       CALL CLinear(  cP, cS, rho, Medium, N1, Freq, Task, ENVFile, PRTFile )
    CASE ( 'S' )  !  Cubic spline profile option 
       CALL CCubic(   cP, cS, rho, Medium, N1, Freq, Task, ENVFile, PRTFile )
    CASE DEFAULT  !  Non-existent profile option 
       WRITE( PRTFile, * ) 'Profile option: ', SSP%Type
       CALL ERROUT( PRTFile, 'F', 'EvaluateSSP', 'Unknown profile option' )
    END SELECT

    RETURN
  END SUBROUTINE EvaluateSSP

  !**********************************************************************!

  SUBROUTINE N2Linear( cP, cS, rho, Medium, N1, Freq, Task, ENVFile, PRTFile )

    ! Tabulate cP, cS, rho for specified Medium
    ! Uses N2-linear segments for P and S-wave speeds
    ! Uses rho-linear segments for density

    INTEGER,           INTENT(IN)  :: ENVFile, PRTFile
    INTEGER,           INTENT(IN)  :: Medium
    REAL     (KIND=8), INTENT(IN)  :: Freq
    REAL     (KIND=8), INTENT(OUT) :: rho( * )
    COMPLEX  (KIND=8), INTENT(OUT) :: cP( * ), cS( * )
    CHARACTER (LEN=8), INTENT(IN)  :: Task
    INTEGER              N, N1
    REAL     (KIND=8) :: h, R, zT
    COMPLEX  (KIND=8) :: N2Bot, N2Top

    ! If Task = 'INIT' then this is the first call and SSP is read.
    ! Any other call is a request for SSP subtabulation.

    IF ( Task( 1 : 4 ) == 'INIT' ) THEN   ! Task 'INIT' for initialization
       CALL ReadSSP( Medium, N1, Freq, ENVFile, PRTFile )
    ELSE   ! Task = 'TABULATE'
       ILoc = SSP%Loc( Medium )
       N    = N1 - 1
       h    = ( SSP%z( ILoc + SSP%NPts( Medium ) ) - SSP%z( ILoc + 1 ) ) / N
       Lay  = 1

       DO iz = 1, N1
          zT = SSP%z( ILoc + 1 ) + ( iz - 1 ) * h
          IF ( iz == N1 ) zT = SSP%z( ILoc + SSP%NPts( Medium ) )   ! Make sure no overshoot

          DO WHILE ( zT > SSP%z( ILoc + Lay + 1 ) )
             Lay = Lay + 1
          END DO

          R = ( zT - SSP%z( ILoc + Lay ) ) / ( SSP%z( ILoc + Lay+1 ) - SSP%z( ILoc + Lay ) )

          ! P-wave
          N2Top    = 1.0 / SSP%cp( ILoc + Lay     )**2
          N2Bot    = 1.0 / SSP%cp( ILoc + Lay + 1 )**2
          cP( iz ) = 1.0 / SQRT( ( 1.0 - R ) * N2Top + R * N2Bot )

          ! S-wave
          IF ( SSP%cs(ILoc + Lay) /= 0.0 ) THEN
             N2Top    = 1.0 / SSP%cs( ILoc + Lay     )**2
             N2Bot    = 1.0 / SSP%cs( ILoc + Lay + 1 )**2
             cS( iz ) = 1.0 / SQRT( ( 1.0 - R ) * N2Top + R * N2Bot )
          ELSE
             cS( iz ) = 0.0
          ENDIF

          rho( iz ) = ( 1.0 - R ) * SSP%rho( ILoc + Lay ) + R * SSP%rho( ILoc + Lay + 1 )
       END DO

    ENDIF

    RETURN
  END SUBROUTINE N2Linear
  
  !**********************************************************************!
  
  SUBROUTINE CLinear( cP, cS, rho, Medium, N1, Freq, Task, ENVFile, PRTFile  )

    ! Tabulate cP, cS, rho for specified Medium

    ! Uses c-linear segments for P and S-wave speeds
    ! Uses rho-linear segments for density
    INTEGER,           INTENT(IN)  :: ENVFile, PRTFile
    INTEGER,           INTENT(IN)  :: Medium
    REAL     (KIND=8), INTENT(IN)  :: Freq
    REAL     (KIND=8), INTENT(OUT) :: rho( * )
    COMPLEX  (KIND=8), INTENT(OUT) :: cP( * ), cS( * )
    CHARACTER (LEN=8), INTENT(IN)  :: Task
    INTEGER           :: N, N1
    REAL     (KIND=8) :: h, R, zT

    ! If Task = 'INIT' then this is the first call and SSP is read.
    ! Any other call is a request for SSP subtabulation.

    IF ( Task( 1 : 4 ) == 'INIT' ) THEN   ! Task 'INIT' for initialization
       CALL ReadSSP( Medium, N1, Freq, ENVFile, PRTFile )
    ELSE   ! Task = 'TABULATE'
       ILoc = SSP%Loc( Medium )
       N    = N1 - 1
       h    = ( SSP%z( ILoc + SSP%NPts( Medium ) ) - SSP%z( ILoc + 1 ) ) / N
       Lay  = 1

       DO iz = 1, N1
          zT = SSP%z( ILoc + 1 ) + ( iz - 1 ) * h
          IF ( iz == N1 ) zT = SSP%z( ILoc + SSP%NPts( Medium ) )   ! Make sure no overshoot

          DO WHILE ( zT > SSP%z( ILoc + Lay + 1 ) )
             Lay = Lay + 1
          END DO

          R = ( zT - SSP%z( ILoc + Lay ) ) / ( SSP%z( ILoc + Lay + 1 ) - SSP%z( ILoc + Lay ) )
          cP(  iz ) = ( 1.0 - R ) * SSP%cp(  ILoc + Lay ) + R * SSP%cp(  ILoc + Lay+1 )
          cS(  iz ) = ( 1.0 - R ) * SSP%cs(  ILoc + Lay ) + R * SSP%cs(  ILoc + Lay+1 )
          rho( iz ) = ( 1.0 - R ) * SSP%rho( ILoc + Lay ) + R * SSP%rho( ILoc + Lay+1 )
       END DO
    ENDIF

    RETURN
  END SUBROUTINE CLinear
  
  !**********************************************************************!
  
  SUBROUTINE CCubic( cP, cS, rho, Medium, N1, Freq, Task, ENVFile, PRTFile  )

    ! Tabulate cP, cS, rho for specified Medium
    ! using cubic spline interpolation

    INTEGER,           INTENT(IN)  :: ENVFile, PRTFile
    INTEGER,           INTENT(IN)  :: Medium
    REAL     (KIND=8), INTENT(IN)  :: Freq
    REAL     (KIND=8), INTENT(OUT) :: rho( * )
    COMPLEX  (KIND=8), INTENT(OUT) :: cP( * ), cS( * )
    CHARACTER (LEN=8), INTENT(IN)  :: Task
    INTEGER              N, N1, IBCBeg, IBCEnd
    REAL     (KIND=8) :: h, zT, HSPLNE
    COMPLEX  (KIND=8) :: SPLINE

    ! If Task = 'INIT' then this is the first call and SSP is read.
    ! Any other call is a request for SSP subtabulation.

    IF ( Task( 1 : 4 ) == 'INIT' ) THEN   ! Task 'INIT' for initialization
       CALL ReadSSP( Medium, N1, Freq, ENVFile, PRTFile )

       ILoc = SSP%Loc( Medium )

       DO iz = 1, SSP%Npts( Medium )
          SSP%cpSpline(  1, ILoc + iz ) = SSP%cp( ILoc + iz )
          SSP%csSpline(  1, ILoc + iz ) = SSP%cs( ILoc + iz )
          SSP%rhoSpline( 1, ILoc + iz ) = SSP%rho( ILoc + iz )
       END DO

       ! Compute spline coefs 
       IBCBEG = 0
       IBCEND = 0
       CALL CSPLINE( SSP%z( ILoc + 1 ), SSP%cpSpline(  1, ILoc + 1 ), SSP%NPts( Medium ), IBCBEG, IBCEND, SSP%NPts( Medium ) )
       CALL CSPLINE( SSP%z( ILoc + 1 ), SSP%csSpline(  1, ILoc + 1 ), SSP%NPts( Medium ), IBCBEG, IBCEND, SSP%NPts( Medium ) )
       CALL CSPLINE( SSP%z( ILoc + 1 ), SSP%rhoSpline( 1, ILoc + 1 ), SSP%NPts( Medium ), IBCBEG, IBCEND, SSP%NPts( Medium ) )
    ELSE   ! Task = 'TABULATE'
       ILoc = SSP%Loc( Medium )
       N    = N1 - 1
       h    = ( SSP%z( ILoc + SSP%NPts( Medium ) ) - SSP%z( ILoc + 1 ) ) / N
       Lay  = 1

       DO iz = 1, N1
          zT = SSP%z( ILoc + 1 ) + ( iz - 1 ) * h
          IF ( iz == N1 ) zT = SSP%z( ILoc + SSP%NPts( Medium ) )   ! Make sure no overshoot
          DO WHILE ( zT > SSP%z( ILoc + Lay + 1 ) )
             Lay = Lay + 1
          END DO

          HSPLNE = zT - SSP%z( ILoc + Lay )

          cP(  iz ) =       SPLINE( SSP%cpSpline(  1, ILoc + Lay ), HSPLNE )
          cS(  iz ) =       SPLINE( SSP%csSpline(  1, ILoc + Lay ), HSPLNE )
          rho( iz ) = DBLE( SPLINE( SSP%rhoSpline( 1, ILoc + Lay ), HSPLNE ) )

       END DO
    ENDIF

    RETURN
  END SUBROUTINE CCubic

!**********************************************************************!

  SUBROUTINE ReadSSP( Medium, N1, Freq, ENVFile, PRTFile )

    ! reads the SSP data from the environmental file
    
    INTEGER, INTENT(IN) :: ENVFile, PRTFile
    INTEGER, INTENT(IN) :: Medium
    REAL     (KIND=8), INTENT(IN)  :: Freq
    INTEGER             :: N1
    COMPLEX  (KIND=8)   :: CRCI

    SSP%NPts( Medium ) = N1

    ! The variable SSP%Loc( Medium ) points to the starting point for the
    ! data in the arrays z, alpha, beta and rho
    IF ( Medium == 1 ) THEN
       SSP%Loc( Medium ) = 0
    ELSE
       SSP%Loc( Medium ) = SSP%Loc( Medium - 1 ) + SSP%NPts( Medium - 1 )
    ENDIF
    ILoc = SSP%Loc( Medium )

    !  Read in data and convert attenuation to Nepers/m 
    N1 = 1
    DO iz = 1, MaxSSP

       READ(  ENVFile, *    ) SSP%z( ILoc + iz ), alphaR, betaR, rhoR, alphaI, betaI
       WRITE( PRTFile, FMT="( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )" ) SSP%z( ILoc + iz ), alphaR, betaR, rhoR, alphaI, betaI

       SSP%cp(  ILoc + iz ) = CRCI( alphaR, alphaI, Freq, SSP%AttenUnit )
       SSP%cs(  ILoc + iz ) = CRCI( betaR,  betaI,  Freq, SSP%AttenUnit )
       SSP%rho( ILoc + iz ) = rhoR

       ! Did we read the last point?
       IF ( ABS( SSP%z( ILoc + iz ) - SSP%Depth( Medium + 1 ) ) < 100. * EPSILON( 1.0e0 ) ) THEN
          SSP%NPts( Medium ) = N1
          IF ( Medium == 1 ) SSP%Depth( 1 ) = SSP%z( 1 )
          IF ( SSP%NPts( Medium ) == 1 ) THEN
              WRITE( PRTFile, * ) '#SSP points: ', SSP%NPts( Medium )
              CALL ERROUT( PRTFile, 'F', 'ReadSSP', 'The SSP must have at least 2 points in each layer' )
          END IF

          RETURN
       ENDIF

       N1 = N1 + 1
    END DO

    ! Fall through means too many points in the profile
    WRITE( PRTFile, * ) 'Max. #SSP points: ', MaxSSP
    CALL ERROUT( PRTFile, 'F', 'ReadSSP', 'Number of SSP points exceeds limit' )

  END SUBROUTINE ReadSSP

END MODULE sspmod
