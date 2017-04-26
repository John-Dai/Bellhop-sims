PROGRAM KRAKENC

  ! Ocean acoustic normal modes.

  ! Copyright (C) 2009 Michael B. Porter

  ! This program is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! This program is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ! GNU General Public License for more details.

  ! You should have received a copy of the GNU General Public License
  ! along with this program.  If not, see <http://www.gnu.org/licenses/>.

  ! Originally developed as part of the author's dissertation under the supervision
  ! of Prof. Edward L. Reiss, Northwestern University

  USE KrakencMod
  USE SdRdRMod
  USE RefCoMod
  USE sspMod
  IMPLICIT NONE
  INTEGER            :: Min_Loc( 1 ), IFirst, ILast, IRec
  REAL               :: zMin, zMax
  REAL      (KIND=8) :: omega, Error
  CHARACTER (LEN=80) :: FileRoot

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  Profile: DO iProf = 1, 9999   ! Loop over a sequence of profiles

     NV( 1 : 5 ) = [ 1, 2, 4, 8, 16 ]

     ! Read in environmental info
     Title = 'KRAKENC-'
     CALL READIN( FileRoot, Title, Freq, MaxMedium, TopOpt, HSTop, NG, BotOpt, HSBot, ENVFile, PRTFile )

     READ(  ENVFile, * ) cLow, cHigh   ! Spectral limits
     WRITE( PRTFile, "( /, ' cLow = ', G12.5, 'm/s      cHigh = ', G12.5, 'm/s' )" ) cLow, cHigh

     READ(  ENVFile, * ) RMax          ! Maximum range for calculations
     WRITE( PRTFile, * ) 'RMax = ', RMax

     ! Read source/receiver depths
     zMin = SNGL( SSP%Depth( 1 ) )
     zMax = SNGL( SSP%Depth( SSP%NMedia + 1 ) )
     CALL ReadSdRd( ENVFile, PRTFile, zMin, zMax )

     omega2 = ( 2.0D0 * pi * Freq ) ** 2
     M      = MaxM

     ! Bot, Top reflection coefficients
     IF ( iProf == 1 ) CALL ReadReflectionCoefficient( FileRoot, HSBot%BC, HSTop%BC, PRTFile )

     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) 'Mesh multiplier   CPU seconds'

     DO ISet = 1, NSets   ! Main loop: solve the problem for a sequence of meshes
        N( 1 : SSP%NMedia ) = NG( 1 : SSP%NMedia ) * NV( ISet )
        h( 1 : SSP%NMedia ) = ( SSP%Depth( 2 : SSP%NMedia + 1 ) - SSP%Depth( 1 : SSP%NMedia ) ) / N( 1 : SSP%NMedia )
        hV( ISet )      = h( 1 )
        CALL Solve( FileRoot, Error )

        IF ( Error * 1000.0 * RMax < 1.0 ) GOTO 3000
     END DO

     ! Fall through indicates failure to converge
     CALL ERROUT( PRTFile, 'W', 'KRAKENC', 'Too many meshes needed: check convergence' )

3000 omega   = SQRT( omega2 )   ! Solution complete: discard modes with phase velocity above cHigh
     Min_Loc = MinLOC( DBLE( Extrap( 1, 1 : M ) ), DBLE( Extrap( 1, 1 : M ) ) > omega2 / cHigh ** 2 )
     M       = Min_Loc( 1 )

     ! Write eigenvalues to PRTFile and MODFile
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) '   I          k             alpha          Phase Speed       Group Speed'

     ! k() contains scatter losses; Extrap contains extrapolated values of wavenumbers, so combine ...
     k( 1 : M ) = SQRT( Extrap( 1, 1 : M ) + k( 1 : M ) )

     DO mode = 1, M, MAX( 1, M / 200 )  ! print every mode unless there are an awful lot
        WRITE( PRTFile, '( I5, 4G18.10 )' ) mode, k( mode ), omega / DBLE( k( mode ) ), VG( mode )
     END DO

     ! Zero out positive imaginary part which would cause growth in range. Should be small.
     WHERE ( AIMAG( k ) > 0.0D0 ) k = REAL( k )

     WRITE( MODFile, REC = IRecProfile + 4 ) M, LRecordLength

     IFirst = 1
     DO IREC = 1, 1 + ( 2 * M - 1 ) / LRecordLength
        ILast  = MIN( M, IFirst + LRecordLength / 2 - 1 )
        WRITE( MODFile, REC = IRecProfile + 5 + M + IREC ) CMPLX( k( IFirst : ILast ) )
        IFirst = ILast + 1
     END DO

     ! set record pointer to beginning of next mode set
     IRecProfile = IRecProfile + 6 + M + 1 + ( 2 * M - 1 ) / LRecordLength

  END DO Profile
  CLOSE( ENVFile )
  CLOSE( MODFile )

END PROGRAM KRAKENC

!**********************************************************************

SUBROUTINE Initialize

  ! Initializes arrays defining difference equations

  USE KrakencMod
  USE sspMod
  IMPLICIT NONE
  LOGICAL           :: ElasticFlag = .FALSE.
  INTEGER           :: IAllocStat, ii, J, Medium, NPoints, N1
  REAL     (KIND=8) :: Two_h
  COMPLEX  (KIND=8) :: cP2, cS2
  COMPLEX  (KIND=8), ALLOCATABLE :: cP( : ), cS( : )
  CHARACTER (LEN=8) :: Task

  cMin          = HUGE( cMin )
  FirstAcoustic = 0
  Loc( 1 )      = 0

  ! Allocate storage for finite-difference coefficients

  NPoints = SUM( N( 1 : SSP%NMedia ) ) + SSP%NMedia

  IF ( ALLOCATED( B1 ) ) DEALLOCATE( B1, B2, B3, B4, rho )
  ALLOCATE ( B1( NPoints ), B2( NPoints ), B3( NPoints ), B4( NPoints ), rho( NPoints ), &
             cP( NPoints ), cS( NPoints ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) &
       CALL ERROUT( PRTFile, 'F', 'KRAKENC - Initialize', 'Insufficient memory to allocate B1, B2, B3, B4 vectors. Reduce mesh.' )

  Media: DO Medium = 1, SSP%NMedia   ! Loop over media

     IF ( Medium /= 1 ) Loc( Medium ) = Loc( Medium - 1 ) + N( Medium - 1 ) + 1
     N1  = N(   Medium ) + 1
     ii  = Loc( Medium ) + 1

     ! EvaluateSSP reads in the data for a medium
     Task = 'TAB'
     CALL EvaluateSSP( cP( ii ), cS( ii ), rho( ii ), Medium, N1, Freq, Task, ENVFile, PRTFile)

     ! Load diagonals of the finite-difference equations

     IF ( cS( ii ) == ( 0.0, 0.0 ) ) THEN  ! Case of an acoustic medium

        SSP%Material( Medium ) = 'ACOUSTIC'
        IF ( FirstAcoustic == 0 ) FirstAcoustic = Medium
        LastAcoustic = Medium

        cMin = MIN( cMin, MINVAL( DBLE( cP( ii : ii + N( Medium ) ) ) ) )

        B1( ii : ii + N( Medium ) ) = -2.0D0 + h( Medium ) ** 2 * omega2 / cP( ii : ii + N( Medium ) ) ** 2

     ELSE                                  ! Case of an elastic medium
        IF ( SSP%sigma( Medium ) /= 0.0D0 ) &
           CALL ERROUT( PRTFile, 'F', 'KRAKENC', 'Rough elastic interfaces are not allowed' )

        SSP%Material( Medium ) = 'ELASTIC'
        ElasticFlag        = .TRUE.
        Two_h              = 2.0D0 * h( Medium )

        DO J = ii, ii + N( Medium )
           cMin = MIN( DBLE( cS( J ) ), cMin )

           cP2 = cP( J ) ** 2
           cS2 = cS( J ) ** 2

           B1(  J ) = Two_h / ( rho( J ) * cS2 )
           B2(  J ) = Two_h / ( rho( J ) * cP2 )
           B3(  J ) = 4.0D0 * Two_h * rho( J ) * cS2 * ( cP2 - cS2 ) / cP2
           B4(  J ) = Two_h * ( cP2 - 2.0D0 * cS2 ) / cP2
           rho( J ) = Two_h * omega2 * rho( J )
        END DO

     END IF

  END DO Media

  ! (cLow, cHigh) = phase speed interval for the mode search
  ! user specified interval is reduced if it exceeds domain
  ! of possible mode phase speeds

  ! Bottom properties
  IF ( HSBot%BC == 'A' ) THEN
     IF ( HSBot%cS /= ( 0.0, 0.0 ) ) THEN    ! Elastic  bottom half-space
        ElasticFlag = .TRUE.
        cMin   = MIN( cMin, DBLE( HSBot%cS ) )
     ELSE                                    ! Acoustic bottom half-space
        cMin   = MIN( cMin, DBLE( HSBot%cP ) )
     END IF
  END IF

  ! Top properties
  IF ( HSTop%BC == 'A' ) THEN
     IF ( HSTop%cS /= ( 0.0, 0.0 ) ) THEN    ! Elastic  top half-space
        ElasticFlag = .TRUE.
        cMin   = MIN( cMin, DBLE( HSTop%cS ) )
     ELSE                                    ! Acoustic top half-space
        cMin   = MIN( cMin, DBLE( HSTop%cP ) )
     END IF
  END IF

  ! If elastic medium then reduce cMin for Scholte wave
  IF ( ElasticFlag ) cMin = 0.85 * cMin
  cLow = MAX( cLow, 0.99 * cMin )

END SUBROUTINE Initialize

!**********************************************************************

SUBROUTINE Solve( FileRoot, Error )

  ! Solves the eigenvalue problem at the current mesh and produces a new extrapolation

  USE KrakencMod
  IMPLICIT NONE
  CHARACTER (LEN=80), INTENT(  IN ) :: FileRoot
  REAL    (KIND=8),   INTENT( OUT ) :: Error
  INTEGER          :: Min_Loc( 1 ), J, Key
  INTEGER          :: NTotal, NTotal1   ! number of mesh points (where eigenvector is sampled)
  REAL    (KIND=8) :: x1, x2, TStart, TEnd
  COMPLEX (KIND=8) :: T1, T2, F1, F2, xTemp( MaxM )


  CALL CPU_TIME( Tstart )
  CALL Initialize           ! set up the finite-difference mesh
  CALL Solve2( FileRoot )   ! solve for the eigenvalues

  Extrap( ISet, 1 : M ) = EVMat( ISet, 1 : M )

  ! Remove any eigenvalues outside the spectral limits
  ! Typically the last 'eigenvalue' results from forcing a zero in funct( x ) for x outside the limits
  ! Inverse iteration would usually fail for that mode

  Min_Loc = MINLOC( REAL( Extrap( 1, 1 : M ) ), REAL( Extrap( 1, 1 : M ) ) > omega2 / cHigh ** 2 )
  M       = Min_Loc( 1 )

  NTotal  = SUM( N( FirstAcoustic : LastAcoustic ) )
  NTotal1 = NTotal + 1
  
  IF ( ISet == 1 ) CALL Vector( FileRoot, NTotal, NTotal1 )   ! If this is the first mesh, compute the eigenvectors

  xTemp( 1 : M ) = Extrap( Iset, 1 : M )
  CALL ORDER( xTemp, M )    ! order eigenvalues by real part
  EVMat(  ISet, 1 : M ) = xTemp( 1 : M )
  Extrap( ISet, 1 : M ) = xTemp( 1 : M )

  ! Richardson extrapolation to improve the accuracy

  Error = 1.0D10          ! initialize error to a large number
  KEY   = 2 * M / 3 + 1   ! index of element used to check convergence

  IF ( ISet > 1 ) THEN
     T1 = Extrap( 1, KEY )

     DO J = ISet - 1, 1, -1
        ModeLoop: DO mode = 1, M
           x1 = NV( J    ) ** 2
           x2 = NV( ISet ) ** 2
           F1 = Extrap( J,     mode )
           F2 = Extrap( J + 1, mode )
           Extrap( J, mode ) = F2 - ( F1 - F2 ) / ( x2 / x1 - 1.0D0 )
        END DO ModeLoop
     END DO

     T2    = Extrap( 1, KEY )
     Error = ABS( T2 - T1 )
  END IF

  CALL CPU_TIME( Tend )   ! check elapsed time
  ET( ISet ) = Tend - Tstart
  WRITE( PRTFile, '( 1X, I8, 6X, G15.3, ''s'' )' ) NV( ISet ), ET( ISet )

END SUBROUTINE Solve

!**********************************************************************

SUBROUTINE Solve2( FileRoot )

  ! Provides initial guess to root finder for each EVMat( I )

  ! A sequence of restarts is done to compute a complete set of modes for the first two sets (ISet = 1, 2).
  ! 'Complete' in this case means that all the modes in the user-specified interval are found.
  ! However, if, for ISet=2, the same number of modes is found as for ISet=1 then it is assumed
  ! that the set is complete and no restarts are performed

  USE KrakencMod
  USE RootFinderSecantMod
  IMPLICIT NONE
  INTEGER            :: ITry, Iteration, MaxIteration, J, MaxTries, NzTab, IFirst, ILast, IRec, ii
  REAL               :: RVAR1, RVAR2, kLow, kHigh
  REAL      (KIND=8) :: Tolerance, x1, x2
  COMPLEX            :: kIG( MaxM )
  COMPLEX   (KIND=8) :: x, P( 10 ), xTemp( MaxM ), cTry, kTry
  CHARACTER (LEN=80) :: ErrorMessage, FileRoot

  x            = ( 1.0D0, 0.0D0 ) * omega2 / cLow ** 2
  MaxIteration = 1000   ! maximum # of iterations in root-finder

  IF ( TopOpt( 5 : 5 ) == '.' .AND. ISet <= 2 ) THEN
     MaxTries = 5000   ! MaxTries = # of restarts for the root finder
  ELSE
     MaxTries = 1
  END IF

  WRITE( PrtFile, * ) 'Max. number of restarts in root finder, MaxTries = ', MaxTries

  ITry = 1   ! counter that keeps track of # of tries

  ! optionally read an existing mode file to get initial guesses
  IF ( ISet == 1 .AND. TopOpt( 4 : 4 ) == 'I' ) THEN
     OPEN ( FILE = 'MODFile', UNIT = MODFile, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED' )

     READ( MODFile, REC = 5 ) M, LRecordLength

     IFirst = 1
     DO IREC = 1, 1 + ( 2 * M - 1 ) / LRecordLength
        ILast  = MIN( M, IFirst + LRecordLength / 2 - 1 )
        READ( MODFile, REC = 6 + M + IREC ) ( kIG( mode ), mode = IFirst, ILast )
        IFirst = ILast + 1
     END DO
  END IF

  ! start looking for each root in succession
  ModeLoop: DO mode = 1, M
     ! For first or second meshes, use a high guess
     ! Otherwise use extrapolation to produce an initial guess

5500 CONTINUE
     x = 1.00001 * x

     IF ( ISet == 1 .AND. TopOpt( 4 : 4 ) == 'I' ) THEN
        x = kIG( mode ) ** 2   ! initial guess from MODFile
     ELSE IF ( ISet >= 2 ) THEN
        P( 1 : ISet - 1 ) = EVMat( 1 : ISet - 1, mode )

        IF ( ISet >= 3 ) THEN
           DO ii = 1, ISet - 2
              DO J = 1, ISet - ii - 1
                 x1 = hV( J      ) ** 2
                 x2 = hV( J + ii ) ** 2

                 P( J ) = ( ( hV( ISet ) ** 2 - x2 ) * P( J     ) - &
                            ( hV( ISet ) ** 2 - x1 ) * P( J + 1 ) ) &
                          / ( x1 - x2 )
              END DO
           END DO
           x = P( 1 )
        END IF

     END IF

     ! Use the secant method to refine the eigenvalue
     ! Tolerance below must be satisfied by 3 iterates in the root finding.
     ! Therefore you will normally get much better accuary on the returned value
     ! 5.0 below was necessary for PGI Fortran on double.env case.
     ! 4.0 was fine for all the other compilers

     Tolerance = ABS( x ) * 10.0D0 ** ( 5.0 - PRECISION( x ) )
     CALL RootFinderSecant( x, Tolerance, Iteration, MaxIteration, ErrorMessage )

     IF ( ErrorMessage /= ' ' ) THEN   ! any problems in root finder?
        CALL ERROUT( PRTFile, 'W', 'KRAKENC-RootFinderSecant', ErrorMessage )
        x = TINY( REAL( x ) )   ! make sure value discarded later
     END IF

     ! Is the mode inside the user specified spectral limits?

     IF ( SQRT( omega2 ) / cHigh > REAL( SQRT( x ) ) ) THEN
        ! Mode outside, restart at a random point
        CALL RANDOM_NUMBER( RVAR1 )
        CALL RANDOM_NUMBER( RVAR2 )

        !cTry  = cLow + RVAR1 * ( cHigh - cLow ) + ( 0.0, 0.01 ) * RVAR2 * cLow
        !x     = omega2 / cTry ** 2

        kLow  = REAL( sqrt( omega2 ) / cHigh )
        kHigh = REAL( sqrt( omega2 ) / cLow  )
        kTry  =  kLow + RVAR1 * ( kHigh - kLow ) - ( 0, 0.01D0 ) * RVAR2 * kLow
        x     = kTry ** 2
        cTry  = sqrt( omega2 ) / kTry

        MaxIteration = 30 ! but don't let root-finder search long

        IF ( ITry < MaxTries ) THEN
           ! MaxTries = MAX( MaxTries, 5 * mode )
           ITry     = ITry + 1
           ! WRITE( *, * ) 'Restart at phase speed = ', ITry, cTry
           GOTO 5500
        ELSE
           M = mode - 1
           EXIT   ! done searching for modes
        END IF
     END IF

     EVMat( ISet, mode ) = x
     xTemp( mode )       = x

  END DO ModeLoop

  ! If no modes, open a dummy MODFile for use by FIELD3D
  IF ( M == 0 ) THEN
     LRecordLength = 32
     NzTab         = 0

     ! Open MODFile and write header
     IF ( iProf == 1 ) THEN
        OPEN ( FILE = TRIM( FileRoot) //'.mod', UNIT = MODFile, ACCESS = 'DIRECT', RECL = 4 * LRecordLength, FORM = 'UNFORMATTED' )
     END IF

     WRITE( MODFile, REC = 1 ) LRecordLength, Title, REAL( Freq ), 1, NzTab, NzTab
     WRITE( MODFile, REC = 5 ) M, LRecordLength
     CALL ERROUT( PRTFile, 'F', 'KRAKENC', 'No modes for given phase speed interval' )
  END IF

  IF ( M == MaxM ) THEN   ! Have we hit the ceiling on max # modes
     WRITE( PRTFile, * ) 'Number of modes = ', M
     WRITE( PRTFile, * ) 'Number of modes allowable = ', MaxM
     CALL ERROUT( PRTFile, 'W', 'KRAKENC', 'Too many modes: Program will compute as many as it can' )
  END IF

  CALL ORDER( xTemp, M )    ! order eigenvalues by real part
  EVMat( ISet, 1 : M ) = xTemp( 1 : M )

END SUBROUTINE Solve2

!**********************************************************************

SUBROUTINE FUNCT( x, Delta, iPower )

  ! FUNCT( x ) = 0 is the dispersion relation

  USE KrakencMod
  IMPLICIT NONE
  INTEGER,           PARAMETER :: iPowerR = 50, iPowerF = -50
  REAL (KIND=8),     PARAMETER :: Roof = 1.0D+50, Floor = 1.0D-50
  INTEGER                      :: iPower, iPowerBot, J
  COMPLEX (KIND=8), INTENT(IN) :: x
  COMPLEX (KIND=8)             :: Delta, f, g, fBot, gBot

!!$  IF ( REAL( x ) <= omega2 / cHigh ** 2 ) THEN    ! For a k below the spectrum limit, force a zero
!!$     Delta  = 0.0D0
!!$     iPower = 0
!!$     RETURN
!!$  END IF

  CALL BCImpedance(    x, 'BOT', HSBot, f,    g,    iPower  )   ! Bottom impedance
  CALL AcousticLayers( x,               f,    g,    iPower  )   ! Shoot through acoustic layers
  CALL BCImpedance(    x, 'TOP', HSTop, fBot, gBot, iPowerBot ) ! Top impedance

  Delta  = f * gBot - g * fBot
  iPower = iPower + iPowerBot

  ! Delta = f /g - f1 / g1
  ! Delta = g / f
  ! Delta = 1 / ( log( f * g1 - g * f1 ) + log( 10.0d0 ) * ( iPower + iPower1 ) )
  ! iPower = 0

  ! Deflate previous roots

  IF ( mode > 1 ) THEN
     DO J = 1, mode - 1
        Delta = Delta / ( x - EVMat( ISet, J ) )

        ! Scale if necessary
        DO WHILE ( ABS( DBLE( Delta ) ) < Floor .AND. ABS( Delta ) > 0.0D0 )
           Delta  = Roof * Delta
           iPower = iPower - iPowerR
        END DO

        DO WHILE ( ABS( DBLE( Delta ) ) > Roof )
           Delta  = Floor * Delta
           iPower = iPower - iPowerF
        END DO

     END DO
  END IF

END SUBROUTINE FUNCT

!**********************************************************************

SUBROUTINE AcousticLayers( x, f, g, iPower )

  ! Shoot through acoustic layers

  USE KrakencMod
  USE sspMod
  IMPLICIT NONE
  INTEGER,           PARAMETER :: iPowerF = -50
  REAL    (KIND=8),  PARAMETER :: Roof = 1.0D+50, Floor = 1.0D-50
  INTEGER                      :: iPower, ii, Medium
  REAL    (KIND=8)             :: rhoMedium
  COMPLEX (KIND=8), INTENT(IN) :: x
  COMPLEX (KIND=8)             :: F, G, p0 = 0, p1, p2, h2k2
  ! for 2-layer code at end: COMPLEX (KIND=8) :: gamma1, gamma2

  IF ( FirstAcoustic == 0 ) RETURN

  Media: DO Medium = LastAcoustic, FirstAcoustic, -1   ! Loop over successive acoustic media

     h2k2      = h( Medium ) ** 2 * x
     ii        = Loc( Medium ) + N( Medium ) + 1
     rhoMedium = rho(  Loc( Medium ) + 1  )   ! density is made homogeneous using value at top of each medium

     p1 = -2.0D0 * g
     p2 = ( B1( ii ) - h2k2 ) * g - 2.0D0 * h( Medium ) * f * rhoMedium

     ! Shoot (towards surface) through a single medium
     DO ii = Loc( Medium ) + N( Medium ), Loc( Medium ) + 1, -1

        p0 = p1
        p1 = p2
        p2 = ( h2k2 - B1( ii ) ) * p1 - p0

        DO WHILE ( ABS(  DBLE( p2 ) ) > Roof )  ! Scale if necessary
           p0     = Floor * p0
           p1     = Floor * p1
           p2     = Floor * p2
           iPower = iPower - iPowerF
        END DO

     END DO

     ! F = P' / rho and G = -P since F P + G P' / rho = 0
     rhoMedium = rho( Loc( Medium ) + 1 )   ! density at top of layer
     f         = -( p2 - p0 ) / ( 2.0D0 * h( Medium ) ) / rhoMedium
     g         = -p1
  END DO Media

  ! here's a little two-layer analytic formula as a test
  !d1 =  5000
  !d2 = 10000
  !gamma1 = SQRT( ( 2 + b1(    1 ) ) / h( 1 ) ** 2 - x )
  !gamma2 = SQRT( ( 2 + b1( 6000 ) ) / h( 2 ) ** 2 - x )
  !rho2 = 1.8
  !f = 0
  !g = SIN( gamma1 * d1 ) * COS( gamma2 * (d1-d2) ) / gamma1 - COS( gamma1 * d1 ) * SIN( gamma2 * ( d1-d2 ) ) / gamma2
  !iPower = 0

END SUBROUTINE AcousticLayers

!**********************************************************************

SUBROUTINE Vector( FileRoot, NTotal, NTotal1 )

  ! Do inverse iteration to compute each of the eigenvectors and write these to the disk file

  USE KrakencMod
  USE SdRdRMod
  USE sspMod
  USE InverseIterationMod
  IMPLICIT NONE
  INTEGER, INTENT( IN ) :: NTotal, NTotal1   ! number of mesh points (where eigenvector is sampled)
  CHARACTER (LEN=80), INTENT( IN ) :: FileRoot
  INTEGER              :: IErr, iPower, ii, ITP, J, jj, L, Medium, NzTab
  INTEGER, ALLOCATABLE :: IzTab( : )
  REAL                 :: zTab( Pos%Nsd + Pos%Nrd )
  REAL                 :: z( NTotal1 )
  REAL,   ALLOCATABLE  :: WTS( : )
  REAL    (KIND=8)     :: h_rho = 0
  COMPLEX (KIND=8)     :: xh2, x, fTop, gTop, fBot, gBot
  COMPLEX (KIND=8)     :: Phi( NTotal1 ), d( NTotal1 ), e( NTotal1 + 1 )
  COMPLEX, ALLOCATABLE :: PhiTab( : )

  ! Tabulate z-coordinates and off-diagonals of matrix

  J      = 1
  z( 1 ) = SNGL( SSP%Depth( FirstAcoustic ) )

  DO Medium = FirstAcoustic, LastAcoustic

     h_rho = h( Medium ) * rho( Loc( Medium ) + 1 )        ! density at the top of each layer

     e( J + 1 : J + N( Medium ) ) = 1.0D0 / h_rho
     z( J + 1 : J + N( Medium ) ) = z( J ) + SNGL( h( Medium ) * [ (jj, jj = 1, N( Medium ) ) ] )

     J = J + N( Medium )
  END DO

  e( NTotal1 + 1 ) = 1.0D0 / h_rho       ! Dummy value; never used

  ! Calculate the indices, weights, ... for mode interpolation
  CALL MergeVectors( Pos%sd, Pos%Nsd, Pos%rd, Pos%Nrd, zTab, NzTab )
  ALLOCATE( WTS( NzTab ), IzTab( NzTab ), PhiTab( NzTab ) )
  CALL Weight( z, NTotal1, zTab, NzTab, WTS, IzTab )

  ! Open MODFile and write header

  IF ( iProf == 1 ) THEN
     ! LRecordLength must not increase between profiles !!!
     LRecordLength = MAX( 2 * NzTab, 32, 3 * ( LastAcoustic - FirstAcoustic + 1 ) )   ! Logical record length in `longwords' (4 bytes)
     OPEN ( FILE = TRIM( FileRoot) //'.mod', UNIT = MODFile, ACCESS = 'DIRECT', RECL = 4 * LRecordLength, FORM = 'UNFORMATTED' )
  END IF

  WRITE( MODFile, REC = IRecProfile     ) LRecordLength, Title, REAL( Freq ), LastAcoustic - FirstAcoustic + 1, NzTab, NzTab
  WRITE( MODFile, REC = IRecProfile + 1 ) ( N( Medium ), SSP%Material( Medium ), Medium = FirstAcoustic, LastAcoustic )
  WRITE( MODFile, REC = IRecProfile + 2 ) HSTop%BC, CMPLX( HSTop%cP ), CMPLX( HSTop%cS ), REAL( HSTop%rho ), &
                                          REAL( SSP%Depth( 1          ) ), &
                                          HSBot%BC, CMPLX( HSBot%cP ), CMPLX( HSBot%cS ), REAL( HSBot%rho ), &
                                          REAL( SSP%Depth( SSP%NMedia + 1 ) )
  WRITE( MODFile, REC = IRecProfile + 3 ) ( REAL( SSP%Depth( Medium ) ), REAL( rho( Loc( Medium ) + 1 ) ), &
                                          Medium = FirstAcoustic, LastAcoustic )
  WRITE( MODFile, REC = IRecProfile + 5 ) zTab( 1 : NzTab )

  ! Main loop: for each eigenvalue call SInverseIteration to get eigenvector

  ModeLoop: DO mode = 1, M
     x = EVMat( 1, mode )

     ! Corner elt requires top impedance
     CALL BCImpedance( x, 'TOP', HSTop, fTop, gTop, iPower )

     IF ( gTop == 0.0D0 ) THEN
        d( 1 ) = 1.0D0
        e( 2 ) = 0.0D0
     ELSE
        L      = Loc( FirstAcoustic ) + 1
        xh2    = x * h( FirstAcoustic ) * h( FirstAcoustic )
        h_rho  = h( FirstAcoustic ) * rho( L )
        d( 1 ) = ( B1( L ) - xh2 ) / h_rho / 2.0D0 + fTop / gTop
     END IF

     ! Set up the diagonal
     ITP = NTotal
     J   = 1
     L   = Loc( FirstAcoustic ) + 1

     Media: DO Medium = FirstAcoustic, LastAcoustic
        xh2   = x * h( Medium ) ** 2
        h_rho = h( Medium ) * rho( Loc( Medium ) + 1 )

        IF ( Medium >= FirstAcoustic + 1 ) THEN
           L      = L + 1
           d( J ) = ( d( J ) + ( B1( L ) - xh2 ) / h_rho ) / 2.0D0
        END IF

        DO ii = 1, N( Medium )
           J      = J + 1
           L      = L + 1
           d( J ) = ( B1( L ) - xh2 ) / h_rho

           IF ( REAL( B1( L ) - xh2 ) + 2.0D0 > 0.0D0 ) THEN   ! Find index of turning point nearest top
              ITP = MIN( J, ITP )
           END IF
        END DO

     END DO Media

     ! Corner elt requires bottom impedance
     CALL BCImpedance( x, 'BOT', HSBot, fBot, gBot, iPower )

     IF ( gBot == 0.0D0 ) THEN
        d( NTotal1 ) = 1.0D0
        e( NTotal1 ) = 0.0D0
     ELSE
        d( NTotal1 ) = d( NTotal1 ) / 2.0D0 - fBot / gBot
     END IF

     CALL InverseIteration( NTotal1, d, e, IERR, Phi )   ! Inverse iteration to compute eigenvector

     IF ( IERR /= 0 ) THEN
        WRITE( PRTFile, * ) 'mode = ', mode
        CALL ERROUT( PRTFile, 'W', 'KRAKENC-SInverseIteration', 'Inverse iteration failed to converge' )
        ! EVMat(  1, mode ) = 0.0D0   ! remove that eigenvalue
        ! Extrap( 1, mode ) = 0.0D0
     ELSE
        CALL Normalize( Phi, ITP, NTotal1, x )   ! Normalize the eigenvector

        ! Tabulate the modes at the source/rcvr depths and write to disk
        PhiTab = CMPLX( Phi( IzTab ) ) + WTS * CMPLX( Phi( IzTab + 1 ) - Phi( IzTab ) )
        WRITE( MODFile, REC = IRecProfile + 5 + mode ) PhiTab
     END IF
  END DO ModeLoop

  DEALLOCATE( WTS, IzTab, PhiTab )

END SUBROUTINE Vector

!**********************************************************************

SUBROUTINE Normalize( Phi, ITP, NTotal1, x )

  ! Normalize the eigenvector:
  ! SqNorm = Integral(Phi ** 2) by the trapezoidal rule:
  ! Integral( F ) = H * ( F(1) + ... + F(N-1) + 0.5 * ( F(0) + F(N) ) )

  ! Compute perturbation due to material absorption
  ! Call ScatterLoss to figure interfacial scatter loss

  USE KrakencMod
  USE sspMod
  IMPLICIT NONE
  INTEGER,           INTENT( IN    ) :: ITP, NTotal1     ! index of turning point, number of points in mesh
  COMPLEX  (KIND=8), INTENT( IN    ) :: x                ! eigenvalue
  COMPLEX  (KIND=8), INTENT( INOUT ) :: Phi( NTotal1 )   ! eigenvector
  INTEGER                            :: iPower, J, J1, L, L1, Medium
  REAL     (KIND=8)                  :: rhoMedium, rho_omega_h2
  COMPLEX  (KIND=8)                  :: x1, x2, DrhoDx, DetaDx, SqNorm, RN, ScaleFactor, Slow
  COMPLEX  (KIND=8)                  :: fTop1, gTop1, fTop2, gTop2, fBot1, gBot1, fBot2, gBot2

  SqNorm = 0.0D0
  Slow   = 0.0D0

  ! Compute contribution from the top half-space
  IF ( HSTop%BC == 'A' ) THEN
     Slow  = Slow + Phi( 1 ) ** 2 / ( 2 * SQRT( x - omega2 / HSTop%cP ** 2 ) ) / ( HSTop%rho * HSTop%cP ** 2 )
  END IF

  ! Compute contribution from the volume
  L = Loc( FirstAcoustic )
  J = 1

  Media: DO Medium = FirstAcoustic, LastAcoustic
     L            = L + 1
     rhoMedium    = rho( L )
     rho_omega_h2 = rhoMedium * omega2 * h( Medium ) ** 2

     ! top interface
     SqNorm = SqNorm + 0.5D0 * h( Medium ) *                    Phi( J ) ** 2 / rhoMedium
     Slow   = Slow   + 0.5D0 * h( Medium ) * ( B1( L ) + 2. ) * Phi( J ) ** 2 / rho_omega_h2

     ! medium
     L1 = L + 1
     L  = L + N( Medium ) - 1
     J1 = J + 1
     J  = J + N( Medium ) - 1

     SqNorm = SqNorm + h( Medium ) * SUM(                         Phi( J1 : J ) ** 2 ) / rhoMedium
     Slow   = Slow   + h( Medium ) * SUM( ( B1( L1 : L ) + 2. ) * Phi( J1 : J ) ** 2 ) / rho_omega_h2

     ! bottom interface
     L = L + 1
     J = J + 1

     SqNorm = SqNorm + 0.5D0 * h( Medium ) *                    Phi( J ) ** 2 / rhoMedium
     Slow   = Slow   + 0.5D0 * h( Medium ) * ( B1( L ) + 2. ) * Phi( J ) ** 2 / rho_omega_h2

  END DO Media

  ! Compute contribution from the bottom half-space

  IF ( HSBot%BC == 'A' ) THEN
     Slow  = Slow + Phi( J ) ** 2 / ( 2 * SQRT( x - omega2 / HSBot%cP ** 2 ) ) / ( HSBot%rho * HSBot%cP ** 2 )
  END IF

  ! Compute derivative of top admitance
  x1 = 0.9999999D0 * x
  x2 = 1.0000001D0 * x

  CALL BCImpedance( x1, 'TOP', HSTop, fTop1, gTop1, iPower )
  CALL BCImpedance( x2, 'TOP', HSTop, fTop2, gTop2, iPower )

  DrhoDx = 0.0D0
  IF ( gTop1 /= 0.0D0 ) DrhoDx = ( fTop2 / gTop2 - fTop1 / gTop1 ) / ( x2 - x1 )

  ! Compute derivative of bottom admitance

  CALL BCImpedance( x1, 'BOT', HSBot, fBot1, gBot1, iPower )
  CALL BCImpedance( x2, 'BOT', HSBot, fBot2, gBot2, iPower )

  DetaDx = 0.0D0
  IF ( gBot1 /= 0.0D0 ) DetaDx = ( fBot2 / gBot2 - fBot1 / gBot1 ) / ( x2 - x1 )

  ! Scale the mode
  RN          = SqNorm - DrhoDx * Phi( 1 ) ** 2 + DetaDx * Phi( NTotal1 ) ** 2
  ScaleFactor = 1.0D0 / SQRT( RN )
  IF ( REAL( Phi( ITP ) ) < 0.0D0 ) ScaleFactor = -ScaleFactor  ! make sign consistent at mode turning point

  Phi        = ScaleFactor * Phi
  Slow       = ScaleFactor ** 2 * Slow * SQRT( omega2 / x )
  VG( mode ) = DBLE( 1 / Slow )

  CALL ScatterLoss( Phi, x )   ! Compute interfacial scatter loss

END SUBROUTINE Normalize

!**********************************************************************

SUBROUTINE ScatterLoss( Phi, x )

  ! Figure scatter loss

  USE KrakencMod
  USE sspMod
  IMPLICIT NONE
  COMPLEX  (KIND=8), INTENT( IN ) :: Phi( * ), x   ! eigenvector, eigenvalue
  INTEGER           :: J, L, Medium
  REAL     (KIND=8) :: rho1, rho2, rhoInside, h2
  COMPLEX  (KIND=8) :: Perturbation_k, eta1Sq, eta2Sq, KupIng, U, PekerisRoot

  Perturbation_k = 0.0D0
  J = 1
  L = Loc( FirstAcoustic )

  Media: DO Medium = FirstAcoustic - 1, LastAcoustic   ! Loop over media

     ! Calculate rho1, eta1Sq, Phi, U

     IF ( Medium == FirstAcoustic - 1 ) THEN   ! Top properties
        SELECT CASE ( HSTop%BC )
        CASE ( 'A' )          ! Acousto-elastic
           rho1      = HSTop%rho
           eta1Sq    = x - omega2 / HSTop%cP ** 2
           U         = PekerisRoot( eta1Sq ) * Phi( 1 ) / HSTop%rho
        CASE ( 'V' )          ! Vacuum
           rho1      = 1.0D-9
           eta1Sq    = 1.0D0
           rhoInside = rho( Loc( FirstAcoustic ) + 1 )
           U         = Phi( 2 ) / h( FirstAcoustic ) / rhoInside
        CASE ( 'R' )          ! Rigid
           rho1      = 1.0D+9
           eta1Sq    = 1.0D0
           U         = 0.0D0
        CASE DEFAULT          ! Tabulated
           rho1      = 0.0D0
           eta1Sq    = 0.0D0
           U         = 0.0D0
        END SELECT
     ELSE
        h2 = h( Medium ) ** 2
        J  = J + N( Medium )
        L  = Loc( Medium ) + N( Medium ) + 1

        rho1   = rho( L )
        eta1Sq = ( 2.0D0 + B1( L ) ) / h2 - x
        U      = ( -Phi( J - 1 ) - 0.5D0 * ( B1( L ) - h2 * x ) * Phi( J ) ) / ( h( Medium ) * rho1 )
     END IF

     ! Calculate rho2, eta2

     IF ( Medium == LastAcoustic ) THEN   ! Bottom properties
        SELECT CASE ( HSBot%BC )
        CASE ( 'A' )          ! Acousto-elastic
           rho2   = HSBot%rho
           eta2Sq = omega2 / HSBot%cP ** 2 - x
        CASE ( 'V' )          ! Vacuum
           rho2   = 1.0D-9
           eta2Sq = 1.0D0
        CASE ( 'R' )          ! Rigid
           rho2   = 1.0D+9
           eta2Sq = 1.0D0
        CASE DEFAULT          ! Tabulated
           rho2   = 0.0D0
           eta2Sq = 0.0D0
        END SELECT
     ELSE
        rho2   = rho( L + 1 )
        eta2Sq = ( 2.0D0 + B1( L + 1 ) ) / h( Medium + 1 ) ** 2 - x
     END IF

     Perturbation_k = Perturbation_k + KupIng( SSP%sigma( Medium + 1 ), eta1Sq, rho1, eta2Sq, rho2, Phi( J ), U )

  END DO Media

  k( mode ) = Perturbation_k

END SUBROUTINE ScatterLoss

!**********************************************************************

SUBROUTINE ORDER( x, N )

  ! Does an insertion sort of the complex vector x in order of decreasing real part

  ! At the Ith step, the first I-1 positions contain a sorted vector.
  ! We shall insert the Ith value into its place in that
  ! vector, shifting up to produce a new vector of length I.

  IMPLICIT NONE
  INTEGER, INTENT( IN ) :: N
  INTEGER               :: ILeft, IMiddle, IRight, I
  COMPLEX (KIND=8)      :: x( N ), xTemp

  IF ( N == 1 ) RETURN

  DO I = 2, N

     xTemp = x( I )

     IF ( REAL( xTemp ) > REAL( x( 1 ) ) ) THEN
        x( 2 : I ) = x( 1 : I - 1 )
        x( 1 )     = xTemp  ! goes in the first position
     ELSE IF ( REAL( xTemp ) > REAL( x( I - 1 ) ) ) THEN

        ! Binary search for its place

        IRight = I - 1
        ILeft  = 1

        DO WHILE ( IRight > ILeft + 1 )
           IMiddle = ( ILeft + IRight ) / 2

           IF ( REAL( xTemp ) > REAL( x( IMiddle ) ) ) THEN
              IRight = IMiddle
           ELSE
              ILeft  = IMiddle
           END IF
        END DO

        x( IRight + 1 : I ) = x( IRight : I - 1 )
        x( IRight ) = xTemp

     END IF

  END DO

END SUBROUTINE ORDER
