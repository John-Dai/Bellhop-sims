SUBROUTINE GetModes( FileRoot, IProf, MaxM, rd, Nrd, Comp, k, PhiR, M, Freq, Title )                                     

  ! Read in modes and extract values at rcvr depths                   

  ! If IProf = 1 then the file is opened and IRecProfile is set to the first record          
 
  IMPLICIT NONE
  INTEGER, PARAMETER     :: PRTFile = 6, MaxN = 16001, MaxMedium = 51
  INTEGER, INTENT( IN  ) :: IProf               ! profile number for the mode set
  INTEGER, INTENT( IN  ) :: MaxM                ! row dimension of PhiR in calling program
  INTEGER, INTENT( IN  ) :: Nrd                 ! number of receiver depths (where modes are sampled)
  REAL,    INTENT( IN  ) :: rd( Nrd )           ! receiver depths
  CHARACTER (LEN=*), INTENT( IN  ) :: FileRoot  ! Name of the mode file
  INTEGER, INTENT( OUT ) :: M                   ! number of modes
  REAL,    INTENT( OUT ) :: freq                ! frequency for which modes were calculated
  COMPLEX, INTENT( OUT ) :: k( MaxM )           ! modal wavenumbers
  COMPLEX, INTENT( OUT ) :: PhiR( MaxM, Nrd )   ! mode shapes
  CHARACTER (LEN=*), INTENT( OUT ) :: Title     ! title in mode file
  CHARACTER (LEN=1), INTENT( OUT ) :: Comp      ! component (vertical, horizontal,...) extracted (ignored if no shear)
  INTEGER            :: ir, iz, IRecProfile, LRecl, Mode, NMat, NMedia, NTot, N( MaxMedium ), ird( Nrd ) 
  REAL               :: Z( MaxN ), W( Nrd ), Depth( MaxMedium ), rho( MaxMedium ), Tol, wT
  REAL               :: rhoT, rhoB, csT, csB, depthT, depthB   ! halfspace density, shear, depth
  COMPLEX            :: PhiT( MaxN ), kTop2, kBot2
  CHARACTER (LEN=8)  :: Material( MaxMedium )
  CHARACTER (LEN=1)  :: BCTop, BCBot            ! boundary condition type for top and bottom
  CHARACTER (LEN=7)  :: Model                   ! KRAKEN or KRAKENC
  SAVE iRecProfile

  ! Read the header data from the mode file
  CALL ReadModeHeader( FileRoot, IProf, IRecProfile, LRecl, Title, Freq, NMedia, NTot, NMat, N, Material, Depth, rho,  &
       BCTop, rhoT, DepthT, & 
       BCBot, rhoB, DepthB, M, MaxM, Z, k, kTop2, kBot2 )

  IF ( M <= 0 ) RETURN  ! no modes? quick return.
  CALL Weight( Z, NTot, rd, Nrd, W, ird )   ! Locate indices of receiver points

  ! Loop over receiver depths to check for safe interpolation
  ! Receivers must be within a fraction of a wavelength
  ! of tabulated pts. We accept one wavelength at 1500 m/s

  Tol = 1500.0 / Freq 
  Model = Title( 1 : 7 )   ! Should be KRAKEN or KRAKENC

  RcvrDepth: DO ir = 1, Nrd 

     iz = ird( ir ) 
     WT = ABS( MIN( W( ir ), 1.0 - W( ir ) ) ) 

     IF ( rd( ir ) < DepthT ) THEN        ! Rcvr in upper halfspace
        cST = 0.0   ! should be passed by ReadModeHeader
        IF ( cST /= 0.0 .OR. BCTop(1:1) /= 'A' ) THEN 
           WRITE( PRTFile, * ) 'Receiver depth: ', rd( ir )
           WRITE( PRTFile, * ) 'Highest valid depth: ', DepthT
           CALL ERROUT( PRTFile, 'W', 'GetMode', 'Rcvr above highest valid depth' )
        ENDIF

     ELSE IF ( rd( ir ) > DepthB ) THEN   ! Rcvr in lower halfspace
        cSB = 0.0   ! should be passed by ReadModeHeader
        IF ( cSB /= 0.0 .OR. BCBot(1:1) /= 'A' ) THEN 
           WRITE( PRTFile, * ) 'Receiver depth: ', rd( ir )
           WRITE( PRTFile, * ) 'Lowest valid depth: ', DepthB
           CALL ERROUT( PRTFile, 'W', 'GetMode', 'Rcvr below lowest valid depth' )
        ENDIF

     ELSE IF ( NTot > 1 ) THEN            ! Rcvr between two grid points or large extrapolation
        IF ( WT * ( Z( iz + 1 ) - Z( iz ) ) > Tol ) THEN
           WRITE( PRTFile, * ) 'Receiver depth: ', rd( ir )
           WRITE( PRTFile, * ) 'Nearest depths: ', Z( iz ), Z( iz + 1 )
           WRITE( PRTFile, * ) 'Tolerance: ', Tol 
           CALL ERROUT( PRTFile, 'W', 'GetMode', 'Modes not tabulated near requested pt.' )
        ENDIF
     ELSE                                 ! Rcvr near a single grid point
        IF ( ABS( rd( ir ) - Z( iz ) ) > Tol ) THEN 
           WRITE( PRTFile, * ) 'Rd, Tabulation depth ', rd( ir), Z( iz )
           WRITE( PRTFile, * ) 'Tolerance: ', Tol 
           CALL ERROUT( PRTFile, 'W', 'GetMode', 'Modes not tabulated near requested pt.' )
        ENDIF
     ENDIF

  END DO RcvrDepth

  ! Read in the modes
  ModeLoop: DO Mode = 1, M
     CALL ReadOneMode( Mode, IRecProfile, NTot, NMat, W, ird, N, Material, NMedia, Comp, &
          kTop2, DepthT, BCTop, kBot2, DepthB, BCBot, rd, Nrd, k, PhiT, Model )
     PhiR( Mode, 1:Nrd ) = PhiT( 1:Nrd )
  END DO ModeLoop

  IRecProfile = IRecProfile + 6 + M + 1 + ( 2 * M - 1 ) / LRecL   ! advance to next profile

END SUBROUTINE GetModes

!**********************************************************************C

SUBROUTINE ReadModeHeader( FileRoot, IProf, IRecProfile, LRecl, Title, Freq, NMedia, NTot, NMat, N, Material, Depth, rho, &
     &   BCTop, rhoT, DepthT, BCBot, rhoB, DepthB, M, MaxM, Z, k, kTop2, kBot2 )

  ! Reads the header information from ModeFile                          
  ! Note T suffix means top
  !      B suffix means bottom                                        

  ! IProf     is a profile number
  ! FileRoot  is the user-provided file name                             
  ! FileNameT is the temporary name we build
  ! These have to be two separate variables to ensure there
  ! is space allocated to construct the file name even when           
  ! the user calls us with FileName = ' ' 

  ! IRecProfile must point to the first record of the profile                             

  IMPLICIT NONE
  INTEGER, PARAMETER :: PRTFile = 6, ModeFile = 30 
  REAL,    PARAMETER :: pi = 3.141592
  LOGICAL            :: OpenFlag
  INTEGER            :: N( * ) , NMedia, NTot, NMat, M, MaxM, IFirst, ILast, iProf, &
                        iostat, IRecProfile, IRec, LRecL, Medium
  REAL               :: Depth( * ), Z( * ), rho( * ), Freq, rhoT, DepthT, rhoB, DepthB 
  COMPLEX            :: k( * ), cPT, cST, cPB, cSB, kTop2, kBot2
  CHARACTER (LEN= *) :: Title
  CHARACTER (LEN=80) :: FileNameT
  CHARACTER (LEN= *) :: FileRoot
  CHARACTER (LEN= 8) :: Material( * )
  CHARACTER (LEN= 1) :: BCTop, BCBot

  ! open ModeFile
  FileNameT = TRIM( FileRoot ) // '.mod'

  INQUIRE( FILE = FileNameT, OPENED = OpenFlag )
  IF ( .NOT. OpenFlag ) THEN
     OPEN( UNIT = ModeFile, FILE = FileNameT, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', &
          RECL = 100, IOSTAT = iostat )
     IF ( IOSTAT /= 0 ) THEN
        WRITE( PRTFile, * ) 'Mode file = ', FileNameT
        CALL ERROUT( PrtFile, 'F', 'GetMode - ReadModeHeader', 'Unable to open the mode file' )
     END IF

     READ( ModeFile, REC = 1 ) LRecL
     CLOSE( UNIT = ModeFile )
     OPEN(  UNIT = ModeFile, FILE = FileNameT, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', &
          RECL = 4 * LRecL, IOSTAT = iostat )
  END IF

  ! If this is the first profile, reset the record counter to the beginning of the file
  IF ( IProf == 1 ) THEN
     IRecProfile = 1   ! set counter pointing to the first record to read
  END IF

  ! Read header info
  READ( ModeFile, REC = IRecProfile     ) LRecL, Title( 1 : 80 ), Freq, NMedia, NTot, NMat  
  READ( ModeFile, REC = IRecProfile + 1 ) ( N( Medium ), Material( Medium ), Medium = 1, NMedia)
  READ( ModeFile, REC = IRecProfile + 2 ) BCTop(1:1), cPT, cST, rhoT, DepthT, BCBot(1:1), cPB, cSB, rhoB, DepthB
  READ( ModeFile, REC = IRecProfile + 3 ) ( Depth( Medium ), rho( Medium ), Medium = 1, NMedia )
  READ( ModeFile, REC = IRecProfile + 4 ) M, LRecL


  IF ( M == 0 ) RETURN 
  READ( ModeFile, REC = IRecProfile + 5 ) Z( 1 : NTot ) 

  ! Read in eigenvalues, k( I )
  ! They come at the end of the ModeFile, because they're calculated after the eigenvectors
  IFirst = 1 
  DO IRec = 1, 1 + ( 2 * MIN( M, MaxM ) - 1 ) / LRecL
     ILast = MIN( M, IFirst + LRecL / 2 - 1 )
     READ( ModeFile, REC = IRecProfile + 5 + M + IRec ) k( IFirst : ILast )      
     IFirst = ILast + 1 
  END DO

  IF ( BCTop( 1 : 1 ) == 'A' ) kTop2 = ( 2.0 * pi * Freq / cPT ) **2 
  IF ( BCBot( 1 : 1 ) == 'A' ) kBot2 = ( 2.0 * pi * Freq / cPB ) **2

  ! Reduce the number of modes, if the storage is inadequate
  IF ( M > MaxM ) THEN 
     WRITE( PRTFile, * ) 'M = ', M, '   MaxM = ', MaxM
     CALL ERROUT( PRTFile, 'W', 'ReadModeHeader', 'Insufficient storage to read all the modes: increase MaxM' )
     M = MaxM
  ENDIF

END SUBROUTINE ReadModeHeader

!**********************************************************************C

SUBROUTINE ReadOneMode( Mode, IRecProfile, NTot, NMat, W, ird, N, Material, NMedia, Comp, &
     &   kTop2, DepthT, BCTop, kBot2, DepthB, BCBot, rd, Nrd, k, PhiR , Model )

  ! Read in a single eigenfunction and extract receiver values
  ! Results are returned in PhiR

  IMPLICIT NONE
  INTEGER, PARAMETER    :: ModeFile = 30
  INTEGER, INTENT( IN ) :: NMedia, Nrd, N( NMedia ), ird( Nrd ), IRecProfile, NTot, NMat
  REAL,    INTENT( IN ) :: rd( Nrd ), W( Nrd ), DepthT, DepthB
  COMPLEX, INTENT( IN ) :: k( * ), kTop2, kBot2
  CHARACTER (LEN=1), INTENT( IN ) :: Comp
  CHARACTER (LEN=8), INTENT( IN ) :: Material( NMedia )
  CHARACTER (LEN=1), INTENT( IN ) :: BCTop, BCBot
  CHARACTER (LEN=7), INTENT( IN ) :: Model                   ! KRAKEN or KRAKENC
  COMPLEX, INTENT( OUT) :: PhiR( Nrd )
  LOGICAL           :: TufLuk 
  INTEGER           :: ir, iz, j, Mode
  COMPLEX           :: Phi( NMat ), gammaT = 0.0, gammaB = 0.0
  COMPLEX  (KIND=8) :: PekerisRoot, gamma2

  READ( ModeFile, REC = IRecProfile + 5 + Mode ) ( Phi( j ), j = 1, NMat ) 

  ! Is there an elastic medium in the problem?
  TufLuk = .FALSE. 
  IF ( ANY( Material( 1 : NMedia ) == 'ELASTIC' ) ) TufLuk = .TRUE.

  ! Extract the component specified by 'Comp'
  IF ( TufLuk ) CALL EXTRACT( Phi, N, Material, NMedia, Comp ) 

  ! Extract values at receiver depths
  ! n.b. should be using real( k(mode) ) for KRAKEN
  IF ( BCTop == 'A' ) THEN
     SELECT CASE( Model )
     CASE ( 'KRAKENC' )
        gamma2 = k( Mode ) ** 2 - kTop2
     CASE DEFAULT   ! KRAKEN
        gamma2 = REAL( k( Mode ) ) ** 2 - kTop2
     END SELECT
     gammaT = CMPLX( PekerisRoot( gamma2 ) )
  END IF

  IF ( BCBot == 'A' ) THEN 
     SELECT CASE( Model )
     CASE ( 'KRAKENC' )
        gamma2 = k( Mode ) ** 2 - kBot2
     CASE DEFAULT   ! KRAKEN
        gamma2 = REAL( k( Mode ) ) ** 2 - kBot2
     END SELECT
     gammaB = CMPLX( PekerisRoot( gamma2 ) )
  END IF

  RcvrDepth: DO ir = 1, Nrd 
     IF ( rd( ir ) < DepthT ) THEN      ! Rcvr in upper halfspace
        PhiR( ir ) = Phi( 1    ) * EXP( -gammaT * ( DepthT - rd( ir  ) ) )

     ELSE IF ( rd( ir ) > DepthB ) THEN ! Rcvr in lower halfspace
        PhiR( ir ) = Phi( NTot ) * EXP( -gammaB * ( rd( ir ) - DepthB ) )

     ELSE IF ( NTot > 1 ) THEN 
        iz         = ird( ir )
        PhiR( ir ) = Phi( iz ) + W( ir ) * ( Phi( iz + 1 ) - Phi( iz ) )

     ELSE                               ! mode is tabulated at only one depth
        iz = ird( ir ) 
        PhiR( ir ) = Phi( iz ) 
     ENDIF

  END DO RcvrDepth

END SUBROUTINE ReadOneMode
!**********************************************************************C
SUBROUTINE EXTRACT( Phi, N, Material, NMedia, Comp ) 

  ! For elastic problems where a stress-displacement vector is output,
  ! extracts the desired component                                    

  IMPLICIT NONE
  INTEGER,           INTENT( IN    ) :: N( * ), NMedia
  CHARACTER (LEN=1), INTENT( IN    ) :: Comp
  CHARACTER (LEN=8), INTENT( IN    ) :: Material( * )
  COMPLEX,           INTENT( INOUT ) :: Phi( * ) 
  INTEGER                            :: i, j = 1, k = 1, Medium

  MediumLoop: DO Medium = 1, NMedia 

     DO i = 1, N( Medium ) + 1
        SELECT CASE ( Material( Medium ) )
        CASE ( 'ACOUSTIC' )
           Phi( j ) = Phi( k )
           k = k + 1
        CASE ( 'ELASTIC' )
           SELECT CASE ( Comp )
           CASE ( 'H' )
              Phi( j ) = Phi( k     )
           CASE ( 'V' )
              Phi( j ) = Phi( k + 1 )
           CASE ( 'T' )
              Phi( j ) = Phi( k + 2 )
           CASE ( 'N' )
              Phi( j ) = Phi( k + 3 )
           END SELECT   ! Comp

           k = k + 4
        END SELECT   ! Material
        j = j + 1 
     END DO   ! next point in medium

  END DO MediumLoop

END SUBROUTINE EXTRACT
