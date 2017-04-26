SUBROUTINE ReadHeader( SHDFile, FileName, Title, Freq, Atten, PlotType )

  ! Read header from disk file

  ! FileName is a SHDFIL for complex pressure or a GRNFIL for a Green's function
  ! Title  arbitrary title
  ! theta  vector of bearing lines,   theta( 1: Ntheta )
  ! sd     vector of source   depths, sd( 1 : Nsd )
  ! rd     vector of receiver depths, rd( 1 : Nrd )
  ! r      vector of receiver ranges, r(  1 : Nr  )
  ! Freq   frequency
  ! Atten  stabilizing attenuation (which is only important for FFP runs. Use zero for other cases.)

  USE SdRdRMod
  IMPLICIT NONE
  INTEGER, PARAMETER                :: PRTFile = 6
  REAL,               INTENT( OUT ) :: Freq, Atten
  INTEGER                           :: IAllocStat, IOStat, LRecL, SHDFile
  CHARACTER (LEN=80)                :: Title, FileName
  CHARACTER (LEN=10), INTENT( OUT ) :: PlotType

  ! Open file, read header
  IF ( SHDFile == 0 ) SHDFile = 25
  IF ( FileName( 1 : 1 ) == ' ' ) FileName = 'SHDFIL'

  ! INQUIRE( FILE = FileName, RECL = IRECL )
  OPEN( UNIT = SHDFile,   FILE = FileName, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', RECL = 4, IOSTAT = IOStat )
  IF ( IOStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadHeader', 'Unable to open shade file' )
  READ( SHDFile, REC = 1 ) LRecl
  CLOSE( UNIT = SHDFile )
  OPEN(  UNIT = SHDFile,   FILE = FileName, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', RECL = 4 * LRecl )

  READ( SHDFile, REC = 1 ) LRecl, Title
  READ( SHDFile, REC = 2 ) PlotType
  READ( SHDFile, REC = 3 ) Freq, Pos%Ntheta, Pos%Nsx, Pos%Nsy, Pos%Nsd, Pos%Nrd, Pos%Nr, atten

  ALLOCATE( Pos%sd( Pos%Nsd ), Pos%rd( Pos%Nrd ), Pos%r( Pos%Nr ), Pos%theta( Pos%Ntheta ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadHeader', 'Too many source/receiver combinations' )

  READ( SHDFile, REC = 4 ) Pos%theta
  READ( SHDFile, REC = 5 ) Pos%sd
  READ( SHDFile, REC = 6 ) Pos%rd
  READ( SHDFile, REC = 7 ) Pos%r

  ! Pos%deltaR = Pos%r( Pos%Nr ) - Pos%r( Pos%Nr - 1 )

END SUBROUTINE ReadHeader
