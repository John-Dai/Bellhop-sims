MODULE ScooterMod

   USE MathConstantsMod

   SAVE
   INTEGER, PARAMETER :: ENVFile = 5, PRTFile = 6, GRNFile = 25, IRCFile = 12, MaxMedium = 500
   INTEGER            :: FirstAcoustic, LastAcoustic, Nk
   REAL      (KIND=8) :: Cmin, Clow, CHigh, omega2, RMax, Atten
   CHARACTER  (LEN=8) :: TopOpt, BotOpt

  ! Halfspace properties
  TYPE HSInfo
     CHARACTER (LEN=1) :: BC                          ! Boundary condition type
     COMPLEX  (KIND=8) :: cP, cS                      ! P-wave, S-wave speeds
     REAL     (KIND=8) :: rho, BumpDensity, eta, xi   ! density, boss parameters
  END TYPE HSInfo

  TYPE( HSInfo )       :: HSTop, HSBot

  ! finite-difference grid
  INTEGER                       :: Loc( MaxMedium ), N( MaxMedium )
  REAL      (KIND=8)            :: H(   MaxMedium )

  ! storage for finite-difference equations
  REAL    (KIND=8),    ALLOCATABLE :: rho( : )
  COMPLEX (KIND=8),    ALLOCATABLE :: B1( : ), B2( : ), B3( : ), B4( : )

END MODULE ScooterMod
