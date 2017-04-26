MODULE FieldSMod

  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER   :: FLPFile = 10, PRTFile = 6, GRNFile = 20, SHDFile = 25
  REAL,    PARAMETER   :: pi = 3.14159265, RadDeg = 180.0 / pi
  COMPLEX (KIND=8), PARAMETER :: i  = ( 0.0, 1.0 )
  INTEGER Nrr
  INTEGER              :: IAllocStat, Nk, NrrLast, Nt, Nt2
  REAL                 :: Atten, AttInt, Freq, kmax, Delta_k, Delta_kInterp, &
                          Rmin, Rmax, RMinKM, RMaxKM, Delta_r
  CHARACTER   (LEN=80) :: PlotTitle, FileRoot
  CHARACTER   (LEN= 4) :: Option
  REAL,    ALLOCATABLE :: k( : ), kInterp( : )
  COMPLEX, ALLOCATABLE :: G( : ), Ginterp( : ), P( : ), Temp( : )

END MODULE FieldSMod
