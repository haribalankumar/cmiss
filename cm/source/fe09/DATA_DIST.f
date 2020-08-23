      REAL*8 FUNCTION DATA_DIST(nd1,nd2,ZD)

C#### Function: DATA_DIST
C###  Type: REAL*8
C###  Description:
C###    Determines the distance between 2 data points from
C###    CMISS data array ZD
C###  See-Also: NODE_DIST

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER nd1,nd2
      REAL*8 ZD(NJM,NDM)
!     Local Variables
      INTEGER nj
C      REAL*8 DATA_DIST

      DATA_DIST=0.0d0
      DO nj=1,NJT
        DATA_DIST=DATA_DIST+ (ZD(nj,nd2)-ZD(nj,nd1))**2
      ENDDO
      DATA_DIST=DSQRT(DATA_DIST)

      RETURN
      END


