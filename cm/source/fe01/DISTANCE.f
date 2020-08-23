      REAL*8 FUNCTION DISTANCE(NDIM,A,B)

C#### Function: DISTANCE
C###  Type: REAL*8
C###  Description:
C###    Returns the distance between two points in NDIM dimensions
C**** Created by Peter Bier, April 2003
      
!     Parameter List
      INTEGER NDIM
      REAL*8 A(NDIM),B(NDIM)
!     local variables
      INTEGER nj
      REAL*8 DSQUARED

      DSQUARED = 0.0d0
      DO nj = 1,NDIM
C GBS 24-Jul-2003 Fixed computation
        DSQUARED = DSQUARED + (A(nj) - B(nj))**2
      ENDDO

      DISTANCE = DSQRT(DSQUARED)

      RETURN
      END


