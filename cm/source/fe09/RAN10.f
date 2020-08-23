      REAL*8 FUNCTION RAN10(i,j,k)
C#### Function: RAN10
C###  Type: REAL*8
C###  Description:
C###    Returns a random number with a gaussian (normal) distribution
C###    and a standard distribution of 1. The seeds are only set the
C###    first time which the function is called.
C###
C###    NOTE : Input Parameters i,j,k refer to the 3 seeds used
C###           generate the random number

      IMPLICIT NONE
!     Parameter List
      INTEGER i,j,k
!     Local Variables
      REAL*8 factor,r
      REAL*8 v1,v2
      REAL*8 RANDSEEDS

      r=2.D0
      DO WHILE (r.GT.1.D0)
        v1=2.D0*RANDSEEDS(i,j,k)-1.0D0
        v2=2.D0*RANDSEEDS(i,j,k)-1.0D0
        r=v1*v1+v2*v2
      ENDDO
      factor=DSQRT(-2.D0*log(r)/r)
      RAN10=factor*v1
      END



C FE30 Functions
C ==============

