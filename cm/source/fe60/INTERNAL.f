      LOGICAL FUNCTION INTERNAL(ELEM,IREGION)

C#### Function: INTERNAL
C###  Description:
C###    Tests whether the tetrahedral is internal. Genmesh routine.
CC JMB 20-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
!     Parameter List
      INTEGER ELEM(LDELEM,4),IREGION(0:N_GM)
!     Local Variables
      INTEGER I

      INTERNAL=.TRUE.
      DO I=1,NPTS
        INTERNAL=INTERNAL.AND.(IREGION(ELEM(1,I)).GE.2)
      ENDDO !i

      RETURN
      END


