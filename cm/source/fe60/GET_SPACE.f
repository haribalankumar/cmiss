      SUBROUTINE GET_SPACE(ELEM,CIRC,FSPACE,NODE)

C#### Subroutine: GET_SPACE
C###  Description:
C###    Genemesh routine.
CC JMB 15-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
!     Parameter List
      INTEGER ELEM(LDELEM,4)
      REAL*8 CIRC(LDCIRC,3),FSPACE(0:N_GM),NODE(LDNODE,3)
!     Local Variables
      INTEGER I
      REAL*8 RESULT(4),SUM,WORK(3)

      CALL DCOPY(3,CIRC(1,1),LDCIRC,WORK,1)
      CALL CONSTRUCT(ELEM,WORK,NODE,RESULT)
      SUM=ZERO
      DO I=1,NPTS
        SUM=SUM+RESULT(I)*FSPACE(ELEM(1,I))
      ENDDO !i
      FSPACE(0)=SUM

      RETURN
      END


