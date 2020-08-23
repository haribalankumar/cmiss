      LOGICAL FUNCTION DEGREE(FACE,REF,CIRC,RSQ)

C#### Function: DEGREE
C###  Description:
C###    Genmesh function
CC JMB 22-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER FACE(LDFACE,NFACES),REF
      REAL*8 CIRC(LDCIRC,3),RSQ(N_GM)
!     Local Variables
      INTEGER I, J
      REAL*8 SUM

      DO I=1,NFACES
        SUM=ZERO
        DO J=1,NJT
          SUM=SUM+(CIRC(REF,J)-CIRC(FACE(REF,I),J))**2.d0
        ENDDO !j
        IF(SUM.LT.RSQ(FACE(REF,I)))THEN
          DEGREE=.TRUE.
          GOTO 10
        ENDIF !sum
      ENDDO !i
      DEGREE = .FALSE.

   10 RETURN
      END


