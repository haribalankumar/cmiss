      LOGICAL FUNCTION ENCASPULATED(ELEM,CIRC,NODE)

C#### Function: ENCAPSULATED
C###  Description: Determines whether the internal point is
C###   encapsulated within the tetrahedral. Genmesh function.
CC JMB 22-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ELEM(LDELEM,4)
      REAL*8 CIRC(LDCIRC,3),NODE(LDNODE,3)
!     Local Variables
      INTEGER I,J
      REAL*8 COORD(3),DTRMNT(4),SUM,WORK(4,4),XI(4)
!     Local Functions
      REAL*8 DETERMINANT,IND

      ENCASPULATED=.TRUE.
      DO I=1,NPTS
        CALL DCOPY(4,NODE(ELEM(1,I),1),LDNODE,WORK(1,I),1)
      ENDDO !i
      DO J=1,NJT
        SUM=ZERO
        DO I=1,NPTS
          SUM=SUM+NODE(ELEM(1,I),J)
        ENDDO !i
        COORD(J)=0.25d0*SUM
      ENDDO !j
      XI(1)=IND(WORK(1,1),WORK(1,2),WORK(1,3),COORD)
      XI(2)=IND(WORK(1,2),WORK(1,4),WORK(1,3),COORD)
      XI(3)=IND(WORK(1,3),WORK(1,4),WORK(1,1),COORD)
      XI(4)=IND(WORK(1,4),WORK(1,2),WORK(1,1),COORD)
      DO I=1,NPTS
        DO J=1,NJT
          WORK(I,J)=WORK(I,J)-CIRC(1,J)
        ENDDO !j
      ENDDO !i
      DTRMNT(1)=DETERMINANT(1,2,3,WORK,4)*XI(1)
      DTRMNT(2)=DETERMINANT(2,4,3,WORK,4)*XI(2)
      DTRMNT(3)=DETERMINANT(3,4,1,WORK,4)*XI(3)
      DTRMNT(4)=DETERMINANT(4,2,1,WORK,4)*XI(4)
      ENCASPULATED=.TRUE.
      DO I=1,NPTS
        ENCASPULATED=ENCASPULATED.AND.(DTRMNT(I).GT.ZERO)
      ENDDO !i
      IF(.NOT.ENCASPULATED)THEN
        ENCASPULATED=.TRUE.
        DO I=1,NPTS
          ENCASPULATED=ENCASPULATED.AND.(DTRMNT(I).LT.ZERO)
        ENDDO !i
      ENDIF !encapsulated

      RETURN
      END


