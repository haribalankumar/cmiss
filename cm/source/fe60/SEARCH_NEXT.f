      INTEGER FUNCTION SEARCH_NEXT(ELEM,INODE,NODE)

C#### Function: SEARCH_NEXT
C###  Description:
C###    Determines the corresponding neighbouring face
C###    tetrahedral closest to the node INODE. Genmesh function
CC JMB 15-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ELEM(LDELEM,4),INODE
      REAL*8 NODE(LDNODE,3)
!     Local Variables
      INTEGER I, INDEX, LWORK
      PARAMETER (LWORK=4*8)
      REAL*8 MINVAL, RESULT(4),WORK(LWORK)
!     Local Functions
      REAL*8 DLAMCH

      CALL DCOPY(NJT,NODE(INODE,1),LDNODE,WORK,1)
      CALL CONSTRUCT(ELEM,WORK,NODE,RESULT)

      MINVAL=RESULT(1)
      INDEX=1
      DO I=2,NPTS
        IF(RESULT(I).LT.MINVAL)THEN
          MINVAL=RESULT(I)
          INDEX=I
        ENDIF !result
      ENDDO !i
      IF(RESULT(INDEX).LT.-SQRT(DLAMCH('E')))THEN
        SEARCH_NEXT=INDEX
      ELSE
        SEARCH_NEXT=0
      ENDIF !result

      RETURN
      END


