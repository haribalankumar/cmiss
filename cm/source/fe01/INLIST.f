      LOGICAL FUNCTION INLIST(i,A_LIST,NTLIST,N1LIST)

C#### Function: INLIST
C###  Type: LOGICAL
C###  Description:
C###    INLIST returns .TRUE. if i belongs to
C###    A_LIST(nolist),nolist=1,NTLIST.
C###    As a side effect, N1LIST returns the position of the first
C###    located incidence of i in A_LIST.

      IMPLICIT NONE
!     Parameter List
      INTEGER i,A_LIST(*),N1LIST,NTLIST
!     Local Variables
      INTEGER nolist

C     CALL ENTERS('INLIST',*9999)
      INLIST=.FALSE.
      DO nolist=1,NTLIST
        IF(i.EQ.A_LIST(nolist)) THEN
          INLIST=.TRUE.
          N1LIST=nolist
          GO TO 900
        ENDIF
      ENDDO
C900  CALL EXITS('INLIST')
 900  RETURN
      END


