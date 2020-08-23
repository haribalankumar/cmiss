      LOGICAL FUNCTION MADEOF(CDATA,CLIST)

C#### Function: MADEOF
C###  Type: LOGICAL
C###  Description:
C###    MADEOF returns .TRUE. if the character string CDATA is made
C###    up of the characters in CLIST and no others.

      IMPLICIT NONE
!     Parameter List
      CHARACTER CDATA*(*),CLIST*(*)
!     Local Variables
      INTEGER IDATA,ilist

C     CALL ENTERS('MADEOF',*9999)
      DO 2 IDATA=1,LEN(CDATA)
        DO 1 ilist=1,LEN(CLIST)
          IF(CDATA(IDATA:IDATA).EQ.CLIST(ilist:ilist)) GOTO 2
    1   CONTINUE
        MADEOF=.FALSE.
C        CALL EXITS('MADEOF')
        RETURN
    2 CONTINUE
      MADEOF=.TRUE.

C     CALL EXITS('MADEOF')
      RETURN
      END

      
