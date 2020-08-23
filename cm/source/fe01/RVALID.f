      LOGICAL FUNCTION RVALID(CDATA)

C#### Function: RVALID
C###  Type: LOGICAL
C###  Description:
C###    RVALID returns .TRUE. if CDATA is a valid character
C###    representation of a REAL*8 number.

      IMPLICIT NONE
!     Parameter List
      CHARACTER CDATA*(*)
!     Local Variables
      INTEGER i,IC,ICDATA
      LOGICAL MADEOF

C     CALL ENTERS('RVALID',*9999)
      ICDATA=LEN(CDATA)
      DO 1 i=1,ICDATA
        IF(CDATA(i:i).NE.' ') THEN
          IC=i
          GOTO 2
        ENDIF
    1 CONTINUE
      RVALID=.FALSE.
C      CALL EXITS('RVALID')
      RETURN
    2 CONTINUE
      IF(MADEOF(CDATA(IC:ICDATA),'0123456789+-.DdEe')) THEN
C****   A check for valid REAL*8 input must be inserted here.
        RVALID=.TRUE.
      ELSE
        RVALID=.FALSE.
      ENDIF

C     CALL EXITS('RVALID')
      RETURN
      END

      
