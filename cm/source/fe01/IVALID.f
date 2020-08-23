      LOGICAL FUNCTION IVALID(CDATA)

C#### Function: IVALID
C###  Type: LOGICAL
C###  Description:
C###    IVALID returns .TRUE. if CDATA is a valid character
C###    representation of an integer.

      IMPLICIT NONE
!     Parameter List
      CHARACTER CDATA*(*)
!     Local Variables
      INTEGER i,IC,ICDATA
      LOGICAL MADEOF

C     CALL ENTERS('IVALID',*9999)
      ICDATA=LEN(CDATA)
      DO i=1,ICDATA
        IF(CDATA(I:I).NE.' ') THEN
          IC=i
          GOTO 2
        ENDIF
      ENDDO !i
      IVALID=.FALSE.
C      CALL EXITS('IVALID')
      RETURN
    2 IF(IC.EQ.ICDATA) THEN
        IF(MADEOF(CDATA(IC:IC),'0123456789')) THEN
          IVALID=.TRUE.
        ELSE
          IVALID=.FALSE.
        ENDIF
      ELSE
        IF(MADEOF(CDATA(IC:IC),'0123456789+-')) THEN
          IF(MADEOF(CDATA(IC+1:ICDATA),'0123456789')) THEN
            IVALID=.TRUE.
          ELSE
            IVALID=.FALSE.
          ENDIF
        ELSE
          IVALID=.FALSE.
        ENDIF
      ENDIF

C     CALL EXITS('IVALID')
      RETURN
      END


