      LOGICAL FUNCTION ABBRV(SHORT,LONG,MNCH)

C#### Function: ABBRV
C###  Type: LOGICAL
C###  Description:
C###    ABBRV returns .TRUE. if the character string SHORT is an
C###    abbreviation of LONG.  SHORT must be at least MNCH
C###    characters long.

      IMPLICIT NONE
!     Parameter List
      INTEGER MNCH
      CHARACTER LONG*(*),SHORT*(*)
!     Local Variables
      INTEGER noch,NXCH
C      INTEGER LEN_LONG,LEN_SHORT
C      CHARACTER NUL*1,TAB*1
C      LOGICAL SAME
C      PARAMETER (NUL=CHAR(0),TAB=CHAR(9))

      ABBRV=.FALSE.
      NXCH=MIN(LEN(LONG),LEN(SHORT))
      DO noch=MNCH,NXCH
        IF(SHORT.EQ.LONG(:noch)) THEN
          ABBRV=.TRUE.
          RETURN
        ENDIF
      ENDDO !noch

C KAT 10Oct98: alternative suggestion
C      DO noch=1,MNCH
C        IF(SHORT(noch:noch).NE.LONG(noch:noch)) THEN
C          ABBRV=.FALSE.
C          RETURN
C        ENDIF
C      ENDDO
C      LEN_LONG=LEN(LONG)
C      LEN_SHORT=LEN(SHORT)
C      NXCH=MIN(LEN_LONG,LEN_SHORT)
C      SAME=.TRUE.
C      noch=MNCH+1
C      DO WHILE(SAME.AND.noch.LE.NXCH)
C        IF(SHORT(noch:noch).NE.LONG(noch:noch)) THEN
C          SAME=.FALSE.
C        ELSE
C          noch=noch+1
C        ENDIF
C      ENDDO
C      IF(noch.LE.LEN_SHORT) THEN
C        ABBRV=SHORT(noch:noch).EQ.' '.
C     '    OR.SHORT(noch:noch).EQ.TAB.
C     '    OR.SHORT(noch:noch).EQ.NUL
C      ELSE
C        ABBRV=.TRUE.
C      ENDIF

      RETURN
      END



