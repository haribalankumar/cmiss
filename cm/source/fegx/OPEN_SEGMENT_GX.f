      SUBROUTINE OPEN_SEGMENT_GX(ISEGNUM,ISEG,ERROR,*)

C#### Subroutine: OPEN_SEGMENT_GX
C###  Description:
C###    Opens graphics segment ISEGNUM

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='OPEN_SEGMENT_GX')
      INCLUDE 'gx.inc'
!     Parameter List
      INTEGER ISEG(*),ISEGNUM
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERR

      CALL ENTERS(ROUTINENAME,*9999)

      CALL OPOBJT(ISEGNUM,ERR)
      IF(err.ne.0) THEN
        ERROR='>>Error from GX'
        GOTO 9999
      ENDIF
      IF(ISEG(ISEGNUM).EQ.1) THEN
        CALL SOBVIS(ISEGNUM,gxNOVIS,ERR)
        IF(err.ne.0) THEN
          ERROR='>>Error from GX'
          GOTO 9999
        ENDIF
      ELSE IF(ISEG(ISEGNUM).EQ.2) THEN
        CALL SOBVIS(ISEGNUM,gxVIS,ERR)
        IF(err.ne.0) THEN
          ERROR='>>Error from GX'
          GOTO 9999
        ENDIF
      ENDIF
      CALL SOBPIK(ISEGNUM,gxNOPIK,ERR)
      IF(err.ne.0) THEN
        ERROR='>>Error from GX'
        GOTO 9999
      ENDIF

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END


