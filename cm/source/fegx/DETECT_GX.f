      SUBROUTINE DETECT_GX(iw,ISEG,ISEGNUM,CLASS,ERROR,*)

C#### Subroutine: DETECT_GX
C###  Description:
C###    DETECT_GX changes ISEGNUM to CLASS='DETECTABLE' or 'UNDETECTABLE'

C**** Change ISEG(ISEGNUM) to 3 if detectable, 2 if not.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='DETECT_GX')
      INCLUDE 'gx.inc'
!     Parameter List
      INTEGER ISEG(*),ISEGNUM,iw
      CHARACTER CLASS*(*),ERROR*(*)
!     Local Variables
      INTEGER ERR

      CALL ENTERS(ROUTINENAME,*9999)

      CALL SETWIN(iw,ERR)
      IF(err.ne.0) THEN
        ERROR='>>Error from GX'
        GOTO 9999
      ENDIF
      IF(CLASS(1:10).EQ.'DETECTABLE') THEN
        CALL SOBPIK(ISEGNUM,gxPIK,ERR)
        IF(err.ne.0) THEN
          ERROR='>>Error from GX'
          GOTO 9999
        ENDIF
        ISEG(ISEGNUM)=3
      ELSE
        CALL SOBPIK(ISEGNUM,gxNOPIK,ERR)
        IF(err.ne.0) THEN
          ERROR='>>Error from GX'
          GOTO 9999
        ENDIF
        ISEG(ISEGNUM)=2
      ENDIF

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END


