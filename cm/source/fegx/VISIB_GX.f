      SUBROUTINE VISIB_GX(ISEG,ISEGNUM,CLASS,ERROR,*)

C#### Subroutine: VISIB_GX
C###  Description:
C###    Changes ISEGNUM to CLASS='VISIBLE' or 'INVISIBLE'.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='VISIB_GX')
      INCLUDE 'gx.inc'
!     Parameter List
      INTEGER ISEG(*),ISEGNUM
      CHARACTER CLASS*(*),ERROR*(*)
!     Local Variables
      INTEGER ERR

      CALL ENTERS(ROUTINENAME,*9999)

      IF(CLASS(1:7).EQ.'VISIBLE') THEN
        CALL SOBVIS(ISEGNUM,gxVIS,ERR)
        IF(err.ne.0) THEN
          ERROR='>>Error from GX'
          GOTO 9999
        ENDIF
        ISEG(ISEGNUM)=2
      ELSE
        CALL SOBVIS(ISEGNUM,gxNOVIS,ERR)
        IF(err.ne.0) THEN
          ERROR='>>Error from GX'
          GOTO 9999
        ENDIF
        ISEG(ISEGNUM)=1
      ENDIF

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END
