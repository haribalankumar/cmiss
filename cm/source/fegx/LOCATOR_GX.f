      SUBROUTINE LOCATOR_GX(INSTAT,XREF,XWC,YREF,YWC,ERROR,*)

C#### Subroutine: LOCATOR_GX
C###  Description:

C**** INSTAT is returned as 1 if locate is successful, 0 otherwise
C**** XREF,YREF are initial coords of echo
C**** XWC,YWC are returned world coords of located point

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='LOCATOR_GX')
      INCLUDE 'gx.inc'
!     Parameter List
      INTEGER INSTAT
      REAL XREF,XWC,YREF,YWC
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CANCEL,ERR
      REAL BNDRCT(4)

      CALL ENTERS(ROUTINENAME,*9999)

      BNDRCT(1)=0.0
      BNDRCT(2)=0.0
      BNDRCT(3)=0.0
      BNDRCT(4)=0.0
      CALL LOCATR(gxLOC_POINT,BNDRCT,XREF,YREF,XWC,YWC,CANCEL,ERR)
      IF(err.ne.0) THEN
        ERROR='>>Error from GX'
        GOTO 9999
      ENDIF
      IF(CANCEL.GT.0) THEN
        INSTAT=0
      ELSE
        INSTAT=1
      ENDIF

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END


