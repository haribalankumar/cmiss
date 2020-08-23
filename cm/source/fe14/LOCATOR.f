      SUBROUTINE LOCATOR(INSTAT,D_XREF,D_XWC,D_YREF,D_YWC,ERROR,*)

C#### Subroutine: LOCATOR
C###  Description:
C###    LOCATOR calls GKS locator.

C**** INIT specifies whether locator is to be initialised (1:yes,2:no)
C**** iw specifies workstation number (which is also transformation no)
C**** INSTAT is returned as 1 if locate is successful, 0 otherwise
C**** MODE is input mode - 'REQUEST','SAMPLE' or 'EVENT'
C**** NECHO is 0 on entry for no echo
C**** NECHO is 3 on entry for default   prompt/echo
C**** NECHO is 4 on entry for line      prompt/echo
C**** NECHO is 5 on entry for rectangle prompt/echo
C**** NECHO is 6 on entry for digital   prompt/echo
C**** XREF,YREF are initial coords of echo
C**** XWC,YWC are returned world coords of located point

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'disp00.cmn'
!     Parameter List
      INTEGER INSTAT
      REAL*8 D_XREF,D_XWC,D_YREF,D_YWC
      CHARACTER ERROR*(*)
!     Local Variables
      REAL XWC,YWC

      CALL ENTERS('LOCATOR',*9999)
      IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
        CALL LOCATOR_GX(INSTAT,REAL(D_XREF),XWC,REAL(D_YREF),YWC,
     '    ERROR,*9999)
      ENDIF
      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' XWC='',E11.3,'' YWC='',E11.3)') XWC,YWC
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      D_XWC=DBLE(XWC)
      D_YWC=DBLE(YWC)

      CALL EXITS('LOCATOR')
      RETURN

 9999 CALL ERRORS('LOCATOR',ERROR)
      CALL EXITS('LOCATOR')
      RETURN 1
      END


