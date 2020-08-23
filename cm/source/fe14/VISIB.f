      SUBROUTINE VISIB(iw,ISEG,ISEGNUM,CLASS,ERROR,*)

C#### Subroutine: VISIB
C###  Description:
C###    VISIB changes ISEGNUM to CLASS='VISIBLE' or 'INVISIBLE'.

C**** change ISEG(ISEGNUM) to 2 if visible, 1 if not.

      IMPLICIT NONE
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'disp00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISEGNUM,iw
      CHARACTER CLASS*(*),ERROR*(*)
!     Local variables

      CALL ENTERS('VISIB',*9999)

      IF(IWKT(iw).EQ.1) THEN      !GKS window
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          CALL VISIB_GX(ISEG,ISEGNUM,CLASS,ERROR,*9999)
        ENDIF
      ELSE IF(IWKT(iw).EQ.2) THEN !Phigs window
      ENDIF

      CALL EXITS('VISIB')
      RETURN
 9999 CALL ERRORS('VISIB',ERROR)
      CALL EXITS('VISIB')
      RETURN 1
      END
