      SUBROUTINE DELETE_SEGMENT(ISEGNUM,ISEG,iw,ERROR,*)

C#### Subroutine: DELETE_SEGMENT
C###  Description:
C###    DELETE_SEGMENT deletes graphics segment ISEGNUM.

      IMPLICIT NONE
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'disp00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISEGNUM,iw
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER ERR

      CALL ENTERS('DELETE_SEGMENT',*9999)
      IF(IWKT(iw).EQ.1) THEN      !GKS
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          CALL DLOBJT(ISEGNUM,ERR)
          IF(err.ne.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
        ENDIF
      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
      ENDIF
      ISEG(ISEGNUM)=0
      ISEGNUM=0

      CALL EXITS('DELETE_SEGMENT')
      RETURN
 9999 CALL ERRORS('DELETE_SEGMENT',ERROR)
      CALL EXITS('DELETE_SEGMENT')
      RETURN 1
      END


