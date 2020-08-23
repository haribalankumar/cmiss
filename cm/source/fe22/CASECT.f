      SUBROUTINE CASECT(ISEG,ISSECT,STRING,ERROR,*)

C#### Subroutine: CASECT
C###  Description:
C###    CASECT cancels section segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISSECT(NGRSEGM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,nosect

      CALL ENTERS('CASECT',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel section
C###  Description:
C###    Cancel section segment.

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CASECT',ERROR,*9999)
      ELSE
C LKC 2-MAY-1998 Using iw
C
C        IF(IWKS(11).GT.0) THEN
C          CALL ACWK(11,1,ERROR,*9999)
        iw=11
        IF(IWKS(iw).GT.0) THEN
          CALL ACWK(iw,1,ERROR,*9999)
          DO nosect=1,NTSECT
            IF(ISSECT(nosect).GT.0) THEN
              CALL DELETE_SEGMENT(ISSECT(nosect),ISEG,iw,ERROR,*9999)
            ENDIF
          ENDDO
          NTSECT=0
C LKC 2-MAY-1998 Using iw
C          CALL DAWK(11,1,ERROR,*9999)
          CALL DAWK(iw,1,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('CASECT')
      RETURN
 9999 CALL ERRORS('CASECT',ERROR)
      CALL EXITS('CASECT')
      RETURN 1
      END


