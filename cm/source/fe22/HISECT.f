      SUBROUTINE HISECT(ISEG,ISSECT,STRING,ERROR,*)

C#### Subroutine: HISECT
C###  Description:
C###    HISECT hides section segments.

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

      CALL ENTERS('HISECT',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide section
C###  Description:
C###    Hide section plot.

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HISECT',ERROR,*9999)
      ELSE
C LKC 2-MAY-1998 use iw
C        IF(IWKS(11).GT.0) THEN
C          CALL ACWK(11,1,ERROR,*9999)
        iw=11
        IF(IWKS(iw).GT.0) THEN
          CALL ACWK(iw,1,ERROR,*9999)
          DO nosect=1,NTSECT
            IF(ISEG(ISSECT(nosect)).EQ.2) THEN
              CALL VISIB(iw,ISEG,ISSECT(nosect),'INVISIBLE',ERROR,*9999)
            ENDIF
          ENDDO
C LKC 2-MAY-1998 use iw
C          CALL DAWK(11,1,ERROR,*9999)
          CALL DAWK(iw,1,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('HISECT')
      RETURN
 9999 CALL ERRORS('HISECT',ERROR)
      CALL EXITS('HISECT')
      RETURN 1
      END


