      SUBROUTINE CAHIST(ISEG,ISHIST,STRING,ERROR,*)

C#### Subroutine: CAHIST
C###  Description:
C###    CAHIST cancels history segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISHIST(0:NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,np

      CALL ENTERS('CAHIST',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel history
C###  Description:
C###    Cancel history display.

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAHIST',ERROR,*9999)
      ELSE
        iw=1 !AJP Hack fix.  What should iw be?
        IF(IWKS(iw).GT.0) THEN
          CALL ACWK(10,1,ERROR,*9999)
          DO np=0,NPT(1)
            IF(ISHIST(np).GT.0) THEN
              CALL DELETE_SEGMENT(ISHIST(np),ISEG,iw,ERROR,*9999)
            ENDIF
          ENDDO
          CALL DAWK(10,1,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('CAHIST')
      RETURN
 9999 CALL ERRORS('CAHIST',ERROR)
      CALL EXITS('CAHIST')
      RETURN 1
      END


