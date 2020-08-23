      SUBROUTINE CAPROF(ISEG,ISPROF,STRING,ERROR,*)

C#### Subroutine: CAPROF
C###  Description:
C###    CAPROF cancels profile segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
!     Parameter List
      INTEGER ISEG(*),ISPROF(2)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,noprof

      CALL ENTERS('CAPROF',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel profile
C###  Description:
C###    Cancel profile display.

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAPROF',ERROR,*9999)
      ELSE
        DO noprof=1,3
          iw=30+noprof
          CALL ACWK(iw,1,ERROR,*9999)
          IF(ISPROF(noprof).GT.0) THEN
            CALL DELETE_SEGMENT(ISPROF(noprof),ISEG,iw,ERROR,*9999)
          ENDIF
          CALL DAWK(iw,1,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('CAPROF')
      RETURN
 9999 CALL ERRORS('CAPROF',ERROR)
      CALL EXITS('CAPROF')
      RETURN 1
      END


