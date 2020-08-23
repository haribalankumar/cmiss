      SUBROUTINE CAAXES(ISAXES,ISEG,STRING,ERROR,*)

C#### Subroutine: CAAXES
C###  Description:
C###    CAAXES cancels axes segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISAXES(NWM),ISEG(*)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw

      CALL ENTERS('CAAXES',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel axes;s
C###  Description:
C###    Cancel axes segment.

        OP_STRING(1)=STRING(1:IEND)//';s'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAAXES',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        DO iw=1,2*NJT-3
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            IF(ISAXES(iw).GT.0) THEN
              CALL DELETE_SEGMENT(ISAXES(iw),ISEG,iw,ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('CAAXES')
      RETURN
 9999 CALL ERRORS('CAAXES',ERROR)
      CALL EXITS('CAAXES')
      RETURN 1
      END


