      SUBROUTINE CACLOC(ISCLOC,ISEG,STRING,ERROR,*)

C#### Subroutine: CACLOC
C###  Description:
C###    CACLOC cancels clock segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISCLOC(NWM),ISEG(*)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw

      CALL ENTERS('CACLOC',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel clock;s
C###  Description:
C###    Cancels clock definition.

        OP_STRING(1)=STRING(1:IEND)//';s'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CACLOC',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        DO iw=1,2*NJT-3
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            IF(ISCLOC(iw).GT.0) THEN
              CALL DELETE_SEGMENT(ISCLOC(iw),ISEG,iw,ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('CACLOC')
      RETURN
 9999 CALL ERRORS('CACLOC',ERROR)
      CALL EXITS('CACLOC')
      RETURN 1
      END


