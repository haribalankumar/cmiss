      SUBROUTINE CACROS(ISCROS,ISEG,STRING,ERROR,*)

C#### Subroutine: CACROS
C###  Description:
C###    CACROS cancels cross-section segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISCROS(NWM,NGRSEGM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),nocros,noiw,NTIW

      CALL ENTERS('CACROS',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel cross-section;s
C###  Parameter:      <on (all/WS#s)[all]>
C###    Specify the worksation (GX window) to cancel the
C###    cross-sections from.
C###  Description:
C###    Cancel crosssection segment on specified workstations.

        OP_STRING(1)=STRING(1:IEND)//';s'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CACROS',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO nocros=1,NTCROS
              IF(ISCROS(iw,nocros).GT.0) THEN
                CALL DELETE_SEGMENT(ISCROS(iw,nocros),ISEG,iw,ERROR,
     '            *9999)
              ENDIF
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
        NTCROS=0
      ENDIF

      CALL EXITS('CACROS')
      RETURN
 9999 CALL ERRORS('CACROS',ERROR)
      CALL EXITS('CACROS')
      RETURN 1
      END


