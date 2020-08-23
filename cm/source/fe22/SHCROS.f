      SUBROUTINE SHCROS(ISCROS,ISEG,STRING,ERROR,*)

C#### Subroutine: SHCROS
C###  Description:
C###    SHCROS shows cross-section segments.

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

      CALL ENTERS('SHCROS',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show cross-section
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation (GX window) to show the
C###    cross-section segments on.
C###  Description:
C###    Make the cross-section segments visible.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHCROS',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO nocros=1,NTCROS
              IF(ISEG(ISCROS(iw,nocros)).EQ.1) THEN
                CALL VISIB(iw,ISEG,ISCROS(iw,nocros),'VISIBLE',ERROR,
     '            *9999)
              ELSE IF(ISEG(ISCROS(iw,nocros)).EQ.0) THEN
                WRITE(OP_STRING,'('' >>Isochrone '',I3,'
     '            //''' is not defined on '',I1)') nocros,iw
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('SHCROS')
      RETURN
 9999 CALL ERRORS('SHCROS',ERROR)
      CALL EXITS('SHCROS')
      RETURN 1
      END


