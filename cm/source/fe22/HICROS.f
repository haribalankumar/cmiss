      SUBROUTINE HICROS(ISCROS,ISEG,STRING,ERROR,*)

C#### Subroutine: HICROS
C###  Description:
C###    HICROS hides cross-section segments.

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

      CALL ENTERS('HICROS',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide cross-section
C###  Parameter:   <on (all/WS#s)[all]>
C###    Specify the worksation (GX window) to hide the
C###    cross-sections from.
C###  Description:
C###    Hide cross-sections on specified workstations.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HICROS',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO nocros=1,NTCROS
              IF(ISEG(ISCROS(iw,nocros)).EQ.2) THEN
                CALL VISIB(iw,ISEG,ISCROS(iw,nocros),'INVISIBLE',
     '            ERROR,*9999)
              ENDIF
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HICROS')
      RETURN
 9999 CALL ERRORS('HICROS',ERROR)
      CALL EXITS('HICROS')
      RETURN 1
      END


