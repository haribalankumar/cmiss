      SUBROUTINE CAPLIN(ISEG,ISPLIN,STRING,ERROR,*)

C#### Subroutine: CAPLIN
C###  Description:
C###    CAPLIN cancels polyline segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'plin00.cmn'
!     Parameter List
      INTEGER ISPLIN(NWM,NGRSEGM),ISEG(*)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),nj,noiw,noplin,nopts,NTIW
      LOGICAL ABBREV,SEGME

      CALL ENTERS('CAPLIN',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel polyline;s
C###  Parameter:      <on WS#[1]>
C###    Specify the workstation (GX window) to cancel the
C###    polyline on.
C###  Description:
C###    Cancel polyline segment on specified workstation.

        OP_STRING(1)=STRING(1:IEND) //';s'
        OP_STRING(2)=BLANK(1:15) //'<on WS#[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAPLIN',ERROR,*9999)
      ELSE
        CALL CHECKQ(' S',noco,1,CO,COQU,STRING,*1)
        SEGME=.FALSE.
        IF(ABBREV(COQU(noco,1),'S',1)) THEN
          SEGME=.TRUE.
          CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        ENDIF

        IF(SEGME) THEN
          DO noiw=1,NTIW
            iw=IWK(noiw)
            IF(IWKS(iw).GT.0) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              DO noplin=1,NTPLIN
                IF(ISPLIN(iw,noplin).GT.0) THEN
                  CALL DELETE_SEGMENT(ISPLIN(iw,noplin),ISEG,iw,ERROR,
     '              *9999)
                ENDIF
              ENDDO
              CALL DAWK(iw,1,ERROR,*9999)
            ENDIF
          ENDDO
          NTPLIN=0

        ELSE
          DO noplin=1,NT_PLIN
            INDEX_PLIN_TYPE(noplin)=0
            NT_PLIN_SECTIONS(noplin)=0
            DO nopts=1,20
              DO nj=1,3
                PLIN_DATA(nj,nopts,noplin)=0.0d0
              ENDDO
            ENDDO
          ENDDO
          NT_PLIN=0
        ENDIF
      ENDIF

      CALL EXITS('CAPLIN')
      RETURN
 9999 CALL ERRORS('CAPLIN',ERROR)
      CALL EXITS('CAPLIN')
      RETURN 1
      END


