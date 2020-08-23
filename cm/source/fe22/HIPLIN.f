      SUBROUTINE HIPLIN(ISEG,ISPLIN,STRING,ERROR,*)

C#### Subroutine: HIPLIN
C###  Description:
C###    HIPLIN hides polyline segment.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISPLIN(NWM,NGRSEGM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),noiw,noplin,NTIW

      CALL ENTERS('HIPLIN',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide polyline
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to hide the
C###    polyline on.
C###  Description:
C###    Hide polylines on specified workstations.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HIPLIN',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO noplin=1,NTPLIN
              IF(ISEG(ISPLIN(iw,noplin)).EQ.2) THEN
                CALL VISIB(iw,ISEG,ISPLIN(iw,noplin),'INVISIBLE',ERROR,
     '            *9999)
              ENDIF
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HIPLIN')
      RETURN
 9999 CALL ERRORS('HIPLIN',ERROR)
      CALL EXITS('HIPLIN')
      RETURN 1
      END


