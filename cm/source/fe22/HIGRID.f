      SUBROUTINE HIGRID(ISEG,ISGRID,STRING,ERROR,*)

C#### Subroutine: HIGRID
C###  Description:
C###    HIGRID hides grid segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISGRID(NWM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),noiw,NTIW

      CALL ENTERS('HIGRID',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide grid
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to hide the
C###    grid on.
C###  Description:
C###    Hide grid on specified workstations.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HIGRID',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            IF(ISEG(ISGRID(iw)).EQ.2) THEN
              CALL VISIB(iw,ISEG,ISGRID(iw),'INVISIBLE',ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HIGRID')
      RETURN
 9999 CALL ERRORS('HIGRID',ERROR)
      CALL EXITS('HIGRID')
      RETURN 1
      END


