      SUBROUTINE SHGRID(ISEG,ISGRID,STRING,ERROR,*)

C#### Subroutine: SHGRID
C###  Description:
C###    SHGRID shows grid segments.

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

      CALL ENTERS('SHGRID',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show grid
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to show the
C###    grid on.
C###  Description:
C###    Make the grid segments visible on the specified workstation.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHGRID',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            IF(ISEG(ISGRID(iw)).EQ.1) THEN
              CALL VISIB(iw,ISEG,ISGRID(iw),'VISIBLE',ERROR,*9999)
            ELSE IF(ISEG(ISGRID(iw)).EQ.0) THEN
              WRITE(OP_STRING,'('' >>Grid is not defined on '',I1)') iw
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO !noiw
      ENDIF

      CALL EXITS('SHGRID')
      RETURN
 9999 CALL ERRORS('SHGRID',ERROR)
      CALL EXITS('SHGRID')
      RETURN 1
      END


