      SUBROUTINE HICLOC(ISCLOC,ISEG,STRING,ERROR,*)

C#### Subroutine: HICLOC
C###  Description:
C###    HICLOC hides clock segment.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISCLOC(NWM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),noiw,NTIW

      CALL ENTERS('HICLOC',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide clock
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation (GX window) to hide the
C###    clock from.
C###  Description:
C###    Hide the clock on a specified workstation.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HICLOC',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            IF(ISEG(ISCLOC(iw)).EQ.2) THEN
              CALL VISIB(iw,ISEG,ISCLOC(iw),'INVISIBLE',ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HICLOC')
      RETURN
 9999 CALL ERRORS('HICLOC',ERROR)
      CALL EXITS('HICLOC')
      RETURN 1
      END


