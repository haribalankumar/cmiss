      SUBROUTINE CAREAC(ISEG,ISREAC,STRING,ERROR,*)

C#### Subroutine: CAREAC
C###  Description:
C###    CAREAC cancels reaction segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISREAC(NWM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),noiw,NTIW

      CALL ENTERS('CAREAC',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel reactions;s
C###  Parameter:      <on (all/WS#s)[all]>
C###    Specify the worksation(s). The all keyword specifies all
C###    currently defined workstations.
C###  Description:
C###    Cancel reaction segment on specified workstations.

        OP_STRING(1)=STRING(1:IEND)//';s'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAREAC',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            IF(ISREAC(iw).GT.0) THEN
              CALL DELETE_SEGMENT(ISREAC(iw),ISEG,iw,ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('CAREAC')
      RETURN
 9999 CALL ERRORS('CAREAC',ERROR)
      CALL EXITS('CAREAC')
      RETURN 1
      END


