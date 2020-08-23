      SUBROUTINE HIINCR(ISEG,ISINCR,STRING,ERROR,*)

C#### Subroutine: HIINCR
C###  Description:
C###    HIINCR hides increment segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISINCR(NWM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),noiw,NTIW

      CALL ENTERS('HIINCR',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide increments
C###  Description:
C###    Hide increments.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s).
C###    The all keyword specifies all currently defined workstations.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HIINCR',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            IF(ISEG(ISINCR(iw)).EQ.2) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              CALL VISIB(iw,ISEG,ISINCR(iw),'INVISIBLE',ERROR,*9999)
              CALL DAWK(iw,1,ERROR,*9999)
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HIINCR')
      RETURN
 9999 CALL ERRORS('HIINCR',ERROR)
      CALL EXITS('HIINCR')
      RETURN 1
      END


