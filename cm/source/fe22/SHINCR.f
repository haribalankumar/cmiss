      SUBROUTINE SHINCR(ISEG,ISINCR,STRING,ERROR,*)

C#### Subroutine: SHINCR
C###  Description:
C###    SHINCR shows increment segments.

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

      CALL ENTERS('SHINCR',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show increments
C###  Description:
C###    Make the increment segments visible.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s).
C###    The all keyword specifies all currently defined workstations.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHINCR',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            IF(ISEG(ISINCR(iw)).EQ.1) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              CALL VISIB(iw,ISEG,ISINCR(iw),'VISIBLE',ERROR,*9999)
              CALL DAWK(iw,1,ERROR,*9999)
            ELSE IF(ISEG(ISINCR(iw)).EQ.0) THEN
              WRITE(OP_STRING,
     '          '('' >>Increment is not defined on '',I1)') iw
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('SHINCR')
      RETURN
 9999 CALL ERRORS('SHINCR',ERROR)
      CALL EXITS('SHINCR')
      RETURN 1
      END


