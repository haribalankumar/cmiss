      SUBROUTINE SHRESI(ISEG,ISRESI,STRING,ERROR,*)

C#### Subroutine: SHRESI
C###  Description:
C###    SHRESI shows residual segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISRESI(NWM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),noiw,NTIW

      CALL ENTERS('SHRESI',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show residuals
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s). The all keyword specifies all
C###    currently defined workstations.
C###  Description:
C###    Make the residual segments visible on the specified
C###    workstations.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHRESI',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            IF(ISEG(ISRESI(iw)).EQ.1) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              CALL VISIB(iw,ISEG,ISRESI(iw),'VISIBLE',ERROR,*9999)
              CALL DAWK(iw,1,ERROR,*9999)
            ELSE IF(ISEG(ISRESI(iw)).EQ.0) THEN
              WRITE(OP_STRING,
     '          '('' >>Residuals are not defined on '',I1)') iw
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('SHRESI')
      RETURN
 9999 CALL ERRORS('SHRESI',ERROR)
      CALL EXITS('SHRESI')
      RETURN 1
      END


