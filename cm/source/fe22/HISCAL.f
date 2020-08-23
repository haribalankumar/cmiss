      SUBROUTINE HISCAL(ISEG,ISSCAL,STRING,ERROR,*)

C#### Subroutine: HISCAL
C###  Description:
C###    HISCAL hides scale segment.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISSCAL(NWM,NGRSEGM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,iw,IWK(6),N3CO,noiw,NOSCAL,NTIW
      LOGICAL CBBREV

      CALL ENTERS('HISCAL',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide scale
C###  Parameter:    <at (last/SCALE_INDEX#)[last]>
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s). The all keyword specifies all
C###    currently defined workstations.
C###  Description:
C###    Hide scale on specified workstations.

        OP_STRING(1)=STRING(1:IEND) //' <at (last/SCALE_INDEX#)[last]>'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HISCAL',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'AT',2,noco+1,noco+3,N3CO)) THEN
          NOSCAL=IFROMC(CO(N3CO+1))
        ELSE
          NOSCAL=NTSCAL
        ENDIF

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            IF(ISEG(ISSCAL(iw,NOSCAL)).EQ.2) THEN
              CALL VISIB(iw,ISEG,ISSCAL(iw,NOSCAL),'INVISIBLE',ERROR,
     '          *9999)
            ELSE IF(ISEG(ISSCAL(iw,NOSCAL)).EQ.0) THEN
              WRITE(OP_STRING,'('' >>Scale at '',I3,'
     '          //''' is not defined on '',I1)') NOSCAL,iw
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HISCAL')
      RETURN
 9999 CALL ERRORS('HISCAL',ERROR)
      CALL EXITS('HISCAL')
      RETURN 1
      END


