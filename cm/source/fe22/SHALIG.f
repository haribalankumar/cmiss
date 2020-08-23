      SUBROUTINE SHALIG(ISALIG,ISEG,STRING,ERROR,*)

C#### Subroutine: SHALIG
C###  Description:
C###    SHALIG shows alignment segment.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISALIG(NWM),ISEG(*)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),noiw,NTIW

      CALL ENTERS('SHALIG',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show alignment
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation (GX window) to show the
C###    alignment segment on.
C###  Description:
C###     Make the alignment segment visible.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHALIG',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            IF(ISEG(ISALIG(iw)).EQ.1) THEN
              CALL VISIB(iw,ISEG,ISALIG(iw),'VISIBLE',ERROR,*9999)
            ELSE IF(ISEG(ISALIG(iw)).EQ.0) THEN
              WRITE(OP_STRING,'('' >>Alignment not '
     '          //'defined on '',I1)') iw
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('SHALIG')
      RETURN
 9999 CALL ERRORS('SHALIG',ERROR)
      CALL EXITS('SHALIG')
      RETURN 1
      END


