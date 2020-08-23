      SUBROUTINE SHLINE(ISEG,ISLINE,ISLINO,STRING,ERROR,*)

C#### Subroutine: SHLINE
C###  Description:
C###    SHLINE shows line segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISLINE(NWM,2*NGRSEGM),ISLINO(NWM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,iw,IWK(6),N3CO,noiw,noline,NTIW
      LOGICAL ABBREV,CBBREV

      CALL ENTERS('SHLINE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show lines
C###  Description:
C###    Make the specified line segments visible.
C###  Parameter:    <at (last/LINE_INDEX)[last]>
C###    Specify either the 'last' line drawn or the line segment number.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s).
C###    The all keyword specifies all currently defined workstations.

        OP_STRING(1)=STRING(1:IEND) //' <at (last/LINE_INDEX)[last]>'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM show line numbers
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to show the
C###    line numbers on.
C###  Description:
C###    Make the line number segments visible on the specified
C###    workstation.

        OP_STRING(1)=STRING(1:IEND) //' numbers'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHLINE',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'AT',2,noco+1,noco+3,N3CO)) THEN
          noline=IFROMC(CO(N3CO+1))
        ELSE
          noline=NTLINE
        ENDIF

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            IF(ISEG(ISLINE(iw,noline)).EQ.1) THEN
              CALL VISIB(iw,ISEG,ISLINE(iw,noline),'VISIBLE',ERROR,
     '          *9999)
            ELSE IF(ISEG(ISLINE(iw,noline)).EQ.0) THEN
              WRITE(OP_STRING,'('' >>Line at '',I3,'
     '          //''' is not defined on '',I1)') noline,iw
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
        IF(ABBREV(CO(noco+1),'NUMBERS',1)) THEN
          DO noiw=1,NTIW
            iw=IWK(noiw)
            CALL ACWK(iw,1,ERROR,*9999)
            IF(ISEG(ISLINO(iw)).EQ.1) THEN
              CALL VISIB(iw,ISEG,ISLINO(iw),'VISIBLE',ERROR,*9999)
            ELSE IF(ISEG(ISLINO(iw)).EQ.0) THEN
              WRITE(OP_STRING,'('' >>Line numbers are not '
     '          //'defined on '',I1)') iw
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('SHLINE')
      RETURN
 9999 CALL ERRORS('SHLINE',ERROR)
      CALL EXITS('SHLINE')
      RETURN 1
      END


