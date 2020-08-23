      SUBROUTINE GRPLIN(STRING,ERROR,*)

C#### Subroutine: GRPLIN
C###  Description:
C###    GRPLIN groups polylines.
C**** NTGRPL is number of polyline groups currently defined.
C**** LAGRPL(nogrpl) is label given to group number NOGRPL.
C**** LIGRPL(0,nogrpl) is number in list for group number NOGRPL.
C**** LIGRPL(1..,nogrpl) is list for group number NOGRPL.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'plin00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG2,IEND,IEND2,N3CO,no_plin
      CHARACTER CHAR2*2
      LOGICAL CBBREV

      CALL ENTERS('GRPLIN',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
C        CHAR2=CFROMI(NTGRPL+1,'(I2)')
        WRITE(CHAR2,'(I2)') NTGRPL+1
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM group polylines
C###  Parameter:     <(POLYLINE_INDEX#s/all)[all]>
C###    Specifies the polyline indices to include.
C###    The 'all' command prompts all currently defined polylines
C###    to be included.
C###  Parameter:     <as LABEL[poly_1]>
C###    Specifies the name for the polyline group.
C###  Description:
C###    Creates a polyline group.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<(POLYLINE_INDEX#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)
     '    //'<as LABEL[poly_'//CHAR2(IBEG2:IEND2)//']>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','GRPLIN',ERROR,*9999)
      ELSE
        NTGRPL=NTGRPL+1
        CALL ASSERT(NTGRPL.LE.GRPL_MAXGRP,
     '    '>>Increase array sizes in grou00.cmn',ERROR,*9999)
        IF(NTCO.EQ.4.OR.NTCO.EQ.6) THEN
          CALL PARSIL(CO(noco+1),NT_PLIN,LIGRPL(0,NTGRPL),
     '  LIGRPL(1,NTGRPL),
     '      ERROR,*9999)
        ELSE
          CALL ASSERT(NT_PLIN.LE.GRPL_MAXDATA,
     '      '>>Increase array sizes in grou00.cmn',ERROR,*9999)
          LIGRPL(0,NTGRPL)=NT_PLIN
          DO no_plin=1,NT_PLIN
            LIGRPL(no_plin,NTGRPL)=no_plin
          ENDDO
        ENDIF
        IF(CBBREV(CO,'AS',1,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          LAGRPL(NTGRPL)=CO(N3CO+1)(IBEG:IEND) !is polyline group label
        ELSE
C          CHAR2=CFROMI(NTGRPL,'(I2)')
          WRITE(CHAR2,'(I2)') NTGRPL
          CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
          LAGRPL(NTGRPL)='poly_'//CHAR2(IBEG2:IEND2) !new label
        ENDIF
      ENDIF

      CALL EXITS('GRPLIN')
      RETURN
 9999 CALL ERRORS('GRPLIN',ERROR)
      CALL EXITS('GRPLIN')
      RETURN 1
      END


