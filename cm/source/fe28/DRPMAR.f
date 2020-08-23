      SUBROUTINE DRPMAR(ISEG,ISPMAR,CSEG,STRING,ERROR,*)

C#### Subroutine: DRPMAR
C###  Description:
C###    DRPMAR draws polymarker segment.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISPMAR(NWM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,INDEX,INDEX_POLYMARKER,iw,IWK(6),
     '  N3CO,noiw,NTIW,NTPTS,NTX,NTY
      REAL*8 XLIST(100),YLIST(100)
      LOGICAL CBBREV

      CALL ENTERS('DRPMAR',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw polymarker x=X_LIST y=Y_LIST
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to draw the
C###    points on.
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Description:
C###    Draws a polymarker.

        OP_STRING(1)=STRING(1:IEND)//' x=X_LIST y=Y_LIST'
        OP_STRING(2)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<rgb=RGB[black]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRPMAR',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        IF(CBBREV(CO,'X',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSRL(CO(N3CO+1),100,NTX,XLIST,ERROR,*9999)
        ENDIF
        IF(CBBREV(CO,'Y',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSRL(CO(N3CO+1),100,NTY,YLIST,ERROR,*9999)
        ENDIF
        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYMARKER(0,'PLUS','SIZE1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYMARKER(0,'PLUS','SIZE1','BLACK')
        ENDIF

        NTPTS=NTX
        IF(NTY.LT.NTPTS) NTPTS=NTY
        DO noiw=1,NTIW
          IW=IWK(noiw)
          CALL ACWK(iw,1,ERROR,*9999)
          CALL SGPMAR(INDEX,ISEG,ISPMAR(iw),iw,NTPTS,CSEG,XLIST,YLIST,
     '      ERROR,*9999)
          CALL DAWK(iw,1,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('DRPMAR')
      RETURN
 9999 CALL ERRORS('DRPMAR',ERROR)
      CALL EXITS('DRPMAR')
      RETURN 1
      END


