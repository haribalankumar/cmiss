      SUBROUTINE DRAXES(ISAXES,ISEG,CSEG,STRING,ERROR,*)

C#### Subroutine: DRAXES
C###  Description:
C###    DRAXES draws axes segment.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'axes00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISAXES(NWM),ISEG(*)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,INDEX,INDEX_POLYLINE,IWK(6),iw,N3CO,noiw,NTIW
      LOGICAL CBBREV

      CALL ENTERS('DRAXES',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw axes
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the window number on which to draw the axes.
C###    The default is to draw axes on all windows.
C###  Parameter:    <x=X_TIC#[0]>
C###    Draw ticks on x axis at specified points (separated by commas)
C###  Parameter:    <y=Y_TIC#[0]>
C###    Draw ticks on y axis at specified points (separated by commas)
C###  Parameter:    <z=Z_TIC#[0]>
C###    Draw ticks on z axis at specified points (separated by commas)
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Description:
C###    Draws a set of axes on the specified workstation, of the
C###    colour given by the rgb value. x, y, and z can be used to
C###    draw tick marks on the axes at the specified points on the
C###    respective axis.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<x=X_TIC#[0]>'
        OP_STRING(4)=BLANK(1:15)//'<y=Y_TIC#[0]>'
        OP_STRING(5)=BLANK(1:15)//'<z=Z_TIC#[0]>'
        OP_STRING(6)=BLANK(1:15)//'<rgb=RGB[black]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRAXES',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'X',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSRL(CO(N3CO+1),50,NTX,XLIST,ERROR,*9999)
        ELSE
          NTX=0
        ENDIF
        IF(CBBREV(CO,'Y',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSRL(CO(N3CO+1),50,NTY,YLIST,ERROR,*9999)
        ELSE
          NTY=0
        ENDIF
        IF(CBBREV(CO,'Z',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSRL(CO(N3CO+1),50,NTZ,ZLIST,ERROR,*9999)
        ELSE
          NTZ=0
        ENDIF

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
        ENDIF

        DO noiw=1,NTIW
          IW=IWK(noiw)
C CPB 10/9/92 Changed workstation update mode
          CALL ACWK(iw,1,ERROR,*9999)
          CALL SGAXES(INDEX,ISAXES(iw),ISEG,iw,CSEG,ERROR,*9999)
          CALL DAWK(iw,1,ERROR,*9999)
        ENDDO

      ENDIF

      CALL EXITS('DRAXES')
      RETURN
 9999 CALL ERRORS('DRAXES',ERROR)
      CALL EXITS('DRAXES')
      RETURN 1
      END


