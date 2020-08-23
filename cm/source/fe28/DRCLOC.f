      SUBROUTINE DRCLOC(ISCLOC,ISEG,CSEG,STRING,ERROR,*)

C#### Subroutine: DRCLOC
C###  Description:
C###    DRCLOC draws clock.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISCLOC(NWM),ISEG(*)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,iw,IWK(6),INDEX,
     '  INDEX_POLYLINE,N3CO,NTIW
      LOGICAL CBBREV

      CALL ENTERS('DRCLOC',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM draw clock
C###  Parameter:    <on WS_ID#[1]>
C###    Specify the worksation (GX window) to draw the
C###    clock on.
C###  Parameter:    <rgb=RGB[red]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Description:
C###    Draws a clock on the specified workstation.  The clock shows
C###    the proportion of the way that the program is through the
C###    solution of the problem.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<on WS_ID#[1]>'
        OP_STRING(3)=BLANK(1:15)//'<rgb=RGB[red]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRCLOC',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IW=IWK(1)
        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','RED')
        ENDIF
        CALL ACWK(iw,1,ERROR,*9999)
        CALL SGCLOC(INDEX,ISCLOC(iw),ISEG,iw,CSEG,0.0D0,ERROR,*9999)
        CALL DAWK(iw,1,ERROR,*9999)
      ENDIF

      CALL EXITS('DRCLOC')
      RETURN
 9999 CALL ERRORS('DRCLOC',ERROR)
      CALL EXITS('DRCLOC')
      RETURN 1
      END


