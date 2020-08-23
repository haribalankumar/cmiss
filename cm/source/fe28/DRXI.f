      SUBROUTINE DRXI(STRING,ERROR,*)

C#### Subroutine: DRXI
C###  Description:
C###    DRXI draws nodal parameters.
C**** NPNODE(0,nr) is the total number of nodes in region nr.
C**** NPNODE(nonode,nr), nonode=1..NPNODE(0,nr) are the node numbers.
C**** NPT(nr) is the highest node number.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,IWK(6),
     '  N3CO,NTIW
      LOGICAL CBBREV

      CALL ENTERS('DRXI',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM draw xi
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to draw.
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Description:
C###

        OP_STRING(1)=BLANK(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<rgb=RGB[black]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRXI',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
C GMH 3/9/95 Unused          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
C GMH 3/9/95 Unused          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','GREEN')
        ENDIF
      ENDIF

      CALL EXITS('DRXI')
      RETURN
 9999 CALL ERRORS('DRXI',ERROR)
      CALL EXITS('DRXI')
      RETURN 1
      END
