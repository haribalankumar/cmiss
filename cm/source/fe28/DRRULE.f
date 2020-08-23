      SUBROUTINE DRRULE(ISEG,ISRULE,CSEG,STRING,ERROR,*)

C#### Subroutine: DRRULE
C###  Description:
C###    <HTML>
C###    DRRULE draws ruled lines for following coordinate systems:
C###    <PRE>
C###    ITYP10(1)=1: Rectangular cartesian
C###        "     2: Cylindrical polar
C###        "     3: Spherical polar
C###        "     4: Prolate spheroidal
C###        "     5: Oblate spheroidal
C###    </PRE> </HTML>

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISRULE(NWM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,INDEX,INDEX_POLYLINE,iw,IWK(6),
     '  N3CO,noiw,NTIW
      LOGICAL CBBREV

      CALL ENTERS('DRRULE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM draw rule
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s). The all keyword specifies all
C###    currently defined workstations.
C###  Parameter:    <rgb=RGB[cyan]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Description:
C###    Draws a ruled grid, for most common coordinate systems, on the
C###    specified workstation, in the specified colour.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<rgb=RGB[cyan]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRRULE',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'DOTTED','WIDTH1','BLACK')
        ENDIF

        DO noiw=1,NTIW
          IW=IWK(noiw)
          CALL ACWK(iw,1,ERROR,*9999)
          CALL SGRULE(INDEX,ISEG,ISRULE(iw),iw,CSEG,ERROR,*9999)
          CALL DAWK(iw,1,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('DRRULE')
      RETURN
 9999 CALL ERRORS('DRRULE',ERROR)
      CALL EXITS('DRRULE')
      RETURN 1
      END


