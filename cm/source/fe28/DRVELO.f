      SUBROUTINE DRVELO(NEELEM,NELIST,STRING,ERROR,*)

C#### Subroutine: DRVELO
C###  Description:
C###    DRVELO draws velocity field segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,N3CO,noelem,nolist,nr
      LOGICAL CBBREV

      CALL ENTERS('DRVELO',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw velocity
C###  Parameter:    <in (all/ELEMENT#s[all]>
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Description:
C###    Draw velocity field segments in the specified elements using
C###    the specified colour.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<in (all/ELEMENT#s)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<rgb=RGB[black]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRVELO',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        IF(CBBREV(CO,'IN',1,noco+1,noco+1,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NE_R_M,NELIST(0),NELIST(1),ERROR,*9999)
        ELSE
          NELIST(0)=0
          DO nr=1,NRT
            DO noelem=NELIST(0)+1,NELIST(0)+NEELEM(0,nr)
              NELIST(noelem)=NEELEM(noelem,nr)
            ENDDO
            NELIST(0)=NELIST(0)+NEELEM(0,nr)
          ENDDO
        ENDIF
        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
C GMH 2/9/95 Unused          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
C GMH 2/9/95 Unused          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
        ENDIF

        IW=2*NJT-3
        NTVELO=NTVELO+1
        CALL ASSERT(NTVELO.LE.NRM,'>>NRM too small',ERROR,*9999)
        CALL ACWK(iw,1,ERROR,*9999)
        DO nolist=1,NELIST(0)
C GMH 2/9/95 Unused          ne=NELIST(nolist)
C         CALL SGVELO(INDEX)
        ENDDO
        CALL DAWK(iw,1,ERROR,*9999)
      ENDIF

      CALL EXITS('DRVELO')
      RETURN
 9999 CALL ERRORS('DRVELO',ERROR)
      CALL EXITS('DRVELO')
      RETURN 1
      END


