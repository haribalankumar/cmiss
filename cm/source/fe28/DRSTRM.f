      SUBROUTINE DRSTRM(ISEG,ISSTRM,NEELEM,NELIST,
     '  CSEG,STRING,ERROR,*)

C#### Subroutine: DRSTRM
C###  Description:
C###    DRSTRM draws streamline segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISSTRM(NEM,NGRSEGM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,INDEX,INDEX_POLYLINE,iw,N3CO,ne,noelem,nolist,nr
      LOGICAL CBBREV

      CALL ENTERS('DRSTRM',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw streamlines
C###  Parameter:    <in (all/ELEMENT#s[all]>
C###    Specify the element(s) in which to draw the streamlines. The
C###    "all" keyword specifies all currently defined elements.
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Description:
C###    Draws streamline segments in the specified elements using the
C###    specified colour

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<in (all/ELEMENT#s)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<rgb=RGB[black]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRSTRM',ERROR,*9999)
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
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
        ENDIF

        CALL STRING_TRIM(FILE03,IBEG,IEND)
 201    CALL OPENF(14,'DISK',FILE03(IBEG:IEND)//'.strmline','OLD',
     '    'DIRECT','FORMATTED',132,ERROR,*203)
 203    IF(ERROR(1:10).EQ.'Iostat= 30') THEN
          WRITE(OP_STRING,'('' Streamline file is locked'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          GO TO 201
        ELSE IF(ERROR(2:2).NE.' ') THEN
          GO TO 9999
        ENDIF

        IW=2*NJT-3
        IF(ADD) THEN
          NTSTRM=NTSTRM+1
        ELSE IF(NTSTRM.EQ.0) THEN
          NTSTRM=1
        ENDIF
        CALL ASSERT(NTSTRM.LE.NRM,'>>NRM too small',ERROR,*9999)
        CALL ACWK(iw,1,ERROR,*9999)
        DO nolist=1,NELIST(0)
          ne=NELIST(nolist)
          CALL SGSTRM(INDEX,ISEG,ISSTRM(ne,NTSTRM),CSEG,ERROR,*9999)
        ENDDO
        CALL DAWK(iw,1,ERROR,*9999)
        CALL CLOSEF(14,ERROR,*9999)
      ENDIF

      CALL EXITS('DRSTRM')
      RETURN
 9999 CALL ERRORS('DRSTRM',ERROR)
      CALL EXITS('DRSTRM')
      RETURN 1
      END


