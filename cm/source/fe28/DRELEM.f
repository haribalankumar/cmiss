      SUBROUTINE DRELEM(IBT,IDO,INP,ISEG,ISELNO,ISERR,MXI,
     '  NAN,NBJ,NEELEM,NELIST,NKJE,NLL,
     '  NPF,NPL,NPNE,NRE,NRLIST,NVJE,DL,NEERR,SE,XA,XE,XP,XW,
     '  CSEG,STRING,ERROR,*)

C#### Subroutine: DRELEM
C###  Description:
C###    DRELEM draws element parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISEG(*),ISELNO(NWM,NEM),ISERR(NWM,NEM),MXI(2,NEM),
     '  NAN(NIM,NAM,NBFM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NPF(9,NFM),NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),
     '  NRE(NEM),NRLIST(0:NRM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 DL(3,NLM),NEERR(NEM,3),SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XP(NKM,NVM,NJM,NPM),XW(NJM,NUM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,INDEX,
     '  INDEX_POLYLINE,iw,IWK(6),
     '  N3CO,nb,ne,nj,njj,noiw,nolist,no_nrlist,
     '  nr,NTIW
      REAL*8 PXI,X(3),XI(3),Z(3)
      LOGICAL ALL_REGIONS,BOUNDARY,BOX,CBBREV,ERR

      CALL ENTERS('DRELEM',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM draw elements
C###  Description:
C###    Draws the specified elements on the specified workstation,
C###    draws only the element numbers. Elements are drawn in the
C###    colour specified by the rgb value.
C###  Parameter:    <elements (all/#s/GROUP)[all]>
C###    Specify the elements which are to be drawn. The default is
C###    to draw all elements in the current region. A list of numbers
C###    or a group name or the 'all' option is allowed.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the number of the window to which the elements will
C###    be drawn.
C###  Parameter:    <rgb=RGB[blue]>
C###    This parameter specifies what colour the boundaries are
C###    to be drawn if the boundary option is selected. The colours
C###    allowed are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY
C###  Parameter:    <boundary>
C###    This parameter draws the boundaries of the elements as well
C###    as the element numbers.
C###  Parameter:    <box>
C###    This parameter puts a box around the element number.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<elements (all/#s/GROUP)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<rgb=RGB[blue]>'
        OP_STRING(5)=BLANK(1:15)//'<boundary>'
        OP_STRING(6)=BLANK(1:15)//'<box>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw elements error
C###  Description:
C###    Draws the relative energy norm error in elements on the
C###    specified workstation. Errors are drawn in the
C###    colour specified by the rgb value.
C###  Parameter:    <elements (all/#s/GROUP)[all]>
C###    Specify the elements which are to be drawn. The default is
C###    to draw all elements in the current region. A list of numbers
C###    or a group name or the 'all' option is allowed.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the number of the window to which the elements will
C###    be drawn.
C###  Parameter:    <rgb=RGB[blue]>
C###    This parameter specifies what colour the errors are
C###    to be drawn. The colours allowed are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.

        OP_STRING(1)=STRING(1:IEND)//' error'
        OP_STRING(2)=BLANK(1:15)//'<elements (all/#s/GROUP)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<rgb=RGB[blue]>'
        OP_STRING(5)=BLANK(1:15)//'<box>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRELEM',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)

        CALL ASSERT(NBT.GT.0,'>>No basis functions are defined',
     '    ERROR,*9999)
        CALL ASSERT(NET(NRLIST(1)).GT.0,'>>No elements defined',
     '    ERROR, *9999)

        IF(CBBREV(CO,'BOUNDARY',1,noco+1,NTCO,N3CO)) THEN
          BOUNDARY=.TRUE.
        ELSE
          BOUNDARY=.FALSE.
        ENDIF

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLUE')
        ENDIF

        IF(CBBREV(CO,'BOX',1,noco+1,NTCO,N3CO)) THEN
          BOX=.TRUE.
        ELSE
          BOX=.FALSE.
        ENDIF

        IF(CBBREV(CO,'ERROR',2,noco+1,NTCO,N3CO)) THEN
          ERR=.TRUE.
        ELSE
          ERR=.FALSE.
        ENDIF


        IF(.NOT.ERR) THEN
          XI(1)=0.5d0
          XI(2)=0.5d0
          XI(3)=0.5d0
        ELSE
          XI(1)=0.5d0
          XI(2)=0.25d0
          XI(3)=0.5d0
        ENDIF
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          DO nolist=1,NELIST(0)
            ne=NELIST(nolist)
            IF(IBT(1,1,NBJ(1,ne)).EQ.3) THEN !Simplex
              IF(.NOT.ERR) THEN
                XI(1)=0.3333d0
                XI(2)=0.3333d0
                XI(3)=0.0d0
              ELSE
                XI(1)=0.1666d0
                XI(2)=0.3333d0
                XI(3)=0.0d0
              ENDIF
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '          NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
              DO nj=1,NJT
                nb=NBJ(nj,ne)
                X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,1,XI,XE(1,nj))
              ENDDO
            ENDIF
            IF(nr.EQ.NRE(ne)) THEN
C ***         Draw element number
              IF(IBT(1,1,NBJ(1,ne)).NE.3) THEN !Not simplex
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)

C new MPN 10-Apr-96: Draw elem # at central Xi location not averaged
C                    global nodal coords
                CALL XEXW(0,IBT,IDO,INP,NAN,NBJ(1,ne),nr,XE,XW,XI,
     '            ERROR,*9999)
                DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj,nr)
                  X(nj)=XW(nj,1)
                ENDDO
              ENDIF
              CALL XZ(ITYP10(1),X,Z)
              DO noiw=1,NTIW
                IW=IWK(noiw)
                CALL ACWK(iw,1,ERROR,*9999)
                IF(.NOT.ERR) THEN
                  CALL SGELEM(INDEX,ISEG,ISELNO(iw,ne),iw,MXI(1,ne),
     '              NBJ(1,ne),ne,NLL(1,ne),NPL,nr,
     '              BOUNDARY,BOX,CSEG,DL,XP,Z,ERROR,*9999)
                ELSE
                  CALL SGERR(INDEX,ISEG,ISERR(iw,ne),iw,MXI(1,ne),
     '              ne,nr,CSEG,NEERR,Z,ERROR,*9999)
                ENDIF
C old
C            DO nj=1,NJT
C              ZCENTR(nj)=0.0D0
C            ENDDO
C            nb=NBJ(1,ne)
C            DO nn=1,NNT(nb)
C              np=NPNE(nn,nb,ne)
C              DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
C                X(nj)=XP(1,1,nj,np)
C              ENDDO
C              CALL XZ(ITYP10(1),X,Z)
C              DO nj=1,NJT
C                ZCENTR(nj)=ZCENTR(nj)+Z(nj)
C              ENDDO
C            ENDDO
C            DO nj=1,NJT
C              ZC(nj,ne)=ZCENTR(nj)/DBLE(NNT(nb))
C            ENDDO
C            CALL SGELEM(INDEX,ISEG,ISELNO(iw,ne),iw,MXI(1,ne),
C     '        NBJ(1,ne),ne,NLL(1,ne),NPL,nr,
C     '        BOUNDARY,BOX,CSEG,DL,XP,ZC(1,ne),ERROR,*9999)
                CALL DAWK(iw,1,ERROR,*9999)
              ENDDO !noiw (iw)
            ENDIF
          ENDDO !nolist (ne)
        ENDDO !no_nrlist (nr)
      ENDIF

      CALL EXITS('DRELEM')
      RETURN
 9999 CALL ERRORS('DRELEM',ERROR)
      CALL EXITS('DRELEM')
      RETURN 1
      END


