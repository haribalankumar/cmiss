      SUBROUTINE DISECT(IBT,IDO,INP,ISEG,ISSECT,NBH,NBJ,
     '  NEELEM,NHE,NHP,NKJE,NKH,NPF,NPLIST,NPNE,NPNODE,
     '  NRE,NVHP,NVJE,NYNE,NYNP,
     '  CSEG,SE,STRING,XA,XE,XP,YP,ZA,ZD,ZP,ERROR,*)

C#### Subroutine: DISECT
C###  Description:
C###    DISECT displays section time profiles, chord pressure plots,
C###    or transmural fibre plots.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'map000.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'post00.cmn'
      INCLUDE 'press00.cmn'
      INCLUDE 'sect00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISEG(*),ISSECT(NGRSEGM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKJE(NKM,NNM,NJM,NEM),
     '  NKH(NHM,NPM,NCM,0:NRM),NPF(9,NFM),NPLIST(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NRE(NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZD(NJM,NDM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,INDEX,INDEX_POLYLINE,INSTAT,IW,IW2,
     '  IWK(6),N3CO,ne,NELEM(10),ng,NOELEM,NONODE,nr,NTELEM,NTIW,nx
      REAL*8 XREF,XWC,YREF,YWC
      REAL*8 RFROMC,XD(3),XI(3)
      CHARACTER CHAR*2,FILE*255,TYPE*10
      LOGICAL ABBREV,CBBREV,MOUSE,RETAIN_POSTSCRIPT

      CALL ENTERS('DISECT',*9999)

      nx=1 ! temporary cpb 22/11/94

      XREF=0.d0 ! Initialisation
      YREF=0.d0 ! Ditto

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE02,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM display section<;FILENAME>
C###  Parameter:       <at NODES[all]>
C###    Specify which nodes to display
C###  Parameter:       <scale FACTOR[1]>
C###    Specify a scaling function.
C###  Parameter:       <rgb=RGB[blue]>
C###    Define colour (eg red,green,blue,cyan)
C###  Description:
C###    Display the solution variable vs position (at specified nodes),
C###    picking up the nodal parameters from FILENAME.history.

        OP_STRING(1)=STRING(1:IEND)
     '    //'<;FILENAME['//FILE02(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<at NODES[all]>'
        OP_STRING(3)=BLANK(1:15)//'<scale FACTOR[1]>'
        OP_STRING(4)=BLANK(1:15)//'<rgb=RGB[blue]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM display section pressure
C###  Parameter:       <at CHORDS[all]>
C###    Specify location to display pressure on a CHORD location.
C###  Parameter:       <scale FACTOR[1]>
C###    Specify a scaling factor
C###  Parameter:       <rgb=RGB[blue]>
C###    Define colour (eg red,green,blue,cyan)
C###  Description:
C###    Display a graph of the pressure distribution over a sail.
C###    Refer to Masters Thesis Fiona McPheat.

        OP_STRING(1)=STRING(1:IEND)//' pressure'
        OP_STRING(2)=BLANK(1:15)//'<at CHORDS[all]>'
        OP_STRING(3)=BLANK(1:15)//'<scale FACTOR[1]>'
        OP_STRING(4)=BLANK(1:15)//'<rgb=RGB[blue]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM display section;m<;FILENAME> fibre
C###  Parameter:       <mu INCREMENT[5deg]>
C###    Specify mu position to display
C###  Parameter:       <theta INCREMENT[5deg]
C###    Specify theta position to display
C###  Parameter:       <on WS[4]>
C###    Specify a workstation number
C###  Parameter:       <rgb=RGB[blue]>
C###    Define colour (eg red,green,blue,cyan)
C###  Description:
C###    Display the fibre variation through the wall of the geometry at
C###    locations selected by the mouse from the specified workstation.

        OP_STRING(1)=STRING(1:IEND) //';m<;FILENAME> fibre'
        OP_STRING(2)=BLANK(1:15)//'<mu INCREMENT[5deg]>'
        OP_STRING(3)=BLANK(1:15)//'<theta INCREMENT[5deg]>'
        OP_STRING(4)=BLANK(1:15)//'<on WS[4]>'
        OP_STRING(5)=BLANK(1:15)//'<rgb=RGB[blue]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','DISECT',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        nr=1 !may need generalizing
        MOUSE=.FALSE.
        IF(ABBREV(COQU(noco,1),'M',1)) THEN
          MOUSE=.TRUE.
        ENDIF

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLUE')
        ENDIF

        IF(MOUSE) THEN
          CALL ASSERT(NJT.EQ.3,'>>Only for 3D',ERROR,*9999)
          IF(CBBREV(CO,'FIBRE',3,noco+1,NTCO,N3CO)) THEN
            TYPE='FIBRE'
          ELSE IF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO)) THEN
            TYPE='FIELD'
          ENDIF
          IF(CBBREV(CO,'MU',1,noco+1,NTCO,N3CO)) THEN
            DMU=RFROMC(CO(N3CO+1))*PI/180.0d0
          ELSE
            DMU=5.0d0*PI/180.0d0
          ENDIF
          IF(CBBREV(CO,'THETA',1,noco+1,NTCO,N3CO)) THEN
            DTHETA=RFROMC(CO(N3CO+1))*PI/180.0d0
          ELSE
            DTHETA=5.0d0*PI/180.0d0
          ENDIF
          IF(NTCOQU(noco).EQ.2) THEN
            TYPE(6:9)='FILE'
            FILE=COQU(noco,2)
          ENDIF
          CALL WS_LIST(IWK,4,NTIW,noco,NTCO,CO,ERROR,*9999)
          IF(POSTSCRIPT) THEN
            RETAIN_POSTSCRIPT=.TRUE.
          ELSE
            RETAIN_POSTSCRIPT=.FALSE.
          ENDIF

        ELSE IF(.NOT.MOUSE) THEN
          IF(CBBREV(CO,'PRESSURE',1,noco+1,NTCO,N3CO)) THEN
            TYPE='PRESSURE'
          ELSE
            TYPE='FILE'
          ENDIF

          IF(TYPE(1:8).EQ.'PRESSURE') THEN
            IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),NPM,NPLIST(0),NPLIST(1),ERROR,
     '          *9999)
            ELSE
              NPLIST(0)=NPNODE(0,NR)
              DO NONODE=1,NPLIST(0)
                NPLIST(NONODE)=NPNODE(NONODE,NR)
              ENDDO
            ENDIF

          ELSE IF(TYPE(1:4).EQ.'FILE') THEN
            IF(NTCOQU(noco).EQ.1) THEN
              FILE=COQU(noco,1)
            ELSE
              FILE=FILE02
            ENDIF

            IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),NPM,NPLIST(0),NPLIST(1),ERROR,
     '          *9999)
            ELSE
              NPLIST(0)=NPNODE(0,NR)
              DO NONODE=1,NPLIST(0)
                NPLIST(NONODE)=NPNODE(NONODE,NR)
              ENDDO
            ENDIF
          ENDIF

          IF(CBBREV(CO,'SCALE',1,noco+1,NTCO,N3CO)) THEN
            FACTOR=RFROMC(CO(N3CO+1))
          ELSE
            FACTOR=1.0D0
          ENDIF

        ENDIF

        IF(TYPE(1:8).EQ.'PRESSURE') THEN
          PRESS_MAX=0.0D0
          DO NOELEM=1,NEELEM(0,NR)
            NE=NEELEM(NOELEM,NR)
            DO ng=1,9
              IF(PRESS(ng,NE).GT.PRESS_MAX) PRESS_MAX=PRESS(ng,NE)
            ENDDO
          ENDDO
          NTSECT=NTSECT+1
          CALL ASSERT(NTSECT.LE.NRM,'>>NRM too small',ERROR,*9999)
          IW2=11
          CALL ACWK(IW2,0,ERROR,*9999)
          CALL SGSECT(INDEX,IBT,IDO,INP,ISEG,ISSECT(NTSECT),IW2,NBH,NBJ,
     '      NEELEM,NELEM,NHE(1,nx),NHP(1,nr,nx),NKJE,NKH(1,1,1,nr),
     '      NPF,NPLIST,NPNE,NPNODE,NRE,NTELEM,NTSECT,
     '      NVHP(1,1,1,nr),NVJE,NYNE,NYNP,CSEG,SE,TYPE,XA,XD,XE,XI,XP,
     '      YP(1,1,nx),ZA,ZD,ZP,ERROR,*9999)
          CALL DAWK(IW2,0,ERROR,*9999)

        ELSE IF(TYPE(1:4).EQ.'FILE') THEN
          NTSECT=NTSECT+1
          CALL ASSERT(NTSECT.LE.NRM,'>>NRM too small',ERROR,*9999)
          CALL STRING_TRIM(FILE,IBEG,IEND)
C!!! same unit as IOOUT for `set output'!
          CALL OPENF(9,'DISK',FILE(IBEG:IEND)//'.history','OLD',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
          IW2=11
          CALL ACWK(IW2,0,ERROR,*9999)
          CALL SGSECT(INDEX,IBT,IDO,INP,ISEG,ISSECT(NTSECT),IW2,NBH,NBJ,
     '      NEELEM,NELEM,NHE(1,nx),NHP(1,nr,nx),NKJE,NKH(1,1,1,nr),
     '      NPF,NPLIST,NPNE,NPNODE,NRE,NTELEM,NTSECT,
     '      NVHP(1,1,1,nr),NVJE,NYNE,NYNP,CSEG,SE,TYPE,XA,XD,XE,XI,XP,
     '      YP(1,1,nx),ZA,ZD,ZP,
     '      ERROR,*9999)
          CALL DAWK(IW2,0,ERROR,*9999)
          CALL CLOSEF(9,ERROR,*9999)

        ELSE IF(TYPE(1:5).EQ.'FIBRE') THEN
          IW=IWK(1)
          IF(TYPE(6:9).EQ.'FILE') THEN
            CALL STRING_TRIM(FILE,IBEG,IEND)
C!!! same unit as IOOUT for `set output'!
            CALL OPENF(9,'DISK',FILE(IBEG:IEND)//'.opfib','NEW',
     '        'SEQUEN','FORMATTED',132,ERROR,*9999)
          ENDIF
          WRITE(OP_STRING,'('' >>Locate sites on map'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          IF(RETAIN_POSTSCRIPT) POSTSCRIPT=.TRUE.
          CALL ACWK(IW,0,ERROR,*9999)
          CALL LOCATOR(INSTAT,XREF,XWC,YREF,YWC,ERROR,
     '      *9999)
          DO WHILE (INSTAT.EQ.1) !valid input received
            NTSECT=NTSECT+1
            WRITE(CHAR,'(I2)') NTSECT
            CALL STRING_TRIM(CHAR,IBEG,IEND)
            XD(1)=XWC
            XD(2)=YWC
            PROJEC='RECTANGULAR'
            CALL TEXT(1,IW,CHAR(IBEG:IEND),XD,ERROR,*9999)
            PROJEC=MAP_PROJEC
            CALL PSCOORD1(XD,XWC,YWC,ERROR,*9999)
            CALL PSCOORD2(NBJ,NEELEM,NELEM,NKJE,NPF,NPNE,
     '        NTELEM,NVJE,SE,XA,XD,XE,XI,XP,ERROR,*9999)
            CALL DAWK(IW,0,ERROR,*9999)
            IF(TYPE(6:9).EQ.'FILE') THEN
              WRITE(9,'(/'' Fibre angles at prolate coordinates'','
     '          //'3F8.2)') XD(1),XD(2)*180.0D0/PI,XD(3)*180.0D0/PI
            ENDIF
            IW2=12
            IF(RETAIN_POSTSCRIPT) POSTSCRIPT=.TRUE.
            CALL ACWK(IW2,0,ERROR,*9999)
            CALL SGSECT(INDEX,IBT,IDO,INP,ISEG,ISSECT(NTSECT),IW2,NBH,
     '        NBJ,NEELEM,NELEM,NHE(1,nx),NHP(1,nr,nx),NKJE,
     '        NKH(1,1,1,nr),NPF,NPLIST,NPNE,NPNODE,NRE,NTELEM,
     '        NTSECT,NVHP(1,1,1,nr),NVJE,
     '        NYNE,NYNP,CSEG,SE,TYPE,XA,XD,XE,XI,XP,YP(1,1,nx),
     '        ZA,ZD,ZP,ERROR,*9999)
            CALL DAWK(IW2,0,ERROR,*9999)
            IF(RETAIN_POSTSCRIPT) POSTSCRIPT=.TRUE.
            CALL ACWK(IW,0,ERROR,*9999)
            CALL LOCATOR(INSTAT,XREF,XWC,YREF,YWC,
     '        ERROR,*9999)
          ENDDO
          CALL DAWK(IW,0,ERROR,*9999)
          CALL CLOSEF(9,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('DISECT')
      RETURN
 9999 CALL ERRORS('DISECT',ERROR)
      CLOSE(UNIT=9)
      CALL EXITS('DISECT')
      RETURN 1
      END


