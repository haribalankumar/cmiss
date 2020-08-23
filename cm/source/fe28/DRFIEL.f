      SUBROUTINE DRFIEL(IBT,IDO,INP,ISEG,ISFIEL,ISSCAL,MXI,NBH,NBJ,
     '  NEELEM,NELIST,NGAP,NHE,NHP,NHQ,NKH,NKHE,NKJ,NKJE,NLATNE,
     '  NPF,
     '  NPNE,NPNODE,NPNY,NQNLAT,
     '  NQNE,NQNY,NQS,NQXI,NRE,NRLIST,NRLIST2,
     '  NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,NYNR,NYQNR,CURVCORRECT,PG,SE,
     '  XA,XE,XG,XP,XQ,YG,YQ,YQS,YP,ZA,ZE,ZG,ZP,CSEG,STRING,ERROR,*)

C#### Subroutine: DRFIEL
C###  Description:
C###    DRFIEL draws a field variable or segment.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'gks001.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'scal00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISEG(*),ISFIEL(NWM,NEM),ISSCAL(NWM,NGRSEGM),MXI(2,NEM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NGAP(NIM,NBM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NHQ(NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),NLATNE(NEQM+1),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NQNLAT(NEQM*NQEM),
     '  NQNE(NEQM,NQEM),NQNY(2,NYQM,0:NRCM,NXM),NQS(NEQM),
     '  NQXI(0:NIM,NQSCM),NRE(NEM),NRLIST(0:NRM),NRLIST2(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),
     '  SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  XQ(NJM,NQM),
     '  YG(NIYGM,NGM,NEM),YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM),
     '  YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IBEG3,IDTYPE,ID_TYPE,IEND,IEND1,IEND2,
     '  IEND3,IFROMC,INDEX,iw,IWK(6),na,N3CO,nb,nb_type,NBJ_TEMP(12),
     '  nc,ne,nh,nhx,NIQLIST(0:1),NIQSLIST(0:40),NIYLIST(0:16),nj,njj1,
     '  njj2,nk,noiw,nolist,nonr,no_time,notime,nr,NTIW,nonrlist,
     '  nonode,np,NUMTIMEDATA,nx,nxc,ny
      REAL DELAY
      REAL*8 DISPLACEMENT,RFROMC,SCALE_FACTOR,TIME,TMAX,TMIN,YPMAX(16),
     '  YPMIN(16),YMAX,YMIN
      CHARACTER CHAR2*9,CHAR3*9,FILEFORMAT*6,HISTORYFNAME*100,TYPE*10,
     '  ERROR1*255
      LOGICAL ALL_REGIONS,ALL_WINDOWS,ANIMATE,BINARYFILE,CBBREV,
     '  COLOUR,DEFORMED,DRAW_SCALE,ENDFILE,YPDATA,YQDATA,YQSDATA

      CALL ENTERS('DRFIEL',*9999)
      nc=1 !temporary MPN 12-Nov-94

      ANIMATE=.FALSE.
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF
        NRLIST(0)=1
        NRLIST(1)=nr
        nx=1 !temporary
        IF(NJ_LOC(njl_fibr,0,nr).GT.0) THEN !field variables defined
          CHAR2='field'
        ELSE IF(ITYP2(nr,nx).GT.0) THEN !dependent variables defined
          CHAR2='dependent'
        ELSE
          CHAR2='Gauss'          !Gauss variables defined
        ENDIF
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
        IF(NINDICES.GE.150) THEN
          CHAR3='colour'
        ELSE
          CHAR3='greyscale'
        ENDIF
        CALL STRING_TRIM(CHAR3,IBEG3,IEND3)

C---------------------------------------------------------------------

C#### Command: FEM draw field
C###  Parameter:    <(field/dependent/Gauss/grid)[Gauss] ID#>
C###  Parameter:    <element (#s/all)[all]>
C###    Define what elements are to be included.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation (GX window) to draw the points on.
C###  Parameter:    <basis (geometry/#)[geometry]>
C###  Parameter:    <(greyscale/colour)[greyscale]>
C###    Draws field in greyscale or colour.
C###  Parameter:    <zmin=(calculated/ZMIN#)[calculated]>
C###    Either specify zmin or calculate it
C###  Parameter:    <zmax=(calculated/ZMAX#)[calculated]>
C###    Either specify zmax or calculate it
C###  Parameter:    <(scale/noscale)[scale]>
C###    Choice of showing the reference scale or not.
C###  Parameter:    <region #[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:    <nc #[1]>

C###  Description:
C###    Draws the specified field, over the specified elements on the
C###    specified workstation, using the given interpolation basis.
C###    Draws in either colour or greyscale, between the specified
C###    maximum (zmax) and minimum (zmin) values.

        OP_STRING(1)=STRING(1:IEND)//' <(field/dependent/Gauss/grid)['
     '    //CHAR2(IBEG2:IEND2)//'] ID#>'
        OP_STRING(2)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<basis (geometric/#)[geometric]>'
        OP_STRING(5)=BLANK(1:15)
     '    //'<(greyscale/colour)['//CHAR3(IBEG3:IEND3)//']>'
        OP_STRING(6)=BLANK(1:15)//
     '    '<zmin=(calculated/ZMIN#)[calculated]>'
        OP_STRING(7)=BLANK(1:15)//
     '    '<zmax=(calculated/ZMAX#)[calculated]>'
        OP_STRING(8)=BLANK(1:15)//'<(scale/noscale)[scale]>'
        OP_STRING(9)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(10)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(11)=BLANK(1:15)//'<nc #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw field deformed
C###  Parameter:    <(field/dependent/gauss/grid)[Gauss] ID#>
C###    Choice of drawing field, dependant variable, gauss point values
C###    or grid point variables.
C###  Parameter:    <element (#s/all)[all]>
C###    Define what elements are to be included.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation (GX window) to draw the points on.
CC###  Parameter:    <basis (geometric/#)[geometric]>
C###  Parameter:    <(greyscale/colour)[greyscale]>
C###    Draws field in greyscale or colour.
C###  Parameter:    <zmin=(calculated/ZMIN#)[calculated]>
C###    Either specify zmin or calculate it
C###  Parameter:    <zmax=(calculated/ZMAX#)[calculated]>
C###    Either specify zmax or calculate it
C###  Parameter:    <(scale/noscale)[scale]>
C###    Choice of showing the reference scale or not.
C###  Parameter:    <scale_factor FACTOR#[1.0]>
C###    Specify the scale factor
C###  Parameter:    <region #[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:    <nc #[1]>
C###  Description:
C###    Draws the specified field, over the specified elements on the
C###    specified workstation, using the given interpolation basis.
C###    Draws in either colour or greyscale, between the specified
C###    maximum (zmax) and minimum (zmin) values.

        OP_STRING(1)=STRING(1:IEND)//' deformed'
        OP_STRING(2)=BLANK(1:15)//' <(field/dependent/Gauss/grid)['
     '    //CHAR2(IBEG2:IEND2)//'] ID#>'
        OP_STRING(3)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(5)=BLANK(1:15)//'<basis (geometric/#)[geometric]>'
        OP_STRING(6)=BLANK(1:15)
     '    //'<(greyscale/colour)>['//CHAR3(IBEG3:IEND3)//']'
        OP_STRING(7)=BLANK(1:15)//
     '    '<zmin=(calculated/ZMIN#)[calculated]>'
        OP_STRING(8)=BLANK(1:15)//
     '    '<zmax=(calculated/ZMAX#)[calculated]>'
        OP_STRING(9)=BLANK(1:15)//'<(scale/noscale)[scale]>'
        OP_STRING(10)=BLANK(1:15)//'<scale_factor FACTOR#[1.0]>'
        OP_STRING(11)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(12)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(13)=BLANK(1:15)//'<nc #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw field animate
C###  Parameter:    <history FILENAME[$current]>
C###    Specify wiether the field history is read from a
C###    file name or is current
C###  Parameter:    <(ascii/binary)[ascii]>
C###    Specify the format of the history file.
C###  Parameter:    <delay #[0.0]>
C###  Parameter:    <element (#s/all)[all]>
C###    Define what elements are to be included.
C###  Parameter:    <window #[1]>
C###  Parameter:    <(scale/noscale)[scale]>
C###    Choice of showing the reference scale or not.
C###  Parameter:    <region #[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <intoregion (#/same)[same]>'
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Draws an animated field from a history file. Delay controls
C###    the length of the pause between successive time frames.

        OP_STRING(1)=STRING(1:IEND)//' animate'
        OP_STRING(2)=BLANK(1:15)//'<history FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(3)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(4)=BLANK(1:15)//'<delay #[0.0]>'
        OP_STRING(5)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<window #[1]>'
        OP_STRING(7)=BLANK(1:15)//'<(scale/noscale)[scale]>'
        OP_STRING(8)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(9)=BLANK(1:15)//'<intoregion (#/same)[same]>'
        OP_STRING(9)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRFIEL',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        nr=NRLIST(1)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'ANIMATE',2,noco+1,NTCO,N3CO)) THEN
          ANIMATE=.TRUE.
        ELSE
          ANIMATE=.FALSE.
        ENDIF

C LKC & MLB - text colour initialization
        INDEX=1

        IF(.NOT.ANIMATE) THEN
          IF(CBBREV(CO,'FIELD',2,noco+1,NTCO,N3CO)) THEN
            TYPE='FIELD'
            IDTYPE=IFROMC(CO(N3CO+1))
            ID_TYPE=IDTYPE
          ELSE IF(CBBREV(CO,'DEPENDENT',2,noco+1,NTCO,N3CO)) THEN
            TYPE='DEPENDENT'
            IDTYPE=IFROMC(CO(N3CO+1))
            ID_TYPE=IDTYPE
          ELSE IF(CBBREV(CO,'GAUSS',2,noco+1,NTCO,N3CO)) THEN
            TYPE='GAUSS'
            IDTYPE=IFROMC(CO(N3CO+1))
            ID_TYPE=IDTYPE
          ELSE IF(CBBREV(CO,'GRID',2,noco+1,NTCO,N3CO)) THEN
            TYPE='GRID'
            IDTYPE=IFROMC(CO(N3CO+1))
            ID_TYPE=IDTYPE
          ELSE
            IF(CBBREV(CO,'CLASS',2,noco+1,NTCO,N3CO)) THEN
              CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
              nxc=NXLIST(1)
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
              CALL ASSERT(nx.GT.0,'>>No nx defined for this class',
     '          ERROR,*9999)
              IDTYPE=NH_LOC(1,nx)
              ID_TYPE=1
            ELSE
              nx=1
              IDTYPE=1
              ID_TYPE=1
            ENDIF

            IF(NJ_LOC(NJL_FIEL,0,nr).GT.0) THEN
              !field variables defined
              TYPE='FIELD'
            ELSE IF(ITYP2(nr,nx).GT.0) THEN !dependent variables defined
              TYPE='DEPENDENT'
            ELSE !Gauss variables defined
              TYPE='GAUSS'
            ENDIF
          ENDIF

          IF((TYPE(1:9).EQ.'DEPENDENT').OR.(TYPE(1:4).EQ.'GRID')) THEN
            CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
            nxc=NXLIST(1)
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,'>>No nx defined for this class',
     '        ERROR,*9999)
          ELSE
            nx=1
          ENDIF

          IF(TYPE(1:4).NE.'GRID') THEN
            IF(ITYP5(nr,nx).EQ.2) THEN !time dep
              IF(ITYP2(nr,nx).EQ.9) THEN !ionic curent activation
                CALL ASSERT(.FALSE.,
     '            '>>Use draw field grid for grid problems',ERROR,*9999)
              ENDIF
            ENDIF
          ENDIF

          IF(CBBREV(CO,'DEFORMED',2,noco+1,NTCO,N3CO)) THEN
            DEFORMED=.TRUE.
          ELSE
            DEFORMED=.FALSE.
          ENDIF

          IF(CBBREV(CO,'BASIS',2,noco+1,NTCO,N3CO)) THEN
            nb_type=IFROMC(CO(N3CO+1))
          ELSE
            nb_type=0
          ENDIF

          IF(CBBREV(CO,'ny',1,noco+1,NTCO,N3CO)) THEN
C GMH 2/9/95 Unused          DEPEND=.TRUE.
            ny=IFROMC(CO(N3CO+1))
            CALL ASSERT(NHM.EQ.NJM,'>>nhm.NE.NJM',ERROR,*9999)
            CALL ASSERT(ny.LE.16,'>>ny.gt.16',ERROR,*9999)
          ELSE
C GMH 2/9/95 Unused          DEPEND=.FALSE.
          ENDIF

          IF(CBBREV(CO,'GREYSCALE',1,noco+1,NTCO,N3CO)) THEN
            COLOUR=.FALSE.
          ELSE
            IF(NINDICES.GE.150) THEN
              COLOUR=.TRUE.
            ELSE
              COLOUR=.FALSE.
            ENDIF
          ENDIF

          IF(CBBREV(CO,'NOSCALE',3,noco+1,NTCO,N3CO)) THEN
            DRAW_SCALE=.FALSE.
          ELSE
            DRAW_SCALE=.TRUE.
          ENDIF

          IF(CBBREV(CO,'SCALE_FACTOR',7,noco+1,NTCO,N3CO)) THEN
            SCALE_FACTOR=RFROMC(CO(N3CO+1))
          ELSE
            SCALE_FACTOR=1.0d0
          ENDIF

          IF(TYPE(1:9).EQ.'DEPENDENT'.OR.
     '      TYPE(1:5).EQ.'FIELD') THEN !find range
            CALL CM_RANGE(IDTYPE,nb_type,NBH,NBJ,NEELEM,NHE(1,nx),
     '        NHP(1,0,nx),NKH,NKHE,NKJE,NPF,NPNE,NPNODE,NRLIST,NVHE,
     '        NVHP,NVJE,NW(1,1,nx),nx,NYNE,NYNP,
     '        CURVCORRECT,PG,SE,TYPE,XA,XE,XG,XP,
     '        YG,YP(1,1,nx),ZA,ZE,ZG,ZP,ERROR,*9999)
          ENDIF

C rgb 14/05/98 - adding nc option for draw field dependent
          IF(TYPE(1:9).EQ.'DEPENDENT') THEN
            IF(CBBREV(CO,'NC',2,noco+1,NTCO,N3CO)) THEN
              nc=IFROMC(CO(N3CO+1))
              CALL ASSERT(nc.LE.NCM,'>>nc > NCM',ERROR,*9999)
            ENDIF
          ELSE
            IF(CBBREV(CO,'NC',2,noco+1,NTCO,N3CO)) THEN
              ERROR='nc only set up for draw field dependent'
              GOTO 9999
            ENDIF
          ENDIF

          IF(CBBREV(CO,'ZMIN',4,noco+1,NTCO,N3CO)) THEN
            ZMINI=RFROMC(CO(N3CO+1))
          ENDIF
          IF(CBBREV(CO,'ZMAX',4,noco+1,NTCO,N3CO)) THEN
            ZMAXI=RFROMC(CO(N3CO+1))
          ENDIF
          ZDIFF=ZMAXI-ZMINI
          IF(DABS(ZDIFF).LT.1.0d-6) ZDIFF=1.0d0

          IF(DEFORMED) THEN
            CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '        NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,
     '        NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
C CS 25/5/98 new
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,nhx,nr)
                nh=NH_LOC(nhx,nx)
                DO nk=1,NKH(nh,np,nc,nr)
                  IF(ITYP6(nr,nx).EQ.1.OR.KTYP58(nr).EQ.2) THEN
                    DISPLACEMENT=ZP(nk,1,nh,np,nc) !linear analysis
                  ELSE IF(ITYP6(nr,nx).EQ.2) THEN !non-linear analy
                    DISPLACEMENT=ZP(nk,1,nh,np,nc)-XP(nk,1,nj,np)
                  ENDIF
                  IF(nk.LE.NKJ(nj,np)) THEN
                    ZP(nk,1,nh,np,nc)=XP(nk,1,nj,np)+
     '                SCALE_FACTOR*DISPLACEMENT
                  ELSE
                    ZP(nk,1,nh,np,nc)=SCALE_FACTOR*DISPLACEMENT
                  ENDIF
                ENDDO !nk
              ENDDO !nhx
            ENDDO !nonode
          ENDIF

        ELSE !Animate
C!!! Setting nx for animation cs 24-2-1997
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this class',
     '      ERROR,*9999)

          IF(CBBREV(CO,'HISTORY',2,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
            HISTORYFNAME=CO(N3CO+1)(IBEG:IEND)
          ELSE
            CALL STRING_TRIM(FILE00,IBEG,IEND)
            HISTORYFNAME=FILE00(IBEG:IEND)
          ENDIF

          IF(CBBREV(CO,'BINARY',1,noco+1,NTCO,N3CO)) THEN
            BINARYFILE=.TRUE.
            FILEFORMAT='BINARY'
          ELSE
            BINARYFILE=.FALSE.
            FILEFORMAT='ASCII'
          ENDIF

          IF(CBBREV(CO,'DELAY',2,noco+1,NTCO,N3CO)) THEN
            DELAY=REAL(RFROMC(CO(N3CO+1)))
          ELSE
            DELAY=0.0
          ENDIF

          CALL PARSE_WINDOWS(IWK,6,noco,NTCO,NTIW,CO,ALL_WINDOWS,
     '      ERROR,*9999)
          iw=IWK(1)

          IF(CBBREV(CO,'NOSCALE',3,noco+1,NTCO,N3CO)) THEN
            DRAW_SCALE=.FALSE.
          ELSE
            DRAW_SCALE=.TRUE.
          ENDIF

          IF(CBBREV(CO,'INTOREGION',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NRM,NRLIST2(0),NRLIST2(1),
     '        ERROR,*9999)
          ELSE
            NRLIST2(0)=NRLIST(0)
            DO nonr=1,NRLIST(0)
              NRLIST2(nonr)=NRLIST(nonr)
            ENDDO
          ENDIF

        ENDIF


C        ERROR_FLAG=.FALSE.
        IF(.NOT.ANIMATE) THEN
          DO noiw=1,NTIW
            IW=IWK(noiw)
            CALL ACWK(iw,1,ERROR,*9999)
C An attempt at parralisation by cst 24/8/96 works on a single thread
C C$DOACROSS local(nolist,ne,njj1,njj2,NBJ_TEMP,XE,ZE,nb,
C C$&             CSEG,nj)
C C$&        share(ERROR_FLAG,NELIST,NJ_LOC)
            DO nolist=1,NELIST(0)
C               IF(.NOT.ERROR_FLAG) THEN
                ne=NELIST(nolist)
                DO njj1=1,3 !geom + fibres + field
                  DO njj2=1,NJ_LOC(njj1,0,nr)
                    nj=NJ_LOC(njj1,njj2,nr)
                    CALL ASSERT(nj.LE.12,'>>ERROR: increase size of '
     '                //'NBJ_temp to NJM',ERROR,*9999)
                    IF(nb_type.EQ.0) THEN
                      NBJ_TEMP(nj)=NBJ(nj,ne)
                    ELSE
                      NBJ_TEMP(nj)=nb_type
                    ENDIF
                  ENDDO !njj2
                ENDDO !njj1
                IF(TYPE(1:4).NE.'GRID') THEN
                  CALL XPXE(NBJ_TEMP,NKJE(1,1,1,ne),NPF(1,1),
     '              NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '              SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
                  IF(DEFORMED) THEN
                    CALL ZPZE(NBJ_TEMP,nc,NJ_LOC(NJL_GEOM,0,nr),
     '                NKHE(1,1,1,ne),NPF(1,1),
     '                NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1,nx),nx,
     '                CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,
     '                ZP,ERROR,*9999)
                  ENDIF
                ENDIF
                nb=NBJ_TEMP(1)
                IF(TYPE(1:5).EQ.'FIELD') THEN
                  IF(DEFORMED) THEN
                    CALL SGFIEL(IBT,IDO,INP,ISEG,ISFIEL(iw,ne),iw,
     '                MXI(1,ne),NBJ(NJ_LOC(NJL_FIEL,ID_TYPE,nr),ne),
     '                NBJ_TEMP,ne,NGAP(1,nb),NLATNE,NQNE,NQNLAT,
     '                NQS,
     '                NQXI,CSEG,TYPE,ZE,XQ,YG(1,1,ne),YQ(1,1,1,nx),
     '                XE(1,NJ_LOC(NJL_FIEL,ID_TYPE,nr)),ERROR,*9999)
                  ELSE
                    CALL SGFIEL(IBT,IDO,INP,ISEG,ISFIEL(iw,ne),iw,
     '                MXI(1,ne),NBJ(NJ_LOC(NJL_FIEL,ID_TYPE,nr),ne),
     '                NBJ_TEMP,ne,NGAP(1,nb),NLATNE,NQNE,NQNLAT,
     '                NQS,
     '                NQXI,CSEG,TYPE,XE,XQ,YG(1,1,ne),YQ(1,1,1,nx),
     '                XE(1,NJ_LOC(NJL_FIEL,ID_TYPE,nr)),ERROR,*9999)
                  ENDIF
                ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
                  CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx),NKHE(1,1,1,ne),
     '              NPF(1,1),NPNE(1,1,ne),NRE(ne),NVHE(1,1,1,ne),
     '              NW(ne,1,nx),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '              ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
                  IF(DEFORMED) THEN
                    CALL SGFIEL(IBT,IDO,INP,ISEG,ISFIEL(iw,ne),iw,
     '                MXI(1,ne),NBH(IDTYPE,1,ne),NBJ_TEMP,ne,
     '                NGAP(1,nb),NLATNE,NQNE,NQNLAT,
     '                NQS,NQXI,CSEG,TYPE,
     '                ZE,XQ,YG(1,1,ne),
     '                YQ(1,1,1,nx),ZE(1,ID_TYPE),ERROR,*9999)
                  ELSE
                    CALL SGFIEL(IBT,IDO,INP,ISEG,ISFIEL(iw,ne),iw,
     '                MXI(1,ne),NBH(IDTYPE,1,ne),NBJ_TEMP,ne,
     '                NGAP(1,nb),NLATNE,NQNE,NQNLAT,
     '                NQS,NQXI,CSEG,TYPE,
     '                XE,XQ,YG(1,1,ne),
     '                YQ(1,1,1,nx),ZE(1,ID_TYPE),ERROR,*9999)
                    ENDIF
                ELSE IF(TYPE(1:5).EQ.'GAUSS') THEN
                  IF(DEFORMED) THEN
                    CALL SGFIEL(IBT,IDO,INP,ISEG,ISFIEL(iw,ne),iw,
     '                MXI(1,ne),nb,NBJ_TEMP,ne,NGAP(1,nb),NLATNE,NQNE,
     '                NQNLAT,
     '                NQS,NQXI,CSEG,TYPE,ZE,XQ,YG(ID_TYPE,1,ne),
     '                YQ(1,1,1,nx),XE(1,1),ERROR,*9999)
                  ELSE
                    CALL SGFIEL(IBT,IDO,INP,ISEG,ISFIEL(iw,ne),iw,
     '                MXI(1,ne),nb,NBJ_TEMP,ne,NGAP(1,nb),NLATNE,
     '                NQNE,
     '                NQNLAT,NQS,NQXI,CSEG,TYPE,XE,XQ,YG(ID_TYPE,1,ne),
     '                YQ(1,1,1,nx),XE(1,1),ERROR,*9999)
                  ENDIF
                ELSE IF(TYPE(1:4).EQ.'GRID') THEN
                  CALL SGFIEL(IBT,IDO,INP,ISEG,ISFIEL(iw,ne),iw,
     '              MXI(1,ne),nb,NBJ_TEMP,ne,NGAP(1,nb),NLATNE,NQNE,
     '              NQNLAT,NQS,NQXI,CSEG,
     '              TYPE,XE,XQ,YG(1,1,ne),YQ(1,ID_TYPE,1,nx),XE(1,1),
     '              ERROR,*9999)
                ENDIF
C              ENDIF
C              IF(ERROR_FLAG) THEN
C C               This statement is designed to be skipped if no error
C C               occur. However if a error occurs within a subroutine
C C              the alternate return points to line 100 to set the flag
C 100            CONTINUE
C C$              call mp_setlock()
C                ERROR_FLAG=.TRUE.
C                WRITE(OP_STRING,'(/'' >>ERROR: An error occurred - '
C     '            //'results may be unreliable!'')')
C                CALL WRITES(IODI,OP_STRING,ERROR,*101)
C C$              call mp_unsetlock()
C 101            CONTINUE
C              ENDIF
            ENDDO

            IF(DRAW_SCALE) THEN
              IF(NTSCAL.EQ.0) NTSCAL=1
              CALL SGSCAL(INDEX,ISEG,ISSCAL(iw,NTSCAL),iw,NTSCAL,COLOUR,
     '          CSEG,ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDDO
        ELSE !animate
          ID_TYPE=1 !temporary
          nb_type=0 !temporary
          COLOUR=.TRUE. !temporary
          TYPE='DEPENDENT'
          NRLIST(0)=1
          NRLIST2(0)=1
          nr=NRLIST(1)
          NIYLIST(0)=1
          NIYLIST(1)=1
          NIQLIST(0)=0
          NIQLIST(1)=0
          NIQSLIST(0)=0
          NIQSLIST(1)=0
          na=1
          CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA,
     '      nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),
     '      YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',FILEFORMAT,HISTORYFNAME,
     '      'OPEN',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)
          YMAX=-RMAX
          YMIN=RMAX
          TMIN=0.0d0
          IF(BINARYFILE) THEN
            DO no_time=1,NUMTIMEDATA
              CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '          NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '          TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',
     '          FILEFORMAT,HISTORYFNAME,'TIME_DATA',ENDFILE,.TRUE.,
     '          YPDATA,YQDATA,YQSDATA,ERROR,*9999)
              CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '          NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),
     '          nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                IF(ZP(1,1,id_type,np,1).GT.YMAX)
     '            YMAX=ZP(1,1,id_type,np,1)
                IF(ZP(1,1,id_type,np,1).LE.YMIN)
     '            YMIN=ZP(1,1,id_type,np,1)
              ENDDO !no_nynr (ny)
            ENDDO !no_time
          ELSE
            ENDFILE=.FALSE.
            NUMTIMEDATA=0
            DO WHILE(.NOT.ENDFILE)
              CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '          NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '          TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',
     '          FILEFORMAT,HISTORYFNAME,'TIME_DATA',ENDFILE,.TRUE.,
     '          YPDATA,YQDATA,YQSDATA,ERROR,*9999)
              IF(.NOT.ENDFILE) THEN
                NUMTIMEDATA=NUMTIMEDATA+1
                CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '            NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),
     '            nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  IF(ZP(1,1,id_type,np,1).GT.YMAX)
     '              YMAX=ZP(1,1,id_type,np,1)
                  IF(ZP(1,1,id_type,np,1).LE.YMIN)
     '              YMIN=ZP(1,1,id_type,np,1)
                ENDDO !no_nynr (ny)
              ENDIF
            ENDDO
          ENDIF
          TMAX=TIME
          ZMINI=YMIN
          ZMAXI=YMAX
          ZDIFF=ZMAXI-ZMINI
          IF(DABS(ZDIFF).LT.1.0D-6) ZDIFF=1.0D0
          IF(DOP) THEN
            WRITE(OP_STRING,'('' YMIN = '',D11.4,'', YMAX = '',D11.4,'
     '        //'/,'' TMIN = '',D11.4,'', TMAX = '',D11.4,'
     '        //'/,'' ZMINI = '',D11.4,'', ZMAXI = '',D11.4,'
     '        //''', ZDIFF = '',D11.4)') YMIN,YMAX,TMIN,TMAX,ZMINI,
     '        ZMAXI,ZDIFF
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL ACWK(iw,1,ERROR,*9999)
          IF(DRAW_SCALE) THEN
            IF(NTSCAL.EQ.0) NTSCAL=1
            CALL SGSCAL(INDEX,ISEG,ISSCAL(iw,NTSCAL),iw,NTSCAL,COLOUR,
     '        CSEG,ERROR,*9999)
          ENDIF
          CALL DAWK(iw,1,ERROR,*9999)
C*** Reset the file to the begining of the time data
          CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA,
     '      nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),
     '      YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',FILEFORMAT,HISTORYFNAME,
     '      'RESET',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)
          DO notime=1,NUMTIMEDATA
            IF(DOP) THEN
              WRITE(OP_STRING,'('' notime = '',I5)') notime
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '        NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA,
     '        nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),
     '        YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',FILEFORMAT,
     '        HISTORYFNAME,'TIME_DATA',ENDFILE,.TRUE.,YPDATA,YQDATA,
     '        YQSDATA,ERROR,*9999)
            DO nonrlist=1,NRLIST(0)
              nr=NRLIST(nonrlist)
              CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '          NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),
     '          nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
            ENDDO !nonrlist
            CALL ACWK(iw,1,ERROR,*9999)
            DO nolist=1,NELIST(0)
              ne=NELIST(nolist)
              DO njj1=1,3 !geom + fibres + field
                DO njj2=1,NJ_LOC(njj1,0,nr)
                  nj=NJ_LOC(njj1,njj2,nr)
                  CALL ASSERT(nj.LE.12,'>>ERROR: increase size of '
     '              //'NBJ_temp to NJM',ERROR,*9999)
                  IF(nb_type.EQ.0) THEN
                    NBJ_TEMP(nj)=NBJ(nj,ne)
                  ELSE
                    NBJ_TEMP(nj)=nb_type
                  ENDIF
                ENDDO !njj2
              ENDDO !njj1
              nb=NBJ_TEMP(1)
              CALL XPXE(NBJ_TEMP,NKJE(1,1,1,ne),NPF(1,1),
     '          NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
              CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx),NKHE(1,1,1,ne),
     '          NPF(1,1),NPNE(1,1,ne),NRE(ne),NVHE(1,1,1,ne),
     '          NW(ne,1,nx),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '          ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
              CALL SGFIEL(IBT,IDO,INP,ISEG,ISFIEL(iw,ne),iw,
     '          MXI(1,ne),NBH(ID_TYPE,1,ne),NBJ_TEMP,ne,NGAP(1,nb),
     '          NLATNE,NQNE,NQNLAT,NQS,NQXI,CSEG,TYPE,XE,XQ,YG(1,1,ne),
     '          YQ(1,1,1,nx),ZE(1,ID_TYPE),ERROR,*9999)
            ENDDO !noelem (ne)
            CALL DAWK(iw,1,ERROR,*9999)
            CALL REFRESH_GRAPHICS(DELAY,ERROR,*9999)
          ENDDO !notime

C*** Close history file
          CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '      NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '      TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',
     '      FILEFORMAT,HISTORYFNAME,' ',ENDFILE,.TRUE.,YPDATA,
     '      YQDATA,YQSDATA,ERROR,*9999)

        ENDIF
      ENDIF

      CALL EXITS('DRFIEL')
      RETURN
 9999 IF(ANIMATE) THEN
        CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '    NPNY(0,1,0,nx),
     '    NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),
     '    NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),
     '    YQS,'CLOSE',FILEFORMAT,HISTORYFNAME,' ',ENDFILE,.TRUE.,YPDATA,
     '    YQDATA,YQSDATA,ERROR1,*1111)
      ENDIF
 1111 CALL ERRORS('DRFIEL',ERROR)
      CALL EXITS('DRFIEL')
      RETURN 1
      END


