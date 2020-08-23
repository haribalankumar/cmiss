      SUBROUTINE DRLINE(IBT,ISEG,ISLINE,ISLINO,NBH,NBJ,NEELEM,NELIST,
     '  NHE,NHP,NKH,NKJ,NLL,NLLINE,NLLIST,NONY,NPL,NPNODE,NPNY,NRLIST,
     '  NVHP,NW,NXLIST,NYNE,NYNP,NYNR,CONY,DL,EIGVEC,XP,YP,ZA,ZP,
     '  CSEG,STRING,ERROR,*)

C#### Subroutine: DRLINE
C###  Description:
C###    DRLINE draws lines.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'eige00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'moti00.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'time00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),ISEG(*),
     '  ISLINE(NWM,2*NGRSEGM),ISLINO(NWM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),NKJ(NJM,NPM),NLL(12,NEM),
     '  NLLINE(0:NL_R_M,0:NRM),NLLIST(0:NLM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),NPL(5,0:3,NLM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM,NXM),
     '  NRLIST(0:NRM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NW(NEM,3,NXM),NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  DL(3,NLM),EIGVEC(NOM,NTM,2),XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,IFROMC,INDEX,INDEX_POLYLINE,iw,
     '  IWK(6),MODE_NUMBER,N3CO,nae,nb,nc,ne,nh,nhx,nj,
     '  nk,NKEF,nks,nl,no,no_coeff,noiw,nolist,no_nrlist,nonode,
     '  no_nynr,noy,np,nr,NTIW,nx,nxc,ny
      REAL*8 co1,DISPLACEMENT,PF1,RFROMC,SCALE_FACTOR
      CHARACTER DEFORM*10,SOLID*8,WIDTH*6
      LOGICAL ALL_REGIONS,CBBREV,MODE,STATIC,Use_line_list

      CALL ENTERS('DRLINE',*9999)
      nc=1 ! Temporary MPN 12-Nov-94
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM draw lines
C###  Description:
C###    Draw the element lines.  The lines are either solid or dotted.
C###  Parameter:    <with (all/LINE#s)[all]>
C###    Specify the line numbers.
C###  Parameter:    <in (all/ELEMENT#s)[all]>
C###    Specify the element numbers.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s).
C###    The all keyword specifies all currently defined workstations.
C###  Parameter:    <(undeformed/deformed/field)[undeformed]>
C###    Specify whether to use the geometry of the 'undeformed' mesh,
C###    'deformed' mesh, or 'field' variables.
C###  Parameter:    <when TIME>
C###    For non-static deformed meshes, specify the time for the mesh.
C###  Parameter:    <scale FACTOR[1.0]>
C###    For static deformed meshes, multiply the displacements by
C###    FACTOR.
C###  Parameter:    <(solid/dotted)[solid]>
C###    Specify the appearance of the lines.
C###  Parameter:    <width (1/2/3/4)[1]>
C###    Specify the pixel width.
C###  Parameter:    <rgb=RGB[red]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING( 1)=STRING(1:IEND)//' <with (all/LINE#s)[all]>'
        OP_STRING( 2)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING( 3)=BLANK(1:15)//'<in (all/ELEMENT#s)[all]>'
        OP_STRING( 4)=BLANK(1:15)
     '    //'<(undeformed/deformed/field)[undeformed]>'
        OP_STRING( 5)=BLANK(1:15)//'<when TIME>'
        OP_STRING( 6)=BLANK(1:15)//'<scale FACTOR[1.0]>'
        OP_STRING( 7)=BLANK(1:15)//'<(solid/dotted)[solid]>'
        OP_STRING( 8)=BLANK(1:15)//'<width (1/2/3/4)[1]>'
        OP_STRING( 9)=BLANK(1:15)//'<rgb=RGB[red]>'
        OP_STRING(10)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(11)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw lines mode MODE#
C###  Description:
C###    Draw the lines of the mesh deformed under mode MODE#.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s).
C###    The all keyword specifies all currently defined workstations.
C###  Parameter:    <scale FACTOR#[1.0]>
C###    Multiply the displacements by FACTOR.
C###  Parameter:    <(solid/dotted)[solid]>
C###    Specify the appearance of the lines.
C###  Parameter:    <width (1/2/3/4)[1]>
C###    Specify the pixel width.
C###  Parameter:    <rgb=RGB[red]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING( 1)=STRING(1:IEND)//' mode MODE#'
        OP_STRING( 2)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING( 3)=BLANK(1:15)//'<scale FACTOR[1.0]>'
        OP_STRING( 4)=BLANK(1:15)//'<(solid/dotted)[solid]>'
        OP_STRING( 5)=BLANK(1:15)//'<width (1/2/3/4)[1]>'
        OP_STRING( 6)=BLANK(1:15)//'<rgb=RGB[red]>'
        OP_STRING( 7)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING( 8)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRLINE',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL ASSERT(CALL_LINE,'>>no lines defined',ERROR,*9999)

        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        Use_line_list=.FALSE.
        IF(CBBREV(CO,'WITH',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NL_R_M,NLLIST(0),NLLIST(1),ERROR,*9999)
        ELSE IF(CBBREV(CO,'IN',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NE_R_M,NELIST(0),NELIST(1),ERROR,*9999)
          NLLIST(0)=0
          DO nolist=1,NELIST(0)
            ne=NELIST(nolist)
            DO nae=1,NLE(NBJ(1,ne))
              nl=NLL(nae,ne)
              IF(nl.NE.0) THEN
                NLLIST(0)=NLLIST(0)+1
                NLLIST(NLLIST(0))=nl
              ENDIF
            ENDDO
          ENDDO
        ELSE IF(CBBREV(CO,'MODE',1,noco+1,NTCO,N3CO)) THEN
          MODE=.TRUE.
          MODE_NUMBER=IFROMC(CO(N3CO+1))
        ELSE
          MODE=.FALSE.
          Use_line_list=.TRUE. !get lines directly from NLLINE
        ENDIF

        IF(CBBREV(CO,'UNDEFORMED',1,noco+1,NTCO,N3CO)) THEN
          DEFORM='UNDEFORMED'
        ELSE IF(CBBREV(CO,'FIELD',1,noco+1,NTCO,N3CO)) THEN
          DEFORM='FIELD'
        ELSE IF(MODE.OR.CBBREV(CO,'DEFORMED',1,noco+1,NTCO,N3CO)) THEN
          DEFORM='DEFORMED'
c old MPN 14-May-96: handled below
c          IF(FIRST_DEF_LINE) THEN !1st drawing of a deformed line
c            ADD=.TRUE.            !so create new line segment
c            FIRST_DEF_LINE=.FALSE.
c          ELSE
c            ADD=.FALSE.
c          ENDIF
        ELSE
          DEFORM='UNDEFORMED'
        ENDIF

        IF(MODE.OR.DEFORM(1:8).EQ.'DEFORMED') THEN
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)
        ENDIF

        IF(MODE) THEN
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            CALL ASSERT(ITYP5(nr,nx).EQ.3,'>>Define modal analysis '
     '        //'first',ERROR,*9999)
          ENDDO !no_nrlist (nr)
          CALL ASSERT(MODE_NUMBER.GE.1.AND.MODE_NUMBER.LE.NUMEIGEN,
     '      '>>Invalid mode number',ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'WHEN',2,noco+1,NTCO,N3CO)) THEN
          STATIC=.FALSE.
          TIME=RFROMC(CO(N3CO+1))
        ELSE
          STATIC=.TRUE.
          TIME=0.0d0
        ENDIF

        IF(CBBREV(CO,'SCALE',1,noco+1,NTCO,N3CO)) THEN
          SCALE_FACTOR=RFROMC(CO(N3CO+1))
        ELSE
          SCALE_FACTOR=1.0d0
        ENDIF

        IF(CBBREV(CO,'SOLID',1,noco+1,NTCO,N3CO)) THEN
          SOLID='SOLID'
        ELSE IF(CBBREV(CO,'DOTTED',1,noco+1,NTCO,N3CO)) THEN
          SOLID='DOTTED'
        ELSE
          SOLID='SOLID'
        ENDIF

        IF(CBBREV(CO,'WIDTH',1,noco+1,NTCO,N3CO)) THEN
          WIDTH='WIDTH'//CO(N3CO+1)(1:1)
        ELSE
          WIDTH='WIDTH1'
        ENDIF

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,SOLID,WIDTH,CO(N3CO+1))
        ELSE
          IF(WIDTH(1:6).EQ.'WIDTH1'.AND.SOLID.EQ.'SOLID') THEN
            INDEX=INDEX_POLYLINE(0,SOLID,WIDTH,'RED')
          ELSE
            INDEX=INDEX_POLYLINE(0,SOLID,WIDTH,'BLACK')
          ENDIF
        ENDIF


        DO no_nrlist=1,NRLIST(0) !loop over regions
          nr=NRLIST(no_nrlist)

C new MPN 14-May-96: adding for deformed case
          NTLINE=no_nrlist
          IF(MODE.OR.DEFORM(1:8).EQ.'DEFORMED') THEN
            NTLINE=NTLINE+NGRSEGM
          ENDIF
          IF(ADD) NTLINE=NTLINE+1
          CALL ASSERT(NTLINE.LE.2*NGRSEGM,'>>Increase NGRSEGM',
     '      ERROR,*9999)
C old
C          IF(ADD) THEN
C            NTLINE=NTLINE+1
C          ELSE
C            NTLINE=no_nrlist
C          ENDIF

          IF(Use_line_list) THEN
            NLLIST(0)=NLLINE(0,nr)
            DO nolist=1,NLLIST(0) !use lines in region nr
              NLLIST(nolist)=NLLINE(nolist,nr)
            ENDDO
          ENDIF !Use_line_list
          DO noiw=1,NTIW
            IW=IWK(noiw)
            CALL ACWK(iw,1,ERROR,*9999)

            IF(.NOT.MODE.AND.(DEFORM(1:10).EQ.'UNDEFORMED'.OR.
     '        DEFORM(1:5).EQ.'FIELD').AND.STATIC) THEN
c PJH 11Oct94 need new segment for each region (temporary?)
c fixed MPN 15Feb96
              CALL SGLINE(INDEX,ISEG,ISLINE(iw,NTLINE),ISLINO(iw),iw,
     '          NLLIST,NTLINE,NPL,nr,nx,CSEG,DEFORM,DL,SOLID,STATIC,
     '          XP,ZP,ERROR,*9999)
c              CALL SGLINE(INDEX,ISEG,ISLINE(iw,nr),ISLINO(iw),iw,
c     '          NLLIST,NTLINE,NPL,nr,nx,CSEG,DEFORM,DL,SOLID,STATIC,
c     '          XP,ZP,ERROR,*9999)

            ELSE IF(MODE.OR.(DEFORM(1:8).EQ.'DEFORMED').OR.
     '          (.NOT.STATIC)) THEN
              IF(STATIC) THEN
                IF(MODE) THEN
                  DO no_nynr=1,NYNR(0,0,1,nr,nx)
                    ny=NYNR(no_nynr,0,1,nr,nx) !is global var number
                    IF(NPNY(0,ny,0,nx).EQ.1) THEN
                      np=NPNY(4,ny,0,nx)
C GMH 8/1/97 Update cmgui link
                      CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                    ENDIF
                    DO noy=1,NONY(0,ny,2,nr,nx)
                      no=NONY(noy,ny,2,nr,nx)
                      co1=CONY(noy,ny,2,nr,nx)
                      YP(ny,1,nx)=EIGVEC(no,MODE_NUMBER,1)*co1
                    ENDDO !noy
                  ENDDO !no_nynr (ny)
                ENDIF
                IF(ITYP6(nr,nx).EQ.1..AND.ITYP1(nr,nx).EQ.4.
     '            AND.NW(NEELEM(1,nr),1,nx).EQ.3) THEN
                  nb=NBJ(2,NEELEM(1,nr))         !beam element
                  CALL ASSERT(NKT(0,nb).EQ.2,
     '              '>>y-coord s.b. cubic Hermite',ERROR,*9999)
                ENDIF
                CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '            NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,
     '            NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
              ELSE !motion
                CALL YPZP(MOTION_IY,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '            NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '            YP(1,1,nx),ZA,ZP,ERROR,*9999)
              ENDIF
              IF(STATIC) THEN
                nc=1 !Temporary AJP 19-12-91
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                    nj=NJ_LOC(NJL_GEOM,nhx,nr)
                    nh=NH_LOC(nhx,nx)
                    DO nk=1,NKH(nh,np,nc,nr)  !note nkh AAY
                      IF(ITYP6(nr,nx).EQ.1.OR.KTYP58(nr).EQ.2) THEN
                        DISPLACEMENT=ZP(nk,1,nh,np,nc) !linear analysis
                      ELSE IF(ITYP6(nr,nx).EQ.2) THEN  !non-linear analy
                        DISPLACEMENT=ZP(nk,1,nh,np,nc)-XP(nk,1,nj,np)
                      ENDIF
                      IF(nk.LE.NKJ(nj,np)) THEN
                        ZP(nk,1,nh,np,nc)=XP(nk,1,nj,np)+
     '                    SCALE_FACTOR*DISPLACEMENT
                      ELSE
                        ZP(nk,1,nh,np,nc)=SCALE_FACTOR*DISPLACEMENT
                      ENDIF
                    ENDDO !nk
                  ENDDO !nhx
                ENDDO !nonode

              ELSE IF(.NOT.STATIC) THEN !interpolate ZP in time first
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                    nh=NH_LOC(nhx,nx)
                    nb=NBH(nh,1,NEELEM(1,nr)) !Fourier basis number
                    DO nks=1,NKJ(nh,np)      !spatial dof number
C                     interpolate node np in time
                      DISPLACEMENT=0.0d0
                      DO no_coeff=1,IBT(2,NIT(nb),nb)
                        NKEF=(no_coeff-1)*NKJ(nh,np)+nks !Fourier dof#
                        DISPLACEMENT=DISPLACEMENT+
     '                    PF1(no_coeff,1,TIME)*ZP(NKEF,1,nh,np,nc)
                      ENDDO !no_coeff
                      IF(KTYP58(nr).EQ.2) THEN
                        ZP(nks,1,nh,np,nc)=XP(nks,1,nh,np)+
     '                    DISPLACEMENT
                      ELSE
                        ZP(nks,1,nh,np,nc)=DISPLACEMENT
                      ENDIF
                    ENDDO !nks
                  ENDDO !nhx
                ENDDO !nonode
              ENDIF

              CALL SGLINE(INDEX,ISEG,ISLINE(iw,NTLINE),ISLINO(iw),iw,
     '          NLLIST,NTLINE,NPL,nr,nx,CSEG,DEFORM,DL,SOLID,
     '          STATIC,XP,ZP,ERROR,*9999)

            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('DRLINE')
      RETURN
 9999 CALL ERRORS('DRLINE',ERROR)
      CALL EXITS('DRLINE')
      RETURN 1
      END


