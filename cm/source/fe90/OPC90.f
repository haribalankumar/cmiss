      SUBROUTINE OPC90(IBT,NBH,NBJ,NEELEM,
     '  NHE,NHP,NKH,NKHE,NKJE,NLL,NP_INTERFACE,NPF,NPNE,NPNODE,
     '  nr,NRE,NVHE,NVHP,NVJE,NW,nx,NYNE,NYNP,
     '  CE,CURVCORRECT,DET,DL,DRDN,PG,RAD,RD,RG,SE,WG,XA,XE,
     '  XG,XG1,XN,XP,XR,YD,YP,ZA,ZE,ZF,ZG,ZP,ERROR,*)

C#### Subroutine: OPC90
C###  Description:
C###    OPC90 prints solutions ZP(nk,nv,nh,np,nc) and, if requested,
C###    the element Gauss point solutions ZG(nh,ni). It also calculates
C###    the solution at any specified point.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'ktyp100.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NHE(NEM),NHP(NPM),
     '  NKH(NHM,NPM,NCM),NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NLL(12,NEM),NP_INTERFACE(0:NPM,0:3),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),CURVCORRECT(2,2,NNM,NEM),
     '  DET(NBFM,0:NNM,NGM,6),DL(3,NLM),DRDN(NGM),
     '  PG(NSM,NUM,NGM,NBM),RAD(NGM),RD(NGM),RG(NGM),
     '  SE(NSM,NBFM,NEM),WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG(NJM,NUM),XG1(NJM,NUM,NGM),
     '  XN(NJM,NGM),XP(NKM,NVM,NJM,NPM),XR(NJM,NGM),
     '  YD(NHM),YP(NYM,NIYM),ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),
     '  ZF(NSM,NHM),ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,ILOOP,INFO,intscheme,iotype_loc,nbbem,nbqh,nbqhp,
     '  nb1j,nb1jp,nc,ne,ng,nh,nhx,ni,nj,nk,nnmin,
     '  noelem,nonode,NOQUES,np,ns,NU1(0:3)
      REAL*8 ARCLEN,DXIX(3,3),DZDN,FLUX,FLUX1,FLUX2,GL(3,3),
     '  GU(3,3),MINDIST,PGG,PGX,RWG,S,SUM,XN_LOCAL(3),XPFP(3)
      LOGICAL FILEIP,INTERFACE
      DATA NU1/1,2,4,7/

      CALL ENTERS('OPC90',*9999)

      NOQUES=0
      ICHAR=999

      CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NVHP,nx,NYNE,NYNP,YP,
     '  ZA,ZP,ERROR,*9999)
      WRITE(OP_STRING,'(/'' Nodal solutions:''/)')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      nc=1
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        WRITE(OP_STRING,'('' Node'',I3,'' ZP(nk,1,1,'',I3,'',1): '','
     '    //'5E11.3,/(23X,5E11.3))') np,np,(ZP(nk,1,1,np,nc),
     '    nk=1,MAX(NKH(NH_LOC(1,nx),np,1)-KTYP93(1,nr),1))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nhx=2,NHP(np)
         nh=NH_LOC(nhx,nx)
         WRITE(OP_STRING,
     '     '(9X,''ZP(nk,'',I1,'','',I1,'','',I3,'',1): '','
     '     //'5E11.3,/(23X,5E11.3))') nh,1,np,(ZP(nk,1,nh,np,nc),
     '      nk=1,MAX(NKH(nh,np,1)-KTYP93(1,nr),1))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDDO

      WRITE(OP_STRING,'(/'' Nodal fluxes:''/)')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      nc=2
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        IF(NVHP(1,np,2).LE.1) THEN !no corner node
          WRITE(OP_STRING,'('' Node'',I3,'' ZP(nk,1,1,'',I3,'',2): '','
     '      //'5E11.3,/(23X,5E11.3))') np,np,(ZP(nk,1,1,np,nc),
     '      nk=1,MAX(NKH(NH_LOC(1,nx),np,2)-KTYP93(2,nr),1))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(NVHP(1,np,nc).EQ.2) THEN !2d corner or 3d edge
          WRITE(OP_STRING,'('' Node'',I3)') np
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
!FIND ne
          WRITE(OP_STRING,
     '      '('' Corner Node'',12X,'' element '',I3,2X,'
     '      //'''element '',I3)')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(9X,''ZP(nk,1,1,'',I3,'',nc): '','
     '      //'5(E11.3,2X,E11.3),/(31X,5(E11.3,2X,E11.3)))') np,
     '      (ZP(nk,1,1,np,nc),ZP(nk,1,1,np,nc),
     '       nk=1,MAX(NKH(NH_LOC(1,nx),np,2)-KTYP93(2,nr),1))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE !3d corner
          WRITE(OP_STRING,'('' Node'',I3)') np
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
!FIND ne
          WRITE(OP_STRING,
     '      '('' Corner Node'',11X,'' element '',I3,2X,'
     '      //'''element '',I3,2X,''element '',I3)')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(9X,''ZP(nk,1,1,'',I3,'',nc): '','
     '      //'5(E11.3,2X,E11.3,2X,E11.3),'
     '      //'/(31X,5(E11.3,2X,E11.3,2X,E11.3)))')
     '      np,(ZP(nk,2,1,np,nc),ZP(nk,3,1,np,nc),ZP(nk,4,1,np,nc),
     '      nk=1,MAX(NKH(NH_LOC(1,nx),np,2)-KTYP93(2,nr),1))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        DO nhx=2,NHP(np)
          nh=NH_LOC(nhx,nx)
          WRITE(OP_STRING,'('' Node'',I3,'' ZP(nk,2,1,'',I3,'',nc): '','
     '      //'5E11.3,/(23X,5E11.3))') nh,2,np,(ZP(nk,2,nh,np,nc),
     '      nk=1,MAX(NKH(nh,np,2)-KTYP93(2,nr),1))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDDO

      ILOOP=1
      iotype_loc=0
 7000 IF(ILOOP.EQ.1) THEN
        FORMAT='(/$,'' Do you want solution at a specified point [N]?'//
     '  ' '',A)'
        CALL GINOUT(iotype_loc,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &    FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      ELSE IF(ILOOP.EQ.2) THEN
        FORMAT='(/$,'' Do you want solution at a specified point [Y]?'//
     '  ' '',A)'
        CALL GINOUT(iotype_loc,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &    FORMAT,1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      ENDIF
      IF(ADATA(1).EQ.'Y') THEN
        ne=NEELEM(1,nr) ! This identifies an element in the BE region.
        CALL EQTYPE(IBT,NBH,NEELEM,nr,NW,nx,ERROR,*9999)
        FORMAT='(/$,'' The Xj-coords of the point are [0,..]? '','
     '    //'3E12.4)'
        DO nj=1,NJT
          RDEFLT(nj)=0.0D0
        ENDDO
        CALL GINOUT(iotype_loc,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &    FORMAT,NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     &    IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        DO nj=1,NJT
          XPFP(nj)=RDATA(nj)
        ENDDO
        IF(ITYP10(nr).GT.1) XPFP(2)=XPFP(2)*PI/180.0D0
        IF(ITYP10(nr).GT.2) XPFP(3)=XPFP(3)*PI/180.0D0
        WRITE(OP_STRING,
     '    '(/'' Solution at Xj-coordinates: '',3E12.4)')
     '    (XPFP(nj),nj=1,NJT)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C MLB Please leave - domsol2 needs more testing
C        CALL DOMSOL(NBH,NBJ,NEELEM,NHE,
C     '    NHP,NKH,NKHE,NKJE,NLL,NPF,NP_INTERFACE,NPNE,NPNODE,
C     '    nr,NRE,NVHE,NVHP,NVJE,NW,nx,NYNE,NYNP,
C     '    CE,CURVCORRECT,DL,PG,RG,SE,WG,XA,XE,XG1,XP,XPFP,
C     '    YD,YP,ZA,ZE,ZF,ZP,ERROR,*9999)
C        IF(.NOT.COMPLEX) THEN
C          WRITE(OP_STRING,'('' First dependent variable'',E12.4)')
C     '      YD(1)
C          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C        ELSE !Complex solution
C        ENDIF

        CALL DOMSOL2(NBH,NBJ,NEELEM,NHE,NHP,NKH,NKHE,NKJE,NLL,NPF,
     '    NP_INTERFACE,NPNE,NPNODE,nr,NRE,NVHE,NVHP,NVJE,NW,nx,NYNE,
     '    NYNP,CE,CURVCORRECT,DET,DL,DRDN,PG,RAD,RD,RG,SE,WG,XA,XE,XG1,
     '    XN,XP,XPFP,XR,YD,YP,ZA,ZE,ZF,ZP,ERROR,*9999)
        IF(.NOT.COMPLEX) THEN
          WRITE(OP_STRING,'('' First dependent variable'',E12.4)')
     '      YD(1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE !Complex solution
        ENDIF
        IF((ITYP2(nr,nx).EQ.4).AND.(ITYP3(nr,nx).EQ.2)) THEN !Yukawa equation
          WRITE(OP_STRING,'('' Moisture Content'',E12.4)')
     '      YD(1)*DEXP(CE(1,1)*XPFP(NJT))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        ILOOP=2
        GO TO 7000
      ENDIF

C Calculation of the total dimensionless flux for a Boundary element
C problems involving the Yukawa equation.
C NOTE: We have used outward normals (i.e. the normal is INTO the
C cavity)
C The flux desired is the flux FROM the cavity i.e. need to use a
C minus sign
      IF ((ITYP2(nr,nx).EQ.4).AND.(ITYP3(nr,nx).EQ.2)) THEN
        IF(ktyp100.NE.2) THEN
          FORMAT='(/$,'' Do you want the flux from the cavity'//
     '  '[N]? '',A)'
        ELSE IF(KTYP100.EQ.2) THEN !Flow around and into cavities
          FORMAT='(/$,'' Do you want the flux into the cavity'//
     '  '[N]? '',A)'
        ENDIF
        CALL GINOUT(iotype_loc,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &    FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(ADATA(1).EQ.'Y') THEN
          S=CE(1,NEELEM(1,nr))
          FLUX=  0.0D0
          FLUX1= 0.0D0
          FLUX2= 0.0D0
          ARCLEN=0.0D0
          CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NVHP,nx,
     '      NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
C ***       Transfer global solution vector to element solution vector
            CALL DIST(intscheme,NBJ(1,ne),NLL(1,ne),nnmin,
     '        NPNE(1,1,ne),nr,NW(ne,2),DL,MINDIST,XP,XPFP,ERROR,*9999)
            CALL QUADBE(intscheme,NBJ(1,ne),nbbem,ERROR,*9999)
C Changed to make sure calcs are using high order scheme
C            nbqhp=NFBASE(1,NBASEF(NBH(NH_LOC(1,nx),2,ne),nbbem))
C            nb1jp=NFBASE(1,NBASEF(NBJ(1,ne),nbbem))
C            !Parent family basis function for geometric (j),
C            !dependent (U) and normal derivative (q).
C            nbqh=NBASEF(NBH(NH_LOC(1,nx),2,ne),nbbem)
C            nb1j=NBASEF(NBJ(1,ne),nbbem)
C
CC AJPs - 191297 - rgb
c            nbqhp=NBJ(1,ne)
c            nb1jp=NBJ(1,ne)
c            nbqh=NBJ(1,ne)
c            nb1j=NBJ(1,ne)
            nbqhp=NBH(1,2,ne)
            nb1jp=NBJ(1,ne)
            nbqh=NBH(1,2,ne)
            nb1j=NBJ(1,ne)
CC AJPe
            CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,2),nx,
     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '        ERROR,*9999)
            CALL ZPZE(NBH(1,1,ne),2,NHE(ne),NKHE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,2),nx,
     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZF,ZP,
     '        ERROR,*9999)
            !ZE(ns,nhx) contains the dependent variable values
            !ZF(ns,nh) contains the normal derivative values
            DO ng=1,NGT(nb1j)
              CALL XEXG(NBJ(1,ne),ng,nr,PG,XE,XG,
     '          ERROR,*9999)
              CALL XGMG(0,NIT(nb1jp),nb1jp,nr,DXIX,GL,GU,RG(ng),XG,
     '          ERROR,*9999)
              !Transfer dependent variable to ZG array
              CALL ZEZG(1,NBH(1,1,ne),ng,NHE(ne),nx,DXIX,PG,ZE,ZG,
     '          ERROR,*9999)
              RWG=RG(ng)*WG(ng,nb1j)
              INTERFACE=.FALSE. !May need changing AJP 21-1-93
              CALL NORMAL(ne,nr,NW,XG,XN_LOCAL,INTERFACE,ERROR,*9999)
              IF(ITYP10(nr).EQ.1) THEN
                IF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN
                  DZDN=XN_LOCAL(2)/DSQRT(XN_LOCAL(1)*XN_LOCAL(1)+
     '              XN_LOCAL(2)*XN_LOCAL(2))
                ELSE
                  DZDN=XN_LOCAL(3)/DSQRT(XN_LOCAL(1)*XN_LOCAL(1)+
     '              XN_LOCAL(2)*XN_LOCAL(2)+XN_LOCAL(3)*XN_LOCAL(3))
                ENDIF
              ELSE
                DZDN=0.0D0 !Needs fixing
              ENDIF
              IF(ktyp100.NE.2) THEN
                FLUX1=FLUX1+S*ZG(1,1)*DEXP(S*XG(2,1))*DZDN*RWG
              ELSE IF(KTYP100.EQ.2) THEN !Flow around and into cavities
                FLUX1=FLUX1+S*(2.0D0+ZG(1,1)*DEXP(S*XG(2,1)))*DZDN*RWG
              ENDIF
              !Transfer normal derivative to ZG array
              DO ni=0,NIT(nbqhp) !(Code extracted from ZEZG).
                SUM=0.0D0
                DO ns=1,NST(nbqhp)
                  IF(ni.EQ.0) THEN
                    PGG=PG(ns,NU1(ni),ng,nbqh)
                  ELSE
                    IF(NIT(nbqhp).EQ.2) THEN
                      PGG=PG(ns,NU1(1),ng,nbqh)*DXIX(1,ni) +
     '                    PG(ns,NU1(2),ng,nbqh)*DXIX(2,ni)
                    ELSE IF(NIT(nbqhp).EQ.3) THEN
                      PGG=PG(ns,NU1(1),ng,nbqh)*DXIX(1,ni) +
     '                    PG(ns,NU1(2),ng,nbqh)*DXIX(2,ni) +
     '                    PG(ns,NU1(3),ng,nbqh)*DXIX(3,ni)
                    ELSE
                      PGG=PGX(nbqh,ni,ns,DXIX,PG(1,1,ng,nbqh))
                    ENDIF
                  ENDIF
                  SUM=SUM+PGG*ZF(ns,NU1(ni))
                ENDDO
CC AJPs 191297 - rgb
C                ZF(ns,NU1(ni))=SUM
                ZG(ns,NU1(ni))=SUM
CC AJPe
              ENDDO
              FLUX2=FLUX2-ZF(1,1)*DEXP(S*XG(2,1))*RWG
              ARCLEN=ARCLEN+RWG
            ENDDO ! End of Gauss Point loop
          ENDDO ! End of loop of elements
          IF(ktyp100.NE.2) THEN !Calculate the flux FROM the cavity
            FLUX=-(FLUX1+FLUX2)
          ENDIF
          WRITE(OP_STRING,'(/'' Dimensionless Flux:'',E14.6)') FLUX
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          FLUX=FLUX/ARCLEN
          WRITE(OP_STRING,
     '      '(/'' Mean Dimensional Flux:'',E14.6)') FLUX
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(/'' Surface area of cavity :'',E14.6)')
     '      ARCLEN
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF ! End of choice to get flux calculated
      ENDIF ! End of ITYP2(nr,nx) and ITYP3(nr,nx) if loop

      CALL EXITS('OPC90')
      RETURN
 9999 CALL ERRORS('OPC90',ERROR)
      CALL EXITS('OPC90')
      RETURN 1
      END


