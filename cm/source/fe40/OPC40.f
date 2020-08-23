      SUBROUTINE OPC40(IBT,IDO,INP,NBH,NBJ,NEELEM,NHE,NHP,NKJE,
     '  NKH,NPF,NPNE,NPNODE,nr,NRE,NVHE,NVHP,NVJE,NW,nx,NYNE,NYNP,CE,
     '  CG,CGE,CP,PG,RG,SE,XA,XE,XG,XP,YP,ZA,ZE,ZP,ERROR,*)

C#### Subroutine: OPC40
C###  Description:
C###    OPC40 prints solutions ZP(nk,nv,nh,np,nc) and ZA(na,nh,nc,ne)
C###    and nodal reactions and, if requested, the element Gauss point
C###    solutions ZG(nh,ni) and the solution at any specified point
C###    in theta coordinates.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b40.cmn'
      INCLUDE 'dx00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NHE(NEM),
     '  NHP(NPM),NKJE(NKM,NNM,NJM,NEM),
     '  NKH(NHM,NPM,NCM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),
     '  CP(NMM,NPM),PG(NSM,NUM,NGM,NBM),
     '  RG(NGM),SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,IREC,j,k,KTOT,nb,nc,ne,ng,nh,nhx,ni,nj,nk,
     '  noelem,nonode,np,nu,nv
      REAL*8 DXIX(3,3),GGM(2,2),GGPX(3,3),GL(3,3),GU(3,3),POIS,PXI,
     '  RRM(2,2),SMX(2,2),SNX(2,2),SSPX(3,3),THIC,TOL_XCOORD,
     '  X3G(4,3,25),XB_LOCAL(10),XI(3),XI3,XPOINT(3),YMOD,
     '  ZGG(10,3)
      CHARACTER ANS

      CALL ENTERS('OPC40',*9999)
      nc=1 !Temporary MPN 12-Nov-94
      nv=1 !Temporary cpb 13/11/94
      nb=1 !Temporary gbs 19/09/95
C GMH 21-Apr-95      nx=1 !Temporary no longer
      TOL_XCOORD=1.0d-06 !MHT 18/11/96 for call to XCOORD

      WRITE(OP_STRING,'(/'' ******OPC40'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      KTOT=KTYP17
C***  temporary
      IF(KTOT.GT.16) KTOT=16
      IF(KTOT.GT.0) THEN
C IOFILE2 replaces unit 3 IOFILE3 replaces unit 4 23/11/92
C 1010     WRITE(3,'(/'' Do you want eigenvector output (Y/N)?'')')
C 101      FORMAT(/' Do you want eigenvector output (Y/N)?')
C          READ(3,501,IOSTAT=IO,END=1011) CDATA
C 1011     CALL DATYP(CDATA,IO,1,40,ANS,LDATA,0,0.,0,0,*1010)
C 1011     CONTINUE
C          INQUIRE(UNIT=1,NEXTREC=IREC)
C          WRITE(1,1012,REC=IREC) ANS
C 1012     FORMAT(/' Do you want eigenvector output (Y/N)?  ',A)
C        IF(ANS.EQ.'N') KTOT=0
      ENDIF
      DO 80 k=1,KTOT
      IF(ITYP5(nr,nx).EQ.3.OR.ITYP5(nr,nx).EQ.5) THEN
        WRITE(OP_STRING,'(//'' Eigenvector no.'',I2,/1X,17(''=''))')
     '    k
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF
      WRITE(OP_STRING,'(/'' Nodal solutions:''/)')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      CALL YPZP(k,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '  nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        WRITE(OP_STRING,'('' Node'',I3,'' ZP(nk,1,1,'',I3,'',1): '','
     '    //'5E16.8,/(23X,5E16.8))')
     '    np,np,(ZP(nk,nv,1,np,nc),nk=1,NKH(NH_LOC(1,nx),np,nc))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nhx=2,NHP(np)
          nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
          WRITE(OP_STRING,
     '      '(9X,''ZP(nk,'',I1,'','',I1,'','',I3,'',1): '','
     '      //'5E16.8,/(23X,5E16.8))')
     '      nh,nv,np,(ZP(nk,nv,nh,np,nc),nk=1,NKH(nh,np,nc))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDDO

C:    DO 30 ne=1,NET
C:      DO 25 nvar=1,NVE(NW(ne,1))
C:        nh=NHV(nvar,NW(ne,1))
C:        nb=NBH(nh,nc,ne)
C:        IF(NAT(nb).GT.0) THEN
C:          WRITE(IO4,250)         nh,ne,(ZA(na,nh,nc,ne),NA=1,NAT(nb))
C250        FORMAT(' ZA(na,1,',I1,',',I2,'): ',6E11.4,/(14X,6E11.4))
C:        ENDIF
C25     CONTINUE
C30   CONTINUE

      IF(ITYP5(nr,nx).EQ.1) THEN
        WRITE(OP_STRING,'(/'' Nodal reactions:''/)')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '    nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
        nc=2 !Temporary AJP 18-12-91
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          WRITE(OP_STRING,'('' Node'',I3,'' ZP(nk,1,1,'',I3,'',1): '','
     '      //'5E16.8,/(23X,5E16.8))')
     '      np,np,(ZP(nk,nv,1,np,nc),nk=1,NKH(NH_LOC(1,nx),np,nc))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nhx=2,NHP(np)
            nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
            WRITE(OP_STRING,
     '        '(9X,''ZP(nk,'',I1,'','',I1,'','',I3,'',1): '','
     '        //'5E16.8,/(23X,5E16.8))')
     '        nh,nv,np,(ZP(nk,nv,nh,np,nc),nk=1,NKH(nh,np,nc))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
      ENDIF

      CALL YPZP(k,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '  nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)

      IF(ITYP5(nr,nx).LE.2) THEN
C IOFILE2 replaces unit 3
C MLB 16/4/97 initialising ANS
        ANS='N'

C 5000   CONTINUE
        WRITE(3,
     '    '(/'' Do you want solution printed out at all element Gauss'
     '    //' points (Y/N)?'')')
C 500    FORMAT(/' Do you want solution printed out at all element Gauss'
C     '         ,' points (Y/N)?')
c        READ(3,501,IOSTAT=IO,END=5001) CDATA
c 501    FORMAT(A)
C 5001   CONTINUE
C        CALL DATYP(CDATA,IO,1,40,ANS,LDATA,0,0.0d0,0,0,*5000)
        INQUIRE(UNIT=1,NEXTREC=IREC)
        WRITE(1,502,REC=IREC) ANS
 502    FORMAT(/' Do you want solution printed out at all element Gauss'
     '         ,' points (Y/N)?  ',A)
      IF(ANS.EQ.'Y') THEN
        WRITE(OP_STRING,'(/'' Gauss point solutions:''/)')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO 60 noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          WRITE(OP_STRING,'(/'' Element '',I2)') ne
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
C ***     CALL ZPZE(NBH(1,nc,ne),nc,
C ***'      NHE(ne),NPNE(1,1,ne),NPF(1,1),nr,NVHE(1,1,1,ne),NW(ne,1),
C ***'      CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
C ***'      ERROR,*9999)
          CALL ZPZE20(NBH(1,nc,ne),NBJ(1,ne),NHE(ne),NPNE(1,1,ne),
     '      NVHE(1,1,1,ne),nx,SE(1,1,ne),ZE,ZP,ERROR,*9999)

C      Rotate the element parameters from global to local system
C<<<
C>>>      CALL ZERT(NBH(1,nc,ne),NBJ(1,ne),NHE(ne),XE,ZE)
*         DO 55 nh=1,NHE(ne)
*           nb=NBH(nh,nc,ne)
*           IF(DOP) THEN
*             WRITE(IO4,5080)         nh,(ZE(ns,nh),ns=1,NST(nb))
*5080         FORMAT(//' ZE(ns,',I1,'): ',6E11.3,(/11X,6E11.3),
*    '              (/11X,6E11.3))
*           ENDIF
*           DO 54 ng=1,NGT(nb)
C             CALL ZEZG(1,NBH(1,nc,ne),ng,NHE(ne),nx,DXIX(1,1,ng),
C ***'          PG,ZE,ZG,ERROR,*9999)
C             WRITE(IO4,519)         NG,nh,(ZGU(nu,nh),nu=1,6)
C519          FORMAT(' Gauss point',I2,' ZGG(nu,',I1,'): ',6E12.4)
C ***         WRITE(IO4,520)         NG,nh,(ZG(nh,ni),NI=0,NIT(nb))
C520          FORMAT(' Gauss point',I2,' ZG(',I1,',ni): ',4E12.4)
*54         CONTINUE
*55       CONTINUE
            DO 540 ng=1,NGT(nb)
              CALL XEXG(NBJ(1,ne),ng,NRE(ne),PG,
     '          XE,XG,ERROR,*9999)
              CALL XGMG(0,NIT(1),NBJ(1,ne),NRE(ne),DXIX,
     '          GL,GU,RG(ng),XG,ERROR,*9999)
              CALL ZEZG20(NBH(1,nc,ne),ng,NHE(ne),nx,PG,ZE,ZGG,
     '          ERROR,*9999)
C             CALL ZGZLG(NBH(NH_LOC(1,nx),nc,ne),NHE(ne),GL,GU,ZGG)
C             CALL ZLZUG(NBH(NH_LOC(1,nx),nc,ne),NHE(ne),GL,GU,ZGG,ZGU)
          nc=1 !Temporary AJP 18-12-91
          DO 550 nhx=1,NHE(ne)
            nh=NH_LOC(nhx,nx)
            nb=NBH(nh,nc,ne)
            WRITE(OP_STRING,'('' Dependent variable nh= '',I1)') nhx
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C              WRITE(OP_STRING,'('' Gauss point'',I2,'
C     '          //''' ZGU(nu,'',I1,''): '',6E12.4)')
C     '          ng,nh,(ZGU(nu,nh),nu=1,6)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
 550      CONTINUE
 540        CONTINUE
 60     CONTINUE
      ENDIF

 65   CONTINUE
      NJT=NJ_LOC(NJL_GEOM,0,nr)
C IOFILE2 replaces unit 3 23/11/92
C 7000 CONTINUE
      WRITE(3,
     '    '(/'' Do you want solution at a specified point (Y/N)?'')')
C 700    FORMAT(/' Do you want solution at a specified point (Y/N)?')
c        READ(3,501,IOSTAT=IO,END=7001) CDATA
C 7001   CONTINUE
C       CALL DATYP(CDATA,IO,1,10,ANS,LDATA,0,0.0d0,0,0,*7000)
        INQUIRE(UNIT=1,NEXTREC=IREC)
        WRITE(1,702,REC=IREC) ANS
 702    FORMAT(/' Do you want solution at a specified point (Y/N)?  ',
     '          A)
        IF(ANS.EQ.'Y') THEN
C IOFILE2 replaces unit 3 23/11/92
          WRITE(3,'(/'' The Xj-coordinates of the point are?'')')
C 704      FORMAT(/' The Xj-coordinates of the point are?')
          READ(3,*)    (XPOINT(nj),nj=1,NJT)
          INQUIRE(UNIT=1,NEXTREC=IREC)
          WRITE(1,706,REC=IREC) (XPOINT(nj),NJ=1,NJT)
 706      FORMAT(/' The Xj-coordinates of the point are?'/,1X,3E12.4)
        ENDIF
      IF(ANS.EQ.'Y') THEN
        IF(ITYP10(1).GT.1) XPOINT(2)=XPOINT(2)*PI/180.d0
        IF(ITYP10(1).GT.2) XPOINT(3)=XPOINT(3)*PI/180.d0
          WRITE(OP_STRING,
     '      '(/'' Solution at Xj-coordinates: '',3E12.4)')
     '      (XPOINT(nj),nj=1,NJT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          CALL XCOORD(IBT,IDO,INP,NBJ,ne,NEELEM,
     '      NKJE,NPF,NPNE,NPNODE,nr,NVJE,
     '      SE,TOL_XCOORD,XA,XE,XI,XP,XPOINT,ERROR,*9999)
        nb=NBJ(1,ne)
        WRITE(OP_STRING,'('' The point is in element '',I3,'
     '    //''' at Xi-coords: '',3E11.3)')
     '    ne,(XI(ni),ni=1,NIT(nb))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C ***   CALL ZPZE(NBH(1,nc,ne),nc,
C ***'    NHE(ne),NPF(1,1),NPNE(1,1,ne,NVHE(1,1,1,ne),
C ***'    CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
C ***'    ZE,ZP,ERROR,*9999)
        CALL ZPZE20(NBH(1,nc,ne),NBJ(1,ne),NHE(ne),NPNE(1,1,ne),
     '    NVHE(1,1,1,ne),nx,SE(1,1,ne),ZE,ZP,ERROR,*9999)
        nc=1 !Temporary AJP 18-12-91
        DO 75 nhx=1,NHE(ne)
          nh=NH_LOC(nhx,nx)
          nb=NBH(nh,nc,ne)
          DO 74 nu=1,NUT(nb)
            XB_LOCAL(nu)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                 nb,nu,XI,ZE(1,nh))
 74       CONTINUE
          WRITE(OP_STRING,
     '      '(/'' Xi-derivatives (up to 2nd order) of dependent'
     '      //' variable '' ,I1,'': ''//,(1X,6E12.4))')
     '      nhx,(XB_LOCAL(nu),nu=1,NUT(nb))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
 75     CONTINUE
        CALL CPCG(NW(ne,1),NBJ(1,ne),NPNE(1,1,ne),nr,nx,
     '    CE(1,ne),CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '    NRE(ne),NVJE(1,1,1,ne),
     '    SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
C       CALL STRESS(IBT,IDO,INP,NBH(1,nc,ne),NBJ(1,ne),NHE(ne),
C    '    NKE(1,1,1,ne),CG(1,1),XE,XI,ZE)
        GO TO 65
      ENDIF
C
C Computation of the stresses and strains for the shell case (NW=7).
C They are computed at each gauss point of each element.
C
C IOFILE2 replaces unit 3 23/11/92
C 8000 CONTINUE
      WRITE(3,
     '    '(/'' Do you want stresses/strains at gauss points (Y/N)?'')')
C 800    FORMAT(/' Do you want stresses/strains at gauss points (Y/N)?')
c        READ(3,501,IOSTAT=IO,END=8001) CDATA
C 8001   CONTINUE
C        CALL DATYP(CDATA,IO,1,10,ANS,LDATA,0,0.0d0,0,0,*8000)
        INQUIRE(UNIT=1,NEXTREC=IREC)
        WRITE(1,802,REC=IREC) ANS
 802    FORMAT(/' Do you want stresses/strains at gauss points (Y/N)? ',
     '          A)
      IF(ANS.EQ.'Y') THEN
C IOFILE2 replaces unit 3 23/11/92
      WRITE(3,'(/,'' Stresses and Strains (in local coordinates)'','
     '  //'/,'' ==========================================='')')
      WRITE(OP_STRING,
     '  '(/,'' Stresses and Strains (in local coordinates)'','
     '  //'/,'' ==========================================='')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO 85 noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        WRITE(OP_STRING,'(/,'' Element Number:'',I3/)')
     '    ne
C IOFILE2 replaces unit 3 23/11/92
        WRITE(3,'(/,'' Element Number:'',I3/)') ne
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '    NRE(ne),NVJE(1,1,1,ne),
     '    SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        CALL CPCG(NW(ne,1),NBJ(1,ne),NPNE(1,1,ne),nr,nx,
     '    CE(1,ne),CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
        CALL ZPZE20(NBH(1,nc,ne),NBJ(1,ne),NHE(ne),NPNE(1,1,ne),
     '    NVHE(1,1,1,ne),nx,SE(1,1,ne),ZE,ZP,ERROR,*9999)
        CALL X3XG(NBJ(1,ne),D3PG,XE,X3G,ERROR,*9999)
C ***   WRITE(3,763) ne
C ***   WRITE(IO4,763) ne
C ***   DO 761 nj=1,NJ_LOC(NJL_GEOM,0,nr)
C ***   DO 761 nd=1,4
C ***     WRITE(3,762)nj,nd,(X3G(nd,nj,ng),ng=1,NGT(2))
C ***     WRITE(IO4,762)nj,nd,(X3G(nd,nj,ng),ng=1,NGT(2))
C761    CONTINUE
C762    FORMAT(' nj=',I2,' nd=',I2,1X,9E10.3)
C763    FORMAT(/,' Element No. ',I2,' X3G(ng,nj,ng),ng=1,2,3,.......',/)
C       WRITE(IO4,*)' Before ZERT***'
C       DO 770 ns=1,NST(NBH(NH_LOC(1,nx),nc,ne))
C         WRITE(IO4,775)ns,(ZE(ns,nh),nh=1,3)
C770    CONTINUE
C<<<
C>>>    CALL ZERT(NBH(1,nc,ne),NBJ(1,ne),NHE(ne),XE,ZE)
C===    CALL ZEPU(NBH(1,nc,ne),NBJ(1,ne),NHE(ne),XE,ZE)
C       WRITE(IO4,*)' After ZERT***'
C       DO 771 ns=1,NST(NBH(NH_LOC(1,nx),nc,ne))
C         WRITE(IO4,775)ns,(ZE(ns,NH_LOC(nhx,nx),nhx=1,3)
C771    CONTINUE
C775    FORMAT(' ns=',I2,' ze(ns,nh)=',3E16.8)
        DO 84 ng=1,NGT(NBH(NH_LOC(1,nx),nc,ne))
          CALL XEXG(NBJ(1,ne),ng,NRE(ne),PG,
     '      XE,XG,ERROR,*9999)
          CALL XGMG(0,NIT(NBJ(1,ne)),NBJ(1,ne),NRE(ne),
     '      DXIX,GL,GU,RG(ng),XG,ERROR,*9999)
          WRITE(OP_STRING,'(/,'' Gauss Point ng='',I3,/)')
     '      ng
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(3,'(/,'' Gauss Point ng='',I3,/)')
     '      ng
C          WRITE(IO4,850)ng
C          WRITE(3,850)ng
C***      DO 83 nh=1,NHE(ne)
C ***       nb=NBH(nh,nc,ne)
C
C  Transfers the ZE array into the Gauss point array ZGG.
C
            CALL ZEZG20(NBH(1,nc,ne),ng,NHE(ne),nx,PG,ZE,ZGG,
     '        ERROR,*9999)
            CALL ZGCHK(NBH(1,nc,ne),NHE(ne),nx,ZGG,ERROR,*9999)
C             CALL ZGZLG(NBH(1,nc,ne),NHE(ne),GL,GU,ZGG)
C83       CONTINUE
C         nc=1 !Temporary AJP 18-12-91
C         DO 832 nh=1,NHE(ne)
C           WRITE(3,833)nh,(ZGG(ns,nh),ns=1,NUT(NBH(nh,nc,ne)))
C           WRITE(IO4,833)nh,(ZGG(ns,nh),ns=1,NUT(NBH(nh,nc,ne)))
C832      CONTINUE
C833      FORMAT(' nh=',I2,' zgg(nu,nh)= ',10E11.3,/)

C  Determines the middle surface Strain and Change in Curvature
C  tensors

          CALL GGRRM(NBJ(1,ne),GGM,GU,RRM,
     '      X3G(1,1,ng),XG,ZGG,ERROR,*9999)
          YMOD=CG(1,ng)
          POIS=CG(2,ng)
          THIC=CG(NMP(NW(ne,1))+1,ng)
C Internal surface......
          XI3=-0.50d0*THIC
C External surface
C         XI3= 0.50d0*THIC
          CALL SSGGX(NBJ(1,ne),GGM,GGPX,GL,GU,
     '      POIS,RRM,SSPX,X3G(1,1,ng),XG,XI3,YMOD,ERROR,*9999)
          CALL MIDMN(NBJ(1,ne),GGM,GL,GU,
     '      POIS,RRM,SMX,SNX,THIC,YMOD,ERROR,*9999)
          DO i=1,3
            WRITE(OP_STRING,
     '        '('' ggpx:'','' i='',I2,3('' j='',I2,E11.3),''  & '
     '        //' sspx:'','' i='',I2,3('' j='',I2,E11.3))')
     '        i,(j,GGPX(i,j),j=1,3),i,(j,SSPX(i,j),j=1,3)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(3,'('' ggpx:'','' i='',I2,3('' j='',I2,E11.3))')
     '        i,(j,GGPX(i,j),j=1,3)
          ENDDO
C 830      CONTINUE
          WRITE(3,*)' '
          DO i=1,3
            WRITE(3,'('' sspx:'','' i='',I2,3('' j='',I2,E11.3))')
     '        i,(j,SSPX(i,j),j=1,3)
          ENDDO
C 831      CONTINUE
          WRITE(3,*)' '
          WRITE(OP_STRING,*)' '
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO i=1,2
            WRITE(OP_STRING,
     '        '('' smx:'','' i='',I2,2('' j='',I2,E11.3),''  & '
     '        //' snx:'','' i='',I2,2('' j='',I2,E11.3))')
     '        i,(j,SMX(i,j),j=1,2),i,(j,SNX(i,j),j=1,2)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(3,'('' smx:'','' i='',I2,2('' j='',I2,E11.3))')
     '        i,(j,SMX(i,j),j=1,2)
          ENDDO
C 833      CONTINUE
          WRITE(3,*)' '
          WRITE(OP_STRING,*)' '
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO i=1,2
            WRITE(3,'('' snx:'','' i='',I2,2('' j='',I2,E11.3))')
     '        i,(j,SNX(i,j),j=1,2)
          ENDDO
C 834      CONTINUE
 84     CONTINUE
 85   CONTINUE
C 850  FORMAT(/,' Gauss Point ng=',I3,/)
C 86   FORMAT(' ggpx:',' i=',I2,3(' j=',I2,E11.3),'  & ',
C     '       ' sspx:',' i=',I2,3(' j=',I2,E11.3))
C 87   FORMAT(' ggpx:',' i=',I2,3(' j=',I2,E11.3))
C 88   FORMAT(' sspx:',' i=',I2,3(' j=',I2,E11.3))
C 860  FORMAT(' smx:',' i=',I2,2(' j=',I2,E11.3),'  & ',
C     '       ' snx:',' i=',I2,2(' j=',I2,E11.3))
C 870  FORMAT(' smx:',' i=',I2,2(' j=',I2,E11.3))
C 880  FORMAT(' snx:',' i=',I2,2(' j=',I2,E11.3))
      ENDIF
      ENDIF
 80   CONTINUE

      CALL EXITS('OPC40')
      RETURN
 9999 CALL ERRORS('OPC40',ERROR)
      CALL EXITS('OPC40')
      RETURN 1
      END


