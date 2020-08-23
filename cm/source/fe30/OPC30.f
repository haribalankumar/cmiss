      SUBROUTINE OPC30(IBT,IDO,INP,NBH,NBJ,NEELEM,NHE,NHP,
     '  NKH,NKHE,NKJE,NPF,NPNE,NPNODE,nr,NRE,NVHE,NVHP,NVJE,
     '  NW,NYNE,NYNP,CE,CG,CGE,CP,CURVCORRECT,PG,RG,SE,WG,XA,XE,XG,
     '  XP,YP,ZA,ZE,ZG,ZP,ERROR,*)

C#### Subroutine: OPC30
C###  Description:
C###    OPC30 prints solns ZP(nk,nv,nh,np,nc) and ZA(na,nh,nc,ne) and
C###    nodal reactions and, if requested, the element Gauss point
C###    solutions ZG(nh,ni) and the solution at any specified point in
C###    theta coordinates. For Boundary Element problems it calculates
C###    a global flux.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NHE(NEM),
     '  NHP(NPM),NKH(NHM,NPM,NCM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),CP(NMM,NPM),
     '  CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),
     '  RG(NGM),SE(NSM,NBFM,NEM),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM),ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,ITMP,na,nb,nc,ne,ng,nh,nhx,ni,nj,nk,
     '  noelem,nonode,NOQUES,np,nu,NU1(0:3),nv,nx
      REAL*8 AB,CG4,CG6,DXIX(3,3),GL(3,3),
     '  GU(3,3),PXI,RAD,RWG,SATE(25),SAT,SATTOT,SUM,SUMTOT,TOL_XCOORD,
     '  VOL,VOLTOT,XB_LOCAL(10),XI(3),XPOINT(3)
      LOGICAL FILEIP
      DATA NU1/1,2,4,7/

      CALL ENTERS('OPC30',*9999)
      nc=1 !Temporary
      nv=1 !Temporary
      nx=1 !Temporary
      TOL_XCOORD=1.0d-06 !MHT 8-11-96 tolerance for calling XCOORD
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      IF(ITYP2(nr,nx).EQ.5.AND.ITYP3(nr,nx).EQ.2) THEN
        !lung gas transport
      ELSE
        WRITE(OP_STRING,'(/'' Nodal solutions:''/)')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
!new MPN 6-Jan-95: current soln now stored in YP(ny,1) for nonlin probs
        CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '    nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
!old
c        IF(ITYP6(nr,nx).EQ.1) THEN
c          CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
c     '      nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
c        ELSE IF(ITYP6(nr,nx).EQ.2) THEN
c          CALL YPZP(4,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
c     '      nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
c        ENDIF
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          WRITE(OP_STRING,
     '      '('' Node'',I3,'' ZP(nk,1,1,'',I3,'',nc): '',5E11.3,'
     '      //'/(23X,5E11.3))') np,np,
     '      (ZP(nk,1,1,np,nc),nk=1,NKH(NH_LOC(1,nx),np,1))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          nv=1 !Temporary AJP 18-12-91
          DO nhx=2,NHP(np)
            nh=nh_loc(nhx,nx)
            WRITE(OP_STRING,
     '        '(9X,''ZP(nk,'',I1,'','',I1,'','',I3,'',nc): '','
     '        //'5E11.3,/(23X,5E11.3))') nh,nc,np,
     '        (ZP(nk,nv,nh,np,nc),nk=1,NKH(nh,np,nc))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO

        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nhx=1,NHE(ne)
            nh=nh_loc(nhx,nx)
            nc=1 !Temporary AJP 18-12-91
            nb=NBH(nh,nc,ne)
            IF(NAT(nb).GT.0) THEN
              WRITE(OP_STRING,
     '          '('' ZA(na,'',I1,'','',I1,'','',I2,''): '',6E11.4,'
     '          //'/(14X,6E11.4))')
     '          nh,nc,ne,(ZA(na,nh,nc,ne),na=1,NAT(nb))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
        ENDDO

        WRITE(OP_STRING,'(/'' Nodal reactions:''/)')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(KTYP90.GT.1)THEN !Coupled FE flux values
          nc=2
        ELSE
          nc=1
        ENDIF
        CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '    nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          WRITE(OP_STRING,
     '      '('' Node'',I3,'' ZP(nk,1,'',I1,'','',I3,'',nc): '',5E11.3,'
     '      //'/(23X,5E11.3))')
     '      np,nv,np,(ZP(nk,nv,1,np,nc),nk=1,NKH(NH_LOC(1,nx),np,nc))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nhx=2,NHP(np)
            nh=nh_loc(nhx,nx)
            WRITE(OP_STRING,
     '        '(9X,''ZP(nk,'',I1,'','',I1,'','',I3,'',nc): '','
     '        //'5E11.3,/(23X,5E11.3))') nh,nv,np,
     '        (ZP(nk,nv,nh,np,nc),nk=1,
     '        NKH(nh,np,nc))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO

!new MPN 6-Jan-95: current soln now stored in YP(ny,1) for nonlin probs
        CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '    nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
!old
c        IF(ITYP6(nr,nx).EQ.1) THEN
c          CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
c     '      nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
c        ELSE IF(ITYP6(nr,nx).EQ.2) THEN
c          CALL YPZP(4,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
c     '      nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
c        ENDIF

        FORMAT='(/$,'' Do you want solution printed out at all'//
     '    ' element Gauss points [N]? '',A)'
        ITMP=0
        CALL GINOUT(ITMP,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,1,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
c       CALL AINOUT(ITMP,IVDU,IFILE,FORMAT,1,ADATA,ANO,INFO,ERROR,*9999)
        IF(ADATA(1).EQ.'Y') THEN
          WRITE(OP_STRING,'(/'' Gauss point solutions:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          VOLTOT=0.0d0
          SUMTOT=0.0d0
          SATTOT=0.0d0
          nc=1 !Temporary AJP 18-12-91
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            WRITE(OP_STRING,'(/'' Element '',I2)') ne
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c GBS IBDRY not used
c            IF(KTYP6 .GT. 0) THEN
c              IBDRY=NW(ne,1)-20
c            ELSE IF(KTYP6 .EQ. 0) THEN
c              IBDRY=NW(ne,1)
c            ENDIF
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '        NRE(ne),NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            CALL ZPZE(NBH(1,1,ne),nc,NHE(ne),NKHE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '        ERROR,*9999)
            CALL CPCG(1,NBJ(1,ne),NPNE(1,1,ne),nr,nx,CE(1,ne),CG,
     '        CGE(1,1,ne),CP,PG,ERROR,*9999)
            DO nhx=1,NHE(ne)
              nh=nh_loc(nhx,nx)
              nb=NBH(nh,nc,ne)
              WRITE(OP_STRING,'('' Dependent variable nh= '',I1) ')
     '          nhx
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              SUM=0.0d0
              VOL=0.0d0
              SAT=0.0d0
              DO ng=1,NGT(nb)
                CALL XEXG(NBJ(1,ne),ng,NRE(ne),PG,
     '            XE,XG,ERROR,*9999)
                CALL XGMG(0,NIT(1),1,NRE(ne),DXIX,GL,GU,RG(ng),
     '            XG,ERROR,*9999)
                RWG=RG(ng)*WG(ng,NBJ(1,ne))
                IF(JTYP4.EQ.2) THEN      !cyl symmetric about x or r
                  IF(ITYP10(1).EQ.1) THEN !rect.cart coords
                    RAD=XG(2,1)
                  ELSE IF(ITYP10(1).EQ.2.OR.ITYP10(1).EQ.3) THEN !polar coords
                    ERROR='Not implemented'
                    GO TO 9999
                  ELSE IF(ITYP10(1).EQ.4) THEN
                    ERROR='Not implemented'
                    GO TO 9999
                  ELSE IF(ITYP10(1).EQ.5) THEN
                    ERROR='Not implemented'
                    GO TO 9999
                  ENDIF
                  RWG=RWG*2.0d0*PI*RAD
                ELSE IF(JTYP4.EQ.3) THEN !cylindrically symmetric about y or z
                  IF(ITYP10(1).EQ.1) THEN !rect.cart coords
                    RAD=XG(1,1)
                  ELSE IF(ITYP10(1).EQ.2.OR.ITYP10(1).EQ.3) THEN !polar coords
                    ERROR='Not implemented'
                    GO TO 9999
                  ELSE IF(ITYP10(1).EQ.4) THEN
                    ERROR='Not implemented'
                    GO TO 9999
                  ELSE IF(ITYP10(1).EQ.5) THEN
                    ERROR='Not implemented'
                    GO TO 9999
                  ENDIF
                  RWG=RWG*2.0d0*PI*RAD
                ELSE IF(JTYP4.EQ.4) THEN !spherically   symmetric
                  RWG=RWG*4.0d0*PI*XG(1,1)**2
                ENDIF
                CALL ZEZG(1,NBH(1,1,ne),ng,NHE(ne),nx,DXIX,PG,ZE,
     '            ZG,ERROR,*9999)
                WRITE(OP_STRING,
     '            '('' Gauss point'',I2,'' ZG('',I1,'',ni): '','
     '            //'4E12.4)') ng,nh,(ZG(nh,NU1(ni)),ni=0,NIT(nb))
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                IF(ITYP5(nr,nx).EQ.1.AND.ITYP2(nr,nx).EQ.8) THEN
                  CG4=DLOG(1.0d0/DABS(1.0d0-0.50d0**(1.0d0/CG(5,ng))))
     '              /CG(4,ng)
                  CG6=DLOG(1.0d0/DABS(1.0d0-0.50d0**(1.0d0/CG(7,ng))))
     '              /CG(6,ng)
                  SUM=SUM+CG(3,ng)*(DABS(1.0d0-DEXP(-CG4*ZG(1,1))))
     '              **CG(5,ng)*RWG
                  VOL=VOL+RWG
                  IF(CG(7,ng).gt.0.0d0) THEN
                    SATE(ng)=(1.0d0-DEXP(-CG6*DABS(ZG(1,1))))**CG(7,ng)
                  ELSE
                    AB=(ZG(1,1)/CG(6,ng))**DABS(CG(7,ng))
                    SATE(ng)=AB/(1.0d0+AB)
                  ENDIF
                  SAT=SAT+SATE(ng)*RWG
                ENDIF
              ENDDO
            ENDDO
            IF(ITYP5(nr,nx).EQ.1.AND.ITYP2(nr,nx).EQ.8) THEN
              IF(DABS(CG(3,1)).gt.1.d-6.AND.DABS(CG(8,1)).lt.1.d-6) THEN
                WRITE(OP_STRING,'('' Element is muscle'')')
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                SATTOT=SATTOT+SAT
              ENDIF
              WRITE(OP_STRING,
     '          '('' M/Hb saturation at Gauss points: '','
     '          //'25F4.1)')(SATE(ng),ng=1,NGT(NBH(NH_LOC(1,nx),1,ne)))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,
     '          '('' Metabolism averaged over element ='',E12.4,'
     '          //''' (W/g)'')') SUM/VOL
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              SUMTOT=SUMTOT+SUM
              VOLTOT=VOLTOT+VOL
              WRITE(OP_STRING,
     '          '('' M/Hb sat.n averaged over element ='',E12.4)')
     '          SAT/VOL
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
          IF(ITYP5(nr,nx).EQ.1.AND.ITYP2(nr,nx).EQ.8) THEN
            WRITE(OP_STRING,
     '        '(//'' Metabolism /unit total volume ='',E12.4,'
     '        //''' (W/g)'')') SUMTOT/VOLTOT
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,
     '        '(//'' Mb satur.n /unit total volume ='',E12.4)')
     '        SATTOT/VOLTOT
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

        NJT=NJ_LOC(NJL_GEOM,0,nr)
 7000   FORMAT=
     '    '(/$,'' Do you want solution at a specified point [N]?'',A)'
        ITMP=0
        CALL GINOUT(ITMP,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,1,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
c       CALL AINOUT(ITMP,IVDU,IFILE,FORMAT,1,ADATA,ANO,INFO,ERROR,*9999)
        IF(ADATA(1).EQ.'Y') THEN
          FORMAT=
     '      '(/$'' The Xj-coords of the point are [0,..]? '',3E12.4)'
          DO nj=1,NJT
            RDEFLT(nj)=0.0d0
          ENDDO
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '      IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,
     '      INFO,ERROR,*9999)
c         CALL RINOUT(0,IVDU,IFILE,FORMAT,NJT,RDATA,RDEFLT,
c    '      -RMAX,RMAX,INFO,ERROR,*9999)
          DO nj=1,NJT
            XPOINT(nj)=RDATA(nj)
          ENDDO
          IF(ITYP10(1).GT.1) XPOINT(2)=XPOINT(2)*PI/180.0d0
          IF(ITYP10(1).GT.2) XPOINT(3)=XPOINT(3)*PI/180.0d0
          WRITE(OP_STRING,
     '      '(/'' Solution at Xj-coordinates: '',3E12.4)')
     '      (XPOINT(nj),nj=1,NJT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          CALL XCOORD(IBT,IDO,INP,NBJ,ne,NEELEM,
     '      NKJE,NPF,NPNE,NPNODE,nr,NVJE,
     '      SE,TOL_XCOORD,XA,XE,XI,XP,XPOINT,ERROR,*9999)
          nb=NBJ(1,ne)
          WRITE(OP_STRING,'('' The point is in element'','
     '      //'I5,'' at Xi-coords:'',3E12.3)')
     '      ne,(XI(ni),ni=1,NIT(nb))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          CALL ZPZE(NBH(1,1,ne),nc,NHE(ne),NKHE(1,1,1,ne),
     '      NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),
     '      nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '      ZE,ZP,ERROR,*9999)
          nc=1 !Temporary AJP 18-12-91
          DO nhx=1,NHE(ne)
            nh=nh_loc(nhx,nx)
            nb=NBH(nh,nc,ne)
            DO nu=1,NUT(nb)
              XB_LOCAL(nu)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '          INP(1,1,nb),nb,nu,XI,ZE(1,nhx))
            ENDDO
            WRITE(OP_STRING,
     '        '('' Xi-derivatives (up to 2nd order) of dependent'
     '        //' variable '' ,I1,'': '',(1X,6E12.4))') nhx,
     '        (XB_LOCAL(nu),nu=1,NUT(nb))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
          GO TO 7000
        ENDIF
      ENDIF

      CALL EXITS('OPC30')
      RETURN
 9999 CALL ERRORS('OPC30',ERROR)
      CALL EXITS('OPC30')
      RETURN 1
      END



