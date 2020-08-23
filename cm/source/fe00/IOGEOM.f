      SUBROUTINE IOGEOM(COMAND,IUNIT,IBT,IDO,INP,ITHRES,NAN,NBH,
     '  NBJ,NBJF,NDET,NEELEM,NEL,NFF,NGAP,NHE,NHP,
     '  NKJE,NKEF,NKH,NKJ,NLF,NLL,NNF,NNL,NONY,NPF,NPL,NPNE,
     '  NPNODE,NPNY,NP_INTERFACE,NRE,NVHE,NVHP,NVJE,NVJP,NW,
     '  NXI,NYNE,NYNO,NYNP,
     '  CE,CONY,CP,DET,DL,FEXT,PE,PF,PG,SE,THRES,
     '  WG,XA,XIG,XP,YP,FIX,FIXP,ERROR,*)

C#### Subroutine: IOGEOM
C###  Description:
C###    IOGEOM reads and writes IOD files.

C**** COMAND determines whether data is to be read from ('READ') or
C**** written to ('WRITE') IUNIT.

      IMPLICIT NONE
      INCLUDE 'b10.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'b40.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'four00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp40.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'ktyp60.cmn'
      INCLUDE 'ktyp70.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxbc00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ITHRES(3,NGM,NEM),NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NDET(NBFM,0:NNM),NEELEM(0:NE_R_M,0:NRM),
     '  NEL(0:NELM,NLM),NFF(6,NEM),NGAP(NIM,NBM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKJE(NKM,NNM,NJM,NEM),
     '  NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),NKJ(NJM,NPM),
     '  NLF(4,NFM),NLL(12,NEM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  NPF(9,NFM),NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM,NXM),
     '  NP_INTERFACE(0:NPM,0:3),NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJP(NJM,NPM),NW(NEM,3,NXM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM,NXM),CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  CP(NMM,NPM,NXM),DET(NBFM,0:NNM,NGM,6),DL(3,NLM),
     '  FEXT(NIFEXTM,NGM,NEM),
     '  PE(2,NEM),PF(2,NEM),PG(NSM,NUM,NGM,NBM),
     '  SE(NSM,NBFM,NEM),THRES(3,NGM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XIG(NIM,NGM,NBM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM)
      CHARACTER COMAND*(*),ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM,NXM),FIXP(2,*)
!     Local Variables
      INTEGER i,id,ie,iee,iface,il,IUNIT,iy,
     '  j,mg,na,nb,nbb,nb2,NBTOP,nc,ncc,ne,nEE,nf,nfee,
     '  ng,nh,NHH,nhx,ni,nj,njj,njj1,njj2,nk,nkk,nl,nm,nn,nnn,
     '  noelem,nonode,NOY,np,NPP,nr,nrc,nrcc,
     '  NRR,ns,nv,nvv,nx,nx1,nxx,ny,NYY
      CHARACTER FMT*500
      LOGICAL DEPENDENT,FIRST

      CALL ENTERS('IOGEOM',*9999)

      nc=1 !Temporary MPN 23/11/94
C      nx=1 !Temporary

!********************** read iod file ********************************

      IF(COMAND.EQ.'READ') THEN
!nx_list
        FMT='('' NX_LIST(0..): '',10I3)'
        READ(IUNIT,FMT) NX_LIST(0),(NX_LIST(nxx),nxx=1,NX_LIST(0))
!n*t params
        FMT='('' NBT='',I6,'' NDT='',I6,'' NFT='',I6,'
     '     //''' NJT='',I6,'' NLT='',I6,'' NQT='',I6,/'
     '     //''' NRT='',I6,'' NBFT='',I5)'
        READ(IUNIT,FMT) NBT,NDT,NFT,NJT,NLT,NQT,NRT,NBFT

        DO nxx=1,NX_LIST(0)
          nx=NX_LIST(nxx)
          FMT='(/'' NZT(1,nx): '',I6,I6)'
          READ(IUNIT,FMT) nx1,NZT(1,nx)
        ENDDO !nxx
!jtypes
        FMT='(/'' JTYP1 ='',I2,'' JTYP2A ='',I2,'' JTYP2B ='',I2,/'
     '      //''' JTYP2C ='',I2,'' JTYP3 ='',I2,/'
     '      //''' JTYP4 ='',I2,'' JTYP5 ='',I2,'' JTYP6 ='',I2,/'
     '      //''' JTYP7 ='',I2,'' JTYP8 ='',I2,'' JTYP9 ='',I2,'
     '      //''' JTYP10='',I2,'' JTYP11='',I2,'' JTYP12='',I2,/'
     '      //''' JTYP13='',I2,'' JTYP14='',I2,'' JTYP15='',I2)'
        READ(IUNIT,FMT) JTYP1,JTYP2A,JTYP2B,JTYP2C,JTYP3,JTYP4,
     '                  JTYP5,JTYP6,
     '                  JTYP7,JTYP8,JTYP9,JTYP10,JTYP11,JTYP12,
     '                  JTYP13,JTYP14,JTYP15
!net,npt
        FMT='(/'' NET(0)='',I6,'' NPT(0)='',I6)'
        READ(IUNIT,FMT) NET(0),NPT(0)
!region params
        NJ_LOC(0,0,0)=0
        DO nr=1,NRT
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            FMT='(/'' nx= '',I6)'
            READ(IUNIT,FMT) nx1
            FMT='(/'' Region '',I1,'': ITYP1 ='',I2,'' ITYP2 ='',I2,'
     '                         //''' ITYP3 ='',I2,'' ITYP4 ='',I2,/10X,'
     '                         //''' ITYP5 ='',I2,'' ITYP6 ='',I2,'
     '                         //''' ITYP7 ='',I2,'' ITYP8 ='',I2,/10X,'
     '                         //''' ITYP9 ='',I2,'' ITYP12 ='',I2)'
            READ(IUNIT,FMT) nrr,ITYP1(nr,nx),ITYP2(nr,nx),ITYP3(nr,nx),
     '        ITYP4(nr,nx),ITYP5(nr,nx),ITYP6(nr,nx),ITYP7(nr,nx),
     '        ITYP8(nr,nx),ITYP9(nr,nx),ITYP12(nr,nx)
          ENDDO !nx
          FMT='(/'' Region '',I1,'': ITYP10='',I2,'
     '                       //''' ITYP11='',I2,/10X,'
     '                       //''' ITYP13='',I2,'' ITYP14='',I2,'
     '                       //''' ITYP15='',I2)'
C KAT 31Oct97: ITYP15 and ITYP16 should be in the nx loop but they are
C not set up properly yet to avoid corrupting present IOD files.
          READ(IUNIT,FMT) nrr,ITYP10(nr),ITYP11(nr),ITYP13(nr),
     &      ITYP14(nr),ITYP15(nr,1)
          FMT='(10X,'' NET('',I2,'')='',I6,'' NPT('',I2,'')='',I6)'
          READ(IUNIT,FMT) NRR,NET(nr),NRR,NPT(nr)
!nj_loc
          FMT='(/'' NJ_LOC(1..3,0,nr): '',3I3)'
          READ(IUNIT,FMT) (NJ_LOC(njj1,0,nr),njj1=1,3)
          FMT='( '' NJ_LOC(1,1..,nr):  '',12I3)'
          READ(IUNIT,FMT) (NJ_LOC(1,njj2,nr),njj2=1,NJ_LOC(1,0,nr))
          FMT='( '' NJ_LOC(2,1..,nr):  '',12I3)'
          READ(IUNIT,FMT) (NJ_LOC(2,njj2,nr),njj2=1,NJ_LOC(2,0,nr))
          FMT='( '' NJ_LOC(3,1..,nr):  '',12I3)'
          READ(IUNIT,FMT) (NJ_LOC(3,njj2,nr),njj2=1,NJ_LOC(3,0,nr))
C ***     Calculate NJ_LOC(0,0) and NJ_TYPE
          NJ_LOC(0,0,nr)=0
          DO njj1=1,3
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              NJ_TYPE(nj,1)=njj1
              NJ_TYPE(nj,2)=njj2
              IF(nj.GT.NJ_LOC(0,0,nr)) NJ_LOC(0,0,nr)=nj
            ENDDO
          ENDDO
          IF(NJ_LOC(0,0,nr).GT.NJ_LOC(0,0,0))
     '      NJ_LOC(0,0,0)=NJ_LOC(0,0,nr)
        ENDDO

!basis function info
        DO nb=1,NBFT
          FMT='(/'' Basis type nb='',I2)'
          READ(IUNIT,FMT) id
          FMT='(/'' NAT(nb)='',I3,''  NBI(nb)='',I3,''  NFE(nb)='',I3,'
     '      //'''  NGT(nb)='',I3,''  NIT(nb)='',I3,/'
     '       //''' NLE(nb)='',I3,''  NNT(nb)='',I3,''  NST(nb)='',I3,'
     '      //'''  NUT(nb)='',I3,''  NBC(nb)='',I3,/'
     '       //''' NABTYP(nb)='',I3,'' NMGT='',I3)'
          READ(IUNIT,FMT) NAT(nb),NBI(nb),NFE(nb),NGT(nb),NIT(nb),
     '      NLE(nb),NNT(nb),NST(nb),NUT(nb),NBC(nb),NABTYP(nb),NMGT
          FMT='('' NKT(nn=0..,nb):'',18I3)'
          READ(IUNIT,FMT) (NKT(nn,nb),nn=0,NNT(nb))
          FMT='('' IBT(i=1..,ni=1..,nb):'',40I2)'
          READ(IUNIT,FMT) ((IBT(i,ni,nb),i=1,3),ni=1,NIT(nb))
          DO nn=1,NNT(nb)
            FMT='('' IDO(nk=1..,nn='',I2,'',ni=1..,nb):'',40I2)'
            READ(IUNIT,FMT) nnn,((IDO(nk,nn,ni,nb),nk=1,NKT(nn,nb)),
     '        ni=0,NIT(nb))
          ENDDO
          FMT='('' INP(nn=1..,ni=1..,nb):'',40I2)'
          READ(IUNIT,FMT) ((INP(nn,ni,nb),nn=1,NNT(nb)),ni=1,NIT(nb))
          FMT='('' NAN(ni=1..,na=1..,nb):'',40I2)'
          READ(IUNIT,FMT) ((NAN(ni,na,nb),ni=1,NIT(nb)),na=1,NAT(nb))
          FMT='('' NGAP(ni=1..,nb):'',20I4)'
          READ(IUNIT,FMT) (NGAP(ni,nb),ni=1,NIT(nb))
          FMT='('' NNL(i=0..,j=1..,nb):'',/(12I4))'
          READ(IUNIT,FMT) ((NNL(i,j,nb),i=0,4),j=1,12)
          DO nf=1,NFE(nb)
            FMT='('' NNF(0,nf='',I1,'',nb):'',I3)'
            READ(IUNIT,FMT) nfee,NNF(0,nf,nb)
            FMT='('' NNF(1..,nf='',I1,'',nb):'',17I3)'
            READ(IUNIT,FMT) nfee,(NNF(i,nf,nb),i=1,NNF(0,nf,nb))
            FMT='('' NKEF(0,nn=1..,nf='',I1,'',nb):'',16I2)'
            READ(IUNIT,FMT) nfee,(NKEF(0,j,nf,nb),j=1,NNF(0,nf,nb))
            FMT='('' NKEF(1..,nn=1..,nf='',I1,'',nb):'',16I2)'
            READ(IUNIT,FMT) nfee,((NKEF(i,j,nf,nb),i=1,NKEF(0,j,nf,nb)),
     '        j=1,NNF(0,nf,nb))
          ENDDO !nfe
          IF(NBC(nb).EQ.4) THEN
            FMT='('' OMEGA='',E25.17)'
            READ(IUNIT,FMT) OMEGA
          ENDIF
          IF((NBC(nb).EQ.5).OR.(NBC(nb).EQ.6)) THEN !BE basis function
            FMT='('' NGLIMITS(ni=1..,nb,1..2):'',/(20I4))'
            READ(IUNIT,FMT) ((NGLIMITS(ni,nb,i),ni=1,NIT(nb)),i=1,2)
          ENDIF
        ENDDO
!focus
        FMT='(/'' Focus:'',E25.17)'
        READ(IUNIT,FMT) FOCUS

!nodes
C LKC 28-JUN-2001
C        READ(IUNIT,'(/30X,''*** Nodes ***'')')
        FMT='(/30X,''*** Nodes ***'')'
        READ(IUNIT,FMT)

        FMT='('' NPNODE(0,0):    '',I5)'
        READ(IUNIT,FMT) NPNODE(0,0)
        DO nr=1,NRT
          FMT='(/'' NPNODE(0,nr='',I2,''):'',I5)'
          READ(IUNIT,FMT) nrr,NPNODE(0,nr)
          FMT='('' NPNODE(nonode=1..,nr='',I2,''):''/,(15I5))'
          READ(IUNIT,FMT) nrr,(NPNODE(nonode,nr),nonode=1,NPNODE(0,nr))
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            FMT='(/'' Node np='',I4)'
            READ(IUNIT,FMT) NPP
            FMT='('' NVJP(nj=1..,np):'',12I3)'
            READ(IUNIT,FMT) ((NVJP(NJ_LOC(njj1,njj2,nr),np),
     '        njj2=1,NJ_LOC(njj1,0,nr)),njj1=1,3)
            FMT='('' NKJ(nj=1..,np):'',20I3)'
            READ(IUNIT,FMT) ((NKJ(NJ_LOC(njj1,njj2,nr),np),
     '        njj2=1,NJ_LOC(njj1,0,nr)),njj1=1,3)
            DO njj1=1,3  !geometry/fibres/field
              DO njj2=1,NJ_LOC(njj1,0,nr)
                nj=NJ_LOC(njj1,njj2,nr)
                DO nv=1,NVJP(nj,np)
                  FMT='('' XP(nk=1..,nv='',I2,'',nj='',I1,'',np):'','
     '              //'3E25.17,:/(3E25.17))'
                  READ(IUNIT,FMT)
     '              nvv,njj,(XP(nk,nv,nj,np),nk=1,NKJ(nj,np))
                ENDDO !nv
              ENDDO !njj2
            ENDDO !njj1
            FMT='('' NP_INTERFACE(np,0): '',I2)'
            READ(IUNIT,FMT) NP_INTERFACE(np,0)
            FMT='('' NP_INTERFACE(np,i=1..): '',3I5)'
            READ(IUNIT,FMT)
     '        (NP_INTERFACE(np,i),i=1,NP_INTERFACE(np,0))
          ENDDO !nonode
        ENDDO !nr
!elements
C LKC 28-JUN-2001
C        READ(IUNIT,'(/30X,''*** Elements ***'')')
        FMT='(/30X,''*** Elements ***'')'
        READ(IUNIT,FMT)

        FMT='('' NEELEM(0,0):    '',I5)'
        READ(IUNIT,FMT) NEELEM(0,0)
        DO nr=1,NRT
          FMT='(/'' NEELEM(0,nr='',I2,''):'',I5)'
          READ(IUNIT,FMT) nrr,NEELEM(0,nr)
          FMT='('' ITYP10(nr)='',I2)'
          READ(IUNIT,FMT) ITYP10(nr)
          FMT='('' NEELEM(noelem=1..,nr='',I2,''):''/,(15I5))'
          READ(IUNIT,FMT) nrr,(NEELEM(noelem,nr),noelem=1,NEELEM(0,nr))
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            FMT='(/'' Element ne='',I4)'
            READ(IUNIT,FMT) NEE
            FMT='('' NRE(ne)='',I2)'
            READ(IUNIT,FMT) NRE(ne)
            FMT='('' NBJ(nj=1..,ne):'',20I3)'
            READ(IUNIT,FMT) ((NBJ(NJ_LOC(njj1,njj2,nr),ne),
     '        njj2=1,NJ_LOC(njj1,0,nr)),njj1=1,3)
            FMT='('' NLL(i=1..,ne):'',20I4)'
            READ(IUNIT,FMT) (NLL(i,ne),i=1,12)
            FMT='('' NFF(i=1..,ne):'',20I4)'
            READ(IUNIT,FMT) (NFF(i,ne),i=1,6)
            DO nb=1,NBFT
              FMT='(/'' NPNE(nn=1..,nb='',I2,'',ne):'',20I4)'
              READ(IUNIT,FMT) NBB,(NPNE(nn,nb,ne),nn=1,NNT(nb))
C              DO nn=1,NNT(nb)
C                FMT='('' NKE(nk=1..,nn='',I2,'',nb='',I2,'',ne):'','
C     '            //'(20I4))'
C                READ(IUNIT,FMT) nnn,nbb,
C     '            (NKE(nk,nn,nb,ne),nk=1,NKT(nn,nb))
C              ENDDO
              FMT='('' SE(ns=1..,nb='',I2,'',ne):'',/(3E25.17))'
              READ(IUNIT,FMT) nbb,(SE(ns,nb,ne),ns=1,NST(nb))
              DO njj1=1,3  !geometry/fibres/field
                DO njj2=1,NJ_LOC(njj1,0,nr)
                  nj=NJ_LOC(njj1,njj2,nr)
                  FMT='('' NVJE(nn=1..,nb='',I2,'',nj='',I2,'
     '              //''',ne):'',20I2)'
                  READ(IUNIT,FMT) NBB,NJJ,
     '              (NVJE(nn,nb,nj,ne),nn=1,NNT(nb))
                ENDDO !njj1
              ENDDO !njj2
            ENDDO !nb
!xa
            FIRST=.TRUE.
            DO njj1=1,3  !geometry/fibres/field
              DO njj2=1,NJ_LOC(njj1,0,nr)
                nj=NJ_LOC(njj1,njj2,nr)
                nb=NBJ(nj,ne)
                DO nn=1,NNT(nb)
                  FMT='('' NKJE(nk=1..,nn='',I2,'',nj='',I2,'',ne):'','
     '              //'(20I4))'
                  READ(IUNIT,FMT) nnn,njj,
     '              (NKJE(nk,nn,nj,ne),nk=1,NKT(nn,nb))
                ENDDO
                IF(NAT(nb).NE.0) THEN
                  IF(FIRST) THEN
                    FMT='()'
                    READ(IUNIT,FMT)
                    FIRST=.FALSE.
                  ENDIF !first
                  FMT='('' XA(na=1..,nj='',I2,'',ne):'',2E25.17'//
     '              ':/(3E25.17))'
                  READ(IUNIT,FMT) NJJ,(XA(na,nj,ne),na=1,NAT(nb))
                ENDIF !NAT.NE.0
              ENDDO !njj1
            ENDDO !njj2

          ENDDO !noelem
        ENDDO !nr
!nxi
C LKC 28-JUN-2001
C        READ(IUNIT,'(/30X,''***  NXI  ***'')')
        FMT='(/30X,''***  NXI  ***'')'
        READ(IUNIT,FMT)

C new MPN 22Apr97: wrong ne index for NXI
        DO nr=1,NRT
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO ni=-NIM,NIM
              FMT='('' NXI(ni='',I4,'',0-4,ne='',I4,''):'',20I4)'
              READ(IUNIT,FMT) id,id,(NXI(ni,i,ne),i=0,4)
            ENDDO !ni
          ENDDO !noelem
        ENDDO !nr
C old
C        DO ne=1,NEELEM(0,0)
C          DO ni=-NIM,NIM
C            FMT='('' NXI(ni='',I4,'',0-4,ne='',I4,''):'',20I4)'
C            READ(IUNIT,FMT) id,id,(NXI(ni,i,ne),i=0,4)
C          ENDDO !ni
C        ENDDO !ne
C!xa
C        FMT='(/'' XA(1,nj=1..,ne=1..):'',/(3E25.17))'
C        READ(IUNIT,FMT) ((XA(1,nj,ne),nj=1,NJM),ne=1,NET(0))


!faces
C LKC 28-JUN-2001
C        READ(IUNIT,'(/30X,''*** Faces ***'')')
        FMT='(/30X,''*** Faces ***'')'
        READ(IUNIT,FMT)
        DO nf=1,NFT
          FMT='(/'' Face nf='',I4)'
          READ(IUNIT,FMT) id
          FMT='('' NPF(i=1..,nf):'',9I6)'
          READ(IUNIT,FMT) (NPF(i,nf),i=1,9)
          FMT='('' NLF(i=1..,nf):'',4I4)'
          READ(IUNIT,FMT) (NLF(i,nf),i=1,4)
          FMT='('' NBJF(nj=1..,nf):'',(3I4))'
          READ(IUNIT,FMT) (NBJF(nj,nf),nj=1,NJT)
        ENDDO !nf

!lines
C LKC 28-JUN-2001
C        READ(IUNIT,'(/30X,''*** Lines ***'')')
        FMT='(/30X,''*** Lines ***'')'
        READ(IUNIT,FMT)
        DO nl=1,NLT
          FMT='(/'' Line nl='',I4)'
          READ(IUNIT,FMT) id
          DO nj=0,3
            FMT='('' NPL(i=1..,nj='',I1,'',nl):'',5(1X,I4))'
            READ(IUNIT,FMT) njj, (NPL(i,nj,nl),i=1,5)
          ENDDO
          FMT='('' NEL(i=0..,nl):'',10(1X,I4),/,15X,10(1X,I4))'
          READ(IUNIT,FMT) NEL(0,nl),(NEL(i,nl),i=1,NEL(0,nl))
          FMT='('' DL(i=1..,nl):'',3E25.17)'
          READ(IUNIT,FMT) (DL(i,nl),i=1,3)
        ENDDO !nl
!call logicals
        FMT='(/'' CALL_BASE='',L1,'' CALL_ELEM='',L1,'
     '    //'  '' CALL_LINE='',L1,'' CALL_NODE='',L1,'
     '    //'  '' CALL_FIBR='',L1,'' CALL_ELFB='',L1)'
        READ(IUNIT,FMT) CALL_BASE,CALL_ELEM,CALL_LINE,CALL_NODE,
     '    CALL_FIBR,CALL_ELFB

            !***********  dependent variables *************

C        nx=1 !Temporary
        DEPENDENT=.FALSE.
        DO nxx=1,NX_LIST(0)
          nx=NX_LIST(nxx)
          DO nr=1,NRT
            IF(ITYP1(nr,nx).GT.1) DEPENDENT=.TRUE.
          ENDDO
        ENDDO

        IF(DEPENDENT) THEN
!ktypes
          FMT='(/'' KTYP1 ='',I2,'' KTYP2 ='',I2,'' KTYP3 ='',I2,'
     '      //'  '' KTYP4 ='',I2,'' KTYP5 ='',I2,'' KTYP6 ='',I2,/'
     '      //'  '' KTYP7 ='',I2,'' KTYP8 ='',I2,'' KTYP9 ='',I2,'
     '      //'  '' KTYP10='',I2,'' KTYP11='',I2,'' KTYP12='',I2,/'
     '      //'  '' KTYP13='',I2,'' KTYP14='',I2,'' KTYP15='',I2,'
     '      //'  '' KTYP16='',I2,'' KTYP17='',I2)'
          READ(IUNIT,FMT) KTYP1,KTYP2,KTYP3,KTYP4,KTYP5,KTYP6,
     '      KTYP7,KTYP8,KTYP9,KTYP10,KTYP11,KTYP12,
     '      KTYP13,KTYP14,KTYP15,KTYP16,KTYP17
          FMT='( '' KTYP19='',I2,'' KTYP1A='',I2,'' KTYP1B='',I2,'
     '      //'  '' KTYP1C='',I2,'' KTYP1D='',I2,'' KTYP1E='',I2,/'
     '      //'  '' KTYP20='',I2,'' KTYP21='',I2,'' KTYP22='',I2,'
     '      //'  '' KTYP23='',I2,'' KTYP24='',I2,'' KTYP25='',I2,/'
     '      //'  '' KTYP26='',I2,'' KTYP27='',I2,'' KTYP28='',I2,'
     '      //'  '' KTYP29='',I2)'
          READ(IUNIT,FMT) KTYP19,KTYP1A,KTYP1B,KTYP1C,KTYP1D,KTYP1E,
     '      KTYP20,KTYP21,KTYP22,KTYP23,KTYP24,KTYP25,
     '      KTYP26,KTYP27,KTYP28,KTYP29
          FMT='(/'' KTYP30='',I2,'' KTYP31='',I2,'' KTYP32='',I2,'
     '      //'  '' KTYP33='',I2,'' KTYP34='',I2,'' KTYP35='',I2)'
          READ(IUNIT,FMT) KTYP30,KTYP31,KTYP32,KTYP33,KTYP34,
     '      KTYP35
          FMT='(/'' KTYP40='',I2,'' KTYP41='',I2,'' KTYP42='',I2,'
     '      //'  '' KTYP43='',I2,'' KTYP44='',I2,'' KTYP45='',I2)'
          READ(IUNIT,FMT) KTYP40,KTYP41,KTYP42,KTYP43,KTYP44,
     '      KTYP45
          DO nr=1,NRT
            FMT='(/'' Region '',I1)'
            READ(IUNIT,FMT) nrr
            FMT='(/'' KTYP50='',I2,'' KTYP51='',I2,'' KTYP52='',I2,'
     '        //'  '' KTYP53='',I2,'' KTYP54='',I2,'' KTYP55='',I2,'
     '        //' /'' KTYP56='',I2,'' KTYP57='',I2,'' KTYP58='',I2,'
     '        //'  '' KTYP59='',I2,'' KTYP5A='',I2,'' KTYP5B='',I2,'
     '        //' /'' KTYP5C='',I2,'' KTYP5D='',I2,'' KTYP5E='',I2)'
            READ(IUNIT,FMT) KTYP50(nr),KTYP51(nr),KTYP52(nr),
     '        KTYP53(nr),KTYP54(nr),KTYP55(nr),KTYP56(nr),
     '        KTYP57(nr),KTYP58(nr),KTYP59(nr),KTYP5A(nr),
     '        KTYP5B(nr),KTYP5C(nr),KTYP5D(nr),KTYP5E(nr)
          ENDDO
          FMT='('' IL_thickness ='',I2,'' IL_sarcomere='',I2,'
     '      //'/'' IL_time_delay='',I2,'' IL_fluid_conductivity='',I2)'
          READ(IUNIT,FMT) IL_thickness,IL_sarcomere,
     '      IL_time_delay,IL_fluid_conductivity
          FMT='(/'' KTYP60='',I2,'' KTYP61='',I2,'' KTYP62='',I2,'
     '      //'  '' KTYP63='',I2,'' KTYP64='',I2,'' KTYP65='',I2)'
          READ(IUNIT,FMT) KTYP60,KTYP61,KTYP62,KTYP63,KTYP64,
     '      KTYP65
          FMT='(/'' KTYP70='',I2,'' KTYP71='',I2,'' KTYP72='',I2,'
     '      //'  '' KTYP73='',I2,'' KTYP74='',I2,'' KTYP75='',I2)'
          READ(IUNIT,FMT) KTYP70,KTYP71,KTYP72,KTYP73,KTYP74,
     '      KTYP75
          FMT='(/'' KTYP90='',I2,'' KTYP91='',I2,'' KTYP92='',I2)'
          READ(IUNIT,FMT) KTYP90,KTYP91,KTYP92
          FMT='( '' KTYP93(1,nr)='',9I2,'' KTYP93(2,nr)='',9I2)'
          READ(IUNIT,FMT) (KTYP93(1,nr),nr=1,NRT),(KTYP93(2,nr),
     '      nr=1,NRT)
          FMT='( '' KTYP94='',I2,'' KTYPMBC='',I2)'
          READ(IUNIT,FMT) KTYP94,KTYPMBC
!iwrit
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            FMT='(/'' nx= '',I6)'
            READ(IUNIT,FMT) nx1
            FMT='(/'' IWRIT1(nr=1..,nx):'',9I2)'
            READ(IUNIT,FMT) (IWRIT1(nr,nx),nr=1,NRT)
            FMT='( '' IWRIT2(nr=1..,nx):'',9I2)'
            READ(IUNIT,FMT) (IWRIT2(nr,nx),nr=1,NRT)
            FMT='( '' IWRIT3(nr=1..,nx):'',9I2)'
            READ(IUNIT,FMT) (IWRIT3(nr,nx),nr=1,NRT)
            FMT='( '' IWRIT4(nr=1..,nx):'',9I2)'
            READ(IUNIT,FMT) (IWRIT4(nr,nx),nr=1,NRT)
          ENDDO
!time params
          FMT='(/'' DT='',E25.17,'' TINCR='',E25.17,'' T0='',E25.17,'
     '      //''' T1='',E25.17,/'' THETA(1)='',E25.17,'' THETA(2)='','
     '      //'E25.17,'//''' THETA(3)='',E25.17)'
          READ(IUNIT,FMT) DT,TINCR,T0,T1,THETA(1),THETA(2),THETA(3)
!logicals
          FMT='(/'' INIT='',L1,'' LUMP='',L1,'' PROMPT='',L1,'
     '      //''' LRESID='',L1,'' LSTEP='',L1,'' RESTAR='',L1)'
          READ(IUNIT,FMT) INIT,LUMP,PROMPT,LRESID,LSTEP,RESTAR

          FMT='('' NX_CLASS: '',9I3)'
          READ(IUNIT,FMT) (NX_CLASS(NX_LIST(nx)),nx=1,NX_LIST(0))
          FMT='('' NX_LOCKS: '',9I3)'
          READ(IUNIT,FMT) (NX_LOCKS(NX_LIST(nx)),nx=1,NX_LIST(0))
          FMT='('' NX_TYPE: '',9I3)'
          READ(IUNIT,FMT) (NX_TYPE(NX_LIST(nx)),nx=1,NX_LIST(0))
          nrc=2 !temporary AJP
          DO nx1=1,NX_LIST(0)
            nx=NX_LIST(nx1)
            FMT='(/'' NYT(nrc,1,'',I1,'')='',I6,'' NOT(nrc,0,'',I1,'
     '        //''')='',I6)'
            READ(IUNIT,FMT) NXX,NYT(nrc,1,nx),NXX,NOT(nrc,1,0,nx)
          ENDDO !nx1
!call dimchk
          DO nx1=1,NX_LIST(0)
            nx=NX_LIST(nx1)
            CALL DIMCHK(IOOP,nx,ERROR,*9999)
          ENDDO !nx

          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              FMT='(/'' Region nr='',I2)'
              READ(IUNIT,FMT) nrr
              FMT='('' NCT(nr='',I2,'',nx='',I2,'')='',I6)'
              READ(IUNIT,FMT) nrr,nx1,NCT(nr,nx)
              DO nc=1,NCT(nr,nx)
                FMT='('' NYT(nrc='',I1,''nc='',I2,'',nx=1..)='',9I6)'
                READ(IUNIT,FMT) NRCC,NCC,NYT(nrc,nc,nx)
              ENDDO
              FMT='('' NOT(nrc,nr='',I2,'',nx=1..)='',9I6)'
              READ(IUNIT,FMT) nrr,NOT(nrc,1,nr,nx)

              IF(ITYP1(nr,nx).EQ.4) THEN !linear elasticity
                DO ie=1,12
                  FMT='(/'' Element type ie='',I3)'
                  READ(IUNIT,FMT) iee  !iee is never used
                  FMT='('' ETYP(ie)='',L1,'' ILT(ie,nr,nx)='',I3,'
     '              //''' IMT(ie)='','
     '              //'I3,'' NGP(ie)='',I3,'' NLP(ie)='','
     '              //'I3,'' NMP(ie)='',I3,'' NVE(ie)='','
     '              //'I3,'' NHV(nv=1..6,ie):'',6I2)'
                  READ(IUNIT,FMT) ETYP(ie),ILT(ie,nr,nx),IMT(ie),
     '              NGP(ie),NLP(ie),NMP(ie),NVE(ie),(NHV(nv,ie),nv=1,6)
                  FMT='('' ILP(il=1..,ie,nr,nx):'',20I2)'
                  READ(IUNIT,FMT) (ILP(il,ie,nr,nx),il=1,ILT(ie,nr,nx))
                  FMT='('' NMB(il=1..,ie,nx):'',20I2)'
                  READ(IUNIT,FMT) (NMB(il,ie,nx),il=1,ILT(ie,nr,nx))
                  FMT='('' IT(1..3,ie):'',3I3)'
                  READ(IUNIT,FMT) IT(1,ie),IT(2,ie),IT(3,ie)
                ENDDO

              ELSE IF(ITYP1(nr,nx).EQ.3.  !pdes
     '          OR.ITYP1(nr,nx).EQ.5.  !finite elasticity
     '          OR.ITYP1(nr,nx).EQ.9) THEN !BEM
                FMT='('' ETYP(nr)='',L1,'' ILT(1,nr,nx)='','
     '            //'I3,'' IMT(nr)='','
     '           //'I3,'' NGP(nr)='',I3,'' NLP(nr)='',I3,'' NMP(nr)='','
     '            //'I3,'' NVE(nr)='',I3,'' NHV(nv=1..6,nr):'',6I2)'
                READ(IUNIT,FMT) ETYP(nr),
     '            ILT(1,nr,nx),IMT(nr),NGP(nr),NLP(nr),NMP(nr),
     '            NVE(nr),(NHV(nv,nr),nv=1,6)
                FMT='('' ILP(il=1..,1,nr,nx):'',20I2)'
                READ(IUNIT,FMT) (ILP(il,1,nr,nx),il=1,ILT(1,nr,nx))
                FMT='('' NMB(il=1..,nr,nx):'',20I2)'
                READ(IUNIT,FMT) (NMB(il,nr,nx),il=1,ILT(1,nr,nx))
                FMT='('' IT(1..3,nr):'',3I3)'
                READ(IUNIT,FMT) IT(1,nr),IT(2,nr),IT(3,nr)
              ENDIF
            ENDDO !nx
          ENDDO !nr

          nc=1 ! Temporary MPN 23-Nov-94
          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                FMT='(/'' Region nr='',I2,''   Node np='',I4)'
                READ(IUNIT,FMT) NRR,NPP
                FMT='('' NHP(np,nr,nx)='',I2)'
                READ(IUNIT,FMT) NHP(np,nr,nx)
                FMT='('' NKH(nh=1..,np,nc,nr):'',12I3)'
                READ(IUNIT,FMT) (NKH(NH_LOC(nhx,nx),np,nc,nr),
     '            nhx=1,NHP(np,nr,nx))
                FMT='('' NVHP(nh=1..,np,nc,nr):'',12I3)'
                READ(IUNIT,FMT) (NVHP(NH_LOC(nhx,nx),np,nc,nr),
     '            nhx=1,NHP(np,nr,nx))
                DO nhx=1,NHP(np,nr,nx)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,nc,nr)
                    DO nk=1,NKH(nh,np,nc,nr)
                      DO nrc=0,2
                        FMT='('' NYNP(nk='',I2,'',nv='',I2,'',nh='',I2,'
     '                  //''',np,nrc='',I1,'',nc='',I1,'',nr='',I1,''):'
     '                  //''',10I6,/(37X,10I6))'
                        READ(IUNIT,FMT) nkk,nvv,nhh,nrcc,ncc,nrr,
     '                    NYNP(nk,nv,nh,np,nrc,nc,nr)
                        ny=NYNP(nk,nv,nh,np,nrc,nc,nr)
                        FMT='('' NPNY(0..6,ny='',I5,'',nrc='',I1,'
     '                    //''',nx='',I1,'') = '',7I6)'
                        READ(IUNIT,FMT) NYY,NRCC,NX1,
     '                    (NPNY(i,ny,nrc,nx),i=0,6)
                      ENDDO !nrc
                    ENDDO
                  ENDDO
                ENDDO
                FMT='('' CP(nm=1..,np,nx):'',/(3E25.17))'
                READ(IUNIT,FMT) (CP(nm,np,nx),NM=1,NMM)
              ENDDO
            ENDDO !nx
          ENDDO !nr

          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                FMT='(/'' Region nr='',I2,''   Element ne='',I4)'
                READ(IUNIT,FMT) NRR,NEE
                FMT='('' NHE(ne,nx)='',I2)'
                READ(IUNIT,FMT) NHE(ne,nx)
                DO nhx=1,NHE(ne,nx)
                  nh=NH_LOC(nhx,nx)
                  FMT='('' NBH(nh='',I2,'',nc,ne):'',20I3)'
                  READ(IUNIT,FMT) nhh,NBH(nh,nc,ne)
                ENDDO !nh

                DO nb=1,NBFT
                  DO nhx=1,NHE(ne,nx)
                    nh=NH_LOC(nhx,nx)
                    FMT='('' NVHE(nn=1..,nb='',I2,'',nh='',I2,'
     '                //''',ne):'',20I2)'
                    READ(IUNIT,FMT) NBB,NHH,
     '                (NVHE(nn,nb,nh,ne),nn=1,NNT(nb))
                  ENDDO !nh
                ENDDO !nb

                FMT='(/'' Solution type nx='',I2,)'
                READ(IUNIT,FMT) NX1
                DO nhx=1,NHE(ne,nx)
                  nh=NH_LOC(nhx,nx)
                  DO nrc=0,2
                    nb=NBH(nh,nc,ne)
                    DO nn=1,NNT(nb)
                      DO nv=1,NVHE(nn,nb,nh,ne)
                        FMT='('' NYNE(na=1..,nh='',I2,'',nrc='',I1,'
     '                    //'''nc='',I2,'',ne):'',10I6)'
                        READ(IUNIT,FMT) NHH,nrcc,ncc,
     '                    (NYNE(na,nh,nrc,nc,ne),na=1,NAT(nb))
                      ENDDO !nv
                    ENDDO !nn
                  ENDDO !nrc
                ENDDO !nh
C news MPN 14Jun2000: added nx index to NW
                FMT='('' NW(ne,1,nx)='',I3)'
                READ(IUNIT,FMT) NW(ne,1,nx)
                FMT='('' CE(nm=1..,ne,nx):'',/(3E25.17))'
                READ(IUNIT,FMT) (CE(nm,ne,nx),nm=1,NMM)
              ENDDO !noelem
            ENDDO !nx
          ENDDO !nr

          nrc=0 !Temporary AJP 1-12-94
          DO nx1=1,NX_LIST(0)
            nx=NX_LIST(nx1)
            DO nr=1,NRT
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO nhx=1,NHP(np,nr,nx)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,nc,nr)
                    DO nk=1,NKH(nh,np,nc,nr)
                      ny=NYNP(nk,nv,nh,np,nrc,nc,nr)
                      FMT='(/'' FIX('',I5,'',iy=1..,nx='',I1,''): '''
     '                  //',5L1)'
                      READ(IUNIT,FMT) NYY,NXX,
     '                  (FIX(ny,iy,nx),iy=1,NIYFIXM)
                      FMT='('' YP('',I5,'',iy=1..,nx='',I1,''): '',/,'
     '                  //'(3E25.17))'
                      READ(IUNIT,FMT) NYY,NXX,(YP(ny,iy,nx),iy=1,NIYM)
C!!! MPN needs fixing
                      nrc=1 !temporary
                      FMT='('' NONY(0,'',I5,'',nr='',I1,'',nx='','
     '                  //'I1,''):'',I3)'
                      READ(IUNIT,FMT) NYY,NRR,NXX,NONY(0,ny,nrc,nr,nx)
                      FMT='('' NONY(noy=1..,'',I5,'',nr='','
     '                  //'I1,'',nx='',I1,''):'',9I3)'
                      READ(IUNIT,FMT) NYY,NRR,NXX,
     '                  (NONY(NOY,ny,nrc,nr,nx),
     '                  NOY=1,NONY(0,ny,nrc,nr,nx))
                      FMT='('' CONY(noy=1..,'',I5,'',nr='','
     '                  //'I1,'',nx='',I1,''):'',(3E25.17))'
                      READ(IUNIT,FMT) NYY,NRR,NXX,
     '                  (CONY(NOY,ny,nrc,nr,nx),
     '                  NOY=1,NONY(0,ny,nrc,nr,nx))
C CPB 13/11/94 Temporarily fixed NYNO will need to be redone when
C NYNO is calculated and used properly
                      FMT='('' NYNO(NONY(1..,ny,nr,nx),nr='',I1,'
     '                  //''',nx='',I1,'') = '',10I5)'
                      READ(IUNIT,FMT) NRR,NXX,
     '                  (NYNO(1,NONY(NOY,ny,nrc,nr,nx),nrc,nr,nx),
     '                  NOY=1,NONY(0,ny,nrc,nr,nx))
                    ENDDO !nk
                  ENDDO !nv
                ENDDO !nh
              ENDDO !nonode

              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                DO nhx=1,NHE(ne,nx)
                  nh=NH_LOC(nhx,nx)
                  nb=NBH(nh,nc,ne)
                  DO nn=1,NNT(nb)
                    DO nv=1,NVHE(nn,nb,nh,ne)
                      DO na=1,NAT(nb)
                  ny=NYNE(na,nh,nrc,nc,ne)
                  FMT='(/'' FIX('',I5,'',iy=1..,nx='',I1,''): '','
     '                    //'5L1)'
                        READ(IUNIT,FMT) NYY,NXX,
     '                    (FIX(ny,iy,nx),iy=1,NIYFIXM)
                        FMT='('' YP('',I5,'',iy=1..,nx='',I1,''): '',/,'
     '                    //'(3E25.17))'
                  READ(IUNIT,FMT) NYY,NXX,(YP(ny,iy,nx),iy=1,NIYM)
C!!! MPN needs fixing
                        nrc=1 !temporary
                        FMT='('' NONY(0,'',I5,'',nr='',I1,'',nx='','
     '                    //'I1,''):'',I3)'
                  READ(IUNIT,FMT) NYY,NRR,NXX,NONY(0,ny,nrc,nr,nx)
                  FMT='('' NONY(noy=1..,'',I5,'',nr='','
     '                    //'I1,'',nx='',I1,''):'',9I3)'
                  READ(IUNIT,FMT) NYY,NRR,NXX,
     '              (NONY(NOY,ny,nrc,nr,nx),
     '                    NOY=1,NONY(0,ny,nrc,nr,nx))
                        FMT='('' CONY(noy=1..,'',I5,'',nr='','
     '                    //'I1,'',nx='',I1,''):'',(3E25.17))'
                  READ(IUNIT,FMT) NYY,NRR,NXX,
     '              (CONY(NOY,ny,nrc,nr,nx),
     '                    NOY=1,NONY(0,ny,nrc,nr,nx))
                      ENDDO !na
                    ENDDO !nv
                  ENDDO !nn
                ENDDO !nh
              ENDDO !noelem
            ENDDO !nr
          ENDDO !nx

          DO nr=1,NRT
            IF(KTYP57(nr).EQ.2.OR.KTYP57(nr).EQ.4.OR.KTYP57(nr).EQ.5)
     '        THEN !Press incr
              FMT='(/'' Region '',I1,'':Pressure Boundary Conditions'')'
              READ(IUNIT,FMT) NRR
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                FMT='(''   FIXP(iface=1..2,ne='',I3,''): '',2L1)'
                READ(IUNIT,FMT) NEE,(FIXP(iface,ne),iface=1,2)
                FMT='(''   PE(iface=1..2,ne='',I3,''): '',2E25.17)'
                READ(IUNIT,FMT) NEE,(PE(iface,ne),iface=1,2)
                FMT='(''   PF(iface=1..2,ne='',I3,''): '',2E25.17)'
                READ(IUNIT,FMT) NEE,(PF(iface,ne),iface=1,2)
              ENDDO !noelem
            ENDIF !ktyp57(nr)
          ENDDO !nr

          FMT='(/'' FEXT(i=1..,ng,ne):'')'
          READ(IUNIT,FMT)
          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(NBH(NH_LOC(1,nx),nc,ne).GT.0) THEN
                  DO ng=1,NGT(NBH(NH_LOC(1,nx),nc,ne))
                    FMT='('' FEXT(i,'',I3,'','',I4,''): '',3E25.17,'
     '                //'/(4E25.17))'
                    READ(IUNIT,FMT) mg,ne,(FEXT(i,ng,ne),i=1,NIFEXTM)
                  ENDDO !ng
                ENDIF
              ENDDO !noelem
              IF(KTYP53(nr).EQ.3)THEN !Active stress included
                FMT='(/'' CALL_ACTI='',L1)'
                READ(IUNIT,FMT) CALL_ACTI
              ENDIF
            ENDDO !nx
          ENDDO !nr

          IF(KTYP1.EQ.8.AND.KTYP3.EQ.2) THEN !activation variables
            DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
                nx=NX_LIST(nxx)
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  FMT='(/'' Element ne='',I4)'
                  READ(IUNIT,FMT) ne
                  nb=NBH(NH_LOC(1,nx),nc,ne)
                  IF(nb.GT.0) THEN
                    FMT='(/'' ITHRES(1,ng,ne):'',/(1X,75I1))'
                    READ(IUNIT,FMT) (ITHRES(1,ng,ne),NG=1,NGT(nb))
                    FMT='(/'' ITHRES(2,ng,ne):'',/(1X,75I1))'
                    READ(IUNIT,FMT) (ITHRES(2,ng,ne),NG=1,NGT(nb))
                    FMT='(/'' THRES(1,ng,ne):'',/(1X,10E12.4))'
                    READ(IUNIT,FMT) (THRES(1,ng,ne),NG=1,NGT(nb))
                    FMT='(/'' THRES(2,ng,ne):'',/(1X,10E12.4))'
                    READ(IUNIT,FMT) (THRES(2,ng,ne),NG=1,NGT(nb))
                    FMT='(/'' THRES(3,ng,ne):'',/(1X,10E12.4))'
                    READ(IUNIT,FMT) (THRES(3,ng,ne),NG=1,NGT(nb))
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF

          FMT='(/'' CALL_EQUA='',L1,'' CALL_INIT='',L1,'
     '      //'  '' CALL_MATE='',L1,'' CALL_SOLV='',L1)'
          READ(IUNIT,FMT) CALL_EQUA,CALL_INIT,CALL_MATE,CALL_SOLV
          FMT='('' CALL_DATA='',L1,'' CALL_FIT ='',L1,'
     '      //' '' CALL_GROW='',L1,'' CALL_MOTI='',L1)'
          READ(IUNIT,FMT) CALL_DATA,CALL_FIT ,CALL_GROW,CALL_MOTI
          IF(KTYP90.GT.0)THEN
            FMT='('' CALL_COUP='',L1)'
            READ(IUNIT,FMT) CALL_COUP
          ENDIF

          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF(CALL_SOLV.AND.ITYP4(nr,nx).EQ.2)THEN !BE
                FMT='(/'' ADAPINT='',L1,'' COMPLEX='',L1,'
     '            //'  '' HERMITE='',L1,'' HYP='',L1)'
                READ(IUNIT,FMT) ADAPINT,COMPLEX,HERMITE,HYP
              ENDIF
            ENDDO
          ENDDO

          IF(KTYP1.EQ.1) THEN !linear elasticity
            IF(ETYP(7).OR.ETYP(8)) THEN !shell elements
              nb=NBJ(1,NEELEM(1,1)) !must be cubic Hermite
              CALL ASSERT(NKT(0,nb).EQ.4,
     '          '>>Geometry must be cubic Hermite',ERROR,*9999)
c D3PG is not passed through - PJH 5-Jan-1991
c             CALL GAUS20(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
c    '          nb,NGAP(1,nb),D3PG(1,1,1,nb))
            ENDIF
          ENDIF

        ENDIF !dependent

            !************ call gauss subroutines ***********

        DO nb=1,NBFT
          IF(NBC(nb).EQ.1) THEN      !Lagrange/Hermite tensor prod basis
            CALL GAUSS1(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '        NGAP(1,nb),PG(1,1,1,nb),WG(1,nb),XIG(1,1,nb),ERROR,*9999)
            NBASEF(nb,0)=1
          ELSE IF(NBC(nb).EQ.2) THEN !Simplex/Serendipity/Lagrange basis
            IF(IBT(2,1,nb).NE.4)THEN !Not a hermite simplex
              CALL GAUSS2(IBT,INP,nb,NGAP(1,nb),PG,WG,XIG,ERROR,*9999)
              NBASEF(nb,0)=1
            ELSE
!             !Hermite simplex
              CALL GAUSS2_HERMITE(IDO,INP,nb,NGAP(1,nb),PG,WG,XIG,
     '          ERROR,*9999)
            ENDIF
          ELSE IF(NBC(nb).EQ.3) THEN !B-spline tensor product basis
            WRITE(OP_STRING,'('' Warning: No call to Gauss3 from '','
     '        //'''IOGEOM!!'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            NBASEF(nb,0)=1
          ELSE IF(NBC(nb).EQ.4) THEN !Fourier Series tensor prod basis
            CALL GAUSS4(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '        NGAP(1,nb),PG(1,1,1,nb),WG(1,nb),XIG(1,1,nb),ERROR,*9999)
            NBASEF(nb,0)=1
          ELSE IF(NBC(nb).EQ.5) THEN !Boundary element basis
            NBTOP=NBFT !Calculate highest basis function number so far
            DO nb2=1,nb-1
              NBTOP=NBTOP+NBASEF(nb2,0)-1
            ENDDO
            CALL GAUSS5(IBT,IDO,INP,nb,NBTOP,NDET,NGAP,
     '        DET,PG,WG,XIG,ERROR,*9999)
          ELSE IF(NBC(nb).EQ.6) THEN !Singular basis function
            NBTOP=NBFT !Calculate highest basis function number so far
            DO nb2=1,NB-1
              NBTOP=NBTOP+NBASEF(nb2,0)-1
            ENDDO
            CALL GAUSS6(IBT,IDO,INP,nb,NBTOP,NDET,NGAP,
     '        DET,PG,WG,XIG,ERROR,*9999)
          ELSE IF(NBC(nb).EQ.7) THEN !Extended Lagrange basis
            CALL GAUSS7(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '        NGAP(1,nb),PG(1,1,1,nb),WG(1,nb),XIG(1,1,nb),ERROR,*9999)
            NBASEF(nb,0)=1
          ENDIF
          IF(NBC(nb).LE.2) THEN !Auxillary basis
            IF(NAT(nb).GT.0) THEN
              CALL GAUSS8(NAN(1,1,nb),nb,NGAP(1,nb),PG(1,1,1,nb),
     '          WG(1,nb),XIG(1,1,nb),ERROR,*9999)
              NBASEF(nb,0)=1
            ENDIF
          ENDIF
        ENDDO

!********************* write iod file ********************************

      ELSE IF(COMAND.EQ.'WRITE') THEN
!nx_list
        FMT='('' NX_LIST(0..): '',10I3)'
        WRITE(IUNIT,FMT) NX_LIST(0),(NX_LIST(nxx),nxx=1,NX_LIST(0))
!n*t params
        FMT='('' NBT='',I6,'' NDT='',I6,'' NFT='',I6,'
     '     //''' NJT='',I6,'' NLT='',I6,'' NQT='',I6,/'
     '     //''' NRT='',I6,'' NBFT='',I5)'
        WRITE(IUNIT,FMT) NBT,NDT,NFT,NJT,NLT,NQT,NRT,NBFT

        DO nxx=1,NX_LIST(0)
          nx=NX_LIST(nxx)
          FMT='(/'' NZT(1,nx): '',I6,I6)'
          WRITE(IUNIT,FMT) nx,NZT(1,nx)
        ENDDO !nxx
!jtypes
        FMT='(/'' JTYP1 ='',I2,'' JTYP2A ='',I2,'' JTYP2B ='',I2,/'
     '      //''' JTYP2C ='',I2,'' JTYP3 ='',I2,/'
     '      //''' JTYP4 ='',I2,'' JTYP5 ='',I2,'' JTYP6 ='',I2,/'
     '      //''' JTYP7 ='',I2,'' JTYP8 ='',I2,'' JTYP9 ='',I2,'
     '      //''' JTYP10='',I2,'' JTYP11='',I2,'' JTYP12='',I2,/'
     '      //''' JTYP13='',I2,'' JTYP14='',I2,'' JTYP15='',I2)'
        WRITE(IUNIT,FMT) JTYP1,JTYP2A,JTYP2B,JTYP2C,JTYP3,
     '                   JTYP4,JTYP5,JTYP6,
     '                   JTYP7,JTYP8,JTYP9,JTYP10,JTYP11,JTYP12,
     '                   JTYP13,JTYP14,JTYP15
!net,npt
        FMT='(/'' NET(0)='',I6,'' NPT(0)='',I6)'
        WRITE(IUNIT,FMT) NET(0),NPT(0)
!region params
        DO nr=1,NRT
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            FMT='(/'' nx= '',I6)'
            WRITE(IUNIT,FMT) nx
            FMT='(/'' Region '',I1,'': ITYP1 ='',I2,'' ITYP2 ='',I2,'
     '                         //''' ITYP3 ='',I2,'' ITYP4 ='',I2,/10X,'
     '                         //''' ITYP5 ='',I2,'' ITYP6 ='',I2,'
     '                         //''' ITYP7 ='',I2,'' ITYP8 ='',I2,/10X,'
     '                         //''' ITYP9 ='',I2,'' ITYP12='',I2)'
            WRITE(IUNIT,FMT) nr,ITYP1(nr,nx),ITYP2(nr,nx),ITYP3(nr,nx),
     '        ITYP4(nr,nx),ITYP5(nr,nx),ITYP6(nr,nx),ITYP7(nr,nx),
     '        ITYP8(nr,nx),ITYP9(nr,nx),ITYP12(nr,nx)
          ENDDO !nx
          FMT='(/'' Region '',I1,'': ITYP10='',I2,'
     '                       //''' ITYP11='',I2,/10X,'
     '                       //''' ITYP13='',I2,'' ITYP14='',I2,'
     '                       //''' ITYP15='',I2)'
C KAT 31Oct97: ITYP15 and ITYP16 should be in the nx loop but they are
C not set up properly yet to avoid corrupting present IOD files.
          WRITE(IUNIT,FMT) nr,ITYP10(nr),ITYP11(nr),ITYP13(nr),
     &      ITYP14(nr),ITYP15(nr,1)
          FMT='(10X,'' NET('',I2,'')='',I6,'' NPT('',I2,'')='',I6)'
          WRITE(IUNIT,FMT) nr,NET(nr),nr,NPT(nr)
!nj_loc
          FMT='(/'' NJ_LOC(1..3,0,nr): '',3I3)'
          WRITE(IUNIT,FMT) (NJ_LOC(njj1,0,nr),njj1=1,3)
          FMT='( '' NJ_LOC(1,1..,nr):  '',12I3)'
          WRITE(IUNIT,FMT) (NJ_LOC(1,njj2,nr),njj2=1,NJ_LOC(1,0,nr))
          FMT='( '' NJ_LOC(2,1..,nr):  '',12I3)'
          WRITE(IUNIT,FMT) (NJ_LOC(2,njj2,nr),njj2=1,NJ_LOC(2,0,nr))
          FMT='( '' NJ_LOC(3,1..,nr):  '',12I3)'
          WRITE(IUNIT,FMT) (NJ_LOC(3,njj2,nr),njj2=1,NJ_LOC(3,0,nr))
        ENDDO
!basis function info
        DO nb=1,NBFT
          FMT='(/'' Basis type nb='',I2)'
          WRITE(IUNIT,FMT) NB
          FMT='(/'' NAT(nb)='',I3,''  NBI(nb)='',I3,''  NFE(nb)='',I3,'
     '      //'''  NGT(nb)='',I3,''  NIT(nb)='',I3,/'
     '       //''' NLE(nb)='',I3,''  NNT(nb)='',I3,''  NST(nb)='',I3,'
     '      //'''  NUT(nb)='',I3,''  NBC(nb)='',I3,/'
     '       //''' NABTYP(nb)='',I3,'' NMGT='',I3)'
          WRITE(IUNIT,FMT) NAT(nb),NBI(nb),NFE(nb),NGT(nb),NIT(nb),
     '      NLE(nb),NNT(nb),NST(nb),NUT(nb),NBC(nb),NABTYP(nb),NMGT
          FMT='('' NKT(nn=0..,nb):'',18I3)'
          WRITE(IUNIT,FMT) (NKT(nn,nb),NN=0,NNT(nb))
          FMT='('' IBT(i=1..,ni=1..,nb):'',40I2)'
          WRITE(IUNIT,FMT) ((IBT(i,ni,nb),i=1,3),NI=1,NIT(nb))
          DO nn=1,NNT(nb)
            FMT='('' IDO(nk=1..,nn='',I2,'',ni=1..,nb):'',40I2)'
            WRITE(IUNIT,FMT) nn,((IDO(nk,nn,ni,nb),nk=1,NKT(nn,nb)),
     '        ni=0,NIT(nb))
          ENDDO
          FMT='('' INP(nn=1..,ni=1..,nb):'',40I2)'
          WRITE(IUNIT,FMT) ((INP(nn,ni,nb),NN=1,NNT(nb)),NI=1,NIT(nb))
          FMT='('' NAN(ni=1..,na=1..,nb):'',40I2)'
          WRITE(IUNIT,FMT) ((NAN(ni,na,nb),NI=1,NIT(nb)),NA=1,NAT(nb))
          FMT='('' NGAP(ni=1..,nb):'',20I4)'
          WRITE(IUNIT,FMT) (NGAP(ni,nb),NI=1,NIT(nb))
          FMT='('' NNL(i=0..,j=1..,nb):'',/(12I4))'
          WRITE(IUNIT,FMT) ((NNL(i,j,nb),i=0,4),j=1,12)
          DO nf=1,NFE(nb)
            FMT='('' NNF(0,nf='',I1,'',nb):'',I3)'
            WRITE(IUNIT,FMT) nf,NNF(0,nf,nb)
            FMT='('' NNF(1..,nf='',I1,'',nb):'',17I3)'
            WRITE(IUNIT,FMT) nf,(NNF(i,nf,nb),i=1,NNF(0,nf,nb))
            FMT='('' NKEF(0,nn=1..,nf='',I1,'',nb):'',16I2)'
            WRITE(IUNIT,FMT) nf,(NKEF(0,j,nf,nb),j=1,NNF(0,nf,nb))
            FMT='('' NKEF(1..,nn=1..,nf='',I1,'',nb):'',16I2)'
            WRITE(IUNIT,FMT) nf,((NKEF(i,j,nf,nb),i=1,NKEF(0,j,nf,nb)),
     '        j=1,NNF(0,nf,nb))
          ENDDO !nfe
          IF(NBC(nb).EQ.4) THEN
            FMT='('' OMEGA='',E25.17)'
            WRITE(IUNIT,FMT) OMEGA
          ENDIF
          IF((NBC(nb).EQ.5).OR.(NBC(nb).EQ.6)) THEN !BE basis function
            FMT='('' NGLIMITS(ni=1..,nb,1..2):'',/(20I4))'
            WRITE(IUNIT,FMT) ((NGLIMITS(ni,nb,i),NI=1,NIT(nb)),i=1,2)
          ENDIF
        ENDDO
!focus
        FMT='(/'' Focus:'',E25.17)'
        WRITE(IUNIT,FMT) FOCUS
!nodes
        WRITE(IUNIT,'(/30X,''*** Nodes ***'')')
        FMT='('' NPNODE(0,0):    '',I5)'
        WRITE(IUNIT,FMT) NPNODE(0,0)
        DO nr=1,NRT
          FMT='(/'' NPNODE(0,nr='',I2,''):'',I5)'
          WRITE(IUNIT,FMT) NR,NPNODE(0,nr)
          FMT='('' NPNODE(nonode=1..,nr='',I2,''):''/,(15I5))'
          WRITE(IUNIT,FMT) NR,(NPNODE(nonode,nr),nonode=1,NPNODE(0,nr))
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            FMT='(/'' Node np='',I4)'
            WRITE(IUNIT,FMT) NP
            FMT='('' NVJP(nj=1..,np):'',12I3)'
            WRITE(IUNIT,FMT) ((NVJP(NJ_LOC(njj1,njj2,nr),np),
     '        njj2=1,NJ_LOC(njj1,0,nr)),njj1=1,3)
            FMT='('' NKJ(nj=1..,np):'',20I3)'
            WRITE(IUNIT,FMT) ((NKJ(NJ_LOC(njj1,njj2,nr),np),
     '        njj2=1,NJ_LOC(njj1,0,nr)),njj1=1,3)
            DO njj1=1,3  !geometry/fibres/field
              DO njj2=1,NJ_LOC(njj1,0,nr)
                nj=NJ_LOC(njj1,njj2,nr)
                DO nv=1,NVJP(nj,np)
                  FMT='('' XP(nk=1..,nv='',I2,'',nj='',I1,'',np):'','
     '              //'3E25.17,:/3E25.17)'
                  WRITE(IUNIT,FMT)
     '              nv,nj,(XP(nk,nv,nj,np),NK=1,NKJ(nj,np))
                ENDDO !nv
              ENDDO !njj2
            ENDDO !njj1
            FMT='('' NP_INTERFACE(np,0): '',I2)'
            WRITE(IUNIT,FMT) NP_INTERFACE(np,0)
            FMT='('' NP_INTERFACE(np,i=1..): '',3I5)'
            WRITE(IUNIT,FMT)
     '        (NP_INTERFACE(np,i),i=1,NP_INTERFACE(np,0))
          ENDDO !nonode
        ENDDO !nr
!elements
        WRITE(IUNIT,'(/30X,''*** Elements ***'')')
        FMT='('' NEELEM(0,0):    '',I5)'
        WRITE(IUNIT,FMT) NEELEM(0,0)
        DO nr=1,NRT
          FMT='(/'' NEELEM(0,nr='',I2,''):'',I5)'
          WRITE(IUNIT,FMT) nr,NEELEM(0,nr)
          FMT='('' ITYP10(nr)='',I2)'
          WRITE(IUNIT,FMT) ITYP10(nr)
          FMT='('' NEELEM(noelem=1..,nr='',I2,''):''/,(15I5))'
          WRITE(IUNIT,FMT) nr,(NEELEM(noelem,nr),noelem=1,NEELEM(0,nr))
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            FMT='(/'' Element ne='',I4)'
            WRITE(IUNIT,FMT) NE
            FMT='('' NRE(ne)='',I2)'
            WRITE(IUNIT,FMT) NRE(ne)
            FMT='('' NBJ(nj=1..,ne):'',20I3)'
            WRITE(IUNIT,FMT) ((NBJ(NJ_LOC(njj1,njj2,nr),ne),
     '        njj2=1,NJ_LOC(njj1,0,nr)),njj1=1,3)
            FMT='('' NLL(i=1..,ne):'',20I4)'
            WRITE(IUNIT,FMT) (NLL(i,ne),i=1,12)
            FMT='('' NFF(i=1..,ne):'',20I4)'
            WRITE(IUNIT,FMT) (NFF(i,ne),i=1,6)
            DO nb=1,NBFT
              FMT='(/'' NPNE(nn=1..,nb='',I2,'',ne):'',20I4)'
              WRITE(IUNIT,FMT) nb,(NPNE(nn,nb,ne),NN=1,NNT(nb))
C              DO nn=1,NNT(nb)
C                FMT='('' NKE(nk=1..,nn='',I2,'',nb='',I2,'',ne):'','
C     '            //'(20I4))'
C                WRITE(IUNIT,FMT) nn,nb,
C     '            (NKE(nk,nn,nb,ne),NK=1,NKT(nn,nb))
C              ENDDO
              FMT='('' SE(ns=1..,nb='',I2,'',ne):'',/(3E25.17))'
              WRITE(IUNIT,FMT) nb,(SE(ns,nb,ne),NS=1,NST(nb))
              DO njj1=1,3  !geometry/fibres/field
                DO njj2=1,NJ_LOC(njj1,0,nr)
                  nj=NJ_LOC(njj1,njj2,nr)
                  FMT='('' NVJE(nn=1..,nb='',I2,'',nj='',I2,'
     '              //''',ne):'',20I2)'
                  WRITE(IUNIT,FMT) nb,nj,
     '              (NVJE(nn,nb,nj,ne),NN=1,NNT(nb))
                ENDDO !njj1
              ENDDO !njj2
            ENDDO !nb
!xa
            FIRST=.TRUE.
            DO njj1=1,3  !geometry/fibres/field
              DO njj2=1,NJ_LOC(njj1,0,nr)
                nj=NJ_LOC(njj1,njj2,nr)
                nb=NBJ(nj,ne)
                DO nn=1,NNT(nb)
                  FMT='('' NKJE(nk=1..,nn='',I2,'',nj='',I2,'',ne):'','
     '              //'(20I4))'
                  WRITE(IUNIT,FMT) nn,nj,
     '              (NKJE(nk,nn,nj,ne),nk=1,NKT(nn,nb))
                ENDDO
                IF(NAT(nb).NE.0) THEN
                  IF(FIRST) THEN
                    FMT='()'
                    WRITE(IUNIT,FMT)
                    FIRST=.FALSE.
                  ENDIF !first
                  FMT='('' XA(na=1..,nj='',I2,'',ne):'',2E25.17'//
     '              ':/(3E25.17))'
                  WRITE(IUNIT,FMT) nj,(XA(na,nj,ne),na=1,NAT(nb))
                ENDIF !NAT.NE.0
              ENDDO !njj1
            ENDDO !njj2

          ENDDO !noelem
        ENDDO !nr
!nxi
        WRITE(IUNIT,'(/30X,''***  NXI  ***'')')
C new MPN 22Apr97: wrong ne index for NXI
        DO nr=1,NRT
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO ni=-NIM,NIM
              FMT='('' NXI(ni='',I4,'',0-4,ne='',I4,''):'',20I4)'
              WRITE(IUNIT,FMT) ni,ne,(NXI(ni,i,ne),i=0,4)
            ENDDO !ni
          ENDDO !noelem
        ENDDO !nr
C old
C        DO ne=1,NEELEM(0,0)
C          DO ni=-NIM,NIM
C            FMT='('' NXI(ni='',I4,'',0-4,ne='',I4,''):'',20I4)'
C            WRITE(IUNIT,FMT) ni,ne,(NXI(ni,i,ne),i=0,4)
C          ENDDO !ni
C        ENDDO !ne
C!xa
C        FMT='(/'' XA(1,nj=1..,ne=1..):'',/(3E25.17))'
C        WRITE(IUNIT,FMT) ((XA(1,nj,ne),nj=1,NJM),ne=1,NET(0))
!faces
        WRITE(IUNIT,'(/30X,''*** Faces ***'')')
        DO nf=1,NFT
          FMT='(/'' Face nf='',I4)'
          WRITE(IUNIT,FMT) NF
          FMT='('' NPF(i=1..,nf):'',9I6)'
          WRITE(IUNIT,FMT) (NPF(i,nf),i=1,9)
          FMT='('' NLF(i=1..,nf):'',4I4)'
          WRITE(IUNIT,FMT) (NLF(i,nf),i=1,4)
          FMT='('' NBJF(nj=1..,nf):'',(3I4))'
          WRITE(IUNIT,FMT) (NBJF(nj,nf),nj=1,NJT)
        ENDDO !nf
!lines
        WRITE(IUNIT,'(/30X,''*** Lines ***'')')
        DO nl=1,NLT
          FMT='(/'' Line nl='',I4)'
          WRITE(IUNIT,FMT) nl
          DO nj=0,3
            FMT='('' NPL(i=1..,nj='',I1,'',nl):'',5(1X,I4))'
            WRITE(IUNIT,FMT) nj, (NPL(i,nj,nl),i=1,5)
          ENDDO
          FMT='('' NEL(i=0..,nl):'',10(1X,I4),/,15X,10(1X,I4))'
          WRITE(IUNIT,FMT) NEL(0,nl),(NEL(i,nl),i=1,NEL(0,nl))
          FMT='('' DL(i=1..,nl):'',3E25.17)'
          WRITE(IUNIT,FMT) (DL(i,nl),i=1,3)
        ENDDO !nl
!call logicals
        FMT='(/'' CALL_BASE='',L1,'' CALL_ELEM='',L1,'
     '    //'  '' CALL_LINE='',L1,'' CALL_NODE='',L1,'
     '    //'  '' CALL_FIBR='',L1,'' CALL_ELFB='',L1)'
        WRITE(IUNIT,FMT) CALL_BASE,CALL_ELEM,CALL_LINE,CALL_NODE,
     '    CALL_FIBR,CALL_ELFB

            !***********  dependent variables *************

C        nx=1 !Temporary
        DEPENDENT=.FALSE.
        DO nxx=1,NX_LIST(0)
          nx=NX_LIST(nxx)
          DO nr=1,NRT
            IF(ITYP1(nr,nx).GT.1) DEPENDENT=.TRUE.
          ENDDO
        ENDDO

        IF(DEPENDENT) THEN
!ktypes
          FMT='(/'' KTYP1 ='',I2,'' KTYP2 ='',I2,'' KTYP3 ='',I2,'
     '      //'  '' KTYP4 ='',I2,'' KTYP5 ='',I2,'' KTYP6 ='',I2,/'
     '      //'  '' KTYP7 ='',I2,'' KTYP8 ='',I2,'' KTYP9 ='',I2,'
     '      //  ''' KTYP10='',I2,'' KTYP11='',I2,'' KTYP12='',I2,/'
     '      //'  '' KTYP13='',I2,'' KTYP14='',I2,'' KTYP15='',I2,'
     '      //'  '' KTYP16='',I2,'' KTYP17='',I2)'
          WRITE(IUNIT,FMT) KTYP1,KTYP2,KTYP3,KTYP4,KTYP5,KTYP6,
     '      KTYP7,KTYP8,KTYP9,KTYP10,KTYP11,KTYP12,
     '      KTYP13,KTYP14,KTYP15,KTYP16,KTYP17
          FMT='( '' KTYP19='',I2,'' KTYP1A='',I2,'' KTYP1B='',I2,'
     '      //'  '' KTYP1C='',I2,'' KTYP1D='',I2,'' KTYP1E='',I2,/'
     '      //'  '' KTYP20='',I2,'' KTYP21='',I2,'' KTYP22='',I2,'
     '      //'  '' KTYP23='',I2,'' KTYP24='',I2,'' KTYP25='',I2,/'
     '      //'  '' KTYP26='',I2,'' KTYP27='',I2,'' KTYP28='',I2,'
     '      //'  '' KTYP29='',I2)'
          WRITE(IUNIT,FMT) KTYP19,KTYP1A,KTYP1B,KTYP1C,KTYP1D,KTYP1E,
     '      KTYP20,KTYP21,KTYP22,KTYP23,KTYP24,KTYP25,
     '      KTYP26,KTYP27,KTYP28,KTYP29
          FMT='(/'' KTYP30='',I2,'' KTYP31='',I2,'' KTYP32='',I2,'
     '      //'  '' KTYP33='',I2,'' KTYP34='',I2,'' KTYP35='',I2)'
          WRITE(IUNIT,FMT) KTYP30,KTYP31,KTYP32,KTYP33,KTYP34,
     '      KTYP35
          FMT='(/'' KTYP40='',I2,'' KTYP41='',I2,'' KTYP42='',I2,'
     '       //' '' KTYP43='',I2,'' KTYP44='',I2,'' KTYP45='',I2)'
          WRITE(IUNIT,FMT) KTYP40,KTYP41,KTYP42,KTYP43,KTYP44,
     '      KTYP45
          DO nr=1,NRT
            FMT='(/'' Region '',I1)'
            WRITE(IUNIT,FMT) nr
            FMT='(/'' KTYP50='',I2,'' KTYP51='',I2,'' KTYP52='',I2,'
     '        //'  '' KTYP53='',I2,'' KTYP54='',I2,'' KTYP55='',I2,'
     '        //' /'' KTYP56='',I2,'' KTYP57='',I2,'' KTYP58='',I2,'
     '        //'  '' KTYP59='',I2,'' KTYP5A='',I2,'' KTYP5B='',I2,'
     '        //' /'' KTYP5C='',I2,'' KTYP5D='',I2,'' KTYP5E='',I2)'
            WRITE(IUNIT,FMT) KTYP50(nr),KTYP51(nr),KTYP52(nr),
     '        KTYP53(nr),KTYP54(nr),KTYP55(nr),KTYP56(nr),
     '        KTYP57(nr),KTYP58(nr),KTYP59(nr),KTYP5A(nr),
     '        KTYP5B(nr),KTYP5C(nr),KTYP5D(nr),KTYP5E(nr)
          ENDDO
          FMT='('' IL_thickness ='',I2,'' IL_sarcomere='',I2,'
     '      //'/'' IL_time_delay='',I2,'' IL_fluid_conductivity='',I2)'
          WRITE(IUNIT,FMT) IL_thickness,IL_sarcomere,
     '      IL_time_delay,IL_fluid_conductivity
          FMT='(/'' KTYP60='',I2,'' KTYP61='',I2,'' KTYP62='',I2,'
     '      //'  '' KTYP63='',I2,'' KTYP64='',I2,'' KTYP65='',I2)'
          WRITE(IUNIT,FMT) KTYP60,KTYP61,KTYP62,KTYP63,KTYP64,
     '      KTYP65
          FMT='(/'' KTYP70='',I2,'' KTYP71='',I2,'' KTYP72='',I2,'
     '      //'  '' KTYP73='',I2,'' KTYP74='',I2,'' KTYP75='',I2)'
          WRITE(IUNIT,FMT) KTYP70,KTYP71,KTYP72,KTYP73,KTYP74,
     '      KTYP75
          FMT='(/'' KTYP90='',I2,'' KTYP91='',I2,'' KTYP92='',I2)'
          WRITE(IUNIT,FMT) KTYP90,KTYP91,KTYP92
          FMT='( '' KTYP93(1,nr)='',9I2,'' KTYP93(2,nr)='',9I2)'
          WRITE(IUNIT,FMT) (KTYP93(1,nr),nr=1,NRT),(KTYP93(2,nr),
     '      nr=1,NRT)
          FMT='( ''KTYP94='',I2,'' KTYPMBC='',I2)'
          WRITE(IUNIT,FMT) KTYP94,KTYPMBC
!iwrit
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            FMT='(/'' nx= '',I6)'
            WRITE(IUNIT,FMT) nx
            FMT='(/'' IWRIT1(nr=1..,nx):'',9I2)'
            WRITE(IUNIT,FMT) (IWRIT1(nr,nx),nr=1,NRT)
            FMT='( '' IWRIT2(nr=1..,nx):'',9I2)'
            WRITE(IUNIT,FMT) (IWRIT2(nr,nx),nr=1,NRT)
            FMT='( '' IWRIT3(nr=1..,nx):'',9I2)'
            WRITE(IUNIT,FMT) (IWRIT3(nr,nx),nr=1,NRT)
            FMT='( '' IWRIT4(nr=1..,nx):'',9I2)'
            WRITE(IUNIT,FMT) (IWRIT4(nr,nx),nr=1,NRT)
          ENDDO
!time params
          FMT='(/'' DT='',E25.17,'' TINCR='',E25.17,'' T0='',E25.17,'
     '      //''' T1='',E25.17,/'' THETA(1)='',E25.17,'' THETA(2)='','
     '      //'E25.17,'//''' THETA(3)='',E25.17)'
          WRITE(IUNIT,FMT) DT,TINCR,T0,T1,THETA(1),THETA(2),THETA(3)
!logicals
          FMT='(/'' INIT='',L1,'' LUMP='',L1,'' PROMPT='',L1,'
     '      //''' LRESID='',L1,'' LSTEP='',L1,'' RESTAR='',L1)'
          WRITE(IUNIT,FMT) INIT,LUMP,PROMPT,LRESID,LSTEP,RESTAR

          FMT='('' NX_CLASS: '',9I3)'
          WRITE(IUNIT,FMT)  (NX_CLASS(NX_LIST(nxx)),nxx=1,NX_LIST(0))
          FMT='('' NX_LOCKS: '',9I3)'
          WRITE(IUNIT,FMT)  (NX_LOCKS(NX_LIST(nxx)),nxx=1,NX_LIST(0))
          FMT='('' NX_TYPE: '',9I3)'
          WRITE(IUNIT,FMT)   (NX_TYPE(NX_LIST(nxx)),nxx=1,NX_LIST(0))
          nrc=2 !temporary AJP
          DO nxx=1,NX_LIST(0)
            nx=NX_LIST(nxx)
            FMT='(/'' NYT(nrc,1,'',I1,'')='',I6,'' NOT(nrc,0,'',I1,'
     '        //''')='',I6)'
            WRITE(IUNIT,FMT) nx,NYT(nrc,1,nx),nx,NOT(nrc,1,0,nx)
          ENDDO !nxx

          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              FMT='(/'' Region nr='',I2)'
              WRITE(IUNIT,FMT) NR
              FMT='('' NCT(nr='',I2,'',nx='',I2,'')='',I6)'
              WRITE(IUNIT,FMT) nr,nx,NCT(nr,nx)
              DO nc=1,NCT(nr,nx)
                FMT='('' NYT(nrc='',I1,''nc='',I2,'',nx=1..)='',9I6)'
                WRITE(IUNIT,FMT) nrc,nc,NYT(nrc,nc,nx)
              ENDDO !nc
              FMT='('' NOT(nrc,nr='',I2,'',nx=1..)='',9I6)'
              WRITE(IUNIT,FMT) nr,NOT(nrc,1,nr,nx)

              IF(ITYP1(nr,nx).EQ.4) THEN !linear elasticity
                DO ie=1,12
                  FMT='(/'' Element type ie='',I3)'
                  WRITE(IUNIT,FMT) ie
                  FMT='('' ETYP(ie)='',L1,'' ILT(ie,nr,nx)='','
     '              //'I3,'' IMT(ie)='','
     '              //'I3,'' NGP(ie)='',I3,'' NLP(ie)='','
     '              //'I3,'' NMP(ie)='',I3,'' NVE(ie)='','
     '              //'I3,'' NHV(nv=1..6,ie):'',6I2)'
                  WRITE(IUNIT,FMT) ETYP(ie),ILT(ie,nr,nx),IMT(ie),
     '              NGP(ie),NLP(ie),NMP(ie),NVE(ie),(NHV(nv,ie),NV=1,6)
                  FMT='('' ILP(il=1..,ie,nr,nx):'',20I2)'
                  WRITE(IUNIT,FMT) (ILP(IL,ie,nr,nx),IL=1,ILT(ie,nr,nx))
                  FMT='('' NMB(il=1..,ie,nx):'',20I2)'
                  WRITE(IUNIT,FMT) (NMB(IL,ie,nx),IL=1,ILT(ie,nr,nx))
                  FMT='('' IT(1..3,ie):'',3I3)'
                  WRITE(IUNIT,FMT) IT(1,ie),IT(2,ie),IT(3,ie)
                ENDDO !ie

              ELSE IF(ITYP1(nr,nx).EQ.3.  !pdes
     '             OR.ITYP1(nr,nx).EQ.5.  !finite elasticity
     '             OR.ITYP1(nr,nx).EQ.9) THEN !BEM
                FMT='('' ETYP(nr)='',L1,'' ILT(1,nr,nx)='','
     '            //'I3,'' IMT(nr)='','
     '           //'I3,'' NGP(nr)='',I3,'' NLP(nr)='',I3,'' NMP(nr)='','
     '            //'I3,'' NVE(nr)='',I3,'' NHV(nv=1..6,nr):'',6I2)'
                WRITE(IUNIT,FMT) ETYP(nr),
     '            ILT(1,nr,nx),IMT(nr),NGP(nr),NLP(nr),NMP(nr),
     '            NVE(nr),(NHV(nv,nr),NV=1,6)
                FMT='('' ILP(il=1..,nr,nx):'',20I2)'
                WRITE(IUNIT,FMT) (ILP(IL,1,nr,nx),IL=1,ILT(1,nr,nx))
                FMT='('' NMB(il=1..,nr,nx):'',20I2)'
                WRITE(IUNIT,FMT) (NMB(IL,nr,nx),IL=1,ILT(1,nr,nx))
                FMT='('' IT(1..3,nr):'',3I3)'
                WRITE(IUNIT,FMT) IT(1,nr),IT(2,nr),IT(3,nr)
              ENDIF
            ENDDO !nx
          ENDDO !nr

          nc=1 ! Temporary MPN 23/11/94
          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                FMT='(/'' Region nr='',I2,''   Node np='',I4)'
                WRITE(IUNIT,FMT) nr,NP
                FMT='('' NHP(np,nr,nx)='',I2)'
                WRITE(IUNIT,FMT) NHP(np,nr,nx)
                FMT='('' NKH(nh=1..,np,nc,nr):'',12I3)'
                WRITE(IUNIT,FMT) (NKH(NH_LOC(nhx,nx),np,nc,nr),
     '            nhx=1,NHP(np,nr,nx))
                FMT='('' NVHP(nh=1..,np,nc,nr):'',12I3)'
                WRITE(IUNIT,FMT) (NVHP(NH_LOC(nhx,nx),np,nc,nr),
     '            nhx=1,NHP(np,nr,nx))
                DO nhx=1,NHP(np,nr,nx)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,nc,nr)
                    DO nk=1,NKH(nh,np,nc,nr)
                      DO nrc=0,2
                        FMT='('' NYNP(nk='',I2,'',nv='',I2,'',nh='',I2,'
     '                  //''',np,nrc='',I1,'',nc='',I1,'',nr='',I1,''):'
     '                  //''',10I6,/(37X,10I6))'
                        WRITE(IUNIT,FMT) nk,nv,nh,nrc,nc,nr,
     '                    NYNP(nk,nv,nh,np,nrc,nc,nr)
                        ny=NYNP(nk,nv,nh,np,nrc,nc,nr)
                        FMT='('' NPNY(0..6,ny='',I5,'',nrc='',I1,'
     '                    //''',nx='',I1,'') = '',7I6)'
                        WRITE(IUNIT,FMT) nyy,NRCC,nx,
     '                    (NPNY(i,ny,nrc,nx),i=0,6)
                      ENDDO !nrc
                    ENDDO
                  ENDDO
                ENDDO
                FMT='('' CP(nm=1..,np,nx):'',/(3E25.17))'
                WRITE(IUNIT,FMT) (CP(nm,np,nx),NM=1,NMM)
              ENDDO !nonode
            ENDDO !nx
          ENDDO !nr

          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                FMT='(/'' Region nr='',I2,''   Element ne='',I4)'
                WRITE(IUNIT,FMT) nr,NE
                FMT='('' NHE(ne,nx)='',I2)'
                WRITE(IUNIT,FMT) NHE(ne,nx)
                DO nhx=1,NHE(ne,nx)
                  nh=NH_LOC(nhx,nx)
                  FMT='('' NBH(nh='',I2,'',nc,ne):'',20I3)'
                  WRITE(IUNIT,FMT) nh,NBH(nh,nc,ne)
                ENDDO !nh

                DO nb=1,NBFT
                  DO nhx=1,NHE(ne,nx)
                    nh=NH_LOC(nhx,nx)
                    FMT='('' NVHE(nn=1..,nb='',I2,'',nh='',I2,'
     '                //''',ne):'',20I2)'
                    WRITE(IUNIT,FMT) nb,nh,
     '                (NVHE(nn,nb,nh,ne),nn=1,NNT(nb))
                  ENDDO !nh
                ENDDO !nb

                FMT='(/'' Solution type nx='',I2,)'
                WRITE(IUNIT,FMT) NX
                DO nhx=1,NHE(ne,nx)
                  nh=NH_LOC(nhx,nx)
                  DO nrc=0,2
                    nb=NBH(nh,nc,ne)
                    DO nn=1,NNT(nb)
                      DO nv=1,NVHE(nn,nb,nh,ne)
                        FMT='('' NYNE(na=1..,nh='',I2,'',nrc='',I1,'
     '                    //'''nc='',I2,'',ne):'',10I6)'
                        WRITE(IUNIT,FMT) NH,nrc,nc,
     '                    (NYNE(na,nh,nrc,nc,ne),na=1,NAT(nb))
                      ENDDO !nv
                    ENDDO !nn
                  ENDDO !nrc
                ENDDO !nh
C news MPN 14Jun2000: added nx index to NW
                FMT='('' NW(ne,1,nx)='',I3)'
                WRITE(IUNIT,FMT) NW(ne,1,nx)
                FMT='('' CE(nm=1..,ne,nx):'',/(3E25.17))'
                WRITE(IUNIT,FMT) (CE(nm,ne,nx),NM=1,NMM)
              ENDDO !noelem
            ENDDO !nx
          ENDDO !nr

          nrc=0 !Temporary AJP 1-12-94
          DO nx1=1,NX_LIST(0)
            nxx=NX_LIST(nx1)
            DO nr=1,NRT
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO nhx=1,NHP(np,nr,nxx)
                  nh=NH_LOC(nhx,nxx)
                  DO nv=1,NVHP(nh,np,nc,nr)
                    DO nk=1,NKH(nh,np,nc,nr)
                      ny=NYNP(nk,nv,nh,np,nrc,nc,nr)
                      FMT='(/'' FIX('',I5,'',iy=1..,nx='',I1,''): '''
     '                  //',5L1)'
                      WRITE(IUNIT,FMT) ny,nxx,
     '                  (FIX(ny,IY,nxx),IY=1,NIYFIXM)
                      FMT='('' YP('',I5,'',iy=1..,nx='',I1,''): '',/,'
     '                  //'(3E25.17))'
                      WRITE(IUNIT,FMT) ny,nxx,(YP(ny,IY,nxx),IY=1,NIYM)
C!!! MPN needs fixing
                      nrc=1 !temporary
                      FMT='('' NONY(0,'',I5,'',nr='',I1,'',nx='','
     '                  //'I1,''):'',I3)'
                      WRITE(IUNIT,FMT) ny,nr,nxx,NONY(0,ny,nrc,nr,nxx)
                      FMT='('' NONY(noy=1..,'',I5,'',nr='','
     '                  //'I1,'',nx='',I1,''):'',9I3)'
                      WRITE(IUNIT,FMT) ny,nr,nxx,
     '                  (NONY(NOY,ny,nrc,nr,nxx),
     '                  NOY=1,NONY(0,ny,nrc,nr,nxx))
                      FMT='('' CONY(noy=1..,'',I5,'',nr='','
     '                  //'I1,'',nx='',I1,''):'',(3E25.17))'
                      WRITE(IUNIT,FMT) ny,nr,nxx,
     '                  (CONY(NOY,ny,nrc,nr,nxx),
     '                  NOY=1,NONY(0,ny,nrc,nr,nxx))
C CPB 13/11/94 Temporarily fixed NYNO will need to be redone when
C NYNO is calculated and used properly
                      FMT='('' NYNO(NONY(1..,ny,nr,nxx),nr='',I1,'
     '                  //''',nx='',I1,'') = '',10I5)'
                      WRITE(IUNIT,FMT) nr,nx,
     '                  (NYNO(1,NONY(NOY,ny,nrc,nr,nxx),nrc,nr,nx),
     '                  NOY=1,NONY(0,ny,nrc,nr,nxx))
                    ENDDO !nk
                  ENDDO !nv
                ENDDO !nh
              ENDDO !nonode

              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                DO nhx=1,NHE(ne,nxx)
                  nh=NH_LOC(nhx,nxx)
                  nb=NBH(nh,nc,ne)
                  DO nn=1,NNT(nb)
                    DO nv=1,NVHE(nn,nb,nh,ne)
                      DO na=1,NAT(nb)
                  ny=NYNE(na,nh,nrc,nc,ne)
                  FMT='(/'' FIX('',I5,'',iy=1..,nx='',I1,''): '','
     '                    //'5L1)'
                        WRITE(IUNIT,FMT) ny,nxx,
     '                    (FIX(ny,IY,nxx),IY=1,NIYFIXM)
                        FMT='('' YP('',I5,'',iy=1..,nx='',I1,''): '',/,'
     '                    //'(3E25.17))'
                  WRITE(IUNIT,FMT) ny,nxx,
     '                    (YP(ny,IY,nxx),IY=1,NIYM)
C!!! MPN needs fixing
                        nrc=1 !temporary
                        FMT='('' NONY(0,'',I5,'',nr='',I1,'',nx='','
     '                    //'I1,''):'',I3)'
                  WRITE(IUNIT,FMT) ny,nr,nxx,NONY(0,ny,nrc,nr,nxx)
                  FMT='('' NONY(noy=1..,'',I5,'',nr='','
     '                    //'I1,'',nx='',I1,''):'',9I3)'
                  WRITE(IUNIT,FMT) ny,nr,nxx,
     '              (NONY(NOY,ny,nrc,nr,nxx),
     '                    NOY=1,NONY(0,ny,nrc,nr,nxx))
                        FMT='('' CONY(noy=1..,'',I5,'',nr='','
     '                    //'I1,'',nx='',I1,''):'',(3E25.17))'
                  WRITE(IUNIT,FMT) ny,nr,nxx,
     '              (CONY(NOY,ny,nrc,nr,nxx),
     '                    NOY=1,NONY(0,ny,nrc,nr,nxx))
                      ENDDO !na
                    ENDDO !nv
                  ENDDO !nn
                ENDDO !nh
              ENDDO !noelem
            ENDDO !nr
          ENDDO !nx1

          DO nr=1,NRT
            IF(KTYP57(nr).EQ.2.OR.KTYP57(nr).EQ.4.OR.KTYP57(nr).EQ.5)
     '        THEN !Press incr
              FMT='(/'' Region '',I1,'':Pressure Boundary Conditions'')'
              WRITE(IUNIT,FMT) NR
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                FMT='(''   FIXP(iface=1..2,ne='',I3,''): '',2L1)'
                WRITE(IUNIT,FMT) ne,(FIXP(iface,ne),iface=1,2)
                FMT='(''   PE(iface=1..2,ne='',I3,''): '',2E25.17)'
                WRITE(IUNIT,FMT) ne,(PE(iface,ne),iface=1,2)
                FMT='(''   PF(iface=1..2,ne='',I3,''): '',2E25.17)'
                WRITE(IUNIT,FMT) ne,(PF(iface,ne),iface=1,2)
              ENDDO
            ENDIF
          ENDDO

          FMT='(/'' FEXT(i=1..,ng,ne):'')'
          WRITE(IUNIT,FMT)
          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(NBH(NH_LOC(1,nx),nc,ne).GT.0) THEN
                  DO ng=1,NGT(NBH(NH_LOC(1,nx),nc,ne))
                    FMT='('' FEXT(i,'',I3,'','',I4,''): '',3E25.17,'
     '                //'/(4E25.17))'
                    WRITE(IUNIT,FMT) ng,ne,(FEXT(i,ng,ne),i=1,NIFEXTM)
                  ENDDO
                ENDIF
              ENDDO
              IF(KTYP53(nr).EQ.3)THEN !Active stress included
                FMT='(/'' CALL_ACTI='',L1)'
                WRITE(IUNIT,FMT) CALL_ACTI
              ENDIF
            ENDDO !nx
          ENDDO !nr

          IF(KTYP1.EQ.8.AND.KTYP3.EQ.2) THEN !activation variables
            DO nr=1,NRT
              DO nxx=1,NX_LIST(0)
                nx=NX_LIST(nxx)
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  FMT='(/'' Element ne='',I4)'
                  WRITE(IUNIT,FMT) NE
                  nb=NBH(NH_LOC(1,nx),nc,ne)
                  IF(nb.GT.0) THEN
                    FMT='(/'' ITHRES(1,ng,ne):'',/(1X,75I1))'
                    WRITE(IUNIT,FMT) (ITHRES(1,ng,ne),ng=1,NGT(nb))
                    FMT='(/'' ITHRES(2,ng,ne):'',/(1X,75I1))'
                    WRITE(IUNIT,FMT) (ITHRES(2,ng,ne),ng=1,NGT(nb))
                    FMT='(/'' THRES(1,ng,ne):'',/(1X,10E12.4))'
                    WRITE(IUNIT,FMT) (THRES(1,ng,ne),ng=1,NGT(nb))
                    FMT='(/'' THRES(2,ng,ne):'',/(1X,10E12.4))'
                    WRITE(IUNIT,FMT) (THRES(2,ng,ne),ng=1,NGT(nb))
                    FMT='(/'' THRES(3,ng,ne):'',/(1X,10E12.4))'
                    WRITE(IUNIT,FMT) (THRES(3,ng,ne),ng=1,NGT(nb))
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF

          FMT='(/'' CALL_EQUA='',L1,'' CALL_INIT='',L1,'
     '      //'  '' CALL_MATE='',L1,'' CALL_SOLV='',L1)'
          WRITE(IUNIT,FMT) CALL_EQUA,CALL_INIT,CALL_MATE,CALL_SOLV
          FMT='('' CALL_DATA='',L1,'' CALL_FIT ='',L1,'
     '      //' '' CALL_GROW='',L1,'' CALL_MOTI='',L1)'
          WRITE(IUNIT,FMT) CALL_DATA,CALL_FIT ,CALL_GROW,CALL_MOTI
          IF(KTYP90.GT.0)THEN
            FMT='('' CALL_COUP='',L1)'
            WRITE(IUNIT,FMT) CALL_COUP
          ENDIF
          DO nr=1,NRT
            DO nxx=1,NX_LIST(0)
              nx=NX_LIST(nxx)
              IF(CALL_SOLV.AND.ITYP4(nr,nx).EQ.2)THEN !BE
                FMT='(/'' ADAPINT='',L1,'' COMPLEX='',L1,'
     '            //'  '' HERMITE='',L1,'' HYP='',L1)'
                WRITE(IUNIT,FMT) ADAPINT,COMPLEX,HERMITE,HYP
              ENDIF
            ENDDO
          ENDDO

        ENDIF !dependent

      ELSE
        ERROR=' Command error: COMAND='//COMAND
        GOTO 9999
      ENDIF

      CALL EXITS('IOGEOM')
      RETURN
 9999 CALL ERRORS('IOGEOM',ERROR)
      CALL EXITS('IOGEOM')
      RETURN 1
      END


