      SUBROUTINE DOMSOLOPTI3D(IBT,IDO,INP,NBH,NBJ,NDET,NEELEM,NGAP,NHE,
     '  NKH,NKHE,NKJE,NLL,NPF,NP_INTERFACE,NPNE,nr,NVHE,NVJE,NW,nx,CE,
     '  CPCUTOFF,CURVCORRECT,DETLOC,DET_ADAPT,DL,DRDN,PG,PG_J,PG_Q,PG_U,
     '  RAD,RD,RG,SE,WG,XE,XG1,XIG,XIG_J,XIG_Q,XIG_U,XN,XP,XPFP,XR,YD,
     '  ZA,ZE,ZF,ZP,ERROR,*)

C#### Subroutine: DOMSOLOPTI3D
C###  Description:
C###    DOMSOLOPTI3D calculates BE solution at the point XPFP but
C###    is based on the optimised 3d laplace code. This version
C###    includes adaptive integration and element splitting
C**** Created by Martin Buist - June 2001

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NDET(NBFM,0:NNM),
     '  NEELEM(0:NE_R_M,0:NRM),NGAP(NIM,NBM),NHE(NEM),NKH(NHM,NPM,NCM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),nr,
     '  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx
      REAL*8 CE(NMM,NEM),CPCUTOFF,CURVCORRECT(2,2,NNM,NEM),
     '  DETLOC(NBFM,0:NNM,NGM,6),DET_ADAPT(NBFM,0:NNM,NGM),DL(3,NLM),
     '  DRDN(NGM),PG(NSM,NUM,NGM,NBM),PG_J(NSM,NUM,NGM),
     '  PG_Q(NSM,NUM,NGM),PG_U(NSM,NUM,NGM),RAD(NGM),RD(NGM),RG(NGM),
     '  SE(NSM,NBFM,NEM),WG(NGM,NBM),XE(NSM,NJM),
     '  XG1(NJM,NUM,NGM),XIG(NIM,NGM,NBM),XIG_J(NIM,NGM),XIG_Q(NIM,NGM),
     '  XIG_U(NIM,NGM),XN(NJM,NGM),XP(NKM,NVM,NJM,NPM),XPFP(*),
     '  XR(NJM,NGM),YD,ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),
     '  ZF(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER intscheme,ISP,mk,nb,nb2,nbuh,nbuhp,nb1j,nb1jp,nbbem,nbqh,
     '  nbqhp,ne,ng,nh,nhx,nk1,nn,NNMIN,nnsp,noelem,np,npnn,ns,nsplit
      REAL*8 D(2),DDOT,DGREEN,DXXI(3,3),FACTOR,FACTOR1,GREEN,MINDIST,
     '  SUMC,SUMU,SUMQ,VLENGTH,XIMIN(2)
      LOGICAL INTERFACE,USEADAPINT

      CALL ENTERS('DOMSOLOPTI3D',*9999)

      YD=0.0d0
      SUMC=0.0d0

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        np=NPNE(1,NBJ(1,ne),ne)
        INTERFACE=(NP_INTERFACE(np,0).GT.1).AND.
     '    (NP_INTERFACE(np,1).NE.nr)

C***    Transfer global parameters to local parameters XPXE
        nb=NBJ(1,ne)
        ns=0
        DO nn=1,NNT(nb)
          np=NPNE(nn,nb,ne)
          DO mk=1,NKT(nn,nb)
            ns=ns+1
            XE(ns,1)=XP(NKJE(mk,nn,1,ne),NVJE(nn,nb,1,ne),1,np)*
     '        SE(ns,nb,ne)
            XE(ns,2)=XP(NKJE(mk,nn,2,ne),NVJE(nn,nb,2,ne),2,np)*
     '        SE(ns,nb,ne)
            XE(ns,3)=XP(NKJE(mk,nn,3,ne),NVJE(nn,nb,3,ne),3,np)*
     '        SE(ns,nb,ne)
          ENDDO
        ENDDO
C***    Transfer global potential vector to element solution vector
        CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,2),nx,
     '    CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '    ERROR,*9999)
C***    Transfer global current vector to element solution vector
        CALL ZPZE(NBH(1,1,ne),2,NHE(ne),NKHE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,2),nx,
     '    CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZF,ZP,
     '    ERROR,*9999)
        !ZE(ns,nhx) contains the dependent variable values
        !ZF(ns,nh) contains the normal derivative values

C***    Determine integration scheme to use for element ne
        CALL DIST_LOOSE(intscheme,NBJ(1,ne),NLL(1,ne),NNMIN,
     '    NPNE(1,1,ne),nr,NW(ne,2),DL,MINDIST,XP,XPFP,ERROR,*9999)

        IF(intscheme.NE.1) THEN
          NNMIN=0
        ENDIF

C***    Decide upon the appropriate basis for the integration scheme
        CALL QUADBE(intscheme,NBJ(1,ne),nbbem,ERROR,*9999)
        nb1jp=NFBASE(1,NBASEF(NBJ(1,ne),nbbem))

C***    Setup adaptive integration information if nescessary.
        USEADAPINT=.FALSE.
        IF(intscheme.EQ.2.AND.NW(ne,2).NE.11.AND.NW(ne,2).NE.13) THEN
C***      Don't use this on distorted cubic linear elements.
C***      Find the minimum distance in the element to the singular point
C***      for Telles' rule.
          USEADAPINT=.TRUE.
          CALL TELLES(IBT,IDO,INP,20,NBJ(1,ne),NLL(1,ne),NW(ne,2),D,DL,
     '      1.0d-6,XE,XIMIN,XPFP,ERROR,*9999)
C***      Set up adaptive basis functions based on Telles' rule.
          nb1j=NBASEF(NBJ(1,ne),nbbem)
          nbuh=NBASEF(NBH(NH_LOC(1,nx),1,ne),nbbem)
          nbqh=NBASEF(NBH(NH_LOC(1,nx),2,ne),nbbem)
          CALL GAUSS11(IBT,IDO,INP,nb1j,nb1jp,nbqh,nbuh,NGAP,D,DETLOC,
     '      DET_ADAPT,PG_J,PG_Q,PG_U,XIG,XIG_J,XIG_Q,XIG_U,XIMIN,
     '      ERROR,*9999)
        ENDIF

C***    Calculate the element integrals for the current node
C***    and dependent variable.
        ISP=0
        IF(intscheme.EQ.1) THEN !Element splitting
C***      Find number of split elements used at nodes < NNSPLIT(nonode)
          DO nnsp=0,NNMIN-1
            ISP=ISP+NDET(nb1jp,nnsp)
          ENDDO
        ENDIF

        nhx=1
        nh=NH_LOC(nhx,nx)
C***    Loop over number of split elements (if any)
        DO nsplit=1,NDET(nb1jp,NNMIN)
          nb2=nbbem+ISP
C***      Decide basis functions from the parent family basis function
C***      for dependent (u) and normal derivative (q).
          nbuhp=NFBASE(1,NBASEF(NBH(nh,1,ne),nbbem))
          nbqhp=NFBASE(1,NBASEF(NBH(nh,2,ne),nbbem))

C***      nbuh,nbqh and nb1j are the global basis functions with the
C***      appropriate quadrature for the dependent (H), normal
C***      derivative (Q) and geometric (J) variables.
          nbuh=NBASEF(NBH(nh,1,ne),nb2)
          nbqh=NBASEF(NBH(nh,2,ne),nb2)
          nb1j=NBASEF(NBJ(1,ne),nb2)
C***      set up the gauss point dependent arrays.
          DO ng=1,NGT(nb1j)
            IF(USEADAPINT) THEN
              XG1(1,1,ng)=DDOT(NST(nb1jp),PG_J(1,1,ng),1,XE(1,1),1)
              XG1(2,1,ng)=DDOT(NST(nb1jp),PG_J(1,1,ng),1,XE(1,2),1)
              XG1(3,1,ng)=DDOT(NST(nb1jp),PG_J(1,1,ng),1,XE(1,3),1)
              XG1(1,2,ng)=DDOT(NST(nb1jp),PG_J(1,2,ng),1,XE(1,1),1)
              XG1(2,2,ng)=DDOT(NST(nb1jp),PG_J(1,2,ng),1,XE(1,2),1)
              XG1(3,2,ng)=DDOT(NST(nb1jp),PG_J(1,2,ng),1,XE(1,3),1)
              XG1(1,4,ng)=DDOT(NST(nb1jp),PG_J(1,4,ng),1,XE(1,1),1)
              XG1(2,4,ng)=DDOT(NST(nb1jp),PG_J(1,4,ng),1,XE(1,2),1)
              XG1(3,4,ng)=DDOT(NST(nb1jp),PG_J(1,4,ng),1,XE(1,3),1)
            ELSE
              XG1(1,1,ng)=DDOT(NST(nb1jp),PG(1,1,ng,nb1j),1,XE(1,1),1)
              XG1(2,1,ng)=DDOT(NST(nb1jp),PG(1,1,ng,nb1j),1,XE(1,2),1)
              XG1(3,1,ng)=DDOT(NST(nb1jp),PG(1,1,ng,nb1j),1,XE(1,3),1)
              XG1(1,2,ng)=DDOT(NST(nb1jp),PG(1,2,ng,nb1j),1,XE(1,1),1)
              XG1(2,2,ng)=DDOT(NST(nb1jp),PG(1,2,ng,nb1j),1,XE(1,2),1)
              XG1(3,2,ng)=DDOT(NST(nb1jp),PG(1,2,ng,nb1j),1,XE(1,3),1)
              XG1(1,4,ng)=DDOT(NST(nb1jp),PG(1,4,ng,nb1j),1,XE(1,1),1)
              XG1(2,4,ng)=DDOT(NST(nb1jp),PG(1,4,ng,nb1j),1,XE(1,2),1)
              XG1(3,4,ng)=DDOT(NST(nb1jp),PG(1,4,ng,nb1j),1,XE(1,3),1)
            ENDIF
C***        RG(ng) stores the Jacobian for a length or area integral
C***        at Gauss point ng.
C***        DXXI(nj,ni) contains the values of dXj/dPSIi at the
C***        Gauss point.
            DXXI(1,1)=XG1(1,2,ng)
            DXXI(2,1)=XG1(2,2,ng)
            DXXI(3,1)=XG1(3,2,ng)
            DXXI(1,2)=XG1(1,4,ng)
            DXXI(2,2)=XG1(2,4,ng)
            DXXI(3,2)=XG1(3,4,ng)
            RG(ng)=DSQRT(
     '        (DXXI(2,1)*DXXI(3,2)-DXXI(2,2)*DXXI(3,1))*
     '        (DXXI(2,1)*DXXI(3,2)-DXXI(2,2)*DXXI(3,1))+
     '        (DXXI(1,2)*DXXI(3,1)-DXXI(1,1)*DXXI(3,2))*
     '        (DXXI(1,2)*DXXI(3,1)-DXXI(1,1)*DXXI(3,2))+
     '        (DXXI(1,1)*DXXI(2,2)-DXXI(1,2)*DXXI(2,1))*
     '        (DXXI(1,1)*DXXI(2,2)-DXXI(1,2)*DXXI(2,1)))
C***        Find the unit outward normal.
            XN(1,ng)=XG1(2,2,ng)*XG1(3,4,ng)-XG1(3,2,ng)*XG1(2,4,ng)
            XN(2,ng)=XG1(3,2,ng)*XG1(1,4,ng)-XG1(1,2,ng)*XG1(3,4,ng)
            XN(3,ng)=XG1(1,2,ng)*XG1(2,4,ng)-XG1(2,2,ng)*XG1(1,4,ng)
            VLENGTH=DSQRT(XN(1,ng)*XN(1,ng)+XN(2,ng)*XN(2,ng)+
     '        XN(3,ng)*XN(3,ng))
            CALL ASSERT(VLENGTH.GT.RDELTA,'>>Zero length normal',
     '        ERROR,*9999)
C           IF((INTERFACE.AND.(NW(ne,3).NE.1)).OR.(NW(ne,3).EQ.1)) THEN
            IF(INTERFACE) THEN
              XN(1,ng)=-XN(1,ng)/VLENGTH
              XN(2,ng)=-XN(2,ng)/VLENGTH
              XN(3,ng)=-XN(3,ng)/VLENGTH
            ELSE
              XN(1,ng)=XN(1,ng)/VLENGTH
              XN(2,ng)=XN(2,ng)/VLENGTH
              XN(3,ng)=XN(3,ng)/VLENGTH
            ENDIF
          ENDDO !End ng loop

C***      XEP bit
          IF(USEADAPINT) THEN
C            IF(IGREN(nr).EQ.2) CE(1,ne)=1.0d0
            FACTOR=4.0d0*PI
            FACTOR1=FACTOR*CE(1,ne)

C***        Calculate gauss point dependent arrays
            DO ng=1,NGT(nb1j)
              XR(1,ng)=XG1(1,1,ng)-XPFP(1)
              XR(2,ng)=XG1(2,1,ng)-XPFP(2)
              XR(3,ng)=XG1(3,1,ng)-XPFP(3)
              RAD(ng)=DSQRT(XR(1,ng)*XR(1,ng)+XR(2,ng)*XR(2,ng)+
     '          XR(3,ng)*XR(3,ng))
              DRDN(ng)=(XR(1,ng)*XN(1,ng)+XR(2,ng)*XN(2,ng)+
     '          XR(3,ng)*XN(3,ng))/RAD(ng)
            ENDDO !ng

C***        GK matrix components
            ns=0
            DO ng=1,NGT(nbuh)
              DGREEN=-DET_ADAPT(nbuhp,NNMIN,ng)*DRDN(ng)/
     '          (RAD(ng)*RAD(ng))
              RD(ng)=DGREEN*RG(ng)*WG(ng,nbuh)
              SUMC=SUMC-RD(ng)/FACTOR
            ENDDO !ng
            DO nn=1,NNT(nbuhp) !Dependent variable loop
              npnn=NPNE(nn,nbuhp,ne)
              DO nk1=1,NKT(nn,nbuhp)
                MK=NKHE(nk1,nn,nh,ne) !MK is zero at a distorted node
                ns=ns+1
                IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,1)-
     '            KTYP93(1,nr),1)) THEN
                  SUMU=0.0d0
                  DO ng=1,NGT(nbuh)
                    SUMU=SUMU+RD(ng)*PG_U(ns,1,ng)
                  ENDDO !End of ng loop
                  YD=YD-SUMU*ZE(ns,nhx)/FACTOR
                ENDIF !End of MK > 0 loop
              ENDDO !End of nk1 loop
            ENDDO !End of nn loop

C***        GQ matrix components
            ns=0
            DO ng=1,NGT(nbqh)
              GREEN=DET_ADAPT(nbuhp,NNMIN,ng)/RAD(ng)
              RD(ng)=GREEN*RG(ng)*WG(ng,nbqh)
            ENDDO !ng
            DO nn=1,NNT(nbqhp) !normal derivative loop
              npnn=NPNE(nn,nbqhp,ne)
              DO nk1=1,NKT(nn,nbqhp)
                MK=NKHE(nk1,nn,nh,ne) !MK is zero at a distorted node
                ns=ns+1
                IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,2)-
     '            KTYP93(2,nr),1)) THEN
                  SUMQ=0.0d0
                  DO ng=1,NGT(nbqh)
                    SUMQ=SUMQ+RD(ng)*PG_Q(ns,1,ng)
                  ENDDO !End ng loop
                  YD=YD+SUMQ*ZF(ns,nh)/FACTOR1
                ENDIF !End of MK > 0 loop
              ENDDO !End of nk1 loop
            ENDDO !End of nn loop

          ELSE !not adaptive int

C            IF(IGREN(nr).EQ.2) CE(1,ne)=1.0d0
            FACTOR=4.0d0*PI
            FACTOR1=FACTOR*CE(1,ne)

C***        Calculate gauss point dependent arrays
            DO ng=1,NGT(nb1j)
              XR(1,ng)=XG1(1,1,ng)-XPFP(1)
              XR(2,ng)=XG1(2,1,ng)-XPFP(2)
              XR(3,ng)=XG1(3,1,ng)-XPFP(3)
              RAD(ng)=DSQRT(XR(1,ng)*XR(1,ng)+XR(2,ng)*XR(2,ng)+
     '          XR(3,ng)*XR(3,ng))
              DRDN(ng)=(XR(1,ng)*XN(1,ng)+XR(2,ng)*XN(2,ng)+
     '          XR(3,ng)*XN(3,ng))/RAD(ng)
            ENDDO !ng

C***        GK matrix components
            ns=0
            DO ng=1,NGT(nbuh)
              DGREEN=-DETLOC(nbuhp,NNMIN,ng,nsplit)*DRDN(ng)/
     '          (RAD(ng)*RAD(ng))
              RD(ng)=DGREEN*RG(ng)*WG(ng,nbuh)
              SUMC=SUMC-RD(ng)/FACTOR
            ENDDO !ng
            DO nn=1,NNT(nbuhp) !Dependent variable loop
              npnn=NPNE(nn,nbuhp,ne)
              DO nk1=1,NKT(nn,nbuhp)
                MK=NKHE(nk1,nn,nh,ne) !MK is zero at a distorted node
                ns=ns+1
                IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,1)-
     '            KTYP93(1,nr),1)) THEN
                  SUMU=0.0d0
                  DO ng=1,NGT(nbuh)
                    SUMU=SUMU+RD(ng)*PG(ns,1,ng,nbuh)
                  ENDDO !End of ng loop
                  YD=YD-SUMU*ZE(ns,nhx)/FACTOR
                ENDIF !End of MK > 0 loop
              ENDDO !End of nk1 loop
            ENDDO !End of nn loop

C***        GQ matrix components
            ns=0
            DO ng=1,NGT(nbqh)
              GREEN=DETLOC(nbuhp,NNMIN,ng,nsplit)/RAD(ng)
              RD(ng)=GREEN*RG(ng)*WG(ng,nbqh)
            ENDDO !ng
            DO nn=1,NNT(nbqhp) !normal derivative loop
              npnn=NPNE(nn,nbqhp,ne)
              DO nk1=1,NKT(nn,nbqhp)
                MK=NKHE(nk1,nn,nh,ne) !MK is zero at a distorted node
                ns=ns+1
                IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,2)-
     '            KTYP93(2,nr),1)) THEN
                  SUMQ=0.0d0
                  DO ng=1,NGT(nbqh)
                    SUMQ=SUMQ+RD(ng)*PG(ns,1,ng,nbqh)
                  ENDDO !End ng loop
                  YD=YD+SUMQ*ZF(ns,nh)/FACTOR1
                ENDIF !End of MK > 0 loop
              ENDDO !End of nk1 loop
            ENDDO !End of nn loop
          ENDIF
          ISP=ISP+nsplit
        ENDDO !End of nsplit loop
      ENDDO !ne

      IF(DABS(SUMC).LT.CPCUTOFF) THEN
        YD=0.0d0
      ELSE
        YD=YD/SUMC
      ENDIF

      CALL EXITS('DOMSOLOPTI3D')
      RETURN
 9999 CALL ERRORS('DOMSOLOPTI3D',ERROR)
      CALL EXITS('DOMSOLOPTI3D')
      RETURN 1
      END


