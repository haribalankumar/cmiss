      SUBROUTINE GRADDOMSOLOPTI3D(NBH,NBJ,NEELEM,NHE,NKH,NKHE,NKJE,
     '  NPF,NP_INTERFACE,NPNE,nr,NVHE,NVJE,NW,nx,
     '  CE,CURVCORRECT,DET,DRDN,GRADPHI,PG,RAD,RG,SE,
     '  WG,XE,XG1,XN,XP,XPFP,XR,ZA,ZE,ZF,ZP,ERROR,*)

C#### Subroutine: GRADDOMSOLOPTI3D
C###  Description:
C###    GRADDOMSOLOPTI3D calculates the gradient of the potential at
C###    the internal point XPFP(nj).
C**** Created by Martin Buist - June 2001

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NHE(NEM),NKH(NHM,NPM,NCM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  nr,NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx
      REAL*8 CE(NMM,NEM),CURVCORRECT(2,2,NNM,NEM),DET(NBFM,0:NNM,NGM,6),
     '  DRDN(NGM),GRADPHI(*),PG(NSM,NUM,NGM,NBM),RAD(NGM),
     '  RG(NGM),SE(NSM,NBFM,NEM),WG(NGM,NBM),
     '  XE(NSM,NJM),XG1(NJM,NUM,NGM),XN(NJM,NGM),XP(NKM,NVM,NJM,NPM),
     '  XPFP(*),XR(NJM,NGM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZF(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER intscheme,MK,nb,nbbt,nbuh,nbuhp,nb1j,nb1jp,nbbem,nbqh,
     '  nbqhp,ne,ng,nh,nhx,nj,nj_grad,nk1,nn,noelem,np,
     '  npnn,ns
      REAL*8 DDOT,DGREEN,DISTANCE,DXXI(3,3),HYPGREEN,MINDIST,SUM,SUMC,
     '  SUMQ,SUMU,VLENGTH,XNO(3),XR_LOCAL(3)
      LOGICAL INTERFACE

      CALL ENTERS('GRADDOMSOLOPTI3D',*9999)

      nhx=1
      nh=NH_LOC(nhx,nx)
      GRADPHI(1)=0.0d0
      GRADPHI(2)=0.0d0
      GRADPHI(3)=0.0d0
      SUMC=0.0d0

C      CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NVHP,nx,NYNE,NYNP,
C     '  YP,ZA,ZP,ERROR,*9999)

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        np=NPNE(1,NBJ(1,ne),ne)
        INTERFACE=(NP_INTERFACE(np,0).GT.1).AND.
     '    (NP_INTERFACE(np,1).NE.nr)
C        IF(IGREN(nr).EQ.1) CE(1,ne)=1.0d0
C        IF(IGREN(nr).EQ.2) CE(1,ne)=1.0d0

C***    Transfer global parameters to local parameters
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

C***    Determine the adaptive integration scheme for ne
C        DUMMY=ADAPINT
C        ADAPINT=.FALSE.
C        CALL DIST(intscheme,NBJ(1,ne),NLL(1,ne),nnmin,NPNE(1,1,ne),nr,
C     '    NW(ne,2),DL,MINDIST,XP,XPFP,ERROR,*9999)
C        ADAPINT=DUMMY

        MINDIST=RMAX
        DO nn=1,NNT(nb)
          SUM=0.0d0
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            XR_LOCAL(nj)=XP(1,1,nj,NPNE(nn,nb,ne))-XPFP(nj)
            SUM=SUM+XR_LOCAL(nj)*XR_LOCAL(nj)
          ENDDO !nj
          DISTANCE=DSQRT(SUM)
          IF(DISTANCE.LT.MINDIST) THEN
            MINDIST=DISTANCE
          ENDIF
        ENDDO !nn
        nbbt=NBASEF(nb,0)
        IF(MINDIST.LE.DLIM(nbbt-3,NW(ne,2))) THEN
          intscheme=3 !high order gauss scheme
        ELSE IF(MINDIST.LE.DLIM(nbbt-2,NW(ne,2))) THEN !Medium scheme
          intscheme=4 !medium order gauss scheme
        ELSE
          intscheme=5 !low order gauss scheme
        ENDIF

C***    No splitting or adaptive int. implemented here
C        IF(intscheme.EQ.1) THEN
C          intscheme=3
C        ELSEIF(intscheme.EQ.2) THEN
C          intscheme=3
C        ENDIF

C***    Decide upon the appropriate basis for the integration scheme
        CALL QUADBE(intscheme,NBJ(1,ne),nbbem,ERROR,*9999)
        nb1jp=NFBASE(1,NBASEF(NBJ(1,ne),nbbem))
        nbuhp=NFBASE(1,NBASEF(NBH(nh,1,ne),nbbem))
        nbqhp=NFBASE(1,NBASEF(NBH(nh,2,ne),nbbem))
C***    nbuh,nbqh and nb1j are the global
C***    basis functions with the appropriate quadrature for the
C***    dependent (H), normal derivative (Q) and geometric (J)
C***    variables.
        nbuh=NBASEF(NBH(nh,1,ne),nbbem)
        nbqh=NBASEF(NBH(nh,2,ne),nbbem)
        nb1j=NBASEF(NBJ(1,ne),nbbem)

        DO ng=1,NGT(nb1j)
          XG1(1,1,ng)=DDOT(NST(nb1jp),PG(1,1,ng,nb1j),1,XE(1,1),1)
          XG1(2,1,ng)=DDOT(NST(nb1jp),PG(1,1,ng,nb1j),1,XE(1,2),1)
          XG1(3,1,ng)=DDOT(NST(nb1jp),PG(1,1,ng,nb1j),1,XE(1,3),1)
          XG1(1,2,ng)=DDOT(NST(nb1jp),PG(1,2,ng,nb1j),1,XE(1,1),1)
          XG1(2,2,ng)=DDOT(NST(nb1jp),PG(1,2,ng,nb1j),1,XE(1,2),1)
          XG1(3,2,ng)=DDOT(NST(nb1jp),PG(1,2,ng,nb1j),1,XE(1,3),1)
          XG1(1,4,ng)=DDOT(NST(nb1jp),PG(1,4,ng,nb1j),1,XE(1,1),1)
          XG1(2,4,ng)=DDOT(NST(nb1jp),PG(1,4,ng,nb1j),1,XE(1,2),1)
          XG1(3,4,ng)=DDOT(NST(nb1jp),PG(1,4,ng,nb1j),1,XE(1,3),1)
C***      RG(ng) stores the Jacobian for a length or area
C***      integral at Gauss point ng.
C***      DXXI(nj,ni) contains the values of dXj/dPSIi at
C***      the Gauss point.
          DXXI(1,1)=XG1(1,2,ng)
          DXXI(2,1)=XG1(2,2,ng)
          DXXI(3,1)=XG1(3,2,ng)
          DXXI(1,2)=XG1(1,4,ng)
          DXXI(2,2)=XG1(2,4,ng)
          DXXI(3,2)=XG1(3,4,ng)
          RG(ng)=DSQRT(
     '      (DXXI(2,1)*DXXI(3,2)-DXXI(2,2)*DXXI(3,1))*
     '      (DXXI(2,1)*DXXI(3,2)-DXXI(2,2)*DXXI(3,1))+
     '      (DXXI(1,2)*DXXI(3,1)-DXXI(1,1)*DXXI(3,2))*
     '      (DXXI(1,2)*DXXI(3,1)-DXXI(1,1)*DXXI(3,2))+
     '      (DXXI(1,1)*DXXI(2,2)-DXXI(1,2)*DXXI(2,1))*
     '      (DXXI(1,1)*DXXI(2,2)-DXXI(1,2)*DXXI(2,1)))
C***      Find the unit outward normal.
          XN(1,ng)=XG1(2,2,ng)*XG1(3,4,ng)-XG1(3,2,ng)*XG1(2,4,ng)
          XN(2,ng)=XG1(3,2,ng)*XG1(1,4,ng)-XG1(1,2,ng)*XG1(3,4,ng)
          XN(3,ng)=XG1(1,2,ng)*XG1(2,4,ng)-XG1(2,2,ng)*XG1(1,4,ng)
          VLENGTH=DSQRT(XN(1,ng)*XN(1,ng)+XN(2,ng)*XN(2,ng)+
     '      XN(3,ng)*XN(3,ng))
C          CALL ASSERT(VLENGTH.GT.RDELTA,'>>Zero length normal',
C     '      ERROR,*9999)
          IF((INTERFACE.AND.(NW(ne,3).NE.1)).OR.(NW(ne,3).EQ.1)) THEN
            XN(1,ng)=-XN(1,ng)/VLENGTH
            XN(2,ng)=-XN(2,ng)/VLENGTH
            XN(3,ng)=-XN(3,ng)/VLENGTH
          ELSE
            XN(1,ng)=XN(1,ng)/VLENGTH
            XN(2,ng)=XN(2,ng)/VLENGTH
            XN(3,ng)=XN(3,ng)/VLENGTH
          ENDIF
          XR(1,ng)=XG1(1,1,ng)-XPFP(1)
          XR(2,ng)=XG1(2,1,ng)-XPFP(2)
          XR(3,ng)=XG1(3,1,ng)-XPFP(3)

C***      Extra calc's to get sumc
          RAD(ng)=DSQRT(XR(1,ng)*XR(1,ng)+XR(2,ng)*XR(2,ng)+
     '      XR(3,ng)*XR(3,ng))
          DRDN(ng)=(XR(1,ng)*XN(1,ng)+XR(2,ng)*XN(2,ng)+
     '      XR(3,ng)*XN(3,ng))/RAD(ng)
          SUMC=SUMC+((DET(nbuhp,0,ng,1)*DRDN(ng)*RG(ng)*
     '      WG(ng,nbuh))/(RAD(ng)*RAD(ng)*4.0d0*PI))
        ENDDO !End of ng loop for variables not depending on which
              !direction the equation was differentiated.

        DO nj_grad=1,NJT !Loop over each derivative direction.
          XNO(1)=0.0d0
          XNO(2)=0.0d0
          XNO(3)=0.0d0
          XNO(nj_grad)=1.0d0
          DO ng=1,NGT(nb1j)
            RAD(ng)=DSQRT(XR(1,ng)*XR(1,ng)+XR(2,ng)*XR(2,ng)+
     '        XR(3,ng)*XR(3,ng))
            DRDN(ng)=-XR(nj_grad,ng)/RAD(ng) !dR/dn0 where n0 is a unit
          ENDDO !ng                !vector in the nj_grad direction.

C***      GK matrix components
          ns=0
          DO nn=1,NNT(nbuhp)     !Dependent variable loop
            npnn=NPNE(nn,nbuhp,ne)
            DO nk1=1,NKT(nn,nbuhp)
              MK=NKHE(nk1,nn,nh,ne)
              ns=ns+1
              IF(mk.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,1)-KTYP93(1,nr),1))
     '          THEN
                SUMU=0.0d0
                DO ng=1,NGT(nbuh)
                  SUMU=SUMU-PG(ns,1,ng,nbuh)*HYPGREEN(IGREN(nr),0,
     '              CE(1,ne),RAD(ng),DET(nbuhp,0,ng,1),XN(1,ng),XNO,
     '              XR(1,ng))*RG(ng)*WG(ng,nbuh)
                ENDDO !End of ng loop
                GRADPHI(nj_grad)=GRADPHI(nj_grad)+SUMU*ZE(ns,nhx)
              ENDIF !End of MK > 0 loop
            ENDDO !End of nk loop
          ENDDO !End of nn loop

C***      GQ matrix components
          ns=0
          DO nn=1,NNT(nbqhp) !normal derivative loop
            npnn=NPNE(nn,nbqhp,ne)
            DO nk1=1,NKT(nn,nbqhp)
              MK=NKHE(nk1,nn,nh,ne)
              ns=ns+1
              IF(mk.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,1)-KTYP93(1,nr),1))
     '          THEN
                SUMQ=0.0d0
                DO ng=1,NGT(nbqh)
                  SUMQ=SUMQ+PG(ns,1,ng,nbqh)*DGREEN(IGREN(nr),1,nh,nh,
     '              CE(1,ne),DRDN(ng),RAD(ng),DET(nbqhp,0,ng,1),
     '              XN(1,ng),XR(1,ng))*RG(ng)*WG(ng,nbqh)
                ENDDO !End of ng loop
                GRADPHI(nj_grad)=GRADPHI(nj_grad)+SUMQ*ZF(ns,nh)
              ENDIF !End of MK > 0 loop.
            ENDDO !End of nk loop
          ENDDO !End of nn loop
        ENDDO !End of nj_grad loop (i.e. derivative direction loop).
      ENDDO !End of element loop

C***  Only want gradients that are definitely internal
      IF(DABS(SUMC).LT.0.9d0) THEN
        DO nj=1,NJT
          GRADPHI(nj)=0.0d0
        ENDDO
      ENDIF

      CALL EXITS('GRADDOMSOLOPTI3D')
      RETURN
9999  CALL ERRORS('GRADDOMSOLOPTI3D',ERROR)
      CALL EXITS('GRADDOMSOLOPTI3D')
      RETURN 1
      END


