      SUBROUTINE GRADDOMSOL(NBH,NBJ,NEELEM,NHE,NHP,NKH,NKHE,NKJE,NLL,
     '  NPF,NP_INTERFACE,NPNE,NPNODE,nr,NRE,NVHE,NVHP,NVJE,NW,nx,NYNE,
     '  NYNP,CE,CURVCORRECT,DET,DL,DRDN,GRADPHI,PG,RAD,RG,SE,
     '  WG,XA,XE,XG1,XN,XP,XPFP,XR,YP,ZA,ZE,ZF,ZP,ERROR,*)

C#### Subroutine: GRADDOMSOL
C###  Description:
C###    GRADDOMSOL calculates the gradient of the potential at the
C###    internal point XPFP(nj).
C**** Created by Martin Buist - June 2001

      IMPLICIT NONE

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
     '  NHE(NEM),NHP(NPM),NKH(NHM,NPM,NCM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),NPF(9,NFM),
     '  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),CURVCORRECT(2,2,NNM,NEM),DET(NBFM,0:NNM,NGM,6),
     '  DL(3,NLM),DRDN(NGM),GRADPHI(*),PG(NSM,NUM,NGM,NBM),RAD(NGM),
     '  RG(NGM),SE(NSM,NBFM,NEM),WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG1(NJM,NUM,NGM),XN(NJM,NGM),XP(NKM,NVM,NJM,NPM),
     '  XPFP(*),XR(NJM,NGM),YP(NYM,NIYM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZF(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER intscheme,MK,nbuh,nbuhp,nb1j,nb1jp,nbbem,nbqh,nbqhp,ne,ng,
     '  nh,nhx,ni,NITB,nj,nj_grad,nk1,nn,NNMIN,noelem,np,npnn,ns,
     '  NU1(0:3)
      REAL*8 DGREEN,DXXI(3,3),HYPGREEN,MINDIST,SUM,SUMQ,SUMR,SUMU,SUMXG,
     ' XNO(3)
      LOGICAL DUMMY,INTERFACE,PROXIM

      DATA NU1/1,2,4,7/

      CALL ENTERS('GRADDOMSOL',*9999)

      nhx=1
      nh=NH_LOC(nhx,nx)
      DO nj=1,NJT
        GRADPHI(nj)=0.0d0
      ENDDO
      PROXIM=.FALSE.

      CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NVHP,nx,NYNE,NYNP,
     '  YP,ZA,ZP,ERROR,*9999)

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        np=NPNE(1,NBJ(1,ne),ne)
        INTERFACE=(NP_INTERFACE(np,0).GT.1).AND.
     '    (NP_INTERFACE(np,1).NE.nr)
        IF(IGREN(nr).EQ.1) CE(1,ne)=1.0d0
        IF(IGREN(nr).EQ.2) CE(1,ne)=1.0d0

C***    Transfer global parameters to local parameters
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '    NRE(ne),NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,
     '    ERROR,*9999)
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
        DUMMY=ADAPINT
        ADAPINT=.FALSE.
        CALL DIST(intscheme,NBJ(1,ne),NLL(1,ne),nnmin,NPNE(1,1,ne),nr,
     '    NW(ne,2),DL,MINDIST,XP,XPFP,ERROR,*9999)
        ADAPINT=DUMMY

C***    Prevent crazy results at nodes
        IF(DABS(MINDIST).LE.LOOSE_TOL) PROXIM=.TRUE.

C***    No splitting or adaptive int. implemented here
        IF(intscheme.EQ.1) THEN
          intscheme=3
        ELSEIF(intscheme.EQ.2) THEN
          intscheme=3
        ENDIF

C***    Decide upon the appropriate basis for the integration scheme
        CALL QUADBE(intscheme,NBJ(1,ne),nbbem,ERROR,*9999)
        nb1jp=NFBASE(1,NBASEF(NBJ(1,ne),nbbem))
        NITB=NIT(nb1jp)
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
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            DO ni=0,NITB
              SUMXG=0.0d0
              DO ns=1,NST(nb1jp)
                SUMXG=SUMXG+PG(ns,NU1(ni),ng,nb1j)*XE(ns,nj)
              ENDDO !ns
              XG1(nj,NU1(ni),ng)=SUMXG
            ENDDO !ni
          ENDDO !nj
C***      RG(ng) stores the Jacobian for a length or area
C***      integral at Gauss point ng.
C***      DXXI(nj,ni) contains the values of dXj/dPSIi at
C***      the Gauss point.
          DO ni=1,NITB
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              DXXI(nj,ni)=XG1(nj,NU1(ni),ng)
            ENDDO
          ENDDO
          IF(NITB.EQ.1) THEN
C***        Jacobian for 1D integral in 2D space
            RG(ng)=DSQRT(DXXI(1,1)*DXXI(1,1)+DXXI(2,1)*DXXI(2,1))
          ELSE !Calculate cross product of 2 tangent vectors
C***        Jacobian for 2D integral in 3D space
            RG(ng)=DSQRT(
     '        (DXXI(2,1)*DXXI(3,2)-DXXI(2,2)*DXXI(3,1))*
     '        (DXXI(2,1)*DXXI(3,2)-DXXI(2,2)*DXXI(3,1))+
     '        (DXXI(1,2)*DXXI(3,1)-DXXI(1,1)*DXXI(3,2))*
     '        (DXXI(1,2)*DXXI(3,1)-DXXI(1,1)*DXXI(3,2))+
     '        (DXXI(1,1)*DXXI(2,2)-DXXI(1,2)*DXXI(2,1))*
     '        (DXXI(1,1)*DXXI(2,2)-DXXI(1,2)*DXXI(2,1)))
          ENDIF
C***      Find the unit outward normal.
          CALL NORMAL(ne,nr,NW,XG1(1,1,ng),XN(1,ng),INTERFACE,
     '      ERROR,*9999)
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            XR(nj,ng)=XG1(nj,1,ng)-XPFP(nj)
          ENDDO
        ENDDO !End of ng loop for variables not depending on which
              !direction the equation was differentiated.

        DO nj_grad=1,NJT !Loop over each derivative direction.
          DO nj=1,NJT
            XNO(nj)=0.0d0
          ENDDO
          XNO(nj_grad)=1.0d0
          DO ng=1,NGT(nb1j)
            SUM=0.0d0
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              SUM=SUM+XR(nj,ng)*XR(nj,ng)
            ENDDO
            SUMR=XR(nj_grad,ng) !the njth component of XR
            RAD(ng)=DSQRT(SUM)
            DRDN(ng)=-SUMR/RAD(ng) !dR/dn0 where n0 is a unit vector
          ENDDO !ng                !in the nj_grad direction.

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

      IF(PROXIM) THEN
        DO nj=1,NJT
          GRADPHI(nj)=0.0d0
        ENDDO
      ENDIF

      CALL EXITS('GRADDOMSOL')
      RETURN
9999  CALL ERRORS('GRADDOMSOL',ERROR)
      CALL EXITS('GRADDOMSOL')
      RETURN 1
      END


