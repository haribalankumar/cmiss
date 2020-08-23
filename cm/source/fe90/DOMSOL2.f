      SUBROUTINE DOMSOL2(NBH,NBJ,NEELEM,NHE,NHP,NKH,NKHE,NKJE,NLL,NPF,
     '  NP_INTERFACE,NPNE,NPNODE,nr,NRE,NVHE,NVHP,NVJE,NW,nx,NYNE,NYNP,
     '  CE,CURVCORRECT,DET,DL,DRDN,PG,RAD,RD,RG,SE,WG,XA,XE,XG1,XN,XP,
     '  XPFP,XR,YD,YP,ZA,ZE,ZF,ZP,ERROR,*)

C#### Subroutine: DOMSOL2
C###  Description:
C###    DOMSOL2 calculates BE solution at the point XPFP (Stored in YD).
C**** Created by Martin Buist - June 2001

      IMPLICIT NONE

      INCLUDE 'b10.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NHE(NEM),NHP(NPM),NKH(NHM,NPM,NCM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),CURVCORRECT(2,2,NNM,NEM),DET(NBFM,0:NNM,NGM,6),
     '  DL(3,NLM),DRDN(NGM),PG(NSM,NUM,NGM,NBM),RAD(NGM),RD(NGM),
     '  RG(NGM),SE(NSM,NBFM,NEM),WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG1(NJM,NUM,NGM),XN(NJM,NGM),XP(NKM,NVM,NJM,NPM),
     '  XPFP(*),XR(NJM,NGM),YD(NHM),YP(NYM,NIYM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZF(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER intscheme,nb2,nbuh,nbuhp,nb1j,nb1jp,nbbem,nbqh,nbqhp,
     '  ne,ng,nh,nhx,ni,NITB,nj,nnmin,
     '  noelem,np,ns,NU1(0:3)
      REAL*8 DXXI(3,3),MINDIST,
     '  SUMXG
      LOGICAL DUMMY,INTERFACE

      DATA NU1/1,2,4,7/

      CALL ENTERS('DOMSOL2',*9999)

      DO nh=1,NHM
        YD(nh)=0.0d0
      ENDDO

C***  Transfer dependent variable to ZP
      CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NVHP,nx,NYNE,NYNP,
     '  YP,ZA,ZP,ERROR,*9999)

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        np=NPNE(1,NBJ(1,ne),ne)
        INTERFACE=(NP_INTERFACE(np,0).GT.1).AND.
     '    (NP_INTERFACE(np,1).NE.nr)

C***    Transfer global geom parameters to local parameters
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

C***    Determine integration scheme to use for element ne
        DUMMY=ADAPINT
        ADAPINT=.FALSE.
        CALL DIST(intscheme,NBJ(1,ne),NLL(1,ne),nnmin,NPNE(1,1,ne),nr,
     '    NW(ne,2),DL,MINDIST,XP,XPFP,ERROR,*9999)
        ADAPINT=DUMMY

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

C***    Loop over dependent variables
        DO nhx=1,NHE(ne)
          nh=NH_LOC(nhx,nx)

C***      Decide basis functions from the parent family basis
C***      function for dependent (u) and normal derivative (q).
          nbuhp=NFBASE(1,NBASEF(NBH(nh,1,ne),nbbem))
          nbqhp=NFBASE(1,NBASEF(NBH(nh,2,ne),nbbem))

C***      nbuh,nbqh and nb1j are the global
C***      basis functions with the appropriate quadrature for the
C***      dependent (H), normal derivative (Q) and geometric (J)
C***      variables.
          nb2=nbbem
          nbuh=NBASEF(NBH(nh,1,ne),nb2)
          nbqh=NBASEF(NBH(nh,2,ne),nb2)
          nb1j=NBASEF(NBJ(1,ne),nb2)

C***      set up the gauss point dependent arrays.
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
C***        RG(ng) stores the Jacobian for a length or area
C***        integral at Gauss point ng.
C***        DXXI(nj,ni) contains the values of dXj/dPSIi at
C***        the Gauss point.
            DO ni=1,NITB
              DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                DXXI(nj,ni)=XG1(nj,NU1(ni),ng)
              ENDDO !nj
            ENDDO !ni
            IF(NITB.EQ.1) THEN
C***          Jacobian for 1D integral in 2D space
              RG(ng)=DSQRT(DXXI(1,1)*DXXI(1,1)+DXXI(2,1)*DXXI(2,1))
            ELSE !Calculate cross product of 2 tangent vectors
C***          Jacobian for 2D integral in 3D space
              RG(ng)=DSQRT(
     '          (DXXI(2,1)*DXXI(3,2)-DXXI(2,2)*DXXI(3,1))*
     '          (DXXI(2,1)*DXXI(3,2)-DXXI(2,2)*DXXI(3,1))+
     '          (DXXI(1,2)*DXXI(3,1)-DXXI(1,1)*DXXI(3,2))*
     '          (DXXI(1,2)*DXXI(3,1)-DXXI(1,1)*DXXI(3,2))+
     '          (DXXI(1,1)*DXXI(2,2)-DXXI(1,2)*DXXI(2,1))*
     '          (DXXI(1,1)*DXXI(2,2)-DXXI(1,2)*DXXI(2,1)))
            ENDIF
C***        Find the unit outward normal.
            CALL NORMAL(ne,nr,NW,XG1(1,1,ng),XN(1,ng),INTERFACE,
     '        ERROR,*9999)
            IF(JTYP4.GT.1) THEN
              IF(JTYP4.EQ.2) THEN
                RD(ng)=XG1(2,1,ng) !rotat symm about x/r axis
              ELSE IF(JTYP4.EQ.3) THEN
                RD(ng)=XG1(1,1,ng) !rotat symm about y/z axis
              ENDIF
            ENDIF
          ENDDO !End ng loop

          CALL XEPGKGQDOMSOL(nb1j,nbqh,nbqhp,nbuh,nbuhp,nh,NHE(ne),nhx,
     '      NKHE(1,1,1,ne),NKH,NPNE(1,1,ne),nr,nx,CE(1,ne),DET,DRDN,
     '      PG(1,1,1,nbqh),PG(1,1,1,nbuh),RAD,RD,RG,WG,XG1,XN,XPFP,XR,
     '      YD,ZE,ZF,ERROR,*9999)

        ENDDO !End of nh loop
      ENDDO !ne

      CALL EXITS('DOMSOL2')
      RETURN
9999  CALL ERRORS('DOMSOL2',ERROR)
      CALL EXITS('DOMSOL2')
      RETURN 1
      END


