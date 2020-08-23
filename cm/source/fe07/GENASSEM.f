      SUBROUTINE GENASSEM(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IBT,IDO,INP,
     '  ISC_GD,ISC_GK,ISC_GKK,ISC_GM,ISC_GQ,ISR_GD,ISR_GK,ISR_GKK,
     '  ISR_GM,ISR_GQ,LGE,NBH,NBJ,NDET,NDIPOLES,NEELEM,NENP,NGAP,NHE,
     '  NHP,NKH,NKHE,NKJE,NLL,NONY,NORD,NPF,NP_INTERFACE,NPB,NPNE,
     '  NPNODE,NPNY,nr,nr_gkk,NRE,NRLIST,NVHE,NVHP,NVJE,NW,nx,NYNE,NYNP,
     '  NYNR,CE,CG,CGE,CONY,CP,CURVCORRECT,DET,DIPOLE_CEN,DIPOLE_DIR,DL,
     '  DRDN,DRDNO,ED,EM,ER,ES,GD,GK,GKK,GM,GQ,GR,PG,RAD,RD,RG,SE,WG,XA,
     '  XE,XG,XG1,XIG,XN,XN_GRAD,XP,XR,XR_GRAD,YG,ZA,ZE,ZG,ZP,ERROR,*)

C#### Subroutine: GENASSEM
C###  Description:
C###    <HTML>
C###    GENASSEM assembles the global unreduced matrices GK,GM,etc.
C###    for class nx and region nr. This routine calls various assemble
C###    routines depending on the solution process and problem type in
C###    the region.  These are:
C###    <UL>
C###    <LI>ASSEMBLE1  Static linear FEM problems
C###    <LI>ASSEMBLE2  Static linear BEM problems
C###    <LI>ASSEMBLE3  Time dependent FEM problems
C###    <LI>ASSEMBLE4  Time dependent BEM problems
C###    <LI>ASSEMBLE5  Nonlinear FEM problems
C###    <LI>ASSEMBLE6  Fourier analysis/Complex problems
C###    <LI>ASSEMBLE7  Modal analysis problems
C###    <LI>ASSEMBLE8  Laplace transform analysis problems
C###    </UL>
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'solv00.cmn'
!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM,NXM),
     '  IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISC_GD(NISC_GDM),ISC_GK(NISC_GKM),
     '  ISC_GKK(NISC_GKKM),ISC_GM(NISC_GMM),ISC_GQ(NISC_GQM),
     '  ISR_GD(NISR_GDM),ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),
     '  ISR_GM(NISR_GMM),ISR_GQ(NISR_GQM),LGE(NHM*NSM,NRCM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NDET(NBFM,0:NNM),
     '  NDIPOLES(NRM,NXM),
     '  NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),NGAP(NIM,NBM),
     '  NHE(NEM),NHP(NPM),NKH(NHM,NPM,NCM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NONY(0:NOYM,NYM,NRCM),NORD(5,NE_R_M),NPF(9,NFM),NPB(0:NP_R_M,5),
     '  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),nr,nr_gkk,NRE(NEM),
     '  NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),
     '  CONY(0:NOYM,NYM,NRCM),CP(NMM,NPM),
     '  CURVCORRECT(2,2,NNM,NEM),DET(NBFM,0:NNM,NGM,6),
     '  DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),DL(3,NLM),
     '  DRDN(NGM),DRDNO(NGM,NKM),ED(NHM*NSM,NHM*NSM),
     '  EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),
     '  GD(NZ_GD_M),GK(NZ_GK_M),GKK(NZ_GKK_M),GM(NZ_GM_M),
     '  GQ(NZ_GQ_M),GR(NYROWM),PG(NSM,NUM,NGM,NBM),RAD(NGM),
     '  RD(NGM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XG1(NJM,NUM,NGM),XIG(NIM,NGM,NBM),XN(NJM,NGM),
     '  XN_GRAD(NJM,NGM),XP(NKM,NVM,NJM,NPM),XR(NJM,NGM),
     '  XR_GRAD(NJM,NGM),YG(NIYGM,NGM,NEM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 TIME
      LOGICAL DYNAM1,DYNAM2,GQ_ASSEM,UPDATE_MATRIX

      CALL ENTERS('GENASSEM',*9999)

      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,'(/'' Region '',I1,'':'')') nr
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      GQ_ASSEM=.FALSE.
      IF(NRLIST(0).GT.1) THEN
        !more than one region being solved for.  Need to determine if
        !a GQ matrix needs to be assembled for any FEM region (i.e.
        !there is some coupling between an FE and a BE region or the
        !FE region is working with point values of fluxes). Currently
        !not checking the FE point values since these are not yet
        !handled correctly in ipini3 (the poitn values are converted
        !to integrated value, and the basis function used is that
        !for the geomtery).
        GQ_ASSEM=COUPLED_BEM(nx) !coupled_bem set up in DESOLV
      ENDIF

      IF(ITYP5(nr,nx).EQ.1.OR.ITYP5(nr,nx).EQ.4) THEN !problem is static
                                                      !or quasistatic.
        IF(ITYP6(nr,nx).EQ.1) THEN !problem is linear

C LKC 15-APR-2002 have the TIME set to T0 for all problems...
C  for a static problem this is normally 0.0 but
C  allows it to be solved at a specific instance in time
C  specified at the command line from solve
C  T0 should be == TSTART anyway.
C          IF(ITYP5(nr,nx).EQ.4) THEN !Quasi-static
C            TIME=T0
C          ELSE
C            TIME=0.0d0
C          ENDIF
          TIME=TSTART


          IF(ITYP4(nr,nx).EQ.1) THEN !problem is f.e.s only
            CALL ASSEMBLE1(IBT,IDO,INP,ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,
     '        ISR_GKK,ISR_GQ,NBH,NBJ,NEELEM,NHE,NHP,NKJE,NONY,NORD,NPF,
     '        NP_INTERFACE,NPNE,NPNY,nr,nr_gkk,NRE,NVHE,NVJE,NW,nx,
     '        NYNE,NYNP,NYNR,CE,CGE,CONY,CP,CURVCORRECT,GK,GKK,GQ,GR,PG,
     '        SE,WG,XA,XP,YG,GQ_ASSEM,.TRUE.,.TRUE.,ERROR,*9999)
          ELSE IF(ITYP4(nr,nx).EQ.2.OR.ITYP4(nr,nx).EQ.3) THEN !BEM
            CALL ASSEMBLE2(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     '        IBT,IDO,INP,ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,
     '        NBH,NBJ,NDET,NDIPOLES,NEELEM,NENP,NGAP,NHE,NHP,NKH,
     '        NKHE,NKJE,NLL,NONY,NPF,NP_INTERFACE,NPNE,NPNODE,NPNY,
     '        nr,nr_gkk,NRE,NVHP,NVJE,NW,nx,NYNE,NYNP,NYNR,CE,CONY,
     '        CURVCORRECT,DET,DIPOLE_CEN,DIPOLE_DIR,DL,GD,GK,
     '        GKK,GQ,PG,SE,TIME,WG,XA,XE,XG,XIG,
     '        XP,.TRUE.,.TRUE.,ERROR,*9999)
          ENDIF
        ENDIF !ITYP6

      ELSE IF(ITYP5(nr,nx).EQ.2) THEN !problem is time-dependent
        IF(ITYP4(nr,nx).EQ.1) THEN !problem is finite elements only
          IF(ITYP2(nr,nx).NE.4) THEN  !finite diffs not characteristics
            IF(ITYP2(nr,nx).NE.8.AND.ITYP2(nr,nx).NE.9) THEN !parabolic FEM

C!!!cpb 9/12/94 These need to be set correctly ?

              DYNAM1=.TRUE.
              DYNAM2=.TRUE.
              UPDATE_MATRIX=.TRUE.

              CALL ASSEMBLE3(IBT,IDO,INP,ISC_GD,ISC_GK,ISC_GM,ISR_GD,
     '          ISR_GK,ISR_GM,LGE,NBH,NBJ,NEELEM,NHE,NKHE,NKJE,NORD,
     '          NPF,NPNE,nr,NVHE,NVJE,NW,nx,NYNE,
     '          NYNP,NYNR,CE,CG,CGE,CP,ED,EM,ER,ES,GD,GK,GM,
     '          GR,PG,RG,SE,WG,XA,XE,XG,XP,YG,ZA,ZE,ZG,ZP,DYNAM1,
     '          DYNAM2,UPDATE_MATRIX,ERROR,*9999)

            ENDIF
          ELSE IF(ITYP4(nr,nx).EQ.2.OR.ITYP4(nr,nx).EQ.3) THEN !BEM
            CALL ASSEMBLE4(ERROR,*9999)
          ENDIF
        ENDIF

      ELSE IF(ITYP5(nr,nx).EQ.3) THEN !problem is modal analysis
        IF(ITYP2(nr,nx).EQ.1) THEN !linear elasticity
          DYNAM1=.FALSE.
          DYNAM2=.TRUE.
          UPDATE_MATRIX=.TRUE.
          CALL ASSEMBLE6(IBT,IDO,INP,ISC_GD,ISC_GK,ISC_GM,ISR_GD,
     '      ISR_GK,ISR_GM,LGE,NBH,NBJ,NEELEM,NHE,NKJE,NORD,NPF,NPNE,
     '      nr,NVHE,NVJE,NW,nx,NYNE,NYNP,NYNR,CE,CG,CGE,CP,ED,EM,ER,ES,
     '      GD,GK,GM,GR,PG,RG,SE,WG,XA,XE,XG,XP,YG,ZE,ZG,DYNAM1,
     '      DYNAM2,UPDATE_MATRIX,ERROR,*9999)
        ELSE IF(ITYP2(nr,nx).EQ.8) THEN !Vocal tract equations
          ERROR='>>Not implemented'
          GOTO 9999
        ELSE
          ERROR='>>Not implemented'
          GOTO 9999
        ENDIF
C cpb 1/5/95 replace Fourier analysis with Quasi-static analysis
C      ELSE IF(ITYP5(nr,nx).EQ.4) THEN !problem is Fourier analysis
C        CALL ASSEMBLE7(ERROR,*9999)

      ELSE IF(ITYP5(nr,nx).EQ.4) THEN !problem is Quasi-static analysis

      ELSE IF(ITYP5(nr,nx).EQ.5) THEN !problem is Laplace transform
C       Not implemented (Would call ASSEMBLE8)
        ERROR='>> Laplace transform analysis is not implemented'
        GOTO 9999

      ELSE IF(ITYP5(nr,nx).EQ.6) THEN !problem is buckling analysis
        CALL ASSEMBLE7(ERROR,*9999)

      ENDIF !ITYP5

      ASSEMBLE_GLOBAL(nr,nx)=.TRUE.

      CALL EXITS('GENASSEM')
      RETURN
 9999 CALL ERRORS('GENASSEM',ERROR)
      CALL EXITS('GENASSEM')
      RETURN 1
      END


