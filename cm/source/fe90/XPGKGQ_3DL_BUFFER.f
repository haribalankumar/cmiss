      SUBROUTINE XPGKGQ_3DL_BUFFER(IBT,IDO,INP,ISC_GK,ISC_GKK,ISC_GQ,
     '  ISR_GK,ISR_GKK,ISR_GQ,LGKE,LGQE,NBH,NBJ,NDET,NEELEM,NENP,NGAP,
     '  NHP,NKH,NKHE,NKJE,NLL,NONY,np,NP_INTERFACE,NPF,NPNE,NPNY,nr,
     '  nr_gkk,NRE,NVJE,NW,nx,NYNE,NYNP,CE,CONY,CURVCORRECT,DET,
     '  DET_ADAPT,DL,DRDN,DRDNO,GK,GKK,GKES,GQ,GQES,PG,PG_J,PG_Q,PG_U,
     '  RAD,RG,SE,WG,XA,XE,XG,XG1,XIG,XIG_J,XIG_Q,XIG_U,XN,XN_GRAD,XP,
     '  XR,XR_GRAD,INTERFACE,FULL_EQUATIONS,ERROR,*)

C#### Subroutine: XPGKGQ_3DL_BUFFER
C###  Description:
C###    XPGKGQ_3DL_BUFFER is the buffered efficient version of XPGKGQ
C###    for 3D Laplace's equation in rectangular cartesian coordinates.
C###  See-Also: XPGKGQ_3DL

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'solv00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISC_GK(NISC_GKM),ISC_GKK(NISC_GKKM),ISC_GQ(NISC_GQM),
     '  ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),
     '  ISR_GQ(NISR_GQM),LGKE(-NKM:NHM*NSM,NKM,2),
     '  LGQE(-NKM:NHM*NSM,NKM,2),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NDET(NBFM,0:NNM),NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),
     '  NGAP(NIM,NBM),NHP(NPM),NKH(NHM,NPM,NCM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),NONY(0:NOYM,NYM,NRCM),np,
     '  NP_INTERFACE(0:NPM,0:3),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNY(0:6,NYM,0:NRCM),nr,nr_gkk,NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),CONY(0:NOYM,NYM,NRCM),CURVCORRECT(2,2,NNM,NEM),
     '  DET(NBFM,0:NNM,NGM,6),DET_ADAPT(NBFM,0:NNM,NGM),DL(3,NLM),
     '  DRDN(NGM),DRDNO(NGM,NKM),GK(NZ_GK_M),GKK(NZ_GKK_M),
     '  GKES(-NKM:NHM*NSM,NKM),GQ(NZ_GQ_M),GQES(-NKM:NHM*NSM,NKM),
     '  PG(NSM,NUM,NGM,NBM),PG_J(NSM,NUM,NGM),PG_Q(NSM,NUM,NGM),
     '  PG_U(NSM,NUM,NGM),RAD(NGM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XG1(NJM,NUM,NGM),XIG(NIM,NGM,NBM),XIG_J(NIM,NGM),
     '  XIG_Q(NIM,NGM),XIG_U(NIM,NGM),XN(NJM,NGM),XN_GRAD(NJM,NGM),
     '  XP(NKM,NVM,NJM,NPM),XR(NJM,NGM),XR_GRAD(NJM,NGM)
      CHARACTER ERROR*(*)
      LOGICAL INTERFACE,FULL_EQUATIONS
!     Local Variables
      INTEGER GETNYR,IADAPTIVE,intscheme,ISP,NB2,nb1j,nb1jp,nbbem,
     '  nbqh,nbqhp,nbuh,nbuhp,NDTOT,ne,ng,nh,nhs,NHSTGK,NHSTGQ,nhx,nkk,
     '  NNMIN,NNSP,no1,no2,noy1,noy2,noelem,NSPLIT,ny1,ny2,ny3,nz,nzz
      REAL*8 co1,co2,D(2),DDOT,DSDX(3,3),DXDS(3,3),DXXI(3,3),MINDIST,
     '  VLENGTH,XIMIN(2),XNO(3,4),XPFP(3)
      LOGICAL GEN_CONVENT_EQUATION,GEN_DERIV_EQUATION,CALC_XR,USEADAPINT

      EXTERNAL DDOT

      CALL ENTERS('XPGKGQ_3DL_BUFFER',*9999)

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(XPGKGQ_3DL_BUFFER_1)
        WRITE(OP_STRING,'(/'' Node number '',I5)') np
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(XPGKGQ_3DL_BUFFER_1)
      ENDIF

C***  Determine location of the singularity
      XPFP(1)=XP(1,1,1,np)
      XPFP(2)=XP(1,1,2,np)
      XPFP(3)=XP(1,1,3,np)

C*** Determine what equations to generate
      IF(KTYP92.EQ.3) THEN
        GEN_CONVENT_EQUATION=.FALSE.
      ELSE
        GEN_CONVENT_EQUATION=.TRUE.
      ENDIF
      GEN_DERIV_EQUATION=HYP.AND.(NKH(NH_LOC(1,nx),np,1).GT.1)

C*** Calculate direction to differentiate BIE equation in if necessary
      IF(GEN_DERIV_EQUATION) THEN
        CALL XPXNO(NBH,NENP,NP_INTERFACE,np,nr,NW,nx,DSDX,DXDS,XG,
     '    XNO,XP,ERROR,*9999)
      ENDIF

C***  Loop over elements in the region
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)

        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(XPGKGQ_3DL_BUFFER_2)
          WRITE(OP_STRING,'(/'' Element number '',I5)') ne
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(XPGKGQ_3DL_BUFFER_2)
        ENDIF

C***    Determine integration scheme to use for this element
        CALL DIST(intscheme,NBJ(1,ne),NLL(1,ne),NNMIN,
     '    NPNE(1,1,ne),nr,NW(ne,2),DL,MINDIST,XP,XPFP,ERROR,*9999)

C***    Transfer global parameters to local parameters

        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '    NRE(ne),NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,
     '    ERROR,*9999)

C***    Decide upon the appropriate basis for the integration scheme
        CALL QUADBE(intscheme,NBJ(1,ne),nbbem,ERROR,*9999)

        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(XPGKGQ_3DL_BUFFER_3)
         WRITE(OP_STRING,'(/'' Integration scheme '',I1,'', nbbem '','
     '      //'I2)') intscheme,nbbem
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(XPGKGQ_3DL_BUFFER_3)
        ENDIF

        IF(intscheme.NE.1) THEN
          NNMIN=0
        ENDIF

        nb1jp=NFBASE(1,NBASEF(NBJ(1,ne),nbbem))
C        NITB=NIT(nb1jp)


C***    Loop over dependent variables
        DO nhx=1,NHP(np)
          nh=NH_LOC(nhx,nx)

C***      Decide basis functions from the parent family basis function
C***      for dependent (u) and normal derivative (q).
          nbuhp=NFBASE(1,NBASEF(NBH(nh,1,ne),nbbem))
          nbqhp=NFBASE(1,NBASEF(NBH(nh,2,ne),nbbem))

C***      Setup adaptive integration information if nescessary.
          USEADAPINT=.FALSE.
          IF(intscheme.EQ.2.AND.NW(ne,2).NE.11.AND.NW(ne,2).NE.13) THEN
C***        Don't use this on distorted cubic linear elements.
C***        Find the minimum distance in the element to the singular
C***        point for Telles' rule.
            USEADAPINT=.TRUE.
            CALL TELLES(IBT,IDO,INP,20,NBJ(1,ne),NLL(1,ne),
     '        NW(ne,2),D,DL,1.0d-6,XE,XIMIN,XPFP,ERROR,*9999)
C***        Set up adaptive basis functions based on Telles' rule.
            nb1j=NBASEF(NBJ(1,ne),nbbem)
            nbuh=NBASEF(NBH(NH_LOC(1,nx),1,ne),nbbem)
            nbqh=NBASEF(NBH(NH_LOC(1,nx),2,ne),nbbem)
            CALL GAUSS11(IBT,IDO,INP,
     '        nb1j,nb1jp,nbqh,nbuh,NGAP,D,DET,DET_ADAPT,PG_J,
     '        PG_Q,PG_U,XIG,XIG_J,XIG_Q,XIG_U,XIMIN,ERROR,*9999)
            IADAPTIVE=1
          ELSE
            IADAPTIVE=0
          ENDIF

          IF(GEN_DERIV_EQUATION) THEN
C*** Calculate the number of derivative equations to generate (NDTOT)
            NDTOT=MAX(NKH(nh,np,1)-1-KTYP93(1,nr),1)
            IF(.NOT.FULL_EQUATIONS.AND.NP_INTERFACE(np,2)
     '        .EQ.nr.AND.KTYP92.EQ.3) THEN
              NDTOT=1
            ENDIF
            IF(KTYP92.EQ.3) THEN !All derivative equations
              NDTOT=NDTOT+1
            ENDIF
          ELSE
            NDTOT=0
          ENDIF

C cpb 18/12/97 Adding BEM stiffness matrices

          CALL MELGEB(ISC_GK,ISR_GK,LGKE,nbuhp,1,NDTOT,nh,NHSTGK,
     '      NHP(np),NKHE(1,1,nh,ne),NKH(1,1,1),np,NPNE(1,nbuhp,ne),
     '      NPNY,nr,nx,NYNE,NYNP,NZ_GK_M,ERROR,*9999)
          CALL MELGEB(ISC_GQ,ISR_GQ,LGQE,nbqhp,2,NDTOT,nh,NHSTGQ,
     '      NHP(np),NKHE(1,1,nh,ne),NKH(1,1,2),np,NPNE(1,nbqhp,ne),
     '      NPNY,nr,nx,NYNE,NYNP,NZ_GQ_M,ERROR,*9999)

          DO nkk=1,NDTOT+1
            DO nhs=-NDTOT-1,NHSTGK
              GKES(nhs,nkk)=0.0d0
            ENDDO !nhs
            DO nhs=-NDTOT-1,NHSTGQ
              GQES(nhs,nkk)=0.0d0
            ENDDO !nhs
          ENDDO !nkk

C***      Calculate the element integrals for the current node and
C***      dependent variable.
          ISP=0
          IF(intscheme.EQ.1) THEN !Element splitting

C***        Find number of split elements used at nodes less than
C***        NNMIN
            DO nnsp=0,NNMIN-1
              ISP=ISP+NDET(nb1jp,nnsp)
            ENDDO !nnsp
          ENDIF

C***      Loop over number of split elements (if any)
          DO nsplit=1,NDET(nb1jp,NNMIN)

            nb2=nbbem+ISP
C***        Determine basis function for the dependent variable for the
C***        split element.
C***        nbuh,nbqh and nb1j are the global basis functions with the
C***        appropriate quadrature for the dependent (H), normal
C***        derivative (Q) and geometric (J) variables.
            nbuh=NBASEF(NBH(nh,1,ne),nb2)
            nbqh=NBASEF(NBH(nh,2,ne),nb2)
C cpb 2/11/98 Shouldn't use nh for NBJ. Use nj=1 instead for now.
C            nb1j=NBASEF(NBJ(nh,ne),nb2)
            nb1j=NBASEF(NBJ(1,ne),nb2)
C***        If element spliting/adaptive integration or the first node
C***        in the integration scheme set up the gauss point dependent
C***        arrays.
            IF(intscheme.EQ.1.OR.intscheme.EQ.2.OR..TRUE.) THEN
              IF(intscheme.EQ.1.OR.intscheme.EQ.2) CALC_XR=.TRUE.
              DO ng=1,NGT(nb1j)
c cpb 9/7/95 reorganising and unrolling loops and using blas
C                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
C                    DO ni=0,NITB
C                      SUMXG=0.0d0
C                      DO ns=1,NST(nb1jp)
C                        SUMXG=SUMXG+PG(ns,NU1(ni),ng,nb1j)*XE(ns,nj)
C                      ENDDO !ns
C                      XG1(nj,NU1(ni),ng)=SUMXG
C                    ENDDO !ni
C                  ENDDO !nj
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
                  XG1(1,1,ng)=DDOT(NST(nb1jp),PG(1,1,ng,nb1j),1,
     '              XE(1,1),1)
                  XG1(2,1,ng)=DDOT(NST(nb1jp),PG(1,1,ng,nb1j),1,
     '              XE(1,2),1)
                  XG1(3,1,ng)=DDOT(NST(nb1jp),PG(1,1,ng,nb1j),1,
     '              XE(1,3),1)
                  XG1(1,2,ng)=DDOT(NST(nb1jp),PG(1,2,ng,nb1j),1,
     '              XE(1,1),1)
                  XG1(2,2,ng)=DDOT(NST(nb1jp),PG(1,2,ng,nb1j),1,
     '              XE(1,2),1)
                  XG1(3,2,ng)=DDOT(NST(nb1jp),PG(1,2,ng,nb1j),1,
     '              XE(1,3),1)
                  XG1(1,4,ng)=DDOT(NST(nb1jp),PG(1,4,ng,nb1j),1,
     '              XE(1,1),1)
                  XG1(2,4,ng)=DDOT(NST(nb1jp),PG(1,4,ng,nb1j),1,
     '              XE(1,2),1)
                  XG1(3,4,ng)=DDOT(NST(nb1jp),PG(1,4,ng,nb1j),1,
     '              XE(1,3),1)
                ENDIF
C***            RG(ng) stores the Jacobian for a length or area
C***            integral at Gauss point ng.
C***            DXXI(nj,ni) contains the values of dXj/dPSIi at the
C***            Gauss point.
c cpb 9/7/95 unrolling loops
C                DO ni=1,NITB
C                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
C                    DXXI(nj,ni)=XG1(nj,NU1(ni),ng)
C                  ENDDO !nj
C                ENDDO !ni
                DXXI(1,1)=XG1(1,2,ng)
                DXXI(2,1)=XG1(2,2,ng)
                DXXI(3,1)=XG1(3,2,ng)
                DXXI(1,2)=XG1(1,4,ng)
                DXXI(2,2)=XG1(2,4,ng)
                DXXI(3,2)=XG1(3,4,ng)
                RG(ng)=DSQRT(
     '            (DXXI(2,1)*DXXI(3,2)-DXXI(2,2)*DXXI(3,1))*
     '            (DXXI(2,1)*DXXI(3,2)-DXXI(2,2)*DXXI(3,1))+
     '            (DXXI(1,2)*DXXI(3,1)-DXXI(1,1)*DXXI(3,2))*
     '            (DXXI(1,2)*DXXI(3,1)-DXXI(1,1)*DXXI(3,2))+
     '            (DXXI(1,1)*DXXI(2,2)-DXXI(1,2)*DXXI(2,1))*
     '            (DXXI(1,1)*DXXI(2,2)-DXXI(1,2)*DXXI(2,1)))
C*** Find the unit outward normal.
c cpb 9/7/95 Do explicit normal for rectangular cartesian
C                CALL NORMAL(ne,nr,NW,XG1(1,1,ng),XN(1,ng),
C     '            INTERFACE,ERROR,*9999)
                XN(1,ng)=XG1(2,2,ng)*XG1(3,4,ng)-XG1(3,2,ng)*
     '            XG1(2,4,ng)
                XN(2,ng)=XG1(3,2,ng)*XG1(1,4,ng)-XG1(1,2,ng)*
     '            XG1(3,4,ng)
                XN(3,ng)=XG1(1,2,ng)*XG1(2,4,ng)-XG1(2,2,ng)*
     '            XG1(1,4,ng)
                VLENGTH=DSQRT(XN(1,ng)*XN(1,ng)+XN(2,ng)*XN(2,ng)+
     '            XN(3,ng)*XN(3,ng))
                CALL ASSERT(VLENGTH.GT.RDELTA,'>>Zero length normal',
     '            ERROR,*9999)
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
            ENDIF
            CALC_XR=.TRUE.

C***        Evaluate standard BIE if not a hermite problem or hermite
C***        with normal BIE used.
            IF(GEN_CONVENT_EQUATION) THEN

C***          XEPGKGQ calculates the standard BIE equations for element
C***          ne and singularity at node np.  It uses
C***          ny1=nynp(1,1,nh,1,np,nrc,nc,nr) i.e. it assumes nk=1.
              IF(USEADAPINT) THEN
                CALL XEPGKGQ_3DL(nb1j,nbqh,
     '            nbqhp,nbuh,nbuhp,NDTOT,nh,NHSTGK,NHSTGQ,
     '            NKHE(1,1,1,ne),NKH,
     '            NNMIN,np,NPNE(1,1,ne),nr,CE(1,ne),
     '            CURVCORRECT(1,1,1,ne),DET_ADAPT,DRDN,GKES,
     '            GQES,PG_Q,PG_U,RAD,RG,SE(1,1,ne),WG,XG1,XN,
     '            XPFP,XR,CALC_XR,ERROR,*9999)
C SMAR009 18/01/99 removed GK,GQ,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,nx,NYNP,
              ELSE
                CALL XEPGKGQ_3DL(nb1j,nbqh,
     '            nbqhp,nbuh,nbuhp,NDTOT,nh,NHSTGK,NHSTGQ,
     '            NKHE(1,1,1,ne),NKH,
     '            NNMIN,np,NPNE(1,1,ne),nr,CE(1,ne),
     '            CURVCORRECT(1,1,1,ne),DET(1,0,1,nsplit+IADAPTIVE),
     '            DRDN,GKES,GQES,PG(1,1,1,nbqh),PG(1,1,1,nbuh),
     '            RAD,RG,SE(1,1,ne),WG,XG1,XN,XPFP,XR,CALC_XR,
     '            ERROR,*9999)
C SMAR009 18/01/99 removed GK,GQ,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,nx,NYNP,
              ENDIF
            ENDIF

C***        Evaluate hypersingular BIE's if a hermite problem.

            IF(GEN_DERIV_EQUATION) THEN

              IF(intscheme.NE.1) THEN !not element splitting

C***            The routine XEPGKGQHYP evaluates the derivative BIE(s)
C***            for the other values of nk (nk=2 only for 2D case) when
C***            np is not in ne.
                IF(USEADAPINT) THEN
                  CALL XEPGKGQHYP_3DL(nb1j,nb1jp,nbqh,nbqhp,
     '              nbuh,nbuhp,NDTOT,ne,nh,
     '              NHSTGK,NHSTGQ,NKHE,NKH,np,NPNE(1,1,ne),nr,
     '              CE(1,ne),CURVCORRECT(1,1,1,ne),DET_ADAPT,
     '              DRDNO,DSDX,GKES,GQES,PG_Q,PG_U,RAD,RG,SE,
     '              WG,XG1,XN,XNO,XN_GRAD,XPFP,XR,XR_GRAD,CALC_XR,
     '              ERROR,*9999)
C SMAR009 19/01/99 removed GK,GQ,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,nx,NYNP,
                ELSE
                  CALL XEPGKGQHYP_3DL(nb1j,nb1jp,nbqh,nbqhp,
     '              nbuh,nbuhp,NDTOT,ne,nh,
     '              NHSTGK,NHSTGQ,NKHE,NKH,np,NPNE(1,1,ne),nr,
     '              CE(1,ne),CURVCORRECT(1,1,1,ne),
     '              DET(1,0,1,nsplit+IADAPTIVE),DRDNO,DSDX,GKES,
     '              GQES,PG(1,1,1,nbqh),PG(1,1,1,nbuh),
     '              RAD,RG,SE,WG,XG1,XN,XNO,XN_GRAD,XPFP,XR,XR_GRAD,
     '              CALC_XR,ERROR,*9999)
C SMAR009 19/01/99 removed GK,GQ,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,nx,NYNP,
                ENDIF

              ELSE !np is contained in element ne (local node NNMIN)

C***            XEPGKGQHYPS is used when np belongs to ne. If ktyp92=3
C***            then derivative equations are calculated for each value
C***            of nk.
                IF(USEADAPINT) THEN
                  CALL XEPGKGQHYPS_3DL(nb1j,nb1jp,nbqh,nbqhp,
     '              nbuh,nbuhp,NDTOT,nh,
     '              NHSTGK,NHSTGQ,NKHE(1,1,1,ne),NKH,NNMIN,np,
     '              NPNE(1,1,ne),nr,CE(1,ne),
     '              CURVCORRECT(1,1,1,ne),DET_ADAPT,DRDNO,DSDX,
     '              GKES,GQES,PG_Q,PG_U,RAD,RG,SE(1,1,ne),WG,
     '              XG1,XN,XNO,XN_GRAD,XPFP,XR,XR_GRAD,CALC_XR,
     '              ERROR,*9999)
C SMAR009 19/01/99 removed GK,GQ,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,nx,NYNP,
                ELSE
                  CALL XEPGKGQHYPS_3DL(nb1j,nb1jp,nbqh,nbqhp,
     '              nbuh,nbuhp,NDTOT,nh,
     '              NHSTGK,NHSTGQ,NKHE(1,1,1,ne),NKH,NNMIN,np,
     '              NPNE(1,1,ne),nr,CE(1,ne),
     '              CURVCORRECT(1,1,1,ne),DET(1,0,1,nsplit+IADAPTIVE),
     '              DRDNO,DSDX,GKES,GQES,PG(1,1,1,nbqh),
     '              PG(1,1,1,nbuh),RAD,RG,SE(1,1,ne),WG,XG1,XN,XNO,
     '              XN_GRAD,XPFP,XR,XR_GRAD,CALC_XR,ERROR,*9999)
C SMAR009 19/01/99 removed GK,GQ,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,nx,NYNP,
                ENDIF
              ENDIF
            ENDIF
            ISP=ISP+nsplit
          ENDDO !End of nsplit loop

C!!! cpb 20/10/98 Should not need locking for XPGKGQ as each node
C!!! *SHOULD* be data independent in GK and GQ.

C CPB 30/10/98 Adding direct assembly of solution matrices

C!!! Note: No locking will be used in order to retain parallel
C!!! scalability. This *WILL* be a problem if multiple global equations
C!!! are mapped to the same solution equation.

          IF(CALC_GLOBAL(nr,nx)) THEN
            DO nkk=1,NDTOT+1
              DO nhs=-NDTOT-1,-1
                nz=LGKE(nhs,nkk,2)
                IF(nz.NE.0) GK(nz)=GK(nz)+GKES(nhs,nkk)
              ENDDO !nhs
              DO nhs=1,NHSTGK
                nz=LGKE(nhs,nkk,2)
                IF(nz.NE.0) GK(nz)=GK(nz)+GKES(nhs,nkk)
              ENDDO !nhs
              DO nhs=-NDTOT-1,-1
                nz=LGQE(nhs,nkk,2)
                IF(nz.NE.0) GQ(nz)=GQ(nz)+GQES(nhs,nkk)
              ENDDO !nhs
              DO nhs=1,NHSTGQ
                nz=LGQE(nhs,nkk,2)
                IF(nz.NE.0) GQ(nz)=GQ(nz)+GQES(nhs,nkk)
              ENDDO !nhs
            ENDDO !nkk
          ELSE
            DO nkk=1,NDTOT+1
              ny1=LGKE(0,nkk,1)
              ny3=GETNYR(1,NPNY,nr,0,1,ny1,NYNE,NYNP) !Equiv glob var.
              DO noy1=1,NONY(0,ny1,1)
                no1=NONY(noy1,ny1,1)
                co1=CONY(noy1,ny1,1)
                DO nhs=-NDTOT-1,NHSTGK
                  IF(nhs.NE.0) THEN
                    ny2=LGKE(nhs,nkk,1)
                    IF(ny2.NE.0) THEN
                      DO noy2=1,NONY(0,ny2,2)
                        no2=NONY(noy2,ny2,2)
                        co2=CONY(noy2,ny2,2)
                        CALL SPARSE(no1,no2,NOT(1,1,nr_gkk,nx),nzz,
     '                    NZ_GKK_M,NZZT(1,nr_gkk,nx),ISC_GKK,ISR_GKK,
     '                    SPARSEGKK(nx),ERROR,*9999)
                        IF(nzz.NE.0) GKK(nzz)=GKK(nzz)+
     '                    GKES(nhs,nkk)*co1*co2
                      ENDDO !noy2
                      IF(NONY(0,ny2,2).EQ.0.OR.NONY(0,ny3,2).EQ.0) THEN
                        nz=LGKE(nhs,nkk,2)
                        IF(nz.NE.0) GK(nz)=GK(nz)+GKES(nhs,nkk)
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO !nhs
              ENDDO !noy1
              IF(NONY(0,ny1,1).EQ.0) THEN
                DO nhs=-NDTOT-1,-1
                  nz=LGKE(nhs,nkk,2)
                  IF(nz.NE.0) GK(nz)=GK(nz)+GKES(nhs,nkk)
                ENDDO !nhs
                DO nhs=1,NHSTGK
                  nz=LGKE(nhs,nkk,2)
                  IF(nz.NE.0) GK(nz)=GK(nz)+GKES(nhs,nkk)
                ENDDO !nhs
              ENDIF
              ny1=LGQE(0,nkk,1)
              ny3=GETNYR(2,NPNY,nr,0,1,ny1,NYNE,NYNP) !Equiv glob var.
              DO noy1=1,NONY(0,ny1,1)
                no1=NONY(noy1,ny1,1)
                co1=CONY(noy1,ny1,1)
                DO nhs=-NDTOT-1,NHSTGQ
                  IF(nhs.NE.0) THEN
                    ny2=LGQE(nhs,nkk,1)
                    IF(ny2.NE.0) THEN
                      DO noy2=1,NONY(0,ny2,2)
                        no2=NONY(noy2,ny2,2)
                        co2=CONY(noy2,ny2,2)
                        CALL SPARSE(no1,no2,NOT(1,1,nr_gkk,nx),nzz,
     '                    NZ_GKK_M,NZZT(1,nr_gkk,nx),ISC_GKK,ISR_GKK,
     '                    SPARSEGKK(nx),ERROR,*9999)
                        IF(nzz.NE.0) GKK(nzz)=GKK(nzz)-
     '                    GQES(nhs,nkk)*co1*co2
                      ENDDO !noy2
                      IF(NONY(0,ny2,2).EQ.0.OR.NONY(0,ny3,2).EQ.0) THEN
                        nz=LGQE(nhs,nkk,2)
                        IF(nz.NE.0) GQ(nz)=GQ(nz)+GQES(nhs,nkk)
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO !nhs
              ENDDO !noy1
              IF(NONY(0,ny1,1).EQ.0) THEN
                DO nhs=-NDTOT-1,-1
                  nz=LGQE(nhs,nkk,2)
                  IF(nz.NE.0) GQ(nz)=GQ(nz)+GQES(nhs,nkk)
                ENDDO !nhs
                DO nhs=1,NHSTGQ
                  nz=LGQE(nhs,nkk,2)
                  IF(nz.NE.0) GQ(nz)=GQ(nz)+GQES(nhs,nkk)
                ENDDO !nhs
              ENDIF
            ENDDO !nkk
          ENDIF

        ENDDO !End of nh loop
      ENDDO !End of ne loop

      CALL EXITS('XPGKGQ_3DL_BUFFER')
      RETURN
 9999 CALL ERRORS('XPGKGQ_3DL_BUFFER',ERROR)
      CALL EXITS('XPGKGQ_3DL_BUFFER')
      RETURN 1
      END


