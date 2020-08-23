      SUBROUTINE XPGKGQ(IBT,IDO,INP,ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,
     '  ISR_GKK,ISR_GQ,LGKE,LGQE,NBH,NBJ,NDET,NEELEM,NENP,NGAP,NHE,
     '  NHP,NKH,NKHE,NKJE,NLL,NONY,np,NP_INTERFACE,NPF,NPNE,NPNY,nr,
     '  nr_gkk,NRE,NVJE,NW,nx,NYNE,NYNP,CE,CONY,CURVCORRECT,DET,
     '  DET_ADAPT,DL,DRDN,DRDNO,GK,GKK,GKES,GQ,GQES,PG,PG_J,PG_Q,
     '  PG_U,RAD,RD,RG,SE,WG,XA,XE,XG,XG1,XIG,XIG_J,XIG_Q,XIG_U,
     '  XN,XN_GRAD,XP,XR,XR_GRAD,INTERFACE,FULL_EQUATIONS,ERROR,*)

C#### Subroutine: XPGKGQ
C###  Description:
C###    XPGKGQ calculates the contributions to global stiffness matrix
C###    from the boundary integral equations for a singularity at node
C###    np.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'jtyp00.cmn'
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
     '  NGAP(NIM,NBM),NHE(NEM),
     '  NHP(NPM),NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NKH(NHM,NPM,NCM),NLL(12,NEM),NONY(0:NOYM,NYM,NRCM),np,
     '  NP_INTERFACE(0:NPM,0:3),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNY(0:6,NYM,0:NRCM),nr,nr_gkk,NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),CONY(0:NOYM,NYM,NRCM),CURVCORRECT(2,2,NNM,NEM),
     '  DET(NBFM,0:NNM,NGM,6),DET_ADAPT(NBFM,0:NNM,NGM),DL(3,NLM),
     '  DRDN(NGM),DRDNO(NGM,NKM),GK(NZ_GK_M),GKK(NZ_GKK_M),
     '  GKES(-NKM:NHM*NSM,NKM),GQ(NZ_GQ_M),GQES(-NKM:NHM*NSM,NKM),
     '  PG(NSM,NUM,NGM,NBM),PG_J(NSM,NUM,NGM),PG_Q(NSM,NUM,NGM),
     '  PG_U(NSM,NUM,NGM),RAD(NGM),RD(NGM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XG1(NJM,NUM,NGM),XIG(NIM,NGM,NBM),XIG_J(NIM,NGM),
     '  XIG_Q(NIM,NGM),XIG_U(NIM,NGM),XN(NJM,NGM),XN_GRAD(NJM,NGM),
     '  XP(NKM,NVM,NJM,NPM),XR(NJM,NGM),XR_GRAD(NJM,NGM)
      CHARACTER ERROR*(*)
      LOGICAL INTERFACE,FULL_EQUATIONS
!     Local Variables
      INTEGER GETNYR,GREENFLAG,IADAPTIVE,intscheme,ISP,NB2,nb1j,
     '  nb1jp,nbbem,nbqh,nbqhp,nbuh,nbuhp,NDTOT,ne,ng,nh,nhs,
     '  NHSTGK,NHSTGQ,nhx,ni,NITB,nj,nkk,NNMIN,NNSP,
     '  no1,no2,noy1,noy2,noelem,ns,NSPLIT,NU1(0:3),numtypes,ny1,
     '  ny2,ny3,nz,nzz,TYPE
      REAL*8 co1,co2,D(2),DSDX(3,3),DXDS(3,3),
     '  DXXI(3,3),MINDIST,SUMXG,XIMIN(2),XNO(3,4),XPFP(3)
      LOGICAL GEN_CONVENT_EQUATION,GEN_DERIV_EQUATION,CALC_XR,
     '  USEADAPINT,USELOGBEM

      DATA NU1/1,2,4,7/

      CALL ENTERS('XPGKGQ',*9999)

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(XPGKGQ_1)
        WRITE(OP_STRING,'(/'' Node number '',I5)') np
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(XPGKGQ_1)
      ENDIF

C***  Determine location of the singularity
      DO nj=1,NJT
        XPFP(nj)=XP(1,1,nj,np)
      ENDDO !nj

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
CC$OMP CRITICAL(XPGKGQ_2)
          WRITE(OP_STRING,'(/'' Element number '',I5)') ne
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(XPGKGQ_2)
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
CC$OMP CRITICAL(XPGKGQ_3)
          WRITE(OP_STRING,'(/'' Integration scheme '',I1,'', nbbem '','
     '      //'I2)') intscheme,nbbem
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(XPGKGQ_3)
        ENDIF

        IF(intscheme.NE.1) THEN
          NNMIN=0
        ENDIF

        nb1jp=NFBASE(1,NBASEF(NBJ(1,ne),nbbem))
        NITB=NIT(nb1jp)

C***    Loop over dependent variables
        DO nhx=1,NHP(np)
          nh=NH_LOC(nhx,nx)

C***      Decide basis functions from the parent family basis function
C***      for dependent (u) and normal derivative (q).
          nbuhp=NFBASE(1,NBASEF(NBH(nh,1,ne),nbbem))
          nbqhp=NFBASE(1,NBASEF(NBH(nh,2,ne),nbbem))

          IF(BEMLOGGAUSS.AND.intscheme.EQ.1) THEN
            USELOGBEM=.TRUE.
            NUMTYPES=2
          ELSE
            USELOGBEM=.FALSE.
            NUMTYPES=1
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
     '      NHE(ne),NKHE(1,1,nh,ne),NKH(1,1,1),np,NPNE(1,nbuhp,ne),
     '      NPNY,nr,nx,NYNE,NYNP,NZ_GK_M,ERROR,*9999)
          CALL MELGEB(ISC_GQ,ISR_GQ,LGQE,nbqhp,2,NDTOT,nh,NHSTGQ,
     '      NHE(ne),NKHE(1,1,nh,ne),NKH(1,1,2),np,NPNE(1,nbqhp,ne),
     '      NPNY,nr,nx,NYNE,NYNP,NZ_GQ_M,ERROR,*9999)

          DO nkk=1,NDTOT+1
            DO nhs=-NDTOT-1,NHSTGK
              GKES(nhs,nkk)=0.0d0
            ENDDO !nhs
            DO nhs=-NDTOT-1,NHSTGQ
              GQES(nhs,nkk)=0.0d0
            ENDDO !nhs
          ENDDO !nkk

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

C cpb 5/11/96 Adding log Gauss quad. Most of the time NUM_TYPES will be
C 1 except for when you are doing 2D laplace in which case you have
C to integrate the convential equation in two parts 1 for the GK and
C 1 for the GQ with logarithmic quadrature.

          DO type=1,NUMTYPES
C***        Calculate the element integrals for the current node and
C***        dependent variable.
            ISP=0
            IF(intscheme.EQ.1) THEN !Element splitting

C***          Find number of split elements used at nodes less than
C***          NNMIN
              DO nnsp=0,NNMIN-1
                ISP=ISP+NDET(nb1jp,nnsp)
              ENDDO
            ENDIF

C***        Loop over number of split elements (if any)
            DO nsplit=1,NDET(nb1jp,NNMIN)

              IF(USELOGBEM) THEN
                IF(type.EQ.1) THEN
                  GREENFLAG=1
                  nb2=nbbem+ISP
                ELSE IF(type.EQ.2) THEN
                  GREENFLAG=2
                  nb2=nbbem+ISP+NNT(nb1jp)
                ENDIF
              ELSE
                GREENFLAG=0
                nb2=nbbem+ISP
              ENDIF

C***          Determine basis function for the dependent variable for
C***          the split element.
C***          nbuh,nbqh and nb1j are the global basis functions with
C***          the appropriate quadrature for the dependent (H), normal
C***          derivative (Q) and geometric (J) variables.

              nbuh=NBASEF(NBH(nh,1,ne),nb2)
              nbqh=NBASEF(NBH(nh,2,ne),nb2)
C cpb 2/11/98 Shouldn't use nh for NBJ. Use nj=1 instead for now.
C              nb1j=NBASEF(NBJ(nh,ne),nb2)
              nb1j=NBASEF(NBJ(1,ne),nb2)

C***          If element spliting/adaptive integration or the first node
C***          in the integration scheme set up the gauss point dependent
C***          arrays.

              IF(intscheme.EQ.1.OR.intscheme.EQ.2.OR..TRUE.) THEN
                IF(intscheme.EQ.1.OR.intscheme.EQ.2) CALC_XR=.TRUE.
                DO ng=1,NGT(nb1j)
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    IF(USEADAPINT) THEN
                      DO ni=0,NITB
                        SUMXG=0.0d0
                        DO ns=1,NST(nb1jp)
                          SUMXG=SUMXG+PG_J(ns,NU1(ni),ng)*XE(ns,nj)
                        ENDDO !ns
                        XG1(nj,NU1(ni),ng)=SUMXG
                      ENDDO !ni
                    ELSE
                      DO ni=0,NITB
                        SUMXG=0.0d0
                        DO ns=1,NST(nb1jp)
                          SUMXG=SUMXG+PG(ns,NU1(ni),ng,nb1j)*XE(ns,nj)
                        ENDDO !ns
                        XG1(nj,NU1(ni),ng)=SUMXG
                      ENDDO !ni
                    ENDIF
                  ENDDO !nj

C***              RG(ng) stores the Jacobian for a length or area
C***              integral at Gauss point ng.
C***              DXXI(nj,ni) contains the values of dXj/dPSIi at the
C***              Gauss point.

                  DO ni=1,NITB
                    DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                      DXXI(nj,ni)=XG1(nj,NU1(ni),ng)
                    ENDDO !nj
                  ENDDO !ni
                  IF(NITB.EQ.1) THEN
C***                Jacobian for 1D integral in 2D space
                    RG(ng)=DSQRT(DXXI(1,1)*DXXI(1,1)+DXXI(2,1)*
     '                DXXI(2,1))
                  ELSE !Calculate cross product of 2 tangent vectors
C***                Jacobian for 2D integral in 3D space
                    RG(ng)=DSQRT(
     '                (DXXI(2,1)*DXXI(3,2)-DXXI(2,2)*DXXI(3,1))*
     '                (DXXI(2,1)*DXXI(3,2)-DXXI(2,2)*DXXI(3,1))+
     '                (DXXI(1,2)*DXXI(3,1)-DXXI(1,1)*DXXI(3,2))*
     '                (DXXI(1,2)*DXXI(3,1)-DXXI(1,1)*DXXI(3,2))+
     '                (DXXI(1,1)*DXXI(2,2)-DXXI(1,2)*DXXI(2,1))*
     '                (DXXI(1,1)*DXXI(2,2)-DXXI(1,2)*DXXI(2,1)))
                  ENDIF
C***              Find the unit outward normal.
                  CALL NORMAL(ne,nr,NW,XG1(1,1,ng),XN(1,ng),
     '              INTERFACE,ERROR,*9999)
                  IF(JTYP4.GT.1) THEN
                    IF(JTYP4.EQ.2) THEN
                      RD(ng)=XG1(2,1,ng) !rotat symm about x/r axis
                    ELSE IF(JTYP4.EQ.3) THEN
                      RD(ng)=XG1(1,1,ng) !rotat symm about y/z axis
                    ENDIF
                  ENDIF
                ENDDO !End ng loop
              ENDIF
              CALC_XR=.TRUE.

C***          Evaluate standard BIE if not a hermite problem or hermite
C***          with normal BIE used.
              IF(GEN_CONVENT_EQUATION) THEN

C***            XEPGKGQ calculates the standard BIE equations for
C***            element ne and singularity at node np.  It uses
C***            ny1=nynp(1,1,nh,1,np,nrc,nc,nr) i.e. it assumes nk=1.

                IF(USEADAPINT) THEN
                  CALL XEPGKGQ(GREENFLAG,
     '              nb1j,nb1jp,nbqh,nbqhp,nbuh,nbuhp,NDTOT,nh,NHE(ne),
     '              NHSTGK,NHSTGQ,NKHE(1,1,1,ne),NKH,NNMIN,np,
     '              NPNE(1,1,ne),nr,nx,CE(1,ne),
     '              CURVCORRECT(1,1,1,ne),DET_ADAPT,DRDN,GKES,
     '              GQES,PG_Q,PG_U,RAD,RD,RG,
     '              SE(1,1,ne),WG,XIG,XG1,XN,XPFP,XR,CALC_XR,
     '              ERROR,*9999)
C  SMAR009 18/01/99 removed GK,GQ,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,NYNP,
                ELSE
                  CALL XEPGKGQ(GREENFLAG,
     '              nb1j,nb1jp,nbqh,nbqhp,nbuh,nbuhp,NDTOT,nh,NHE(ne),
     '              NHSTGK,NHSTGQ,NKHE(1,1,1,ne),NKH,NNMIN,np,
     '              NPNE(1,1,ne),nr,nx,CE(1,ne),
     '              CURVCORRECT(1,1,1,ne),DET(1,0,1,nsplit+IADAPTIVE),
     '              DRDN,GKES,GQES,PG(1,1,1,nbqh),PG(1,1,1,nbuh),
     '              RAD,RD,RG,SE(1,1,ne),WG,XIG,XG1,XN,XPFP,XR,CALC_XR,
     '              ERROR,*9999)
C  SMAR009 18/01/99 removed GK,GQ,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,NYNP,
                ENDIF

              ENDIF

C***          Evaluate hypersingular BIE's if a hermite problem.

              IF(GEN_DERIV_EQUATION.AND.TYPE.NE.2) THEN

                IF(intscheme.NE.1) THEN !not element splitting

C***               The routine XEPGKGQHYP evaluates the derivative
C***               BIE(s) for the other values of nk (nk=2 only for
C***               2D case) when np is not in ne.

                  IF(USEADAPINT) THEN
                    CALL XEPGKGQHYP(nb1j,nb1jp,nbqh,nbqhp,
     '                nbuh,nbuhp,NDTOT,ne,
     '                nh,NHSTGK,NHSTGQ,NKHE,NKH,NNMIN,np,NPNE(1,1,ne),
     '                nr,CE(1,ne),CURVCORRECT(1,1,1,ne),
     '                DET_ADAPT,DRDNO,DSDX,GKES,GQES,PG_Q,
     '                PG_U,RAD,RG,SE,WG,XG1,XN,XNO,
     '                XN_GRAD,XPFP,XR,XR_GRAD,CALC_XR,ERROR,*9999)
C SMAR009 19/01/99 removed GK,GQ,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,nx,NYNP,
                  ELSE
                    CALL XEPGKGQHYP(nb1j,nb1jp,nbqh,nbqhp,nbuh,
     '                nbuhp,NDTOT,ne,nh,
     '                NHSTGK,NHSTGQ,NKHE,NKH,NNMIN,np,NPNE(1,1,ne),nr,
     '                CE(1,ne),CURVCORRECT(1,1,1,ne),
     '                DET(1,0,1,nsplit+IADAPTIVE),DRDNO,DSDX,
     '                GKES,GQES,PG(1,1,1,nbqh),PG(1,1,1,nbuh),
     '                RAD,RG,SE,WG,XG1,XN,XNO,XN_GRAD,XPFP,XR,XR_GRAD,
     '                CALC_XR,ERROR,*9999)
C SMAR009 19/01/99 removed GK,GQ,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,nx,NYNP,
                  ENDIF

                ELSE !np is contained in element ne (local node NNMIN)

C***              XEPGKGQHYPS is used when np belongs to ne. If
C***              ktyp92=3 then derivative equations are calculated for
C***              each value of nk.

                  IF(USEADAPINT) THEN
                    CALL XEPGKGQHYPS(nb1j,nb1jp,nbqh,nbqhp,
     '                nbuh,nbuhp,NDTOT,nh,
     '                NHSTGK,NHSTGQ,NKHE(1,1,1,ne),NKH,NNMIN,np,
     '                NPNE(1,1,ne),nr,CE(1,ne),
     '                CURVCORRECT(1,1,1,ne),DET_ADAPT,DRDNO,DSDX,
     '                GKES,GQES,PG_Q,PG_U,RAD,RG,SE(1,1,ne),
     '                WG,XG1,XN,XNO,XN_GRAD,XPFP,XR,XR_GRAD,CALC_XR,
     '                ERROR,*9999)
C SMAR009 19/01/99 removed GK,GQ,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,nx,NYNP,
                  ELSE
                    CALL XEPGKGQHYPS(nb1j,nb1jp,nbqh,nbqhp,
     '                nbuh,nbuhp,NDTOT,nh,
     '                NHSTGK,NHSTGQ,NKHE(1,1,1,ne),NKH,NNMIN,np,
     '                NPNE(1,1,ne),nr,CE(1,ne),
     '                CURVCORRECT(1,1,1,ne),
     '                DET(1,0,1,nsplit+IADAPTIVE),DRDNO,DSDX,
     '                GKES,GQES,PG(1,1,1,nbqh),PG(1,1,1,nbuh),
     '                RAD,RG,SE(1,1,ne),WG,XG1,XN,XNO,XN_GRAD,XPFP,XR,
     '                XR_GRAD,CALC_XR,ERROR,*9999)
C SMAR009 19/01/99 removed GK,GQ,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,nx,NYNP,
                  ENDIF
                ENDIF
              ENDIF
              ISP=ISP+nsplit
            ENDDO !End of nsplit loop
          ENDDO !num_types loop
C cpb 2/1/98 Adding BEM element stiffness matrices

C!!! cpb 20/10/98 Should not need locking for XPGKGQ as each node
C!!! *SHOULD* be data independent in GK and GQ.

C CPB 30/10/98 Adding direct assembly of solution matrices

C!!! Note: No locking will be used in order to retain parallel
C!!! scalability. This *WILL* be a problem if multiple global equations
C!!! are mapped to the same solution equation.

          IF(NDTOT.EQ.0) THEN
            IF(CALC_GLOBAL(nr,nx)) THEN
              DO nhs=1,NHSTGK
                nz=LGKE(nhs,1,2)
                IF(nz.NE.0) GK(nz)=GK(nz)+GKES(nhs,1)
              ENDDO !nhs
              DO nhs=1,NHSTGQ
                nz=LGQE(nhs,1,2)
                IF(nz.NE.0) GQ(nz)=GQ(nz)+GQES(nhs,1)
              ENDDO !nhs
            ELSE
              ny1=LGKE(0,1,1)
              ny3=GETNYR(1,NPNY,nr,0,1,ny1,NYNE,NYNP) !Equiv glob var.
              DO noy1=1,NONY(0,ny1,1)
                no1=NONY(noy1,ny1,1)
                co1=CONY(noy1,ny1,1)
                DO nhs=1,NHSTGK
                  ny2=LGKE(nhs,1,1)
                  IF(ny2.NE.0) THEN
                    DO noy2=1,NONY(0,ny2,2)
                      no2=NONY(noy2,ny2,2)
                      co2=CONY(noy2,ny2,2)
                      CALL SPARSE(no1,no2,NOT(1,1,nr_gkk,nx),nzz,
     '                  NZ_GKK_M,NZZT(1,nr_gkk,nx),ISC_GKK,ISR_GKK,
     '                  SPARSEGKK(nx),ERROR,*9999)
                      IF(nzz.NE.0) GKK(nzz)=GKK(nzz)+
     '                  GKES(nhs,1)*co1*co2
                    ENDDO !noy2
                    IF(NONY(0,ny2,2).EQ.0.OR.NONY(0,ny3,2).EQ.0) THEN
                      nz=LGKE(nhs,1,2)
                      IF(nz.NE.0) GK(nz)=GK(nz)+GKES(nhs,1)
                    ENDIF
                  ENDIF
                ENDDO !nhs
              ENDDO !noy1
              IF(NONY(0,ny1,1).EQ.0) THEN
                DO nhs=1,NHSTGK
                  nz=LGKE(nhs,1,2)
                  IF(nz.NE.0) GK(nz)=GK(nz)+GKES(nhs,1)
                ENDDO !nhs
              ENDIF
              ny1=LGQE(0,1,1)
              ny3=GETNYR(2,NPNY,nr,0,1,ny1,NYNE,NYNP) !Equiv glob var.
              DO noy1=1,NONY(0,ny1,1)
                no1=NONY(noy1,ny1,1)
                co1=CONY(noy1,ny1,1)
                DO nhs=1,NHSTGQ
                  ny2=LGQE(nhs,1,1)
                  IF(ny2.NE.0) THEN
                    DO noy2=1,NONY(0,ny2,2)
                      no2=NONY(noy2,ny2,2)
                      co2=CONY(noy2,ny2,2)
                      CALL SPARSE(no1,no2,NOT(1,1,nr_gkk,nx),nzz,
     '                  NZ_GKK_M,NZZT(1,nr_gkk,nx),ISC_GKK,ISR_GKK,
     '                  SPARSEGKK(nx),ERROR,*9999)
                      IF(nzz.NE.0) GKK(nzz)=GKK(nzz)-
     '                  GQES(nhs,1)*co1*co2
                    ENDDO !noy2
                    IF(NONY(0,ny2,2).EQ.0.OR.NONY(0,ny3,2).EQ.0) THEN
                      nz=LGQE(nhs,1,2)
                      IF(nz.NE.0) GQ(nz)=GQ(nz)+GQES(nhs,1)
                    ENDIF
                  ENDIF
                ENDDO !nhs
              ENDDO !noy1
              IF(NONY(0,ny1,1).EQ.0) THEN
                DO nhs=1,NHSTGQ
                  nz=LGQE(nhs,1,2)
                  IF(nz.NE.0) GQ(nz)=GQ(nz)+GQES(nhs,1)
                ENDDO !nhs
              ENDIF
            ENDIF
          ELSE
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
                ny3=GETNYR(1,NPNY,nr,0,1,ny1,NYNE,NYNP) !Equiv glo var
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
     '                      NZ_GKK_M,NZZT(1,nr_gkk,nx),ISC_GKK,ISR_GKK,
     '                      SPARSEGKK(nx),ERROR,*9999)
                          IF(nzz.NE.0) GKK(nzz)=GKK(nzz)+
     '                      GKES(nhs,nkk)*co1*co2
                        ENDDO !noy2
                        IF(NONY(0,ny2,2).EQ.0.OR.NONY(0,ny3,2).EQ.0)
     '                    THEN
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
     '                      NZ_GKK_M,NZZT(1,nr_gkk,nx),ISC_GKK,ISR_GKK,
     '                      SPARSEGKK(nx),ERROR,*9999)
                          IF(nzz.NE.0) GKK(nzz)=GKK(nzz)-
     '                      GQES(nhs,nkk)*co1*co2
                        ENDDO !noy2
                        IF(NONY(0,ny2,2).EQ.0.OR.NONY(0,ny3,2).EQ.0)
     '                    THEN
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
          ENDIF
        ENDDO !End of nh loop
      ENDDO !End of ne loop

      CALL EXITS('XPGKGQ')
      RETURN
 9999 CALL ERRORS('XPGKGQ',ERROR)
      CALL EXITS('XPGKGQ')
      RETURN 1
      END


