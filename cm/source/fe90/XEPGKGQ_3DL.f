      SUBROUTINE XEPGKGQ_3DL(nb1j,nbqh,
     '  nbqhp,nbuh,nbuhp,NDTOT,nh,NHSTGK,NHSTGQ,NKHE,NKH,NNMIN,np,
     '  NPNE,nr,CE,CURVCORRECT,DET,DRDN,GKES,GQES,
     '  PG_Q,PG_U,RAD,RG,SE,WG,XG1,XN,XPFP,XR,CALC_XR,ERROR,*)
C SMAR009 18/01/99 removed GK,GQ,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,nx,NYNP,
C#### Subroutine: XEPGKGQ_3DL
C###  Description:
C###    XEPGKGQ_3DL is an efficient version of XEPGKGQ for 3D Laplace
C###    equation in rectangular cartesian coordinates.
C###  See-Also: XEPGKGQ

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER
     '  nb1j,nbqh,nbqhp,nbuh,nbuhp,NDTOT,nh,
     '  NHSTGK,NHSTGQ,NKHE(NKM,NNM,NHM),NKH(NHM,NPM,NCM),
     '  NNMIN,np,NPNE(NNM,NBFM),nr
C SMAR009 18/01/99 removed ISC_GK(NISC_GKM),ISC_GQ(NISC_GQM),
C  ISR_GK(NISR_GKM),ISR_GQ(NISR_GQM),nx,
C  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM),CURVCORRECT(2,2,NNM),DET(NBFM,0:NNM,NGM),
     '  DRDN(NGM),GKES(-NKM:NHM*NSM,NKM),
     '  GQES(-NKM:NHM*NSM,NKM),PG_Q(NSM,NUM,NGM),
     '  PG_U(NSM,NUM,NGM),RAD(NGM),RG(NGM),SE(NSM,NBFM),
     '  XG1(NJM,NUM,NGM),XN(NJM,NGM),XPFP(*),XR(NJM,NGM),WG(NGM,NBM)
C SMAR009 18/01/99 removed GK(NZ_GK_M),GQ(NZ_GQ_M),
      CHARACTER ERROR*(*)
      LOGICAL CALC_XR
!     Local Variables
      INTEGER MK,ng,nhs,nk1,nkk,nn,npnn,ns,nsnn
       !SMAR009 22/12/98,ny1,ny2,nz
      REAL*8 DGREEN,FACTOR,FACTOR1,FACTOR2,GREEN,SUMU,SUMQ,TERM(999)

      CALL ENTERS('XEPGKGQ_3DL',*9999)

      FACTOR=4.0d0*PI
      IF(IGREN(nr).EQ.8) THEN !Generalised Laplace (3D)
        IF(DABS(CE(1)).GT.ZERO_TOL) THEN
          FACTOR1=FACTOR*CE(1)
        ELSE
          ERROR='>>Error in material property of equation'
          GOTO 9999
        ENDIF
        FACTOR2=1.0d0
      ELSE IF(IGREN(nr).EQ.9.OR.IGREN(nr).EQ.10) THEN
        !Poisson equation with a special source term - for use in
        !constructing uniform-dipole layer transfer matrices
        FACTOR1=1.0d0
        IF(DABS(CE(1)).GT.ZERO_TOL) THEN
          FACTOR2=CE(1)
        ELSE
          ERROR='>>Error in material property of Poisson equation'
          GOTO 9999
        ENDIF
      ELSE
        FACTOR1=FACTOR
        FACTOR2=1.0d0
      ENDIF

C*** Calculate gauss point dependent arrays
      IF(CALC_XR) THEN
        DO ng=1,NGT(nb1j)
c cpb 9/7/95 unrolling loops
C            SUM= 0.0d0
C            SUMR=0.0d0
C            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
C              XR(nj,ng)=XG1(nj,1,ng)-XPFP(nj)
C              SUM=SUM+XR(nj,ng)*XR(nj,ng)
C              SUMR=SUMR+XR(nj,ng)*XN(nj,ng)
C            ENDDO !nj
C            RAD(ng)=DSQRT(SUM)
C            DRDN(ng)=SUMR/RAD(ng)
          XR(1,ng)=XG1(1,1,ng)-XPFP(1)
          XR(2,ng)=XG1(2,1,ng)-XPFP(2)
          XR(3,ng)=XG1(3,1,ng)-XPFP(3)
          RAD(ng)=DSQRT(XR(1,ng)*XR(1,ng)+XR(2,ng)*XR(2,ng)+
     '      XR(3,ng)*XR(3,ng))
          DRDN(ng)=(XR(1,ng)*XN(1,ng)+XR(2,ng)*XN(2,ng)+
     '      XR(3,ng)*XN(3,ng))/RAD(ng)
        ENDDO !ng
        CALC_XR=.TRUE.
      ENDIF

C**** GK matrix components

C cpb 2/1/98 Adding BEM element stiffness matrices
      ns=0
      nhs=0
C      ny1=NYNP(1,1,nh,np,1,1,nr) !eqtn number (row number) of GK.
c cpb 10/7/95 Precalc DGREEN*RG*WG
      DO ng=1,NGT(nbuh)
        DGREEN=-DET(nbuhp,NNMIN,ng)*DRDN(ng)/(RAD(ng)*RAD(ng))
        TERM(ng)=DGREEN*RG(ng)*WG(ng,nbuh)
      ENDDO !ng
      DO nn=1,NNT(nbuhp) !Dependent variable loop
        npnn=NPNE(nn,nbuhp)
        DO nk1=1,NKT(nn,nbuhp)
          MK=NKHE(nk1,nn,nh)  !MK is zero at a distorted node
          ns=ns+1
          nhs=nhs+1
C cpb 3/1/98 I think this should test on NKH(nh,npnn,1)
C          IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,np,1)-KTYP93(1,nr),1)) THEN
          IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,1)-KTYP93(1,nr),1)) THEN
C            ny2=NYNP(MK,1,nh,npnn,2,1,nr) !local column number of GK
            SUMU=0.0d0
            DO ng=1,NGT(nbuh)
C***          Evaluate element integrals
c cpb 10/7/95 Inlining Green's function
C              SUMU=SUMU+PG(ns,1,ng,nbuh)*
C     '          DGREEN(IGREN(nr),0,nh,nh2,CE,DRDN(ng),RAD(ng),
C     '          DET(nbuhp,NNMIN,ng),XN(1,ng),XR(1,ng))*
C     '          RG(ng)*WG(ng,nbuh)
              SUMU=SUMU+TERM(ng)*PG_U(ns,1,ng)
            ENDDO !End of ng loop
C            CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,NZT(1,nx),
C     '        ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C            GK(nz)=GK(nz)+SUMU*SE(ns,nbuhp)/FACTOR
            GKES(nhs,1)=GKES(nhs,1)+SUMU*SE(ns,nbuhp)/FACTOR
          ENDIF !End of MK > 0 loop
        ENDDO !End of nk1 loop
      ENDDO !End of nn loop
C***  GQ matrix components
      ns=0
      nhs=0
C      ny1=NYNP(1,1,nh,np,1,2,nr) !eqtn number (row number) of GQ.
c cpb 10/7/95 Precalc GREEN*RG*WG
      DO ng=1,NGT(nbqh)
        GREEN=DET(nbuhp,NNMIN,ng)/RAD(ng)
        TERM(ng)=GREEN*RG(ng)*WG(ng,nbqh)
      ENDDO !ng
      DO nn=1,NNT(nbqhp) !normal derivative loop
        npnn=NPNE(nn,nbqhp)
        nsnn=ns
        DO nk1=1,NKT(nn,nbqhp)
          MK=NKHE(nk1,nn,nh)  !MK is zero at a distorted node
          ns=ns+1
          nhs=nhs+1
C cpb 3/1/98 I think this should test on NKH(nh,npnn,2)
C          IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,np,2)-KTYP93(2,nr),1)) THEN
          IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,2)-KTYP93(2,nr),1)) THEN
C            ny2=NYNP(MK,1,nh,npnn,2,2,nr) !local column number of GQ
            SUMQ=0.0d0
            DO ng=1,NGT(nbqh)
C***            Evaluate element integrals
c cpb 9/7/95 inlining Greens function
C                SUMQ=SUMQ+PG(ns,1,ng,nbqh)*
C     '            GREEN(IGREN(nr),nh,nh2,CE,RAD(ng),
C     '            DET(nbqhp,NNMIN,ng),XR(1,ng))*
C     '            RG(ng)*WG(ng,nbqh)
              SUMQ=SUMQ+TERM(ng)*PG_Q(ns,1,ng)
            ENDDO !End ng loop
C            CALL SPARSE(ny1,ny2,NYT(1,2,nx),nz,NZ_GQ_M,NZT(2,nx),
C     '        ISC_GQ,ISR_GQ,KTYP24,ERROR,*9999)
C            GQ(nz)=GQ(nz)+SUMQ*SE(ns,nbqhp)/(FACTOR1*FACTOR2)
            GQES(nhs,1)=GQES(nhs,1)+SUMQ*SE(ns,nbqhp)/(FACTOR1*FACTOR2)
C cpb 20/11/96 Adding bem curvature correction
            IF(BEMCURVATURECORRECTION.AND.(nk1.EQ.2.OR.nk1.EQ.3)) THEN
C              ny2=NYNP(2,1,nh,npnn,2,1,nr) !loc col # of GK
C              CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
C     '          NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C              GK(nz)=GK(nz)-SUMQ*SE(ns,nbqhp)/FACTOR1*
C     '          CURVCORRECT(1,nk1-1,nn)
C              ny2=NYNP(3,1,nh,npnn,2,1,nr) !loc col # of GK
C              CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
C     '          NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C              GK(nz)=GK(nz)-SUMQ*SE(ns,nbqhp)/FACTOR1*
C     '          CURVCORRECT(2,nk1-1,nn)
              GKES(nsnn+2,1)=GKES(nsnn+2,1)-SUMQ*SE(ns,nbqhp)/FACTOR1*
     '          CURVCORRECT(1,nk1-1,nn)
              GKES(nsnn+3,1)=GKES(nsnn+3,1)-SUMQ*SE(ns,nbqhp)/FACTOR1*
     '          CURVCORRECT(2,nk1-1,nn)
            ENDIF
         ENDIF !End of MK > 0 loop
        ENDDO !End of nk1 loop
      ENDDO !End of nn loop

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(XEPGKQ_3DL_1)
        WRITE(OP_STRING,'('' Element contribution to'//
     '   ' Global matrices'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Node '',I5)') np
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        !GK matrix
C        ny1=NYNP(1,1,nh,np,1,1,nr) !eqtn number (row number) of GK.
C        DO nn=1,NNT(nbuhp)
C          npnn=NPNE(nn,nbuhp)
C          DO nk1=1,NKT(nn,nbuhp)
C            MK=NKE(nk1,nn,nbuhp)
C            IF(MK.GT.0.AND.
C     '        nk1.LE.MAX(NKH(nh,np,1)-KTYP93(1,nr),1)) THEN
C              ny2=NYNP(MK,1,nh,npnn,2,1,nr) !local col number of GK
C              CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
C     '          NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C              WRITE(OP_STRING,'('' GK('',I5,'','',I5,'')='',D12.4)')
C     '          ny1,ny2,GK(nz)
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C            ENDIF
C          ENDDO !nk1
C        ENDDO !nn
        DO nkk=1,NDTOT+1
          WRITE(OP_STRING,'('' GKES(nhs,'',I1,'')='',4(1X,D12.5),'
     '      //'/:(13X,4(1X,D12.5)))') nkk,(GKES(nhs,nkk),
     '      nhs=-NDTOT-1,NHSTGK)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nkk
        !GQ matrix
C        ny1=NYNP(1,1,nh,np,1,2,nr) !eqtn number (row number) of GQ.
C        DO nn=1,NNT(nbqhp)
C          npnn=NPNE(nn,nbqhp)
C          DO nk1=1,NKT(nn,nbqhp)
C            MK=NKE(nk1,nn,nbqhp)
C            IF(MK.GT.0.AND.
C     '        nk1.LE.MAX(NKH(nh,np,2)-KTYP93(2,nr),1)) THEN
C              ny2=NYNP(MK,1,nh,npnn,2,2,nr) !local col number of GQ
C              CALL SPARSE(ny1,ny2,NYT(1,2,nx),nz,NZ_GQ_M,
C     '          NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,ERROR,*9999)
C              WRITE(OP_STRING,'('' GQ('',I5,'','',I5,'')='',D12.4)')
C     '          ny1,ny2,GQ(nz)
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C            ENDIF
C          ENDDO !nk1
C        ENDDO !nn
        DO nkk=1,NDTOT+1
          WRITE(OP_STRING,'('' GQES(nhs,'',I1,'')='',4(1X,D12.5),'
     '      //'/:(13X,4(1X,D12.5)))') nkk,(GQES(nhs,nkk),
     '      nhs=-NDTOT-1,NHSTGQ)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nkk
CC$OMP END CRITICAL(XEPGKQ_3DL_1)
      ENDIF

      CALL EXITS('XEPGKGQ_3DL')
      RETURN
 9999 CALL ERRORS('XEPGKGQ_3DL',ERROR)
      CALL EXITS('XEPGKGQ_3DL')
      RETURN 1
      END


