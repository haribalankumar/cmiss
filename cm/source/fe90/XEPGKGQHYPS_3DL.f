      SUBROUTINE XEPGKGQHYPS_3DL(nb1j,nb1jp,nbqh,nbqhp,nbuh,
     '  nbuhp,NDTOT,nh,NHSTGK,NHSTGQ,
     '  NKHE,NKH,NNMIN,np,NPNE,nr,CE,CURVCORRECT,DET,DRDNO,
     '  DSDX,GKES,GQES,PG_Q,PG_U,RAD,RG,SE,WG,XG1,XN,XNO,
     '  XN_GRAD,XPFP,XR,XR_GRAD,CALC_XR,ERROR,*)
C SMAR009 19/01/99 removed GK,GQ,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,nx,NYNP,

C#### Subroutine: XEPGKGQHYPS_3DL
C###  Description:
C###    XEPGKGQHYPS_3DL is an efficient version of XEPGKGQHYPS for 3D
C###    Laplace's equation in rectangular cartesian coordinates.
C###  See-Also: XEPGKGQHYPS

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER nb1j,nb1jp,nbqh,nbqhp,nbuh,nbuhp,NDTOT,
     '  nh,NHSTGK,NHSTGQ,NKHE(NKM,NNM,NHM),NKH(NHM,NPM,NCM),NNMIN,
     '  np,NPNE(NNM,NBFM),nr
C SMAR009 19/01/99 removed ISC_GK(NISC_GKM),ISC_GQ(NISC_GQM),
C  ISR_GK(NISR_GKM),ISR_GQ(NISR_GQM),nx,
C  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM),CURVCORRECT(2,2,NNM),
     '  DET(NBFM,0:NNM,NGM),DRDNO(NGM,NKM),DSDX(3,*),
     '  GKES(-NKM:NHM*NSM,NKM),
     '  GQES(-NKM:NHM*NSM,NKM),PG_Q(NSM,NUM,NGM),PG_U(NSM,NUM,NGM),
     '  RAD(NGM),RG(NGM),SE(NSM,NBFM),WG(NGM,NBM),XG1(NJM,NUM,NGM),
     '  XN(NJM,NGM),XNO(3,*),XN_GRAD(NJM,NGM),XPFP(*),XR(NJM,NGM),
     '  XR_GRAD(NJM,NGM)
C SMAR009 19/01/99 removed GK(NZ_GK_M),GQ(NZ_GQ_M),
      CHARACTER ERROR*(*)
      LOGICAL CALC_XR
!     Local Variables
      INTEGER MK,nd,ng,nhs,nk,nk1,nkk,nn,npnn,ns,nsnn
      !SMAR009 22/12/98 ,ny1,ny2,NYP2,nz
      REAL*8 DGREEN,HYPGREEN,FACTOR,FACTOR1,FACTOR2,NDOTNO,RDOTN,
     '  RDOTNO,SUM1,SUM2,SUM3(3),SUMQ,TERM(999),TERM1

      CALL ENTERS('XEPGKGQHYPS_3DL',*9999)

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
        FACTOR1=FACTOR
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
c cpb 10/7/95 unrolling loops
C          SUM=0.0d0
C          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
C            XR(nj,ng)=XG1(nj,1,ng)-XPFP(nj)
C            SUM=SUM+XR(nj,ng)*XR(nj,ng)
C          ENDDO !nj
C          RAD(ng)=DSQRT(SUM)
          XR(1,ng)=XG1(1,1,ng)-XPFP(1)
          XR(2,ng)=XG1(2,1,ng)-XPFP(2)
          XR(3,ng)=XG1(3,1,ng)-XPFP(3)
          RAD(ng)=DSQRT(XR(1,ng)*XR(1,ng)+XR(2,ng)*XR(2,ng)+
     '      XR(3,ng)*XR(3,ng))
        ENDDO !ng
        CALC_XR=.TRUE.
      ENDIF
      DO ng=1,NGT(nb1j)
c cpb 10/7/95 unrolling loops
C        DO nd=1,NDTOT
C          SUMR=0.0d0
C          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
C            SUMR=SUMR+XR(nj,ng)*XNO(nj,nd) !dR/dn0 (note minus sign)
C          ENDDO !nj
C          DRDNO(ng,nd)=-SUMR/RAD(ng)
C        ENDDO !nd
        DRDNO(ng,1)=-(XR(1,ng)*XNO(1,1)+XR(2,ng)*XNO(2,1)+
     '    XR(3,ng)*XNO(3,1))/RAD(ng)
        DRDNO(ng,2)=-(XR(1,ng)*XNO(1,2)+XR(2,ng)*XNO(2,2)+
     '    XR(3,ng)*XNO(3,2))/RAD(ng)
C*** Evaluate the dot products of XR and XN with the gradients of
C*** psi and the outward normal.
c cpb 10/7/95 unrolling loops
C        DO nj1=1,NJ_LOC(NJL_GEOM,0,nr)
C          XR_GRAD(nj1,ng)=0.0d0
C          XN_GRAD(nj1,ng)=0.0d0
C          DO nj2=1,NJ_LOC(NJL_GEOM,0,nr)
C            XR_GRAD(nj1,ng)=XR_GRAD(nj1,ng)+DSDX(nj1,nj2)*XR(nj2,ng)
C            XN_GRAD(nj1,ng)=XN_GRAD(nj1,ng)+DSDX(nj1,nj2)*XN(nj2,ng)
C          ENDDO !nj2
C        ENDDO !nj1
        XR_GRAD(1,ng)=DSDX(1,1)*XR(1,ng)+DSDX(1,2)*XR(2,ng)+
     '    DSDX(1,3)*XR(3,ng)
        XR_GRAD(2,ng)=DSDX(2,1)*XR(1,ng)+DSDX(2,2)*XR(2,ng)+
     '    DSDX(2,3)*XR(3,ng)
        XR_GRAD(3,ng)=DSDX(3,1)*XR(1,ng)+DSDX(3,2)*XR(2,ng)+
     '    DSDX(3,3)*XR(3,ng)
        XN_GRAD(1,ng)=DSDX(1,1)*XN(1,ng)+DSDX(1,2)*XN(2,ng)+
     '    DSDX(1,3)*XN(3,ng)
        XN_GRAD(2,ng)=DSDX(2,1)*XN(1,ng)+DSDX(2,2)*XN(2,ng)+
     '    DSDX(2,3)*XN(3,ng)
        XN_GRAD(3,ng)=DSDX(3,1)*XN(1,ng)+DSDX(3,2)*XN(2,ng)+
     '    DSDX(3,3)*XN(3,ng)
      ENDDO !ng

C cpb 2/1/98 Adding BEM elemental stiffness matrices

C*** Loop over derivative directions
      DO nd=1,2

C*** Set derivative number nk for NYNP array i.e. the 1st, 2nd, 3rd
C*** derivative equations are for the 2nd, 3rd ,4th nk values.
        nk=nd+1
        IF((KTYP92.EQ.3).AND.(nd.EQ.2)) THEN
C*** A derivative equation replaces the conventional BIE equation.
C*** Assemble this in the nk=1 place.
          nk=1
        ENDIF
        IF(NDTOT.GT.NKT(0,nbuhp)) nk=1

C*** Start of first integral calculation
        SUM1=0.0d0
        DO ng=1,NGT(nb1j)
c cpb 10/7/95 Inlining Green's functions
C          SUM1=SUM1+RG(ng)*WG(ng,nb1j)*
C     '      HYPGREEN(IGREN(nr),1,CE,RAD(ng),DET(nb1jp,NNMIN,ng),
C     '      XN(1,ng),XNO(1,nd),XR(1,ng))*XR_GRAD(NJ_LOC(NJL_GEOM,0,nr),ng)
          NDOTNO=XN(1,ng)*XNO(1,nd)+XN(2,ng)*XNO(2,nd)+XN(3,ng)*
     '      XNO(3,nd)
          RDOTN=XR(1,ng)*XN(1,ng)+XR(2,ng)*XN(2,ng)+XR(3,ng)*XN(3,ng)
          RDOTNO=XR(1,ng)*XNO(1,nd)+XR(2,ng)*XNO(2,nd)+XR(3,ng)*
     '      XNO(3,nd)
          HYPGREEN=(NDOTNO-3.0d0*RDOTN*RDOTNO/(RAD(ng)*RAD(ng)))*
     '      DET(nb1jp,NNMIN,ng)/(RAD(ng)*RAD(ng)*RAD(ng))
          SUM1=SUM1+RG(ng)*WG(ng,nb1j)*HYPGREEN*XR_GRAD(3,ng)
        ENDDO !ng
C        ny1=NYNP(nk,1,nh,np,1,2,nr) !eqtn number (row number) of GQ.
C        NYP2=NYNP(1,1,nh,np,2,2,nr) !local column number of GQ.
C        CALL SPARSE(ny1,NYP2,NYT(1,2,nx),nz,NZ_GQ_M,NZT(2,nx),
C     '    ISC_GQ,ISR_GQ,KTYP24,ERROR,*9999)
C        GQ(nz)=GQ(nz)+SUM1/(FACTOR1*FACTOR2)
        GQES(-1,nk)=GQES(-1,nk)+SUM1/(FACTOR1*FACTOR2)
C*** GK matrix components
        ns=0
        nhs=0
C        ny1=NYNP(nk,1,nh,np,1,1,nr) !eqtn number (row number) of GK.
c CPB 10/7/95 Precalculate HYPGREEN*RG*WG
        DO ng=1,NGT(nbuh)
          NDOTNO=XN(1,ng)*XNO(1,nd)+XN(2,ng)*XNO(2,nd)+XN(3,ng)*
     '      XNO(3,nd)
          RDOTN=XR(1,ng)*XN(1,ng)+XR(2,ng)*XN(2,ng)+XR(3,ng)*XN(3,ng)
          RDOTNO=XR(1,ng)*XNO(1,nd)+XR(2,ng)*XNO(2,nd)+XR(3,ng)*
     '      XNO(3,nd)
          HYPGREEN=(NDOTNO-3.0d0*RDOTN*RDOTNO/(RAD(ng)*RAD(ng)))*
     '      DET(nb1jp,NNMIN,ng)/(RAD(ng)*RAD(ng)*RAD(ng))
          TERM(ng)=HYPGREEN*RG(ng)*WG(ng,nbuh)
        ENDDO !ng
        DO nn=1,NNT(nbuhp) !Dependent variable loop
          npnn=NPNE(nn,nbuhp)
          DO nk1=1,NKT(nn,nbuhp)
            MK=NKHE(nk1,nn,nh)
            ns=ns+1
            nhs=nhs+1
C cpb 3/1/98 I think this should test on NKH(nh,npnn,1)
C            IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,np,1)-KTYP93(1,nr),1)) THEN
            IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,1)-KTYP93(1,nr),1))
     '        THEN
C              ny2=NYNP(MK,1,nh,npnn,2,1,nr) !local column number of GK
              SUM2=0.0d0
              IF(nn.EQ.NNMIN) THEN
                IF(nk1.EQ.1) THEN
                  DO ng=1,NGT(nbuh)
c cpb 10/7/95 Inlining Green's functions
C                    SUM2=SUM2+RG(ng)*
C     '                 HYPGREEN(IGREN(nr),0,CE,RAD(ng),
C     '                 DET(nbuhp,NNMIN,ng),
C     '                 XN(1,ng),XNO(1,nd),XR(1,ng))*
C     '                 (PG(ns,1,ng,nbuh)-1.0d0)*WG(ng,nbuh)
                    SUM2=SUM2+TERM(ng)*(PG_U(ns,1,ng)-1.0d0)
                  ENDDO !End of ng loop
                ELSE IF(nk1.LT.4) THEN
                  DO ng=1,NGT(nbuh)
c cpb 10/7/95 Inlining Green's functions
C                    SUM2=SUM2+RG(ng)*
C     '                 HYPGREEN(IGREN(nr),0,CE,RAD(ng),
C     '                 DET(nbuhp,NNMIN,ng),XN(1,ng),
C     '                 XNO(1,nd),XR(1,ng))*WG(ng,nbuh)*
C     '                 (SE(ns,nbuhp)*PG(ns,1,ng,nbuh)-
C     '                  XR_GRAD(nk1-1,ng))
                    SUM2=SUM2+TERM(ng)*(SE(ns,nbuhp)*PG_U(ns,1,ng)-
     '                XR_GRAD(nk1-1,ng))
                  ENDDO !End of ng loop
                ELSE
                  DO ng=1,NGT(nbuh)
c cpb 10/7/95 Inlining Green's functions
C                    SUM2=SUM2+RG(ng)*
C     '               HYPGREEN(IGREN(nr),0,CE,RAD(ng),
C     '               DET(nbuhp,NNMIN,ng),
C     '               XN(1,ng),XNO(1,nd),XR(1,ng))*WG(ng,nbuh)*
C     '               SE(ns,nbuhp)*PG(ns,1,ng,nbuh)
                    SUM2=SUM2+TERM(ng)*PG_U(ns,1,ng)
                  ENDDO !End of ng loop
                  SUM2=SUM2*SE(ns,nbuhp)
                ENDIF
              ELSE
                DO ng=1,NGT(nbuh)
c cpb 10/7/95 Inlining Green's functions
C                  SUM2=SUM2+RG(ng)*HYPGREEN(IGREN(nr),0,
C     '               CE,RAD(ng),DET(nbuhp,NNMIN,ng),
C     '               XN(1,ng),XNO(1,nd),XR(1,ng))*PG(ns,1,ng,nbuh)*
C     '               WG(ng,nbuh)
                  SUM2=SUM2+TERM(ng)*PG_U(ns,1,ng)
                ENDDO !End of ng loop
                SUM2=SUM2*SE(ns,nbuhp)
              ENDIF
C              CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
C     '          NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C              GK(nz)=GK(nz)+SUM2/FACTOR
              GKES(nhs,nk)=GKES(nhs,nk)+SUM2/FACTOR
            ENDIF
          ENDDO !End of nk1 loop
        ENDDO !End of nn loop
C*** End of first integral calculation
C*** Start of second integral calculation
c cpb 10/7/95 Unrolling loops
C        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)-1
C          SUM3(nj)=0.0d0
C        ENDDO
        SUM3(1)=0.0d0
        SUM3(2)=0.0d0
        DO ng=1,NGT(nb1j)
c cpb 10/7/95 Inlining Green's functions
C          TERM=DGREEN(IGREN(nr),0,nh,nh,CE,DRDNO(ng,nd),RAD(ng),
C     '         DET(nb1jp,NNMIN,ng),XN(1,ng),XR(1,ng))
C     '         *RG(ng)*WG(ng,nb1j)
          DGREEN=-DET(nb1jp,NNMIN,ng)*DRDNO(ng,nd)/(RAD(ng)*RAD(ng))
          TERM1=DGREEN*RG(ng)*WG(ng,nb1j)
c cpb 10/7/95 unrolling loops
C          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)-1
C            SUM3(nj)=SUM3(nj)+TERM*XN_GRAD(nj,ng)
C          ENDDO !nj
          SUM3(1)=SUM3(1)+TERM1*XN_GRAD(1,ng)
          SUM3(2)=SUM3(2)+TERM1*XN_GRAD(2,ng)
        ENDDO !End of ng loop
C cpb 10/7/95 unrolling loops
C        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)-1
C          nk1=nj+1
C          NYP2=NYNP(nk1,1,nh,np,2,1,nr) !local column number of GK
C          CALL SPARSE(ny1,NYP2,NYT(1,1,nx),nz,NZ_GK_M,NZT(1,nx),
C     '      ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C          GK(nz)=GK(nz)+SUM3(nj)
C        ENDDO
C        NYP2=NYNP(2,1,nh,np,2,1,nr) !local column number of GK
C        CALL SPARSE(ny1,NYP2,NYT(1,1,nx),nz,NZ_GK_M,NZT(1,nx),
C     '    ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C        GK(nz)=GK(nz)+SUM3(1)/FACTOR
        GKES(-2,nk)=GKES(-2,nk)+SUM3(1)/FACTOR
C        NYP2=NYNP(3,1,nh,np,2,1,nr) !local column number of GK
C        CALL SPARSE(ny1,NYP2,NYT(1,1,nx),nz,NZ_GK_M,NZT(1,nx),
C     '    ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C        GK(nz)=GK(nz)+SUM3(2)/FACTOR
        GKES(-3,nk)=GKES(-3,nk)+SUM3(2)/FACTOR
        ns=0
        nhs=0
C*** GQ matrix
C        ny1=NYNP(nk,1,nh,np,1,2,nr) !eqtn number (row number) of GQ.
c CPB 10/7/95 Precalculate DGREEN*RG*WG
        DO ng=1,NGT(nbqh)
          DGREEN=-DET(nbqhp,NNMIN,ng)*DRDNO(ng,nd)/(RAD(ng)*RAD(ng))
          TERM(ng)=DGREEN*RG(ng)*WG(ng,nbqh)
        ENDDO !ng
        DO nn=1,NNT(nbqhp)
          npnn=NPNE(nn,nbqhp)
          nsnn=ns
          DO nk1=1,NKT(nn,nbqhp)
            MK=NKHE(nk1,nn,nh)
            ns=ns+1
            nhs=nhs+1
C cpb 3/1/98 I think this should test on NKH(nh,npnn,2)
C            IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,np,2)-KTYP93(2,nr),1))
C     '        THEN
            IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,2)-KTYP93(2,nr),1))
     '        THEN
C              ny2=NYNP(MK,1,nh,npnn,2,2,nr) !local column number of GQ
              SUMQ=0.0d0
              IF((nn.EQ.NNMIN).AND.(nk1.EQ.1)) THEN
                DO ng=1,NGT(nbqh)
c cpb 10/7/95 Inlining Greens's functions
C                  SUMQ=SUMQ+(PG(ns,1,ng,nbqh)*SE(ns,nbqhp)-
C     '               XN_GRAD(NJ_LOC(NJL_GEOM,0,nr),ng))*
C     '               DGREEN(IGREN(nr),1,nh,nh,CE,DRDNO(ng,nd),RAD(ng),
C     '               DET(nbqhp,NNMIN,ng),XN(1,ng),XR(1,ng))*
C     '               RG(ng)*WG(ng,nbqh)
                  SUMQ=SUMQ+TERM(ng)*(SE(ns,nbqhp)*PG_Q(ns,1,ng)-
     '              XN_GRAD(3,ng))
                ENDDO !End ng loop
              ELSE
                DO ng=1,NGT(nbqh)
c cpb 10/7/95 Inlining Green's functions
C                  SUMQ=SUMQ+PG(ns,1,ng,nbqh)*SE(ns,nbqhp)*
C     '               DGREEN(IGREN(nr),1,nh,nh,CE,DRDNO(ng,nd),RAD(ng),
C     '               DET(nbqhp,NNMIN,ng),XN(1,ng),XR(1,ng))*
C     '               RG(ng)*WG(ng,nbqh)
                  SUMQ=SUMQ+TERM(ng)*PG_Q(ns,1,ng)
                ENDDO !End ng loop
                SUMQ=SUMQ*SE(ns,nbqhp)
              ENDIF
C              CALL SPARSE(ny1,ny2,NYT(1,2,nx),nz,NZ_GQ_M,
C     '          NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,ERROR,*9999)
C              GQ(nz)=GQ(nz)+SUMQ/(FACTOR1*FACTOR2)
              GQES(nhs,nk)=GQES(nhs,nk)+SUMQ/(FACTOR1*FACTOR2)
C cpb 20/11/96 Adding bem curvature correction
              IF(BEMCURVATURECORRECTION.AND.(nk1.EQ.2.OR.nk1.EQ.3)) THEN
C                ny2=NYNP(2,1,1,npnn,2,1,nr) !loc col # of GK
C                CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
C     '            NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C                GK(nz)=GK(nz)-SUMQ/FACTOR1*CURVCORRECT(1,nk1-1,nn)
C                ny2=NYNP(3,1,1,npnn,2,1,nr) !loc col # of GK
C                CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
C     '            NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C                GK(nz)=GK(nz)-SUMQ/FACTOR1*CURVCORRECT(2,nk1-1,nn)
                GKES(nsnn+2,nk)=GKES(nsnn+2,nk)-SUMQ/FACTOR1*
     '            CURVCORRECT(1,nk1-1,nn)
                GKES(nsnn+3,nk)=GKES(nsnn+3,nk)-SUMQ/FACTOR1*
     '            CURVCORRECT(2,nk1-1,nn)
              ENDIF
            ENDIF
          ENDDO !End of nk1 loop
        ENDDO !End of nn loop
C*** End of second integral calculation

        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(XEPGKGQHYPS_3DL_1)
          WRITE(OP_STRING,'('' Element contribution to'//
     '     ' Global matrices'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Node '',I5)') np
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          !GK matrix
C          ny1=NYNP(nk,1,nh,np,1,1,nr) !eqtn number (row number) of GK.
C          DO nn=1,NNT(nbuhp)
C            npnn=NPNE(nn,nbuhp)
C            DO nk1=1,NKT(nn,nbuhp)
C              MK=NKE(nk1,nn,nbuhp)
C              IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,np,1)-KTYP93(1,nr),1))
C     '          THEN
C                ny2=NYNP(MK,1,nh,npnn,2,1,nr)!local column number of GK
C                CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
C     '            NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C                WRITE(OP_STRING,'('' GK('',I5,'','',I5,'')='',D12.4)')
C     '            ny1,ny2,GK(nz)
C                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C              ENDIF
C            ENDDO
C          ENDDO
          DO nkk=1,NDTOT+1
            WRITE(OP_STRING,'('' GKES(nhs,'',I1,'')='',4(1X,D12.5),'
     '        //'/:(13X,4(1X,D12.5)))') nkk,(GKES(nhs,nkk),
     '        nhs=-NDTOT-1,NHSTGK)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO !nkk
          !GQ matrix
C          ny1=NYNP(nk,1,nh,np,1,2,nr) !eqtn number (row number) of GQ.
C          DO nn=1,NNT(nbqhp)
C            npnn=NPNE(nn,nbqhp)
C            DO nk1=1,NKT(nn,nbqhp)
C              MK=NKE(nk1,nn,nbqhp)
C              IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,np,2)-KTYP93(2,nr),1))
C     '          THEN
C                ny2=NYNP(MK,1,nh,npnn,2,2,nr) !local col number of GQ
C                CALL SPARSE(ny1,ny2,NYT(1,2,nx),nz,NZ_GQ_M,
C     '            NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,ERROR,*9999)
C                WRITE(OP_STRING,'('' GQ('',I5,'','',I5,'')='',D12.4)')
C     '            ny1,ny2,GQ(nz)
C                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C              ENDIF
C            ENDDO !nk1
C          ENDDO !nn
          DO nkk=1,NDTOT+1
            WRITE(OP_STRING,'('' GQES(nhs,'',I1,'')='',4(1X,D12.5),'
     '        //'/:(13X,4(1X,D12.5)))') nkk,(GQES(nhs,nkk),
     '        nhs=-NDTOT-1,NHSTGQ)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO !nkk
CC$OMP END CRITICAL(XEPGKGQHYPS_3DL_1)
        ENDIF !dop
      ENDDO !End of nk loop

      CALL EXITS('XEPGKGQHYPS_3DL')
      RETURN
 9999 CALL ERRORS('XEPGKGQHYPS_3DL',ERROR)
      CALL EXITS('XEPGKGQHYPS_3DL')
      RETURN 1
      END


