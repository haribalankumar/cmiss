      SUBROUTINE XEPGKGQHYP(nb1j,nb1jp,
     '  nbqh,nbqhp,nbuh,nbuhp,NDTOT,ne,nh,NHSTGK,NHSTGQ,NKHE,NKH,NNMIN,
     '  np,NPNE,nr,CE,CURVCORRECT,DET,DRDNO,DSDX,GKES,
     '  GQES,PG_Q,PG_U,RAD,RG,SE,WG,XG1,XN,XNO,XN_GRAD,XPFP,XR,XR_GRAD,
     '  CALC_XR,ERROR,*)
C SMAR009 19/01/99 removed GK,GQ,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,nx,NYNP,
C#### Subroutine: XEPGKGQHYP
C###  Description:
C###    XEPGKGQHYP calculates the derivative BEM integrals for an
C###    element ne and a singularity at node np for the case
C###    when np is not contained in ne.

C**** This routine is required when hermite
C**** interpolation is used for the dependent variable.
C**** DSDX(i,j) contains dsi/dXj and dn/dXj evaluated at the point np.
C**** XNO contains the direction(s) in which the BIE was
C**** differentiated.  If KTYP92=3 then the conventional BIE has
C**** been replaced by a derivative BIE.
C**** NOTE no need to calculate RG etc in XEGKGQHYP if the same basis
C**** is used in XEGKGQ and XEGKGQHYP and XEGKGQ has been called
C**** (i.e.ktyp92 < 3).

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER  nb1j,nb1jp,nbqh,nbqhp,nbuh,nbuhp,
     '  NDTOT,ne,nh,NHSTGK,NHSTGQ,NKHE(NKM,NNM,NHM,NEM),
     '  NKH(NHM,NPM,NCM),NNMIN,np,NPNE(NNM,NBFM),nr
C SMAR009 19/01/99 removed ISR_GK(NISR_GKM),ISC_GK(NISC_GKM),
C  ISC_GQ(NISC_GQM),ISR_GQ(NISR_GQM),nx,
C  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM),CURVCORRECT(2,2,NNM),DET(NBFM,0:NNM,NGM),
     '  DRDNO(NGM,NKM),DSDX(3,3),GKES(-NKM:NHM*NSM,NKM),
     '  GQES(-NKM:NHM*NSM,NKM),PG_Q(NSM,NUM,NGM),
     '  PG_U(NSM,NUM,NGM),RAD(NGM),RG(NGM),SE(NSM,NBFM,NEM),
     '  XG1(NJM,NUM,NGM),XN(NJM,NGM),XNO(3,*),XN_GRAD(NJM,NGM),
     '  XPFP(*),XR(NJM,NGM),XR_GRAD(NJM,NGM),WG(NGM,NBM)
C SMAR009 19/01/99 removed GK(NZ_GK_M),GQ(NZ_GQ_M),
      CHARACTER ERROR*(*)
      LOGICAL CALC_XR
!     Local Variables
      INTEGER MK,nd,ng,nhs,nj,nj1,nj2,nk,nkk,nk1,nn,npnn,ns,nsnn
      !SMAR009 22/12/98  ny1,ny2,NYP1,NYP2,nz
      REAL*8 DGREEN,FACTOR2,HYPGREEN,SUM,SUM1,SUM2,SUM3(3),SUMQ,SUMR,
     '  TERM1,TERM2,XGC(3),XPC(3)
C      CHARACTER FORMAT*500

      CALL ENTERS('XEPGKGQHYP',*9999)

      IF(IGREN(nr).EQ.9.OR.IGREN(nr).EQ.10) THEN
        !Poisson equation with a special source term - for use in
        !constructing uniform-dipole layer transfer matrices
        IF(DABS(CE(1)).GT.ZERO_TOL) THEN
          FACTOR2=CE(1)
        ELSE
          ERROR='>>Error in material property of Poisson equation'
          GOTO 9999
        ENDIF
      ELSE
        FACTOR2=1.0d0
      ENDIF

C*** Calculate gauss point dependent arrays
      IF(CALC_XR) THEN
        IF(ITYP10(nr).EQ.1) THEN
          DO ng=1,NGT(nb1j)
            SUM=0.0d0
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              XR(nj,ng)=XG1(nj,1,ng)-XPFP(nj)
              SUM=SUM+XR(nj,ng)*XR(nj,ng)
            ENDDO !nj
            RAD(ng)=DSQRT(SUM)
          ENDDO !ng
        ELSE
          DO NG=1,NGT(nb1j)
C*** Transform into cartesian to find distance (RAD).
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              XR(nj,ng)=XPFP(nj)
            ENDDO !nj
            CALL COORD(ITYP10(nr),1,XR(1,ng),XPC,ERROR,*9999)
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              XR(nj,ng)=XG1(nj,1,ng)
            ENDDO !nj
            CALL COORD(ITYP10(nr),1,XR(1,ng),XGC,ERROR,*9999)
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              XR(nj,ng)=XGC(nj)-XPC(nj)
              SUM=SUM+XR(nj,ng)*XR(nj,ng)
            ENDDO !nj
            RAD(ng)=DSQRT(SUM)
          ENDDO !ng
        ENDIF
        CALC_XR=.TRUE.
      ENDIF
      DO ng=1,NGT(nb1j)
        DO nd=1,NDTOT
          SUMR=0.0d0
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            SUMR=SUMR+XR(nj,ng)*XNO(nj,nd) !dR/dn0 (note minus sign)
          ENDDO !nj
          DRDNO(ng,nd)=-SUMR/RAD(ng)
        ENDDO !nd
C*** Evaluate the dot products of XR and XN with the gradients of
C*** psi and the outward normal.
        DO nj1=1,NJ_LOC(NJL_GEOM,0,nr)
          XR_GRAD(nj1,ng)=0.0d0
          XN_GRAD(nj1,ng)=0.0d0
          DO nj2=1,NJ_LOC(NJL_GEOM,0,nr)
            XR_GRAD(nj1,ng)=XR_GRAD(nj1,ng)+DSDX(nj1,nj2)*XR(nj2,ng)
            XN_GRAD(nj1,ng)=XN_GRAD(nj1,ng)+DSDX(nj1,nj2)*XN(nj2,ng)
          ENDDO !nj2
        ENDDO !nj1
      ENDDO !ng

C*** Loop over derivative directions
      DO nd=1,NDTOT

C*** Set derivative number nk for NYNP array i.e. the 1st, 2nd, 3rd
C*** derivative equations are for the 2nd, 3rd ,4th nk values.
        nk=nd+1
        IF((KTYP92.EQ.3).AND.(nd.EQ.NDTOT)) THEN
C*** A derivative equation replaces the conventional BIE equation.
C*** Assemble this in the nk=1 place.
          nk=1
        ENDIF
        IF(NDTOT.GT.NKT(0,nbuhp)) nk=1  !???? Compare to nkh(nh,np,1)

C*** Start of first integral calculations
        SUM1=0.0d0
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          SUM3(nj)=0.0d0
        ENDDO
        DO ng=1,NGT(nb1j)
C*** Term1 has no conductivity term in it, Term2 has.
          TERM1=RG(ng)*WG(ng,nb1j)*HYPGREEN(IGREN(nr),0,CE,RAD(ng),
     '         DET(nb1jp,NNMIN,ng),XN(1,ng),XNO(1,nd),XR(1,ng))
          TERM2=RG(ng)*WG(ng,nb1j)*HYPGREEN(IGREN(nr),1,CE,RAD(ng),
     '         DET(nb1jp,NNMIN,ng),XN(1,ng),XNO(1,nd),XR(1,ng))
          SUM1=SUM1+TERM1
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)-1
            SUM3(nj)=SUM3(nj)+TERM1*XR_GRAD(nj,ng)
          ENDDO
          SUM3(NJ_LOC(NJL_GEOM,0,nr))
     '      =SUM3(NJ_LOC(NJL_GEOM,0,nr))+
     '      TERM2*XR_GRAD(NJ_LOC(NJL_GEOM,0,nr),ng)
        ENDDO !End of ng loop
C cpb 2/1/98 adding bem element stiffness matrices
C        ny1=NYNP(nk,1,nh,np,1,1,nr) !eqtn number (row number) of GK.
C        NYP1=NYNP(1,1,nh,np,2,1,nr) !local column number of GK
C        CALL SPARSE(ny1,NYP1,NYT(1,1,nx),nz,NZ_GK_M,NZT(1,nx),
C     '    ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C        GK(nz)=GK(nz)-SUM1
        GKES(-1,nk)=GKES(-1,nk)-SUM1
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)-1
          nk1=nj+1
C          NYP2=NYNP(nk1,1,nh,np,2,1,nr) !local column number of GK
C          CALL SPARSE(ny1,NYP2,NYT(1,1,nx),nz,NZ_GK_M,NZT(1,nx),
C     '      ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C          GK(nz)=GK(nz)-SUM3(nj)
          GKES(-nk1,nk)=GKES(-nk1,nk)-SUM3(nj)
        ENDDO
C        ny1=NYNP(nk,1,nh,np,1,2,nr) !eqtn number (row number) of GQ.
C        NYP2=NYNP(1,1,nh,np,2,2,nr) !local column number of GQ
C        CALL SPARSE(ny1,NYP2,NYT(1,2,nx),nz,NZ_GQ_M,NZT(2,nx),
C     '    ISC_GQ,ISR_GQ,KTYP24,ERROR,*9999)
C        GQ(nz)=GQ(nz)+SUM3(NJ_LOC(NJL_GEOM,0,nr))/FACTOR2
        GQES(-1,nk)=GQES(-1,nk)+SUM3(NJ_LOC(NJL_GEOM,0,nr))/FACTOR2
        !GK matrix components
        ns=0
        nhs=0
C        nv=1
C        ny1=NYNP(nk,1,nh,np,1,1,nr) !eqtn number (row number) of GK.
        DO nn=1,NNT(nbuhp) !Dependent variable loop
          npnn=NPNE(nn,nbuhp)
          DO nk1=1,NKT(nn,nbuhp)
            MK=NKHE(nk1,nn,nh,ne)
            ns=ns+1
            nhs=nhs+1
C cpb 3/1/98 I think this should test on NKH(nh,npnn,1)
C            IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,np,1)-KTYP93(1,nr),1)) THEN
            IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,1)-KTYP93(1,nr),1))
     '        THEN
C              ny2=NYNP(MK,nv,nh,npnn,2,1,nr) !local column number of GK
              SUM2=0.0d0
              DO ng=1,NGT(nbuh)
                SUM2=SUM2+RG(ng)*HYPGREEN(IGREN(nr),0,CE,RAD(ng),
     '               DET(nbuhp,NNMIN,ng),XN(1,ng),XNO(1,nd),XR(1,ng))*
     '               PG_U(ns,1,ng)*WG(ng,nbuh)
              ENDDO !End of ng loop
C              CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
C     '          NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C              GK(nz)=GK(nz)+SUM2*SE(ns,nbuhp,ne)
              GKES(nhs,nk)=GKES(nhs,nk)+SUM2*SE(ns,nbuhp,ne)
            ENDIF
          ENDDO !End of nk1 loop
        ENDDO !End of nn loop
        !End of first integral calculations
        !Start of second integral calculations
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          SUM3(nj)=0.0d0
        ENDDO
        DO ng=1,NGT(nb1j)
          TERM1=DGREEN(IGREN(nr),0,nh,nh,CE,DRDNO(ng,nd),RAD(ng),
     '         DET(nb1jp,NNMIN,ng),XN(1,ng),XR(1,ng))*
     '         RG(ng)*WG(ng,nb1j)
          TERM2=DGREEN(IGREN(nr),1,nh,nh,CE,DRDNO(ng,nd),RAD(ng),
     '         DET(nb1jp,NNMIN,ng),XN(1,ng),XR(1,ng))*
     '         RG(ng)*WG(ng,nb1j)
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)-1
            SUM3(nj)=SUM3(nj)+TERM1*XN_GRAD(nj,ng)
          ENDDO
          SUM3(NJ_LOC(NJL_GEOM,0,nr))=SUM3(NJ_LOC(NJL_GEOM,0,nr))+
     '      TERM2*XN_GRAD(NJ_LOC(NJL_GEOM,0,nr),ng)
        ENDDO !End of ng loop
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)-1
          nk1=nj+1
C          NYP2=NYNP(nk1,nv,nh,np,2,1,nr) !local column number of GK
C          CALL SPARSE(ny1,NYP2,NYT(1,1,nx),nz,NZ_GK_M,NZT(1,nx),
C     '      ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C          GK(nz)=GK(nz)+SUM3(nj)
          GKES(-nk1,nk)=GKES(-nk1,nk)+SUM3(nj)
        ENDDO
        !GQ matrix
C        ny1=NYNP(nk,1,nh,np,1,2,nr) !eqtn number (row number) of GQ.
C        NYP2=NYNP(1,nv,nh,np,2,2,nr) !local column number of GQ
C        CALL SPARSE(ny1,NYP2,NYT(1,2,nx),nz,NZ_GQ_M,NZT(2,nx),
C     '    ISC_GQ,ISR_GQ,KTYP24,ERROR,*9999)
C        GQ(nz)=GQ(nz)-SUM3(NJ_LOC(NJL_GEOM,0,nr))/FACTOR2
        GQES(-1,nk)=GQES(-1,nk)-SUM3(NJ_LOC(NJL_GEOM,0,nr))/FACTOR2
        ns=0
        nhs=0
C        nv=1 !temporary
        DO nn=1,NNT(nbqhp)
          npnn=NPNE(nn,nbqhp)
          nsnn=ns
          DO nk1=1,NKT(nn,nbqhp)
            MK=NKHE(nk1,nn,nh,ne)
            ns=ns+1
            nhs=nhs+1
C cpb 3/1/98 I think this should test on NKH(nh,npnn,2)
C            IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,np,2)-KTYP93(2,nr),1)) THEN
            IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,2)-KTYP93(2,nr),1))
     '        THEN
C              ny2=NYNP(MK,nv,nh,npnn,2,2,nr) !local column number of GQ
              SUMQ=0.0d0
              DO ng=1,NGT(nbqh)
                SUMQ=SUMQ+PG_Q(ns,1,ng)*
     '              DGREEN(IGREN(nr),1,nh,nh,CE,DRDNO(ng,nd),RAD(ng),
     '              DET(nbqhp,NNMIN,ng),XN(1,ng),XR(1,ng))*
     '              RG(ng)*WG(ng,nbqh)
              ENDDO !End ng loop
C              CALL SPARSE(ny1,ny2,NYT(1,2,nx),nz,NZ_GQ_M,
C     '          NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,ERROR,*9999)
C              GQ(nz)=GQ(nz)+SUMQ*SE(ns,nbqhp,ne)/FACTOR2
              GQES(nhs,nk)=GQES(nhs,nk)+SUMQ*SE(ns,nbqhp,ne)/FACTOR2
C cpb 20/11/96 Adding bem curvature correction
              IF(BEMCURVATURECORRECTION.AND.(nk1.EQ.2.OR.nk1.EQ.3)) THEN
                IF(NIT(nb1jp).EQ.1) THEN
C                  ny2=NYNP(MK,nv,1,npnn,2,1,nr) !local column # of GK
C                  CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
C     '              NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C                  GK(nz)=GK(nz)-SUMQ*SE(ns,nbqhp,ne)*CURVCORRECT(1,1,nn)
                  GKES(nhs,nk)=GKES(nhs,nk)-SUMQ*SE(ns,nbqhp,ne)*
     '              CURVCORRECT(1,1,nn)
                ELSE
                  DO nkk=1,2
C                    ny2=NYNP(nkk+1,nv,1,npnn,2,1,nr) !loc col # of GK
C                    CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
C     '                NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C                    GK(nz)=GK(nz)-SUMQ*SE(ns,nbqhp,ne)*
C     '                CURVCORRECT(nkk,nk1-1,nn)
                    GKES(nsnn+nkk+1,nk)=GKES(nsnn+nkk+1,nk)-SUMQ*
     '                SE(ns,nbqhp,ne)*CURVCORRECT(nkk,nk1-1,nn)
                  ENDDO !nkk
                ENDIF
              ENDIF
            ENDIF
          ENDDO !End of nk1 loop
        ENDDO !End of nn loop
        !End of second integral calculations.

        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(XEPGKGQHYP_1)
          WRITE(OP_STRING,'('' Element contribution to'//
     '     ' Global matrices'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Node '',I5)')np
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          !GK matrix
C          ny1=NYNP(nk,1,nh,np,1,1,nr) !eqtn number (row number) of GK.
C          DO nn=1,NNT(nbuhp)
C            npnn=NPNE(nn,nbuhp)
C            FORMAT='('' GK('',I3,'','',I3,'')='',D12.4)'
C            DO nk1=1,NKT(nn,nbuhp)
C              MK=NKE(nk1,nn,nbuhp,ne)
C              IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,np,1)-KTYP93(1,nr),1))
C     '          THEN
C                ny2=NYNP(MK,1,nh,npnn,2,1,nr)!local column number of GK
C                CALL SPARSE(ny1,NYP2,NYT(1,1,nx),nz,NZ_GK_M,
C     '            NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C                WRITE(OP_STRING,FORMAT)ny1,ny2,GK(nz)
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
C          nv=1 !temporary
C          DO nn=1,NNT(nbqhp)
C            npnn=NPNE(nn,nbqhp)
C            FORMAT='('' GQ('',I3,'','',I3,'')='',D12.4)'
C            DO nk1=1,NKT(nn,nbqhp)
C              MK=NKE(nk1,nn,nbqhp,ne)
C              IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,np,2)-KTYP93(2,nr),1))
C     '          THEN
C                ny2=NYNP(MK,nv,nh,npnn,2,2,nr) !local col number of GQ
C                CALL SPARSE(ny1,ny2,NYT(1,2,nx),nz,NZ_GQ_M,
C     '            NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,ERROR,*9999)
C                WRITE(OP_STRING,FORMAT)ny1,ny2,GQ(nz)
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
CC$OMP END CRITICAL(XEPGKGQHYP_1)
        ENDIF !dop
      ENDDO !End of nk loop

      CALL EXITS('XEPGKGQHYP')
      RETURN
 9999 CALL ERRORS('XEPGKGQHYP',ERROR)
      CALL EXITS('XEPGKGQHYP')
      RETURN 1
      END


