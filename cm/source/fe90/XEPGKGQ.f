      SUBROUTINE XEPGKGQ(GREENFLAG,
     '  nb1j,nb1jp,nbqh,nbqhp,nbuh,nbuhp,NDTOT,nh,NHE,NHSTGK,NHSTGQ,
     '  NKHE,NKH,NNMIN,np,NPNE,nr,nx,CE,CURVCORRECT,DET,
     '  DRDN,GKES,GQES,PG_Q,PG_U,RAD,RD,RG,SE,WG,XIG,XG1,
     '  XN,XPFP,XR,CALC_XR,ERROR,*)
C  SMAR009 18/01/99 removed GK,GQ,ISC_GK,ISC_GQ,ISR_GK,ISR_GQ,NYNP,
C#### Subroutine: XEPGKGQ
C###  Description:
C###    XEPGKGQ contains the main Boundary Element loop. Calculates the
C###    element integrals for element ne for a given node np and given
C###    dependent variable number nh.

C**** For hypersingular BIE the routine XEGKGQHYP is used.
C**** nbuh is the basis function number of the dependent variable nh
C**** nbqh is the basis function number of its normal derivative
C**** nb1j is the basis function number of the geometric variable
C**** nbuhp is the parent basis function number of the dependent var nh
C**** nbqhp is the parent basis function number of its normal derivative
C**** nb1jp is the parent basis function number of the geometric var nh
C**** DET(nbpf,0:nn,ng) is the determinant of the mapping from
C**** the ith part of the split element in (psi1,psi2) space to a
C**** unit element in (psi1',psi2') space.  The splitting is detemined
C**** by the location nn of the singularity np.  If element is not
C**** split, nn=0 and the identity map is used.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER GREENFLAG,
     '  nb1j,nb1jp,nbqh,nbqhp,nbuh,nbuhp,NDTOT,nh,NHE,
     '  NHSTGK,NHSTGQ,NKHE(NKM,NNM,NHM),NKH(NHM,NPM,NCM),NNMIN,np,
     '  NPNE(NNM,NBFM),nr,nx
C  SMAR009 18/01/99 removed ISC_GK(NISC_GKM),ISC_GQ(NISC_GQM),
C   ISR_GK(NISR_GKM),ISR_GQ(NISR_GQM),
C   NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM),CURVCORRECT(2,2,NNM),DET(NBFM,0:NNM,NGM),DRDN(NGM),
     '  GKES(-NKM:NHM*NSM,NKM),
     '  GQES(-NKM:NHM*NSM,NKM),PG_Q(NSM,NUM,NGM),PG_U(NSM,NUM,NGM),
     '  RAD(NGM),RD(NGM),RG(NGM),SE(NSM,NBFM),XIG(NIM,NGM,NBM),
     '  XG1(NJM,NUM,NGM),XN(NJM,NGM),XPFP(*),XR(NJM,NGM),WG(NGM,NBM)
C  SMAR009 18/01/99 removed GK(NZ_GK_M),GQ(NZ_GQ_M),
      CHARACTER ERROR*(*)
      LOGICAL CALC_XR
!     Local Variables
      INTEGER MK,ng,nh2,nh2x,nhs,nj,nk1,nkk,nn,npnn,ns,nsnn
       !SMAR009 22/12/98 ,ny1,ny2,nz
      REAL*8 DGREEN,DGREENA,FACTOR2,GREEN,GREENA,SUM,SUMR,SUMU,
     '  SUMQ,XGC(3),XPC(3)
c      CHARACTER FORMAT*500

      CALL ENTERS('XEPGKGQ',*9999)

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
            SUM= 0.0d0
            SUMR=0.0d0
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              XR(nj,ng)=XG1(nj,1,ng)-XPFP(nj)
              SUM=SUM+XR(nj,ng)*XR(nj,ng)
              SUMR=SUMR+XR(nj,ng)*XN(nj,ng)
            ENDDO !nj
            RAD(ng)=DSQRT(SUM)
            DRDN(ng)=SUMR/RAD(ng)
          ENDDO !ng
        ELSE
C*** Transform into cartesian to find distance (RAD).
          DO ng=1,NGT(nb1j)
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
              SUMR=SUMR+XR(nj,ng)*XN(nj,ng)
            ENDDO !nj
            RAD(ng)=DSQRT(SUM)
            DRDN(ng)=SUMR/RAD(ng)
          ENDDO !ng
        ENDIF
        CALC_XR=.TRUE.
      ENDIF

C**** GK matrix components
      IF(GREENFLAG.NE.2) THEN
        nhs=0
C        nv=1 !temporary
C        ny1=NYNP(1,1,nh,np,1,1,nr) !eqtn number (row number) of GK.
        DO nh2x=1,NHE !Loop over dependent variable columns
          nh2=NH_LOC(nh2x,nx)
          ns=0
          DO nn=1,NNT(nbuhp) !Dependent variable loop
            npnn=NPNE(nn,nbuhp)
            DO nk1=1,NKT(nn,nbuhp)
              MK=NKHE(nk1,nn,nh) !MK is zero at a distorted node
              ns=ns+1
              nhs=nhs+1
C cbp 3/1/98 I think this should test on NKH(nh,npnn,1)
C              IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,np,1)-KTYP93(1,nr),1))
C     '          THEN
              IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,1)-KTYP93(1,nr),1))
     '          THEN
C                ny2=NYNP(MK,nv,nh2,npnn,2,1,nr) !local column # of GK
                SUMU=0.0d0
                DO ng=1,NGT(nbuh)
C                 Evaluate element integrals
                  IF(JTYP4.EQ.1) THEN !unsymmetric
                    SUMU=SUMU+PG_U(ns,1,ng)*
     '                DGREEN(IGREN(nr),0,nh,nh2,CE,DRDN(ng),RAD(ng),
     '                DET(nbuhp,NNMIN,ng),XN(1,ng),XR(1,ng))*
     '                RG(ng)*WG(ng,nbuh)
                  ELSE IF((JTYP4.EQ.2.OR.JTYP4.EQ.3).AND.
     '                ITYP10(nr).EQ.1) THEN
C                   cartesians,cylindrically sym. The Integral
C                   calculated after trans to polar coords. The
C                   additional trans introduces a polar "r" into the
C                   integrand which is given by
C                     RD = XG1(2,1,ng) if JTYP4=2 .
C                          XG1(1,1,ng) if JTYP4=3
                    SUMU=SUMU+PG_U(ns,1,ng)*RD(ng)*
     '                DGREENA(IGREN(nr),CE,XN(1,ng),XPFP(1),XPFP(2),
     '                XG1(1,1,ng),XG1(2,1,ng))*RG(ng)*WG(ng,nbuh)
                  ENDIF
                ENDDO !End of ng loop
C cpb 2/1/98 Adding BEM stiffness matrices
C                CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
C     '            NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C                GK(nz)=GK(nz)+SUMU*SE(ns,nbuhp)
                GKES(nhs,1)=GKES(nhs,1)+SUMU*SE(ns,nbuhp)
              ENDIF !End of MK > 0 loop
            ENDDO !End of nk1 loop
          ENDDO !End of nn loop
        ENDDO !End of nh2 loop
      ENDIF
C     GQ Loop
C      ny1=NYNP(1,1,nh,np,1,2,nr) !eqtn number (row number) of GQ.
      nhs=0
C      nv=1 !temporary
C     New meaning of nv 31-Oct-94
C     If NVJP(nv=0) > 1 Then we have a corner node. Considering the
C     list of surrounding ne's. If the current ne is the smallest in
C     in the list then nv=2. If it is the next highest then nv=3,
C     and if it is the highest (3D only) then nv=4.
C     Can avoid using NVJE for this case -> NVJE to be removed.
      DO nh2x=1,NHE !Loop over dependent variable columns
        nh2=NH_LOC(nh2x,nx)
        ns=0
        DO nn=1,NNT(nbqhp) !normal derivative loop
C Old        nv=NVJE(nn,nbqh,1,ne) !nv=1 unless at a corner or edge
          npnn=NPNE(nn,nbqhp)
          nsnn=ns
          DO nk1=1,NKT(nn,nbqhp)
            MK=NKHE(nk1,nn,nh) !MK is zero at a distorted node
            ns=ns+1
            nhs=nhs+1
C cpb 3/1/98 I think this should test on NKH(nh,npnn,2)
C            IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,np,2)-KTYP93(2,nr),1)) THEN
            IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,2)-KTYP93(2,nr),1))
     '        THEN
C              ny2=NYNP(MK,nv,nh2,npnn,2,2,nr) !local column # of GQ
              SUMQ=0.0D0
              DO ng=1,NGT(nbqh) !Evaluate element integrals
                IF(JTYP4.EQ.1) THEN !unsymmetric
                  IF(GREENFLAG.EQ.1) THEN
                    IF(NNMIN.EQ.1) THEN
                      SUMQ=SUMQ+PG_Q(ns,1,ng)*GREEN(IGREN(nr),
     '                  GREENFLAG,nh,nh2,CE,RAD(ng)/XIG(1,ng,nbqh),
     '                  DET(nbqhp,NNMIN,ng),XR(1,ng))*RG(ng)*WG(ng,nbqh)
                    ELSE IF(NNMIN.EQ.2) THEN
                      SUMQ=SUMQ+PG_Q(ns,1,ng)*GREEN(IGREN(nr),
     '                  GREENFLAG,nh,nh2,CE,RAD(ng)/(1.0d0-
     '                  XIG(1,ng,nbqh)),DET(nbqhp,NNMIN,ng),XR(1,ng))*
     '                  RG(ng)*WG(ng,nbqh)
                    ENDIF
                  ELSE
                    SUMQ=SUMQ+PG_Q(ns,1,ng)*GREEN(IGREN(nr),
     '                GREENFLAG,nh,nh2,CE,RAD(ng),DET(nbqhp,NNMIN,ng),
     '                XR(1,ng))*RG(ng)*WG(ng,nbqh)
                  ENDIF
                ELSE IF((JTYP4.EQ.2.OR.JTYP4.EQ.3).AND.ITYP10(nr).EQ.1)
     '            THEN
                  !cartesians,cylindrically sym.
                  !The Integral calculated after trans to polar coords
                  !The additional trans introduces a polar "r" into the
                  !integrand which is given by
                  !                   RD = XG1(2,1,ng) if JTYP4=2 .
                  !                        XG1(1,1,ng) if JTYP4=3
                  SUMQ=SUMQ+PG_Q(ns,1,ng)*RD(ng)*
     '                      GREENA(IGREN(nr),CE,XPFP(1),XPFP(2),
     '                      XG1(1,1,ng),XG1(2,1,ng))*RG(ng)*WG(ng,nbqh)
                ENDIF
              ENDDO !End ng loop
C cpb 2/1/98 Adding BEM stiffness matrices
C              CALL SPARSE(ny1,ny2,NYT(1,2,nx),nz,NZ_GQ_M,
C     '          NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,ERROR,*9999)
C              GQ(nz)=GQ(nz)+SUMQ*SE(ns,nbqhp)/FACTOR2
              GQES(nhs,1)=GQES(nhs,1)+SUMQ*SE(ns,nbqhp)/FACTOR2
C cpb 20/11/96 Adding bem curvature correction
              IF(BEMCURVATURECORRECTION.AND.(nk1.EQ.2.OR.nk1.EQ.3)) THEN
                IF(NIT(nb1jp).EQ.1) THEN
C                  ny2=NYNP(MK,nv,nh2,npnn,2,1,nr) !local column # of GK
C                  CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
C     '              NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C                  GK(nz)=GK(nz)-SUMQ*SE(ns,nbqhp)*CURVCORRECT(1,1,nn)
                  GKES(nhs,1)=GKES(nhs,1)-SUMQ*SE(ns,nbqhp)*
     '              CURVCORRECT(1,1,nn)
                ELSE
                  DO nkk=1,2
C                    ny2=NYNP(nkk+1,nv,nh2,npnn,2,1,nr) !loc col # of GK
C                    CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
C     '                NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C                    GK(nz)=GK(nz)-SUMQ*SE(ns,nbqhp)*
C     '                CURVCORRECT(nkk,nk1-1,nn)
                    GKES(nsnn+nkk+1,1)=GKES(nsnn+nkk+1,1)-SUMQ*
     '                SE(ns,nbqhp)*CURVCORRECT(nkk,nk1-1,nn)
                  ENDDO !nkk
                ENDIF
              ENDIF
            ENDIF !End of MK > 0 loop
          ENDDO !End of nk1 loop
        ENDDO !End of nn loop
      ENDDO !End of nh2 loop

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(XEPGKGQ_1)
        WRITE(OP_STRING,'('' Element contribution to'//
     '   ' Global matrices'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Node '',I5)') np
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        !GK matrix
C cpb 2/1/98 Adding BEM stiffness matrices
C        ny1=NYNP(1,1,nh,np,1,1,nr) !eqtn number (row number) of GK.
C        DO nn=1,NNT(nbuhp)
C          npnn=NPNE(nn,nbuhp)
C          FORMAT='('' GK('',I3,'','',I3,'')='',D12.4)'
C          DO nh2x=1,NHE !Loop over dependent variable columns
C            nh2=NH_LOC(nh2x,nx)
C            DO nk1=1,NKT(nn,nbuhp)
C              MK=NKE(nk1,nn,nbuhp)
C              IF(MK.GT.0.AND.
C     '          nk1.LE.MAX(NKH(nh,np,1)-KTYP93(1,nr),1)) THEN
C                ny2=NYNP(MK,1,nh2,npnn,2,1,nr) !local col number of GK
C                CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
C     '            NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
C                WRITE(OP_STRING,FORMAT)ny1,ny2,GK(nz)
C                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C              ENDIF
C            ENDDO
C          ENDDO !End of nh2 loop
C        ENDDO
        DO nkk=1,NDTOT+1
          WRITE(OP_STRING,'('' GKES(nhs,'',I1,'')='',4(1X,D12.5),'
     '      //'/:(13X,4(1X,D12.5)))') nkk,(GKES(nhs,nkk),
     '      nhs=-NDTOT-1,NHSTGK)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nkk
        !GQ matrix
C        ny1=NYNP(1,1,nh,np,1,2,nr) !eqtn number (row number) of GQ.
C        nv=1 !temporary
C        DO nn=1,NNT(nbqhp)
C          npnn=NPNE(nn,nbqhp)
C          FORMAT='('' GQ('',I3,'','',I3,'')='',D12.4)'
C          DO nh2x=1,NHE !Loop over dependent variable columns
C            nh2=NH_LOC(nh2x,nx)
C            DO nk1=1,NKT(nn,nbqhp)
C              MK=NKE(nk1,nn,nbqhp)
C              IF(MK.GT.0.AND.
C     '          nk1.LE.MAX(NKH(nh,np,2)-KTYP93(2,nr),1)) THEN
C                ny2=NYNP(MK,nv,nh2,npnn,2,2,nr)!local col number of GQ
C                CALL SPARSE(ny1,ny2,NYT(1,2,nx),nz,NZ_GQ_M,
C     '            NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,ERROR,*9999)
C                WRITE(OP_STRING,FORMAT)ny1,ny2,GQ(nz)
C                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C              ENDIF
C            ENDDO
C          ENDDO !End of nh2 loop
C        ENDDO
        DO nkk=1,NDTOT+1
          WRITE(OP_STRING,'('' GQES(nhs,'',I1,'')='',4(1X,D12.5),'
     '      //'/:(13X,4(1X,D12.5)))') nkk,(GQES(nhs,nkk),
     '      nhs=-NDTOT-1,NHSTGQ)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nkk
CC$OMP END CRITICAL(XEPGKGQ_1)
      ENDIF

      CALL EXITS('XEPGKGQ')
      RETURN
 9999 CALL ERRORS('XEPGKGQ',ERROR)
      CALL EXITS('XEPGKGQ')
      RETURN 1
      END


