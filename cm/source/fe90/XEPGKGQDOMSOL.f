      SUBROUTINE XEPGKGQDOMSOL(nb1j,nbqh,nbqhp,nbuh,nbuhp,nh,NHE,nhx,
     '  NKHE,NKH,NPNE,nr,nx,CE,DET,DRDN,PG_Q,PG_U,RAD,RD,RG,WG,XG1,XN,
     '  XPFP,XR,YD,ZE,ZF,ERROR,*)

C#### Subroutine: XEPGKGQDOMSOL
C###  Description:
C###    XEPGKGQDOMSOL contains code for finding a domain solution.
C###    It calculates the
C###    element integrals for element ne and given
C###    dependent variable number nh.
C**** Created by Martin Buist - June 2001

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
      INTEGER nb1j,nbqh,nbqhp,nbuh,nbuhp,nh,NHE,nhx,NKHE(NKM,NNM,NHM),
     '  NKH(NHM,NPM,NCM),NPNE(NNM,NBFM),nr,nx
      REAL*8 CE(NMM),DET(NBFM,0:NNM,NGM),DRDN(NGM),PG_Q(NSM,NUM,NGM),
     '  PG_U(NSM,NUM,NGM),RAD(NGM),RD(NGM),RG(NGM),XG1(NJM,NUM,NGM),
     '  XN(NJM,NGM),XPFP(*),XR(NJM,NGM),WG(NGM,NBM),YD(NHM),
     '  ZE(NSM,NHM),ZF(NSM,NHM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER MK,ng,nh2,nh2x,nhs,nj,nk1,nn,npnn,ns
      REAL*8 DGREEN,DGREENA,GREEN,GREENA,SUM,SUMR,SUMU,SUMQ,XGC(3),
     '  XPC(3)

      CALL ENTERS('XEPGKGQDOMSOL',*9999)

      IF(IGREN(nr).EQ.1) CE(1)=1.0d0
      IF(IGREN(nr).EQ.2) CE(1)=1.0d0

C***  Calculate gauss point dependent arrays
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
C***    Transform into cartesian to find distance (RAD).
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

C**** GK matrix components
      nhs=0
      DO nh2x=1,NHE !Loop over dependent variable columns
        nh2=NH_LOC(nh2x,nx)
        ns=0
        DO nn=1,NNT(nbuhp) !Dependent variable loop
          npnn=NPNE(nn,nbuhp)
          DO nk1=1,NKT(nn,nbuhp)
            MK=NKHE(nk1,nn,nh) !MK is zero at a distorted node
            ns=ns+1
            nhs=nhs+1
            IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,1)-KTYP93(1,nr),1))
     '        THEN
              SUMU=0.0d0
              DO ng=1,NGT(nbuh)
C               Evaluate element integrals
                IF(JTYP4.EQ.1) THEN !unsymmetric
                  SUMU=SUMU-PG_U(ns,1,ng)*DGREEN(IGREN(nr),0,nh,nh2,CE,
     '              DRDN(ng),RAD(ng),DET(nbuhp,0,ng),XN(1,ng),
     '              XR(1,ng))*RG(ng)*WG(ng,nbuh)
                ELSE IF((JTYP4.EQ.2.OR.JTYP4.EQ.3).AND.ITYP10(nr).EQ.1)
     '            THEN
C                 cartesians,cylindrically sym. The Integral
C                 calculated after trans to polar coords. The
C                 additional trans introduces a polar "r" into the
C                 integrand which is given by
C                   RD = XG1(2,1,ng) if JTYP4=2 .
C                        XG1(1,1,ng) if JTYP4=3
                  SUMU=SUMU-PG_U(ns,1,ng)*RD(ng)*DGREENA(IGREN(nr),CE,
     '              XN(1,ng),XPFP(1),XPFP(2),XG1(1,1,ng),XG1(2,1,ng))*
     '              RG(ng)*WG(ng,nbuh)
                ENDIF
              ENDDO !End of ng loop
              YD(nh)=YD(nh)+SUMU*ZE(ns,nhx)
            ENDIF !End of MK > 0 loop
          ENDDO !End of nk1 loop
        ENDDO !End of nn loop
      ENDDO !End of nh2 loop

C**** GQ matrix components
      nhs=0
      DO nh2x=1,NHE !Loop over dependent variable columns
        nh2=NH_LOC(nh2x,nx)
        ns=0
        DO nn=1,NNT(nbqhp) !normal derivative loop
          npnn=NPNE(nn,nbqhp)
          DO nk1=1,NKT(nn,nbqhp)
            MK=NKHE(nk1,nn,nh) !MK is zero at a distorted node
            ns=ns+1
            nhs=nhs+1
            IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,2)-KTYP93(2,nr),1))
     '        THEN
              SUMQ=0.0D0
              DO ng=1,NGT(nbqh) !Evaluate element integrals
                IF(JTYP4.EQ.1) THEN !unsymmetric
                  SUMQ=SUMQ+PG_Q(ns,1,ng)*GREEN(IGREN(nr),0,nh,nh2,CE,
     '              RAD(ng),DET(nbqhp,0,ng),XR(1,ng))*RG(ng)*
     '              WG(ng,nbqh)
                ELSE IF((JTYP4.EQ.2.OR.JTYP4.EQ.3).AND.ITYP10(nr).EQ.1)
     '            THEN
                  !cartesians,cylindrically sym.
                  !The Integral calculated after trans to polar coords
                  !The additional trans introduces a polar "r" into the
                  !integrand which is given by
                  !                   RD = XG1(2,1,ng) if JTYP4=2 .
                  !                        XG1(1,1,ng) if JTYP4=3
                  SUMQ=SUMQ+PG_Q(ns,1,ng)*RD(ng)*GREENA(IGREN(nr),CE,
     '              XPFP(1),XPFP(2),XG1(1,1,ng),XG1(2,1,ng))*RG(ng)*
     '              WG(ng,nbqh)
                ENDIF
              ENDDO !End ng loop
              YD(nh)=YD(nh)+SUMQ*ZF(ns,nh)
            ENDIF !End of MK > 0 loop
          ENDDO !End of nk1 loop
        ENDDO !End of nn loop
      ENDDO !End of nh2 loop

      CALL EXITS('XEPGKGQDOMSOL')
      RETURN
 9999 CALL ERRORS('XEPGKGQDOMSOL',ERROR)
      CALL EXITS('XEPGKGQDOMSOL')
      RETURN 1
      END


