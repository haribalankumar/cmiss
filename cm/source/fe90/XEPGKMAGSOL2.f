      SUBROUTINE XEPGKMAGSOL2(nb1j,nbuh,nbuhp,nh,NHE,nhx,NKHE,NKH,NPNE,
     '  nr,DSIGMA,PG_U,RAD,RG,WG,XG1,XPFP,XR,YD_LOCAL,ZE,ERROR,*)

C#### Subroutine: XEPGKMAGSOL2
C###  Description:
C###    XEPGKMAGSOL2 contains code for finding a magnetic vector
C###    potential domain solution.
C###    It calculates the
C###    element integrals for element ne and given
C###    dependent variable number nh.
C**** Created by Martin Buist - December 2001

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER nb1j,nbuh,nbuhp,nh,NHE,nhx,NKHE(NKM,NNM,NHM),
     '  NKH(NHM,NPM,NCM),NPNE(NNM,NBFM),nr
      REAL*8 DSIGMA,PG_U(NSM,NUM,NGM),RAD(NGM),RG(NGM),WG(NGM,NBM),
     '  XG1(NJM,NUM,NGM),XPFP(*),XR(NJM,NGM),
     '  YD_LOCAL(9,9),ZE(NSM,NHM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER MK,ng,nhs,nh2x,nj,nk1,nn,npnn,ns
      REAL*8 GREEN_LOC,GREENSTATIC,SUM,XGC(3),XPC(3)

      CALL ENTERS('XEPGKMAGSOL2',*9999)

C***  Calculate gauss point dependent arrays
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
          ENDDO !nj
          RAD(ng)=DSQRT(SUM)
        ENDDO !ng
      ENDIF

C**** Potential field integration
      nhs=0
      GREENSTATIC=DSIGMA/(4.0d0*PI)
      DO nh2x=1,NHE !Loop over dependent variable columns
        ns=0
        DO nn=1,NNT(nbuhp) !Dependent variable loop
          npnn=NPNE(nn,nbuhp)
          DO nk1=1,NKT(nn,nbuhp)
            MK=NKHE(nk1,nn,nh) !MK is zero at a distorted node
            ns=ns+1
            nhs=nhs+1
            IF(MK.GT.0.AND.nk1.LE.MAX(NKH(nh,npnn,1)-KTYP93(1,nr),1))
     '        THEN
              SUM=0.0d0
              DO ng=1,NGT(nbuh)
C***            Evaluate element integrals
                IF(JTYP4.EQ.1) THEN !unsymmetric
                  !GREEN gives dsigma/(4*pi*r)
                  GREEN_LOC=GREENSTATIC/RAD(ng)
                  SUM=SUM+PG_U(ns,1,ng)*GREEN_LOC*RG(ng)*WG(ng,nbuh)
                ELSE IF(JTYP4.EQ.2.OR.JTYP4.EQ.3) THEN
                  ERROR='>>Not implemented'
                  GOTO 9999
                ENDIF
              ENDDO !End of ng loop
              YD_LOCAL(nh,1)=YD_LOCAL(nh,1)-SUM*ZE(ns,nhx)
              YD_LOCAL(nh,2)=YD_LOCAL(nh,2)-SUM*ZE(ns,nhx)
              YD_LOCAL(nh,3)=YD_LOCAL(nh,3)-SUM*ZE(ns,nhx)
            ENDIF !End of MK > 0 loop
          ENDDO !End of nk1 loop
        ENDDO !End of nn loop
      ENDDO !End of nh2 loop

      CALL EXITS('XEPGKMAGSOL2')
      RETURN
 9999 CALL ERRORS('XEPGKMAGSOL2',ERROR)
      CALL EXITS('XEPGKMAGSOL2')
      RETURN 1
      END


