      SUBROUTINE XEPGKMAGSOL(nb1j,nbuh,nbuhp,ne,nh,NHE,nhx,NKHE,NKH,
     &  nnmin,NPNE,nr,NW,DETLOC,DSIGMA,PG_U,RAD,RG,WG,XG1,XN,XPFP,XR,
     &  YD_LOCAL,ZE,ERROR,*)

C#### Subroutine: XEPGKMAGSOL
C###  Description:
C###    XEPGKMAGSOL contains code for finding a magnetic
C###    field intensity domain solution.
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
      INTEGER nb1j,nbuh,nbuhp,ne,nh,NHE,nhx,NKHE(NKM,NNM,NHM),
     '  NKH(NHM,NPM,NCM),nnmin,NPNE(NNM,NBFM),nr,NW(NEM,3)
      REAL*8 DETLOC(NBFM,0:NNM,NGM),DSIGMA,PG_U(NSM,NUM,NGM),RAD(NGM),
     '  RG(NGM),WG(NGM,NBM),XG1(NJM,NUM,NGM),XN(NJM,NGM),XPFP(*),
     '  XR(NJM,NGM),YD_LOCAL(9,9),ZE(NSM,NHM)
      CHARACTER ERROR*(*)
      
!     Local Variables
      INTEGER MK,ng,nhs,nh2x,nj,nk1,nn,npnn,ns
      REAL*8 DGREEN_LOC,DGREENSTATIC,SUM,SUMU1,SUMU2,SUMU3,
     '  TERM1,TERM2,TERM3,XGC(3),XPC(3)

      CALL ENTERS('XEPGKMAGSOL',*9999)

C***  Calculate gauss point dependent arrays
      IF(ITYP10(nr).EQ.1) THEN
        DO ng=1,NGT(nb1j)
          SUM=0.0d0
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)

C LKC Committed 14-JAN-2011
C
C LKC 27-JAN-2009 calculate the distance between sensor and g.p.
C position. Previously it was calculating the distance between g.p. and
C sensor (the opposite sign).            
C            XR(nj,ng)=XG1(nj,1,ng)-XPFP(nj)
            XR(nj,ng)=XPFP(nj)-XG1(nj,1,ng)            

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
      DGREENSTATIC=DSIGMA/(4.0d0*PI)
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
              SUMU1=0.0d0
              SUMU2=0.0d0
              SUMU3=0.0d0
              DO ng=1,NGT(nbuh)
C               Evaluate element integrals
                IF(JTYP4.EQ.1) THEN !unsymmetric
                  !DGREEN gives -dsigma/(4*pi*r^2) * (1/r)
                  DGREEN_LOC=DGREENSTATIC*DETLOC(nbuhp,NNMIN,ng)/
     '              (RAD(ng)*RAD(ng)*RAD(ng))

C LKC Committed 14-JAN-2011
C
C LKC 27-JAN-2009 - was previously calculating XR x XN
C                  (should really be XN x NR)
C                  
C                  TERM1=DGREEN_LOC*(XR(2,ng)*XN(3)-XR(3,ng)*XN(2)
c                  TERM2=DGREEN_LOC*(XR(3,ng)*XN(1)-XR(1,ng)*XN(3)
C                  TERM3=DGREEN_LOC*(XR(1,ng)*XN(2)-XR(2,ng)*XN(1)

C LKC Committed 14-JAN-2011
C                  
C LKC 27-JAN-2009 -- Previously XN was calculated (incorrectly) in
C MAGSOL. Calculate XN here.
                  CALL NORMAL(ne,nr,NW,XG1(1,1,ng),XN(1,ng),.FALSE.,
     '              ERROR,*9999)
C LKC 27-JAN-2009 -- (sigma_in - sigma_out)/4pi * (XN x XR)
                  TERM1=DGREEN_LOC*(XN(2,ng)*XR(3,ng)-XN(3,ng)*XR(2,ng))
                  TERM2=DGREEN_LOC*(XN(3,ng)*XR(1,ng)-XN(1,ng)*XR(3,ng))
                  TERM3=DGREEN_LOC*(XN(1,ng)*XR(2,ng)-XN(2,ng)*XR(1,ng))

                  SUMU1=SUMU1+(PG_U(ns,1,ng)*RG(ng)*WG(ng,nbuh)*TERM1)
                  SUMU2=SUMU2+(PG_U(ns,1,ng)*RG(ng)*WG(ng,nbuh)*TERM2)
                  SUMU3=SUMU3+(PG_U(ns,1,ng)*RG(ng)*WG(ng,nbuh)*TERM3)
                ELSE IF(JTYP4.EQ.2.OR.JTYP4.EQ.3) THEN
                  ERROR='>>Not implemented'
                  GOTO 9999
                ENDIF
              ENDDO !End of ng loop
              YD_LOCAL(nh,1)=YD_LOCAL(nh,1)+SUMU1*ZE(ns,nhx)
              YD_LOCAL(nh,2)=YD_LOCAL(nh,2)+SUMU2*ZE(ns,nhx)
              YD_LOCAL(nh,3)=YD_LOCAL(nh,3)+SUMU3*ZE(ns,nhx)
            ENDIF !End of MK > 0 loop
          ENDDO !End of nk1 loop
        ENDDO !End of nn loop
      ENDDO !End of nh2 loop

      CALL EXITS('XEPGKMAGSOL')
      RETURN
 9999 CALL ERRORS('XEPGKMAGSOL',ERROR)
      CALL EXITS('XEPGKMAGSOL')
      RETURN 1
      END


