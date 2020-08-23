      REAL*8 FUNCTION ULTRAQUICK(nfv,NFVC,nh,NPNODE,nr,nvc,NVCNODE,
     '  NYNP,TC,TD,VC,XNFV,YP)

C#### Function: ULTRAQUICK
C###  Type: REAL*8
C###  Description: Computes the Quadratic upstream interpolation for
C###    convective kinematics with an ULTRA (Universal limiter for
C###    tight resolution and accuracy) limiter. It assumes that flow
C###    is from velocity TC to TD. (Upwinding on velocity TC).
C###    See Leonard & Mokhtari, Int J Num Meth Eng, Vol 30,
C###    pp 729-766 (1990).

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'voro00.inc'
!     Parameter List
      REAL*8 BIG
      PARAMETER (BIG=100.d0)
      INTEGER nfv,NFVC(2,0:NFVCM,NVCM),nh,NPNODE(0:NP_R_M,0:NRM),
     '  nr,nvc,NVCNODE(2,NP_R_M),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 TC,TD,VC(0:NVCM),XNFV(-(NJM+1):NJM,NFVM),YP(NYM,NIYM)
!     Local Variables
      INTEGER nfvl,cnp,nfv2,cny,cnonode !,cnv,ny
      REAL*8 DDOT,CURVN,DIST,TEMP,TU,T3,T4,T5,DMEDIAN,CF

      IF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN
        CF=0.125d0
      ELSEIF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN
        CF=0.166666666666666666666667d0
      ENDIF

      CURVN=0.d0
      DO nfvl=1,NFVC(1,0,nvc)
        cnonode=NFVC(1,nfvl,nvc)

        IF(NVCNODE(TYPE,cnonode).NE.BOUNDARY) THEN
          cnp=NPNODE(cnonode,nr)
          cny=NYNP(1,1,nh,cnp,0,1,nr)
          nfv2=NFVC(2,nfvl,nvc)

C         ..Compute normal component of curvature
          TEMP=(YP(cny,8)-TC)*XNFV(FAREA,nfv2)*XNFV(IDIST,nfv2)
          CURVN=CURVN+TEMP*DABS(DDOT(NJ_LOC(NJL_GEOM,0,nr),
     '      XNFV(1,nfv2),1,XNFV(1,nfv),1))

        ENDIF
      ENDDO

      CURVN=CURVN/VC(nvc)
      DIST=1.d0/XNFV(IDIST,nfv)
      CURVN=CURVN*DIST**2

      TU=CURVN+2.d0*TC-TD
      T3=TU+BIG*(TC-TU)
      T4=DMEDIAN(TC,TD,T3)
      T5=0.5d0*(TC+TD)-CF*CURVN

      ULTRAQUICK=DMEDIAN(TC,T4,T5)

      RETURN
      END

