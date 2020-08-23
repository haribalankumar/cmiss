      REAL*8 FUNCTION QUICKEST(nfv,NFVC,nh,NPNODE,nr,nvc,NVCNODE,
     '  NYNP,COUR,TC,TD,VC,XNFV,YP)

C#### Function: QUICKEST
C###  Type: REAL*8
C###  Description: Computes the Quadratic upstream interpolation for
C###    convective kinematics with estimated streaming terms.
C###    It assumes that flow is from velocity TC to TD.
C###    (Upwinding on velocity TC).
C###    See Leonard & Mokhtari, Int J Num Meth Eng, Vol 30,
C###    pp 729-766 (1990).

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'voro00.inc'
      REAL*8 CF
      PARAMETER (CF=0.166666666666666667D0)
!     Parameter List
      INTEGER nfv,NFVC(2,0:NFVCM,NVCM),nh,NPNODE(0:NP_R_M,0:NRM),
     '  nr,nvc,NVCNODE(2,NP_R_M),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 COUR,TC,TD,VC(0:NVCM),XNFV(-(NJM+1):NJM,NFVM),YP(NYM,NIYM)
!     Local Variables
      INTEGER nfvl,nfv2,cny,cnp,cnonode  !,cnv,ny
      REAL*8 DDOT,CURVN,DIST,TEMP        !,DMEDIAN

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

      QUICKEST=0.5d0*(TC+TD)-0.5d0*COUR*(TD-TC)-CF*CURVN

      RETURN
      END


