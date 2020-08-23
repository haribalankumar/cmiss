      REAL*8 FUNCTION RADIUS_MEAN(NBJ,ne,NPNE,NVJE,XP)

C#### Function: RADIUS_MEAN
C###  Description:
C###    RADIUS_MEAN evaluates the radius of 1D airway tree by
C###    integrating between Xi= 0.25 to Xi=0.75.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),ne,NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
!     Local variables
      INTEGER nb,np1,np2,nv1,nv2

      nb=NBJ(nj_radius,ne)
      np1=NPNE(1,nb,ne)
      np2=NPNE(2,nb,ne)
      nv1=NVJE(1,nb,nj_radius,ne)
      nv2=NVJE(2,nb,nj_radius,ne)

      RADIUS_MEAN=0.5d0*(XP(1,nv1,nj_radius,np1)+
     &  XP(1,nv2,nj_radius,np2))+11.d0/96.d0*(XP(2,nv1,nj_radius,np1)
     &  -XP(2,nv2,nj_radius,np2))
      
      RETURN
      END


