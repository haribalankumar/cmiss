      REAL*8 FUNCTION RADIUS_1D(NBJ,ne,NPNE,NVJE,XP)

C#### Function: RADIUS_1D
C###  Description:
C###    RADIUS_1D calculates the radius of a 1D linear element

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),ne,NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
!     Local variables
      INTEGER nb,nn,npnn(2),nvnn(2)

      nb=NBJ(nj_radius,ne) !basis function for radii
      DO nn=1,NNT(nb) !versions for radii at nodes
        nvnn(nn)=NVJE(nn,nb,nj_radius,ne)
        npnn(nn)=NPNE(nn,nb,ne)
      ENDDO !nn
      RADIUS_1D=0.5d0*(XP(1,nvnn(1),nj_radius,npnn(1))+
     &  XP(1,nvnn(2),nj_radius,npnn(2)))
      
      RETURN
      END


