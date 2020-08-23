      REAL*8 FUNCTION LENGTH_1D(NBJ,ne,NPNE,NVJE,XP)

C#### Function: LENGTH_1D
C###  Description:
C###    LENGTH_1D calculates the length of a 1D linear element

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),ne,NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
!     Local variables
      INTEGER nb,nj,nn,npnn(2),nvnn(2)
      REAL*8 length

      nb=NBJ(1,ne) !for geometry
      DO nn=1,NNT(nb)
        npnn(nn)=NPNE(nn,nb,ne)
        nvnn(nn)=NVJE(nn,nb,1,ne)
      ENDDO !nn
      length=0.d0
      DO nj=1,NJT
        length=length+(XP(1,nvnn(2),nj,npnn(2))-XP(1,nvnn(1),nj,
     &    npnn(1)))**2
      ENDDO !nj
      LENGTH_1D=DSQRT(length)
      
      RETURN
      END


