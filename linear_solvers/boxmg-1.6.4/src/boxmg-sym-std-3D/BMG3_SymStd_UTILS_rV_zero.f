
      subroutine BMG3_SymStd_UTILS_rV_zero(q,ii,jj,kk)

C***  BEGIN PROLOGUE  BMG3_SymStd_UTILS_rV_zero
C***  DESCRIPTION
C     BMG3_SymStd_UTILS_rV_zero zeroes the Q array on grid KG.
C***  PARAMETERS
C***  INPUT
C     KG        KG is the grid number.
C     II        Number of grid ponts in the x direction,
C     including two fictitious points.
C     JJ        Number of grid points in the y direction,
C     including two fictitious points.
C     KK        Number of grid points in the z direction,
C     including two fictitious points.
C***  OUTPUT
C     Q         Refer to BMG3D.
C***  ROUTINES CALLED (NONE)
C***  CHANGES
C  991207       code indented, do loops end on enddo, by M.Berndt
C***  END PROLOGUE  BMG3_SymStd_UTILS_rV_zero
      IMPLICIT NONE
C     CALLING ARGUMENTS
      integer ii, jj, kk, kg
      real*8 q(ii,jj,kk)
C     LOCAL VARIABLES
      integer i, j, k, kl, maxIJ, maxIJK
      real*8 rZERO
C***  FIRST EXECUTABLE STATEMENT  MGOP3
      rZERO = 0

      maxIJK = ii*jj*kk - 1
      maxIJ  = ii*jj

C$OMP PARALLEL DO 
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k,kl)
C$OMP& SHARED(ii,jj,kk,maxIJ,maxIJK,q,rZERO)
      DO kl = 0,maxIJK-1
         i = mod(kl,ii)+1
         j = mod(kl/ii,jj)+1
         k = (kl/maxIJ)+1
         q(i,j,k)= rZERO
      ENDDO
C$OMP END PARALLEL DO

      RETURN
      END
