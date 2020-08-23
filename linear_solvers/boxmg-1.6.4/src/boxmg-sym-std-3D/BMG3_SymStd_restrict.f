      SUBROUTINE BMG3_SymStd_restrict(
     &                KFG, KCG, 
     &                Q, QC, CI, Nx, Ny, Nz, Nxc, Nyc, Nzc 
     &                )

C ==========================================================================
C
C   BMG3_SymStd_restrict.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_restrict computes the restriction of a vector on the
C   fine grid, Q, to a vector on the coarse grid, QC.  The weights 
C   involve the transpose of the interpolation operator from the coarse 
C   grid to the fine grid.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C
C   2000/03/06  - written (JDM)
C               - the core computation was cut from mgrcl3.f
C
C ==================================================================
C   INPUT:
C ========================
C
C
C
C ==================================================================
C   OUTPUT:
C ===========================
C
C
C
C ==================================================================
C   LOCAL:
C ========================
C
C
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  Nx, Nxc, Ny, Nyc, Nz, Nzc

      INTEGER  KCG, KFG
      REAL*8   CI(Nxc,Nyc,Nzc,26), Q(Nx,Ny,Nz), QC(Nxc,Nyc,Nzc)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, ic, j, jc, k, kc, kl, maxXY, maxXYZ, Nxc2, Nyc2, Nzc2

C ==========================================================================

C ---------------------------------------------
C     Restrict the vector Q -> QC:
C ---------------------------------------------

      Nxc2 = Nxc-2
      Nyc2 = Nyc-2
      Nzc2 = Nzc-2

      maxXYZ = Nxc2*Nyc2*Nzc2
      maxXY  = Nxc2*Nyc2

C$OMP PARALLEL DO 
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,ic,j,jc,k,kc)
C$OMP& SHARED(maxXY,maxXYZ,Nxc2,Nyc2,Nzc2,CI,Q,QC)
      DO kl = 0,maxXYZ-1
         
         ic = MOD(kl,Nxc2)+2
         i = 2*(ic-1)
         jc = MOD(kl/Nxc2,Nyc2)+2
         j = 2*(jc-1)
         kc = (kl/maxXY)+2
         k = 2*(kc-1)

         QC(ic,jc,kc) = CI(ic,jc,kc,lxyne)*Q(i-1,j-1,k)
     &        + CI(ic,jc,kc,lxya)*Q(i,j-1,k)
     &        + CI(ic+1,jc,kc,lxynw)*Q(i+1,j-1,k)
     &        + CI(ic,jc,kc,lxyr)*Q(i-1,j,k)
     &        + Q(i,j,k)
     &        + CI(ic+1,jc,kc,lxyl)*Q(i+1,j,k)
     &        + CI(ic,jc+1,kc,lxyse)*Q(i-1,j+1,k)
     &        + CI(ic,jc+1,kc,lxyb)*Q(i,j+1,k)
     &        + CI(ic+1,jc+1,kc,lxysw)*Q(i+1,j+1,k)
     &        + CI(ic,jc,kc,ltne)*Q(i-1,j-1,k-1)
     &        + CI(ic,jc,kc,lyznw)*Q(i,j-1,k-1)
     &        + CI(ic+1,jc,kc,ltnw)*Q(i+1,j-1,k-1)
     &        + CI(ic,jc,kc,lxzne)*Q(i-1,j,k-1)
     &        + CI(ic,jc,kc,lxza)*Q(i,j,k-1)
     &        + CI(ic+1,jc,kc,lxznw)*Q(i+1,j,k-1)
     &        + CI(ic,jc+1,kc,ltse)*Q(i-1,j+1,k-1)
     &        + CI(ic,jc+1,kc,lyzne)*Q(i,j+1,k-1)
     &        + CI(ic+1,jc+1,kc,ltsw)*Q(i+1,j+1,k-1)
     &        + CI(ic,jc,kc+1,lbne)*Q(i-1,j-1,k+1)
     &        + CI(ic,jc,kc+1,lyzsw)*Q(i,j-1,k+1)
     &        + CI(ic+1,jc,kc+1,lbnw)*Q(i+1,j-1,k+1)
     &        + CI(ic,jc,kc+1,lxzse)*Q(i-1,j,k+1)
     &        + CI(ic,jc,kc+1,lxzb)*Q(i,j,k+1)
     &        + CI(ic+1,jc,kc+1,lxzsw)*Q(i+1,j,k+1)
     &        + CI(ic,jc+1,kc+1,lbse)*Q(i-1,j+1,k+1)
     &        + CI(ic,jc+1,kc+1,lyzse)*Q(i,j+1,k+1)
     &        + CI(ic+1,jc+1,kc+1,lbsw)*Q(i+1,j+1,k+1)
              
      ENDDO
C$OMP END PARALLEL DO
         
C ==========================================================================

      return
      end






