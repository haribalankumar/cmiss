      SUBROUTINE BMG2_SymStd_restrict( 
     &                       KF, KC, Q, QC, CI, 
     &                       Nx, Ny, Nxc, Nyc, IBC 
     &                       )

C ==========================================================================
C
C   BMG2_SymStd_restrict.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_restrict computes the restriction of a vector on the
C   fine grid, Q, to a vector on the coarse grid, QC.  The weights 
C   involve the transpose of the interpolation operator from the coarse 
C   grid to the fine grid.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2000/02/23 (JDM)
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
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  Nx, Nxc, Ny, Nyc

      INTEGER  IBC, KC, KF
      REAL*8   CI(Nxc,Nyc,8), Q(Nx,Ny), QC(Nxc,Nyc)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, ic, j, jc, kk, LOOPEND

C ==========================================================================

C ---------------------------------------------
C     Restrict the vector Q -> QC:
C ---------------------------------------------

C$OMP PARALLEL DO 
C$OMP& DEFAULT(NONE)
C$OMP& FIRSTPRIVATE(Nx,Ny,Nxc,Nyc)
C$OMP& PRIVATE(i,ic,j,jc,kk)
C$OMP& SHARED(CI,Q,QC)
	DO kk = 1,(Nxc-2)*(Nyc-2)
	   ic = mod(kk-1,Nxc-2)+2
	   i = 2*(ic-1)
	   jc = (kk-1)/(Nxc-2)+2
	   j = 2*(jc-1)
	   QC(ic,jc) = CI(ic,jc,LNE)*Q(i-1,j-1)
     &             + CI(ic,jc,LA)*Q(i,j-1)
     &             + CI(ic+1,jc,LNW)*Q(i+1,j-1)
     &             + CI(ic,jc,LR)*Q(i-1,j)
     &             + Q(i,j)
     &             + CI(ic+1,jc,LL)*Q(i+1,j)
     &             + CI(ic,jc+1,LSE)*Q(i-1,j+1)
     &             + CI(ic,jc+1,LB)*Q(i,j+1)
     &             + CI(ic+1,jc+1,LSW)*Q(i+1,j+1)
      ENDDO
C$OMP END PARALLEL DO

C ===========================================================================

       RETURN
       END
