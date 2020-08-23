      SUBROUTINE BMG2_SymStd_SETUP_recip( 
     &                SO, SOR, Nx, Ny, NStncl, NSOR_v
     &                )

C ==========================================================================
C
C   BMG2_SymStd_cons_interp.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_cons_recip.f computes the reciprocal of the 
C   central stencil coefficient.  
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2000/02/28 (JDM)
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

C ---------------------------
C    Argument Declarations:
C
      INTEGER   NSOR_v, NStncl, Nx, Ny
      REAL*8    SO(Nx,Ny,NStncl), SOR(Nx,Ny,NSOR_v)

C --------------------------
C     Local Declarations:
C
      INTEGER   i, j, k

	INTEGER omp_get_thread_num
	EXTERNAL omp_get_thread_num

C ==========================================================================

C$OMP PARALLEL DO
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k),
C$OMP& SHARED(Nx,Ny,SO,SOR)
	DO k = 1, (Nx-2)*(Ny-2)
	   i = mod(k-1,Nx-2)+2
	   j = (k-1)/(Nx-2)+2
	   SOR(i,j,msor)=rONE/SO(i,j,ko)
      ENDDO
C$OMP END PARALLEL DO

C ==========================================================================

C =======================

         RETURN
         END
