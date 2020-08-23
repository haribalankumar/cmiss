      SUBROUTINE BMG2_SymStd_SETUP_lines_x(
     &                SO, SOR, Nx, Ny, NStncl 
     &                )

C ==========================================================================
C
C   BMG2_SymStd_SETUP_lines_x.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_SETUP_lines_x.f performs a factorization of the tridiagonal
C   matrix that arises in x-line relaxation.  It assumes that the system
C   is diagonally dominant and it works directly with the stencil.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2000/02/23 (JDM)
C   - cut and paste from mgcoef.f
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
      INTEGER   Nx, Ny, NStncl
      REAL*8    SO(Nx,Ny,NStncl), SOR(Nx,Ny,2)

C ----------------------------
C     Local Declarations
C
      INTEGER i, j, k, INFO

C ==========================================================================

C$OMP PARALLEL DO
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k),
C$OMP& SHARED(Nx,Ny,SO,SOR)
	DO k = 1, (Nx-2)*(Ny-2)
	   i = mod(k-1,Nx-2)+2
	   j = (k-1)/(Nx-2)+2
         SOR(i,j,2) = -SO(i,j,KW)  ! off diagonal
         SOR(i,j,1) =  SO(i,j,KO)  ! diagonal
	ENDDO
C$OMP END PARALLEL DO

C     calculate the L*D*L' factorizations for the lines in x-direction

C$OMP PARALLEL DO
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(INFO,j),
C$OMP& SHARED(Nx,Ny,SOR)
      DO j=2,Ny-1
         CALL DPTTRF (Nx-2, SOR(2,j,1), SOR(3,j,2), INFO)
         IF (INFO .NE. 0) THEN
            WRITE(*,*) 'SETUP lines - x INFO = ',INFO
         ENDIF
      ENDDO
C$OMP END PARALLEL DO

c      DO j=2, Ny-1
c         SOR(2,j,MSOR)=rONE/SO(2,j,KO)
c         DO i=3, Nx-1
c            SOR(i,j,MSOR)=rONE
c     &                   /( SO(i,j,KO) - SOR(i-1,j,MSOR)*SO(i,j,KW)**2 )
c         ENDDO
c      ENDDO

C ==========================================================================

      RETURN
      END

