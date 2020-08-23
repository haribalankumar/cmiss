      SUBROUTINE BMG2_SymStd_SETUP_lines_y(
     &                SO, SOR, Nx, Ny, NStncl
     &                )

C ==========================================================================
C
C   BMG2_SymStd_SETUP_lines_y.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_SETUP_lines_y.f performs a factorization of the tridiagonal
C   matrix that arises in y-line relaxation.  It assumes that the system
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
      REAL*8    SO(Nx,Ny,NStncl), SOR(Ny,Nx,2)

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
         SOR(j,i,2) = -SO(i,j,KS)  ! off diagonal
         SOR(j,i,1) =  SO(i,j,KO)  ! diagonal
      ENDDO
C$OMP END PARALLEL DO

C     calculate the L*D*L' factorizations for the lines in y-direction

C$OMP PARALLEL DO
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(INFO,j),
C$OMP& SHARED(Nx,Ny,SOR)
      DO i=2, Nx-1
         CALL DPTTRF (Ny-2, SOR(2,i,1), SOR(3,i,2), INFO)
         IF (INFO .NE. 0) THEN
            write(*,*) 'SETUP lines - y INFO = ',INFO
         ENDIF
      ENDDO
C$OMP END PARALLEL DO

c      DO i=2, Nx-1
c         SOR(i,2,MSOS)=rONE/SO(i,2,KO)
c         DO j=3, Ny-1
c            SOR(i,j,MSOS)=rONE
c     &                   /( SO(i,j,KO) - SOR(i,j-1,MSOS)*SO(i,j,KS)**2 )
c         ENDDO
c      ENDDO

C ==========================================================================

      RETURN
      END

