      SUBROUTINE BMG_CONVERT_csr_matrix(
     &              Nd, Nx, Ny, Nz, 
     &              BMG_iPARMS, BMG_rWORK, BMG_pWORK, BMG_BDY,  
     &              BMG_TO_CM, AA, IA, JA, NStncl
     &              )

C =======================================================================
C 
C   BMG_CONVERT_csr_matrix
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   This subroutine converts a matrix stored in CSR format to a matrix
C   stored in the Dendy-Cell-Symmetric format for structured grids.  It
C   is assumed that the matrix has been generated from a lexicographical
C   construction so that the unknowns may be accessed via (ijk) entries.
C
C   AUTHOR  AUSTIN, TRAVIS
C           BIOENGINEERING INSTITUTE
C           UNIVERSITY OF AUCKLAND
C           NEW ZEALAND
C           E-MAIL: t.austin@auckland.ac.nz
C
C   Written:    2005/01/22 (TMA)
C
C =======================================================================
C   INPUT:
C ========================
C
C   ----------------------
C   Fine Grid Dimensions:
C   ---------------------- 
C
C   Nd       Dimension of problem
C   Nx       Number of points in x-direction (excluding ghost points)
C   Ny       Number of points in y-direction (excluding ghost points)
C   Nz       Number of points in z-direction (excluding ghost points)
C
C   ----------------------
C   CSR Matrix
C   ---------------------- 
C
C   AA      Nonzeros in matrix stored in CSR format.
C   IA      Pointers to rows in matrix A.
C   JA      Column number for corresponding entry in AA 
C
C   Note:  Size(AA) = Size(JA)
C
C   -----------------------
C   BOXMG Parameter Arrays:
C   -----------------------
C
C   BMG_pWORK      Array containing pointers to locations in rWORK and iWORK.
C
C =======================================================================
C   OUTPUT:
C ========================
C
C   -----------------------
C   BOXMG Parameter Arrays:
C   -----------------------
C
C   BMG_rWORK      Arrays containing solution, rhs, matrix, residual, etc.
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_stencils.h'

C ========================================================================== 

C ---------------- 
C     Arguments  
C
      INTEGER  Nd, Nx, Ny, Nz, NStncl 
      INTEGER  BMG_BDY(*), BMG_iPARMS(*), BMG_pWORK(*), BMG_TO_CM(*), 
     &         IA(*), JA(*)
      REAL*8   AA(*), BMG_rWORK(*)

C ---------------- 
C     Local Arguments  
C
      INTEGER IFD

C ========================================================================== 


      IF( Nd.EQ.2 ) THEN

         CALL BMG2_CONVERT_csr_matrix( 
     &           Nx,Ny,
     &           BMG_rWORK(BMG_pWORK(ip_SO)),BMG_BDY,BMG_iPARMS,
     &           BMG_TO_CM,AA,IA,JA,NStncl
     &           )

      ELSEIF( Nd.EQ.3 ) THEN

         CALL BMG3_CONVERT_csr_matrix( 
     &           Nx,Ny,Nz,
     &           BMG_rWORK(BMG_pWORK(ip_SO)),BMG_BDY,BMG_iPARMS,
     &           BMG_TO_CM,AA,IA,JA,NStncl
     &           )

      ENDIF


      RETURN
      END
