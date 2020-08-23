      SUBROUTINE BMG3_SymStd_READ_matrix( 
     &                BMG_IOFLAG, SO, Nx, Ny, Nz, IFD, NStncl, kg, NOG
     &                )

C ==========================================================================
C
C   BMG3_SymStd_DUMP_stencil.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_DUMP_stencil.f outputs the stencil interactively
C   to the screen or dumps the entire stencil into an ASCII text
C   file. Thus it is intended for debugging the stencil generation
C   on reasonable size grids.
C
C   - assumes that SO(p_SO) has been used in the call, where
C     p_SO is returned from BMG3key.
C   - assumes that you have provided the value of NStncl
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2003/01/30 (JDM)
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

C ---------------------------
C     Argument Declarations:
C
      INTEGER IFD, kg, NOG, NStncl, Nx, Ny, Nz
      REAL*8  SO(Nx*Ny*Nz*Nstncl)
      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)

C --------------------------
C     Local Declarations:
C
      INTEGER   i, j, k, lstncl
      CHARACTER ANS*1, FILEo*30

C ==========================================================================


      FILEo = 'output/'
     &     //'Matrix.'//CHAR(ICHAR('0')+kg)      

      OPEN (10,FILE=FILEo,STATUS='UNKNOWN',FORM='UNFORMATTED') 

      READ(10) Nx, Ny, Nz, NStncl
      READ(10) (SO(i),i=1,NStncl*Nx*Ny*Nz)

      CLOSE(10)

      RETURN
      END
         
