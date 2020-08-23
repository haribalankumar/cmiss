      SUBROUTINE BMG3_SymStd_GET_pointers( 
     &               k, IGRD, NOG,
     &               p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz )

C ===========================================================================

C***BEGIN PROLOGUE BMG3_SymStd_GET_pointers
C
C***PURPOSE  
C
C     BMG3_SymStd_GET_pointers is a subroutine provided to access 
C     the data structures used in BMG3D.
C
C***DESCRIPTION
C
C     This routine fetches the parameters for grid K. This
C     includes the number of points in the x, y, and directions,
C     the x, y, and z grid spacings, and pointers
C     to the start of storage of arrays in BMG3D.
C
C***PARAMETERS
C
C***INPUT
C
C   k       The grid number.
C   NOG     The number of grids and leading dimension of the
C           workspace array IGRD.
C
C***OUTPUT
C
C   Nx      The number of grid points in the x direction for grid k, 
C           including two fictitious points.
C   Ny      The number of grid points in the y directiom for grid K, 
C           including two fictitious points.
C   Nz      The number of grid points in the z direction for grid K,
C           including two fictitious points
C
C   p_CI    Pointer to the interpolation on grid k, in CI.
C   p_SO    Pointer to the stencil on grid k, in SO.
C   p_SOR   Pointer to the stencil workspace on grid k, in SOR.
C   p_U     Pointer to vectors on grid k, in Q and QF
C
C WORK ARRAYS
C
C   IGRD      An integer work array.
C             It contains indexing information for the different
C             grids. IGRD is set up by MGGRD3.
C
C***REFERENCES  "the multigrid method for the diffusion equation
C                 with strongly discontinuous coefficients"
C                 r. e. alcouffe, achi brandt, j. e. dendy, jr.
C                 and j. w. painter, siam j. sci. stat. comput.,
C                 vol. 2, no. 4, dec. 1981, pp. 430-454.
C               "black box multigrid", j. e. dendy, jr., j. comp. phys.,
C                 vol. 48, no. 3, dec. 1982, pp. 366-386.
C
C***ROUTINES CALLED (NONE)
C
C***REVISION HISTORY  (YYYYMMDD)
C
C   1986/12/09  - WRITTEN (JED)
C               - based on bmgk3.f from the 2D black box code.
C   2000/03/04  - Cleaned up declarations, IMPLICIT NONE (JDM)
C   2000/03/05  - Updated comments (JDM)
C
C***END PROLOGUE BMG3_SymStd_GET_pointers
C
C ===========================================================================

      IMPLICIT NONE

C ---------------------------
C    Includes:
C
      INCLUDE 'BMG_workspace.h'

C ---------------------------
C    Argument Declarations:
C     
      INTEGER k, NOG, Nx, Ny, Nz

      INTEGER IGRD(NOG,NBMG_pIGRD), p_CI, p_SO, p_SOR, p_U

C ===========================================================================

      Nx = igrd(k,idL_BMG_Nx)
      Ny = igrd(k,idL_BMG_Ny)
      Nz = igrd(k,idL_BMG_Nz)

      p_U   = IGRD(k,ipL_BMG_U)
      p_SO  = IGRD(k,ipL_BMG_SO)
      p_SOR = IGRD(k,ipL_BMG_SOR)
      p_CI  = IGRD(k,ipL_BMG_CI)

C ===========================================================================

      RETURN
      END

