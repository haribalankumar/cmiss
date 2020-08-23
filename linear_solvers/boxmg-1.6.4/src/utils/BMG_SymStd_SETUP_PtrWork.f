      SUBROUTINE BMG_SymStd_SETUP_PtrWork( 
     &                Nd, Nx, Ny, Nz, 
     &                BMG_iPARMS, 
     &                NOGm, NFm, NSOm, 
     &                NBMG_iWORKm, NBMG_rWORKm,
     &                NBMG_iWORK_PLm, NBMG_rWORK_PLm,
     &                BMG_pWORK, BMG_InWORK, pSR, pSI,
     &                BMG_IOFLAG
     &                )

C =======================================================================
C 
C   BMG_SymStd_SETUP_PtrWork
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG_SymStd_SETUP_PtrWork calls either BMG2_SymStd_SETUP_PtrWork or it
C   calls BMG3_SymStd_SETUP_PtrWork depending on ND.
C
C   Written:    2004/12/20 (TMA)
C
C =======================================================================
C   PARAMETERS:
C ========================
C
C   pNx   - index of x-direction dimensions
C   pNy   - index of y-direction dimensions
C   pNz   - index of z-direction dimensions
C   pN3   - 3D dimensions
C
C =======================================================================
C   INPUT:
C ========================
C
C   ----------------------
C   Fine Grid Dimensions:
C   ---------------------- 
C
C   ND       Dimension of problem
C   NX       Number of points in x-direction (excluding ghost points)
C   NY       Number of points in y-direction (excluding ghost points)
C   NZ       Number of points in z-direction (excluding ghost points)
C
C  ------------------------
C   Discretization:
C  ------------------------
C
C   IFD      Discrete operator index:
C            iFD .EQ. 1 => a 7  point discretization
C            iFD .NE. 1 => a 27 point discretization
C
C   ------------------------
C   Smoothing:
C   ------------------------   
C
C   IRELAX   Relaxation index (refer to BMG3D).
C
C   --------------------------
C   Work
C   --------------------------
C
C   IWORK     Integer work array (4*NOG)
C
C   -------------------------
C   Multigrid:
C   -------------------------
C
C   NOGm      Maximum number of grids.
C   NSOm      Maximum allowed space for 3D stencils (all grids).
C   NFm       Maximum allowed space for grid functions (all grids).
C
C   * note: NSOm and NFm may be zero if these arrays are to be stored
C           in the real work array RWORK.
C
C =======================================================================
C   OUTPUT:
C ===========================
C
C   NOG       Number of grids needed for the given (Nx,Ny,Nz)
C   NF        Storage for arrays Q and QF on all grids
C   NC        Storage for Q and QF on coarser grids
C             - Q the solution vector on all 3D grids
C             - QF the source vector on all 3D grids
C
C   NSO       Storage for the array SO
C             - SO holds the stencil on all 3D grids.
C   NSOR      Storage for the array SOR: 
C             - SOR holds the current residual, the reciprocal of 
C               the central stencil weight, and the LU decomposition
C               if iRELAX > 1.
C   NCI       Storage for the array CI
C             - CI holds the 3D interpolation weights on all 3D grids
C
C   NCBW      First dimension of ABD, which is set to the bandwidth
C             of the coarsest grid stencil (including the diagonal).
C   NCU       Second dimension of ABD, the dimension of the coarse grid.
C
C =======================================================================
C   LOCAL:
C ========================
C
C   Nxc       Number of points in the x-direction on a coarser grid
C   Nyc       Number of points in the y-direction on a coarser grid
C   Nzc       Number of points in the z-direction on a coarser grid
C
C =======================================================================

      IMPLICIT   NONE

      INCLUDE    'BMG_workspace.h'
      INCLUDE    'BMG_constants.h'
      INCLUDE    'BMG_parameters.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER  NFm, NBMG_iWORKm, NBMG_iWORK_PLm, NOGm,
     &         NBMG_rWORKm, NBMG_rWORK_PLm, NSOm

      INTEGER  BMG_iPARMS(NBMG_iPARMS), BMG_pWORK(NBMG_pWORK), 
     &         Nd, Nx, Ny, Nz, pSI, pSR

      LOGICAL  BMG_InWORK(NBMG_InWORK)

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C ========================================================================


      IF( Nd.EQ.2 ) THEN

         CALL  BMG2_SymStd_SETUP_PtrWork( 
     &                Nx, Ny, BMG_iPARMS, 
     &                NOGm, NFm, NSOm, NBMG_IWORKm, NBMG_rWORKm, 
     &                BMG_pWORK, BMG_InWORK, pSR, pSI,
     &                BMG_IOFLAG
     &                )

      ELSEIF( Nd.EQ.3 ) THEN

         CALL  BMG3_SymStd_SETUP_PtrWork( 
     &                Nx, Ny, Nz, BMG_iPARMS, 
     &                NOGm, NFm, NSOm, 
     &                NBMG_iWORKm, NBMG_rWORKm,
     &                NBMG_iWORK_PLm, NBMG_rWORK_PLm,
     &                BMG_pWORK, BMG_InWORK, pSR, pSI,
     &                BMG_IOFLAG
     &                )

      ENDIF

      RETURN
      END


