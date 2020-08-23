      SUBROUTINE BMG_INIT_RES(Flag, Nd, Nx, Ny, Nz,
     &              BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &              BMG_rWORK, BMG_iWORK, BMG_pWORK,
     &              BMG_iWORK_PL, NBMG_iWORK_PL, 
     &              BMG_rWORK_PL, NBMG_rWORK_PL, RES_L2)

C =======================================================================
C 
C   BMG_INIT_RES
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG_SOLVE calls either 2D BOXMG or 3D BOXMG-PCG solve or it calls 
c   3D BOXMG or 3D BOXMG-PCG solve depending on ND and SOLVER. 
C
C   Written:    2005/01/19 (TMA)
C
C =======================================================================
C   INPUT:
C ========================
C
C   SOLVER   Specifies either BOXMG or BOXMG-PCG solver.
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
C   -----------------------
C   BOXMG Parameter Arrays:
C   -----------------------
C
C   BMG_iPARMS     See BMG2_SymStd_SOLVE_boxmg or BMG3_SymStd_SOLVE_boxmg.
C   BMG_rPARMS                  "              "              "
C   BMG_IOFLAG                  "              "              "
C   BMG_rWORK      Arrays containing solution, rhs, matrix, residual, etc.
C   BMG_iWORK      Array containing any workspace needed for integer work.
C   BMG_pWORK      Array containing pointers to locations in rWORK and iWORK.
C 
c   BMG_iWORK_PL   Array containing INT workspace for 3D plane relaxation.
C   BMG_rWORK_PL   Array containing REAL "        "        "         "
C
C   NBMG_iWORK_PL  Size of BMG_iWORK_PL.
C   NBMG_rWORK_PL  Size of BMG_rWORK_PL.
C
C =======================================================================
C   LOCAL:
C ========================
C
C   NOG    Number of grids.
C   NF     Number of fine grid points.
C   NC     Total number of all coarse grid points.
C   NCI    Size of array for interpolation in BMG_rWORK.
C   NSO    Size of array for operator on all grids in BMG_rWORK.
C   NSOR   Size of array SOR in BMG_rWORK.
C   NCBW   Sizes related to coarsest grid equations.
C   NCU    Sizes related to coarsest grid equations
C   NOGc   Number of coarsest grid.
C
C =======================================================================

      IMPLICIT   NONE

      INCLUDE    'BMG_workspace.h'
      INCLUDE    'BMG_constants.h'
      INCLUDE    'BMG_parameters.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER  FLag, Nd, Nx, Ny, Nz
      INTEGER  NBMG_iWORK_PL, NBMG_rWORK_PL 

      INTEGER  BMG_iPARMS(NBMG_iPARMS), 
     &         BMG_pWORK(NBMG_pWORK), 
     &         BMG_iWORK(*)

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      REAL*8   BMG_rPARMS(NBMG_rPARMS)

      REAL*8   BMG_rWORK(*), RES_L2

      INTEGER  BMG_iWORK_PL(NBMG_iWORK_PL)
      REAL*8   BMG_rWORK_PL(NBMG_rWORK_PL)

C ---------------------------
C     Local Variables:
C
      INTEGER  NOG, NF, NC, NCI, NSO, NSOR, NCBW, NCU, NOGc
      INTEGER  IFD, NStncl, I

      INTEGER p_U, p_SO, p_SOR, p_CI, II, JJ, KK


C ========================================================================

      IF( Nd.EQ.2 ) THEN

         NOG  = BMG_iPARMS(id_BMG2_DIM_NOG)
         IFD  = BMG_iPARMS(id_BMG2_STENCIL)

         IF ( IFD.NE.1 ) THEN
            NStncl=5
         ELSE
            NStncl=3
         ENDIF

         CALL BMG2_SymStd_GET_pointers( 
     &        NOG, BMG_iWORK(BMG_pWORK(ip_iG)), NOG, 
     &        II, JJ, p_U, p_SO, p_SOR, p_CI
     &        )

         CALL BMG2_SymStd_residual( 
     &        NOG, 
     &        BMG_rWORK(BMG_pWORK(ip_SO)),
     &        BMG_rWORK(BMG_pWORK(ip_Q)), 
     &        BMG_rWORK(BMG_pWORK(ip_U)), 
     &        BMG_rWORK(BMG_pWORK(ip_RES)),
     &        II, JJ, NOG, IFD, NStncl
     &        )

         CALL BMG2_SymStd_UTILS_norm_l2(
     &        BMG_rWORK(BMG_pWORK(ip_RES)), 
     &        II, JJ, KK, RES_L2  
     &        )

         IF( Flag.EQ.0 ) THEN
            WRITE (*,405) '*** THE INITIAL RESIDUAL (l2-NORM) = ',RES_L2
         ELSEIF( Flag.EQ.1 ) THEN
            WRITE (*,405) '*** THE FINAL RESIDUAL (l2-NORM) = ', RES_L2
         ENDIF

      ELSEIF( Nd.EQ.3 ) THEN

         NOG  = BMG_iPARMS(id_BMG3_DIM_NOG)
         IFD  = BMG_iPARMS(id_BMG3_STENCIL)
         
         IF ( IFD.NE.1 ) THEN
            NStncl=14
         ELSE
            NStncl=4
         ENDIF

         CALL BMG3_SymStd_GET_pointers( 
     &        NOG, BMG_iWORK(BMG_pWORK(ip_iG)), NOG, 
     &        p_U, p_SO, p_SOR, p_CI, II, JJ, KK 
     &        )

         CALL BMG3_SymStd_residual( 
     &        NOG, NOG, IFD, 
     &        BMG_rWORK(BMG_pWORK(ip_U)), 
     &        BMG_rWORK(BMG_pWORK(ip_Q)), 
     &        BMG_rWORK(BMG_pWORK(ip_SO)),
     &        BMG_rWORK(BMG_pWORK(ip_RES)),
     &        II, JJ, KK, NStncl
     &        )

         CALL BMG3_SymStd_UTILS_norm_l2(
     &        BMG_rWORK(BMG_pWORK(ip_RES)), 
     &        II, JJ, KK, RES_L2  
     &        )
                    
      ENDIF


C ==========================================================================

 405  FORMAT (/,1X,A,1X,1P,E16.9,/)

C ==========================

      RETURN
      END
