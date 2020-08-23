      SUBROUTINE BMG_SymStd_SOLVE_boxmg(
     &                Nd, Nx, Ny, Nz,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                SO, Q, QF, BMG_rWORK, BMG_iWORK, BMG_pWORK,
     &                BMG_iWORK_PL, NBMG_iWORK_PL, 
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )

C =======================================================================
C 
C   BMG_SymStd_SOLVE_boxmg
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG_SymStd_SOLVE_boxmg calls either BMG2_SymStd_SOLVE_boxmg or it
C   calls BMG3_SymStd_SOLVE_boxmg depending on ND.
C
C   Written:    2004/12/20 (TMA)
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
      INTEGER  Nd, Nx, Ny, Nz
      INTEGER  NBMG_iWORK_PL, NBMG_rWORK_PL 

      INTEGER  BMG_iPARMS(NBMG_iPARMS), BMG_pWORK(NBMG_pWORK), 
     &         BMG_iWORK(*)

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      REAL*8   BMG_rPARMS(NBMG_rPARMS)

      REAL*8   SO(*), Q(*), QF(*), BMG_rWORK(*)

      INTEGER  BMG_iWORK_PL(NBMG_iWORK_PL)
      REAL*8   BMG_rWORK_PL(NBMG_rWORK_PL)

C ---------------------------
C     Local Variables:
C
      INTEGER  NOG, NF, NC, NCI, NSO, NSOR, NCBW, NCU, NOGc

C ========================================================================



      IF( Nd.EQ.2 ) THEN

         NOG  = BMG_iPARMS(id_BMG2_DIM_NOG)
         NF   = BMG_iPARMS(id_BMG2_DIM_NF)
         NC   = BMG_iPARMS(id_BMG2_DIM_NC)
         NCI  = BMG_iPARMS(id_BMG2_DIM_NCI)
         NSO  = BMG_iPARMS(id_BMG2_DIM_NSO)
         NSOR = BMG_iPARMS(id_BMG2_DIM_NSOR)
         NCBW = BMG_iPARMS(id_BMG2_DIM_NCBW)
         NCU  = BMG_iPARMS(id_BMG2_DIM_NCU)

         CALL BMG2_SymStd_SOLVE_boxmg(
     &             Nx, Ny,
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &             Q, QF, BMG_rWORK(BMG_pWORK(ip_RES)), NF, NC, 
     &             SO, NSO, BMG_rWORK(BMG_pWORK(ip_SOR)), NSOR, 
     &             BMG_rWORK(BMG_pWORK(ip_CI)), NCI,
     &             BMG_rWORK(BMG_pWORK(ip_CSO)), 
     &             BMG_rWORK(BMG_pWORK(ip_CU)),
     &             NCBW, NCU, BMG_iWORK(BMG_pWORK(ip_iG)), 
     &             NOG, NOGc
     &             )         

      ELSEIF( Nd.EQ.3 ) THEN

         NOG  = BMG_iPARMS(id_BMG3_DIM_NOG)
         NF   = BMG_iPARMS(id_BMG3_DIM_NF)
         NC   = BMG_iPARMS(id_BMG3_DIM_NC)
         NCI  = BMG_iPARMS(id_BMG3_DIM_NCI)
         NSO  = BMG_iPARMS(id_BMG3_DIM_NSO)
         NSOR = BMG_iPARMS(id_BMG3_DIM_NSOR)
         NCBW = BMG_iPARMS(id_BMG3_DIM_NCBW)
         NCU  = BMG_iPARMS(id_BMG3_DIM_NCU)

         CALL BMG3_SymStd_SOLVE_boxmg(
     &             Nx, Ny, Nz, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &             Q, QF, BMG_rWORK(BMG_pWORK(ip_RES)), NF, NC, SO, NSO,
     &             BMG_rWORK(BMG_pWORK(ip_SOR)), NSOR, 
     &             BMG_rWORK(BMG_pWORK(ip_CI)), NCI, 
     &             BMG_rWORK(BMG_pWORK(ip_CSO)),
     &             BMG_rWORK(BMG_pWORK(ip_CU)), NCBW, NCU,
     &             BMG_iWORK(BMG_pWORK(ip_iG)), NOG, NOGc,
     &             BMG_iWORK_PL, NBMG_iWORK_PL, 
     &             BMG_rWORK_PL, NBMG_rWORK_PL
     &             )

      ENDIF


      RETURN
      END
