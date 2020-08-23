      SUBROUTINE BMG2_SymStd_PRECON_boxmg( 
     &                       Nx, Ny, 
     &                       BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                       Q, QF, NFm, SO, NSOm, NOGm,
     &                       BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                       BMG_rWORK, NBMG_rWORKm
     &                       )

C ==========================================================================
C 
C   BMG2_SymStd_PRECON_boxmg
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_PRECON_boxmg is used to wrap the call to 
C   the BoxMG solver, BMG2_SymStd_SOLVE_boxmg, for 
C   preconditioning in BMG2_SymStd_SOLVE_pcg.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2003/07/10 (JDM)
C
C =======================================================================
C   INPUT:
C ========================



C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  NBMG_iWORKm, NBMG_rWORKm, NOGm, NFm, NSOm, Nx, Ny

      INTEGER  BMG_iPARMS(NBMG_iPARMS), BMG_iWORK(NBMG_iWORKm),
     &         BMG_pWORK(NBMG_pWORK)
      REAL*8   BMG_rPARMS(NBMG_rPARMS), BMG_rWORK(NBMG_rWORKm),
     &         Q(NFm), QF(NFm), SO(NSOm)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C ----------------------------
C     Local Declarations
C
      INTEGER  NOG, NF, NC, NSO, NCI, NSOR, NCBW, NCU,
     &         p_CI, p_CSO, p_CU, p_iG, p_RES, p_SOR

C =========================================================================

      !
      !  BoxMG dimensional parameters
      !
      NOG  = BMG_iPARMS(id_BMG2_DIM_NOG)
      NF   = BMG_iPARMS(id_BMG2_DIM_NF)
      NC   = BMG_iPARMS(id_BMG2_DIM_NC)
      NSO  = BMG_iPARMS(id_BMG2_DIM_NSO)
      NCI  = BMG_iPARMS(id_BMG2_DIM_NCI)
      NSOR = BMG_iPARMS(id_BMG2_DIM_NSOR)
      NCBW = BMG_iPARMS(id_BMG2_DIM_NCBW)
      NCU  = BMG_iPARMS(id_BMG2_DIM_NCU)

      !
      !  BoxMG pointers
      !
      p_RES = BMG_pWORK(ip_RES)
      p_SOR = BMG_pWORK(ip_SOR)
      p_CI  = BMG_pWORK(ip_CI)
      p_CSO = BMG_pWORK(ip_CSO)
      p_CU  = BMG_pWORK(ip_CU)
      p_iG  = BMG_pWORK(ip_iG)
      
      CALL BMG2_SymStd_SOLVE_boxmg( 
     &                 Nx, Ny, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                 Q, QF, BMG_rWORK(p_RES),
     &                 NF, NC, SO, NSO, BMG_rWORK(p_SOR), NSOR,
     &                 BMG_rWORK(p_CI), NCI,
     &                 BMG_rWORK(p_CSO), BMG_rWORK(p_CU), NCBW, NCU,
     &                 BMG_iWORK(p_iG), NOGm, NOG
     &                 )
      IF (BMG_iPARMS(id_BMG2_Err_Code) .ne. iZERO) THEN
         RETURN
      ENDIF

C =========================================================================

      RETURN
      END

