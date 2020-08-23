      SUBROUTINE BMG3_SymStd_relax_CG( 
     &                KG, SO, QF, Q, RES, SOR, II, JJ, KK,
     &                RES_L2, BMG_IOFLAG, NOG, IFD, NStncl, 
     &                BMG_IJK_MAP
     &                )


C ==========================================================================
C
C***BEGIN PROLOGUE  BMG3_SymStd_relax_CG
C
C***PURPOSE  
C 
C     This subroutine is conjugate gradients used as a relaxation method
C     for Black Box Multigrid.
C
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ------------------------------------------------
C     Argument Declarations
C
      INTEGER II, JJ, KK, NStncl

      INTEGER IFD, KG, NOG, BMG_IJK_MAP(*)
      REAL*8  Q(II,JJ,KK), QF(II,JJ,KK), RES(II,JJ,KK), RES_L2,
     &        SO(II,JJ,KK,NStncl), SOR(4*II*JJ*KK)
      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)

C ----------------------------
C     Local Declarations
C
      INTEGER  J, p_D, p_R, p_Z, p_P
      REAL*8   malpha, alpha, delta0, delta1, beta

C =========================================================================

      p_D = 1
      p_R = p_D + II*JJ*KK
      p_P = p_R + II*JJ*KK
      p_Z = p_P + II*JJ*KK

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>>>>> BEGIN: ZERO ARRAYS <<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      CALL BMG3_SymStd_UTILS_rV_zero( 
     &          SOR(p_R), II, JJ, KK
     &          )

      CALL BMG3_SymStd_UTILS_rV_zero( 
     &          SOR(p_P), II, JJ, KK
     &          )

      CALL BMG3_SymStd_UTILS_rV_zero( 
     &          SOR(p_Z), II, JJ, KK
     &          )

      CALL BMG3_SymStd_residual(
     &          KG, NOG, IFD, 
     &          Q, QF, SO, SOR(p_R), II, JJ, KK, NStncl,
     &          BMG_IJK_MAP
     &          )

      CALL BMG3_SymStd_UTILS_dcopy(
     &          SOR(p_R), SOR(p_Z), II, JJ, KK
     &          )

      CALL BMG3_SymStd_UTILS_diag_scal(
     &          SOR(p_D), SOR(p_Z), II, JJ, KK
     &          )

      CALL BMG3_SymStd_UTILS_dcopy(
     &          SOR(p_Z), SOR(p_P), II, JJ, KK
     &          )

      CALL BMG3_SymStd_UTILS_dot_l2(
     &          SOR(p_R), SOR(p_Z), II, JJ, KK, delta0
     &          )

C      CALL BMG3_SymStd_UTILS_matvec(
C     &           KG, SOR(p_P), SOR(p_Z), SO, 
C     &           II, JJ, KK, NStncl, BMG_IJK_MAP 
C     &           )

C      CALL BMG3_SymStd_UTILS_dot_l2(
C     &           SOR(p_Z), SOR(p_P), II, JJ, KK, delta1
C     &           )

C      alpha = delta0/delta1

C ----------------------------------
C        Calculate q <- q + alpha*p
C ----------------------------------
      
C      CALL BMG3_SymStd_UTILS_daxpy( 
C     &           alpha, SOR(p_P), Q, II, JJ, KK
C     &           )
      
C ----------------------------------
C        Calculate r <- r - alpha*z
C ----------------------------------
         
C         malpha = -1.0*alpha

C         CALL BMG3_SymStd_UTILS_daxpy(
C     &           malpha, SOR(p_Z), SOR(p_R), II, JJ, KK
C     &           )


C ================================
C     Start main loop.
C ================================


      DO J=1,2

C ---------------------------------------
C
C        Calculate delta0/(Apj,pj) by
C           (2) Zj     = Apj
C           (3) delta1 = (Zj,Pj)
C           (4) alpha  = delta0/delta1
C
C ---------------------------------------

         CALL BMG3_SymStd_UTILS_matvec(
     &             KG, SOR(p_P), SOR(p_Z), SO, 
     &             II, JJ, KK, NStncl, BMG_IJK_MAP 
     &             )

         CALL BMG3_SymStd_UTILS_dot_l2(
     &             SOR(p_Z), SOR(p_P), II, JJ, KK, delta1
     &             )

         alpha = delta0/delta1

C ----------------------------------
C        Calculate q <- q + alpha*p
C ----------------------------------

         CALL BMG3_SymStd_UTILS_daxpy( 
     &             alpha, SOR(p_P), Q, II, JJ, KK
     &             )
         
C ----------------------------------
C        Calculate r <- r - alpha*z
C ----------------------------------
         
         malpha = -1.0*alpha

         CALL BMG3_SymStd_UTILS_daxpy(
     &             malpha, SOR(p_Z), SOR(p_R), II, JJ, KK
     &             )
         
         CALL BMG3_SymStd_UTILS_dcopy(
     &             SOR(p_R), SOR(p_Z), II, JJ, KK
     &             )

         CALL BMG3_SymStd_UTILS_diag_scal(
     &             SOR(p_D), SOR(p_Z), II, JJ, KK
     &             )

C ------------------------------------------------------
C        Calculate delta1 = <R,Z> = <R,inv(M)R>
C ------------------------------------------------------

         CALL BMG3_SymStd_UTILS_dot_l2(
     &             SOR(p_R), SOR(p_Z), II, JJ, KK, delta1
     &             )

C ----------------------------------------------
C        Calculate beta. 
C ---------------------------------------
     
         beta   = delta1/delta0
         delta0 = delta1

C ----------------------------------
C        Calculate p <- z + beta*p
C ----------------------------------

         CALL BMG3_SymStd_UTILS_dxpby(
     &             beta, SOR(p_R), SOR(p_P), II, JJ, KK
     &             )


      ENDDO

C =======================================================================

      IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) ) THEN
         CALL BMG3_SymStd_residual( 
     &             KG, NOG, IFD, Q, QF, SO, RES, II, JJ, KK, NStncl,
     &             BMG_IJK_MAP
     &             )
         CALL BMG3_SymStd_UTILS_norm_l2( 
     &             RES, II, JJ, KK, RES_L2
     &             )
      ENDIF

C ============================

      RETURN
      END
