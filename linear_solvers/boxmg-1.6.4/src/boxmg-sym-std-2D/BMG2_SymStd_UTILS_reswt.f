      REAL*8 FUNCTION BMG2_SymStd_UTILS_reswt(
     &                     CI, R, IIC, JJC, ICL, JCL, IL, JL 
     &                     )

C ==========================================================================
C
C   BMG2_SymStd_UTILS_reswt.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_UTILS_reswt.f local residual weighting routing used
C   by BMG2_SymStd_cg_ITLI_bl.f.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2004/01/17 (JED)
C
C ==================================================================
C
C
C ==========================================================================

      IMPLICIT NONE

C ----------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'

C ----------------------------
C     Argument Declarations
C
      INTEGER ICL, IIC, IL, JCL, JJC, JL
      REAL*8  CI(IIC,JJC,8), R(7,7)

C ----------------------------
C     Local Variables
C     
      INTEGER ICP, JCP

C ==========================================================================

      ICP = MIN(ICL+1,IIC)
      JCP = MIN(JCL+1,JJC)
  
      BMG2_SymStd_UTILS_reswt = CI(ICL,JCL,LNE)*R(IL-1,JL-1)
     &                        + CI(ICL,JCL,LA)*R(IL,JL-1)
     &                        + CI(ICP,JCL,LNW)*R(IL+1,JL-1)
     &                        + CI(ICL,JCL,LR)*R(IL-1,JL)
     &                        + R(IL,JL)
     &                        + CI(ICP,JCL,LL)*R(IL+1,JL)
     &                        + CI(ICL,JCP,LSE)*R(IL-1,JL+1)
     &                        + CI(ICL,JCP,LB)*R(IL,JL+1)
     &                        + CI(ICP,JCP,LSW)*R(IL+1,JL+1)

C ==========================================================================

      RETURN
      END
