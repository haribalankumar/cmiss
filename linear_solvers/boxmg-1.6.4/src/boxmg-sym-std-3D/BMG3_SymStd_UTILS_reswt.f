      REAL*8 FUNCTION BMG3_SymStd_UTILS_reswt(
     &                CI, R, IIC, JJC, KKC, ICL, JCL, KCL, IL, JL, KL 
     &                )

      IMPLICIT NONE

C ----------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'

C ----------------------------
C     Argument Declarations
C
      
      INTEGER ICL, IIC, IL, JCL, JJC, JL, KCL, KKC, KL
      REAL*8  CI(IIC,JJC,KKC,26), R(7,7,7)

C ----------------------------
C     Local Variables
C     
      INTEGER ICP, JCP, KCP

      ICP = MIN0(ICL+1,IIC)
      JCP = MIN0(JCL+1,JJC)
      KCP = MIN0(KCL+1,KKC)
  
      BMG3_SymStd_UTILS_reswt = 
     &       CI(ICL,JCL,KCL,LXYNE)*R(IL-1,JL-1,KL)
     &     + CI(ICL,JCL,KCL,LXYA)*R(IL,JL-1,KL)
     &     + CI(ICP,JCL,KCL,LXYNW)*R(IL+1,JL-1,KL)
     &     + CI(ICL,JCL,KCL,LXYR)*R(IL-1,JL,KL)
     &     + R(IL,JL,KL)
     &     + CI(ICP,JCL,KCL,LXYL)*R(IL+1,JL,KL)
     &     + CI(ICL,JCP,KCL,LXYSE)*R(IL-1,JL+1,KL)
     &     + CI(ICL,JCP,KCL,LXYB)*R(IL,JL+1,KL)
     &     + CI(ICP,JCP,KCL,LXYSW)*R(IL+1,JL+1,KL)
     &     + CI(ICL,JCL,KCL,LTNE)*R(IL-1,JL-1,KL-1)
     &     + CI(ICL,JCL,KCL,LYZNW)*R(IL,JL-1,KL-1)
     &     + CI(ICP,JCL,KCL,LTNW)*R(IL+1,JL-1,KL-1)
     &     + CI(ICL,JCL,KCL,LXZNE)*R(IL-1,JL,KL-1)
     &     + CI(ICL,JCL,KCL,LXZA)*R(IL,JL,KL-1)
     &     + CI(ICP,JCL,KCL,LXZNW)*R(IL+1,JL,KL-1)
     &     + CI(ICL,JCP,KCL,LTSE)*R(IL-1,JL+1,KL-1)
     &     + CI(ICL,JCP,KCL,LYZNE)*R(IL,JL+1,KL-1)
     &     + CI(ICP,JCP,KCL,LTSW)*R(IL+1,JL+1,KL-1)
     &     + CI(ICL,JCL,KCP,LBNE)*R(IL-1,JL-1,KL+1)
     &     + CI(ICL,JCL,KCP,LYZSW)*R(IL,JL-1,KL+1)
     &     + CI(ICP,JCL,KCP,LBNW)*R(IL+1,JL-1,KL+1)
     &     + CI(ICL,JCL,KCP,LXZSE)*R(IL-1,JL,KL+1)
     &     + CI(ICL,JCL,KCP,LXZB)*R(IL,JL,KL+1)
     &     + CI(ICP,JCL,KCP,LXZSW)*R(IL+1,JL,KL+1)
     &     + CI(ICL,JCP,KCP,LBSE)*R(IL-1,JL+1,KL+1)
     &     + CI(ICL,JCP,KCP,LYZSE)*R(IL,JL+1,KL+1)
     &     + CI(ICP,JCP,KCP,LBSW)*R(IL+1,JL+1,KL+1)

      RETURN
      END
