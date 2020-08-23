      SUBROUTINE BMG3_SymStd_SETUP_ITLI27_bl(
     &                SO, SOC, CI, IIF, JJF, KKF, IIC, JJC, KKC
     &                )

C ========================================================================

      IMPLICIT NONE

C ----------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'

C ----------------------------
C     Argument Declarations
C
      INTEGER IIC, IIF, JJC, JJF, KKC, KKF
      REAL*8  CI(IIC,JJC,KKC,26), SO(IIF,JJF,KKF,14),
     &        SOC(IIC,JJC,KKC,14)

C ----------------------------
C     Externals
C
      REAL*8   BMG3_SymStd_UTILS_reswt
      EXTERNAL BMG3_SymStd_UTILS_reswt

C ----------------------------
C     Local Declarations
C
      INTEGER I, IC, IIC1, IIF1, IL, ILG, ILGP, J, JC, JJC1, JJF1, JL,
     &        JLG, JLGP, K, KC, KK, KKC1, KKF1, KL, KLG, KLGP
      REAL*8  Q(7,7,7), R(7,7,7)

C ========================================================================

      IIC1 = IIC-1
      JJC1 = JJC-1
      KKC1 = KKC-1
      IIF1 = IIF-1
      JJF1 = JJF-1
      KKF1 = KKF-1

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,IC,IL,ILG,ILGP,J,JC,JL,JLG,JLGP,
C$OMP&         K,KC,KK,KL,KLG,KLGP,Q,R)
C$OMP& SHARED(CI,IIC,IIC1,IIF,JJC,JJC1,JJF,KKC,KKC1,KKF,SO,SOC)
      DO KK = 1, (IIC1-1)*(JJC1-1)*(KKC1-1)
         IC = MOD((KK-1),IIC1-1)+2
         I  = 2*(IC-1)
         JC = MOD((KK-1)/(IIC1-1),JJC1-1)+2
         J  = 2*(JC-1)
         KC = (KK-1)/((IIC1-1)*(JJC1-1))+2
         K  = 2*(KC-1)

         DO KL = 1,7
            DO JL = 1,7
               DO IL = 1,7
                  Q(IL,JL,KL) = rZERO
                  R(IL,JL,KL) = rZERO
               ENDDO
            ENDDO
         ENDDO

         !     XY-PLANE
         Q(4,4,4) = rONE
         Q(3,4,4) = CI(IC,JC,KC,LXYR)
         Q(5,4,4) = CI(IC+1,JC,KC,LXYL)
         Q(4,3,4) = CI(IC,JC,KC,LXYA)
         Q(4,5,4) = CI(IC,JC+1,KC,LXYB)
         Q(3,3,4) = CI(IC,JC,KC,LXYNE)
         Q(3,5,4) = CI(IC,JC+1,KC,LXYSE)
         Q(5,3,4) = CI(IC+1,JC,KC,LXYNW)
         Q(5,5,4) = CI(IC+1,JC+1,KC,LXYSW)

         !    XZ-PLANE
         Q(4,4,3) = CI(IC,JC,KC,LXZA)
         Q(4,4,5) = CI(IC,JC,KC+1,LXZB)
         Q(3,4,3) = CI(IC,JC,KC,LXZNE)
         Q(3,4,5) = CI(IC,JC,KC+1,LXZSE)
         Q(5,4,3) = CI(IC+1,JC,KC,LXZNW)
         Q(5,4,5) = CI(IC+1,JC,KC+1,LXZSW)

         !     YZ-PLANE
         Q(4,3,3) = CI(IC,JC,KC,LYZNW)
         Q(4,3,5) = CI(IC,JC,KC+1,LYZSW)
         Q(4,5,3) = CI(IC,JC+1,KC,LYZNE)
         Q(4,5,5) = CI(IC,JC+1,KC+1,LYZSE)
         
         !     CORNERS
         Q(3,3,3) = CI(IC,JC,KC,LTNE)
         Q(3,5,3) = CI(IC,JC+1,KC,LTSE)
         Q(5,3,3) = CI(IC+1,JC,KC,LTNW)
         Q(5,5,3) = CI(IC+1,JC+1,KC,LTSW)
         Q(3,3,5) = CI(IC,JC,KC+1,LBNE)
         Q(3,5,5) = CI(IC,JC+1,KC+1,LBSE)
         Q(5,3,5) = CI(IC+1,JC,KC+1,LBNW)
         Q(5,5,5) = CI(IC+1,JC+1,KC+1,LBSW)
         
         !   APPLY FINE-GRID OPERATOR

         DO KL = 2,6
            DO JL = 2,6
               DO IL = 2,6
                  
                  ILG = I + IL - 4
                  JLG = J + JL - 4
                  KLG = K + KL - 4
                  
                  ILGP = MIN(ILG+1,IIF)
                  JLGP = MIN(JLG+1,JJF)
                  KLGP = MIN(KLG+1,KKF)
                  
                  ILG = MIN(ILG,IIF)
                  JLG = MIN(JLG,JJF)
                  KLG = MIN(KLG,KKF)
                  
                  ILG = MAX(ILG,1)
                  JLG = MAX(JLG,1)
                  KLG = MAX(KLG,1)
                  
                  R(IL,JL,KL) 
     &                 =   SO(ILG,JLG,KLG,KP)*Q(IL,JL,KL)
     &                 - SO(ILG,JLG,KLG,KPW)*Q(IL-1,JL,KL)
     &                 - SO(ILG,JLGP,KLG,KPS)*Q(IL,JL+1,KL)
     &                 - SO(ILGP,JLG,KLG,KPW)*Q(IL+1,JL,KL)
     &                 - SO(ILG,JLG,KLG,KPS)*Q(IL,JL-1,KL)
     &                 - SO(ILG,JLG,KLG,KB)*Q(IL,JL,KL-1)
     &                 - SO(ILG,JLG,KLGP,KB)*Q(IL,JL,KL+1)
     &                 - SO(ILG,JLGP,KLG,KPNW)*Q(IL-1,JL+1,KL)
     &                 - SO(ILGP,JLGP,KLG,KPSW)*Q(IL+1,JL+1,KL)
     &                 - SO(ILGP,JLG,KLG,KPNW)*Q(IL+1,JL-1,KL)
     &                 - SO(ILG,JLG,KLG,KPSW)*Q(IL-1,JL-1,KL)
     &                 - SO(ILG,JLG,KLG,KBW)*Q(IL-1,JL,KL-1)
     &                 - SO(ILG,JLGP,KLG,KBNW)*Q(IL-1,JL+1,KL-1)
     &                 - SO(ILG,JLGP,KLG,KBN)*Q(IL,JL+1,KL-1)
     &                 - SO(ILGP,JLGP,KLG,KBNE)*Q(IL+1,JL+1,KL-1)
     &                 - SO(ILGP,JLG,KLG,KBE)*Q(IL+1,JL,KL-1)
     &                 - SO(ILGP,JLG,KLG,KBSE)*Q(IL+1,JL-1,KL-1)
     &                 - SO(ILG,JLG,KLG,KBS)*Q(IL,JL-1,KL-1)
     &                 - SO(ILG,JLG,KLG,KBSW)*Q(IL-1,JL-1,KL-1)
     &                 - SO(ILG,JLG,KLGP,KBE)*Q(IL-1,JL,KL+1)
     &                 - SO(ILG,JLGP,KLGP,KBSE)*Q(IL-1,JL+1,KL+1)
     &                 - SO(ILG,JLGP,KLGP,KBS)*Q(IL,JL+1,KL+1)
     &                 - SO(ILGP,JLGP,KLGP,KBSW)*Q(IL+1,JL+1,KL+1)
     &                 - SO(ILGP,JLG,KLGP,KBW)*Q(IL+1,JL,KL+1)
     &                 - SO(ILGP,JLG,KLGP,KBNW)*Q(IL+1,JL-1,KL+1)
     &                 - SO(ILG,JLG,KLGP,KBN)*Q(IL,JL-1,KL+1)
     &                 - SO(ILG,JLG,KLGP,KBNE)*Q(IL-1,JL-1,KL+1)
                  
               ENDDO
            ENDDO
         ENDDO  

         !     APPLY RESIDUAL WEIGHTING

         SOC(IC,JC,KC,KP) 
     &        =   BMG3_SymStd_UTILS_reswt(
     &             CI, R, IIC, JJC, KKC, IC, JC, KC, 4, 4, 4
     &            )

         SOC(IC,JC,KC,KPW) 
     &        = - BMG3_SymStd_UTILS_reswt(
     &             CI, R, IIC, JJC, KKC, IC-1, JC, KC, 2, 4, 4
     &            )

         SOC(IC,JC,KC,KPS) 
     &        = - BMG3_SymStd_UTILS_reswt(
     &             CI, R, IIC, JJC, KKC, IC, JC-1, KC, 4, 2, 4
     &            )

         SOC(IC,JC,KC,KPSW) 
     &        = - BMG3_SymStd_UTILS_reswt(
     &             CI, R, IIC, JJC, KKC, IC-1, JC-1, KC, 2, 2, 4
     &            )

         SOC(IC,JC+1,KC,KPNW)
     &        = - BMG3_SymStd_UTILS_reswt(
     &             CI, R, IIC, JJC, KKC, IC-1, JC+1, KC, 2, 6, 4 
     &            )

         SOC(IC,JC,KC,KB) 
     &        = - BMG3_SymStd_UTILS_reswt(
     &             CI, R, IIC, JJC, KKC, IC, JC, KC-1, 4, 4, 2 
     &            )

         SOC(IC,JC,KC,KBSW)
     &        = - BMG3_SymStd_UTILS_reswt(
     &             CI, R, IIC, JJC, KKC, IC-1, JC-1, KC-1, 2, 2, 2 
     &            )

         SOC(IC,JC,KC,KBW) 
     &        = - BMG3_SymStd_UTILS_reswt(
     &             CI, R, IIC, JJC, KKC, IC-1, JC, KC-1, 2, 4, 2 
     &            )
         
         SOC(IC,JC+1,KC,KBNW)
     &        = - BMG3_SymStd_UTILS_reswt(
     &             CI, R, IIC, JJC, KKC, IC-1, JC+1, KC-1, 2, 6, 2
     &            )
         SOC(IC,JC+1,KC,KBN)
     &        = - BMG3_SymStd_UTILS_reswt(
     &             CI, R, IIC, JJC, KKC, IC, JC+1, KC-1, 4, 6, 2
     &            )
         SOC(IC+1,JC+1,KC,KBNE)
     &        = - BMG3_SymStd_UTILS_reswt(
     &             CI, R, IIC, JJC, KKC, IC+1, JC+1, KC-1, 6, 6, 2 
     &            )
         SOC(IC+1,JC,KC,KBE)
     &        = - BMG3_SymStd_UTILS_reswt(
     &             CI, R, IIC, JJC, KKC, IC+1, JC, KC-1, 6, 4, 2 
     &            )
         SOC(IC+1,JC,KC,KBSE) 
     &        = - BMG3_SymStd_UTILS_reswt(
     &             CI, R, IIC, JJC, KKC, IC+1, JC-1, KC-1, 6, 2, 2 
     &            )
         SOC(IC,JC,KC,KBS) 
     &        = - BMG3_SymStd_UTILS_reswt(
     &             CI, R, IIC, JJC, KKC, IC, JC-1, KC-1, 4, 2, 2
     &            )
              
      ENDDO
C$OMP END PARALLEL DO

      RETURN         
      END
