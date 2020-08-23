      SUBROUTINE BMG3_SymStd_interp_add(
     &                KCG, KFG, Q ,QC, SO, RES, CI,
     &                IIC, JJC, KKC, IIF, JJF, KKF, NStncl
     &                )

C ==========================================================================
C
C   BMG3_SymStd_interp_add.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_interp_add.f interpolates Q from the coarse mesh, KCG, to 
C   the fine mesh, KFG, and adds the result to Q on fine mesh.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   2000/03/06   - written (JDM)
C                - the core computation was cut from mgia3.f 
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
      INCLUDE 'BMG_stencils.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER IIC, IIF, JJC, JJF, KCG, KFG, KKC, KKF, KL, NStncl
      REAL*8  CI(IIC,JJC,KKC,26), Q(IIF,JJF,KKF), QC(IIC,JJC,KKC), 
     &        RES(IIF,JJF,KKF), SO(IIF,JJF,KKF,NStncl)

C ----------------------------
C     Local Declarations
C
      INTEGER IC, I, IICF, IICF2, IICF3, IIC2, IIF2,
     &        JC, J, JJCF, JJCF2, JJCF3, JJC2, JJF2,
     &        KC, K, KKCF, KKCF2, KKCF3, KKC2, KKF2,
     &        MAXIJ, MAXIJK
      REAL*8  A, AQ

C ==========================================================================

C -------------------------------------------------
C     Useful index bounds:
C -------------------------------------------------

      IIC2=IIC-2
      JJC2=JJC-2
      KKC2=KKC-2

      IIF2=IIF-2
      JJF2=JJF-2
      KKF2=KKF-2

      IICF=(IIF2/2)+3
      JJCF=(JJF2/2)+3
      KKCF=(KKF2/2)+3

      IICF2=IICF-2
      JJCF2=JJCF-2
      KKCF2=KKCF-2

      IICF3=IICF-3
      JJCF3=JJCF-3
      KKCF3=KKCF-3

C --------------------------------------------------
C     NB: division is only possible in the interior
C --------------------------------------------------
      
      MAXIJK = IIF2*JJF2*KKF2
      MAXIJ  = IIF2*JJF2

C$OMP PARALLEL DO 
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,J,K,KL)
C$OMP& SHARED(IIF2,JJF2,MAXIJ,MAXIJK,RES,SO)
      DO KL = 0,MAXIJK-1
         I = MOD(KL,IIF2)+2
         J = MOD(KL/IIF2,JJF2)+2
         K = (KL/MAXIJ)+2
         RES(I,J,K)=RES(I,J,K)/SO(I,J,K,kp)
      ENDDO
C$OMP END PARALLEL DO

C --------------------------------------------------
C   interpolate answers from coarse to fine mesh 
C   and add to answers on fine mesh.
C --------------------------------------------------

      MAXIJK = IICF2*JJCF2*KKC2
      MAXIJ  = IICF2*JJCF2

C$OMP PARALLEL DO 
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(A,AQ,I,IC,J,JC,K,KC,KL)
C$OMP& SHARED(CI,IICF2,JJCF2,MAXIJ,MAXIJK,Q,QC,RES)
      DO KL = 0,MAXIJK-1

         IC = MOD(KL,IICF2)+2
         I = 2*(IC-1)
         JC = MOD(KL/IICF2,JJCF2)+2
         J = 2*(JC-1)
         KC = (KL/MAXIJ)+2
         K = 2*(KC-1)

         Q(I,J,K) = Q(I,J,K) + QC(IC,JC,KC)

         IF( IC.NE.2 ) THEN
            A = CI(IC,JC,KC,lxyr)*QC(IC,JC,KC)
     &        + CI(IC,JC,KC,lxyl)*QC(IC-1,JC,KC)
            Q(I-1,J,K) = Q(I-1,J,K) + A + RES(I-1,J,K)
         ENDIF

         IF( JC.NE.2 ) THEN
            AQ = CI(IC,JC,KC,lxya)*qc(IC,JC,KC)
     &         + CI(IC,JC,KC,lxyb)*qc(IC,JC-1,KC)
            Q(I,J-1,K) = Q(I,J-1,K) + AQ + RES(I,J-1,K)
         ENDIF

         IF( IC.NE.2 .AND. JC.NE.2 ) THEN
            A = CI(IC,JC,KC,lxysw)*QC(IC-1,JC-1,KC)
     &        + CI(IC,JC,KC,lxynw)*QC(IC-1,JC,KC)
     &        + CI(IC,JC,KC,lxyne)*QC(IC,JC,KC)
     &        + CI(IC,JC,KC,lxyse)*QC(IC,JC-1,KC)
            Q(I-1,J-1,K) = Q(I-1,J-1,K) + A + RES(I-1,J-1,K)
         ENDIF

      ENDDO
C$OMP END PARALLEL DO

      MAXIJK = IIC2*JJCF2*KKCF3
      MAXIJ  = IIC2*JJCF2

C$OMP PARALLEL DO 
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,IC,J,JC,K,KC,KL)
C$OMP& SHARED(CI,IIC2,JJCF2,MAXIJ,MAXIJK,Q,QC,RES)
      DO KL = 0, MAXIJK-1

         IC = MOD(KL,IIC2)+2
         I = 2*(IC-1)
         JC = MOD(KL/IIC2,JJCF2)+2
         J = 2*(JC-1)
         KC = (KL/MAXIJ)+3
         K = 2*(KC-1)-1

         IF( JC.EQ.2 ) THEN

            Q(i,j,k) = Q(i,j,k) + CI(IC,JC,KC,lxza)*QC(IC,JC,KC)
     &           + CI(IC,JC,KC,lxzb)*QC(IC,JC,KC-1) + RES(i,j,k)

         ELSE
            
            Q(i,j,k) = Q(i,j,k) 
     &               + CI(IC,JC,KC,lxza)*QC(IC,JC,KC)
     &               + CI(IC,JC,KC,lxzb)*QC(IC,JC,KC-1)
     &               + RES(i,j,k)

            Q(i,j-1,k) = Q(i,j-1,k)
     &                 + CI(IC,JC,KC,lyznw)*QC(IC,JC,KC)
     &                 + CI(IC,JC,KC,lyzne)*QC(IC,JC-1,KC)
     &                 + CI(IC,JC,KC,lyzsw)*QC(IC,JC,KC-1)
     &                 + CI(IC,JC,KC,lyzse)*QC(IC,JC-1,KC-1)
     &                 + RES(i,j-1,k)
         
         ENDIF

      ENDDO
C$OMP END PARALLEL DO
      
      MAXIJK = IICF3*JJCF2*KKCF3
      MAXIJ  = IICF3*JJCF2

C$OMP PARALLEL DO 
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,IC,J,JC,K,KC,KL)
C$OMP& SHARED(CI,IICF3,JJCF2,MAXIJ,MAXIJK,Q,QC,RES)
      DO KL = 0, MAXIJK-1
         
         IC = MOD(KL,IICF3)+3
         I = 2*(IC-1)-1
         JC = MOD(KL/IICF3,JJCF2)+2
         J = 2*(JC-1)
         KC = (KL/MAXIJ)+3
         K = 2*(KC-1)-1
         
         IF(JC.EQ.2) THEN

            Q(i,j,k) = Q(i,j,k)
     &           + CI(IC,JC,KC,lxznw)*QC(IC-1,JC,KC)
     &           + CI(IC,JC,KC,lxzne)*QC(IC,JC,KC)
     &           + CI(IC,JC,KC,lxzsw)*QC(IC-1,JC,KC-1)
     &           + CI(IC,JC,KC,lxzse)*QC(IC,JC,KC-1)
     &           + RES(i,j,k)

         ELSE

            Q(i,j,k) = Q(i,j,k) 
     &           + CI(IC,JC,KC,lxznw)*QC(IC-1,JC,KC)
     &           + CI(IC,JC,KC,lxzne)*QC(IC,JC,KC)
     &           + CI(IC,JC,KC,lxzsw)*QC(IC-1,JC,KC-1)
     &           + CI(IC,JC,KC,lxzse)*QC(IC,JC,KC-1)
     &           + RES(i,j,k)

            Q(i,j-1,k) = Q(i,j-1,k) 
     &           + CI(IC,JC,KC,ltnw)*QC(IC-1,JC,KC)
     &           + CI(IC,JC,KC,ltne)*QC(IC,JC,KC)
     &           + CI(IC,JC,KC,ltsw)*QC(IC-1,JC-1,KC)
     &           + CI(IC,JC,KC,ltse)*QC(IC,JC-1,KC)
     &           + CI(IC,JC,KC,lbnw)*QC(IC-1,JC,KC-1)
     &           + CI(IC,JC,KC,lbne)*QC(IC,JC,KC-1)
     &           + CI(IC,JC,KC,lbsw)*QC(IC-1,JC-1,KC-1)
     &           + CI(IC,JC,KC,lbse)*QC(IC,JC-1,KC-1)
     &           + RES(i,j-1,k)

         ENDIF

      ENDDO
COMP END DO PARALLEL

C ==========================================================================

      return
      end
