      SUBROUTINE BMG2_SymStd_SETUP_ITLI_bl( 
     &                KF, KC, SO, SOC, CI, 
     &                IIF, JJF, IIC, JJC, NOG, IFD, NStncl,
     &                BMG_IOFLAG, BMG_iPARMS
     &                )

C ==========================================================================
C
C   BMG2_SymStd_SETUP_cg_ITLI.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_SETUP_cg_ITLI.f constructs the Galerkin (variational)
C   coarse-grid operator on the coasre grid, KC, given the fine-grid
C   stencil, SO, and interpolation stencil, CI, on the fine-grid, KF.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2004/01/17 (JED)
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
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'
      
C ---------------------------
C    Argument Declarations:
C
      INTEGER   IIC, IIF, IFD, JJC, JJF, KC, KF, NOG, NStncl
      REAL*8    CI(IIC,JJC,8), SO(IIF,JJF,NStncl), SOC(IIC,JJC,5) 

      LOGICAL   BMG_IOFLAG(NBMG_IOFLAG)
      INTEGER   BMG_iPARMS(NBMG_iPARMS)

C --------------------------
C     Externals
C      
      REAL*8   BMG2_SymStd_UTILS_reswt
      EXTERNAL BMG2_SymStd_UTILS_reswt

C --------------------------
C     Local Declarations:
C
      INTEGER   IC, I, IIC1, IIF1, IL, ILG, ILGP, 
     &          JC, J, JJC1, JJF1, JL, JLG, JLGP
     	INTEGER   KK
      REAL*8    Q(7,7), R(7,7) 


C ==========================================================================

C ----------------------------------
C     Sanity Check:
C ----------------------------------
      
      IF (KF-1.NE.KC ) THEN
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,*) 'ERROR: BMG2_SymStd_SETUP_cg_ITLI   .... '
            WRITE(*,*) '*****  KC = ', KC
            WRITE(*,*) '*****  KF = ', KF
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,16)
         RETURN

      ENDIF

C -----------------------------------
C     Useful indexing bounds:
C -----------------------------------

      IIC1=IIC-1
      JJC1=JJC-1

      IIF1=IIF-1
      JJF1=JJF-1


      IF ( IFD.NE.1 .OR. KF.LT.NOG ) THEN

C******************************
C     begin computation of grid kc difference operator when kf difference
C     operator is nine point unless kc. ge. NOG
C     

C$OMP PARALLEL DO
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,IC,IL,ILG,ILGP,J,JC,JL,JLG,JLGP,KK,Q,R)
C$OMP& SHARED(CI,IIC,IIF,JJC,JJF,SO,SOC)
	   DO KK=1,(IIC-2)*(JJC-2)

		IC=MOD(KK-1,IIC-2)+2
		I=2*(IC-1)

	      JC=(KK-1)/(IIC-2)+2
		J=2*(JC-1)

            DO JL = 1,7
               DO IL = 1,7
                  Q(IL,JL) = rZERO
                  R(IL,JL) = rZERO
               ENDDO
            ENDDO   

            Q(4,4) = rONE
            Q(3,4) = CI(IC,JC,LR)
            Q(5,4) = CI(IC+1,JC,LL)
            Q(4,3) = CI(IC,JC,LA)
            Q(4,5) = CI(IC,JC+1,LB)
            Q(3,3) = CI(IC,JC,LNE)
            Q(3,5) = CI(IC,JC+1,LSE)
            Q(5,3) = CI(IC+1,JC,LNW)
            Q(5,5) = CI(IC+1,JC+1,LSW)

            !
            ! Apply fine-grid operator
            !
            DO JL = 2,6
               DO IL = 2,6

                  ILG = I + IL - 4
                  JLG = J + JL - 4

                  ILGP = MIN(ILG+1,IIF)
                  JLGP = MIN(JLG+1,JJF)

                  ILG = MIN(ILG,IIF)
                  JLG = MIN(JLG,JJF)

                  ILG = MAX(ILG,1)
                  JLG = MAX(JLG,1)

                  R(IL,JL) = SO(ILG,JLG,KO)*Q(IL,JL)
     &                     - SO(ILG,JLG,KW)*Q(IL-1,JL)
     &                     - SO(ILG,JLGP,KS)*Q(IL,JL+1)
     &                     - SO(ILGP,JLG,KW)*Q(IL+1,JL)
     &                     - SO(ILG,JLG,KS)*Q(IL,JL-1)
     &                     - SO(ILG,JLGP,KNW)*Q(IL-1,JL+1)
     &                     - SO(ILGP,JLGP,KSW)*Q(IL+1,JL+1)
     &                     - SO(ILGP,JLG,KNW)*Q(IL+1,JL-1)
     &                     - SO(ILG,JLG,KSW)*Q(IL-1,JL-1)


               ENDDO   
            ENDDO   
            !
            ! Apply residual weighting
            !
            SOC(IC,JC,KO) 
     &         =   BMG2_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, IC, JC, 4, 4 
     &                  )
            SOC(IC,JC,KW) 
     &         = - BMG2_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, IC-1, JC, 2, 4 
     &                  )
            SOC(IC,JC,KS) 
     &         = - BMG2_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, IC, JC-1, 4, 2 
     &                  )
            SOC(IC,JC,KSW) 
     &         = - BMG2_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, IC-1, JC-1, 2, 2 
     &                  )
            SOC(IC,JC+1,KNW) 
     &         = - BMG2_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, IC-1, JC+1, 2, 6
     &                  )
            !
	   ENDDO
C$OMP END PARALLEL DO

C     end of computation of kc difference operator when kf difference
C     operator is nine point
C******************************

         ELSE

C******************************
C   begin computation of kc difference operator when kf difference
C   operator is five point unless kc.ge.NOG
C

C$OMP PARALLEL DO
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,IC,IL,ILG,ILGP,J,JC,JL,JLG,JLGP,KK,Q,R)
C$OMP& SHARED(CI,IIC,IIF,JJC,JJF,SO,SOC)
	   DO KK=1,(IIC-2)*(JJC-2)

		IC=MOD(KK-1,IIC-2)+2
		I=2*(IC-1)

	      JC=(KK-1)/(IIC-2)+2
		J=2*(JC-1)

            DO JL = 1,7
               DO IL = 1,7
                  Q(IL,JL) = rZERO
                  R(IL,JL) = rZERO
               ENDDO
            ENDDO   
            !
            Q(4,4) = rONE
            Q(3,4) = CI(IC,JC,LR)
            Q(5,4) = CI(IC+1,JC,LL)
            Q(4,3) = CI(IC,JC,LA)
            Q(4,5) = CI(IC,JC+1,LB)
            Q(3,3) = CI(IC,JC,LNE)
            Q(3,5) = CI(IC,JC+1,LSE)
            Q(5,3) = CI(IC+1,JC,LNW)
            Q(5,5) = CI(IC+1,JC+1,LSW)
            !
            ! Apply fine-grid operator
            !
            DO JL = 2, 6
               DO IL = 2, 6

                  ILG = I + IL - 4
                  JLG = J + JL - 4

                  ILGP = MIN(ILG+1,IIF)
                  JLGP = MIN(JLG+1,JJF)

                  ILG = MIN(ILG,IIF)
                  JLG = MIN(JLG,JJF)

                  ILG = MAX(ILG,1)
                  JLG = MAX(JLG,1)

                  R(IL,JL) = SO(ILG,JLG,KO)*Q(IL,JL)
     &                     - SO(ILG,JLG,KW)*Q(IL-1,JL)
     &                     - SO(ILG,JLGP,KS)*Q(IL,JL+1)
     &                     - SO(ILGP,JLG,KW)*Q(IL+1,JL)
     &                     - SO(ILG,JLG,KS)*Q(IL,JL-1)

              ENDDO   
           ENDDO   
           !
           !  Apply residual weighting
           !
           SOC(IC,JC,KO) 
     &        =   BMG2_SymStd_UTILS_reswt(
     &                 CI, R, IIC, JJC, IC, JC, 4, 4 
     &                 )
           SOC(IC,JC,KW) 
     &        = - BMG2_SymStd_UTILS_reswt(
     &                 CI, R, IIC, JJC, IC-1, JC, 2, 4 
     &                 )
           SOC(IC,JC,KS) 
     &        = - BMG2_SymStd_UTILS_reswt(
     &                 CI, R, IIC, JJC, IC, JC-1, 4, 2
     &                 )
           SOC(IC,JC,KSW) 
     &        = - BMG2_SymStd_UTILS_reswt(
     &                 CI, R, IIC, JJC, IC-1, JC-1, 2, 2
     &                 )
           SOC(IC,JC+1,KNW) 
     &        = - BMG2_SymStd_UTILS_reswt(
     &                 CI, R, IIC, JJC, IC-1, JC+1, 2, 6 
     &                 )
           !
	   ENDDO
C$OMP END PARALLEL DO

      ENDIF

C   end of computation of grid kc difference operator, when kf
C   difference operator is five point
C******************************

      RETURN
      END
