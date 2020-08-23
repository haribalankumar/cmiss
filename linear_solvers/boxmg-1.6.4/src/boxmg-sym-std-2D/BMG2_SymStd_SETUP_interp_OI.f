      SUBROUTINE BMG2_SymStd_SETUP_interp_OI( 
     &                KF, KC, SO, SOC, CI, 
     &                IIF, JJF, IIC, JJC, NOG, IFD, NStncl, IBC, 
     &                IRELAX, BMG_IOFLAG, BMG_iPARMS
     &                )

C ==========================================================================
C
C   BMG2_SymStd_SETUP_interp_OI.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_SETUP_interp_OI.f constructs the operator-induced 
C   interpolation operator CI from the fine-grid stencil, SO.
C   The operator CI interpolates a vector from the coarse grid
C   KC, to the fine grid, KF. 
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Rewritten:  20045/01/11 (TMA)
C   - changes to loops to make use of OMP directives
C   Written:    2000/02/28 (JDM)
C   - taken, almost entirely from mgcoef.f 
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
      INTEGER   IBC, IIC, IIF, IFD, IRELAX, 
     &          JJC, JJF, KC, KF, NOG, NStncl
      REAL*8    CI(IIC,JJC,8), SO(IIF,JJF,NStncl), SOC(IIC,JJC,5) 

      LOGICAL   BMG_IOFLAG(NBMG_IOFLAG)
      INTEGER   BMG_iPARMS(NBMG_iPARMS)

C --------------------------
C     OpenMP Declarations:
C
	INTEGER   my_thread, my_ithread, my_jthread, nthreads, 
     &	    nithreads, njthreads
	INTEGER   IC_STEP, JC_STEP, IC_BEG, IC_END, JC_BEG, JC_END

C --------------------------
C     Local Declarations:
C
      INTEGER   IC, I, IIC1, IICF, IICF1, IICF2, IIF1, 
     &          JC, J, JJC1, JJCF, JJCF1, JJCF2, JJF1, KK
      REAL*8    A, B, D1MACH, EP, EPSILON, SUM, S

C ==========================================================================

C ----------------------------------
C     Sanity Check:
C ----------------------------------
      
      IF (KF-1.NE.KC ) THEN
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,*) 'ERROR: BMG2_SymStd_SETUP_interp_OI   .... '
            WRITE(*,*) '*****  KC = ', KC
            WRITE(*,*) '*****  KF = ', KF
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,19)
         RETURN

      ENDIF

C ----------------------------------
C     Useful Constants:
C ----------------------------------

      EPSILON = D1MACH(3)

C -----------------------------------
C     Useful indexing bounds:
C -----------------------------------

      IIC1=IIC-1
      JJC1=JJC-1

      IIF1=IIF-1
      JJF1=JJF-1

      IICF=(IIF-2)/2+3
      JJCF=(JJF-2)/2+3
      IICF1=IICF-1
      JJCF1=JJCF-1
      IICF2=IICF-2
      JJCF2=JJCF-2

C******************************
C   begin computation of i when kf difference operator is nine point
C

      IF ( IFD.NE.1 .OR. KF.LT.NOG ) THEN

C$OMP PARALLEL DO
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(A,B,I,IC,J,JC,SUM,EP),
C$OMP& SHARED(JJC1,IICF1,CI,SO,EPSILON)
         DO KK=1,(JJC1-1)*(IICF1-2)
            IC=MOD(KK-1,IICF1-2)+3
            I=2*(IC-1)
            JC=(KK-1)/(IICF1-2)+2
            J=2*(JC-1) 
            A=SO(I,J,KW)+SO(I,J,KNW)+SO(I,J+1,KSW)
            B=SO(I-1,J,KW)+SO(I-1,J,KSW)+SO(I-1,J+1,KNW)
            EP=MIN(ABS(A/SO(I-1,J,KO)),ABS(B/SO(I-1,J,KO)))
            SUM=A+B+SO(I-1,J,KS)+SO(I-1,J+1,KS)
            SUM=A+B+(SO(I-1,J,KO)-SUM)*MAX(SO(I-1,J,KO)
     &           -(rONE+EP)*SUM,rZERO)/(ABS(SO(I-1,J,KO)
     &           -(rONE+EP)*SUM)+EPSILON)
            SUM=rONE/SUM
            CI(IC,JC,LR)=A*SUM
            CI(IC,JC,LL)=B*SUM
         ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(A,B,I,IC,J,JC,SUM,EP),
C$OMP& SHARED(JJC1,IICF1,CI,SO,EPSILON)
         DO KK=1,(JJCF1-2)*(IIC1-1)
            IC=MOD(KK-1,IIC1-1)+2
            I=2*(IC-1)
            JC=(KK-1)/(IIC1-1)+3
            J=2*(JC-1) 
            A=SO(I,J,KS)+SO(I,J,KNW)+SO(I+1,J,KSW)
            B=SO(I,J-1,KS)+SO(I,J-1,KSW)+SO(I+1,J-1,KNW)
            EP=MIN(ABS(A/SO(I,J-1,KO)),ABS(B/SO(I,J-1,KO)))
            SUM=A+B+SO(I,J-1,KW)+SO(I+1,J-1,KW)
            SUM=A+B+(SO(I,J-1,KO)-SUM)*MAX(SO(I,J-1,KO)
     &           -(rONE+EP)*SUM,rZERO)/(ABS(SO(I,J-1,KO)
     &           -(rONE+EP)*SUM)+EPSILON)
            SUM=rONE/SUM
            CI(IC,JC,LA)=A*SUM
            CI(IC,JC,LB)=B*SUM
         ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(A,B,I,IC,J,JC,SUM,EP,S),
C$OMP& SHARED(JJCF1,IICF1,CI,SO,EPSILON)
         DO KK=1,(JJCF1-2)*(IICF1-2)
            IC=MOD(KK-1,IICF1-2)+3
            I=2*(IC-1)
            JC=(KK-1)/(IICF1-2)+3
            J=2*(JC-1) 
            SUM=SO(I-1,J-1,KW)+SO(I-1,J,KNW)+SO(I-1,J,KS)
     &           +SO(I,J,KSW)+SO(I,J-1,KW)+SO(I,J-1,KNW)
     &           +SO(I-1,J-1,KS)+SO(I-1,J-1,KSW)
            EP=MIN(ABS((SO(I-1,J-1,KSW)+SO(I-1,J-1,KW)
     &           +SO(I-1,J,KNW))/SO(I-1,J-1,KO)),
     &           ABS((SO(I-1,J,KNW)+SO(I-1,J,KS)
     &           +SO(I,J,KSW))/SO(I-1,J-1,KO)),
     &           ABS((SO(I,J,KSW)+SO(I,J-1,KW)
     &           +SO(I,J-1,KNW))/SO(I-1,J-1,KO)),
     &           ABS((SO(I,J-1,KNW)+SO(I-1,J-1,KS)
     &           +SO(I-1,J-1,KSW))/SO(I-1,J-1,KO))
     &           )
            SUM=SUM+(SO(I-1,J-1,KO)-SUM)*MAX(SO(I-1,J-1,KO)
     &           -(rONE+EP)*SUM,rZERO)/(ABS(SO(I-1,J-1,KO)
     &           -(rONE+EP)*SUM)+EPSILON)
            S=rONE/SUM
            CI(IC,JC,LSW)=(SO(I-1,J-1,KS)*CI(IC,JC-1,LL)
     &           +SO(I-1,J-1,KW)*CI(IC-1,JC,LB)
     &           +SO(I-1,J-1,KSW))*S
            CI(IC,JC,LSE)=(SO(I-1,J-1,KS)*CI(IC,JC-1,LR)
     &           +SO(I,J-1,KW)*CI(IC,JC,LB)
     &           +SO(I,J-1,KNW))*S
            CI(IC,JC,LNW)=(SO(I-1,J-1,KW)*CI(IC-1,JC,LA)
     &           +SO(I-1,J,KS)*CI(IC,JC,LL) 
     &           +SO(I-1,J,KNW))*S
            CI(IC,JC,LNE)=(SO(I-1,J,KS)*CI(IC,JC,LR)
     &           +SO(I,J-1,KW)*CI(IC,JC,LA)
     &           +SO(I,J,KSW))*S
         ENDDO
C$OMP END PARALLEL DO
	
C     end of computation of i when kf difference operator is nine point
C******************************

      ELSE

C******************************
C     begin computation of i when kf difference operator is five point
C
     
C$OMP PARALLEL DO
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(A,B,I,IC,J,JC,SUM,EP)
C$OMP& SHARED(JJC1,IICF1,CI,SO,EPSILON)
         DO KK=1,(JJC1-1)*(IICF1-2)
            IC=MOD(KK-1,IICF1-2)+3
            I=2*(IC-1)
            JC=(KK-1)/(IICF1-2)+2
            J=2*(JC-1) 
            A=SO(I,J,KW)
            B=SO(I-1,J,KW)
            EP=MIN(ABS(A/SO(I-1,J,KO)),ABS(B/SO(I-1,J,KO)))
            SUM=A+B+SO(I-1,J,KS)+SO(I-1,J+1,KS)
            SUM=A+B+(SO(I-1,J,KO)-SUM)*MAX(SO(I-1,J,KO)
     &           -(rONE+EP)*SUM,rZERO)/(ABS(SO(I-1,J,KO)
     &           -(rONE+EP)*SUM)+EPSILON)
            SUM=rONE/SUM
            CI(IC,JC,LR)=A*SUM
            CI(IC,JC,LL)=B*SUM
         ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(A,B,I,IC,J,JC,SUM,EP)
C$OMP& SHARED(JJC1,IICF1,CI,SO,EPSILON)
         DO KK=1,(JJCF1-2)*(IIC1-1)
            IC=MOD(KK-1,IIC1-1)+2
            I=2*(IC-1)
            JC=(KK-1)/(IIC1-1)+3
            J=2*(JC-1) 
            A=SO(I,J,KS)
            B=SO(I,J-1,KS)
            EP=MIN(ABS(A/SO(I,J-1,KO)),ABS(B/SO(I,J-1,KO)))
            SUM=A+B+SO(I,J-1,KW)+SO(I+1,J-1,KW)
            SUM=A+B+(SO(I,J-1,KO)-SUM)*MAX(SO(I,J-1,KO)
     &           -(rONE+EP)*SUM,rZERO)/(ABS(SO(I,J-1,KO)
     &           -(rONE+EP)*SUM)+EPSILON)
            SUM=rONE/SUM
            CI(IC,JC,LA)=A*SUM
            CI(IC,JC,LB)=B*SUM
         ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(A,B,I,IC,J,JC,SUM,EP,S)
C$OMP& SHARED(JJC1,IICF1,CI,SO,EPSILON)
         DO KK=1,(JJCF1-2)*(IICF1-2)
            IC=MOD(KK-1,IICF1-2)+3
            I=2*(IC-1)
            JC=(KK-1)/(IICF1-2)+3
            J=2*(JC-1) 
            SUM=SO(I-1,J-1,KW)+SO(I-1,J,KS)+SO(I,J-1,KW)
     &           +SO(I-1,J-1,KS)
            EP=MIN(ABS(SO(I-1,J-1,KW)/SO(I-1,J-1,KO)),
     &           ABS(SO(I-1,J,KS)/SO(I-1,J-1,KO)),
     &           ABS(SO(I,J-1,KW)/SO(I-1,J-1,KO)),
     &           ABS(SO(I-1,J-1,KS)/SO(I-1,J-1,KO))
     &           )
            SUM=SUM+(SO(I-1,J-1,KO)-SUM)*MAX(SO(I-1,J-1,KO)
     &           -(rONE+EP)*SUM,rZERO)/(ABS(SO(I-1,J-1,KO)
     &           -(rONE+EP)*SUM)+EPSILON)
            S=rONE/SUM
            CI(IC,JC,LSW)=(SO(I-1,J-1,KS)*CI(IC,JC-1,LL)
     &           + SO(I-1,J-1,KW)*CI(IC-1,JC,LB))*S
            CI(IC,JC,LSE)=(SO(I-1,J-1,KS)*CI(IC,JC-1,LR)
     &           + SO(I,J-1,KW)*CI(IC,JC,LB))*S
            CI(IC,JC,LNW)=(SO(I-1,J-1,KW)*CI(IC-1,JC,LA)
     &           + SO(I-1,J,KS)*CI(IC,JC,LL))*S
            CI(IC,JC,LNE)=(SO(I-1,J,KS)*CI(IC,JC,LR)
     &           + SO(I,J-1,KW)*CI(IC,JC,LA))*S
         ENDDO
C$OMP END PARALLEL DO

C     end computation of i when kf difference operator is five point
C*****************************
*
      ENDIF

      RETURN
      END
