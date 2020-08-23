      SUBROUTINE BMG2_SymStd_SETUP_interp_BI( 
     &                KF, KC, SO, SOC, CI, 
     &                IIF, JJF, IIC, JJC, NOG, IFD, NStncl, IBC, 
     &                IRELAX, BMG_IOFLAG, BMG_iPARMS
     &                )

C ==========================================================================
C
C   BMG2_SymStd_SETUP_interp_BI.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_SETUP_interp_BI.f constructs the standard bilinear
C   interpolation operator CI.  It is only being used at the time
C   for problems for which the opertors are 9-pt operators.  The 
C   operator CI interpolates a vector from the coarse grid KC, to 
C   the fine grid, KF. 
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Rewritten:  2005/01/11 (TMA)
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
      INTEGER   IC_STEP, JC_STEP, IC_BEG, IC_END, JC_BEG, JC_END

C --------------------------
C     Local Declarations:
C
      INTEGER   IC, I, IIC1, IICF, IICF1, IIF1, IIF2, JJF2,
     &          JC, J, JJC1, JJCF, JJCF1, JJF1, KK
      REAL*8    A, B, D1MACH, EP, EPSILON, SUM, S, OFFSUM

      REAL*8    rHALF, rFOURTH
      PARAMETER (rHALF=0.5D0, rFOURTH=0.25D0)

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

      IIF2=IIF-2
      JJF2=JJF-2

      IICF=(IIF-2)/2+3
      JJCF=(JJF-2)/2+3
      IICF1=IICF-1
      JJCF1=JJCF-1


C******************************
C   begin computation of i when kf difference operator is nine point
C

      IF ( IFD.NE.1 .OR. KF.LT.NOG ) THEN

         DO KK=1,(JJC1-1)*(IICF1-2)
            IC=MOD(KK-1,IICF1-2)+3
            I=2*(IC-1)
            JC=(KK-1)/(IICF1-2)+2
            J=2*(JC-1) 
            CI(IC,JC,LR)=rHALF
            CI(IC,JC,LL)=rHALF
         ENDDO

         DO KK=1,(JJCF1-2)*(IIC1-1)
            IC=MOD(KK-1,IIC1-1)+2
            I=2*(IC-1)
            JC=(KK-1)/(IIC1-1)+3
            J=2*(JC-1) 
            CI(IC,JC,LA)=rHALF
            CI(IC,JC,LB)=rHALF
         ENDDO

         DO KK=1,(JJCF1-2)*(IICF1-2)
            IC=MOD(KK-1,IICF1-2)+3
            I=2*(IC-1)
            JC=(KK-1)/(IICF1-2)+3
            J=2*(JC-1)
            CI(IC,JC,LSW)=rFOURTH
            CI(IC,JC,LSE)=rFOURTH
            CI(IC,JC,LNW)=rFOURTH
            CI(IC,JC,LNE)=rFOURTH
         ENDDO
	
C     end of computation of i when kf difference operator is nine point
C******************************

      ELSE

         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,*) 'ERROR: BMG2_SymStd_SETUP_interp_OI   .... '
            WRITE(*,*) '***** Bilinear interp not designed for 5-pt ops'
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,21)

      ENDIF

      I=2
      DO J=2,JJF1
         OFFSUM = DABS(SO(I  ,J  ,KW )) 
     &          + DABS(SO(I+1,J  ,KW )) 
     &          + DABS(SO(I  ,J  ,KS )) 
     &          + DABS(SO(I  ,J+1,KS )) 
     &          + DABS(SO(I  ,J  ,KSW)) 
     &          + DABS(SO(I+1,J  ,KNW)) 
     &          + DABS(SO(I  ,J+1,KNW)) 
     &          + DABS(SO(I+1,J+1,KSW))
         IF(OFFSUM.LE.1E-10) THEN
            print *,' and what is this?'
            IF(MOD(J,2).EQ.0) THEN
               IC=I/2+1
               JC=J/2+1
               CI(IC,JC,LA) = rZERO
               CI(IC,JC+1,LB) = rZERO
               CI(IC,JC,LR) = rZERO
               CI(IC+1,JC,LL) = rZERO
               CI(IC+1,JC,LNW) = rZERO
               CI(IC+1,JC+1,LSW) = rZERO
               CI(IC,JC+1,LSE) = rZERO
               CI(IC,JC,LNE) = rZERO
            ELSE
               IC=I/2+1
               JC=J/2+1
               CI(IC,JC,LA) = rZERO
               CI(IC,JC+1,LB) = rZERO
            ENDIF
         ENDIF
      ENDDO

      I=IIF1
      DO J=2,JJF1
         OFFSUM = DABS(SO(I  ,J  ,KW )) 
     &        + DABS(SO(I+1,J  ,KW )) 
     &        + DABS(SO(I  ,J  ,KS )) 
     &        + DABS(SO(I  ,J+1,KS )) 
     &        + DABS(SO(I  ,J  ,KSW)) 
     &        + DABS(SO(I+1,J  ,KNW)) 
     &        + DABS(SO(I  ,J+1,KNW)) 
     &        + DABS(SO(I+1,J+1,KSW)) 
         IF(OFFSUM.LE.1E-10) THEN
            print *,' and what is this?'
            IF(MOD(J,2).EQ.0) THEN
               IC=I/2+1
               JC=J/2+1
               CI(IC,JC,LA) = rZERO
               CI(IC,JC+1,LB) = rZERO
               CI(IC,JC,LR) = rZERO
               CI(IC+1,JC,LL) = rZERO
               CI(IC+1,JC,LNW) = rZERO
               CI(IC+1,JC+1,LSW) = rZERO
               CI(IC,JC+1,LSE) = rZERO
               CI(IC,JC,LNE) = rZERO
            ELSE
               IC=I/2+1
               JC=J/2+1
               CI(IC,JC,LA) = rZERO
               CI(IC,JC+1,LB) = rZERO
            ENDIF
         ENDIF
      ENDDO  

      J=2
      DO I=2,IIF1
         OFFSUM = DABS(SO(I  ,J  ,KW )) 
     &          + DABS(SO(I+1,J  ,KW )) 
     &          + DABS(SO(I  ,J  ,KS )) 
     &          + DABS(SO(I  ,J+1,KS )) 
     &          + DABS(SO(I  ,J  ,KSW)) 
     &          + DABS(SO(I+1,J  ,KNW)) 
     &          + DABS(SO(I  ,J+1,KNW)) 
     &          + DABS(SO(I+1,J+1,KSW)) 
         IF(OFFSUM.LE.1E-10) THEN
            print *,' and what is this?'
            IF(MOD(I,2).EQ.0) THEN
               IC=I/2+1
               JC=J/2+1
               CI(IC,JC,LR) = rZERO
               CI(IC+1,JC,LL) = rZERO
               CI(IC,JC,LA) = rZERO
               CI(IC,JC+1,LB) = rZERO
               CI(IC+1,JC,LNW) = rZERO
               CI(IC+1,JC+1,LSW) = rZERO
               CI(IC,JC+1,LSE) = rZERO
               CI(IC,JC,LNE) = rZERO
            ELSE
               IC=I/2+1
               JC=J/2+1
               CI(IC,JC,LR) = rZERO
               CI(IC+1,JC,LL) = rZERO
            ENDIF
         ENDIF
      ENDDO

      IF( (JJF/2).EQ.0 ) THEN
         print *,' what is this?'
         J=JJF1
         DO I=2,IIF1
            OFFSUM = DABS(SO(I  ,J  ,KW )) 
     &          + DABS(SO(I+1,J  ,KW )) 
     &          + DABS(SO(I  ,J  ,KS )) 
     &          + DABS(SO(I  ,J+1,KS )) 
     &          + DABS(SO(I  ,J  ,KSW)) 
     &          + DABS(SO(I+1,J  ,KNW)) 
     &          + DABS(SO(I  ,J+1,KNW)) 
     &          + DABS(SO(I+1,J+1,KSW)) 
            IF(OFFSUM.LE.1E-10) THEN
               IF(MOD(I,2).EQ.0) THEN
                  IC=I/2+1
                  JC=J/2+1
                  CI(IC,JC,LR) = rZERO
                  CI(IC+1,JC,LL) = rZERO
                  CI(IC,JC,LA) = rZERO
                  CI(IC,JC+1,LB) = rZERO
                  CI(IC+1,JC,LNW) = rZERO
                  CI(IC+1,JC+1,LSW) = rZERO
                  CI(IC,JC+1,LSE) = rZERO
                  CI(IC,JC,LNE) = rZERO
               ELSE
                  IC=I/2+1
                  JC=J/2+1
                  CI(IC,JC,LR) = rZERO
                  CI(IC+1,JC,LL) = rZERO
               ENDIF
            ENDIF
         ENDDO      
      ENDIF

      RETURN 
      END 
