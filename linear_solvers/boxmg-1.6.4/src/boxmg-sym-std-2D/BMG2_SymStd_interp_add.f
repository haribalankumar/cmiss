      SUBROUTINE BMG2_SymStd_interp_add(
     &                       KC, KF, Q ,QC, RES, SO, CI,
     &                       IIC, JJC, IIF, JJF, NStncl, IBC 
     &                       )

C ==========================================================================
C
C   BMG2_SymStd_interp_add.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_interp_add.f interpolates Q from the coarse mesh, KC, to 
C   the fine mesh, KF, and adds the result to Q on fine mesh.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2000/02/26 (JDM)
C   - taken, almost entirely from mgiad.f 
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
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER IBC, IIC, IIF, JJC, JJF, KC, KF, NStncl
      REAL*8  CI(IIC,JJC,8), Q(IIF,JJF), QC(IIC,JJC), 
     &        SO(IIF,JJF,NStncl), RES(IIF,JJF)

C ----------------------------
C     Local Declarations
C
      INTEGER IC, I, IICF, IICF1, IIC1, IIF1,
     &        JC, J, JJCF, JJCF1, JJC1, JJF1
	INTEGER kk
      REAL*8  A, AQ

	REAL*8 T1, T2

C ==========================================================================

C -------------------------------------------------
C     Useful index bounds:
C -------------------------------------------------

      IIF1=IIF-1
      JJF1=JJF-1

      IIC1=IIC-1
      JJC1=JJC-1

      IICF=(IIF-2)/2+3
      JJCF=(JJF-2)/2+3
      IICF1=IICF-1
      JJCF1=JJCF-1

C      DO i=1,JJC
C         DO j=1,IIC
C            IF( (i-1)*(i-IIC).EQ.0 .OR. (j-1)*(JJC-j).EQ.0 ) THEN 
C               QC(i,j) = 0.D0
C            ENDIF
C         ENDDO
C      ENDDO

C --------------------------------------------------
C     NB: division is only possible in the interior
C --------------------------------------------------

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,J,KK)
C$OMP& SHARED(IIF1,JJF1,RES,SO)
	DO KK = 1, (IIF1-1)*(JJF1-1)
	   I = mod(KK-1,IIF1-1)+2
	   J = (KK-1)/(IIF1-1)+2
	   RES(I,J)=RES(I,J)/SO(I,J,ko)
	ENDDO	
C$OMP END PARALLEL DO

C --------------------------------------------------
C   interpolate answers from coarse to fine mesh 
C   and add to answers on fine mesh.
C --------------------------------------------------

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(A,AQ,I,IC,J,JC,KK)
C$OMP& SHARED(CI,IICF1,JJCF1,Q,QC,RES)
	DO KK = 1, (IICF1-1)*(JJCF1-1)
	
	   IC = mod(KK-1,IICF1-1)+2
	   I = 2*(IC-1)
	   JC = (KK-1)/(IICF1-1)+2
	   J = 2*(JC-1)

	   Q(I,J) = Q(I,J) + QC(IC,JC)

	   IF( IC.NE.2 ) THEN	   
              A = CI(IC,JC,LR)*QC( IC ,JC) 
     &          + CI(IC,JC,LL)*QC(IC-1,JC)
              Q(I-1,J) = Q(I-1,J) + A + RES(I-1,J)
	   ENDIF

	   IF( JC.NE.2 ) THEN
              AQ = CI(IC,JC,LA)*QC(IC, JC ) 
     &           + CI(IC,JC,LB)*QC(IC,JC-1)
              Q(I,J-1) = Q(I,J-1) + AQ + RES(I,J-1)
	   ENDIF

	   IF( IC.NE.2 .AND. JC.NE.2 ) THEN
              A = CI(IC,JC,LSW)*QC(IC-1,JC-1)
     &          + CI(IC,JC,LNW)*QC(IC-1, JC )
     &          + CI(IC,JC,LNE)*QC( IC , JC )
     &          + CI(IC,JC,LSE)*QC( IC ,JC-1)
              Q(I-1,J-1) = Q(I-1,J-1) + A + RES(I-1,J-1)
	   ENDIF	

	ENDDO	
C$OMP END PARALLEL DO		  

C ---------------------------------------
C     Periodicity copying:
C ---------------------------------------

      IF ( ABS(IBC).EQ.1 .OR. ABS(IBC).EQ.3 ) THEN 
         DO i=1, IIF
            Q(i,1)   = Q(i,JJF1)
            Q(i,JJF) = Q(i,2)
         ENDDO
      ENDIF

      IF( ABS(IBC).EQ.2 .OR. ABS(IBC).EQ.3 ) THEN
         DO j=1, JJF
            Q(1,j)   = Q(IIF1,j)
            Q(IIF,j) = Q(2,j)
         ENDDO
      ENDIF



C ==========================================================================

C =======================

      RETURN
      END
