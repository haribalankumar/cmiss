      SUBROUTINE BMG2_PerSymStd_residual(
     &                KF, KC, IMULT, SO, SOC, QF, QFC,
     &                Q, QC, CI, SOR, IIC, JJC, IIF, JJF,
     &                TAUERR, NStncl, NSORv, JPN 
     &                )

C
C***BEGIN PROLOGUE  BMG2_PerSymStd_residual
C***SUBSIDIARY
C***PURPOSE  BMG2_PerSymStd_residual computes weighted averages of 
C            residuals on the fine mesh to lay down on the coarse mesh. 
C            (The static residuals are computed in BMG2_PerSymStd_relax.)
C            The weights involve the transpose of the interpolation operator
C            from the coarse grid to the fine grid. In one mode 
C            BMG2_PerSymStd_residual computes a quantity to be used in 
C            estimating truncation error. BMG2_PerSymStd_residual also does
C            some prliminary setup for the interpolation routine 
C            BMG2_PerSymStd_interp_add.
C
C
C***LIBRARY   SLATEC
C***AUTHOR  DENDY, J. E. JR.
C             LOS ALAMOS NATIONAL LABORATORY
C           VOYTKO, M. H.
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C
C          Copyright, 1988. The Regents of the University of California.
C          This software was produced under a U.S. Government contract
C          (W-7405-ENG-36) by Los Alamso National Laboratory, which is
C          operated by the University of California for the U.S. Depart-
C          ment of Energy. The U.S. Government is licensed to use,
C          reproduce, and distribute this software. Permission is
C          granted to the public to copy and use this software without
C          charge, provided that this Notice and any statement of
C          authorship are reproduced on all copies. Neither the
C          Government nor the University makes any warranty, expressed
C          or implied, or assumes any liability or responsibility for
C          the use of this software.
C
C***PARAMETERS
C***INPUT
C   KF        KF is the fine grid number.
C   KC        KC is the coarse grid number.
C   IMULT     If IMULT.ne.0, MGRCAP computes a weighted average of
C             residuals on the fine grid to lay down on the
C             coarse grid. If IMULT.eq.0 MGRCAP computes
C             a quantity TAUERR used in estimating
C             truncation error.
C   SO        Refer to BOXMG.
C   SOC       SO on grid KC.
C   QF        Refer to BOXMG.
C   Q         Refer to BOXMG.
C   QC        Q on grid KC.
C   CI        Refer to BOXMG.
C   IIC       Number of points in x direction on coarse grid,
C             including two fictitious points.
C   JJC       Number of points in y direction on coarse grid,
C             including two fictitious points.
C   IIF       Number of grid points in x direction on fine grid,
C             including two fictitious points.
C   JJF       Number of grid points in y direction on fine grid,
C             including two fictitious points.
C   IPN       Refer to BOXMG.
C***INPUT/OUTPUT
C   SOR       Refer to BOXMG.
C***OUTPUT
C   QFC       QF on coarse grid.
C   TAUERR    TAUERR is a quantity used in estimating truncation
C             error.
C***ROUTINES CALLED  (NONE)
C
C***REVISION HISTORY  (YYMMDD)
C   1983/09/25  DATE WRITTEN
C   1990/06/27  - modified to conform to the 4/10/90 "Guide to the SLATEC
C                 Common Mathematical Library" by Victor A. Bandy
C   1999/12/06  - made do-loops end on enddo, changes gotos to if statements
C                 indented the code, replace (IF,JF)<--(I,J) by M.Berndt
C   2000/01/10  - updated declarations, mostly cosmetic 
C               - added INCLUDE to initialize constants
C               - fixed initialization of TAUERR so if no calculation
C                 is performed a rZERO value is returned.
C
C***END PROLOGUE  MGRCAP

C ==========================================================================

      IMPLICIT NONE

C ---------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'

C ----------------------------
C     Argument Declarations
C
      INTEGER IIC, IIF, IMULT, JJC, JJF, JPN, KC, KF, NSORv, NStncl
      REAL*8  CI(IIC,JJC,8), Q(IIF,JJF), QC(IIC,JJC), 
     &        QF(IIF,JJF), QFC(IIC,JJC), SO(IIF,JJF,NStncl), 
     &        SOC(IIC,JJC,5), SOR(IIF,JJF,NSORv), TAUERR

C ----------------------------
C     Local Declarations
C
      INTEGER I, IC, IIC1, IIC2, IIF1, IIF2, IPN, 
     &        J, JC, JJC1, JJC2, JJF1, JJF2
      REAL*8  XITAU, TEMP

C ==========================================================================

C
C***FIRST EXECUTABLE STATEMENT  MGRCAP
C
      IPN=IABS(JPN)
      XITAU=IMULT
      IIC1=IIC-1
      JJC1=JJC-1
      IIF1=IIF-1
      JJF1=JJF-1
      IIC2=IIC-2
      JJC2=JJC-2
      IIF2=IIF-2
      JJF2=JJF-2

      TAUERR=RZERO

      IF (IMULT.EQ.0) THEN
         DO J=2,JJF1
            DO I=2,IIF1
               SOR(I,J,MTOT)=QF(I,J)
            ENDDO
         ENDDO
      ENDIF
C
C   compute weighted averages of residuals on fine mesh to lay
C   down on coarse mesh. the weights are i(transpose) where i is the
C   the interpolation operator from grid kc to grid kf.
C
      IF(IPN.EQ.1.OR.IPN.EQ.3) THEN
         DO I=1,IIF
            SOR(I,1,MTOT)=SOR(I,JJF1,MTOT)
            SOR(I,JJF,MTOT)=SOR(I,2,MTOT)
         ENDDO
      ENDIF
      IF(IPN.EQ.2.OR.IPN.EQ.3) THEN
         DO J=1,JJF
            SOR(1,J,MTOT)=SOR(IIF1,J,MTOT)
            SOR(IIF,J,MTOT)=SOR(2,J,MTOT)
         ENDDO
      ENDIF
      J=0
      DO JC=2,JJC1
         J=J+2
         I=0
         DO IC=2,IIC1
            I=I+2
            QFC(IC,JC) = CI(IC-1,JC-1,LNE)*SOR(I-1,J-1,MTOT)
     &                 + CI(IC,JC,LA)*SOR(I,J-1,MTOT)
     &                 + CI(IC,JC-1,LNW)*SOR(I+1,J-1,MTOT)
     &                 + CI(IC,JC,LR)*SOR(I-1,J,MTOT)
     &                 + SOR(I,J,MTOT)
     &                 + CI(IC+1,JC,LL)*SOR(I+1,J,MTOT)
     &                 + CI(IC-1,JC,LSE)*SOR(I-1,J+1,MTOT)
     &                 + CI(IC,JC+1,LB)*SOR(I,J+1,MTOT)
     &                 + CI(IC,JC,LSW)*SOR(I+1,J+1,MTOT)
            
         ENDDO
      ENDDO
C
C   if imult=0, estimate truncation errors. otherwise, set up
C   quantities for mginadp.
C
      IF (IMULT.NE.0) THEN
         DO J=2,JJF1
            DO I=2,IIF1
               SOR(I,J,MTOT)=SOR(I,J,MTOT)/SO(I,J,KO)
            ENDDO
         ENDDO
      ELSE
         J=0
         DO JC=2,JJC1
            J=J+2
            I=0
            DO IC=2,IIC1
               I=I+2
               QC(IC,JC)=Q(I,J)
            ENDDO
         ENDDO
         IF(IPN.EQ.1.OR.IPN.EQ.3) THEN
            DO IC=1,IIC
               QC(IC,1)=QC(IC,JJC1)
               QC(IC,JJC)=QC(IC,2)
            ENDDO
         ENDIF
         IF(IPN.EQ.2.OR.IPN.EQ.3) THEN
            DO JC=1,JJC
               QC(IIC,JC)=QC(2,JC)
               QC(1,JC)=QC(IIC1,JC)
            ENDDO
         ENDIF
         DO JC=2,JJC1
            DO IC=2,IIC1
               TEMP=SOC(IC,JC,KW)*QC(IC-1,JC)+SOC(IC+1,JC,KW)
     &              *QC(IC+1,JC)+SOC(IC,JC,KS)*QC(IC,JC-1)
     &              +SOC(IC,JC+1,KS)*QC(IC,JC+1)+SOC(IC,JC,KSW)
     &              *QC(IC-1,JC-1)+SOC(IC+1,JC,KNW)*QC(IC+1,JC-1)
     &              +SOC(IC,JC+1,KNW)*QC(IC-1,JC+1)+SOC(IC+1,JC+1,KSW)
     &              *QC(IC+1,JC+1)-SOC(IC,JC,KO)*QC(IC,JC)
               QFC(IC,JC)=QFC(IC,JC)+TEMP
               TAUERR=TAUERR+QFC(IC,JC)**2
            ENDDO
         ENDDO
         TAUERR=SQRT(TAUERR)
      ENDIF
C     
      RETURN
      END
