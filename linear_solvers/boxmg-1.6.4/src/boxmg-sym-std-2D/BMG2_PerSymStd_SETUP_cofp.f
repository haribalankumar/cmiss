*DECK MGCOFP
      SUBROUTINE BMG2_PerSymStd_SETUP_cofp(
     &                KF, KC, SO, SOC, QF, QFC, SOR, SORC, CI,
     &                IIF, JJF, IIC, JJC, ICOEF, IFD, NStncl,
     &                NSORv, IBC, IRELAX, ISTRT, ISKIP 
     &                )

C
C***BEGIN PROLOGUE  BMG2_PerSymStd_SETUP_cofp
C
C***PURPOSE  BMG2_PerSymStd_SETUP_cofp computes the interpolation operator
C            from the coarse grid to the fine grid. It also calculates
C            the difference operator on the coarse grid and computes
C            the right hand side for the coarse grid if required. Coef
C            will also compute lu decompositions if line relaxation is
C            specified.
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
C   KF        fine grid number
C   KC        coarse grid number
C   IIF       Number of grid points in x direction on fine grid,
C             including two fictitious points.
C   JJF       Number of grid points in y direction on fine grid,
C             including two fictitious points.
C   IIC       Number of grid points in x direction on coarse grid,
C             including two fictitious points.
C   JJC       Number of grid points in y direction on coarse grid,
C             including two fictitious points.
C   ICOEF     Refer to BOXMG.
C   IFD       Refer to BOXMG.
C   IRELAX    Refer to BOXMG.
C   ISTRT     Refer to BOXMG.
C   SO        Refer to BOXMG.
C   SOR       Refer to BOXMG.
C   CI        Refer to BOXMG.
C   ISKIP     Refer to BOXMG.
C   IPN       Refer to BOXMG.
C***OUTPUT
C   SOC       SO for coarse grid
C   QFC       QF for coarse grid
C   SOR       SOR for coarse grid
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830925  Joel E. Dendy
C           - DATE WRITTEN
C   900627  Victor A. Bandy
C           - modified to conform to the 4/10/90 "Guide to the SLATEC
C             Common Mathematical Library"
C   931022  J. David Moulton
C           - Bug Fixes:
C               1) The special case (Nxcf.NE.Nxc. AND. IPN.EQ.3) requires
C                  additional copying of the interpolation operator
C               2) With periodicity in x (IPN.EQ.2 .OR. IPN.EQ.3) a special
C                  case had been defined in the construction of SO(ic+1,jc,KNW)
C                  This was incorrect (refer also to bug fix #3)
C               3) The indexing in the computation of SOC(ic+1,jc,KNW) 
C                  caused an array index bound error. The calculation 
C                  has simply been shifted to SOC(ic,jc,KNW). As a result of
C                  this shift the copying of SOC(ic,jc,KNW) that is required
C                  for periodicity has also been modified. 
C
C***END PROLOGUE  MGCOFP
C
     
      IMPLICIT NONE

C ----------------------------
C     Includes
C ----------------------------

      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'

C ---------------------------
C    Argument Declarations:
C ---------------------------

      INTEGER   ICOEF, IIC, IIF, IFD, IRELAX, ISTRT, ISKIP, 
     &          JJC, JJF, KC, KF, NSORv, NStncl

      REAL*8    CI(IIC,JJC,8), QF(IIF,JJF), QFC(IIC,JJC),
     &          SO(IIF,JJF,NStncl), SOC(IIC,JJC,5), 
     &          SOR(IIF,JJF,NSORv), SORC(IIC,JJC,NSORv)

C --------------------------
C     Local Declarations:
C --------------------------

      INTEGER   IBEGC, IC, I, IIC1, IICF, IICF1, IICF2, IIF1, IIF2, 
     &          INDX, IPN, JBEGC, JC, J, JJC1, JJCF, JJCF1, JJCF2, 
     &          JJF1, JJF2, IBC, k
      REAL*8    A, B, CE, CEA, CENW, CFNW, CN, CNE, CNW, CO, COA, CONW, 
     &          COSW, CS, CSA, CSE, CSEA, CSENW, CSNW, CSWA, CSSW, CSW,
     &          CSWSW, CW, CWA, CWSW, D1MACH, EP, ONE, SUM, S, ZEPS,
     &          ZERO

C
C***FIRST EXECUTABLE STATEMENT  MGCOFP
C
      ZERO = 0
      ONE  = 1
      ZEPS = D1MACH(3)
C
      IPN=IABS(IBC)
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
      IICF2=IICF-2
      JJCF2=JJCF-2
C
C   for finest grid compute reciprocal of so and store in
C   sor.
C
      DO J=2,JJF1
         DO I=2,IIF1
            SOR(I,J,MSOR)=ONE/SO(I,J,KO)
            SOR(I,J,MTOT)=QF(I,J)*SOR(I,J,MSOR)
         ENDDO
      ENDDO
      IF( ISKIP.EQ.2 )GO TO 230
C
C   if kf difference operator is five point and kf.ge.icoef, go to 1600
C
      IF (IFD.EQ.1.AND.KF.GE.ICOEF) GO TO 120

C******************************
C   begin computation of i when kf difference operator is nine point
C
      J=0
      DO JC=2,JJC1
         J=J+2
         I=2
         IBEGC=3
         IF( IABS(IBC).NE.2 .AND. IABS(IBC).NE.3 ) GO TO 39
         I=0
         IBEGC=2
 39      DO IC=IBEGC,IICF1
            I=I+2
            A=SO(I,J,KW)+SO(I,J,KNW)+SO(I,J+1,KSW)
            B=SO(I-1,J,KW)+SO(I-1,J,KSW)+SO(I-1,J+1,KNW)
            EP = MIN(ABS(A),ABS(B),ONE)
            SUM=A+B+SO(I-1,J,KS)+SO(I-1,J+1,KS)
            SUM=A+B+(SO(I-1,J,KO)-SUM)*MAX(SO(I-1,J,KO)
     &           -(ONE+EP)*SUM,ZERO)/(ABS(SO(I-1,J,KO)
     &           -(ONE+EP)*SUM)+ZEPS)
            SUM=ONE/SUM
            CI(IC,JC,LR)=A*SUM
            CI(IC,JC,LL)=B*SUM
         ENDDO
      ENDDO

      IF( IIC.NE.IICF .OR. ( IABS(IBC).NE.2 .AND. IABS(IBC).NE.3) )THEN
         GO TO 45
      ENDIF

      I=3
      IC=IIC
      J=0
      DO JC=2,JJC1
         J=J+2
         A=SO(I,J,KW)+SO(I,J,KNW)+SO(I,J+1,KSW)
         B=SO(I-1,J,KW)+SO(I-1,J,KSW)+SO(I-1,J+1,KNW)
         EP=MIN(ABS(A),ABS(B),ONE)
         SUM=A+B+SO(I-1,J,KS)+SO(I-1,J+1,KS)
         SUM=A+B+(SO(I-1,J,KO)-SUM)*MAX(SO(I-1,J,KO)-(ONE+EP)*SUM,ZERO)
     &        /(ABS(SO(I-1,J,KO)-(ONE+EP)*SUM)+ZEPS)
         SUM=ONE/SUM
         CI(IC,JC,LR)=A*SUM
         CI(IC,JC,LL)=B*SUM
      ENDDO

   45 IF( IABS(IBC).NE.1 .AND. IABS(IBC).NE.3 ) GO TO 70

      IF( JJCF.EQ.JJC ) GO TO 53

      DO IC=2,IIC
         CI(IC,1,LL)=CI(IC,JJC1,LL)
         CI(IC,1,LR)=CI(IC,JJC1,LR)
         CI(IC,JJC,LL)=CI(IC,2,LL)
         CI(IC,JJC,LR)=CI(IC,2,LR)
      ENDDO
      GO TO 70

   53 J=JJF2
      JC=1
      INDX=0
      GO TO 55
   54 J=3
      JC=JJC
   55 IBEGC=3
      I=2

      IF( IABS(IBC).NE.2 .AND. IABS(IBC).NE.3 ) GO TO 57

      IBEGC=2
      I=0
   57 DO 56 IC=IBEGC,IICF1
      I=I+2
      A=SO(I,J,KW)+SO(I,J,KNW)+SO(I,J+1,KSW)
      B=SO(I-1,J,KW)+SO(I-1,J,KSW)+SO(I-1,J+1,KNW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(I-1,J,KS)+SO(I-1,J+1,KS)
      SUM=A+B+(SO(I-1,J,KO)-SUM)*MAX(SO(I-1,J,KO)-(ONE+EP)*SUM,ZERO)
     1/(ABS(SO(I-1,J,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LR)=A*SUM
      CI(IC,JC,LL)=B*SUM
   56 CONTINUE
      IF(IICF.NE.IIC.OR.(IABS(IBC).NE.2.AND.IABS(IBC).NE.3))GO TO 58
      I=3
      IC=IIC
      A=SO(I,J,KW)+SO(I,J,KNW)+SO(I,J+1,KSW)
      B=SO(I-1,J,KW)+SO(I-1,J,KSW)+SO(I-1,J+1,KNW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(I-1,J,KS)+SO(I-1,J+1,KS)
      SUM=A+B+(SO(I-1,J,KO)-SUM)*MAX(SO(I-1,J,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(I-1,J,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LR)=A*SUM
      CI(IC,JC,LL)=B*SUM
   58 IF(INDX.NE.0)GO TO 70
      INDX=1
      GO TO 54

   70 J=2
      JBEGC=3
      IF( IABS(IBC).NE.1 .AND. IABS(IBC).NE.3 ) GO TO 75

      J=0
      JBEGC=2
   75 DO 80 JC=JBEGC,JJCF1
         J=J+2
         I=0
         DO 79 IC=2,IIC1
            I=I+2
            A=SO(I,J,KS)+SO(I,J,KNW)+SO(I+1,J,KSW)
            B=SO(I,J-1,KS)+SO(I,J-1,KSW)+SO(I+1,J-1,KNW)
            SUM=A+B+SO(I,J-1,KW)+SO(I+1,J-1,KW)
            EP=MIN(ABS(A),ABS(B),ONE)
            SUM=A+B+(SO(I,J-1,KO)-SUM)*MAX(SO(I,J-1,KO)
     &           -(ONE+EP)*SUM,ZERO)
     &           /(ABS(SO(I,J-1,KO)-(ONE+EP)*SUM)+ZEPS)
            SUM=ONE/SUM
            CI(IC,JC,LA)=A*SUM
            CI(IC,JC,LB)=B*SUM
   79 CONTINUE

      IF(IICF.NE.IIC.OR.(IABS(IBC).NE.2.AND.IABS(IBC).NE.3))GO TO 59
      INDX=0
      I=IIF2
      IC=1
      GO TO 78
   77 I=3
      IC=IIC
   78 A=SO(I,J,KS)+SO(I,J,KNW)+SO(I+1,J,KSW)
      B=SO(I,J-1,KS)+SO(I,J-1,KSW)+SO(I+1,J-1,KNW)
      SUM=A+B+SO(I,J-1,KW)+SO(I+1,J-1,KW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+(SO(I,J-1,KO)-SUM)*MAX(SO(I,J-1,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(I,J-1,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LA)=A*SUM
      CI(IC,JC,LB)=B*SUM
      IF(INDX.NE.0)GO TO 80
      INDX=1
      GO TO 77
   59 IF(IICF.EQ.IIC.OR.(IABS(IBC).NE.2.AND.IABS(IBC).NE.3))GO TO 80
      CI(1,JC,LA)=CI(IIC1,JC,LA)
      CI(1,JC,LB)=CI(IIC1,JC,LB)
      CI(IIC,JC,LA)=CI(2,JC,LA)
      CI(IIC,JC,LB)=CI(2,JC,LB)
   80 CONTINUE
      IF(JJC.NE.JJCF.OR.(IABS(IBC).NE.1.AND.IABS(IBC).NE.3))GO TO 84
      J=3
      I=0
      DO 81 IC=2,IIC1
      I=I+2
      A=SO(I,J,KS)+SO(I,J,KNW)+SO(I+1,J,KSW)
      B=SO(I,J-1,KS)+SO(I,J-1,KSW)+SO(I+1,J-1,KNW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(I,J-1,KW)+SO(I+1,J-1,KW)
      SUM=A+B+(SO(I,J-1,KO)-SUM)*MAX(SO(I,J-1,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(I,J-1,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JJC,LA)=A*SUM
      CI(IC,JJC,LB)=B*SUM
   81 CONTINUE
      IF(IICF.NE.IIC.OR.IABS(IBC).NE.3)GO TO 835
      INDX=0
      I=IIF2
      IC=1
      GO TO 83
   82 I=3
      IC=IIC
   83 A=SO(I,J,KS)+SO(I,J,KNW)+SO(I+1,J,KSW)
      B=SO(I,J-1,KS)+SO(I,J-1,KSW)+SO(I+1,J-1,KNW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(I,J-1,KW)+SO(I+1,J-1,KW)
      SUM=A+B+(SO(I,J-1,KO)-SUM)*MAX(SO(I,J-1,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(I,J-1,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JJC,LA)=A*SUM
      CI(IC,JJC,LB)=B*SUM
      IF(INDX.NE.0)GO TO 84
      INDX=1
      GO TO 82
C---------------------------------------------------------------------------
C   bug fix #1(a) special case was missed
C
 835  IF (IICF.EQ.IIC .OR. IABS(IBC).NE.3) GOTO 84
      CI(1,JJC,LA)=CI(IIC1,JJC,LA)
      CI(1,JJC,LB)=CI(IIC1,JJC,LB)
      CI(IIC,JJC,LA)=CI(2,JJC,LA)
      CI(IIC,JJC,LB)=CI(2,JJC,LB)      
C----------------------------------------------------------------------------
   84 J=0
      JBEGC=2
      IF(IABS(IBC).NE.1.AND.IABS(IBC).NE.3)GO TO 85
      JBEGC=1
      J=-2
   85 DO 100 JC=JBEGC,JJCF2
      J=J+2
      I=0
      IBEGC=2
      IF(IABS(IBC).NE.2.AND.IABS(IBC).NE.3)GO TO 89
      IBEGC=1
      I=-2
   89 DO 90 IC=IBEGC,IICF2
      I=I+2
      SUM=SO(I+1,J+1,KW)+SO(I+1,J+2,KNW)+SO(I+1,J+2,KS)
     1+SO(I+2,J+2,KSW)+SO(I+2,J+1,KW)+SO(I+2,J+1,KNW)
     2+SO(I+1,J+1,KS)+SO(I+1,J+1,KSW)
      EP=MIN(ABS(SO(I+1,J+1,KSW)+SO(I+1,J+1,KW)+SO(I+1,J+2,KNW)),
     &       ABS(SO(I+1,J+2,KNW)+SO(I+1,J+2,KS)+SO(I+2,J+2,KSW)),
     &       ABS(SO(I+2,J+2,KSW)+SO(I+2,J+1,KW)+SO(I+2,J+1,KNW)),
     &       ABS(SO(I+2,J+1,KNW)+SO(I+1,J+1,KS)+SO(I+1,J+1,KSW)),
     &       ONE)
      SUM=SUM+(SO(I+1,J+1,KO)-SUM)*MAX(SO(I+1,J+1,KO)-(ONE+EP)*SUM,
     &    ZERO)/(ABS(SO(I+1,J+1,KO)-(ONE+EP)*SUM)+ZEPS)
      S=ONE/SUM
      CI(IC,JC,LSW)=(SO(I+1,J+1,KS)*CI(IC+1,JC,LL)+SO(I+1,J+1,KW)*
     & CI(IC,JC+1,LB)+SO(I+1,J+1,KSW))*S
      CI(IC,JC,LSE)=(SO(I+1,J+1,KS)*CI(IC+1,JC,LR)+SO(I+2,J+1,KW)*
     & CI(IC+1,JC+1,LB)+SO(I+2,J+1,KNW))*S
      CI(IC,JC,LNW)=(SO(I+1,J+1,KW)*CI(IC,JC+1,LA)+SO(I+1,J+2,KS)*
     &CI(IC+1,JC+1,LL) +SO(I+1,J+2,KNW))*S
      CI(IC,JC,LNE)=(SO(I+1,J+2,KS)*CI(IC+1,JC+1,LR)+SO(I+2,J+1,KW)*
     & CI(IC+1,JC+1,LA)+SO(I+2,J+2,KSW))*S
   90 CONTINUE
      IF(IIC.NE.IICF.OR.(IABS(IBC).NE.2.AND.IABS(IBC).NE.3))GO TO 100
      I=1
      IC=IIC1
      SUM=SO(I+1,J+1,KW)+SO(I+1,J+2,KNW)+SO(I+1,J+2,KS)
     &+SO(I+2,J+2,KSW)+SO(I+2,J+1,KW)+SO(I+2,J+1,KNW)
     &+SO(I+1,J+1,KS)+SO(I+1,J+1,KSW)
      EP=MIN(ABS(SO(I+1,J+1,KSW)+SO(I+1,J+1,KW)+SO(I+1,J+2,KNW)),
     &       ABS(SO(I+1,J+2,KNW)+SO(I+1,J+2,KS)+SO(I+2,J+2,KSW)),
     &       ABS(SO(I+2,J+2,KSW)+SO(I+2,J+1,KW)+SO(I+2,J+1,KNW)),
     &       ABS(SO(I+2,J+1,KNW)+SO(I+1,J+1,KS)+SO(I+1,J+1,KSW)),
     &       ONE)
      SUM=SUM+(SO(I+1,J+1,KO)-SUM)*MAX(SO(I+1,J+1,KO)-(ONE+EP)*SUM,
     &    ZERO)/(ABS(SO(I+1,J+1,KO)-(ONE+EP)*SUM)+ZEPS)
      S=ONE/SUM
      CI(IC,JC,LSW)=(SO(I+1,J+1,KS)*CI(IC+1,JC,LL)+SO(I+1,J+1,KW)*
     & CI(IC,JC+1,LB)+SO(I+1,J+1,KSW))*S
      CI(IC,JC,LSE)=(SO(I+1,J+1,KS)*CI(IC+1,JC,LR)+SO(I+2,J+1,KW)*
     & CI(IC+1,JC+1,LB)+SO(I+2,J+1,KNW))*S
      CI(IC,JC,LNW)=(SO(I+1,J+1,KW)*CI(IC,JC+1,LA)+SO(I+1,J+2,KS)*
     &CI(IC+1,JC+1,LL) +SO(I+1,J+2,KNW))*S
      CI(IC,JC,LNE)=(SO(I+1,J+2,KS)*CI(IC+1,JC+1,LR)+SO(I+2,J+1,KW)*
     & CI(IC+1,JC+1,LA)+SO(I+2,J+2,KSW))*S
  100 CONTINUE

      IF( JJC.NE.JJCF .OR. ( IABS(IBC).NE.1 .AND. IABS(IBC).NE.3 )) THEN
         GO TO 105
      ENDIF

      J=1
      JC=JJC1
      I=0
      IBEGC=2
      IF(IABS(IBC).NE.2.AND.IABS(IBC).NE.3)GO TO 103
      IBEGC=1
      I=-2
  103 DO 101 IC=IBEGC,IICF2
      I=I+2
      SUM=SO(I+1,J+1,KW)+SO(I+1,J+2,KNW)+SO(I+1,J+2,KS)
     &+SO(I+2,J+2,KSW)+SO(I+2,J+1,KW)+SO(I+2,J+1,KNW)
     &+SO(I+1,J+1,KS)+SO(I+1,J+1,KSW)
      EP=MIN(ABS(SO(I+1,J+1,KSW)+SO(I+1,J+1,KW)+SO(I+1,J+2,KNW)),
     &       ABS(SO(I+1,J+2,KNW)+SO(I+1,J+2,KS)+SO(I+2,J+2,KSW)),
     &       ABS(SO(I+2,J+2,KSW)+SO(I+2,J+1,KW)+SO(I+2,J+1,KNW)),
     &       ABS(SO(I+2,J+1,KNW)+SO(I+1,J+1,KS)+SO(I+1,J+1,KSW)),
     &       ONE)
      SUM=SUM+(SO(I+1,J+1,KO)-SUM)*MAX(SO(I+1,J+1,KO)-(ONE+EP)*SUM,
     &    ZERO)/(ABS(SO(I+1,J+1,KO)-(ONE+EP)*SUM)+ZEPS)
      S=ONE/SUM
      CI(IC,JC,LSW)=(SO(I+1,J+1,KS)*CI(IC+1,JC,LL)+SO(I+1,J+1,KW)*
     & CI(IC,JC+1,LB)+SO(I+1,J+1,KSW))*S
      CI(IC,JC,LSE)=(SO(I+1,J+1,KS)*CI(IC+1,JC,LR)+SO(I+2,J+1,KW)*
     & CI(IC+1,JC+1,LB)+SO(I+2,J+1,KNW))*S
      CI(IC,JC,LNW)=(SO(I+1,J+1,KW)*CI(IC,JC+1,LA)+SO(I+1,J+2,KS)*
     &CI(IC+1,JC+1,LL) +SO(I+1,J+2,KNW))*S
      CI(IC,JC,LNE)=(SO(I+1,J+2,KS)*CI(IC+1,JC+1,LR)+SO(I+2,J+1,KW)*
     & CI(IC+1,JC+1,LA)+SO(I+2,J+2,KSW))*S
  101 CONTINUE
      IF(IIC.NE.IICF.OR.IABS(IBC).NE.3)GO TO 105
      I=1
      IC=IIC1
      SUM=SO(I+1,J+1,KW)+SO(I+1,J+2,KNW)+SO(I+1,J+2,KS)
     &+SO(I+2,J+2,KSW)+SO(I+2,J+1,KW)+SO(I+2,J+1,KNW)
     &+SO(I+1,J+1,KS)+SO(I+1,J+1,KSW)
      EP=MIN(ABS(SO(I+1,J+1,KSW)+SO(I+1,J+1,KW)+SO(I+1,J+2,KNW)),
     &       ABS(SO(I+1,J+2,KNW)+SO(I+1,J+2,KS)+SO(I+2,J+2,KSW)),
     &       ABS(SO(I+2,J+2,KSW)+SO(I+2,J+1,KW)+SO(I+2,J+1,KNW)),
     &       ABS(SO(I+2,J+1,KNW)+SO(I+1,J+1,KS)+SO(I+1,J+1,KSW)),
     &       ONE)
      SUM=SUM+(SO(I+1,J+1,KO)-SUM)*MAX(SO(I+1,J+1,KO)-(ONE+EP)*SUM,
     &    ZERO)/(ABS(SO(I+1,J+1,KO)-(ONE+EP)*SUM)+ZEPS)
      S=ONE/SUM
      CI(IC,JC,LSW)=(SO(I+1,J+1,KS)*CI(IC+1,JC,LL)+SO(I+1,J+1,KW)*
     & CI(IC,JC+1,LB)+SO(I+1,J+1,KSW))*S
      CI(IC,JC,LSE)=(SO(I+1,J+1,KS)*CI(IC+1,JC,LR)+SO(I+2,J+1,KW)*
     & CI(IC+1,JC+1,LB)+SO(I+2,J+1,KNW))*S
      CI(IC,JC,LNW)=(SO(I+1,J+1,KW)*CI(IC,JC+1,LA)+SO(I+1,J+2,KS)*
     & CI(IC+1,JC+1,LL) +SO(I+1,J+2,KNW))*S
      CI(IC,JC,LNE)=(SO(I+1,J+2,KS)*CI(IC+1,JC+1,LR)+SO(I+2,J+1,KW)*
     & CI(IC+1,JC+1,LA)+SO(I+2,J+2,KSW))*S
  105 CONTINUE
C   end of computation of i when kf difference operator is nine point
C******************************
C   begin computation of grid kc difference operator when kf difference
C   operator is nine point unless kc. ge. icoef
C
      IF (KC.GE.ICOEF) GO TO 230

C-------------------------------------------------------------------------
C   Bug fix #2(a): With periodicity in x (IABS(IBC).EQ.2 .OR. IABS(IBC).EQ.3) a 
C                  special case had been defined in the construction 
C                  of SO(ic+1,jc,KNW).  This was incorrect (refer also 
C                  to bug fix #3)
C-------------------------------------------------------------------------

C-------------------------------------------------------------------------
C   Bug Fix #3(a): The computation of SOC(ic+1,jc,KNW) -> SOC(ic,jc,KNW)
C                  to prevent an array indexing bound error
C-------------------------------------------------------------------------

      J=0
      DO 110 JC=2,JJC1
      J=J+2
      I=0
      DO 109 IC=2,IIC1
      I=I+2
      CO=SO(I,J+1,KNW)*CI(IC-1,JC,LSW)+SO(I,J,KW)*CI(IC,JC,LL)+
     & SO(I,J,KSW)*CI(IC-1,JC-1,LNW)
      CS=SO(I,J-1,KW)*CI(IC-1,JC-1,LNW)+SO(I,J,KNW)*CI(IC,JC,LL)
      CSW=-SO(I-1,J-1,KO)*CI(IC-1,JC-1,LNW)+SO(I-1,J-1,KW)*
     &CI(IC-1,JC,LA)+SO(I-1,J,KNW)+SO(I-1,J,KS)*CI(IC,JC,LL)
      CW=-SO(I-1,J,KO)*CI(IC,JC,LL)+SO(I-1,J,KS)*CI(IC-1,JC-1,LNW)+
     & SO(I-1,J,KSW)*CI(IC-1,JC,LA)+SO(I-1,J,KW)+SO(I-1,J+1,KNW)*
     & CI(IC-1,JC+1,LB)+SO(I-1,J+1,KS)*CI(IC-1,JC,LSW)
      CNW=-SO(I-1,J+1,KO)*CI(IC-1,JC,LSW)+
     & SO(I-1,J+1,KS)*CI(IC,JC,LL)+SO(I-1,J+1,KSW)+
     & SO(I-1,J+1,KW)*CI(IC-1,JC+1,LB)
      CN=SO(I,J+1,KSW)*CI(IC,JC,LL)+SO(I,J+1,KW)*CI(IC-1,JC,LSW)
      SOC(IC,JC,KW)=CO+CI(IC,JC,LA)*CS+CI(IC-1,JC-1,LNE)*CSW+
     & CI(IC,JC,LR)*CW+CI(IC-1,JC,LSE)*CNW+CI(IC,JC+1,LB)*CN
      COSW=SO(I,J,KSW)*CI(IC-1,JC-1,LSW)
      CSSW=SO(I,J-1,KSW)*CI(IC,JC-1,LL)+
     & SO(I,J-1,KW)*CI(IC-1,JC-1,LSW)
      CSWSW=-SO(I-1,J-1,KO)*CI(IC-1,JC-1,LSW)+SO(I-1,J-1,KS)*
     & CI(IC,JC-1,LL)+SO(I-1,J-1,KSW)+SO(I-1,J-1,KW)*CI(IC-1,JC,LB)
      CWSW=SO(I-1,J,KS)*CI(IC-1,JC-1,LSW)+
     & SO(I-1,J,KSW)*CI(IC-1,JC,LB)
      SOC(IC,JC,KSW)=COSW+CI(IC,JC,LA)*CSSW+CI(IC-1,JC-1,LNE)*CSWSW+
     & CI(IC,JC,LR)*CWSW
      COA=SO(I,J,KSW)*CI(IC-1,JC-1,LSE)+SO(I,J,KS)*CI(IC,JC,LB)+
     &SO(I+1,J,KNW)*CI(IC,JC-1,LSW)
      CEA=SO(I+1,J,KS)*CI(IC,JC-1,LSW)+SO(I+1,J,KSW)*CI(IC,JC,LB)
      CSEA=-SO(I+1,J-1,KO)*CI(IC,JC-1,LSW)+SO(I+1,J-1,KS)*
     & CI(IC+1,JC-1,LL)+SO(I+1,J-1,KSW)+SO(I+1,J-1,KW)*CI(IC,JC,LB)
      CSA=-SO(I,J-1,KO)*CI(IC,JC,LB)+SO(I+1,J-1,KW)*CI(IC,JC-1,LSW)
     &+SO(I+1,J-1,LSW)*CI(IC+1,JC-1,LL)+SO(I,J-1,KS)+SO(I,J-1,KSW)
     &*CI(IC,JC-1,LR)+SO(I,J-1,KW)*CI(IC-1,JC-1,LSE)
      CSWA=-SO(I-1,J-1,KO)*CI(IC-1,JC-1,LSE)+SO(I,J-1,KW)*
     & CI(IC,JC,LB)+SO(I,J-1,KNW)+SO(I-1,J-1,KS)*CI(IC,JC-1,LR)
      CWA=SO(I,J,KNW)*CI(IC,JC,LB)+SO(I-1,J,KS)*CI(IC-1,JC-1,LSE)
      SOC(IC,JC,KS)=COA+CI(IC+1,JC,LL)*CEA+CI(IC,JC-1,LNW)*CSEA+
     & CI(IC,JC,LA)*CSA+CI(IC-1,JC-1,LNE)*CSWA+CI(IC,JC,LR)*CWA
      CONW=SO(I-1,J,KNW)*CI(IC-1,JC-1,LSE)
      CENW=SO(I,J,KNW)*CI(IC,JC,LB)+SO(I-1,J,KS)*CI(IC-1,JC-1,LSE)
      CSENW=-SO(I-1,J-1,KO)*CI(IC-1,JC-1,LSE)+SO(I,J-1,KW)*
     & CI(IC,JC,LB)+SO(I,J-1,KNW)+SO(I-1,J-1,KS)*CI(IC,JC-1,LR)
      CFNW=SO(I-1,J-1,KW)*CI(IC-1,JC-1,LSE)+SO(I-1,J-1,KNW)*
     & CI(IC,JC-1,LR)
      SOC(IC,JC,KNW)=CONW+CI(IC,JC,LL)*CENW+CI(IC-1,JC-1,LNW)*CSENW+
     & CI(IC-1,JC,LA)*CFNW
  109 CONTINUE
  110 CONTINUE

  118 J=0
      DO 112 JC=2,JJC1
      J=J+2
      I=0
      DO 111 IC=2,IIC1
      I=I+2
      CO=SO(I,J,KW)*CI(IC,JC,LR)+SO(I,J+1,KNW)*CI(IC-1,JC,LSE)+
     &SO(I,J+1,KS)*CI(IC,JC+1,LB)+SO(I+1,J+1,KSW)*CI(IC,JC,LSW)
     &+SO(I+1,J,KW)*CI(IC+1,JC,LL)+SO(I+1,J,KNW)*CI(IC,JC-1,LNW)
     &+SO(I,J,KS)*CI(IC,JC,LA)+SO(I,J,KSW)*CI(IC-1,JC-1,LNE)-
     & SO(I,J,KO)
      CW=-SO(I-1,J,KO)*CI(IC,JC,LR)+SO(I-1,J+1,KS)*CI(IC-1,JC,LSE)
     &+SO(I,J+1,KSW)*CI(IC,JC+1,LB)+SO(I,J,KW)+SO(I,J,KNW)*
     &CI(IC,JC,LA)+SO(I-1,J,KS)*CI(IC-1,JC-1,LNE)
      CNW=-SO(I-1,J+1,KO)*CI(IC-1,JC,LSE)+
     & SO(I,J+1,KW)*CI(IC,JC+1,LB)+SO(I,J+1,KNW)+
     & SO(I-1,J+1,KS)*CI(IC,JC,LR)
      CN=-SO(I,J+1,KO)*CI(IC,JC+1,LB)+SO(I+1,J+1,KNW)*CI(IC+1,JC,LL)
     & +SO(I,J+1,KS)+SO(I,J+1,KSW)*CI(IC,JC,LR)+SO(I,J+1,KW)*
     & CI(IC-1,JC,LSE)+SO(I+1,J+1,KW)*CI(IC,JC,LSW)
      CNE=-SO(I+1,J+1,KO)*CI(IC,JC,LSW)+
     & SO(I+1,J+1,KS)*CI(IC+1,JC,LL)+SO(I+1,J+1,KSW)+
     & SO(I+1,J+1,KW)*CI(IC,JC+1,LB)
      CE=-SO(I+1,J,KO)*CI(IC+1,JC,LL)+SO(I+1,J,KS)*CI(IC,JC-1,LNW)+
     & SO(I+1,J,KSW)*CI(IC,JC,LA)+SO(I+1,J,KW)+SO(I+1,J+1,KNW)*
     & CI(IC,JC+1,LB)+SO(I+1,J+1,KS)*CI(IC,JC,LSW)
      CSE=-SO(I+1,J-1,KO)*CI(IC,JC-1,LNW)+
     & SO(I+1,J-1,KW)*CI(IC,JC,LA)+SO(I+1,J,KNW)+
     & SO(I+1,J,KS)*CI(IC+1,JC,LL)
      CS=-SO(I,J-1,KO)*CI(IC,JC,LA)+SO(I,J-1,KW)*CI(IC-1,JC-1,LNE)+
     &SO(I,J,KNW)*CI(IC,JC,LR)+SO(I,J,KS)+SO(I+1,J,KSW)*
     &CI(IC+1,JC,LL)+SO(I+1,J-1,KW)*CI(IC,JC-1,LNW)
      CSW=-SO(I-1,J-1,KO)*CI(IC-1,JC-1,LNE)+
     & SO(I-1,J,KS)*CI(IC,JC,LR)+SO(I,J,KSW)+
     & SO(I,J-1,KW)*CI(IC,JC,LA)
      SOC(IC,JC,KO)=-CI(IC-1,JC,LSE)*CNW-CI(IC,JC+1,LB)*CN-CI(IC,JC,LSW)
     & *CNE-CI(IC,JC,LR)*CW-CO-CI(IC+1,JC,LL)*CE-CI(IC-1,JC-1,LNE)*CSW-
     & CI(IC,JC,LA)*CS-CI(IC,JC-1,LNW)*CSE
      SORC(IC,JC,MSOR)=ONE/SOC(IC,JC,KO)
  111 CONTINUE
  112 CONTINUE

C-------------------------------------------------------------------------
C     Bug Fix #3(aa): Copying of SOC(ic,jc,KNW) modified
C-------------------------------------------------------------------------

      IF(IABS(IBC).NE.1.AND.IABS(IBC).NE.3)GO TO 114
      DO 113 IC=1,IIC
      SOC(IC,JJC,KS)=SOC(IC,2,KS)
      SOC(IC,JJC,KW)=SOC(IC,2,KW)
      SOC(IC,JJC,KNW)=SOC(IC,2,KNW)
      SOC(IC,JJC,KSW)=SOC(IC,2,KSW)
      SOC(IC,JJC,KO)=SOC(IC,2,KO)
      SOC(IC,1,KO)=SOC(IC,JJC1,KO)
      SOC(IC,1,KW)=SOC(IC,JJC1,KW)
      SOC(IC,1,KS)=SOC(IC,JJC1,KS)
      SOC(IC,1,KSW)=SOC(IC,JJC1,KSW)
      SOC(IC,1,KNW)=SOC(IC,JJC1,KNW)
  113 CONTINUE
  114 IF(IABS(IBC).NE.2.AND.IABS(IBC).NE.3)GO TO 116
      DO 115 JC=1,JJC
      SOC(IIC,JC,KW)=SOC(2,JC,KW)
      SOC(IIC,JC,KS)=SOC(2,JC,KS)
      SOC(IIC,JC,KNW)=SOC(2,JC,KNW)
      SOC(IIC,JC,KSW)=SOC(2,JC,KSW)
      SOC(IIC,JC,KO)=SOC(2,JC,KO)
      SOC(1,JC,KW)=SOC(IIC1,JC,KW)
      SOC(1,JC,KS)=SOC(IIC1,JC,KS)
      SOC(1,JC,KSW)=SOC(IIC1,JC,KSW)
      SOC(1,JC,KNW)=SOC(IIC1,JC,KNW)
      SOC(1,JC,KO)=SOC(IIC1,JC,KO)
  115 CONTINUE
  116 CONTINUE
      GO TO 230
C   end of computation of kc difference operator when kf difference
C   operator is nine point
C******************************
C   begin computation of i when kf difference operator is five point
C
  120 J=0
      DO 160 JC=2,JJC1
      J=J+2
      I=2
      IBEGC=3
      IF(IABS(IBC).NE.2.AND.IABS(IBC).NE.3)GO TO 149
      I=0
      IBEGC=2
  149 DO 150 IC=IBEGC,IICF1
      I=I+2
      A=SO(I,J,KW)
      B=SO(I-1,J,KW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(I-1,J,KS)+SO(I-1,J+1,KS)
      SUM=A+B+(SO(I-1,J,KO)-SUM)*MAX(SO(I-1,J,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(I-1,J,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LR)=A*SUM
      CI(IC,JC,LL)=B*SUM
  150 CONTINUE
  160 CONTINUE
      IF(IIC.NE.IICF.OR.(IABS(IBC).NE.2.AND.IABS(IBC).NE.3))GO TO 155
      I=3
      IC=IIC
      J=0
      DO 153 JC=2,JJC1
      J=J+2
      A=SO(I,J,KW)
      B=SO(I-1,J,KW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(I-1,J,KS)+SO(I-1,J+1,KS)
      SUM=A+B+(SO(I-1,J,KO)-SUM)*MAX(SO(I-1,J,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(I-1,J,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LR)=A*SUM
      CI(IC,JC,LL)=B*SUM
  153 CONTINUE
  155 IF(IABS(IBC).NE.1 .AND. IABS(IBC).NE.3)GO TO 180
      IF(JJCF.EQ.JJC)GO TO 163
      DO 161 IC=2,IIC
      CI(IC,1,LL)=CI(IC,JJC1,LL)
      CI(IC,1,LR)=CI(IC,JJC1,LR)
      CI(IC,JJC,LL)=CI(IC,2,LL)
      CI(IC,JJC,LR)=CI(IC,2,LR)
  161 CONTINUE
      GO TO 180
  163 J=JJF2
      JC=1
      INDX=0
      GO TO 165
  164 J=3
      JC=JJC
  165 IBEGC=3
      I=2
      IF(IABS(IBC).NE.2.AND.IABS(IBC).NE.3)GO TO 167
      IBEGC=2
      I=0
  167 DO 166 IC=IBEGC,IICF1
      I=I+2
      A=SO(I,J,KW)
      B=SO(I-1,J,KW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(I-1,J,KS)+SO(I-1,J+1,KS)
      SUM=A+B+(SO(I-1,J,KO)-SUM)*MAX(SO(I-1,J,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(I-1,J,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LR)=A*SUM
      CI(IC,JC,LL)=B*SUM
  166 CONTINUE
      IF(IICF.NE.IIC.OR.(IABS(IBC).NE.2.AND.IABS(IBC).NE.3))GO TO 168
      I=3
      IC=IIC
      A=SO(I,J,KW)
      B=SO(I-1,J,KW)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(I-1,J,KS)+SO(I-1,J+1,KS)
      SUM=A+B+(SO(I-1,J,KO)-SUM)*MAX(SO(I-1,J,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(I-1,J,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LR)=A*SUM
      CI(IC,JC,LL)=B*SUM
  168 IF(INDX.NE.0)GO TO 180
      INDX=1
      GO TO 164

  180 J=2
      JBEGC=3
      IF( IABS(IBC).NE.1 .AND. IABS(IBC).NE.3 ) GO TO 175
      J=0
      JBEGC=2
  175 DO 190 JC=JBEGC,JJCF1
      J=J+2
      I=0

      DO 189 IC=2,IIC1
      I=I+2
      A=SO(I,J,KS)
      B=SO(I,J-1,KS)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(I,J-1,KW)+SO(I+1,J-1,KW)
      SUM=A+B+(SO(I,J-1,KO)-SUM)*MAX(SO(I,J-1,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(I,J-1,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LA)=A*SUM
      CI(IC,JC,LB)=B*SUM
  189 CONTINUE
      IF(IICF.NE.IIC.OR.(IABS(IBC).NE.2.AND.IABS(IBC).NE.3))GO TO 169
      INDX=0
      I=IIF2
      IC=1
      GO TO 188
  187 I=3
      IC=IIC
  188 A=SO(I,J,KS)
      B=SO(I,J-1,KS)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(I,J-1,KW)+SO(I+1,J-1,KW)
      SUM=A+B+(SO(I,J-1,KO)-SUM)*MAX(SO(I,J-1,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(I,J-1,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JC,LA)=A*SUM
      CI(IC,JC,LB)=B*SUM
      IF(INDX.NE.0)GO TO 190
      INDX=1
      GO TO 187
  169 IF(IICF.EQ.IIC.OR.(IABS(IBC).NE.2.AND.IABS(IBC).NE.3))GO TO 190
      CI(1,JC,LA)=CI(IIC1,JC,LA)
      CI(1,JC,LB)=CI(IIC1,JC,LB)
      CI(IIC,JC,LA)=CI(2,JC,LA)
      CI(IIC,JC,LB)=CI(2,JC,LB)
  190 CONTINUE
      IF(JJC.NE.JJCF.OR.(IABS(IBC).NE.1.AND.IABS(IBC).NE.3))GO TO 194
      J=3
      I=0
      DO 191 IC=2,IIC1
      I=I+2
      A=SO(I,J,KS)
      B=SO(I,J-1,KS)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(I,J-1,KW)+SO(I+1,J-1,KW)
      SUM=A+B+(SO(I,J-1,KO)-SUM)*MAX(SO(I,J-1,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(I,J-1,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JJC,LA)=A*SUM
      CI(IC,JJC,LB)=B*SUM
  191 CONTINUE
      IF(IICF.NE.IIC.OR.IABS(IBC).NE.3)GO TO 1935
      INDX=0
      I=IIF2
      IC=1
      GO TO 193
  192 I=3
      IC=IIC
  193 A=SO(I,J,KS)
      B=SO(I,J-1,KS)
      EP=MIN(ABS(A),ABS(B),ONE)
      SUM=A+B+SO(I,J-1,KW)+SO(I+1,J-1,KW)
      SUM=A+B+(SO(I,J-1,KO)-SUM)*MAX(SO(I,J-1,KO)-(ONE+EP)*SUM,ZERO)
     &   /(ABS(SO(I,J-1,KO)-(ONE+EP)*SUM)+ZEPS)
      SUM=ONE/SUM
      CI(IC,JJC,LA)=A*SUM
      CI(IC,JJC,LB)=B*SUM
      IF(INDX.NE.0)GO TO 194
      INDX=1
      GO TO 192

C------------------------------------------------------------------------
C     Bug fix #1(b): special was case missed
C
 1935 IF (IICF.EQ.IIC .OR. IABS(IBC).NE.3) GOTO 194
      CI(1,JJC,LA)=CI(IIC1,JJC,LA)
      CI(1,JJC,LB)=CI(IIC1,JJC,LB)
      CI(IIC,JJC,LA)=CI(2,JJC,LA)
      CI(IIC,JJC,LB)=CI(2,JJC,LB)      
C------------------------------------------------------------------------

  194 J=0
      JBEGC=2
      IF(IABS(IBC).NE.1.AND.IABS(IBC).NE.3)GO TO 195
      JBEGC=1
      J=-2
  195 DO 210 JC=JBEGC,JJCF2
      J=J+2
      I=0
      IBEGC=2
      IF(IABS(IBC).NE.2.AND.IABS(IBC).NE.3)GO TO 209
      IBEGC=1
      I=-2
  209 DO 200 IC=IBEGC,IICF2
         I=I+2
         SUM = SO(I+1,J+1,KW)
     &       + SO(I+1,J+2,KS)
     &       + SO(I+2,J+1,KW)
     &       + SO(I+1,J+1,KS)
         EP=MIN(ABS(SO(I+1,J+1,KW)),ABS(SO(I+1,J+2,KS)),
     &        ABS(SO(I+2,J+1,KW)),ABS(SO(I+1,J+1,KS)),ONE)
         SUM=SUM+(SO(I+1,J+1,KO)-SUM)*MAX(SO(I+1,J+1,KO)
     &        - (ONE+EP)*SUM,ZERO)/(ABS(SO(I+1,J+1,KO)
     &        - (ONE+EP)*SUM)+ZEPS)
         S=ONE/SUM
         !PRINT *, IC,JC, ' : ', I+1,J+1,IC+1,JC
         CI(IC,JC,LSW)=(SO(I+1,J+1,KS)*CI(IC+1,JC,LL)
     &                + SO(I+1,J+1,KW)*CI(IC,JC+1,LB))*S
         CI(IC,JC,LSE)=(SO(I+1,J+1,KS)*CI(IC+1,JC,LR)
     &                + SO(I+2,J+1,KW)*CI(IC+1,JC+1,LB))*S
         CI(IC,JC,LNW)=(SO(I+1,J+1,KW)*CI(IC,JC+1,LA)
     &                + SO(I+1,J+2,KS)*CI(IC+1,JC+1,LL))*S
         CI(IC,JC,LNE)=(SO(I+1,J+2,KS)*CI(IC+1,JC+1,LR)
     &                + SO(I+2,J+1,KW)*CI(IC+1,JC+1,LA))*S
  200 CONTINUE

      IF(IIC.NE.IICF.OR.(IABS(IBC).NE.2.AND.IABS(IBC).NE.3))GO TO 210
      I=1
      IC=IIC1
      SUM=SO(I+1,J+1,KW)+SO(I+1,J+2,KS)+SO(I+2,J+1,KW)
     &+SO(I+1,J+1,KS)
      EP=MIN(ABS(SO(I+1,J+1,KW)),ABS(SO(I+1,J+2,KS)),
     &       ABS(SO(I+2,J+1,KW)),ABS(SO(I+1,J+1,KS)),ONE)
      SUM=SUM+(SO(I+1,J+1,KO)-SUM)*MAX(SO(I+1,J+1,KO)-(ONE+EP)*SUM,
     &    ZERO)/(ABS(SO(I+1,J+1,KO)-(ONE+EP)*SUM)+ZEPS)
      S=ONE/SUM
      CI(IC,JC,LSW)=(SO(I+1,J+1,KS)*CI(IC+1,JC,LL)+SO(I+1,J+1,KW)*
     & CI(IC,JC+1,LB))*S
      CI(IC,JC,LSE)=(SO(I+1,J+1,KS)*CI(IC+1,JC,LR)+SO(I+2,J+1,KW)*
     & CI(IC+1,JC+1,LB))*S
      CI(IC,JC,LNW)=(SO(I+1,J+1,KW)*CI(IC,JC+1,LA)+SO(I+1,J+2,KS)*
     & CI(IC+1,JC+1,LL))*S
      CI(IC,JC,LNE)=(SO(I+1,J+2,KS)*CI(IC+1,JC+1,LR)+SO(I+2,J+1,KW)*
     & CI(IC+1,JC+1,LA))*S
  210 CONTINUE
      IF(JJC.NE.JJCF.OR.(IABS(IBC).NE.1.AND.IABS(IBC).NE.3))GO TO 215
      J=1
      JC=JJC1
      I=0
      IBEGC=2
      IF(IABS(IBC).NE.2.AND.IABS(IBC).NE.3)GO TO 205
      IBEGC=1
      I=-2
  205 DO 211 IC=IBEGC,IICF2
      I=I+2
      SUM=SO(I+1,J+1,KW)+SO(I+1,J+2,KS)+SO(I+2,J+1,KW)
     &+SO(I+1,J+1,KS)
      EP=MIN(ABS(SO(I+1,J+1,KW)),ABS(SO(I+1,J+2,KS)),
     &       ABS(SO(I+2,J+1,KW)),ABS(SO(I+1,J+1,KS)),ONE)
      SUM=SUM+(SO(I+1,J+1,KO)-SUM)*MAX(SO(I+1,J+1,KO)-(ONE+EP)*SUM,
     &    ZERO)/(ABS(SO(I+1,J+1,KO)-(ONE+EP)*SUM)+ZEPS)
      S=ONE/SUM
      CI(IC,JC,LSW)=(SO(I+1,J+1,KS)*CI(IC+1,JC,LL)+SO(I+1,J+1,KW)*
     & CI(IC,JC+1,LB))*S
      CI(IC,JC,LSE)=(SO(I+1,J+1,KS)*CI(IC+1,JC,LR)+SO(I+2,J+1,KW)*
     & CI(IC+1,JC+1,LB))*S
      CI(IC,JC,LNW)=(SO(I+1,J+1,KW)*CI(IC,JC+1,LA)+SO(I+1,J+2,KS)*
     & CI(IC+1,JC+1,LL))*S
      CI(IC,JC,LNE)=(SO(I+1,J+2,KS)*CI(IC+1,JC+1,LR)+SO(I+2,J+1,KW)*
     & CI(IC+1,JC+1,LA))*S
  211 CONTINUE
      IF(IIC.NE.IICF.OR.IABS(IBC).NE.3)GO TO 215
      I=1
      IC=IIC1
      SUM=SO(I+1,J+1,KW)+SO(I+1,J+2,KS)+SO(I+2,J+1,KW)
     &+SO(I+1,J+1,KS)
      EP=MIN(ABS(SO(I+1,J+1,KW)),ABS(SO(I+1,J+2,KS)),
     &       ABS(SO(I+2,J+1,KW)),ABS(SO(I+1,J+1,KS)),ONE)
      SUM=SUM+(SO(I+1,J+1,KO)-SUM)*MAX(SO(I+1,J+1,KO)-(ONE+EP)*SUM,
     &    ZERO)/(ABS(SO(I+1,J+1,KO)-(ONE+EP)*SUM)+ZEPS)
      S=ONE/SUM
      CI(IC,JC,LSW)=(SO(I+1,J+1,KS)*CI(IC+1,JC,LL)+SO(I+1,J+1,KW)*
     1 CI(IC,JC+1,LB))*S
      CI(IC,JC,LSE)=(SO(I+1,J+1,KS)*CI(IC+1,JC,LR)+SO(I+2,J+1,KW)*
     1 CI(IC+1,JC+1,LB))*S
      CI(IC,JC,LNW)=(SO(I+1,J+1,KW)*CI(IC,JC+1,LA)+SO(I+1,J+2,KS)*
     1 CI(IC+1,JC+1,LL))*S
      CI(IC,JC,LNE)=(SO(I+1,J+2,KS)*CI(IC+1,JC+1,LR)+SO(I+2,J+1,KW)*
     1 CI(IC+1,JC+1,LA))*S
  215 CONTINUE
C   end computation of i when kf difference operator is five point

C******************************
C   begin computation of kc difference operator when kf difference
C   operator is five point unless kc.ge.icoef
C
      IF (KC.GE.ICOEF) GO TO 230

C-------------------------------------------------------------------------------
C     Bug Fix 2(b): With periodicity in x 
C                   (IABS(IBC).EQ.2 .OR. IABS(IBC).EQ.3) a special
C                   case had been defined in the construction of SO(ic+1,jc,KNW)
C                   This was incorrect (refer also to bug fix #3)
C-------------------------------------------------------------------------------

C-------------------------------------------------------------------------------
C     Bug Fix 3(b): The computation of SOC(ic+1,jc,KNW) -> SOC(ic,jc,KNW)
C                   to prevent an array indexing bound error
C-------------------------------------------------------------------------------

      J=0
      DO JC=2,JJC1
         J=J+2
         I=0
         DO IC=2,IIC1
            I=I+2

            CO=SO(I,J,KW)*CI(IC,JC,LL)
            CS=SO(I,J-1,KW)*CI(IC-1,JC-1,LNW)
            CSW=-SO(I-1,J-1,KO)*CI(IC-1,JC-1,LNW)+SO(I-1,J-1,KW)
     &           *CI(IC-1,JC,LA)+SO(I-1,J,KS)*CI(IC,JC,LL)
            CW=-SO(I-1,J,KO)*CI(IC,JC,LL)+SO(I-1,J,KS)
     &           *CI(IC-1,JC-1,LNW)+SO(I-1,J,KW)+SO(I-1,J+1,KS)
     &           *CI(IC-1,JC,LSW)
            CNW=-SO(I-1,J+1,KO)*CI(IC-1,JC,LSW)+SO(I-1,J+1,KS)
     &           *CI(IC,JC,LL)+SO(I-1,J+1,KW)*CI(IC-1,JC+1,LB)
            CN=SO(I,J+1,KW)*CI(IC-1,JC,LSW)
            SOC(IC,JC,KW)=CO+CI(IC,JC,LA)*CS+CI(IC-1,JC-1,LNE)*CSW
     &           +CI(IC,JC,LR)*CW+CI(IC-1,JC,LSE)*CNW
     &           +CI(IC,JC+1,LB)*CN

            CSSW=SO(I,J-1,KW)*CI(IC-1,JC-1,LSW)
            CSWSW=-SO(I-1,J-1,KO)*CI(IC-1,JC-1,LSW)
     &           +SO(I-1,J-1,KS)*CI(IC,JC-1,LL)+SO(I-1,J-1,KW)
     &           *CI(IC-1,JC,LB)
            CWSW=SO(I-1,J,KS)*CI(IC-1,JC-1,LSW)
            SOC(IC,JC,KSW)=CI(IC,JC,LA)*CSSW+CI(IC-1,JC-1,LNE)
     &           *CSWSW+CI(IC,JC,LR)*CWSW

            COA=SO(I,J,KS)*CI(IC,JC,LB)
            CEA=SO(I+1,J,KS)*CI(IC,JC-1,LSW)
            CSEA=-SO(I+1,J-1,KO)*CI(IC,JC-1,LSW)+SO(I+1,J-1,KS)
     &           *CI(IC+1,JC-1,LL)+SO(I+1,J-1,KW)*CI(IC,JC,LB)
            CSA=-SO(I,J-1,KO)*CI(IC,JC,LB)+SO(I+1,J-1,KW)
     &           *CI(IC,JC-1,LSW)+SO(I,J-1,KS)+SO(I,J-1,KW)
     &           *CI(IC-1,JC-1,LSE)
            CSWA=-SO(I-1,J-1,KO)*CI(IC-1,JC-1,LSE)+SO(I,J-1,KW)
     &           *CI(IC,JC,LB)+SO(I-1,J-1,KS)*CI(IC,JC-1,LR)
            CWA=SO(I-1,J,KS)*CI(IC-1,JC-1,LSE)
            SOC(IC,JC,KS)=COA+CI(IC+1,JC,LL)*CEA+CI(IC,JC-1,LNW)
     &           *CSEA+CI(IC,JC,LA)*CSA+CI(IC-1,JC-1,LNE)*CSWA
     &           +CI(IC,JC,LR)*CWA

            CENW=SO(I-1,J,KS)*CI(IC-1,JC-1,LSE)
            CSENW=-SO(I-1,J-1,KO)*CI(IC-1,JC-1,LSE)+SO(I,J-1,KW)*
     &           CI(IC,JC,LB)+SO(I-1,J-1,KS)*CI(IC,JC-1,LR)
            CSNW=SO(I-1,J-1,KW)*CI(IC-1,JC-1,LSE)
            SOC(IC,JC,KNW)=CI(IC,JC,LL)*CENW+CI(IC-1,JC-1,LNW)*CSENW+
     &           CI(IC-1,JC,LA)*CSNW

            CO=SO(I,J,KW)*CI(IC,JC,LR)+SO(I,J+1,KS)
     &           *CI(IC,JC+1,LB)+SO(I+1,J,KW)*CI(IC+1,JC,LL)
     &           +SO(I,J,KS)*CI(IC,JC,LA)-SO(I,J,KO)
            CW=-SO(I-1,J,KO)*CI(IC,JC,LR)+SO(I-1,J+1,KS)
     &           *CI(IC-1,JC,LSE)+SO(I,J,KW)+SO(I-1,J,KS)
     &           *CI(IC-1,JC-1,LNE)
            CNW=-SO(I-1,J+1,KO)*CI(IC-1,JC,LSE)+SO(I,J+1,KW)
     &           *CI(IC,JC+1,LB)+SO(I-1,J+1,KS)*CI(IC,JC,LR)
            CN=-SO(I,J+1,KO)*CI(IC,JC+1,LB)+SO(I,J+1,KS)+
     &           SO(I,J+1,KW)*CI(IC-1,JC,LSE)+SO(I+1,J+1,KW)
     &           *CI(IC,JC,LSW)
            CNE=-SO(I+1,J+1,KO)*CI(IC,JC,LSW)+SO(I+1,J+1,KS)
     &           *CI(IC+1,JC,LL)+SO(I+1,J+1,KW)*CI(IC,JC+1,LB)
            CE=-SO(I+1,J,KO)*CI(IC+1,JC,LL)+SO(I+1,J,KS)
     &           *CI(IC,JC-1,LNW)+SO(I+1,J,KW)+SO(I+1,J+1,KS)
     &           *CI(IC,JC,LSW)
            CSE=-SO(I+1,J-1,KO)*CI(IC,JC-1,LNW)+
     &           SO(I+1,J-1,KW)*CI(IC,JC,LA)+SO(I+1,J,KS)*CI(IC+1,JC,LL)
            CS=-SO(I,J-1,KO)*CI(IC,JC,LA)+SO(I,J-1,KW)
     &           *CI(IC-1,JC-1,LNE)+SO(I,J,KS)+SO(I+1,J-1,KW)
     &           *CI(IC,JC-1,LNW)
            CSW=-SO(I-1,J-1,KO)*CI(IC-1,JC-1,LNE)+
     &           SO(I-1,J,KS)*CI(IC,JC,LR)+SO(I,J-1,KW)*CI(IC,JC,LA)
            SOC(IC,JC,KO)=-CI(IC-1,JC,LSE)*CNW-CI(IC,JC+1,LB)*CN
     &           -CI(IC,JC,LSW)*CNE-CI(IC,JC,LR)*CW-CO
     &           -CI(IC+1,JC,LL)*CE-CI(IC-1,JC-1,LNE)*CSW
     &           -CI(IC,JC,LA)*CS-CI(IC,JC-1,LNW)*CSE

            SORC(IC,JC,MSOR)=ONE/SOC(IC,JC,KO)

         ENDDO
      ENDDO
      
C------------------------------------------------------------------------------
C     Bug Fix 3(bb): Copying of SOC(ic,jc,KNW) modified
C------------------------------------------------------------------------------

  226 IF(IABS(IBC).NE.1.AND.IABS(IBC).NE.3)GO TO 222
      DO 221 IC=1,IIC
      SOC(IC,JJC,KS)=SOC(IC,2,KS)
      SOC(IC,JJC,KW)=SOC(IC,2,KW)
      SOC(IC,JJC,KNW)=SOC(IC,2,KNW)
      SOC(IC,JJC,KSW)=SOC(IC,2,KSW)
      SOC(IC,JJC,KO)=SOC(IC,2,KO)
      SOC(IC,1,KO)=SOC(IC,JJC1,KO)
      SOC(IC,1,KW)=SOC(IC,JJC1,KW)
      SOC(IC,1,KS)=SOC(IC,JJC1,KS)
      SOC(IC,1,KSW)=SOC(IC,JJC1,KSW)
      SOC(IC,1,KNW)=SOC(IC,JJC1,KNW)
  221 CONTINUE
  222 IF(IABS(IBC).NE.2.AND.IABS(IBC).NE.3)GO TO 224
      DO 223 JC=1,JJC
      SOC(IIC,JC,KW)=SOC(2,JC,KW)
      SOC(IIC,JC,KS)=SOC(2,JC,KS)
      SOC(IIC,JC,KNW)=SOC(2,JC,KNW)
      SOC(IIC,JC,KSW)=SOC(2,JC,KSW)
      SOC(IIC,JC,KO)=SOC(2,JC,KO)
      SOC(1,JC,KW)=SOC(IIC1,JC,KW)
      SOC(1,JC,KS)=SOC(IIC1,JC,KS)
      SOC(1,JC,KSW)=SOC(IIC1,JC,KSW)
      SOC(1,JC,KNW)=SOC(IIC1,JC,KNW)
      SOC(1,JC,KO)=SOC(IIC1,JC,KO)
  223 CONTINUE
  224 CONTINUE
C   end of computation of grid kc difference operator, when kf
C   difference operator is five point
C******************************

 230  IF (ISTRT.LT.0) GO TO 250 ! IF FMG cycles then restrict QF

C
C   unless istrt.lt.0 form right hand side for grid kc
C   by weighting right hand side of grid kf with i(transpose).
C

C     Ensure periodicity prior to restriction 
C     (consistent with mgrcap.f)

      IF ( IABS(IBC).EQ.2 .OR. IABS(IBC).EQ.3 ) THEN
         DO J=2,JJF1
            QF(1,J)  = QF(IIF1,J)
            QF(IIF,J)= QF(2,J)
         ENDDO
      ENDIF
         
      IF ( IABS(IBC).EQ.1 .OR. IABS(IBC).EQ.3 ) THEN
         DO I=2,IIF1
            QF(I,1)  = QF(I,JJF1)
            QF(I,JJF)= QF(I,2)
         ENDDO
      ENDIF

      J=0
      DO 241 JC=2,JJC1
         !
         J=J+2
         I=0
         !
         DO 240 IC=2,IIC1
            !
            I=I+2
            !
            QFC(IC,JC) = CI(IC-1,JC-1,LNE)*QF(I-1,J-1)
     &                 + CI(IC,JC,LA)*QF(I,J-1)
     &                 + CI(IC,JC-1,LNW)*QF(I+1,J-1)
     &                 + CI(IC,JC,LR)*QF(I-1,J)
     &                 + QF(I,J)
     &                 + CI(IC+1,JC,LL)*QF(I+1,J)
     &                 + CI(IC-1,JC,LSE)*QF(I-1,J+1)
     &                 + CI(IC,JC+1,LB)*QF(I,J+1)
     &                 + CI(IC,JC,LSW)*QF(I+1,J+1)
            !
            SORC(IC,JC,MTOT) = QFC(IC,JC)
            !
 240     CONTINUE
 241  CONTINUE

      IF( ISKIP.EQ.2 ) RETURN   ! Skip setup of relaxation 
C
C   if irelax=1, point relaxation. hence, return.
C   if irelax=2, relaxation by lines in x. form lu decomposition
C   of tridiagonal matrices along y=const. lines. if irelax=3,
C   relaxation by lines in y. form lu decomposition of tridiagonal
C   matrices along x=const. lines. if irelax=4, relaxation by lines in x
C   and y. hence need lu decomposition of both sets of tridiagonal
C   matrices.
C
 250  GO TO (260,270,290,270), IRELAX

 260  RETURN                    ! point relaxation so return

      ! Factorization for lines in x
 270  DO 281 J=2,JJF1
         SOR(2,J,MSOR)=ONE/SO(2,J,KO)
         SOR(2,J,MTOT)=QF(2,J)*SOR(2,J,MSOR)
         DO 280 I=3,IIF1
            SOR(I,J,MSOR)=ONE
     &           /(SO(I,J,KO)-SOR(I-1,J,MSOR)*SO(I,J,KW)**2)
            SOR(I,J,MTOT)=QF(I,J)/SO(I,J,KO)
 280     CONTINUE
 281  CONTINUE

      IF(IABS(IBC).NE.2.AND.IABS(IBC).NE.3)GO TO 296  ! not periodic in x jump to y

      ! Periodic fixup of factorization (use ghost part of QF for workspace)
      DO 294 J=2,JJF1
         QF(1,J)=SO(IIF,J,KW)
         QF(IIF,J)=SO(2,J,KW)*SOR(2,J,MSOR)
         SOR(IIF1,J,MSOR)=SO(IIF1,J,KO)-QF(1,J)*QF(IIF,J)
         DO 293 I=3,IIF2
            QF(IIF,J)=QF(IIF,J)*SO(I,J,KW)*SOR(I,J,MSOR)
            QF(1,J)=QF(1,J)*SO(I,J,KW)*SOR(I-1,J,MSOR)
            SOR(IIF1,J,MSOR)=SOR(IIF1,J,MSOR)-QF(IIF,J)*QF(1,J)
 293     CONTINUE
         SOR(IIF1,J,MSOR)=ONE/( SOR(IIF1,J,MSOR) 
     &        -(SO(IIF1,J,KW)+QF(1,J))
     &        *SOR(IIF2,J,MSOR)*SO(IIF1,J,KW)
     &        -SO(IIF1,J,KW)*QF(IIF,J) )
 294  CONTINUE
      ! Restore periodicity in QF workspace
      DO 295 J=2,JJF1
         QF(1,J)  = QF(IIF1,J)
         QF(IIF,J)= QF(2,J)
 295  CONTINUE

 296  IF (IRELAX.NE.4) RETURN  ! just line in x so return

      ! Factorization for lines in y
 290  DO 301 I=2,IIF1
         SOR(I,2,MSOS)=ONE/SO(I,2,KO)
         SOR(I,2,MTOT)=QF(I,2)*SOR(I,2,MSOS)
         DO 300 J=3,JJF1
            SOR(I,J,MSOS)=ONE
     &           /(SO(I,J,KO)-SOR(I,J-1,MSOS)*SO(I,J,KS)**2)
            SOR(I,J,MTOT)=QF(I,J)/SO(I,J,KO)
 300     CONTINUE
 301  CONTINUE

      IF (IABS(IBC).NE.1 .AND. IABS(IBC).NE.3) RETURN ! not periodic in y

      ! Periodic fixup of factorization (use ghost part of QF for workspace)
      DO 310 I=2,IIF1
         QF(I,1)=SO(I,JJF,KS)
         QF(I,JJF)=SO(I,2,KS)*SOR(I,2,MSOS)
         SOR(I,JJF1,MSOS)=SO(I,JJF1,KO)-QF(I,1)*QF(I,JJF)
         DO 305 J=3,JJF2
            QF(I,JJF)=QF(I,JJF)*SO(I,J,KS)*SOR(I,J,MSOS)
            QF(I,1)=QF(I,1)*SO(I,J,KS)*SOR(I,J-1,MSOS)
            SOR(I,JJF1,MSOS)=SOR(I,JJF1,MSOS)-QF(I,JJF)*QF(I,1)
 305     CONTINUE
         SOR(I,JJF1,MSOS)=ONE/( SOR(I,JJF1,MSOS)
     &        -(SO(I,JJF1,KS)+QF(I,1))
     &        *SOR(I,JJF2,MSOS)*SO(I,JJF1,KS)
     &        -SO(I,JJF1,KS)*QF(I,JJF) )
 310  CONTINUE
      ! Restore periodicity in QF workspace
      DO 315 I=2,IIF1
         QF(I,1)  = QF(I,JJF1)
         QF(I,JJF)= QF(I,2)
 315  CONTINUE
C     
      RETURN
      END
