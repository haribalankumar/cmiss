      SUBROUTINE CHARAC(AP,AQ,AR,BP,BQ,BR,CP,CQ,CR,FP,FQ,FR,
     '  XP,XQ,XR,TP,TQ,TR,PP,PQ,PR,QP,QQ,QR,UP,UQ,UR)

C#### Subroutine: CHARAC
C###  Description:
C###    CHARAC calculates solution at point R from known values
C###    at points P and  Q on previous solution front.

C**** Calls Functions ALFA and BETA to calc characteristic gradients
C**** through points P & Q.
C**** Calls Functions A_COEFF,B_COEFF,C_COEFF,F_COEFF to calculate PDE
C**** coefficients A,B,C,F at each solution point.

      IMPLICIT NONE
!     Parameter List
      REAL*8 AP,AQ,AR,BP,BQ,BR,CP,CQ,CR,FP,FQ,FR,PP,PQ,PR,QP,QQ,QR,
     '  TP,TQ,TR,UP,UQ,UR,XP,XQ,XR
!     Local Variables
      INTEGER NIT
      REAL*8 A_COEFF,AALFAC_RP,ABETAC_RQ,ALFA,ALFA_P,ALFA_R,ALFA_RP,
     '  B_COEFF,BETA,BETA_Q,BETA_R,BETA_RQ,BOT,C_COEFF,F_COEFF,FTC_RP,
     '  FTC_RQ,P_RP,P_RQ,Q_RP,Q_RQ,TOP
      LOGICAL left_to_right

      left_to_right=.TRUE.
      NIT=1

 10   ALFA_P=ALFA(AP,BP,CP)

      IF(NIT.EQ.1) THEN
        ALFA_R=ALFA_P
      ELSE
        ALFA_R=ALFA(AR,BR,CR)
      ENDIF

      BETA_Q=BETA(AQ,BQ,CQ)

      IF(NIT.EQ.1) THEN
        BETA_R=BETA_Q
      ELSE
        BETA_R=BETA(AR,BR,CR)
      ENDIF

      ALFA_RP=0.5D0*(ALFA_R+ALFA_P)
      BETA_RQ=0.5D0*(BETA_R+BETA_Q)


      XR=(TQ-TP+ALFA_RP*XP-BETA_RQ*XQ)/(ALFA_RP-BETA_RQ)
      TR=TP+ALFA_RP*(XR-XP)


C *** Evaluate A,B,C,F at R for iteration if they are f(X,T).

C *** AR=
C *** BR=
C *** CR=
C *** FR=

      IF(NIT.EQ.1) THEN
        AR=AQ
        CR=CQ
        FR=FQ
      ENDIF

      FTC_RQ=(FR+FQ)*(TR-TQ)/(CR+CQ)
      ABETAC_RQ=(AR*BETA_R+AQ*BETA_Q)/(CR+CQ)

      IF(NIT.EQ.1) THEN
        AR=AP
        CR=CP
        FR=FP
      ENDIF

      FTC_RP=(FR+FP)*(TR-TP)/(CR+CP)
      AALFAC_RP=(AR*ALFA_R+AP*ALFA_P)/(CR+CP)

C *** ABETAC_RQ=0.5D0*(AR+AQ)*(BETA_R+BETA_Q)/(CR+CQ)
C *** AALFAC_RP=0.5D0*(AR+AP)*(ALFA_R+ALFA_P)/(CR+CP)

      TOP=(QQ-QP+FTC_RQ-FTC_RP+PQ*ABETAC_RQ-PP*AALFAC_RP)
      BOT=ABETAC_RQ-AALFAC_RP
      PR=TOP/BOT
      QR=QP+FTC_RP-AALFAC_RP*(PR-PP)

C *** Evaluate A,B,C,F at R for iteration if they are f(X,T,P,Q)

C *** AR=
C *** BR=
C *** CR=
C *** FR=

      IF(left_to_right) THEN
        P_RP=0.5D0*(PR+PP)
        Q_RP=0.5D0*(QR+QP)
        UR=UP+P_RP*(XR-XP)+Q_RP*(TR-TP)
      ELSE
        P_RQ=0.5D0*(PR+PQ)
        Q_RQ=0.5D0*(QR+QQ)
        UR=UQ+P_RQ*(XR-XQ)+Q_RQ*(TR-TQ)
      ENDIF

C *** Evaluate A,B,C,F at R for iteration if they are f(X,T,P,Q,U)

      AR=A_COEFF()
      BR=B_COEFF(XR)
      CR=C_COEFF(XR)
      FR=F_COEFF()

      NIT=NIT+1

      IF(NIT.LT.3) GO TO 10

C *** WRITE (*,*)'ALFA_RP=', ALFA_RP
C *** WRITE (*,*)'BETA_RQ=', BETA_RQ

      RETURN
      END


