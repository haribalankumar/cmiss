      SUBROUTINE AZTG45(nw,AXL,AZL,CG,EG,TG,ERROR,*)

C#### Subroutine: AZTG45
C###  Description:
C###    <HTML> <PRE>
C###    AZTG45 is for  membrane problems only. Stresses are in kN/m.
C###    For isotropic materials (IMT(ie)=1):
C###      CG(1) is Young's modulus in GPa;
C###      CG(2) is Poisson's ratio;
C###      CG(3) is membrane thickness in mm;
C###      YMH   is Young's modulus * membrane thickness in kN/m.
C###    For orthotropic materials (IMT(ie)=4):
C###      CG(1) is Young's modulus in fibre direction in GPa;
C###      CG(2) is Young's modulus in transverse direction in GPa;
C###      CG(3) is in-plane shear modulus in GPa;
C###      CG(4) is major Poisson's ratio Nu12;
C###      CG(5) is membrane thickness in mm.
C###      EH1,EH2,GH12 are modulus * membrane thickness in kN/m.
C###    </PRE> </HTML>

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER nw
      REAL*8 AXL(3,*),AZL(3,*),CG(NMM),EG(3,*),TG(3,*)
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 DENOM,EH1,EH2,EKK,EVAL(2),GH12,NU12,NU21,RLA,RMU,
     '  RPR,YMH

      CALL ENTERS('AZTG45',*9999)
      EG(1,1)=0.5D0*(AZL(1,1)-AXL(1,1))
      EG(2,2)=0.5D0*(AZL(2,2)-AXL(2,2))
      EG(2,1)=0.5D0*(AZL(2,1)-AXL(2,1))
      EG(1,2)=EG(2,1)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' EG: '',3E11.3)') EG(1,1),EG(2,2),EG(1,2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      IF(IMT(nw).EQ.1) THEN
C ***   Isotropic materials
        YMH=CG(1)*CG(3)*1.D3
        RPR=CG(2)
        RLA=YMH*RPR/(1.D0+RPR)/(1.D0-2.D0*RPR)
        RMU=YMH/2.D0/(1.D0+RPR)
        EG(3,3)=-RPR/(1.D0-RPR)*(EG(1,1)+EG(2,2))
        EKK=EG(1,1)+EG(2,2)+EG(3,3)
        TG(1,1)=RLA*EKK+2.D0*RMU*EG(1,1)
        TG(2,2)=RLA*EKK+2.D0*RMU*EG(2,2)
        TG(2,1)=        2.D0*RMU*EG(2,1)
        TG(1,2)=TG(2,1)
        IF(TG(1,1).LT.0.D0) TG(1,1)=0.D0
        IF(TG(2,2).LT.0.D0) TG(2,2)=0.D0

      ELSE IF(IMT(nw).EQ.2) THEN

      ELSE IF(IMT(nw).EQ.3) THEN

      ELSE IF(IMT(nw).EQ.4) THEN
C ***   Orthotropic materials. Use 1D relation if princ strain is LE 0.
        CALL EVALUE(2,EG,EVAL,ERROR,*9999)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' E.values:'',2E12.4)') EVAL(1),EVAL(2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
c        IF(EVAL(1).GT.0.D0.AND.EVAL(2).GT.0.D0) THEN
          EH1 =CG(1)*CG(5)*1.D3
          EH2 =CG(2)*CG(5)*1.D3
          GH12=CG(3)*CG(5)*1.D3
          NU12=CG(4)
          NU21=NU12*EH2/EH1
          DENOM=1.D0-NU12*NU21
          TG(1,1)=(     EH1*EG(1,1)+NU21*EH1*EG(2,2))/DENOM
          TG(2,2)=(NU12*EH2*EG(1,1)+     EH2*EG(2,2))/DENOM
          TG(2,1)= GH12*EG(2,1)
c        ELSE IF(EVAL(1).GT.0.D0.AND.EVAL(2).LE.0.D0) THEN
c          CALL EVECTR(2,EG,EVAL(1),EVEC,ERROR,*9999)
c          S11=EVEC(1)
c          S12=EVEC(2)
c          EH1=CG(1)*CG(5)*1.E3
c          PRINC_STRESS_1=EH1*EVAL(1)
c          TG(1,1)=PRINC_STRESS_1*S11*S11
c          TG(2,1)=PRINC_STRESS_1*S11*S12
c          TG(2,2)=PRINC_STRESS_1*S12*S12
c        ELSE IF(EVAL(2).GT.0.D0.AND.EVAL(1).LE.0.D0) THEN
c          CALL EVECTR(2,EG,EVAL(2),EVEC,ERROR,*9999)
c          S21=EVEC(1)
c          S22=EVEC(2)
c          EH2=CG(2)*CG(5)*1.E3
c          PRINC_STRESS_2=EH2*EVAL(2)
c          TG(1,1)=PRINC_STRESS_2*S21*S21
c          TG(2,1)=PRINC_STRESS_2*S21*S22
c          TG(2,2)=PRINC_STRESS_2*S22*S22
c        ENDIF
        TG(1,2)=TG(2,1)
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' TG: '',3E11.3)') TG(1,1),TG(2,2),TG(1,2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('AZTG45')
      RETURN
 9999 CALL ERRORS('AZTG45',ERROR)
      CALL EXITS('AZTG45')
      RETURN 1
      END


