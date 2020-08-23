      SUBROUTINE CALC_O2_KELMAN(content,Hb,pH,PCO2,PO2,SHbO2,Temp,
     &  ERROR,*)

C#### Subroutine: CALC_O2_KELMAN
C###  Description:
C###     Calculates the oxyhemoglobin saturation and oxygen content 
C###     from the partial pressure. Based on Kelman 1966 equations.
C###     Not accurate at low PO2 values (PO2 < 25 mmHg).
C###  Inputs:
C###     Hb -> Hemoglobin content molar (standard value = 2.33E-3 = 0.15g/ml)
C###     pH -> pH of plasma (standard value = 7.39)
C###     PCO2 -> CO2 partial pressure mmHg (standard value = 40)
C###     PO2 -> O2 partial pressure mmHg (standard value = 100)
C###     Temp -> temperature degrees celcius (standard value = 37)
C###  Outputs:
C###     content -> O2 content ml O2 / ml blood
C###     SHbO2 -> fractional oxyhemoglobin saturation
C**** Created by AJS, Nov 2009

     
      IMPLICIT NONE

      REAL*8 content,Hb,pH,PCO2,PO2,SHbO2,Temp
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 A1,A2,A3,A4,A5,A6,A7,alphaO2,O2molarVol,Wbl,X

      CALL ENTERS('CALC_O2_KELMAN',*9999)

C     Parameters
      alphaO2=1.46d-6 ! O2 solubilitiy in water at T=37, molar/mmHg
      O2molarVol=22.256d0  ! ml/mol
      Wbl=0.81d0 !fractional water content of blood

!       write(*,*) 'PO2=',PO2,' content=',content
      IF(PO2.EQ.0.0d0) THEN
        SHbO2=0.0d0
        content=0.0d0
      ELSE

C     Calculate Hb-O2 saturation
      A1=-8.538889d3
      A2=2.121401d3
      A3=-6.707399d1
      A4=9.359609d5
      A5=-3.134626d4
      A6=2.396167d3
      A7=-6.710441d1
      X=PO2*10.0d0**(0.024d0*(37.0d0-Temp)+0.4d0*(pH-7.4d0)+
     &  0.06d0*(DLOG10(DBLE(40.0d0))-DLOG10(DBLE(PCO2))))
      SHbO2=(X*(X*(X*(X+A3)+A2)+A1))/(X*(X*(X*(X+A7)+A6)+A5)+A4)

C     Calculate O2 content (convert from molar to ml O2 per ml blood)
      content=(Wbl*alphaO2*PO2+4.0d0*Hb*SHbO2)*O2molarVol

      ENDIF

      IF(content.LT.0.0d0) content=0.0d0 !curve fit behaves poorly at low PO2
      IF(SHbO2.LT.0.0d0) SHbO2=0.0d0

      CALL EXITS('CALC_O2_KELMAN')
      RETURN
 9999 CALL ERRORS('CALC_O2_KELMAN',ERROR)
      CALL EXITS('CALC_O2_KELMAN')
      RETURN 1
      END

