      SUBROUTINE ENERGY(nr,CG,DW,P1,P2,P3,P4,P5,P6,YG,ERROR,*)

C#### Subroutine: ENERGY
C###  Description:
C###    ENERGY calculates derivatives of strain energy function wrt
C###    either principal strain invariants (KTYP55(nr)=1), or
C###    principal extensions (KTYP55(nr)=2), or physical strains
C###    (KTYP55(nr)=3) at current Gauss point.
C**** For KTYP55(nr)=1:
C****      P1 is the First  principal invariant RI1;
C****      P2 is the Second principal invariant RI2;
C****      P3 is the Third  principal invariant RI3;
C****      P4 is the First  transverse isotropic invariant RK1;
C****      P5 is the Second transverse isotropic invariant RK2;
C****      DW(1..5) are dW/dI1,dW/dI2,dW/dI3,dW/dK1,dW/dK2.
C**** For KTYP55(nr)=2:
C****      P1 is the First  principal extension ratio RL1;
C****      P2 is the Second principal extension ratio RL2;
C****      P3 is the Third  principal extension ratio RL3;
C****      DW(1..3) are dW/dL1,dW/dL2,dW/dL3.
C**** For KTYP55(nr)=3:
C****      P1..P6 are the physical components of Green's strain wrt
C****      theta (KTYP53(nr)=1) or nu coords (KTYP53(nr)>1);
C****      E(1,1),E(2,2),E(3,3),E(1,2),E(1,3) and E(2,3);
C****      DW(1..6) are dW/dE(1,1) .. dW/dE(2,3).

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'lung00.cmn'
!     Parameter List
      INTEGER nr
      REAL*8 CG(NMM),DW(6),P1,P2,P3,P4,P5,P6,YG(NIYGM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,i1,i2,i3,il,iterm
      REAL*8 alpha1,alpha2,beta,coeff1,coeff2,CEXPQ,DENOM1,DENOM2,
     '  DF1,DF2,DF3,DtoAlpha1,DtoAlpha2,DWslope1,DWslope2,
     '  exp_term,FI1,FI2,FI3,pole1,pole2,
     '  R1,R2,R3,RI1,RI2,RI3,RK1,RL1,RL2,RL3,
     '  strain,strain0,strain1,strain2,TOL,ZERO,mytol
      PARAMETER(ZERO=0.0d0)

      CALL ENTERS('ENERGY',*9999)


C news VJ 10Dec2003
      IF(KTYP54(nr).EQ.3) THEN !gauss point stresses (grid based)
C Refer to Holger Schmid.et.al technical note on definition of stress in terms of 
C derivatives of the strain energy function. "Ambiguities in Hyperelastic Constitutive Law Formulation
C The formula for the stress components of a material depend on the formulation of the strain energy
C density function. If the strain energy function may be expressed in terms of 6 or 9 different strain
C components - W(E11,E22,E33,E12,E13,E21,E23) or W(E11,E22,E33,E12,E13,E23,E21,E32,E31). If expressed in terms
C of 6 components then 
C TGij deviatoric = dWdEij for i=j i=1...3, j=1...3
C TGij deviatoric = 0.5dWdEij for i<>j i=1...3, j=1...3
C else if expressed in terms of 9 components
C TGij deviatoric = dWdEij for i,j..i=1...3, j=1...3
C Hence it has been decided, in order for people to have the freedom to express W in their own way, DW array
C will contain the deviatoric component of stress calculated in the CellML file and in the YG array. 
        DW(1) = YG(1)
        DW(2) = YG(2)
        DW(3) = YG(3)
        DW(4) = YG(4)
        DW(5) = YG(5)
        DW(6) = YG(6)
C newe VJ 10Dec2003
      ELSE IF(KTYP54(nr).NE.3) THEN !Other forms of constitutive laws
        IF(KTYP55(nr).EQ.1) THEN      !strain invariants
          RI1=P1
          RI2=P2
          RI3=P3
          RK1=P4
          R1=RI1-3.0d0
          R2=RI2-3.0d0
          R3=RI3-1.0d0
        ELSE IF(KTYP55(nr).EQ.2) THEN !principal extensions
          RL1=P1
          RL2=P2
          RL3=P3
          R1=RL1-1.0d0
          R2=RL2-1.0d0
          R3=RL3-1.0d0
        ELSE IF(KTYP55(nr).EQ.3) THEN !fibre and transverse strains
          R1=P1
          R2=P2
          R3=P3
        ENDIF
  
        IF(KTYP56(nr).EQ.1) THEN !W is polynomial
          il=0
          DW(1)=0.0d0
          DW(2)=0.0d0
          DW(3)=0.0d0
          DW(4)=0.0d0
          DW(5)=0.0d0
          DO i3=0,IT(3,nr)
            FI3=1.0d0
            DF3=1.0d0
            IF(i3.GT.0) FI3=R3**i3
            IF(i3.GT.1) DF3=DBLE(i3)*R3**(i3-1)
            IF(KTYP52(nr).EQ.4) THEN !compressible + fluid
              DF3=DF3*(1.d0-CG(IL_porosity))
            ENDIF
            DO i2=0,IT(2,nr)
              FI2=1.0d0
              DF2=1.0d0
              IF(i2.GT.0) FI2=R2**i2
              IF(i2.GT.1) DF2=DBLE(i2)*R2**(i2-1)
              DO i1=0,IT(1,nr)
                FI1=1.0d0
                DF1=1.0d0
                IF(i1.GT.0) FI1=R1**i1
                IF(i1.GT.1) DF1=DBLE(i1)*R1**(i1-1)
                il=il+1
                IF(i1.GT.0) DW(1)=DW(1)+CG(il)*DF1*FI2*FI3
                IF(i2.GT.0) DW(2)=DW(2)+CG(il)*FI1*DF2*FI3
                IF(i3.GT.0) DW(3)=DW(3)+CG(il)*FI1*FI2*DF3
              ENDDO !i1
            ENDDO !i2
          ENDDO !i3


        ELSE IF(KTYP56(nr).EQ.2) THEN !W is special function\
          DO i=1,6
            DW(i)=0.d0
          ENDDO !i
          IF(KTYP55(nr).EQ.1) THEN      !W is a func of I1,I2,I3
c            IF(KTYP52(nr).EQ.1) THEN !W is Blatz-Ko material
c              DW(1)=CG(1)/RI3
c              DW(2)=0.0d0
c              DW(3)=-CG(1)*(RI1+2.0d0)/(RI3*RI3)
  
!   Alex 19-Nov-02: replaced old Blatz-Ko function
C   W=0.5*u*a*((J1-3)+B*(J3**(-2/B)-1))=0.5*(1-a)*((J2-3)+B*(J3**(2/B)
C     -1))

             IF(KTYP52(nr).EQ.1) THEN      !W is Blatz-Ko material
 
!   Alex 21-Nov-02: calculate beta from nu

               beta=(1.0d0-2.0d0*CG(3))/CG(3)
               DW(1)=CG(1)*CG(2)/2.0d0
               DW(2)=CG(1)*(1.0d0-CG(2))/(2.0d0*RI3)
               DW(3)=-0.5d0*CG(1)*(1.0d0-CG(2))*RI2/(RI3*RI3)
     '           -0.5d0*CG(1)*CG(2)*(RI3**(-(1.0d0+beta)/beta))
     '           +0.5d0*CG(1)*(1.0d0-CG(2))*(RI3**((1.0d0-beta)/beta))
            ELSE                      !W is Exponential in I1,I2 (Demiray)
              DW(1)=CG(1)*CG(2)*DEXP(CG(2)*(RI1-3.0d0))
              DW(2)=0.0d0
            ENDIF
          ELSE IF(KTYP55(nr).EQ.2) THEN !func of princ extn ratios (Ogden)
            IF(KTYP52(nr).EQ.1) THEN
              DW(1)=CG(1)*P1**(CG(2)-1)
              DW(2)=CG(1)*P2**(CG(2)-1)
              DW(3)=CG(1)*P3**(CG(2)-1)
            ELSE
              DW(1)=CG(1)*P1**(CG(2)-1)
              DW(2)=CG(1)*P2**(CG(2)-1)
            ENDIF
          ELSE IF(KTYP55(nr).EQ.3) THEN !W=func of fibre and trans strains
            TOL=1.0D-08         !trans isotrop expon law W=Cexp(Q) (Fung)
            IF((DABS(P1).LT.TOL).AND.(DABS(P2).LT.TOL).AND.
     '        (DABS(P3).LT.TOL).AND.(DABS(P4).LT.TOL).AND.
     '        (DABS(P5).LT.TOL).AND.(DABS(P6).LT.TOL)) THEN
              IF(KTYP55a(nr).EQ.3) THEN ! Tong & Fung
                CEXPQ=0.5d0*CG(5)
              ELSE ! Fung or Holmes
                CEXPQ=CG(1)
              ENDIF
            ELSE
C CS 18/7/2001 Allowing multiple types of this law
C and adding Jeffery Holmes's distributed fibre formulation
              IF(KTYP55a(nr).EQ.1) THEN
C               Transversely isotropic law from
C               Guccione 1991 J. Biomech. Eng. 113: 42-55
                CEXPQ=CG(1)*DEXP(2.0d0*CG(2)*(P1+P2+P3)
     '            +     CG(3)* P1*P1
     '            +     CG(4)*(P2*P2+P3*P3+2.0d0*P6*P6)
     '            + 2.0d0*CG(5)*(P4*P4+P5*P5))
C MPN 4-Feb-1994: replaced Guccione exponential law for D. Bloomgarden
c           CEXPQ1=CG(1)*DEXP(CG(2)*P1*P1+CG(3)*(P4*P4+P5*P5))
c           CEXPQ2=CG(4)*DEXP(CG(5)*(P2+P3)**2+CG(6)*(P2*P3-P6*P6))
 
              ELSE IF(KTYP55a(nr).EQ.2) THEN
C               Jeffery Holme's distributed fibre formulation
                CEXPQ=CG(1)*DEXP(CG(2)*2.0d0*(P1+P2+P3)
     '            +    CG(3)*(CG(4)*P1*P1
     '            +    CG(5)*P2*P2
     '            +    CG(6)*(2.0d0*P1*P2+4.0d0*P4)
     '            +    CG(7)*4.0d0*P1*P4
     '            +    CG(8)*4.0d0*P4*P2))
              ELSE IF(KTYP55a(nr).EQ.3) THEN
C               Tong & Fung skin equation
                CEXPQ=0.5d0*(CG(1)*P1*P1+CG(2)*P2*P2+CG(3)*P4*P4+
     '            2*CG(4)*P1*P2)
     '            + 0.5d0*CG(5)*DEXP(CG(6)*P1*P1+CG(7)*P2*P2
     '            + CG(8)*P4*P4+2.0d0*CG(9)*P1*P2
     '            + CG(10)*P1*P1*P1
     '            + CG(11)*P2*P2*P2
     '            + CG(12)*P1*P1*P2
     '            + CG(13)*P1*P2*P2)
              ELSE
C               Orthotropic law from Costa, 2001, Philosophical Transactions
C               of the Royal Society, 359(1783): 1233-1250.
                CEXPQ=CG(1)*DEXP( CG(2)*P1*P1 + CG(3)*P2*P2
     '             + CG(4)*P3*P3 + 2.0d0*CG(5)*P4*P4
     '             + 2.0d0*CG(6)*P5*P5
     '             + 2.0d0*CG(7)*P6*P6 )
              ENDIF
            ENDIF

            IF(KTYP55a(nr).EQ.1) THEN
C             Transversely isotropic law from
C             Guccione 1991 J. Biomech. Eng. 113: 42-55
              DW(1)= CEXPQ*(CG(2)+CG(3)*P1)
              DW(2)= CEXPQ*(CG(2)+CG(4)*P2)
              DW(3)= CEXPQ*(CG(2)+CG(4)*P3)
C KFA 2003-09-26: recorrected errors in derivatives of W wrt shear strains
C MPN 4Apr2003: corrected errors in derivatives of W  wrt shear strains
C              DW(4)= CEXPQ*2.0d0*CG(5)*P4
C              DW(5)= CEXPQ*2.0d0*CG(5)*P5
C              DW(6)= CEXPQ*2.0d0*CG(4)*P6
              DW(4)= CEXPQ*CG(5)*P4
              DW(5)= CEXPQ*CG(5)*P5
              DW(6)= CEXPQ*CG(4)*P6
C MPN 4-Feb-1994: replaced Guccione exponential law for D. Bloomgarden
c            DW(1)= CEXPQ1*2.0d0*CG(2)*P1
c            DW(2)= CEXPQ2*(2.0d0*CG(5)*(P2+P3)+CG(6)*P3)
c            DW(3)= CEXPQ2*(2.0d0*CG(5)*(P2+P3)+CG(6)*P2)
c            DW(4)= CEXPQ1*CG(3)*2.0d0*P4
c            DW(5)= CEXPQ1*CG(3)*2.0d0*P5
c            DW(6)=-CEXPQ2*CG(6)*2.0d0*P6

            ELSE IF(KTYP55a(nr).EQ.2) THEN
C             Jeffery Holme's distributed fibre formulation
              DW(1)= CEXPQ*(CG(2)+CG(3)*(CG(4)*P1+
     '                 CG(6)*P2+2.0d0*CG(7)*P4))
              DW(2)= CEXPQ*(CG(2)+CG(3)*(CG(5)*P2+
     '                 CG(6)*P1+2.0d0*CG(8)*P4))
              DW(3)= CEXPQ*CG(2)
C MPN 4Apr2003: fixed bug in deriv wrt E12 (P4)
              DW(4)= CEXPQ*CG(3)*
     '                 (2.0d0*CG(6)*P4+CG(7)*2.0d0*P1+CG(8)*2.0d0*P2)
C              DW(4)= CEXPQ*CG(3)*
C     '                 (2.0d0*CG(6)*P4+CG(7)*P1+CG(8)*P2)
              DW(5)= 0.0d0
              DW(6)= 0.0d0
C DAH 12-Nov-2002: Adding in equations for skin.
            ELSE IF(KTYP55a(nr).EQ.3) THEN
C             Tong & Fung skin equation
              DW(1)=CG(1)*P1+CG(4)*P2+CG(5)*
     '          (CG(6)*P1+CG(9)*P2+(3.0d0/2.0d0)*CG(10)*P1*P1+
     '          CG(12)*P1*P2+0.5d0*CG(13)*P2*P2)*
     '          DEXP(CG(6)*P1*P1+CG(7)*P2*P2
     '            + CG(8)*P4*P4+2.0d0*CG(9)*P1*P2
     '            + CG(10)*P1*P1*P1
     '            + CG(11)*P2*P2*P2
     '            + CG(12)*P1*P1*P2
     '          + CG(13)*P1*P2*P2)
              DW(2)=CG(4)*P1+CG(2)*P2+CG(5)*
     '          (CG(9)*P1+CG(7)*P2+(3.0d0/2.0d0)*CG(11)*P1*P1+
     '          CG(13)*P1*P2+0.5d0*CG(12)*P1*P1)*
     '          DEXP(CG(6)*P1*P1+CG(7)*P2*P2
     '            + CG(8)*P4*P4+2.0d0*CG(9)*P1*P2
     '            + CG(10)*P1*P1*P1
     '            + CG(11)*P2*P2*P2
     '            + CG(12)*P1*P1*P2
     '          + CG(13)*P1*P2*P2)
              DW(3)=0.0d0
              DW(4)=CG(3)*P4+CG(5)*CG(8)*P4*
     '          DEXP(CG(6)*P1*P1+CG(7)*P2*P2
     '            + CG(8)*P4*P4+2.0d0*CG(9)*P1*P2
     '            + CG(10)*P1*P1*P1
     '            + CG(11)*P2*P2*P2
     '            + CG(12)*P1*P1*P2
     '          + CG(13)*P1*P2*P2)
              DW(5)=0.0d0
              DW(6)=0.0d0
            ELSE
C             Costa orthotropic
              DW(1)=CEXPQ*CG(2)*P1
              DW(2)=CEXPQ*CG(3)*P2
              DW(3)=CEXPQ*CG(4)*P3
              DW(4)=CEXPQ*CG(5)*P4
              DW(5)=CEXPQ*CG(6)*P5
              DW(6)=CEXPQ*CG(7)*P6
            ENDIF
          ENDIF
  
        ELSE IF(KTYP56(nr).EQ.3) THEN   !W is special function
          IF(KTYP55(nr).EQ.1) THEN        !W is a func of I1,I2,I3
            IF(KTYP52(nr).EQ.1) THEN        !W is
              DW(1)=0.0d0
              DW(2)=0.0d0
              DW(3)=0.0d0
            ELSE IF(KTYP52(nr).EQ.6)THEN !W is exponential from Fung
              exp_term=CG(2)/4.d0*R1**2+CG(3)/4.d0*(R2-2.d0*R1)
              DW(1)=CG(1)/4.d0*DEXP(exp_term)*(CG(2)*R1-CG(3))
              DW(2)=CG(1)*CG(3)/8.d0*DEXP(exp_term)
              DW(3)=0.0
            ELSE                    !W is
              DW(1)=0.0d0
              DW(2)=0.0d0
            ENDIF
          ELSE IF(KTYP55(nr).EQ.2) THEN !func of princ extn ratios (Ogden)
            IF(KTYP52(nr).EQ.1) THEN
              DW(1)=0.0d0
              DW(2)=0.0d0
              DW(3)=0.0d0
            ELSE
              DW(1)=0.0d0
              DW(2)=0.0d0
            ENDIF
          ELSE IF(KTYP55(nr).EQ.3) THEN !func of fibre and transv strains
C           pole-zero law - orthotropic model (incl shear terms)
C 25/2/97 LC archived section : init extns handled by 'growth' defm tensor
C                    see ZGTG53.

C           Calculate the 2nd deriv of the pole-zero law
C           for strains between 0 and 90% of the pole for each term.
C           For strains beyond 90% of their pole the stress/strain is
C           flat with magnitude equal to the stress at a strain of 90%
C           of the pole.
C           For compressive strains use linear stress/strain relationship
C           with stiffness equal to the 2nd deriv of W wrt e
C           evaluated at e=0
            TOL=0.90d0
  
            DO iterm=1,6
              coeff1=CG((iterm-1)*3+1) !coefficient for current term
              pole1 =CG((iterm-1)*3+2) !pole for current term
              alpha1=CG((iterm-1)*3+3) !curvature for current term
              IF(iterm.EQ.1.OR.iterm.EQ.2.OR.iterm.EQ.3) THEN  !Axial
                IF(iterm.EQ.1) THEN
                  strain=P1 !=EG(1,1)
                ELSE IF(iterm.EQ.2) THEN
                  strain=P2 !=EG(2,2)
                ELSE IF(iterm.EQ.3) THEN
                  strain=P3 !=EG(3,3)
                ENDIF
                coeff2=coeff1
                pole2 =pole1
                alpha2=alpha1
              ELSE IF(iterm.EQ.4.OR.iterm.EQ.5.OR.iterm.EQ.6) THEN  !Shear
                IF(iterm.EQ.4) THEN
                  strain=P4 !=EG(1,2)
                ELSE IF(iterm.EQ.5) THEN
                  strain=P5 !=EG(1,3)
                ELSE IF(iterm.EQ.6) THEN
                  strain=P6 !=EG(2,3)
                ENDIF
                coeff2=CG((iterm-1)*3+10)
                pole2 =CG((iterm-1)*3+11)
                alpha2=CG((iterm-1)*3+12)
              ENDIF
C EWR 02Nov03: Making negative shear strain symmetric to the postive 
C                shear strain about the origin. Calculate strain energy
C                as it was positve, then set DW=-DW at the end.
              strain0=strain
              IF(iterm.EQ.4.OR.iterm.EQ.5.OR.iterm.EQ.6) THEN  !Shear
                IF(strain.LE.ZERO) THEN !compressive strain
                  strain=-strain
                ENDIF
              ENDIF
  
              strain1=strain
              strain2=strain
              IF(strain.LE.ZERO) THEN !compressive strain
                strain1=ZERO
                strain2=ZERO
                IF(DOP.AND.strain.LT.ZERO) THEN
C KAT 14May01:   Can't branch out of critical section.
C                Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,
     '              '('' >>Compressive strain in ENERGY'')')
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
              ELSE IF(strain.GT.(TOL*pole1).OR.strain.GT.(TOL*pole2))
     '          THEN
                IF(strain1.GT.(TOL*pole1)) strain1=TOL*pole1 !yielded
                IF(strain2.GT.(TOL*pole2)) strain2=TOL*pole2 !yielded
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C                Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,
     '              '('' Pole strain exceeded in ENERGY'')')
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
              ENDIF
              DENOM1=pole1-strain1
              DtoAlpha1=1.0d0
              IF(alpha1.GT.ZERO) DtoAlpha1=DENOM1**alpha1
              DENOM2=pole2-strain2
              DtoAlpha2=1.0d0
              IF(alpha2.GT.ZERO) DtoAlpha2=DENOM2**alpha2
C             If outside pole-zero range - calc slopes for linear relns
              IF(strain.LE.ZERO.OR.strain.GT.(TOL*pole1)) THEN
C               Calc slope for first set of params
                DWslope1=(2.0d0*coeff1*DENOM1*(DENOM1
     '            +2.0d0*alpha1*strain1)
     '            +coeff1*alpha1*(alpha1+1.0d0)*strain1*strain1)
     '            /(DtoAlpha1*DENOM1*DENOM1)
              ENDIF
              IF(strain.LE.ZERO.OR.strain.GT.(TOL*pole2)) THEN
C               Calc slope for second set of params
                DWslope2=(2.0d0*coeff2*DENOM2*(DENOM2
     '            +2.0d0*alpha2*strain2)
     '            +coeff2*alpha2*(alpha2+1.0d0)*strain2*strain2)
     '            /(DtoAlpha2*DENOM2*DENOM2)
              ENDIF
              IF(strain.GT.ZERO) THEN !Tensile strain
                DW(iterm)=0.5d0*
     '            (coeff1*strain1*
     '            (2.0d0+alpha1*strain1/DENOM1)/DtoAlpha1
     '            +coeff2*strain2*
     '            (2.0d0+alpha2*strain2/DENOM2)/DtoAlpha2)
                IF(strain.GT.(TOL*pole1))
     '            DW(iterm)=DW(iterm)+(strain-strain1)*DWslope1/2.0d0
                IF(strain.GT.(TOL*pole2))
     '            DW(iterm)=DW(iterm)+(strain-strain2)*DWslope2/2.0d0
C EWR 02Nov03: Making negative shear strain symmetric about the origin 
                IF(iterm.EQ.4.OR.iterm.EQ.5.OR.iterm.EQ.6) THEN  !Shear
                  IF(strain0.LE.ZERO) THEN !compressive strain
                    DW(iterm)=-DW(iterm)
                  ENDIF
                ENDIF
              ELSE !Compressive strain
C               Use deriv of DW wrt strain at strain=0 for compress. slope
                DW(iterm)=strain*(DWslope1+DWslope2)/2.0d0
              ENDIF
            ENDDO !iterm
          ENDIF
  
        ELSE IF(KTYP56(nr).EQ.4) THEN   !W is special function
          IF(KTYP55(nr).EQ.1) THEN      !W is a func of I1,I2,I3
            IF(KTYP52(nr).EQ.1) THEN    !W is
              DW(1)=0.0d0
              DW(2)=0.0d0
              DW(3)=0.0d0
            ELSE                    !W is
              DW(1)=0.0d0
              DW(2)=0.0d0
            ENDIF
          ELSE IF(KTYP55(nr).EQ.2) THEN !func of princ extn ratios (Ogden)
            IF(KTYP52(nr).EQ.1) THEN
              DW(1)=0.0d0
              DW(2)=0.0d0
              DW(3)=0.0d0
            ELSE
              DW(1)=0.0d0
              DW(2)=0.0d0
            ENDIF
          ELSE IF(KTYP55(nr).EQ.3) THEN !W=func of fibre and trans strains
            DW(1)=0.0d0
            DW(2)=0.0d0
            DW(3)=0.0d0
            DW(4)=0.0d0
            DW(5)=0.0d0
            DW(6)=0.0d0
          ENDIF

        ELSE IF(KTYP56(nr).EQ.5) THEN   !W is user defined function
          IF(KTYP55(nr).EQ.1) THEN      !in princ strain invariants
            CALL USER51(CG,DW,RK1,ERROR,*9999)
          ELSE IF(KTYP55(nr).EQ.2) THEN !in princ extension ratios
            CALL USER52(ERROR,*9999)
          ELSE IF(KTYP55(nr).EQ.3) THEN !in fibre & transverse strains
            CALL USER53(CG,DW,P1,P2,P3,P4,P5,P6,ERROR,*9999)
          ENDIF

        ELSE IF(KTYP56(nr).EQ.6) THEN   !linear viscous relation
!         This case is dealt with in ZGTG53/4 etc
        ENDIF !ktyp56
      ENDIF !ktyp54
  
      CALL EXITS('ENERGY')
      RETURN
 9999 CALL ERRORS('ENERGY',ERROR)
      CALL EXITS('ENERGY')
      RETURN 1
      END


