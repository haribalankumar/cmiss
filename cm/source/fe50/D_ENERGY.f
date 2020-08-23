      SUBROUTINE D_ENERGY(PARAMTYPE,NMNO,nr,CG,D_DW,P1,P2,P3,P4,P5,P6,
     '  ERROR,*)

C#### Subroutine: D_ENERGY
C###  Description:
C###    D_ENERGY calculates 2nd derivatives of strain energy function
C###    wrt either the NMNO'th material parameter or the
C###    physical strains (KTYP55(nr)=3) at current Gauss point

C**** P1..P6 are the physical components of Green's strain wrt nu
C**** coordinates
C**** (KTYP53(nr)>1):E(1,1),E(2,2),E(3,3),E(1,2),E(1,3) and E(2,3)
C**** D_DW(1..6) are d2W/dE2(1,1) .. d2W/dE2(2,3)

      IMPLICIT NONE
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NMNO,nr
      REAL*8 CG(NMM),D_DW(6),P1,P2,P3,P4,P5,P6
      CHARACTER ERROR*(*),PARAMTYPE*(*)
!     Local Variables
      INTEGER iterm,iterm1,iterm2,pole1num,pole2num
      REAL*8 alpha,alpha1,alpha2,coeff,coeff1,coeff2,
     '  D_DWslope,D_DWslope1,D_DWslope2,DENOM,DENOM1,DENOM2,
     '  Dpole,DtoAlpha,DtoAlpha1,DtoAlpha2,pole,pole1,pole2,
     '  strain,strain1,strain11,strain2,strain21,TOL,ZERO
      PARAMETER(ZERO=0.0d0)

      CALL ENTERS('D_ENERGY',*9999)

      CALL ASSERT(KTYP55(nr).EQ.3,'Strain energy fn must be given as a'
     '  //' function of fibre and transverse strains',ERROR,*9999)
      CALL ASSERT(KTYP56(nr).EQ.3,'Strain energy fn must be'
     '  //' pole-zero in fibre and transverse strains',ERROR,*9999)

C old
C      L0_fibre=CG(28)           !initial fibre ext ratio
C      L0_sheet=CG(29)           !initial sheet ext ratio
C      L0_sheetnormal=CG(30)     !initial sheet-normal ext ratio
C      E0_fibre=0.5d0*(L0_fibre*L0_fibre-1.0d0) !init fibre strain
C      E0_sheet=0.5d0*(L0_sheet*L0_sheet-1.0d0) !init sheet strain
C      E0_sheetnormal=0.5d0*(L0_sheetnormal*L0_sheetnormal-1.0d0) !init sheet-normal strain
C      PP1=P1+E0_fibre
C      PP2=P2+E0_sheet
C      PP3=P3+E0_sheetnormal

      DO iterm=1,6
        D_DW(iterm)=0.0d0
      ENDDO

C      Calculate the 2nd deriv of the pole-zero law
C      for strains between 0 and 90% of the pole for each term.
C      For strains beyond 90% of their pole the stress/strain is flat
C      with magnitude equal to the stress at a strain of 90% of the pole
C      For compressive strains use linear stress/strain relationship
C      with stiffness equal to the 2nd deriv of W wrt e
C      evaluated at e=0
      TOL=0.90d0

      IF(PARAMTYPE.EQ.'MATERIAL_PARAMETERS') THEN

        CALL ASSERT(NMNO.LE.9,' >>Analytic derivs only defined'
     '    //' wrt material param no.s 1-9',ERROR,*9999)

C       Determine term in constitutive law the depends on material param
C       to be optimised
        iterm=(NMNO-1)/3+1    !NOTE: need integer division here
        IF(iterm.EQ.1) THEN
          strain=P1 !=EG(1,1)
        ELSE IF(iterm.EQ.2) THEN
          strain=P2 !=EG(2,2)
        ELSE IF(iterm.EQ.3) THEN
          strain=P3 !=EG(3,3)
        ENDIF

        IF(NMNO.EQ.1.OR.NMNO.EQ.4.OR.NMNO.EQ.7) THEN
C         Parameter is a coefficient.
C !!!     NOTE: if fibre distribution model is altered to express shear
C !!!     coefficients in terms of axial coefficients then this will
C !!!     need updating
          coeff=CG(NMNO)
          pole =CG(NMNO+1)
          alpha=CG(NMNO+2)
          strain1=strain
          IF(strain.LE.ZERO) THEN !compressive strain
            strain1=ZERO
          ELSE IF(strain.GT.(TOL*pole)) THEN
            strain1=TOL*pole !yielded
          ENDIF
          DENOM=pole-strain1
          DtoAlpha=1.0d0
          IF(alpha.GT.ZERO) DtoAlpha=DENOM**alpha
C         If outside pole-zero range - calc slopes for linear relns
          IF(strain.LE.ZERO.OR.strain.GT.(TOL*pole)) THEN
C           Calc slope
            D_DWslope=(2.0d0*DENOM*(DENOM+2.0d0*alpha*strain1)
     '        +alpha*(alpha+1.0d0)*strain1*strain1)
     '        /(DtoAlpha*DENOM*DENOM)
          ENDIF
          IF(strain.GT.ZERO) THEN !Tensile strain
            D_DW(iterm)=strain1*(2.0d0+alpha*strain1/DENOM)/DtoAlpha
            IF(strain.GT.(TOL*pole))
     '        D_DW(iterm)=D_DW(iterm)+(strain-strain1)*D_DWslope
          ELSE !Compressive strain
C           Use deriv of DW at strain=0 for compresive slope
            D_DW(iterm)=strain*D_DWslope
          ENDIF

        ELSE IF(NMNO.EQ.2.OR.NMNO.EQ.5.OR.NMNO.EQ.8) THEN
C         Parameter is a pole - set params for main term in
C         constitutive law
C !!!     NOTE: if fibre distribution model is altered for the shear
C !!!     poles then this will need updating
          coeff=CG(NMNO-1)
          pole =CG(NMNO)
          alpha=CG(NMNO+1)

C         Calculate deriv wrt axial pole
          strain1=strain
          IF(strain.LE.ZERO) THEN !compressive strain
            strain1=ZERO
          ELSE IF(strain.GT.(TOL*pole)) THEN
            strain1=TOL*pole !yielded
          ENDIF
          DENOM=pole-strain1
          DtoAlpha=1.0d0
          IF(alpha.GT.ZERO) DtoAlpha=DENOM**alpha
C         If outside pole-zero range - calc slopes for linear relns
          IF(strain.LE.ZERO.OR.strain.GT.(TOL*pole)) THEN
C           Calc slope
            D_DWslope=(4.0d0*coeff*DENOM*(DENOM+alpha*strain1)
     '        -coeff*(alpha+2.0d0)*(2.0d0*DENOM*(DENOM
     '        +2.0d0*alpha*strain1)
     '        +alpha*(alpha+1.0d0)*strain1*strain1))
     '        /(DtoAlpha*DENOM*DENOM*DENOM)
          ENDIF
          IF(strain.GT.ZERO) THEN !Tensile strain
            D_DW(iterm)=-alpha*coeff*strain1
     '        *(2.0d0+(alpha+1.0d0)*strain1/DENOM)/(DtoAlpha*DENOM)
            IF(strain1.GT.(TOL*pole))
     '        D_DW(iterm)=D_DW(iterm)+(strain-strain1)*D_DWslope
          ELSE !Compressive strain
            D_DW(iterm)=strain*D_DWslope
          ENDIF

C         Set params for shear terms in constitutive law.
          IF(NMNO.EQ.2) THEN !1 dirn
C           params associated with 1-2 deformation mode
            strain1=P4 !EG(1,2)
            iterm1=4
            pole1num=11
C           params associated with 1-3 deformation mode
            strain2=P5 !EG(1,3)
            iterm2=5
            pole2num=14
          ELSE IF(NMNO.EQ.5) THEN  !2 dirn
C           params associated with 2-1 deformation mode
            strain1=P4 !EG(2,1)=EG(1,2)
            iterm1=4
            pole1num=17
C           params associated with 2-3 deformation mode
            strain2=P6 !EG(2,3)
            iterm2=6
            pole2num=20
          ELSE IF(NMNO.EQ.8) THEN  !3 dirn
C           params associated with 3-1 deformation mode
            strain1=P5 !EG(3,1)=EG(1,3)
            iterm1=5
            pole1num=23
C           params associated with 3-2 deformation mode
            strain2=P6 !EG(3,2)=EG(2,3)
            iterm2=6
            pole2num=26
          ENDIF

C         Set parameters appropriate for selected terms of the
C         constitutive law
          coeff1=CG(pole1num-1)
          pole1 =CG(pole1num)
          alpha1=CG(pole1num+1)
          coeff2=CG(pole2num-1)
          pole2 =CG(pole2num)
          alpha2=CG(pole2num+1)

C         Calc deriv of shear poles wrt axial pole to be optimised
          Dpole=2.0d0*(1.0d0+pole)/((1.0d0+2.0d0*pole)**1.5d0)

C         Calc d(DW(iterm1))/d(pole) using chain rule. ie)
C         d(DW(iterm1))/d(pole)=d(DW(iterm1))/d(pole1)*d(pole1)/d(pole)
          strain11=strain1
          IF(strain1.LE.ZERO) THEN !compressive strain
            strain11=ZERO
          ELSE IF(strain1.GT.(TOL*pole1)) THEN
            strain11=TOL*pole1 !yielded
          ENDIF
          DENOM1=pole1-strain11
          DtoAlpha1=1.0d0
          IF(alpha1.GT.ZERO) DtoAlpha1=DENOM1**alpha1
C         If outside pole-zero range - calc slopes for linear relns
          IF(strain1.LE.ZERO.OR.strain1.GT.(TOL*pole1)) THEN
C           Calc slope
            D_DWslope1=Dpole*(4.0d0*coeff1*DENOM1*(DENOM1
     '        +alpha1*strain11)
     '        -coeff1*(alpha1+2.0d0)*(2.0d0*DENOM1*(DENOM1
     '        +2.0d0*alpha1*strain11)
     '        +alpha1*(alpha1+1.0d0)*strain11*strain11))
     '        /(DtoAlpha1*DENOM1*DENOM1*DENOM1)
          ENDIF
          IF(strain1.GT.ZERO) THEN !Tensile strain
            D_DW(iterm1)=0.5d0*Dpole*(-alpha1*coeff1*strain11
     '        *(2.0d0+(alpha1+1.d0)*strain11/DENOM1)/(DtoAlpha1*DENOM1))
            IF(strain1.GT.(TOL*pole1)) D_DW(iterm1)=D_DW(iterm1)
     '        +(strain1-strain11)*D_DWslope1/2.0d0
          ELSE !Compressive strain
            D_DW(iterm1)=strain1*D_DWslope1/2.0d0
          ENDIF

C         Calc d(DW(iterm2))/d(pole) using chain rule. ie)
C         d(DW(iterm2))/d(pole)=d(DW(iterm2))/d(pole2)*d(pole2)/d(pole)
          strain21=strain2
          IF(strain2.LE.ZERO) THEN !compressive strain
            strain21=ZERO
          ELSE IF(strain2.GT.(TOL*pole2)) THEN
            strain21=TOL*pole2 !yielded
          ENDIF
          DENOM2=pole2-strain21
          DtoAlpha2=1.0d0
          IF(alpha2.GT.ZERO) DtoAlpha2=DENOM2**alpha2
C         If outside pole-zero range - calc slopes for linear relns
          IF(strain2.LE.ZERO.OR.strain2.GT.(TOL*pole2)) THEN
C           Calc slope
            D_DWslope2=Dpole*(4.0d0*coeff2*DENOM2*(DENOM2
     '        +alpha2*strain21)
     '        -coeff2*(alpha2+2.0d0)*(2.0d0*DENOM2*(DENOM2
     '        +2.0d0*alpha2*strain21)
     '        +alpha2*(alpha2+1.0d0)*strain21*strain21))
     '        /(DtoAlpha2*DENOM2*DENOM2*DENOM2)
          ENDIF
          IF(strain2.GT.ZERO) THEN !Tensile strain
            D_DW(iterm2)=0.5d0*Dpole*(-alpha2*coeff2*strain21
     '        *(2.0d0+(alpha2+1.d0)*strain21/DENOM2)/(DtoAlpha2*DENOM2))
            IF(strain2.GT.(TOL*pole2)) D_DW(iterm2)=D_DW(iterm2)
     '        +(strain2-strain21)*D_DWslope2/2.0d0
          ELSE !Compressive strain
            D_DW(iterm2)=strain2*D_DWslope2/2.0d0
          ENDIF

        ELSE IF(NMNO.EQ.3.OR.NMNO.EQ.6.OR.NMNO.EQ.9) THEN
C         Parameter is a curvature.
C !!!     NOTE: if fibre distribution model is altered to express shear
C !!!     curvatures in terms of axial curvatures then this will
C !!!     need updating
          coeff=CG(NMNO-2)
          pole =CG(NMNO-1)
          alpha=CG(NMNO)
          IF(strain.LE.ZERO) THEN !compressive strain
            strain1=ZERO
          ELSE IF(strain.GT.(TOL*pole)) THEN
            strain1=TOL*pole !yielded
          ENDIF
          DENOM=pole-strain1
          DtoAlpha=1.0d0
          IF(alpha.GT.ZERO) DtoAlpha=DENOM**alpha
C         If outside pole-zero range - calc slopes for linear relns
          IF(strain.LE.ZERO.OR.strain.GT.(TOL*pole)) THEN
C           Calc slope
            D_DWslope=(coeff*strain1*(4.0d0*DENOM
     '        +strain1*(2.0d0*alpha+1.0d0))
     '        +coeff*(2.0d0*DENOM*(DENOM+2.0d0*alpha*strain1)
     '        +alpha*(alpha+1.0d0)*strain1*strain1)*DLOG(DENOM))
     '        /(DtoAlpha*DENOM*DENOM)
          ENDIF
          IF(strain.GT.ZERO) THEN !Tensile strain
            D_DW(iterm)=coeff*strain1*(strain1/DENOM
     '        -(2.0d0+alpha*strain1/DENOM)*DLOG(DENOM))/DtoAlpha
            IF(strain.GT.(TOL*pole))
     '        D_DW(iterm)=D_DW(iterm)+(strain-strain1)*D_DWslope
          ELSE !Compressive strain
            D_DW(iterm)=strain*D_DWslope
          ENDIF
        ENDIF

      ELSE IF(PARAMTYPE.EQ.'GEOMETRIC_PARAMETERS') THEN

        DO iterm=1,6
          coeff1=CG((iterm-1)*3+1)
          pole1 =CG((iterm-1)*3+2)
          alpha1=CG((iterm-1)*3+3)
          IF(iterm.EQ.1.OR.iterm.EQ.2.OR.iterm.EQ.3) THEN  !Axial
            IF(iterm.EQ.1) THEN
              strain=P1
            ELSE IF(iterm.EQ.2) THEN
              strain=P2
            ELSE IF(iterm.EQ.3) THEN
              strain=P3
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

          strain1=strain
          strain2=strain
          IF(strain.LE.ZERO) THEN !compressive strain
            strain1=ZERO
            strain2=ZERO
          ELSE IF(strain.GT.(TOL*pole1).OR.strain.GT.(TOL*pole2)) THEN
            IF(strain1.GT.(TOL*pole1)) strain1=TOL*pole1 !yielded
            IF(strain1.GT.(TOL*pole2)) strain2=TOL*pole2 !yielded
          ENDIF
          DENOM1=pole1-strain1
          DtoAlpha1=1.0d0
          IF(alpha1.GT.ZERO) DtoAlpha1=DENOM1**alpha1
          DENOM2=pole2-strain2
          DtoAlpha2=1.0d0
          IF(alpha2.GT.ZERO) DtoAlpha2=DENOM2**alpha2
C         If outside pole-zero range - calc slopes for linear relns
          IF(strain.LE.ZERO.OR.strain.GT.(TOL*pole1).OR.
     '      strain.GT.(TOL*pole2)) THEN
C           Calc slope for first set of params
            D_DWslope1=(2.0d0*coeff1*DENOM1*(DENOM1
     '        +2.0d0*alpha1*strain1)
     '        +coeff1*alpha1*(alpha1+1.0d0)*strain1*strain1)
     '        /(DtoAlpha1*DENOM1*DENOM1)
C           Calc slope for second set of params
            D_DWslope2=(2.0d0*coeff2*DENOM2*(DENOM2
     '        +2.0d0*alpha2*strain2)
     '        +coeff2*alpha2*(alpha2+1.0d0)*strain2*strain2)
     '        /(DtoAlpha2*DENOM2*DENOM2)
          ENDIF
          IF(strain.GT.ZERO) THEN !Tensile strain
            D_DW(iterm)=0.5d0*(
     '        2.0d0*coeff1/DtoAlpha1
     '        +4.0d0*coeff1*alpha1*strain1/(DtoAlpha1*DENOM1)
     '        +coeff1*alpha1*(alpha1+1.0d0)*strain1*strain1
     '        /(DtoAlpha1*DENOM1*DENOM1)
     '        +2.0d0*coeff2/DtoAlpha2
     '        +4.0d0*coeff2*alpha2*strain2/(DtoAlpha2*DENOM2)
     '        +coeff2*alpha2*(alpha2+1.0d0)*strain2*strain2
     '        /(DtoAlpha2*DENOM2*DENOM2))
            IF(strain.GT.(TOL*pole1).OR.strain.GT.(TOL*pole2)) THEN
              D_DW(iterm)=ZERO !term above is a const for these cases
              IF(strain.GT.(TOL*pole1))
     '          D_DW(iterm)=D_DW(iterm)+D_DWslope1/2.0d0
              IF(strain.GT.(TOL*pole2))
     '          D_DW(iterm)=D_DW(iterm)+D_DWslope2/2.0d0
            ENDIF
          ELSE !Compressive strain
C           Use deriv of D_DW wrt strain at strain=0 for compress. slope
            D_DW(iterm)=(D_DWslope1+D_DWslope2)/2.0d0
          ENDIF
        ENDDO !iterm
      ENDIF

      CALL EXITS('D_ENERGY')
      RETURN
 9999 CALL ERRORS('D_ENERGY',ERROR)
      CALL EXITS('D_ENERGY')
      RETURN 1
      END


