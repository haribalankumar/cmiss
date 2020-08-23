      SUBROUTINE CALC_DIFF_CAPS(DCO2,Dh,DO2,HCT,SA,SO2,tau_h,
     &  VC,ERROR,*)
        
C#### Subroutine: CALC_DIFF_CAPS
C###  Description:
C###    CALC_DIFF_CAPS: Calculates the oxygen and carbon dioxide
C###     diffusing capacities (diffusive conductances). 
C###     Calculation of DO2 from Weibel 1997.
C###     Calculation of DCO2 from Wagner 1972.
C**** Created by AJS, June 2007.
      
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'pulm00.cmn'    
            
!     Parameter List
      REAL*8 DCO2,Dh,DO2,HCT,SA,SO2,tau_h,VC
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 alpha,CC_Kelman,CO2_INITIAL,CO2_FINAL,CP_Kelman,D_Kelman,
     '  DeCO2,DeO2,DmCO2,DmO2,DpO2,DtO2,Hb_conc,K,kc_O2,length,phi,
     '  tau_hp,tau_ht,thetaO2,thetaCO2
      LOGICAL SEPARATE_LAYERS

      CALL ENTERS('CALC_DIFF_CAPS',*9999)

C...  Description of parameters and variables: 
C      DO2,DCO2 = Diffusing capacities for O2 and CO2, l/s/mmHg
C      K = permeation coeffient oxygen in tissue, m^2/min/mmHg
C      SA = contact surface area between blood and air, m^2
C      tau_h = air-blood barrier thickness, m
C      Vc = volume of capillary blood, litre
C      tau_ht = tissue barrier thickness, m
C      tau_hp = effective plasma barrier thickness, m (Note: This assumes that RBC is cylindrical and occupies the centre of the cylindrical capillary segment)
C      length = effective length of capillary, m (for calculating diff cap of a capillary segment)
C      Dh = hydraulic diameter of capillary, m (for calculating diff cap of a capillary segment)
C      thetaO2,thetaCO2 = reaction rates between gas and blood, ml(gas)/min/mmHg/ml(blood)
C      Hb_conc = concentration of hemoglobin, g(Hb)/100ml(blood)
C      alpha = Bunsen solubility coefficient for O2 in plasma, ml(O2)/ml/atm (Holland:1977, pg 306)
C      phi = correction factor for effective diffusion surface area		
C      kc_O2 = inital reaction velocity of O2 in blood, 1/mM/sec (mM=mmol/litre) (Weibel:1997, pg 1152)
C      HCT = hematocrit, fractional
C      SO2 = oxyhemoglobin saturation fraction
C      thalf = halftime for reaction of CO2 with blood, sec
C      DmO2,DmCO2 = membrane component of diffusing capacities, l/s/mmHg 
C      DeO2,DeCO2 = ethrocyte component of diffusing capacities, l/s/mmHg

      SEPARATE_LAYERS=.FALSE.

      IF(VC.EQ.0.d0.OR.SA.EQ.0.d0) THEN
!        write(*,*) 'CALC_DIFF_CAPS: VC.EQ.0.d0.OR.SA.EQ.0.d0'
        DO2=0.d0
        DCO2=0.d0
      ELSE

C     Calculate O2 diffusing capacity
      K=3.3d-12 
      phi=1.0d0!0.8d0  !removing phi (see more recent articles: Weibel 1997 or 2005)
      alpha=0.0227d0 
      Hb_conc=HCT*MCH/(MCV)/10.0d0 
      kc_O2=440.0d0 ! 1/mM/sec (mM=mmol/litre)(Weibel:1997, pg 1152)
      thetaO2=kc_O2*0.0587d0*alpha*(1.0d0-SO2)*0.01333d0*
     &  Hb_conc*60.d0  !reaction rate, ml O2/min/mmHg/ml
      DeO2=(VC*thetaO2)/60.d0
      IF(SEPARATE_LAYERS)THEN ! consider tissue and plasma separately (appropriate for small scale only)
        length=2.0d0*SQRT(VC/PI)
        tau_ht=tau_h !value set in *.ipmate
        tau_hp=0.5d0*(Dh-2.d0*SQRT(HCT*VC/(PI*length))) !assumes RBC is spherical
        DtO2=(K*phi*SA/tau_ht)
        DpO2=(K*phi*SA/tau_hp)
        DmO2=1.d0/(1.d0/DtO2+1.d0/DpO2)*1.0d3/60.0d0
      ELSE ! consider tissue and plasma as one layer
        DmO2=(K*phi*SA/tau_h)*1.0d3/60.0d0
      ENDIF
      DO2=1.0d0/(1.0d0/DmO2+1.0d0/DeO2)
!       IF(DO2.EQ.0.d0)THEN
!         write(*,*) 'kc_O2',kc_O2,' alpha',alpha,' SO2',SO2
!         write(*,*) ' Hb_conc',Hb_conc,' HCT',HCT,' MCH',MCH,' MCV',MCV
!         write(*,*) 'DmO2=',DmO2,' DeO2=',DeO2,' thetaO2=',thetaO2
!         write(*,*) 'SO2=',SO2,' tau_h=',tau_h,' SA=',SA,' VC=',VC
!       ENDIF

C     Convert PCO2 to CO2 content Kelman:1967 (assume pH=7.4, temp=37 degrees)        
      D_Kelman=0.6617d0
C     Mixed venous blood 
      CP_Kelman=0.6565d0*INITIAL_PCB
      CC_Kelman=D_Kelman*CP_Kelman
      CO2_INITIAL=(HCT*CC_Kelman+(1.d0-HCT)*CP_Kelman)*2.22d0
C     Arterial blood (assume blood equilibrates with alveolar air)
      CP_Kelman=0.6565d0*INITIAL_PCA
      CC_Kelman=D_Kelman*CP_Kelman
      CO2_FINAL=(HCT*CC_Kelman+(1.d0-HCT)*CP_Kelman)*2.22d0
 
C     Calculate CO2 diffusing capacity (Wagner:1972)
      DmCO2=20.d0*DmO2
      thetaCO2=5.4d0 !27.d0 !ml(CO2)/ml(blood)/min/mmHg (mean value from Wagner:1972)
      !thalf=0.150 !sec
      !tau=thalf/log(2.d0)
      !thetaCO2=-(CO2_INITIAL-CO2_FINAL)/tau*exp(-time/tau) !ml(CO2)/ml(blood)/min
      !thetaCO2=-thetaCO2/(INITIAL_PCB-INITIAL_PCA); !ml(CO2)/ml(blood)/min/mmHg
      DeCO2=(VC*thetaCO2)/60.d0
      DCO2=1.0d0/(1.0d0/DmCO2+1.0d0/DeCO2)

      ENDIF !VC.EQ.0.d0.OR.SA.EQ.0.d0

      CALL EXITS('CALC_DIFF_CAPS')
      RETURN
 9999 CALL ERRORS('CALC_DIFF_CAPS',ERROR)
      CALL EXITS('CALC_DIFF_CAPS')
      RETURN 1
      END
