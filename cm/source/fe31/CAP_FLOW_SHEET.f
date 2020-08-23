      SUBROUTINE CAP_FLOW_SHEET(ne,SHEET_RES,Q_c,Hart,Hven,RBC_TT,
     & zone,Pin_sheet,Pout_sheet,Ppl,area,area_new,recruited,
     & ERROR,*) 

C#### Subroutine CAP_FLOW_SHEET
C###  CREATED: JUNE 2009 by ARC

C###  Description:
C###  This subroutine calculates resistance across a capillary sheet.
C###  It can be used as a stand-alone code or in combination with
C###  'CAP_FLOW_LADDER.f'
C###  It uses Fung's sheet flow model to calculate flow and resistance 

C###  Called from:
C###  /fe30/CAP_FLOW_SIMPLE - for non-ladder acini
C###  /fe30/CAP_FLOW_LADDER - for ladder acini
C###  /fe07/SOLVE11 - for multibranching acini

C###  Input:
C###  The input to this subroutine is: 
C###  ne=element number
C###  Pin_sheet= Pressure into the sheet
C###  Pout_sheet=Pressure out of the sheet
C###  Ppl= Pleural pressure
C###  area= Sheet area - this is an input as it depends on the 
C###        and venule structure.

C###  Output:
C###  Important output to large vessel models 
C###  SHEET_RES= Resistance across the sheet
  
C###  
C###  In addition this subroutine outputs:
C###  Flow and RBC transit times through the sheet. 

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'pulm00.cmn'

      !INPUT AND OUTPUT PARAMETER LIST
      INTEGER ne,zone
      REAL*8 SHEET_RES,Q_c,Hart,Hven,Pin_sheet,Pout_sheet,Ppl,area    
      CHARACTER ERROR*(*)

      !Local variables
      REAL*8 area_new,C,Hmax_art,Hmax_ven,Hub,L_new,P_a,Palv_pa,P_ave,
     &    P_v,RBC_tt,recruited,waterfall_scale

 
C###  ENTER SUBROUTINE 'CAP_FLOW_SHEET'
      CALL ENTERS('CAP_FLOW_SHEET',*9999)
C...  area of sheet is an input. Now area and pathlength need to be scaled. 
      area_new=area*area_scale;
      L_new=L_c*L_scale
C... Maximum sheet heights
      Hmax_art=H0+alpha_c*Pub_c 
      Hmax_ven=Hmax_art
C...  Blood pressure in and out of capillary sheet
      P_a=Pin_sheet           
      P_v=Pout_sheet
      Palv_pa=Palv*98.06d0
C...   CAPILLARY RECRUITMENT
C...   CALCULATE AVERAGE CAPILLARY BLOOD PRESSURE
       P_ave=(P_a+P_v)/2.d0
C      calculate recruitment proportion
       recruited=1.d0-F_rec*EXP(-(P_ave**2.d0)/(sigma_rec**2.d0))
C      multiply flows and transit times
       area_new=area_new*recruited

       IF((P_v-Palv_pa).LE.Plb_c)THEN
C ### Zone 2 (waterfaall) scaling 
         waterfall_scale=1-F_sheet+F_sheet*
     &   EXP(-(P_v-Palv_pa)**2.d0/(2.d0*sigma_cap**2.d0))
         area_new=area_new*waterfall_scale
       ENDIF
C###  -----CALCULATING SHEET FLOW----
C...  Sheet flow model constant, C
      C=area_new/(4.d0*mu_c*K_cap*F_cap*L_new**2*alpha_c)
C###  Determine what zone we are in and calculate flow and TT
      IF((P_a-Palv_pa).LT.Plb_c)THEN 
C...   ZONE 1:
C...   Arteriole and venous pressure both less than alveolar pressure - the
C...   sheet is collapsed
         !write(*,*) 'sheet collapsed'
         zone=1
         Hart=0.d0
         Hven=0.d0            
         Q_c=0.d0
         RBC_tt=0.d0 !TEMPORARY
        ELSEIF((P_a-Palv_pa).LE.Pub_c.AND.(P_v-Palv_pa).LE.Plb_c)THEN
C...    ZONE 2:
           zone=2
           Hart=H0+alpha_c*(P_a-Palv_pa)
           Hven=H0
           Q_c=C*(Hart**4.d0-H0**4.d0)
           RBC_tt=(3.d0*mu_c*F_cap*K_cap*L_new**2.d0*alpha_c)
     &     /(Hart**3.d0-H0**3.d0)
        ELSEIF((P_a-Palv_pa).GT.Pub_c.AND.(P_v-Palv_pa).LE.Plb_c)THEN
C...       ZONE 2:
           zone=2
           Hart=Hmax_art
           Hven=H0
           Q_c=4.d0*C*alpha_c*(Hmax_art**3.d0*(P_a-Palv_pa-Pub_c)
     &       +(Hmax_art**4.d0-H0**4.d0)/(4.d0*alpha_c))
           RBC_tt=(3.d0*mu_c*F_cap*K_cap*L_new**2.d0*alpha_c)/
     &       ((3.d0*alpha_c*Hmax_art**2.d0*
     &         (P_a-Palv_pa-Pub_c)+(Hmax_art**3.d0-H0**3.d0))) 
        ELSEIF((P_a-Palv_pa).LE.Pub_c.AND.(P_v-Palv_pa).LE.Pub_c)THEN
C...       ZONE3 or 4: The boundary between zone 3 and 4 is not clearcut.
           zone=3 !tmp!!should = 3
           Hart=H0+alpha_c*(P_a-Palv_pa)
           Hven=H0+alpha_c*(P_v-Palv_pa)
           Q_c=C*(Hart**4.d0-Hven**4.d0)
           RBC_tt=(3.d0*mu_c*F_cap*K_cap*L_new**2.d0*alpha_c)
     &              /(Hart**3.d0-Hven**3.d0)
        ELSEIF((P_a-Palv_pa).GT.Pub_c.AND.(P_v-Palv_pa).LE.Pub_c)THEN
           zone=3
           Hart=Hmax_art
           Hven=H0+alpha_c*(P_v-Palv_pa)
           Hub=H0+alpha_c*Pub_c
           Q_c=4.d0*C*alpha_c*(Hmax_art**3.d0*(P_a-Palv_pa-Pub_c)+
     &          (Hub**4.d0-Hven**4.d0)/(4*alpha_c))
           RBC_tt=(3.d0*mu_c*F_cap*K_cap*L_new**2*alpha_c)
     &      /(3.d0*alpha_c*Hmax_art**2.d0*(P_a-Palv_pa-Pub_c)+
     &            (Hmax_art**3.d0-Hven**3.d0))
        ELSEIF((P_a-Palv_pa).GT.Pub_c.AND.(P_v-Palv_pa).GT.Pub_c)THEN
           zone=4 !!!tmp should = 3 
           Hart=Hmax_art
           Hven=Hmax_ven
           Q_c=4.d0*C*alpha_c*Hart**3.d0*(P_a-P_v)
           RBC_tt=(mu_c*F_cap*K_cap*L_new**2.d0)/(Hart**2.d0*(P_a-P_v))
        ELSE
C....   SOMETHING HAS GONE TERRIBLY WRONG AS THIS SHOULD BE ALL THE OPTIONS!!
           Q_c=1.d-8
           RBC_tt=0.d0
       write(*,*) 'aaaaarrrrgh, something has gone terribly wrong!'
       ENDIF

C...   RBC transit time (1.4 times faster than blood)
       RBC_tt=RBC_tt/1.4d0 ;
      IF(alpha_c.eq.0d0)THEN
      C=area_new/(4.d0*mu_c*K_cap*F_cap*L_new**2)
           Hart=Hmax_art
           Hven=Hmax_ven
           Q_c=4.d0*C*Hart**3.d0*(P_a-P_v)
           RBC_tt=(mu_c*F_cap*K_cap*L_new**2.d0)/(Hart**2.d0*(P_a-P_v))
      ENDIF

      IF(Q_c.le.0.d0)THEN ! to acount for small negative flows when pressure in and out of a capillary are the same. 
         zone=5
       Q_c=1.d-15
      ENDIF
    
      IF(zone.eq.1.or.zone.eq.5)THEN 
C...   ZONE 1:
       SHEET_RES=1d12
      ELSE
C...  Resistance through a single capillary sheet, Pa.s/m^3
       SHEET_RES=(Pin_sheet-Pout_sheet)/Q_c
      ENDIF

C###  EXIT SUBROUTINE 'CAP_FLOW_SHEET'
      CALL EXITS('CAP_FLOW_SHEET')
      RETURN
 9999 CALL ERRORS('CAP_FLOW_SHEET',ERROR)
      CALL EXITS('CAP_FLOW_SHEET')
      RETURN 1
      END
