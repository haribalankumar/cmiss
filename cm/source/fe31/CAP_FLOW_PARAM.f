      SUBROUTINE CAP_FLOW_PARAM(ne,L_in,L_out,Ppl,
     &  R_in,R_out,stretch,ERROR,*)
C#### Subroutine CAP_FLOW_PARAM
C###  CREATED: JAN 2010 by ARC
C###  Description:
C###  Assigns values to parameters used in acinar blood flow models that
C###  'step down' from large vessels through the acinus or are 
C###  calculated from whole lung measures. Parameters are either obtained 
C###  from the .ipvalu file for pulmonary perfusion problems or from the 
C###  adjacent extra-acinar vessel. All parameters are defined in common 
C###   block PULM00.cmn.
C###  Called from: SOLVE11
C###  Input:
C###  ne->element number, Ppl->Pleural pressure in that element.
C###  L_in->length of upstream artery,L_out -> length of downstream vein,
C###  R_in->radius of upstream artery
C###  R_out -> radius of downstream vein ,stretch -> element stretch factor

      IMPLICIT NONE
      INCLUDE 'pulm00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'valu00.cmn' 
      !INPUT AND OUTPUT PARAMETER LIST
      INTEGER ne
      REAL*8 L_in,L_out,Ppl,R_in,R_out,stretch
      CHARACTER ERROR*(*)   

      !Local variables
      INTEGER i
      REAL*8 Ptp

      CALL ENTERS('CAP_FLOW_PARAM',*9999)

C... 09/12/10 KSB: CHECK PAREMETERS ARE DEFINED IN .ipvalu FILE. 
      CALL ASSERT(PARAMETERS_DEFINED,
     &  'Define parameter values first (fem define values)',ERROR,*9999)
  
C...  --CAPILLARY SHEET PROPERTIES--
C...  RELATIONSHIP BETWEEN ALVEOLAR VOLUME AND TRANSPULMONARY PRESURE
C...  Transpulmonary pressure should be defined in Pa so convert to cmH20
C...  for use in this relationship
      Ptp=Palv-(-Ppl/98.06d0) !Pa -> cmH2O
C...  This assumes that zero pressure reference volume is 20% of TLC
      IF(Ptp.LT.0.d0) THEN
         L_scale=(1.d0-0.8d0*exp(-0.1d0*0.d0))**(1.d0/3.d0) 
         area_scale=(1.d0-0.8d0*exp(-0.1d0*0.d0))**(2.d0/3.d0) 
      ELSE
         L_scale=(1.d0-0.8d0*exp(-0.1d0*Ptp))**(1.d0/3.d0)
         area_scale=(1.d0-0.8d0*exp(-0.1d0*Ptp))**(2.d0/3.d0) 
      ENDIF

C...   KSB 18/08/09: Including effect of lung inflation on sheet compliance
C...   Below value is scaled from dog lung measurements. This is valid in the 
C...   human only. When we consider other animals we may need to add an 
C...   option which defines alpha_c in each species.
      alpha_c=1.26d0*((-2.04762d-09*Ptp+1.3019d-07)/98.06d0) !Units m/cmH2O -> m/Pa
      CALL ASSERT(alpha_c.GT.0.d0,'>>Capillary alpha < 0!',ERROR,*9999)

C    --ARTERIOLE AND VENULE PROPERTIES--
C###  APPARENT BLOOD VISCOSITY: 
C...  Stepping down linearly with each generation from Fung's
C...  estimate at 45% hematocrit (4e-3Pa.s) to his estimate at 30% hematocrit
C...  (1.92e-3Pa.s). (Biomechanics: Circulation)
      DO i=1,num_symm_gen
         mu_app(i)=(4d-3-(i-1)*(4d-3-1.92d-3)/(num_symm_gen-1));
      ENDDO

C!!!  ALTERNATIVE OPTION FOR APPARENT VISCOSITY
C!!!      u_plasma=1.2d0*1.d-3 !units (cP->Pa.s)
C!!!      Dm=2.70d-6              !(mm)=diam of smallest vessel RBC can pass through
C!!!     Hd_norm=0.4d0 !currently setting this = normal inlet value
C!!!      u_cyto=exp(0.48d0+2.35d0*Hd_Norm)*1.0d-3 !cP->Pa.s (RBC)
C!!!      D_star=(2.03d0-2.0d0*Hd_Norm)*1.0d-6 !um -> m

C...  --- ARTERIOLE and VENULE LENGTHS AND RADII (m) ---
C***  OPTION 1: STEP DOWN LINEARLY FROM INLET TO CAPILLARY RADIUS
C.... ARC Feb 2010: Read in the Radii and the length of the inlet arteries and veins 
C...  (NB length is the same in both the arteries and veins (trees are identical) 
      R_in=R_in/1000.d0    !radius of input artery mm->m
      R_out=R_out/1000.d0  ! radius of outlet vein mm->m
      L_in=L_in/1000.d0
      L_out=L_out/1000.d0
C...  Note that stretch is included here so we don't have to add it later
      DO i=1,num_symm_gen
         L_a(i)=L_in-i*(L_in-L_art_terminal*stretch)/
     &    num_symm_gen
         L_v(i)=L_out-i*(L_in-L_vein_terminal*stretch)/
     &    num_symm_gen
      ENDDO     
      DO i=1,num_symm_gen
        rad_a(i)=(R_art_terminal*sqrt(1.d0/stretch)-R_in)*i
     &    /DBLE(num_symm_gen)+R_in
        rad_v(i)=(R_vein_terminal*sqrt(1.d0/stretch)-R_out)*i
     &    /DBLE(num_symm_gen)+R_out
      ENDDO

C...  -- ARTERIOLE AND VENULE DISTENSIBILITY--
C...  Same values for arterioles and venules based on inlet radius
C...  -Vessel elasticity constant (/Pa)-
      alpha_a=R_in/(PULMAT(3)*1000.d0)   
      alpha_v=R_out/(PULMAT(3)*1000.d0)


      CALL EXITS('CAP_FLOW_PARAM')
      RETURN
 9999 CALL ERRORS('CAP_FLOW_PARAM',ERROR)
      CALL EXITS('CAP_FLOW_PARAM')
      RETURN 1
      END
