      SUBROUTINE METABOLISM_CELL(TIME,Y,DY,CONTROL,MODEL,SIZES,
     '  VARIANT,DERIVED,PARAM,PROTOCOL,AII,AIO,ARI,ARO,ERR_CODE)

C#### Subroutine: METABOLISM_CELL
C###  Description:
C###    Solve using the cardiac myocyte metabolic equations developed by
C###    Peter J Mulquiney
CC***  Created by Nicolas Smith Aug 2000

      IMPLICIT NONE

      INCLUDE 'cell_metab.inc'
      INCLUDE 'cell_reserved.inc'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER SIZES(11)
      INTEGER CONTROL(*),MODEL(*),VARIANT,AII(*),AIO(*),ERR_CODE
      REAL*8 TIME(*),Y(*),PARAM(*),PROTOCOL(*),ARI(*),DY(*),
     '  DERIVED(*),ARO(*)
!     Local variables
      REAL*8 voxphos1,voxphos2,vglycl,vglyc2,ff,katpase1, katpase2,
     '  vatpase1,vatpase2, K_napump_nai,K_napump_ki,K_napump_ke,
     '  K_napump_ne,K_napump_mgatp,kapp_napump1f, kapp_napump1b,
     '  kapp_napump2f,kapp_napump2b,K_napump_nae,kapp_napump3f,
     '  kapp_napump3b,kapp_napump4f,kapp_napump4b,
     '  N1napump,N2napump,D1napump,D2napump,D3napump,D4napump,
     '  D5napump,D6napump,D7napump,D8napump,D9napump,D10napump,
     '  D11napump,D12napump,D13napump,D14napump,D15napump,D16napump,
     '  denom,v_napump,kapp_mct1,kapp_mct1b,k_lactrans2,k_lactrans2b,
     '  vmct,vlactrans,vphostrans,vsarccampump,N1_nhe,N2_nhe,D1_nhe,
     '  D2_nhe,D3_nhe,D4_nhe,D5_nhe,D6_nhe,D7_nhe,D8_nhe,denom_nhe,
     '  vco2trans,vna_influx,vke_flux,vnce,vnhe,knbc_2b, N1_nbc,N2_nbc,
     '  D1_nbc,D2_nbc,D3_nbc,D4_nbc,D5_nbc,D6_nbc,D7_nbc,D8_nbc,
     '  denom_nbc,v_nbc,kae_2b,N1_ae,N2_ae,D1_ae,D2_ae,D3_ae,
     '  D4_ae,D5_ae, D6_ae,D7_ae,D8_ae,denom_ae,vae,N1_che,N2_che,
     '  D1_che,D2_che, D3_che,D4_che,D5_che,D6_che,D7_che,D8_che,
     '  denom_che,v_che,v_fh2co3,v_fh2co3e,fatp,fadp,famp,k_ck_f,
     '  k_ck_r,v_ck,kak_f,kak_r,v_ak,Kmgatp,Kapp_mgatp,kfmgatp_f,
     '  vfmgatp,Kmgadp,Kapp_mgadp,kfmgadp_f,vfmgadp,Kapp_mghphos,
     '  Kmghphos,fhp,
     '  kfmghphos_f,vfmghphos,Kapp_mgPCr,k_fmgPCr_f,vfmgPCr,Kapp_hatp,
     '  kfhatp_f,vfhatp,Kapp_hadp,kfhadp_f,vfhadp,Kapp_hamp,kfhamp_f,
     '  vfhamp,Kapp_hphos,kfh2phos_f,vfh2phos,vfh2phose,
     '  vkintrinsicbuffer1,vkintrinsicbuffer2,vexternalproteinbuffer,
     '  v_trop,v_cal,test1,test2,test3,test4,test5,test6
C *** DPN 20 October 2000 - Removing h2phos from above declaration since
C     it is given in cell_metab.inc as an integer
      LOGICAL DEBUG

      ERR_CODE=0

      DEBUG=.FALSE.


c Rate equations

c-----------------------------oxydative phosphorylation ----------------

      voxphos1=PARAM(Voxphos)*Y(MgADP)/(Y(MgADP)+PARAM(Koxphos)) !oxydative phosphorylation
      voxphos2=0.0d0


c-----------------------------gycolgenolysis----------------

      vglycl=0.0d0
      vglyc2=PARAM(Vglyco)*Y(Glyc)/(Y(Glyc)+PARAM(Kglyco))


c--------------------------ATPases-------------------

      ff=1-((TIME(TCELL)/1000.0d0)**3.0d0) !fudge factor

      katpase1=0.00026d0/(Y(MgATP)+0.0005)
      katpase2=0.00026*ff/(Y(MgATP)+0.0005)

      vatpase1=katpase1*Y(MgATP)
      vatpase2=katpase2*Y(MgATP)

c--------------------------Sodium Pump-------------------

      PARAM(knapump4b)=894.328d0/(1+(5.01187e-07/Y(h)))



      K_napump_nai=2e-3*exp(-0.2d0/3.0d0*PARAM(F)*Y(Vm)/
     '  (PARAM(R)*PARAM(Temp)))

      K_napump_ki=6e-3*exp(-0.0d0/2.0d0*PARAM(F)*Y(Vm)/!**zeros and exponent
     '  (PARAM(R)*PARAM(Temp)))

      K_napump_ke=0.8e-3*exp(0.0d0/2.0d0*PARAM(F)*Y(Vm)/!**zeros and exponent
     '  (PARAM(R)*PARAM(Temp)))

      K_napump_nae=200e-3*exp(0.8d0/3.0d0*PARAM(F)*Y(Vm)/!**zeros and exponent
     '  (PARAM(R)*PARAM(Temp)))

      K_napump_mgatp=0.6e-03


      PARAM(knapump3f)=(PARAM(knapump4b)*PARAM(knapump3b)*
     '  PARAM(knapump1f)*PARAM(knapump2f)*(K_napump_ki**2.0d0)*
     '  (K_napump_nae**3.0d0))/(PARAM(knapump2b)*PARAM(knapump1b)*
     '  PARAM(knapump4f)*(K_napump_ke**2.0d0)*K_napump_mgatp*
     '  (K_napump_nai**3.0d0)*0.011d0)


      kapp_napump1f=PARAM(knapump1f)/(K_napump_nai**3.0)/(  !**brackets??
     '  1+((Y(K)/K_napump_ki)**2.0d0)+(2.0d0*Y(K)/K_napump_ki)+
     '  (3.0d0*Y(Na)/K_napump_nai)+
     '  (3.0d0*(Y(Na)**2.0d0)/(K_napump_nai**2.0d0))+
     '  ((Y(Na)**3.0d0)/(K_napump_nai**3.0d0)))

       kapp_napump1b=-PARAM(knapump1b)
       kapp_napump2f=PARAM(knapump2f)


       kapp_napump2b=PARAM(knapump2b)/(K_napump_nae**3.0)/(  !**brackets??
     '  1+((Y(Ke)/K_napump_ke)**2.0d0)+(2.0d0*Y(Ke)/K_napump_ke)+
     '  (3.0d0*Y(Nae)/K_napump_nae)+
     '  (3.0d0*(Y(Nae)**2.0d0)/(K_napump_nae**2.0d0))+
     '  ((Y(Nae)**3.0d0)/(K_napump_nae**3.0d0)))

       kapp_napump3f=PARAM(knapump3f)/(K_napump_ki**2.0)/(  !**brackets??
     '  1+((Y(K)/K_napump_ki)**2.0d0)+(2.0d0*Y(K)/K_napump_ki)+
     '  (3.0d0*Y(Na)/K_napump_nai)+
     '  (3.0d0*(Y(Na)**2.0d0)/(K_napump_nai**2.0d0))+
     '  ((Y(Na)**3.0d0)/(K_napump_nai**3.0d0)))

       kapp_napump3b=(PARAM(knapump3b)/K_napump_mgatp)/
     '   (1+(Y(MgATP)/K_napump_mgatp))

       kapp_napump4f=PARAM(knapump4f)/
     '   (1+(Y(MgATP)/K_napump_mgatp))

      kapp_napump4b=PARAM(knapump4b)/(K_napump_ke**2.0)/(  !**brackets??
     '  1+((Y(Ke)/K_napump_ke)**2.0d0)+(2.0d0*Y(Ke)/K_napump_ke)+
     '  (3.0d0*Y(Nae)/K_napump_nae)+
     '  (3.0d0*(Y(Nae)**2.0d0)/(K_napump_nae**2.0d0))+
     '  ((Y(Nae)**3.0d0)/(K_napump_nae**3.0d0)))


C note should be able to simplify above with equations for kapp_napump2b

      N1napump=kapp_napump4b*kapp_napump3b*kapp_napump1f*kapp_napump2f
      N2napump=kapp_napump2b*kapp_napump1b*kapp_napump3f*kapp_napump4f
      D1napump=kapp_napump2f*kapp_napump3f*kapp_napump4f
      D2napump=kapp_napump4b*kapp_napump2f*kapp_napump3f
      D3napump=kapp_napump1b*kapp_napump3f*kapp_napump4f
      D4napump=kapp_napump4b*kapp_napump1b*kapp_napump3f
      D5napump=kapp_napump4b*kapp_napump3b*kapp_napump2f
      D6napump=kapp_napump4b*kapp_napump3b*kapp_napump1b
      D7napump=kapp_napump1f*kapp_napump2f*kapp_napump4f
      D8napump=kapp_napump4b*kapp_napump1f*kapp_napump2f
      D9napump=kapp_napump3b*kapp_napump1f*kapp_napump2f

      D10napump=kapp_napump4b*kapp_napump3b*kapp_napump1f
      D11napump=kapp_napump2b*kapp_napump3f*kapp_napump4f
      D12napump=kapp_napump2b*kapp_napump1b*kapp_napump3f
      D13napump=kapp_napump2b*kapp_napump1b*kapp_napump4f
      D14napump=kapp_napump3b*kapp_napump2b*kapp_napump1b
      D15napump=kapp_napump2b*kapp_napump1f*kapp_napump4f
      D16napump=kapp_napump3b*kapp_napump2b*kapp_napump1f


      denom=D1napump*Y(hphos)*Y(h)*(Y(K)**2.0)+D2napump*(Y(K)**2.0)*
     '  (Y(Ke)**2.0)+
     '  D3napump*Y(hphos)*Y(h)*(Y(K)**2.0)*Y(MgATP)+
     '  D4napump*(Y(K)**2.0d0)*(Y(Ke)**2.0)*Y(MgADP)+
     '  D5napump*(Y(Ke)**2.0d0)*Y(MgADP)+
     '  D6napump*(Y(Ke)**2.0)*Y(MgATP)*Y(MgADP)+
     '  D7napump*Y(hphos)*Y(h)*(Y(Na)**3.0d0)+
     '  D8napump*(Y(Ke)**2.0d0)*(Y(Na)**3.0d0)+
     '  D9napump*Y(MgATP)*(Y(Na)**3.0d0)+
     '  D10napump*(Y(Ke)**2.0d0)*Y(MgATP)*(Y(Na)**3.0d0)+
     '  D11napump*(Y(Na)**3.0d0)*Y(hphos)*Y(h)*(Y(K)**2.0d0)+
     '  D12napump*(Y(K)**2.0d0)*Y(MgADP)*(Y(Nae)**3.0d0)+
     '  D13napump*Y(hphos)*Y(h)*Y(MgADP)*(Y(Nae)**3.0d0)+
     '  D14napump*Y(MgADP)*Y(MgATP)*(Y(Nae)**3.0d0)+
     '  D15napump*Y(hphos)*Y(h)*(Y(Na)**3.0d0)*(Y(Nae)**3.0d0)+
     '  D16napump*Y(MgATP)*(Y(Na)**3.0d0)*(Y(Nae)**3.0d0)


      v_napump=(PARAM(conapump)/denom)*(
     '  (N1napump*(Y(Ke)**2.0d0)*Y(MgATP)*(Y(Na)**3.0d0))-
     '  (N2napump*Y(hphos)*Y(h)*(Y(K)**2.0d0)*Y(MgADP)*
     '  (Y(Nae)**3.0d0)))

c--------------------Metabolite Transport--------------

      kapp_mct1=PARAM(kmct1)/(1.0d0+
     '  (Y(h)/PARAM(Kmct_h))+
     '  (Y(Lac)*Y(h)/(PARAM(Kmct_h)*PARAM(Kmct_lac))))

      kapp_mct1b=PARAM(kmct1b)/(1.0d0+
     '  (Y(he)/PARAM(Kmct_h))+
     '  (Y(Lace)*Y(he)/(PARAM(Kmct_h)*PARAM(Kmct_lac))))

      k_lactrans2=(PARAM(kmct2)*Y(Lac)*Y(h)/(PARAM(Kmct_h)*
     '  PARAM(Kmct_lac)))/
     '  (1.0d0+(Y(h)/PARAM(Kmct_h))+
     '  (Y(Lac)*Y(h)/(PARAM(Kmct_h)*PARAM(Kmct_lac))))

      k_lactrans2b=(PARAM(kmct2b)*Y(Lac)*Y(he)/(PARAM(Kmct_h)*
     '  PARAM(Kmct_lac)))/
     '  (1.0d0+(Y(he)/PARAM(Kmct_h))+
     '  (Y(Lace)*Y(he)/(PARAM(Kmct_h)*PARAM(Kmct_lac))))

      vmct=(kapp_mct1+k_lactrans2)*Y(MCT)-
     '  (kapp_mct1b+k_lactrans2b)*Y(MCTe)

      vlactrans=k_lactrans2*Y(MCT)-k_lactrans2b*Y(MCTe)


C--------------------------phosphate transport----------


      vphostrans=-(PARAM(V_NaPi_symport)*Y(h2phos))/
     '   (Y(h2phos)+PARAM(K_NaPi_symport))+
     '   (PARAM(V_Piefflux)*Y(h2phos))/
     '   (Y(h2phos)+PARAM(K_Piefflux))


C--------------------------C02 Transport-------------------

      vco2trans=PARAM(rho)*PARAM(P_co2trans)*(Y(CO2)-Y(CO2e))


C--------------------------ion tansport-------------------

      vna_influx=0.40e-04*(Y(Nae)-Y(Na))
      vke_flux=0.52e-04*(Y(K)-Y(Ke))

C--------------------------NCE-------------------

      test1=(Y(Na)**3.0d0)*Y(Cae)
      test2=exp(0.35d0*Y(Vm)*PARAM(F)/(PARAM(R)*PARAM(Temp)))
      test3=exp(-0.65d0*Y(Vm)*PARAM(F)/(PARAM(R)*PARAM(Temp)))
      test4=(Y(Nae)**3.0d0)*Y(Ca)
      test5=((PARAM(Knce_nae)**3.0d0)+(Y(Nae)**3.0d0))
     '  *(PARAM(Knce_cae)+Y(Cae))
      test6=(1.0d0+0.1d0*exp(-0.65d0*Y(Vm)*PARAM(F)
     '  /(PARAM(R)*PARAM(Temp))))

      vnce=(PARAM(Vncef)*
     '  (exp(0.35d0*Y(Vm)*PARAM(F)/(PARAM(R)*PARAM(Temp)))*
     '  (Y(Na)**3.0d0)*Y(Cae)-
     '  exp(-0.65d0*Y(Vm)*PARAM(F)/(PARAM(R)*PARAM(Temp)))*
     '  (Y(Nae)**3.0d0)*Y(Ca)))/
     '  (((PARAM(Knce_nae)**3.0d0)+(Y(Nae)**3.0d0))*
     '  (PARAM(Knce_cae)+Y(Cae))*
     '  (1.0d0+0.1d0*exp(-0.65d0*Y(Vm)*PARAM(F)
     '  /(PARAM(R)*PARAM(Temp)))))

C----------------------Sarcolemal Ca Pump---------------

      vsarccampump=(7.0d0**(-5.0))*Y(Ca)/(0.5e-06+Y(Ca))

C----------------------Protein Transport--------------


      PARAM(Knhe_Na)=PARAM(Knhe_h)*PARAM(Knhe_Nae)/PARAM(Knhe_he)
      PARAM(knhe2b)=PARAM(knhe1b)*PARAM(knhe2)/PARAM(knhe1)

      test6= PARAM(Knhe_Na)
      N1_nhe=PARAM(Conhe)*PARAM(Knhe_Na)*PARAM(Knhe_he)*
     '    PARAM(knhe1b)*PARAM(knhe2)
      N2_nhe=PARAM(Conhe)*PARAM(Knhe_h)*PARAM(Knhe_Nae)*
     '    PARAM(knhe2b)*PARAM(knhe1)

      D1_nhe=PARAM(Knhe_h)*PARAM(Knhe_Na)*PARAM(Knhe_Nae)*
     '    PARAM(knhe2b)
      D2_nhe=PARAM(Knhe_he)*PARAM(Knhe_Na)*PARAM(Knhe_Nae)*
     '    PARAM(knhe2)
      D3_nhe=PARAM(Knhe_h)*PARAM(Knhe_Na)*PARAM(Knhe_he)*
     '    PARAM(knhe1b)
      D4_nhe=PARAM(Knhe_h)*PARAM(Knhe_he)*PARAM(Knhe_Nae)*
     '    PARAM(knhe1)
      D5_nhe=PARAM(Knhe_Na)*PARAM(Knhe_Nae)*(PARAM(knhe2b)+
     '    PARAM(knhe2))
      D6_nhe=PARAM(Knhe_h)*PARAM(Knhe_Nae)*(PARAM(knhe2b)+
     '    PARAM(knhe1))
      D7_nhe=PARAM(Knhe_Na)*PARAM(Knhe_he)*(PARAM(knhe1b)+
     '    PARAM(knhe2))
      D8_nhe=PARAM(Knhe_h)*PARAM(Knhe_he)*(PARAM(knhe1)+
     '    PARAM(knhe1b))

      denom_nhe=D1_nhe*Y(he)+D2_nhe*Y(h)+D3_nhe*Y(Nae)+
     '    D4_nhe*Y(Na)+D5_nhe*Y(h)*Y(he)+
     '    D6_nhe*Y(he)*Y(Na)+D7_nhe*Y(h)*Y(Nae)+
     '    D8_nhe*Y(Na)*Y(Nae)

      vnhe=((1.0d0+(PARAM(Knhe)/Y(h))**PARAM(nnhe))**(-1.0d0))*
     '    ((N1_nhe*Y(h)*Y(Nae)-N2_nhe*Y(he)*PARAM(Na))
     '    /denom_nhe)

C------------------NBC (sodium bicarbonate) tansport---------------

      knbc_2b=PARAM(knbc_1b)*PARAM(knbc_2)/PARAM(knbc_1)
      PARAM(Knbc_hco3)=0.00354813**PARAM(nnbc)


      N1_nbc=PARAM(Co_nbc)*PARAM(Knbc_Na)*PARAM(Knbc_hco3)*
     '    PARAM(knbc_1b)*PARAM(knbc_2)
      N2_nbc=PARAM(Co_nbc)*PARAM(Knbc_Na)*PARAM(Knbc_hco3)*
     '    PARAM(knbc_1)*knbc_2b

      D1_nbc=(PARAM(Knbc_hco3)**2.0d0)*(PARAM(Knbc_Na)**2.0d0)*
     '    (PARAM(knbc_1b)+PARAM(knbc_1))

      D2_nbc=PARAM(Knbc_hco3)*(PARAM(Knbc_Na)**2.0d0)*
     '    PARAM(knbc_1b)

      D3_nbc=PARAM(Knbc_hco3)*(PARAM(Knbc_Na)**2.0d0)*
     '    PARAM(knbc_1)

      D4_nbc=PARAM(Knbc_Na)*PARAM(Knbc_hco3)*
     '    (PARAM(knbc_1b)+PARAM(knbc_2))

      D5_nbc=PARAM(Knbc_Na)*PARAM(Knbc_hco3)*
     '    (knbc_2b+PARAM(knbc_1))

      D6_nbc=PARAM(Knbc_Na)*PARAM(knbc_2)

      D7_nbc=PARAM(Knbc_Na)*knbc_2b

      D8_nbc=knbc_2b+PARAM(knbc_2)


      denom_nbc=D1_nbc+D2_nbc*(Y(hco3)**PARAM(nnbc))+
     '     D3_nbc*(Y(hco3e)**PARAM(nnbc))+
     '     D4_nbc*Y(Na)*(Y(hco3)**PARAM(nnbc))+
     '     D5_nbc*Y(Nae)*(Y(hco3e)**PARAM(nnbc))+
     '     D6_nbc*Y(Na)*(Y(hco3)**PARAM(nnbc))*(Y(hco3e)**PARAM(nnbc))+
     '     D7_nbc*(Y(hco3)**PARAM(nnbc))*(Y(hco3e)**PARAM(nnbc))*Y(Nae)+
     '     D8_nbc*(Y(hco3)**PARAM(nnbc))*(Y(hco3e)**PARAM(nnbc))
     '     *Y(Nae)*Y(Na)


      v_nbc=(-1.0d0/denom_nbc)*(N1_nbc*(Y(hco3)**PARAM(nnbc))*Y(Na)-
     '     N2_nbc*(Y(hco3e)**PARAM(nnbc))*Y(Nae))


c-------------------------Ac-----------------------------------

      kae_2b=PARAM(kae_2)*PARAM(kae_1b)/PARAM(kae_1)

      N1_ae=PARAM(Coae)*PARAM(Kae_Cl)*PARAM(Kae_hco3e)*
     '     PARAM(kae_1b)*PARAM(kae_2)
      N2_ae=PARAM(Coae)*PARAM(Kae_Cle)*PARAM(Kae_hco3)*
     '     kae_2b*PARAM(kae_1)

      D1_ae=PARAM(Kae_hco3)*PARAM(Kae_Cl)*PARAM(Kae_Cle)*kae_2b
      D2_ae=PARAM(Kae_Cl)*PARAM(Kae_hco3e)*PARAM(Kae_Cle)*PARAM(kae_2)
      D3_ae=PARAM(Kae_hco3)*PARAM(Kae_Cl)*PARAM(Kae_hco3e)*PARAM(kae_1b)
      D4_ae=PARAM(Kae_hco3)*PARAM(Kae_hco3e)*PARAM(Kae_Cle)*PARAM(kae_1)
      D5_ae=PARAM(Kae_Cl)*PARAM(Kae_Cle)*(kae_2b+PARAM(kae_2))
      D6_ae=PARAM(Kae_hco3)*PARAM(Kae_Cle)*(kae_2b+PARAM(kae_1))
      D7_ae=PARAM(Kae_Cl)*PARAM(Kae_hco3e)*(PARAM(kae_1b)+PARAM(kae_2))
      D8_ae=PARAM(Kae_hco3)*PARAM(Kae_hco3e)*
     '     (PARAM(kae_1)+PARAM(kae_1b))

      denom_ae=D1_ae*Y(hco3e)+D2_ae*Y(hco3)+D3_ae*Y(Nae)+
     '     D4_ae*Y(Na)+D5_ae*Y(hco3)*Y(hco3e)+
     '     D6_ae*Y(hco3e)*Y(Na)+D7_ae*Y(hco3)*Y(Nae)+
     '     D8_ae*Y(Na)*Y(Nae)


      vae=((1.0d0+(PARAM(Kae)/Y(h))**(-PARAM(n_ae)))**(-1.0d0))*
     '     ((N1_ae*Y(hco3)*Y(Cle)-N2_ae*Y(hco3e)*Y(Cl))/denom_ae)


c----------------------CHE-----------------------------


      N1_che=PARAM(Co_che)*PARAM(Kche_Cl)*PARAM(Kche_ohe)*
     '     PARAM(kche_1b)*PARAM(kche_2)
      N2_che=PARAM(Co_che)*PARAM(Kche_oh)*PARAM(Kche_Cle)*
     '     PARAM(kche_2b)*PARAM(kche_1)

      D1_che=PARAM(Kche_oh)*PARAM(Kche_Cl)*PARAM(Kche_Cle)*
     '     PARAM(kche_2b)
      D2_che=PARAM(Kche_Cl)*PARAM(Kche_ohe)*PARAM(Kche_Cle)*
     '     PARAM(kche_2)
      D3_che=PARAM(Kche_oh)*PARAM(Kche_Cl)*PARAM(Kche_ohe)*
     '     PARAM(kche_1b)
      D4_che=PARAM(Kche_oh)*PARAM(Kche_ohe)*PARAM(Kche_Cle)*
     '     PARAM(kche_1)
      D5_che=PARAM(Kche_Cl)*PARAM(Kche_Cle)*(PARAM(kche_2b)+
     '     PARAM(kche_2))
      D6_che=PARAM(Kche_oh)*PARAM(Kche_Cle)*(PARAM(kche_2b)+
     '     PARAM(kche_1))
      D7_che=PARAM(Kche_Cl)*PARAM(Kche_ohe)*(PARAM(kche_1b)+
     '     PARAM(kche_2))
      D8_che=PARAM(Kche_oh)*PARAM(Kche_ohe)*(PARAM(kche_1)+
     '     PARAM(kche_1b))


      denom_che=D1_che*1e-14/Y(he)+D2_che*1e-14/Y(h)+
     '     D3_che*Y(Nae)+D4_che*Y(Na)+D5_che*1e-28/Y(h)/Y(he)+
     '     D6_che*1e-14/Y(he)*Y(Na)+D7_che*1e-14/Y(h)*Y(Nae)+
     '     D8_che*Y(Na)*Y(Nae)

      v_che=(1.0d0/denom_che)*(N1_che*1e-14/Y(h)*Y(Cle)-
     '     N2_che*1e-14/Y(he)*Y(Cl))

c-------------------proton buffering---------------------
c-------------------caribuc abydrase---------------------

      v_fh2co3=PARAM(kfh2co3_f)*Y(h)*Y(hco3)-PARAM(kfh2co3_r)*Y(CO2)


c-------------------hydration of CO2---------------------

      v_fh2co3e=PARAM(kfh2co3e_f)*Y(he)*Y(hco3e)-PARAM(kfh2co3e_r)
     '  *Y(CO2e)


c-------------------Fast reactions-----------------------
c-------------------Energy Metabolism--------------------
c-------------------creatine kinase----------------------

      fatp=1.0d0+PARAM(k1)*PARAM(Kkatp) ! ???
      fadp=1.0d0+PARAM(k1)*PARAM(Kkadp) ! ???
      famp=1.0d0+PARAM(k1)*PARAM(Kkamp) ! ???

      k_ck_f=3.88e08*1.2e03*fatp
      k_ck_r=1.2e03*fadp

      v_ck=k_ck_f*Y(PCr)*Y(MgADP)*Y(h)-k_ck_r*Y(Cr)*Y(MgATP)

c-----------------------adenylate kinase-------------------

      kak_f=3.8e03*famp
      kak_r=1.0e03*fadp
      v_ak=kak_f*Y(MgADP)*Y(ADP)-kak_r*Y(MgATP)*Y(AMP)

c------------------------Mg buffering---------------------------

      Kmgatp=2.6e04*(1+0.15*PARAM(Kkatp)+(10**(-7.2))*PARAM(Khatp))
      Kapp_mgatp=Kmgatp/fatp
      kfmgatp_f=Kapp_mgatp*1.2e03
      vfmgatp=kfmgatp_f*Y(Mg)*Y(ATP)-PARAM(kfmgatp_r)*Y(MgATP)

c------------------------MgADP---------------------------------

      Kmgadp=
     '  2.3e03*(1.0d0+0.15d0*PARAM(Kkadp)+(10**(-7.2))*PARAM(Khadp))
      Kapp_mgadp=Kmgadp/fadp
      kfmgadp_f=Kapp_mgadp*1.2e03
      vfmgadp=kfmgadp_f*Y(Mg)*Y(ADP)-PARAM(kfmgadp_r)*Y(MgADP)

c-----------------------Mghphos---------------------------------

      Kmghphos=34.0d0*(1.0d0+0.15d0*PARAM(Kkhphos)+
     '      (10**(-7.2))*PARAM(Khphos))

      fhp=1.0d0+PARAM(k1)*PARAM(Kkhphos)
      Kapp_mghphos=Kmghphos/fhp
      kfmghphos_f=Kapp_mghphos*1.2e03
      vfmghphos=kfmghphos_f*Y(Mg)*Y(hphos)-PARAM(kfmghphos_r)*Y(Mghphos)

c----------------------MgPCr-------------------------------

      Kapp_mgPCr=PARAM(KmgPCr)
      k_fmgPCr_f=Kapp_mgPCr*1.2e03
      vfmgPCr=k_fmgPCr_f*Y(Mg)*Y(PCr)-PARAM(kfmgPCr_r)*Y(MgPCr)

c------------------------Proton Buffering-----------------------
c HATP

      Kapp_hatp=PARAM(Khatp)/fatp
      kfhatp_f=Kapp_hatp*1.2e03
      vfhatp=kfhatp_f*Y(h)*Y(ATP)-PARAM(kfhatp_r)*Y(hATP)
c HADP

      Kapp_hadp=PARAM(Khadp)/fadp
      kfhadp_f=Kapp_hadp*1.2e03
      vfhadp=kfhadp_f*Y(h)*Y(ADP)-PARAM(kfhadp_r)*Y(hADP)

c HAMP

      Kapp_hamp=PARAM(Khamp)/famp
      kfhamp_f=Kapp_hamp*1.2e03
      vfhamp=kfhamp_f*Y(h)*Y(AMP)-PARAM(kfhamp_r)*Y(hAMP)

c H2Phos

      Kapp_hphos=PARAM(Khphos)/fhp
      kfh2phos_f=Kapp_hphos*1.2e03
      vfh2phos=kfh2phos_f*Y(h)*Y(hphos)-PARAM(kfh2phos_r)*Y(h2phos)
      vfh2phose=kfh2phos_f*Y(he)*Y(hphose)-PARAM(kfh2phos_r)*Y(h2phose)

c intrinsic buffers

      vkintrinsicbuffer1=PARAM(kintrinsicbuffer1_f)*Y(B1)*Y(h)-
     '    PARAM(kintrinsicbuffer1_r)*Y(B1h)

      vkintrinsicbuffer2=PARAM(kintrinsicbuffer2_f)*Y(B2)*Y(h)-
     '    PARAM(kintrinsicbuffer2_r)*Y(B2h)

c external protein buffers

      vexternalproteinbuffer=PARAM(kexternalproteinbuffer_f)*Y(Be)*Y(He)
     '    -PARAM(kexternalproteinbuffer_r)*Y(Behe)

c Ca Buffers

      v_trop=PARAM(ktrop_f)*Y(Ca)*Y(Trop)-PARAM(ktrop_r)*Y(CaTrop)
      v_cal=PARAM(kcal_f)*Y(Ca)*Y(Cal)-PARAM(kcal_r)*Y(CaCal)


c-------------------------------time derivatives---------------------
C      DY(Vm)=0.0d0 !need to check what to do with Vm
C      DY(ADP)=-v_ak-vfhadp-vfmgadp
C      DY(AMP)=v_ak-vfhamp
C      DY(ATP)=-vfhatp-vfmgatp
C      DY(B1)=-vkintrinsicbuffer1
C      DY(B1h)=vkintrinsicbuffer1
C      DY(B2)=-vkintrinsicbuffer2
C      DY(B2h)=vkintrinsicbuffer2
C      DY(Be)=-vexternalproteinbuffer
C      DY(Behe)=vexternalproteinbuffer
C      DY(Ca)=0!-v_cal+vicab+vnce-v_trop couldn't find the rate variable check with PJM
C      DY(CaCal)=v_cal
C      DY(Cae)=0!-vicab-vnce couldn't find the rate variable check with PJM
C      DY(Cal)=-v_cal
C      DY(CaTrop)=v_trop
C      DY(CO2)=-vco2trans+v_fh2co3
C      DY(CO2e)=vco2trans+v_fh2co3e
C      DY(Cr)=v_ck
C      DY(Glyc)=-vglycl
C
C      IF(MODEL(USE_ClOSED_SYS).EQ.1) THEN
C        DY(h)=vatpase2+v_che-v_ck-v_fh2CO3-
C     '  vfh2phos-vfhadp-vfhamp-vfhatp-
C     '    vglyc2-vkintrinsicbuffer1-vkintrinsicbuffer2-
C     '    vlactrans+v_napump-vnhe-voxphos2
C        DY(hphos)=vatpase2-vfh2phos-vfmghphos-3.0d0*vglyc2+v_napump-
C     '    voxphos2
C        DY(Lac)=2.0d0*vglyc2-vlactrans
C        DY(MgADP)=-v_ak+vatpase2-v_ck+vfmgadp-(3.0d0*vglyc2)+
C     '    v_napump-voxphos2
C        DY(MgATP)=v_ak-vatpase2+v_ck+vfmgadp+(3.0d0*vglyc2)-
C     '    v_napump+voxphos2
C      ELSE !open system
C        DY(h)=vatpase1+v_che-v_ck-v_fh2CO3-
C     '    vfh2phos-vfhadp-vfhamp-vfhatp-
C     '    vglycl-vkintrinsicbuffer1-vkintrinsicbuffer2-
C     '    vlactrans+v_napump-vnhe-voxphos1
C        DY(hphos)=vatpase1-vfh2phos-vfmghphos-3.0d0*vglycl+v_napump-
C     '    voxphos1
C        DY(Lac)=2.0d0*vglycl-vlactrans
C        DY(MgADP)=-v_ak+vatpase1-v_ck+vfmgadp-(3.0d0*vglycl)+
C     '    v_napump-voxphos1
C        DY(MgATP)=v_ak-vatpase1+v_ck+vfmgadp+(3.0d0*vglycl)-
C     '    v_napump+voxphos1
C      ENDIF
C
C      DY(h2phos)=vfh2phos-vphostrans
C      DY(h2phose)=vfh2phose+vphostrans
C      DY(hADP)=vfhadp
C      DY(hAMP)=vfhamp
C      DY(hATP)=vfhatp
C      DY(hco3)=-vae-v_fh2co3+v_nbc
C      DY(hco3e)=vae-v_fh2co3e-v_nbc   !should be zero for open system ?? CHECK
C      DY(he)=-v_che-vexternalproteinbuffer-v_fh2co3e-vfh2phose+
C     '  vlactrans+vnhe
C      DY(hphose)=vfh2phose
C      DY(K)=-vke_flux+2.0d0*v_napump
C      DY(Ke)=vke_flux-2.0d0*v_napump !should be zero for open system ?? CHECK
C      DY(Lace)=vlactrans
C      DY(MCT)=-vmct
C      DY(MCTe)=vmct
C      DY(Mg)=-vfmgadp-vfmgatp-vfmghphos-vfmgPCr
C      DY(Mghphos)=vfmghphos
C      DY(MgPCr)=vfmgPCr
C      DY(Na)=vna_influx-3.0d0*v_napump+v_nbc+vnhe
C      DY(Nae)=-vna_influx+3.0d0*v_napump-v_nbc+3.0d0*vnce-vnhe
C      DY(PCr)=-v_ck-vfmgPCr
C      DY(R1)=0!-vr1 couldn't find the rate variable check with PJM
C      DY(R1A)=0!vr1 couldn't find the rate variabl check with PJMe
C      DY(R2)=-0!vr2 couldn't find the rate variable check with PJM
C      DY(R2A)=0!vr2 couldn't find the rate variable check with PJM
C      DY(Trop)=-v_trop
C      DY(Cl)=0.0d0    !could not find rate equation check with PJM
C      DY(Cle)=0.0d0   !could not find rate equation check with PJM

c---------------------------------------------------------------------------

c---------------------implimenting simple open model--------------------


      DY(Vm)=0.0d0 !need to check what to do with Vm
      DY(ADP)=-v_ak-vfhadp-vfmgadp
      DY(AMP)=v_ak-vfhamp
      DY(ATP)=-vfhatp-vfmgatp
      DY(B1)=-vkintrinsicbuffer1
      DY(B1h)=vkintrinsicbuffer1
      DY(B2)=-vkintrinsicbuffer2
      DY(B2h)=vkintrinsicbuffer2
      DY(Be)=0.0d0
      DY(Behe)=0.0d0
      DY(Ca)=0.0d0
      DY(CaCal)=0.0d0
      DY(Cae)=0.0d0
      DY(Cal)=0.0d0
      DY(CaTrop)=0.0d0
      DY(CO2)=-vco2trans+v_fh2co3
      DY(CO2e)=0.0d0
      DY(Cr)=v_ck
      DY(Glyc)=-vglycl
      DY(h)=vatpase1+v_che-v_ck-v_fh2CO3-
     '  vfh2phos-vfhadp-vfhamp-vfhatp-
     '  vglycl-vkintrinsicbuffer1-vkintrinsicbuffer2-
     '  vlactrans-vnhe-voxphos1
      DY(h2phos)=vfh2phos-vphostrans
      DY(h2phose)=0.0d0
      DY(hADP)=vfhadp
      DY(hAMP)=vfhamp
      DY(hATP)=vfhatp
      DY(hco3)=-vae-v_fh2co3+v_nbc
      DY(hco3e)=((vae-v_nbc)/2.0d0)-v_fh2co3e
      DY(he)=0.0d0
      DY(hphos)=vatpase1-vfh2phos-vfmghphos-3.0d0*vglycl-
     '  voxphos1
      DY(Lac)=2.0d0*vglycl-vlactrans
      DY(MgADP)=-v_ak+vatpase1-v_ck+vfmgadp-(3.0d0*vglycl)-
     '  voxphos1
      DY(MgATP)=v_ak-vatpase1+v_ck+vfmgatp+(3.0d0*vglycl)+
     '  voxphos1
      DY(hphose)=0.0d0
      DY(K)=0.0d0
      DY(Ke)=0.0d0
      DY(Lace)=0.0d0
      DY(MCT)=-vmct
      DY(MCTe)=vmct
      DY(Mg)=-vfmgadp-vfmgatp-vfmghphos-vfmgPCr
      DY(Mghphos)=vfmghphos
      DY(MgPCr)=vfmgPCr
      DY(Na)=0.0d0
      DY(Nae)=0.0d0
      DY(PCr)=-v_ck-vfmgPCr
      DY(R1)=0.0d0
      DY(R1A)=0.0d0
      DY(R2)=0.0d0
      DY(R2A)=0.0d0
      DY(Trop)=0.0d0
      DY(Cl)=0.0d0
      DY(Cle)=0.0d0


        DERIVED(DVOXPHOS1)=0.0D0
        DERIVED(rv_ak)=v_ak
        DERIVED(rvfhadp)=vfhadp
        DERIVED(rvfmgadp)=vfmgadp
        DERIVED(rvfhamp)=vfhamp
        DERIVED(rvfhatp)=vfhatp
        DERIVED(rvfmgatp)=vfmgatp


      RETURN
      END


C      SUBROUTINE LR_INIT_CELL(ICQS,RCQS,YQS,ERROR,*)
C
CC#### Subroutine: LR_INIT_CELL
CC###  Description:
CC###    Initialise the Luo-Rudy model
CC***  Created by David Nickerson, May 1999
C
C      IMPLICIT NONE
C
C      INCLUDE 'cmiss$reference:cell00.cmn'
C      INCLUDE 'cmiss$reference:cell02.cmn'
C      INCLUDE 'cmiss$reference:cell_lr.inc'
C      INCLUDE 'cmiss$reference:cell_reserved.inc'
C      INCLUDE 'cmiss$reference:geom00.cmn'
C
C      !parameter list
C      INTEGER ICQS(NQIM)
C      REAL*8 RCQS(NQRM),YQS(NIQSM,NQM)
C      CHARACTER ERROR*(*)
C      !local variables
C      INTEGER nd
C
C      CALL ENTERS('LR_INIT_CELL',*9999)
C
C      !initialise stuff
Cc      RCQS(CELL_PARAMETERS_OFFSET-1+Cai_on)=0.0d0
Cc      RCQS(CELL_PARAMETERS_OFFSET-1+Cai2)=0.0d0
Cc      ICQS(CELL_CONTROL_OFFSET-1+NUM_APPLIED_STIMULII)=0
C      ICQS(1)=ICQS(1)
C      RCQS(1)=RCQS(1)
C      DO nd=0,CELL_NUM_DERIVED-1
C        YQS(CELL_DERIVED_OFFSET+nd,1)=0.0d0
C      ENDDO
C
C      CALL EXITS('LR_INIT_CELL')
C      RETURN
C 9999 CALL ERRORS('LR_INIT_CELL',ERROR)
C      CALL EXITS('LR_INIT_CELL')
C      RETURN 1
C      END


C *** DPN 14 March 2000 - fixing up units
C      SUBROUTINE NOBLE98_CELL(T,Y,DY,CONTROL,MODEL,SIZES,VARIANT,
C     '  DERIVED,PARAMETERS,PROTOCOL,AII,AIO,ARI,ARO,ERR_CODE)
C
CC#### Subroutine: NOBLE98_CELL
CC###  Description:
CC###    Calculates the RHS of the user defined cell 1 odes.
C
C      IMPLICIT NONE
C
C      INCLUDE 'cmiss$reference:cell_reserved.inc'
C      INCLUDE 'cmiss$reference:cell_n98.inc'
C
C      INTEGER SIZES(10)
C      INTEGER CONTROL(*),MODEL(*),VARIANT,AII(*),AIO(*),ERR_CODE
C      REAL*8 T,Y(*),DY(*),DERIVED(*),PARAM(*),PROTOCOL(*),
C     '  ARI(*),ARO(*)
C
C      !local variables
C      REAL*8 IK1, !Time-independent (background) K+ current (nA.mm^-2)
C     '  Ito, !Transient outward K+ current (nA.mm^-2)
C     '  IKr1, !Time-dependent (delayed) fast K+ current (nA.mm^-2)
C     '  IKr2, !Time-dependent (delayed) fast K+ current (nA.mm^-2)
C     '  IKs, !Time-dependent (delayed) slow K+ current (nA.mm^-2)
C     '  IKNa, !Na+-dependent K+ current (nA.mm^-2)
C     '  IbK, !Background K+ current (nA.mm^-2)
C     '  IKATP, !ATP-dependent K+ current (nA.mm^-2)
C     '  IKACh,  !ACh-dependent K+ current (nA.mm^-2)
C     '  INa, !Fast Na+ current (nA.mm^-2)
C     '  IbNa, !Background Na+ current (nA.mm^-2)
C     '  IpNa, !Persistent Na+ current (nA.mm^-2)
C     '  ICaLK, !L-type K+ current (nA.mm^-2)
C     '  ICaLNa, !L-type Na+ current (nA.mm^-2)
C     '  ICaLCa, !L-type Ca++ current (nA.mm^-2)
C     '  ICaLKDS, !Diadic Space L-type K+ current (nA.mm^-2)
C     '  ICaLNaDS, !Diadic Space L-type Na+ current (nA.mm^-2)
C     '  ICaLCaDS, !Diadic Space L-type Ca++ current (nA.mm^-2)
C     '  IbCa, !Background Ca++ current (nA.mm^-2)
C     '  INaK, !Na+/K+ exchange current (nA.mm^-2)
C     '  INaCa, !Na+/Ca++ exchange current (nA.mm^-2)
C     '  INaCaDS, !Na+/Ca++ Diadic space exchange current (nA.mm^-2)
C     '  IK_stretch, !K+ stretch activated current (nA.mm^-2)
C     '  INa_stretch, !Na+ stretch activated current (nA.mm^-2)
C     '  ICa_stretch, !Ca++ stretch activated current (nA.mm^-2)
C     '  INs_stretch, !Ns stretch activated current (nA.mm^-2)
C     '  IAn_stretch, !An stretch activated current (nA.mm^-2)
C     '  Istim, !Stimulus current (nA.mm^-2)
C     '  IK_tot, !Total K+ current (nA.mm^-2)
C     '  INa_tot, !Total Na+ current (nA.mm^-2)
C     '  ICa_tot, !Total Ca+ current (nA.mm^-2)
C     '  ICaDS_tot !Total DS Ca+ current (nA.mm^-2)
C      REAL*8 jup, !SR Ca++ uptake flux (mM.L^-1.ms^-1)
C     '  jtr, !SR Ca++ translocation flux (mM.L^-1.ms^-1)
C     '  jrel, !SR Ca++ release flux (mM.L^-1.ms^-1)
C     '  jleak, !SR Ca++ leakage flux (mM.L^-1.ms^-1)
C     '  jdecay !DS Ca++ decay flux (mM.L^-1.ms^-1)
C      REAL*8 RTF, !RT/F (mV)
C     '  ENa, !Na reversal potential (mV)
C     '  EK, !K reversal potential (mV)
C     '  EKs, !K (IKs) reversal potential (mV)
C     '  ECa, !Ca reversal potential (mV)
C     '  Emh !Na (INa) reversal potential (mV)
C      REAL*8 alpha_xr1, !IKr x1 alpha rate (ms^-1)
C     '  beta_xr1, !IKr x1 beta rate (ms^-1)
C     '  alpha_xr2, !IKr x2 alpha rate (ms^-1)
C     '  beta_xr2, !IKr x2 beta rate (ms^-1)
C     '  alpha_xs, !IKs x alpha rate (ms^-1)
C     '  beta_xs, !IKs x beta rate (ms^-1)
C     '  alpha_sto, !Ito s alpha rate (ms^-1)
C     '  beta_sto, !Ito s beta rate (ms^-1)
C     '  alpha_mNa, !INa m alpha rate (ms^-1)
C     '  beta_mNa, !INa m beta rate (ms^-1)
C     '  alpha_hNa, !INa h alpha rate (ms^-1)
C     '  beta_hNa, !INa h beta rate (ms^-1)
C     '  alpha_dCa, !ILCa d alpha rate (ms^-1)
C     '  beta_dCa, !ILCa d beta rate (ms^-1)
C     '  alpha_fCa, !ILCa f alpha rate (ms^-1)
C     '  beta_fCa, !ILCa f beta rate (ms^-1)
C     '  alpha_xACh1, !IKACh x1 alpha rate (ms^-1)
C     '  beta_xACh1, !IKACh x1 beta rate (ms^-1)
C     '  alpha_xACh2, !IKACh x2 alpha rate (ms^-1)
C     '  beta_xACh2, !IKACh x2 beta rate (ms^-1)
C     '  alpha_act, !ICaL activation rate (ms^-1)
C     '  alpha_inact !ICaL inactivation rate (ms^-1)
C      REAL*8 c1,c2,FractBackSRSites,FractCaUptakeSites,f_stretch,
C     '  OpenReleaseChannelFract,PrecursorFract,RegulatoryBindingSite,
C     '  VRTF,z1,z2
C
CC     Physical constants
C      REAL*8 R,F,nNaK,nNaCa
C      PARAMETER(R=8314.41d0) !Gas constant (mJ.Mol^-1.K^-1)
C      PARAMETER(F=96485.0d0) !Faradays constant (C.Mol^-1)
C      PARAMETER(nNaK=1.5d0) !Na+/K+ exchange stoichometry ratio ()
C      PARAMETER(nNaCa=3.0d0) !Na+/Ca++ exchange stoichometry ratio ()
C
C      IF(SIZES(5).LT.117) THEN
C        ERR_CODE=1
C      ELSE IF(SIZES(3).LT.1) THEN
C        ERR_CODE=2
C      ELSE IF(SIZES(1).LT.1) THEN
C        ERR_CODE=3
C      ELSE IF(SIZES(1).LT.25) THEN
C        ERR_CODE=4
C      ELSE IF(SIZES(6).LT.10) THEN
C        ERR_CODE=5
C      ELSE
C        IF(CONTROL(RETURN_CURRENT).NE.0.AND.SIZES(4).LT.37) THEN
C          ERR_CODE=6
C        ELSE
C          ERR_CODE=0
C
C          Istim=PROTOCOL(PSTIMULUS)
C          IF(T.GE.PROTOCOL(STIMULUS1_ON).AND.
C     '      T.LT.PROTOCOL(STIMULUS1_OFF))
C     '      Istim=Istim+PROTOCOL(STIMULUS1_MAG)
C          IF(T.GE.PROTOCOL(STIMULUS2_ON).AND.
C     '      T.LT.PROTOCOL(STIMULUS2_OFF))
C     '      Istim=Istim+PROTOCOL(STIMULUS2_MAG)
C          IF(PROTOCOL(STIMULUS_FREQ_PERIOD).GT.1.0d-6) THEN
C            IF(DMOD(T,PROTOCOL(STIMULUS_FREQ_PERIOD)).LT.
C     '        PROTOCOL(STIMULUS_FREQ_DURATION))
C     '        Istim=Istim+PROTOCOL(STIMULUS_FREQ_MAG)
C          ENDIF
C          Istim=1000.0d0*Istim/PARAM(Am) !Convert to nA.mm^-2
C
C          RTF = R*PARAM(TEMP)/F
C
CC         Equilibrium potentials
C          ENa = RTF*DLOG(PARAM(Nao)/Y(Nai)) !(mV) Equation 1
C          EK = RTF*DLOG(PARAM(Ko)/Y(Ki)) !(mV) Equation 2
C          EKs = RTF*DLOG((PARAM(Ko)+PARAM(PKNa)*
C     '      PARAM(Nao))/(Y(Ki)+PARAM(PKNa)*Y(Nai))) !(mV) Equation 3
C          ECa = 0.5d0*RTF*DLOG(PARAM(Cao)/Y(Cai)) !(mV) Equation 4
C          Emh = RTF*DLOG((PARAM(Nao)+PARAM(PNaK)*
C     '      PARAM(Ko))/(Y(Nai)+PARAM(PNaK)*Y(Ki))) !(mV) Equation 5
C
CC         Time-independent (background) K+ current, IK1
C          IK1 = PARAM(gK1)*PARAM(Ko)/(PARAM(Ko)+
C     '      PARAM(KmK1))*(Y(Vm)-EK)/(1.0d0+
C     '      DEXP(PARAM(STEEP_IK1)*(Y(Vm)-EK+10.0d0-
C     '      PARAM(SHIFT_IK1))/RTF)) !(nA.mm^-2) Equation 6
C
CC         Transient outward current, Ito
C          alpha_sto = PARAM(epsilon_sto)*0.000033d0*
C     '      DEXP(-(Y(Vm)-PARAM(SHIFT_sto))/
C     '      (2.125d0*PARAM(STEEP_sto))) !(ms^-1) Equation 10
C          beta_sto =  PARAM(epsilon_sto)*0.033d0/
C     '      (1.0d0+DEXP(-(Y(Vm)+10.0d0-PARAM(SHIFT_sto))/
C     '      PARAM(STEEP_sto))) !(ms^-1) Equation 11
C          Ito = PARAM(gto)*(PARAM(gtos)+Y(sto)*
C     '      (1.0d0-PARAM(gtos)))*Y(rto)*(Y(Vm)-EK) !(nA.mm^-2) Equation 7
C          DY(rto) = PARAM(epsilon_rto)*PARAM(alpha_rto)*
C     '      (1.0d0/(1.0d0+DEXP(-(Y(Vm)+4.0d0-PARAM(SHIFT_rto))/
C     '      PARAM(STEEP_rto)))-Y(rto)) !(ms^-1) Equation 8
C          DY(sto) = alpha_sto*(1.0d0-Y(sto))-beta_sto*Y(sto) !(ms^-1) Equation 9
C
CC         Time-dependent (delayed) rapid K+ current, IKr1
C          alpha_xr1 = PARAM(epsilon_alpha_xr1)*0.05d0/
C     '      (1.0d0+DEXP(-(Y(Vm)-5.0d0-PARAM(SHIFT_alpha_xr1))/
C     '      9.0d0)) !(ms^-1) Equation 14
C          beta_xr1 = PARAM(epsilon_beta_xr1)*0.00005d0*
C     '      DEXP(-(Y(Vm)-20.0d0-PARAM(SHIFT_beta_xr1))/
C     '      15.0d0) !(ms^-1) Equation 15
C          IKr1 = PARAM(gKr1)*Y(xr1)*(Y(Vm)-EK)/(1.0d0+
C     '      DEXP((Y(Vm)+9.0d0-PARAM(SHIFT_IKr))/
C     '      22.4d0)) !(nA.mm^-2) Equation 12
C          DY(xr1) = alpha_xr1*(1.0d0-Y(xr1))-beta_xr1*Y(xr1) !(ms^-1) Equation 13
C
CC         Time-dependent (delayed) rapid K+ current, IKr2
C          alpha_xr2 = PARAM(epsilon_alpha_xr2)*0.05d0/
C     '      (1.0d0+DEXP(-(Y(Vm)-5.0d0-PARAM(SHIFT_alpha_xr2))/
C     '      9.0d0)) !(ms^-1) Equation 18
C          beta_xr2 = PARAM(epsilon_beta_xr2)*0.0004d0*
C     '      DEXP(-((Y(Vm)+30.0d0-PARAM(SHIFT_beta_xr2))/
C     '      30.0d0)**3) !(ms^-1) Equation 19
C          IKr2 = PARAM(gKr2)*Y(xr2)*(Y(Vm)-EK)/(1.0d0+
C     '      DEXP((Y(Vm)+9.0d0+PARAM(SHIFT_IKr))/
C     '      22.4d0)) !(nA.mm^-2) Equation 16
C          DY(xr2) = alpha_xr2*(1.0d0-Y(xr2))-beta_xr2*Y(xr2) !(ms^-1) Equation 17
C
CC         Time-dependent (delayed) slow K+ current, IKs
C          alpha_xs = PARAM(epsilon_alpha_xs)*0.014d0/
C     '      (1.0d0+DEXP(-(Y(Vm)-40.0d0-PARAM(SHIFT_alpha_xs))/
C     '      9.0d0)) !(ms^-1) Equation 22
C          beta_xs = PARAM(epsilon_beta_xs)*0.001d0*
C     '      DEXP(-(Y(Vm)-PARAM(SHIFT_beta_xs))/45.0d0) !(ms^-1) Equation 23
C          IKs = PARAM(gKs)*Y(xs)**2*(Y(Vm)-EKs) !(nA.mm^-2) Equation 20
C          DY(xs) = alpha_xs*(1.0d0-Y(xs))-beta_xs*Y(xs) !(ms^-1) Equation 21
C
CC         Na+ dependent K+ current, IKNa
C          IKNa = PARAM(gKNa)*Y(Nai)/(Y(Nai)+PARAM(KmgKNa))*
C     '      (Y(Vm)-EK) !(nA.mm^-2) Equation 24
C
CC         Background Potassium Current, IbK
C          IbK = PARAM(gbK)*(Y(Vm)-EK) !(nA.mm^-2) Equation 25
C
CC         ATP-dependent K+ current, IKATP
C          IKATP = PARAM(gKATP)*(Y(Vm)-PARAM(EKATPrev))/
C     '      (1.0d0+(PARAM(ATP)/PARAM(KmATP))**2) !(nA.mm^-2) Equation 26
C
CC         ACh-dependent K+ current (Mark Boyett Formulation), IKACh
C          alpha_xACh1 = 0.003684211d0 !(ms^-1) Equation 30
C          beta_xACh1 = 0.00582d0/(1.0d0+DEXP(-(Y(Vm)+
C     '      PARAM(SHIFT_ACh1))/15.0d0)) !(ms^-1) Equation 31
C          alpha_xACh2 = 0.07309924d0 !(ms^-1) Equation 32
C          beta_xACh2 = 0.12d0/(1.0d0+DEXP(-(Y(Vm)+
C     '      PARAM(SHIFT_ACh2))/15.0d0)) !(ms^-1) Equation 33
C          IKACh = PARAM(gKACh)*(PARAM(ACh)**1.4969d0/
C     '      (PARAM(ACh)**1.4969d0+PARAM(KmACh)**1.4969d0))*
C     '      (PARAM(Ko)/(PARAM(Ko)+PARAM(KmKACh)))*
C     '      Y(xACh1)*Y(xACh2)*(Y(Vm)-EK)/(1.0d0+DEXP(0.4d0*(Y(Vm)-
C     '      EK-140.0d0)/RTF)) !(nA.mm^-2) Equation 27
C          DY(xACh1) = alpha_xACh1*(1.0d0-Y(xACh1))-
C     '      beta_xACh1*Y(xACh1) !(ms^-1) Equation 28
C          DY(xACh2) = alpha_xACh2*(1.0d0-Y(xACh2))-
C     '      beta_xACh2*Y(xACh2) !(ms^-1) Equation 29
C
CC         Fast Na+ current, INa
C          IF(DABS(Y(Vm)+41.0d0-PARAM(SHIFT_mNa)).LT.0.00001d0) THEN
C            alpha_mNa = 2.0d0 !(ms^-1) Equation 37
C          ELSE
C            alpha_mNa = 0.2d0*(Y(Vm)+41.d0-PARAM(SHIFT_mNa))/
C     '        (1.0d0-DEXP(-(Y(Vm)+41.0d0-PARAM(SHIFT_mNa))/
C     '        10.0d0)) !(ms^-1) Equation 37
C          ENDIF
C          beta_mNa = 8.0d0*DEXP(-(Y(Vm)+66.d0-PARAM(SHIFT_mNa))/
C     '      18.0d0) !(ms^-1) Equation 38
C          alpha_hNa = 0.02d0*DEXP(-(Y(Vm)+75.d0-
C     '      PARAM(SHIFT_hNa))/8.0d0) !(ms^-1) Equation 39
C          beta_hNa = 2.0d0/(1.d0+320.0d0*DEXP(-(Y(Vm)+
C     '      75.0d0-PARAM(SHIFT_hNa))/10.0d0)) !(ms^-1) Equation 40
C          INa = PARAM(gNa)*Y(mNa)**3*Y(hNa)*(Y(Vm)-Emh) !(nA.mm^-2) Equation 34
C          DY(mNa) = alpha_mNa*(1.0d0-Y(mNa))-beta_mNa*Y(mNa) !(ms^-1) Equation 35
C
C          DY(hNa) = alpha_hNa*(1.0d0-Y(hNa))-beta_hNa*Y(hNa) !(ms^-1) Equation 36
C
CC         Background Na+ current, IbNa
C          IbNa = PARAM(gbNa)*(Y(Vm)-ENa) !(nA.mm^-2) Equation 41
C
CC         Persistent Na+ current, IpNa
C          IpNa = PARAM(gpNa)*(Y(Vm)-ENa)/(1.0d0+DEXP(-(Y(Vm)+
C     '      52.0d0)/8.0d0)) !(nA.mm^-2) Equation 42
C
CC         Background Ca++ current, IbCa
C          IbCa  = PARAM(gbCa)*(Y(Vm)-ECa) !(nA.mm^-2) Equation 43
C
CC         L-type Ca++ channel, ICaL
C          IF(DABS(Y(Vm)+24.0d0-PARAM(SHIFT_dCa)).LT.0.001d0) THEN
C            alpha_dCa = 0.12d0*PARAM(epsilon_d) !(ms^-1) Equation 54
C            beta_dCa = 0.12d0*PARAM(epsilon_d) !(ms^-1) Equation 55
C          ELSE
C            alpha_dCa = 0.03d0*PARAM(epsilon_d)*(Y(Vm)+24.0d0-
C     '        PARAM(SHIFT_dCa))/(1.0d0-DEXP(-(Y(Vm)+24.d0-
C     '        PARAM(SHIFT_dCa))/
C     '        PARAM(STEEP_dCa))) !(ms^-1) Equation 54
C            beta_dCa = -0.012d0*PARAM(epsilon_d)*
C     '        (Y(Vm)+24.0d0-PARAM(SHIFT_dCa))/
C     '        (1.0d0-DEXP((Y(Vm)+24.0d0-PARAM(SHIFT_dCa))/
C     '        (2.5d0*PARAM(STEEP_dCa)))) !(ms^-1) Equation 54
C          ENDIF
C          IF(DABS(Y(Vm)+34.0d0-PARAM(SHIFT_fCa)).LT.0.001d0) THEN
C            alpha_fCa = 0.00625d0*PARAM(epsilon_f)*
C     '        PARAM(STEEP_fCa) !(ms^-1) Equation 56
C          ELSE
C            alpha_fCa = -0.00625d0*PARAM(epsilon_f)*(Y(Vm)+
C     '        34.0d0-PARAM(SHIFT_fCa))/(1.0d0-DEXP((Y(Vm)+
C     '        34.0d0-PARAM(SHIFT_fCa))/
C     '        PARAM(STEEP_fCa))) !(ms^-1) Equation 56
C          ENDIF
C          beta_fCa = 0.012d0*PARAM(epsilon_f)/(1.0d0+
C     '      DEXP(-(Y(Vm)+34.0d0-PARAM(SHIFT_fCa))/
C     '      PARAM(STEEP_fCa))) !(ms^-1) Equation 57
C
C          z1 = PARAM(PCa)*(1.0d0-PARAM(ACh)/
C     '      (PARAM(ACh)+PARAM(KmAChICaL)))
C          z2 = DEXP(PARAM(Esurf)/RTF)*Y(Ki)-DEXP(-(Y(Vm)-
C     '      PARAM(Esurf))/RTF)*PARAM(Ko)
C          ICaLK = (1.0d0-PARAM(ICaLFract))*PARAM(PCaK)*
C     '      Y(dCa)*Y(fCa)*Y(f2Ca)*z1*z2 !(nA.mm^-2) Equation 44
C          ICaLKDS = PARAM(ICaLFract)*PARAM(PCaK)*Y(dCa)*
C     '      Y(fCa)*Y(f2CaDS)*z1*z2 !(nA.mm^-2) Equation 47
C          z2 = DEXP(PARAM(Esurf)/RTF)*Y(Nai)-DEXP(-(Y(Vm)-
C     '      PARAM(Esurf))/RTF)*PARAM(Nao)
C
C          ICaLNa = (1.0d0-PARAM(ICaLFract))*PARAM(PCaNa)*
C     '      Y(dCa)*Y(fCa)*Y(f2Ca)*z1*z2 !(nA.mm^-2) Equation 45
C          ICaLNaDS = PARAM(ICaLFract)*PARAM(PCaNa)*Y(dCa)*
C     '      Y(fCa)*Y(f2CaDS)*z1*z2 !(nA.mm^-2) Equation 48
C          z1 = PARAM(PCa)*(1.0d0-PARAM(ACh)/
C     '      (PARAM(ACh)+PARAM(KmAChICaL)))*
C     '      (Y(Vm)-PARAM(Esurf))/(RTF*(1.0d0-
C     '      DEXP(-2.0d0*(Y(Vm)-PARAM(ESurf))/RTF)))
C          z2 = DEXP(2.0d0*PARAM(Esurf)/RTF)*Y(Cai)-
C     '      DEXP(-2.0d0*(Y(Vm)-PARAM(Esurf))/RTF)*PARAM(Cao)
C          ICaLCa = (1.0d0-PARAM(ICaLFract))*4.0d0*
C     '      Y(dCa)*Y(fCa)*Y(f2Ca)*z1*z2 !(nA.mm^-2) Equation 46
C          z2 = DEXP(2.0d0*PARAM(Esurf)/RTF)*Y(CaDS)-
C     '      DEXP(-2.0d0*(Y(Vm)-PARAM(Esurf))/RTF)*PARAM(Cao)
C          ICaLCaDS = PARAM(ICaLFract)*4.0d0*
C     '      Y(dCa)*Y(fCa)*Y(f2CaDS)*z1*z2 !(nA.mm^-2) Equation 49
C
C          DY(dCa) = alpha_dCa*(1.0d0-Y(dCa))-beta_dCa*Y(dCa) !(ms^-1) Equation 50
C
C          DY(fCa) = alpha_fCa*(1.0d0-Y(fCa))-beta_fCa*Y(fCa) !(ms^-1) Equation 51
C
C          DY(f2Ca) = PARAM(alpha_CaInact)*(1.0d0-
C     '      Y(Cai)/(Y(Cai)+PARAM(KmCaInact))-Y(f2Ca)) !(ms^-1) Equation 52
C
C          DY(f2CaDS) = PARAM(alpha_CaDSInact)*(1.0d0-
C     '      Y(CaDS)/(Y(CaDS)+PARAM(KmCaDSInact))-
C     '      Y(f2CaDS)) !(ms^-1) Equation 53
C
CC         Na+/K+ Exchange, INaK (Ip)
C          INaK = PARAM(INaKmax)*PARAM(Ko)/
C     '      (PARAM(Ko)+PARAM(KmK))*
C     '      Y(Nai)/(Y(Nai)+PARAM(KmNa)) !(nA.mm^-2) Equation 58
C
CC         Na+/Ca++ exchange current, INaCa
C          VRTF = Y(Vm)/RTF
C          c1 = DEXP(PARAM(gamma_NaCa)*(nNaCa-2.0d0)*VRTF)
C          c2 = DEXP((PARAM(gamma_NaCa)-1.0d0)*(nNaCa-2.0d0)*VRTF)
C          z1 = c1*Y(Nai)**nNaCa*PARAM(Cao)-
C     '      c2*PARAM(Nao)**nNaCa*Y(Cai)
C          z2 = 1.0d0+PARAM(dNaCa)*(PARAM(Nao)**nNaCa*
C     '      Y(Cai)+Y(Nai)**nNaCa*PARAM(Cao))
C          INaCa = (1.0d0-PARAM(INaCaFract))*PARAM(INaCamax)*
C     '      (z1/z2)*1.0d0/(1.0d0+Y(Cai)/0.0069d0) !(nA.mm^-2) Equation 59
C
CC         Diadic-space Na+/Ca++ exchange current, INaCaDS
C          z1 = c1*Y(Nai)**nNaCa*PARAM(Cao)-
C     '      c2*PARAM(Nao)**nNaCa*Y(CaDS)
C          z2 = 1.0d0+PARAM(dNaCa)*(PARAM(Nao)**nNaCa*
C     '      Y(CaDS)+Y(Nai)**nNaCa*PARAM(Cao))
C          INaCaDS = PARAM(INaCaFract)*PARAM(INaCamax)*
C     '      (z1/z2)*1.0d0/(1.0d0+Y(CaDS)/0.0069d0) !(nA.mm^-2) Equation 60
C
C          IF(MODEL(USE_SAC).EQ.1) THEN
CC           Stretch activated channels
C            f_stretch = 1.0d0/(1.0d0+DEXP(-2.0d0*
C     '        PARAM(gamma_SACSL)*(PARAM(SL)-
C     '        PARAM(SLRef)))) !() Equation 61
C
CC           K+ stretch activated channel
C            IK_stretch = PARAM(gK_stretch)*f_stretch*
C     '        (Y(Vm)-EK) !(nA.mm^-2) Equation 62
C
CC           Na+ stretch activated channel
C            INa_stretch = PARAM(gNa_stretch)*f_stretch*
C     '        (Y(Vm)-ENa) !(nA.mm^-2) Equation 63
C
CC           Ca++ stret ch activated channel
C            ICa_stretch = PARAM(gCa_stretch)*f_stretch*
C     '        (Y(Vm)-ECa) !(nA.mm^-2) Equation 64
C
CC           INs stretch activated channel
C            INs_stretch = PARAM(gNs_stretch)*f_stretch*
C     '        (Y(Vm)-PARAM(ENs_stretch)) !(nA.mm^-2) Equation 65
C
CC           IAn stretch activated channel
C            IAn_stretch = PARAM(gAn_stretch)*f_stretch*
C     '        (Y(Vm)-PARAM(EAn_stretch)) !(nA.mm^-2) Equation 66
C          ELSE
C            INa_stretch = 0.0d0
C            IK_stretch = 0.0d0
C            ICa_stretch = 0.0d0
C            INs_stretch = 0.0d0
C            IAn_stretch = 0.0d0
C          ENDIF
C
CC         SR Ca++ uptake flux, jup
C          c1 = PARAM(KcyCa)*PARAM(Kxcs)/PARAM(KSRCa)
C          c2 = Y(Cai)+Y(Caup)*c1+PARAM(KcyCa)*PARAM(Kxcs)+
C     '      PARAM(KcyCa)
C          FractCaUptakeSites = Y(Cai)/c2 !() Equation 68
C          FractBackSRSites = Y(Caup)*c1/c2 !() Equation 69
C          jup = PARAM(alpha_up)*FractCaUptakeSites-
C     '      PARAM(beta_up)*FractBackSRSites !(mM.L^-1.ms^-1) Equation 67
C
CC         SR Ca++ translocation flux, jtr
C          jtr = PARAM(alpha_tr)*(Y(Caup)-Y(Carel)) !(mM.L^-1.ms^-1) Equation 70
C
CC         SR Ca++ leakage flux, jleak
C          jleak = DEXP(PARAM(gamma_SRSL)*PARAM(SL))*
C     '      PARAM(alpha_SRLeak)*Y(Carel) !(mM.L^-1.ms^-1) Equation 71
C
CC         SR Ca++ release flux, jrel
C          OpenReleaseChannelFract=(Y(ActivatorFract)/
C     '      (Y(ActivatorFract)+0.25d0))**2 !() Equation 72
C          jrel = PARAM(alpha_rel)*OpenReleaseChannelFract*
C     '      Y(Carel) !(mM.L^-1.ms^-1) Equation 71
C          RegulatoryBindingSite=(Y(Cai)/(Y(Cai)+
C     '      PARAM(KmCa))+(1.0d0-Y(Cai)/(Y(Cai)+
C     '      PARAM(KmCa)))*Y(CaDS)/(Y(CaDS)+
C     '      PARAM(KmCaDS)))**2 !() Equation 78
C          PrecursorFract=1.0d0-Y(ActivatorFract)-Y(ProductFract) !() Equation 75
C          alpha_act = 0.5d0*RegulatoryBindingSite !(ms^-1) Equation 76
C          alpha_inact = 0.06d0+0.5d0*RegulatoryBindingSite !(ms^-1) Equation 77
C          IF(Y(Vm).LT.-50.0d0) THEN
C            DY(ActivatorFract) = PARAM(epsilon_sr)*
C     '        (alpha_act*PrecursorFract-alpha_inact*Y(ActivatorFract)) !(ms^-1) Equation 73
C            DY(ProductFract) = PARAM(epsilon_sr)*
C     '        (alpha_inact*Y(ActivatorFract)-0.001d0*Y(ProductFract)) !(ms^-1) Equation 74
C          ELSE
C            DY(ActivatorFract) = alpha_act*PrecursorFract-
C     '        alpha_inact*Y(ActivatorFract) !(ms^-1) Equation 73
C            DY(ProductFract) = alpha_inact*Y(ActivatorFract)-
C     '        0.001d0*Y(ProductFract) !(ms^-1) Equation 74
C          ENDIF
C
CC         DS Ca++ decay flux, jdecay
C          jdecay = PARAM(alpha_CaDSdecay)*Y(CaDS)
C
CC         Membrane potential change
CC         Units are nA.mm^-2/uF.mm^-2 or mV.s^-1 therefore divide by
CC         1000 to get mV.ms^-1
C          DY(Vm) = (IStim-IK1-Ito-IKr1-IKr2-IKs-IKNa-IbK-IKATP-
C     '      IKACh-INa-IbNa-IpNa-ICaLK-ICaLNa-ICaLCa-
C     '      ICaLKDS-ICaLNaDS-ICaLCaDS-IbCa-INaK-INaCa-INaCaDS-
C     '      IK_stretch-INa_stretch-ICa_stretch-INs_stretch-
C     '      IAn_stretch)/(1000.0d0*PARAM(Cm)) !(mV.ms^-1) Equation 79
C
CC         Concentration changes
CC         Units are nA.mm^-2.mm^-1/C.Mol^-1 or mMol.L^-1.s^-1 therefore
CC         divide by 1000 to get mMol.L^-1.ms^-1
CC         Intracellular K+, Ki
C          DY(Ki) = -PARAM(Am)/(1000.0d0*PARAM(Vi)*F)*
C     '      (IK1+Ito+IKr1+IKr2+IKs+IKNa+IbK+IKATP+IKACh+ICaLK+
C     '      ICaLKDS-1.0d0/(nNaK-1.0d0)*INaK+IK_stretch) !(mM.L^-1.ms^-1) Equation 80
C
CC         Intracellular Na+, Nai
C          DY(Nai) = -PARAM(Am)/(1000.0d0*PARAM(Vi)*F)*
C     '      (INa+IbNa*PARAM(Nao)/140.0d0+IpNa+ICaLNa+ICaLNaDS+
C     '      nNaK/(nNaK-1.0d0)*INaK+
C     '      nNaCa/(nNaCa-2.0d0)*(INaCa+INaCaDS)) !(mM.L^-1.ms^-1) Equation 81
C
CC         Ca++ bound to calmodulin, CaCALM
C          DY(CaCALM) = PARAM(alpha_CALM)*(PARAM(CALM)-
C     '      Y(CaCALM))*Y(Cai)-PARAM(beta_CALM)*
C     '      Y(CaCALM) !(mM.L^-1.ms^-1) Equation 86
C
CC         Ca++ bound to troponin, CaTROP
C          DY(CaTROP) = PARAM(alpha_TROP)*(PARAM(TROP)-
C     '      Y(CaTROP))*Y(Cai)-PARAM(beta_TROP)*
C     '      Y(CaTROP) !(mM.L^-1.ms^-1) Equation 87
C
CC         Ca++ bound to indicator, CaIND
C          DY(CaIND) = PARAM(alpha_IND)*(PARAM(IND)-
C     '      Y(CaIND))*Y(Cai)-PARAM(beta_IND)*
C     '      Y(CaIND) !(mM.L^-1.ms^-1) Equation 88
C
CC         Intracellular Ca++, Cai
C          DY(Cai) = -PARAM(Am)/(2000.0d0*PARAM(Vi)*F)*
C     '      (ICaLCa+IbCa-2.0d0/(nNaCa-2.0d0)*INaCa+ICa_stretch)+
C     '      jrel*PARAM(Vrel)/PARAM(Vi)+
C     '      jleak*PARAM(Vrel)/PARAM(Vi)-jup+
C     '      jdecay*PARAM(VDS)/PARAM(Vi)-DY(CaCALM)-
C     '      DY(CaTROP)-DY(CaIND) !(mM.L^-1.ms^-1) Equation 82
C
CC         Diadic space Ca++, CaDS
C          DY(CaDS) = -PARAM(Am)/(2000.0d0*PARAM(VDS)*F)*
C     '      (ICaLCaDS-2.0d0/(nNaCa-2.0d0)*INaCaDS)-jdecay !(mM.L^-1.ms^-1) Equation 83
C
CC         SR uptake Ca++, Caup
C          DY(Caup) = jup*PARAM(Vi)/PARAM(Vup)-
C     '      jtr !(mM.L^-1.ms^-1) Equation 84
C
CC         SR release Ca++, Carel
C          DY(Carel) = jtr*PARAM(Vup)/PARAM(Vrel)-
C     '      jrel-jleak !(mM.L^-1.ms^-1) Equation 85
C
C          IF(CONTROL(RETURN_CURRENT).NE.0) THEN
C             DERIVED(DIK1) = IK1
C             DERIVED(DIto) = Ito
C             DERIVED(DIKr1) = IKr1
C             DERIVED(DIKr2) = IKr2
C             DERIVED(DIKs) = IKs
C             DERIVED(DIKNa) = IKNa
C             DERIVED(DIbK) = IbK
C             DERIVED(DIKATP) = IKATP
C             DERIVED(DIKACh) = IKACh
C             DERIVED(DINa) = INa
C             DERIVED(DIbNa) = IbNa
C             DERIVED(DIpNa) = IpNa
C             DERIVED(DICaLK) = ICaLK
C             DERIVED(DICaLNa) = ICaLNa
C             DERIVED(DICaLCa) = ICaLCa
C             DERIVED(DICaLKDS) = ICaLKDS
C             DERIVED(DICaLNaDS) = ICaLNaDS
C             DERIVED(DICaLCaDS) = ICaLCaDS
C             DERIVED(DIbCa) = IbCa
C             DERIVED(DINaK) = INaK
C             DERIVED(DINaCa) = INaCa
C             DERIVED(DINaCaDS) = INaCaDS
C             DERIVED(DIK_stretch) = IK_stretch
C             DERIVED(DINa_stretch) = INa_stretch
C             DERIVED(DICa_stretch) = ICa_stretch
C             DERIVED(DINs_stretch) = INs_stretch
C             DERIVED(DIAn_stretch) = IAn_stretch
C             DERIVED(Djup) = jup
C             DERIVED(Djtr) = jtr
C             DERIVED(Djrel) = jrel
C             DERIVED(Djleak) = jleak
C             DERIVED(Djdecay) = jdecay
C             DERIVED(DIK) = IKr1+IKr2+IKs
C             IK_tot = IK1+Ito+IKr1+IKr2+IKs+IKNa+IbK+IKATP+IKACh+ICaLK+
C     '         ICaLKDS+INaK+IK_stretch
C             DERIVED(DIK_tot) = IK_tot
C             INa_tot = INa+IbNa+IpNa+ICaLNa+ICaLNaDS+INaK+INaCa+INaCaDS
C             DERIVED(DINa_tot) = INa_tot
C             ICa_tot = ICaLCa+IbCa+INaCa+ICa_stretch
C             DERIVED(DICa_tot) = ICa_tot
C             ICaDS_tot = ICaLCaDS+INaCaDS
C             DERIVED(DICaDS_tot) = ICaDS_tot
C         ENDIF
C        ENDIF
C      ENDIF
C
C      RETURN
C      END
