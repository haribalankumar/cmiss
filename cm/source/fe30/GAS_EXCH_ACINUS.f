      SUBROUTINE GAS_EXCH_ACINUS(NBJ,NEELEM,NENP,nh,NPNE,nr,
     &  NVJE,NXI,NYNP,ACINUS,BBM,CE,T,XAB,XP,YP,ERROR,*)

C#### Subroutine: GAS_EXCH_ACINUS
C###  Description:
C###     Calculates the partial pressures (mmHg) for a BBM 
C###     ("black box model") which is a lumped acinar unit
C###      BBM(1,ne) = alveolar volume mm3
C###      BBM(2,ne) = concentration (fractional)
C###      ACINUS(nm_time,ne)=time elapsed for RBC in capillaries
C###      ACINUS(nm_avPAO2,ne)=breath-averaged alveolar air PO2
C###      ACINUS(nm_avPcO2,ne)=breath-averaged end-capillary PO2
C###      ACINUS(nm_PcO2,ne)=end-capillary PO2
C###  For each BBM, update capillary blood volume and air-blood surface 
C###  area, based on volume change. Calculate new PAO2, PBO2.
C**** Created by AJS, Feb 2010
     
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'     
      INCLUDE 'lung00.cmn'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'odes00.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),NENP(NPM,0:NEPM),nh,
     &  NPNE(NNM,NBFM,NEM),nr,NVJE(NNM,NBFM,NJM,NEM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNP(NKM,NVM,NHM,NPM)
      REAL*8 ACINUS(4,NEM),BBM(2,NEM),CE(NMM,NEM),T,XAB(NORM,NEM),
     &  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IPAR(0:MAX_PAR),gen,nb,ne,ne0,noelem,nn,np,nv,nv2,
     &  ny,kount
      REAL*8 conc_new,conc_old,DCO2,DO2,flux,O2,O2CONTENT_VENOUS,pH,
     &  PO2,PCO2,QTotal,RPAR(0:MAX_PAR),SO2,SO2_VENOUS,tbeg,Temp,
     &  tend,Th,volCapBlood,Y(1)
      CHARACTER ERRMSG*50
      EXTERNAL DPO2DT,JAC,MAS,SOLOUT

      CALL ENTERS('GAS_EXCH_ACINUS',*9999)
!       nk=1  !no derivatives

C Parameters - ACINUS_SURFACE_AREA, ACINUS_VBLOOD set in command file
      Th=0.45d0*MCH/(MCV*MW) !mol/l assumes Hct=0.45
!       BARRIER_THICKNESS=1.11d-6 !air-blood barrier thickness (tissue+plasma), m
      pH=7.40d0 !pH of plasma
      Temp=37.0d0 !body temperature
      PCO2=40.0d0 !assume PCO2 is constant 
      O2CONTENT=0.0d0
      QTotal=0.0d0
      DO2_TOTAL=0.0d0

C Calculate airside PO2 for each conducting airway using concentrations
C Assumes air total pressure in conducting airways = atmospheric pressure
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        nb=NBJ(nj_source,ne)
        DO nn=1,NNT(nb)
          nv=NVJE(nn,nb,nj_source,ne)
          np=NPNE(nn,nb,ne)
          ny=NYNP(1,1,nh,np)
          XP(1,nv,nj_poa,np)=YP(ny,1)*(O2molVol*1.d3)*(P_ATM-p_water) 
        ENDDO
      ENDDO

C Calculate blood PO2 and O2 flux
      kount=0
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        nb=NBJ(nj_flow,ne)
        nv=NVJE(2,nb,nj_flow,ne)
        np=NPNE(2,nb,ne) !terminal node
        IF(NXI(1,0,ne).EQ.0)THEN ! terminal
          kount=kount+1
C Initial PO2 and saturations
          PO2=XP(1,nv,nj_pob,np) 
          CALL CALC_O2_KELMAN(O2,Th,pH,PCO2,PO2,SO2,Temp,ERROR,*9999)
C Diffusing capacities l/s/mmHg
          volCapBlood=XP(1,nv,nj_Vc,np)*1.0d-6 !convert mm3 to litre
          CALL CALC_DIFF_CAPS(DCO2,7.d-6,DO2,0.45d0,ACINUS_SURFACE_AREA,
     &      SO2,BARRIER_THICKNESS,volCapBlood,ERROR,*9999) !Vblood set in ACINUS array for each ne
!           CALL CALC_DIFF_CAPS(DCO2,7.d-6,DO2,0.45d0,ACINUS_SURFACE_AREA,
!      &      SO2,BARRIER_THICKNESS,ACINUS_VBLOOD,ERROR,*9999)!set Dh=7micron (diameter of capillary)
          DO2_TOTAL=DO2_TOTAL+DO2 !Sum diffusing conductance for all elements
C Save parameters into parameter array 
          RPAR(1)=P_ATM !temp
          RPAR(2)=XP(1,nv,nj_Vc,np)*1.0d-6 !capillary blood volume, litre
!           RPAR(2)=ACINUS_VBLOOD !capillary blood volume, litre
          RPAR(3)=DO2 !O2 diffusing conductance, l/s/mmHg
          RPAR(4)=BBM(2,ne)*(O2molVol*1.d3)! fractional 
!           RPAR(4)=YP(ny,1)*(O2molVol*1.d3) !air O2 concentration, convert from mmol/L to fractional
          RPAR(5)=Th !hemoglobin concentration
C Update acinar airside PO2=conc*(PA-pw) (where conc is fractional)
          XP(1,nv,nj_poa,np)=RPAR(4)*(RPAR(1)-p_water) 
C Reset blood PO2 to mixed venous partial pressures at start of transit
          IF(ACINUS(nm_time,ne).EQ.0.d0)THEN
            XP(1,nv,nj_pob,np)=MIXED_VENOUS_PO2 !40.0d0
!             CALL CALC_O2_KELMAN(O2,Th,pH,PCO2,PBO2,SO2,Temp,
!      &        ERROR,*9999)
          ELSE 
C Initial conditions for odes 
            Y(1)=XP(1,nv,nj_pob,np)  !blood PO2
C Call ode solver. Solve from t=blood_transit to t=blood_transit+tstep.
            tend=ACINUS(nm_time,ne)+DT
            tbeg=ACINUS(nm_time,ne)
            CALL PREODES(1,ERROR,*9999) !number of variables=1
            CALL RADAU5(1,DPO2DT,tbeg,Y,tend,H,RTOL,ATOL,ITOL,JAC,IJAC,
     &        MLJAC,MUJAC,MAS,IMAS,MLMAS,MUMAS,SOLOUT,IOUT,WORK,LWORK,
     &        IWORK,LIWORK,RPAR,IPAR,IDID)
C Error conditions
            ERRMSG=' '
            IF(IDID.EQ.2) ERRMSG=' >> COMPUT. SUCCESSFUL '// 
     &        '(INTERRUPTED BY SOLOUT)'
            IF(IDID.EQ.-1) ERRMSG=' >> INPUT IS NOT CONSISTENT'
            IF(IDID.EQ.-2) ERRMSG=' >> LARGER NMAX IS NEEDED'
            IF(IDID.EQ.-3) ERRMSG=' >> STEP SIZE BECOMES TOO SMALL'
            IF(IDID.EQ.-4) ERRMSG=' >> MATRIX IS REPEATEDLY SINGULAR' 
            CALL ASSERT(IDID.EQ.1,ERRMSG,ERROR,*9999)
            IF(IDID.NE.1)WRITE(*,*) '*** ERROR for node ',np,ERRMSG
C Store new solution for np in XP
            XP(1,nv,nj_pob,np)=Y(1)
            CALL CALC_O2_KELMAN(O2,Th,pH,PCO2,XP(1,nv,nj_pob,np),SO2,
     &        Temp,ERROR,*9999)

          ENDIF !tt_elapsed

C Flux in mmol/s 
          flux=-DO2/(standard_molar_vol*1.d-3)*(XP(1,nv,nj_poa,np)-
     &      XP(1,nv,nj_pob,np)) !convert DO2 to mmol(O2)/s/mmHg
          XP(1,nv,nj_source,np)=flux
!
C Reset elapsed time if blood transit is greater than total RBC transit time
C Update end capillary blood PO2 for acinus if tt_elapsed>ttRBC
          ACINUS(nm_time,ne)=ACINUS(nm_time,ne)+DT
          IF(ACINUS(nm_time,ne).GT.XP(1,nv,nj_tt,np)) THEN
            ACINUS(nm_time,ne)=0.0d0 
            ACINUS(nm_PcO2,ne)=XP(1,nv,nj_pob,np) !Record "end-capillary" PbO2
            WRITE(OP_STRING,'(''>>Error: End-capillary PO2 less than '',
     &        ''zero: element='',I5,'', PO2='',F8.2)')
     &        ne,ACINUS(nm_PcO2,ne)
            CALL ASSERT(ACINUS(nm_PcO2,ne).GE.0.0d0,OP_STRING,
     &        ERROR,*9999)
          ENDIF !tt_elapsed
          CALL CALC_O2_KELMAN(O2,Th,pH,PCO2,ACINUS(nm_PcO2,ne),SO2,Temp,
     &      ERROR,*9999) !end-capillary O2 content
          nv2=NVJE(2,nb,nj_Qdot,ne)
          QTotal=QTotal+ABS(XP(1,nv2,nj_Qdot,np)) !mm3/s
          O2CONTENT=O2CONTENT+(O2*ABS(XP(1,nv2,nj_Qdot,np))) !flow-weighted

        ENDIF !terminal
      ENDDO !noelem

C Mixed end-capillary blood (arterial)
!       CALL ASSERT(QTotal.GT.0.0d0,' >> Total perfusion <= 0',
!      &  ERROR,*9999)

C MHT 06-04-11 commented out for now. needs to do correct summation
c      CALL ASSERT(QTotal.GT.0.999d0.AND.QTotal.LT.1.001d0,
c     &  ' >> Total perfusion not equal to 1. Acinar perfusion values'//
c     &  ' should be fractional.',ERROR,*9999)
      O2CONTENT=O2CONTENT/QTotal !ml O2 / ml blood
      CALL CALC_PO2_FROM_O2(O2CONTENT,Th,pH,PCO2,PO2_arterial,
     &  SO2_arterial,Temp,ERROR,*9999)

C Calculate O2 uptake = bloodflowx(CaO2 - CvO2)
      CALL CALC_O2_KELMAN(O2CONTENT_VENOUS,Th,pH,PCO2,MIXED_VENOUS_PO2,
     &  SO2_VENOUS,Temp,ERROR,*9999) 
! CARDIAC_OUTPUT is in units of mm3/s (convert to ml per min)
      O2UPTAKE=CARDIAC_OUTPUT*1.d-3*60.d0*(O2CONTENT-O2CONTENT_VENOUS)

      CALL EXITS('GAS_EXCH_ACINUS')
      RETURN
 9999 CALL ERRORS('GAS_EXCH_ACINUS',ERROR)
      CALL EXITS('GAS_EXCH_ACINUS')
      RETURN 1
      END
