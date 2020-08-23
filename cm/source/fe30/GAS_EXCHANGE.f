      SUBROUTINE GAS_EXCHANGE(NPLIST,nr,nx,XP,ERROR,*)

C#### Subroutine: GAS_EXCHANGE
C###  Description:
C###     Calculates the partial pressures (mmHg) in a unit (typically  
C###     an acinus). Air-side and blood-side partial pressures are 
C###     evaluated.
C###  Inputs required for each unit: CAP_BLOOD_VOL, AIR_BLOOD_SURFACE,
C###    HCT_INITIAL, ALV_AIR_VOL 
C###  Pre-defined transport fields: Vdot (ventilation), Qdot (blood 
C###    flow), PA (alveolar pressure)
C**** Created by AJS, August 2007

     
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'     
      INCLUDE 'lungex00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'odes00.cmn'
      
!     Parameter List
      INTEGER NPLIST(0:NP_R_M),nr
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER n,nk,nonode,np,nv,nx
      REAL*8 DCO2,DO2,dt,HCT,IPAR(0:MAX_PAR),PO2,SO2,
     &  RPAR(0:MAX_PAR),t,tend,Th,Y(1:MAX_ODE_DIM)
      CHARACTER CHAR1*5,CHAR2*12,CHAR3*12,CHAR4*12,CHAR5*12,CHAR6*12,
     &  CHAR7*12,CHAR8*12,ERRMSG*50
      EXTERNAL BENTAL_2006,JAC,MAS,SOLOUT
          
      CALL ENTERS('GAS_EXCHANGE',*9999)

      PI=3.14159265358979d0 !set value of PI (declared in b00.cmn)

      CALL ASSERT((ITYP7(nr,nx).EQ.2.OR.ITYP7(nr,nx).EQ.3),
     &  '>> No gas exchange model specified.',ERROR,*9999)

C     Field indices
      nj_pob=NJ_LOC(NJL_FIEL,1,nr) !=4 
      nj_pcb=NJ_LOC(NJL_FIEL,2,nr) !=5 
      nj_poa=NJ_LOC(NJL_FIEL,3,nr) !=6 
      nj_pca=NJ_LOC(NJL_FIEL,4,nr) !=7 
      nv=1  !one value at each node
      nk=1  !no derivatives
      n=6   !Number of solution variables in array Y

C     Hemoglobin concentration (node-independent)
      Th=HCT_INITIAL*MCH/(MCV*MW) !mol/l
     
C     For each node
      DO nonode=1,NPLIST(0) 

C       Initialise arrays
        np=NPLIST(nonode)
        XP(nk,nv,nj_pob,np)=INITIAL_POB
        XP(nk,nv,nj_pcb,np)=INITIAL_PCB
        XP(nk,nv,nj_poa,np)=INITIAL_POA
        XP(nk,nv,nj_pca,np)=INITIAL_PCA

C       Air volume
        ALV_AIR_VOL=ALV_AIR_VOL+XP(nk,nv,nj_Vdot,np) !initial vol + vol change

C       Initial PO2 and saturations
        PO2=XP(nk,nv,nj_pob,np) 
        SO2=(L*kt*sigma_o*PO2*(1.0d0+kt*sigma_o*PO2)**3.0d0+kr*
     &    sigma_o*PO2*(1.0d0+kr*sigma_o*PO2)**3)/(L*(1.0d0+kt*
     &    sigma_o*PO2)**4.0d0+(1.0d0+kr*sigma_o*PO2)**4.0d0)

        CALL CALC_DIFF_CAPS(DCO2,DO2,HCT_INITIAL,AIR_BLOOD_SURFACE,
     &    SO2,0.0d0,CAP_BLOOD_VOL,ERROR,*9999)

        IF(ITYP7(nr,nx).EQ.2)THEN !Hill 1973
         ! to be coded

        ELSEIF(ITYP7(nr,nx).EQ.3)THEN !Ben-Tal 2006

C         Save parameters into parameter array for BENTAL_2006
          GASPAR(1)=P_ATM 
          GASPAR(2)=0.21d0            !O2 conc in the mouth (fom), mmHg  
          GASPAR(3)=0.0d0             !CO2 conc in the mouth (fcm), mmHg 
          GASPAR(4)=LUNG_ELASTANCE    !lung elastance, mmHg/l 
          GASPAR(5)=ALV_AIR_VOL       !alveolar air volume, litre
          GASPAR(6)=CAP_BLOOD_VOL     !capillary blood volume, litre
          GASPAR(7)=Th                !hemoglobin concentration, mol/l
          GASPAR(8)=DO2               !O2 diffusing conductance, l/s/mmHg
          GASPAR(9)=DCO2              !CO2 diffusing conductance, l/s/mmHg    

C         Initial conditions for odes Y = [Pa, fo, fc, po, pc, z]
          Y(1)=XP(nk,nv,nj_PA,np)
          Y(2)=XP(nk,nv,nj_poa,np)/(XP(nk,nv,nj_PA,np)-p_water)
          Y(3)=XP(nk,nv,nj_pca,np)/(XP(nk,nv,nj_PA,np)-p_water)
          Y(4)=XP(nk,nv,nj_pob,np)  
          Y(5)=XP(nk,nv,nj_pcb,np)
          Y(6)=46.0d0*sigma_c*rtwo/(ltwo*hconc)
          t=0.0d0

C         Time that blood is in capillary bed (integration time)
          IF(XP(nk,nv,nj_Qdot,np).EQ.0.0d0)THEN
            WRITE(OP_STRING,'('' Blood flow is zero for node '',I5)')np
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            GOTO 100
          ENDIF
          tend=CAP_BLOOD_VOL/XP(nk,nv,nj_Qdot,np)

C         Call ode solver
          CALL PREODES(n,ERROR,*9999)
          CALL RADAU5(n,BENTAL_2006,t,Y,tend,H,RTOL,ATOL,ITOL,JAC,
     &      IJAC,MLJAC,MUJAC,MAS,IMAS,MLMAS,MUMAS,SOLOUT,IOUT,WORK,
     &      LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)

C         Error conditions
          ERRMSG=' '
          IF(IDID.EQ.2) ERRMSG=' >> COMPUT. SUCCESSFUL '// 
     &      '(INTERRUPTED BY SOLOUT)'
          IF(IDID.EQ.-1) ERRMSG=' >> INPUT IS NOT CONSISTENT'
          IF(IDID.EQ.-2) ERRMSG=' >> LARGER NMAX IS NEEDED'
          IF(IDID.EQ.-3) ERRMSG=' >> STEP SIZE BECOMES TOO SMALL'
          IF(IDID.EQ.-3) ERRMSG=' >> MATRIX IS REPEATEDLY SINGULAR' 
          CALL ASSERT(IDID.EQ.1,ERRMSG,ERROR,*9999)

C         Store solution for np in XP	    
          XP(nk,nv,nj_pob,np)=Y(4)
          XP(nk,nv,nj_pcb,np)=Y(5) 
          XP(nk,nv,nj_poa,np)=Y(2)*(Y(1)-p_water)
          XP(nk,nv,nj_pca,np)=Y(3)*(Y(1)-p_water)

        ENDIF !ITYP7

100   ENDDO !nonode

                           
      CALL EXITS('GAS_EXCHANGE')
      RETURN
 9999 CALL ERRORS('GAS_EXCHANGE',ERROR)
      CALL EXITS('GAS_EXCHANGE')
      RETURN 1
      END
