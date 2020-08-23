      SUBROUTINE OPLUNG(nb,NEELEM,NPNE,NPNODE,nr,nx,NVJE,NYNP,
     &  BBM,CE,CP,T,XAB,XP,YP,GAS,ERROR,*)

C#### Subroutine: OPLUNG
C###  Description:
C###    OPLUNG outputs information during simulations of gas mixing
C###    or water and heat transfer in lung models.
C***  Created by Merryn Howatson Tawhai, August 1999

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
c      INCLUDE 'mach00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER nb,NEELEM(0:NE_R_M),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M),
     &  nr,nx,NVJE(NNM,NBFM,NJM,NEM),NYNP(NKM,NVM,NHM,NPM)
      REAL*8 BBM(2,NEM),CE(NMM,NEM),CP(NMM,NPM),
     &  T,XAB(NORM,NEM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL GAS
!     Local Variables
      INTEGER j,ne,ne_max,ne_min,nh,nlpm,nonode,np,nv,ny_top,ny1,ny2
      REAL*8 MAX_CONC,MEAN_CONC,MIN_CONC,VALUES(3,10),ZMAX,ZMIN
      REAL*8 CONC_EVAL,MassErr,VolErr
      LOGICAL PRESSURE

      CALL ENTERS('OPLUNG',*9999)
C Thermofluid model output:
C temperature at entrance of model, and other specified nodes
C absolute humidity at entrance of model, and other specified nodes
C relative humidity at entrance of model, and other specified nodes
C snapshot of temperature and humidity everywhere
      
      nh=NH_LOC(1,nx)
      ny_top=NYNP(1,1,nh,NPNODE(1))

C....Temperature & humidity simulations
      IF(.NOT.GAS)THEN
        IF(N_SOLN.EQ.IWRIT6)THEN
C write output to screen          
          WRITE(OP_STRING,'('' Time='',F12.4,''   DV='',
     &      F7.4,''   T_inlet='',F7.2,'' degC  C_inlet='',F7.2,
     &      '' mg/L  Flow='',F5.2,'' L/s'')') T,DV_TOTAL(nx)/1.d6,
     &      YP(ny_top,1)-273.15d0,YP(ny_top+1,1),INLET_FLOW(nx)/1.d6
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

C write output to history file
        IF(LUNG_OP.GE.1)THEN !write to .history file
          np=NPNODE(1) !first node
          ny1=NYNP(1,1,nh,np) !for temperature in K
          ny2=NYNP(1,1,nh+1,np) !for humidity in g/mm^3
          MAX_CONC=CONC_EVAL(YP(ny1,1))
          VALUES(1,1)=YP(ny1,1)-273.15d0 !temperature
          VALUES(2,1)=YP(ny2,1) !absolute humidity
          VALUES(3,1)=YP(ny2,1)/MAX_CONC*100.d0 !relative humidity
          DO nonode=1,NP_OP_LIST(0)
            np=NP_OP_LIST(nonode)
            ny1=NYNP(1,1,nh,np) !for temperature in K
            ny2=NYNP(1,1,nh+1,np) !for humidity in g/mm^3
            MAX_CONC=CONC_EVAL(YP(ny1,1))
            VALUES(1,1+nonode)=YP(ny1,1)-273.15d0 !temperature
            VALUES(2,1+nonode)=YP(ny2,1) !absolute humidity
            VALUES(3,1+nonode)=YP(ny2,1)/MAX_CONC*100.d0 !relative humidity
          ENDDO !nonode
          WRITE(IOFILE2,'(F10.4,30(F8.2))')  T,(VALUES(1,j),j=1,
     &      NP_OP_LIST(0)+1),(VALUES(2,j),j=1,NP_OP_LIST(0)+1),
     &      (VALUES(3,j),j=1,NP_OP_LIST(0)+1)
        ENDIF
        
        IF(LUNG_OP.GE.2)THEN
          np=NPNODE(1) !first node
          VALUES(1,1)=CP(4,np)-273.15d0 !temperature
          VALUES(2,1)=CP(8,np) !depth
          DO nonode=1,NP_OP_LIST(0)
            np=NP_OP_LIST(nonode)
            VALUES(1,1+nonode)=CP(4,np)-273.15d0 !temperature
            VALUES(2,1+nonode)=CP(8,np) !depth
          ENDDO !nonode
          WRITE(IOFILE3,'(F10.4,30(F12.6))')  T,(VALUES(1,j),j=1,
     &      NP_OP_LIST(0)+1),(VALUES(2,j),j=1,NP_OP_LIST(0)+1)
        ENDIF
      ENDIF !.NOT. GAS
      
C....Gas transport simulations
      IF(GAS)THEN
        IF(N_SOLN.EQ.IWRIT6)THEN

C Volume and mass errors
          VolErr=1.d2*(CURRENT_VOLUME(nx)-(INITIAL_VOLUME
     &      +DV_TOTAL(nx)))/(INITIAL_VOLUME+DV_TOTAL(nx))
          MassErr=1.d2*(MASS_CURRENT(nx)-IDEAL_MASS(nx))/
     &      IDEAL_MASS(nx)

C write output to screen        
         IF(ITYP7(nr,nx).EQ.1)THEN !not gas exchange
            VolErr=1.d2*(CURRENT_VOLUME(nx)-(INITIAL_VOLUME
     &        +DV_TOTAL(nx)))/(INITIAL_VOLUME+DV_TOTAL(nx))
            MassErr=1.d2*(MASS_CURRENT(nx)-IDEAL_MASS(nx))/
     &        IDEAL_MASS(nx)
            IF(DV_TOTAL(nx).LT.1.d3)THEN !small values
              WRITE(OP_STRING,'('' Time='',F8.3,''   DV='',E12.4,
     &        ''   Mass='',E12.4,''   Flow='',E12.4,
     &        ''   %VolERR='',F10.3,''   %MassERR='',F10.3,
     &        ''   Inlet='',E12.4)')
     &        T,DV_TOTAL(nx),MASS_CURRENT(nx),INLET_FLOW(nx),VolErr,
     &        MassErr,YP(ny_top,1)
            ELSE
              WRITE(OP_STRING,'('' Time='',F8.3,''   DV='',E12.4,
     &        ''   Mass='',E12.4,''   Flow='',E12.4,
     &        ''   %VolERR='',F10.3,''   %MassERR='',F10.3,
     &        ''   Inlet='',E12.4)')
     &        T,DV_TOTAL(nx)/1.d6,MASS_CURRENT(nx)/1.d6,
     &        INLET_FLOW(nx)/1.d6,VolErr,MassErr,YP(ny_top,1)
            ENDIF
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

          ELSE !gas exchange

C AJS output partial pressures
            WRITE(OP_STRING,'(X,F7.3,X,F7.3,X,F7.3,X,F12.3,X,
     &        F12.3,X,F12.3,X,E14.4,X,E16.4)') 
     &        T,VolErr,MassErr,YP(1,1)*(O2molVol*1.0d3)*(P_ATM-p_water),
     &        PAO2_AVERAGE,PO2_arterial,DO2_TOTAL*1.d3,O2UPTAKE
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

            IF(nx.EQ.2.AND.(TFINISH-T).LT.DT)THEN !end of expiration
              WRITE(OP_STRING,'('' Spatial and temporal average PAO2 '//
     &          '(averaged for all acini and averaged over a breath) '//
     &          '= '',F8.3,'' mmHg'')')PAO2_BREATH
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Temporal average arterial PO2 '//
     &          '(averaged over a breath) = '',F8.3,'' mmHg'')')
     &          PartO2_BREATH
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Average P(A-a)O2 over a breath = '',
     &          F8.3,'' mmHg'')')PAO2_BREATH-PartO2_BREATH
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          
          ENDIF!gas exchange

C write output to snhist file        
!           WRITE(IOFILE3,'(E14.4)') MassErr

        ENDIF !N_SOLN

C write gas exchange output to files
        IF(ITYP7(nr,nx).GT.1)THEN
          IF(LUNG_OP.GE.1)THEN !.history file
            WRITE(IOFILE2,'(4(F10.3))') T,PAO2_AVERAGE,
     &        PO2_arterial,YP(1,1)*(O2molVol*1.0d3)*(P_ATM-p_water)
          ENDIF
          IF(LUNG_OP.GE.2)THEN!.error file
            WRITE(IOFILE3,'(E14.4)') MassErr
          ENDIF
          IF(LUNG_OP.GE.3.AND.NP_OP_LIST(0).GT.0)THEN !.nodes file
            DO nonode=1,NP_OP_LIST(0)
              np=NP_OP_LIST(nonode)
              VALUES(1,nonode)=XP(1,1,nj_poa,np) !airside PO2
            ENDDO !nonode
            WRITE(IOFILE4,'(F10.4,30(F8.3))')  T,(VALUES(1,j),j=1,
     &        NP_OP_LIST(0))
          ENDIF
        ELSE

C write output to snhist file        
          WRITE(IOFILE3,'(E14.4)') MassErr

C write output to history file 
          IF(IWRIT1(nr,nx).NE.0)THEN !write to .history file
            WRITE(IOFILE2,'(2(E15.5))') SUM_EXP_VOL/1.d6,YP(1,1)
          ENDIF

        ENDIF !gas exchange

      ENDIF !GAS

      IF(N_SOLN.EQ.IWRIT6) N_SOLN=0
      
      CALL EXITS('OPLUNG')
      RETURN

 9999 CALL ERRORS('OPLUNG',ERROR)
      CALL EXITS('OPLUNG')
      RETURN 1
      END



