      SUBROUTINE EVFLOW_DYNAM(Gdirn,ITER_USER,MECHANICS_FILETYPE,NBJ,
     &  NBREATHS,NDP,NEELEM,NENP,NORD,NPLIST,NPNE,NPNODE,NVJE,NVJP,NXI,
     &  COV,CW,dt,dt_init,dt_max,ERR_USER,FlowIN,FlowOUT,I_TO_E_RATIO,
     &  MeanCompliance,Pmax,Pmin,
     &  Ppl_step,PressureIN,T_interval,refvol,RMaxMean,RMinMean,
     &  volume_target,BBM,CE,
     &  undef,XAB,XP,ZD,COMPLIANCE_BC,DIAG_OP,FIRST_ORDER,
     &  INITIAL,NORMALISE,PATHLENGTHS,PEAK,PRINT,
     &  READ_VOLUMES,SETUP,UNIFORM,EXTEND_FRC,filename,P_TYPE,ERROR,*)

C#### Subroutine: EVFLOW_DYNAM
C###  Description:
C###    EVFLOW_DYNAM solves pulmonary 1D flow distributions. Loop over 
C###    solution for given number of breaths, or until tidal volume target
C###    is reached.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn' 
      INCLUDE 'geom00.cmn' 
      INCLUDE 'lung00.cmn' 
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'b00.cmn' 
      INCLUDE 'valu00.cmn'
!     Parameter List
      INTEGER Gdirn,ITER_USER,MECHANICS_FILETYPE,NBJ(NJM,NEM),NBREATHS,
     &  NDP(NDM),NEELEM(0:NE_R_M),NENP(NPM,0:NEPM),NORD(5,NE_R_M),
     &  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M),
     &  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 COV,CW,dPl,dt,dt_init,dt_max,ERR_USER,volume_target,
     &  BBM(2,NEM),CE(NMM,NEM),FlowIN,FlowOUT,I_TO_E_RATIO,
     &  MeanCompliance,Pmax,Pmin,
     &  Ppl_step,PressureIN,undef,refvol,RMaxMean,
     &  RMinMean,T_interval,XAB(NORM,NEM),XP(NKM,NVM,NJM,NPM),
     &  ZD(NJM,NDM)
      LOGICAL COMPLIANCE_BC,DIAG_OP,CONSTANT,CONTINUE,FIRST_ORDER,
     &  INITIAL,NORMALISE,PEAK,PRINT,
     &  READ_VOLUMES,SETUP,UNIFORM,PATHLENGTHS
      CHARACTER EXTEND_FRC*1,filename*200,FILENAME2*200,
     &  P_TYPE*50,ERROR*(*)
!     Local variables
      REAL*8 endtime,FRC,InitialVolume,maxPmus,maxPpl,
     &  MeanVolume,Pmus,Ppl_current,
     &  Ppl_factor,Ppl_inc,PressureIN_step,previous_flow,sum_exhaled,
     &  sum_tidal,sum_total,time,sumvolume,ttime,Texpn,Tinsp
      INTEGER count_exports,N,nb,ncount,ne,ne0,ngen,noelem,nonode,np,
     &  i,np1,np2,nv,nv1,nv2,Nterms,Nend,IBEG,IEND,step
      LOGICAL FIRST,FAST,FIRST_READ

      CALL ENTERS('EVFLOW_DYNAM',*9999)   

c      FAST = .TRUE.
      FAST=.FALSE.
      CONSTANT=.FALSE.
      previous_flow=1.d0 ! initialise to a positive value
      
      IF(UNIFORM)THEN
        IF(nj_flow.EQ.0)THEN
          nj_flow=4 !AJS - should not be hard-coded
        ENDIF
C.....Calculate the mesh volume
        CALL VOLUMEOFMESH(NBJ,NEELEM,NORD,NPNE,NVJE,NXI,MeanVolume,
     &    BBM,CE,XP,ERROR,*9999)
        INITIAL_VOLUME=CE(nm_vol_bel,1)
c        call SetInitialVolume(Gdirn,NBJ,NEELEM,NPNE,NVJE,NXI,BBM,CE,
c     &    COV,MeanVolume,refvol,RMaxMean,RMinMean,sumvolume,undef,
c     &    XP,ERROR,*9999)
        
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          ne0=NXI(-1,1,ne) !parent element
          nb=NBJ(1,ne)
          np1=NPNE(1,nb,ne) !start node
          np2=NPNE(2,nb,ne) !end node
          nv1=NVJE(1,nb,nj_flow,ne)
          nv2=NVJE(2,nb,nj_flow,ne)

          XP(1,nv1,nj_flow,np1)=(CE(nm_vol_bel,ne)-CE(nm_volumes,ne))
     &      /(CE(nm_vol_bel,1)-CE(nm_volumes,1))
          XP(1,nv2,nj_flow,np2)=(CE(nm_vol_bel,ne)-CE(nm_volumes,ne))
     &      /(CE(nm_vol_bel,1)-CE(nm_volumes,1))
          XAB(1,ne)=0.d0 !proportion of volume change
        ENDDO
        
      ELSE !not uniform

! AJS 16/2/11 moved error check so that it doesn't occur for uniform flow option
        CALL ASSERT(PARAMETERS_DEFINED,
     &    'Define parameter values first (fem define values)',
     &     ERROR,*9999) 

        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          CE(nm_Ppl,ne)=0.d0 !initialise the pleural pressure at an acinus
          XAB(1,ne)=0.d0 !initialise the adjacent pressure
          XAB(3,ne)=0.d0 !initialise the adjacent pressure
        ENDDO
        
        CALL VOLUMEOFMESH(NBJ,NEELEM,NORD,NPNE,NVJE,NXI,MeanVolume,
     &   BBM,CE,XP,ERROR,*9999)
        CALL SET_USER_DOUBLE('deadspace',CE(nm_volumes,1)/1.d3,ERROR)
        IF(CE(nm_volumes,1)/1.d3.GT.1.d3)THEN
          WRITE(OP_STRING,'(''Anatomical deadspace = '',F8.3,'' L'')')
     &      CE(nm_volumes,1)/1.d6
        ELSE
          WRITE(OP_STRING,'(''Anatomical deadspace = '',F8.3,'' ml'')')
     &      CE(nm_volumes,1)/1.d3
        ENDIF
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        Nterms=0
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          IF(NXI(1,0,ne).EQ.0)THEN
            Nterms=Nterms+1
          ENDIF
        ENDDO

        IF(READ_VOLUMES)THEN
          FIRST_READ=.TRUE.
          CALL READINITIALVOLUME(MECHANICS_FILETYPE,NBJ,NDP,NEELEM,
     &      NENP,NPNE,NXI,BBM,CE,sumvolume,undef,XP,ZD,FIRST_READ,
     &      ERROR,*9999)
          FIRST_READ=.FALSE.
        ELSE
          CALL SETINITIALVOLUME(Gdirn,NBJ,NEELEM,NPNE,NVJE,NXI,BBM,CE,
     &      COV,MeanVolume,refvol,RMaxMean,RMinMean,sumvolume,undef,
     &      XP,filename,ERROR,*9999)
        ENDIF

        IF(sumvolume/1.d3.GT.1.d3)THEN
          WRITE(OP_STRING,'(''Respiratory volume   = '',F8.3,'' L'')')
     &      sumvolume/1.d6
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,'(''Respiratory volume   = '',F8.3,'' ml'')')
     &      sumvolume/1.d3
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

        IF(.NOT.SETUP)THEN !solve the flow problem
C.......Set the initial pleural pressure balanced by elastic recoil        
          CALL TISSUECOMPLIANCE(NEELEM,NPNE,NXI,BBM,CE,CW,dPl,
     &      undef,.TRUE.,ERROR,*9999)
C MHT. Calculate the initial pleural pressure as -Pel.
C i.e. CE(nm_Ppl) contains Pel
          Ppl_current=0.d0
          Nend=0
          DO noelem=1,NEELEM(0)
            ne=NEELEM(noelem)
            IF(NXI(1,0,ne).EQ.0)THEN
              Ppl_current=Ppl_current-CE(nm_Ppl,ne)
              Nend=Nend+1
            ENDIF
          ENDDO
          Ppl_current=Ppl_current/DBLE(Nend)

          maxPmus=0.d0
          maxPpl=0.d0

          FRC=sumvolume
          InitialVolume = FRC + CE(nm_volumes,1)
          initial_volume=InitialVolume !vol at start of simulation
          
          WRITE(OP_STRING,
     &      '(''   Time        Inflow        V_tidal       '//
     &      ' Raw          Rtotal      Compliance    Palv-Pmouth      '
     &      //' Ppl           Ptp       Vol Lung        Pmusc'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''   (s)         (mL/sec)        (mL)   '
     &      //'   (cmH2O/L.s)   (cmH2O/L.s)    (L/cmH2O)      (cmH2O)'
     &      //'       (cmH2O)        (cmH2O)      (mL)      (cmH2O)'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          
          PressureIN_step=0.d0
          time = 0.d0
          sum_tidal=0.d0
          sum_total=0.d0
          sum_exhaled=0.d0
          
          CALL BRANCHRESISTANCE(NBJ,NEELEM,NPNE,NVJE,CE,XP,
     &      ERROR,*9999)
          CALL SET_DPPL_GRAD(Gdirn,NBJ,NEELEM,NPNE,NVJE,NXI,CE,
     &      Pmax,Pmin,Ppl_step,XP,ERROR,*9999)

          sum_tidal=0.d0
          sum_exhaled=0.d0
          n=0

          DO noelem=1,NEELEM(0)
             ne=NEELEM(noelem)
             BBM(2,ne)=0.d0
          ENDDO

          Ppl_factor=1.d0
          CONTINUE=.TRUE.
          
cc          sum_tidal=volume_target
            sum_tidal=0.d0
          DO WHILE(CONTINUE)
            n=n+1
            write(*,*) 'breath',n,ttime,(Tinsp+Texpn),endtime
            PRINT=.TRUE.
C     Modify the driving pressure to get close to the target volume          
cc            IF(n.GT.1.AND.DABS(volume_target).GT.1.d-5)THEN
cc!     modify driving pressure by volume_target/sum_tidal
cc!     this increases Ppl for volume_target>sum_tidal, and
cc!     decreases Ppl for volume_target<sum_tidal
cc              Ppl_factor=Ppl_factor*DABS(volume_target/sum_tidal)
cc            ENDIF
             
            ttime=0.d0
c            sum_tidal=0.d0
            Pmus = 0.d0         !initialise the muscle pressure to zero
            FIRST=.TRUE.
            endtime=T_interval*DBLE(n)-0.5d0*dt

            Texpn = T_interval / (1.d0+I_TO_E_RATIO)
            Tinsp = T_interval - Texpn

c            DO noelem=1,NEELEM(0)
c              ne=NEELEM(noelem)
c              BBM(2,ne)=0.d0
c            ENDDO
            count_exports=0
            ncount=0
            IF(READ_VOLUMES)THEN
              CALL READINITIALVOLUME(MECHANICS_FILETYPE,NBJ,NDP,NEELEM,
     &          NENP,NPNE,NXI,BBM,CE,sumvolume,undef,XP,ZD,
     &          FIRST_READ,ERROR,*9999)
            ELSE
c              CALL RESETINITIALVOLUME(NEELEM,NXI,BBM,CE,sumvolume,
c     &          ERROR,*9999) !resets BBM with gravitational effect
            ENDIF
            CALL VOLUMEOFMESH(NBJ,NEELEM,NORD,NPNE,NVJE,NXI,MeanVolume,
     &        BBM,CE,XP,ERROR,*9999)
c            InitialVolume=CE(nm_vol_bel,1)

            step=0
            DO WHILE (time.LT.endtime)
              step=step+1
              time = time + dt
              ttime=ttime+dt
              PressureIN_step=PressureIN_step+PressureIN/T_interval*dt

              IF(P_TYPE(1:8).EQ.'SINUSOID')THEN
                 IF(ttime.LT.Tinsp)THEN
                    dpl=Ppl_factor*pi*sin(2.d0*pi/(2.d0*Tinsp)*ttime)/
     &                   (2.d0*Tinsp)*dt
                 ELSEIF(ttime.LE.Tinsp+Texpn)THEN
                    dpl=Ppl_factor*pi*sin(2.d0*pi*(0.5d0+(ttime-Tinsp)
     &                  /(2.d0*Texpn)))/(2.d0*Texpn)*dt
                 ENDIF

              ELSEIF(P_TYPE(1:6).EQ.'LINEAR')THEN
                 IF(ttime.LT.Tinsp)THEN
                    dpl=2.d0*dt/(2.d0*Tinsp)
                 ELSEIF(ttime.LE.Tinsp+Texpn)THEN
                    dpl=-2.d0*dt/(2.d0*Texpn)
                 ENDIF
              ENDIF

C Pl_current is the current total pleural pressure, including the initial
C pressure at FRC (verified MHT 09-05-11)

              CALL BRANCHRESISTANCE(NBJ,NEELEM,NPNE,NVJE,CE,
     &          XP,ERROR,*9999)
              CALL TISSUECOMPLIANCE(NEELEM,NPNE,NXI,BBM,CE,CW,dPl,
     &          undef,.FALSE.,ERROR,*9999)
C calculate the mean pleural pressure based on current Pel (=Ptp) and Palv
C i.e. Ppl(unit) = -Pel(unit)+Palv(unit)
              ppl_inc = 0.d0
              Ppl_current = 0.d0
              Nend=0
              DO noelem=1,NEELEM(0)
                ne=NEELEM(noelem)
                IF(NXI(1,0,ne).EQ.0)THEN
                  nb=NBJ(1,ne)
                  np=NPNE(2,nb,ne)
                  ppl_inc=ppl_inc+dPl*CE(nm_dpl,ne)
                  ppl_current=ppl_current-CE(nm_Ppl,ne)+
     &              XP(1,1,nj_pressure,np)
                  Nend=Nend+1
                ENDIF
              ENDDO 

              ppl_inc=ppl_inc/DBLE(Nend)
              Pmus = Pmus + ppl_inc
              Ppl_current = Ppl_current/DBLE(Nend)

              IF(DABS(Pmus).GT.DABS(maxPmus)) maxPmus = Pmus
              IF(Ppl_current.LT.maxPpl) maxPpl = -Ppl_current

c              CALL UPDATEPRESSUREDT(NBJ,NEELEM,NPNE,NXI,CE,dPl,
c     &          dt,time,ttime,T_interval,XAB,XP,.FALSE.,ERROR,
c     &          *9999)
              previous_flow = XP(1,1,nj_flow,1)
              CALL DIRECTSOLUTION(ITER_USER,NBJ,NEELEM,NORD,NPLIST,NPNE,
     &          NVJE,NVJP,NXI,BBM,CE,CW,dPl,dt,ERR_USER,InitialVolume,
     &          MeanVolume,Ppl_current,Pmus,time,ttime,T_interval,
     &          undef,XAB,XP,DIAG_OP,FAST,ERROR,*9999)

              IF(ttime.LE.Tinsp)THEN
C MHT change this so that 'tidal' is summation of any positive flow
C              IF(XP(1,1,nj_flow,1).GT.0.d0)THEN
                sum_tidal=sum_tidal+XP(1,1,nj_flow,1)*dt
              ELSE
                sum_exhaled=sum_exhaled-XP(1,1,nj_flow,1)*dt
              ENDIF
                
              sum_total=sum_total+XP(1,1,nj_flow,1)*dt
              sumvolume=sumvolume+XP(1,1,nj_flow,1)*dt !mm^3

              IF(count_exports.EQ.0.AND.ttime/T_interval.GE.0.25d0)THEN
                DO nonode=1,NPNODE(0)
                  np=NPNODE(nonode)
                  DO nv=1,3
                    XP(1,nv,nj_radius+1,np)=XP(1,1,nj_pressure,np)
                  ENDDO
                ENDDO
              ENDIF

              IF(EXTEND_FRC(1:1).eq.'Y')THEN !option to return to FRC
                 IF(ttime.GE.(Tinsp+Texpn))THEN !check at end expn
                    IF(CE(nm_vol_bel,1).GT.initial_volume)THEN
                       endtime=endtime+dt  !extend simulation by dt
                    ENDIF
                 ENDIF
              ENDIF
C... output 'end of breath' information when changing from -ve to +ve flow
!              IF(XP(1,1,nj_flow,1).GT.0.d0.
!     &          AND.previous_flow.LT.0.d0)THEN
             !write(*,*) ttime,(Tinsp+Texpn),n,nbreaths
             IF(ttime.GE.(Tinsp+Texpn-dt*0.5d0)
     &         .AND.n.LT.nbreaths)THEN
		!write(*,*) ttime,(Tinsp+Texpn),n,nbreaths
                 CALL SET_USER_DOUBLE('tidalvol',sum_tidal/1.d3,ERROR)
                 WRITE(OP_STRING,'(''##################'')')
                 CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                 WRITE(OP_STRING,'(''For breath'',I3,'' V_T ='',F8.2,'
     &                //''' L'')') n,sum_tidal/1.d6
                 CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                 WRITE(OP_STRING,'('' Difference from target = '',F8.2,'
     &                //''' %'')') 100.d0*(volume_target-sum_tidal)
     &                /volume_target
                 CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                 WRITE(OP_STRING,'('' Max Pmus = '',F8.2,'
     &                //''' cmH2O'')') maxPmus/98.0655d0
                 CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                 WRITE(OP_STRING,'('' Max Ppl = '',F8.2,'
     &                //''' cmH2O'')') -maxPpl/98.0655d0
                 CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                 
                 WRITE(OP_STRING,'('' Change in FRC = '',F8.2,'
     &                //''' % '')') 100.d0*(CE(nm_vol_bel,1)
     &                -Initial_Volume)/Initial_Volume
                 CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                 WRITE(OP_STRING,'(''##################'')')
                 CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

                 IF(n.GT.1.AND.DABS(volume_target).GT.1.d-5)THEN
!     modify driving pressure by volume_target/sum_tidal
!     this increases Ppl for volume_target>sum_tidal, and
!     decreases Ppl for volume_target<sum_tidal
                    Ppl_factor=Ppl_factor*DABS(volume_target/sum_tidal)
                    write(*,*) 'Ppl_factor', Ppl_factor
                    write(*,*) n,'length',endtime
                 ENDIF

                 sum_tidal=0.d0 !reset the tidal volume
                 DO noelem=1,NEELEM(0)
                    ne=NEELEM(noelem)
                    BBM(2,ne)=0.d0 !reset acinar tidal volume
!                    BBM(3,ne)=BBM(1,ne) !store initial acinar volume
                 ENDDO

              ENDIF             !output

             ENDDO !while ttime<endtime

c...  CHECK WHETHER TO CONTINUE          
             IF(n.GE.nbreaths)THEN
               CONTINUE=.FALSE.
             ELSE IF(DABS(volume_target).GT.1.d-3)THEN
                IF(DABS(100.d0*(volume_target-sum_tidal)
     &               /volume_target).GT.0.1d0.OR.(n.LT.2))THEN
                   CONTINUE=.TRUE.
                ELSE
                   CONTINUE=.FALSE.
                ENDIF
             ENDIF 
 
             IF(.NOT.CONTINUE)THEN
             CALL SET_USER_DOUBLE('tidalvol',sum_tidal/1.d3,ERROR)             
             WRITE(OP_STRING,'(''##################'')')
             CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
             WRITE(OP_STRING,'(''For breath'',I3,'' V_T ='',F8.2,'
     &          //''' L'')') n,sum_tidal/1.d6
             CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
             WRITE(OP_STRING,'(''   Difference from target = '',F8.2,'
     &          //''' %'')') 100.d0*(volume_target-sum_tidal)
     &            /volume_target
             CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
             WRITE(OP_STRING,'(''   Max Pmus = '',F8.2,'
     &          //''' cmH2O'')') maxPmus/98.0655d0
             CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
             WRITE(OP_STRING,'(''   Max Ppl = '',F8.2,'
     &          //''' cmH2O'')') -maxPpl/98.0655d0
             CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

             WRITE(OP_STRING,'(''   Change in FRC = '',F8.2,'
     &          //''' % '')') 100.d0*(CE(nm_vol_bel,1)
     &          -Initial_Volume)/Initial_Volume
             CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
             WRITE(OP_STRING,'(''##################'')')
             CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
         ENDDO
          
c...END OF THE FLOW/BREATH LOOP        
          
        CALL TREERESISTANCE(NEELEM,NORD,CE,filename,ERROR,*9999)
          
C.......Set the flows to be the time averaged flow over the last inspiration
          DO noelem=1,NEELEM(0)
            ne=NEELEM(noelem)
            IF(NXI(1,0,ne).EQ.0)THEN !terminal
              np=NPNE(2,1,ne)
              XP(1,1,nj_flow,np)=BBM(2,ne)/Tinsp !in L/s
              XP(2,1,nj_flow,np)=BBM(2,ne)/Tinsp
              XP(1,1,nj_vent,np)=BBM(2,ne)/(Tinsp+Texpn) !in L/s
              XP(2,1,nj_vent,np)=BBM(2,ne)/(Tinsp+Texpn)
!              XP(1,1,nj_sV,np)=BBM(2,ne)/BBM(3,ne)
C.......For later use in gas mixing        
              BBM(2,ne)=0.d0 !stores the concentration
            ENDIF
            XAB(1,ne)=0.d0 !proportion of volume change
          ENDDO
          CALL FLOW_SUMMATION(NBJ,NEELEM,NPNE,NVJE,NVJP,NXI,XP,ERROR,
     &      *9999)
c          IF(.NOT.FAST)THEN
C.......Calculate the shear stress using the time-averaged flows        
c          CALL TREESHEAR(NBJ,ncount,NEELEM,NORD,NPNE,NVJE,
c     &      torr_in_order(1,1,ncount),XP,filename,ERROR,
c     &      *9999)
c            IF(READ_VOLUMES)THEN
c              CALL READINITIALVOLUME(MECHANICS_FILETYPE,NBJ,NDP,NEELEM,
c     &          NENP,NPNE,NXI,BBM,CE,sumvolume,undef,XP,ZD,
c     &          FIRST_READ,ERROR,*9999)
c            ELSE
c              CALL RESETINITIALVOLUME(NEELEM,NXI,BBM,CE,sumvolume,
c     &          ERROR,*9999) !resets BBM with gravitational effect
c            ENDIF
            
C.......Calculate the specific ventilation following volume reset
c            DO noelem=1,NEELEM(0)
c              ne=NEELEM(noelem)
c              np=NPNE(2,1,ne)
c              IF(NXI(1,0,ne).EQ.0)THEN !terminal
c                XP(1,1,nj_radius+1,np)=XP(1,1,nj_flow,np)*T_interval
c     &            /BBM(1,ne)
c              ELSE
c                XP(1,1,nj_radius+1,np)=0.d0
c              ENDIF
c            ENDDO
c          ENDIF !FAST
          
          IF(.NOT.FAST)THEN
            CALL VOLUMEOFMESH(NBJ,NEELEM,NORD,NPNE,NVJE,NXI,MeanVolume,
     &        BBM,CE,XP,ERROR,*9999)
            INITIAL_VOLUME=CE(nm_vol_bel,1)
          ENDIF
c        CALL SHEARSTRESS(NBJ,NEELEM,NPNE,NVJE,XP,ERROR,*9999)
          
          IF(NORMALISE) CALL NORMALISEFLOWS(NBJ,NEELEM,NPNE,NXI,NVJE,XP,
     &      ERROR,*9999)
          
        ENDIF !setup
      ENDIF !uniform

      CALL EXITS('EVFLOW_DYNAM')
      RETURN
 9999 CALL ERRORS('EVFLOW_DYNAM',ERROR)
      CALL EXITS('EVFLOW_DYNAM')
      RETURN 1
      END
      
      
