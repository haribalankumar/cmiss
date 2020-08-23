      SUBROUTINE IPMOTI_LUNG_PV(NBJ,NEELEM,NELIST,NELIST2,NENP,
     &  NJJ_COEFF,NJJ_FLOW,NKJ,NKJE,NORD,NPNE,NPNODE,nr,NVJE,NVJP,NXI,
     &  CE,FRC,MINPCNT,TLC,XAB,XP,BELOW,LUMPED_PARAMETER,SET_FRC,
     &  SET_TLC,SCALE,VMIN_STATE,ERROR,*)

C#### Subroutine: IPMOTI_LUNG
C###  Description:
C###    IPMOTI_LUNG inputs motion parameters for pulmonary problems in
C###    region nr.


      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),NELIST(0:NEM),
     &  NELIST2(0:NEM),NENP(NPM,0:NEPM),NJJ_COEFF,NJJ_FLOW,
     &  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),NORD(5,NE_R_M),
     &  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M),nr,NVJE(NNM,NBFM,NJM,NEM),
     &  NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 CE(NMM,NEM),FRC,MINPCNT,TLC,XAB(NORM,NEM),
     &  XP(NKM,NVM,NJM,NPM)
      LOGICAL BELOW,LUMPED_PARAMETER,SET_FRC,SET_TLC,VMIN_STATE
      CHARACTER SCALE*(10),ERROR*(*)
!     Local Variables
      INTEGER i,M,N,nb,ne,ne0,ne1,ne2,
     &  NE_OLD(1000),NE_TEMP(1000),N_EFF,nj,NJJ_START,nk,nlpm,nn,noelem,
     &  nonode,np,nrr,NT_BNS,NTB_ACINI,
     &  NTB_ACINI_EFF,NUM_NODES,nv,new,cn,no_in_base,c1,
     &  no_constricted,ntop,nbottom
      REAL*8 FLOWPROP,
     &  length,pxi,radius,RATIO_SCALE,specific_ventilation,
     &  sum_pxi,SMAX,SMIN,SUM_ACINI,
     &  TOTAL_CHANGING_VOLUME,VMAX,VMEAN,VMIN,volume_below(NE_R_M),
     &  volume_dv(NE_R_M),volumes(NE_R_M),z,z_max,z_min,sum_newvol,
     &  sum_ce,a,b,c,d,P_nlpm,Vmaxi_acin, Percent_Vmax_acin,V_acin,
     &  P_max0,P_max1,P_max2,Vol_0,Vol_1,Ans0,Ans1,Ans2,FIXED_VOLUME,
     &  P_ng0,P_ng1,P_ng2,top,bottom,max_vol_top,max_vol_bottom
      LOGICAL SET_ZED
!     Functions
      REAL*8 LENGTH_1D

      CALL ENTERS('IPMOTI_LUNG_PV',*9999)

      CALC_INIT_VOL=.TRUE. !Use this in Pressure/Volume stuff
      TOTAL_LUNG_CAP=TLC
      FUNC_RES_CAP=FRC !set as a global variable

      IF(PRESSURE_DISTRIBUTION)THEN
        c=C_VEN
        d=D_VEN

        IF(ASET)THEN
          A_VEN=A_VEN
        ELSE
          A_VEN=5.0d0 !fixed percentage volume for min acinus
        ENDIF
        
        IF(BSET)THEN
          B_VEN=B_VEN
        ELSE
          B_VEN=95.0d0 !fixed percentage volume for min acinus
        ENDIF

        a=A_VEN
        b=B_VEN
        D_SPREAD=0.0d0 !modified 13/7/06
        CULM_DELP=0.0d0

      ENDIF!PRESSURE_DISTRIBUTION
      
      NELIST2(0)=0
      IF(BELOW)THEN
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          NT_BNS=1
          NE_OLD(1)=ne
          DO WHILE(NT_BNS.NE.0)
            NUM_NODES=NT_BNS
            NT_BNS=0
            DO M=1,NUM_NODES
              ne1=NE_OLD(M)
              DO N=1,NXI(1,0,ne1)
                NT_BNS=NT_BNS+1
                ne2=NXI(1,N,ne1)
                NE_TEMP(NT_BNS)=ne2
                NELIST2(0)=NELIST2(0)+1
                NELIST2(NELIST2(0))=ne2
              ENDDO !N
            ENDDO !M
            DO N=1,NT_BNS
              NE_OLD(N)=NE_TEMP(N)
            ENDDO !N
          ENDDO !WHILE
        ENDDO !noelem
      ELSE
        CALL ILIST_COPY(NELIST(0),NELIST(1),NELIST2(1))
        NELIST2(0)=NELIST(0)
      END IF !BELOW
      
      TERMINAL_NO=NELIST(0)
      
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        CE(1,ne)=0.d0 !initialise the percentage volume change
        XAB(22,ne)=0.0d0!initilise the percentage volume change
        volume_dv(ne)=0.d0
      ENDDO !noelem
   

      IF(FIRST_SET)THEN !initialise the number in base
        no_in_base=0
        SET_ZED=.TRUE.
        DO ne=1,NEELEM(0)
          XAB(17,ne)=42.0d0!cant be between 0,1 anything else is ok
          XAB(16,ne)=0.0d0
          XAB(18,ne)=0.0d0
        ENDDO !Initiliase XAB(16,nlpm) to 0.0
        FIRST_SET=.FALSE.
        ENDIF!FIRST SET

        IF(PRESSURE_DISTRIBUTION)THEN
c......Need to calculate this first time through
        z_min=1.0d6
        z_max=-1.0d6
        nlpm=0
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          nlpm=nlpm+1
          XAB(1,nlpm)=DBLE(ne)
          XAB(21,nlpm)=XAB(1,nlpm) !use because XAB(1,nlpm)is overwritten in PRF SOLN
          nb=NBJ(1,ne)
          np=NPNE(2,nb,ne)
          z=XP(1,1,3,np) !z position
c          write(*,*)'nlpm ',nlpm,' z ',z
          z_min=MIN(z_min,z)
          z_max=max(z_max,z)
        ENDDO !noelem
c      ENDIF !FIRST_SET
       
        DO nlpm=1,TERMINAL_NO
          ne=NINT(XAB(1,nlpm))
c          write(*,*)'XAB(19,ne) in IPMOTI_LUNG.f',XAB(19,ne)
          np=NPNE(2,nb,ne)
          z=XP(1,1,3,np) !z position
          pxi=(z-z_min)/(z_max-z_min)
          XAB(4,nlpm)=pxi
          XAB(14,nlpm)=pxi
        ENDDO !nlpm
        ENDIF
        
      IF(USE_ATELEM)THEN
        DO cn=1,CNC !For each node chosen
          new=ATELEM(cn)
C.........Loop through each of the subtended elements
          NT_BNS=1
          NE_OLD(1)=new
          DO WHILE(NT_BNS.NE.0)
            NUM_NODES=NT_BNS
            NT_BNS=0
            DO M=1,NUM_NODES
              ne1=NE_OLD(M)
              DO N=1,NXI(1,0,ne1)
                NT_BNS=NT_BNS+1
                ne2=NXI(1,N,ne1)
                NE_TEMP(NT_BNS)=ne2
                DO nlpm=1,TERMINAL_NO
                  IF(ne2.EQ.NINT(XAB(1,nlpm)))THEN
                    XAB(16,nlpm)=1.0d0! sets y/n condition
                  ENDIF
                ENDDO !TERMINAL_NO
              ENDDO !N
            ENDDO !M
            DO N=1,NT_BNS
              NE_OLD(N)=NE_TEMP(N)
            ENDDO !N
          ENDDO !WHILE
        ENDDO !cn
      ENDIF !(USE_ATELEM)


      no_constricted=0
      DO nlpm=1,TERMINAL_NO
        IF(XAB(16,nlpm).GT.0.0d0)THEN
          no_constricted=no_constricted+1
        ENDIF
      ENDDO

      IF(no_constricted.GT.0)THEN!have already set up XAB(16,nlpm) for ZED
        SET_ZED=.FALSE.
      ELSE
        SET_ZED=.TRUE.
      ENDIF
      
        IF(ZED_DISTRIBUTION.AND.SET_ZED)THEN
            no_in_base=0
c...........Go through all Terminal Elements and select if z< ZED
            DO nlpm=1,TERMINAL_NO
              IF(XAB(4,nlpm).LE.ZED)THEN
                no_in_base=no_in_base+1
                XAB(17,nlpm)=XAB(4,nlpm)!identifies the acinii
              ENDIF !If LE ZED HEIGHT IN LUNG
            ENDDO !TERMINAL NO
            
            c1=1
            DO nlpm=1,TERMINAL_NO
              IF(XAB(17,nlpm).NE.42.0d0)THEN !is in base
                IF(PATCHY_PERCENT(c1).EQ.1)THEN
                  XAB(16,nlpm)=1.0d0
                ENDIF
                c1=c1+1
                IF(c1.GT.10)THEN
                  c1=1
                ENDIF
              ENDIF
            ENDDO
            
          ENDIF !ZED_DISTRIBUTION (The nlpm is constricted if not 42.0 


          IF(RNDCOMP)THEN
            c1=1
            DO nlpm=1,TERMINAL_NO
              XAB(18,nlpm)=R_VALUE(c1)
              c1=c1+1
              IF(c1.GT.50)THEN!hardcoded to need an index of 50.
                c1=1
              ENDIF
            ENDDO
          ENDIF
          
   
C.....Calculate the current mesh volume      
      CALL MESH_VOLUME(NBJ,NEELEM,NORD,NPNE,NVJE,NXI,volumes,
     &  volume_below,XAB,XP,ERROR,*9999)
      WRITE(OP_STRING,'(''Initial volume is'',F8.3,''L'')')
     &  volume_below(1)/1.d6
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)


C.....Initialise changing volume
      DO noelem=1,NELIST2(0)
        ne=NELIST2(noelem)
        volume_dv(ne)=volumes(ne)
      ENDDO !noelem



C.....Add lumped parameter volume to volume of terminal elements
      DO nlpm=1,NTB
        ne=NINT(XAB(1,nlpm))
        volume_dv(ne)=volume_dv(ne)+XAB(2,nlpm)
      ENDDO



C.....Calculate the total changing volume      
      DO noelem=NEELEM(0),1,-1
        ne=NEELEM(noelem)
        ne0=NXI(-1,1,ne)
        IF(ne0.NE.0)THEN !not stem branch
          IF(NORD(5,ne).EQ.1)THEN !start of a 'half' branch
            volume_dv(ne0)=volume_dv(ne0)+2.0d0*volume_dv(ne)
          ELSE !within a tube branch
            volume_dv(ne0)=volume_dv(ne0)+volume_dv(ne)
          ENDIF !NORD(5)
        ENDIF !ne0
      ENDDO !noelem
      TOTAL_CHANGING_VOLUME=volume_dv(NEELEM(1))



C.....Update the proportion of tidal volume received
      DO noelem=1,NELIST2(0)
        ne=NELIST2(noelem)
        CE(1,ne)=volumes(ne)/TOTAL_CHANGING_VOLUME
      ENDDO !noelem

      
      DO nlpm=1,NTB
        XAB(7,nlpm)=XAB(2,nlpm)/TOTAL_CHANGING_VOLUME
      ENDDO

      INITIAL_VOLUME=volume_below(NEELEM(1))
      FIXED_VOLUME=INITIAL_VOLUME-TOTAL_CHANGING_VOLUME
      FIXED_VOLUMES=FIXED_VOLUME
     
      
      IF(SET_FRC)THEN
        CALL MESH_SET_FRC(NBJ,NEELEM,NELIST2,NORD,NPNE,NVJE,NXI,
     &    FIXED_VOLUME,FRC,volume_below,volumes,XAB,XP,LUMPED_PARAMETER,
     &    SCALE,ERROR,*9999)
      ELSE
        FRC=volume_below(NEELEM(1))
      ENDIF

c      write(*,*)'if this is last output then volume not converging'
   
     
      IF(VMIN_STATE)THEN
c        CALL MESH_SET_DIST(NBJ,NELIST,NE_OLD,NE_TEMP,NORD,NPNE,
c     &    NTB_ACINI,NTB_ACINI_EFF,NVJE,NXI,CE,FIXED_VOLUME,FRC,MINPCNT,
c     &    volume_below,volumes,XAB,XP,BELOW,ERROR,*9999)

      IF(NTB.GT.0)THEN
        VMEAN=(FRC-FIXED_VOLUME)/NTB
        VMIN=VMEAN*(MINPCNT/100.d0)
        sum_pxi=0.0d0
        z_min=1.0d6
        z_max=-1.0d6
        DO nlpm=1,NTB !calculate the vertical range
          ne=NINT(XAB(1,nlpm))
          nb=NBJ(1,ne)
          np=NPNE(2,nb,ne)
          z=XP(1,1,3,np) !z position
          z_min=MIN(z_min,z)
          z_max=MAX(z_max,z)
        ENDDO !nlpm

        DO nlpm=1,NTB !calculate normalized acini positions, and sum
          ne=NINT(XAB(1,nlpm)) !supplying element
          nb=NBJ(1,ne)
          np=NPNE(2,nb,ne) !end node of supplying element
          z=XP(1,1,3,np) !z position
          pxi=(z-z_min)/(z_max-z_min) !normalize the position
          XAB(4,nlpm)=pxi !stores normalized position of acinus
          XAB(14,nlpm)=pxi !also stores (and doesn't get changed in fe07)
          sum_pxi=sum_pxi+pxi !sum the normalized positions
        ENDDO !nlpm
      
C.......FRC = sum of initial volumes
C.......    = sum(1-pxi).VM<<<<<<< IPMOTI_LUNG.f
        VMAX=(FRC-FIXED_VOLUME-NTB*VMIN+sum_pxi*VMIN)/sum_pxi
        VMEAN=0.0d0
        DO nlpm=1,NTB
          ne=NINT(XAB(1,nlpm))
          nb=NBJ(1,ne)
          np=NPNE(2,nb,ne)
          pxi=XAB(4,nlpm)
          XAB(2,nlpm)=(1.d0-pxi)*VMIN+pxi*VMAX !calculate volume
          XAB(3,nlpm)=XAB(2,nlpm)
          VMEAN=VMEAN+XAB(2,nlpm)
        ENDDO !nlpm
        VMEAN=VMEAN/NTB

c The pressure volume stuff only works for alveolated airways

        
      ELSE !do for alveolated airways
        CALL ASSERT(BELOW,'>>Use BELOW terminal element list',ERROR,
     &    *9999)

 !nelist must contain the terminal conducting elements
 !proportion of TLC = element volume/total changing volume

c    Calculate the normalised position of the acinii       
        z_min=1.0d6
        z_max=-1.0d6
        VMEAN=0.0d0
        
        nlpm=0
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          nlpm=nlpm+1
          XAB(1,nlpm)=DBLE(ne)
          XAB(3,nlpm)=volume_below(ne)-volumes(ne) !store below parent branch = acinus
          nb=NBJ(1,ne)
          np=NPNE(2,nb,ne)
          z=XP(1,1,3,np) !z position
          z_min=MIN(z_min,z)
          z_max=max(z_max,z)
        ENDDO !noelem
          TERMINAL_NO=nlpm

c      Calculate the effective number of acinii          
        N_EFF=1
        NTB_ACINI=nlpm
        NTB_ACINI_EFF=0
        SUM_ACINI=0
        sum_pxi=0.0d0
        DO nlpm=1,NTB_ACINI
          ne=NINT(XAB(1,nlpm))
          np=NPNE(2,nb,ne)
          z=XP(1,1,3,np) !z position
          pxi=(z-z_min)/(z_max-z_min)
          XAB(4,nlpm)=pxi
 !Use this to determine where any given acinar unit is in a normalised lung
 !As is now determines the acinii closest to the isovolume point         
c          IF(pxi.LT.0.57d0)THEN
c            IF(pxi.GT.0.56d0)THEN
c              write(*,*)'nlpm (0.56-0.57 of pxi',nlpm,' pxi',pxi
c            ENDIF
c          ENDIF
          XAB(14,nlpm)=pxi
          ne0=ne
          N_EFF=1
          DO WHILE (ne0.NE.0)
            IF(NORD(5,ne0).EQ.1)N_EFF=N_EFF*2
            ne0=NXI(-1,1,ne0)
          ENDDO
          NTB_ACINI_EFF=NTB_ACINI_EFF+N_EFF
          sum_pxi=sum_pxi+pxi*N_EFF
        ENDDO !nlpm

        
        IF(PRESSURE_DISTRIBUTION) THEN
          Vmaxi_acin=(TLC-FIXED_VOLUME)/NTB_ACINI_EFF
          VMAX_ACIN=Vmaxi_acin
          VMIN_ACIN=0.10d0*VMAX_ACIN
          VMAX=VMAX_ACIN
          VMIN=VMIN_ACIN
          VMEAN=(FRC-FIXED_VOLUME)/NTB_ACINI_EFF
        ELSE  
          VMEAN=(FRC-FIXED_VOLUME)/NTB_ACINI_EFF
          VMIN=VMEAN*(MINPCNT/100.d0)
          VMAX=(FRC-FIXED_VOLUME-NTB_ACINI_EFF*VMIN+sum_pxi*VMIN)
     '      /sum_pxi
        ENDIF

       
        
        sum_newvol=0.0d0

        IF(PRESSURE_DISTRIBUTION)THEN
c.......This is sometimes failing to be bounded
          P_max0=P_MAX-3.5d0*(P_MAX-P_MIN)
          P_max1=P_MAX+3.5d0*(P_MAX-P_MIN)
          P_max2=P_MAX
          P_ng0=P_MIN-10.0d0
          P_ng1=P_MIN+10.0d0
          P_ng2=P_MIN
        ENDIF
      
c.......Cacluate the volumes for P_max0,P_max1
c        Pmax0
         SUM_ACINI=0.0d0 !set the sum of acinii to be zero at each iteration
         DO nlpm=1,NTB_ACINI
           
          IF(pxi.GT.0.5d0)THEN
            XAB(15,nlpm)=VMAX+0.0d0*(1.0d0-pxi) !can use this to distribute volumes
          ELSE
            XAB(15,nlpm)=VMAX-0.0d0*(0.5d0-pxi) !can use this to distribute volumes
          ENDIF
            ne=NINT(XAB(1,nlpm))
            np=NPNE(2,nb,ne)
            z=XP(1,1,3,np) !z position
            pxi=XAB(4,nlpm)

           
            IF(PRESSURE_DISTRIBUTION) THEN
              IF(NO_GRAV)THEN
                P_nlpm=P_ng0
              ELSE
                P_nlpm=(1.0d0-pxi)*P_MIN+pxi*P_max0
              ENDIF
              c=C_VEN+(1.0d0-pxi)*C_SPREAD
              c=c+c*XAB(16,nlpm)*CHANGE_COMPLIANCE+XAB(18,nlpm)*C_RND
              d=D_VEN+(1.0d0-pxi)*D_SPREAD
              Percent_Vmax_acin=a+(b/(1+EXP(-(P_nlpm-c)/d)))
              V_acin=Percent_Vmax_acin*XAB(15,nlpm)/100.0d0 !looks ok
c              write(*,*)'V_acin is',V_acin
              XAB(2,nlpm)=V_acin-V_acin*XAB(16,nlpm)*CONSTRICT/100.0d0
c              write(*,*)'XAB(2,',nlpm,') is ',XAB(2,nlpm)
              XAB(11,nlpm)=XAB(2,nlpm)
            ELSE!UNIFORM
              XAB(2,nlpm)=(1.d0-pxi)*VMIN+pxi*VMAX
              XAB(11,nlpm)= XAB(2,nlpm)
            ENDIF

c        write(*,*)XAB(2,15),'  ',XAB(2,22),'  ',XAB(2,31),'  ',XAB(2,83)
            
            sum_newvol=sum_newvol+XAB(2,nlpm)
           
            ne0=ne
            N_EFF=1
            DO WHILE (ne0.NE.0)
              IF(NORD(5,ne0).EQ.1) N_EFF=N_EFF*2
              ne0=NXI(-1,1,ne0)
            ENDDO
            SUM_ACINI=SUM_ACINI+XAB(2,nlpm)*N_EFF
          ENDDO

          
          
c........begins Bisection method for establishement of FRC in Pressure Distribution
          IF(PRESSURE_DISTRIBUTION)THEN
c.........Store Vol_0,Ans0          
          Vol_0=SUM_ACINI+FIXED_VOLUME
          Ans0=FRC-Vol_0
c          write(*,*)'Vol0',Vol_0
c.........For P_Max1          
          SUM_ACINI=0.0d0 !set the sum of acinii to be zero at each iteration
          DO nlpm=1,NTB_ACINI
            ne=NINT(XAB(1,nlpm))
            np=NPNE(2,nb,ne)
            z=XP(1,1,3,np) !z position
            pxi=XAB(4,nlpm)
              IF(NO_GRAV)THEN
                P_nlpm=P_ng1
              ELSE
                P_nlpm=(1.0d0-pxi)*P_MIN+pxi*P_max1
              ENDIF
              c=C_VEN+(1.0d0-pxi)*C_SPREAD
              c=c+c*XAB(16,nlpm)*CHANGE_COMPLIANCE+XAB(18,nlpm)*C_RND
              d=D_VEN+(1.0d0-pxi)*D_SPREAD
              Percent_Vmax_acin=a+(b/(1+EXP(-(P_nlpm-c)/d)))
              V_acin=Percent_Vmax_acin*XAB(15,nlpm)/100.0d0 !looks ok
              XAB(2,nlpm)=V_acin-V_acin*XAB(16,nlpm)*CONSTRICT/100.0d0
              XAB(11,nlpm)=XAB(2,nlpm)
            sum_newvol=sum_newvol+XAB(2,nlpm)
            ne0=ne
            N_EFF=1
            DO WHILE (ne0.NE.0)
              IF(NORD(5,ne0).EQ.1) N_EFF=N_EFF*2
              ne0=NXI(-1,1,ne0)
            ENDDO
            SUM_ACINI=SUM_ACINI+XAB(2,nlpm)*N_EFF
          ENDDO
c.........Store Vol1,Ans1          
          Vol_1=SUM_ACINI+FIXED_VOLUME
          Ans1=FRC-Vol_1
          SUM_ACINI=0.0d0
       
         
          DO WHILE(ABS(FRC-(SUM_ACINI+FIXED_VOLUME)).GT.1.0d-6)
            SUM_ACINI=0.0d0 !set the sum of acinii to be zero at each iteration
            DO nlpm=1,NTB_ACINI
              ne=NINT(XAB(1,nlpm))
              np=NPNE(2,nb,ne)
              z=XP(1,1,3,np) !z position
              pxi=XAB(4,nlpm)
                IF(NO_GRAV)THEN
                  XAB(13,nlpm)=P_ng2
                ELSE
                  XAB(13,nlpm)=(1.0d0-pxi)*P_MIN+pxi*P_max2
                ENDIF
                c=C_VEN+(1.0d0-pxi)*C_SPREAD
                c=c+C*XAB(16,nlpm)*CHANGE_COMPLIANCE+XAB(18,nlpm)*C_RND
                d=D_VEN+(1.0d0-pxi)*D_SPREAD
                Percent_Vmax_acin=a+(b/(1+EXP(-(XAB(13,nlpm)-c)/d)))
                V_acin=Percent_Vmax_acin*XAB(15,nlpm)/100.0d0 !looks ok
                XAB(2,nlpm)=V_acin-V_acin*XAB(16,nlpm)*CONSTRICT/100.0d0
                XAB(11,nlpm)=XAB(2,nlpm)
              sum_newvol=sum_newvol+XAB(2,nlpm)
              ne0=ne
              N_EFF=1
              DO WHILE (ne0.NE.0)
                IF(NORD(5,ne0).EQ.1) N_EFF=N_EFF*2
                ne0=NXI(-1,1,ne0)
              ENDDO
              SUM_ACINI=SUM_ACINI+XAB(2,nlpm)*N_EFF
            ENDDO
            Ans2=FRC-(SUM_ACINI+FIXED_VOLUME)
            IF(Ans1*Ans2.GT.0.0d0)THEN
              IF(NO_GRAV)THEN
                P_ng1=P_ng2
                P_ng0=P_ng0
              ELSE
                P_max1=P_max2
                P_max0=P_max0
              ENDIF
              Ans1=Ans2
              Ans0=Ans0
            ELSEIF(Ans2*Ans0.GT.0.0d0)THEN
              IF(NO_GRAV)THEN
                P_ng0=P_ng2
                P_ng1=P_ng1
              ELSE
                P_max0=P_max2
                P_max1=P_max1
              ENDIF
              Ans0=Ans2
              Ans1=Ans1
            ENDIF
            IF(NO_GRAV)THEN
              P_ng2=P_ng0+0.5d0*(P_ng1-P_ng0)
            ELSE
              P_max2=P_max0+0.5d0*(P_max1-P_max0)
            ENDIF
c          write(*,*)'Ans2',Ans2
          ENDDO !While loop for volume convergence
        IF(NO_GRAV)THEN
          write(*,*)'The lung pressure is',P_ng2
        ELSE
          P_MAX=P_max2
          write(*,*)'PMAX',P_MAX
        ENDIF
        ENDIF!Ends bisection method for finding Volume
    
          
          sum_newvol=0.0d0
          DO nlpm=1,NTB_ACINI
            ne=NINT(XAB(1,nlpm))
            np=NPNE(2,nb,ne)
            RATIO_SCALE=XAB(2,nlpm)/XAB(3,nlpm)
            length=LENGTH_1D(NBJ,ne,NPNE,NVJE,XP)
            volumes(ne)=volumes(ne)*RATIO_SCALE
            CE(1,ne)=volumes(ne)/SUM_ACINI
            XAB(20,ne)=volumes(ne)/SUM_ACINI
            radius=DSQRT(volumes(ne)/(length*PI))
            CALL RADIUS_1D_CHANGE(NBJ,ne,NPNE,NVJE,radius,XP,ERROR,
     &        *9999)
            
C...........Loop through each of the subtended elements in the acinus
            NT_BNS=1
            NE_OLD(1)=ne
            DO WHILE(NT_BNS.NE.0)
              NUM_NODES=NT_BNS
              NT_BNS=0
              DO M=1,NUM_NODES
                ne1=NE_OLD(M)
                DO N=1,NXI(1,0,ne1)
                  NT_BNS=NT_BNS+1
                  ne2=NXI(1,N,ne1)
                  NE_TEMP(NT_BNS)=ne2
                  length=LENGTH_1D(NBJ,ne2,NPNE,NVJE,XP)
                  volumes(ne2)=volumes(ne2)*RATIO_SCALE
                  radius=DSQRT(volumes(ne2)/(length*PI))
                  CALL RADIUS_1D_CHANGE(NBJ,ne2,NPNE,NVJE,radius,XP,
     &              ERROR,*9999)
                  CE(1,ne2)=volumes(ne2)/SUM_ACINI
                  XAB(20,ne2)=volumes(ne2)/SUM_ACINI
c                  CE(1,ne2)=RATIO_SCALE
c                  sum_newvol=sum_newvol+CE(1,ne2)
                  sum_newvol=sum_newvol+volumes(ne2)
                ENDDO !N
              ENDDO !M
              DO N=1,NT_BNS
                NE_OLD(N)=NE_TEMP(N)
              ENDDO !N
            ENDDO !WHILE
          ENDDO !nlpm  
        ENDIF !
      ENDIF !VMINSTATE
      

      CALL MESH_VOLUME(NBJ,NEELEM,NORD,NPNE,NVJE,NXI,volumes,
     &  volume_below,XAB,XP,ERROR,*9999)
      write(*,*) 'Initial_volume (check)',INITIAL_VOLUME

c     For uniform model output Acinar volumes here
c      write(*,*)'Initial volumes nlpm 15,22,31,83'
c      write(*,*)XAB(2,15),'  ',XAB(2,22),'  ',XAB(2,31),'  ',XAB(2,83)

      IF(SET_TLC)THEN
c        CALL MESH_SET_TLC(NELIST,NTB_ACINI,NTB_ACINI_EFF,NXI,CE,
c     &    FIXED_VOLUME,FRC,TLC,volume_below,volumes,XAB,BELOW,ERROR,
c     &    *9999)
      SMAX=-1.0d6
      SMIN= 1.0d6
      
      IF(NTB.GT.0)THEN
        DO nlpm=1,NTB
          XAB(7,nlpm)=((TLC-FIXED_VOLUME)/NTB-XAB(2,nlpm))/(TLC-FRC)
          specific_ventilation=XAB(7,nlpm)/XAB(2,nlpm)
          SMAX=MAX(SMAX,specific_ventilation)
          SMIN=MIN(SMIN,specific_ventilation)
        ENDDO !nlpm
      ELSE
        CALL ASSERT(BELOW,'>>Use BELOW terminal element list',ERROR,
     &    *9999)
C.......nelist must contain the terminal conducting elements
C.......proportion of TLC = element volume/total changing volume
C.......Calculate the real and effective number of acini
        nlpm=0
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          nlpm=nlpm+1
          XAB(1,nlpm)=DBLE(ne)
        ENDDO
        NTB_ACINI=NELIST(0)
        NTB_ACINI_EFF=0
        DO nlpm=1,NTB_ACINI
          ne=NINT(XAB(1,nlpm))
          ne0=ne
          N_EFF=1
          DO WHILE (ne0.NE.0)
            IF(NORD(5,ne0).EQ.1) N_EFF=N_EFF*2
            ne0=NXI(-1,1,ne0)
          ENDDO
          NTB_ACINI_EFF=NTB_ACINI_EFF+N_EFF
        ENDDO !nlpm
        nlpm=0
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          nlpm=nlpm+1
C.........store below parent branch = acinus
          XAB(3,nlpm)=volume_below(ne)-volumes(ne)
C.........proportion of TLC-FRC to this acinus
c          XAB(2,nlpm)=(TLC/NTB_ACINI_EFF-XAB(3,nlpm))/(TLC-FRC)
          XAB(2,nlpm)=((TLC-FIXED_VOLUME)/NTB_ACINI_EFF-XAB(3,nlpm))
     &      /(TLC-FRC)
        ENDDO !noelem
        sum_ce=0.0d0
        DO nlpm=1,NTB_ACINI
          ne=NINT(XAB(1,nlpm))
C.........Loop through each of the subtended elements in the acinus
          NT_BNS=1
          NE_OLD(1)=ne
          DO WHILE(NT_BNS.NE.0)
            NUM_NODES=NT_BNS
            NT_BNS=0
            DO M=1,NUM_NODES
              ne1=NE_OLD(M)
              DO N=1,NXI(1,0,ne1)
                NT_BNS=NT_BNS+1
                ne2=NXI(1,N,ne1)
                NE_TEMP(NT_BNS)=ne2
                IF(UNIFORM_DISTRIBUTION)THEN
                  CE(1,ne2)=volumes(ne2)/XAB(3,nlpm)*XAB(2,nlpm)
c                  write(*,*)CE(1,ne2)
                ENDIF
                specific_ventilation=CE(1,ne2)/volumes(ne2)
                SMAX=MAX(SMAX,specific_ventilation)
                SMIN=MIN(SMIN,specific_ventilation)
                sum_ce=sum_ce+CE(1,ne2)
              ENDDO !N
            ENDDO !M
            DO N=1,NT_BNS
              NE_OLD(N)=NE_TEMP(N)
            ENDDO !N
          ENDDO !WHILE
        ENDDO !nlpm
      ENDIF
      ENDIF

C     Set up a new field to store the flow
      IF(NJJ_FLOW.NE.0.OR.NJJ_COEFF.NE.0)THEN
        IF(NJJ_FLOW.NE.0)THEN !still to put in njj_flow correctly
          IF(NJ_LOC(NJL_FIEL,NJJ_FLOW,nr).GT.0)THEN
            nj_flow=NJ_LOC(NJL_FIEL,NJJ_FLOW,nr)
          ELSE
            NJJ_START=NJ_LOC(NJL_FIEL,0,nr)
C nj_flow        
            nj=NJ_LOC(NJL_FIEL,NJJ_START,nr)+1
            CALL ASSERT(nj.LE.NJM,' >>Increase NJM',ERROR,*9999)
            NJ_LOC(NJL_FIEL,NJJ_START+1,nr)=nj
            NJ_TYPE(nj,1)=NJL_FIEL
            NJ_TYPE(nj,2)=NJJ_START+1
            IF(nj.GT.NJ_LOC(0,0,nr)) NJ_LOC(0,0,nr)=nj
            IF(nj.GT.NJ_LOC(NJL_FIEL,0,0)) NJ_LOC(NJL_FIEL,0,0)=nj
            nj_flow=nj
            NJ_LOC(NJL_FIEL,0,nr)=NJJ_START+1
            DO nrr=1,NRT
              IF(NJ_LOC(0,0,nrr).GT.NJ_LOC(0,0,0))
     &          NJ_LOC(0,0,0)=NJ_LOC(0,0,nrr)
            ENDDO !nrr
          ENDIF
        ENDIF
        IF(NJJ_COEFF.NE.0)THEN !still to put in njj_coeff correctly
C nj_coefficient
          NJJ_START=NJ_LOC(NJL_FIEL,0,nr)
          nj=NJ_LOC(NJL_FIEL,NJJ_START,nr)+1
          CALL ASSERT(nj.LE.NJM,' >>Increase NJM',ERROR,*9999)
          NJ_LOC(NJL_FIEL,NJJ_START+1,nr)=nj
          NJ_TYPE(nj,1)=NJL_FIEL
          NJ_TYPE(nj,2)=NJJ_START+1
          IF(nj.GT.NJ_LOC(0,0,nr)) NJ_LOC(0,0,nr)=nj
          IF(nj.GT.NJ_LOC(NJL_FIEL,0,0)) NJ_LOC(NJL_FIEL,0,0)=nj
          nj_coeff=nj
          NJ_LOC(NJL_FIEL,0,nr)=NJJ_START+1
          DO nrr=1,NRT
            IF(NJ_LOC(0,0,nrr).GT.NJ_LOC(0,0,0))
     &        NJ_LOC(0,0,0)=NJ_LOC(0,0,nrr)
          ENDDO !nrr
          CALL ASSERT(NJ_LOC(NJL_FIEL,0,nr).LE.NJ_LOC_MX,
     &      '>>Increase dimension of PROMPT_NV',ERROR,*9999)
        
          DO noelem=1,NEELEM(0)
            ne=NEELEM(noelem)
            nb=NBJ(1,ne)
            NBJ(nj_flow,ne)=nb !basis function the same as for geometry
            NBJ(nj_coeff,ne)=nb
            DO nn=1,2
              np=NPNE(nn,nb,ne)
              DO i=1,NENP(np,0)
                ne2=NENP(np,i)
                IF(ne2.EQ.ne)THEN
                  NVJE(nn,nb,nj_flow,ne)=i !the version of node at nn'th position
                  NVJE(nn,nb,nj_coeff,ne)=i
                ENDIF
              ENDDO !i
              DO nk=1,NKT(nn,nb)
                NKJE(nk,nn,nj_flow,ne)=1
                NKJE(nk,nn,nj_coeff,ne)=1
              ENDDO !nk
            ENDDO !nn
          ENDDO !noelem
          DO nonode=1,NPNODE(0)
            np=NPNODE(nonode)
            NVJP(nj_flow,np)=NENP(np,0) !# of versions = # of adjacent elements
            NKJ(nj_flow,np)=1 !# of derivatives = 0, vlaue=1
            NVJP(nj_coeff,np)=NENP(np,0) !# of versions = # of adjacent elements
            NKJ(nj_coeff,np)=1 !# of derivatives = 0, vlaue=1
          ENDDO !np
        ENDIF
    
        DO nonode=1,NPNODE(0)
          np=NPNODE(nonode)
          IF(NJJ_FLOW.NE.0) NVJP(nj_flow,np)=NENP(np,0) !#vers = #adj elems
          IF(NJJ_COEFF.NE.0) NVJP(nj_coeff,np)=NENP(np,0)
          IF(NJJ_FLOW.NE.0) NKJ(nj_flow,np)=0 !# of derivatives = 0
          IF(NJJ_COEFF.NE.0) NKJ(nj_coeff,np)=0
        ENDDO !np
      
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          nb=NBJ(1,ne)
          IF(NJJ_FLOW.NE.0) NBJ(nj_flow,ne)=nb !bf same as geometry
          IF(NJJ_COEFF.NE.0) NBJ(nj_coeff,ne)=nb
          DO nn=1,NNT(nb)
            np=NPNE(nn,nb,ne)
            DO i=1,NENP(np,0)
              ne2=NENP(np,i)
              IF(ne2.EQ.ne)THEN
                IF(NJJ_FLOW.NE.0) NVJE(nn,nb,nj_flow,ne)=i !vers node at nn'th pos
                IF(NJJ_COEFF.NE.0) NVJE(nn,nb,nj_coeff,ne)=i
              ENDIF
            ENDDO !i
            IF(NJJ_FLOW.NE.0)  NKJE(1,nn,nj_flow,ne)=1
            IF(NJJ_COEFF.NE.0) NKJE(1,nn,nj_coeff,ne)=1
          ENDDO !nn
          DO nn=1,NNT(nb)
            np=NPNE(nn,nb,ne)
            nv=NVJE(nn,nb,nj_flow,ne)
            IF(NJJ_FLOW.NE.0)  XP(1,nv,nj_flow,np)=0.d0
            IF(NJJ_COEFF.NE.0) XP(1,nv,nj_coeff,np)=0.d0
          ENDDO !nn
        ENDDO !noelem
        DO nonode=1,NPNODE(0)
          np=NPNODE(nonode)
          IF(NJJ_FLOW.NE.0)THEN
            NVJP(nj_flow,np)=NENP(np,0) 
            NKJ(nj_flow,np)=1
          ENDIF
          IF(NJJ_COEFF.NE.0) NVJP(nj_coeff,np)=NENP(np,0)
          IF(NJJ_COEFF.NE.0) NKJ(nj_coeff,np)=1
        ENDDO !np
      ENDIF
C.....Recalculate the mesh volume
      CALL MESH_VOLUME(NBJ,NEELEM,NORD,NPNE,NVJE,NXI,volumes,
     &  volume_below,XAB,XP,ERROR,*9999)
 

C.....Sum the flow proportions to check = 1      
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        XAB(3,ne)=CE(1,ne)
      ENDDO !noelem
      DO nlpm=1,NTB
        ne=NINT(XAB(1,nlpm))
        XAB(3,ne)=XAB(7,nlpm)
      ENDDO
    
      DO noelem=NEELEM(0),1,-1
        ne=NEELEM(noelem)
        ne0=NXI(-1,1,ne)
        IF(ne0.NE.0)THEN !not stem branch
          IF(NORD(5,ne).EQ.1)THEN !start of a 'half' branch
            XAB(3,ne0)=XAB(3,ne0)+2.d0*XAB(3,ne)
          ELSE !within a tube branch
            XAB(3,ne0)=XAB(3,ne0)+XAB(3,ne)
          ENDIF !NORD(5)
        ENDIF !ne0
      ENDDO !noelem

      FLOWPROP=XAB(3,1) !proportion of flow below stem
      
c      write(*,*)'Initial volumes(Pressure Distribution)nlpm 15,22,31,83'
c       write(*,*)XAB(2,15),'  ',XAB(2,22),'  ',XAB(2,31),'  ',XAB(2,83)
c     write out all the acinii valomues 
c       DO nlpm=1,TERMINAL_NO
c         write(*,*)XAB(2,nlpm)
c       ENDDO
      IF(PRESSURE_DISTRIBUTION)THEN
        top=0.0d0
        bottom=0.0d0
        ntop=0
        nbottom=0
        DO nlpm=1,TERMINAL_NO
          IF(XAB(4,nlpm).GT.0.80d0)THEN
            top=top+XAB(11,nlpm)
            ntop=ntop+1
          ENDIF
          IF(XAB(4,nlpm).LT.0.2d0)THEN
            bottom=bottom+XAB(11,nlpm)
            nbottom=nbottom+1
          ENDIF
        ENDDO
        max_vol_top=ntop*VMAX_ACIN
        max_vol_bottom=nbottom*VMAX_ACIN
        
        write(*,*)'volume top 20% ',top,'   vol bottom 20% ',bottom
        write(*,*)'ntop',ntop,' nbottom',nbottom
        write(*,*)'max_vol_top',max_vol_top,' max_bottom',
     &    max_vol_bottom    
      ENDIF !PRESSURE_DISTIBUTION
       
c      write(*,*)'FLOWPRP IS',FLOWPROP
c     update VOLD to use in MESH_FLOW
      IF(PRESSURE_DISTRIBUTION)THEN
        V_OLD=INITIAL_VOLUME-FIXED_VOLUME
      ELSE
        V_OLD=FRC !do so that example_91 runs with PULMAT(2).NE.0.0d0
      ENDIF
      
c      write(*,*)'INITIAL_VOLUME',INITIAL_VOLUME
c      write(*,*)'FIXED_VOLUME',FIXED_VOLUME
c      write(*,*)'V_OLD',V_OLD all ok
c        DO nlpm=1,TERMINAL_NO
c        write(*,*)'XAB(2,nlpm)',XAB(2,nlpm),' nlpm',nlpm
c        ENDDO !seems to be working for ATELEM XAB(16,nlpm) 1 if True 0 if False
c      This is now used to set up the Pressure _Volume stuff
  

      
      IF(BELOW)THEN
        WRITE(OP_STRING,
     &    '(''Scaled volume is'',F8.3,''L for '',I7,''acini'')'
     &    )INITIAL_VOLUME/1.d6,MAX(NTB_ACINI,NTB)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''Equivalent to'',I7,'' acini'')')
     &    MAX(NTB_ACINI_EFF,NTB)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ELSE
        WRITE(OP_STRING,'(''Scaled volume is '',F8.3,''L'')')
     &    INITIAL_VOLUME/1.d6
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
c      IF(NTB.GT.0)THEN
      IF(SET_FRC)THEN
        write(*,*)'VMEAN',VMEAN
          WRITE(OP_STRING,'(''Lumped model volumes (mm^3): min ='',F8.3,
     &      '' max = '',F8.3,'' mean ='',F8.3)') VMIN,VMAX,VMEAN
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(SET_TLC)THEN
          WRITE(OP_STRING,'(''Specific ventilation to TLC: min ='',D9.3,
     &      '' max = '',D9.3)') SMIN,SMAX
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
c      ENDIF
      WRITE(OP_STRING,'(''Flow percentage is'',F12.3)') FLOWPROP*100.d0
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL EXITS('IPMOTI_LUNG_PV')
      RETURN
 9999 CALL ERRORS('IPMOTI_LUNG_PV',ERROR)
      CALL EXITS('IPMOTI_LUNG_PV')
      RETURN 1
      END


