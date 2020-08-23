      SUBROUTINE INIT_2_11(NBJ,NEELEM,NELIST,NENP,NPLIST,NPNE,NPNODE,nr,
     &  NTIME_POINTS,NVJE,nx,NXI,NYNE,NYNP,BBM,CE,CP,XP,YP,BC,FIX,
     &  FLOW_FIXED,PRESSURE_FIXED,ERROR,*)

C#### Subroutine: INIT_2_11
C###  Description:
C###    INIT_2_11 sets initial conditions and boundary conditions for
C###    ITYP5(nr,nx)=2 (time integration) and ITYP2(nr,nx)=11 (pulmonary
C###    transport).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn'

!     Parameter List      
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),NELIST(0:NEM),
     &  NENP(NPM,0:NEPM,0:NRM),NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),
     &  NPNODE(0:NP_R_M),nr,NTIME_POINTS(NTIMEVARSM),
     &  NVJE(NNM,NBFM,NJM,NEM),nx,NXI(-NIM:NIM,0:NEIM,0:NEM),
     &  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM)
      REAL*8 BBM(2,NEM),CE(NMM,NEM),CP(NMM,NPM),XP(NKM,NVM,NJM,NPM),
     &  YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL BC,FIX(NYM,NIYFIXM),FLOW_FIXED,PRESSURE_FIXED
!     Local Variables
      INTEGER i,nb,ne,nh,nn,noelem,nonode,np,nv,ny,ny1,ny2,N3CO
      REAL*8 Cwall,Twall
!     Functions      
      REAL*8 CONC_EVAL,RFROMC
      LOGICAL CBBREV

      CALL ENTERS('INIT_2_11',*9999)

      NTIME_POINTS(1)=0
      nh=NH_LOC(1,nx)

      IF(.NOT.ADD)THEN
        DO ny=1,NYM
          FIX(ny,1)=.FALSE.
        ENDDO !ny
      ENDIF
      
      IF(ITYP3(nr,nx).EQ.1)THEN !inert gas mixing
        MAX_SMOOTH_CONC = MAX(CONC_INIT,CONC_IN)
        FIRST_BREATH=.TRUE.
        IF(.NOT.BC)THEN
          IF(.NOT.ADD)THEN
C       Initial concentrations
            DO nonode=1,NPNODE(0)
              np=NPNODE(nonode)
              ny=NYNP(1,1,nh,np,1,1)
              YP(ny,1)=CONC_INIT
              YP(ny,3)=CONC_INIT
              YP(ny,4)=CONC_INIT
C AJS 02/2011: Uncommenting the following lines to assign the initial
C concentration to lumped acini (BBM). This shouldn't affect mass
C calculations for geometries w/o BBMs b/c they will have zero volume.
              DO nn=1,NENP(np,0,nr)
                ne=NENP(np,nn,nr)
                IF(NXI(1,0,ne).EQ.0.AND.ne.NE.INLET_ELEMENT.AND.
     &            ne.NE.0)THEN
                  BBM(2,ne)=CONC_INIT
                  nb=NBJ(nj_radius,ne)
                  nv=NVJE(nn,nb,nj_radius,ne)
C AJS 02/2011 Calculate partial pressure at terminal (acinar) nodes
                  IF(ITYP7(nr,nx).GT.1.AND.nv.NE.0)THEN
                    XP(1,nv,nj_poa,np)=CONC_INIT*(O2molVol*1.d3)*
     &                (P_ATM-p_water)
                  ENDIF
                ENDIF
              ENDDO !nn
            ENDDO !noelem
          ELSE
            DO nonode=1,NPLIST(0)
              np=NPLIST(nonode)
              ny=NYNP(1,1,nh,np,1,1)
              YP(ny,1)=CONC_INIT
              YP(ny,3)=CONC_INIT
              YP(ny,4)=CONC_INIT
            ENDDO !noelem
          ENDIF
          IF(ITYP7(nr,nx).GT.1)
     &       WRITE(*,'('' >> NOTE: Airway concentrations must be '//
     &       'consistent with gas exchange units (i.e. mmol/mm3)'')')
          IDEAL_MASS(nx)=CONC_INIT*INITIAL_VOLUME

        ELSE IF(BC)THEN
          IF(.NOT.FLUX_BC)THEN
C       Fixed boundary conditions        
            DO nonode=1,NPLIST(0)
              np=NPLIST(nonode) !fixed node
              
              ny=NYNP(1,1,nh,np,1,1) !for fixed BC
              FIX(ny,1)=.TRUE. !set fixed BC
              YP(ny,1)=CONC_IN 
              YP(ny,3)=CONC_IN
              YP(ny,4)=CONC_IN
              
              ny=NYNP(1,1,nh,np,0,2) !for flux BC
              FIX(ny,1)=.FALSE. !set fixed BC
            ENDDO
          ELSE IF(FLUX_BC)THEN
C       Flux boundary conditions        
            DO nonode=1,NPLIST(0)
              np=NPLIST(nonode)
              ny=NYNP(1,1,nh,np,1,1)
              FIX(ny,1)=.FALSE.
              ny=NYNP(1,1,nh,np,0,2)
              FIX(ny,1)=.TRUE.
              YP(ny,3)=0.d0
            ENDDO
          ENDIF
        ENDIF
        
      ELSE IF(ITYP3(nr,nx).EQ.2)THEN !water and heat transfer
        IF(.NOT.ADD)THEN
 !do for first node
          ne=NEELEM(1)
          nb=NBJ(1,ne)
          np=NPNE(1,nb,ne)
          ny=NYNP(1,1,nh,np,1,1)
          Twall=CORE_TEMP
          Cwall=CONC_EVAL(Twall)
          YP(ny,1)=Twall
          YP(ny+1,1)=Cwall
          CP(1,np)=XP(1,1,nj_radius,np)
          CP(2,np)=0.007d0 !initial airway surface liquid depth
          CP(3,np)=Cwall
          DO i=4,7
            CP(i,np)=Twall
          ENDDO
          DO i=2,7
            CP(i+6,np)=CP(i,np)
          ENDDO
          
          DO noelem=1,NEELEM(0)
            ne=NEELEM(noelem)
            nb=NBJ(1,ne)
            np=NPNE(2,nb,ne)
            ny=NYNP(1,1,nh,np,1,1)
            Twall=CORE_TEMP
            Cwall=CONC_EVAL(Twall)
            YP(ny,1)=Twall
            YP(ny+1,1)=Cwall
            CP(1,np)=XP(1,1,nj_radius,np)
            CP(2,np)=0.007d0 !initial airway surface liquid depth
            CP(3,np)=Cwall
            DO i=4,7
              CP(i,np)=Twall
            ENDDO
            DO i=2,7
              CP(i+6,np)=CP(i,np)
            ENDDO
          ENDDO !noelem(ne)
        ENDIF !ADD
        
        IF(.NOT.FLUX_BC)THEN
C       Fixed boundary conditions        
          DO nonode=1,NPLIST(0)
            np=NPLIST(nonode)
            ny=NYNP(1,1,nh,np,1,1)
            YP(ny,1)=TEMP_IN
            YP(ny+1,1)=AH_IN
            FIX(ny,1)=.TRUE. !set fixed temp/humidity at mouth
            FIX(ny+1,1)=.TRUE. !set fixed temp/humidity at mouth
          ENDDO !nonode
        ELSE IF(FLUX_BC)THEN
C       Flux boundary conditions        
          DO nonode=1,NPLIST(0)
            np=NPLIST(nonode)
            ny=NYNP(1,1,nh,np,1,1)
            FIX(ny,1)=.FALSE.
            FIX(ny+1,1)=.FALSE.
            ny=NYNP(1,1,nh,np,0,2)
            FIX(ny,1)=.TRUE.
            FIX(ny+1,1)=.TRUE.
            YP(ny,3)=0.d0
          ENDDO
        ENDIF

      ELSE IF(ITYP3(nr,nx).EQ.4)THEN !simple P-R-F
        IF(.NOT.ADD)THEN
          DO ny=1,NYM
            FIX(ny,1)=.FALSE.
          ENDDO !ny
C       Initial values, default everything to zero currently
          DO ny=1,NYM
            YP(ny,1)=0.d0
          ENDDO !ny
          IF(CBBREV(CO,'TRANSPULMONARY',4,noco+1,NTCO,N3CO))THEN
            ptrans=RFROMC(CO(N3CO+1))
          ELSE
            ptrans=-5.d0
C MHT new
            DO noelem=1,NEELEM(0)
              ne=NEELEM(noelem)
              CE(20,ne)=ptrans
            ENDDO
          ENDIF
        ENDIF
        
        IF(.NOT.FLUX_BC)THEN
C       Fixed boundary conditions        
          IF(PRESSURE_FIXED)THEN
            DO nonode=1,NPLIST(0)
              np=NPLIST(nonode) !fixed node
              ny1=NYNP(1,1,nh,np,1,1) !for fixed pressure BC
c              ny2=NYNP(1,1,nh,np,0,2) !for flux BC
              FIX(ny1,1)=.TRUE. !set fixed
c              WRITE(*,*) ny1,ny2
c              FIX(ny2,1)=.FALSE. !set flux
              YP(ny1,1)=0.d0
            ENDDO !nonode
          ELSE IF(FLOW_FIXED)THEN
            DO noelem=1,NELIST(0)
              ne=NELIST(noelem)
              ny1=NYNE(1,nh,0,1,ne) !fixed
              ny2=NYNE(1,nh,0,2,ne) !flux
              FIX(ny1,1)=.TRUE. !set fixed
              FIX(ny2,1)=.FALSE. !set flux
              YP(ny1,1)=0.d0
            ENDDO
          ENDIF
        ELSE IF(FLUX_BC)THEN
C       Flux boundary conditions
C haven't done this yet          
        ENDIF
C...    If pulmonary circulation coupled via LPM capillary model 
C...    (Define: Palv, kc, Hd) 
        IF((COMPLIANCE.EQ.3.OR.COMPLIANCE.EQ.6).
     &    AND.COUPLE_VIA_LPM.EQ.'Y') THEN 
          IF(CBBREV(CO,'PALV',2,noco+1,NTCO,N3CO))THEN             
            Palv=RFROMC(CO(N3CO+1))
          ELSE
            Palv=0.d0
          ENDIF
          IF(CBBREV(CO,'KC',2,noco+1,NTCO,N3CO))THEN
            kc=RFROMC(CO(N3CO+1))
          ELSE
            kc=22.5d0
          ENDIF
          IF(CBBREV(CO,'HD',2,noco+1,NTCO,N3CO))THEN
            Hd=RFROMC(CO(N3CO+1))
          ELSE
            Hd=0.4d0
          ENDIF          
        ENDIF !COMPLIANCE.EQ.3
      ENDIF !problem type
      
      CALL EXITS('INIT_2_11')
      RETURN
 9999 CALL ERRORS('INIT_2_11',ERROR)
      CALL EXITS('INIT_2_11')
      RETURN 1
      END


