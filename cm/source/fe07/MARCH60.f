      SUBROUTINE MARCH60(ISC_GKK,ISR_GKK,NBH,NBJ,NEELEM,NELIST,NENFVC,
     '  NENP,NFVC,NHE,NHP,NHQ,NKJE,NKH,NODENVC,NODENVCB,NPLIST,NPNE,
     '  NPNODE,NPNY,NQNY,nr,NRE,NRLIST,NRLIST2,NVCB,NVCNODE,NVHP,NVJE,
     '  nx,NXI,NYNE,NYNP,NYNR,NYQNR,GKK,GRR,SE,VC,VC_INIT,XNFV,XO,XP,YP,
     '  YQ,YQS,ZA,ZNFV,ZP,ERROR,*)

C#### Subroutine: MARCH60
C###  Description:
C###    MARCH60 is used to solve fe60 fluid dynamics problems.
C###    This method is based on an adaption of MAC-like explicit
C###    methods to a collocated unstructured Voronoi mesh. Rhie-Chow
C###    like interpolation is used to avoid the problem of spurious
C###    pressure oscillations. CJ Were is the original implementor
C###    of the algorithm, now adapted by RG Boyes for moving mesh
C###    flow problems in CMISS.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'marc00.cmn'
      INCLUDE 'time01.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER ISC_GKK(NISC_GKKM),ISR_GKK(NISR_GKKM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NENFVC(0:NFVCM,NFVM),NENP(NPM,0:NEPM,0:NRM),
     '  NFVC(2,0:NFVCM,NVCM),NHE(NEM),NHP(NPM),
     '  NHQ(NRM),NKJE(NKM,NNM,NJM,NEM),NKH(NHM,NPM,NCM),NODENVC(NVCM),
     '  NODENVCB(NVCBM),NPLIST(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM),NQNY(2,NYQM,0:NRCM),nr,NRE(NEM),
     '  NRLIST(0:NRM),NRLIST2(0:NRM),NVCB(-1:3,NVCBM),NVCNODE(2,NP_R_M),
     '  NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),
     '  nx,NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),NYQNR(0:NYQM,0:NRCM,NCM,0:NRM)
      REAL*8 GKK(NZ_GKK_M),GRR(NOM),SE(NSM,NBFM,NEM),VC(0:NVCM),
     '  VC_INIT(2,NVCM),XNFV(-(NJM+1):NJM,NFVM),XO(NOM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM),YQ(NYQM,NIQM,NAM),
     '  YQS(NIQSM,NQM),ZA(NAM,NHM,NCM,NEM),ZNFV(NFVM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER HIST_file_counter,na,NIQLIST(0:16),NIQSLIST(0:16),
     '  NIYLIST(0:16),nh,nhx,
     '  nonode,np,NSTEP,NUMTIMEDATA,ny,!SMAR009 22/12/98 ,nfv
     '  nvc,iy,nj,njj,N_OUTLET,MAXLOC_FACES_OUTLET,
     '  FACETOT_OUTLET,noelem,ne,nb,nn
      INTEGER*4 NVCB_OUTLET_PTR,XNFV_OUTLET_PTR,NFVC_OUTLET_PTR,
     '  NVCBNODE_OUTLET_PTR,OUTLET_MATRIX_PTR,IPIVOT_OUTLET_PTR,
     '  OUTLET_RHS_PTR,OLDVC_PTR,ZNFVMSH_PTR
      REAL*8 T,YPMAX(16),YPMIN(16),AZERO(1)
      REAL ELAPSED_TIME,TIME_START1(1),TIME_START2(1),TIME_STOP(1)
      CHARACTER FILEFORMAT*6,CHAR3*10
      LOGICAL CONTINUE,ENDFILE,OUTPUT,YPDATA,YQDATA,YQSDATA,FIRST_A,
     '  UPDATE_MATRIX,x_INIT,MEM_INIT,REFACT_OUTLET,VALVES_SWITCHED

      SAVE NSTEP
      CALL ENTERS('MARCH60',*9999)

      DT=TINCR
      AZERO(1)=0.d0
      IF(BINTIMEFILE.GT.0) THEN
        FILEFORMAT='BINARY'
      ELSE
        FILEFORMAT='ASCII'
      ENDIF !bintimefile
      IF(RESTART) THEN
        ERROR='>>Not implemented yet'
        GOTO 9999
      ELSE IF(.NOT.RESTART) THEN !perform initial tasks

C       ..Timing information
        CALL CPU_TIMER(CPU_USER,TIME_START1)
        CALL CPU_TIMER(CPU_USER,TIME_START2)

C       ..Require storage of old Voronoi cell volumes if mesh moving
C         and mesh flux storage
        IF(.NOT.MESHFIXD) THEN
          OLDVC_PTR=0
          ZNFVMSH_PTR=0
          MEM_INIT=.TRUE.
          CALL ALLOCATE_MEMORY((NVCM+1),1,DPTYPE,OLDVC_PTR,
     '      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(NFVM,1,DPTYPE,ZNFVMSH_PTR,
     '      MEM_INIT,ERROR,*9999)
        ELSE
          OLDVC_PTR=0
          ZNFVMSH_PTR=0
          MEM_INIT=.TRUE.
          CALL ALLOCATE_MEMORY(1,1,DPTYPE,OLDVC_PTR,
     '      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(1,1,DPTYPE,ZNFVMSH_PTR,
     '      MEM_INIT,ERROR,*9999)
        ENDIF

C       ..Some checking of array sizes
        CALL ASSERT(NOM.GE.NVCT,'>>Increase NOM >= NVCM',ERROR,*9999)
        IF(.NOT.MESHFIXD) THEN
          CALL ASSERT(NVM.GE.2,'>>Increase NVM >= 2',ERROR,*9999)
        ENDIF

C       ..Form Poisson matrix if mesh is static
        CALL FORMPPEM(ISC_GKK,ISR_GKK,NFVC,nr,NVCNODE,nx,
     '    GKK,XNFV,ERROR,*9999)
        IF(.NOT.CONFINED.AND.DUCTFLOW) THEN
          NVCB_OUTLET_PTR=0
          NFVC_OUTLET_PTR=0
          XNFV_OUTLET_PTR=0
          NVCBNODE_OUTLET_PTR=0
          OUTLET_MATRIX_PTR=0
          IPIVOT_OUTLET_PTR=0
          OUTLET_RHS_PTR=0
          MEM_INIT=.TRUE.
          CALL OUTLET_MEM(FACETOT_OUTLET,MAXLOC_FACES_OUTLET,NBJ,
     '      NENP,NPLIST,NPNE,NPNODE,nr,NVCB,NVCNODE,N_OUTLET,
     '      ERROR,*9999)
          CALL ALLOCATE_MEMORY(N_OUTLET,1,INTTYPE,NVCB_OUTLET_PTR,
     '      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(2*(MAXLOC_FACES_OUTLET+1)*N_OUTLET,1,
     '      INTTYPE,NFVC_OUTLET_PTR,MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(5*FACETOT_OUTLET,1,DPTYPE,
     '      XNFV_OUTLET_PTR,MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(NP_R_M,1,INTTYPE,NVCBNODE_OUTLET_PTR,
     '      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(N_OUTLET**2,1,DPTYPE,OUTLET_MATRIX_PTR,
     '      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(N_OUTLET,1,DPTYPE,OUTLET_RHS_PTR,
     '      MEM_INIT,ERROR,*9999)
          CALL ALLOCATE_MEMORY(N_OUTLET,1,INTTYPE,IPIVOT_OUTLET_PTR,
     '      MEM_INIT,ERROR,*9999)
          CALL FORMPPEM_OUTLET(FACETOT_OUTLET,
     '      MAXLOC_FACES_OUTLET,NBJ,NENP,%VAL(NFVC_OUTLET_PTR),NPLIST,
     '      NPNE,NPNODE,nr,NVCB,%VAL(NVCBNODE_OUTLET_PTR),
     '      %VAL(NVCB_OUTLET_PTR),NVCNODE,N_OUTLET,
     '      %VAL(OUTLET_MATRIX_PTR),%VAL(XNFV_OUTLET_PTR),
     '      XP,ZA,ERROR,*9999)
        ENDIF

C       ..Initialize matrix factorization parameters
        FIRST_A=.TRUE.
        UPDATE_MATRIX=.TRUE.
        x_INIT=.TRUE.
        REFACT_OUTLET=.TRUE.

C       ..Initialize history file stuff
        YPDATA=.TRUE.
        YQDATA=.FALSE.
        YQsDATA=.FALSE.
        IF(.NOT.MESHFIXD) THEN
          NIYLIST(0)=2
          NIYLIST(1)=1
          NIYLIST(2)=5
        ELSE
          NIYLIST(0)=1
          NIYLIST(1)=1
        ENDIF
        NIQLIST(0)=0
        NIQLIST(1)=0
        NIQSLIST(0)=0
        NRLIST(0)=1
        NRLIST(1)=nr
        NRLIST2(0)=1
        NRLIST2(1)=nr
        na=1
        IF(HIST_file_intervals.GT.0) THEN !history file output
          CALL IOHIST(IOFILE1,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,
     '      NQNY,NRLIST,NRLIST2,NUMTIMEDATA,nx,NYNR,NYQNR,T,YP,YPMAX,
     '      YPMIN,YQ,YQS,'WRITE',FILEFORMAT,FILE02,'OPEN',ENDFILE,
     '      .TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)
          HIST_file_counter=0 !initialise counter
        ENDIF !history file output
      ENDIF !.not.restart

C     ..Initialize time loop and arrays
      T=TSTART
      NSTEP=0
      CONTINUE=.TRUE.
      CALL DCOPY(NFVT,AZERO,0,ZNFV,1)
      DV=1.d0

C     ..Initialize moving mesh variables
      IF(.NOT.MESHFIXD) THEN
        SYSTOLE=.FALSE.
        VALVES_SWITCHED=.FALSE.
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            XP(1,2,nj,np)=0.d0
          ENDDO
        ENDDO
        DO ny=1,NYM
          YP(ny,5)=0.d0
        ENDDO
        ! Initialise meshflux
        CALL DCOPY(NFVT,AZERO,0,%VAL(ZNFVMSH_PTR),1)
      ENDIF

C     ..Enter the main time loop
      DO WHILE(CONTINUE)

        NSTEP=NSTEP+1

C       ..Calculate the timestep if automatic stepping
        IF(KTYP23.EQ.2) THEN
          CALL TIMESTEP(NFVC,NODENVC,NPNODE,nr,NVCNODE,nx,NYNP,XNFV,
     '      YP,ERROR,*9999)
C         ..This step is to create a very small time step if there
C           are no specified velocities anywhere
          IF(.NOT.MESHFIXD.AND.NSTEP.EQ.1) DT=1.d-4
        ENDIF

        T=T+DT

C       ..check if current time is < finish time, and field change is
C         greater than a tolerance
        IF(T.LE.TFINISH.AND.DV.GT.DVTOL) THEN
          IF(IWRIT1(nr,nx).EQ.0) THEN
            OUTPUT=.FALSE.
          ELSE
            IF(MOD(NSTEP,IWRIT1(nr,nx)).EQ.0)
     '        THEN
              OUTPUT=.TRUE.
            ELSE
              OUTPUT=.FALSE.
            ENDIF
          ENDIF !iwrit1

C ..Do all The Voronoi stuff in this section

C         ..Move nodes and calculate mesh flux if mesh deforms
          IF(.NOT.MESHFIXD) THEN

C           Time dependent boundary conditions
            IF(KTYP3_init(nx).EQ.2) THEN
C             ..User_10 operates the valves - see feuser.f
              CALL USER_10(NPLIST,NPNODE,nr,NVCB,NVCNODE,T,
     '          VALVES_SWITCHED,ERROR,*9999)
            ENDIF
            IF(VALVES_SWITCHED) THEN

              IF(.NOT.CONFINED.AND.DUCTFLOW) THEN
                CALL FREE_MEMORY(NVCB_OUTLET_PTR,ERROR,*9999)
                CALL FREE_MEMORY(NFVC_OUTLET_PTR,ERROR,*9999)
                CALL FREE_MEMORY(XNFV_OUTLET_PTR,ERROR,*9999)
                CALL FREE_MEMORY(NVCBNODE_OUTLET_PTR,ERROR,*9999)
                CALL FREE_MEMORY(IPIVOT_OUTLET_PTR,ERROR,*9999)
                CALL FREE_MEMORY(OUTLET_MATRIX_PTR,ERROR,*9999)
                CALL FREE_MEMORY(OUTLET_RHS_PTR,ERROR,*9999)
                NVCB_OUTLET_PTR=0
                NFVC_OUTLET_PTR=0
                XNFV_OUTLET_PTR=0
                NVCBNODE_OUTLET_PTR=0
                OUTLET_MATRIX_PTR=0
                IPIVOT_OUTLET_PTR=0
                OUTLET_RHS_PTR=0
                MEM_INIT=.TRUE.
                CALL OUTLET_MEM(FACETOT_OUTLET,MAXLOC_FACES_OUTLET,NBJ,
     '            NENP,NPLIST,NPNE,NPNODE,nr,NVCB,NVCNODE,
     '            N_OUTLET,ERROR,*9999)
C     SMAR009 12/01/99 removed ,NODENVC from list
                CALL ALLOCATE_MEMORY(N_OUTLET,1,INTTYPE,NVCB_OUTLET_PTR,
     '            MEM_INIT,ERROR,*9999)
                CALL ALLOCATE_MEMORY(2*(MAXLOC_FACES_OUTLET+1)*N_OUTLET,
     '            1,INTTYPE,NFVC_OUTLET_PTR,MEM_INIT,ERROR,*9999)
                CALL ALLOCATE_MEMORY(5*FACETOT_OUTLET,1,DPTYPE,
     '            XNFV_OUTLET_PTR,MEM_INIT,ERROR,*9999)
                CALL ALLOCATE_MEMORY(NP_R_M,1,INTTYPE,
     '            NVCBNODE_OUTLET_PTR,MEM_INIT,ERROR,*9999)
                CALL ALLOCATE_MEMORY(N_OUTLET**2,1,DPTYPE,
     '            OUTLET_MATRIX_PTR,MEM_INIT,ERROR,*9999)
                CALL ALLOCATE_MEMORY(N_OUTLET,1,DPTYPE,OUTLET_RHS_PTR,
     '            MEM_INIT,ERROR,*9999)
                CALL ALLOCATE_MEMORY(N_OUTLET,1,INTTYPE,
     '            IPIVOT_OUTLET_PTR,MEM_INIT,ERROR,*9999)
                CALL FORMPPEM_OUTLET(FACETOT_OUTLET,
     '            MAXLOC_FACES_OUTLET,NBJ,NENP,%VAL(NFVC_OUTLET_PTR),
     '            NPLIST,NPNE,NPNODE,nr,NVCB,%VAL(NVCBNODE_OUTLET_PTR),
     '            %VAL(NVCB_OUTLET_PTR),NVCNODE,N_OUTLET,
     '            %VAL(OUTLET_MATRIX_PTR),%VAL(XNFV_OUTLET_PTR),
     '            XP,ZA,ERROR,*9999)
                REFACT_OUTLET=.TRUE.
              ENDIF
              VALVES_SWITCHED=.FALSE.
            ENDIF

C           ..Initialize node velocities
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,njj,nr)
                XP(1,2,nj,np)=0.d0
              ENDDO
            ENDDO

C           ..User_11 moves the boundaries - see feuser.f

            CALL USER_11(NODENVCB,NPNODE,nr,NVCB,NVCNODE,T,XP,
     '        ERROR,*9999)
            CALL MVINTNOD(NFVC,NODENVC,NPNODE,nr,NVCNODE,XNFV,XP,
     '        ERROR,*9999)
            CALL MESHFLUX(NFVC,NODENVC,NPNODE,nr,XNFV,XP,
     '        %VAL(ZNFVMSH_PTR),ERROR,*9999)
          ENDIF


C         ..Calculate net flux of momentum into each Voronoi cell
          CALL NETTFLUX(NFVC,NODENVC,NPNODE,nr,NVCB,NVCNODE,nx,
     '      NYNP,VC,XNFV,YP,ZNFV,%VAL(ZNFVMSH_PTR),ERROR,*9999)

          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(MARCH60_1)
            WRITE(OP_STRING,'(/$,'' Nett Fluxes:'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' ############'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            iy=1
            CALL YPZP(iy,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NVHP,nx,NYNE,
     '        NYNP,YP,ZA,ZP,ERROR,*9999)
            DO nvc=1,NVCT
              nonode=NODENVC(nvc)
              np=NPNODE(nonode,nr)
              WRITE(OP_STRING,'('' Voronoi cell: '',I7,
     '          ''  (u,v,w): '',3d16.6)')
     '          nvc,(ZP(1,1,nh_loc(nhx,nx),np,1),nhx=1,nh_loc(0,nx)-1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !nvc
CC$OMP       END CRITICAL(MARCH60_1)
          ENDIF !dop
C
C         ..Reconnect mesh in the event of deformation
          IF(.NOT.MESHFIXD) THEN
            CALL RECONNECT(NBJ,NEELEM,NELIST,NENP,NKJE,NPLIST,NPNE,
     '        NPNODE,nr,NRE,NVJE,NXI,SE,XP,ZA,ERROR,*9999)

C           ..Check the new mesh for boundary violations and
C             mesh quality
            IF(MESHCHECK) THEN
              CALL CHECKMSH(NBJ,NEELEM,NFVC,NPLIST,NPNE,NPNODE,nr,
     '          NVCNODE,VC,ERROR,*9999)
            ENDIF
            IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP         CRITICAL(MARCH60_2)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                nb=NBJ(1,ne)
                WRITE(OP_STRING,'(/,''Delaunay simplex:  '',I12)') ne
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(''Circumcentre:      '',3D12.4)')
     '            (ZA(1,NJ_LOC(NJL_GEOM,nj,nr),1,ne),
     '            nj=1,NJ_LOC(NJL_GEOM,0,nr))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(''Connections:'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(''============'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(''  Node  Vers   Opp      Coord'//
     '            'inates'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                DO nn=1,NNT(nb)
                  WRITE(OP_STRING,'(I6,I6,I6,''    '',3D12.4)')
     '              NPNE(nn,nb,ne),NVJE(nn,nb,1,ne),
     '              NXI(0,nn,ne),
     '              (XP(1,NVJE(nn,nb,1,ne),NJ_LOC(NJL_GEOM,nj,nr),
     '              NPNE(nn,nb,ne)),nj=1,NJ_LOC(NJL_GEOM,0,nr))
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDDO
              ENDDO !noelem
CC$OMP         END CRITICAL(MARCH60_2)
            ENDIF

C           ..Swap the volume storage
            CALL DCOPY(NVCT,VC,1,%VAL(OLDVC_PTR),1)

            CALL CALC_VORO(NBJ,NENFVC,NENP,NFVC,NODENVC,NPLIST,
     '        NPNE,NPNODE,nr,NXI,VC,VC_INIT,XNFV,XP,ZA,ERROR,*9999)
            CALL DCOPY(NFVT,AZERO,0,ZNFV,1)
            CALL DCOPY(NFVT,AZERO,0,%VAL(ZNFVMSH_PTR),1)
          ENDIF

C         ..Calculate tentative velocities
          CALL TENTVELO(NODENVC,NPNODE,nr,nx,NYNP,%VAL(OLDVC_PTR),
     '      VC,YP,ERROR,*9999)

          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(MARCH60_3)
            WRITE(OP_STRING,'(/$,'' Tentative velocities:'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' #####################'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            iy=1
            CALL YPZP(iy,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NVHP,nx,NYNE,
     '        NYNP,YP,ZA,ZP,ERROR,*9999)
            DO nvc=1,NVCT
              nonode=NODENVC(nvc)
              np=NPNODE(nonode,nr)
              WRITE(OP_STRING,'('' Voronoi cell: '',I7,''  (u,v,w): '',
     '          3d16.6)') nvc,
     '          (ZP(1,1,nh_loc(nhx,nx),np,1),nhx=1,nh_loc(0,nx)-1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !nvc
CC$OMP       END CRITICAL(MARCH60_3)
          ENDIF !dop

C         ..Need to recompute Poisson equation matrix if mesh changes
          IF(.NOT.MESHFIXD) THEN
            CALL FORMPPEM(ISC_GKK,ISR_GKK,NFVC,nr,NVCNODE,nx,
     '        GKK,XNFV,ERROR,*9999)
c     SMAR009 12/01/99 removed  ,NODENVC from list
            UPDATE_MATRIX=.TRUE.
            FIRST_A=.TRUE.
          ENDIF

C         ..Calculate RHS = divergence of u divided by time step
          IF(.NOT.MESHFIXD) THEN
C           ..If mesh is moving, then need to calculate the divergence
C             of the fluid out of each cell for the *NEW* mesh.
            CALL MESHFLUX(NFVC,NODENVC,NPNODE,nr,XNFV,XP,
     '        %VAL(ZNFVMSH_PTR),ERROR,*9999)
          ENDIF
          CALL CALCDIVU(NFVC,NODENVC,NPNODE,nr,NVCB,NVCNODE,nx,
     '      NYNP,GRR,XNFV,YP,%VAL(ZNFVMSH_PTR),ERROR,*9999)
          IF(OUTPUT) THEN
            WRITE(OP_STRING,'(/,'' Divergence magnitude before '//
     '        'pressure correction:'',D16.8)') DIVMAG
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

C         ..Solve Pressure Correction Equation
          CALL SOLVPPEQ(FACETOT_OUTLET,%VAL(IPIVOT_OUTLET_PTR),ISC_GKK,
     '      ISR_GKK,MAXLOC_FACES_OUTLET,NFVC,%VAL(NFVC_OUTLET_PTR),
     '      NODENVC,NODENVCB,NPLIST,NPNODE,nr,NVCB,
     '      %VAL(NVCBNODE_OUTLET_PTR),%VAL(NVCB_OUTLET_PTR),
     '      NVCNODE,nx,NYNP,N_OUTLET,GKK,GRR,%VAL(OUTLET_MATRIX_PTR),
     '      %VAL(OUTLET_RHS_PTR),XNFV,
     '      %VAL(XNFV_OUTLET_PTR),XO,YP,FIRST_A,REFACT_OUTLET,
     '      UPDATE_MATRIX,x_INIT,ERROR,*9999)

C         ..If only the first step, & the matrix is fixed ,
C           then the matrix only has to be factorised once.
          IF(NSTEP.EQ.1) THEN
            IF(MESHFIXD) THEN
              FIRST_A=.FALSE.
              UPDATE_MATRIX=.FALSE.
            ENDIF
          ENDIF

C         ..Calculate Rhie Chow fluxes
          CALL RHIECHOW(NFVC,NODENVC,NPNODE,nr,NVCB,
     '      NVCNODE,nx,NYNP,XNFV,YP,ZNFV,%VAL(ZNFVMSH_PTR),
     '      ERROR,*9999)
          IF(OUTPUT) THEN
            WRITE(OP_STRING,'(/,'' Divergence magnitude after '//
     '        'pressure correction:'',D16.8)') DIVMAG
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

C         ..Correct the velocities to get the new velocity field
          CALL CORRVELO(NFVC,NODENVC,NPNODE,nr,NVCB,
     '      %VAL(NVCBNODE_OUTLET_PTR),NVCNODE,nx,NYNP,N_OUTLET,
     '      %VAL(OUTLET_RHS_PTR),VC,XNFV,YP,ERROR,*9999)

C ..End of Voronoi computations

C         ..write history file output
          IF(HIST_file_intervals.GT.0) THEN !history file output
            HIST_file_counter=HIST_file_counter+1
            IF(HIST_file_counter.EQ.HIST_file_intervals) THEN

C ..This section is used to write nodal positions (time-dependent)
C ..to the history file. I have used iy = 5, and assumed that
C ..dependent variables have the same indexing as geometric. Bit
C ..of a hack, but it should suit the purpose well. Ends justifying
C ..means ??
              IF(.NOT.MESHFIXD) THEN
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
C                 need some inter-relation between nj and nh ??
                  DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                    nj=NJ_LOC(NJL_GEOM,njj,nr)
                    ny=NYNP(1,1,nj,np,0,1,nr)
                    YP(ny,5)=XP(1,1,nj,np)
                  ENDDO
                ENDDO
              ENDIF
              CALL IOHIST(IOFILE1,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY,NQNY,NRLIST,NRLIST2,NUMTIMEDATA,nx,NYNR,NYQNR,
     '          T,YP,YPMAX,YPMIN,YQ,YQS,'WRITE',FILEFORMAT,FILE02,
     '          'TIME_DATA',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '          ERROR,*9999)
              HIST_file_counter=0
            ENDIF !counter
          ENDIF !history file output

C         ..write output
          IF(OUTPUT) THEN
            WRITE(OP_STRING,'(/,'' I = '',I7,'' T'//
     '        ' = '',D13.6,'' DT = '',D13.6,'' DV = '',D13.6)')
     '        NSTEP,T,DT,DV
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            IF(.NOT.CONFINED) THEN
              WRITE(OP_STRING,'('' Influx = '',D11.4,'' Outflux'//
     '          ' = '',D11.4,'' Volume swept = '',D11.4,'' Error'//
     '          '= '',D11.4)') INFLUX,OUTFLUX,VSWEPT,
     '          (INFLUX+OUTFLUX+VSWEPT)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            IF(.NOT.MESHFIXD) THEN
              WRITE(OP_STRING,'('' Total Volume = '',D11.4)')
     '          VC(0)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP         CRITICAL(MARCH60_4)
              iy=8
              CALL YPZP(iy,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NVHP,nx,
     '          NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
              DO nvc=1,NVCT
                nonode=NODENVC(nvc)
                np=NPNODE(nonode,nr)
                WRITE(OP_STRING,'(/$,''Voronoi cell: '',I6,'' (np = '',
     '            I6,'')'')') nvc,np
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(''============================='//
     '            '====='')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(''Cell centre position: '',3D16.8)')
     '            (XP(1,1,NJ_LOC(NJL_GEOM,njj,nr),np),
     '            njj=1,NJ_LOC(NJL_GEOM,0,nr))
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                DO nhx=1,nh_loc(0,nx)
                  IF(nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN
                    WRITE(CHAR3,'(A10)') 'pressure p'
                  ELSE
                    IF(nhx.EQ.1) THEN
                      WRITE(CHAR3,'(A10)') 'velocity u'
                    ELSEIF(nhx.EQ.2) THEN
                      WRITE(CHAR3,'(A10)') 'velocity v'
                    ELSEIF(nhx.EQ.3) THEN
                      WRITE(CHAR3,'(A10)') 'velocity w'
                    ENDIF
                  ENDIF
                  nh=nh_loc(nhx,nx)
                  WRITE(OP_STRING,'('''//CHAR3(1:10)//' = '',D16.8)')
     '              ZP(1,1,nh,np,1)
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDDO !nhx
              ENDDO !nvc
CC$OMP         END CRITICAL(MARCH60_4)
            ENDIF !dop
          ENDIF !output
        ELSE
          CONTINUE=.FALSE.
        ENDIF !T
      ENDDO !time loop

C ..End main time loop

C     ..Calculate elapsed time
      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)

C     ..Do some output
      WRITE(OP_STRING,'(/$,'' Solution information'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' ===================='')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C     ..Check to see if solution reached steady state
      IF(T.LE.TFINISH) THEN
        WRITE(OP_STRING,'('' Solution reached steady state.'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Time taken to reach steady '//
     '    'state..... '',D14.6,'' secs.'')') T
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Number of iterations.......'//
     '    '.......... '',I14)') NSTEP
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Magnitude change in velocit'//
     '    'y field... '',D14.6)') DV
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' CPU time taken.............'//
     '    '.......... '',D14.6,'' secs.'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ELSE
        WRITE(OP_STRING,'('' Solution didnt reach steady state.'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Flow still unsteady at time'//
     '    '.......... '',D14.6,'' secs.'')') T
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Number of iterations.......'//
     '    '.......... '',I14)') NSTEP
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/$,'' CPU time taken.............'//
     '    '.......... '',D14.6,'' secs.'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

C     .. Free the memory used for the outlet solver ..
      IF(.NOT.CONFINED.AND.DUCTFLOW) THEN
        CALL FREE_MEMORY(NVCB_OUTLET_PTR,ERROR,*9999)
        CALL FREE_MEMORY(NFVC_OUTLET_PTR,ERROR,*9999)
        CALL FREE_MEMORY(XNFV_OUTLET_PTR,ERROR,*9999)
        CALL FREE_MEMORY(NVCBNODE_OUTLET_PTR,ERROR,*9999)
        CALL FREE_MEMORY(IPIVOT_OUTLET_PTR,ERROR,*9999)
        CALL FREE_MEMORY(OUTLET_MATRIX_PTR,ERROR,*9999)
        CALL FREE_MEMORY(OUTLET_RHS_PTR,ERROR,*9999)
      ENDIF
      CALL FREE_MEMORY(OLDVC_PTR,ERROR,*9999)
      CALL FREE_MEMORY(ZNFVMSH_PTR,ERROR,*9999)

C     ..Free the tetrahedral linked list
      IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN
        CALL ALE_destroy()
      ENDIF

      IF(HIST_file_intervals.GT.0) THEN !close history file
        CALL IOHIST(IOFILE1,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,
     '    NQNY,NRLIST,NRLIST2,NUMTIMEDATA,nx,NYNR,NYQNR,T,YP,YPMAX,
     '    YPMIN,YQ,YQS,'CLOSE',FILEFORMAT,FILE02,'  ',ENDFILE,
     '    .TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)
      ENDIF !history file output

      CALL CLOSEF(IOFILE2,ERROR,*9999)
      CALL CLOSEF(IOFILE4,ERROR,*9999)
      CALL EXITS('MARCH60')
      RETURN

 9999 IF(HIST_file_intervals.GT.0) THEN !close history file
        CALL IOHIST(IOFILE1,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,
     '    NQNY,NRLIST,NRLIST2,NUMTIMEDATA,nx,NYNR,NYQNR,T,YP,YPMAX,
     '    YPMIN,YQ,YQS,'CLOSE',FILEFORMAT,FILE02,'  ',ENDFILE,
     '    .TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*1111)
      ENDIF !history file output
 1111 CALL ERRORS('MARCH60',ERROR)
      CALL EXITS('MARCH60')
      RETURN 1
      END


