      SUBROUTINE GNBDMESH_CONTINUOUS(LD,LD_NP,nb,NBJ,NEELEM,NELIST2,
     &  NENP,NE_OLD,NE_TEMP,NKJ,NKJE,NNB,NP_INTERFACE,NPNE,NPNODE,nr,
     &  NRE,NVJE,NVJP,NXI,POINT_LIMIT,SCALE_DIST_LIMIT,angle_max,
     &  fraction,length_limit,min_length,SE,WD,XP,ZD,ERROR,*)

C####  Subroutine: GNBDMESH
C###   Description:
C###     GNBDMESH uses a bifurcating distributive algorithm to fill a host
C###     volume with a tree-like mesh. Details of the algorithm can be
C###     found in Tawhai et al. Ann. Biomed. Eng. 28(7):793-802 (2000). 

C#### Variable: LD_NP(NDM)
C###  Type: INTEGER
C###  Set_up: GNBDMESH
C###  Description:
C###  LD_NP(NDM) stores the single mesh node that is associated with a
C###  unique data point at the end of a volume-filling tree generation.
C###  This is used when coupling parenchymal soft tissue mechanics to
C###  airflow in the tree mesh (see UPVERTEX for application).

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
!     Parameter values
      INTEGER LD(NDM),LD_NP(NDM),nb,NBJ(NJM,NEM),
     &  NEELEM(0:NE_R_M,0:NRM),NELIST2(0:NEM),NENP(NPM,0:NEPM,0:NRM),
     &  NE_OLD(NE_R_M),NE_TEMP(NE_R_M),NKJ(NJM,NPM),
     &  NKJE(NKM,NNM,NJM,NEM),NNB(4,4,4,NBFM),NP_INTERFACE(0:NPM,0:3),
     &  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),
     &  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     &  POINT_LIMIT
      REAL*8 angle_max,fraction,length_limit,min_length,
     &  SCALE_DIST_LIMIT,SE(NSM,NBFM,NEM),WD(NJM,NDM),XP(NKM,NVM,NJM,
     &  NPM),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER i,LD_NUM,LD_TEMP(NDM),M,N,nd,nde1,nde2,
     '  nd_min,ne,NE_REACTIVATE(0:NEM),N_ELM,N_ELM_TEMP,
     &  ne_min,ne1,ne2,neb,nes,ngen,nj,nk,nn,noelem,
     &  noelem_gen,noelem_initial,nonode,np,np1,np2,np3,
     &  NSTEM(NE_R_M,2),NT_BNS,numtb,numzero
      REAL*8 COFM(3),length,LENGTH_PARENT,DIST,
     &  DISTANCE_LIMIT,MIN_DIST,XP1(3)
      LOGICAL BRANCH,SS(2)
      
      CALL ENTERS('GNBDMESH_CONTINUOUS',*9999)

      CALL ASSERT(NDT.GT.0,'>>Define data as seed points',ERROR,*9999)

C     Initialise NE_OLD (stores parent elements for current generation)
C     to the list of parent elements in NELIST2, and NT_BNS (current
C     number of terminal branches) to the number of parent branches.
      DO N=1,NELIST2(0) !for each parent element in list
        NE_OLD(N)=NELIST2(N)
        NSTEM(NELIST2(N),1)=NELIST2(N)
        NSTEM(NELIST2(N),2)=0
      ENDDO !N
      NT_BNS=NELIST2(0) !initial number of 'terminal' branches
      N_ELM=NT_BNS
      
C     Initialise LD array to map each seed point to a parent branch.
C     For a single parent, all seed points will initially be mapped to
C     it; for multiple parents MESH_REPOINT is called to calculate
C     the closest parent end-point to each seed point.
      DO nd=1,NDT !for every data point
        LD_NP(nd)=LD(nd) !temp storage of LD 
        LD(nd)=NELIST2(1) !assigns data points to the single parent
        DO nj=1,NJT !set weighting factor all=1.0
          WD(nj,NDT)=1.0d0
        ENDDO !nj
      ENDDO !nd

      IF(NELIST2(0).GT.1)THEN
        CALL MESH_REPOINT(LD,nb,N_ELM,NE_OLD,NE_REACTIVATE,
     '    NPNE,500.d0,XP,ZD,.TRUE.,ERROR,*9999)
        DO nd=1,NDT
          IF(LD(nd).NE.0)THEN
            ne_min=LD(nd)
            NSTEM(ne_min,2)=NSTEM(ne_min,2)+1
            LD_TEMP(nd)=LD(nd)
          ENDIF !LD
        ENDDO !nd

C       to get things right if an initial branch has no seed points
        NT_BNS=0
        DO N=1,NELIST2(0) !for each parent element in list
          IF(NSTEM(NELIST2(N),2).NE.0)THEN
            NT_BNS=NT_BNS+1
            NE_OLD(NT_BNS)=NELIST2(N)
          ENDIF
        ENDDO !N
        N_ELM=NT_BNS
        N_ELM_TEMP=N_ELM

C...    KSB 23/06/05 removing terminal elements (only 1 data pt
C...    associated) and corresponding data point from the group
        DO N=1,N_ELM
          ne_min=NE_OLD(N)
          IF(NSTEM(ne_min,2).LT.2)THEN  !terminal branch -remove data pt
            DO nd=1,NDT
              IF(LD(nd).EQ.ne_min)THEN
                LD_TEMP(nd)=LD(nd)
                LD(nd)=0
                NE_OLD(N)=0
                N_ELM_TEMP=N_ELM_TEMP-1
              ENDIF
            ENDDO !nd
          ENDIF
        ENDDO
        numtb=0
        DO N=1,N_ELM
          IF(NE_OLD(N).EQ.0)THEN
            I=0
            DO WHILE((N+I.LT.N_ELM).AND.(NE_OLD(N+I).EQ.0))
              I=I+1
            ENDDO
            DO M=N,N_ELM-I
              NE_OLD(M)=NE_OLD(M+I)
            ENDDO !M
            numtb=numtb+1
          ENDIF !NE_OLD
        ENDDO !N
        N_ELM=N_ELM_TEMP !removes terminal elements from loop
        
c        DO N=1,N_ELM
c          ne_min=NE_OLD(N)
c          IF(NSTEM(ne_min,2).LT.10)THEN !find closest 10 points
c            np2=NPNE(2,nb,ne_min)
c            DO i=1,10
c              ND_C(i)=0
c            ENDDO !i
c            DO i=1,10
c              MIN_DIST=1.d10
c              DO nd_count=1,NDT
c                nd=LD_TEMP(nd_count)
c                IF(nd.NE.0)THEN
c                  DIST=DSQRT((ZD(1,nd)-XP(1,1,1,np2))**2.d0+(ZD(2,nd)
c     '              -XP(1,1,2,np2))**2.d0+(ZD(3,nd)-XP(1,1,3,np2))
c     '              **2.d0)
c                  IF(DIST.LT.MIN_DIST)THEN
c                    TEST=.TRUE.
c                    DO j=1,10
c                      IF(nd.EQ.ND_C(j)) TEST=.FALSE.
c                    ENDDO !i
c                    IF(TEST)THEN
c                      ND_C(i)=nd
c                      MIN_DIST=dist
c                    ENDIF
c                  ENDIF !DIST
c                ENDIF !nd.NE.0
c              ENDDO !nd_count
c            ENDDO !i
c            DO i=1,4
c              LD(ND_C(i))=ne_min
c              LD_TEMP(ND_C(i))=0 !can't be assigned to another branch
c            ENDDO !i

c            N_ELM=N_ELM+1
c            NE_OLD(N_ELM)=ne_min
c          ENDIF !NSTEM
c        ENDDO !N
      ENDIF !NELIST2.GT.1
      
      LD_NUM=0
      DO nd=1,NDT
        IF(LD(nd).GT.0) LD_NUM=LD_NUM+1
      ENDDO !nd

      WRITE(OP_STRING,'('' gen   #brn  total# #term  #data'')')
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      
C     Set initial values for local and global nodes and elements
      ne=NET(0) !initialise mesh global element #
      np=NPT(0) !initialise mesh global node #
      noelem=NEELEM(0,nr) !initialise local element #
      nonode=NPNODE(0,nr) !initialise local node #
      noelem_initial=noelem
      
      ne1=NELIST2(1) !first parent element
      ngen=3
      ne_min=ne
      
C.... BIFURCATING DISTRIBUTIVE ALGORITHM
      DO WHILE(N_ELM.NE.0) !while still some parent branches
        ngen=ngen+1 !increment the generation from the parent
        NT_BNS=N_ELM 
        N_ELM=0
        numtb=0
        numzero=0
        noelem_gen=0
        DO M=1,NT_BNS
          ne1=NE_OLD(M) !parent global element #
C........ Calculate COFM of each seed point set
          CALL MESH_COFM(LD,ne1,COFM,ZD,ERROR,*9999)
          ne2=NXI(-1,1,ne1) !grandparent global element #
          CALL ASSERT(ne2.GT.0,
     &      '>>Require more than 2 parent generations',ERROR,*9999)
          np1=NPNE(2,nb,ne1) !parent global end node #
          np2=NPNE(1,nb,ne1) !parent global start node #
          np3=NPNE(1,nb,ne2) !grandparent global start node #
          LENGTH_PARENT=0.d0
          DO nj=1,NJT
            LENGTH_PARENT=LENGTH_PARENT+(XP(1,1,nj,np1)-XP(1,1,nj,np2))
     '        **2
          ENDDO !nj
          LENGTH_PARENT=DSQRT(LENGTH_PARENT)
C........ Split each set of seed points using the plane defined by the
C........ parent branch and the centre of mass.
          CALL MESH_SPLIT(LD,nde1,nde2,POINT_LIMIT,ne1,
     &      ne,np1,np2,np3,COFM,XP,ZD,SS,ERROR,*9999)
          
          IF(SS(1).AND.SS(2))THEN
            DO N=1,2
C............ Calculate the centre of mass of each new set of seed points              
              CALL MESH_COFM(LD,ne+N,COFM,ZD,ERROR,*9999)
              DO nj=1,3
                XP1(nj)=XP(1,1,nj,np2) !start node
              ENDDO
            ENDDO
            DO N=1,2
C............ Set up arrays for new element and node              
              CALL GN1DNE(nb,NBJ,ne,NEELEM,NENP,ne1,NKJ,NKJE,noelem,
     '          nonode,np,NP_INTERFACE,np1,NPNE,NPNODE,nr,NRE,NVJE,NVJP,
     &          NXI,SE,.TRUE.,ERROR,*9999)
              noelem_gen=noelem_gen+1
              CALL MESH_COFM(LD,ne,COFM,ZD,ERROR,*9999)
C............ Generate a branch directed towards the centre of mass              
              CALL MESH_BRANCH_CONTINUOUS(LD,LD_TEMP,numzero,ne,np1,np2,
     &          POINT_LIMIT,COFM,fraction,ANGLE_MAX,LENGTH_LIMIT,
     '          LENGTH_PARENT,min_length,XP,XP1,ZD,BRANCH,SS(N),ERROR,
     &          *9999)
              CALL ASSERT(np.LE.NPM,'>>NPM too small',ERROR,*9999)
              length=0.0d0
              DO nj=1,NJT
                XP(1,1,nj,np)=XP1(nj)
                length=length+(XP(1,1,nj,np)-XP(1,1,nj,np1))**2.0d0
              ENDDO !nj
              length=DSQRT(length)
              CALL ASSERT(NVM.GE.2,'>>Increase NVM to 2',ERROR,*9999)              
              DO nj=1,NJT !unit vector for branch direction
                XP(1,2,nj,np)=(XP(1,1,nj,np)-XP(1,1,nj,np1))/length
              ENDDO !nj
              NSTEM(ne,1)=NSTEM(ne1,1)
              IF(BRANCH.AND.SS(N))THEN
                N_ELM=N_ELM+1
                NE_TEMP(N_ELM)=ne
              ELSE !terminal
                numtb=numtb+1
              ENDIF !BRANCH
            ENDDO !N
            IF(NXI(1,0,ne2).EQ.2)THEN
              IF(ne1.EQ.NXI(1,1,ne2))THEN
                nes=NXI(1,2,ne2)
              ELSE
                nes=NXI(1,1,ne2)
              ENDIF
            ELSE IF(NXI(1,0,ne2).EQ.1)THEN
              nes=ne2
              DO WHILE(NXI(1,0,nes).NE.2)
                neb=nes
                nes=NXI(-1,1,nes)
              ENDDO
              IF(neb.EQ.NXI(1,1,nes))THEN
                nes=NXI(1,2,nes)
              ELSE
                nes=NXI(1,1,nes)
              ENDIF
            ENDIF
          ELSE
            numtb=numtb+1
            np2=NPNE(2,nb,ne1)
            MIN_DIST=1.d10
            DO nd=1,NDT
              IF(LD(nd).NE.0)THEN
                IF(LD(nd).EQ.nde1) LD(nd)=ne1
                IF(LD(nd).EQ.nde2) LD(nd)=ne1
                DIST=DSQRT((ZD(1,nd)-XP(1,1,1,np2))**2.d0+(ZD(2,nd)
     '            -XP(1,1,2,np2))**2.d0+(ZD(3,nd)-XP(1,1,3,np2))**2.d0)
                IF(DIST.LT.MIN_DIST)THEN
                  nd_min=nd
                  MIN_DIST=DIST
                ENDIF !DIST
              ENDIF !LD
            ENDDO !nd
            LD_TEMP(nd_min)=ne1
            LD(nd_min)=0
            
            LD_NUM=0
            DO nd=1,NDT
              IF(LD(nd).GT.0) LD_NUM=LD_NUM+1
            ENDDO !nd

          ENDIF
        ENDDO !M

        DO N=1,N_ELM
          NE_OLD(N)=NE_TEMP(N)
          NSTEM(NE_OLD(N),2)=0 !initialise the count of nd
        ENDDO !N

        DISTANCE_LIMIT=SCALE_DIST_LIMIT
c        !height of lung = 300.0, s.t. at ngen=30 will = 0.0
c        DISTANCE_LIMIT=MAX(300.d0-ngen*20.d0,5.d0)
        CALL MESH_REPOINT(LD,nb,N_ELM,NE_OLD,NE_REACTIVATE,
     '    NPNE,DISTANCE_LIMIT,XP,ZD,.FALSE.,ERROR,*9999)
        
        DO N=1,NE_REACTIVATE(0)
          DO M=N+1,NE_REACTIVATE(0)
            IF(NE_REACTIVATE(M).EQ.NE_REACTIVATE(N)) NE_REACTIVATE(M)=0
          ENDDO
          IF(NE_REACTIVATE(N).NE.0)THEN
            N_ELM=N_ELM+1
            NE_OLD(N_ELM)=NE_REACTIVATE(N)
          ENDIF
        ENDDO

        DO nd=1,NDT
          IF(LD(nd).NE.0)THEN
            ne_min=LD(nd)
            NSTEM(ne_min,2)=NSTEM(ne_min,2)+1
          ENDIF !LD
        ENDDO !nd

C...... Remove any terminal elements
        N_ELM_TEMP=N_ELM
        DO N=1,N_ELM
          ne_min=NE_OLD(N)
          IF(NSTEM(ne_min,2).EQ.0)THEN !find closest point
            np2=NPNE(2,nb,ne_min)
            MIN_DIST=1.d10
            DO nd=1,NDT
              IF(LD(nd).NE.0)THEN
                DIST=DSQRT((ZD(1,nd)-XP(1,1,1,np2))**2.d0+(ZD(2,nd)
     '            -XP(1,1,2,np2))**2.d0+(ZD(3,nd)-XP(1,1,3,np2))**2.d0)
                IF(DIST.LT.MIN_DIST)THEN
                  nd_min=nd
                  MIN_DIST=DIST
                ENDIF !DIST
              ENDIF !LD
            ENDDO !nd
            LD_TEMP(nd_min)=ne_min
            LD(nd_min)=0
            N_ELM_TEMP=N_ELM_TEMP-1
            NE_OLD(N)=0
            
          ELSE IF(NSTEM(ne_min,2).EQ.1)THEN
            DO nd=1,NDT
              IF(LD(nd).EQ.ne_min)THEN
                LD_TEMP(nd)=LD(nd)
                LD(nd)=0
                NE_OLD(N)=0
                N_ELM_TEMP=N_ELM_TEMP-1
              ENDIF
            ENDDO !nd
          ENDIF !NSTEM
        ENDDO !N
        DO N=1,N_ELM
          IF(NE_OLD(N).EQ.0)THEN
            I=0
            DO WHILE((N+I.LT.N_ELM).AND.(NE_OLD(N+I).EQ.0))
              I=I+1
            ENDDO
            DO M=N,N_ELM-I
              NE_OLD(M)=NE_OLD(M+I)
            ENDDO !M
            numtb=numtb+1
          ENDIF !NE_OLD
        ENDDO !N
        N_ELM=N_ELM_TEMP
        LD_NUM=0
        DO nd=1,NDT
          IF(LD(nd).GT.0) LD_NUM=LD_NUM+1
        ENDDO !nd

        WRITE(OP_STRING,'(I4,4(I7))') ngen,noelem_gen,(noelem
     '    -noelem_initial),numtb,LD_NUM
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)

      ENDDO !WHILE NT_BNS

C     Put end node associated with data point into LD
C      IF(LSP(0,1).EQ.NDT)THEN
C        NSP=0
C      ELSE
C        NSP=LSP(0,1)
C      ENDIF

      DO nd=1,NDT
        ne2=LD_TEMP(nd)
        nb=NBJ(1,ne2)
        np1=NPNE(2,nb,ne2) !end node
        LD(nd)=LD_NP(nd) !replace with temporary storage from start
        LD_NP(nd)=np1 !records closest branch node
c        LD_NP(nd)=NPNE(2,nb,ne2) !records closest branch node
      ENDDO !nd
              
      NPT(nr)=np !highest node # in nr
      NET(nr)=ne !highest element # in nr
      NPT(0)=NPT(nr) !highest node# in all regions
      NET(0)=NET(nr) !highest element# in all regions
      NEELEM(0,nr)=noelem !# of elements in nr
      NPNODE(0,nr)=nonode !# of nodes in nr
      NPNODE(0,0)=NPNODE(0,0)+noelem !number nodes in all regions
      NEELEM(0,0)=NEELEM(0,0)+nonode !number elements in all regions
      CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
      CALL CALC_NENP_1D(nb,NEELEM,NENP,NPNE,nr,ERROR,*9999)
      CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)
      
      CALL_TREE=.TRUE. !logical to specify tree generation has been performed

      CALL EXITS('GNBDMESH_CONTINUOUS')
      RETURN
 9999 CALL ERRORS('GNBDMESH_CONTINUOUS',ERROR)
      CALL EXITS('GNBDMESH_CONTINUOUS')
      RETURN 1
      END

      

