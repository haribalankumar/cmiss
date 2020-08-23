      SUBROUTINE GNMESH1(nb,NBJ,N_DELTA,NEELEM,NELIST2,NENP,NKJ,NKJE,
     &  NORD,NP_INTERFACE,NPNE,NPNODE,nr,NRE,Nrefine,NVJE,NVJP,
     &  NXI,ACINUS_VOLUME,BBM,SE,XP,MESH_TYPE,ERROR,*)
      
      
C#### Subroutine: GNMESH1
C###  Description:
C###    GNMESH1 generates conducting airway, pulmonary venous,
C###    and pulmonary arterial trees

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER nb,NBJ(NJM,NEM),N_DELTA,NEELEM(0:NE_R_M,0:NRM),
     &  NELIST2(0:NEM),NENP(NPM,0:NEPM,0:NRM),
     &  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),NORD(5,NE_R_M),
     &  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     &  NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),Nrefine,
     &  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 ACINUS_VOLUME,BBM(2,NEM),SE(NSM,NBFM,NEM),
     &  XP(NKM,NVM,NJM,NPM)
      CHARACTER MESH_TYPE*30,ERROR*(*)
!     Local Variables
      INTEGER i,DELTA(2),Gdirn,j,M,N,N_DGTHR,ne,ne2,
     &  NE_OLD(NE_R_M),N_ELEM,NELIST_PARENT(0:NE_R_M),ne0,ne00,
     &  ne_parent,ngen,ngen_parent,nj,nlpm,noelem,np00,noelem_parent,
     &  nonode,np,np0,NT_BNS,NUM_NODES,nv,INUM,NE_TEMP(NE_R_M),nx
      REAL*8 ANG_Y,BRANCH_L,DIST,u(3),v(3),NRML(3),
     &  NRML2(3),A(3,3),B(3),temp,NORM_RAND_NUM,ANG_Y_MN,max_z,
     &  min_z,range_z,Xi
      
      CALL ENTERS('GNMESH1',*9999)
      
      CALL ASSERT(NVM.GE.2,'>> Increase NVM to 2 ',ERROR,*9999)

      IF(ADD.OR.MESH_TYPE.EQ.'DEFAULT')THEN

        NELIST_PARENT(0)=0
        IF(NELIST2(0).EQ.0) THEN
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            IF(NXI(1,0,ne).EQ.0)THEN !APPEND MODEL
              NELIST_PARENT(0)=NELIST_PARENT(0)+1
              NELIST_PARENT(NELIST_PARENT(0))=ne
            ENDIF
          ENDDO !noelem
        ELSE
          DO noelem=1,NELIST2(0)
            ne=NELIST2(noelem)
            NELIST_PARENT(0)=NELIST_PARENT(0)+1
            NELIST_PARENT(NELIST_PARENT(0))=ne
          ENDDO
        ENDIF
        
        np=NPT(0) !current highest global node #
        ne=NET(0) !current highest global element #
        noelem=NEELEM(0,nr) !current # of elements in nr
        nonode=NPNODE(0,nr) !current # of nodes in nr

      ELSE
C     Overwrite the current elements

        SYMMETRIC_NE=.FALSE.
        DO N=1,NE_R_M
          DO i=1,5
            NORD(i,N)=0
          ENDDO
        ENDDO
        
        NELIST_PARENT(0)=1
        NELIST_PARENT(1)=0
        
        IF(NEELEM(0,nr).GT.0)THEN
          WRITE(OP_STRING,'('' Overwriting existing elements''
     &      /''Use FEM DEFINE;ADD MESH to append to existing mesh'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

        ne=0
        noelem=0
        np=1
        nonode=1
        NPNODE(nonode,nr)=np
        NPNE(1,nb,1)=0 !to start the mesh
        NPNE(2,nb,1)=1 !to start the mesh
        NP_INTERFACE(np,0)=1 !set up for initial node
        NP_INTERFACE(np,1)=nr
        DO nj=1,NJT
          XP(1,1,nj,np)=0.d0
          XP(1,2,nj,np)=0.d0
        ENDDO !nj
        XP(1,2,3,np)=-1.d0
      ENDIF

      IF(MESH_TYPE.EQ.'SYMMETRIC')THEN
        SYMMETRIC_NE=.TRUE.
C  SYMMETRIC AIRWAYS
C  A single branch represents all branches from its generation.
C  The default values for length are 1.d0.
C  NORD(5,ne)=1 tags element ne as at the start of the next
C  generation.  This is used in the solution to effectively double
C  the cross-sectional area at this point. 

        DO noelem_parent=1,NELIST_PARENT(0)
          ne0=NELIST_PARENT(noelem_parent)
          IF(ne0.NE.0)THEN
            np0=NPNE(2,nb,ne0)
          ELSE
            np0=1 !ne0=0 for overwriting mesh
          ENDIF
          DO ngen=1,NTT_GEN !for each symmetric generation
            DIST=L_BRANCH_GEN(ngen)/DBLE(Nrefine)
c            DO N=1,B_BRANCH_GEN(ngen) !for each element in the airway
            DO N=1,Nrefine !for each element in the airway
              CALL GN1DNE(nb,NBJ,ne,NEELEM,NENP,ne0,NKJ,NKJE,
     '          noelem,nonode,np,NP_INTERFACE,np0,NPNE,NPNODE,nr,NRE,
     &          NVJE,NVJP,NXI,SE,.TRUE.,ERROR,*9999)
              NORD(1,ne)=ngen !record branch generation
              IF(N.EQ.1)THEN
                IF((ngen.NE.1).OR.(ngen.EQ.1.AND.ADD))THEN
                  NORD(5,ne)=1 !start of 'half' branch
                ENDIF
              ELSE
                NORD(5,ne)=0
              ENDIF
              DO nj=1,2
                DO nv=1,2
                  XP(1,nv,nj,np)=0.d0
                ENDDO !nv
              ENDDO !nj
C             XP(1,1,3,np)=XP(1,1,3,np-1)-DIST !position of new node
C!! KSB - think the line above should really be this:
              XP(1,1,3,np)=XP(1,1,3,np0)-DIST !position of new node
              XP(1,2,3,np)=-1.d0
              np0=np
              ne0=ne
            ENDDO !N
          ENDDO !ngen
        ENDDO !noelem_cond

      ELSE IF(MESH_TYPE.EQ.'HORSFIELD')THEN

        IF(ADD)THEN
          ngen=NTT_GEN
          NORD(2,ne)=NORD(2,1) !NORD(2,1) will be overwritten later
          NT_BNS=1 !create the rest of the airways, by generation
          NELIST2(1)=ne !start from last symmetric airway element
          DO WHILE(NT_BNS.NE.0)
            NUM_NODES=NT_BNS
            NT_BNS=0
            ngen=ngen+1
            DO M=1,NUM_NODES
              ne0=NELIST2(M) !parent global element #
              np0=NPNE(2,nb,ne0) !parent global end node #
              DELTA(1)=NORD(2,ne0)-1 !highest daughter order
              IF(NORD(2,ne0).EQ.2)THEN
                DELTA(2)=1 !order 1, will be terminal branch
                N_DGTHR=1
                NORD(5,ne+1)=1
              ELSE
                IF(NORD(2,ne0).EQ.3.OR.NORD(2,ne0).EQ.4.OR.
     '            NORD(2,ne0).EQ.5)THEN
                  DELTA(2)=1
                ELSEIF(NORD(2,ne0).EQ.6)THEN
                  DELTA(2)=3
                ELSE
                  DELTA(2)=NORD(2,ne0)-(N_DELTA+1)
                ENDIF !NORD(2,ne0)
                N_DGTHR=2
              ENDIF

C             ELSEIF(NORD(2,ne0).EQ.3)THEN
C             DELTA(2)=2
C             ELSEIF(NORD(2,ne0).LE.5)THEN
C             DELTA(2)=3
C             ELSEIF(NORD(2,ne0).LE.7)THEN
C             DELTA(2)=NORD(2,ne0)-3
C             ELSEIF(NORD(2,ne0).GE.8)THEN !high enough to use Delta
C             DELTA(2)=NORD(2,ne0)-(N_DELTA+1)
C             ENDIF !NORD(2,ne0)
              DO N=1,N_DGTHR !for each new daughter branch
                CALL GN1DNE(nb,NBJ,ne,NEELEM,NENP,ne0,NKJ,NKJE,
     '            noelem,nonode,np,NP_INTERFACE,np0,NPNE,NPNODE,nr,NRE,
     &            NVJE,NVJP,NXI,SE,.TRUE.,ERROR,*9999)
                NORD(1,ne)=ngen
                NORD(2,ne)=DELTA(N)
                IF(DELTA(N).GT.1)THEN
                  NT_BNS=NT_BNS+1
                  NELIST2(NUM_NODES+NT_BNS)=ne
                ENDIF
                XP(1,1,3,np)=XP(1,1,3,np0)
     &            -HORSFIELD_AIRWAY_LENGTH(DELTA(N))
                XP(1,2,3,np)=-1.d0 !vector
              ENDDO !N
            ENDDO !M
            DO N=1,NT_BNS
              NELIST2(N)=NELIST2(NUM_NODES+N) !updates list of prev gen elem#s
            ENDDO !N
          ENDDO !WHILE
          NTT_GEN=NORD(1,ne)
          
        ELSEIF(.NOT.ADD)THEN
          
C***  HORSFIELD DELTA MODEL
C***  See Horsfield et al.(1971) JAP 31:207-217. Delta is the
C***  difference in order of daughter branches. i.e. if parent is order
C***  N, then daughter of highest order will be N-1 and other daughter
C***  will be N-1-delta.  For the lower orders, the scheme departs
C***  from Delta to ensure that all terminal branches are order 1.

          CALL ASSERT(NKM.GE.2,'>>Increase NKM to 2',ERROR,*9999)
          DO nj=1,NJT !position top of trachea at [0,0,0]
            XP(1,1,nj,np)=0.d0
            XP(2,1,nj,np)=0.d0
          ENDDO !nj
          NPNODE(nonode,nr)=np
          !Trachea
          ngen=1
          DIST=HORSFIELD_AIRWAY_LENGTH(28)
          np0=np
          ne0=ne
          CALL GN1DNE(nb,NBJ,ne,NEELEM,NENP,ne0,NKJ,NKJE,noelem,
     '      nonode,np,NP_INTERFACE,np0,NPNE,NPNODE,nr,NRE,NVJE,NVJP,
     &      NXI,SE,.TRUE.,ERROR,*9999)
          NORD(1,ne)=ngen !record branch generation
          NORD(2,ne)=28 !set the trachea order
          XP(1,1,3,np)=XP(1,1,3,np-1)-DIST !position of new node
          XP(1,2,3,np)=-1.d0

          NT_BNS=1 !create the rest of the airways, by generation
          NELIST2(1)=ne
          DO WHILE(NT_BNS.NE.0)
            NUM_NODES=NT_BNS
            NT_BNS=0
            ngen=ngen+1
            DO M=1,NUM_NODES
              ne0=NELIST2(M) !parent global element #
              np0=NPNE(2,nb,ne0) !parent global end node #
              DELTA(1)=NORD(2,ne0)-1 !highest daughter order
              IF(NORD(2,ne0).EQ.2)THEN !other daughter order is determined
                DELTA(2)=1 !order 1, will be terminal branch
              ELSEIF(NORD(2,ne0).EQ.3)THEN
                DELTA(2)=2
              ELSEIF(NORD(2,ne0).LE.5)THEN
                DELTA(2)=3
              ELSEIF(NORD(2,ne0).LE.7)THEN
                DELTA(2)=NORD(2,ne0)-3
              ELSEIF(NORD(2,ne0).GE.8)THEN !high enough to use Delta
                DELTA(2)=NORD(2,ne0)-(N_DELTA+1)
              ENDIF !NORD(2,ne0)
              DO N=1,2 !for each new daughter branch
                CALL GN1DNE(nb,NBJ,ne,NEELEM,NENP,ne0,NKJ,NKJE,
     '            noelem,nonode,np,NP_INTERFACE,np0,NPNE,NPNODE,nr,NRE,
     &            NVJE,NVJP,NXI,SE,.TRUE.,ERROR,*9999)
                NORD(1,ne)=ngen
                NORD(2,ne)=DELTA(N)
                IF(DELTA(N).GT.1)THEN
                  NT_BNS=NT_BNS+1
                  NELIST2(NUM_NODES+NT_BNS)=ne
                ENDIF
                XP(1,1,3,np)=XP(1,1,3,np0)
     &            -HORSFIELD_AIRWAY_LENGTH(DELTA(N))
                XP(1,2,3,np)=-1.d0 !vector
              ENDDO !N
            ENDDO !M
            DO N=1,NT_BNS
              NELIST2(N)=NELIST2(NUM_NODES+N) !updates list of prev gen elem#s
            ENDDO !N
          ENDDO !WHILE

        ENDIF!ADD  

      ELSE IF(MESH_TYPE.EQ.'MULTI_HBW')THEN

        IF(.NOT.ADD)THEN !must create the first branch
          ne0=ne
          np0=np
          CALL GN1DNE(nb,NBJ,ne,NEELEM,NENP,ne0,NKJ,NKJE,noelem,nonode,
     &      np,NP_INTERFACE,np0,NPNE,NPNODE,nr,NRE,NVJE,NVJP,NXI,SE,
     &      .TRUE.,ERROR,*9999) !makes new node and element
          NORD(1,ne)=1 !record generation #
          NT_BNS=1
          NE_TEMP(NT_BNS)=ne !records elem #s at this gen
          DO nj=1,NJT
            XP(1,1,nj,np)=0.d0
            XP(1,2,nj,np)=0.d0
          ENDDO !nj
          XP(1,1,3,np)=-HBW_LENGTH(1)
          XP(1,2,3,np)=-1.d0
          NELIST_PARENT(0)=1
          NELIST_PARENT(1)=ne
        ENDIF

        NT_GEN=11 !11 gens for HBW
        DO ngen=1,NT_GEN 
          W_ANGLE_Y(ngen)=COND_ANGLE(ngen)
        ENDDO !ngen
        DO ngen=1,NT_GEN !for each respiratory generation
c          W_ANGLE_Y(ngen)=ACINUS_ANGLE(ngen+1)
          MNB(ngen)=HBW_BRANCH_MEAN(ngen+1)
          SDNB(ngen)=HBW_BRANCH_SD(ngen+1)
          B_LENGTH(ngen)=HBW_LENGTH(ngen+1)
        ENDDO !ngen
        W_ANGLE_SD(1)=0.0d0
        B_ANGLE_XY(1)=0.d0
        DO ngen=2,NT_GEN !Calc std devs and initial angles
          W_ANGLE_SD(ngen)=0.1d0*W_ANGLE_Y(ngen)
          B_ANGLE_XY(ngen)=B_ANGLE_XY(2)
        ENDDO !ngen

        DO noelem_parent=1,NELIST_PARENT(0)
          ne_parent=NELIST_PARENT(noelem_parent)
          IF(ADD)THEN
            ngen_parent=NORD(1,ne_parent)
          ELSE
            ngen_parent=0
          ENDIF
          NE_OLD(1)=ne_parent !stores previous generation node #s
          NT_BNS=1 !Total number of branches in generation

          DO ngen=1,NT_GEN !for each respiratory generation
            N_ELEM=NT_BNS !# of elements in previous generation
            NT_BNS=0 !initialise # elements in new generation
            DO M=1,N_ELEM !for each parent element
              ne0=NE_OLD(M) !parent elem #
              IF(ne0.NE.0)THEN
                np0=NPNE(2,nb,ne0) !start node for new branch
                np00=NPNE(1,nb,ne0) !parent start node
              ELSE
                np0=np !initialised to 1
                np00=0
              ENDIF
              IF(np00.NE.0)THEN
                IF(Nrefine.GT.1)THEN !use from previous generation
                  ne2=NXI(-1,1,ne0)
                  DO WHILE(NORD(1,ne0).NE.NORD(1,ne2))
                    ne2=NXI(-1,1,ne2)
                  ENDDO
                  DO nj=1,NJT !get directions
                    u(nj)=XP(1,2,nj,NPNE(2,nb,ne2)) !parent unit vector
                    v(nj)=XP(1,2,nj,NPNE(1,nb,ne2)) !grandparent unit vector
                  ENDDO !nj
                ELSE
                  DO nj=1,NJT !get directions
                    u(nj)=XP(1,2,nj,np0) !parent unit vector
                    v(nj)=XP(1,2,nj,np00) !grandparent unit vector
                  ENDDO !nj
                ENDIF
                CALL CROSS(u,v,NRML) !get normal to parent-grandparent
                IF(DABS(NRML(1)).LT.ZERO_TOL.AND.DABS(NRML(2)).LT.
     &            ZERO_TOL.AND.DABS(NRML(3)).LT.ZERO_TOL)THEN
                  NRML(1)=0.d0
                  NRML(2)=-1.d0/DSQRT(2.d0)
                  NRML(3)=-1.d0/DSQRT(2.d0)
                ELSE
                  CALL NORMALISE(NJT,NRML,ERROR,*9999) !unit vector
                ENDIF
              ELSE
                NRML(1)=0.d0
                NRML(2)=-1.d0/DSQRT(2.d0)
                NRML(3)=-1.d0/DSQRT(2.d0)
              ENDIF
              temp=NORM_RAND_NUM(ISEED_BRANCH,MNB(ngen),SDNB(ngen))
              INUM=NINT(temp) !integer # of daughter branches
              IF(INUM.LE.1.AND.ngen.GT.1) INUM=0 !# of daughter branches
              IF(INUM.GT.2) INUM=2 !always 0 or 2
              DO N=1,INUM !for each daughter branch
                NT_BNS=NT_BNS+1 !count # of elems in new generation

                ANG_Y_MN=W_ANGLE_Y(ngen)*PI/180.d0 !mean angle
                ANG_Y=NORM_RAND_NUM(ISEED_GEOM,ANG_Y_MN,
     '            W_ANGLE_SD(ngen)*PI/180.d0) !generate an angle
                BRANCH_L=B_LENGTH(ngen)
                CALL CROSS(u,NRML,NRML2) !calculate branching plane
                CALL NORMALISE(NJT,NRML2,ERROR,*9999) !unit vector
                DO nj=1,NJT !set up system of linear equations
                  A(1,nj)=u(nj) !dotprod parent and new element
                  A(2,nj)=NRML(nj) !dotprod normal and new element
                  A(3,nj)=NRML2(nj) !dotprod plane and new element
                ENDDO !nj
                B(1)=DCOS(ANG_Y) !angle btwn parent & element=ang_y
                IF(N.EQ.1)THEN !for first daughter
                  B(2)=DCOS(PI/2.d0-ANG_Y) !angle btwn normal & element
                ELSE !for second daughter
                  B(2)=DCOS(PI/2.d0+ANG_Y) !angle btwn normal & element
                ENDIF !N.EQ.1
                B(3)=0.d0 !on plane:(w-p).nrml=const;nrml.p=const
                CALL FSOLV3V(A,B,ERROR,*9999) !solve ax=b
                CALL NORMALISE(NJT,B,ERROR,*9999) !unit

                ne00=ne0 !temporary store
                np00=np0 !temporary store
                DO j=1,Nrefine
                  CALL GN1DNE(nb,NBJ,ne,NEELEM,NENP,ne0,NKJ,NKJE,
     &              noelem,nonode,np,NP_INTERFACE,np0,NPNE,NPNODE,nr,
     &              NRE,NVJE,NVJP,NXI,SE,.TRUE.,ERROR,*9999) !makes new node and element
                  NORD(1,ne)=ngen+ngen_parent !record generation #
                  NE_TEMP(NT_BNS)=ne !records elem #s at this gen
                  DO nj=1,NJT !set up XP for new point
                    IF(j.EQ.1)THEN
                      XP(1,2,nj,np)=B(nj) !direction stored in nv=2, very bad
                    ELSE
                      XP(1,2,nj,np)=XP(1,2,nj,np-1)
                    ENDIF
                    XP(1,1,nj,np)=XP(1,1,nj,np0)+BRANCH_L/DBLE(Nrefine)
     &                *XP(1,2,nj,np)
                  ENDDO !nj
                  ne0=ne
                  np0=np
                ENDDO !j
                ne0=ne00
                np0=np00
              ENDDO !N
            ENDDO !M
            DO N=1,NT_BNS !for each element in new generation
              NE_OLD(N)=NE_TEMP(N) !updates list of prev gen elem#s
            ENDDO !N
          ENDDO !ngen
        ENDDO !NELIST_PARENT
        NEELEM(0,nr)=noelem

      ELSE IF(MESH_TYPE.EQ.'LUMPED_PARAMETER')THEN
        CALL ASSERT(ADD,'>> Use DEFINE;ADD MESH;C AIRWAY',ERROR,*9999)
        nlpm=0
        DO noelem_parent=1,NELIST_PARENT(0)
          ne_parent=NELIST_PARENT(noelem_parent)
c          nlpm=nlpm+1
c          XAB(1,nlpm)=ne_parent
c          XAB(2,nlpm)=ACINUS_VOLUME
c          XAB(3,nlpm)=ACINUS_VOLUME
          DO nx=1,NXM
            BBM(1,ne_parent)=ACINUS_VOLUME !initial volume
            BBM(2,ne_parent)=0.d0 !initial concentration
          ENDDO !nx
        ENDDO !noelem_parent
        NTB=NELIST_PARENT(0)
c        IF(Spread.GT.0.d0)THEN
c          Gdirn=3 !for z direction
c          max_z=-1.d6
c          min_z=1.d6
c          DO noelem=1,NELIST_PARENT(0)
c            ne=NELIST_PARENT(noelem)
c            nb=NBJ(1,ne)
c            np=NPNE(2,nb,ne)
c            max_z=MAX(max_z,XP(1,1,Gdirn,np))
c            min_z=MIN(min_z,XP(1,1,Gdirn,np))
c          ENDDO
c          range_z=DABS(max_z-min_z)
c          IF(DABS(range_z).LE.1d-5) range_z=1.d0
c          DO noelem=1,NELIST_PARENT(0)
c            ne=NELIST_PARENT(noelem)
c            nb=NBJ(1,ne)
c            np=NPNE(2,nb,ne)
c            Xi=(XP(1,1,Gdirn,np)-min_z)/range_z
c            BBM(1,ne_parent)=(ACINUS_VOLUME-0.5d0*Spread)*Xi
c     &        +(ACINUS_VOLUME+0.5d0*Spread)*(1.d0-Xi)
c          ENDDO
c        ENDIF !spread
      ENDIF !mesh_types

      NPT(nr)=np
      NET(nr)=ne
      NPT(0)=NPT(nr)
      NET(0)=NET(nr)
      NEELEM(0,nr)=noelem
      NPNODE(0,nr)=nonode
      
      nb=NBJ(1,NEELEM(1,nr))
      
      CALL CALC_NENP_1D(nb,NEELEM,NENP,NPNE,nr,ERROR,*9999)
      CALL NENXI_1D(nb,NEELEM,NENP,NPNE,nr,NXI,ERROR,*9999)
      
      CALL EXITS('GNMESH1')
      RETURN
 9999 CALL ERRORS('GNMESH1',ERROR)
      CALL EXITS('GNMESH1')
      RETURN 1
      END



