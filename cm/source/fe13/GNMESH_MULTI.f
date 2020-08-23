      SUBROUTINE GNMESH_MULTI(nb,NBJ,ne,NEELEM,NELIST2,NENP,NKJ,NKJE,
     &  noelem,nonode,NORD,np,NP_INTERFACE,NPNE,NPNODE,nr,NRE,Nrefine,
     &  NVJE,NVJP,NXI,SE,XP,MESH_TYPE,ERROR,*)
      
      
C#### Subroutine: GNMESH_MULTI
C###  Description:
C###    GNMESH_MULTI 

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER nb,NBJ(NJM,NEM),ne,NEELEM(0:NE_R_M,0:NRM),
     &  NELIST2(0:NEM),NENP(NPM,0:NEPM,0:NRM),
     &  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),noelem,nonode,NORD(5,NE_R_M),
     &  np,NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     &  NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),Nrefine,
     &  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER MESH_TYPE*30,ERROR*(*)
!     Local Variables
      INTEGER j,M,N,ne2,
     &  NE_OLD(NE_R_M),N_ELEM,NELIST_PARENT(0:NE_R_M),ne0,ne00,
     &  ne_parent,ngen,ngen_parent,nj,np00,noelem_parent,
     &  np0,NT_BNS,INUM,NE_TEMP(NE_R_M)
      REAL*8 ANG_Y,BRANCH_L,u(3),v(3),NRML(3),
     &  NRML2(3),A(3,3),B(3),temp,NORM_RAND_NUM,ANG_Y_MN
      
      CALL ENTERS('GNMESH_MULTI',*9999)
      
      IF(.NOT.ADD)THEN !must create the first branch
        CALL GN1DNE(nb,NBJ,ne,NEELEM,NENP,0,NKJ,NKJE,noelem,nonode,
     &    np,NP_INTERFACE,1,NPNE,NPNODE,nr,NRE,NVJE,NVJP,NXI,SE,
     &    .TRUE.,ERROR,*9999) !makes new node and element
        NORD(1,ne)=1 !record generation #
        DO nj=1,NJT
          XP(1,1,nj,np)=0.d0
          XP(1,2,nj,np)=0.d0
        ENDDO !nj
        XP(1,1,3,np)=-L_BRANCH_GEN(1)*0.5d0
        XP(1,2,3,np)=-1.d0
        
        ne0=ne
        np0=np
        CALL GN1DNE(nb,NBJ,ne,NEELEM,NENP,ne0,NKJ,NKJE,noelem,nonode,
     &    np,NP_INTERFACE,np0,NPNE,NPNODE,nr,NRE,NVJE,NVJP,NXI,SE,
     &    .TRUE.,ERROR,*9999) !makes new node and element
        NORD(1,ne)=1 !record generation #
        DO nj=1,NJT
          XP(1,1,nj,np)=0.d0
          XP(1,2,nj,np)=0.d0
        ENDDO !nj
        XP(1,1,3,np)=-L_BRANCH_GEN(1)
        XP(1,2,3,np)=-1.d0

        NT_BNS=1
        NE_TEMP(NT_BNS)=ne !records elem #s at this gen
        NELIST_PARENT(0)=1
        NELIST_PARENT(1)=ne
      ENDIF

      DO ngen=1,NT_GEN 
        W_ANGLE_Y(ngen)=COND_ANGLE(ngen)
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
        
        DO ngen=2,NT_GEN !for each respiratory generation
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
     &          ZERO_TOL.AND.DABS(NRML(3)).LT.ZERO_TOL)THEN
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
     '          W_ANGLE_SD(ngen)*PI/180.d0) !generate an angle
              BRANCH_L=L_BRANCH_GEN(ngen)
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
     &            noelem,nonode,np,NP_INTERFACE,np0,NPNE,NPNODE,nr,
     &            NRE,NVJE,NVJP,NXI,SE,.TRUE.,ERROR,*9999) !makes new node and element
                NORD(1,ne)=ngen+ngen_parent !record generation #
                NE_TEMP(NT_BNS)=ne !records elem #s at this gen
                DO nj=1,NJT !set up XP for new point
                  IF(j.EQ.1)THEN
                    XP(1,2,nj,np)=B(nj) !direction stored in nv=2, very bad
                  ELSE
                    XP(1,2,nj,np)=XP(1,2,nj,np-1)
                  ENDIF
                  XP(1,1,nj,np)=XP(1,1,nj,np0)+BRANCH_L/DBLE(Nrefine)
     &              *XP(1,2,nj,np)
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
      
      CALL EXITS('GNMESH_MULTI')
      RETURN
 9999 CALL ERRORS('GNMESH_MULTI',ERROR)
      CALL EXITS('GNMESH_MULTI')
      RETURN 1
      END



