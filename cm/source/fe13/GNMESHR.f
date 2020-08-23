      SUBROUTINE GNMESHR(nb,NBJ,NEELEM,NELIST_PARENT,NENP,NKJ,NKJE,NORD,
     &  NP_INTERFACE,NPNE,NPNODE,nr,NRE,Nrefine,NSTORE,NVJE,NVJP,
     &  NXI,CE,RSTORE,SE,Spread,XP,XSTORE,MESH_TYPE,ERROR,*)
      
      
C#### Subroutine: GNMESHR
C###  Description:
C###    GNMESHR generates conducting airway, pulmonary venous,
C###    and pulmonary arterial trees for READING IN from file

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER nb,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     &  NELIST_PARENT(0:NE_R_M),NENP(NPM,0:NEPM,0:NRM),
     &  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),NORD(5,NE_R_M),
     &  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     &  NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),Nrefine,
     &  NSTORE(NE_R_M,3),NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 CE(NMM,NEM),RSTORE(NE_R_M,3),
     &  SE(NSM,NBFM,NEM),Spread,XP(NKM,NVM,NJM,NPM),XSTORE(NE_R_M,3)
      CHARACTER MESH_TYPE*30,ERROR*(*)
!     Local Variables
      INTEGER i,N,ne,ne0,ne2,ne_parent,ngen,ngen_parent,nj,noelem,
     &  noelem2,noelem_parent,nonode,np,np0,np2,nv,nx
      REAL*8 direction(3),DIST,max_z,min_z,range_z,Xi
      
      CALL ENTERS('GNMESHR',*9999)
      
      CALL ASSERT(NVM.GE.2,'>> Increase NVM to 2 ',ERROR,*9999)
      
      IF(ADD.OR.MESH_TYPE.EQ.'DEFAULT')THEN

        IF(NELIST_PARENT(0).EQ.0) THEN
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            IF(NXI(1,0,ne).EQ.0)THEN !APPEND MODEL
              NELIST_PARENT(0)=NELIST_PARENT(0)+1
              NELIST_PARENT(NELIST_PARENT(0))=ne
            ENDIF
          ENDDO !noelem
        ELSE
          DO noelem=1,NELIST_PARENT(0)
            ne=NELIST_PARENT(noelem)
c            NELIST_PARENT(0)=NELIST_PARENT(0)+1
c            NELIST_PARENT(NELIST_PARENT(0))=ne
c            DO nx=1,NXM
c              BBM(1,ne,nx)=0.d0 !reset the LPM volume to zero
c            ENDDO
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
c            DIST=L_BRANCH_GEN(ngen)
            DIST=L_BRANCH_GEN(ngen)/DBLE(Nrefine)
            DO N=1,B_BRANCH_GEN(ngen)*Nrefine !for each element in the airway
c            DO N=1,Nrefine !for each element in the airway
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
              DO nj=1,NJT
                XP(1,2,nj,np)=XP(1,2,nj,np0) !same direction as parent
                XP(1,1,nj,np)=XP(1,1,nj,np0)+DIST*XP(1,2,nj,np)
              ENDDO !nj
              np0=np
              ne0=ne
            ENDDO !N
          ENDDO !ngen
        ENDDO !noelem_cond

      ELSE IF(MESH_TYPE.EQ.'MULTI_BRANCH')THEN
        CALL GNMESH_MULTI(nb,NBJ,ne,NEELEM,NELIST_PARENT,NENP,NKJ,NKJE,
     &    noelem,nonode,NORD,np,NP_INTERFACE,NPNE,NPNODE,nr,NRE,Nrefine,
     &    NVJE,NVJP,NXI,SE,XP,MESH_TYPE,ERROR,*9999)
      ELSE IF(MESH_TYPE.EQ.'MULTI_READ')THEN

c        IF(.NOT.ADD)THEN !must create the first branch
c          CALL GN1DNE(nb,NBJ,ne,NEELEM,NENP,0,NKJ,NKJE,noelem,nonode,
c     &      np,NP_INTERFACE,1,NPNE,NPNODE,nr,NRE,NVJE,NVJP,NXI,SE,
c     &      .TRUE.,ERROR,*9999) !makes new node and element
c          NORD(1,ne)=1 !record generation #
c          DO nj=1,NJT
c            XP(1,1,nj,np)=0.d0
c            XP(1,2,nj,np)=0.d0
c          ENDDO !nj
c          XP(1,1,3,np)=-HBW_LENGTH(1)
c          XP(1,2,3,np)=-1.d0
c          NELIST_PARENT(0)=1
c          NELIST_PARENT(1)=ne
c        ENDIF

        IF(Spread.GT.0.d0)THEN
          max_z=-1.d6
          min_z=1.d6
          DO noelem2=1,NEELEM(0,nr)
            ne2=NEELEM(noelem2,nr)
            IF(NXI(1,0,ne2).EQ.0)THEN !only for the terminals
              nb=NBJ(1,ne2)
              np2=NPNE(2,nb,ne2)
              max_z=MAX(max_z,XP(1,1,3,np2))
              min_z=MIN(min_z,XP(1,1,3,np2))
            ENDIF
          ENDDO
          range_z=DABS(max_z-min_z)
          IF(DABS(range_z).LE.1d-5) range_z=1.d0
        ENDIF
        
        DO noelem_parent=1,NELIST_PARENT(0)
          ne_parent=NELIST_PARENT(noelem_parent)
          IF(ADD)THEN
            ngen_parent=NORD(1,ne_parent)
          ELSE
            ngen_parent=0
          ENDIF

          DO i=1,N_RESP
            IF(NSTORE(i,1).EQ.0)THEN
              ne0=ne_parent
            ELSE
              ne0=NET(0)+NSTORE(i,1)*Nrefine
     &          +(noelem_parent-1)*N_RESP*Nrefine
            ENDIF
            IF(ne0.NE.0)THEN
              np0=NPNE(2,nb,ne0)
              IF(Spread.GT.0.d0) Xi=(XP(1,1,3,np0)-min_z)/range_z
C...........Calculate the direction of the branch
              IF(NSTORE(i,1).EQ.0)THEN
                DO nj=1,NJT
                  direction(nj)=XSTORE(i,nj)-0.d0
                ENDDO !nj
                CALL NORMALISE(NJT,direction,ERROR,*9999)
              ELSE
                DO nj=1,NJT
                  direction(nj)=XSTORE(i,nj)-XSTORE(NSTORE(i,1),nj)
                ENDDO !nj
                CALL NORMALISE(NJT,direction,ERROR,*9999)
              ENDIF
            ELSE
              !make a starting node at 0,0,0
              np0=1
              ne0=0
              np=1
              ne=0
              nonode=1
              noelem=0
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
              direction(1)=0.d0
              direction(2)=0.d0
              direction(3)=-1.d0
            ENDIF


            DO N=1,Nrefine !for each element in the airway
              CALL GN1DNE(nb,NBJ,ne,NEELEM,NENP,ne0,NKJ,NKJE,
     &          noelem,nonode,np,NP_INTERFACE,np0,NPNE,NPNODE,nr,
     &          NRE,NVJE,NVJP,NXI,SE,.TRUE.,ERROR,*9999)
              NORD(1,ne)=NSTORE(i,2)+ngen_parent !record generation #
              IF(N.EQ.1)THEN
                NORD(5,ne)=NSTORE(i,3) !symmetry when = 1
              ELSE
                NORD(5,ne)=0
              ENDIF
              IF(Spread.GT.0.d0)THEN
                CE(1,ne)=DSQRT(1.d0+Spread)*RSTORE(i,2)*Xi+
     &            DSQRT(1.d0-Spread)*RSTORE(i,2)*(1.d0-Xi)
              ELSE
                CE(1,ne)=RSTORE(i,2) !temporary store for area
              ENDIF
              CE(2,ne)=RSTORE(i,3) !temporary store for a/A
              DO nj=1,NJT
                XP(1,2,nj,np)=direction(nj)
              ENDDO !nj
              DO nj=1,NJT
                XP(1,1,nj,np)=XP(1,1,nj,np0)+XP(1,2,nj,np)*
     &            RSTORE(i,1)/DBLE(Nrefine)
              ENDDO !nj
c              XAB(1,ne)=1.d0
              np0=np
              ne0=ne
            ENDDO !N
          ENDDO !noelem
        ENDDO !NELIST_PARENT
        NEELEM(0,nr)=noelem

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
      
      CALL EXITS('GNMESHR')
      RETURN
 9999 CALL ERRORS('GNMESHR',ERROR)
      CALL EXITS('GNMESHR')
      RETURN 1
      END



