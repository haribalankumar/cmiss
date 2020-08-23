      SUBROUTINE GENPCAP(INLET,nb,nb_host,NBJ,NEELEM,NELIST,NELIST2,
     '  NENP,NKJ,NKJE,NP_INTERFACE,NPLIST,NPNE,NPNODE,NRE,nr,nr_feed,
     '  nr_host,NVJE,NVJP,NXI,OUTLET,C,SE,XP,REGULAR,ERROR,*)

C#### Subroutine: GENPCAP
C###  Description:
C###    GENPCAP creates a Voronoi capillary mesh over the surface of
C###    a unit sphere. Given the elements making up each alveolar
C###    group, and the centre of each alveolus, this voronoi mesh is
C###    then projected onto the surface of an alveolar mesh.

C***   Regular points are generated on the surface of a unit sphere.
C***   These points are then projected onto alveolar surface & the
C***   capillary mesh nodes are grouped according to alveolar elements.
C***   This ensures only 1 capillary sheet is created between adjacent
C***   alveoli. The capillary nodes are then projected back onto surface
C***   of unit half-sphere & delaunay triangulation & voronoi diagram
C***   is created. This mesh is then projected onto alveolar surface.
C***   Adjacent alveolar-capillary meshes are made continuous.

C###    Created by KSB, September 2001.


      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn' !MEM_INIT
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'mesh00.cmn'
      INCLUDE 'pulm00.cmn'
      INCLUDE 'mach00.inc' !DP_TYPE

!     Parameter List
      INTEGER INLET,nb,nb_host,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NELIST2(0:NEM),NENP(NPM,0:NEPM,0:NRM),
     '  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     '  NP_INTERFACE(0:NPM,0:3),NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,nr_feed,nr_host,NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 C(MAX_ALVEOLI,NJT),SE(NSM,NBFM,NEM),
     '  XP(NKM,NVM,NJM,NPM)
      LOGICAL REGULAR
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ne,ne2,ne_cap,nj,nn,noelem,noelem2,
     '  noelem3,nonode,nonode2,nonode3,np,np1,np2,np_cap,OUTLET
      REAL*8 radius,tol
      LOGICAL ELEM_PRESENT
      INTEGER*4 BOUNDARY_PTR,BOUNDARY2_PTR,BOUNDARY_TEMP_PTR,
     &  DNP_VNP_PTR,NE_REPEAT_PTR,NPI_PTR,NPLIST2_PTR,
     &  NPLIST3_PTR,NPNE_ALV_PTR,VERTICES_XYZ_PTR

      CALL ENTERS('GENPCAP',*9999)

      np_cap=NPT(0)+2 !Adding nodes & elements after max # present
      ne_cap=NET(0)+2 !Inlet node & element 1st # & in 1st position
      noelem2=2 !loop variable for positioning element in NEELEM array
      nonode=2
      IF(CAP_MESH.EQ.2) THEN !simple test mesh on regular mesh
        noelem2=1
        ne_cap=NET(0)+1
        DO nonode=1,NPNODE(0,nr_host) !copying regular mesh nodes into
          np=NPNODE(nonode,nr_host) !capillary region NPNODE
          NPNODE(nonode+1,nr)=np
        ENDDO !nonode
        DO noelem=1,NEELEM(0,nr_host)
          ne=NEELEM(noelem,nr_host)
          DO nn=2,3 !loops over all nodes for 2D element
            ELEM_PRESENT=.FALSE.
            np1=NPNE(1,nb_host,ne) !element start node #
            np2=NPNE(nn,nb_host,ne)
C... loop ensures elements not duplicated
            DO noelem3=2,NEELEM(0,nr)
              ne2=NEELEM(noelem3,nr)
              IF((NPNE(1,nb,ne2).EQ.np1.AND.NPNE(2,nb,ne2).EQ.np2).OR.
     '          (NPNE(2,nb,ne2).EQ.np1.AND.NPNE(1,nb,ne2).EQ.np2))THEN
                ELEM_PRESENT=.TRUE.
              ENDIF
            ENDDO !noelem3
            IF(.NOT.ELEM_PRESENT) THEN
              ne_cap=ne_cap+1
              noelem2=noelem2+1
              NEELEM(noelem2,nr)=ne_cap
              NPNE(1,nb,ne_cap)=np1
              NPNE(2,nb,ne_cap)=np2
              CALL GN1DNEJ(nb,NBJ,ne_cap,NKJ,NKJE,NPNE(2,nb,ne_cap),
     '          NPNE(1,nb,ne_cap),nr,NRE,NVJE,NVJP,SE,ERROR,*9999)
              NEELEM(0,nr)=noelem2
            ENDIF !.NOT.ELEM_PRESENT
          ENDDO !nn
          DO nn=2,3
            ELEM_PRESENT=.FALSE.
            np1=NPNE(nn,nb_host,ne)
            np2=NPNE(4,nb_host,ne)
C... loop ensures elements not duplicated
            DO noelem3=2,NEELEM(0,nr)
              ne2=NEELEM(noelem3,nr)
              IF((NPNE(1,nb,ne2).EQ.np1.AND.NPNE(2,nb,ne2).EQ.np2).OR.
     '          (NPNE(2,nb,ne2).EQ.np1.AND.NPNE(1,nb,ne2).EQ.np2))THEN
                ELEM_PRESENT=.TRUE.
              ENDIF
            ENDDO !noelem3
            IF(.NOT.ELEM_PRESENT) THEN
              ne_cap=ne_cap+1
              noelem2=noelem2+1
              NEELEM(noelem2,nr)=ne_cap
              NPNE(1,nb,ne_cap)=np1
              NPNE(2,nb,ne_cap)=np2
              CALL GN1DNEJ(nb,NBJ,ne_cap,NKJ,NKJE,NPNE(2,nb,ne_cap),
     '          NPNE(1,nb,ne_cap),nr,NRE,NVJE,NVJP,SE,ERROR,*9999)
              NEELEM(0,nr)=noelem2
            ENDIF !.NOT.ELEM_PRESENT
          ENDDO !nn
        ENDDO !noelem
      ELSEIF(CAP_MESH.EQ.3) THEN !Voronoi mesh used to create network
        radius=1.d0 !unit sphere, simplifies geometry
        BOUNDARY_PTR=0
        BOUNDARY2_PTR=0
        BOUNDARY_TEMP_PTR=0
        DNP_VNP_PTR=0
        NE_REPEAT_PTR=0
        NPI_PTR=0
        NPLIST2_PTR=0
        NPLIST3_PTR=0
        NPNE_ALV_PTR=0
        VERTICES_XYZ_PTR=0
        CALL ALLOCATE_MEMORY((N_BOUNDARY+1)*N_ALVEOLI*2*2,1,INTTYPE,
     &    BOUNDARY_PTR,MEM_INIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(N_BOUNDARY*2*2,1,INTTYPE,BOUNDARY2_PTR,
     &    MEM_INIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(N_BOUNDARY*2*2,1,INTTYPE,BOUNDARY_TEMP_PTR,
     &    MEM_INIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(N_ALVEOLI*MAX_TRIANGLES*4,1,INTTYPE,
     &    DNP_VNP_PTR,MEM_INIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(NEREPM+1,1,INTTYPE,NE_REPEAT_PTR,MEM_INIT,
     &    ERROR,*9999)
        CALL ALLOCATE_MEMORY(N_VERT_EST,1,INTTYPE,NPI_PTR,MEM_INIT,
     &    ERROR,*9999)
        CALL ALLOCATE_MEMORY(NPM+1,1,INTTYPE,NPLIST2_PTR,MEM_INIT,ERROR,
     &    *9999)
        CALL ALLOCATE_MEMORY(NPM+1,1,INTTYPE,NPLIST3_PTR,MEM_INIT,ERROR,
     &    *9999)
        CALL ALLOCATE_MEMORY((NP_NE+1)*NE_R_M,1,INTTYPE,NPNE_ALV_PTR,
     &    MEM_INIT,ERROR,*9999)
        CALL ALLOCATE_MEMORY(3*N_VERT_EST,1,DPTYPE,VERTICES_XYZ_PTR,
     '    MEM_INIT,ERROR,*9999)
        
        CALL GENPCAP_SUB(%VAL(BOUNDARY_PTR),%VAL(BOUNDARY2_PTR),
     &    %VAL(BOUNDARY_TEMP_PTR),%VAL(DNP_VNP_PTR),nb,NBJ,NEELEM,
     &    NELIST,NELIST2,%VAL(NE_REPEAT_PTR),NKJ,NKJE,%VAL(NPI_PTR),
     &    NP_INTERFACE,NPLIST,%VAL(NPLIST2_PTR),%VAL(NPLIST3_PTR),NPNE,
     &    %VAL(NPNE_ALV_PTR),NPNODE,nr,NRE,NVJE,NVJP,C,radius,SE,
     &    %VAL(VERTICES_XYZ_PTR),XP,REGULAR,ERROR,*9999)

        CALL FREE_MEMORY(BOUNDARY_PTR,ERROR,*9999)
        CALL FREE_MEMORY(BOUNDARY2_PTR,ERROR,*9999)
        CALL FREE_MEMORY(BOUNDARY_TEMP_PTR,ERROR,*9999)
        CALL FREE_MEMORY(DNP_VNP_PTR,ERROR,*9999)
        CALL FREE_MEMORY(NE_REPEAT_PTR,ERROR,*9999)
        CALL FREE_MEMORY(NPI_PTR,ERROR,*9999)
        CALL FREE_MEMORY(NPLIST2_PTR,ERROR,*9999)
        CALL FREE_MEMORY(NPLIST3_PTR,ERROR,*9999)
        CALL FREE_MEMORY(NPNE_ALV_PTR,ERROR,*9999)
        CALL FREE_MEMORY(VERTICES_XYZ_PTR,ERROR,*9999)
        
        IF(DOP) THEN
          WRITE(OP_STRING,'(''Done GENPCAP_SUB'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF !creates & projects the voronoi capillary mesh
      ENDIF !CAP_MESH
      IF(CAP_MESH.EQ.2) THEN !create arteriole & venule
        IF(OUTLET.EQ.3) THEN !outlet venule(else join to gencirc mesh)
          DO nj=1,NJT
            XP(1,1,nj,np_cap)=1.20d0
          ENDDO !nj
          nonode=nonode+1
          np1=np_cap-2
          NPNODE(nonode,nr)=np_cap
          np2=np_cap
          CAP_OUTLET=ne_cap+1
          NEELEM(noelem2+1,nr)=CAP_OUTLET
          NPNE(1,nb,CAP_OUTLET)=np1 !Creates a venule outlet element
          NPNE(2,nb,CAP_OUTLET)=np2
          DO nj=1,NJT
            NBJ(nj,CAP_OUTLET)=nb !Basis function in each direction
          ENDDO !nj
          NRE(CAP_OUTLET)=nr !new elem in nr
        ENDIF !OUTLET
        IF(INLET.EQ.3) THEN !inlet arteriole(else join to gencirc mesh)
          DO nj=1,NJT
            XP(1,1,nj,(NPT(0)+1))=-0.2d0 !1st node in region
          ENDDO !nj
          NPNODE(1,nr)=NPT(0)+1 !1st node
          np1=NPT(0)+1
          np2=NPNODE(1,nr_host)
          CAP_INLET=NET(0)+1
          NEELEM(1,nr)=CAP_INLET !This = 1st element in region
          NPNE(1,nb,CAP_INLET)=np1 !Creates arteriole inlet element
          NPNE(2,nb,CAP_INLET)=np2
          DO nj=1,NJT
            NBJ(nj,CAP_INLET)=nb !Basis function in each direction
          ENDDO  !nj
          NRE(CAP_INLET)=nr !new elem in nr
          IF(CAP_MESH.EQ.3) THEN !create inlet & outlet ne
            ne_cap=CAP_INLET-1 !temp fix so correct ne # used below
            noelem2=noelem2-1
          ENDIF !CAP_MESH.EQ.3
        ENDIF !INLET
      ENDIF !CAP_MESH.EQ.2.OR.CAP_MESH.EQ.3
      IF(CAP_MESH.EQ.2) THEN
        NPNODE(0,nr)=nonode !# nodes in region nr
        NPNODE(0,0)=np_cap !Highest node number in all regions
        NPT(nr)=np_cap !Highest node # in nr
        NPT(0)=np_cap !Highest node in all regions
        NET(nr)=ne_cap+1 !highest element # in region nr
        NET(0)=ne_cap+1 !highest element # in any region
        NEELEM(0,nr)=noelem2+1
        NEELEM(0,0)=ne_cap+1
      ENDIF !CAP_MESH
      IF(CAP_MESH.EQ.3) THEN !voronoi
C... also check whether there is already a node at these co-ordinates
C... mirror image delauanay triangles, even though not shared triangles
C... will create voronoi point with the same co-ords
        tol=1.d-6
C... called here because NENP required below. This is called again
C... after the final mesh is made
        CALL CALC_NENP_1D(nb,NEELEM,NENP,NPNE,nr,ERROR,*9999)
        nonode=0
        DO WHILE(nonode.LT.NPNODE(0,nr))
          nonode=nonode+1
          np=NPNODE(nonode,nr) !all nodes in voronoi mesh
          DO nonode2=1,NPNODE(0,nr)
            np2=NPNODE(nonode2,nr)
            IF(np.NE.np2) THEN
              IF(XP(1,1,1,np).GT.(XP(1,1,1,np2)-tol).AND.
     '          XP(1,1,1,np).LT.(XP(1,1,1,np2)+tol).AND.
     '          XP(1,1,2,np).GT.(XP(1,1,2,np2)-tol).AND.
     '          XP(1,1,2,np).LT.(XP(1,1,2,np2)+tol).AND.
     '          XP(1,1,3,np).GT.(XP(1,1,3,np2)-tol).AND.
     '          XP(1,1,3,np).LT.(XP(1,1,3,np2)+tol)) THEN
                NPNODE(0,nr)=NPNODE(0,nr)-1 !remove np2
                DO nonode3=nonode2,NPNODE(0,nr)
                  NPNODE(nonode3,nr)=NPNODE(nonode3+1,nr)
                ENDDO !nonode3
C... replace np2 with np in NPNE array
                DO noelem=1,NENP(np2,0,nr)
                  ne=NENP(np2,noelem,nr) !any ne np2 is in
                  DO nn=1,NNT(nb)
                    IF(NPNE(nn,nb,ne).EQ.np2) NPNE(nn,nb,ne)=np !replace
                  ENDDO !nn
                  IF(NPNE(1,nb,ne).EQ.NPNE(2,nb,ne)) THEN !remove elem
                    noelem2=0
                    ne2=0 !initialise
                    DO WHILE(ne.NE.ne2.AND.noelem2.LT.NEELEM(0,nr))
                      noelem2=noelem2+1 !finds position of ne in NEELEM
                      ne2=NEELEM(noelem2,nr)
                    ENDDO
                    DO noelem3=noelem2,NEELEM(0,nr) !shifts ne's up 1
                      NEELEM(noelem3,nr)=NEELEM(noelem3+1,nr)
                    ENDDO !noelem3
                    NEELEM(0,nr)=NEELEM(0,nr)-1
C... check if this element is now duplicated i.e same nodes
                  ELSE
                    DO noelem2=1,NEELEM(0,nr) !loop over every element
                      ne2=NEELEM(noelem2,nr)
                      IF(ne.NE.ne2) THEN
                        IF(NPNE(1,nb,ne).EQ.NPNE(1,nb,ne2).AND.
     '                    NPNE(2,nb,ne).EQ.NPNE(2,nb,ne2).OR.
     '                    NPNE(1,nb,ne).EQ.NPNE(2,nb,ne2).AND.
     '                    NPNE(2,nb,ne).EQ.NPNE(1,nb,ne2)) THEN !remove
                          NEELEM(0,nr)=NEELEM(0,nr)-1
                          DO noelem3=noelem2,NEELEM(0,nr) !moves ne's up
                            NEELEM(noelem3,nr)=NEELEM(noelem3+1,nr)
                          ENDDO !noelem3
                        ENDIF
                      ENDIF
                    ENDDO !noelem3
                  ENDIF
                ENDDO !noelem
              ENDIF
            ENDIF
          ENDDO !nonode2
        ENDDO !WHILE nonode
        radius=0.5d0 !for test problems
      ENDIF
      CALL CAP_IO(nb,NBJ,NEELEM,NKJ,NKJE,NPNE,NPNODE,nr,nr_feed,
     &  NRE,NVJE,NVJP,SE,XP,ERROR,*9999)
      WRITE(OP_STRING,'('' Total # elements in mesh '',I6)')
     '  NEELEM(0,nr)
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      CALL CALC_NENP_1D(nb,NEELEM,NENP,NPNE,nr,ERROR,*9999)
      CALL NENXI_1D(nb,NEELEM,NENP,NPNE,nr,NXI,ERROR,
     '  *9999)
C      CALL CAP_NE(NBJ,NEELEM(0,nr),NENP,NPNE,nr,CE,XP,
C     &  ERROR,*9999)
!28/06/04 - Diameters etc now set up via command:- fem update
!            mesh geometry capillary
      
      CALL EXITS('GENPCAP')
      RETURN
 9999 CALL ERRORS('GENPCAP',ERROR)
      CALL EXITS('GENPCAP')
      RETURN 1
      END


