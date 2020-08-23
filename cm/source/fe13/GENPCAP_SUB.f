      SUBROUTINE GENPCAP_SUB(BOUNDARY,BOUNDARY2,BOUNDARY_TEMP,DNP_VNP,
     &  nb,NBJ,NEELEM,NELIST,NELIST2,NE_REPEAT,NKJ,NKJE,NPI,
     &  NP_INTERFACE,NPLIST,NPLIST2,NPLIST3,NPNE,NPNE_ALV,NPNODE,nr,NRE,
     &  NVJE,NVJP,C,radius,SE,VERTICES_XYZ,XP,REGULAR,ERROR,*)
      
C#### Subroutine: GENPCAP_SUB
C###  Description:
C###    GENPCAP_SUB generates the Voronoi capillary mesh over
C###    alveolar surfaces

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'genmesh.cmn'
      INCLUDE 'pulm00.cmn'
!     Parameter list
      INTEGER nb,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     &  NELIST2(0:NEM),NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     &  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     &  NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),
     &  NVJP(NJM,NPM)
      REAL*8 C(MAX_ALVEOLI,NJT),radius,SE(NSM,NBFM,NEM),
     '  XP(NKM,NVM,NJM,NPM)
      LOGICAL REGULAR
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER alveoli,B,BOUNDARY(0:N_BOUNDARY*N_ALVEOLI*2,2),
     '  BOUNDARY2(N_BOUNDARY*2,2),BOUNDARY_TEMP(N_BOUNDARY*2,2),counter,
     '  DNP_VNP(N_ALVEOLI*MAX_TRIANGLES,4),I,IB,IBEG,IEND,
     '  j,k,group,ne,ne2,NE_REPEAT(0:NEREPM),nj,noelem,noelem2,
     '  nonode,nonode2,np,NPI(N_VERT_EST),NPLIST(0:NPM),NPLIST2(0:NPM),
     '  NPLIST3(0:NPM),NPNE_ALV(0:NP_NE,NE_R_M)
      REAL*8 alpha1,alpha2,BNODE,centre(NJT),IBNODE,TOL,
     '  VERTICES_XYZ(3*N_VERT_EST),XP_point(NJT),X_TRANS(NJT)
      LOGICAL PROJECT1,REPEAT
      CHARACTER CDATA(1)*100,LABEL*30,STRING*255

      CALL ENTERS('GENPCAP_SUB',*9999)
      write(*,*) 'enter sub'

      counter=0
      DO noelem=1,NE_R_M
        NPNE_ALV(0,noelem)=0 !initialise totals in array
      ENDDO !noelem
      IF(NRM.LT.nr+2) THEN
        WRITE(OP_STRING,'('' Increase NRM (NRM>=nr+2)'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
C...   Generate points on surface of (unit) sphere
      CALL SPHERE_POINTS(BOUNDARY,NKJ,NPNODE,nr+1,NVJP,VERTICES_XYZ,
     '  XP,REGULAR,ERROR,*9999)
      NPLIST(0)=0
      DO nonode=1,NPNODE(0,nr+1) !currently always all nodes in region
        np=NPNODE(nonode,nr+1) !groups nodes created in SPHERE_POINTS
        NPLIST(0)=NPLIST(0)+1
        NPLIST(NPLIST(0))=np
      ENDDO !nonode  !groups nodes
      STRING='sphere' !should have this input on command line
      CALL GRNODE_SUB(NPLIST,STRING,.TRUE.,ERROR,*9999) !AS=.TRUE.=1
      BOUNDARY(0,1)=0
      BOUNDARY(0,2)=0

      DO group=1,N_ALVEOLI !each alveolar group
        PROJECT1=.TRUE.
        CALL CUPPER(ALVEOLI_NAMES(group),STRING)
        CDATA(1)='ELEMENTS' !puts element group list in NELIST
        CALL PARSILG(NELIST,NEM,CDATA(1),STRING,ERROR,*9999)
        NE_REPEAT(0)=0 !initialise
        alveoli=0                  
        IF(.NOT.FIRST_PROJECT) THEN !find elements common to groups
C... to find elements which are shared by alveoli
          LABEL='start' 
          DO WHILE(LABEL.NE.STRING.AND.alveoli.LT.N_ALVEOLI)
C... for now only do previous groups, but could use the fact that if a
C... alveolar group has been projected a group "alveolar group_cap"
C... will be created
            alveoli=alveoli+1
            CALL CUPPER(ALVEOLI_NAMES(alveoli),LABEL)
            CALL STRING_TRIM(LABEL,IBEG,IEND)
            IF(LABEL.NE.STRING) THEN
              CDATA(1)='ELEMENTS'
              CALL PARSILG(NELIST2,NEM,CDATA(1),LABEL,ERROR,*9999)
              DO noelem=1,NELIST(0) !element group currently using
                ne=NELIST(noelem)
                noelem2=0
                DO WHILE(noelem2.LT.NELIST2(0).AND.
     '            NE_REPEAT(NE_REPEAT(0)).NE.ne)
                  noelem2=noelem2+1
                  ne2=NELIST2(noelem2)
                  IF(ne.EQ.ne2) THEN !ne shared by another alveolus
                    NE_REPEAT(0)=NE_REPEAT(0)+1
                    NE_REPEAT(NE_REPEAT(0))=ne
                  ENDIF
                ENDDO !noelem2
              ENDDO !noelem
            ENDIF !LABEL.NE.STRING
          ENDDO !WHILE
        ENDIF !.NOT.FIRST
        DO nj=1,NJT
          centre(nj)=C(group,nj) !centre of alveolar elem group
        ENDDO
        CALL CAP_PROJ(group,NBJ,NELIST,NE_REPEAT,NKJ,NPNE,
     '    NPNE_ALV,NPLIST,NPNODE,nr+2,NVJP,centre,XP,PROJECT1,
     '    ERROR,*9999)
C... first projection for group creates mesh nodes
        PROJECT1=.FALSE.
C... CAP_PROJ2 projects points back out to sphere surface
        CALL CAP_PROJ2(NELIST,NPNE_ALV,centre,XP,ERROR,*9999)
C... move delaunay point co-ords from XP to VERTICES_XYZ & set
C... up new BOUNDARY array BOUNDARY2
        NPLIST2(0)=0
        I=0
        j=0 !indexing for BOUNDARY2 array
        k=0
C... firstly check if any shared nodes are already in BOUNDARY2
        DO noelem=1,NE_REPEAT(0)
          ne=NE_REPEAT(noelem)
          DO nonode=1,NPNE_ALV(0,ne)
            np=NPNE_ALV(nonode,ne) !nodes used again
            DO nonode2=1,BOUNDARY(0,1)
              IB=BOUNDARY(nonode2,1)
              IF(np.EQ.IB) THEN
                j=j+1
                BOUNDARY_TEMP(j,1)=np !temporary storage
              ENDIF
            ENDDO
            DO nonode2=1,BOUNDARY(0,2)
              B=BOUNDARY(nonode2,2)
              IF(np.EQ.B) THEN
                k=k+1
                BOUNDARY_TEMP(k,2)=np
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        DO nonode=1,j !move into correct array
          BOUNDARY2(nonode,1)=BOUNDARY_TEMP(nonode,1)
        ENDDO
        DO nonode=1,k
          BOUNDARY2(nonode,2)=BOUNDARY_TEMP(nonode,2)
        ENDDO
        TOL=1.d-6
        alpha1=(PI-low_limit)-z_bound*PI !angle for IB nodes
        alpha2=(PI-low_limit)+z_bound*PI !angle for B nodes
        IBNODE=radius*DCOS(alpha1) !IB node position in z
        BNODE=radius*DCOS(alpha2) !B node position in z
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          REPEAT=.FALSE.
          DO noelem2=1,NE_REPEAT(0)
            ne2=NE_REPEAT(noelem2)
            IF(ne.EQ.ne2) REPEAT=.TRUE.
          ENDDO
          DO nonode=1,NPNE_ALV(0,ne)
            np=NPNE_ALV(nonode,ne)
            NPLIST2(0)=NPLIST2(0)+1
            NPLIST2(NPLIST2(0))=np
            IF(VORONOI_ALV) THEN
              DO nj=1,NJT
                XP_point(nj)=XP(1,1,nj,np)
              ENDDO
              CALL COORD_TRANS_INV(ANGLE_X,ANGLE_Y,ANGLE_Z,X_TRANS,
     '          centre,XP_point,ERROR,*9999)
            ENDIF
            DO nj=1,NJT
              I=I+1
C... sphere_delaunay assumes unit sphere at origin therefore undo
C... translation & rotation of points done in CAP_PROJ, re-done in CAP_PROJ
              IF(VORONOI_ALV) THEN !undo translation & rotation
C... NB/ XP co-ords are changed in tri_ne_np
                VERTICES_XYZ(I)=X_TRANS(nj)
              ELSE
                VERTICES_XYZ(I)=XP(1,1,nj,np)-centre(nj) !co-ord of np
              ENDIF
            ENDDO !nj
            IF(.NOT.REPEAT) THEN !repeated B & IB nodes already done
              IF(VERTICES_XYZ(I).GT.(BNODE-TOL).AND.
     '          VERTICES_XYZ(I).LT.(BNODE+TOL)) THEN
                k=k+1
                BOUNDARY2(k,2)=np !B node
              ELSE IF(VERTICES_XYZ(I).GT.(IBNODE-TOL).AND.
     '            VERTICES_XYZ(I).LT.(IBNODE+TOL)) THEN
                j=j+1
                BOUNDARY2(j,1)=np !IB node
              ENDIF
            ENDIF
            NPI(I/3)=np !stores np # for co-ords in VERTICES_XYZ(I)
          ENDDO !nonode
        ENDDO !noelem
        N_BOUND=k !#N_BOUNDARY may increase due to shared faces
        N_IBOUND=j !# internal boundary points
        DO nonode=1,N_IBOUND
          BOUNDARY(0,1)=BOUNDARY(0,1)+1
          BOUNDARY(BOUNDARY(0,1),1)=BOUNDARY2(nonode,1)
        ENDDO
        DO nonode=1,N_BOUND
          BOUNDARY(0,2)=BOUNDARY(0,2)+1
          BOUNDARY(BOUNDARY(0,2),2)=BOUNDARY2(nonode,2)
        ENDDO
        N_VERTICES=NPLIST2(0) !final # of vertices to triangulate
        FIRST_PROJECT=.FALSE.
        DO J=1,3
          DNP_VNP(1,J)=1 !just to start off
        ENDDO
        CALL CALC_VORO_CAPS(BOUNDARY2,counter,DNP_VNP,nb,NBJ,
     '    NEELEM,NE_REPEAT,NKJ,NKJE,NPI,NPLIST3,NPNE,NPNE_ALV,
     '    NPNODE,nr,NRE,NVJE,NVJP,SE,VERTICES_XYZ,XP,ERROR,*9999)
C... projects points back onto alveolar surface, after triangulation,
C... without making new nodes, only changes XP for nodes
        IF(NRM.GE.nr+2) THEN
          CALL CAP_PROJ(group,NBJ,NELIST,NE_REPEAT,
     '      NKJ,NPNE,NPNE_ALV,NPLIST2,NPNODE,nr+2,NVJP,centre,XP,
     '      PROJECT1,ERROR,*9999) !region = nr+2
        ENDIF
C... above projects delaunay nodes (required for adjacent spheres)
        CALL CAP_PROJ(group,NBJ,NELIST,NE_REPEAT,
     '    NKJ,NPNE,NPNE_ALV,NPLIST3,NPNODE,nr,NVJP,centre,XP,
     '    PROJECT1,ERROR,*9999) !projects Voronoi nodes (region = nr)
        CALL STRING_TRIM(STRING,IBEG,IEND) !group name for capillaries
        CALL APPENDC(IEND,'_cap',STRING) ! is "(host group name)_cap"
C... group nodes & elems for an alveolus group
C       CALL GRNODE_SUB(NPLIST2,STRING,.TRUE.,ERROR,*9999) !AS=.TRUE.
C       CALL GRELEM_SUB(NELIST3,STRING,.TRUE.,ERROR,*9999)
      ENDDO !group
      CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)

      CALL EXITS('GENPCAP_SUB')
      RETURN
 9999 CALL ERRORS('GENPCAP_SUB',ERROR)
      CALL EXITS('GENPCAP_SUB')
      RETURN 1
      END


