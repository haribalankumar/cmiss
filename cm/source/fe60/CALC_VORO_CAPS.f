      SUBROUTINE CALC_VORO_CAPS(BOUNDARY2,counter,DNP_VNP,
     '  nb,NBJ,NEELEM,NE_REPEAT,NKJ,NKJE,NPI,NPLIST3,NPNE,
     '  NPNE_ALV,NPNODE,nr,NRE,NVJE,NVJP,SE,VERTICES_XYZ,XP,ERROR,*)

C###  Subroutine: CALC_VORO_CAPS
C###  Description:

C###    Generates randomly spaced points on the surface of a sphere,
C###    calls SPHERE_DELAUNAY to triangulate the points, calculates the
C###    Voronoi dual for the Delaunay triangulation, and creates a
C###    capillary segment along each Voronoi cell edge.

C***  Created by: Kelly Burrowes, February 2002
C***  Last modified: 6th June 2002

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'pulm00.cmn'
!     Parameter List
      INTEGER BOUNDARY2(N_BOUNDARY*2,2),counter,
     '  DNP_VNP(N_ALVEOLI*MAX_TRIANGLES,4),nb,NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NE_REPEAT(0:NEREPM),
     '  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     &  NPI(N_VERT_EST),NPLIST3(0:NPM),NPNE(NNM,NBFM,NEM),
     &  NPNE_ALV(0:NP_NE,NE_R_M),NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),
     &  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM)
      REAL*8 SE(NSM,NBFM,NEM),VERTICES_XYZ(3*N_VERT_EST),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER I,ne,noelem,nonode,np,N_REPEAT,NP_REPEAT,ntri,ntri2,
     '  TRI(3),TRIANGLES(3*MAX_TRIANGLES),TRI_ADJACENT(3,MAX_TRIANGLES),
     '  TRI_REPEATED(MAX_TRIANGLES,2)
      REAL*8 RWORKING(4*MAX_TRIANGLES)
      LOGICAL VNP_FOUND,ONBOUND(MAX_TRIANGLES)
      INTEGER IERROR

      CALL ENTERS('CALC_VORO_CAPS',*9999)

      
C...  Calculate Delaunay triangulation of surface points
      CALL SPHERE_DELAUNAY(VERTICES_XYZ,TRIANGLES,RWORKING,
     '  IERROR,ERROR,*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'(''Done SPHERE_DELAUNAY'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
C...  Check that the triangulation has worked, using IERROR
      CALL ASSERT(IERROR.LE.0,'>>Error in Delaunay triangulation',
     '  ERROR,*9999)
C...  Retriangulate triangles involving boundary nodes.
      CALL REJECT_TRIANGLES(BOUNDARY2,NPI,TRIANGLES,ONBOUND,
     '  ERROR,*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'(''Done REJECT_TRIANGLES'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
C...  Create nodes and elements (for display) if enough regions set up
      IF(NRM.GE.nr+2)THEN
        CALL TRI_NE_NP(nb,NBJ,NEELEM,NKJ,NKJE,NPI,NPNE,
     '    nr+2,NRE,NVJE,NVJP,TRIANGLES,SE,VERTICES_XYZ,XP,ERROR,*9999)
      ENDIF
      N_REPEAT=0
      IF(NE_REPEAT(0).GT.0) THEN !only do if elements shared
        DO ntri=1,N_TRIANGLES !each delaunay triangle
          NP_REPEAT=0
          DO I=1,3 !each triangle vertex
            TRI(I)=NPI(TRIANGLES(3*(ntri-1)+I)) !node # in triangle
            noelem=0
            DO WHILE(NP_REPEAT.NE.3.AND.noelem.LT.NE_REPEAT(0))
              noelem=noelem+1
              ne=NE_REPEAT(noelem) !repeated alveolar element
              nonode=0
              DO WHILE(NP_REPEAT.NE.3.AND.nonode.LT.NPNE_ALV(0,ne))
                nonode=nonode+1
                np=NPNE_ALV(nonode,ne) !shared node #
                IF(TRI(I).EQ.np) NP_REPEAT=NP_REPEAT+1 !# nps repeated
              ENDDO !WHILE
            ENDDO !WHILE
          ENDDO !I
          IF(NP_REPEAT.EQ.3) THEN !all 3 nodes of triangle shared, need
C... to use voronoi node # of original triangle in place of the repeated 1
            N_REPEAT=N_REPEAT+1
            TRI_REPEATED(N_REPEAT,1)=ntri
            ntri2=0
            VNP_FOUND=.FALSE.
            DO WHILE(.NOT.VNP_FOUND.AND.ntri2.LT.counter)
              ntri2=ntri2+1
              IF(TRI(1).EQ.DNP_VNP(ntri2,1).OR.TRI(1).EQ.
     '          DNP_VNP(ntri2,2).OR.TRI(1).EQ.DNP_VNP(ntri2,3)) THEN
                IF(TRI(2).EQ.DNP_VNP(ntri2,1).OR.TRI(2).EQ.
     '            DNP_VNP(ntri2,2).OR.TRI(2).EQ.DNP_VNP(ntri2,3)) THEN
                  IF(TRI(3).EQ.DNP_VNP(ntri2,1).OR.TRI(3).EQ.
     '              DNP_VNP(ntri2,2).OR.TRI(3).EQ.DNP_VNP(ntri2,3)) THEN
                    TRI_REPEATED(N_REPEAT,2)=DNP_VNP(ntri2,4) !voro np #
                    VNP_FOUND=.TRUE.
                  ENDIF
                ENDIF
              ENDIF
            ENDDO !WHILE
            IF(.NOT.VNP_FOUND) THEN
C... must be different triangle, find triangle with same 2 points
              ntri2=0
              DO WHILE(.NOT.VNP_FOUND.AND.ntri2.LT.counter)
                ntri2=ntri2+1
                IF(TRI(1).EQ.DNP_VNP(ntri2,1).OR.TRI(1).EQ.
     '            DNP_VNP(ntri2,2).OR.TRI(1).EQ.DNP_VNP(ntri2,3)) THEN
                  IF(TRI(2).EQ.DNP_VNP(ntri2,1).OR.TRI(2).EQ.
     '              DNP_VNP(ntri2,2).OR.TRI(2).EQ.DNP_VNP(ntri2,3)) THEN
                    TRI_REPEATED(N_REPEAT,2)=DNP_VNP(ntri2,4) !voro np #
                    VNP_FOUND=.TRUE.
                  ENDIF
                ENDIF
              ENDDO !WHILE
            ENDIF
            IF(.NOT.VNP_FOUND) THEN !else can only look for 1 np common
              ntri2=0
              DO WHILE(.NOT.VNP_FOUND.AND.ntri2.LT.counter)
                ntri2=ntri2+1
                IF(TRI(1).EQ.DNP_VNP(ntri2,1).OR.TRI(1).EQ.
     '            DNP_VNP(ntri2,2).OR.TRI(1).EQ.DNP_VNP(ntri2,3)) THEN
                  TRI_REPEATED(N_REPEAT,2)=DNP_VNP(ntri2,4) !voro np #
                  VNP_FOUND=.TRUE.
                ENDIF
              ENDDO !WHILE
            ENDIF
            IF(.NOT.VNP_FOUND) THEN
              WRITE(OP_STRING,'(''WARNING: voronoi node on shared '//
     '          'face not found '')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDDO !tri
      ENDIF
C...  Create nodes and elements in region nr, where nodes are placed
C...  at the triangle circumcentres, and 1D elements join neighbouring
C...  circumcentres.
      CALL CREATE_CAPILLARY(counter,DNP_VNP,nb,NBJ,NEELEM,NKJ,NKJE,
     '  NPI,NPLIST3,NPNE,NPNODE,nr,NRE,N_REPEAT,NVJE,NVJP,TRIANGLES,
     '  TRI_ADJACENT,TRI_REPEATED,SE,VERTICES_XYZ,XP,ONBOUND,ERROR,
     &  *9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'(''Done CREATE_CAPILLARY'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      
      WRITE(OP_STRING,'('' Number of capillaries in mesh:'',I5)')
     '  NEELEM(0,nr)
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      
      
      CALL EXITS('CALC_VORO_CAPS')
      RETURN
 9999 CALL ERRORS('CALC_VORO_CAPS',ERROR)
      CALL EXITS('CALC_VORO_CAPS')
      RETURN 1
      END


