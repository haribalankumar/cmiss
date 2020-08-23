      SUBROUTINE CREATE_CAPILLARY(counter,DNP_VNP,nb,NBJ,
     '  NEELEM,NKJ,NKJE,NPI,NPLIST3,NPNE,NPNODE,nr,NRE,N_REPEAT,NVJE,
     &  NVJP,TRIANGLES,TRI_ADJACENT,TRI_REPEATED,SE,VERTICES_XYZ,XP,
     '  ONBOUND,ERROR,*)

C###  Subroutine: CREATE_CAPILLARY
C###  Description:
C###    Creates a capillary mesh by joinng the circumcentres of Voronoi
C###    cells.

C***  Created by: Kelly Burrowes, February 2002
C***  Last modified: June 2002

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'pulm00.cmn'
!     Parameter List
      INTEGER counter,DNP_VNP(N_ALVEOLI*MAX_TRIANGLES,4),nb,
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NKJ(NJM,NPM),
     &  NKJE(NKM,NNM,NJM,NEM),NPI(N_VERT_EST),NPLIST3(0:NPM),
     &  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),N_REPEAT,
     &  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),TRIANGLES(3*MAX_TRIANGLES),
     &  TRI_ADJACENT(3,MAX_TRIANGLES),TRI_REPEATED(MAX_TRIANGLES,2)
      REAL*8 SE(NSM,NBFM,NEM),VERTICES_XYZ(3*N_VERT_EST),
     &  XP(NKM,NVM,NJM,NPM)
      LOGICAL ONBOUND(MAX_TRIANGLES)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER I,IERROR,J,ne,nj,noelem,nonode,nonode2,np,
     '  NPLIST(0:NPM),ntri,ntri2,ntri3,ntri4,tri1,tri2,TRI_A(3),
     '  TRI_B(3)
      REAL*8 XYZ_CENTRE1(3),XYZ_TRI(3,3)
      LOGICAL REPEAT

      CALL ENTERS('CREATE_CAPILLARY',*9999)

C     Initialise adjacency array
      DO ntri=1,N_TRIANGLES
        DO I=1,3
          TRI_ADJACENT(I,ntri)=0 !initialise
        ENDDO !I
      ENDDO !ntri

C     Create nodes at circumcentres of (non-B) triangles, and set up
C     adjacency array TRI-ADJACENT
      nonode=NPNODE(0,nr) !highest local node in capillary region
      np=NPT(0) !highest global node number
      I=1
      DO ntri=1,N_TRIANGLES
        REPEAT=.FALSE.
        ntri3=0
        DO WHILE(ntri3.LT.N_REPEAT.AND..NOT.REPEAT)
          ntri3=ntri3+1
          IF(ntri.EQ.TRI_REPEATED(ntri3,1)) REPEAT=.TRUE.
        ENDDO !WHILE
        IF(.NOT.REPEAT) THEN
          counter=counter+1 !counts total # triangles for all alveoli
          DO J=1,3 !for each triangle vertex
            TRI_A(J)=TRIANGLES(I) !the vertex #
            DNP_VNP(counter,J)=NPI(TRI_A(J)) !delaunay node #
            XYZ_TRI(1,J)=VERTICES_XYZ(3*TRI_A(J)-2) !x
            XYZ_TRI(2,J)=VERTICES_XYZ(3*TRI_A(J)-1) !y
            XYZ_TRI(3,J)=VERTICES_XYZ(3*TRI_A(J)) !z
            I=I+1 !next vertex
          ENDDO !J

C... calculates circumcenter of center for boundary triangles
          CALL CIRCUM2(IERROR,ntri,XYZ_CENTRE1,XYZ_TRI,ONBOUND,
     '      ERROR,*9999)
          np=np+1
          nonode=nonode+1
          CALL ASSERT(nonode.LE.NP_R_M,'>>Increase NP_R_M',
     '      ERROR,*9999)
          CALL ASSERT(np.LE.NPM,'>>Increase NPM',ERROR,*9999)
          NPNODE(nonode,nr)=np
          NPLIST(ntri)=np
          DNP_VNP(counter,4)=np
C... DNP_VNP stores 3 delaunay nodes of triangles & corresponding
C... voronoi node for that triangle

 !put circumcentre coordinates into XP
          DO nj=1,3 !assuming 3D, must be for sphere
            XP(1,1,nj,np)=XYZ_CENTRE1(nj)
          ENDDO !nj
        ELSE IF(REPEAT) THEN
          DO J=1,3 !for each triangle vertex
            TRI_A(J)=TRIANGLES(I) !the vertex #
            XYZ_TRI(1,J)=VERTICES_XYZ(3*TRI_A(J)-2) !x
            XYZ_TRI(2,J)=VERTICES_XYZ(3*TRI_A(J)-1) !y
            XYZ_TRI(3,J)=VERTICES_XYZ(3*TRI_A(J)) !z
            I=I+1 !next vertex
          ENDDO !J
          NPLIST(ntri)=TRI_REPEATED(ntri3,2) !voronoi np # already there
        ENDIF !.NOT.REPEAT !only create new np if not one already there
 !set up triangle adjacency array
        DO ntri2=ntri+1,N_TRIANGLES
          TRI_B(1)=TRIANGLES(3*ntri2-2) !vertex #s of next triangle
          TRI_B(2)=TRIANGLES(3*ntri2-1)
          TRI_B(3)=TRIANGLES(3*ntri2)
          IF(TRI_B(1).EQ.TRI_A(1))THEN
            IF(TRI_B(2).EQ.TRI_A(2).OR.TRI_B(2).EQ.TRI_A(3))THEN
              TRI_ADJACENT(1,ntri2)=ntri !1-2 edge
              IF(TRI_B(2).EQ.TRI_A(2))THEN
                TRI_ADJACENT(1,ntri)=ntri2 !1-2 edge
              ELSE
                TRI_ADJACENT(3,ntri)=ntri2 !1-3 edge
              ENDIF
            ELSE IF(TRI_B(3).EQ.TRI_A(2).OR.TRI_B(3)
     '          .EQ.TRI_A(3))THEN
              TRI_ADJACENT(3,ntri2)=ntri !1-3 edge
              IF(TRI_B(3).EQ.TRI_A(2))THEN
                TRI_ADJACENT(1,ntri)=ntri2 !1-2 edge
              ELSE
                TRI_ADJACENT(3,ntri)=ntri2 !1-3 edge
              ENDIF
            ENDIF
          ELSE IF(TRI_B(2).EQ.TRI_A(1))THEN
            IF(TRI_B(1).EQ.TRI_A(2).OR.TRI_B(1).EQ.TRI_A(3))THEN
              TRI_ADJACENT(1,ntri2)=ntri !1-2 edge
              IF(TRI_B(1).EQ.TRI_A(2))THEN
                TRI_ADJACENT(1,ntri)=ntri2 !1-2 edge
              ELSE
                TRI_ADJACENT(3,ntri)=ntri2 !1-3 edge
              ENDIF
            ELSE IF(TRI_B(3).EQ.TRI_A(2).OR.TRI_B(3)
     '          .EQ.TRI_A(3))THEN
              TRI_ADJACENT(2,ntri2)=ntri !2-3 edge
              IF(TRI_B(3).EQ.TRI_A(2))THEN
                TRI_ADJACENT(1,ntri)=ntri2 !1-2 edge
              ELSE
                TRI_ADJACENT(3,ntri)=ntri2 !1-3 edge
              ENDIF
            ENDIF
          ELSE IF(TRI_B(3).EQ.TRI_A(1))THEN
            IF(TRI_B(1).EQ.TRI_A(2).OR.TRI_B(1).EQ.TRI_A(3))THEN
              TRI_ADJACENT(3,ntri2)=ntri !1-3 edge
              IF(TRI_B(1).EQ.TRI_A(2))THEN
                TRI_ADJACENT(1,ntri)=ntri2 !1-2 edge
              ELSE
                TRI_ADJACENT(3,ntri)=ntri2 !1-3 edge
              ENDIF
            ELSE IF(TRI_B(2).EQ.TRI_A(2).OR.TRI_B(2)
     '          .EQ.TRI_A(3))THEN
              TRI_ADJACENT(2,ntri2)=ntri !2-3 edge
              IF(TRI_B(2).EQ.TRI_A(2))THEN
                TRI_ADJACENT(1,ntri)=ntri2 !1-2 edge
              ELSE
                TRI_ADJACENT(3,ntri)=ntri2 !1-3 edge
              ENDIF
            ENDIF
          ELSE IF(TRI_B(1).EQ.TRI_A(2))THEN
            IF(TRI_B(2).EQ.TRI_A(3).OR.TRI_B(3).EQ.TRI_A(3))THEN
              TRI_ADJACENT(2,ntri)=ntri2 !2-3 edge
              IF(TRI_B(2).EQ.TRI_A(3))THEN
                TRI_ADJACENT(1,ntri2)=ntri !1-2 edge
              ELSE
                TRI_ADJACENT(3,ntri2)=ntri !1-3 edge
              ENDIF
            ENDIF
          ELSE IF(TRI_B(2).EQ.TRI_A(2))THEN
            IF(TRI_B(1).EQ.TRI_A(3).OR.TRI_B(3).EQ.TRI_A(3))THEN
              TRI_ADJACENT(2,ntri)=ntri2 !2-3 edge
              IF(TRI_B(1).EQ.TRI_A(3))THEN
                TRI_ADJACENT(1,ntri2)=ntri !1-2 edge
              ELSE
                TRI_ADJACENT(2,ntri2)=ntri !2-3 edge
              ENDIF
            ENDIF
          ELSE IF(TRI_B(3).EQ.TRI_A(2))THEN
            IF(TRI_B(1).EQ.TRI_A(3).OR.TRI_B(2).EQ.TRI_A(3))THEN
              TRI_ADJACENT(2,ntri)=ntri2 !2-3 edge
              IF(TRI_B(1).EQ.TRI_A(3))THEN
                TRI_ADJACENT(1,ntri2)=ntri !1-2 edge
              ELSE
                TRI_ADJACENT(2,ntri2)=ntri !2-3 edge
              ENDIF
            ENDIF
          ENDIF
        ENDDO !ntri2
      ENDDO !ntri
C     Create elements joining the neighbouring circumcentres
      noelem=NEELEM(0,nr) !highest element # in capillary region
      ne=NET(0) !highest global element #
      DO ntri=1,N_TRIANGLES !for each Delaunay triangle
        DO I=1,3 !for each possible adjacent triangle
          REPEAT=.FALSE.
          ntri2=TRI_ADJACENT(I,ntri) !adjacent triangle #
          DO ntri3=1,N_REPEAT !all repeated triangles
            tri1=TRI_REPEATED(ntri3,1) !repeated triangle
            IF(ntri.EQ.tri1) THEN
              DO ntri4=1,N_REPEAT
                tri2=TRI_REPEATED(ntri4,1)
                IF(ntri2.EQ.tri2) REPEAT=.TRUE. !both tri's repeated
              ENDDO !ntri3
            ENDIF
          ENDDO !ntri3
          IF(ntri2.GT.ntri.AND..NOT.REPEAT) THEN
C...only create element if adj bigger & element not already created
C... i.e triangle not shared by adjacent alveoli
            noelem=noelem+1
            ne=ne+1
            CALL ASSERT(noelem.LE.NE_R_M,'>>Increase NE_R_M',
     '        ERROR,*9999)
            CALL ASSERT(ne.LE.NEM,'>>Increase NEM',ERROR,*9999)
            NEELEM(noelem,nr)=ne
            NPNE(1,nb,ne)=NPLIST(ntri) !first node
            NPNE(2,nb,ne)=NPLIST(ntri2) !second node
            CALL GN1DNEJ(nb,NBJ,ne,NKJ,NKJE,NPNE(1,nb,ne),NPNE(2,nb,ne),
     &        nr,NRE,NVJE,NVJP,SE,ERROR,*9999)
          ENDIF !ntri2.GT.ntri
        ENDDO !I
      ENDDO !ntri
C... group voronoi nodes & elems just created for alveolar group
      NPLIST3(0)=0
      DO nonode2=NPT(0)+1,np
        NPLIST3(0)=NPLIST3(0)+1
        NPLIST3(NPLIST3(0))=nonode2
      ENDDO
      NEELEM(0,nr)=noelem
      NPNODE(0,nr)=nonode
      NEELEM(0,0)=ne
      NPNODE(0,0)=np
      NET(nr)=noelem
      NPT(nr)=nonode
      NET(0)=ne
      NPT(0)=np

      CALL EXITS('CREATE_CAPILLARY')
      RETURN
 9999 CALL ERRORS('CREATE_CAPILLARY',ERROR)
      CALL EXITS('CREATE_CAPILLARY')
      RETURN 1
      END


