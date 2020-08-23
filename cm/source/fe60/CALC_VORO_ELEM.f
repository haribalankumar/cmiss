      SUBROUTINE CALC_VORO_ELEM(nb,NBJ,ne,NEELEM,NELEM_LIST,NELEM_LIST3,
     '  NENFVC,NENP,NFVC,NODE_CONVERT_NVC,noelem,nonode,np,NPNE,NPNODE,
     '  nr,nvc,NVC_CLASSIFY,NVCNODE,XP,ZA,ERROR,*)

C#### Subroutine: CALC_VORO_ELEM
C###  Description:
C###    CALC_VORO_ELEM converts Voronoi vertices and face
C###    boundaries to nodes and elements.  These are written into the
C###    same region that the original Delaunay elements were in - the
C###    Delaunay nodes and elements are discarded.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER nb,NBJ(NJM,NEM),ne,NEELEM(0:NE_R_M,0:NRM),
     '  NELEM_LIST(0:NVCM*5),NELEM_LIST3(0:NVCM),NENFVC(0:NFVCM,NFVM),
     '  NENP(NPM,0:NEPM,0:NRM),NFVC(2,0:NFVCM,NVCM),
     '  NODE_CONVERT_NVC(NE_R_M),noelem,nonode,np,NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,nvc,NVC_CLASSIFY(NVCM),
     '  NVCNODE(2,NP_R_M)
      REAL*8 XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER ncount,ned,nf,nfvl,nj,nofv,NPLIST_LOCAL(20),
     '  np_start,np1,nvc_adjacent,nvc_node,N_VERT

      CALL ENTERS('CALC_VORO_ELEM',*9999)

      DO nfvl=1,NFVC(1,0,nvc) !for each face in cell
        nvc_node=NFVC(1,nfvl,nvc) !local node number at centre of adjacent cell
        nf=NFVC(2,nfvl,nvc) !global face #
        nvc_adjacent=NVCNODE(2,nvc_node) !cell number of adjacent cell
        IF((NVC_CLASSIFY(nvc_adjacent).EQ.1).OR.
     '    (NVC_CLASSIFY(nvc_adjacent).EQ.3))THEN
C         Won't create elements on face adjacent to the duct cell, or on
C         shared face that has already been done. Will=3 if the alveolus
C         has already been done.
        ELSE
C         Adjacent cell is an alveolus (or unclassified), so create a
C         single element -covered 'face' between them (the alveolar
C         septum). 
          
C         Make nodes at vertices.
          np_start=np
          DO nofv=1,NENFVC(0,nf)!max # new nodes = #vertices
C           Make sure that the vertices have not already been converted
C           to nodes.  Do this by checking NODE_CONVERT_NVC
            ned=NENFVC(nofv,nf) !Delaunay element # for vertex
            IF(NODE_CONVERT_NVC(ned).EQ.0)THEN
C             Create a new node, and store in a local node list.
              nonode=nonode+1
              np=np+1
              CALL ASSERT(np.LE.NPM,'Increase NPM',
     '          ERROR,*9999)
              NPNODE(nonode,nr)=np
              NENP(np,0,nr)=0
              NPLIST_LOCAL(nofv)=np
            ELSE
C             Has already been set up, so find out which np to use.
              NPLIST_LOCAL(nofv)=NODE_CONVERT_NVC(ned)
            ENDIF
          ENDDO !nofv

          noelem=noelem+1
          ne=ne+1
          CALL ASSERT(ne.LE.NEM,'Increase NEM',ERROR,*9999)
          NEELEM(noelem,nr)=ne
          NELEM_LIST(0)=NELEM_LIST(0)+1
          NELEM_LIST(NELEM_LIST(0))=ne
          NELEM_LIST3(0)=NELEM_LIST3(0)+1
          NELEM_LIST3(NELEM_LIST3(0))=ne

          np1=NPLIST_LOCAL(1)
          NPNE(1,nb,ne)=np1
          NENP(np1,0,nr)=NENP(np1,0,nr)+1
          NENP(np1,NENP(np1,0,nr),nr)=ne

          np1=NPLIST_LOCAL(2)
          NPNE(2,nb,ne)=np1
          NENP(np1,0,nr)=NENP(np1,0,nr)+1
          NENP(np1,NENP(np1,0,nr),nr)=ne

          IF(NENFVC(0,nf).EQ.3)THEN
            np1=NPLIST_LOCAL(3)
          ELSE
            np1=NPLIST_LOCAL(4)
          ENDIF
          NPNE(3,nb,ne)=np1
          NENP(np1,0,nr)=NENP(np1,0,nr)+1
          NENP(np1,NENP(np1,0,nr),nr)=ne

          np1=NPLIST_LOCAL(3)
          NPNE(4,nb,ne)=np1
          NENP(np1,0,nr)=NENP(np1,0,nr)+1
          NENP(np1,NENP(np1,0,nr),nr)=ne
          
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            NBJ(nj,ne)=nb
          ENDDO !nj

          N_VERT=4
          DO WHILE(N_VERT.LT.NENFVC(0,nf))
            noelem=noelem+1
            ne=ne+1
            CALL ASSERT(ne.LE.NEM,'Increase NEM',ERROR,*9999)
            NEELEM(noelem,nr)=ne
            NELEM_LIST(0)=NELEM_LIST(0)+1
            NELEM_LIST(NELEM_LIST(0))=ne
            NELEM_LIST3(0)=NELEM_LIST3(0)+1
            NELEM_LIST3(NELEM_LIST3(0))=ne

            np1=NPLIST_LOCAL(N_VERT)
            NPNE(2,nb,ne)=np1
            NENP(np1,0,nr)=NENP(np1,0,nr)+1
            NENP(np1,NENP(np1,0,nr),nr)=ne
            N_VERT=N_VERT+1
            
            np1=NPLIST_LOCAL(N_VERT)
            NPNE(1,nb,ne)=np1
            NENP(np1,0,nr)=NENP(np1,0,nr)+1
            NENP(np1,NENP(np1,0,nr),nr)=ne
            N_VERT=N_VERT+1

            IF(N_VERT.LE.NENFVC(0,nf))THEN
              np1=NPLIST_LOCAL(N_VERT)
            ELSE
              np1=NPLIST_LOCAL(1)
            ENDIF
            NPNE(3,nb,ne)=np1
            NENP(np1,0,nr)=NENP(np1,0,nr)+1
            NENP(np1,NENP(np1,0,nr),nr)=ne

            np1=NPLIST_LOCAL(1)
            NPNE(4,nb,ne)=np1
            NENP(np1,0,nr)=NENP(np1,0,nr)+1
            NENP(np1,NENP(np1,0,nr),nr)=ne

            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              NBJ(nj,ne)=nb
            ENDDO !nj
          ENDDO!WHILE
          
C         Copy vertex coordinates into XP.  The circum-centres of
C         the Delaunay elements are the Voronoi vertices, and these
C         become the nodes.
          ncount=0
          DO nofv=1,NENFVC(0,nf) !for each new face vertex
            ned=NENFVC(nofv,nf) !Delaunay element # for vertex
            IF(NODE_CONVERT_NVC(ned).EQ.0)THEN
              ncount=ncount+1
              DO nj=1,NJT
                XP(1,1,nj,np_start+ncount)=ZA(1,nj,1,ned) !coordinates of vertex
              ENDDO !nj
              NODE_CONVERT_NVC(ned)=np_start+ncount !so that we know it has been done
            ENDIF !NODE_CONVERT
          ENDDO !nofv
          
        ENDIF !NVC_CLASSIFY
      ENDDO !nfvl

      CALL EXITS('CALC_VORO_ELEM')
      RETURN
 9999 CALL ERRORS('CALC_VORO_ELEM',ERROR)
      CALL EXITS('CALC_VORO_ELEM')
      RETURN 1
      END


