      SUBROUTINE EDGEFLIP(NBJ,ne1,ne2,NPNE,nr,NVJE,NXI,XP,ZA,ERROR,
     '  *)

C#### Subroutine: EDGEFLIP
C###  Description:
C###    Given a complex of two triangles with an adjoining edge, that is
C###    not Delaunay, the adjoining edge is deleted and the new edge
C###    becomes the opposite diagonal

C***    IE The basic thing that happens is:
C***
C***         /\                 /|\
C***        /  \               / | \
C***       /    \             /  |  \
C***       ------      =>        |
C***       \    /             \  |  /
C***        \  /               \ | /
C***         \/                 \|/

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),ne1,ne2,NPNE(NNM,NBFM,NEM),nr,
     '  NVJE(NNM,NBFM,NJM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb1,nn,old_ne1_nodes(3),old_ne1_nbors(3),
     '  old_ne2_nodes(3),old_ne2_nbors(3),new_ne1_nodes(3),
     '  new_ne1_nbors(3),new_ne2_nodes(3),new_ne2_nbors(3),
     '  e1_local_node_points_to_e2,e1_global_node_points_to_e2,
     '  e2_local_node_points_to_e1,e2_global_node_points_to_e1,
     '  e1_local_node_dropped,e1_global_node_dropped,
     '  e1_local_node_kept,e1_global_node_kept,
     '  e2_local_node_dropped,
     '  e2_local_node_kept

      CALL ENTERS('EDGEFLIP',*9999)

      nb1=NBJ(1,ne1)
c      nb2=NBJ(1,ne2)

C     ..First make a local copy of the local nodes
      DO nn=1,NNT(nb1)
        old_ne1_nodes(nn)=NPNE(nn,nb1,ne1)
        old_ne1_nbors(nn)=NXI(0,nn,ne1)
        CALL ASSERT(old_ne1_nbors(nn).NE.0,'>>Cannot flip boundary '//
     '    'nodes',ERROR,*9999)
        old_ne2_nodes(nn)=NPNE(nn,nb1,ne2)
        old_ne2_nbors(nn)=NXI(0,nn,ne2)
        CALL ASSERT(old_ne2_nbors(nn).NE.0,'>>Cannot flip boundary '//
     '    'nodes',ERROR,*9999)
      ENDDO

C     ..Find each elements nodes that point to the adjacent element
      DO nn=1,NNT(nb1)
        IF(old_ne1_nbors(nn).EQ.ne2) THEN
          e1_local_node_points_to_e2=nn
          e1_global_node_points_to_e2=old_ne1_nodes(nn)
        ENDIF
        IF(old_ne2_nbors(nn).EQ.ne1) THEN
          e2_local_node_points_to_e1=nn
          e2_global_node_points_to_e1=old_ne2_nodes(nn)
        ENDIF
      ENDDO

C     ..From the other 2 nodes of e1, pick a node to be dropped from e1
      DO nn=1,NNT(nb1)
        IF(nn.NE.e1_local_node_points_to_e2) THEN
          e1_local_node_dropped=nn
          e1_global_node_dropped=old_ne1_nodes(nn)
          GOTO 10
        ENDIF
      ENDDO
 10   CONTINUE

C     ..The remaining node will be the node kept
      DO nn=1,NNT(nb1)
        IF(nn.NE.e1_local_node_points_to_e2.AND.
     '    nn.NE.e1_local_node_dropped) THEN
          e1_local_node_kept=nn
          e1_global_node_kept=old_ne1_nodes(nn)
          GOTO 20
        ENDIF
      ENDDO
 20   CONTINUE

C     ..Search for nodes of ne2 to equal kept & dropped nodes of ne1
      DO nn=1,NNT(nb1)
        IF(old_ne2_nodes(nn).EQ.e1_global_node_dropped) THEN
          e2_local_node_kept=nn
c          e2_global_node_kept=old_ne2_nodes(nn)
        ELSEIF(old_ne2_nodes(nn).EQ.e1_global_node_kept) THEN
          e2_local_node_dropped=nn
c          e2_global_node_dropped=old_ne2_nodes(nn)
        ENDIF
      ENDDO

C     ..Copy the old mesh stuff into a new copy
      DO nn=1,NNT(nb1)
        new_ne1_nodes(nn)=old_ne1_nodes(nn)
        new_ne1_nbors(nn)=old_ne1_nbors(nn)
        new_ne2_nodes(nn)=old_ne2_nodes(nn)
        new_ne2_nbors(nn)=old_ne2_nbors(nn)
      ENDDO

C     ..Reassign the nodes
      new_ne1_nodes(e1_local_node_dropped)=e2_global_node_points_to_e1
      new_ne2_nodes(e2_local_node_dropped)=e1_global_node_points_to_e2
      new_ne1_nbors(e1_local_node_kept)=ne2
      new_ne2_nbors(e2_local_node_kept)=ne1
      new_ne1_nbors(e1_local_node_points_to_e2)=
     '  old_ne2_nbors(e2_local_node_kept)
      new_ne2_nbors(e2_local_node_points_to_e1)=
     '  old_ne1_nbors(e1_local_node_kept)

C     ..Next 2 statements are not needed, but kept for clarity
      new_ne1_nbors(e1_local_node_dropped)=
     '  old_ne1_nbors(e1_local_node_dropped)
      new_ne2_nbors(e2_local_node_dropped)=
     '  old_ne2_nbors(e2_local_node_dropped)

C     ..Update neighbour's neighbours
      DO nn=1,NNT(nb1)
        IF(NXI(0,nn,old_ne1_nbors(e1_local_node_kept)).EQ.ne1) THEN
          NXI(0,nn,old_ne1_nbors(e1_local_node_kept))=ne2
          GOTO 30
        ENDIF
      ENDDO
      WRITE(ERROR,'('' Unable to find element '',I6,'' opposite '//
     '  'element'',I6)') ne1,old_ne1_nbors(e1_local_node_kept)
      GOTO 9999
 30   CONTINUE
      DO nn=1,NNT(nb1)
        IF(NXI(0,nn,old_ne2_nbors(e2_local_node_kept)).EQ.ne2) THEN
          NXI(0,nn,old_ne2_nbors(e2_local_node_kept))=ne1
          GOTO 40
        ENDIF
      ENDDO
      WRITE(ERROR,'('' Unable to find element '',I6,'' opposite '//
     '  'element'',I6)') ne2,old_ne2_nbors(e2_local_node_kept)
      GOTO 9999
 40   CONTINUE

C     ..Copy the new information into the main arrays
      DO nn=1,NNT(nb1)
        NPNE(nn,nb1,ne1)=new_ne1_nodes(nn)
        NXI(0,nn,ne1)=new_ne1_nbors(nn)
        NPNE(nn,nb1,ne2)=new_ne2_nodes(nn)
        NXI(0,nn,ne2)=new_ne2_nbors(nn)
      ENDDO

C     ..Have two new elements, update circumcircle stuff
      CALL CIRCUM(NBJ,ne1,NPNE,nr,NVJE,XP,ZA,ERROR,*9999)
      CALL CIRCUM(NBJ,ne2,NPNE,nr,NVJE,XP,ZA,ERROR,*9999)

      CALL EXITS('EDGEFLIP')
      RETURN
 9999 CALL ERRORS('EDGEFLIP',ERROR)
      CALL EXITS('EDGEFLIP')
      RETURN 1
      END


