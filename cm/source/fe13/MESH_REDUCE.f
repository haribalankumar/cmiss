      SUBROUTINE MESH_REDUCE(NBJ,NEELEM,NENP,NNB,NORD,NP_INTERFACE,NPNE,
     '  NPNODE,nr,NXI,TERMINAL_ORDER,ORDER,ERROR,*)


C#### Subroutine: MESH_REDUCE
C###  Description:
C###    MESH_REDUCE reduces a 1D tree to a specified generation or order.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'mesh00.cmn'
      INCLUDE 'pulm00.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NNB(4,4,4,NBFM),NORD(5,NE_R_M),
     '  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  TERMINAL_ORDER
      LOGICAL ORDER
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ne,ne0,ne2,NEELEM_INVERSE(NE_R_M),noelem,noelem0,
     '  nonode,norder,norder_parent,norder_sibling,np,np0,
     &  NPNODE_INVERSE(NP_R_M),N_ELEM,N_NODE

      CALL ENTERS('MESH_REDUCE',*9999)

      nb=NBJ(1,NEELEM(1,nr))
      N_ELEM=0
      N_NODE=0
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        NEELEM_INVERSE(ne)=noelem
      ENDDO!noelem
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        NPNODE_INVERSE(np)=nonode
      ENDDO!nonode

      IF(.NOT.ORDER)THEN
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(ne.NE.0)THEN
            norder=NORD(1,ne) !Weibel generation
            IF(norder.GT.TERMINAL_ORDER)THEN
              noelem0=NEELEM_INVERSE(ne)
              NEELEM(noelem0,nr)=0 !deleting ne0
              N_ELEM=N_ELEM+1
              nonode=NPNODE_INVERSE(NPNE(2,nb,ne)) !end node
              N_NODE=N_NODE+1
              NPNODE(nonode,nr)=0 !delete end node
            ENDIF !norder
          ENDIF !ne
        ENDDO !noelem
      ELSE IF(ORDER)THEN
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          ne0=NXI(-1,1,ne) !parent element
          IF(ne.NE.0.AND.ne0.NE.0)THEN
            norder=NORD(2,ne) !Horsfield order
            norder_parent=NORD(2,ne0)
            ! get the sibling branch
            ne2=NXI(1,1,ne0)
            IF(ne2.EQ.ne) ne2=NXI(2,1,ne0)
            norder_sibling=NORD(2,ne2)
            
            IF(norder.LE.TERMINAL_ORDER)THEN
              IF(norder_sibling.GE.TERMINAL_ORDER)THEN
                !keep the element
              ELSEIF((norder_parent-norder.EQ.1).OR.
     &          (norder_parent-norder.EQ.0).OR.
     '          ((norder_parent-norder.GT.1).AND.
     '          (norder_parent.LE.TERMINAL_ORDER)))THEN
 !if order of parent more than 1 higher than branch, then sister branch is higher than norder. Keep the branch
                noelem0=NEELEM_INVERSE(ne)
                NEELEM(noelem0,nr)=0 !deleting ne0
                N_ELEM=N_ELEM+1
                nonode=NPNODE_INVERSE(NPNE(2,nb,ne)) !end node
                N_NODE=N_NODE+1
                NPNODE(nonode,nr)=0 !delete end node
              ENDIF
            ENDIF !norder
          ENDIF !ne
        ENDDO !noelem
      ENDIF
      ne0=0
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        IF(ne.NE.0)THEN
          ne0=ne0+1
          NEELEM(ne0,nr)=ne
        ENDIF !ne
      ENDDO !noelem
      np0=0
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        IF(np.NE.0)THEN
          np0=np0+1
          NPNODE(np0,nr)=np
        ENDIF !np
      ENDDO !nonode
      NEELEM(0,nr)=NEELEM(0,nr)-N_ELEM
      NPNODE(0,nr)=NPNODE(0,nr)-N_NODE
      NPNODE(0,0)=NPNODE(0,0)-N_NODE !number nodes in all regions
      NEELEM(0,0)=NEELEM(0,0)-N_ELEM !number elements in all regions
      CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
      CALL CALC_NENP_1D(nb,NEELEM,NENP,NPNE,nr,ERROR,*9999)
      CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)
      LTYP_R(4)=4 !defining each volume independently

      CALL EXITS('MESH_REDUCE')
      RETURN
 9999 CALL ERRORS('MESH_REDUCE',ERROR)
      CALL EXITS('MESH_REDUCE')
      RETURN 1
      END



