      SUBROUTINE CALC_ELEM_SHAR_VERT(NBJ,ne,NEAM,NEES,
     '  NENP,NEVS,NNB,NPNE,NVERTNE,ERROR,*)

C#### Subroutine: CALC_ELEM_SHAR_VERT
C###  Description:
C###    CALC_ELEM_SHAR_VERT calculates the vertex and element
C###    sharing arrays NEVS and NEES for element
C###    ne, used by CREATE_LATTICE for generating the lattice
C###    grid scheme.

C#### Variable: NEES(nn,ne_adj)
C###  Type: INTEGER
C###  Set_up: CALC_ELEM_SHAR_VERT
C###  Description:
C###    NEES(nn,ne_adj) stores the elements that share vertices
C###    with an element. NEES(nn,1) stores vertices of this element.      
      
C#### Variable: NEVS(12,ne_adj)      
C###  Type: INTEGER
C###  Set_up: CALC_ELEM_SHAR_VERT
C###  Description:
C###    NEVS(12,ne_adj) is a temporary array that stores how a
C###    given element shares its
C###    vertices with adjacent elements. NEVS(1,ne_adj)
C###    stores the global element number for the ne_adj
C###    described by that column. NEVS(2...9,ne_adj) stores the
C###    vertices of ne_adj that are shared with the vertices of an
C###    element.
C###    NEVS(10..12,ne_adj) stores the xi 1,2 and 3 direction mappings
C###    of ne_adj to the element of interest. These mappings store
C###    the xi directions of element ne_adj that are shared with
C###    with the element of interest. For instance if a line is
C###    shared between the element of interest and ne_adj, and this
C###    line is in the xi 1 direction of the element, and the xi 2
C###    direction of ne_adj, NEVS(10,ne_adj)=2. If the line were
C###    in the xi 3 direction of the element, NEVS(12,ne_adj)=2.      
C###    NEVS(:,1) stores the information about the current element.      
      
      IMPLICIT NONE
      
      INCLUDE 'geom00.cmn' 
!     Parameter List
      INTEGER NBJ(NJM,NEM),ne,NEAM,NEES(8,NEAM),NENP(NPM,0:NEPM),
     '  NEVS(0:12,NEAM),NNB(4,4,4,NBFM),NPNE(NNM,NBFM,NEM),
     '  NVERTNE(0:8,2)
      CHARACTER ERROR*(*)
      
!     Local Variables
      INTEGER adjacent_element,adjacent_node,eindex,ETOIJK(12),
     '  ETOOE(3,12),EVLIST(2),FTOE(4,6),FTOOF(6),i,lce,lcf,
     '  lcv,maxadjelem,nae,nae_adj,nb,nbestunique,ne_adj,ni,nf,nn,
     '  NVERTNE_ADJ(0:8,2),none_adj,np,nunique,product,VLIST(2),
     '  VTOE(2,12),VTOF(4,6),VTOOV(7,8)
      LOGICAL NOTFOUND

      DATA ETOIJK /1,1,1,1,2,2,2,2,3,3,3,3/   
      DATA ETOOE /2,3,4,
     '            1,3,4,
     '            1,2,4,
     '            1,2,3,
     '            6,7,8,
     '            5,7,8,
     '            5,6,8,
     '            5,6,7,
     '            10,11,12,
     '            9,11,12,
     '            9,10,12,
     '            9,10,11/
      DATA FTOE /5,6,9,11,
     '           7,8,10,12,
     '           1,3,9,10,
     '           2,4,11,12,
     '           1,2,5,7,
     '           3,4,6,8/
      DATA FTOOF /2,1,4,3,6,5/
      DATA VTOE /1,2,
     '           3,4,
     '           5,6,
     '           7,8,
     '           1,3,
     '           5,7,
     '           2,4,
     '           6,8,
     '           1,5,
     '           2,6,
     '           3,7,
     '           4,8/ 
      DATA VTOF /1,3,5,7,
     '           2,4,6,8,   
     '           1,5,2,6, 
     '           3,7,4,8,
     '           1,2,3,4, 
     '           5,6,7,8/   
      DATA VTOOV /2,3,4,5,6,7,8,
     '            1,3,4,5,6,7,8,
     '            1,2,4,5,6,7,8,
     '            1,2,3,5,6,7,8,
     '            1,2,3,4,6,7,8,
     '            1,2,3,4,5,7,8,
     '            1,2,3,4,5,6,8,
     '            1,2,3,4,5,6,7/ 

      CALL ENTERS('CALC_ELEM_SHAR_VERT',*9999)

!     Initialise arrays and variables

      NOTFOUND=.TRUE.
      maxadjelem=0
      
      DO nn=1,8
        DO ne_adj=1,NEAM
          NEES(nn,ne_adj)=0
        ENDDO
      ENDDO
      DO i=1,12
        DO ne_adj=1,NEAM
          NEVS(i,ne_adj)=0
        ENDDO
      ENDDO
      
C     Set up first column of NEVS and NEES
      NEVS(1,1)=ne
      DO nn=1,8
        NEVS(nn+1,1)=NVERTNE(nn,1)
        NEES(nn,1)=ne
      ENDDO
      DO ni=9,11
        NEVS(ni+1,1)=ni-8
      ENDDO

C     First store the shared vertices in NEVS
      
      DO nn=1,8 !loop over the eight vertices
        np=NVERTNE(nn,2) !get the global node that matches nn
C       Find the number of elements adjacent to vertice nn
        none_adj=NENP(np,0) 
        DO ne_adj=1,none_adj !loop over the adjacent elements
C         Get the global element number for adjacent element
          adjacent_element=NENP(np,ne_adj)
          eindex=2 !eindex keeps track of where we are in NEVS
          NOTFOUND=.TRUE.
C         NENP(np,0) will give all the elements adjacent to global
C         node np. Check we are not in ne.
          IF(adjacent_element.EQ.ne) THEN !in ne
            eindex=1
          ELSE !valid adjacent element
            DO WHILE(NOTFOUND)
C             Look for a zero entry in NEVS(1,:) and store
C             the global number of the adjacent element there.
              IF(NEVS(1,eindex).EQ.0) THEN
                NEVS(1,eindex)=adjacent_element
                NOTFOUND=.FALSE.
              ELSEIF(NEVS(1,eindex).NE.adjacent_element) THEN
C               NEVS(1,eindex) was not zero, so increment eindex.
                eindex=eindex+1
              ELSE
                NOTFOUND=.FALSE.           
              ENDIF                
            ENDDO
            adjacent_node=1
            nb=NBJ(1,adjacent_element)
            DO WHILE(NPNE(adjacent_node,nb,adjacent_element).NE.np)
              adjacent_node=adjacent_node+1
            ENDDO
            NEVS(nn+1,eindex)=adjacent_node
            NEES(nn,eindex)=adjacent_element
          ENDIF
          IF(eindex.GT.maxadjelem) THEN !increase maxadjelem
            maxadjelem=eindex-1
          ENDIF        
        ENDDO !adjacent element loop
      ENDDO !vertices loop

C     Look at the collapsing of the current and adjacent elements 
      
      DO ne_adj=1,maxadjelem !loop over adjacent elements
C       Sort shared vertices into faces. lcf is least common
C       face.
        lcf=0 
        nbestunique=0
        DO nf=1,6 !loop over the six faces
C         This product is zero if the four vertices for the
C         current face are not shared with the adjacent element.
          product=NEVS(VTOF(1,nf)+1,ne_adj+1)*
     '      NEVS(VTOF(2,nf)+1,ne_adj+1)*NEVS(VTOF(3,nf)+1,ne_adj+1)
     '      *NEVS(VTOF(4,nf)+1,ne_adj+1)
          IF(product.NE.0) THEN !find the number of unique vertices
            nunique=1
C           Look at each vertex, comparing it to the other four
C           and counting the number which are unique.
            IF(NEVS(VTOF(1,nf)+1,ne_adj+1).NE.
     '        NEVS(VTOF(2,nf)+1,ne_adj+1)) nunique=nunique+1
            IF(NEVS(VTOF(3,nf)+1,ne_adj+1).NE.
     '        NEVS(VTOF(1,nf)+1,ne_adj+1).AND.
     '        NEVS(VTOF(3,nf)+1,ne_adj+1).NE.
     '        NEVS(VTOF(2,nf)+1,ne_adj+1)) nunique=nunique+1
            IF(NEVS(VTOF(4,nf)+1,ne_adj+1).NE.
     '        NEVS(VTOF(3,nf)+1,ne_adj+1).AND.
     '        NEVS(VTOF(4,nf)+1,ne_adj+1).NE.
     '        NEVS(VTOF(2,nf)+1,ne_adj+1).AND.
     '        NEVS(VTOF(4,nf)+1,ne_adj+1).NE.
     '        NEVS(VTOF(1,nf)+1,ne_adj+1)) nunique=nunique+1           
            IF(nunique.GT.nbestunique) THEN
              lcf=nf
              nbestunique=nunique
            ENDIF
          ENDIF
        ENDDO !end face loop
        
        IF(lcf.EQ.0) THEN !no face found, search for least common edges
C         lce is lease common edge
          lce=0
          nbestunique=0
          DO nae=1,12 !loop over edges
C           product is zero if the vertices of nae for the adjacent
C           element are not shared with ne.
            product=NEVS(VTOE(1,nae)+1,ne_adj+1)*
     '        NEVS(VTOE(2,nae)+1,ne_adj+1)
            IF(product.NE.0) THEN
              nunique=0
              IF(NEVS(VTOE(2,nae)+1,ne_adj+1).NE.
     '          NEVS(VTOE(1,nae)+1,ne_adj+1)) nunique=1
              IF(nunique.GT.nbestunique) THEN
                lce=nae
                nbestunique=nunique
              ENDIF
            ENDIF
          ENDDO !end edge loop
          IF(lce.EQ.0) THEN !work with vertices
C           lcv is least common vertex
            lcv=0
            DO nn=1,8
C             If this is true, then this the shared vertex.
              IF(NEVS(nn+1,ne_adj+1).NE.0) lcv=nn
            ENDDO
            DO nn=1,7
              IF(NEVS(VTOOV(nn,lcv)+1,ne_adj+1).EQ.0) THEN
                NEVS(VTOOV(nn,lcv)+1,ne_adj+1)=NEVS(lcv+1,ne_adj+1)
              ENDIF
            ENDDO
          ELSE !work with edges
            DO nae=1,3
              DO nn=1,2
                IF(NEVS(VTOE(nn,ETOOE(nae,lce))+1,ne_adj+1).EQ.0) THEN
                  NEVS(VTOE(nn,ETOOE(nae,lce))+1,ne_adj+1)=
     '              NEVS(VTOE(nn,lce)+1,ne_adj+1)
                ENDIF
              ENDDO
            ENDDO          
          ENDIF
        ELSE !work with faces
          DO nn=1,4
            IF(NEVS(VTOF(nn,FTOOF(lcf))+1,ne_adj+1).EQ.0) THEN
              NEVS(VTOF(nn,FTOOF(lcf))+1,ne_adj+1)=
     '          NEVS(VTOF(nn,lcf)+1,ne_adj+1)
            ENDIF
          ENDDO
        ENDIF
C       Now find the mappings between the xi directions in ne
C       and the xi directions in the adjacent element.

C       GET_ELEMENT_VERT is called to retrieve the vertices
C       for the adjacent elements. These are returned in
C       NVERTNE_ADJ.
        
        CALL GET_ELEMENT_VERT(NBJ(1,NEVS(1,ne_adj+1)),NNB,
     '    NPNE(1,1,NEVS(1,ne_adj+1)),NVERTNE_ADJ,
     '    ERROR,*9999)

        IF(lcf.NE.0) THEN
          DO nae=1,3,2 !only need to look at edges 1 and 3
            nae_adj=0
            NOTFOUND=.TRUE.
            DO WHILE(NOTFOUND.AND.nae_adj.LT.12)
              nae_adj=nae_adj+1
              IF(NEVS(VTOE(1,FTOE(nae,lcf))+1,ne_adj+1).EQ.
     '          NVERTNE_ADJ(VTOE(1,nae_adj),1).AND.
     '          NEVS(VTOE(2,FTOE(nae,lcf))+1,ne_adj+1).EQ.
     '          NVERTNE_ADJ(VTOE(2,nae_adj),1).AND.
     '          NVERTNE_ADJ(VTOE(1,nae_adj),1).NE.
     '          NVERTNE_ADJ(VTOE(2,nae_adj),1)) THEN
                NEVS(9+ETOIJK(FTOE(nae,lcf)),ne_adj+1)=ETOIJK(nae_adj)
                NOTFOUND=.FALSE.
              ELSEIF(NEVS(VTOE(1,FTOE(nae,lcf))+1,ne_adj+1).EQ.
     '            NVERTNE_ADJ(VTOE(2,nae_adj),1).AND.
     '            NEVS(VTOE(2,FTOE(nae,lcf))+1,ne_adj+1).EQ.
     '            NVERTNE_ADJ(VTOE(1,nae_adj),1).AND.
     '            NVERTNE_ADJ(VTOE(1,nae_adj),1).NE.
     '            NVERTNE_ADJ(VTOE(2,nae_adj),1)) THEN
                NEVS(9+ETOIJK(FTOE(nae,lcf)),ne_adj+1)=-ETOIJK(nae_adj)
                NOTFOUND=.FALSE.                
              ENDIF
            ENDDO
          ENDDO
        ELSEIF(lce.NE.0) THEN
          
C         Sort the vertices of the lce. 
          IF(NEVS(VTOE(1,lce)+1,ne_adj+1).LT.
     '      NEVS(VTOE(2,lce)+1,ne_adj+1)) THEN
            VLIST(1)=NEVS(VTOE(1,lce)+1,ne_adj+1)
            VLIST(2)=NEVS(VTOE(2,lce)+1,ne_adj+1)
          ELSE
            VLIST(1)=NEVS(VTOE(2,lce)+1,ne_adj+1)
            VLIST(2)=NEVS(VTOE(1,lce)+1,ne_adj+1)            
          ENDIF
          DO nae=1,12
C           Also need to sort the vertices of nae.
            IF(NVERTNE_ADJ(VTOE(1,nae),1).LT.
     '        NVERTNE_ADJ(VTOE(2,nae),1)) THEN
              EVLIST(1)=NVERTNE_ADJ(VTOE(1,nae),1)
              EVLIST(2)=NVERTNE_ADJ(VTOE(2,nae),1)
            ELSE
              EVLIST(1)=NVERTNE_ADJ(VTOE(2,nae),1)
              EVLIST(2)=NVERTNE_ADJ(VTOE(1,nae),1)            
            ENDIF
            IF(EVLIST(1).EQ.VLIST(1).AND.EVLIST(2).EQ.VLIST(2))
C           IF nae has the same vertices as lce.
     '        nae_adj=nae
          ENDDO
C         Now determine whether nae_adj is in a positive or
C         negative xi orienation to lce.
          IF(NEVS(VTOE(1,lce)+1,ne_adj+1).EQ.
     '      NVERTNE_ADJ(VTOE(1,nae_adj),1).AND.
     '      NEVS(VTOE(2,lce)+1,ne_adj+1).EQ.
     '      NVERTNE_ADJ(VTOE(2,nae_adj),1)
     '      .AND.NVERTNE_ADJ(VTOE(1,nae_adj),1).NE.
     '      NVERTNE_ADJ(VTOE(2,nae_adj),1)) THEN
            NEVS(9+ETOIJK(nae_adj),ne_adj+1)=ETOIJK(nae_adj)
          ELSEIF(NEVS(VTOE(1,lce)+1,ne_adj+1).EQ.
     '        NVERTNE_ADJ(VTOE(2,nae_adj),1).AND.
     '        NEVS(VTOE(2,lce)+1,ne_adj+1)
     '        .EQ.NVERTNE_ADJ(VTOE(1,nae_adj),1).AND.
     '        NVERTNE_ADJ(VTOE(1,nae_adj),1).NE.
     '        NVERTNE_ADJ(VTOE(2,nae_adj),1)) THEN
            NEVS(9+ETOIJK(nae_adj),ne_adj+1)=-ETOIJK(nae_adj)
          ENDIF
        ENDIF
      ENDDO !end loop over adjacent elements

      NEVS(0,1)=maxadjelem+1 !store the number of adjacent elements

      DO ne_adj=1,maxadjelem !write remainder of NEES
        DO nn=1,8 !loop over local nodes
          IF(NEES(nn,ne_adj+1).EQ.0) THEN
C           Flag that this elem node is not shared by setting its
C           entry in NEES to NEM+1.
            NEES(nn,ne_adj+1)=NEM+1 
          ENDIF
        ENDDO
      ENDDO
            
      CALL EXITS('CALC_ELEM_SHAR_VERT')
      RETURN
 9999 CALL ERRORS('CALC_ELEM_SHAR_VERT',ERROR)
      CALL EXITS('CALC_ELEM_SHAR_VERT')
      RETURN 1      
      END


