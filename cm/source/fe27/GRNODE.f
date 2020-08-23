      SUBROUTINE GRNODE(NBJ,NEELEM,NEL,NELIST,NPL,NPNE,NPNODE,
     '  NPLIST,NPLIST2,NRLIST,NVCNODE,NXI,CE,XP,STRING,ERROR,*)

C#### Subroutine: GRNODE
C###  Description:
C###    Group nodes under a label.

C**** NTGRNO is number of node groups currently defined.
C**** LAGRNO(nogrno) is label given to group number NOGRNO.
C**** LIGRNO(0,nogrno) is number in list for group number NOGRNO.
C**** LIGRNO(1..,nogrno) is list for group number NOGRNO.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),
     '  NELIST(0:NEM),NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPLIST(0:NPM),NPLIST2(0:NPM),
     &  NRLIST(0:NRM),NVCNODE(2,NP_R_M),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 CE(NMM,NEM,NXM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER ERR,IBEG,IBEG2,IEND,IEND2,N3CO,
     '  n1,nb,ne,ni,nl,node1,node2,nonode,nn,
     '  NN_TOT,NNlist1(0:4),NNlist2(0:4),NNlist3(0:4),
     '  noelem,no_nrlist,np,NPEXCL(0:NPM),nplag,NPLIST3(0:500),
     &  NPN(NPM,-3:3),nr,NVERSIONS,Xi1_value,Xi2_value,Xi3_value,
     &  XIDIRN(6)
      REAL*8 radius,xyz(3)
      CHARACTER CHAR2*2,STRING2*255,TYPE*8
      LOGICAL ALL_REGIONS,AS,CBBREV,CENTRE,FOUND,INLIST,
     '  STRING_SET,UNSORTED,Xi1_fixed,Xi2_fixed,Xi3_fixed,XIEND
!     Functions
      REAL*8 RFROMC

      CALL ENTERS('GRNODE',*9999)
      IF(noco.EQ.NTCO) THEN
        CO(noco+1)='?'
        NTCO=NTCO+1
        CALL STRING_TRIM(STRING,IBEG,IEND)
        STRING=STRING(IBEG:IEND)//' nodes'
      ENDIF        
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
C        CHAR2=CFROMI(NTGRNO+1,'(I2)')
        WRITE(CHAR2,'(I2)') NTGRNO+1
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM group nodes NODE#s/NODE_GROUP
C###  Parameter:     <as LABEL[node_1]>
C###    Specifies the name of the node group.
C###  Parameter:     <(sorted/unsorted)[sorted]>
C###    The sorted command sorts the node numbers in ascending order
C###    and removes duplicates.
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###    Group nodes under a label. By default the list is sorted
C###    by element numbers into ascending order and duplicates are
C###    removed - the UNSORTED option prevents this.

        OP_STRING(1)=STRING(1:IEND)//' NODE#s/NODE_GROUP'
        OP_STRING(2)=BLANK(1:15)
     '    //'<as LABEL[node_'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(3)=BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group nodes in
C###  Parameter:     <elements (#s/all)[all]>
C###    Group all nodes that are in the list of elements.
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the region for grouping instead of elements.
C###  Parameter:     <as LABEL[node_1]>
C###    Specifies the name of the node group.
C###  Description:
C###    Group nodes in a list of elements or region under a label. The
C###    nodes are sorted and duplicates removed.

        OP_STRING(1)= STRING(1:IEND)//' in element (#s/all)[all]'
        OP_STRING(2)=BLANK(1:15)
     '    //'<as LABEL[node_'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(3)=BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group nodes between NODE1#,NODE2#
C###  Parameter:     <exclude NODE#s/NODE_GROUP>
C###    Specifies the nodes to exclude from list.
C###  Parameter:     <as LABEL[node_1]>
C###    Specifies the name of the node group.
C###  Parameter:     <(sorted/unsorted)[sorted]>
C###    The sorted command sorts the node numbers in ascending order
C###    and removes duplicates.
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###    Group nodes under a label. By default the list is sorted
C###    by element numbers into ascending order and duplicates are
C###    removed - the UNSORTED option prevents this.
C###    This command should be extended to cope with four nodes that 
C###    will group nodes on a face.

        OP_STRING(1)=STRING(1:IEND)
     &    //' between NODE1#,NODE2#'
        OP_STRING(2)=BLANK(1:15)
     &    //'<as LABEL[node_'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(3)=BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group nodes
C###  Parameter:     <xi1=(0/1)>
C###    Specifies which nodes in each element in the Xi_1 direction
C###    are to be included in the group.
C###  Parameter:     <xi2=(0/1)>
C###    Specifies which nodes in each element in the Xi_2 direction
C###    are to be included in the group.
C###  Parameter:     <xi3=(0/1)>
C###    Specifies which nodes in each element in the Xi_3 direction
C###    are to be included in the group.
C###  Parameter:     <(external/all)[all]>
C###    Specifies that only outer nodes are to be included.
C###  Parameter:     <element (#s/all)[all]>
C###    Only includes the nodes which make up the specified elements.
C###    The 'all' command prompts all currently defined elements
C###    to be included.
C###  Parameter:     <(sorted/unsorted)[sorted]>
C###    The sorted command sorts the node numbers in ascending order
C###    and removes duplicates.
C###  Parameter:     <as LABEL[node_1]>
C###    Specifies the name of the node group.
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###    Creates a node group with name 'LABEL'. Restrictions can be
C###    imposed to only include the nodes at certain positions within
C###    the elements.

        OP_STRING(1)=STRING(1:IEND)//' <xi1=(0/1)'
        OP_STRING(2)=  BLANK(1:15)//'<xi2=(0/1)>'
        OP_STRING(3)=  BLANK(1:15)//'<xi3=(0/1)>'
        OP_STRING(4)=  BLANK(1:15)//'<(external/all)[all]>'
        OP_STRING(5)=  BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(6)=  BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(7)=  BLANK(1:15)
     '    //'<as LABEL[node_'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(8)=  BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C---------------------------------------------------------------------

C#### Command: FEM group nodes voronoi_internal
C###  Description:
C###    Groups all internal voronoi nodes together. Only set up to
C###    handle one region.
C###  Parameter:     <region (#)[1]>
C###    Specify the element file region number to be defined.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Parameter:     <as LABEL[node_1]>
C###    Specifies the name of the node group.

        OP_STRING(1)=STRING(1:IEND)//' VORONOI_INTERNAL'
        OP_STRING(2)=BLANK(1:15)
     '    //'<as LABEL[node_'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(3)=BLANK(1:15)//'<region (#)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C---------------------------------------------------------------------

C#### Command: FEM group nodes common
C###  Description:
C###    Groups all nodes that are common to two specified groups.
C###  Parameter:     <(#s/group name) [1]>
C###    Specify the name of the first group of nodes
C###  Parameter:     <and (#s/group name)[1]>
C###    Specify the name of the second group of nodes
C###  Parameter:     <as LABEL[node_1]>
C###    Specifies the name of the new node group.

        OP_STRING(1)=STRING(1:IEND)//' VORONOI_INTERNAL'
        OP_STRING(2)=BLANK(1:15)
     '    //'<as LABEL[node_'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(3)=BLANK(1:15)//'<region (#)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','GRNODE',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)

        STRING_SET=.FALSE. !default string set from command line
        IF(CBBREV(CO,'UNSORTED',1,noco+1,NTCO,N3CO)) THEN
          UNSORTED=.TRUE.
        ELSE
          UNSORTED=.FALSE.
        ENDIF
        
C DMAL 20 SEP 2003
C  Nodes between NODE1 and NODE2
        IF(CBBREV(CO,'BETWEEN',3,noco+1,NTCO,N3CO)) THEN
          CDATA(1)='NODES'
          CALL PARSILG(NPLIST,NPM,CDATA(1),CO(N3CO+1),ERROR,*9999)
          IF(NPLIST(0).EQ.2) THEN
            node1=NPLIST(1)
            node2=NPLIST(2)
          ELSEIF(NPLIST(0).EQ.4)THEN
            ! To be implemented for grouping nodes on a 2d face
            CALL ASSERT(.FALSE.,
     &        '>> No support for grouping nodes on faces',
     &        ERROR,*9999)
          ELSE
            CALL ASSERT(.FALSE.,
     &      '>> Incorrect number of nodes to group between',
     &      ERROR,*9999)
          ENDIF
        ENDIF !xi1

C  Nodes on Xi# boundary
        IF(CBBREV(CO,'XI1',3,noco+1,NTCO,N3CO)) THEN
          Xi1_fixed=.TRUE.
          CALL INTFROMCHAR(Xi1_value,CO(N3CO+1),ERR)
          IF(ERR.NE.0) GOTO 9998
        ELSE
          Xi1_fixed=.FALSE.
        ENDIF !xi1

        IF(CBBREV(CO,'XI2',3,noco+1,NTCO,N3CO)) THEN
          Xi2_fixed=.TRUE.
          CALL INTFROMCHAR(Xi2_value,CO(N3CO+1),ERR)
          IF(ERR.NE.0) GOTO 9998
        ELSE
          Xi2_fixed=.FALSE.
        ENDIF !xi2

        IF(CBBREV(CO,'XI3',3,noco+1,NTCO,N3CO)) THEN
          Xi3_fixed=.TRUE.
          CALL INTFROMCHAR(Xi3_value,CO(N3CO+1),ERR)
          IF(ERR.NE.0) GOTO 9998
        ELSE
          Xi3_fixed=.FALSE.
        ENDIF !xi3

        IF(Xi1_fixed.OR.Xi2_fixed.OR.Xi3_fixed) THEN
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Xi1_fixed='',L1,'
     '                     //''' Xi2_fixed='',L1,'
     '                     //''' Xi3_fixed='',L1)')
     '      Xi1_fixed,Xi2_fixed,Xi3_fixed
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF !dop
          IF(CBBREV(CO,'EXTERNAL',3,noco+1,NTCO,N3CO)) THEN
            TYPE='EXTERNAL'
          ELSE
            TYPE='ALL'
          ENDIF
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Type = '',A)') TYPE
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF !dop

          NPLIST(0)=0
          DO no_nrlist=1,NRLIST(0) !region list
            nr=NRLIST(no_nrlist)     !region#
            DO noelem=1,NELIST(0) !element list
              ne=NELIST(noelem)
              nb=NBJ(1,ne)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' ne='',I4,'' nb='',I4)') ne,nb
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF

! 1D elements
              IF(NIT(nb).EQ.1) THEN !1D elements
!             Compile lists of element nodes to include in group
                IF(Xi1_fixed) THEN
                  IF(Xi1_value.EQ.0) THEN
                    np=NPNE(1,nb,ne)
                  ELSE IF(Xi1_value.EQ.1) THEN
                    np=NPNE(2,nb,ne)
                  ENDIF
                  NPLIST(0)=NPLIST(0)+1
                  IF(NPLIST(0).LE.NPM) NPLIST(NPLIST(0))=np
                ENDIF !xi1

! 2D elements
              ELSE IF(NIT(nb).EQ.2) THEN !2D elements
!             Compile lists of element nodes to include in group
                IF(Xi1_fixed) THEN
                  IF(Xi1_value.EQ.0) THEN
                    NNlist1(1)=1 !element node 1
                    NNlist1(2)=3 !element node 3
C                    ne_adjacent=NXI(-1,1,ne)
                  ELSE IF(Xi1_value.EQ.1) THEN
                    NNlist1(1)=2 !element node 2
                    NNlist1(2)=4 !element node 4
C                    ne_adjacent=NXI( 1,1,ne)
                  ENDIF
                  NNlist1(0)=2
                ELSE
                  NNlist1(0)=0
                ENDIF !xi1

                IF(Xi2_fixed) THEN
                  IF(Xi2_value.EQ.0) THEN
                    NNlist2(1)=1 !element node 1
                    NNlist2(2)=2 !element node 2
C                    ne_adjacent=NXI(-1,1,ne)
                  ELSE IF(Xi2_value.EQ.1) THEN
                    NNlist2(1)=3 !element node 3
                    NNlist2(2)=4 !element node 4
C                    ne_adjacent=NXI( 1,1,ne)
                  ENDIF
                  NNlist2(0)=2
                ELSE
                  NNlist2(0)=0
                ENDIF !xi2

                DO nn=1,NNT(nb) !loop through all elements nodes
                  IF(   (.NOT.Xi1_fixed.
     .                   OR.INLIST(nn,NNlist1(1),NNlist1(0),n1)).
     '              AND.(.NOT.Xi2_fixed.
     '                   OR.INLIST(nn,NNlist2(1),NNlist2(0),n1)))
     '              THEN
                    np=NPNE(nn,nb,ne) !global node#
                    NPLIST(0)=NPLIST(0)+1
                    IF(NPLIST(0).LE.NPM) NPLIST(NPLIST(0))=np
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' nn='',I4,'' np='',I6)')
     '                  nn,np
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
c                   IF((TYPE(1:8).EQ.'EXTERNAL'.AND.ne_adjacent.EQ.0).OR
c    '                .TYPE(1:3).EQ.'ALL') THEN
c                   ENDIF !type
                  ENDIF !inlist
                ENDDO !nn

! 3D elements
              ELSE IF(NIT(nb).EQ.3) THEN !3D elements
!             Compile lists of element nodes to include in group
                IF(Xi1_fixed) THEN
                  IF(Xi1_value.EQ.0) THEN
                    NNlist1(1)=1 !element node 1
                    NNlist1(2)=3 !element node 3
                    NNlist1(3)=5 !element node 5
                    NNlist1(4)=7 !element node 7
C                    ne_adjacent=NXI(-1,1,ne)
                  ELSE IF(Xi1_value.EQ.1) THEN
                    NNlist1(1)=2 !element node 2
                    NNlist1(2)=4 !element node 4
                    NNlist1(3)=6 !element node 6
                    NNlist1(4)=8 !element node 8
C                    ne_adjacent=NXI( 1,1,ne)
                  ENDIF
                  NNlist1(0)=4
                ELSE
                  NNlist1(0)=0
                ENDIF !xi1
 
                IF(Xi2_fixed) THEN
                  IF(Xi2_value.EQ.0) THEN
                    NNlist2(1)=1 !element node 1
                    NNlist2(2)=2 !element node 2
                    NNlist2(3)=5 !element node 5
                    NNlist2(4)=6 !element node 6
C                    ne_adjacent=NXI(-1,1,ne)
                  ELSE IF(Xi2_value.EQ.1) THEN
                    NNlist2(1)=3 !element node 3
                    NNlist2(2)=4 !element node 4
                    NNlist2(3)=7 !element node 7
                    NNlist2(4)=8 !element node 8
C                    ne_adjacent=NXI( 1,1,ne)
                  ENDIF
                  NNlist2(0)=4
                ELSE
                  NNlist2(0)=0
                ENDIF !xi2

                IF(Xi3_fixed) THEN
                  IF(Xi3_value.EQ.0) THEN
                    NNlist3(1)=1 !element node 1
                    NNlist3(2)=2 !element node 2
                    NNlist3(3)=3 !element node 3
                    NNlist3(4)=4 !element node 4
C                    ne_adjacent=NXI(-1,1,ne)
                  ELSE IF(Xi3_value.EQ.1) THEN
                    NNlist3(1)=5 !element node 5
                    NNlist3(2)=6 !element node 6
                    NNlist3(3)=7 !element node 7
                    NNlist3(4)=8 !element node 8
C                    ne_adjacent=NXI( 1,1,ne)
                  ENDIF
                  NNlist3(0)=4
                ELSE
                  NNlist3(0)=0
                ENDIF !xi3

                DO nn=1,NNT(nb) !loop through all elements nodes
                  IF(   (.NOT.Xi1_fixed.
     .                   OR.INLIST(nn,NNlist1(1),NNlist1(0),n1)).
     '              AND.(.NOT.Xi2_fixed.
     '                   OR.INLIST(nn,NNlist2(1),NNlist2(0),n1)).
     '              AND.(.NOT.Xi3_fixed.
     '                   OR.INLIST(nn,NNlist3(1),NNlist3(0),n1)))
     '              THEN
                    np=NPNE(nn,nb,ne) !global node#
                    NPLIST(0)=NPLIST(0)+1
                    IF(NPLIST(0).LE.NPM) NPLIST(NPLIST(0))=np
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' nn='',I4,'' np='',I6)') nn,np
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDIF !inlist
                ENDDO !nn

c               nfelem=NF_elem(Xi_coord,Xi_value) !local element face#
c               nf    =NFF(nfelem,ne)             !global face#
c               IF(DOP) WRITE(*,'('' nfelem='',I4,'' nf='',I4)')
c    '            nfelem,nf
c               IF((TYPE(1:8).EQ.'EXTERNAL'.AND.NPF(5,nf).EQ.1).OR
c    '             .TYPE(1:3).EQ.'ALL') THEN
c                 DO nne=1,NNF(0,nfelem,nb) !loop over face nodes
c                   nn=NNF(1+nne,nfelem,nb) !element node
c                   np=NPNE(nn,nb,ne)    !global node
c                   IF(DOP) WRITE(*,'('' nn='',I4,'' np='',I4)') nn,np
c                   IF(.NOT.INLIST(np,NPLIST(1),MIN(NPLIST(0),NPM),n1))
c    '                THEN
c                     NPLIST(0)=NPLIST(0)+1
c                     IF(NPLIST(0).LE.NPM) NPLIST(NPLIST(0))=np
c                   ENDIF
c                 ENDDO !nne
c               ENDIF !type

              ENDIF !nit
            ENDDO !noelem
          ENDDO !no_nrlist
          IF(NPLIST(0).GT.NPM) THEN
            CALL WRITE_CHAR(IOER,'Increase NPM to at least ',ERR)
            CALL WRITE_INT(IOER,NPLIST(0),ERR)
            CALL WRITE_CHAR(IOER,NEWLINE,ERR)
          ENDIF
          CALL ASSERT(NPLIST(0).LE.NPM,'>>Increase NPM',ERROR,*9999)

C  All nodes between NODE1 and NODE2
        ELSE IF(CBBREV(CO,'BETWEEN',3,noco+1,NTCO,N3CO)) THEN
          CALL NPLINK(NPN,NPL,ERROR,*9999)
          FOUND=.FALSE.
          XIDIRN(1)=1
          XIDIRN(2)=-1
          XIDIRN(3)=2
          XIDIRN(4)=-2
          XIDIRN(5)=3
          XIDIRN(6)=-3
          DO ni=1,6
            IF(.NOT.FOUND)THEN
              XIEND=.FALSE.
              NPLIST(0)=1
              NPLIST(1)=node1
            ENDIF
            DO WHILE(.NOT.XIEND)
              IF(.NOT.FOUND)THEN
                IF(NPN(NPLIST(NPLIST(0)),XIDIRN(ni)).NE.0)THEN
                  IF(NPN(NPLIST(NPLIST(0)),XIDIRN(ni)).EQ.node2)THEN
                    FOUND=.TRUE.
                    XIEND=.TRUE.
                  ELSE IF(NPN(NPLIST(NPLIST(0)),
     &                XIDIRN(ni)).EQ.node1)THEN !circular
                    FOUND=.FALSE.
                    XIEND=.TRUE.
                  ENDIF
                  NPLIST(NPLIST(0)+1)=NPN(NPLIST(NPLIST(0)),XIDIRN(ni))
                  NPLIST(0)=NPLIST(0)+1
                ELSE
                  XIEND=.TRUE.
                ENDIF
              ENDIF
            ENDDO
          ENDDO
          
          IF(CBBREV(CO,'EXCLUDE',3,noco+1,NTCO,N3CO)) THEN
            CDATA(1)='NODES'
            CALL PARSILG(NPEXCL,NPM,CDATA(1),CO(N3CO+1),ERROR,*9999)
            nplag=0
            DO np=1,NPLIST(0)
              IF(.NOT.INLIST(NPLIST(np),NPEXCL(1),NPEXCL(0),n1))THEN
                nplag=nplag+1
                NPLIST(nplag)=NPLIST(np)
              ENDIF
            ENDDO
            NPLIST(0)=nplag
          ENDIF

C  All external node points
        ELSE IF(CBBREV(CO,'EXTERNAL',3,noco+1,NTCO,N3CO)) THEN
          nr=1 !temporary?
          IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.2) THEN !2D geometry
            NPLIST(0)=0
            DO nl=1,NLT
              IF(NEL(0,nl).EQ.1) THEN !line nl is on boundary
                IF(NPL(1,1,nl).EQ.1) THEN      !linear Lagrange basis
                  NN_TOT=2
                ELSE IF(NPL(1,1,nl).EQ.2) THEN !quadratic Lagrange
                  NN_TOT=3
                ELSE IF(NPL(1,1,nl).EQ.3) THEN !cubic Lagrange
                  NN_TOT=4
                ELSE IF(NPL(1,1,nl).EQ.4) THEN !cubic Hermite
                  NN_TOT=2
                ENDIF
                DO nn=1,NN_TOT
                  IF(DOP) THEN
                    WRITE(*,'('' nl='',I4,'' nn='',I2,'' np='',I5)')
     '                nl,nn,NPL(1+nn,1,nl)
                  ENDIF
                  IF(.NOT.INLIST(NPL(1+nn,1,nl),
     '              NPLIST(1),MIN(NPLIST(0),NPM),n1)) THEN
                    NPLIST(0)=NPLIST(0)+1
                    IF(NPLIST(0).LE.NPM) NPLIST(NPLIST(0))=
     '                NPL(1+nn,1,nl)
                  ENDIF
                ENDDO !nn
              ENDIF !NEL
            ENDDO !nl
            IF(NPLIST(0).GT.NPM) THEN
              CALL WRITE_CHAR(IOER,'Increase NPM to at least ',ERR)
              CALL WRITE_INT(IOER,NPLIST(0),ERR)
              CALL WRITE_CHAR(IOER,NEWLINE,ERR)
            ENDIF
            CALL ASSERT(NPLIST(0).LE.NPM,'>>Increase NPM',ERROR,*9999)

          ELSE
! Need to use faces to check for bdry nodes
            CALL ASSERT(.true.,' >>Not implemented for 1D or 3D',
     '        ERROR,*9999)
          ENDIF !2D
        ELSE IF(CBBREV(CO,'IN',2,noco+1,NTCO,N3CO)) THEN
          !group all nodes in an element or element group
          UNSORTED=.FALSE. !nodes sorted and duplicates removed
          NPLIST(0)=0
          DO noelem=1,NELIST(0) !already called PARSE_ELEMENTS
            ne=NELIST(noelem)
            nb=NBJ(1,ne)
            DO nn=1,NNT(nb)
              NPLIST(0)=NPLIST(0)+1
              IF(NPLIST(0).GT.NPM) THEN
                CALL WRITE_CHAR(IOER,'Increase NPM to at least ',ERR)
                CALL WRITE_INT(IOER,NPLIST(0),ERR)
                CALL WRITE_CHAR(IOER,NEWLINE,ERR)
              ENDIF
              CALL ASSERT(NPLIST(0).LE.NPM,' >>Increase NPM',
     '          ERROR,*9999)
              NPLIST(NPLIST(0))=NPNE(nn,nb,ne)
            ENDDO
          ENDDO
        ELSE IF(CBBREV(CO,'COMMON',2,noco+1,NTCO,N3CO)) THEN
          !group all nodes that are common between two existing groups
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          CALL CUPPER(CO(N3CO+1)(IBEG:IEND),STRING)
          !put first node group into NPLIST
          CALL PARSILG(NPLIST,NPM,'NODES',STRING,ERROR,*9999)
          IF(CBBREV(CO,'AND',3,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
            CALL CUPPER(CO(N3CO+1)(IBEG:IEND),STRING)
            !put second node group into NPLIST2
            CALL PARSILG(NPLIST2,NPM,'NODES',STRING,ERROR,*9999)
          ELSE
            NPLIST2(0)=0
          ENDIF
          NPLIST3(0)=0
          DO nonode=1,NPLIST2(0)
            np=NPLIST2(nonode)
            IF(INLIST(np,NPLIST(1),NPLIST(0),n1))THEN
              NPLIST3(0)=NPLIST3(0)+1
              NPLIST3(NPLIST3(0))=np
              CALL ASSERT(NPLIST3(0).LE.500,
     &          '>>Increase NPLIST3 in GRNODE',ERROR,*9999)
            ENDIF
          ENDDO
          NPLIST(0)=NPLIST3(0)
          DO nonode=1,NPLIST3(0)
            NPLIST(nonode)=NPLIST3(nonode)
          ENDDO
          
C TVK 10/04/2000 Voronoi Internal Nodes Group
        ELSE IF(CBBREV(CO,'VORONOI_INTERNAL',3,noco+1,NTCO,N3CO)) THEN
          nr=NRLIST(1) !Only one region possible
          NPLIST(0)=0
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            IF(NVCNODE(1,nonode).EQ.3) THEN
              NPLIST(0)=NPLIST(0)+1
              NPLIST(NPLIST(0))=np
            ENDIF
          ENDDO
        ELSE IF(CBBREV(CO,'TERMINAL',3,noco+1,NTCO,N3CO)) THEN
c          STRING_SET=.TRUE. !default string set in GRNODE_BY
c          STRING2='TERMINAL' !DEFAULT NAME

          nr=NRLIST(1)
          CENTRE=.FALSE.
          IF(CBBREV(CO,'CENTRE',3,noco+1,NTCO,N3CO)) THEN
            CENTRE=.TRUE.
            xyz(1)=RFROMC(CO(N3CO+1))
            xyz(2)=RFROMC(CO(N3CO+2))
            xyz(3)=RFROMC(CO(N3CO+3))
            
            IF(CBBREV(CO,'RADIUS',3,noco+1,NTCO,N3CO)) THEN
              radius=RFROMC(CO(N3CO+1))
            ELSE
              radius=1.d0
            ENDIF
          ENDIF
          CALL GRNODE_BY(NBJ,NEELEM,NPLIST,NPNE,nr,NXI,CE(1,1,1),
     '      radius,XP,xyz,STRING2,CENTRE,ERROR,*9999)
        ELSE
          CALL PARSE_NODES(NPNODE,NPLIST,noco-1,NRLIST,NTCO,CO,
     '      ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'EXCLUDE',3,noco+1,NTCO,N3CO)) THEN
          CDATA(1)='NODES'
          CALL PARSILG(NPEXCL,NPM,CDATA(1),CO(N3CO+1),ERROR,*9999)
          nplag=0
          DO np=1,NPLIST(0)
            IF(.NOT.INLIST(NPLIST(np),NPEXCL(1),NPEXCL(0),n1))THEN
              nplag=nplag+1
              NPLIST(nplag)=NPLIST(np)
            ENDIF
          ENDDO
          NPLIST(0)=nplag
        ENDIF

        IF(.NOT.UNSORTED) THEN
C         Sort and remove duplicates from the node list
          CALL ILISTRMDUP(NPLIST(0),NPLIST(1),ERROR,*9999)
        ENDIF

C        CALL ASSERT(NPLIST(0).LE.GRNO_MAXDATA,
C     '    '>>Increase array sizes in grou00.cmn',ERROR,*9999)

        IF(CBBREV(CO,'AS',1,noco+1,NTCO,N3CO)) THEN
          AS=.TRUE.
        ELSE
          AS=.FALSE.
        ENDIF
        IF(.NOT.STRING_SET)THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          CALL CUPPER(CO(N3CO+1)(IBEG:IEND),STRING2)
        ELSE
          AS=.TRUE.
        ENDIF
        CALL GRNODE_SUB(NPLIST,STRING2,AS,ERROR,*9999)
      ENDIF

      CALL EXITS('GRNODE')
      RETURN
 9998 ERROR=' '
 9999 CALL ERRORS('GRNODE',ERROR)
      CALL EXITS('GRNODE')
      RETURN 1
      END


