      SUBROUTINE CALC_VORO_ALVEOLI(nb,NBJ,NEELEM,NENFVC,NENP,NFVC,NKJ,
     '  NKJE,NMAX_GENERATION,NODENVC,NP_INTERFACE,NPNE,NPNODE,nr,NRE,
     '  nr_host,nr_target,NVCNODE,NVJE,NVJP,SE,XP,ZA,ERROR,*)

C#### Subroutine: CALC_VORO_ALVEOLI
C###  Description:
C###    CALC_VORO_ALVEOLI converts Voronoi vertices and face
C###    boundaries to nodes and elements.  These are written into the
C###    same region that the original Delaunay elements were in - the
C###    Delaunay nodes and elements are discarded.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER nb,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENFVC(0:NFVCM,NFVM),NENP(NPM,0:NEPM,0:NRM),
     '  NFVC(2,0:NFVCM,NVCM),NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     '  NMAX_GENERATION,NODENVC(NVCM),NP_INTERFACE(0:NPM,0:3),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),nr_host,
     '  nr_target,NVCNODE(2,NP_R_M),NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM)
      REAL*8 SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER ERR,IBEG,IBEG1,IEND,IEND1,NADJACENT_DUCT(0:3,NVCM),
     '  NCELL_ADJACENT(0:1000),ne,ne_alv,NELEM_LIST1(0:NVCM*5),
     '  NELEM_LIST2(0:NVCM*5),NELEM_LIST3(0:NVCM),nfvl,N_GENERATION,nj,
     '  nk,nn,noelem,noelem_alv,nonode,nonvc,nonvc_adjacent,nonvc_temp,
     &  no_nvc,np,npnvc,npnvc_adjacent,ns,num_alveoli,nvc,nvc_adjacent,
     '  NVC_CLASSIFY(NVCM),NVC_DUCT(1000),NVC_DUCT_TEMP(1000),nvc_max,
     '  nvc_max2,nvc_min,nvc_min2,NVC_TERMINAL(0:1000),
     '  NVC_TERMINAL_TEMP(1000),NODE_CONVERT_NVC(NE_R_M)
      REAL*8 ANGLE,DIST,MAX_ANGLE,MIN_ANGLE,MIN_DIST,U(3),
     '  U_ADJACENT(3),U_TERMINAL(3,1000),U_TERMINAL_TEMP(3,1000),
     '  TOL,X_NVC(3),X_NVC_ADJACENT(3),XP_HOST(3)
      LOGICAL BIFURCATING,CONTINUE
      CHARACTER CHAR*1,CHAR1*5,STRING*(255)
!     Functions
      REAL*8 ANG2V
!     External functions
      INTEGER IDIGITS

C NCELL_ADJACENT stores list of adjacent cells that haven't yet been
C     classified
C     NVC_DUCT stores the ith cell number from start of duct. e.g. =4 if
C     the cell is 4th from start.
C     NFVC(1,nfvl,nvc) gives local Delaunay node number (in region nr)
C     for face nfvl of cell nvc
C     NVCNODE(2,nonode) gives cell number for local (in nr) node nonode
      
      CALL ENTERS('CALC_VORO_ALVEOLI',*9999)
      
      CALL OPENF(IOFILE2,'DISK','alveolar_groups.com','NEW',
     '  'SEQUEN','FORMATTED',132,ERROR,*9999)
      
      TOL=1.d-8
C     Starting cell and direction.  Make the start cell be the first
C     cell with an internal node at its centre.
      DO nonvc=1,NVCM
        NVC_CLASSIFY(nonvc)=0
        NADJACENT_DUCT(0,nonvc)=0
      ENDDO
      DO noelem=1,NE_R_M
        NODE_CONVERT_NVC(noelem)=0
      ENDDO
      N_GENERATION=0
      num_alveoli=0
C     Get global coordinates of host node
      DO nj=1,NJT
        U_TERMINAL(nj,1)=DSQRT(1.d0/3.d0)
      ENDDO!nj
      np=NPNODE(1,nr_host)
      DO nj=1,NJT
        XP_HOST(nj)=XP(1,1,nj,np)+0.25d0*U_TERMINAL(nj,1)
      ENDDO !nj

      MIN_DIST=1.d6
      DO nvc=1,NVCT !for each Voronoi cell
        nonode=NODENVC(nvc)
        np=NPNODE(nonode,nr)
        IF(NVCNODE(TYPE,nonode).EQ.INTERNAL)THEN
          DIST=0.d0
          DO nj=1,NJT
            DIST=DIST+(XP(1,1,nj,np)-XP_HOST(nj))**2.d0
          ENDDO !nj
          DIST=DSQRT(DIST)
          IF(DIST.LT.MIN_DIST+TOL)THEN
            MIN_DIST=DIST
            nvc_min=nvc
          ENDIF
        ENDIF
      ENDDO !nvc
      NVC_TERMINAL(0)=1
      NVC_TERMINAL(1)=nvc_min
        

      NVC_CLASSIFY(NVC_TERMINAL(1))=1 !duct
      NVC_DUCT(1)=1
      
C     Do this for a specified number of generations, else until all
C     possible adjacent cells have been classified into duct or
C     alveolus.
      CONTINUE=.TRUE.
      DO WHILE(CONTINUE)
        
        IF(NVC_DUCT(1).EQ.1)THEN
          N_GENERATION=N_GENERATION+1
        ENDIF
C       For each 'terminal' cell (i.e. each cell from last iteration)
C       find adjacent cells, and select next one (duct) or two
C       (bifurcation).
        nonvc_temp=0
        DO nonvc=1,NVC_TERMINAL(0)
C         Set counter for 'terminal' cells to zero
          nvc=NVC_TERMINAL(nonvc) !current cell
          no_nvc=NODENVC(nvc)
          npnvc=NPNODE(no_nvc,nr) !node at centre of current cell
          DO nj=1,NJT
            X_NVC(nj)=XP(1,1,nj,npnvc)
            U(nj)=U_TERMINAL(nj,nonvc) !get the unit direction vector
          ENDDO !nj
          
C         Make a list of all potential adjacent cells.  i.e. all
C         cells that have not been labelled as alveoli or ducts.
          NCELL_ADJACENT(0)=0
 !CHECK use of nonode, np etc.
          DO nfvl=1,NFVC(1,0,nvc) !for each face in current cell
            nonode=NFVC(1,nfvl,nvc) !local node (nr) centre of adjacent cell
            IF(NVCNODE(TYPE,nonode).EQ.INTERNAL)THEN
              nvc_adjacent=NVCNODE(2,nonode) !cell number of adjacent cell
              IF(NVC_CLASSIFY(nvc_adjacent).EQ.0)THEN !hasn't been classified
                NCELL_ADJACENT(0)=NCELL_ADJACENT(0)+1
                NCELL_ADJACENT(NCELL_ADJACENT(0))=nvc_adjacent
              ENDIF !CELL_CLASSIFY
            ENDIF !is it internal?
          ENDDO !nfvl
          IF(NCELL_ADJACENT(0).EQ.0) CONTINUE=.FALSE. !no more cells to classify
C         Check to see whether the duct should bifurcate.  Will do so if
C         the parent cell is the 4th in the duct.
          IF(CONTINUE) THEN
            IF(N_GENERATION.GT.1.AND.NVC_DUCT(nonvc).EQ.3)THEN
              BIFURCATING=.TRUE.
            ELSE IF(N_GENERATION.EQ.1.AND.NVC_DUCT(nonvc).EQ.4)THEN
              BIFURCATING=.TRUE.
            ELSE
              BIFURCATING=.FALSE.
            ENDIF !NCELL_DUCT
            
            IF(BIFURCATING)THEN !will divide, needs two daughter cells
              
C             Find the cell that has the greatest angle (less than
C             90deg) to the parent cell direction.  This will become one of the 
C             'daughter' cells.  The other daughter cell will be the one
C             that has the greatest angle to the first daughter.
              MAX_ANGLE=0.d0 !initialise maximum angle
              DO nonvc_adjacent=1,NCELL_ADJACENT(0) !for each adjacent cell
                nvc_adjacent=NCELL_ADJACENT(nonvc_adjacent) !cell number
                no_nvc=NODENVC(nvc_adjacent) !node at cell centre
                npnvc_adjacent=NPNODE(no_nvc,nr) !node at cell centre
                DO nj=1,NJT !get coordinates of cell centre
                  X_NVC_ADJACENT(nj)=XP(1,1,nj,npnvc_adjacent)
                  U_ADJACENT(nj)=X_NVC_ADJACENT(nj)-X_NVC(nj) !calc direction vector
                ENDDO !nj
                CALL NORMALISE(NJT,U_ADJACENT,ERROR,*9999) !normalize U_ADJACENT
                
C               calculate angle between previous direction and direction
C               to this cell; cos^(-1) of dot product of the two
C               direction
C               vectors.
                ANGLE=ANG2V(U,U_ADJACENT)
 !!! POSSIBLE REASON FOR DIFFERENT RESULTS ON DIFF PLATFORMS.IF(ANGLE.LT.PI/2.d0)
                
               IF(ANGLE.LT.PI/2.d0+TOL)THEN !candidate for daughter cell
                  IF(ANGLE.GT.MAX_ANGLE+TOL)THEN !new maximum angle
                    MAX_ANGLE=ANGLE
                    nvc_max=nvc_adjacent !store the cell number for minimum angle
                    DO nj=1,NJT
                      U_TERMINAL_TEMP(nj,nonvc_temp+1)=U_ADJACENT(nj)
                    ENDDO !nj
                  ENDIF !ANGLE
                ENDIF !ANGLE
                NVC_CLASSIFY(nvc_adjacent)=2 !initially classify as an alveolus
              ENDDO !nonvc_adjacent
              
C             Classify the nvc_max cell as a duct, and store for next
C             iteration
              NVC_CLASSIFY(nvc_max)=1 !duct
              nonvc_temp=nonvc_temp+1 !increment the number of 'terminal' cells
              NVC_TERMINAL_TEMP(nonvc_temp)=nvc_max
              NVC_DUCT_TEMP(nonvc_temp)=1 !start of a new duct, so re-set to 1
              NADJACENT_DUCT(0,nvc)=NADJACENT_DUCT(0,nvc)+1
              NADJACENT_DUCT(NADJACENT_DUCT(0,nvc),nvc)=nvc_max
              
C             Find the other daughter cell.
              MAX_ANGLE=0.d0 !initialise maximum angle
              MIN_ANGLE=2.d0*PI !initialise minimum angle
              nvc_max2=0
              DO nonvc_adjacent=1,NCELL_ADJACENT(0)
                nvc_adjacent=NCELL_ADJACENT(nonvc_adjacent) !cell number
                IF(nvc_adjacent.NE.nvc_max)THEN !continue if not the same cell
                  no_nvc=NODENVC(nvc_adjacent) !node at cell centre
                  npnvc_adjacent=NPNODE(no_nvc,nr)
                  DO nj=1,NJT !get coordinates of cell centre
                    X_NVC_ADJACENT(nj)=XP(1,1,nj,npnvc_adjacent)
                    U_ADJACENT(nj)=X_NVC_ADJACENT(nj)-X_NVC(nj) !calc direction vector
                  ENDDO !nj
                  CALL NORMALISE(NJT,U_ADJACENT,ERROR,*9999) !normalize U_ADJACENT
C                 Calculate angle with parent cell.
                  ANGLE=ANG2V(U,U_ADJACENT)
                  IF(ANGLE.LT.PI/2.d0+TOL)THEN !candidate for second daughter cell
C                   Calculate angle between daughter cells.
                    ANGLE=ANG2V(U_TERMINAL_TEMP(1,nonvc_temp),
     '                U_ADJACENT)
                    IF(ANGLE.GT.MAX_ANGLE+TOL)THEN !new maximum angle
                      MAX_ANGLE=ANGLE
                      nvc_max2=nvc_adjacent !store the cell # for minimum angle
                      DO nj=1,NJT
                        U_TERMINAL_TEMP(nj,nonvc_temp+1)=U_ADJACENT(nj)
                      ENDDO !nj
                    ENDIF !ANGLE.GT.MAX_ANGLE
                  ELSE !just in-case none of them are less than pi/2
                    ANGLE=ANG2V(U_TERMINAL_TEMP(1,nonvc),U_ADJACENT)
                    IF(ANGLE.LT.MIN_ANGLE+TOL)THEN !new minimum angle
                      MIN_ANGLE=ANGLE
                      nvc_min2=nvc_adjacent !store the cell # for minimum angle
                      DO nj=1,NJT
                        U_TERMINAL_TEMP(nj,nonvc_temp+1)=U_ADJACENT(nj)
                      ENDDO !nj
                    ENDIF !ANGLE.GT.MAX_ANGLE
                    
                  ENDIF !ANGLE.LT.PI/2.d0
                ENDIF !nvc_adjacent
              ENDDO !nonvc_adjacent
C             Classify the nvc_max2 cell as a duct, and store for next
C             iteration
              IF(nvc_max2.NE.0)THEN
                NVC_CLASSIFY(nvc_max2)=1 !duct
                nonvc_temp=nonvc_temp+1 !increment the number of 'terminal' cells
                NVC_TERMINAL_TEMP(nonvc_temp)=nvc_max2
                NVC_DUCT_TEMP(nonvc_temp)=1 !start of a new duct, so re-set to 1
                NADJACENT_DUCT(0,nvc)=NADJACENT_DUCT(0,nvc)+1
                NADJACENT_DUCT(NADJACENT_DUCT(0,nvc),nvc)=nvc_max2
              ELSE !was not a cell with angle l.t. pi/2, so use smallest
                NVC_CLASSIFY(nvc_min2)=1 !duct
                nonvc_temp=nonvc_temp+1 !increment the number of 'terminal' cells
                NVC_TERMINAL_TEMP(nonvc_temp)=nvc_min2
                NVC_DUCT_TEMP(nonvc_temp)=1 !start of a new duct, so re-set to 1
                NADJACENT_DUCT(0,nvc)=NADJACENT_DUCT(0,nvc)+1
                NADJACENT_DUCT(NADJACENT_DUCT(0,nvc),nvc)=nvc_min2
              ENDIF !nvc_max2
            ELSE IF(.NOT.BIFURCATING)THEN !continue as straight as possible
              
C             Determine next duct cell.  Will be the cell with direction
C             vector (calculated as direction from cell-to-cell centres)
C             closest to previous direction vector.
              MIN_ANGLE=2.d0*PI !initialise minimum angle
              DO nonvc_adjacent=1,NCELL_ADJACENT(0) !for each adjacent cell
                nvc_adjacent=NCELL_ADJACENT(nonvc_adjacent) !cell number
                no_nvc=NODENVC(nvc_adjacent) !node at cell centre
                npnvc_adjacent=NPNODE(no_nvc,nr)
                DO nj=1,NJT !get coordinates of cell centre
                  X_NVC_ADJACENT(nj)=XP(1,1,nj,npnvc_adjacent)
                  U_ADJACENT(nj)=X_NVC_ADJACENT(nj)-X_NVC(nj) !calc direction vector
                ENDDO !nj
                CALL NORMALISE(NJT,U_ADJACENT,ERROR,*9999) !normalize U_ADJACENT
                
C               calculate angle between previous direction and direction
C               to this cell; cos^(-1) of dot product of the two
C               direction
C               vectors.
                ANGLE=ANG2V(U,U_ADJACENT)
                IF(ANGLE.LT.MIN_ANGLE+TOL)THEN !new minimum angle
                  MIN_ANGLE=ANGLE
                  nvc_min=nvc_adjacent !store the cell number for minimum angle
                  DO nj=1,NJT
                    U_TERMINAL_TEMP(nj,nonvc_temp+1)=U_ADJACENT(nj)
                  ENDDO !nj
                ENDIF !ANGLE
                IF(N_GENERATION.EQ.1.AND.NVC_DUCT(1).EQ.1)THEN
 !don't classify surrounding alveoli
                ELSE
                  NVC_CLASSIFY(nvc_adjacent)=2 !initially classify as an alveolus
                ENDIF !N_GENERATION
              ENDDO !nonvc_adjacent
              
C             Classify the nvc_min cell as a duct, and store for next
C             iteration
              NVC_CLASSIFY(nvc_min)=1 !duct
              nonvc_temp=nonvc_temp+1 !increment the number of 'terminal' cells
              NVC_TERMINAL_TEMP(nonvc_temp)=nvc_min
              NVC_DUCT_TEMP(nonvc_temp)=NVC_DUCT(nonvc)+1
              NADJACENT_DUCT(0,nvc)=NADJACENT_DUCT(0,nvc)+1
              NADJACENT_DUCT(NADJACENT_DUCT(0,nvc),nvc)=nvc_min
              
            ENDIF !BIFURCATING
          ENDIF !CONTINUE
        ENDDO !nonvc
        IF(DOP)THEN
          WRITE(OP_STRING,'('' Cell iteration:'',I6)') N_GENERATION
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        
C       Put 'temp' information into 'current' information, so that the
C       next iteration can be done.
        NVC_TERMINAL(0)=nonvc_temp !the new number of 'terminal' cells
        DO nonvc=1,NVC_TERMINAL(0)
          NVC_TERMINAL(nonvc)=NVC_TERMINAL_TEMP(nonvc)
          NVC_DUCT(nonvc)=NVC_DUCT_TEMP(nonvc)
          DO nj=1,NJT
            U_TERMINAL(nj,nonvc)=U_TERMINAL_TEMP(nj,nonvc)
          ENDDO !nj
        ENDDO !nonvc
        
C       Test to see whether the algorithm should continue.  Check the
C       number of generations (i.e. number of bifurcations) and whether
C       there are any 'terminal' cells to check. CONTINUE will be set to
C       false if the number of specified generations is exceeded, or
C       there are no more adjacent cells to classify.
        IF(NVC_TERMINAL(0).EQ.0)THEN
          CONTINUE=.FALSE.
        ELSEIF(N_GENERATION.EQ.NMAX_GENERATION.AND.NVC_DUCT(1).EQ.3)THEN
          CONTINUE=.FALSE.
        ELSE
          CONTINUE=.TRUE.
        ENDIF
      ENDDO !CONTINUE

C     'cap off' any terminal cells with surrounding alveoli.

      IF(DOP)THEN
        WRITE(OP_STRING,'('' Cap end cells:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      DO nonvc=1,NVC_TERMINAL(0)
        nvc=NVC_TERMINAL(nonvc) !current cell
        DO nfvl=1,NFVC(1,0,nvc) !for each face in current cell
          nonode=NFVC(1,nfvl,nvc) !local node (nr) centre of adjacent cell
          IF(NVCNODE(TYPE,nonode).EQ.INTERNAL)THEN
            nvc_adjacent=NVCNODE(2,nonode) !cell number of adjacent cell
            IF(NVC_CLASSIFY(nvc_adjacent).EQ.0)THEN !hasn't been classified
              NVC_CLASSIFY(nvc_adjacent)=2 !alveolus
            ENDIF !CELL_CLASSIFY
          ENDIF
        ENDDO !nfvl
      ENDDO !nonvc

      
C     Convert the Voronoi cells into elements.  The only cells that will
C     be converted are the 'alveoli'.  Duct cells will be discarded, and
C     alveolar cells will have surface elements created over each face.
C     Faces that adjoin duct cells will not have elements created over
C     them.
      
      IF(DOP)THEN
        WRITE(OP_STRING,'('' Convert cells to elements'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      nonode=0
      np=NPNODE(0,0)
      noelem=0
      ne=NEELEM(0,0)

      NELEM_LIST1(0)=0 !will use this for list of ducts
      NELEM_LIST2(0)=0 !will use this for list of alveoli
      
      DO nvc=1,NVCT
        IF(NVC_CLASSIFY(nvc).EQ.1)THEN
C         Call CALC_VORO_ELEM to create nodes at each face centre and
C         vertex for the duct cells. Put into a list for creating
C         element group 'duct'.
          NELEM_LIST3(0)=0
          CALL CALC_VORO_ELEM(nb,NBJ,ne,NEELEM,NELEM_LIST1,NELEM_LIST3,
     '      NENFVC,NENP,NFVC,NODE_CONVERT_NVC,noelem,nonode,np,NPNE,
     '      NPNODE,nr_target,nvc,NVC_CLASSIFY,NVCNODE,XP,ZA,ERROR,*9999)
c          CALL CALC_VORO_ELEM_1D(NADJACENT_DUCT(0,nvc),NBJ,ne,NEELEM,
c     '      NELEM_LIST1,NENP,nodenvc,noelem,nonode,np,NPNE,NPNODE,nr,
c     '      nr_target,nvc,XP,ERROR,*9999)
        ELSE IF(NVC_CLASSIFY(nvc).EQ.2)THEN !alveolus, make elements
C         Call CALC_VORO_ELEM to create nodes at each face centre and
C         vertex for the alveolar cells. Put into a list for creating
C         element group 'alveoli'.
          NELEM_LIST3(0)=0
          CALL CALC_VORO_ELEM(nb,NBJ,ne,NEELEM,NELEM_LIST2,NELEM_LIST3,
     '      NENFVC,NENP,NFVC,NODE_CONVERT_NVC,noelem,nonode,np,NPNE,
     '      NPNODE,nr_target,nvc,NVC_CLASSIFY,NVCNODE,XP,ZA,ERROR,*9999)
          IF(NELEM_LIST3(0).NE.0) THEN
          !groups elems in each alveolar "cell"
            num_alveoli=num_alveoli+1
            STRING='alveolus'
            WRITE(CHAR1,'(I3)') num_alveoli
            CALL STRING_TRIM(STRING,IBEG,IEND)
            CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
            CALL APPENDC(IEND,CHAR1(IBEG1:IEND1),STRING)
            CALL GRELEM_SUB(NELEM_LIST3,STRING,.TRUE.,ERROR,*9999)
C.. Writing out groups to file so they can be read back in later when
C.. creating capillary mesh:
            CALL WRITE_CHAR(IOFILE2,'fem group elem ',ERR)
            DO noelem_alv=1,NELEM_LIST3(0)
              ne_alv=NELEM_LIST3(noelem_alv)
              WRITE(CHAR,'($,I1)') IDIGITS(ne_alv)
              CALL WRITE_INT(IOFILE2,ne_alv,ERR)
              IF(noelem_alv.LT.NELEM_LIST3(0)) THEN
                CALL WRITE_CHAR(IOFILE2,',',ERR)
              ENDIF            
            ENDDO
            WRITE(OP_STRING,'($'' as '',A46)') STRING
            CALL WRITES(IOFILE2,OP_STRING,ERROR,*9999)            
          ENDIF
        ENDIF !NVC_CLASSIFY(nvc).EQ.2
        NVC_CLASSIFY(nvc_adjacent)=3 !so that we know this cell has been done
      ENDDO !nvc
      
C     Create element groups for the duct cells and alveolar cells.
      IF(NELEM_LIST1(0).NE.0) THEN
        STRING='duct'
        CALL GRELEM_SUB(NELEM_LIST1,STRING,.TRUE.,ERROR,*9999)
      ENDIF
      IF(NELEM_LIST2(0).NE.0) THEN
        STRING='alveoli'
        CALL GRELEM_SUB(NELEM_LIST2,STRING,.TRUE.,ERROR,*9999)
      ENDIF
      
      NPNODE(0,0)=NPNODE(0,0)+nonode
      NPNODE(0,nr_target)=nonode
      NPT(nr_target)=np
      NPT(0)=NPNODE(0,0)
      DO nonode=1,NPNODE(0,nr_target) !for each new node
        np=NPNODE(nonode,nr_target)
        DO nj=1,NJT
          NKJ(nj,np)=1
          NVJP(nj,np)=1
        ENDDO !nj
      ENDDO !nonode(np)
      NEELEM(0,0)=NEELEM(0,0)+noelem
      NEELEM(0,nr_target)=noelem
      NET(nr_target)=ne
      NET(0)=NEELEM(0,0)
      DO noelem=1,NEELEM(0,nr_target) !for each new elem
        ne=NEELEM(noelem,nr_target)
        NRE(ne)=nr_target
c        DO nj=1,NJ_LOC(NJL_GEOM,0,nr_target)
c          NBJ(nj,ne)=nb
c        ENDDO !nj
        DO ns=1,NST(nb)+NAT(nb)
          SE(ns,nb,ne)=1.0d0 !1D elements with unit scale factors
        ENDDO !ns
        DO nn=1,NNT(nb)
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr_target)
            NVJE(nn,nb,nj,ne)=1
            DO nk=1,NKT(nn,nb)
              NKJE(nk,nn,nj,ne)=1
            ENDDO !nk
          ENDDO !nj
        ENDDO !nn
      ENDDO !noelem(ne)

C     Remove the Delaunay region, becuase carrying around the triangles
C     is just TOO big.
      NPNODE(0,0)=NPNODE(0,0)-NPNODE(0,nr)
      NEELEM(0,0)=NEELEM(0,0)-NEELEM(0,nr)
      NPT(0)=NPNODE(0,0)
      NET(0)=NEELEM(0,0)
      NPNODE(0,nr)=0
      NEELEM(0,nr)=0
      
      CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
c      CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)

      CALL CLOSEF(IOFILE2,ERROR,*9999)

      CALL EXITS('CALC_VORO_ALVEOLI')
      RETURN
 9999 CALL ERRORS('CALC_VORO_ALVEOLI',ERROR)
      CALL EXITS('CALC_VORO_ALVEOLI')
      RETURN 1
      END


