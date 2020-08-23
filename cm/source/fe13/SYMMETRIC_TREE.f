      SUBROUTINE SYMMETRIC_TREE(NBJ,NEELEM,NELIST,NENP,NKJ,NKJE,
     &  NP_INTERFACE,NPNE,NPNODE,NORD,nr,num_gen,NVJE,NVJP,NXI,
     &  branch_angle,length_ratio,BRANCH_LENGTH,SE,XP,FRACTAL,ERROR,*)


C#### Subroutine: SYMMETRIC_TREE
C###  Description:
C###  SYMMETRIC_TREE generates a symmetric branching tree 1D mesh. Each
C###  new branch is created at a plane rotation angle of 90 degrees to
C###  the parent branch and with a constant input branching angle
C###  (branch_angle). 


      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      
!     Parameter list
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     &  NENP(NPM,0:NEPM,0:NRM),NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     &  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     &  NPNODE(0:NP_R_M,0:NRM),NORD(5,NE_R_M),nr,NRE(NEM),num_gen,
     &  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 branch_angle,BRANCH_LENGTH(*),length_ratio,
     &  SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL FRACTAL
!     Local variables
      INTEGER gen,i,j,nb,ne,ne0,NELIST2(0:NEM),
     &  ne_parent,ne_parent2,nj,no_branch,noelem,noelem_parent,nonode,
     &  np,np0,np_parent,np_parent2,TMP(0:NEM)
      REAL*8 A(NJT),ANGLE,B(NJT),FACTOR(2),length,length_new,
     &  parent_length,PLANE_ANGLE,N(NJT),NEW(NJT),NEW2(njt),Rx(3,3),
     &  Rv(3,4),V(NJT),R(NJT)

      CALL ENTERS('SYMMETRIC_TREE',*9999)

      FACTOR(1)=1.d0
      FACTOR(2)=-1.d0
      PLANE_ANGLE=1.5708d0
      branch_angle=branch_angle/180.d0*PI !converting degrees to radians
      nb=1 !TEMPORARY
      ne=NET(0)
      noelem=NEELEM(0,nr)
      np=NPT(0)
      nonode=NPNODE(0,nr)
C.. if start element not defined - create start node at 0,0,0
      IF(NELIST(0).EQ.0) THEN
        nonode=nonode+1
        np=np+1
        CALL ASSERT(nonode+1.LE.NP_R_M,'>>Increase NP_R_M',
     '    ERROR,*9999)
        CALL ASSERT(np+1.LE.NPM,'>>Increase NPM',ERROR,*9999)
        NPNODE(nonode,nr)=np
        DO nj=1,NJT
          XP(1,1,nj,np)=0.d0
        ENDDO
        nonode=nonode+1
        np=np+1
        NPNODE(nonode,nr)=np
        XP(1,1,1,np)=0.d0
        XP(1,1,2,np)=0.d0
        length=BRANCH_LENGTH(1) !1st generation branch length
        XP(1,1,3,np)=-length
        ne=ne+1
        noelem=noelem+1
        CALL ASSERT(noelem.LE.NE_R_M,'>>Increase NE_R_M',
     '    ERROR,*9999)
        CALL ASSERT(ne.LE.NEM,'>>Increase NEM',ERROR,*9999)
        NEELEM(noelem,nr)=ne
        NPNE(1,nb,ne)=np-1
        NPNE(2,nb,ne)=np
        CALL GN1DNEJ(nb,NBJ,ne,NKJ,NKJE,NPNE(2,nb,ne),
     '    NPNE(1,nb,ne),nr,NRE,NVJE,NVJP,SE,ERROR,*9999)
        NELIST(0)=1
        NELIST(1)=ne
        NORD(1,ne)=1
      ELSE !check tree has been ordered
        IF(NORD(1,NELIST(1)).EQ.0) THEN
          WRITE(OP_STRING,
     &      '(''*** Make sure orders have been evaluated ***'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)          
        ENDIF
      ENDIF
      DO noelem_parent=1,NELIST(0) !for each pre-defined parent
        ne_parent=NELIST(noelem_parent)
        np_parent=NPNE(2,nb,ne_parent)
        IF(NXI(1,0,ne_parent).EQ.0) THEN !only if terminal node
          gen=NORD(1,ne_parent) !initialise
          NELIST2(0)=1
          NELIST2(1)=ne_parent
          parent_length=0.d0
          IF(FRACTAL) THEN
            DO nj=1,3
              parent_length=parent_length+(XP(1,1,nj,np_parent)-
     &          XP(1,1,nj,NPNE(1,nb,ne_parent)))**2.d0
            ENDDO
            parent_length=DSQRT(parent_length)
          ENDIF !fractal
          DO i=1,num_gen !number of branch generations to create
            TMP(0)=0
            gen=gen+1
            IF(FRACTAL) THEN !If fractal use length_ratio instead of morphometric data
              length=parent_length*length_ratio**i
            ELSE
              length=BRANCH_LENGTH(gen)
            ENDIF
C.. create two new nodes and element
            DO j=1,NELIST2(0) !for each terminal
              ne_parent=NELIST2(j)
              np_parent=NPNE(2,nb,ne_parent) !terminal node
              np0=NPNE(1,nb,ne_parent)
              ne0=NXI(-1,1,ne_parent) !grandparent
              IF(NXI(1,0,ne0).GE.2) THEN
                ne_parent2=NXI(1,2,ne0) !sister branch to ne_parent
                IF(ne_parent2.EQ.ne_parent) ne_parent2=NXI(1,1,ne0)
                np_parent2=NPNE(2,nb,ne_parent2)
                DO nj=1,NJT
                  A(nj)=XP(1,1,nj,np_parent)-XP(1,1,nj,np0)
                  B(nj)=XP(1,1,nj,np_parent2)-XP(1,1,nj,np0)
                ENDDO !nj
                CALL CROSS(A,B,N) !normal vector to 2 parent vectors
                CALL NORMALISE(NJT,N,ERROR,*9999)
              ELSE
                N(1)=1.d0
                N(2)=0.d0 !x-axis
                N(3)=0.d0
              ENDIF
              DO no_branch=1,2
                nonode=nonode+1
                np=np+1
                CALL ASSERT(nonode.LE.NP_R_M,'>>Increase NP_R_M',
     '            ERROR,*9999)
                CALL ASSERT(np.LE.NPM,'>>Increase NPM',ERROR,*9999)
                NPNODE(nonode,nr)=np
                DO nj=1,NJT
                  V(nj)=XP(1,1,nj,np_parent)-XP(1,1,nj,np0)
                  R(nj)=V(nj)
                ENDDO
                !Rotation about the x-axis by branching angle
                ANGLE=FACTOR(no_branch)*branch_angle
                CALL ROTATE_V(ANGLE,N,Rx,ERROR,*9999)
                CALL RESET(Rv) !sets up identity matrix
                IF(MOD(i,2).EQ.0) THEN
                  CALL NORMALISE(NJT,R,ERROR,*9999)
                  CALL ROTATE_V(PLANE_ANGLE,R,Rv,ERROR,*9999)
                ENDIF
                length_new=0.d0
                DO nj=1,NJT
                  NEW(nj)=Rx(nj,1)*V(1)+Rx(nj,2)*V(2)+Rx(nj,3)*V(3)
                ENDDO
                DO nj=1,NJT                  
                  NEW2(nj)=Rv(nj,1)*NEW(1)+Rv(nj,2)*NEW(2)+Rv(nj,3)
     &              *NEW(3)
                ENDDO
                DO nj=1,NJT
                  NEW(nj)=NEW2(nj)
                  length_new=length_new+NEW(nj)**2.d0
                ENDDO
                length_new=DSQRT(length_new)
                DO nj=1,NJT
                  XP(1,1,nj,np)=XP(1,1,nj,np_parent)+NEW(nj)*
     &              (length/length_new)
                ENDDO
                ne=ne+1
                noelem=noelem+1
                CALL ASSERT(noelem.LE.NE_R_M,'>>Increase NE_R_M',
     '            ERROR,*9999)
                CALL ASSERT(ne.LE.NEM,'>>Increase NEM',ERROR,*9999)
                NEELEM(noelem,nr)=ne
                NXI(1,0,ne_parent)=no_branch
                CALL ASSERT(no_branch.LE.NEIM,'>>Increase NEIM',ERROR,
     &            *9999)
                NXI(1,no_branch,ne_parent)=ne
                NXI(-1,0,ne)=1
                NXI(-1,1,ne)=ne_parent
                IF(i.LE.num_gen) THEN
                  TMP(0)=TMP(0)+1
                  TMP(TMP(0))=ne !new parent/terminal node
                ENDIF
                NPNE(1,nb,ne)=np_parent
                NPNE(2,nb,ne)=np
                CALL GN1DNEJ(nb,NBJ,ne,NKJ,NKJE,NPNE(2,nb,ne),
     '            NPNE(1,nb,ne),nr,NRE,NVJE,NVJP,SE,ERROR,*9999)
              ENDDO !no_branch
            ENDDO !j
            NELIST2(0)=0
            DO j=1,TMP(0)
              NELIST2(0)=NELIST2(0)+1 
              NELIST2(NELIST2(0))=TMP(j) !create new terminal list
            ENDDO
          ENDDO !i
        ENDIF !TERMINAL
      ENDDO !noelem_parent
      NEELEM(0,nr)=noelem
      NPNODE(0,nr)=nonode
      NPNODE(0,0)=np
      NEELEM(0,0)=ne
      NPT(nr)=nonode
      NET(nr)=noelem
      NPT(0)=np
      NET(0)=ne
      CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
      CALL CALC_NENP_1D(nb,NEELEM,NENP,NPNE,nr,ERROR,*9999)
      CALL NENXI_1D(nb,NEELEM,NENP,NPNE,nr,NXI,ERROR,*9999)


      CALL EXITS('SYMMETRIC_TREE')
      RETURN
 9999 CALL ERRORS('SYMMETRIC_TREE',ERROR)
      CALL EXITS('SYMMETRIC_TREE')
      RETURN 1
      END


