      SUBROUTINE GNBDMESH(IBT,IDO,INP,LD,LD_NP,LDTMP1,LDTMP2,nb,NBJ,
     &  NDLIST,NEELEM,NELIST,NELIST2,NENP,NE_OLD,NEP,NE_REACTIVATE,
     &  NE_START,NE_TEMP,NHOST,NKJ,NKJE,NNB,NORD,NPF,NP_INTERFACE,
     &  NPNE,NPNODE,NP_START,nr_host,nr,NRE,NSTEM,NVJE,
     &  NVJP,NXI,POINT_LIMIT,SCALE_DIST_LIMIT,angle_major,angle_minor,
     &  fraction,length_limit,min_length,RATIO_LENGTHS,ROTATION_LIMIT,
     &  SE,WD,XA,XE,XIP,XP,ZD,CALC_XI,REDUCING_DISTANCE_LIMIT,
     &  RESTRICT_BRANCH_ROTATION,RESTRICT_BRANCH_PLANE,ERROR,*)

C####  Subroutine: GNBDMESH
C###   Description:
C###     GNBDMESH uses a bifurcating distributive algorithm to fill a host
C###     volume with a tree-like mesh. Details of the algorithm can be
C###     found in Tawhai et al. Ann. Biomed. Eng. 28(7):793-802 (2000). 

C#### Variable: LD_NP(NDM)
C###  Type: INTEGER
C###  Set_up: GNBDMESH
C###  Description:
C###  LD_NP(NDM) stores the single mesh node that is associated with a
C###  unique data point at the end of a volume-filling tree generation.
C###  This is used when coupling parenchymal soft tissue mechanics to
C###  airflow in the tree mesh (see UPVERTEX for application).

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter values
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     &  LD(NDM),LD_NP(NDM),LDTMP1(NDM),LDTMP2(NDM),nb,NBJ(NJM,NEM),
     &  NDLIST(0:NDM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     &  NELIST2(0:NEM),NENP(NPM,0:NEPM,0:NRM),NE_OLD(NE_R_M),NEP(NPM),
     &  NE_REACTIVATE(0:NEM),NE_START,NE_TEMP(NE_R_M),NHOST(NP_R_M),
     &  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),NNB(4,4,4,NBFM),
     &  NORD(5,NE_R_M),NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),
     &  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NP_START,nr_host,nr,
     &  NRE(NEM),NSTEM(NE_R_M,2),NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM),POINT_LIMIT,SCALE_DIST_LIMIT
      REAL*8 angle_major,angle_minor,fraction,length_limit,min_length,
     &  RATIO_LENGTHS,ROTATION_LIMIT,SE(NSM,NBFM,NEM),WD(NJM,NDM),
     &  XA(NAM,NJM,NEM),XE(NSM,NJM),XIP(NIM,NPM),XP(NKM,NVM,NJM,NPM),
     &  ZD(NJM,NDM)
      LOGICAL CALC_XI,REDUCING_DISTANCE_LIMIT,
     &  RESTRICT_BRANCH_ROTATION,RESTRICT_BRANCH_PLANE
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER COUNT,i,itemp,LD_NUM,M,N,ncount,nd,nd_closest,nde1,nde2,j,
     &  nd_min,ne,ne0,N_ELM,N_ELM_TEMP,ne_min,ne1,ne2,neb,ne_parent,
     &  nes,ngen,nj,nk,nn,noelem,noelem_gen,noelem_initial,
     *  noelem_parent,nonode,nonode_start,np,np1,np2,np3,nld,
     &  np4,Nremoved,NT_BNS,numtb,numzero,offset,numtbsum,nclos1,
     &  nclos2,ldnum2,countLD
      REAL*8 angle1,angle2,COFM(3),distance,distance2,length,
     &  LENGTH_PARENT,DIST,DISTANCE_LIMIT,max_dist,MIN_DIST,
     &  VECTOR_LENGTH,XP1(3),rot_angle
      LOGICAL AIRWAY_BRANCHING,BRANCH,INLIST,INTERNAL,IN_USE,SS(2)
      CHARACTER FMT*100

      REAL*8 angle_uw,angle_uv,angle_vw,NRML(3),SCALAR,u(3),
     &  v(3),w(3),A(3,3),B(3),length_w,length_u,ANGLE

      CALL ENTERS('GNBDMESH',*9999)

      AIRWAY_BRANCHING=.false.

      IF(NDT.EQ.0)THEN !JUST GET PARENT XI COORDS
        nr=1
        nr_host=2
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          np1=np
          INTERNAL=.FALSE.
          DO nj=1,3
            XP1(nj)=XP(1,1,nj,np)
          ENDDO
          CALL CHECK_POINT_INTERNAL(IBT,IDO,INP,NBJ,NDLIST,NEELEM,NEP,
     &      NHOST,NKJE,np,np1,NPF,NPNE,nr_host,NVJE,NXI,SE,XA,
     &      XE,XIP,XP,XP1,ZD,INTERNAL,.TRUE.,ERROR,*9999)
        ENDDO
      ENDIF
      CALL ASSERT(NDT.GT.0,'>>Define data as seed points',ERROR,*9999)
      
      numtbsum=0
C...  Calculate the generations for the initial branches
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        IF(NXI(-1,0,ne).EQ.0)THEN !stem
          NORD(1,ne)=1
        ELSE
          ne0=NXI(-1,1,ne)
          NORD(1,ne)=NORD(1,ne0)+1
        ENDIF
        nb=NBJ(1,ne)
        np=NPNE(2,nb,ne)
        NHOST(np)=0
      ENDDO !noelem

c      IF(CALC_XI)THEN
C.......Adapt this so that a larger set of nodes can be calculated, but
C.......not the whole tree        
c        DO nonode=1,NPNODE(0,nr)
c          np=NPNODE(nonode,nr)
c          np1=np
c          INTERNAL=.FALSE.
c          DO nj=1,3
c            XP1(nj)=XP(1,1,nj,np)
c          ENDDO
c          CALL CHECK_POINT_INTERNAL(IBT,IDO,INP,NBJ,NDLIST,NEELEM,NEP,
c     &      NHOST,NKJE,np,np1,NPF,NPNE,nr_host,NVJE,NXI,SE,XA,
c     &      XE,XIP,XP,XP1,ZD,INTERNAL,.TRUE.,ERROR,*9999)
c        ENDDO
c      ENDIF

C     Initialise NE_OLD (stores parent elements for current generation)
C     to the list of parent elements in NELIST2, and NT_BNS (current
C     number of terminal branches) to the number of parent branches.
      DO N=1,NELIST2(0) !for each parent element in list
        NE_OLD(N)=NELIST2(N)
        NSTEM(NELIST2(N),1)=NELIST2(N)
        NSTEM(NELIST2(N),2)=0
      ENDDO !N
      NT_BNS=NELIST2(0) !initial number of 'terminal' branches
      N_ELM=NT_BNS
      
C     Initialise LD array to map each seed point to a parent branch.
C     For a single parent, all seed points will initially be mapped to
C     it; for multiple parents MESH_REPOINT is called to calculate
C     the closest parent end-point to each seed point.
      DO nd=1,NDT !for every data point
c        LD_NP(nd)=LD(nd) !temp storage of LD 
c        LD(nd)=NELIST2(1) !assigns data points to the single parent
        LDTMP2(nd)=NELIST2(1)
        DO nj=1,NJT !set weighting factor all=1.0
          WD(nj,NDT)=1.0d0
        ENDDO !nj
      ENDDO !nd

      IF(NELIST2(0).GT.1)THEN
        CALL ASSERT(NORD(1,NELIST2(1)).GT.0,
     &    '>> Use FEM EVALUATE ORDERING first ',ERROR,*9999)
        CALL MESH_REPOINT(LDTMP2,nb,N_ELM,NE_OLD,NE_REACTIVATE,
     '    NPNE,500.d0,XP,ZD,.TRUE.,ERROR,*9999)
c        CALL MESH_REPOINT_DIAMETER(LDTMP2,nb,N_ELM,NE_OLD,
c     '    NPNE,500.d0,XP,ZD,.TRUE.,ERROR,*9999)
        DO nd=1,NDT
          IF(LDTMP2(nd).NE.0)THEN
            ne_min=LDTMP2(nd)
            NSTEM(ne_min,2)=NSTEM(ne_min,2)+1
            LDTMP1(nd)=LDTMP2(nd)
          ENDIF !LD
        ENDDO !nd

        ldnum2=0
        DO nd=1,NDT
          IF(LDTMP2(nd).GT.0) ldnum2=ldnum2+1
        ENDDO !nd
C       to get things right if an initial branch has no seed points
        NT_BNS=0
        DO N=1,NELIST2(0) !for each parent element in list
          IF(NSTEM(NELIST2(N),2).NE.0)THEN
            NT_BNS=NT_BNS+1
            NE_OLD(NT_BNS)=NELIST2(N)
          ENDIF
        ENDDO !N
        N_ELM=NT_BNS
        N_ELM_TEMP=N_ELM

C...    KSB 23/06/05 removing terminal elements (only 1 data pt
C...    associated) and corresponding data point from the group
        DO N=1,N_ELM
          ne_min=NE_OLD(N)
          IF(NSTEM(ne_min,2).LT.2)THEN  !terminal branch -remove data pt
            DO nd=1,NDT
              IF(LDTMP2(nd).EQ.ne_min)THEN
                LDTMP1(nd)=LDTMP2(nd)
                LDTMP2(nd)=0
                NE_OLD(N)=0
                N_ELM_TEMP=N_ELM_TEMP-1
              ENDIF
            ENDDO !nd
          ENDIF
        ENDDO
        numtb=0
        DO N=1,N_ELM
          IF(NE_OLD(N).EQ.0)THEN
            I=0
            DO WHILE((N+I.LT.N_ELM).AND.(NE_OLD(N+I).EQ.0))
              I=I+1
            ENDDO
            DO M=N,N_ELM-I
              NE_OLD(M)=NE_OLD(M+I)
            ENDDO !M
            numtb=numtb+1
          ENDIF !NE_OLD
        ENDDO !N
        N_ELM=N_ELM_TEMP !removes terminal elements from loop
        
      ENDIF !NELIST2.GT.1

      LD_NUM=0
      DO nd=1,NDT
        IF(LDTMP2(nd).GT.0) LD_NUM=LD_NUM+1
      ENDDO !nd
      ldnum2=LD_NUM

C      WRITE(OP_STRING,'('' gen   #brn  total# #term  #data'')')
      WRITE(OP_STRING,'(''  parent  #seeds  #terminal'')')
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      
C     Set initial values for local and global nodes and elements
      ne=NE_START !initialise mesh global element #
      np=NP_START !initialise mesh global node #
      noelem=NEELEM(0,nr) !initialise local element #
      nonode=NPNODE(0,nr) !initialise local node #
      nonode_start=NPNODE(0,nr)
      noelem_initial=noelem
      
      DO noelem_parent=1,NELIST2(0)
         ne_parent = NELIST2(noelem_parent)
         N_ELM=1
         ngen=3
C store data point numbers associated with this parent in LD
         DO nd=1,NDM
           LD(nd)=0
         ENDDO
         nld=0
         DO nd=1,NDM
           IF(LDTMP2(nd).EQ.ne_parent)THEN
             LD(nd)=ne_parent
             nld=nld+1
           ENDIF
         ENDDO
         N_ELM=1
         NE_OLD(1)=ne_parent
         numtb=0
      
C.... BIFURCATING DISTRIBUTIVE ALGORITHM
         DO WHILE(N_ELM.NE.0)   !while still some parent branches
C shouldn't need to use ngen. arbitrary for new algorithm structure
            ngen=ngen+1         !increment the generation from the parent
            NT_BNS=N_ELM 
            N_ELM=0
            numzero=0
            noelem_gen=0
            max_dist=0.d0
            Nremoved=0

            DO M=1,NT_BNS
               ne1=NE_OLD(M)    !parent global element #
C........Calculate COFM of each seed point set
               CALL MESH_COFM(LD,ne1,COFM,ZD,ERROR,*9999)
               ne2=NXI(-1,1,ne1) !grandparent global element #
               CALL ASSERT(ne2.GT.0,
     &              '>>Require more than 2 parent generations',ERROR,
     &              *9999)
               np1=NPNE(2,nb,ne1) !parent global end node #
               np2=NPNE(1,nb,ne1) !parent global start node #
               np3=NPNE(1,nb,ne2) !grandparent global start node #
               IF(NXI(1,1,ne2).EQ.ne1)THEN
                  IF(NXI(1,2,ne2).EQ.0)THEN
                     np4=NPNE(1,nb,NXI(-1,1,ne2))
                  ELSE
                     np4=NPNE(2,nb,NXI(1,2,ne2))
                  ENDIF
               ELSE
                  np4=NPNE(2,nb,NXI(1,1,ne2))
               ENDIF
               LENGTH_PARENT=0.d0
               DO nj=1,NJT
                  LENGTH_PARENT=LENGTH_PARENT+(XP(1,1,nj,np1)
     &              -XP(1,1,nj,np2))**2
               ENDDO            !nj
               LENGTH_PARENT=DSQRT(LENGTH_PARENT)
C........Split each set of seed points using the plane defined by the
C........parent branch and the centre of mass.

               CALL MESH_SPLIT(LD,nde1,nde2,POINT_LIMIT,ne1,
     &           ne,np1,np2,np3,COFM,XP,ZD,SS,ERROR,*9999)

               IF(SS(1).AND.SS(2))THEN
                  DO N=1,2
C............Calculate the centre of mass of each new set of seed
C............points
C                     CALL MESH_COFM(LD,ne+N,COFM,ZD,ERROR,*9999)
                     DO nj=1,3
                        XP1(nj)=XP(1,1,nj,np2) !start node
                     ENDDO
                  ENDDO
                  DO N=1,2
C............Set up arrays for new element and node
                     CALL GN1DNE(nb,NBJ,ne,NEELEM,NENP,ne1,NKJ,NKJE,
     &                 noelem,nonode,np,NP_INTERFACE,np1,NPNE,NPNODE,
     &                 nr,NRE,NVJE,NVJP,NXI,SE,.TRUE.,ERROR,*9999)
                     NORD(1,ne)=NORD(1,ne1)+1
                     noelem_gen=noelem_gen+1
                     CALL MESH_COFM(LD,ne,COFM,ZD,ERROR,*9999)
C............Generate a branch directed towards the centre of mass              
                     CALL MESH_BRANCH(LD,LD_NUM,LDTMP1,numzero,ne,np1,
     &                    np2,POINT_LIMIT,COFM,fraction,ANGLE_MINOR,
     '                    LENGTH_LIMIT,LENGTH_PARENT,min_length,XP,XP1,
     &                    ZD,BRANCH,SS(N),ERROR,*9999)
! at this point XP1 contains the location of the end point of the new branch
                     CALL ASSERT(np.LE.NPM,'>>NPM too small',ERROR,
     &                 *9999)
                     ldnum2=0
                     DO nd=1,NDT
                        IF(LD(nd).GT.0) ldnum2=ldnum2+1
                     ENDDO      !nd

                     IF(RESTRICT_BRANCH_PLANE)THEN
c                     IF(NORD(1,ne).LE.500)THEN !dummy, not used
                        IF(N.EQ.1)THEN
                           length=0.0d0
                           DO nj=1,NJT
                              XP(1,1,nj,np)=XP1(nj)
                              length=length+(XP(1,1,nj,np)
     &                          -XP(1,1,nj,np1))**2.0d0
                           ENDDO !nj
                           length=DSQRT(length)
                           CALL ASSERT(NVM.GE.2,'>>Increase NVM to 2',
     &                          ERROR,*9999)              
                           DO nj=1,NJT !unit vector for branch direction
                              XP(1,2,nj,np)=(XP(1,1,nj,np)
     &                          -XP(1,1,nj,np1))/length
                           ENDDO !nj
C.............Calculate the branching plane
                           DO nj=1,3
                              u(nj)=XP(1,2,nj,np1) !parent branch direction
                              v(nj)=XP(1,2,nj,np) !new branch direction
                           ENDDO !nj
                           IF(DABS(u(1)-v(1)).LT.1.d-4.AND.DABS(u(2)
     &                          -v(2)).LT.1.d-4.AND.DABS(u(3)
     &                          -v(3)).LT.1.d-4)THEN !perturb
C     Find the closest data point
                              distance=1.d6
                              DO nd=1,NDT
                                 IF(LD(nd).GT.0)THEN
                                    distance2=0.d0
                                    DO nj=1,NJT
                                       distance2=distance2+(ZD(1,nd)
     &                                   -XP1(1))**2
                                    ENDDO
                                    distance2=DSQRT(distance2)
                                    IF(distance2.LT.distance)THEN
                                       distance=distance2
                                       nd_closest=nd
                                    ENDIF
                                 ENDIF !LD(nd).GT.0
                              ENDDO !nd
C...................Direct the branch towards the closest data point
                              DO nj=1,NJT
                                 v(nj)=ZD(nj,nd_closest)-XP(1,1,nj,np1)
                              ENDDO
                              CALL NORMALISE(NJT,v,ERROR,*9999) !unit
                              DO nj=1,NJT
                                 XP(1,1,nj,np)=XP(1,1,nj,np1)+v(nj)
     &                             *length
                                 XP(1,2,nj,np)=v(nj)
                                 XP1(nj)=XP(1,1,nj,np)
                              ENDDO
                           ENDIF
                           CALL MESH_ANGLE_CHECK(N,np1,np2,np3,np,
     &                          angle_minor,LENGTH,length_u,XP,XP1,
     &                          ERROR,*9999)
                           
                           DO nj=1,NJT !unit vector for branch direction
                              XP(1,1,nj,np)=XP1(nj)
                              XP(1,2,nj,np)=(XP1(nj)-XP(1,1,nj,np1))
     &                          /length
                              v(nj)=XP(1,2,nj,np) !new branch direction
                           ENDDO !nj
                           length_w=dsqrt((xp1(1)-xp(1,1,1,np1))**2
     &                       +(xp1(2)-xp(1,1,2,np1))**2+(xp1(3)
     &                       -xp(1,1,3,np1))**2)
                           
                           CALL NORMALISE(NJT,v,ERROR,*9999) !unit vector
                           CALL CROSS(u,v,NRML) !calculate branching plane
                           CALL NORMALISE(NJT,NRML,ERROR,*9999) !unit vector
                           
                        ELSE
C.....Adjust the direction of the second branch s.t. in plane
                           length_w = 0.d0
                           DO nj=1,3
                              w(nj)=XP1(nj)-XP(1,1,nj,np1) !current vector
                              length_w = length_w + w(nj)*w(nj)
                           ENDDO !nj
                           length_w = DSQRT(length_w)
                           CALL NORMALISE(NJT,w,ERROR,*9999) !unit
                           
                           angle_uw=MAX(-1.d0,SCALAR(3,u,w))
                           angle_uw=MIN(1.d0,angle_uw)
                           angle_uw=DACOS(angle_uw)
                           
                           angle_uv=MAX(-1.d0,SCALAR(3,u,v))
                           angle_uv=MIN(1.d0,angle_uv)
                           angle_uv=DACOS(angle_uv)
                           
                           angle_vw=angle_uw+angle_uv
                           DO nj=1,NJT !set up system of linear equations
                              A(1,nj)=u(nj) !dotprod parent and new element
                              A(2,nj)=NRML(nj) !dotprod normal and new element
                              A(3,nj)=v(nj) !dotprod sibling and new element
                           ENDDO !nj
                           B(1)=DCOS(angle_uw)
                           B(2)=0.d0
                           B(3)=DCOS(angle_vw)
                           CALL MESH_A_X_EQ_B(A,B,w,ERROR,*9999)
                           CALL NORMALISE(NJT,w,ERROR,*9999) !unit
                           DO nj=1,NJT !unit vector for branch direction
                              XP1(nj)=XP(1,1,nj,np1)+length_w*w(nj)
     &                          *RATIO_LENGTHS
                           ENDDO !nj
                           CALL MESH_ANGLE_CHECK(N,np1,np2,np3,np,
     &                          angle_major,
     &                          LENGTH,length_u,XP,XP1,ERROR,*9999)
                           DO nj=1,NJT !unit vector for branch direction
                              XP(1,1,nj,np)=XP1(nj)
                              XP(1,2,nj,np)=(XP1(nj)-XP(1,1,nj,np1))
     &                          /length
                           ENDDO !nj
                           CALL ASSERT(NVM.GE.2,'>>Increase NVM to 2',
     &                       ERROR,*9999)
                           
                        ENDIF
                        
                     ELSE
                        length=0.0d0
                        DO nj=1,NJT
                           XP(1,1,nj,np)=XP1(nj)
                           length=length+(XP(1,1,nj,np)
     &                       -XP(1,1,nj,np1))**2.0d0
                        ENDDO   !nj
                        length=DSQRT(length)
                        CALL ASSERT(NVM.GE.2,'>>Increase NVM to 2',
     &                    ERROR,*9999)              
                        DO nj=1,NJT !unit vector for branch direction
                           XP(1,2,nj,np)=(XP(1,1,nj,np)-XP(1,1,nj,np1))
     &                       /length
                        ENDDO   !nj
                     ENDIF !RESTRICT_BRANCH_PLANE
                     
                     NSTEM(ne,1)=NSTEM(ne1,1)
                     IF(BRANCH.AND.SS(N))THEN
                        N_ELM=N_ELM+1
                        NE_TEMP(N_ELM)=ne !records parent elements 
                     ELSE       !terminal 
                        numtb=numtb+1
                        Nremoved=Nremoved+1
                     ENDIF      !BRANCH
                     
                  ENDDO         !N
C...........Make sure the branches aren't too much longer than the parent
                  CALL CHECK_BRANCH_LENGTH(np1,np2,np,XP,ERROR,*9999)
                  
C...  make sure that the branches have the correct rotation angle from
C...  the parent branching plane.
                  
                  IF(RESTRICT_BRANCH_ROTATION)THEN
C                  IF(NORD(1,ne).LE.500)THEN !dummy, not used
                     ROT_ANGLE=ROTATION_LIMIT*PI/180.d0
                     CALL CHECK_ROTATION_ANGLE(ne1,np4,np3,np2,np1,np,
     &                    np-1,ROT_ANGLE,XP,ERROR,*9999)
                     CALL ROTATION_ANGLE(np2,np1,np4,np,np-1,ANGLE,XP,
     &                    ERROR,*9999)
                     INTERNAL=.FALSE.
                     COUNT=0
                     DO WHILE(.NOT.INTERNAL.AND.COUNT.LT.2)
                        DO nj=1,3
                           XP1(nj)=XP(1,1,nj,np-1)
                        ENDDO
                        CALL CHECK_POINT_INTERNAL(IBT,IDO,INP,NBJ,
     &                       NDLIST,NEELEM,NEP,NHOST,NKJE,np-1,np1,NPF,
     &                       NPNE,nr_host,NVJE,NXI,SE,XA,XE,XIP,XP,XP1,
     &                       ZD,INTERNAL,.FALSE.,ERROR,*9999)
                        if(.not.internal)THEN
                           length=0.d0
                           DO nj=1,NJT
                              length=length+(XP(1,1,nj,np1)
     &                          -XP(1,1,nj,np-1))**2
                           ENDDO
                           length=DSQRT(length)
                           DO nj=1,NJT
                              XP(1,1,nj,np-1)=XP(1,1,nj,np1)+
     &                             0.5d0*length*XP(1,2,nj,np-1)
                              XP1(nj)=XP(1,1,nj,np-1)
                           ENDDO
                           CALL LIMIT_ANGLE(np1,np2,np3,np-1,XP,ERROR,
     &                       *9999)
                        ENDIF   !internal
                        COUNT=COUNT+1
                     ENDDO
                     IF(.NOT.INTERNAL)THEN
C...............Remove the branch from the list of next generation parents
                        offset=0
                        DO nes=1,N_ELM
                           IF(NE_TEMP(nes).EQ.ne-1) offset=1
                           NE_TEMP(nes)=NE_TEMP(nes+offset)
                        ENDDO
                        N_ELM=N_ELM-1
                        numtb=numtb+1
                        Nremoved=Nremoved+1
                        
C...............Remove the closest data point to the end of the branch
                        MIN_DIST=1.d10
                        DO nd=1,NDT
                           IF(LD(nd).NE.0)THEN
                              DIST=DSQRT((ZD(1,nd)-XP(1,1,1,np-1))
     &                          **2.d0+(ZD(2,nd)-XP(1,1,2,np-1))**2.d0
     &                          +(ZD(3,nd)-XP(1,1,3,np-1))**2.d0)
                              IF(DIST.LT.MIN_DIST)THEN
                                 nd_min=nd
                                 MIN_DIST=DIST
                              ENDIF !DIST
                           ENDIF !LD
                        ENDDO   !nd
                        LDTMP1(nd_min)=ne-1
                        LD(nd_min)=0
                        LD_NUM=LD_NUM-1
                        ldnum2=0
                        DO nd=1,NDT
                           IF(LD(nd).GT.0) ldnum2=ldnum2+1
                        ENDDO   !nd
                     ENDIF
                     
                     INTERNAL=.FALSE.
                     COUNT=0
                     DO WHILE(.NOT.INTERNAL.AND.COUNT.LT.2)
                        DO nj=1,3
                           XP1(nj)=XP(1,1,nj,np)
                        ENDDO
                        CALL CHECK_POINT_INTERNAL(IBT,IDO,INP,NBJ,
     &                     NDLIST,NEELEM,NEP,NHOST,NKJE,np,np1,NPF,NPNE,
     &                     nr_host,NVJE,NXI,SE,XA,XE,XIP,XP,XP1,ZD,
     &                     INTERNAL,.FALSE.,ERROR,*9999)
                        if(.not.internal)THEN
                           length=0.d0
                           DO nj=1,NJT
                              length=length+(XP(1,1,nj,np1)
     &                          -XP(1,1,nj,np))**2
                           ENDDO
                           length=DSQRT(length)
                           DO nj=1,NJT
                              XP(1,1,nj,np)=XP(1,1,nj,np1)+
     &                             0.5d0*length*XP(1,2,nj,np)
                           ENDDO
                           CALL LIMIT_ANGLE(np1,np2,np3,np,XP,ERROR,
     &                       *9999)
                        ENDIF
                        COUNT=COUNT+1
                     ENDDO      !while
                     IF(.NOT.INTERNAL)THEN
C...............Remove the branch from the list of next generation parents
                        offset=0
                        DO nes=1,N_ELM
                           IF(NE_TEMP(nes).EQ.ne) offset=1
                           NE_TEMP(nes)=NE_TEMP(nes+offset)
                        ENDDO
                        N_ELM=N_ELM-1
                        numtb=numtb+1
                        Nremoved=Nremoved+1
C...............Remove the closest data point to the end of the branch
                        MIN_DIST=1.d10
                        DO nd=1,NDT
                           IF(LD(nd).NE.0)THEN
                              DIST=DSQRT((ZD(1,nd)-XP(1,1,1,np-1))**2.d0
     &                          +(ZD(2,nd)-XP(1,1,2,np-1))**2.d0+
     &                          (ZD(3,nd)-XP(1,1,3,np-1))**2.d0)
                              IF(DIST.LT.MIN_DIST)THEN
                                 nd_min=nd
                                 MIN_DIST=DIST
                              ENDIF !DIST
                           ENDIF !LD
                        ENDDO   !nd
                        LDTMP1(nd_min)=ne
                        LD(nd_min)=0
                        LD_NUM=LD_NUM-1
                        ldnum2=0
                        DO nd=1,NDT
                           IF(LD(nd).GT.0) ldnum2=ldnum2+1
                        ENDDO   !nd
                     ENDIF
                     
                  ENDIF !RESTRICT_BRANCH_ROTATION
                  
                  IF(NXI(1,0,ne2).EQ.2)THEN
                     IF(ne1.EQ.NXI(1,1,ne2))THEN
                        nes=NXI(1,2,ne2)
                     ELSE
                        nes=NXI(1,1,ne2)
                     ENDIF
                  ELSE IF(NXI(1,0,ne2).EQ.1)THEN
                     nes=ne2
                     DO WHILE(NXI(1,0,nes).NE.2)
                        neb=nes
                        nes=NXI(-1,1,nes)
                     ENDDO
                     IF(neb.EQ.NXI(1,1,nes))THEN
                        nes=NXI(1,2,nes)
                     ELSE
                        nes=NXI(1,1,nes)
                     ENDIF
                  ENDIF
               ELSE
C.........Not enough seed points in the set during the split. Set
C.........parent(ne1) as a terminal branch. Find the closest data point
C.........and remove. This should be distinct from the later similar
C.........removal of elements with too few seed points because that step
C.........follows a mesh_repoint
                  Nremoved=Nremoved+1
                  numtb=numtb+1
                  np2=NPNE(2,nb,ne1)
                  MIN_DIST=1.d10
                  DO nd=1,NDT
                     IF(LD(nd).NE.0)THEN
c                        IF(LD(nd).EQ.nde1) LD(nd)=ne1
c                        IF(LD(nd).EQ.nde2) LD(nd)=ne1
                        DIST=DSQRT((ZD(1,nd)-XP(1,1,1,np2))**2.d0+
     &                     (ZD(2,nd)-XP(1,1,2,np2))**2.d0+(ZD(3,nd)
     &                     -XP(1,1,3,np2))**2.d0)
                        IF(DIST.LT.MIN_DIST)THEN
                           nd_min=nd
                           MIN_DIST=DIST
                        ENDIF   !DIST
                     ENDIF      !LD
                  ENDDO         !nd
                  LDTMP1(nd_min)=ne1
                  LD(nd_min)=0
                  LD_NUM=LD_NUM-1
                  ldnum2=0
                  DO nd=1,NDT
                     IF(LD(nd).GT.0) ldnum2=ldnum2+1
                  ENDDO         !nd
               ENDIF
c     CALL MINIMUM_ANGLE(NPNE(1,nb,ne1),NPNE(2,nb,ne1),np,
c     &      LENGTH_LIMIT,XP,ERROR,*9999)
               
            ENDDO               !M
C.......Copy the temporary list of branches to NE_OLD. These become the
C.......parent elements for the next branching        
            DO N=1,N_ELM
               NE_OLD(N)=NE_TEMP(N)
               NSTEM(NE_OLD(N),2)=0 !initialise the count of nd
            ENDDO               !N
            
!height of lung = 300.0, s.t. at ngen=30 will = 0.0
!increasing the distance limit increases the number of generations
            IF(REDUCING_DISTANCE_LIMIT)THEN
               DISTANCE_LIMIT=MAX(300.d0-ngen*SCALE_DIST_LIMIT,5.d0)
            ELSE
               DISTANCE_LIMIT=SCALE_DIST_LIMIT
            ENDIF

            NE_REACTIVATE(0)=0
            CALL MESH_REPOINT(LD,nb,N_ELM,NE_OLD,NE_REACTIVATE,
     &        NPNE,DISTANCE_LIMIT,XP,ZD,.FALSE.,ERROR,*9999)

            DO nd=1,NDT
               IF(LD(nd).NE.0)THEN
                  ne_min=LD(nd)
                  NSTEM(ne_min,2)=NSTEM(ne_min,2)+1
                  LDTMP1(nd)=LD(nd)
               ENDIF            !LD
            ENDDO               !nd

            nclos1=0
            nclos2=0
            
C......If there is 0 or 1 seed point for an element then set it as a
C......terminal and remove a single seed point from a set with >1 seed points
            N_ELM_TEMP=N_ELM
            DO N=1,N_ELM
               ne_min=NE_OLD(N)
               IF(NSTEM(ne_min,2).EQ.0)THEN !find closest point
                  nclos1=nclos1+1
                  np2=NPNE(2,nb,ne_min)
                  MIN_DIST=1.d10
                  DO nd=1,NDT
                     IF(LD(nd).NE.0)THEN
C check that this data point isn't the only one assigned to another terminal
                        IN_USE=.FALSE.
                        DO M=1,N_ELM
                           IF(NE_OLD(M).NE.0)THEN
                             IF(NSTEM(NE_OLD(M),2).EQ.1.AND.
     &                         LD(nd).EQ.NE_OLD(M))THEN
                                IN_USE=.TRUE.
                             ENDIF
                          ENDIF
                        ENDDO
                        IF(.NOT.IN_USE)THEN
                           DIST=DSQRT((ZD(1,nd)-XP(1,1,1,np2))**2.d0+
     &                          (ZD(2,nd)-XP(1,1,2,np2))**2.d0+(ZD(3,nd)
     &                          -XP(1,1,3,np2))**2.d0)
                           IF(DIST.LT.MIN_DIST)THEN
                              nd_min=nd
                              MIN_DIST=DIST
                           ENDIF !DIST
                        ENDIF
                     ENDIF      !LD
                  ENDDO         !nd
                  LDTMP1(nd_min)=ne_min
                  LD(nd_min)=0
                  LD_NUM=LD_NUM-1
                  N_ELM_TEMP=N_ELM_TEMP-1
                  NE_OLD(N)=0
                  numtb=numtb+1
                  Nremoved=Nremoved+1
                  
                  ldnum2=0
                  DO nd=1,NDT
                     IF(LD(nd).GT.0) ldnum2=ldnum2+1
                  ENDDO         !nd
                  
               ELSE IF(NSTEM(ne_min,2).EQ.1)THEN
                  nclos2=nclos2+1
                  DO nd=1,NDT
                     IF(LD(nd).EQ.ne_min)THEN
                        LDTMP1(nd)=LD(nd)
                        LD(nd)=0
                        NE_OLD(N)=0
                        N_ELM_TEMP=N_ELM_TEMP-1
                        LD_NUM=LD_NUM-1
                        numtb=numtb+1
                        Nremoved=Nremoved+1
                     ENDIF
                  ENDDO         !nd
                  
                  ldnum2=0
                  DO nd=1,NDT
                     IF(LD(nd).GT.0) ldnum2=ldnum2+1
                  ENDDO         !nd
                  
               ENDIF            !NSTEM
            ENDDO               !N

            DO N=1,N_ELM
               IF(NE_OLD(N).EQ.0)THEN
                  I=0
                  DO WHILE((N+I.LT.N_ELM).AND.(NE_OLD(N+I).EQ.0))
                     I=I+1
                  ENDDO
                  DO M=N,N_ELM-I
                     NE_OLD(M)=NE_OLD(M+I)
                  ENDDO         !M
               ENDIF            !NE_OLD
            ENDDO               !N
            N_ELM=N_ELM_TEMP
            
            CALL ILISTRMDUP(N_ELM,NE_OLD,ERROR,*9999)

            ldnum2=0
            DO nd=1,NDT
               IF(LD(nd).GT.0) ldnum2=ldnum2+1
            ENDDO               !nd
            
            numtbsum=numtbsum+numtb
         ENDDO                  !WHILE NT_BNS
         WRITE(OP_STRING,'(I7,I8,I9)') ne_parent,nld,numtb
         CALL WRITES(IODI,OP_STRING,ERROR,*9999)

      ENDDO !for each initial parent element

      DO nd=1,NDT
        ne2=LDTMP1(nd)
        IF(ne2.NE.0.AND.NXI(1,0,ne2).EQ.0)THEN
          nb=NBJ(1,ne2)
          np1=NPNE(2,nb,ne2) !end node
          LD(nd)=LD_NP(nd) !replace with temporary storage from start
          LD_NP(nd)=np1 !records closest branch node
        ELSE
          LD(nd)=0
          LDTMP1(nd)=0
        ENDIF
      ENDDO !nd

      NPT(nr)=np !highest node # in nr
      NET(nr)=ne !highest element # in nr
      NPT(0)=NPT(nr) !highest node# in all regions
      NET(0)=NET(nr) !highest element# in all regions
      NEELEM(0,nr)=noelem !# of elements in nr
      NPNODE(0,nr)=nonode !# of nodes in nr
      NPNODE(0,0)=NPNODE(0,0)+noelem !number nodes in all regions
      NEELEM(0,0)=NEELEM(0,0)+nonode !number elements in all regions
      CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
      CALL CALC_NENP_1D(nb,NEELEM,NENP,NPNE,nr,ERROR,*9999)
      CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)

C.....Make the length of all terminal branches sensible
C      DO noelem=1,NEELEM(0,nr)
C        ne=NEELEM(noelem,nr)
C        IF(NXI(1,0,ne).EQ.0)THEN
C          IF(.NOT.INLIST(ne,NELIST2,NELIST2(0),itemp))THEN
C            nb=NBJ(1,ne)
C            np1=NPNE(1,nb,ne)
C            np2=NPNE(2,nb,ne)
C            length=0.d0
C            DO nj=1,NJT
C              length=length+(XP(1,1,nj,np1)-XP(1,1,nj,np2))**2.0d0
C            ENDDO !nj
C            length=DSQRT(length)
C            length=MIN(1.4d0,length) !shorter than 1.4
C            length=MAX(1.0d0,length) !longer than 1.0
C            DO nj=1,NJT
C              XP(1,1,nj,np2)=XP(1,1,nj,np1)+XP(1,2,nj,np2)*length
C            ENDDO
C          ENDIF
C        ENDIF
C      ENDDO
              
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        NHOST(np)=0
      ENDDO
      nclos1=0
      nclos2=0
C count the number of data points associated with end node np1 of
C element ne2. Store in LDTMP2(np1-NP_START)
      DO nd=1,NDT
        ne2=LDTMP1(nd)
        IF(ne2.NE.0)THEN
          nb=NBJ(1,ne2)
          np1=NPNE(2,nb,ne2) !end node
          NHOST(np1)=NHOST(np1)+1
          nclos1=nclos1+1
        ENDIF
      ENDDO !nd
C Check that only 1 data point is associated with the terminal elements
      ncount=0
      DO noelem=NE_START+1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        nb=NBJ(1,ne)
        IF(NXI(1,0,ne).EQ.0)THEN
          ncount=ncount+1
          np=NPNE(2,nb,ne)
          IF(NHOST(np).EQ.0)THEN
            WRITE(OP_STRING,'('' WARNING no data point for node'',I8)')
     &       np
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ELSE
            nd=NHOST(np)
            DO nj=1,3
              XP(1,2,nj,np)=ZD(nj,nd)
              NVJP(nj,np)=2 !set 2 versions at end node
            ENDDO
          ENDIF
          IF(NHOST(np).GT.1)THEN
            WRITE(OP_STRING,'('' WARNING >1 data point for node'',I8)')
     &        np
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ELSE
            nclos2=nclos2+1
          ENDIF
        ENDIF
      ENDDO

      IF(CALL_FIEL.AND.nj_radius.GT.0)THEN
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          NVJP(nj_radius,np)=NENP(np,0,nr)
          NKJ(nj_radius,np)=1
        ENDDO !nonode
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          NBJ(nj_radius,ne)=nb
          DO nn=1,NNT(nb)
            np=NPNE(nn,nb,ne)
            DO i=1,NENP(np,0,nr)
              ne2=NENP(np,i,nr)
              IF(ne2.EQ.ne)THEN
                NVJE(nn,nb,nj_radius,ne)=i
              ENDIF
            ENDDO !i
          ENDDO !nn
          DO nk=1,NKT(nn,nb)
            NKJE(nk,nn,nj_radius,ne)=1
          ENDDO !nk
        ENDDO !noelem
      ENDIF !CALL_FIEL
      
      CALL_TREE=.TRUE. !logical to specify tree generation has been performed

      CALL OPENF(IOFILE2,'DISK','xicoords.ipxi',
     '  'NEW','SEQUEN','FORMATTED',160,ERROR,*9999)

      FMT='(2I7,4E25.16)'
      DO nonode=nonode_start,NPNODE(0,nr)
         np=NPNODE(nonode,nr)
         ne=NEP(np)
         WRITE(IOFILE2,FMT) np,ne,(XIP(nj,np),nj=1,3)
      ENDDO

      CLOSE(IOFILE2)

      CALL EXITS('GNBDMESH')
      RETURN
 9999 CALL ERRORS('GNBDMESH',ERROR)
      CALL EXITS('GNBDMESH')
      RETURN 1
      END

      

