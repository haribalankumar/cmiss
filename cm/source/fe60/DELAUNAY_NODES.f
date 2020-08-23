      SUBROUTINE DELAUNAY_NODES(IBT,IDO,INP,NBJ,NDIVISION,NEELEM,
     '  NELIST,NFF,NFFACE,NFLIST,NFNP,NKJ,NKJE,NNF,NPF,NPLIST,NPNE,
     '  NPNODE,NP_INTERFACE,nr,nr_target,nr_tree,NVCNP,NVJE,NVJP,
     '  NXI,RADIUS,SE,XE,XP,XP_B,XP_IB,XPNP,INTERNAL_NODES,REGULAR,
     '  ERROR,*)

C#### Subroutine: DELAUNAY_NODES
C###  Description:
C###    DELAUNAY_NODES calculates boundary (B) and internal boundary
C###  (IB) nodes on the surface of host elements.  Internal (IN) nodes
C###  are calculated only if INTERNAL_NODES is .TRUE.

C***  This routine and ones that it calls is currently under
C***  development - a lot of tidying up is still to be done.  Please
C***  ignore any unused variables etc at the moment.

C***  Created by: Merryn Howatson Tawhai, August 2001

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NBJ(NJM,NEM),NDIVISION,NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NFF(6,NEM),
     '  NFFACE(0:NF_R_M,NRM),NFLIST(0:NFM),NFNP(0:10,NP_R_M),
     '  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),NNF(0:17,6,NBFM),NPF(9,NFM),
     '  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NP_INTERFACE(0:NPM,0:3),nr,nr_target,nr_tree,
     '  NVCNP(NP_R_M),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJP(NJM,NPM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 RADIUS,SE(NSM,NBFM,NEM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM),
     &  XP_B(3,NE_R_M),XP_IB(3,NE_R_M),XPNP(3,NP_R_M)
      LOGICAL INTERNAL_NODES,REGULAR
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER BOUND,i,in1,in2,IT,j,k,nb,N_BDRY,N_IBDRY,ne,ne2,nef,nef2,
     '  nf,nf_count,nf2,nj,nnode,nobdry,noelem,noface,node,nonode,
     &  nonode2,n_al_ax,no_b,no_ib,no_in,no_in2,N_INTNL,np,np1,np2,
     '  np3,np4,NUM_ALVEOLI1,NUM_ALVEOLI2,N_XI,XI1,XI2,XI3,nnob,nob,
     &  XI_LIST1(4),XI_LIST2(4),fact,XID1,XID2,XI_USE
      INTEGER NXF(-3:3,0:1,NFM)
      REAL*8 ANG_INCREMENT1,ANG_INCREMENT2,ANG1,ANG2,B(3),b_inc,DOT,
     '  L_W,PXI,R,R_NEW,R0,R2,R2_NEW,ROTATE,
     '  U(3),USER_TOL,X(3),XE2(NSM,NJM),XI(3),XI_DIST,XILM(6),XIF(3),
     '  XIS(3),XIS_USE,XI_HALF(3),XI_INCR,XI_LIMIT(3),XP_I(3),XNORM(3),
     '  X1(3),Z(3),Z_LEN,TEMP_B(3),XI_TEST(3),XP_IN(NJM,NVCM)
      LOGICAL ADD_NF,BOUNDARY,COLLAPSE(4),FOUND,FOUND2,KEEP,KEEP2
      CHARACTER STRING*255
      DATA XILM/0.d0,1.d0,0.d0,1.d0,0.d0,1.d0/


      CALL ENTERS('DELAUNAY_NODES',*9999)

      IF(REGULAR)THEN
        XI_DIST=0.5d0*1.d0/NDIVISION !will use this when new alveolar mesh is complete
      ELSE
        XI_DIST=1.d0/NDIVISION
      ENDIF
      USER_TOL=LOOSE_TOL
      XI_INCR=0.2d0*XI_DIST
      nobdry=0 !B node indicator (will be same for B nodes with same IB)
      no_b=0 !B node counter
      no_ib=0 !IB node counter
      nb=NBJ(1,NEELEM(1,nr)) !nb for 1st ne of host region
      DO nj=1,3
        DO np=1,NP_R_M
          NFNP(0,np)=0 !initialising
          XPNP(nj,np)=0.d0
        ENDDO !np
        DO nf=1,NFM
          NXF(-nj,0,nf)=0
          NXF(nj,0,nf)=0
        ENDDO !nf
        DO ne=1,NE_R_M
          XP_B(nj,ne)=0.d0
        ENDDO !ne
      ENDDO !nj
C*** For each face in host region (nr) check whether it is a boundary
C*** face (will only be in 1 element).  Record in NFLIST.
      nf_count=0 !counts the number of boundary faces/lines
      DO noface=1,NFFACE(0,nr) !loop over faces/lines
        nf=NFFACE(noface,nr) !global face number
        IF(NPF(5,nf).EQ.1)THEN !single element adjoining face/line nf
          nf_count=nf_count+1 !increment number of boundary faces/lines
          NFLIST(nf_count)=NFFACE(noface,nr) !store boundary face/line
          ne=NPF(6,nf) !element to which face associated
C***    Put host element information into XE.
          CALL VORO_XPXE(NBJ,ne,NKJE,NPNE,NVJE,SE,XE,XP,ERROR,*9999)
          nef=NPF(8,nf) !local face #
          XI1=NNF(1,nef,nb)
          XI(XI1)=XILM(nef) !0.d0 or 1.d0
          nonode2=1
          XI2=NPF(1,nf) !first Xi direction for face nf
          XI3=NPF(3,nf) !second Xi direction for face nf
          in1=XI2
          in2=XI3
          np1=NPNE(NNF(2,nef,nb),nb,ne)
          np2=NPNE(NNF(3,nef,nb),nb,ne)
          np3=NPNE(NNF(4,nef,nb),nb,ne)
          np4=NPNE(NNF(5,nef,nb),nb,ne)
          XI(XI2)=XI_INCR
          XI(XI3)=XI_INCR
          DO i=1,4
            COLLAPSE(i)=.FALSE.
          ENDDO !i
          ADD_NF=.TRUE.
          IF(np1.EQ.np2)THEN
            XI(XI2)=0.5d0
            COLLAPSE(3)=.TRUE. !-XI3 edge collapsed
          ELSE IF(np1.EQ.np3)THEN
            XI(XI3)=0.5d0
            COLLAPSE(1)=.TRUE. !-XI2 edge collapsed
          ENDIF
          CALL NPNFNP(IBT,IDO,INP,in1,in2,NBJ(1,ne),nf,NFNP,np1,XE,
     '      XI,XPNP,ADD_NF,ERROR,*9999)
          XI(XI2)=1.d0-XI_INCR
          XI(XI3)=XI_INCR
          IF(np2.EQ.np1)THEN
            ADD_NF=.FALSE.
          ELSE IF(np2.EQ.np4)THEN
            XI(XI3)=0.5d0
            COLLAPSE(2)=.TRUE. !XI2 edge collapsed
          ENDIF
          CALL NPNFNP(IBT,IDO,INP,XI2,XI3,NBJ(1,ne),nf,NFNP,np2,XE,
     '      XI,XPNP,ADD_NF,ERROR,*9999)
          XI(XI2)=XI_INCR
          XI(XI3)=1.d0-XI_INCR
          IF(np3.EQ.np1)THEN
            ADD_NF=.FALSE.
          ELSE IF(np3.EQ.np4)THEN
            XI(XI2)=0.5d0
            COLLAPSE(4)=.TRUE. !XI3 edge collapsed
          ENDIF
          CALL NPNFNP(IBT,IDO,INP,in1,in2,NBJ(1,ne),nf,NFNP,np3,XE,
     '      XI,XPNP,ADD_NF,ERROR,*9999)
          XI(XI2)=1.d0-XI_INCR
          IF(np4.EQ.np2.OR.np4.EQ.np3) ADD_NF=.FALSE.
          CALL NPNFNP(IBT,IDO,INP,in1,in2,NBJ(1,ne),nf,NFNP,np4,XE,
     '      XI,XPNP,ADD_NF,ERROR,*9999)
C*** set up NXF
          CALL NFNXF(nb,ne,nf,NFF,NNF,NPF,NXF,NXI,XI2,XI3,COLLAPSE,
     '      ERROR,*9999)
        ENDIF !NPF
      ENDDO !noface (nf)
      NFLIST(0)=nf_count !total number of boundary faces/lines
      DO nonode=1,NPNODE(0,nr) !calculate mean IB projections
        np=NPNODE(nonode,nr)
        IF(NFNP(0,np).NE.0)THEN !has surrounding boundary faces
          ne=NPF(6,NFNP(1,np)) !element for first face
C***    Put host element information into XE.
          CALL VORO_XPXE(NBJ,ne,NKJE,NPNE,NVJE,SE,XE,XP,ERROR,*9999)
          BOUNDARY=.FALSE. !determine whether IB should be on face
          DO noface=2,NFNP(0,np)
            ne2=NPF(6,NFNP(noface,np))
            IF(ne2.NE.ne) BOUNDARY=.TRUE. !more than one ne for faces
          ENDDO !noface
C*** Estimate the IB node position (initial guess)
          DO nj=1,3
            XI_TEST(nj)=0.5d0
          ENDDO !nj
          CALL XNORMXI(IBT,IDO,INP,in1,in2,NBJ(1,ne),X,XE,
     '      XI_TEST,XNORM,ERROR,*9999) !get X from XI
          DO nj=1,3
            Z(nj)=X(nj)-XP(1,1,nj,np) !direction to element centre
          ENDDO
          CALL NORMALISE(3,Z,ERROR,*9999) !normal vector for IB
          DO nj=1,3
            XPNP(nj,np)=XP(1,1,nj,np)+Z(nj)*RADIUS !IB node position
            TEMP_B(nj)=XP(1,1,nj,np)-Z(nj)*RADIUS !first B node
          ENDDO !nj
          IF(BOUNDARY)THEN !make IB node sit on the boundary face
            nf=NFNP(1,np) !global face #
            ne=NPF(6,nf) !XE info should already be correct for ne
            CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),Z_LEN,USER_TOL,XE,
     '        XI,XPNP(1,np),FOUND,.TRUE.,ERROR,*9999) !XI from XPNP
            nef=NPF(8,nf) !local face # in element ne
            XI1=NNF(1,nef,nb) !Xi direction normal to face nf
            IF(nef.EQ.1.OR.nef.EQ.3.OR.nef.EQ.5)THEN
              XIS(XI1)=XI_INCR !1st, 3rd, or 5th face -> at Xi=0
              XI_TEST(XI1)=0.d0
            ELSE
              XIS(XI1)=1.d0-XI_INCR !2nd, 4th, or 6th face -> at Xi=1
              XI_TEST(XI1)=1.d0
            ENDIF !NPF
            XI2=NPF(1,nf) !first Xi direction for face
            XI3=NPF(3,nf) !second Xi direction for face
            DO node=1,4 !# of nodes in the face
              nonode2=NNF(node+1,nef,nb)
              np2=NPNE(nonode2,nb,ne) !global node #
              IF(np2.EQ.np) nnode=node
            ENDDO !nonode
            IF(nnode.EQ.1)THEN
              XIS(XI3)=0.d0
              XIS(XI2)=XI_INCR
            ELSE IF(nnode.EQ.2)THEN
              XIS(XI3)=0.d0
              XIS(XI2)=1.d0-XI_INCR
            ELSE IF(nnode.EQ.3)THEN
              XIS(XI3)=1.d0
              XIS(XI2)=XI_INCR
            ELSE IF(nnode.EQ.4)THEN
              XIS(XI3)=1.d0
              XIS(XI2)=1.d0-XI_INCR
            ENDIF
            DO i=1,3
              XI(i)=XIS(i)
            ENDDO !i
            FOUND2=.FALSE.
            DO WHILE(.NOT.FOUND2)
              CALL XNORMXI(IBT,IDO,INP,XI2,XI3,NBJ(1,ne),X,XE,
     '          XI,XNORM,ERROR,*9999) !get X from XI
              DO nj=1,3
                Z(nj)=X(nj)-XP(1,1,nj,np)
              ENDDO !nj
              CALL NORMALISE(3,Z,ERROR,*9999)
              DO nj=1,3
                XPNP(nj,np)=XP(1,1,nj,np)+RADIUS*Z(nj)
              ENDDO !nj
C*** Checking that XI is on the boundary face
              CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),Z_LEN,USER_TOL,XE,
     '          XI,XPNP(1,np),FOUND,.TRUE.,ERROR,*9999) !XI from XPNP
              IF(DABS(XI(XI3)-XIS(XI3)).LT.LOOSE_TOL)THEN
                FOUND2=.TRUE.
              ELSE
                XI(XI3)=XIS(XI3)
              ENDIF
            ENDDO !WHILE
          ENDIF !BOUNDARY
          no_ib=no_ib+1
          nobdry=nobdry+1
c          no_b=no_b+1
          nnob=0
          DO nj=1,3
            XP_IB(nj,no_ib)=XPNP(nj,np) !store IB coordinates
c            XP_B(nj,no_b)=TEMP_B(nj)
          ENDDO !nj
          DO noface=1,NFNP(0,np) !for each surrounding face
            nf=NFNP(noface,np)
C***  Find point on surface that is normal to IB
            ne2=NPF(6,nf) !host element
C***    Put host element information into XE2.
            CALL VORO_XPXE(NBJ,ne2,NKJE,NPNE,NVJE,SE,XE2,XP,ERROR,*9999)
            CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne2),Z_LEN,USER_TOL,XE2,
     '        XI,XP_IB(1,no_ib),FOUND,.TRUE.,ERROR,*9999) !XI from XPNP
            nef=NPF(8,nf) !local face number
            XI1=NNF(1,nef,nb) !normal to face
            XIS(XI1)=XILM(nef)
            XI2=NPF(1,nf) !first Xi direction for face nf
            XI3=NPF(3,nf) !second Xi direction for face nf
            in1=XI2
            in2=XI3
            KEEP=.TRUE.
            CALL XINORM(IBT,IDO,INP,XI2,XI3,NBJ(1,ne),XI1,RADIUS,
     '        X,XE2,XI,XIS(XI1),XNORM,XP_IB(1,no_ib),KEEP,ERROR,*9999)
            IF(KEEP)THEN
              DOT=-XNORM(1)*Z(1)-XNORM(2)*Z(2)-XNORM(3)*Z(3)
              IF(DOT.EQ.1.d0)THEN
                L_W=RADIUS
              ELSE
                L_W=0.5d0*DSIN(PI-2.d0*DACOS(DOT))
     '            *RADIUS/DSIN(DACOS(DOT))
              ENDIF
              DO nj=1,3
                TEMP_B(nj)=X(nj)+XNORM(nj)*L_W
              ENDDO !nj
            ELSE
              CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),USER_TOL,Z_LEN,
     '          XE,XI,XP_IB(1,no_ib),FOUND,.TRUE.,ERROR,*9999)
              XI(XI1)=XIS(XI1)
              CALL XINORM1D(IBT,IDO,INP,NBJ(1,ne),XI3,X,XE,XI,
     '          XP_IB(1,no_ib),ERROR,*9999)
              DO nj=1,3
                TEMP_B(nj)=2.d0*X(nj)-XP_IB(nj,no_ib)
              ENDDO !nj
            ENDIF !KEEP
            KEEP2=.TRUE.
            DO nob=no_b+1,no_b+nnob
              IF(DABS(XP_B(1,nob)-TEMP_B(1)).LE.LOOSE_TOL.AND.
     '          DABS(XP_B(2,nob)-TEMP_B(2)).LE.LOOSE_TOL.AND.
     '          DABS(XP_B(3,nob)-TEMP_B(3))
     '          .LE.LOOSE_TOL) KEEP2=.FALSE.
            ENDDO !nob
            IF(KEEP2)THEN
              nnob=nnob+1
              CALL ASSERT(no_b+nnob.LE.NE_R_M,'>>Increase NE_R_M',ERROR,
     '          *9999)
              DO nj=1,3
                XP_B(nj,no_b+nnob)=TEMP_B(nj)
              ENDDO !nj
c              BDRY(no_b+nnob)=nobdry
            ENDIF !KEEP2
          ENDDO !noface
          no_b=no_b+nnob
        ENDIF !NFNP(0,np).NE.0
      ENDDO !nonode

C*** Find B and IB for lines
      DO noface=1,NFLIST(0)
        nf=NFLIST(noface)
        ne=NPF(6,nf)
        nef=NPF(8,nf) !local face #
C***    Put host element information into XE.
        CALL VORO_XPXE(NBJ,ne,NKJE,NPNE,NVJE,SE,XE,XP,ERROR,*9999)
        XI1=NNF(1,nef,nb) !normal to face
        IF(nef.EQ.1.OR.nef.EQ.3.OR.nef.EQ.5)THEN
          XIS(XI1)=0.d0 !1st, 3rd, or 5th face -> at Xi=0
          XIF(XI1)=XI_INCR
        ELSE
          XIS(XI1)=1.d0 !2nd, 4th, or 6th face -> at Xi=1
          XIF(XI1)=-XI_INCR
        ENDIF !NPF
        XI2=NPF(1,nf)
        XI3=NPF(3,nf)
        in1=XI2
        in2=XI3
        XI_LIST1(1)=XI2
        XI_LIST1(2)=XI2
        XI_LIST1(3)=XI3
        XI_LIST1(4)=XI3
        XI_LIST2(1)=XI3
        XI_LIST2(2)=XI3
        XI_LIST2(3)=XI2
        XI_LIST2(4)=XI2
        fact=1
        DO i=1,4
          fact=fact*(-1)
          XID1=XI_LIST1(i) !XI2,XI2,XI3,XI3
          XID2=XI_LIST2(i) !XI3,XI3,XI2,XI2
          IF(fact.EQ.1)THEN
            XIS(XID1)=1.d0 !along the +XID1 edge
            XIF(XID1)=-XI_INCR
          ELSE IF(fact.EQ.-1)THEN
            XIS(XID1)=0.d0 !along the -XID1 edge
            XIF(XID1)=XI_INCR
          ENDIF
          nf2=NXF(fact*XID1,1,nf) !global face # of face 'next door'
          IF(nf2.GT.nf)THEN !hasn't been done yet
            XI(XID2)=XI_DIST !first placement of a centre point
            DO WHILE(XI(XID2).LT.1.d0-0.5d0*XI_DIST)
C*** Calculate the 'centre' point on edge, common to both faces
              DO nj=1,3
                Z(nj)=0.d0
              ENDDO !nj
              XI(XI1)=XIS(XI1) !0 or 1
              XI(XID1)=XIS(XID1) !0 or 1
              in1=NPF(1,nf)
              in2=NPF(3,nf)
              CALL XNORMXI(IBT,IDO,INP,in1,in2,NBJ(1,ne),X,XE,
     '          XI,XNORM,ERROR,*9999) !get X from XI
              DO nj=1,3
                XP_I(nj)=X(nj) !record circumcentre
              ENDDO !nj
C*** Decide whether neighbouring face is in same element, or other
              ne2=NPF(6,nf2) !neighbouring element
              nef2=NPF(8,nf2) !local face #
              IF(ne2.EQ.ne)THEN !same element
C*** Put IB position to be at circumcentre + radius x direction
                no_ib=no_ib+1
                nobdry=nobdry+1
                nnob=0
                DO nj=1,3
                  XI_TEST(nj)=0.5d0
                ENDDO
                CALL XNORMXI(IBT,IDO,INP,in1,in2,NBJ(1,ne),X,XE,
     '            XI_TEST,XNORM,ERROR,*9999) !get X from XI
                DO nj=1,3
                  Z(nj)=X(nj)-XP_I(nj)
                ENDDO
                CALL NORMALISE(3,Z,ERROR,*9999) !normal vector for IB
                DO nj=1,3
                  XP_IB(nj,no_ib)=XP_I(nj)+RADIUS*Z(nj) !IB position
                ENDDO !nj
C*** CALCULATE B NODE FOR FIRST FACE
                in1=NPF(1,nf)
                in2=NPF(3,nf)
C*** Call XINORM to find the point on the surface that is orthogonal
C***  to the XP_IB node.
                KEEP=.TRUE.
                CALL XINORM(IBT,IDO,INP,in1,in2,NBJ(1,ne),XI1,
     '            RADIUS,X,XE,XI,XIS(XI1),XNORM,XP_IB(1,no_ib),KEEP,
     '            ERROR,*9999)
                IF(KEEP)THEN
                  DOT=-XNORM(1)*Z(1)-XNORM(2)*Z(2)-XNORM(3)*Z(3)
                  IF(DOT.EQ.1.d0)THEN
                    L_W=RADIUS
                  ELSE
C*** Outward projection distance calculated to give B point at RADIUS
C***  from the circumcentre
                    L_W=0.5d0*DSIN(PI-2.d0*DACOS(DOT))
     '                *RADIUS/DSIN(DACOS(DOT))
                  ENDIF
                  DO nj=1,3
                    TEMP_B(nj)=X(nj)+XNORM(nj)*L_W !B node position
                  ENDDO !nj
C*** Check for identical B nodes.  Will occur when faces have same
C***  normal.  i.e. for faces from different elements
                  KEEP2=.TRUE.
                  DO nob=no_b+1,no_b+nnob
                    IF(DABS(XP_B(1,nob)-TEMP_B(1)).LE.LOOSE_TOL.AND.
     '                DABS(XP_B(2,nob)-TEMP_B(2)).LE.LOOSE_TOL.AND.
     '                DABS(XP_B(3,nob)-TEMP_B(3)).LE.LOOSE_TOL)
     '                KEEP2=.FALSE.
                  ENDDO !nob
                  IF(KEEP2)THEN
                    nnob=nnob+1
                    CALL ASSERT(no_b+nnob.LE.NE_R_M,'>>Increase NE_R_M',
     '                ERROR,*9999)
                    DO nj=1,3
                      XP_B(nj,no_b+nnob)=TEMP_B(nj) !store B position
                    ENDDO !nj
C                    BDRY(no_b+nnob)=nobdry !same nobdry # as for IB node
                  ENDIF !KEEP2
                ENDIF !KEEP

C*** CALCULATE B NODE FOR SECOND FACE
                in1=NPF(1,nf2)
                in2=NPF(3,nf2)
                XIS_USE=XILM(nef2)
                XI_USE=6-in1-in2
                KEEP=.TRUE.
                CALL XINORM(IBT,IDO,INP,in1,in2,NBJ(1,ne),XI_USE,
     '            RADIUS,X,XE,XI,XIS_USE,XNORM,XP_IB(1,no_ib),KEEP,
     '            ERROR,*9999)
                IF(KEEP)THEN
                  DOT=-XNORM(1)*Z(1)-XNORM(2)*Z(2)-XNORM(3)*Z(3)
                  IF(DOT.EQ.1.d0)THEN
                    L_W=RADIUS
                  ELSE
C*** Outward projection distance calculated to give B point at RADIUS
C***  from the circumcentre
                    L_W=0.5d0*DSIN(PI-2.d0*DACOS(DOT))
     '                *RADIUS/DSIN(DACOS(DOT))
                  ENDIF
                  DO nj=1,3
                    TEMP_B(nj)=X(nj)+XNORM(nj)*L_W !B node position
                  ENDDO !nj
C*** Check for identical B nodes.  Will occur when faces have same
C***  normal.  i.e. for faces from different elements
                  KEEP2=.TRUE.
                  DO nob=no_b+1,no_b+nnob
                    IF(DABS(XP_B(1,nob)-TEMP_B(1)).LE.LOOSE_TOL.AND.
     '                DABS(XP_B(2,nob)-TEMP_B(2)).LE.LOOSE_TOL.AND.
     '                DABS(XP_B(3,nob)-TEMP_B(3)).LE.LOOSE_TOL)
     '                KEEP2=.FALSE.
                  ENDDO !nob
                  IF(KEEP2)THEN
                    nnob=nnob+1
                    CALL ASSERT(no_b+nnob.LE.NE_R_M,'>>Increase NE_R_M',
     '                ERROR,*9999)
                    DO nj=1,3
                      XP_B(nj,no_b+nnob)=TEMP_B(nj) !store B position
                    ENDDO !nj
C                    BDRY(no_b+nnob)=nobdry !same nobdry # as for IB node
                  ENDIF !KEEP2
                ENDIF !KEEP
              ELSE IF(ne2.NE.ne)THEN !from different elements
C*** Calculate IB node to be on common face of ne and ne2
                XI(XI1)=XIS(XI1)+XIF(XI1) !XI1 face + offset
                XI(XID1)=XIS(XID1) !at XID1 edge
                in1=NPF(1,nf)
                in2=NPF(3,nf)
                CALL XNORMXI(IBT,IDO,INP,in1,in2,NBJ(1,ne),X,XE,
     '            XI,XNORM,ERROR,*9999) !get X from XI
                DO nj=1,3
                  Z(nj)=X(nj)-XP_I(nj)
                ENDDO !nj
                CALL NORMALISE(3,Z,ERROR,*9999) !get Z direction
                no_ib=no_ib+1
                nobdry=nobdry+1
                nnob=0
                DO nj=1,3
                  XP_IB(nj,no_ib)=XP_I(nj)+RADIUS*Z(nj)
                ENDDO !nj
C*** CALCULATE B NODE FOR FIRST FACE
C*** Calculate XI coordinates of IB point
                CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),USER_TOL,Z_LEN,
     '            XE,XI,XP_IB(1,no_ib),FOUND,.TRUE.,ERROR,*9999)
C*** Find closest orthogonal point on face nf
                XI(XI1)=XIS(XI1) !on XI1 surface, XI(XID1) will change
                in1=NPF(1,nf)
                in2=NPF(3,nf)
                KEEP=.TRUE.
                CALL XINORM(IBT,IDO,INP,in1,in2,NBJ(1,ne),XI1,
     '            RADIUS,X,XE,XI,XIS(XI1),XNORM,XP_IB(1,no_ib),KEEP,
     '            ERROR,*9999)
                IF(.NOT.KEEP)THEN !check ne2
                  CALL VORO_XPXE(NBJ,ne2,NKJE,NPNE,NVJE,SE,XE2,XP,
     '              ERROR,*9999)
                  in1=NPF(1,nf2)
                  in2=NPF(3,nf2)
                  KEEP=.TRUE.
                  CALL XINORM(IBT,IDO,INP,in1,in2,NBJ(1,ne2),XI1,
     '              RADIUS,X,XE2,XI,XIS(XI1),XNORM,XP_IB(1,no_ib),KEEP,
     '              ERROR,*9999)
                  IF(.NOT.KEEP)THEN
                    CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),USER_TOL,
     '                Z_LEN,XE,XI,XP_IB(1,no_ib),FOUND,.TRUE.,ERROR,
     '                *9999)
                    XI(XI1)=XIS(XI1)
                    CALL XINORM1D(IBT,IDO,INP,NBJ(1,ne),XI3,X,XE,XI,
     '                XP_IB(1,no_ib),ERROR,*9999)
                    DO nj=1,3
                      TEMP_B(nj)=2.d0*X(nj)-XP_IB(nj,no_ib)
                    ENDDO !nj
                  ENDIF !.NOT.KEEP
                ENDIF
                IF(KEEP)THEN
                  DOT=-XNORM(1)*Z(1)-XNORM(2)*Z(2)-XNORM(3)*Z(3)
                  IF(DOT.EQ.1.d0)THEN
                    L_W=RADIUS
                  ELSE
                    L_W=0.5d0*DSIN(PI-2.d0*DACOS(DOT))
     '                *RADIUS/DSIN(DACOS(DOT))
                  ENDIF
                  DO nj=1,3
                    TEMP_B(nj)=X(nj)+XNORM(nj)*L_W
                  ENDDO !nj
                ENDIF
                KEEP2=.TRUE.
                DO nob=no_b+1,no_b+nnob
                  IF(DABS(XP_B(1,nob)-TEMP_B(1)).LE.LOOSE_TOL.AND.
     '              DABS(XP_B(2,nob)-TEMP_B(2)).LE.LOOSE_TOL.AND.
     '              DABS(XP_B(3,nob)-TEMP_B(3))
     '              .LE.LOOSE_TOL) KEEP2=.FALSE.
                ENDDO !nob
                IF(KEEP2)THEN
                  nnob=nnob+1
                  CALL ASSERT(no_b+nnob.LE.NE_R_M,'>>Increase NE_R_M',
     '              ERROR,*9999)
                  DO nj=1,3
                    XP_B(nj,no_b+nnob)=TEMP_B(nj)
                  ENDDO !nj
C                  BDRY(no_b+nnob)=nobdry
                ENDIF !KEEP2

C*** CALCULATE B NODE FOR SECOND FACE
C***    Put host element information into XE2.
                CALL VORO_XPXE(NBJ,ne2,NKJE,NPNE,NVJE,SE,XE2,XP,
     '            ERROR,*9999)
C*** Calculate XI coordinates of IB point in ne2
                CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne2),Z_LEN,USER_TOL,
     '            XE2,XI,XP_IB(1,no_ib),FOUND,.TRUE.,ERROR,*9999)
C*** Find closest orthogonal point on face nf2
                XI(XI1)=XIS(XI1) !on XI1 surface, XI(XID1) will change
                in1=NPF(1,nf2)
                in2=NPF(3,nf2)
                KEEP=.TRUE.
                CALL XINORM(IBT,IDO,INP,in1,in2,NBJ(1,ne2),XI1,
     '            RADIUS,X,XE2,XI,XIS(XI1),XNORM,XP_IB(1,no_ib),KEEP,
     '            ERROR,*9999)
                IF(.NOT.KEEP)THEN !check ne
                  in1=NPF(1,nf)
                  in2=NPF(3,nf)
                  KEEP=.TRUE.
                  CALL XINORM(IBT,IDO,INP,in1,in2,NBJ(1,ne),XI1,RADIUS,
     '              X,XE,XI,XIS(XI1),XNORM,XP_IB(1,no_ib),KEEP,ERROR,
     '              *9999)
                  IF(.NOT.KEEP)THEN
                    CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),USER_TOL,
     '                Z_LEN,XE,XI,XP_IB(1,no_ib),FOUND,.TRUE.,ERROR,
     '                *9999)
                    XI(XI1)=XIS(XI1)
                    CALL XINORM1D(IBT,IDO,INP,NBJ(1,ne),XI3,X,XE,XI,
     '                XP_IB(1,no_ib),ERROR,*9999)
                    DO nj=1,3
                      TEMP_B(nj)=2.d0*X(nj)-XP_IB(nj,no_ib)
                    ENDDO !nj
                  ENDIF !.NOT.KEEP
                ENDIF
                IF(KEEP)THEN
                  DOT=-XNORM(1)*Z(1)-XNORM(2)*Z(2)-XNORM(3)*Z(3)
                  IF(DOT.EQ.1.d0)THEN
                    L_W=RADIUS
                  ELSE
                    L_W=0.5d0*DSIN(PI-2.d0*DACOS(DOT))
     '                *RADIUS/DSIN(DACOS(DOT))
                  ENDIF
                  DO nj=1,3
                    TEMP_B(nj)=X(nj)+XNORM(nj)*L_W
                  ENDDO !nj
                  KEEP2=.TRUE.
                ENDIF !KEEP
                DO nob=no_b+1,no_b+nnob
                  IF(DABS(XP_B(1,nob)-TEMP_B(1)).LE.LOOSE_TOL.AND.
     '              DABS(XP_B(2,nob)-TEMP_B(2)).LE.LOOSE_TOL.AND.
     '              DABS(XP_B(3,nob)-TEMP_B(3))
     '              .LE.LOOSE_TOL) KEEP2=.FALSE.
                ENDDO !nob
                IF(KEEP2)THEN
                  nnob=nnob+1
                  CALL ASSERT(no_b+nnob.LE.NE_R_M,'>>Increase NE_R_M',
     '              ERROR,*9999)
                  DO nj=1,3
                    XP_B(nj,no_b+nnob)=TEMP_B(nj)
                  ENDDO !nj
C                  BDRY(no_b+nnob)=nobdry
                ENDIF !KEEP2
              ENDIF !ne2.EQ.ne
              no_b=no_b+nnob
              XI(XID2)=XI(XID2)+XI_DIST
            ENDDO !WHILE
          ENDIF !nf2.GT.nf
        ENDDO !i
      ENDDO !noface

      DO noface=1,NFLIST(0)
        nf=NFLIST(noface) !global face #
        ne=NPF(6,nf)
        CALL VORO_XPXE(NBJ,ne,NKJE,NPNE,NVJE,SE,XE,XP,ERROR,*9999)
        nef=NPF(8,nf) !local face # in element ne
        XI1=NNF(1,nef,nb) !Xi direction normal to face nf
        IF(nef.EQ.1.OR.nef.EQ.3.OR.nef.EQ.5)THEN
          XIS(XI1)=0.d0 !1st, 3rd, or 5th face -> at Xi=0
          XIF(XI1)=XI_INCR
        ELSE
          XIS(XI1)=1.d0 !2nd, 4th, or 6th face -> at Xi=1
          XIF(XI1)=-XI_INCR
        ENDIF !NPF
C***  Calculate B and IB nodes central to the face
        XI2=NPF(1,nf) !first Xi direction for face
        XI3=NPF(3,nf) !second Xi direction for face
        in1=XI2
        in2=XI3
        XIF(XI2)=0.d0
        XIF(XI3)=0.d0
        XIS(XI2)=XI_DIST !initial XI2 location
        XI_LIMIT(XI2)=1.d0-XI_DIST !final XI2 location
        DO WHILE(XIS(XI2).LE.XI_LIMIT(XI2)) !within the element
          XIS(XI3)=XI_DIST !initial XI3 location
          XI_LIMIT(XI3)=1.d0-XI_DIST !final XI3 location
          DO WHILE(XIS(XI3).LE.XI_LIMIT(XI3)) !within the element
            no_b=no_b+1 !one more B node
            no_ib=no_ib+1 !one more IB node
            CALL ASSERT(no_b.LE.NE_R_M,'>>Increase NE_R_M',ERROR,
     '        *9999)
C***          Calculate the global coordinates and derivatives
C***          for surface point
            DO i=1,3 !3 Xi directions assumed
              XI(i)=XIS(i)
            ENDDO !i
            CALL XNORMXI(IBT,IDO,INP,in1,in2,NBJ(1,ne),X,XE,
     '        XI,XNORM,ERROR,*9999)
            DO nj=1,3
              XP_B(nj,no_b)=X(nj)+XNORM(nj)*RADIUS !store B
              XP_IB(nj,no_ib)=X(nj)-XNORM(nj)*RADIUS !store IB
              XP_I(nj)=X(nj)
            ENDDO !nj
            nobdry=nobdry+1
C            BDRY(no_b)=nobdry
            XIS(XI3)=XIS(XI3)+XI_DIST !increment XI3 location
          ENDDO !j
          XIS(XI2)=XIS(XI2)+XI_DIST !increment XI2 location
        ENDDO !i
      ENDDO !noface (nf)

      N_BDRY=no_b !# of B nodes
      N_IBDRY=no_ib !# of IB nodes
      DO no_b=1,N_BDRY !write the B points
        NVCNP(no_b)=1
      ENDDO !no_b
      DO no_ib=1,N_IBDRY !write the IB points
        NVCNP(N_BDRY+no_ib)=1
      ENDDO !no_ib

      no_in=0 !initialise number of internal nodes
      IF(INTERNAL_NODES)THEN !calculate internal nodes
        IF(REGULAR)THEN !regular Xi spacing
C Place IN nodes inside the host element list, with regular Xi spacing
          N_XI=NDIVISION-1
          DO noelem=1,NELIST(0)
            ne=NELIST(noelem)
            CALL VORO_XPXE(NBJ,ne,NKJE,NPNE,NVJE,SE,XE,XP,ERROR,*9999)
            XI(1)=XI_DIST*2.d0
            DO i=1,N_XI
              XI(2)=XI_DIST*2.d0
              DO j=1,N_XI
                XI(3)=XI_DIST*2.d0
                DO k=1,N_XI
                  DO nj=1,3 !note that 3D elements assumed throughout
                    nb=NBJ(nj,ne) !basis fn
                    X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '                1,XI,XE(1,nj)) !global coordinates of surface point
                  ENDDO !nj
                  no_in=no_in+1
                  XP_IN(1,no_in)=X(1)
                  XP_IN(2,no_in)=X(2)
                  XP_IN(3,no_in)=X(3)
                  NVCNP(no_in+N_BDRY+N_IBDRY)=2
                  XI(3)=XI(3)+XI_DIST*2.d0
                ENDDO !k
                XI(2)=XI(2)+XI_DIST*2.d0
              ENDDO !j
              XI(1)=XI(1)+XI_DIST*2.d0
            ENDDO !i
          ENDDO !noelem
          N_XI=NDIVISION
          DO noelem=1,NELIST(0)
            ne=NELIST(noelem)
            CALL VORO_XPXE(NBJ,ne,NKJE,NPNE,NVJE,SE,XE,XP,ERROR,*9999)
            XI(1)=XI_DIST
            DO i=1,N_XI
              XI(2)=XI_DIST
              DO j=1,N_XI
                XI(3)=XI_DIST
                DO k=1,N_XI
                  DO nj=1,3 !note that 3D elements assumed throughout
                    nb=NBJ(nj,ne) !basis fn
                    X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '                1,XI,XE(1,nj)) !global coordinates of surface point
                  ENDDO !nj
                  no_in=no_in+1
                  XP_IN(1,no_in)=X(1)
                  XP_IN(2,no_in)=X(2)
                  XP_IN(3,no_in)=X(3)
                  NVCNP(no_in+N_BDRY+N_IBDRY)=2
                  XI(3)=XI(3)+XI_DIST*2.d0
                ENDDO !k
                XI(2)=XI(2)+XI_DIST*2.d0
              ENDDO !j
              XI(1)=XI(1)+XI_DIST*2.d0
            ENDDO !i
          ENDDO !noelem
          
        ELSE !alveolar spacing
C IN nodes are placed along the centrelines of
          no_in=0
          n_al_ax=4 !# of alveoli axially
          IF(NEELEM(0,nr_tree).EQ.1) THEN
            b_inc=0.784d0/n_al_ax !for single branch
          ELSE
            b_inc=1.d0/n_al_ax
          ENDIF
          ne=NEELEM(1,nr_tree) !first point in region
          nb=NBJ(1,ne) !1D basis function
          np=NPNE(1,nb,ne) !first node in region nr_target
C         no_in=no_in+1
C         DO nj=1,NJT
C         XP_IN(nj,no_in)=XP(1,1,nj,np)
C         XI_HALF(nj)=0.d0
C         ENDDO !nj
C         NVCNP(no_in+N_BDRY+N_IBDRY)=0
          DO noelem=1,NEELEM(0,nr_tree) !for current airway paths in nr_target
            ne=NEELEM(noelem,nr_tree)
            CALL VORO_XPXE(NBJ,ne,NKJE,NPNE,NVJE,SE,XE,XP,ERROR,*9999)
            DO nj=1,3
              XI(nj)=0.d0
            ENDDO !nj
            np=NPNE(1,nb,ne) !first node in element ne
            DO i=1,n_al_ax
              NUM_ALVEOLI1=5
              NUM_ALVEOLI2=20
              R0=0.0935d0
              R=0.09351d0
              IF(NEELEM(0,nr_tree).EQ.1) THEN !for single branch
                R2=0.4725d0
              ELSE
                R2=0.15d0
              ENDIF
              ANG_INCREMENT1=2.d0*PI/NUM_ALVEOLI1
              ANG_INCREMENT2=2.d0*PI/NUM_ALVEOLI2
              XI(1)=XI(1)+b_inc !position along ne of seed point
              no_in=no_in+1
              DO nj=1,3
                nb=NBJ(nj,ne) !basis fn, must be 1D, check with ASSERT
                XP_IN(nj,no_in)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '            INP(1,1,nb),nb,1,XI,XE(1,nj)) !coords of airway centre
              ENDDO !nj
              NVCNP(no_in+N_BDRY+N_IBDRY)=0 !records as airway centre
              XI_HALF(1)=XI(1)-0.5d0*b_inc !axial position of outer seed
              DO nj=1,3
                nb=NBJ(nj,ne) !basis fn, must be 1D, check with ASSERT
                B(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,1,XI_HALF,XE(1,nj)) !translation to origin
                X1(nj)=XP(1,1,nj,np)-B(nj)
              ENDDO !nj
              IF(X1(1).NE.0.d0)THEN
                ANG1=DATAN(X1(2)/X1(1))
                ANG2=DATAN(X1(3)/X1(1))
              ELSE
                ANG1=0.d0
                ANG2=0.d0
              ENDIF

              ROTATE=0.d0
              DO WHILE(ROTATE.LT.2.d0*PI-LOOSE_TOL)
                no_in=no_in+1
                U(1)=R0*DCOS(ROTATE)
                U(2)=R0*DSIN(ROTATE)
                U(3)=0.d0
                CALL COORD_TRANS(ANG2,ANG1,0.d0,XP_IN(1,no_in),
     '            B,U,ERROR,*9999)
                NVCNP(no_in+N_BDRY+N_IBDRY)=1
                ROTATE=ROTATE+ANG_INCREMENT1
              ENDDO

              XI_HALF(1)=XI(1)-0.5d0*b_inc !axial position of outer seed
C             XI_HALF(1)=XI(1)
              DO nj=1,3
                nb=NBJ(nj,ne) !basis fn, must be 1D, check with ASSERT
                B(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,1,XI_HALF,XE(1,nj)) !translation to origin
                X1(nj)=XP(1,1,nj,np)-B(nj)
              ENDDO !nj
              IF(X1(1).NE.0.d0)THEN
                ANG1=DATAN(X1(2)/X1(1))
                ANG2=DATAN(X1(3)/X1(1))
              ELSE
                ANG1=0.d0
                ANG2=0.d0
              ENDIF

              ROTATE=0.d0
              DO WHILE(ROTATE.LT.2.d0*PI-LOOSE_TOL)
                no_in=no_in+1
                U(1)=R*DCOS(ROTATE)
                U(2)=R*DSIN(ROTATE)
                U(3)=0.d0
                CALL COORD_TRANS(ANG2,ANG1,0.d0,XP_IN(1,no_in),
     '            B,U,ERROR,*9999)
                NVCNP(no_in+N_BDRY+N_IBDRY)=2
                ROTATE=ROTATE+ANG_INCREMENT1
              ENDDO
              XI_HALF(1)=XI(1)-0.5d0*b_inc
              DO nj=1,3
                nb=NBJ(nj,ne) !basis fn, must be 1D, check with ASSERT
                B(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,1,XI_HALF,XE(1,nj)) !translation to origin
                X1(nj)=XP(1,1,nj,np)-B(nj)
              ENDDO !nj
              IF(X1(1).NE.0.d0)THEN
                ANG1=DATAN(X1(2)/X1(1))
                ANG2=DATAN(X1(3)/X1(1))
              ELSE
                ANG1=0.d0
                ANG2=0.d0
              ENDIF

              ROTATE=0.d0
c             ROTATE=0.5d0*ANG_INCREMENT
c             c            ROTATE=0.25d0*ANG_INCREMENT
              no_in2=NUM_ALVEOLI2/NUM_ALVEOLI1+1
C             no_in3=NUM_ALVEOLI2/NUM_ALVEOLI1-1
              DO WHILE(ROTATE.LT.2.d0*PI-LOOSE_TOL)
                no_in=no_in+1
                no_in2=no_in2+1
C               no_in3=no_in3+1
                IF(MOD(no_in2,NUM_ALVEOLI2/NUM_ALVEOLI1).EQ.0) THEN
                  R2_NEW=0.95d0*R2
C                 ELSE IF(MOD(no_in3,NUM_ALVEOLI2/NUM_ALVEOLI1).EQ.0)
C                 THEN
C                 R2_NEW=1.05d0*R2
                ELSE
                  R2_NEW=R2
                ENDIF
                U(1)=R2_NEW*DCOS(ROTATE)
                U(2)=R2_NEW*DSIN(ROTATE)
                U(3)=0.d0
                CALL COORD_TRANS(ANG2,ANG1,0.d0,XP_IN(1,no_in),
     '            B,U,ERROR,*9999)
                NVCNP(no_in+N_BDRY+N_IBDRY)=1
                ROTATE=ROTATE+ANG_INCREMENT2
c               c              ROTATE=ROTATE+0.5d0*ANG_INCREMENT
              ENDDO
            ENDDO !i

C... KSB: To create end alveoli
            IF(NEELEM(0,nr_tree).EQ.1) THEN !currently do only for single branch
              NUM_ALVEOLI2=3
              NUM_ALVEOLI1=20
              ANG_INCREMENT2=2.d0*PI/NUM_ALVEOLI2
              ANG_INCREMENT1=PI/3.7d0 !3.5d0 ! angle in z plane
              DO i=1,3 !3 rows of points
                IF(i.EQ.1) THEN
                  R2_NEW=0.0d0
                  BOUND=1
                ELSE IF(i.EQ.2) THEN
                  R2_NEW=R*0.93d0
                  BOUND=2
                ELSE IF(i.EQ.3) THEN
                  R2_NEW=R2*0.6d0 !0.65d0
                  BOUND=1
                  ANG_INCREMENT2=2.d0*PI/NUM_ALVEOLI1
                ENDIF
                DO nj=1,3
                  nb=NBJ(nj,ne)
                  B(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '              nb,1,XI,XE(1,nj)) !translation to origin
                  X1(nj)=XP(1,1,nj,np)-B(nj)
                ENDDO !nj
                IF(X1(1).NE.0.d0)THEN
                  ANG1=DATAN(X1(2)/X1(1))
                  ANG2=DATAN(X1(3)/X1(1))
                ELSE
                  ANG1=0.d0
                  ANG2=0.d0
                ENDIF

                R_NEW=R2_NEW*DSIN(ANG_INCREMENT1)
                ROTATE=0.d0
                DO WHILE(ROTATE.LT.2.d0*PI-LOOSE_TOL)
                  no_in=no_in+1
                  U(1)=R_NEW*DCOS(ROTATE)
                  U(2)=R_NEW*DSIN(ROTATE)
                  U(3)=R2_NEW*DCOS(ANG_INCREMENT1)
                  CALL COORD_TRANS(ANG2,ANG1,0.d0,XP_IN(1,no_in),
     '              B,U,ERROR,*9999)
                  NVCNP(no_in+N_BDRY+N_IBDRY)=BOUND
                  ROTATE=ROTATE+ANG_INCREMENT2
                ENDDO
              ENDDO

              XI(1)=XI(1)+b_inc
              no_in=no_in+1
              DO nj=1,3
                nb=NBJ(nj,ne) !basis fn, must be 1D, check with ASSERT
                XP_IN(nj,no_in)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '            INP(1,1,nb),nb,1,XI,XE(1,nj)) !coords of airway centre
              ENDDO !nj
              NVCNP(no_in+N_BDRY+N_IBDRY)=1
            ENDIF
C... end creating end alveoli
          ENDDO !noelem
        ENDIF
      ENDIF !INTERNAL_NODES
      N_INTNL=no_in !total number of internal nodes

      nonode=0
      np=NPT(0)
C Put XP_B, XP_IB, XP_IN into XP
      DO no_b=1,N_BDRY
        nonode=nonode+1
        np=np+1
        CALL ASSERT(nonode.LE.NPM,'>>Increase NPM',ERROR,*9999)
        CALL ASSERT(np.LE.NP_R_M,'>>Increase NP_R_M',ERROR,*9999)
        NPNODE(nonode,nr_target)=np
        NPLIST(no_b)=np
        DO nj=1,NJT
          XP(1,1,nj,np)=XP_B(nj,no_b)
          NKJ(nj,np)=1 !# derivatives for nj at np
          NVJP(nj,np)=1 !# of versions of nj at np
        ENDDO !j
      ENDDO !N_BDRY
      NPLIST(0)=N_BDRY
      STRING='BOUNDARY'
      CALL GRNODE_SUB(NPLIST,STRING,.TRUE.,ERROR,*9999)

      DO no_ib=1,N_IBDRY
        nonode=nonode+1
        np=np+1
        NPNODE(nonode,nr_target)=np
        NPLIST(no_ib)=np
        DO nj=1,NJT
          XP(1,1,nj,np)=XP_IB(nj,no_ib)
          NKJ(nj,np)=1 !# derivatives for nj at np
          NVJP(nj,np)=1 !# of versions of nj at np
        ENDDO !j
      ENDDO !N_IBDRY
      NPLIST(0)=N_IBDRY
      STRING='INTERNAL_BOUNDARY'
      CALL GRNODE_SUB(NPLIST,STRING,.TRUE.,ERROR,*9999)

      DO no_in=1,N_INTNL
        nonode=nonode+1
        np=np+1
        NPNODE(nonode,nr_target)=np
        NPLIST(no_in)=np
        DO nj=1,NJT
          XP(1,1,nj,np)=XP_IN(nj,no_in)
          NKJ(nj,np)=1 !# derivatives for nj at np
          NVJP(nj,np)=1 !# of versions of nj at np
        ENDDO !j
      ENDDO !N_BDRY
      NPLIST(0)=N_INTNL
      STRING='INTERNAL'
      IF(NPLIST(0).NE.0)THEN
        CALL GRNODE_SUB(NPLIST,STRING,.TRUE.,ERROR,*9999)
      ENDIF

      NPNODE(0,nr_target)=N_BDRY+N_IBDRY+N_INTNL
      NPNODE(0,0)=NPNODE(0,0)+NPNODE(0,nr_target)
      NPT(nr_target)=np
      NPT(0)=np
      CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)

C Put B, IB, and IN nodes into node groups
c      DO no_b=1,N_BDRY !write the B points
c        NPLIST(no_b)=NPT(0)+no_b
c      ENDDO !no_b
c      NPLIST(0)=N_BDRY
c      STRING='BOUNDARY'
c      CALL GRNODE_SUB(NPLIST,STRING,.TRUE.,ERROR,*9999)
c      DO no_ib=1,N_IBDRY !write the IB points
c        nobdry=N_BDRY+no_ib
c        NPLIST(no_ib)=NPT(0)+N_BDRY+no_ib
c      ENDDO !no_ib
c      NPLIST(0)=N_IBDRY
c      STRING='INTERNAL_BOUNDARY'
c      CALL GRNODE_SUB(NPLIST,STRING,.TRUE.,ERROR,*9999)
c      DO no_in=1,N_INTNL !write the IN points
c        nobdry=N_BDRY+N_IBDRY+no_in
c        NPLIST(no_in)=NPT(0)+N_BDRY+N_IBDRY+no_in
c      ENDDO !no_in
c      NPLIST(0)=N_INTNL
c      STRING='INTERNAL'
c      IF(NPLIST(0).NE.0)THEN
c        CALL GRNODE_SUB(NPLIST,STRING,.TRUE.,ERROR,*9999)
c      ENDIF

      CALL EXITS('DELAUNAY_NODES')
      RETURN
 9999 CALL ERRORS('DELAUNAY_NODES',ERROR)
      CALL EXITS('DELAUNAY_NODES')
      RETURN 1
      END


