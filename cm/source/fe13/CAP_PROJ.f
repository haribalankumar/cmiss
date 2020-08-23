      SUBROUTINE CAP_PROJ(group,NBJ,NELIST,NE_REPEAT,NKJ,
     '  NPNE,NPNE_ALV,NPLIST,NPNODE,nr,NVJP,centre,XP,PROJECT1,
     '  ERROR,*)

C#### Subroutine: CAP_PROJ
C###  Description:
C###    CAP_PROJ projects the pulmonary capillary mesh created on the
C###    surface of a unit sphere onto the surface of a pulmonary
C###    alveolar mesh. Delaunay elements from this alveolar mesh are
C###    grouped to each individual alveolus and the spherical capillary
C###    mesh projected onto the surface.

C*** Created by Kelly Burrowes, May 2002.
C*** Last modified 10/10/02 - still under development

C... Equation of plane for node, np1, is (x-XP(np1)).n=0
C... Eqn of line from centre of sphere (c) to np1 is x=c+tXP(np1)
C... point of intersection (x) of these 2 gives distance from c = t

C... array NPNE_ALV(0:NP_NE,NE_R_M) is used to store which
C... capillary nodes are projected onto which alveolar element.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'genmesh.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'pulm00.cmn'
!     Parameter list
      INTEGER group,NBJ(NJM,NEM),NELIST(0:NEM),
     '  NE_REPEAT(0:NEREPM),NKJ(NJM,NPM),NPNE(NNM,NBFM,NEM),
     '  NPNE_ALV(0:NP_NE,NE_R_M),NPLIST(0:NPM),NPNODE(0:NP_R_M,0:NRM),
     '  nr,NVJP(NJM,NPM)
      REAL*8 centre(NJT),XP(NKM,NVM,NJM,NPM)
      LOGICAL PROJECT1
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER i,nb2,ne,ne_i,nj,nolist,nonode,nonode2,np,np1,
     '  np2,np3,np4,np5
      REAL*8 A(NJT),AREA,B(NJT),C(NJT),D(NJT),DOT_PROD,E(NJT),F(NJT),
     '  N(NJT),P_A1(NJT),P_A2(NJT),P_A3(NJT),P_A4(NJT),P_21(NJT),
     '  P_21_CHECK(NJT),P_24(NJT),P_31(NJT),
     '  P_31_CHECK(NJT),SUM,t,tol,XP_np1(NJT),
     '  XP_point(NJT),XP_temp(NJT),X_TRANS(NJT)
      LOGICAL FOUND,NE_REPEATED
      
      CALL ENTERS('CAP_PROJ',*9999)

      np5=NPT(0)
      nonode2=NPNODE(0,nr)
      tol=1.d-2
C... translates XP co-ords of nodes being projected from slave mesh
C... by co-ords centre(x,y,z), input in IPMESH2 & rotates points

      IF(VORONOI_ALV) THEN !if projecting onto voronoi alveoli
        ANGLE_X=ALV_ANGLE(group,1)
        ANGLE_Y=ALV_ANGLE(group,2)
        ANGLE_Z=ALV_ANGLE(group,3)
        DO nonode=1,NPLIST(0)
          np=NPLIST(nonode)
          DO nj=1,NJT
            XP_point(nj)=XP(1,1,nj,np)
          ENDDO
          CALL COORD_TRANS(ANGLE_X,ANGLE_Y,ANGLE_Z,X_TRANS,
     '      centre,XP_point,ERROR,*9999)
          DO nj=1,NJT
            XP(1,1,nj,np)=X_TRANS(nj)
          ENDDO
        ENDDO
      ELSE
        DO nonode=1,NPLIST(0)
          np=NPLIST(nonode)
          DO nj=1,NJT !adds offset, will be subtracted at end to retain
            XP(1,1,nj,np)=XP(1,1,nj,np)+centre(nj) !sphere centre:0,0,0
          ENDDO
        ENDDO
      ENDIF
      DO nonode=1,NPLIST(0)
        np=NPLIST(nonode)
        nolist=1
        FOUND=.FALSE.
        DO WHILE(.NOT.FOUND.AND.nolist.LE.NELIST(0))
          ne=NELIST(nolist) !temporarily using NELIST
          NE_REPEATED=.FALSE. !initialise for ne
          nb2=NBJ(1,ne) !let nj=1
          np1=NPNE(2,nb2,ne) !3 nodes for an element
          np2=NPNE(3,nb2,ne) !only 3 nodes required to calculate plane
          np3=NPNE(4,nb2,ne)
          np4=NPNE(1,nb2,ne)
          DO nj=1,NJT !each direction
            P_21(nj)=XP(1,1,nj,np2)-XP(1,1,nj,np1)
            P_31(nj)=XP(1,1,nj,np3)-XP(1,1,nj,np1)
            P_21_CHECK(nj)=XP(1,1,nj,np2)-XP(1,1,nj,np1)
            P_31_CHECK(nj)=XP(1,1,nj,np3)-XP(1,1,nj,np1)
            IF(np1.NE.np2) THEN !not triangular element
              P_24(nj)=XP(1,1,nj,np2)-XP(1,1,nj,np4)
            ENDIF
          ENDDO !nj
          CALL CROSS(P_21,P_31,N) !calculates normal vector (N)
          DO nj=1,NJT
            XP_np1(nj)=XP(1,1,nj,np1) !first vertex of element
            XP_temp(nj)=XP(1,1,nj,np)-centre(nj) !node being projected
          ENDDO
          IF(DABS(N(1)).LT.tol.AND.DABS(N(2)).LT.tol.AND.
     '        DABS(N(3)).LT.tol) THEN
            IF(np1.NE.np4) THEN
              CALL NORMALISE(NJT,P_21_CHECK,ERROR,*9999)
              CALL NORMALISE(NJT,P_31_CHECK,ERROR,*9999)
              IF(DABS(DACOS(DOT_PROD(P_21_CHECK,P_31_CHECK))).LT.tol)
     '          CALL CROSS(P_21,P_24,N)
            ENDIF
            IF(N(1).EQ.0.d0.AND.N(2).EQ.0.d0.AND.N(3).EQ.0.d0) THEN
              t=0.d0
            ELSE
              t=(DOT_PROD(XP_np1,N)-DOT_PROD(centre,N))/
     '          DOT_PROD(XP_temp,N)
            ENDIF
          ELSE !triangular
            t=(DOT_PROD(XP_np1,N)-DOT_PROD(centre,N))/
     '        DOT_PROD(XP_temp,N)
          ENDIF
          IF(t.GT.0.d0.AND.t.LE.1.d0) THEN
C... plane/line don't intersect face
            DO nj=1,NJT    !new position of np on
              XP_temp(nj)=centre(nj)+t*XP_temp(nj) !alveolar face
            ENDDO !nj
C... check that point is actually on face for each pair of vertices
            DO nj=1,NJT
              P_A1(nj)=XP_temp(nj)-XP(1,1,nj,np1) !vector from np1 to A
              P_A2(nj)=XP_temp(nj)-XP(1,1,nj,np2)
              P_A3(nj)=XP_temp(nj)-XP(1,1,nj,np3)
              P_A4(nj)=XP_temp(nj)-XP(1,1,nj,np4)
            ENDDO
            IF(np1.EQ.np4) THEN !triangular element
              CALL CROSS(P_21,P_31,A)
              CALL CROSS(P_A1,P_A2,B)
              CALL CROSS(P_A2,P_A3,D)
              CALL CROSS(P_A3,P_A1,E)
            ELSE
              CALL CROSS(P_21,P_31,A)
              CALL CROSS(P_21,P_24,B)
              CALL CROSS(P_A1,P_A3,C)
              CALL CROSS(P_A3,P_A2,D)
              CALL CROSS(P_A2,P_A4,E)
              CALL CROSS(P_A4,P_A1,F)
            ENDIF
            AREA=0.d0
            SUM=0.d0
            IF(np1.EQ.np4) THEN
              AREA=0.5d0*DSQRT(A(1)**2.d0+A(2)**2.d0+A(3)**2.d0)
            ELSE
              AREA=0.5d0*(DSQRT(A(1)**2.d0+A(2)**2.d0+A(3)**2.d0)
     '          +DSQRT(B(1)**2.d0+B(2)**2.d0+B(3)**2.d0))
            ENDIF
            IF(np1.EQ.np4) THEN
              SUM=0.5d0*(DSQRT(B(1)**2.d0+B(2)**2.d0+B(3)**2.d0)+
     '          DSQRT(D(1)**2.d0+D(2)**2.d0+D(3)**2.d0)+
     '          DSQRT(E(1)**2.d0+E(2)**2.d0+E(3)**2.d0))
            ELSE
              SUM=0.5d0*(DSQRT(C(1)**2.d0+C(2)**2.d0+C(3)**2.d0)
     '          +DSQRT(D(1)**2.d0+D(2)**2.d0+D(3)**2.d0)+DSQRT(E(1)
     '          **2.d0+E(2)**2.d0+E(3)**2.d0)+DSQRT(F(1)**2.d0+F(2)
     '          **2.d0+F(3)**2.d0))
            ENDIF
            IF(SUM/AREA.GT.1.d0-tol.AND.
     '        SUM/AREA.LT.1.d0+tol) THEN
              FOUND=.TRUE. !still TRUE but node not created
C... check if this alveolar element is being repeated & nodes therefore
C... already projected onto this face (element), i.e new nodes not made
              DO i=1,NE_REPEAT(0)
                ne_i=NE_REPEAT(i)
                IF(ne_i.EQ.ne) NE_REPEATED=.TRUE.
              ENDDO !i
              IF(.NOT.NE_REPEATED.AND.PROJECT1) THEN
C... else node not created, if 2nd time through nodes already created
C... just need to change XP co-ordinates
                np5=np5+1
                nonode2=nonode2+1 !create new node
                NPNODE(nonode2,nr)=np5
                NPNE_ALV(0,ne)=NPNE_ALV(0,ne)+1
                NPNE_ALV(NPNE_ALV(0,ne),ne)=np5
                DO nj=1,NJT
                  XP(1,1,nj,np5)=XP_temp(nj) !stores new nodal co-ords
                  NKJ(nj,np5)=1 !# derivatives - 1D
                  NVJP(nj,np5)=1 !# versions
                ENDDO
              ELSE IF(.NOT.PROJECT1) THEN
                DO nj=1,NJT
                  XP(1,1,nj,np)=XP_temp(nj) !changes nodal co-ords
                ENDDO
              ENDIF
            ENDIF
          ENDIF !t.GT.0
          nolist=nolist+1
        ENDDO !WHILE
        IF(.NOT.FOUND.AND.nolist.GE.NELIST(0)) THEN
          WRITE(OP_STRING,'('' WARNING:node not projected:'',I5)') np
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDDO !nonode
      NPNODE(0,nr)=nonode2
      NPNODE(0,0)=np5
      NPT(0)=np5
      NPT(nr)=nonode2
      IF(VORONOI_ALV) THEN
C... undo translation of points, finds inverse of coordinate transform
        IF(PROJECT1) THEN
          DO nonode=1,NPLIST(0)
            np=NPLIST(nonode)
            DO nj=1,NJT
              XP_point(nj)=XP(1,1,nj,np)
            ENDDO
            CALL COORD_TRANS_INV(ANGLE_X,ANGLE_Y,ANGLE_Z,X_TRANS,
     '        centre,XP_point,ERROR,*9999)
            DO nj=1,NJT
              XP(1,1,nj,np)=X_TRANS(nj)
            ENDDO
          ENDDO !nonode
        ENDIF !PROJECT1
      ELSE
        IF(PROJECT1) THEN !leave centre offset if 2nd projection
          DO nonode=1,NPLIST(0) !because these are the actual mesh nodes
            np=NPLIST(nonode)
            DO nj=1,NJT !subtracts offset to retain original centre(000)
              XP(1,1,nj,np)=XP(1,1,nj,np)-centre(nj)
            ENDDO
          ENDDO
        ENDIF
      ENDIF !VORONOI_ALV

      CALL EXITS('CAP_PROJ')
      RETURN
 9999 CALL ERRORS('CAP_PROJ',ERROR)
      CALL EXITS('CAP_PROJ')
      RETURN 1
      END

