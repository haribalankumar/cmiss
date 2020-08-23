      SUBROUTINE GNNEZD(IBT,IDO,INP,LD,NBJ,NDDATA,NELIST,
     &  NKJE,NPF,NPNE,nr,NVJE,PG,RG,SE,WD,WG,XA,X_DISCRETISATION,
     &  XE,XG,XID,XP,ZD,SPREAD_TYPE,ERROR,*)

C#### Subroutine: GNNEZD
C###  Description:
C###    GNNEZD generates points (coordinates stored in data array ZD)
C###    into a host element.  The points are uniformly spaced in either
C###    xi space or global coordinates.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     &  LD(NDM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),
     &  NELIST(0:NEM),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     &  NPNE(NNM,NBFM,NEM),nr,NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),
     &  WD(NJM,NDM),WG(NGM,NBM),XA(NAM,NJM,NEM),
     &  X_DISCRETISATION,XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),
     &  XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      CHARACTER SPREAD_TYPE*7,ERROR*(*)
!     Local Variables
      INTEGER i,IT,j,k,m,Mii,nb,nbb,ne,nd,ndd,nd_start,NDT2,ng,
     &  ni,n,Nii,NITB,nj,noelem,
     &  nonode,np,number_of_points_xi,p,Pii
      REAL*8 DET,DXRCXI(3,3),MAX(3),MAX_NE(400,3),MAX_NE2(3),MIN(3),
     &  MIN_NE(400,3),MIN_NE2(3),PXI,START_X(3),USER_TOL,VOL,X(3),
     &  XI(3),XI_DISCRETISATION,Z_LEN
      CHARACTER STRING*(155)
      LOGICAL FOUND,FOUND_NE,IDENTICAL,INELEM,POSITIVE

      CALL ENTERS('GNNEZD',*9999)

      nb=NBJ(1,NELIST(1)) !basis function number for 1st element
      NDT2=0
      NDT=0

C MHT 11.02.11 changing the process for generating the grid of points
C for the purpose of easily being able to determine which points are adjacent
C to each other. The main change is that the points are generated to fill a
C rectangular box that encloses the entire domain, and then they are 'pruned'
C by checking whether they are inside or outside the volume mesh.

      IF(SPREAD_TYPE(1:7).EQ.'REGULAR')THEN !regular spacing in global coordinates
        USER_TOL=LOOSE_TOL !tolerance for calculating Xi coordinate
        NDT=0 !count the number of data points
        ne=NELIST(1) !first element in list
        np=NPNE(1,nb,ne) !first node in element (at Xi1,Xi2,Xi3 = 0,0,0 by definition)

C find the bounding coordinates for the volume mesh (total and each element) 
C using max and min coordinates for each element

        DO nj=1,NJT
          MIN(nj)=XP(1,1,nj,np) !initial minimum coordinates
          MAX(nj)=XP(1,1,nj,np) !initial maximum coordinates
        ENDDO !nj
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          WRITE(STRING,'('' Increase array size in GNNEZD'')')
          CALL ASSERT(ne.LE.400,STRING,ERROR,*9999)
          DO nj=1,NJT !initial local min and max
            MIN_NE(ne,nj)=XP(1,1,nj,NPNE(1,nb,ne))
            MAX_NE(ne,nj)=XP(1,1,nj,NPNE(1,nb,ne))
          ENDDO !nj
          DO nonode=1,NNT(nb)
            np=NPNE(nonode,nb,ne)
            DO nj=1,NJT
              IF(XP(1,1,nj,np).LT.MIN(nj)) MIN(nj)=XP(1,1,nj,np)
              IF(XP(1,1,nj,np).GT.MAX(nj)) MAX(nj)=XP(1,1,nj,np)
              IF(XP(1,1,nj,np).LT.MIN_NE(ne,nj)) 
     &             MIN_NE(ne,nj)=XP(1,1,nj,np)
              IF(XP(1,1,nj,np).GT.MAX_NE(ne,nj)) 
     &             MAX_NE(ne,nj)=XP(1,1,nj,np)
            ENDDO !nj
          ENDDO !nonode
        ENDDO !noelem

C Set the x,y,z coordinates for the first data point, at the minimum position
C for the bounding box, offset by the data point spacing

        DO nj=1,NJT
          START_X(nj)=MIN(nj)-0.5d0*X_DISCRETISATION
        ENDDO !nj

C create a grid of data points that fills the entire bounding box
C use the numbering pattern to record data point 'connectivity' to other points

        Mii = INT((MAX(1)-MIN(1))/X_DISCRETISATION)+1
        Nii = INT((MAX(2)-MIN(2))/X_DISCRETISATION)+1
        Pii = INT((MAX(3)-MIN(3))/X_DISCRETISATION)+1

C check that the data arrays will be large enough to store the grid

        WRITE(STRING,'('' Increase NDM in .ippara file to'',
     &      I6)') Mii*Nii*Pii
        CALL ASSERT(NDM.GT.Mii*Nii*Pii,STRING,ERROR,*9999)

        X(1)=START_X(1)
        DO m=1,Mii
          X(1)=X(1)+X_DISCRETISATION
          X(2)=START_X(2)
          DO n=1,Nii
            X(2)=X(2)+X_DISCRETISATION
            X(3)=START_X(3)
            DO p=1,Pii
              NDT=NDT+1
C record the connectivity pattern
c              NDADJ(1,NDT)=NDT-1
c              NDADJ(2,NDT)=NDT+1
c              NDADJ(3,NDT)=NDT-Pii
c              NDADJ(4,NDT)=NDT+Pii
c              NDADJ(5,NDT)=NDT-Pii*Nii
c              NDADJ(6,NDT)=NDT+Pii*Nii
c              IF(p.EQ.1)   NDADJ(1,NDT)=0
c              IF(p.EQ.Pii) NDADJ(2,NDT)=0
c              IF(n.EQ.1)   NDADJ(3,NDT)=0
c              IF(n.EQ.Nii) NDADJ(4,NDT)=0
c              IF(m.EQ.1)   NDADJ(5,NDT)=0
c              IF(m.EQ.Mii) NDADJ(6,NDT)=0

              X(3)=X(3)+X_DISCRETISATION
              DO nj=1,NJT
                ZD(nj,NDT)=X(nj)
              ENDDO !nj
            ENDDO               !Pii
          ENDDO                 !Nii
        ENDDO                   !Mii

C prune the data points by calculating whether they are inside or outside
C the volume. 
        ndd=0
        NDDATA(0,nr)=0
        DO nd=1,NDT
          FOUND_NE=.FALSE. !logical for whether the host element has been found
          DO noelem=1,NELIST(0)
            ne=NELIST(noelem)
            IF(.NOT.FOUND_NE)THEN !only do this if not found already
C check whether data point is in the element only if it is determined
C to be inside the bounding box

              IF(ZD(1,nd).GE.MIN_NE(ne,1).AND.
     &          ZD(1,nd).LE.MAX_NE(ne,1).AND.
     &          ZD(2,nd).GE.MIN_NE(ne,2).AND.
     &          ZD(2,nd).LE.MAX_NE(ne,2).AND.
     &          ZD(3,nd).GE.MIN_NE(ne,3).AND.
     &          ZD(3,nd).LE.MAX_NE(ne,3))THEN !is inside bounding box
C transfer the nodal coordinates to an element array
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     &            NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     &            XA(1,1,ne),XE,XP,ERROR,*9999)
C initialise variables for calculating XI point
                FOUND=.FALSE.
                IT=0
                DO ni=1,3
                  XI(ni)=0.5d0
                ENDDO
C calculate Xi coordinates. if 0<=Xi<=1, then is inside element
                CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),Z_LEN,
     &            USER_TOL,XE,XI,ZD(1,nd),FOUND,.FALSE.,ERROR,*9999)
                IF(FOUND)THEN
                  ndd=ndd+1
                  NDDATA(ndd,nr)=ndd
                  NDDATA(0,nr)=NDDATA(0,nr)+1
                  LD(ndd)=ne           !store the host element number
                  DO nj=1,NJT
                    WD(nj,ndd)=1.0d0   !set the weighting
                    ZD(nj,ndd)=ZD(nj,nd)
                    XID(nj,ndd)=XI(nj) !store the element Xi location
                  ENDDO !nj
c                  NMAP(ndd)=nd
                  FOUND_NE=.TRUE.
                ENDIF !FOUND
              ENDIF !ZD
            ENDIF !FOUND_NE
          ENDDO !noelem

        ENDDO !NDT

        NDT=NDDATA(0,nr) !number of data points remaining

C update the adjacency array to reflect mapping of nd to ndd
C NDADJ already contains the adjacency for the initial cuboid. This can
C be overwritten safely because the data points are numbered sequentially. 
C i.e. won't be overwriting the wrong parts of the array

c        DO nd=1,NDT !for each of the new data points (the subset)
c          ndd=NDMAP(nd) !data number in cuboid
c          DO i=1,6
c            n_adj=NADJ(i,ndd,)
c            NDADJ(i,nd)=NDMAP(n_adj)
c          ENDDO
c        ENDDO !nd

      ELSE IF(SPREAD_TYPE(1:7).EQ.'BYELMNT')THEN !regular spacing in global coordinates
        USER_TOL=LOOSE_TOL
        NDT=0
        ne=NELIST(1)
        np=NPNE(1,nb,ne)
        DO nj=1,NJT
          MIN(nj)=XP(1,1,nj,np) !initial minimum coordinates
          MAX(nj)=XP(1,1,nj,np) !initial maximum coordinates
        ENDDO !nj
 !find the host volume bounding coordinates
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          DO nonode=1,NNT(nb)
            np=NPNE(nonode,nb,ne)
            DO nj=1,NJT
              IF(XP(1,1,nj,np).LT.MIN(nj)) MIN(nj)=XP(1,1,nj,np)
              IF(XP(1,1,nj,np).GT.MAX(nj)) MAX(nj)=XP(1,1,nj,np)
            ENDDO !nj
          ENDDO !nonode
        ENDDO !noelem

        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          NDT2=0 !number of seed points in element
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '      XA(1,1,ne),XE,XP,ERROR,*9999)

          DO nj=1,NJT !initial local min and max
            MIN_NE2(nj)=XP(1,1,nj,NPNE(1,nb,ne))
            MAX_NE2(nj)=XP(1,1,nj,NPNE(1,nb,ne))
          ENDDO !nj
          
          DO nonode=1,NNT(nb)
            np=NPNE(nonode,nb,ne)
            DO nj=1,NJT
              IF(XP(1,1,nj,np).LT.MIN_NE2(nj)) MIN_NE2(nj)=XP(1,1,nj,np)
              IF(XP(1,1,nj,np).GT.MAX_NE2(nj)) MAX_NE2(nj)=XP(1,1,nj,np)
            ENDDO !nj
          ENDDO !nonode

          DO nj=1,NJT
            START_X(nj)=MIN(nj)+DBLE(X_DISCRETISATION*INT((MIN_NE2(nj)
     '        -MIN(nj))/X_DISCRETISATION+0.5d0))+0.5d0*X_DISCRETISATION
          ENDDO !nj
        
          X(1)=START_X(1)
          DO WHILE(X(1).LE.MAX_NE2(1))
            X(2)=START_X(2)
            DO WHILE(X(2).LE.MAX_NE2(2))
              X(3)=START_X(3)
              DO WHILE(X(3).LE.MAX_NE2(3))
                FOUND=.FALSE.
                DO ni=1,3
                  XI(ni)=0.5d0
                ENDDO
                INELEM=.FALSE.
                IT=0
                CALL CLOS31(IBT,IDO,INP,IT,20,NBJ(1,ne),Z_LEN,
     '            USER_TOL,XE,XI,X,FOUND,INELEM,ERROR,*9999)
                IF(FOUND)THEN
                  NDT=NDT+1 !counts RPs for whole lobe
                  NDT2=NDT2+1
                  CALL ASSERT(NDM.GT.NDT,
     &              '>>Increase NDM in .ippara file',ERROR,*9999)
                  LD(NDT)=ne
                  DO nj=1,NJT
                    ZD(nj,NDT)=X(nj)
                    WD(nj,NDT)=1.0d0
                    XID(nj,NDT)=XI(nj) !store xi location
                  ENDDO !nj
                  NDDATA(NDT,nr)=NDT
                  NDDATA(0,nr)=NDDATA(0,nr)+1
                ENDIF !FOUND
                X(3)=X(3)+X_DISCRETISATION
              ENDDO              !WHILE
              X(2)=X(2)+X_DISCRETISATION
            ENDDO                 !WHILE
            X(1)=X(1)+X_DISCRETISATION
          ENDDO !WHILE
          IF(NDT2.EQ.0)THEN !make sure each element contains a point
            NDT2=NDT2+1
            NDT=NDT+1
            DO nj=1,3
              XI(nj)=0.5d0
            ENDDO
            DO nj=1,3
              ZD(nj,NDT)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     &          XI,XE(1,nj))
            ENDDO
            LD(NDT)=ne
            DO nj=1,NJT
              ZD(nj,NDT)=X(nj)
              WD(nj,NDT)=1.0d0
              XID(nj,NDT)=XI(nj) !store xi location
            ENDDO !nj
            NDDATA(NDT,nr)=NDT
            NDDATA(0,nr)=NDDATA(0,nr)+1
          ENDIF
          WRITE(OP_STRING,'('' ne = '',I6,'' #points(ne) = '',
     '      I6,'' total points = '',I6)') ne,NDT2,NDT
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !noelem
        
      ELSE IF(SPREAD_TYPE(1:2).EQ.'XI')THEN !Xi-spacing
        NDT=0
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '      nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
          NDT2=0
C Generate the xi-spaced seed points
          xi_discretisation=x_discretisation
          number_of_points_xi=INT(1.d0/xi_discretisation)-1
          
          XI(1)=xi_discretisation
          DO i=1,number_of_points_xi
            XI(2)=xi_discretisation
            DO j=1,number_of_points_xi
              XI(3)=xi_discretisation
              DO k=1,number_of_points_xi
                DO nj=1,NJT !3
                  X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '              nb,1,XI,XE(1,nj))
                ENDDO !nj
                NDT=NDT+1 !counts RPs for whole lobe
                NDT2=NDT2+1
                CALL ASSERT(NDM.GT.NDT,
     &            '>>Increase NDM in .ippara file',ERROR,*9999)
                LD(NDT)=ne
                DO nj=1,NJT
                  ZD(nj,NDT)=X(nj) !stores random point coords
                  WD(nj,NDT)=1.0d0
                  XID(nj,NDT)=XI(nj) !store xi location
                ENDDO !nj
                NDDATA(NDT,nr)=NDT
                NDDATA(0,nr)=NDDATA(0,nr)+1
                XI(3)=XI(3)+xi_discretisation
              ENDDO !k
              XI(2)=XI(2)+xi_discretisation
            ENDDO !j
            XI(1)=XI(1)+xi_discretisation
          ENDDO !i
          IF(DOP)THEN
            WRITE(OP_STRING,'('' ne = '',I6,'' density = '',
     '        D12.4)') ne,NDT2/VOL
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO !noelem (ne)
      ELSE
        WRITE(OP_STRING,'('' This type not currently implemented '')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      WRITE(OP_STRING,'('' Number of seed points = '',I6)') NDT
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)

      WRITE(OP_STRING,
     &  '('' Note: use FEM DEFINE DATA;W;FILE to write data points'')')
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,
     &  '(''       and FEM DEFINE;XI;W;FILE to write host elements'')')
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)

      CALC_XI=.TRUE.
      CALL_XI=.TRUE.
      CALL_DATA=.TRUE.
      MAKE_DATA=.TRUE.
      
      CALL EXITS('GNNEZD')
      RETURN
 9999 CALL ERRORS('GNNEZD',ERROR)
      CALL EXITS('GNNEZD')
      RETURN 1
      END



