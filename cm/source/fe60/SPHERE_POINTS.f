      SUBROUTINE SPHERE_POINTS(BOUNDARY,NKJ,NPNODE,nr,NVJP,
     '  VERTICES_XYZ,XP,REGULAR,ERROR,*)

C#### Subroutine: SPHERE_POINTS
C###  Description:
C###    SPHERE_POINTS generates regular points on a unit sphere, with a
C###    'hole' defined by an angle limit.

C*** Created by KSB January 2002. Last modified 30th May 2002.
C*** Currently only working for regular points, not random points.
C*** Regular point generation based on icosahedron structure, this
C*** then refined.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'pulm00.cmn'

!     Parameter List
      INTEGER BOUNDARY(0:N_BOUNDARY*N_ALVEOLI*2,2),NKJ(NJM,NPM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NVJP(NJM,NPM)
      REAL*8 VERTICES_XYZ(3*N_VERT_EST),XP(NKM,NVM,NJM,NPM)
      LOGICAL REGULAR
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,IERROR,j,k,l,N,n_b,ncheck,nj,
     '  no,nonode,notri,np,np2,np3,nreject,NTRI,ntri2,
     '  point,point2,TOTAL_POINTS,TRI(3,3000)
      REAL*8 alpha1,alpha2,CM_RANDOM_NUMBER,
     '  LENGTH,LENGTH2,LLIMIT,radius,REG_ALPHA(12),REG_THETA(12),
     '  temp,theta,tri_temp(3),tri_temp2(3),X(3),z_limit
      LOGICAL CONTINU,NP_REPEAT,REJECT

      CALL ENTERS('SPHERE_POINTS',*9999)

      low_limit=low_limit/180.d0*PI !changing to radians
      radius=1.d0
      I=1
      np=NPT(0) !node #
      nonode=NPNODE(0,nr)
      IF(REGULAR)THEN
        !set spherical coordinates for first 12 points
        REG_ALPHA(1)=0.d0
        REG_ALPHA(2)=DACOS((3.d0-DSQRT(3.d0))/3.d0)
        REG_ALPHA(7)=(REG_ALPHA(2)+PI)/2.d0
        REG_ALPHA(12)=PI
        REG_THETA(1)=0.d0
        REG_THETA(2)=0.d0
        REG_THETA(7)=PI/5.d0
        REG_THETA(12)=0.d0
        DO k=1,4
          REG_ALPHA(2+k)=REG_ALPHA(2)
          REG_ALPHA(7+k)=REG_ALPHA(7)
          REG_THETA(2+k)=REG_THETA(1+k)+2.d0*PI/5.d0
          REG_THETA(7+k)=REG_THETA(6+k)+2.d0*PI/5.d0
        ENDDO !k
        I=I-1
        DO k=1,12 !for each of the 12 initial points
          np=np+1
          nonode=nonode+1
          NPNODE(nonode,nr)=np !creates new node in region nr
          VERTICES_XYZ(I+1)=radius*DCOS(REG_THETA(k))*
     '      DSIN(REG_ALPHA(k)) !x
          VERTICES_XYZ(I+2)=radius*DSIN(REG_THETA(k))*
     '      DSIN(REG_ALPHA(k)) !y
          VERTICES_XYZ(I+3)=radius*DCOS(REG_ALPHA(k)) !z
          DO nj=1,3
            XP(1,1,nj,np)=VERTICES_XYZ(I+nj) !nodal co-ordinates
          ENDDO !j
          I=I+3
          n_b=0
        ENDDO !k
        I=I+1
        DO j=1,5  !sets up 1st 5 triangle connectivity
          TRI(1,j)=1
          TRI(2,j)=j+1
          TRI(3,j)=j+2
        ENDDO !j
        N=6
        DO j=6,14,2
          N=N+1
          TRI(1,j)=j/2-1
          TRI(2,j)=N
          TRI(3,j)=j/2+(5*MOD(j,2)) !if even = j/2, if odd = j/2+5
          TRI(1,j+1)=TRI(1,j)+1
          TRI(2,j+1)=N
          TRI(3,j+1)=j/2+(5*MOD(j+1,2))
        ENDDO !j
        DO k=7,11
          TRI(1,j)=k
          TRI(2,j)=k+1
          TRI(3,j)=12
          j=j+1
        ENDDO !k
        TRI(3,5)=2 !correcting connectivity array
        TRI(2,20)=7
        TRI(3,14)=2
        TRI(3,15)=7
        TRI(1,15)=2
        NTRI=20 !number of triangles
        ntri2=NTRI
        DO no=1,N_REFINE
          DO notri=1,NTRI
            DO k=1,3
              l=k+1
              IF(l.GT.3) l=1
              nonode=nonode+1
              np=np+1
              NPNODE(nonode,nr)=np !new node
              np2=NPNODE(NPNODE(0,nr)+TRI(k,notri),nr)
              np3=NPNODE(NPNODE(0,nr)+TRI(l,notri),nr)
              DO nj=1,3 !new pt btwn each 2 pts of tri
                TRI_TEMP(nj)=0.5d0*(XP(1,1,nj,np2)+XP(1,1,nj,np3))
              ENDDO
              CALL SPHERE_PROJECT(TRI_TEMP,TRI_TEMP2,IERROR,ERROR,*9999)
              DO nj=1,3
                VERTICES_XYZ(I)=TRI_TEMP2(nj)
                XP(1,1,nj,np)=VERTICES_XYZ(I)
                I=I+1
              ENDDO
            ENDDO !k
C... make the 4 new triangles
            ntri2=ntri2+1 !subtract NPNODE(0,nr) in case nodes already
            TRI(1,ntri2)=nonode-NPNODE(0,nr) !in region nr
            TRI(2,ntri2)=nonode-NPNODE(0,nr)-1
            TRI(3,ntri2)=TRI(3,notri)
            ntri2=ntri2+1
            TRI(1,ntri2)=nonode-NPNODE(0,nr)-2
            TRI(2,ntri2)=TRI(2,notri)
            TRI(3,ntri2)=nonode-NPNODE(0,nr)-1
            ntri2=ntri2+1
            TRI(1,ntri2)=nonode-NPNODE(0,nr)-2
            TRI(2,ntri2)=nonode-NPNODE(0,nr)-1
            TRI(3,ntri2)=nonode-NPNODE(0,nr)
            TRI(2,notri)=nonode-NPNODE(0,nr)-2
            TRI(3,notri)=nonode-NPNODE(0,nr)
          ENDDO !notri
          NTRI=ntri2
        ENDDO !no
        TOTAL_POINTS=nonode
C... the following is to remove points from the hole in the sphere
C... & to remove any duplicated points (i.e nodes with same XP co-ords)
        np2=NPT(0)
        z_limit=-DCOS(low_limit+z_bound*PI)*radius
        IF(low_limit.EQ.0.d0) z_limit=-1.0d0-z_bound
        DO WHILE(np2.LT.NPT(0)+TOTAL_POINTS)
          np2=np2+1
          np=NPT(0)
          NP_REPEAT=.FALSE.
          DO WHILE(.NOT.NP_REPEAT.AND.np.LT.NPT(0)+TOTAL_POINTS)
            np=np+1
            IF(np.NE.np2) THEN
              IF(XP(1,1,1,np).EQ.XP(1,1,1,np2).AND.XP(1,1,2,np).EQ.
     '          XP(1,1,2,np2).AND.XP(1,1,3,np).EQ.XP(1,1,3,np2))
     '          NP_REPEAT=.TRUE.
            ENDIF
          ENDDO !WHILE
          IF(XP(1,1,3,np2).LE.z_limit.OR.NP_REPEAT) THEN!remove np2
            DO point2=np2,NPT(0)+TOTAL_POINTS-1
              DO nj=1,3
                XP(1,1,nj,point2)=XP(1,1,nj,point2+1) !moves up 1 place
                VERTICES_XYZ((point2-NPT(0)-1)*3+nj)=
     '            VERTICES_XYZ((point2-NPT(0))*3+nj)
              ENDDO !nj
            ENDDO !point2
            np2=np2-1 !point removed
            np=np-1
            nonode=nonode-1
            TOTAL_POINTS=TOTAL_POINTS-1
            I=I-3
          ENDIF !XP
        ENDDO !WHILE
        point=TOTAL_POINTS
      ELSE
C...  generating random points on unit sphere
C...  Random x,y,z, project to surface
        nreject=0
        z_limit=-DCOS(low_limit+0.1d0*PI)*radius!0.1 is arbitrary value
        I=1             ! to ensure that random points sit away from bndry
        LLIMIT=4.d0*DSQRT(PI/(N_VERT_EST-2*N_BOUNDARY)
     '    *DSQRT(3.d0)/6.d0)*(2.d0-z_limit)/2.d0*0.45d0
        point=0
        np=NPT(0)
        nonode=NPNODE(0,nr)
        CONTINU=.TRUE.
        DO WHILE(CONTINU)
          LENGTH=0.d0
          DO J=1,3 !3 coordinates
            temp=CM_RANDOM_NUMBER(0) !random coordinate
            X(J)=2.d0*(temp-0.5d0) !makes the points btwn -1 and 1
            LENGTH=LENGTH+X(J)**2.d0
          ENDDO !J
          LENGTH=DSQRT(LENGTH)
          IF(X(3)/LENGTH.GE.z_limit)THEN !MIGHT accept as a point
            REJECT=.FALSE.
            DO ncheck=NPT(0)+1,np
              LENGTH2=0.d0
              DO nj=1,3
                LENGTH2=LENGTH2+(X(nj)/LENGTH-XP(1,1,nj,ncheck))**2.d0
              ENDDO
              LENGTH2=DSQRT(LENGTH2)
              IF(LENGTH2.LT.LLIMIT)THEN
                REJECT=.TRUE. !reject
              ENDIF
            ENDDO !ncheck
            IF(.NOT.REJECT)THEN
              point=point+1
              nonode=nonode+1
              np=np+1
              NPNODE(nonode,nr)=np !new node
              DO nj=1,3
                VERTICES_XYZ(I)=X(nj)/LENGTH
                XP(1,1,nj,np)=VERTICES_XYZ(I)
                I=I+1
              ENDDO !nj
            ELSE
              nreject=nreject+1
            ENDIF
          ENDIF
          IF(point.GE.N_VERT_EST-2*N_BOUNDARY)THEN
            CONTINU=.FALSE.
          ELSE IF(nreject.GE.2000)THEN
            CONTINU=.FALSE.
            WRITE(OP_STRING,'(''Warning: too many rejected '//
     '        ' points,# made= '',I5)') point
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO !np
      ENDIF
C     Create B & IB nodes for Delaunay trianglulation
      alpha1=(PI-low_limit)-z_bound*PI !angle for IB nodes
      alpha2=(PI-low_limit)+z_bound*PI !angle for B nodes
      IF(N_BOUNDARY.EQ.0) THEN
        theta=2.d0*PI
      ELSE
        theta=PI/N_BOUNDARY !initial IB/B node location
      ENDIF
      radius=1.d0
      n_b=0
      DO WHILE(theta.LT.2.d0*PI)
        nonode=nonode+1
        np=np+1
        NPNODE(nonode,nr)=np
        n_b=n_b+1 !counts # of IB and B points
        BOUNDARY(n_b,1)=np !IB node
        BOUNDARY(n_b,2)=np+1 !B node
        VERTICES_XYZ(I)=radius*dcos(theta)*dsin(alpha1) !x
        XP(1,1,1,np)=VERTICES_XYZ(I)
        I=I+1
        VERTICES_XYZ(I)=radius*dsin(theta)*dsin(alpha1) !y
        XP(1,1,2,np)=VERTICES_XYZ(I)
        I=I+1
        VERTICES_XYZ(I)=radius*dcos(alpha1) !z
        XP(1,1,3,np)=VERTICES_XYZ(I)
        I=I+1
        nonode=nonode+1
        np=np+1
        NPNODE(nonode,nr)=np !new node
        VERTICES_XYZ(I)=radius*dcos(theta)*dsin(alpha2) !x
        XP(1,1,1,np)=VERTICES_XYZ(I)
        I=I+1
        VERTICES_XYZ(I)=radius*dsin(theta)*dsin(alpha2) !y
        XP(1,1,2,np)=VERTICES_XYZ(I)
        I=I+1
        VERTICES_XYZ(I)=radius*dcos(alpha2) !z
        XP(1,1,3,np)=VERTICES_XYZ(I)
        I=I+1
        theta=theta+2.d0*PI/N_BOUNDARY
      ENDDO !while

      DO np2=NPT(0)+1,np
        DO nj=1,3
          NKJ(nj,np2)=1
          NVJP(nj,np2)=1
        ENDDO !nj
      ENDDO !nonde2

      TOTAL_POINTS=np-NPT(0)
      NPNODE(0,nr)=nonode
      NPNODE(0,0)=np
      NPT(0)=np
      IF(N_VERT_EST.NE.TOTAL_POINTS) THEN
        N_VERTICES=TOTAL_POINTS !IPMESH2 calculates incorrect #
C.. for regular mesh N_VERTICES=no. points + no. boundary points
C        write(*,*) "** Number of vertices changed in SPHERE_POINTS"
      ENDIF
      NDT=N_VERTICES

      CALL EXITS('SPHERE_POINTS')
      RETURN
 9999 CALL ERRORS('SPHERE_POINTS',ERROR)
      CALL EXITS('SPHERE_POINTS')
      RETURN 1
      END


