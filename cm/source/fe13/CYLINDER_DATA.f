      SUBROUTINE CYLINDER_DATA(NBJ,NDDATA,NDP,NENP,NPLIST,NPNE,nr,NVJE,
     &  density,WD,XP,ZD,FIRST,ERROR,*)
      
C#### Subroutine: CYLINDER_DATA
C###  Description:
C###  CYLINDER_DATA2 calculates data positions over the surface of a
C###  cylinder-based 1D tree.


      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),NDP(NDM),
     &  NENP(NPM,0:NEPM,0:NRM),NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),nr,
     &  NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 WD(NJM,NDM),XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
      LOGICAL FIRST
!     Local Variables
      INTEGER i,N,nb,nb_radius,ne,nj,noelem,nonode,nonode_1,np,nv
      REAL*8 DOT_PROD
      INTEGER j,k,M,NUM_POINTS,N_CIRCLE,np0,NE_LOCAL(36,4),
     '  num_in_1,num_in_2,num_points_generated,num_points_to_generate,
     &  NPOINT(36)
      REAL*8 angle,AREA_ELEMENT,centroid(3),COS_T,
     &  density,distance,length,mult1,mult2,radius(4),N1(3),
     &  new_p(3,1000),NML(3),NML_0(4),PATH(3,4),phi(3),R(3,3),SIN_T,
     &  THETT,temp(3),temp1,THETT_INC,u(3),v(3),V_I(3,4),W(3),XI1,
     &  XI1_inc,XI2,XI2_inc,XP_BP(3),XP_EP(3,3)
C CM_RANDOM_NUMBER,
      LOGICAL REGULAR

C XP_BP stores the coordinates of the bifurcation point
C XP_EP stores the coordinates of the points half-way along the elements
C PATH stores the direction of each element from the bifurcation point
      
      CALL ENTERS('CYLINDER_DATA',*9999)
      REGULAR=.TRUE.

C     Treat each node in NPLIST as a bifurcation point. Find the 3
C     1D elements that form the bifurcation around the point. Make a
C     ring of points half-way along each of the 3 elements, and 3 semi
C     -rings of points that join over the bifurcation. Join up the
C     points to make psuedo elements. Generate random data points over
C     each 'element'.

      
      N_CIRCLE=12 !number of points == nodes around the cylinder end
      THETT_INC=2.d0*PI/N_CIRCLE !angle increment, for N_CIRCLE points
      NDT=0 !initial number of data points to zero

      IF(FIRST)THEN
        nonode_1=1
      ELSE
        nonode_1=2
      ENDIF
      
      DO nonode=nonode_1,NPLIST(0) !for each specified bifurcation point
        np=NPLIST(nonode)

        DO nj=1,NJT  !bifurcation point is at node np, store in XP_BP
          XP_BP(nj)=XP(1,1,nj,np)
        ENDDO !nj

        NUM_POINTS=0 !initialise the number of points that are created
        i=0
        IF(NENP(np,0,nr).NE.1)THEN !more than one adjacent element
          IF(NENP(np,0,nr).GT.3)THEN
            WRITE(OP_STRING,
     &        '('' Warning: node'',I6,'' adjoins > 3 elements'')') np
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE
            DO noelem=1,NENP(np,0,nr) !for each adjoining element
              ne=NENP(np,noelem,nr)
              nb=NBJ(1,ne)
              i=i+1
              nb_radius=NBJ(nj_radius,ne)
              IF(i.EQ.1)THEN
                np0=NPNE(1,nb,ne) !first node of element
                nv=NVJE(2,nb_radius,nj_radius,ne)
              ELSE
                np0=NPNE(2,nb,ne) !second node of element
                nv=NVJE(1,nb_radius,nj_radius,ne)
              ENDIF
              radius(i)=XP(1,nv,nj_radius,np)
              DO nj=1,NJT
                XP_EP(nj,i)=0.5d0*(XP(1,1,nj,np0)+XP_BP(nj))
                PATH(nj,i)=XP(1,1,nj,np0)-XP_BP(nj) !direction
              ENDDO !nj
              CALL NORMALISE(3,PATH(1,i),ERROR,*9999)

            ENDDO !noelem
            DO nj=1,NJT
              PATH(nj,4)=PATH(nj,1)
            ENDDO !nj
            angle=DACOS(DOT_PROD(PATH(1,1),PATH(1,2)))
            IF(angle.GT.(PI-0.01d0)) write(*,*)'angle too big 12',np
            angle=DACOS(DOT_PROD(PATH(1,2),PATH(1,3)))
            IF(angle.GT.(PI-0.01d0)) write(*,*)'angle too big 23',np
            angle=DACOS(DOT_PROD(PATH(1,3),PATH(1,3)))
            IF(angle.GT.(PI-0.01d0)) write(*,*)'angle too big 31',np
            radius(4)=radius(1)
            
C         Set up a circle of new points at each bifurcation ending. The
C         first point on the circle is the projection of the normal to
C         the plane containing all 3 mid-element points onto the plane
C         containing
C         the circle of points (normal vector is PATH).
            
C         Calculate the normal to the plane containing all XP_EP nodes
            CALL PLANE_FROM_3_PTS(NML_0,2,XP_EP(1,1),XP_EP(1,2),
     '        XP_EP(1,3),ERROR,*9999) !normal to plane is in NML_0
            DO i=1,3 !parent, daughter1, daughter2
C           Calculate the position of the first point, which gives a
C           rotation vector V to define the location of the other
C           points.
              DO nj=1,3
                NML(nj)=PATH(nj,i)
              ENDDO
              THETT=dacos(DOT_PROD(NML,NML_0)) !angle between two normals
              length=0.d0
              IF(THETT.LT.PI/2.d0)THEN
                DO nj=1,3
                  V(nj)=-PATH(nj,i)*radius(i)/dtan(THETT)+NML_0(nj)
     '              *radius(i)/dsin(THETT) !direction of first new point
                  length=length+V(nj)**2
                ENDDO
              ELSE IF(THETT.GT.PI/2.d0)THEN
                THETT=PI-THETT
                DO nj=1,3
                  V(nj)=PATH(nj,i)*radius(i)/dtan(THETT)+NML_0(nj)
     '              *radius(i)/dsin(THETT) !direction of first new point
                  length=length+V(nj)**2
                ENDDO
              ELSE
                DO nj=1,3
                  V(nj)=NML_0(nj)*radius(i)
                  length=length+V(nj)**2
                ENDDO
              ENDIF

              length=dsqrt(length)
              IF(DABS(length).LT.LOOSE_TOL)THEN
                write(*,*) 'zero length for V rotation vector', i
              ENDIF
              
              CALL NORMALISE(3,V,ERROR,*9999) !V is a unit vector
C           Rotate V about the element to calculate the point positions
              THETT=0.d0 !first location
              DO N=1,N_CIRCLE
C             The rotation matrix is used as follows: (from http:
C             //research.microsoft.com/~hollasch/cgindex/math
C             /rotvec.html)
C             let [v] = [vx,vy, vz] be the vector to be rotated
C             [l] = [lx, ly, lz] is unit vector about which to rotate
C             [L] = [0 lz -ly; -lz 0 lz; ly -lx 0]
C             then the rotation matrix [R] is given by:
C             [R] = I + sin(T)*[L] + (1-cos(T))*([L]*[L])
C             rotation of vector is given by [v]' = [R]*[v]
                
                SIN_T=DSIN(THETT)
                COS_T=1.d0-DCOS(THETT)
                
                R(1,1)=1.d0-COS_T*(NML(3)**2+NML(2)**2)
                R(1,2)=NML(3)*SIN_T+COS_T*NML(1)*NML(2)
                R(1,3)=-NML(2)*SIN_T+COS_T*NML(1)*NML(3)
                
                R(2,1)=-NML(3)*SIN_T+COS_T*NML(1)*NML(2)
                R(2,2)=1.d0-COS_T*(NML(3)*NML(3)+NML(1)*NML(1))
                R(2,3)=NML(1)*SIN_T+COS_T*NML(2)*NML(3)
                
                R(3,1)=NML(2)*SIN_T+COS_T*NML(1)*NML(3)
                R(3,2)=-NML(1)*SIN_T+COS_T*NML(2)*NML(3)
                R(3,3)=1.d0-COS_T*(NML(1)*NML(1)+NML(2)*NML(2))
                NUM_POINTS=NUM_POINTS+1
                
C             Apply the rotation to V
                length=0.d0
                DO nj=1,NJT
                  TEMP(nj)=R(nj,1)*V(1)+R(nj,2)*V(2)+R(nj,3)*V(3)
                  length=length+temp(nj)**2
                ENDDO
                length=dsqrt(length)
                IF(DABS(length).LT.LOOSE_TOL)THEN
                  write(*,*) 'zero length for applied rotation',i
                ENDIF
                
                CALL NORMALISE(3,TEMP,ERROR,*9999)
                DO nj=1,3
                  new_p(nj,NUM_POINTS)=TEMP(nj)*radius(i)+XP_EP(nj,i)
                ENDDO !nj
                
                IF(i.EQ.1)THEN
                  THETT=THETT-THETT_INC !increment the angle for the next rotation
                ELSE
                  THETT=THETT+THETT_INC !increment the angle for the next rotation
                ENDIF
              ENDDO !N
            ENDDO !i
            
            DO i=1,3 !each pair of points
              DO nj=1,3
                U(nj)=PATH(nj,i)
                V(nj)=PATH(nj,i+1)
              ENDDO
              CALL CROSS(U,V,N1)

              length=0.d0
              DO nj=1,3
                length=length+N1(nj)**2
              ENDDO
              length=dsqrt(length)
              IF(DABS(length).LT.LOOSE_TOL)THEN
                write(*,*)'zero length for cross-product',i,np
                angle=DACOS(DOT_PROD(U,V))
                write(*,*) 'angle between',angle
                
              ENDIF

              CALL NORMALISE(3,N1,ERROR,*9999)
              phi(1)=DACOS(DOT_PROD(U,V))
              phi(2)=DABS(datan((RADIUS(I)+RADIUS(I+1)*DCOS(PHI(1)))
     '          /(RADIUS(I+1)*DSIN(PHI(1)))))
              phi(3)=DABS(phi(2)-PI/2.d0)
              W(3)=U(1)
              W(3)=W(3)*(DCOS(phi(1)-phi(3))*(N1(1)*U(2)-N1(2)*U(1))
     '          +DCOS(phi(3))*(V(1)*N1(2)-V(2)*N1(1)))
              temp1=((U(1)*V(2)-V(1)*U(2))*(U(1)*N1(3)-N1(1)*U(3))-(U(1)
     '          *N1(2)-N1(1)*U(2))*(U(1)*V(3)-V(1)*U(3)))
              IF(temp1.EQ.0.d0)THEN
                W(3)=0.d0
              ELSE
                W(3)=W(3)/
     '            ((U(1)*V(2)-V(1)*U(2))*(U(1)*N1(3)-N1(1)*U(3))-(U(1)
     '            *N1(2)-N1(1)*U(2))*(U(1)*V(3)-V(1)*U(3)))
              ENDIF
              IF((U(1)*N1(2)-N1(1)*U(2)).NE.0.d0)THEN
                W(2)=(-N1(1)*DCOS(phi(3))-(U(1)*N1(3)-N1(1)*U(3))*W(3))
     '            /(U(1)*N1(2)-N1(1)*U(2))
              ELSE
                W(2)=0.d0
              ENDIF
              IF(U(1).NE.0.d0)THEN
                W(1)=(DCOS(phi(3))-U(2)*W(2)-U(3)*W(3))/U(1)
              ELSE
                U(1)=0.d0
              ENDIF
              
              CALL NORMALISE(3,W,ERROR,*9999)
c              IF(DABS(phi(1)-PI).LT.0.01d0) phi(1)=0.d0
c              IF(DABS(phi(2)-PI).LT.0.01d0) phi(2)=0.d0
              distance=radius(i)/DCOS(PI-phi(1)-phi(2))
              distance=0.5d0*(radius(i)+radius(i+1))
              DO nj=1,3
                V_I(i,nj)=W(nj)*distance
              ENDDO
c              IF(np.EQ.138881)THEN
c                write(*,*) i,distance,radius(i),phi(1),phi(2)
c              ENDIF
            ENDDO
            DO nj=1,3
              NEW_P(nj,NUM_POINTS+1)=XP_BP(nj)+NML_0(nj)*(radius(1)
     &          +radius(2)+radius(3))/3.d0
              NEW_P(nj,NUM_POINTS+7)=XP_BP(nj)-NML_0(nj)*(radius(1)
     &          +radius(2)+radius(3))/3.d0
              NEW_P(nj,NUM_POINTS+4)=XP_BP(nj)+V_I(1,nj)
c              NEW_P(nj,NUM_POINTS+15)=XP_BP(nj)+V_I(2,nj)
              NEW_P(nj,NUM_POINTS+15)=XP_BP(nj)-PATH(nj,1)*(radius(2)
     &          +radius(3))/2.d0
              NEW_P(nj,NUM_POINTS+10)=XP_BP(nj)+V_I(3,nj)
            ENDDO !nj
c            IF(np.EQ.52641)THEN
c              write(*,*) v_i(1,1),v_i(1,2),v_i(1,3)
c              write(*,*) num_points+1,num_points+7,num_points+4,
c     &          num_points+15,num_points+10
c            ENDIF
            
            NUM_POINTS=NUM_POINTS+5
            
            NPOINT(1)=3*N_CIRCLE+1
            NPOINT(2)=3*N_CIRCLE+4
            NPOINT(3)=3*N_CIRCLE+7
            NPOINT(4)=3*N_CIRCLE+10
            NPOINT(5)=3*N_CIRCLE+1
            NPOINT(6)=3*N_CIRCLE+15
            NPOINT(7)=3*N_CIRCLE+7

            DO M=1,6
              mult1=1.d0
              mult2=0.d0
              DO i=1,2
                mult1=mult1-1.d0/3.d0
                mult2=mult2+1.d0/3.d0
                DO nj=1,3
                  temp(nj)=NEW_P(nj,NPOINT(M))*mult1+NEW_P(nj,NPOINT(M
     &              +1))*mult2-XP_BP(nj)
                ENDDO !nj
                CALL NORMALISE(3,temp,ERROR,*9999)
                IF(M.NE.5)THEN
                  DO nj=1,3
                    NEW_P(nj,NPOINT(M)+i)=XP_BP(nj)+temp(nj)*radius(1)
                  ENDDO !nj
                ELSE
                  DO nj=1,3
                    NEW_P(nj,NPOINT(M+1)-(3-i))=XP_BP(nj)+temp(nj)
     '                *radius(1)
                  ENDDO !nj
                ENDIF                
              ENDDO !i
            ENDDO !M
            
            NUM_POINTS=NUM_POINTS+12

!around the parent branch
            DO M=1,N_CIRCLE
              NE_LOCAL(M,1)=3*N_CIRCLE+M
              NE_LOCAL(M,3)=M
              NE_LOCAL(N_CIRCLE+M,1)=N_CIRCLE+M
              NE_LOCAL(N_CIRCLE+M,3)=3*N_CIRCLE+M
              NE_LOCAL(2*N_CIRCLE+M,1)=2*N_CIRCLE+M
              NE_LOCAL(2*N_CIRCLE+M,3)=3*N_CIRCLE+M
            ENDDO
            DO M=1,5
              NE_LOCAL(N_CIRCLE+7+M,3)=54-M
              NE_LOCAL(2*N_CIRCLE+1+M,3)=48+M
            ENDDO !M
            
            DO M=1,3*N_CIRCLE-1
              NE_LOCAL(M,2)=NE_LOCAL(M+1,1)
              NE_LOCAL(M,4)=NE_LOCAL(M+1,3)
            ENDDO
            DO M=1,N_CIRCLE-1
              NE_LOCAL(2*N_CIRCLE+M,4)=NE_LOCAL(2*N_CIRCLE+M+1,3)
            ENDDO !M
            NE_LOCAL(3*N_CIRCLE,4)=NE_LOCAL(2*N_CIRCLE+1,3)
            DO M=1,3
              NE_LOCAL(M*N_CIRCLE,2)=NE_LOCAL(1+(M-1)*N_CIRCLE,1)
              NE_LOCAL(M*N_CIRCLE,4)=NE_LOCAL(1+(M-1)*N_CIRCLE,3)
            ENDDO !M

c            IF(np.EQ.52641)THEN
c              write(*,*) np,new_p(1,40),new_p(2,40),new_p(3,40)
c            ENDIF
            
            num_points_generated=0
            DO i=1,36
              DO nj=1,3
                centroid(nj)=0.d0
              ENDDO !nj
              DO M=1,4
                NPOINT(M)=NE_LOCAL(i,M)
                DO nj=1,3
                  centroid(nj)=centroid(nj)+0.25d0*NEW_P(nj,NPOINT(M))
                ENDDO !nj
              ENDDO !M
              NPOINT(5)=NE_LOCAL(i,1)
C           Assume that the area of the element can be approximated by
C           summing the triangles made by the centroid and each pair of
C           bounding nodes.
              AREA_ELEMENT=0.d0
              DO M=1,4 !take each pair of adjacent nodes in turn
                DO nj=1,3 !make the vectors from the centroid
                  u(nj)=NEW_P(nj,NPOINT(M))-centroid(nj)
                  v(nj)=NEW_P(nj,NPOINT(M+1))-centroid(nj)
                ENDDO !nj
C             Area of triangle = half * length of cross-product              
                AREA_ELEMENT=AREA_ELEMENT+0.25d0*dsqrt((u(2)*v(3)-u(3)
     '            *v(2))**2+(u(3)*v(1)-u(1)*v(3))**2+(u(1)*v(2)-u(2)
     &            *v(1))**2)/2.d0
              ENDDO !M
              
              IF(.NOT.REGULAR)THEN
C           Make data points over the 'surfaces', by generating random
C           xi_1 and xi_2 coordinates (between 0 to 1 inclusive) then
C           scaling linearly back to global coordinates. The estimated
C           area of the element is used to determine the number of
C           points that should be generated.
                num_points_to_generate=NINT(density*AREA_ELEMENT)
                
                CALL ASSERT(num_points_to_generate.GT.0,'Increase
     &            density',ERROR,*9999)
                DO M=1,num_points_to_generate
                  num_points_generated=num_points_generated+1
c              XI1=CM_RANDOM_NUMBER(0)
c              XI2=CM_RANDOM_NUMBER(0)
                  XI1=0.5d0
                  XI2=0.5d0
                  NDT=NDT+1
                  DO nj=1,3
                    ZD(nj,NDT)=
     '                (1.d0-XI1)*(1.d0-XI2)*NEW_P(nj,NPOINT(1))+
     '                XI1*(1.d0-XI2)*NEW_P(nj,NPOINT(2))+
     '                (1.d0-XI1)*XI2*NEW_P(nj,NPOINT(3))+
     '                XI1*XI2*NEW_P(nj,NPOINT(4))
                    WD(nj,NDT)=1.0d0
                  ENDDO !nj
                  NDP(NDT)=NDT
                  NDDATA(NDT,nr)=NDT
                  NDDATA(0,nr)=NDDATA(0,nr)+1
                ENDDO !M
              ELSE !.NOT.RANDOM, do regular spacing in Xi
                num_in_1=4
                num_in_2=6
                XI1_inc=DBLE(1.d0/num_in_1)
                XI2_inc=DBLE(1.d0/num_in_2)
                XI1=0.d0
                DO j=1,num_in_1 !from 0 to < 1
                  XI2=0.d0
                  DO k=1,num_in_2 !from 0 to < 1
                    NDT=NDT+1
                    CALL ASSERT(NDT.LE.NDM,'Increase NDM',ERROR,*9999)
                    DO nj=1,3
                      ZD(nj,NDT)=
     '                  (1.d0-XI1)*(1.d0-XI2)*NEW_P(nj,NPOINT(1))+
     '                  XI1*(1.d0-XI2)*NEW_P(nj,NPOINT(2))+
     '                  (1.d0-XI1)*XI2*NEW_P(nj,NPOINT(3))+
     '                  XI1*XI2*NEW_P(nj,NPOINT(4))
                      WD(nj,NDT)=1.0d0
                    ENDDO !nj
                    NDP(NDT)=NDT
                    NDDATA(NDT,nr)=NDT
                    NDDATA(0,nr)=NDDATA(0,nr)+1
                    XI2=XI2+XI2_inc
                  ENDDO !k
                  XI1=XI1+XI1_inc
                ENDDO !j
              ENDIF
            ENDDO !i
          ENDIF
        ELSE
 !terminal branch
        ENDIF
        IF(nonode.EQ.2)THEN
        CALL OPENF(IOFILE4,'DISK','test_out.ipnode',
     '    'NEW','SEQUEN','FORMATTED',160,ERROR,*9999)
        IOFI=IOFILE4
        
        WRITE(OP_STRING,'('' CMISS Version 1.21 ipnode File Version ''
     '    ''2'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Heading:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' The number of nodes is [    1]: '',I5)'
     '    )NUM_POINTS
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Number of coordinates [3]: 3'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' Do you want prompting for different versions''
     '    '' of nj=1 [N]? N'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' Do you want prompting for different versions''
     '    '' of nj=2 [N]? N'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' Do you want prompting for different versions''
     '    '' of nj=3 [N]? N'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' The number of derivatives for coordinate 1 ''
     '    ''is [0]:0'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' The number of derivatives for coordinate 2 ''
     '    ''is [0]:0'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' The number of derivatives for coordinate 3 ''
     '    ''is [0]:0'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        
        DO i=1,NUM_POINTS
          WRITE(OP_STRING,'('' Node number [    1]: '',I5)') i
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '('' The Xj(1) coordinate is [ 0.00000E+00]: '',D13.5)')
     '      NEW_P(1,i)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '('' The Xj(2) coordinate is [ 0.00000E+00]: '',D13.5)')
     '      NEW_P(2,i)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '('' The Xj(3) coordinate is [ 0.00000E+00]: '',D13.5)')
     '      NEW_P(3,i)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !i
        CLOSE(IOFILE4)
        
        CALL OPENF(IOFILE5,'DISK','test_out.ipelem',
     '    'NEW','SEQUEN','FORMATTED',160,ERROR,*9999)
        IOFI=IOFILE5
        WRITE(OP_STRING,'('' CMISS Version 1.21 ipelem File Version ''
     '    ''2'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Heading:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' The number of elements is [1]: '',I5)')
     '    36
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO i=1,36
          WRITE(OP_STRING,'('' Element number [    1]: '',I5)') i
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '('' The number of geometric Xj-coordinates is [3]: 3'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '('' The basis function type for geometric variable 1 ''
     '      ''is [1]: 1'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '('' The basis function type for geometric variable 2 ''
     '      ''is [1]: 1'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '('' The basis function type for geometric variable 3 ''
     '      ''is [1]: 1'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '('' Enter the 4 global numbers for basis 1: '',4(I5))')
     '      NE_LOCAL(i,1),NE_LOCAL(i,2),NE_LOCAL(i,3),
     '      NE_LOCAL(i,4)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !i
      ENDIF

      ENDDO !nonode
          

c          CALL OPENF(IOFILE6,'DISK','test.exdata',
c     '      'NEW','SEQUEN','FORMATTED',160,ERROR,*9999)
c          IOFI=IOFILE6
c          WRITE(OP_STRING,'('' Group name: surface_data'')')
c          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c          WRITE(OP_STRING,'('' #Fields=1'')')
c          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c          WRITE(OP_STRING,
c     '      '('' 1) coordinates, coordinate, rectangular cartesian,''
c     '      ''#Components=3'')')
c          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c          WRITE(OP_STRING,'('' x.  Value index= 1, #Derivatives=0'')')
c          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c          WRITE(OP_STRING,'('' y.  Value index= 2, #Derivatives=0'')')
c          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c          WRITE(OP_STRING,'('' z.  Value index= 3, #Derivatives=0'')')
c          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          
c          DO i=1,num_points_generated
c            WRITE(OP_STRING,'('' Node:'',I5)') 1000+i
c            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c            WRITE(OP_STRING,'(3(F10.4))') X_DATA(i,1),X_DATA(i,2),
c     '        X_DATA(i,3)
c            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c          enddo
c          CLOSE(IOFILE6)

        
      CALL EXITS('CYLINDER_DATA')
      RETURN
 9999 CALL ERRORS('CYLINDER_DATA',ERROR)
      CALL EXITS('CYLINDER_DATA')
      RETURN 1
      END

