      SUBROUTINE MESH_ANGLE_CHECK(Nth,np1,np2,np3,np,A_LIM,LENGTH,LU,XP,
     &  XP1,ERROR,*)

C#### Subroutine: MESH_ANGLE_CHECK
C###  Description:
C###    MESH_ANGLE_CHECK checks the branching angle between a parent
C###    and daughter branch.  If the branch angle is greater than the
C###    branch angle limit, then the branch angle is reduced to the
C###    limit value, such that the daughter branch remains in the
C###    original branching plane.

C called from GNBDMESH, MESH_BRANCH

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
      INTEGER Nth,np1,np2,np3,np
      REAL*8 A_LIM,LENGTH,LU,XP(NKM,NVM,NJM,NPM),XP1(3)
      CHARACTER ERROR*(*)
      !Local variables
      INTEGER i,j,N,nj,pivot_row
      REAL*8 A(3,4),A_MIN,angle,angle0,LV,max,N_UV(3),
     &  pivot_value,SCALAR,TEMP,U(3),UU(3),V(3),VECTOR(3),W(3),
     &  XP2(3),XP3(3)
c      REAL*8 DOT_PROD

      CALL ENTERS('MESH_ANGLE_CHECK',*9999)

      A_MIN=25.d0*PI/180.d0
      
C***  np1=end of parent; np2=start of parent
C*** Check branch angle to parent
      DO nj=1,3
        XP2(nj)=XP(1,1,nj,np2) !start node
        XP3(nj)=XP(1,1,nj,np1) !end node
        U(nj)=XP(1,1,nj,np1)-XP(1,1,nj,np2) !direction of parent
        V(nj)=XP1(nj)-XP(1,1,nj,np1) !direction of this branch
        W(nj)=V(nj) !will store direction if no angle change
      ENDDO !nj
      LU=DSQRT(SCALAR(3,U,U)) !length of the parent
      LV=DSQRT(SCALAR(3,V,V)) !length of new branch
      CALL ASSERT(LU.GT.LOOSE_TOL,'>>Zero length vector LU',ERROR,*9999)
      CALL ASSERT(LV.GT.LOOSE_TOL,'>>Zero length vector LV',ERROR,*9999)
      CALL MESH_ANGLE(angle,XP2,XP3,XP1,ERROR,*9999)
      CALL NORMALISE2(3,U,TEMP,ERROR,*9999)
      CALL NORMALISE2(3,V,TEMP,ERROR,*9999)
      CALL NORMALISE2(3,W,TEMP,ERROR,*9999)

      IF(DABS(angle).GT.A_LIM)THEN !reduce angle
        angle0=angle !original angle (too large)
        angle=angle-A_LIM !amount of angle to 'remove'
        CALL CROSS(U,V,N_UV)
        CALL NORMALISE2(3,N_UV,TEMP,ERROR,*9999)
C...... n.w = 0
C...... u.w = cos(angle_u)
C...... v.w = cos(angle_v)
        DO nj=1,3
          A(1,nj)=N_UV(nj)
          A(2,nj)=U(nj)
          A(3,nj)=V(nj)
        ENDDO
        VECTOR(1)=0.d0
        VECTOR(2)=DCOS(A_LIM)
        VECTOR(3)=DCOS(angle)

        CALL MESH_A_X_EQ_B(A,VECTOR,W,ERROR,*9999)
        CALL NORMALISE2(3,W,TEMP,ERROR,*9999)

        LV=MIN(LV,LU*1.25d0)
        
        DO nj=1,NJT
          XP1(nj)=XP(1,1,nj,np1)+W(nj)*LV !*0.5d0
          XP2(nj)=XP(1,1,nj,np2)
          XP3(nj)=XP(1,1,nj,np1)
        ENDDO !nj
        CALL MESH_ANGLE(angle,XP2,XP3,XP1,ERROR,*9999)

        DO i=1,3
          V(i)=W(i)
        ENDDO

      ELSEIF(DABS(angle).LT.A_MIN)THEN
        IF(DABS(angle).GT.LOOSE_TOL)THEN
C increase the branch angle to the limit value
          angle0=angle !original angle
          angle=A_MIN-angle0 !amount of angle to 'add'
          CALL CROSS(U,V,N_UV)
          CALL NORMALISE(3,N_UV,ERROR,*9999)
C...... n.w = 0
C...... u.w = cos(angle_u)
C...... v.w = cos(angle_v)
          DO nj=1,3
            A(1,nj)=N_UV(nj)
            A(2,nj)=U(nj)
            A(3,nj)=V(nj)
          ENDDO
          
          VECTOR(1)=0.d0
          VECTOR(2)=DCOS(A_MIN)
          VECTOR(3)=DCOS(angle)
          
          CALL MESH_A_X_EQ_B(A,VECTOR,W,ERROR,*9999)
          CALL NORMALISE(3,W,ERROR,*9999)
          
          LV=MIN(LV,LU*1.25d0)
          
          DO nj=1,NJT
            XP1(nj)=XP(1,1,nj,np1)+W(nj)*LV !*0.5d0
            XP2(nj)=XP(1,1,nj,np2)
            XP3(nj)=XP(1,1,nj,np1)
          ENDDO !nj
          CALL MESH_ANGLE(angle,XP2,XP3,XP1,ERROR,*9999)
c          write(*,*) 'increase angle',np1,np2,angle0*180.d0/PI,angle
c     &      *180.d0/PI
c        pause
          DO i=1,3
            V(i)=W(i)
          ENDDO
          
        ELSE IF(Nth.EQ.2)THEN
          DO nj=1,3
            U(nj)=XP(1,2,nj,np-1) !direction of sibling
          ENDDO !nj
C increase the branch angle to the limit value
          ANGLE=SCALAR(3,U,V) !angle between branch and sibling
          ANGLE=MIN(1.d0,MAX(-1.d0,ANGLE))
          ANGLE=DACOS(ANGLE)

          CALL CROSS(U,V,N_UV)
          CALL NORMALISE(3,N_UV,ERROR,*9999)
C...... n.w = 0
C...... u.w = cos(angle_u+angle_min)
C...... v.w = cos(angle_min)
          DO nj=1,3
            A(1,nj)=N_UV(nj) !normal
            A(2,nj)=U(nj) !sibling
            A(3,nj)=V(nj) !itself
          ENDDO
          
          VECTOR(1)=0.d0
          VECTOR(2)=DCOS(ANGLE+A_MIN)
          VECTOR(3)=DCOS(A_MIN)
          
          CALL MESH_A_X_EQ_B(A,VECTOR,W,ERROR,*9999)
          CALL NORMALISE(3,W,ERROR,*9999)
          
          DO nj=1,NJT
            XP1(nj)=XP(1,1,nj,np1)+W(nj)*LV !*0.5d0
            XP2(nj)=XP(1,1,nj,np2)
            XP3(nj)=XP(1,1,nj,np1)
          ENDDO !nj
          CALL MESH_ANGLE(angle,XP2,XP3,XP1,ERROR,*9999)
c          write(*,*) 'increase angle',np1,np2,angle0*180.d0/PI,angle
c     &      *180.d0/PI
c        pause
          DO i=1,3
            V(i)=W(i)
          ENDDO
          
        ELSE IF(Nth.EQ.1)THEN
          WRITE(*,*) 'WARNING!!!!! Zero branch angle for ',np
        ENDIF
        
      ENDIF !DABS(angle).GT.A_LIM

      length=LV
C.....Check that not less than 0.75 of parent length
      IF(LV.LT.0.5d0*LU)THEN
        DO nj=1,NJT
          XP1(nj)=XP(1,1,nj,np1)+W(nj)*LU*0.5d0
        ENDDO
        length=LU*0.5d0
c        write(*,*) 'increased length from ',LV,' to ',length
c        pause
      ENDIF
      

      CALL EXITS('MESH_ANGLE_CHECK')
      RETURN
 9999 CALL ERRORS('MESH_ANGLE_CHECK',ERROR)
      CALL EXITS('MESH_ANGLE_CHECK')
      RETURN 1
      END



