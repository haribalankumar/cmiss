      SUBROUTINE LIMIT_ANGLE(np1,np2,np3,np,XP,ERROR,*)

C#### Subroutine: LIMIT_ANGLE
C###  Description:

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
      INTEGER np1,np2,np3,np
      REAL*8 LENGTH,LU,XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      !Local variables
      INTEGER i,j,N,nj
      REAL*8 A(3,4),A_LIM,angle,angle0,LV,max,N_UV(3),
     &  SCALAR,TEMP,U(3),UU(3),V(3),VECTOR(3),W(3),XP1(3),
     &  XP2(3),XP3(3)

      CALL ENTERS('LIMIT_ANGLE',*9999)

C***  np1=end of parent; np2=start of parent
C*** Check branch angle to parent
      DO nj=1,3
        XP1(nj)=XP(1,1,nj,np)
        XP2(nj)=XP(1,1,nj,np2) !start node
        XP3(nj)=XP(1,1,nj,np1) !end node
        U(nj)=XP(1,1,nj,np1)-XP(1,1,nj,np2) !direction of parent
        V(nj)=XP1(nj)-XP(1,1,nj,np1) !direction of this branch
        W(nj)=V(nj) !will store direction if no angle change
      ENDDO !nj
      LV=DSQRT(SCALAR(3,V,V)) !length of new branch
      CALL ASSERT(LV.GT.LOOSE_TOL,'>>Zero length vector LV',ERROR,*9999)
      CALL MESH_ANGLE(angle,XP2,XP3,XP1,ERROR,*9999)
      CALL NORMALISE2(3,U,TEMP,ERROR,*9999)
      CALL NORMALISE2(3,V,TEMP,ERROR,*9999)
      CALL NORMALISE2(3,W,TEMP,ERROR,*9999)

      A_LIM=angle*0.5d0
      angle0=angle !original angle (too large)
      angle=0.5d0*angle !amount of angle to 'remove'
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
      
      DO nj=1,NJT
        XP(1,1,nj,np)=XP(1,1,nj,np1)+W(nj)*LV !*0.5d0
        XP1(nj)=XP(1,1,nj,np)
      ENDDO !nj
      CALL MESH_ANGLE(angle,XP2,XP3,XP1,ERROR,*9999)

c      write(*,*)'reduced angle from ',angle0*180.d0/PI,' to ',angle
c     &  *180.d0/PI

      CALL EXITS('LIMIT_ANGLE')
      RETURN
 9999 CALL ERRORS('LIMIT_ANGLE',ERROR)
      CALL EXITS('LIMIT_ANGLE')
      RETURN 1
      END



