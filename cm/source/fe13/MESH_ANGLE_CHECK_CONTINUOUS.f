      SUBROUTINE MESH_ANGLE_CHECK_CONTINUOUS(np1,np2,A_LIM,LENGTH,LU,XP,
     &  XP1,ERROR,*)

C#### Subroutine: MESH_ANGLE_CHECK
C###  Description:
C###    MESH_ANGLE_CHECK checks the branching angle between a parent
C###    and daughter branch.  If the branch angle is greater than the
C###    branch angle limit, then the branch angle is reduced to the
C###    limit value, such that the daughter branch remains in the
C###    original branching plane.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
      INTEGER np1,np2
      REAL*8 A_LIM,LENGTH,LU,XP(NKM,NVM,NJM,NPM),XP1(3)
      CHARACTER ERROR*(*)
      !Local variables
      INTEGER i,j,N,nj,pivot_row
      REAL*8 A(3,4),A_MIN,angle,LV,max,N_UV(3),
     &  pivot_value,SCALAR,TEMP(4),U(3),V(3),W(3),XP2(3),XP3(3)
c      REAL*8 DOT_PROD

      CALL ENTERS('MESH_ANGLE_CHECK_CONTINUOUS',*9999)

      A_MIN=0.d0*PI/180.d0
      
C***  np1=end of parent; np2=start of parent
C*** Check branch angle to parent
      DO nj=1,3
        XP2(nj)=XP(1,1,nj,np2) !start node
        XP3(nj)=XP(1,1,nj,np1) !end node
        U(nj)=XP(1,1,nj,np1)-XP(1,1,nj,np2)
        V(nj)=XP1(nj)-XP(1,1,nj,np1)
      ENDDO !nj
      LU=DSQRT(SCALAR(3,U,U))
      LV=DSQRT(SCALAR(3,V,V)) !length of new branch
      CALL ASSERT(LU.GT.LOOSE_TOL,'>>Zero length vector',ERROR,*9999)
      CALL ASSERT(LV.GT.LOOSE_TOL,'>>Zero length vector',ERROR,*9999)
      CALL MESH_ANGLE(angle,XP2,XP3,XP1,ERROR,*9999)
      CALL NORMALISE(3,U,ERROR,*9999)
      CALL NORMALISE(3,V,ERROR,*9999)

      IF(DABS(angle).GT.A_LIM)THEN !reduce angle
C...... Solve system: n_uv.w=0,v.w=cos(angle),u.w=cos(a_lim)
        angle=angle-A_LIM
        CALL CROSS(U,V,N_UV)
        CALL NORMALISE(3,N_UV,ERROR,*9999)

        DO j=1,3
          A(1,j)=N_UV(j)
          A(2,j)=U(j)
          A(3,j)=V(j)
        ENDDO !j
        A(1,4)=0.d0
        A(2,4)=DCOS(a_lim)
        A(3,4)=DCOS(angle)

        DO N=1,2
          max=0.d0
          DO i=N,3
            IF(DABS(A(i,N)).GT.max)THEN
              max=DABS(A(i,N))
              pivot_row=i
            ENDIF
          ENDDO !i
          IF(pivot_row.NE.N)THEN !swap pivot row with row N
            DO j=1,4
              TEMP(j)=A(N,j)
              A(N,j)=A(pivot_row,j)
              A(pivot_row,j)=TEMP(j)
            ENDDO !j
          ENDIF
          pivot_value=A(N,N)
          DO j=1,4
            A(N,j)=A(N,j)/pivot_value !divide by diagonal entry
          ENDDO !j
          DO i=N+1,3
            DO j=N+1,4
              A(i,j)=A(i,j)-A(i,N)*A(N,j)
            ENDDO
            A(i,N)=0.d0
          ENDDO
        ENDDO !N

        A(3,4)=A(3,4)/A(3,3)
        A(2,4)=A(2,4)-A(3,4)*A(2,3)
        A(1,4)=A(1,4)-A(3,4)*A(1,3)-A(2,4)*A(1,2)

        DO i=1,3
          W(i)=A(i,4)
        ENDDO
        CALL NORMALISE(3,W,ERROR,*9999)
        
        DO nj=1,NJT
          XP1(nj)=XP(1,1,nj,np1)+W(nj)*LV !*0.5d0
          XP2(nj)=XP(1,1,nj,np2)
          XP3(nj)=XP(1,1,nj,np1)
        ENDDO !nj
        CALL MESH_ANGLE(angle,XP2,XP3,XP1,ERROR,*9999)

        DO i=1,3
          V(i)=W(i)
        ENDDO

      ELSE IF(DABS(angle).LT.A_MIN)THEN !increase angle
        angle=A_MIN-angle
        CALL CROSS(U,V,N_UV)
        CALL NORMALISE(3,N_UV,ERROR,*9999)
        DO j=1,3
          A(1,j)=N_UV(j)
          A(2,j)=U(j)
          A(3,j)=V(j)
        ENDDO !j
        A(1,4)=0.d0
        A(2,4)=DCOS(A_MIN)
        A(3,4)=DCOS(angle)
        DO N=1,2
          max=0.d0
          DO i=N,3
            IF(DABS(A(i,N)).GT.max)THEN
              max=DABS(A(i,N))
              pivot_row=i
            ENDIF
          ENDDO !i
          IF(pivot_row.NE.N)THEN
            DO j=1,4
              TEMP(j)=A(N,j)
              A(N,j)=A(pivot_row,j)
              A(pivot_row,j)=TEMP(j)
            ENDDO !j
          ENDIF
          pivot_value=A(N,N)
          DO j=1,4
            A(N,j)=A(N,j)/pivot_value
          ENDDO !j
          DO i=N+1,3
            DO j=N+1,4
              A(i,j)=A(i,j)-A(i,N)*A(N,j)
            ENDDO
            A(i,N)=0.d0
          ENDDO
        ENDDO !N
        A(3,4)=A(3,4)/A(3,3)
        A(2,4)=A(2,4)-A(3,4)*A(2,3)
        A(1,4)=A(1,4)-A(3,4)*A(1,3)-A(2,4)*A(1,2)
        DO i=1,3
          W(i)=A(i,4)
        ENDDO
        CALL NORMALISE(3,W,ERROR,*9999)
        DO nj=1,NJT
          XP1(nj)=XP(1,1,nj,np1)+W(nj)*LV !*0.5d0
          XP2(nj)=XP(1,1,nj,np2)
          XP3(nj)=XP(1,1,nj,np1)
        ENDDO !nj
        CALL MESH_ANGLE(angle,XP2,XP3,XP1,ERROR,*9999)

      ENDIF !DABS(angle).GT.A_LIM
      
      length=lv

      CALL EXITS('MESH_ANGLE_CHECK_CONTINUOUS')
      RETURN
 9999 CALL ERRORS('MESH_ANGLE_CHECK_CONTINUOUS',ERROR)
      CALL EXITS('MESH_ANGLE_CHECK_CONTINUOUS')
      RETURN 1
      END



