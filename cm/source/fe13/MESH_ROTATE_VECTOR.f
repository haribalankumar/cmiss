      SUBROUTINE MESH_ROTATE_VECTOR(N,NVP,angle,length,radius,radius0,
     &  Rx,Ry,Txyz,XPB,XPV,ERROR,*)

C#### Subroutine: MESH_ROTATE_VECTOR
C###  Description:
C###    MESH_ROTATE_VECTOR rotates a vector (starting at (0,0,0)) about
C###    an arbitrary axis. Rotation matrices are in Rx and Ry. angle
C###    is the amount to rotate about z-axis. 

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER N,NVP(5)
      REAL*8 angle,length,radius,radius0,Rx(3,3),Ry(3,3),Txyz(3),
     &  XPB(NKM,NVM,NJM,5),XPV(NKM,NVM,NJM,5)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,nj,nk,nv
      REAL*8 X1(3),X2(3)

      CALL ENTERS('MESH_ROTATE_VECTOR',*9999)
      
      DO i=1,N !each node in ring, plus crux if at bifurcation
        DO nk=1,4
          DO nv=1,NVP(i)
            DO nj=1,3
              XPV(nk,nv,nj,i)=XPB(nk,nv,nj,i)
            ENDDO !nj
            X1(1)=DCOS(angle)*XPV(nk,nv,1,i)+DSIN(angle)*XPV(nk,nv,2,i)
            X1(2)=-DSIN(angle)*XPV(nk,nv,1,i)+DCOS(angle)*XPV(nk,nv,2,i)
            XPV(nk,nv,1,i)=X1(1)
            XPV(nk,nv,2,i)=X1(2)
          ENDDO !nv
        ENDDO !nk
        DO nv=1,NVP(i)
          IF(N.EQ.5)THEN
            IF(i.EQ.2.OR.i.EQ.4)THEN
c              DO nk=1,4
              nk=1
              XPV(nk,nv,1,i)=XPV(nk,nv,1,i)*radius0
              XPV(nk,nv,2,i)=XPV(nk,nv,2,i)*radius0
c              ENDDO !nk
            ELSE
c              DO nk=1,4
              nk=1
              XPV(nk,nv,1,i)=XPV(nk,nv,1,i)*radius*1.05d0
              XPV(nk,nv,2,i)=XPV(nk,nv,2,i)*radius*1.05d0
c              ENDDO
            ENDIF
          ELSE
c            DO nk=1,4
            nk=1
            XPV(nk,nv,1,i)=XPV(nk,nv,1,i)*radius
            XPV(nk,nv,2,i)=XPV(nk,nv,2,i)*radius
c            ENDDO
          ENDIF
          IF(N.LT.5)THEN
            XPV(1,nv,3,i)=length
          ELSE IF(N.EQ.5)THEN
            IF(i.EQ.3.OR.i.EQ.1)THEN
c              XPV(1,nv,3,i)=length-radius
              XPV(1,nv,3,i)=length-radius*0.5d0
            ELSE
              XPV(1,nv,3,i)=length
            ENDIF
          ENDIF
        ENDDO !nv
        DO nk=1,4
          DO nv=1,NVP(i)
            DO nj=1,3
              X1(nj)=0.d0
              X2(nj)=0.d0
            ENDDO !nj
            DO nj=1,3
              DO j=1,3
                X1(nj)=X1(nj)+Ry(nj,j)*XPV(nk,nv,j,i)
              ENDDO !j
            ENDDO !nj
            DO nj=1,3
              DO j=1,3
                X2(nj)=X2(nj)+Rx(nj,j)*X1(j)
              ENDDO !j
            ENDDO !nj
            DO nj=1,3
              IF(nk.EQ.1)THEN !geometry
                XPV(nk,nv,nj,i)=X2(nj)-Txyz(nj)
              ELSE !derivatives
                XPV(nk,nv,nj,i)=X2(nj)
              ENDIF
            ENDDO !nj
          ENDDO !nv
        ENDDO !nk
      ENDDO !i

        
      CALL EXITS('MESH_ROTATE_VECTOR')
      RETURN
 9999 CALL ERRORS('MESH_ROTATE_VECTOR',ERROR)
      CALL EXITS('MESH_ROTATE_VECTOR')
      RETURN 1
      END


