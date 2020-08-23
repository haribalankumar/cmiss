      SUBROUTINE MESH_ROTATION(nb,ne,NML,NPNE,NXI,angle_z,Rx,Ry,
     &  Txyz,XP,XPB,ERROR,*)

C#### Subroutine: MESH_ROTATION
C###  Description:
C###    MESH_ROTATION calculates the rotation matrices and z-angle for
C###    rotation of a vector about an arbitrary axis.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER nb,ne,NPNE(NNM,NBFM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 angle_z,Rx(3,3),Ry(3,3),Txyz(3),
     &  XP(NKM,NVM,NJM,NPM),XPB(NKM,NVM,NJM,5)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,ne0,nj,np0,np1,np2,np3,np4
      REAL*8 angle_z2,length,length_23,NML(4),NMLXP(3),V(3),X1(3),X2(3),
     &  X3(3)
      LOGICAL FOUND
      REAL*8 DOT_PROD

      CALL ENTERS('MESH_ROTATION',*9999)

      np1=NPNE(1,nb,ne) !start of element
      np2=NPNE(2,nb,ne) !end of element
C.....Find next bifurcation nodes, for calculating rotation about the 
C.....z-axis
      IF(NXI(1,0,ne).GE.2)THEN !get adjacent nodes
        np0=np1
        np3=NPNE(2,nb,NXI(1,1,ne)) !end node of first child
        np4=NPNE(2,nb,NXI(1,2,ne)) !end node of second child
      ELSE !find next bifurcation
        FOUND=.FALSE.
        ne0=ne
        DO WHILE(.NOT.FOUND)
          IF(NXI(1,0,ne0).EQ.1)THEN
            ne0=NXI(1,1,ne0)
            np0=NPNE(1,nb,ne0)
          ELSE IF(NXI(1,0,ne0).GE.2)THEN
            FOUND=.TRUE.
            np3=NPNE(2,nb,NXI(1,1,ne0))
            np4=NPNE(2,nb,NXI(1,2,ne0))
          ELSE IF(NXI(1,0,ne0).EQ.0)THEN
            FOUND=.TRUE.
            np3=0
            np4=0
          ENDIF
        ENDDO
      ENDIF
C.....Calculate the rotation matrices      
      length=0.d0
      DO nj=1,NJT
        Txyz(nj)=-XP(1,1,nj,np1)
        NMLXP(nj)=XP(1,1,nj,np2)-XP(1,1,nj,np1)
        length=length+NMLXP(nj)**2
      ENDDO !nj
      CALL NORMALISE(3,NMLXP,ERROR,*9999)
      length=DSQRT(length)
      DO i=1,3
        Rx(1,i)=0.d0 !x row = 0
        Rx(i,1)=0.d0 !x column = 0
        Ry(2,i)=0.d0 !x row = 0
        Ry(i,2)=0.d0 !x column = 0
      ENDDO !i
      Rx(1,1)=1.d0 !x-x = 1
      Ry(2,2)=1.d0 !x-x = 1
      length_23=DSQRT(NMLXP(2)**2+NMLXP(3)**2)
      IF(DABS(length_23).LT.ZERO_TOL) length_23=1.d0
      Rx(2,2)= NMLXP(3)/length_23
      Rx(2,3)= NMLXP(2)/length_23
      Rx(3,2)=-Rx(2,3)
      Rx(3,3)= Rx(2,2)

      Ry(1,1)= length_23
      Ry(1,3)= NMLXP(1)
      Ry(3,1)=-Ry(1,3)
      Ry(3,3)= Ry(1,1)

C.....The angle for rotation about the z-axis is equal to the angle
C.....between the normal to the plane containing bifurcation nodes and
C.....the theta=0 direction
      IF(np3.NE.0)THEN
        DO nj=1,3
          X1(nj)=XP(1,1,nj,np0)
          X2(nj)=XP(1,1,nj,np3)
          X3(nj)=XP(1,1,nj,np4)
        ENDDO !nj
C Note that had error here due to calculation of nml(4)
        CALL PLANE_FROM_3_PTS(NML,2,X1,X2,X3,ERROR,*9999)
        DO nj=1,3
          V(nj)=XPB(1,1,nj,2)
        ENDDO !nj
        CALL NORMALISE(3,V,ERROR,*9999)
C.......Calculate location of V if rotation about x- and y- was applied
        DO nj=1,3
          X1(nj)=0.d0
          X2(nj)=0.d0
        ENDDO !nj
        DO nj=1,3
          DO j=1,3
            X1(nj)=X1(nj)+Ry(nj,j)*V(j)
          ENDDO !j
        ENDDO !nj
        DO nj=1,3
          DO j=1,3
            X2(nj)=X2(nj)+Rx(nj,j)*X1(j)
          ENDDO !j
        ENDDO !nj
        DO nj=1,3
          V(nj)=X2(nj)!-Txyz(nj)
        ENDDO
        CALL NORMALISE(3,V,ERROR,*9999)
        angle_z=MIN(DOT_PROD(V,NML),1.d0)
        angle_z=MAX(DOT_PROD(V,NML),-1.d0)
        IF(angle_z.GE.1.d0)THEN
          angle_z=0.d0
        ELSE IF(angle_z.LE.-1.d0)THEN
          angle_z=PI
        ELSE
          angle_z=DACOS(angle_z) !between normal and 2nd node
        ENDIF
        DO nj=1,3
          V(nj)=XPB(1,1,nj,3)
        ENDDO !nj
        CALL NORMALISE(3,V,ERROR,*9999)
C.......Calculate location of V if rotation about x- and y- was applied
        DO nj=1,3
          X1(nj)=0.d0
          X2(nj)=0.d0
        ENDDO !nj
        DO nj=1,3
          DO j=1,3
            X1(nj)=X1(nj)+Ry(nj,j)*V(j)
          ENDDO !j
        ENDDO !nj
        DO nj=1,3
          DO j=1,3
            X2(nj)=X2(nj)+Rx(nj,j)*X1(j)
          ENDDO !j
        ENDDO !nj
        DO nj=1,3
          V(nj)=X2(nj)!-Txyz(nj)
        ENDDO
        CALL NORMALISE(3,V,ERROR,*9999)
        angle_z2=DACOS(DOT_PROD(V,NML)) !between normal and 2nd node
        IF(angle_z2.LT.PI/2.d0)THEN
          angle_z=-angle_z
        ENDIF
      ELSE
        angle_z=0.d0
      ENDIF
      
      CALL EXITS('MESH_ROTATION')
      RETURN
 9999 CALL ERRORS('MESH_ROTATION',ERROR)
      CALL EXITS('MESH_ROTATION')
      RETURN 1
      END
