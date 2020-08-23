      SUBROUTINE MESH_SMOOTH(NBJ,NELIST,NPNE,NVJE,NXI,ANGLE_LIMIT,
     &  LDMinimum,LENGTH_LIMIT,LENGTH_MAX,XP,ERROR,*)

C#### Subroutine: MESH_SMOOTH
C###  Description:

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'

      INTEGER NBJ(NJM,NEM),NELIST(0:NEM),NPNE(NNM,NBFM,NEM),
     &  NVJE(NNM,NBFM,NJM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 ANGLE_LIMIT,LDMinimum,LENGTH_LIMIT,LENGTH_MAX,
     &  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      !Local variables
      INTEGER M,N,nb,ne,ne0,ne1,ne2,nj,noelem,np,np0,np02,np1,np2,
     &  NE_OLD(1000),NE_TEMP(1000),NT_BNS,N_ELM
      REAL*8 A(3,3),ANGLE0,ANGLE,ANGLE2,diameter,length,N_UV(3),u(3),
     &  v(3),VECTOR(3),w(3),Y(3)
      REAL*8 LENGTH_1D,RADIUS_1D,SCALAR

      CALL ENTERS('MESH_SMOOTH',*9999)

      DO noelem=1,NELIST(0)
        ne=NELIST(noelem)
        nb=NBJ(1,ne)
        IF(NXI(1,0,ne).EQ.0)THEN !terminal
          np1= NPNE(1,nb,ne)
          np2= NPNE(2,nb,ne)
          ne0=NXI(-1,1,ne) !parent
          np0= NPNE(1,nb,ne0)
          DO nj=1,NJT
            u(nj)=XP(1,1,nj,np1)-XP(1,1,nj,np0) !parent vector
            v(nj)=XP(1,1,nj,np2) -XP(1,1,nj,np1) !child vector
          ENDDO !nj
          length=DSQRT(SCALAR(3,V,V)) !length of branch
          CALL NORMALISE(3,U,ERROR,*9999)
          CALL NORMALISE(3,V,ERROR,*9999)
C.........Make length be a minimum limit          
          IF(length.LT.LENGTH_LIMIT+LOOSE_TOL)THEN
            length=LENGTH_LIMIT
            DO nj=1,NJT
              XP(1,1,nj,np2)=XP(1,1,nj,np1)+v(nj)*length
            ENDDO
          ENDIF
C.........Make length be a maximum limit          
          IF(length.GT.LENGTH_MAX+LOOSE_TOL)THEN
            length=LENGTH_MAX
            DO nj=1,NJT
              XP(1,1,nj,np2)=XP(1,1,nj,np1)+v(nj)*length
            ENDDO
          ENDIF

          ANGLE=SCALAR(3,U,V) !branch angle
          ANGLE=MAX(-1.d0,ANGLE)
          ANGLE=MIN(1.d0,ANGLE)
          ANGLE=DACOS(DABS(ANGLE))
c          IF(ANGLE.LT.ANGLE_LIMIT-0.1d0)THEN !increase the angle
c            write(*,*) 'WARNING!!! small terminal angle = ',angle*180.d0
c     &        /PI,' at branch ',ne
c          ENDIF
        ENDIF
      ENDDO !noelem
        
      CALL EXITS('MESH_SMOOTH')
      RETURN
 9999 CALL ERRORS('MESH_SMOOTH',ERROR)
      CALL EXITS('MESH_SMOOTH')
      RETURN 1
      END

