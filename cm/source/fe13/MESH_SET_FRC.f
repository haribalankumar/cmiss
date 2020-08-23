      SUBROUTINE MESH_SET_FRC(NBJ,NEELEM,NELIST2,NORD,NPNE,NVJE,NXI,
     &  FIXED_VOLUME,FRC,volume_below,volumes,XAB,XP,LUMPED_PARAMETER,
     &  SCALE,ERROR,*)

C#### Subroutine: MESH_SET_FRC
C###  Description:
C###    MESH_SET_FRC 

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),NELIST2(0:NEM),
     &  NORD(5,NE_R_M),NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 FIXED_VOLUME,FRC,volume_below(NE_R_M),volumes(NE_R_M),
     &  XAB(NORM,NEM),XP(NKM,NVM,NJM,NPM)
      LOGICAL LUMPED_PARAMETER
      CHARACTER SCALE*(10),ERROR*(*)
!     Local Variables
      INTEGER nb,ne,nlpm,nj,nn,noelem,npnn(2),nvnn(2)
      REAL*8 direction(3),length,radius,RATIO_SCALE

      CALL ENTERS('MESH_SET_FRC',*9999)


      CALL ASSERT(FRC.GT.FIXED_VOLUME,
     &  '>>FRC is less than the fixed volume',ERROR,*9999)
      RATIO_SCALE=(FRC-FIXED_VOLUME)/(volume_below(NEELEM(1))
     &  -FIXED_VOLUME)
      IF(LUMPED_PARAMETER)THEN !only scale the lumped parameter models
        DO nlpm=1,NTB
          XAB(2,nlpm)=XAB(2,nlpm)*RATIO_SCALE
          XAB(3,nlpm)=XAB(2,nlpm)
        ENDDO !nlpm 
      ELSE !scale the lumped parameter models and the list of elements
        DO nlpm=1,NTB
          XAB(2,nlpm)=XAB(2,nlpm)*RATIO_SCALE
          XAB(3,nlpm)=XAB(2,nlpm)
        ENDDO !nlpm 
        DO noelem=1,NELIST2(0)
          ne=NELIST2(noelem)
          nb=NBJ(1,ne)
C.........Calculate the initial element length
          DO nn=1,NNT(nb)
            npnn(nn)=NPNE(nn,nb,ne)
            nvnn(nn)=NVJE(nn,nb,1,ne)
          ENDDO !nn
          length=0.d0
          DO nj=1,NJT
            length=length+(XP(1,nvnn(2),nj,npnn(2))-XP(1,nvnn(1),nj,
     &        npnn(1)))**2
          ENDDO !nj
          length=DSQRT(length)
C.........Calculate the initial element radius
          nb=NBJ(nj_radius,ne) !basis function for radii
          DO nn=1,NNT(nb) !versions for radii at nodes
            nvnn(nn)=NVJE(nn,nb,nj_radius,ne)
          ENDDO !nn
          radius=0.5d0*(XP(1,nvnn(1),nj_radius,npnn(1))+
     &      XP(1,nvnn(2),nj_radius,npnn(2)))
C.........Calculate the scaled element volume
          volumes(ne)=volumes(ne)*RATIO_SCALE
          
          IF(SCALE(1:8).EQ.'DIAMETER')THEN !scale only the radii
c              write(*,*)'scale radii only'
            radius=DSQRT(volumes(ne)/(length*PI))
            nb=NBJ(nj_radius,ne) !basis function for radii
            DO nn=1,NNT(nb) !versions for radii at nodes
              nvnn(nn)=NVJE(nn,nb,nj_radius,ne)
            ENDDO !nn
            XP(1,nvnn(1),nj_radius,npnn(1))=radius
            XP(1,nvnn(2),nj_radius,npnn(2))=radius
          ELSE IF(SCALE(1:6).EQ.'LENGTH')THEN !scale only the lengths
c              write(*,*)'scale length only'
            length=volumes(ne)/(PI*radius**2)
            DO nj=1,NJT
              direction(nj)=XP(1,nvnn(2),nj,npnn(2))-
     &          XP(1,nvnn(1),nj,npnn(1))
            ENDDO !nj
            CALL NORMALISE(3,direction,ERROR,*9999)
            DO nj=1,NJT
              XP(1,nvnn(2),nj,npnn(2))=XP(1,nvnn(1),nj,npnn(1))+
     &          length*direction(nj)
            ENDDO !nj
          ELSE IF(SCALE(1:6).EQ.'VOLUME')THEN !scale lengths and radii
c              write(*,*)'scale volume only'
            length=(volumes(ne)*(1.d0-1.d0/RATIO_SCALE))**(1.d0/3.d0)
            DO nj=1,NJT
              direction(nj)=XP(1,nvnn(2),nj,npnn(2))-
     &          XP(1,nvnn(1),nj,npnn(1))
            ENDDO !nj
            CALL NORMALISE(3,direction,ERROR,*9999)
            DO nj=1,NJT
              XP(1,nvnn(2),nj,npnn(2))=XP(1,nvnn(1),nj,npnn(1))+
     &          length*direction(nj)
            ENDDO !nj
            radius=DSQRT(volumes(ne)/(length*PI))
            nb=NBJ(nj_radius,ne) !basis function for radii
            DO nn=1,NNT(nb) !versions for radii at nodes
              nvnn(nn)=NVJE(nn,nb,nj_radius,ne)
            ENDDO !nn
            XP(1,nvnn(1),nj_radius,npnn(1))=radius
            XP(1,nvnn(2),nj_radius,npnn(2))=radius
          ENDIF
        ENDDO
      ENDIF
C.....Recalculate volumes
      CALL MESH_VOLUME(NBJ,NEELEM,NORD,NPNE,NVJE,NXI,volumes,
     &  volume_below,XAB,XP,ERROR,*9999)
      
      CALL EXITS('MESH_SET_FRC')
      RETURN
 9999 CALL ERRORS('MESH_SET_FRC',ERROR)
      CALL EXITS('MESH_SET_FRC')
      RETURN 1
      END


