      SUBROUTINE MESH_VOLUME(NBJ,NEELEM,NORD,NPNE,NVJE,NXI,
     &  BBM,CE,XP,ERROR,*)

C#### Subroutine: MESH_VOLUME
C###  Description:
C###    MESH_VOLUME calculates the volume of a 1D mesh.

C CE(2) is volume_dv
C CE(3) is volume of element      
C CE(4) is volume below the element

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),NORD(5,NE_R_M),
     &  NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 BBM(2,NEM),CE(NMM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ne,ne0,nj,nlpm,nn,noelem,np1,npnn(2),nvnn(2)
      REAL*8 length,radius
      REAL*8 LENGTH_1D,RADIUS_1D
      REAL*8 sum_volume

      CALL ENTERS('MESH_VOLUME',*9999)

      sum_volume=0.d0
      
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        length=LENGTH_1D(NBJ,ne,NPNE,NVJE,XP)
c        radius=RADIUS_1D(NBJ,ne,NPNE,NVJE,XP)
        np1=NPNE(2,NBJ(nj_radius,ne),ne)
        radius=XP(1,1,nj_radius,np1)
         write(*,*) 'noelem, radius',noelem,radius
        CE(3,ne)=length*PI*radius**2
        CE(4,ne)=CE(3,ne)
        IF(NXI(1,0,ne).EQ.0)THEN !terminal
          CE(4,ne)=CE(4,ne)+BBM(1,ne) !add BBM volume
          sum_volume=sum_volume+BBM(1,ne)
        ENDIF
      ENDDO !noelem

      DO noelem=NEELEM(0),1,-1
        ne=NEELEM(noelem)
        ne0=NXI(-1,1,ne)
        IF(ne0.NE.0)THEN !not stem branch
          IF(NORD(5,ne).EQ.1)THEN !start of a 'half' branch
            CE(4,ne0)=CE(4,ne0)+2.d0*CE(4,ne)
          ELSE !within a tube branch
            CE(4,ne0)=CE(4,ne0)+CE(4,ne)
          ENDIF !NORD(5)
        ELSE
C          CURRENT_VOLUME=CE(4,ne) !volume below the stem branch
        ENDIF !ne0
      ENDDO !noelem

         write(*,*) 'noelem, radius END SUB'

      CALL EXITS('MESH_VOLUME')
      RETURN
 9999 CALL ERRORS('MESH_VOLUME',ERROR)
      CALL EXITS('MESH_VOLUME')
      RETURN 1
      END


