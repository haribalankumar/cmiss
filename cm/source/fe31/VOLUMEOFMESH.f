      SUBROUTINE VOLUMEOFMESH(NBJ,NEELEM,NORD,NPNE,NVJE,NXI,MeanVolume,
     &  BBM,CE,XP,ERROR,*)

C#### Subroutine: MESH_VOLUME
C###  Description:
C###    MESH_VOLUME calculates the volume of a 1D mesh.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),NORD(5,NE_R_M),
     &  NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 MeanVolume,
     &  BBM(2,NEM),CE(NMM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nbbm,ne,ne0,ngen,nn,noelem,np1,np2
      REAL*8 generation_sum(40),length,radius
      REAL*8 LENGTH_1D,RADIUS_1D

      CALL ENTERS('VOLUMEOFMESH',*9999)
      MeanVolume=0.d0
      nbbm=0
      
      DO nn=1,40
        generation_sum(nn)=0.d0
      ENDDO
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        nb=NBJ(1,ne)
        np1=NPNE(1,nb,ne)
        np2=NPNE(2,nb,ne)
        length=LENGTH_1D(NBJ,ne,NPNE,NVJE,XP)
        radius=RADIUS_1D(NBJ,ne,NPNE,NVJE,XP)
        CE(nm_volumes,ne)=length*PI*radius**2
        CE(nm_vol_bel,ne)=CE(nm_volumes,ne)
        ngen=NORD(1,ne) !generation
        generation_sum(ngen)=generation_sum(ngen)+CE(nm_volumes,ne)
      ENDDO !noelem
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        IF(NXI(1,0,ne).EQ.0)THEN
          CE(nm_vol_bel,ne)=CE(nm_vol_bel,ne)+BBM(1,ne)
          MeanVolume=MeanVolume+BBM(1,ne)
          nbbm=nbbm+1
        ENDIF
      ENDDO

      MeanVolume=MeanVolume/DBLE(nbbm)
      
      DO noelem=NEELEM(0),1,-1
        ne=NEELEM(noelem)
        ne0=NXI(-1,1,ne)
        IF(ne0.NE.0)THEN !not stem branch
          IF(NORD(5,ne).EQ.1)THEN !start of a 'half' branch
            CE(nm_vol_bel,ne0)=CE(nm_vol_bel,ne0)+2.d0*CE(nm_vol_bel,ne)
            CE(nm_volumes,ne0)=CE(nm_volumes,ne0)+2.d0*CE(nm_volumes,ne)
          ELSE !within a tube branch
            CE(nm_vol_bel,ne0)=CE(nm_vol_bel,ne0)+CE(nm_vol_bel,ne)
            CE(nm_volumes,ne0)=CE(nm_volumes,ne0)+CE(nm_volumes,ne)
          ENDIF !NORD(5)
        ENDIF !ne0
      ENDDO !noelem

      CALL EXITS('VOLUMEOFMESH')
      RETURN
 9999 CALL ERRORS('VOLUMEOFMESH',ERROR)
      CALL EXITS('VOLUMEOFMESH')
      RETURN 1
      END
