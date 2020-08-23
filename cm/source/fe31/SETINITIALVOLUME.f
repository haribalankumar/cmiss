      SUBROUTINE SetInitialVolume(Gdirn,NBJ,NEELEM,NPNE,NVJE,NXI,BBM,CE,
     &  COV,MeanVolume,refvol,RMaxMean,RMinMean,sumvolume,undef,
     &  XP,filename,ERROR,*)
      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn' 
      INCLUDE 'lung_nej00.cmn' 
      INCLUDE 'pulm00.cmn' 
!     Parameter List
      INTEGER Gdirn,NBJ(NJM,NEM),NEELEM(0:NE_R_M),NPNE(NNM,NBFM,NEM),
     &  NVJE(NNM,NBFM,NJM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 BBM(2,NEM),CE(NMM,NEM),COV,MeanVolume,refvol,RMaxMean,
     &  RMinMean,sumvolume,undef,XP(NKM,NVM,NJM,NPM)
      CHARACTER filename*200,ERROR*(*)

      INTEGER IBEG,IEND,nb,ne,noelem,np1,np2,nv1,nv2
      REAL*8 max_z,min_z,RadiusMax,RadiusMin,random_number,range_z,Vmax,
     &  Vmin,Xi
      CHARACTER FILENAME2*200
        
C      FILENAME2=filename
C      CALL STRING_TRIM(FILENAME2,IBEG,IEND)
C      CALL APPENDC(IEND,'.vinit',FILENAME2)
C      CALL OPENF(IOFILE3,'DISK',FILENAME2,
C     ' 'NEW','SEQUEN','FORMATTED',160,ERROR,*9999)

      undef=MeanVolume*refvol
      sumvolume=0.d0
      
      Vmin = RMinMean*MeanVolume
      Vmax = RMaxMean*MeanVolume

      RadiusMax=(RMaxMean**0.33d0)*1.d0
      RadiusMin=(RMinMean**0.33d0)*1.d0
      
      max_z=-1.d6
      min_z=1.d6
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        IF(NXI(1,0,ne).EQ.0)THEN
          np2=NPNE(2,1,ne)
          max_z=MAX(max_z,XP(1,1,Gdirn,np2))
          min_z=MIN(min_z,XP(1,1,Gdirn,np2))
        ENDIF
      ENDDO
      range_z=DABS(max_z-min_z)
      IF(DABS(range_z).LE.1d-5) range_z=1.d0

      random_number=-1.1d0
      
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        nb=NBJ(1,ne)
        np1=NPNE(1,1,ne)
        np2=NPNE(2,1,ne)
        nv1=NVJE(1,nb,nj_radius,ne)
        nv2=NVJE(2,nb,nj_radius,ne)
        Xi=(XP(1,1,Gdirn,np2)-min_z)/range_z
        random_number=random_number+0.1d0
        IF(random_number.GT.1.d0) random_number=-1.1d0
        IF(NXI(1,0,ne).EQ.0)THEN
          BBM(1,ne)=Vmax*Xi+Vmin*(1.d0-Xi) !in mm^3
          BBM(1,ne)=BBM(1,ne)*(1.d0+COV*random_number)
          CE(nm_vinit,ne)=BBM(1,ne)
          BBM(2,ne)=0.d0
          sumvolume=sumvolume+BBM(1,ne)
C          WRITE(IOFILE3,'(I7,D12.4,D12.4,D12.4)') ne,BBM(1,ne),
C     &      XP(1,1,Gdirn,np2),Xi
        ENDIF
c        IF(ne.GT.7)THEN
c         CE(nm_lambda,ne)=((Vmax*Xi+Vmin*(1.d0-Xi))/undef)
c     &      **(1.d0/3.d0) !extn ratio
c          XP(2,nv1,nj_radius,np1)=XP(1,nv1,nj_radius,np1)/
c     &      ((1.d0/refvol)**(1.d0/3.d0))
c          XP(1,nv1,nj_radius,np1)=XP(1,nv1,nj_radius,np1)
c     &      *CE(nm_lambda,ne)/((1.d0/refvol)**(1.d0/3.d0))
c          XP(2,nv2,nj_radius,np2)=XP(1,nv2,nj_radius,np2)/
c     &      ((1.d0/refvol)**(1.d0/3.d0))
c          XP(1,nv2,nj_radius,np2)=XP(1,nv2,nj_radius,np2)
c     &      *CE(nm_lambda,ne)/((1.d0/refvol)**(1.d0/3.d0))
c        ENDIF

      ENDDO

C      CLOSE(IOFILE3)

      RETURN
 9999 CALL ERRORS('SetInitialVolume',ERROR)
      RETURN
      END

