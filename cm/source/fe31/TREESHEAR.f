      SUBROUTINE TREESHEAR(NBJ,ncount,NEELEM,NORD,NPNE,NVJE,
     &  torr_in_order,XP,filename,ERROR,*)

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn' 
      INCLUDE 'geom00.cmn' 
      INCLUDE 'lung_nej00.cmn' 

      INTEGER NBJ(NJM,NEM),ncount,NEELEM(0:NE_R_M),NORD(5,NE_R_M),
     &  NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      CHARACTER filename*200,FILENAME2*200

      INTEGER IBEG,IEND,n,nb,ne,ngen,NMAXGEN,noelem,np1,np2,nv1,
     &  nv2,N_in_order(3,40)
      REAL*8 torr,torr_in_order(3,40)

      CALL SHEARSTRESS(NBJ,NEELEM,NPNE,NVJE,XP,ERROR,*9999)
      nj_mu=10
      
      DO ngen=1,40
        DO n=1,3
          N_in_order(n,ngen)=0
          torr_in_order(n,ngen)=0.d0
        ENDDO
      ENDDO
      NMAXGEN=0
      
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        nb=NBJ(1,ne)
        np1=NPNE(1,nb,ne) !end node
        np2=NPNE(2,nb,ne) !end node
        nv1=NVJE(1,nb,nj_radius,ne) !should be using nj_flow
        nv2=NVJE(2,nb,nj_radius,ne)
        torr=0.5d0*(XP(1,nv1,nj_mu,np1)+XP(1,nv2,nj_mu,np2))
        DO n=1,3
          ngen=NORD(n,ne) !Weibel generation, Horsfield order, Strahler order
          N_in_order(n,ngen)=N_in_order(n,ngen)+1
          torr_in_order(n,ngen)=torr_in_order(n,ngen)+torr
          IF(ngen.GT.NMAXGEN) NMAXGEN=ngen
        ENDDO
      ENDDO

      DO n=1,3
        DO ngen=1,40
          IF(N_in_order(n,ngen).GT.0)THEN
            torr_in_order(n,ngen)=torr_in_order(n,ngen)
     &        /DBLE(N_in_order(n,ngen))
          ENDIF
        ENDDO
      ENDDO

      FILENAME2=filename
      CALL STRING_TRIM(FILENAME2,IBEG,IEND)
c      CALL APPENDC(IEND,'_',FILENAME2)
c      CALL APPENDI(IEND,ncount,FILENAME2)
      CALL APPENDC(IEND,'.shear',FILENAME2)
 
      CALL OPENF(IOFILE3,'DISK',FILENAME2,
     '  'NEW','SEQUEN','FORMATTED',160,ERROR,*9999)

      DO ngen=1,NMAXGEN
        WRITE(IOFILE3,'(I4,D14.5)') ngen,torr_in_order(2,ngen)
      ENDDO
      CLOSE(IOFILE3)

      RETURN
 9999 CALL ERRORS('TREESHEAR',ERROR)
      END
      
