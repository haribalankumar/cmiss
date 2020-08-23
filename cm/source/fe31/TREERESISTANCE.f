      SUBROUTINE TREERESISTANCE(NEELEM,NORD,CE,filename,ERROR,*)

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn' 
      INCLUDE 'geom00.cmn' 
      INCLUDE 'lung00.cmn' 
      INCLUDE 'lung_nej00.cmn' 

      INTEGER NEELEM(0:NE_R_M),NORD(5,NE_R_M)
      REAL*8 CE(NMM,NEM)
      CHARACTER ERROR*(*)
      CHARACTER filename*200,FILENAME2*200

      INTEGER IBEG,IEND,n,ne,ngen,NMAXGEN,noelem,NN_O(40)
      REAL*8 OO_O(40),RR_O(40)

      CALL ENTERS('TREERESISTANCE',*9999)

      DO ngen=1,40
        NN_O(ngen)=0
        RR_O(ngen)=0.d0
        OO_O(ngen)=0.d0
      ENDDO
      NMAXGEN=0
      
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
c        DO n=1,3
          ngen=NORD(2,ne) !Weibel generation
          NN_O(ngen)=NN_O(ngen)+1
          RR_O(ngen)=RR_O(ngen)+CE(nm_Rt,ne)
          OO_O(ngen)=OO_O(ngen)+1.d0/CE(nm_Rt,ne) !like a parallel sum
          IF(ngen.GT.NMAXGEN) NMAXGEN=ngen
c        ENDDO
      ENDDO

c      DO n=1,3
        DO ngen=1,40
          IF(NN_O(ngen).GT.0)THEN
            RR_O(ngen)=RR_O(ngen)/DBLE(NN_O(ngen))
          ENDIF
        ENDDO
c      ENDDO

      FILENAME2=filename
      CALL STRING_TRIM(FILENAME2,IBEG,IEND)
      CALL APPENDC(IEND,'.resist',FILENAME2)
      CALL OPENF(IOFILE3,'DISK',FILENAME2,
     '  'NEW','SEQUEN','FORMATTED',160,ERROR,*9999)

c      WRITE(IOFILE3,*) 'Order  MeanRbyGen SumRbyGen MeanRbyHord SumRbyHord'
      DO ngen=1,NMAXGEN
        WRITE(IOFILE3,'(I4,2(D14.5))') ngen,
     &    RR_O(ngen),1.d0/OO_O(ngen)
      ENDDO
      CLOSE(IOFILE3)

      CALL EXITS('TREERESISTANCE')
      RETURN
 9999 CALL ERRORS('TREERESISTANCE',ERROR)
      CALL EXITS('TREERESISTANCE')
      RETURN 1
      END
      
