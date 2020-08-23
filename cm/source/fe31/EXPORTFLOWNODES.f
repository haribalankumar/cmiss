      SUBROUTINE EXPORTFLOWNODES(ncount,NEELEM,NXI,BBM,filename,ERROR,*)

      IMPLICIT NONE
      INCLUDE 'b00.cmn' 
      INCLUDE 'cbdi02.cmn' 
      INCLUDE 'geom00.cmn' 
      INCLUDE 'lung00.cmn' 
!     Parameter List
      INTEGER ncount,NEELEM(0:NE_R_M),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 BBM(2,NEM)
      CHARACTER ERROR*(*)

      INTEGER IBEG,IEND,ne,noelem
      REAL*8 sum_volume
      CHARACTER filename*200,FILENAME2*200

      CALL ENTERS('EXPORTFLOWNODES',*9999)

      FILENAME2=filename
      CALL STRING_TRIM(FILENAME2,IBEG,IEND)
      CALL APPENDC(IEND,'_',FILENAME2)
      CALL APPENDI(IEND,ncount,FILENAME2)
      CALL APPENDC(IEND,'.extime',FILENAME2)
 
      CALL OPENF(IOFILE3,'DISK',FILENAME2,
     '  'NEW','SEQUEN','FORMATTED',160,ERROR,*9999)

      sum_volume=0.d0
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        IF(NXI(1,0,ne).EQ.0)THEN !terminal
          sum_volume=sum_volume+BBM(1,ne)
        ENDIF
      ENDDO

      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        IF(NXI(1,0,ne).EQ.0)THEN !terminal
          WRITE(IOFILE3,'(I7,D12.4,D12.4)') ne,BBM(1,ne),
     &      sum_volume
        ENDIF
      ENDDO

      CLOSE(IOFILE3)
      
      CALL EXITS('EXPORTFLOWNODES')
      RETURN
 9999 CALL ERRORS('EXPORTFLOWNODES',ERROR)
      CALL EXITS('EXPORTFLOWNODES')
      RETURN 1
      END

