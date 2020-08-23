      SUBROUTINE ReSetInitialVolume(NEELEM,NXI,BBM,CE,sumvolume,
     &  ERROR,*)
      IMPLICIT NONE
      INCLUDE 'geom00.cmn' 
      INCLUDE 'lung00.cmn' 
      INCLUDE 'lung_nej00.cmn'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 BBM(2,NEM),CE(NMM,NEM),sumvolume
      CHARACTER ERROR*(*)

      INTEGER ne,noelem

      sumvolume=0.d0
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        IF(NXI(1,0,ne).EQ.0)THEN
          BBM(1,ne)=CE(nm_vinit,ne)
          BBM(2,ne)=0.d0
          sumvolume=sumvolume+BBM(1,ne)
        ENDIF
      ENDDO

      RETURN
 9999 CALL ERRORS('ResetInitialVolume',ERROR)
      RETURN
      END

