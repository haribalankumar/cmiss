      SUBROUTINE FLIP2D(NBJ,NEELEM,NPNE,nr,NVJE,NXI,XP,ZA,ERROR,*)

C#### Subroutine: FLIP2D
C###  Description:
C###    Performs flipping of the 2 dimensional delaunay triangulation

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter list
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NPNE(NNM,NBFM,NEM),nr,
     '  NVJE(NNM,NBFM,NJM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 XP(NKM,NVM,NJM,NPM),ZA(NAM,NHM,NCM,NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER FLIPS,ne1,ne2,nf,noelem1,nb
      LOGICAL FIN,DELAUNAY_FACE

      CALL ENTERS('FLIP2D',*9999)

      FIN=.FALSE.
      DO WHILE(.NOT.FIN)
        FLIPS=0
        DO noelem1=1,NEELEM(0,nr)
          ne1=NEELEM(noelem1,nr)
          nb=NBJ(1,ne1)
          DO nf=1,NNT(nb)
            ne2=NXI(0,nf,ne1)
            IF(ne2.GT.ne1) THEN
              IF(.NOT.
     '          DELAUNAY_FACE(NBJ,nf,ne1,ne2,NPNE,nr,NVJE,XP,ZA))
     '          THEN
                FLIPS=FLIPS+1
                CALL EDGEFLIP(NBJ,ne1,ne2,NPNE,nr,NVJE,NXI,XP,ZA,
     '            ERROR,*9999)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        IF(FLIPS.EQ.0) FIN=.TRUE.
      ENDDO

      CALL EXITS('FLIP2D')
      RETURN
 9999 CALL ERRORS('FLIP2D',ERROR)
      CALL EXITS('FLIP2D')
      RETURN 1
      END


