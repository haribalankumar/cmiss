      SUBROUTINE CALC_FACE_INFORMATION_DEP(NBH,NBHF,nef,NHE,NKHE,
     '  NKEF,NKHF,NNF,NPNE,NPNF,NVHE,NVHF,nx,SE,SF,ERROR,*)

C#### Subroutine: CALC_FACE_INFORMATION_DEP
C###  Description:
C###    CALC_FACE_INFORMATION_DEP calculates the dependent variable
C###    face information NPNF, NVHF, NKHF, and SF for local face nef of
C###    global element ne.

C#### Variable: NKHF(nk,nn,nh)
C###  Type: INTEGER
C###  Set_up: CALC_FACE_INFORMATION_DEP
C###  Description:
C###    NKHF(nk,nn,nh) is the global derivative number of local
C###    derivative nk of node nn for dependent variable nh.  This is a
C###    temporary array set up for one face at a time.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
!     Parameter List
      INTEGER NBH(NHM),NBHF(NHM),nef,NHE,NKHE(NKM,NNM,NHM),
     '  NKEF(0:4,16,6,NBFM),NKHF(NKM,NNM,NHM),NNF(0:17,6,NBFM),
     '  NPNE(NNM,NBFM),NPNF(NNM,NBFM),NVHE(NNM,NBFM,NHM),
     '  NVHF(NNM,NBFM,NHM),nx
      REAL*8 SE(NSM,NBFM),SF(NSM,NBFM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nbe,nbface,nh,nhx,nke,nkface,nn2,nne,nnface,nsface,nse

      CALL ENTERS('CALC_FACE_INFORMATION_DEP',*9999)

      DO nhx=1,NHE
        nh=NH_LOC(nhx,nx)
        nbe=NBH(nh)
        nbface=NBHF(nh)
        nsface=0
C!!! NNF is not always initialized in FACSEG.
C!!! e.g. if the basis function has no nodes.
        DO nnface=1,NNF(0,nef,nbe)
          nne=NNF(1+nnface,nef,nbe)
          NPNF(nnface,nbface)=NPNE(nne,nbe)
          nse=0
          DO nn2=1,nne-1
            nse=nse+NKT(nn2,nbe)
          ENDDO !nne
          DO nkface=1,NKEF(0,nnface,nef,nbe)
            nke=NKEF(nkface,nnface,nef,nbe)
            nsface=nsface+1
            NKHF(nkface,nnface,nh)=NKHE(nke,nne,nh)
            SF(nsface,nbface)=SE(nse+nke,nbe)
          ENDDO !nk
          NVHF(nnface,nbface,nh)=NVHE(nne,nbe,nh)
        ENDDO !nn
      ENDDO !nhx

      CALL EXITS('CALC_FACE_INFORMATION_DEP')
      RETURN
 9999 CALL ERRORS('CALC_FACE_INFORMATION_DEP',ERROR)
      CALL EXITS('CALC_FACE_INFORMATION_DEP')
      RETURN 1
      END


