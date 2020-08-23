      SUBROUTINE CALC_FACE_INFORMATION_IND(NBJ,NBJF,nef,NKJE,NKEF,
     '  NKJF,NNF,NPNE,NPNF,nr,NVJE,NVJF,SE,SF,ERROR,*)

C#### Subroutine: CALC_FACE_INFORMATION_IND
C###  Description:
C###    CALC_FACE_INFORMATION_IND calculates the independent variable
C###    face information for local face nef of global element ne.

C#### Variable: NKJF(nk,nn,nj)
C###  Type: INTEGER
C###  Set_up: CALC_FACE_INFORMATION_IND
C###  Description:
C###    NKJF(nk,nn,nj) is the global derivative number of local
C###    derivative nk of node nn for independent variable nj.  This is a
C###    temporary array set up for one face at a time.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
!     Parameter List
      INTEGER NBJ(NJM),NBJF(NJM),nef,NKJE(NKM,NNM,NJM),
     '  NKEF(0:4,16,6,NBFM),NKJF(NKM,NNM,NJM),NNF(0:17,6,NBFM),
     '  NPNE(NNM,NBFM),NPNF(NNM,NBFM),nr,NVJE(NNM,NBFM,NJM),
     '  NVJF(NNM,NBFM,NJM)
      REAL*8 SE(NSM,NBFM),SF(NSM,NBFM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb_e,nb_f,nj,njj1,njj2,nk_f,nk_e,nn2,nn_e,nn_f,ns_f,ns_e

      CALL ENTERS('CALC_FACE_INFORMATION_IND',*9999)

      DO njj1=1,3
        DO njj2=1,NJ_LOC(njj1,0,nr)
          nj=NJ_LOC(njj1,njj2,nr)
          nb_e=NBJ(nj)
          nb_f=NBJF(nj)
          ns_f=0
          DO nn_f=1,NNF(0,nef,nb_e)
            nn_e=NNF(1+nn_f,nef,nb_e)
            NPNF(nn_f,nb_f)=NPNE(nn_e,nb_e)
            ns_e=0
            DO nn2=1,nn_e-1
              ns_e=ns_e+NKT(nn2,nb_e)
            ENDDO !nn2
            DO nk_f=1,NKEF(0,nn_f,nef,nb_e)
              nk_e=NKEF(nk_f,nn_f,nef,nb_e)
              ns_f=ns_f+1
              NKJF(nk_f,nn_f,nj)=NKJE(nk_e,nn_e,nj)
              SF(ns_f,nb_f)=SE(ns_e+nk_e,nb_e)
            ENDDO !nk
            NVJF(nn_f,nb_f,nj)=NVJE(nn_e,nb_e,nj)
          ENDDO !nnf
        ENDDO !njj2
      ENDDO !njj1

      CALL EXITS('CALC_FACE_INFORMATION_IND')
      RETURN
 9999 CALL ERRORS('CALC_FACE_INFORMATION_IND',ERROR)
      CALL EXITS('CALC_FACE_INFORMATION_IND')
      RETURN 1
      END


