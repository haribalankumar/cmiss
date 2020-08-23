      SUBROUTINE XPXF(IBT,IDO,NBJ,NBJF,NEA,NEAXT,NF_EA,NIEF,
     '  NKB,NKJE,NNF,NPNE,nr,NSB,NVJE,SE,XDF,XP,ERROR,*)

C#### Subroutine: XPXF
C###  Description:
C###    XPXF transfers global parameters XP to face parameters
C###    XDF for interpolation of geometry within the face.
C**** Needs to be modified for any non-standard interpolations

C#### Variable: XDF(ns_f,jdoxf,neax,nj)
C###  Type: REAL*8
C###  Set_up: XPXF
C###  Description:
C###    XDF(ns_f,jdoxf,neax,nj) are the parameters ns_f for
C###    interpolation over a face of out-of-face xi derivative jdoxf of
C###    geometry/fibre/field variable nj using adjacent element neax.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  NBJ(NJM,NEM),NBJF(NJM),NEA(2),NEAXT,NF_EA(2),
     '  NIEF(0:2),NKB(2,2,2,NNM,NBFM),NKJE(NKM,NNM,NJM,NEM),
     '  NNF(0:17,6,NBFM),NPNE(NNM,NBFM,NEM),
     '  nr,NSB(NKM,NNM,NBFM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 SE(NSM,NBFM,NEM),XDF(NSFM,2,2,NJM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IDOI(3),jdoxf,nb_e,nb_f,ne,neax,nf_e,ni,nj,njj,njl,
     '  nk,nk_e,nk_f,nn_f,nn_e,np,ns_e,ns_f,nv

      CALL ENTERS('XPXF',*9999)

C     With Hermite interpolation in out-of-face direction, facial nodal
C     values of out-of-face first derivs are found from corresponding
C     element nodal values.

      DO neax=1,NEAXT !loop over adjacent elements
        ne=NEA(neax)
        nf_e=NF_EA(neax)

C***    Geometric variable and derivative.
C!!!    Should probably check that appropriate components are continuous.

        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          nb_e=NBJ(nj,ne)
          nb_f=NBJF(nj)
          CALL ASSERT(NST(nb_f).LE.NSFM,'>>NSFM too small',ERROR,*9999)
          ni=NIEF(0)
          IF(IBT(1,ni,nb_e).EQ.2.OR
     '      .(IBT(1,ni,nb_e).EQ.5.AND.IBT(2,ni,nb_e).EQ.4).OR
     '      .(IBT(1,ni,nb_e).EQ.6.AND.IBT(2,ni,nb_e).EQ.4)) THEN
C           Hermite in out-of-face direction
            ns_f=0
            DO nn_f=1,NNT(nb_f)
              nn_e=NNF(1+nn_f,nf_e,nb_e)
              np=NPNE(nn_e,nb_e,ne)
              nv=NVJE(nn_e,nb_e,nj,ne)
              DO nk_f=1,NKT(nn_f,nb_f)
                ns_f=ns_f+1
                IDOI(NIEF(1))=IDO(nk_f,nn_f,1,nb_f)
                IDOI(NIEF(2))=IDO(nk_f,nn_f,2,nb_f)
                DO jdoxf=1,2
                  IDOI(NIEF(0))=jdoxf !partial deriv in cross dirn
                  nk_e=NKB(IDOI(1),IDOI(2),IDOI(3),nn_e,nb_e)
                  IF(nk_e.NE.0) THEN
                    nk=NKJE(nk_e,nn_e,nj,ne)
                    ns_e=NSB(nk_e,nn_e,nb_e)
                    XDF(ns_f,jdoxf,neax,nj)=
     '                XP(nk,nv,nj,np)*SE(ns_e,nb_e,ne)
                  ELSE
                    XDF(ns_f,jdoxf,neax,nj)=0.0d0
                  ENDIF !nk_e
                ENDDO
              ENDDO !nk_f
            ENDDO !nn_f
          ELSE !Not Hermite
            ERROR='>>Only Hermite face integrals are implemented'
            GO TO 9999
          ENDIF !Hermite
        ENDDO !njj

C***    Fibre directions and field variables.

        DO njl=NJL_FIBR,NJL_FIEL,NJL_FIEL-NJL_FIBR
          DO njj=1,NJ_LOC(njl,0,nr)
            nj=NJ_LOC(njl,njj,nr)
            nb_e=NBJ(nj,ne)
            nb_f=NBJF(nj)
            CALL ASSERT(NST(nb_f).LE.NSFM,'>>NSFM too small',
     '        ERROR,*9999)
            ns_f=0
            DO nn_f=1,NNT(nb_f)
              nn_e=NNF(1+nn_f,nf_e,nb_e)
              np=NPNE(nn_e,nb_e,ne)
              nv=NVJE(nn_e,nb_e,nj,ne)
              DO nk_f=1,NKT(nn_f,nb_f)
                ns_f=ns_f+1
                IDOI(NIEF(1))=IDO(nk_f,nn_f,1,nb_f)
                IDOI(NIEF(2))=IDO(nk_f,nn_f,2,nb_f)
                IDOI(NIEF(0))=1 !no calculation of deriv in cross dirn
                nk_e=NKB(IDOI(1),IDOI(2),IDOI(3),nn_e,nb_e)
C               nk_e=NKEF(nk_f,nn_f,nf_e,nb_e) !this does not extend to lines
                nk=NKJE(nk_e,nn_e,nj,ne)
                ns_e=NSB(nk_e,nn_e,nb_e)
                XDF(ns_f,1,neax,nj)=XP(nk,nv,nj,np)*SE(ns_e,nb_e,ne)
              ENDDO !nk_f
            ENDDO !nn_f
          ENDDO !njj
        ENDDO !njl

      ENDDO !neax

      CALL EXITS('XPXF')
      RETURN
 9999 CALL ERRORS('XPXF',ERROR)
      CALL EXITS('XPXF')
      RETURN 1
      END


