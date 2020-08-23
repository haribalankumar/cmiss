      SUBROUTINE MELGEF_FACE(LGE_FACE,NBHF,NHST_FACE,njj,NKH,
     '  NKJF,NPNF,NPNODE,nr,NVHF,NVHP,nx,NYNP,ERROR,*)

C#### Subroutine: MELGEF_FACE
C###  Description:
C###    MELGEF_FACE  determines arrays LGE_FACE and NHST_FACE.
C###    LGE_FACE contains the positions of the components in the
C###    element stifness matrix and RHS vector in the global stifness
C###    matrix and global RHS vector.
C###    NHST_FACE contains total number of face variables.
C**** Written by Kumar Mithraratne, Aug. 2002.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER LGE_FACE(48),NBHF(NHM,NCM),NHST_FACE(NRCM),
     '  njj,NKH(NHM,NPM,NCM,0:NRM),NKJF(NKM,NNM,NJM),NPNF(NNM,NBFM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NVHF(NNM,NBFM,NHM),NVHP(NHM,NPM,NCM),
     '  nx,NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nbf,nh,nhj,nhx,nk,nkf,nkface,nn_f,nonode,np,npface,
     '  nrc,nv,nvface


      CALL ENTERS('MELGEF_FACE',*9999)

      DO nrc=1,2
        NHST_FACE(nrc)=0
        DO nhj=1,NUM_FIT(njj)
          nhx=NLH_FIT(nhj,3,njj)
          nh=NH_LOC(nhx,nx)
          nbf=NBHF(nh,1)
          DO nn_f=1,NNT(nbf)
            npface=NPNF(nn_f,nbf)
            nvface=NVHF(nn_f,nbf,nh)
            DO nkf=1,NKT(nn_f,nbf)
              nkface=NKJF(nkf,nn_f,nh)
              NHST_FACE(nrc)=NHST_FACE(nrc)+1
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO nv=1,NVHP(nh,np,1)
                  DO nk=1,NKH(nh,np,1,nr)
                    IF((npface.EQ.np).AND.(nkface.EQ.nk).AND.
     '                (nvface.EQ.nv)) THEN
                      LGE_FACE(NHST_FACE(nrc))=
     '                  NYNP(nk,nv,nh,np,nrc,1,nr)
                    ENDIF
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nonode
            ENDDO !nkf
          ENDDO !nn_f
        ENDDO !nh
      ENDDO !nrc

      CALL EXITS('MELGEF_FACE')
      RETURN
 9999 CALL ERRORS('MELGEF_FACE',ERROR)
      CALL EXITS('MELGEF_FACE')
      RETURN 1
      END


