      SUBROUTINE CALC_FACE_BASIS_DEP(IBT,NBH,NBHF,NEELEM,NFF,NHE,NNF,
     '  nr,nx,ERROR,*)

C#### Subroutine: CALC_FACE_BASIS_DEP
C###  Description:
C###    CALC_FACE_BASIS_DEP calculates the dependent variable
C###    face bases numbers.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NFF(6,NEM),NHE(NEM),NNF(0:17,6,NBFM),
     '  nr,nx
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nbf,nc,ne,nef,nf,nh,nhx,noelem

      CALL ENTERS('CALC_FACE_BASIS_DEP',*9999)

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        DO nc=1,NCT(nr,nx)
          DO nhx=1,NHE(ne)
            nh=NH_LOC(nhx,nx)
            nb=NBH(nh,nc,ne)
            IF(nb.NE.0) THEN
C??? CS 21/1/98 Should be all NIT(nb).EQ.2 ?
C              IF(NIT(nb).EQ.2.AND.(NBC(nb).EQ.5.OR.NBC(nb).EQ.6)) THEN
CC             Boundary element basis
              IF(NIT(nb).EQ.2) THEN
                nf=NFF(1,ne)
                IF(nf.NE.0) NBHF(nh,nc,nf)=nb
              ELSE
                DO nef=1,NFE(nb)
                  nf=NFF(nef,ne)
                  IF(nf.NE.0) THEN
                    CALL FIND_FACE_BASIS(IBT,nb,nbf,nef,NNF,ERROR,*9999)
                    IF(nbf.NE.0) THEN
                      NBHF(nh,nc,nf)=nbf
                    ELSE
                      ERROR='>>Could not find face basis'
                      GOTO 9999
                    ENDIF
                  ENDIF
                ENDDO !nef
              ENDIF
            ELSE
              WRITE(ERROR,'('' >>No basis specified for nh='',I2,'
     '          //''', nc='',I1)') nh,nc
              GOTO 9999
            ENDIF
          ENDDO !nhx
        ENDDO !nc
      ENDDO !noelem

      CALL EXITS('CALC_FACE_BASIS_DEP')
      RETURN
 9999 CALL ERRORS('CALC_FACE_BASIS_DEP',ERROR)
      CALL EXITS('CALC_FACE_BASIS_DEP')
      RETURN 1
      END


