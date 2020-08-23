      SUBROUTINE CPCGF(nb,NEA,NF_EA,NNF,NPNE,nr,nx,CE,CG,CP,PG,ERROR,*)

C#### Subroutine: CPCGF
C###  Description:
C###    CPCGF transfers element parameters CE(il,ne) or interpolates
C###    nodal parameters CP(il,np) with linear basis functions
C###    nb=NMB(il,2,nx), for il=1,ILT(1,nr,nx), in a face to return
C###    Gauss pt values CG(il,ng). ILP(il,1,nr,nx) indicates whether
C###    parameter IL is constant(1),piecewise constant(2) or piecewise
C###    linear(3).

      IMPLICIT NONE
      INCLUDE 'b13.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER nb,nb_e,nb_f,NEA(2),NF_EA(2),NNF(0:17,6,NBFM),
     '  NPNE(NNM,NBFM,NEM),nr,nx
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),
     '  CP(NMM,NPM),PG(NSM,NUM,NGM,NBM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER il,ne,nf_e,ng,nn_e,nn_f,np
      REAL*8 CPE(64),SUM

      CALL ENTERS('CPCGF',*9999)

      ne=NEA(1) !using first face
      nf_e=NF_EA(1) !using first face

      DO il=1,ILT(1,nr,nx)
        IF(ILP(il,1,nr,nx).EQ.1.OR.ILP(il,1,nr,nx).EQ.2) THEN
C         Constant spatially or may work for element based
          DO ng=1,NGT(nb)
            CG(il,ng)=CE(il,ne)
          ENDDO

        ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN
C         Piecewise linear (nodal parameters)
          nb_e=NMB(il,1,nx)
          nb_f=NMB(il,2,nx)
          DO nn_f=1,NNT(nb_f)
            nn_e=NNF(1+nn_f,nf_e,nb_e)
            np=NPNE(nn_e,nb_e,ne)
            CPE(nn_f)=CP(il,np)
          ENDDO !nn
          DO ng=1,NGT(nb_f)
            SUM=0.0D0
            DO nn_f=1,NNT(nb_f)
              SUM=SUM+PG(nn_f,1,ng,nb_f)*CPE(nn_f)
            ENDDO
            CG(il,ng)=SUM
          ENDDO !ng

        ELSE
          ERROR='Material parameters must be constant spatially '//
     '      'or defined by nodes'
          GO TO 9999
        ENDIF
      ENDDO !il

      CALL EXITS('CPCGF')
      RETURN
 9999 CALL ERRORS('CPCGF',ERROR)
      CALL EXITS('CPCGF')
      RETURN 1
      END


