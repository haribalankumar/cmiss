      SUBROUTINE YGER(NBH,njj,nx,ER,PG,SE,YG,ERROR,*)

C#### Subroutine: YGER
C###  Description:
C###    YGER evaluates element rhs, ER(ns), in calculation of least
C###    squares fit of linear field variables, defined by nodal values
C###    ZP(nk,nv,nj,np,1), to the set of data values YG(niyg,ng)
C####   defined at the Gauss points of basis nb.

      IMPLICIT NONE
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM),njj,nx
      REAL*8 ER(NHM*NSM),PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM),YG(NIYGM,NGM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ng,nh,nhj1,nhs1,nhx,nk1,nn1,ns1
      REAL*8 SUM1

      CALL ENTERS('YGER',*9999)

      nhs1=0
      DO nhj1=1,NUM_FIT(njj)
        nhx=NLH_FIT(nhj1,3,njj)
        nh=NH_LOC(nhx,nx)
        nb=NBH(nh)
        ns1=0
        DO nn1=1,NNT(nb)
          DO nk1=1,NKT(nn1,nb)
            nhs1=nhs1+1
            ns1=ns1+1
            SUM1=0.d0
            DO ng=1,NGT(nb)
              SUM1=SUM1+PG(ns1,1,ng,nb)*YG(NG_FIT(nhj1,njj),ng)
            ENDDO !ng
            ER(nhs1)=SUM1*SE(ns1,nb)
          ENDDO !nk1
        ENDDO !nn1
      ENDDO !nhj1

      CALL EXITS('YGER')
      RETURN
 9999 CALL ERRORS('YGER',ERROR)
      CALL EXITS('YGER')
      RETURN 1
      END


