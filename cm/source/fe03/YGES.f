      SUBROUTINE YGES(NBH,njj,nx,ES,PG,SE,WG,WU,ERROR,*)

C#### Subroutine: YGES
C###  Description:
C###    YGES evaluates element stiffness matrix ES(ms,ns) for Gauss
C###    point data in calculation of least squares fit of linear field
C###    variables, defined by nodal values ZP(nk,nv,nj,np,1), to the
C###    set of data values YG(nh0,ng,ne) defined at the Gauss points of
C###    basis nb.

      IMPLICIT NONE
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM),njj,nx
      REAL*8 ES(NHM*NSM,NHM*NSM),PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM),
     '  WG(NGM,NBM),WU(0:NUM+1)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ng,nh1,nhj1,nhj2,nhs1,nhs2,nhx,nk1,nk2,nn1,nn2,
     '  NNTB,ns1,ns2
      REAL*8 SUM1,SUM2

      CALL ENTERS('YGES',*9999)

      nhs1=0
      DO nhj1=1,NUM_FIT(njj)
        nhx=NLH_FIT(nhj1,3,njj)
        nh1=NH_LOC(nhx,nx)
        nb=NBH(nh1)
        NNTB=NNT(nb)
        ns1=0
        DO nn1=1,NNTB
          DO nk1=1,NKT(nn1,nb)
            nhs1=nhs1+1
            ns1=ns1+1
            nhs2=0
            DO nhj2=1,NUM_FIT(njj)
              nhx=NLH_FIT(nhj2,3,njj)
              ns2=0
              DO nn2=1,NNTB
                DO nk2=1,NKT(nn2,nb)
                  nhs2=nhs2+1
                  ns2=ns2+1
                  SUM1=0.0d0
                  SUM1=0.d0
                  DO ng=1,NGT(nb)
                    SUM1=SUM1+PG(ns1,1,ng,nb)*PG(ns2,1,ng,nb)
                  ENDDO !ng
                  SUM2=0.d0
                  IF(KTYP12.EQ.1) THEN
                    DO ng=1,NGT(nb)
                      SUM2=SUM2+
     '                  (PG(ns1,2,ng,nb)*PG(ns2,2,ng,nb)*WU(2)+
     '                  PG(ns1,3,ng,nb)*PG(ns2,3,ng,nb)*WU(3)+
     '                  PG(ns1,4,ng,nb)*PG(ns2,4,ng,nb)*WU(4)+
     '                  PG(ns1,5,ng,nb)*PG(ns2,5,ng,nb)*WU(5)+
     '                  PG(ns1,6,ng,nb)*PG(ns2,6,ng,nb)*WU(6)+
     '                  PG(ns1,7,ng,nb)*PG(ns2,7,ng,nb)*WU(7)+
     '                  PG(ns1,8,ng,nb)*PG(ns2,8,ng,nb)*WU(8)+
     '                  PG(ns1,9,ng,nb)*PG(ns2,9,ng,nb)*WU(9)+
     '                  PG(ns1,10,ng,nb)*PG(ns2,10,ng,nb)*WU(10))*
     '                  WG(ng,nb)
                    ENDDO !ng
                  ENDIF
                  ES(nhs1,nhs2)=SUM1*SE(ns1,nb)*SE(ns2,nb)+SUM2*WU(0)
                ENDDO !nk2
              ENDDO !nn2
            ENDDO !nhj2
          ENDDO !nk1
        ENDDO !nn1
      ENDDO !nhj1

      CALL EXITS('YGES')
      RETURN
 9999 CALL ERRORS('YGES',ERROR)
      CALL EXITS('YGES')
      RETURN 1
      END


