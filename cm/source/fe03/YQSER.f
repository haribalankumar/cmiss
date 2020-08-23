      SUBROUTINE YQSER(IBT,IDO,INP,NBH,ne,NENQ,njj,NQNE,NQS,
     '  NQXI,nx,ER,SE,YQS,ERROR,*)

C#### Subroutine: YQSER
C###  Description:
C###    YQSER evaluates element rhs, ER(ns), in calculation of least
C###    squares fit of linear field variables, defined by nodal values
C###    ZP(nk,nv,nj,np,1), to the set of data values YQS(niyq,nq)
C####   defined at grid points.

      IMPLICIT NONE
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NBH(NHM),ne,NENQ(0:8,NQM),
     '  njj,NQNE(NEQM,NQEM),
     '  NQS(NEQM),NQXI(0:NIM,NQSCM),nx
      REAL*8 ER(NHM*NSM),SE(NSM,NBFM),YQS(NIQSM,NQM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER II,IJ,IK,nb,neq,ng,nh,nhj1,
     '  nhs1,nhx,nik,nij,nii,nk1,nn1,nq,ns1,SCHEME
      REAL*8 PSI1,SUM1,XI(3)

      CALL ENTERS('YQSER',*9999)

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

            SCHEME=NQS(ne)
            II=MAX(1,NQXI(1,SCHEME))
            IJ=1
            IK=1
            IF(NQXI(0,SCHEME).GT.1) IJ=MAX(1,NQXI(2,SCHEME))
            IF(NQXI(0,SCHEME).GT.2) IK=MAX(1,NQXI(3,SCHEME))

C           Loop over the grid points in each element
            DO nik=1,IK
              DO nij=1,IJ
                DO nii=1,II
                  neq=nii+((nij-1)*NQXI(1,SCHEME)) !local grid pt #
                  IF(NQXI(0,SCHEME).GT.1) neq=neq+((nik-1)*
     '              NQXI(1,SCHEME)*NQXI(2,SCHEME))
                  nq=NQNE(ne,neq)                  !global grid pt #
C                 Evaluate each grid point only once
                  IF(NENQ(1,nq).EQ.ne) THEN
C                   Local Xi coords of grid point nq in element ne
                    IF(II.NE.1) XI(1)=DBLE(nii-1)/DBLE(II-1)
                    IF(IJ.NE.1) XI(2)=DBLE(nij-1)/DBLE(IJ-1)
                    IF(IK.NE.1) XI(3)=DBLE(nik-1)/DBLE(IK-1)

                    SUM1=SUM1+PSI1(IBT(1,1,nb),IDO(1,1,0,nb),
     '                INP(1,1,nb),nb,1,nk1,nn1,XI)*
     '                YQS(NQ_FIT(nhj1,njj),nq)

                  ENDIF !NENQ(1,nq).EQ.ne
                ENDDO !nii
              ENDDO !nij
            ENDDO !nik

            ER(nhs1)=SUM1*SE(ns1,nb)
          ENDDO !nk1
        ENDDO !nn1
      ENDDO !nhj1

      CALL EXITS('YQSER')
      RETURN
 9999 CALL ERRORS('YQSER',ERROR)
      CALL EXITS('YQSER')
      RETURN 1
      END


