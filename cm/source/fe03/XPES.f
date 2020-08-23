      SUBROUTINE XPES(IBT,IDO,INP,NBH,NBJ,njj,NKJE,NPF,NPLIST2,
     '  NPNE,NRE,NVJE,nx,ES,PG,RG,SE,WDL,WG,WU,XA,XE,XG,XIDL,XP,ERROR,*)

C#### Subroutine: XPES
C###  Description:
C###    XPES evaluates element stiffness matrix ES(ms,ns) in calculation
C###    of least squares fit of linear field variables, defined by nodal
C###    values XP(nk,nv,nj,np), to a group of selected nodes. This is
C###    a modification of ZDES

C**** NDLT         is  number of data points within current element ne
C**** ZDL(nh,npe)   "  r.c. coords & value of element node pt npe
C**** WDL(nh,npe)   "  weighting factor for     "      "      "
C**** XIDL(ni,npe)  "  Xi-coordinate of         "      "      "

      IMPLICIT NONE
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBH(NHM),NBJ(NJM),njj,NKJE(NKM,NNM,NJM),NPF(9,NFM),
     '  NPNE(NNM,NBFM),NRE,NVJE(NNM,NBFM,NJM),nx,NPLIST2(0:NPM)
      REAL*8 ES(NHM*NSM,NHM*NSM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM),WDL(NHM,NDEM),WG(NGM,NBM),WU(0:NUM+1),
     '  XA(NAM,NJM),XE(NSM,NJM),XG(NJM,NUM),XIDL(NIM,NDEM),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,npe,ng,nh1,nhj1,nhj2,nhs1,nhs2,nhx,ni,nk1,nk2,nn1,nn2,
     '  NNTB,ns1,ns2,nu,NUTB
      REAL*8 AG(6),DXIX(3,3),G11,G12,G13,G21,G22,G23,G31,G32,G33,
     '  GL(3,3),GU(3,3),PSI1,PSI2,PSI5,PSI2_HERMITE,SUM1,SUM2,SUM3
      LOGICAL SECTOR

      CALL ENTERS('XPES',*9999)

      nhs1=0
      DO nhj1=1,NUM_FIT(njj) !nhj are vars for the fit problem njj
        nhx=NLH_FIT(nhj1,3,njj)
        nh1=NH_LOC(nhx,nx)
        nb=NBH(nh1)
        SECTOR=.FALSE.
        DO ni=1,NIT(nb)
          IF(IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6) SECTOR=.TRUE.
        ENDDO !ni
C SAB 13Sep00 Shifted the evaluation of the Jacobian to the start of the
C     routine as it wasn't changing and was making the calculation very slow.
        IF(KTYP12.EQ.2) THEN !smooth field parameters
          CALL XPXE(NBJ,NKJE,NPF(1,1),NPNE,NRE,NVJE,
     '      SE,XA,XE,XP,ERROR,*9999)
          DO ng=1,NGT(nb)
            CALL XEXG(NBJ(1),ng,NRE,PG,XE,XG,ERROR,*9999)
            CALL XGMG(0,NIT(nb),nb,NRE,DXIX,GL,GU,RG(ng),
     '        XG,ERROR,*9999)
          ENDDO !ng
        ENDIF
        NNTB=NNT(nb)
        ns1=0
        DO nn1=1,NNTB
          DO nk1=1,NKT(nn1,nb)
            nhs1=nhs1+1
            ns1=ns1+1
            nhs2=0
            DO nhj2=1,NUM_FIT(njj) !columns
              ns2=0
              DO nn2=1,NNTB
                DO nk2=1,NKT(nn2,nb)
                  nhs2=nhs2+1
                  ns2=ns2+1
                  IF(nhj2.EQ.nhj1) THEN !to avoid coupling for now
                    SUM1=0.0d0
C CPB 3/4/94 Adding Simplex elements
                    IF(IBT(1,1,nb).EQ.3) THEN !Simplex
                      IF(IBT(1,2,nb).EQ.3) THEN !Hermite-Simplex
                        DO npe=1,nplist2(0)
                          SUM1=SUM1+PSI2_HERMITE(IDO(1,1,0,nb),
     '                      INP(1,1,nb),nb,1,nk1,nn1,XIDL(1,npe))*
     '                      PSI2_HERMITE(IDO(1,1,0,nb),INP(1,1,nb),
     '                      nb,1,nk2,nn2,XIDL(1,npe))*WDL(nh1,npe)
                        ENDDO !npe
                      ELSE !Normal Simplex
                        DO npe=1,nplist2(0)
                          SUM1=SUM1+PSI2(IBT(1,1,nb),INP(1,1,nb),nb,
     '                      1,nn1,XIDL(1,npe)) *
     '                      PSI2(IBT(1,1,nb),
     '                      INP(1,1,nb),nb,1,nn2,XIDL(1,npe))
     '                      *WDL(nh1,npe)
                        ENDDO !npe
                      ENDIF
                    ELSE !Lagrange/Hermite Tensor product
C cpb 2/2/97 Adding sector elements
                      IF(SECTOR) THEN
                        DO npe=1,nplist2(0)
                          SUM1=SUM1+PSI5(IBT(1,1,nb),IDO(1,1,0,nb),
     '                      INP(1,1,nb),nb,1,nk1,nn1,XIDL(1,npe)) *
     '                      PSI5(IBT(1,1,nb),IDO(1,1,0,nb),
     '                      INP(1,1,nb),nb,1,nk2,nn2,XIDL(1,npe))
     '                      *WDL(nh1,npe)
                        ENDDO !npe
                      ELSE
                        DO npe=1,nplist2(0)
                          SUM1=SUM1+PSI1(IBT(1,1,nb),IDO(1,1,0,nb),
     '                      INP(1,1,nb),nb,1,nk1,nn1,XIDL(1,npe)) *
     '                      PSI1(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                      nb,1,nk2,nn2,XIDL(1,npe))
     '                      *WDL(nh1,npe)
                        ENDDO !npe
                      ENDIF
                    ENDIF
                    SUM2=0.0d0

                    IF(KTYP12.EQ.1.OR.KTYP12.EQ.2) THEN !Sobolev smoothing term
C LKC 15-OCT-97 added assert statement
C KAT 23Nov98: modifying assert for different NIT
C                        CALL ASSERT(NUM.GE.6,'>>Increase NUM to 6',
C     '                    ERROR,*9999)
                      IF(NIT(nb).EQ.1) THEN
                        NUTB=3
                      ELSE IF(NIT(nb).EQ.2) THEN
                        NUTB=6
                      ELSE
                        NUTB=11
                      ENDIF
                      IF(NUTB.GT.NUM) THEN
                        WRITE(ERROR,'(''>>Increase NUM to '',I2)') NUTB
                        GO TO 9999
                      ENDIF
                      IF(KTYP12.EQ.1) THEN !smooth field parameters
C KAT 23Nov98: No Jacobian as this makes the problem nonlinear
                        DO ng=1,NGT(nb)
                          SUM3=0.0d0
                          DO nu=2,NUTB
                            SUM3=SUM3+
     '                        PG(ns1,nu,ng,nb)*PG(ns2,nu,ng,nb)*WU(nu)
                          ENDDO !nu
                          SUM2=SUM2+SUM3*WG(ng,nb)
                        ENDDO !ng
                      ELSE !smooth deviation from initial field
C cpb 13/5/98 Adding the Jacobian to the Sobolev integration
C SAB 13Sep00 Shifted the calculation of the Jacobian to the top
                        DO ng=1,NGT(nb)
                          SUM3=0.0d0
                          DO nu=2,NUTB
                            SUM3=SUM3+
     '                        PG(ns1,nu,ng,nb)*PG(ns2,nu,ng,nb)*WU(nu)
                          ENDDO !nu
                          SUM2=SUM2+SUM3*WG(ng,nb)*RG(ng)
                        ENDDO !ng
                      ENDIF

                    ELSE IF(KTYP12.EQ.3) THEN !StrEnergy smoothing term
                      CALL XPXE(NBJ,NKJE,NPF(1,1),NPNE,NRE,NVJE,
     '                  SE,XA,XE,XP,ERROR,*9999)

C cpb 13/5/98 There is no Jacobian in this integration either????

                      DO ng=1,NGT(nb)

                        CALL XEXG(NBJ(1),ng,NRE,PG,XE,XG,ERROR,*9999)
                        CALL XGMG(0,NIT(nb),nb,NRE,DXIX,GL,GU,RG(ng),XG,
     '                    ERROR,*9999)

C***                    Gij are compon of the matrx del(Xi-ij)/del(X-ij)
                        G11=DXIX(1,1)
                        G21=DXIX(2,1)
                        G31=DXIX(3,1)
                        G12=DXIX(1,2)
                        G22=DXIX(2,2)
                        G32=DXIX(3,2)
                        G13=DXIX(1,3)
                        G23=DXIX(2,3)
                        G33=DXIX(3,3)

                        AG(1)=G11**2+G12**2+G13**2
                        AG(2)=G21**2+G22**2+G23**2
                        AG(3)=G31**2+G32**2+G33**2
                        AG(4)=G21*G31+G22*G32+G23*G33
                        AG(5)=G11*G21+G12*G22+G13*G23
                        AG(6)=G11*G31+G12*G32+G13*G33

C***                    Note: Contrib to the stiffn. matrix due to s.e.
C**                           smoothing are the same in all 3 dirns.
                        SUM2=SUM2+WU(2)*
     '                    (AG(1)* PG(ns1,2,ng,nb)*PG(ns2,2,ng,nb)
     '                    +AG(2)* PG(ns1,4,ng,nb)*PG(ns2,4,ng,nb)
     '                    +AG(3)* PG(ns1,7,ng,nb)*PG(ns2,7,ng,nb)
     '                    +AG(4)*(PG(ns1,4,ng,nb)*PG(ns2,7,ng,nb)
     '                    +PG(ns1,7,ng,nb)*PG(ns2,4,ng,nb))
     '                    +AG(5)*(PG(ns1,2,ng,nb)*PG(ns2,4,ng,nb)
     '                    +PG(ns1,4,ng,nb)*PG(ns2,2,ng,nb))
     '                    +AG(6)*(PG(ns1,2,ng,nb)*PG(ns2,7,ng,nb)
     '                    +PG(ns1,7,ng,nb)*PG(ns2,2,ng,nb)))*WG(ng,nb)

                      ENDDO !ng
                    ENDIF

                    ES(nhs1,nhs2)=ES(nhs1,nhs2)+(SUM1+SUM2*WU(0))*
     '                SE(ns1,nb)*SE(ns2,nb)

                  ELSE
                    ES(nhs1,nhs2)=0.0d0
                  ENDIF !nhj2=nhj1
                ENDDO !nk2
              ENDDO !nn2
            ENDDO !nhj2
          ENDDO !nk1
        ENDDO !nn1
      ENDDO !nhj1

      CALL EXITS('XPES')
      RETURN
 9999 CALL ERRORS('XPES',ERROR)
      CALL EXITS('XPES')
      RETURN 1
      END



