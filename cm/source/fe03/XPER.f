      SUBROUTINE XPER(IBT,IDO,INP,NBH,NBJ,
     '  njj,NKJE,NPF,NPLIST2,NPNE,NRE,NVJE,nx,ER,PG,RG,SE,WDL,
     '  WG,WU,XA,XE,XG,XIP,XIDL,XP,ZDL,ZE,ZP,ERROR,*)

C#### Subroutine: XPER
C###  Description:
C###    Evaluates element rhs, ER(ns), in calculation of least squares
C###    fit of linear field variables, defined by nodal values
C###    XP(nk,nv,nj,np). This is a modification of ZDER.

C**** NDLT         is  no. data points within current element ne
C**** NDDL(ne,nde)  "  global data pt no. of element data pt nde
C**** ZDL(nh,nde)   "  r.c. coords & value of   "      "      "
C**** WDL(nh,nde)   "  weighting factor for     "      "      "
C**** XIDL(ni,nde)  "  Xi-coordinate of         "      "      "

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBH(NHM),NBJ(NJM),NPLIST2(0:NPM),njj,
     '  NKJE(NKM,NNM,NJM),NPF(9,NFM),NPNE(NNM,NBFM),NRE,
     '  NVJE(NNM,NBFM,NJM),nx
      REAL*8 ER(NHM*NSM),PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM),
     '  WDL(NHM,NDEM),WG(NGM,NBM),WU(0:NUM+1),
     '  XA(NAM,NJM),XE(NSM,NJM),XG(NJM,NUM),XIP(NIM,NPM),
     '  XIDL(NIM,NDEM),XP(NKM,NVM,NJM,NPM),ZDL(NHM,NDEM),
     '  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,npe,ng,nh,nhj1,nhs1,nhx,ni,nj,nk1,nn1,ns1,ns2,nr,
     '  nu,NUTB
      REAL*8 BCQ(9),DXIX(3,3),F(3,3),G11,G12,G13,G21,G22,G23,G31,G32,
     '  G33,GL(3,3),GU(3,3),PSI1,PSI2,PSI2_HERMITE,PSI5,PXI,R(3,3),
     '  R11,R12,R13,R21,R22,R23,R31,R32,R33,SUM1,SUM2,SUM3,SUM4,U(3,3),X
      LOGICAL SECTOR

      CALL ENTERS('XPER',*9999)

      CALL XPZDL(1,NBH,NBJ,nplist2,nx,WDL,XIP,XIDL,XP,ZDL,ERROR,*9999)
      nhs1=0
      nr=NRE
      DO nhj1=1,NUM_FIT(njj)
        nhx=NLH_FIT(nhj1,3,njj)
        nh=NH_LOC(nhx,nx)
        nb=NBH(nh)
        SECTOR=.FALSE.
        DO ni=1,NIT(nb)
          IF(IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6) SECTOR=.TRUE.
        ENDDO !ni
        DO npe=1,nplist2(0)
          X=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XIDL(1,npe),
     '      ZE(1,nhx))
          ZDL(nh,npe)=ZDL(nh,npe)-X
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' npe='',I5,'' nh='',I1,'' x='',D12.4,'
     '        //''' ZDL='',D12.4)') npe,nh,X,ZDL(nh,npe)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDDO !npe
        ns1=0
        DO nn1=1,NNT(nb)
          DO nk1=1,NKT(nn1,nb)
            nhs1=nhs1+1
            ns1=ns1+1
            SUM1=0.0d0
C CPB 3/4/94 Adding Simplex elements
            IF(IBT(1,1,nb).EQ.3) THEN ! Simplex
              IF(IBT(1,2,nb).EQ.3) THEN ! Hermite-Simplex
                DO npe=1,nplist2(0)
                  SUM1=SUM1+PSI2_HERMITE(IDO(1,1,0,nb),INP(1,1,nb),nb,
     '              1,nk1,nn1,XIDL(1,npe))*ZDL(nh,npe)*WDL(nh,npe)
                ENDDO !npe
              ELSE ! Normal Simplex
                DO npe=1,nplist2(0)
                  SUM1=SUM1+PSI2(IBT(1,1,nb),INP(1,1,nb),nb,1,nn1,
     '              XIDL(1,npe))*ZDL(nh,npe)*WDL(nh,npe)
                ENDDO !npe
              ENDIF
            ELSE ! Lagrange/Hermite Tensor product
C cpb 2/2/97 Adding sector elements
              IF(SECTOR) THEN
                DO npe=1,nplist2(0)
                  SUM1=SUM1+PSI5(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '              nb,1,nk1,nn1,XIDL(1,npe))*ZDL(nh,npe)*WDL(nh,npe)
                ENDDO !npe
              ELSE
                DO npe=1,nplist2(0)
                  SUM1=SUM1+PSI1(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '              nb,1,nk1,nn1,XIDL(1,npe))*ZDL(nh,npe)*WDL(nh,npe)
                ENDDO !npe
              ENDIF
            ENDIF

            SUM2=0.0d0
            IF(KTYP12.EQ.1) THEN !Sobolev smoothing on field paramters
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
              DO ng=1,NGT(nb)
                SUM3=0.0d0
                DO nu=2,NUTB
                  SUM4=0.0d0
                  DO ns2=1,NST(nb)
                    SUM4=SUM4+ZE(ns2,nhx)*PG(ns2,nu,ng,nb)
                  ENDDO !ns2
                  SUM3=SUM3+SUM4*PG(ns1,nu,ng,nb)*WU(nu)
                ENDDO !nu
                SUM2=SUM2-SUM3*WG(ng,nb) !*RG(ng)
              ENDDO !ng

            ELSE IF(KTYP12.EQ.3) THEN !Strain energy smoothing term

              DO ng=1,NGT(nb)

C***            Compute dXI-i/dX-J based on undef configuration (XP)
                CALL XPXE(NBJ,NKJE,NPF(1,1),NPNE,NR,NVJE,
     '            SE,XA,XE,XP,ERROR,*9999)
                CALL XEXG(NBJ,ng,NR,PG,XE,XG,ERROR,*9999)
                CALL XGMG(0,NIT(nb),nb,NR,DXIX,GL,GU,RG(ng),XG,
     '            ERROR,*9999)
C***            Note: DXIX should be referred to reference coords.

C***            Compute dx-i/dXI-k based on deformed configuration (ZP)
                CALL XPXE(NBJ,NKJE,NPF(1,1),NPNE,nr,NVJE,
     '            SE,XA,XE,ZP,ERROR,*9999)
                CALL XEXG(NBJ,ng,NR,PG,XE,XG,ERROR,*9999)

                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  DO ni=1,NIT(nb)
                    F(nj,ni)=XG(nj,2)*DXIX(1,ni)+XG(nj,4)*DXIX(2,ni)
                    IF(NIT(nb).EQ.3) THEN
                      F(nj,ni)=F(nj,ni)+XG(nj,7)*DXIX(3,ni)
                    ENDIF
                  ENDDO !ni
                ENDDO !nj
                IF (NIT(nb).EQ.2) THEN
                  F(3,3)=1.0d0
                  F(3,2)=0.0d0
                  F(3,1)=0.0d0
                  F(1,3)=0.0d0
                  F(2,3)=0.0d0
                ENDIF

                CALL POLAR(2,F,R,U,ERROR,*9999) !RU decomposition

                IF (NIT(nb).EQ.2)THEN
C                 POLAR didn't affect 3rd row or column
                  R(3,1)=0.0d0
                  R(3,2)=0.0d0
                  R(1,3)=0.0d0
                  R(2,3)=0.0d0
                  R(3,3)=1.0d0
                ENDIF

C***            Rij are components of the rotation matrix computed at
C***            the last iteration.

                R11=R(1,1)
                R21=R(2,1)
                R31=R(3,1)
                R12=R(1,2)
                R22=R(2,2)
                R32=R(3,2)
                R13=R(1,3)
                R23=R(2,3)
                R33=R(3,3)

C***            Gij are components of the matrix del(Xi-ij)/del(X-ij)

                G11=DXIX(1,1)
                G21=DXIX(2,1)
                G31=DXIX(3,1)
                G12=DXIX(1,2)
                G22=DXIX(2,2)
                G32=DXIX(3,2)
                G13=DXIX(1,3)
                G23=DXIX(2,3)
                G33=DXIX(3,3)

C***            BCQ(i) are the coefficients on the right hand side terms

                BCQ(1)=G11*(1.0d0-R11)-G12*R12-G13*R13
                BCQ(2)=-G11*R21+G12*(1.0d0-R22)-G13*R23
                BCQ(3)=-G11*R31-G12*R32+G13*(1.0d0-R33)
                BCQ(4)=G21*(1.0d0-R11)-G22*R12-G23*R13
                BCQ(5)=-G21*R21+G22*(1.0d0-R22)-G23*R23
                BCQ(6)=-G21*R31-G22*R32+G23*(1.0d0-R33)
                BCQ(7)=G31*(1.0d0-R11)-G32*R12-G33*R13
                BCQ(8)=-G31*R21+G32*(1.0d0-R22)-G33*R23
                BCQ(9)=-G31*R31-G32*R32+G33*(1.0d0-R33)

                IF(NLH_FIT(nhj1,2,njj).EQ.1) THEN !x-direction term
                  SUM2=SUM2-WU(2)*
     '              (BCQ(1)*PG(ns1,2,ng,nb)
     '              +BCQ(4)*PG(ns1,4,ng,nb)
     '              +BCQ(7)*PG(ns1,7,ng,nb))*WG(ng,nb)
                ELSE IF(NLH_FIT(nhj1,2,njj).EQ.2) THEN !y-direction term
                  SUM2=SUM2-WU(2)*
     '              (BCQ(2)*PG(ns1,2,ng,nb)
     '              +BCQ(5)*PG(ns1,4,ng,nb)
     '              +BCQ(8)*PG(ns1,7,ng,nb))*WG(ng,nb)
                ELSE IF(NLH_FIT(nhj1,2,njj).EQ.3) THEN !z-direction term
                  SUM2=SUM2-WU(2)*
     '              (BCQ(3)*PG(ns1,2,ng,nb)
     '              +BCQ(6)*PG(ns1,4,ng,nb)
     '              +BCQ(9)*PG(ns1,7,ng,nb))*WG(ng,nb)
                ENDIF

              ENDDO !ng
            ENDIF
            ER(nhs1)=ER(nhs1)+(SUM1+SUM2*WU(0))*SE(ns1,nb)
          ENDDO !nk1
        ENDDO !nn1
      ENDDO !nhj1

      CALL EXITS('XPER')
      RETURN
 9999 CALL ERRORS('XPER',ERROR)
      CALL EXITS('XPER')
      RETURN 1
      END


