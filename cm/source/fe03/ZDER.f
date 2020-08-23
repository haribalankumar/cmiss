      SUBROUTINE ZDER(IBT,IDO,INP,NBH,NBJ,NDDL,NDLT,ne,njj,NKJE,
     '  NPF,NPNE,NRE,NVJE,nx,ER,PG,RG,SE,SF,WD,WDL,WG,WU,XA,XE,
     '  XG,XID,XIDL,XP,ZD,ZDL,ZE,ZP,IN_PLANE,ERROR,*)

C#### Subroutine: ZDER
C###  Description:
C###    Evaluates element rhs, ER(ns), in calculation of least squares
C###    fit of linear field variables, defined by nodal values
C###    XP(nk,nv,nj,np), to the set of data values XD(nj,nd) with
C###    weights WD(nj,nd) at local coordinate values XID(ni,nd).

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
     '  NBH(NHM),NBJ(NJM),NDDL(NEM,NDEM),NDLT,ne,njj,
     '  NKJE(NKM,NNM,NJM),NPF(9,NFM),NPNE(NNM,NBFM),NRE,
     '  NVJE(NNM,NBFM,NJM),nx
      REAL*8 ER(NHM*NSM),PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM),
     '  SF(NSM,NBFM),WD(NJM,NDM),WDL(NHM,NDEM),WG(NGM,NBM),WU(0:NUM+1),
     '  XA(NAM,NJM),XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),
     '  XIDL(NIM,NDEM),XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM),ZDL(NHM,NDEM),
     '  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      LOGICAL IN_PLANE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,nb,nde,ng,nh,nh2,nhj1,nhj2,nhs1,nhx,nhx2,ni,nj,nk1,
     '  nn1,ns1,ns2,nr,nu,NUTB
      REAL*8 BCQ(9),dx,DXIX(3,3),dy,dz,F(3,3),G11,G12,G13,G21,G22,
     '  G23,G31,G32,G33,GL(3,3),GU(3,3),JAC(3,3),JAC1,JAC2,JAC3,
     '  JACDET,JACORIG,PSI,PSI1,PSI2,PSI2_HERMITE,PSI5,PXI,R(3,3),
     '  R11,R12,R13,R21,R22,R23,R31,R32,R33,SUM1,SUM2,SUM3,SUM4,
     '  SUMV,U(3,3),VAR,VOL_CON,X,X3,X0,Y0,Z0,ZDL2(NHM,NDEM)
      LOGICAL SECTOR

      CALL ENTERS('ZDER',*9999)

C PM 14Aug02 : Added for face fitting
      IF (KTYP11.EQ.1) THEN        ! element fitting
        CALL ZDZDL(1,NBH,NBJ,NDDL,NDLT,ne,nx,WD,WDL,XID,XIDL,ZD,ZDL,
     '    ERROR,*9999)
      ELSEIF (KTYP11.EQ.2) THEN    ! face fitting
        CALL ZDZDFACE(NDDL,NDLT,ne,NPF,nx,WD,WDL,XID,XIDL,ZD,ZDL,
     '    ERROR,*9999)
      ENDIF

      nhs1=0
      nr=NRE

C EWR 1Apr03 : Added the functionality of only calculating error
C in a plane given by a plane-normal-vector. Only relevant to 3D fitting.
      IF(IN_PLANE.AND.NUM_FIT(njj).EQ.3) THEN
        IN_PLANE=.TRUE.
      ELSE
        IN_PLANE=.FALSE.
      ENDIF

      IF(IN_PLANE) THEN
        DO nhj1=1,NUM_FIT(njj)
          nhx=NLH_FIT(nhj1,3,njj)
          nh=NH_LOC(nhx,nx)
          DO nde=1,NDLT
            ZDL2(nh,nde)=ZDL(nh,nde)
          ENDDO
        ENDDO
      ENDIF

      DO nhj1=1,NUM_FIT(njj)
        nhx=NLH_FIT(nhj1,3,njj)
        nh=NH_LOC(nhx,nx)
        nb=NBH(nh)
        SECTOR=.FALSE.
        DO ni=1,NIT(nb)
          IF(IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6) SECTOR=.TRUE.
        ENDDO !ni
        DO nde=1,NDLT
          X=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XIDL(1,nde),
     '      ZE(1,nhx))
          X3=0
          IF(IN_PLANE) THEN
            nhx2=NLH_FIT(1,3,njj)
            nh2=NH_LOC(nhx2,nx)
            X0=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '          XIDL(1,nde),ZE(1,nhx2))
            dx=(ZDL2(nh2,nde)-X0)*WDL(nh2,nde)
            nhx2=NLH_FIT(2,3,njj)
            nh2=NH_LOC(nhx2,nx)
            Y0=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '          XIDL(1,nde),ZE(1,nhx2))
            dy=(ZDL2(nh2,nde)-Y0)*WDL(nh2,nde)
            nhx2=NLH_FIT(3,3,njj)
            nh2=NH_LOC(nhx2,nx)
            Z0=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '          XIDL(1,nde),ZE(1,nhx2))
            dz=(ZDL2(nh2,nde)-Z0)*WDL(nh2,nde)
            X3=(dx+dy+dz)*WDL(nh,nde)
          ENDIF
          ZDL(nh,nde)=ZDL(nh,nde)-X-X3

          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' nde='',I5,'' nh='',I1,'' x='',D12.4,'
     '        //''' ZDL='',D12.4)') nde,nh,X,ZDL(nh,nde)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDDO !nde

        ns1=0
        DO nn1=1,NNT(nb)
          DO nk1=1,NKT(nn1,nb)
            nhs1=nhs1+1
            ns1=ns1+1
            SUM1=0.0d0
            IF(TAG_FIT.EQV..TRUE.) THEN
              DO nhj2=1,NUM_FIT(njj)
                nh2=NH_LOC(NLH_FIT(nhj2,3,njj),nx)
                DO nde=1,NDLT
                  SUM1=SUM1+PSI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '              nb,nn1,nk1,1,XIDL(1,nde))*ZDL(nh,nde)*WDL(nh2,nde)
     '              *WDL(nh,nde)
                ENDDO
              ENDDO
            ELSE IF(IN_PLANE) THEN
              IF(IBT(1,1,nb).EQ.3) THEN ! Simplex
                IF(IBT(1,2,nb).EQ.3) THEN ! Hermite-Simplex
                  DO nde=1,NDLT
                    SUM1=SUM1+PSI2_HERMITE(IDO(1,1,0,nb),INP(1,1,nb),nb,
     '                1,nk1,nn1,XIDL(1,nde))*ZDL(nh,nde)
                  ENDDO !nde
                ELSE ! Normal Simplex
                  DO nde=1,NDLT
                    SUM1=SUM1+PSI2(IBT(1,1,nb),INP(1,1,nb),nb,1,nn1,
     '                XIDL(1,nde))*ZDL(nh,nde)
                  ENDDO !nde
                ENDIF
              ELSE ! Lagrange/Hermite Tensor product
C cpb 2/2/97 Adding sector elements
                IF(SECTOR) THEN
                  DO nde=1,NDLT
                    SUM1=SUM1+PSI5(IBT(1,1,nb),IDO(1,1,0,nb),
     '                INP(1,1,nb),nb,1,nk1,nn1,XIDL(1,nde))*ZDL(nh,nde)
                  ENDDO !nde
                ELSE
                  DO nde=1,NDLT
                    SUM1=SUM1+PSI1(IBT(1,1,nb),IDO(1,1,0,nb),
     '                INP(1,1,nb),nb,1,nk1,nn1,XIDL(1,nde))*ZDL(nh,nde)
                  ENDDO !nde
                ENDIF
              ENDIF
            ELSE ! not IN_PLANE
C CPB 3/4/94 Adding Simplex elements
              IF(IBT(1,1,nb).EQ.3) THEN ! Simplex
                IF(IBT(1,2,nb).EQ.3) THEN ! Hermite-Simplex
                  DO nde=1,NDLT
                    SUM1=SUM1+PSI2_HERMITE(IDO(1,1,0,nb),INP(1,1,nb),nb,
     '                1,nk1,nn1,XIDL(1,nde))*ZDL(nh,nde)*WDL(nh,nde)
                  ENDDO !nde
                ELSE ! Normal Simplex
                  DO nde=1,NDLT
                    SUM1=SUM1+PSI2(IBT(1,1,nb),INP(1,1,nb),nb,1,nn1,
     '                XIDL(1,nde))*ZDL(nh,nde)*WDL(nh,nde)
                  ENDDO !nde
                ENDIF
              ELSE ! Lagrange/Hermite Tensor product
C cpb 2/2/97 Adding sector elements
                IF(SECTOR) THEN
                  DO nde=1,NDLT
                    SUM1=SUM1+PSI5(IBT(1,1,nb),IDO(1,1,0,nb),
     '                INP(1,1,nb),nb,1,nk1,nn1,XIDL(1,nde))*ZDL(nh,nde)*
     '                WDL(nh,nde)
                  ENDDO !nde
                ELSE
                  DO nde=1,NDLT
                    SUM1=SUM1+PSI1(IBT(1,1,nb),IDO(1,1,0,nb),
     '                INP(1,1,nb),nb,1,nk1,nn1,XIDL(1,nde))*ZDL(nh,nde)*
     '                WDL(nh,nde)
                  ENDDO !nde
                ENDIF
              ENDIF
            ENDIF ! IN_PLANE
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

C JWF 10/8/04 Adding host-mesh volume constraint residual term
            SUMV=0.0d0 ! initialise
            IF(KTYP1V.EQ.1) THEN ! Add volume constraint

              CALL XPXE(NBJ,NKJE,NPF(1,1),NPNE,NR,NVJE,
     &          SE,XA,XE,XP,ERROR,*9999)
              DO ng=1,NGT(nb)
  
                CALL XEXG(NBJ,ng,NR,PG,XE,XG,ERROR,*9999)
                CALL XGMG(0,NIT(nb),nb,NR,DXIX,GL,GU,RG(ng),XG,
     &            ERROR,*9999)

                DO i=1,3
                  nu=1+i*(1+i)/2
                  DO j=1,3
                    JAC(i,j)=XG(j,nu)
                  ENDDO
                ENDDO

                JAC1=(JAC(3,3)*JAC(2,2))-(JAC(3,2)*JAC(2,3))
                JAC2=(JAC(3,3)*JAC(2,1))-(JAC(3,1)*JAC(2,3))
                JAC3=(JAC(3,2)*JAC(2,1))-(JAC(3,1)*JAC(2,2))
                JACDET=dabs((JAC(1,1)*JAC1)-(JAC(1,2)*JAC2)+
     &            (JAC(1,3)*JAC3))

                DO i=1,3
                  nu=1+i*(1+i)/2
                  DO j=1,3
                    JAC(i,j)=XG(j+6,nu)
                  ENDDO
                ENDDO

                JAC1=(JAC(3,3)*JAC(2,2))-(JAC(3,2)*JAC(2,3))
                JAC2=(JAC(3,3)*JAC(2,1))-(JAC(3,1)*JAC(2,3))
                JAC3=(JAC(3,2)*JAC(2,1))-(JAC(3,1)*JAC(2,2))
                JACORIG=dabs((JAC(1,1)*JAC1)-(JAC(1,2)*JAC2)+
     &            (JAC(1,3)*JAC3))

                VOL_CON=VOL_WT*(JACDET-JACORIG)                                

C Jacobian has 6 terms, each term consisting of the product of three 1st derivatives.
C Each term gives 3 terms in the variational statement. ie total of 18 terms

                VAR=(PG(ns1,2,ng,nb)*XG(2,4)*XG(3,7))
     &            +(PG(ns1,4,ng,nb)*XG(1,2)*XG(3,7))
     &            +(PG(ns1,7,ng,nb)*XG(1,2)*XG(2,4))
     &           -(PG(ns1,2,ng,nb)*XG(2,7)*XG(3,4))
     &           -(PG(ns1,4,ng,nb)*XG(1,2)*XG(2,7))
     &           -(PG(ns1,7,ng,nb)*XG(1,2)*XG(3,2))
     &           -(PG(ns1,2,ng,nb)*XG(3,7)*XG(1,4))
     &           -(PG(ns1,4,ng,nb)*XG(2,2)*XG(3,7))
     &           -(PG(ns1,7,ng,nb)*XG(2,2)*XG(1,4))
     &            +(PG(ns1,2,ng,nb)*XG(1,7)*XG(3,4))
     &            +(PG(ns1,4,ng,nb)*XG(2,2)*XG(1,7))
     &            +(PG(ns1,7,ng,nb)*XG(2,2)*XG(3,4))
     &            +(PG(ns1,2,ng,nb)*XG(2,7)*XG(1,4))
     &            +(PG(ns1,4,ng,nb)*XG(3,2)*XG(2,7))
     &            +(PG(ns1,7,ng,nb)*XG(3,2)*XG(1,4))
     &           -(PG(ns1,2,ng,nb)*XG(1,7)*XG(2,4))
     &           -(PG(ns1,4,ng,nb)*XG(3,2)*XG(1,7))
     &           -(PG(ns1,7,ng,nb)*XG(3,2)*XG(2,4))

                SUMV=SUMV-(VAR*VOL_CON*WG(ng,nb))

              ENDDO !ng

            ENDIF !KTYP1V

            IF(KTYP11.EQ.1) THEN !elements
              ER(nhs1)=ER(nhs1)+(SUM1+SUM2*WU(0)+SUMV)*SE(ns1,nb)
            ELSE !faces
              ER(nhs1)=ER(nhs1)+(SUM1+SUM2*WU(0))*SF(ns1,nb)
            ENDIF
          ENDDO !nk1
        ENDDO !nn1
      ENDDO !nhj1

      CALL EXITS('ZDER')
      RETURN
 9999 CALL ERRORS('ZDER',ERROR)
      CALL EXITS('ZDER')
      RETURN 1
      END


