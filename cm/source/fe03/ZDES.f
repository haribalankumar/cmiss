      SUBROUTINE ZDES(IBT,IDO,INP,NBH,NBJ,NDLT,njj,NKJE,NPF,NPNE,NRE,
     '  NVJE,nx,ES,PG,RG,SE,SF,WDL,WG,WU,XA,XE,XG,XIDL,XP,IN_PLANE,
     '  ERROR,*)

C#### Subroutine: ZDES
C###  Description:
C###    ZDES evaluates element stiffness matrix ES(ms,ns) in calculation
C###    of least squares fit of linear field variables, defined by nodal
C###    values XP(nk,nv,nj,np), to the set of data values XD(nj,nd) with
C###    weights WD(nj,nd) at local coordinate values XID(ni,nd), where
C###    nj=NJO.

C**** NDLT         is  number of data points within current element ne
C**** NDDL(ne,nde)  "  global data pt no. of element data pt nde
C**** ZDL(nh,nde)   "  r.c. coords & value of   "      "      "
C**** WDL(nh,nde)   "  weighting factor for     "      "      "
C**** XIDL(ni,nde)  "  Xi-coordinate of         "      "      "

      IMPLICIT NONE
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBH(NHM),NBJ(NJM),NDLT,njj,NKJE(NKM,NNM,NJM),NPF(9,NFM),
     '  NPNE(NNM,NBFM),NRE,NVJE(NNM,NBFM,NJM),nx
      REAL*8 ES(NHM*NSM,NHM*NSM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM),SF(NSM,NBFM),WDL(NHM,NDEM),WG(NGM,NBM),WU(0:NUM+1),
     '  XA(NAM,NJM),XE(NSM,NJM),XG(NJM,NUM),XIDL(NIM,NDEM),
     '  XP(NKM,NVM,NJM,NPM)
      LOGICAL IN_PLANE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nde,ng,nh1,nh2,nhj1,nhj2,nhs1,nhs1_for_nhj1,nhs2,
     '  nhx,nhx2,nk1,nk2,nn1,nn2,NNTB,ns1,ns2,nu,NUTB
      REAL*8 AG(6),DXIX(3,3),G11,G12,G13,G21,G22,G23,G31,G32,G33,
     '  GL(3,3),GU(3,3),INC,PD(NSM),PSI,SUM2,SUM3,SUMV,
     '  TEMP_PG(NGM,NUM,NSM),VAR,XGT(NJM,NUM,NGM)

      CALL ENTERS('ZDES',*9999)

C SAB 25Jan01 Rearranging this routine to speed it up.
C     Separating the data point contributions to the element stiffness
C     matrix from the smoothing terms.  This enables the expensive
C     evaluation of the basis functions at the gauss points to be
C     precalculated and just combined together in the correct combinations.

C EWR 1Apr03 : Added the functionality of only calculating error
C in a plane given by a plane-normal-vector. Only relevant to 3D fitting.
      IF(IN_PLANE.AND.NUM_FIT(njj).EQ.3) THEN
        IN_PLANE=.TRUE.
      ELSE
        IN_PLANE=.FALSE.
      ENDIF

      nhs1=0
      DO nhj1=1,NUM_FIT(njj) !nhj are vars for the fit problem njj
        nhx=NLH_FIT(nhj1,3,njj)
        nh1=NH_LOC(nhx,nx)
        nb=NBH(nh1)
C        SECTOR=.FALSE.
C        DO ni=1,NIT(nb)
C          IF(IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6) SECTOR=.TRUE.
C        ENDDO !ni
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

C JWF 10/08/04 Pre-determination of XG values for host-mesh fitting 
        IF(KTYP1V.EQ.1) THEN 
          CALL XPXE(NBJ,NKJE,NPF(1,1),NPNE,NRE,NVJE,
     &      SE,XA,XE,XP,ERROR,*9999)
          DO ng=1,NGT(nb)
            CALL XEXG(NBJ(1),ng,NRE,PG,XE,XG,ERROR,*9999)
            CALL DCOPY(NJM*NUM,XG(1,1),1,XGT(1,1,ng),1)
            CALL XGMG(0,NIT(nb),nb,NRE,DXIX,GL,GU,RG(ng),
     &        XG,ERROR,*9999)         
          ENDDO !ng
        ENDIF !KTYP1V

C SAB 21Jan01 Zero the element stiffness matrix
        nhs1_for_nhj1 = nhs1
        DO ns1=1,NST(nb)
          nhs1=nhs1+1
          nhs2=0
          DO nhj2=1,NUM_FIT(njj) !columns
            DO ns2=1,NST(nb)
              nhs2=nhs2+1
              ES(nhs1,nhs2)=0.0d0
            ENDDO !ns2
          ENDDO !nhj2

C SAB     Fill in a rearranged PG matrix as this made my problems run much faster.
          IF(KTYP12.EQ.2) THEN !Sobolev smoothing term
C LKC 15-OCT-97 added assert statement
C KAT 23Nov98: modifying assert for different NIT
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
            DO nu=1,NUTB
              DO ng=1,NGT(nb)
                TEMP_PG(ng,nu,ns1) = PG(ns1,nu,ng,nb)
              ENDDO
            ENDDO
          ENDIF !KTYP12

        ENDDO !ns1

C SAB 21Jan01 Loop over each data point and add in its contribution
        NNTB=NNT(nb)
        DO nde=1,NDLT
C SAB 21Jan01 Evaluate the basis functions at gauss points
          nhs1 = nhs1_for_nhj1
          ns1=0
          DO nn1=1,NNTB
            DO nk1=1,NKT(nn1,nb)
              nhs1=nhs1+1
              ns1=ns1+1
              PD(ns1)=PSI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '          nb,nn1,nk1,1,XIDL(1,nde))
C CPB 3/4/94 Adding Simplex elements
C              IF(IBT(1,1,nb).EQ.3) THEN !Simplex
C                IF(IBT(1,2,nb).EQ.3) THEN !Hermite-Simplex
C                  PD(ns1) = PSI2_HERMITE(IDO(1,1,0,nb),
C     '              INP(1,1,nb),nb,1,nk1,nn1,XIDL(1,nde))
C                ELSE !Normal Simplex
C                  PD(ns1) = PSI2(IBT(1,1,nb),INP(1,1,nb),nb,
C     '              1,nn1,XIDL(1,nde))
C                ENDIF
C              ELSE !Lagrange/Hermite Tensor product
C cpb 2/2/97 Adding sector elements
C                IF(SECTOR) THEN
C                  PD(ns1) = PSI5(IBT(1,1,nb),IDO(1,1,0,nb),
C     '              INP(1,1,nb),nb,1,nk1,nn1,XIDL(1,nde))
C                ELSE
C                  PD(ns1) = PSI1(IBT(1,1,nb),IDO(1,1,0,nb),
C     '              INP(1,1,nb),nb,1,nk1,nn1,XIDL(1,nde))
C                ENDIF
C              ENDIF
              
            ENDDO !nk1
          ENDDO !nn1
C SAB 21Jan01 Add in the products of the basis functions
          nhs1 = nhs1_for_nhj1
          DO ns1=1,NST(nb)
            nhs1=nhs1+1
            nhs2=0
            DO nhj2=1,NUM_FIT(njj) !columns
              nhx2 = NLH_FIT(nhj2,3,njj)
              nh2=NH_LOC(nhx2,nx)
              DO ns2=1,NST(nb)
                nhs2=nhs2+1
C KFA 17-DEC-2001 Added tag fit flag, I'm sure there is a better
C  way to do it but this is simple, but might change.
                IF(TAG_FIT.EQV..TRUE.) THEN
                  ES(nhs1,nhs2)=ES(nhs1,nhs2)+PD(ns1)*PD(ns2)*
     '              WDL(nh1,nde)*WDL(nh2,nde)*
     '              SE(ns1,nb)*SE(ns2,nb)
                ELSE
                  IF(nhj2.EQ.nhj1) THEN !to avoid coupling for now
                    IF(IN_PLANE) THEN
                      IF(KTYP11.EQ.1) THEN
                        ES(nhs1,nhs2)=ES(nhs1,nhs2)+PD(ns1)*PD(ns2)
     '                    *SE(ns1,nb)*SE(ns2,nb)
                      ELSE
                        ES(nhs1,nhs2)=ES(nhs1,nhs2)+PD(ns1)*PD(ns2)
     '                    *SF(ns1,nb)*SF(ns2,nb)
                      ENDIF
                    ELSE
                      IF(KTYP11.EQ.1) THEN
                        ES(nhs1,nhs2)=ES(nhs1,nhs2)+PD(ns1)*PD(ns2)
     '                    *WDL(nh1,nde)*SE(ns1,nb)*SE(ns2,nb)
                      ELSE
                        ES(nhs1,nhs2)=ES(nhs1,nhs2)+PD(ns1)*PD(ns2)
     '                    *WDL(nh1,nde)*SF(ns1,nb)*SF(ns2,nb)
                      ENDIF
                    ENDIF !IN_PLANE
                  ENDIF !nhj2=nhj1
                ENDIF
              ENDDO !ns2
            ENDDO !nhj2
          ENDDO !ns1
        ENDDO !nde

C SAB 21Jan01 Add in the smoothing contributions
        ns1=0
        nhs1 = nhs1_for_nhj1
        NNTB=NNT(nb)
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
                      ELSE !KTYP12==2 smooth deviation from initial field
C cpb 13/5/98 Adding the Jacobian to the Sobolev integration
C SAB 13Sep00 Shifted the calculation of the Jacobian to the top
                        DO ng=1,NGT(nb)
                          SUM3=0.0d0
                          DO nu=2,NUTB
                            SUM3=SUM3+TEMP_PG(ng,nu,ns1)*
     '                        TEMP_PG(ng,nu,ns2)*WU(nu)
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

C JWF 10/08/04 Adding host-mesh volume constraint stiffness term.
                    SUMV=0.0d0 ! initialiase                  
                    IF(KTYP1V.EQ.1) THEN ! Add volume constraint

                      DO ng=1,NGT(nb)
                 
                        VAR=(PG(ns1,2,ng,nb)*XGT(2,4,ng)*XGT(3,7,ng))
     &                    +(PG(ns1,4,ng,nb)*XGT(1,2,ng)*XGT(3,7,ng))
     &                    +(PG(ns1,7,ng,nb)*XGT(1,2,ng)*XGT(2,4,ng))
     &                    -(PG(ns1,2,ng,nb)*XGT(2,7,ng)*XGT(3,4,ng))
     &                    -(PG(ns1,4,ng,nb)*XGT(1,2,ng)*XGT(2,7,ng))
     &                    -(PG(ns1,7,ng,nb)*XGT(1,2,ng)*XGT(3,2,ng))
     &                    -(PG(ns1,2,ng,nb)*XGT(3,7,ng)*XGT(1,4,ng))
     &                    -(PG(ns1,4,ng,nb)*XGT(2,2,ng)*XGT(3,7,ng))
     &                    -(PG(ns1,7,ng,nb)*XGT(2,2,ng)*XGT(1,4,ng))
     &                    +(PG(ns1,2,ng,nb)*XGT(1,7,ng)*XGT(3,4,ng))
     &                    +(PG(ns1,4,ng,nb)*XGT(2,2,ng)*XGT(1,7,ng))
     &                    +(PG(ns1,7,ng,nb)*XGT(2,2,ng)*XGT(3,4,ng))
     &                    +(PG(ns1,2,ng,nb)*XGT(2,7,ng)*XGT(1,4,ng))
     &                    +(PG(ns1,4,ng,nb)*XGT(3,2,ng)*XGT(2,7,ng))
     &                    +(PG(ns1,7,ng,nb)*XGT(3,2,ng)*XGT(1,4,ng))
     &                    -(PG(ns1,2,ng,nb)*XGT(1,7,ng)*XGT(2,4,ng))
     &                    -(PG(ns1,4,ng,nb)*XGT(3,2,ng)*XGT(1,7,ng))
     &                    -(PG(ns1,7,ng,nb)*XGT(3,2,ng)*XGT(2,4,ng))

                        INC=(PG(ns2,2,ng,nb)*XGT(2,4,ng)*XGT(3,7,ng))
     &                    +(PG(ns2,4,ng,nb)*XGT(1,2,ng)*XGT(3,7,ng))
     &                    +(PG(ns2,7,ng,nb)*XGT(1,2,ng)*XGT(2,4,ng))
     &                    -(PG(ns2,2,ng,nb)*XGT(2,7,ng)*XGT(3,4,ng))
     &                    -(PG(ns2,4,ng,nb)*XGT(1,2,ng)*XGT(2,7,ng))
     &                    -(PG(ns2,7,ng,nb)*XGT(1,2,ng)*XGT(3,2,ng))
     &                    -(PG(ns2,2,ng,nb)*XGT(3,7,ng)*XGT(1,4,ng))
     &                    -(PG(ns2,4,ng,nb)*XGT(2,2,ng)*XGT(3,7,ng))
     &                    -(PG(ns2,7,ng,nb)*XGT(2,2,ng)*XGT(1,4,ng))
     &                    +(PG(ns2,2,ng,nb)*XGT(1,7,ng)*XGT(3,4,ng))
     &                    +(PG(ns2,4,ng,nb)*XGT(2,2,ng)*XGT(1,7,ng))
     &                    +(PG(ns2,7,ng,nb)*XGT(2,2,ng)*XGT(3,4,ng))
     &                    +(PG(ns2,2,ng,nb)*XGT(2,7,ng)*XGT(1,4,ng))
     &                    +(PG(ns2,4,ng,nb)*XGT(3,2,ng)*XGT(2,7,ng))
     &                    +(PG(ns2,7,ng,nb)*XGT(3,2,ng)*XGT(1,4,ng))
     &                    -(PG(ns2,2,ng,nb)*XGT(1,7,ng)*XGT(2,4,ng))
     &                    -(PG(ns2,4,ng,nb)*XGT(3,2,ng)*XGT(1,7,ng))
     &                    -(PG(ns2,7,ng,nb)*XGT(3,2,ng)*XGT(2,4,ng))

                        SUMV=SUMV+(VOL_WT*VAR*INC*WG(ng,nb))

                      ENDDO !ng

                    ENDIF !KTYP1V

                    IF(KTYP11.EQ.1) THEN
                      ES(nhs1,nhs2)=ES(nhs1,nhs2)+(SUM2*WU(0)+SUMV)*
     '                  SE(ns1,nb)*SE(ns2,nb)
                    ELSE
                      ES(nhs1,nhs2)=ES(nhs1,nhs2)+SUM2*WU(0)*
     '                  SF(ns1,nb)*SF(ns2,nb)
                    ENDIF
                  ENDIF !nhj2=nhj1
                ENDDO !nk2
              ENDDO !nn2
            ENDDO !nhj2
          ENDDO !nk1
        ENDDO !nn1
      ENDDO !nhj1

      CALL EXITS('ZDES')
      RETURN
 9999 CALL ERRORS('ZDES',ERROR)
      CALL EXITS('ZDES')
      RETURN 1
      END

