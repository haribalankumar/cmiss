      SUBROUTINE RESFUN(COORDS,DEFORMED,PARAMTYPE,IBT,IDO,INP,
     '  ISIZE_PHI,ISIZE_TBH,LD,LDR,LGE,LIST_RESID,
     '  MODE,NAN,NBH,NBHF,NBJ,NBJF,nd0,nd1,NEELEM,NEL,NELIST,NENP,
     '  NFF,NFFACE,NGAP,NGLIST,NHE,NHP,NKB,NKEF,NKH,
     '  NKHE,NKJE,NLL,NLNO,NNB,NNF,NNL,NONL,NONY,
     '  NP_INTERFACE,NP1OPT,NPF,NPL,NPLIST3,NPNE,NPNODE,
     '  NPNY,nr,NRE,NRLIST,NSB,NVHE,NVHP,NVJE,NVJL,NW,nx,
     '  NXI,NYNE,NYNO,NYNP,NYNR,PAOPTY,Z_CONT_LIST,AQ,CE,
     '  CG,CGE,CONY,CP,CURVCORRECT,DL,FEXT,FGRAD,
     '  LAPL,LAPLSQR,
     '  PAOPTI,PG,PHI,PHI_H,RE1,RESID,RESJAC,RG,SE,T_BH,WG,
     '  WK1_INV,WU,XA,XE,XG,XID,XIG,XN,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,
     '  ZD,ZE,ZE1,ZG,ZP,ZP1,COUPLED_RES,FIX,ERROR,*)

C#### Subroutine: RESFUN
C###  Description:
C###    RESFUN returns residuals for optimising material parameters etc.

C cpb 28/3/96 This comment needs to be updated for new iy scheme.
C**** For KTYP27=2 problems there are KTYP28 exptl measurements in fit
C****   & first set of coords & reactions are held in YP(ny,4) &
C****   YP(ny,5) and YP(ny,6/8/..) & YP(ny,7/9/..) hold
C****   remaining coords & reactions, respec.
C**** Note: Computed reactions from ZPRP are stored in YP(ny,1).
C**** For KTYP27=3 problems there are NTOPTI integrated fluxes which
C****   will all be zero when the saturated bulbs' surface has been
C****   located.
C****   For this case the solution for both regions is required.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbst02.cmn'
      INCLUDE 'chmesh0.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'opti00.cmn'
      INCLUDE 'stab00.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'trsf00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISIZE_PHI(2),ISIZE_TBH(2),
     '  LD(NDM),LDR(0:NDM),LGE(NHM*NSM,NRCM),LIST_RESID,
     '  MODE,NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  nd0,nd1,
     '  NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),
     '  NFF(6,NEM),NFFACE(0:NF_R_M,NRM),NGAP(NIM,NBM),NHE(NEM),
     '  NGLIST(0:NGM),NHP(NPM,0:NRM),
     '  NKB(2,2,2,NNM,NBFM),NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NLL(12,NEM),NLNO(NOPM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONL(NLM),NONY(0:NOYM,NYM,NRCM,0:NRM),
     '  NP_INTERFACE(0:NPM,0:3),NP1OPT(NOPM),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPLIST3(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),
     '  nr,NRE(NEM),NRLIST(0:NRM),NSB(NKM,NNM,NBFM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJL(4,NJM,NLM),NW(NEM,3),
     '  nx,NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),PAOPTY(NOPM),
     '  Z_CONT_LIST(NDM,2,7)
      REAL*8 AQ(NMAQM,NQM),CE(NMM,NEM),
     '  CG(NMM,NGM),CGE(NMM,NGM,NEM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM),
     '  CP(NMM,NPM),CURVCORRECT(2,2,NNM,NEM),DL(3,NLM),
     '  FEXT(NIFEXTM,NGM,NEM),FGRAD(*),
     '  LAPL(NY_TRANSFER_M,NY_TRANSFER_M),
     '  LAPLSQR(NY_TRANSFER_M,NY_TRANSFER_M),
     '  PAOPTI(*),PG(NSM,NUM,NGM,NBM),
     '  PHI(NY_TRANSFER_M,NTSM),PHI_H(NY_TRANSFER_M,NTSM),RE1(NSM,NHM),
     '  RESID(*),RESJAC(NREM,*),RG(NGM),SE(NSM,NBFM,NEM),
     '  T_BH(NY_TRANSFER_M,NY_TRANSFER_M),WG(NGM,NBM),
     '  WK1_INV(NY_TRANSFER_M,NY_TRANSFER_M),WU(0:NUM+1,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XID(NIM,NDM),
     '  XIG(NIM,NGM,NBM),XN(NJM,NGM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM),
     '  ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),Z_CONT(NDM,2,67),
     '  ZD(NJM,NDM),ZE(NSM,NHM),ZE1(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM)
      CHARACTER COORDS*9,PARAMTYPE*(*),ERROR*(*)
      LOGICAL COUPLED_RES,DEFORMED,FIX(NYM,NIYFIXM)

!     Local Variables
      INTEGER atime,DERIVATIVE_TYPE,DERIV_PARAM_NUM,ERR,HEART_ny,
     '  i,ig,IY,k,nb,nc,nd,ndr,ne,NEIGHBOURS,NEIGHBOUR1,
     '  NEIGHBOUR2,NE1,NE2,ng,nh,nhx,nj,njj2,
     '  nk,nn,nn1,nl,MAP(4,4),no,no2,noelem,no_interface,nol,nolist,
     '  nonrlist,no_nynr,no_nynr1,noopti,nores,noy,np,npp,
     '  NP1,NP2,NP3,nu,
     '  nr_res,ns,nts,nv,ny,ny1,nyo,
     '  count,resid_index
      REAL*8 CALC_TIME_FROM_SAMPLE,CALC_FORWARD_ACTN,co,closest,
     '  CROSS_OVER_PENALTY,
     '  DEPTH,dist,DIFF,DZDX(3,3),EG(3,3),ENERGY,SMOM(2,6),RHO1,RHO2,
     '  VELOC1,VELOC2,Z(3),DZ(3),SUM1(3),SUM2,SUM3,TOTAL,TOT_LAPL,
     '  TOT_RES,PH3,PL2S1,PL2S3,PST(3),u,VOL1,VOL2,WIDTH1,ZVAL,
     '  PHI_MEAN,PHI_H_MEAN,PHI_NORM,PHI_H_NORM,R(3,3),
     '  SS_RESID,CC_RESID,
     '  PHI_PHI_H_SUM,SS_SUM,CC_SUM,D_PHI_H_MEAN,
     '  RG2D,RGZ,RGZ2D,RM(3,3),RI1,RI2,RI3,
     '  TC(3,3),TG(3,3),TGNG(81),TIME,TN(3,3),TNA,WEIGHT(81),XI(3)
      CHARACTER STRESSTYPE*17
      LOGICAL CALCJAC,CONNECTED,FOUND

!     Functions
      REAL*8 DDOT,PSI

      DATA MAP/3,4,0,0,0,0,3,4,1,0,2,0,0,1,0,2/
      DATA STRESSTYPE/'Total'/

      CALL ENTERS('RESFUN',*9999)

      TOT_RES=0.D0
      nc=1 !Temporary GMH 12/9/95

C!!! LKC 29-NOV-1999 This initialisation is not required for
C!!!   the 'potential' stuff. Bypass if speedup is required.
C!!! KFA  8-MAR-2004 This init stuff is detrimental to problems
C!!!   which progressively fill the residual vector (ie example 723)
C!!!   and possibly Carey's-Holmes stuff aswell.
C!!!   In addition the code isn't general as RESID is not nesscesarily
C!!!   NOT(...) long.

C     Initialise RESID
      IF (.NOT.((KTYP26.EQ.1.AND.KTYP27.EQ.5) .OR. 
     '  (KTYP26.EQ.1.AND.KTYP27.EQ.6))) THEN
        DO no=1,NOT(1,1,nr,nx)
          RESID(no)=0.0d0
        ENDDO                  !no
      ENDIF

      IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
        IF(KTYP27.EQ.1) THEN !Obj func is max princ stress difference
          CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),NPNODE,
     '      nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
          PRSTMAX=0.0d0
          PRSTMIN=1.0d6
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),
     '        NW(ne,1),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '        ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
            IF(ITYP1(nr,nx).EQ.4) THEN
              CALL CPCG(NW(ne,1),NBH(1,1,ne),NPNE(1,1,ne),nr,nx,
     '          CE(1,ne),CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
              nb=NBH(NH_LOC(1,nx),1,ne)
              CALL OPST40(nb,NBH(1,1,ne),NBJ(1,ne),ne,NHE(ne),
     '          NW(ne,1),nx,CG,ENERGY,PG,'BLANKS',WG,
     '          XE,XG,YG(1,1,ne),ZE,ZG,.TRUE.,COORDS,ERROR,*9999)
            ELSE IF(ITYP1(nr,nx).EQ.5) THEN
              nb=NBH(NH_LOC(1,nx),1,ne)
              CALL CPCG(1,nb,NPNE(1,1,ne),nr,nx,
     '          CE(1,ne),CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
              CALL OPC50(IBT,IDO,INP,NAN,NBH(1,1,ne),NBJ(1,ne),ne,
     '          NHE(ne),NPNE(1,1,ne),nr,nx,
     '          CE(1,ne),CG,CP,FEXT(1,1,ne),PG,RG,
     '          XE,XG,XIG(1,1,nb),YG(1,1,ne),ZE,ZG,ERROR,*9999)
            ENDIF
            FUNC=(PRSTMAX-PRSTMIN)**2
          ENDDO
          WRITE(OP_STRING,'('' Max principle stress = '',D11.3)')
     '      PRSTMAX
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Min principle stress = '',D11.3)')
     '      PRSTMIN
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Princ. stress diff.^2= '',D11.3)')
     '      FUNC
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        ELSE IF(KTYP27.EQ.2) THEN !Obj func is sum of squared residuals
          CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,atime,MAQ_ACTIV_TIME,
     '      ERROR,*9999)
          FUNC=0.0d0
          IF(KTYP28.EQ.0) THEN      !residuals calc.d during fit
            IF(ITYP5(nr,nx).EQ.2.AND.ITYP19(nr,nx).EQ.1.AND.
     '        ITYP2(nr,nx).EQ.9) THEN !activation model
              DO nd=1,NDT    !electrode loop
                RESID(nd)=ZD(NJT+1,nd)-AQ(atime,LDR(nd))
                IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP             CRITICAL(RESFUN_1)
                  WRITE(OP_STRING,'('' nd='',I5,'' time='',F8.2,'
     '              //'''  :   nq='',I7,'' time='',F8.2)')
     '              nd,ZD(NJT+1,nd),LDR(nd),AQ(atime,LDR(nd))
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP             END CRITICAL(RESFUN_1)
                ENDIF
              ENDDO !nd
            ELSE !all others
              CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),NPNODE,
     '          nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
              IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const V
C               Put reference state for cavity from YP(ny,10) into
C               ZA1,ZP1 for ZERE55
                CALL YPZP(10,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),
     '            NPNODE,nr,
     '            NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA1,ZP1,ERROR,*9999)
              ENDIF
              CALL ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '          NEELEM,NFF,NFFACE,NGAP,NHE,NHP(1,nr),NKB,NKEF,
     '          NKH(1,1,1,nr),NKHE,NKJE,NNB,NNF,NPF,NPNE,
     '          NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP(1,1,1,nr),NVJE,NW,
     '          nx,NXI,NYNE,NYNP,NYNR(0,0,1,nr),Z_CONT_LIST,
     '          CE,CG,CGE,CP,CURVCORRECT,FEXT,PG,RE1,RG,
     '          SE,WG,
     '          XA,XG,XP,YG,YGF,YP,ZA,ZA1,%VAL(0),Z_CONT,ZE,ZE1,
     '          ZP,ZP1,
     '          %VAL(0),FIX,ERROR,*9999)
              DO no_nynr=1,NYNR(0,1,1,nr) !loop over rows
                ny1=NYNR(no_nynr,1,1,nr) !is row number
                DO noy=1,NONY(0,ny1,1,nr) !loop soln rows assoc with ny1
                  no=NONY(noy,ny1,1,nr) !is row number for ny1
                  co=CONY(noy,ny1,1,nr) !is coupling coeff for ny1
                  RESID(no)=RESID(no)+YP(ny1,4)*co
                ENDDO !noy
              ENDDO !no_nynr
            ENDIF !ityp5(nr)

          ELSE IF(KTYP28.GT.0) THEN !use existing residuals

            IF(ITYP5(nr,nx).EQ.2.AND.ITYP19(nr,nx).EQ.1.AND.
     '        ITYP2(nr,nx).EQ.9) THEN !activation model
            ELSE !all others
              DO k=1,KTYP28
                IF(k.EQ.1) THEN
                  IY=1
                ELSE
                  CALL ASSERT(.FALSE.,'>>Code needs updating for '
     '              //'new YP iy locations',ERROR,*9999)
C!!! This code needs checking for IY locations
                  IY=4+2*k
                ENDIF
                CALL YPZP(IY,NBH,NEELEM,NHE,NHP(1,nr),
     '            NKH(1,1,1,nr),NPNODE,
     '            nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
                IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !cnst V
C                 Put reference state for cavity from YP(ny,10) into
C                 ZA1,ZP1 for ZERE55
                  CALL YPZP(10,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),
     '              NPNODE,nr,
     '              NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA1,ZP1,ERROR,*9999)
                ENDIF
                CALL ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '            NEELEM,NFF,NFFACE,NGAP,NHE,NHP(1,nr),NKB,NKEF,
     '            NKH(1,1,1,nr),NKHE,NKJE,NNB,NNF,NPF,NPNE,
     '            NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP(1,1,1,nr),NVJE,NW,
     '            nx,NXI,NYNE,NYNP,NYNR(0,0,1,nr),Z_CONT_LIST,
     '            CE,CG,CGE,CP,CURVCORRECT,FEXT,PG,RE1,RG,
     '            SE,WG,
     '            XA,XG,XP,YG,YGF,YP,ZA,ZA1,%VAL(0),Z_CONT,ZE,ZE1,
     '            ZP,ZP1,
     '            %VAL(0),FIX,ERROR,*9999)
                DO no_nynr=1,NYNR(0,1,1,nr) !loop over rows
                  ny1=NYNR(no_nynr,1,1,nr) !is row number
                  DO noy=1,NONY(0,ny1,1,nr) !loop sol rows assc with ny1
                    no=NONY(noy,ny1,1,nr) !is row number for ny1
                    co=CONY(noy,ny1,1,nr) !is coupling coeff for ny1
                    RESID(no)=RESID(no)+YP(ny1,4)*co
                  ENDDO !noy
                ENDDO !no_nynr
              ENDDO !k
            ENDIF !ktyp28
          ENDIF !ityp5(nr) and ityp2(nr)

        ELSE IF(KTYP27.EQ.3) THEN
C ***     Objective function is difference between FE flux and BE flux
C ***     nr refers to the BE region.

          ERROR='>> Not implemented. See backup version.'
          GOTO 9999

        ELSE IF(KTYP27.EQ.4) THEN
C ***     Objective function is diff betw calc.d height & bdry position
          CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),NPNODE,
     '      nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
          FUNC=0.0d0
          DO noopti=1,NTOPTI
            np1=NP1OPT(noopti)
            FUNC=FUNC+(ZP(1,1,1,np1,nc)*40.0d0+XP(1,1,2,np1))**2
          ENDDO

        ELSE IF(KTYP27.EQ.6) THEN
C ***     Objective function is fluid interface dynamic bdry condition
          DO nr=1,2 !!! Why does this go up to 2 only ?
            CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),NPNODE,
     '        nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
          ENDDO
          NE1=NEELEM(1,1) !is first element of region 1
          NE2=NEELEM(1,2) !is first element of region 2
          RHO1=CE(1,NE1)  !is density for region 1
          RHO2=CE(1,NE2)  !is density for region 2
          VELOC1=CE(2,NE1)  !is velocity for region 1
          VELOC2=CE(2,NE2)  !is velocity for region 2
          WRITE(OP_STRING,'('' rho1='',e12.3,'' rho2='',e12.3,'
     '      //''' veloc1='',e12.3,'' veloc2='',e12.3)')
     '      rho1,rho2,veloc1,veloc2
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          FUNC=0.0D0
          DO no_interface=1,NP_INTERFACE(0,3)
            NP1=NP_INTERFACE(no_interface,1)
            NP2=NP_INTERFACE(no_interface,2)
            NP3=NP_INTERFACE(no_interface,3)
            FUNC=FUNC+(RHO1*((ZP(1,1,1,NP1,nc)-ZP(1,1,2,NP1,nc))/TINCR
     '        +VELOC1*ZP(2,1,1,NP1,nc)+G_ACCEL*ZP(1,1,1,NP3,nc))
     '        -RHO2*((ZP(1,1,1,NP2,nc)-ZP(1,1,2,NP2,nc))/TINCR
     '        +VELOC2*ZP(2,1,1,NP2,nc)+G_ACCEL*ZP(1,1,1,NP3,nc)))**2
            WRITE(OP_STRING,'('' interface node '',i2,'' func= '','
     '        //'e12.3)') no_interface,func
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF

      ELSE IF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
        IF(COUPLED_RES) THEN
          DO nonrlist=1,NRLIST(0)
            nr_res=NRLIST(nonrlist)
            CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr_res),NKH(1,1,1,nr_res),
     '        NPNODE,nr_res,NVHP(1,1,1,nr_res),nx,NYNE,NYNP,YP,ZA,ZP,
     '        ERROR,*9999)
C new MPN 24Dec1999: bug - used wrong nr here
            IF(ITYP2(nr_res,nx).EQ.8.AND.ITYP3(nr_res,nx).EQ.3) THEN
C             constant volume constraint
C             Put reference state for cavity from YP(ny,10) into
C             ZA1,ZP1 for ZERE55
              CALL YPZP(10,NBH,NEELEM,NHE,NHP(1,nr_res),
     '          NKH(1,1,1,nr_res),NPNODE,nr_res,
     '          NVHP(1,1,1,nr_res),nx,NYNE,NYNP,YP,ZA1,ZP1,ERROR,*9999)
C end new
C OLD
C            IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const V
CC             Put reference state for cavity from YP(ny,10) into
CC             ZA1,ZP1 for ZERE55
C              CALL YPZP(10,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),
C     '          NPNODE,nr,
C     '          NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA1,ZP1,ERROR,*9999)
C end OLD
            ENDIF
            CALL ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '        NEELEM,NFF,
     '        NFFACE,NGAP,NHE,NHP(1,nr_res),NKB,NKEF,
     '        NKH(1,1,1,nr_res),NKHE,NKJE,NNB,NNF,NPF,NPNE,
     '        NPNODE,NPNY,
     '        nr_res,NRE,NSB,NVHE,NVHP(1,1,1,nr_res),NVJE,NW,nx,NXI,
     '        NYNE,NYNP,NYNR(0,0,1,nr_res),Z_CONT_LIST,CE,CG,
     '        CGE,
     '        CP,CURVCORRECT,FEXT,PG,RE1,RG,SE,WG,XA,XG,
     '        XP,YG,YGF,YP,
     '        ZA,ZA1,%VAL(0),Z_CONT,ZE,ZE1,ZP,ZP1,%VAL(0),FIX,
     '        ERROR,*9999)
          ENDDO !nonrlist (nr)
        ELSE
          CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),
     '      NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,
     '      ERROR,*9999)
          IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const V
C           Put reference state for cavity from YP(ny,10) into
C           ZA1,ZP1 for ZERE55
            CALL YPZP(10,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),
     '        NPNODE,nr,
     '        NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA1,ZP1,ERROR,*9999)
          ENDIF
          CALL ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,NFF,
     '      NFFACE,NGAP,NHE,NHP(1,nr),NKB,NKEF,NKH(1,1,1,nr),
     '      NKHE,NKJE,NNB,NNF,NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,
     '      NVHP(1,1,1,nr),NVJE,NW,nx,NXI,NYNE,NYNP,NYNR(0,0,1,nr),
     '      Z_CONT_LIST,
     '      CE,CG,CGE,CP,CURVCORRECT,FEXT,PG,RE1,RG,SE,WG,
     '      XA,XG,XP,YG,YGF,YP,ZA,ZA1,%VAL(0),Z_CONT,ZE,ZE1,ZP,ZP1,
     '      %VAL(0),FIX,
     '      ERROR,*9999)
        ENDIF !COUPLED_RES
        IF(KTYP27.EQ.7) THEN !Aero wake press diff & sail stress
          nr=2
C         DO noopti=NL_WAKE(0,1)+2,NTOPTI
          DO noopti=N_OPTI(1)+1,NTOPTI
            ny=NYNO(1,noopti,2,nr)
            RESID(noopti)=YP(ny,4)
          ENDDO

        ELSE
          DO no_nynr=1,NYNR(0,1,1,nr) !loop over rows
            ny1=NYNR(no_nynr,1,1,nr) !is row number
            DO noy=1,NONY(0,ny1,1,nr) !loop soln rows assoc with ny1
              no=NONY(noy,ny1,1,nr) !is row number for ny1
              co=CONY(noy,ny1,1,nr) !is coupling coeff for ny1
              RESID(no)=RESID(no)+YP(ny1,4)*co
            ENDDO !noy
          ENDDO !no_nynr
        ENDIF

      ELSE IF(PARAMTYPE(1:12).EQ.'DATA_FITTING') THEN
C       Copy in the current estimates of XP and DL from PAOPTI
        DO noopti=1,NTOPTI
          IF(PAOPTY(noopti).EQ.1) THEN !Parameter is geometric dof
            DO nyo=1,NYNO(0,noopti,2,nr)
              ny=NYNO(nyo,noopti,2,nr)
              nk=NPNY(1,ny,0)
              nv=NPNY(2,ny,0)
              nj=NPNY(3,ny,0)
              np=NPNY(4,ny,0)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              IF((nk.EQ.2.OR.nk.EQ.3).AND.KTYP1B.EQ.2) THEN !derivative as angle
C               Only calculate if this is the first of the angles
                IF(NONY(1,ny,2,nr).EQ.noopti) THEN
                  IF(NYNO(0,noopti,2,nr).EQ.3) THEN !two angles
                    no2=NONY(2,ny,2,nr) !get the second angle
                    IF(nyo.EQ.1) THEN !x
                      XP(nk,nv,nj,np)=DSIN(PAOPTI(no2))*
     '                  DCOS(PAOPTI(noopti))
                    ELSEIF(nyo.EQ.2) THEN !y
                      XP(nk,nv,nj,np)=DSIN(PAOPTI(no2))*
     '                  DSIN(PAOPTI(noopti))
                    ELSEIF(nyo.EQ.3) THEN !z
                      XP(nk,nv,nj,np)=DCOS(PAOPTI(no2))
                    ELSE
                      ERROR='>>Only three components for two angles'
                      GOTO 9999
                    ENDIF
                  ELSEIF(NYNO(0,noopti,2,nr).EQ.2) THEN !one angle
                    IF(nyo.EQ.1) THEN !equivalent of x
                      XP(nk,nv,nj,np)=DCOS(PAOPTI(noopti))
                    ELSEIF(nyo.EQ.2) THEN !equivalent of y
                      XP(nk,nv,nj,np)=DSIN(PAOPTI(noopti))
                    ELSE
                      ERROR='>>Only two components for one angle'
                      GOTO 9999
                    ENDIF
                  ELSE
                    ERROR='>>Must have at least 2 components for angles'
                    GOTO 9999
                  ENDIF
                ENDIF !first angle
              ELSE
                XP(nk,nv,nj,np)=PAOPTI(noopti)
              ENDIF
            ENDDO !nyo
          ELSE IF(PAOPTY(noopti).EQ.2) THEN !Parameter is a line
            nl=NLNO(noopti)
            DL(1,nl)=PAOPTI(noopti)
            DL(2,nl)=DL(1,nl)
            DL(3,nl)=DL(1,nl)
            IF(JTYP2B.EQ.1.AND.NPL(4,0,nl).GT.0) THEN
              DL(1,NPL(4,0,nl))=-1.0d0*DL(2,nl)
              DL(2,NPL(4,0,nl))=-1.0d0*DL(1,nl)
              DL(3,NPL(4,0,nl))=DL(3,nl)
            ENDIF
          ELSE IF(PAOPTY(noopti).EQ.3) THEN !Parameter is mate
          ELSE
            WRITE(ERROR,'(''>>Invalid PAOPTY('',I3,'')='',I3,'
     '        //''' for data fitting'')') noopti,PAOPTY(noopti)
C           ERROR=' >>Invalid PAOPTY for data fitting'
            GOTO 9999
          ENDIF
        ENDDO
        IF(G1SCALING) THEN
          DO nl=1,NLT
            CALL ARCLEN(IDO,NBJ,NEL(0,nl),nl,
     '        NPL(1,0,nl),NPNE,NVJL(1,1,nl),DL,XP,ERROR,*9999)
          ENDDO
C         Update scale factor arrays
          DO nb=1,NBFT
            CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,DL,
     '        SE,ERROR,*9999)
          ENDDO
        ENDIF
C       Evaluate the residuals
        IF(KTYP29B.EQ.2) THEN !one residual
          nores=1
          RESID(nores)=0.0d0
          IF(MODE.EQ.1.OR.MODE.EQ.2) THEN
            DO noopti=1,NTOPTI
              RESJAC(nores,noopti)=0.0d0
            ENDDO
          ENDIF !mode
        ELSE
          nores=0 !initializing for multiple residuals
        ENDIF !residual type
        DO ndr=1,LDR(0)
          nd=LDR(ndr)
          IF(nd.LT.nd0.OR.nd.GT.nd1) THEN
            IF(KTYP29B.EQ.1) THEN
              nores=nores+NJ_LOC(NJL_GEOM,0,nr)
            ELSE
              nores=nores+1
            ENDIF
          ELSE
            IF(LD(nd).NE.0) THEN
              ne=LD(nd)
              IF(KTYP29B.EQ.3) THEN !residual for each projection
                nores=nores+1
                RESID(nores)=0.0d0
                IF(MODE.EQ.1.OR.MODE.EQ.2) THEN
C KAT 12Nov98:    Initialize RESJAC as there may be more than one element
C                 dof affecting a data point projection mapping to the same
C                 global dof. e.g. a collapsed node
                  DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                    nj=NJ_LOC(NJL_GEOM,njj2,nr)
                    nb=NBJ(nj,ne)
                    DO nn=1,NNT(nb)
                      np=NPNE(nn,nb,ne)
                      nv=NVJE(nn,nb,nj,ne)
                      DO nk=1,NKT(nn,nb)
                        ny=NYNP(nk,nv,nj,np,0,1,nr)
                        DO noy=1,NONY(0,ny,2,nr)
                          noopti=NONY(noy,ny,2,nr)
                          RESJAC(nores,noopti)=0.0d0
                        ENDDO !noy
                      ENDDO !nk
                    ENDDO !nn
                  ENDDO !nj
                ENDIF !MODE
              ENDIF !residual for each projection
              IF(DEFORMED) THEN
                CALL ZPZE(NBH(1,1,ne),1,NHE(ne),
     '            NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     '            NVHE(1,1,1,ne),NW(ne,1),nx,
     '            CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '            ZA(1,1,1,ne),XE,ZP,ERROR,*9999)
              ELSE
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA,XE,XP,ERROR,*9999)
              ENDIF
              DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,njj2,nr)
                nb=NBJ(nj,ne)
                ns=0
C               Store basis functions in RE1
                DO nn=1,NNT(nb)
                  DO nk=1,NKT(nn,nb)
                    ns=ns+1
                    RE1(ns,1)=PSI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                nb,nn,nk,1,XID(1,nd))
                  ENDDO !nk
                ENDDO !nn
                Z(nj)=DDOT(ns,RE1(1,1),1,XE(1,nj),1)
                DZ(nj)=Z(nj)-ZD(nj,nd)
                IF(KTYP29B.EQ.1) THEN !residual for each component
                  nores=nores+1
                  RESID(nores)=DZ(nj)
                ELSE
                  RESID(nores)=RESID(nores)+DZ(nj)**2
                ENDIF
                IF(MODE.EQ.1.OR.MODE.EQ.2) THEN
                  IF(KTYP29B.EQ.1) THEN !residual for each component
C KAT 12Nov98:    Initialize RESJAC as there may be more than one element
C                 dof affecting a data point projection mapping to the same
C                 global dof. e.g. a collapsed node
                    DO nn=1,NNT(nb)
                      np=NPNE(nn,nb,ne)
                      nv=NVJE(nn,nb,nj,ne)
                      DO nk=1,NKT(nn,nb)
                        ny=NYNP(nk,nv,nj,np,0,1,nr)
                        DO noy=1,NONY(0,ny,2,nr)
                          noopti=NONY(noy,ny,2,nr)
                          RESJAC(nores,noopti)=0.0d0
                        ENDDO !noy
                      ENDDO !nk
                    ENDDO !nn
                  ENDIF !residual for each component
C                 Evaluate the objective jacobian
                  ns=0
                  DO nn=1,NNT(nb)
                    np=NPNE(nn,nb,ne)
                    nv=NVJE(nn,nb,nj,ne)
                    DO nk=1,NKT(nn,nb)
                      ns=ns+1
                      ny=NYNP(nk,nv,nj,np,0,1,nr)
                      SUM3=RE1(ns,1)*SE(ns,nb,ne)
C                      SUM3=PSI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
C     '                  nn,nk,1,XID(1,nd))*SE(ns,nb,ne)
                      DO noy=1,NONY(0,ny,2,nr)
                        noopti=NONY(noy,ny,2,nr)
                        IF(KTYP1B.EQ.2.AND.(nk.EQ.2.OR.nk.EQ.3))THEN
C                         derivative as angle
                          IF(NYNO(0,noopti,2,nr).EQ.3) THEN
C                           two angles
                            IF(noy.EQ.1) THEN !del/del theta
                              no2=NONY(2,ny,2,nr) !get phi
                              IF(NYNO(1,noopti,2,nr).EQ.ny) THEN !x
                                SUM3=SUM3*DSIN(PAOPTI(no2))*
     '                            DSIN(PAOPTI(noopti))
                              ELSEIF(NYNO(2,noopti,2,nr).EQ.ny) THEN !y
                                SUM3=SUM3*DSIN(PAOPTI(no2))*
     '                            DCOS(PAOPTI(noopti))
                              ELSEIF(NYNO(3,noopti,2,nr).EQ.ny) THEN !z
                                SUM3=0.0d0 !add nothing
                              ELSE
                                ERROR='>>Invalid component'
                                GOTO 9999
                              ENDIF
                            ELSEIF(noy.EQ.2) THEN !del/del phi
                              no2=NONY(1,ny,2,nr) !get theta
                              IF(NYNO(1,noopti,2,nr).EQ.ny) THEN !x
                                SUM3=SUM3*DCOS(PAOPTI(noopti))*
     '                            DCOS(PAOPTI(no2))
                              ELSEIF(NYNO(2,noopti,2,nr).EQ.ny) THEN !y
                                SUM3=SUM3*DCOS(PAOPTI(noopti))*
     '                            DSIN(PAOPTI(no2))
                              ELSEIF(NYNO(3,noopti,2,nr).EQ.ny) THEN !z
                                SUM3=SUM3*DSIN(PAOPTI(noopti))
                              ELSE
                                ERROR='>>Invalid component'
                                GOTO 9999
                              ENDIF
                            ELSE
                              ERROR='>>Only two angles '
     '                          //'for three components'
                              GOTO 9999
                            ENDIF
                          ELSEIF(NYNO(0,noopti,2,nr).EQ.2) THEN !one angle
                            IF(noy.EQ.1) THEN !del/del theta
                              IF(NYNO(1,noopti,2,nr).EQ.ny) THEN !equiv of x
                                SUM3=SUM3*DSIN(PAOPTI(noopti))
                              ELSEIF(NYNO(2,noopti,2,nr).EQ.ny) THEN !equiv of y
                                SUM3=SUM3*DCOS(PAOPTI(noopti))
                              ELSE
                                ERROR='>>Invalid component'
                                GOTO 9999
                              ENDIF
                            ELSE
                              ERROR='>>Only one angle '
     '                          //'for two components'
                              GOTO 9999
                            ENDIF
                          ENDIF !first angle
                        ENDIF !deriv as angle
                        IF(KTYP29B.NE.1) SUM3=DZ(nj)*SUM3
                        IF(KTYP29B.EQ.3) SUM3=2.0d0*SUM3
                        RESJAC(nores,noopti)=
     '                    RESJAC(nores,noopti)+SUM3
                      ENDDO !noy
                    ENDDO !nk
                  ENDDO !nn
                ENDIF !MODE
              ENDDO !nj
              IF(G1SCALING.AND.(MODE.EQ.1.OR.MODE.EQ.2)) THEN
C CPB 6/6/93 WARNING : Assuming same basis function in each geometric
C                      direction
C             Calculate the jacobian wrt to the line length
                nb=NBJ(1,ne)
                IF(NIT(nb).EQ.1) THEN !1D Element
                  nl=NLL(1,ne)
                  IF(nl.NE.0) THEN
                    noopti=NONL(nl)
                    IF(noopti.GT.0) THEN !Line is to be optimised
                      SUM2=0.0d0
                      DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                        nj=NJ_LOC(NJL_GEOM,njj2,nr)
                        SUM1(nj)=0.0d0
                        DO nn=1,NNT(nb)
                          SUM1(nj)=SUM1(nj)+
     '                      PH3(nn,2,1,XID(1,nd))*
     '                      XP(2,1,nj,NPL(1+nn,1,nl))
                        ENDDO !nn
                        SUM1(nj)=SUM1(nj)*DZ(nj)
                        SUM2=SUM2+SUM1(nj)
                      ENDDO !nj
                      RESJAC(nores,noopti)=2.0d0*SUM2
                    ENDIF
                  ENDIF
                ELSE IF(NIT(nb).EQ.2) THEN !2D Element
                  IF(IBT(1,1,nb).EQ.2.AND.IBT(1,2,nb).EQ.2) THEN !BiHermite
                    DO nol=1,4
                      nl=NLL(nol,ne)
                      IF(nl.NE.0) THEN
                        noopti=NONL(nl)
                        IF(noopti.GT.0) THEN !Line is optimised
                          SUM2=0.0d0
                          DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                            nj=NJ_LOC(NJL_GEOM,njj2,nr)
                            SUM1(nj)=0.0d0
                            DO nn1=1,2
                              np=NPL(1+nn1,1,nl)
                              nk=NPL(4,1,nl)
C***                            Find local node nn of global node np on
C***                            element ne
                              nn=1
                              FOUND=.FALSE.
                              DO WHILE((nn.LE.NNT(nb)).AND.(.NOT.FOUND))
                                IF(np.EQ.NPNE(nn,nb,ne)) THEN
                                  FOUND=.TRUE.
                                ELSE
                                  nn=nn+1
                                ENDIF
                              ENDDO
                              IF(.NOT.FOUND)THEN
                                ERROR='>>Could not find local node'
                                GOTO 9999
                              ENDIF
                              SUM1(nj)=SUM1(nj)+
     '                          PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,
     '                          XID(1,nd))*
     '                          PH3(INP(nn,2,nb),IDO(nk,nn,2,nb),1,
     '                          XID(2,nd))*XP(nk,1,nj,np)
                              nk=4
                              SUM1(nj)=SUM1(nj)+
     '                          PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,
     '                          XID(1,nd))*
     '                          PH3(INP(nn,2,nb),IDO(nk,nn,2,nb),1,
     '                          XID(2,nd))*
     '                          XP(nk,1,nj,np)*DL(1,NLL(MAP(nn,nol),ne))
                            ENDDO !nn1
                            SUM1(nj)=SUM1(nj)*DZ(nj)
                            SUM2=SUM2+SUM1(nj)
                          ENDDO !nj
                          RESJAC(nores,noopti)=2.0d0*SUM2
                        ENDIF
                      ENDIF
                    ENDDO !nol
                  ELSE IF(IBT(1,1,nb).EQ.3.AND.IBT(1,2,nb).EQ.3) THEN !Hermite-Simplex
                    DO nol=1,3
                      nl=NLL(nol,ne)
                      IF(nl.NE.0) THEN
                        noopti=NONL(nl)
                        IF(noopti.GT.0) THEN !Line is optimised
                          SUM2=0.0d0
                          DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                            nj=NJ_LOC(NJL_GEOM,njj2,nr)
                            SUM1(nj)=0.0d0
                            IF(NPL(1,nj,nl).EQ.4) THEN !Cubic Hermite line
                              IF(NKT(1,nb).EQ.1) THEN !Apex at node 1
                                DO nn=2,3
                                  np=NPL(1+nn-1,1,nl)
                                  nk=NPL(5,1,nl)
                                  SUM1(nj)=SUM1(nj)+PH3(INP(nn,1,nb),
     '                              IDO(nk,nn,1,nb),1,XID(1,nd))*
     '                              PL2S1(INP(nn,2,nb),IDO(nk,nn,2,nb),
     '                              1,XID(2,nd))*XP(nk,1,nj,np)
                                  nk=4
                                  SUM1(nj)=SUM1(nj)+PH3(INP(nn,1,nb),
     '                              IDO(nk,nn,1,nb),1,XID(1,nd))*
     '                              PL2S1(INP(nn,2,nb),IDO(nk,nn,2,nb),
     '                              1,XID(2,nd))*XP(nk,1,nj,np)*
     '                              DL(1,NLL(nn,ne))
                                ENDDO
                              ELSE IF(NKT(3,nb).EQ.1) THEN !Apex at node 3
                                DO nn=1,2
                                  np=NPL(1+nn,1,nl)
                                  nk=NPL(4,1,nl)
                                  SUM1(nj)=SUM1(nj)+PH3(INP(nn,1,nb),
     '                              IDO(nk,nn,1,nb),1,XID(1,nd))*
     '                              PL2S3(INP(nn,2,nb),IDO(nk,nn,2,nb),
     '                              1,XID(2,nd))*XP(nk,1,nj,np)
                                  nk=4
                                  SUM1(nj)=SUM1(nj)+PH3(INP(nn,1,nb),
     '                              IDO(nk,nn,1,nb),1,XID(1,nd))*
     '                              PL2S3(INP(nn,2,nb),IDO(nk,nn,2,nb),
     '                              1,XID(2,nd))*XP(nk,1,nj,np)*
     '                              DL(1,NLL(nn+1,ne))
                                ENDDO
                              ELSE
                                ERROR='>>NKT incorrect for Hermite-'
     '                            //'Simplex element'
                                GOTO 9999
                              ENDIF
                            ELSE IF(NPL(1,nj,nl).EQ.6) THEN !Apex at node 1
                              np=NPL(3,1,nl)
                              nk=NPL(5,1,nl)
                              IF(NPNE(2,nb,ne).eq.np) THEN
                                nn=2
                              ELSE
                                nn=3
                              ENDIF
                              SUM1(nj)=SUM1(nj)+
     '                          PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,
     '                          XID(1,nd))*
     '                          PL2S1(INP(nn,2,nb),IDO(nk,nn,2,nb),1,
     '                          XID(2,nd))*XP(nk,1,nj,np)
                              nk=4
                              SUM1(nj)=SUM1(nj)+
     '                          PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,
     '                          XID(1,nd))*
     '                          PL2S1(INP(nn,2,nb),IDO(nk,nn,2,nb),1,
     '                          XID(2,nd))*XP(nk,1,nj,np)*DL(1,NLL(3,
     '                          ne))
                            ELSE IF(NPL(1,nj,nl).EQ.7) THEN !Apex at node 3
                              np=NPL(2,1,nl)
                              nk=NPL(4,1,nl)
                              IF(NPNE(1,nb,ne).eq.np) THEN
                                nn=1
                              ELSE
                                nn=2
                              ENDIF
                              SUM1(nj)=SUM1(nj)+PH3(INP(nn,1,nb),
     '                          IDO(nk,nn,1,nb),1,XID(1,nd))*
     '                          PL2S3(INP(nn,2,nb),IDO(nk,nn,2,nb),1,
     '                          XID(2,nd))*XP(nk,1,nj,np)
                              nk=4
                              SUM1(nj)=SUM1(nj)+PH3(INP(nn,1,nb),
     '                          IDO(nk,nn,1,nb),1,XID(1,nd))*
     '                          PL2S3(INP(nn,2,nb),IDO(nk,nn,2,nb),1,
     '                          XID(2,nd))*XP(nk,1,nj,np)*DL(1,NLL(1,
     '                          ne))
                            ELSE
                              ERROR='>>Unknown line type'
                              GOTO 9999
                            ENDIF
                            SUM1(nj)=SUM1(nj)*DZ(nj)
                            SUM2=SUM2+SUM1(nj)
                          ENDDO !nj
                          RESJAC(nores,noopti)=2.0d0*SUM2
                        ENDIF
                      ENDIF
                    ENDDO !nol
                  ELSE
                    ERROR='>>Element type not yet implemented'
                    GOTO 9999
                  ENDIF
                ELSE
                  ERROR='>>Current value of NITB not yet implemented'
                  GOTO 9999
                ENDIF
              ENDIF !G1SCALING, mode
            ELSE
              ERROR='>>LD(nd)=0 redefine optimisation problem'
              GOTO 9999
            ENDIF
          ENDIF
        ENDDO
        IF(KTYP29B.EQ.2) THEN !one residual
          RESID(nores)=DSQRT(RESID(nores))
          IF(MODE.EQ.1.OR.MODE.EQ.2) THEN
            DO noopti=1,NTOPTI
              RESJAC(nores,noopti)=RESJAC(nores,noopti)/RESID(nores)
            ENDDO
          ENDIF !mode
        ENDIF !residual type

        IF(KTYP12.EQ.1) THEN !Sobolev Smoothing
          CALCJAC=.TRUE.
          CALL SOBOLEV(IBT,IDO,INP,NBJ,NEELEM,NKJE,NLL,NONL,
     '      NONY(0,1,1,nr),NPF,NPNE,nr,NRE,NVJE,
     '      NYNO(0,1,1,nr),NYNP,DL,RESJAC,
     '      FGRAD,PAOPTI,
     '      PG,RG,SE,RESID(NT_RES),WG,WU,XA,XE,XG,XIG,XP,CALCJAC,
     '      ERROR,*9999)
        ENDIF
      ELSE IF(PARAMTYPE(1:15).EQ.'TORSO_CUSTOMISE') THEN

        IF(CUSTOMISATION_TYPE.EQ.1.OR.CUSTOMISATION_TYPE.EQ.3) THEN
          DO nores=1,NT_RES
            ZVAL=CIRMEASURE(2,nores)
            CALL CHMESH_CALC_ARCLENGTHS(IBT,IDO,INP,NBJ,NEELEM,NGAP,
     '        NPNE,1,DL,SE,TOTAL,WG,XIG,XP,ZVAL,ERROR,*9999)
            RESID(nores)=TOTAL-CIRMEASURE(1,nores)
          ENDDO !nores
        ELSEIF(CUSTOMISATION_TYPE.EQ.2) THEN
          DO noelem=1,NEELEM(0,vreg1)
            NELIST(noelem)=NEELEM(noelem,vreg1)
          ENDDO
          NELIST(0)=NEELEM(0,vreg1)
          CALL VOLUME(NBJ,NELIST,NKJE,NPF,NPNE,vreg1,NVJE,NW,
     '      PG,RG,SE,VOL1,WG,XA,XE,XG,XN,XP,ERROR,*9999)
          DO noelem=1,NEELEM(0,vreg2)
            NELIST(noelem)=NEELEM(noelem,vreg2)
          ENDDO
          NELIST(0)=NEELEM(0,vreg2)
          CALL VOLUME(NBJ,NELIST,NKJE,NPF,NPNE,vreg2,NVJE,NW,
     '      PG,RG,SE,VOL2,WG,XA,XE,XG,XN,XP,ERROR,*9999)
          CALL CHK_INTERCEPT(NPNODE,vreg1,vreg2,CROSS_OVER_PENALTY,XP,
     '      ERROR,*9999)
          TOTAL=(VOL1-VOL2)
          TOTAL=100*(TOTAL/VOL1)
          RESID(1)=PERCENT_FAT-TOTAL
          RESID(2)=CROSS_OVER_PENALTY
        ELSEIF(CUSTOMISATION_TYPE.EQ.4) THEN
          DO nores=1,INT(NT_RES*0.5)
            ZVAL=WDMEASURE(3,nores)
            CALL MESHXY(IBT,IDO,INP,NBJ,1,NEELEM,NPNE,DEPTH,SE,XP,
     '        WIDTH1,ZVAL,ERROR,*9999)
            RESID(2*nores-1)=WIDTH1-WDMEASURE(1,nores)
            RESID(2*nores)=DEPTH-WDMEASURE(2,nores)
          ENDDO !nores
        ENDIF

      ELSE IF(PARAMTYPE(1:7).EQ.'MOMENTS') THEN

        CALL MOMENTS(NRLIST,SMOM,XP,ZD,ERROR,*9999)
        DO nores=1,NJT+2*NJT-3
          RESID(nores)=DABS(SMOM(1,nores)-SMOM(2,nores))
        ENDDO

      ELSE IF(PARAMTYPE(1:9).EQ.'POTENTIAL') THEN

        ERR=0

        nts=SOPTI_END-SOPTI_START+1
        IF(ISIZE_PHI(2).LT.nts) THEN
          WRITE(OP_STRING(1),'(''>>PHI has temporal size'',I5)')
     '      ISIZE_PHI(2)
          WRITE(OP_STRING(2),'(''>> Optimising '',I5,'' time steps'')')
     '      nts
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ERROR='>>Temporal size of PHI is less than the optim range'
          GOTO 9999
        ENDIF

        IF(MODE.EQ.1.OR.MODE.EQ.2) THEN !Calculate Jacobian

C!!! LKC 29-NOV-1999 These asserts (slow) should only be done the
C!!!  first time. maybe before the first call to resfun?

          CALL ASSERT(SOPTI_END-SOPTI_START.LE.NY_TRANSFER_M,
     '      '>>WK1_INV is to small for temporal info.'
     '      //' Increase NY_TRANSFER_M',
     '      ERROR,*9999)
          CALL ASSERT(NT_RES.LE.NY_TRANSFER_M,
     '      '>>WK1_INV is too small for resids. Increase NY_TRANSFER_M',
     '      ERROR,*9999)

          IF(ACTN_IREGULARISE.LE.1) THEN
            CALL ASSERT(NT_RES.EQ.ISIZE_TBH(1),
     '        '>>Number of resids <> # Torso nodes',ERROR,*9999)
          ELSEIF(ACTN_IREGULARISE.EQ.2) THEN
            CALL ASSERT(NT_RES-1.EQ.ISIZE_TBH(1),
     '        '>>Number of resids-1 <> # Torso nodes',ERROR,*9999)
          ELSE
            ERROR='>>Update Code: Unknown Regularisation'
            GOTO 9999
          ENDIF


          DO nores=1,NT_RES
            DO nts=SOPTI_START,SOPTI_END
              WK1_INV(nores,nts)=0.0d0
            ENDDO !nts
          ENDDO !nores

        ENDIF !mode=1 or 2

        DO nores=1,NT_RES-ACTN_IREGULARISE+1
!         these should correspond to the surface nodes
!         PHI(nytr,nts) is the measured signal at torso dof
!                       nytr (=nores if linear elements) and time nts
!         PHI_H(nytr,nts) is the estimated signal at torso dof
!                       nytr (=nores if linear elements) and time nts

!         Need to calculate squared difference and integrate over time
!         (or just use sum of squares of differences over time for now)
!         i.e. resid function = sum_over_time{phi-phi_h}^2

!         Need to calculate the correlation coefficient between PHI and
!         PHI_H over time.

          PHI_MEAN=0.0d0
          PHI_H_MEAN=0.0d0
          count=0
          PHI_NORM=0.0d0
          PHI_H_NORM=0.0d0
          IF(DABS(CC_OBJ_WEIGHT).GT.ZERO_TOL) THEN

            DO nts=SOPTI_START,SOPTI_END
              count=count+1
              PHI_MEAN=PHI_MEAN+PHI(nores,nts)
              PHI_H_MEAN=PHI_H_MEAN+PHI_H(nores,nts)
            ENDDO
            IF(COUNT.GT.0) THEN
              PHI_MEAN=PHI_MEAN/DBLE(count)
              PHI_H_MEAN=PHI_H_MEAN/DBLE(count)
            ENDIF

            DO nts=SOPTI_START,SOPTI_END
              PHI_NORM=PHI_NORM+(PHI(nores,nts)-PHI_MEAN)**2
              PHI_H_NORM=PHI_H_NORM+(PHI_H(nores,nts)-PHI_H_MEAN)**2
            ENDDO
            PHI_NORM=DSQRT(PHI_NORM)
            PHI_H_NORM=DSQRT(PHI_H_NORM)
          ENDIF !CC calculations


C*** Evaluate residuals (sum of squares)
          SS_RESID=0.0d0
          CC_RESID=0.0d0
          DO nts=SOPTI_START,SOPTI_END
            IF(DABS(SS_OBJ_WEIGHT).GT.ZERO_TOL) THEN
              SS_RESID=SS_RESID+(PHI(nores,nts)-PHI_H(nores,nts))**2
            ENDIF
            IF(DABS(CC_OBJ_WEIGHT).GT.ZERO_TOL) THEN
              CC_RESID=CC_RESID+(PHI(nores,nts)-PHI_MEAN)*
     '          (PHI_H(nores,nts)-PHI_H_MEAN)
            ENDIF
          ENDDO !nts
          IF(DABS(CC_OBJ_WEIGHT).GT.ZERO_TOL) THEN
            IF(DABS(PHI_NORM*PHI_H_NORM).GT.ZERO_TOL)
     '        CC_RESID=CC_RESID/(PHI_NORM*PHI_H_NORM)
            CC_RESID=DSQRT(1.0d0-CC_RESID)
          ENDIF
          RESID(nores)=SS_OBJ_WEIGHT*SS_RESID+CC_OBJ_WEIGHT*CC_RESID
          TOT_RES=TOT_RES+RESID(nores)

          IF(MODE.EQ.1.OR.MODE.EQ.2) THEN !Calculate Jacobian

            PHI_PHI_H_SUM=0.0d0
            IF(DABS(CC_OBJ_WEIGHT).GT.ZERO_TOL) THEN
              DO nts=SOPTI_START,SOPTI_END
                PHI_PHI_H_SUM=PHI_PHI_H_SUM+(PHI(nores,nts)-PHI_MEAN)*
     '            (PHI_H(nores,nts)-PHI_H_MEAN)
              ENDDO !nts
            ENDIF

            DO noopti=1,NTOPTI

              IF(PAOPTY(noopti).EQ.1) THEN
                DERIVATIVE_TYPE=2
                HEART_ny=NYNO(1,noopti,2,nr)
                nk=NPNY(1,HEART_ny,0)
                nv=NPNY(2,HEART_ny,0)
                nh=NPNY(3,HEART_ny,0)
                np=NPNY(4,HEART_ny,0)
                FOUND=.FALSE.
                i=0
                nolist=1
                DO WHILE(nolist.LE.NPLIST3(0).AND..NOT.FOUND)
                  npp=NPLIST3(nolist)
                  IF(np.EQ.npp) THEN
                    FOUND=.TRUE.
                    i=i+1
                  ELSE
                    DO nhx=1,NHP(npp,TRSF_NR_FIRST)
                      nh=NH_LOC(nhx,nx)
                      DO nv=1,NVHP(nh,npp,1,TRSF_NR_FIRST)
                        DO nk=1,MAX(NKH(nh,npp,1,TRSF_NR_FIRST)-
     '                    KTYP93(1,TRSF_NR_FIRST),1)
                          i=i+1 !Appropriate column of T_BH
                        ENDDO !nk
                      ENDDO !nv
                    ENDDO !nhx
                  ENDIF
                  nolist=nolist+1
                ENDDO !nolist
                IF(.NOT.FOUND) THEN
                  ERROR='>>ERROR: Could not find np in list of heart '
     '              //'nodes'
                  GOTO 9999
                ENDIF
              ELSE IF(PAOPTY(noopti).EQ.2) THEN
                DERIVATIVE_TYPE=3
                DERIV_PARAM_NUM=2
                i=0
                HEART_ny=0
              ELSE IF(PAOPTY(noopti).EQ.3) THEN
                DERIVATIVE_TYPE=3
                DERIV_PARAM_NUM=3
                i=0
                HEART_ny=0
              ELSE
                ERROR='>>Invalid PAOPTY for activation optimisation'
                GOTO 9999
              ENDIF

              D_PHI_H_MEAN=0.0d0

              DO nts=SOPTI_START,SOPTI_END
                TIME=CALC_TIME_FROM_SAMPLE(nts)
                WK1_INV(1,nts)=CALC_FORWARD_ACTN(DERIVATIVE_TYPE,
     '            DERIV_PARAM_NUM,HEART_ny,i,nores,NHP,NKH,
     '            NPLIST3,NVHP,nx,NYNP,0.5d0/TRSF_FREQUENCY,T_BH,
     '            TIME,YP,TRSF_ACTN_WAVE_INTERPOLATE,ERR,ERROR)
                IF(ERR.NE.0) GOTO 9999
                D_PHI_H_MEAN=D_PHI_H_MEAN+WK1_INV(1,nts)
              ENDDO !nts
              IF(count.GT.0)
     '          D_PHI_H_MEAN=D_PHI_H_MEAN/DBLE(count)
              SS_SUM=0.0d0
              CC_SUM=0.0d0
              DO nts=SOPTI_START,SOPTI_END
                IF(DABS(SS_OBJ_WEIGHT).GT.ZERO_TOL) THEN
                  SS_SUM=SS_SUM-2.0d0*(PHI(nores,nts)-
     '              PHI_H(nores,nts))*WK1_INV(1,nts)
                ENDIF
                IF(DABS(CC_OBJ_WEIGHT).GT.ZERO_TOL) THEN
                  CC_SUM=CC_SUM+(PHI_H_NORM**2*(PHI(nores,nts)-
     '              PHI_MEAN)-PHI_PHI_H_SUM*(PHI_H(nores,nts)-
     '              PHI_H_MEAN))*(WK1_INV(1,nts)-D_PHI_H_MEAN)
                ENDIF
              ENDDO
              IF(DABS(CC_OBJ_WEIGHT).GT.ZERO_TOL) THEN
                CC_SUM=-CC_SUM/(2.0d0*CC_RESID)
                IF(DABS(PHI_NORM*PHI_H_NORM).GT.ZERO_TOL) THEN
                  CC_SUM=CC_SUM/(PHI_NORM*PHI_H_NORM**3)
                ENDIF
              ENDIF
              RESJAC(nores,noopti)=SS_OBJ_WEIGHT*SS_SUM+
     '          CC_OBJ_WEIGHT*CC_SUM

            ENDDO !noopti

          ENDIF !mode=1 or 2 (calculate Jacobian)

        ENDDO !nores

C*** Add an additional constraint as nores=NT_RES
        IF(ACTN_IREGULARISE.EQ.1) THEN
          IF(LIST_RESID.GE.2) THEN
            WRITE(OP_STRING(1),'('' '')')
            WRITE(OP_STRING(2),'('' *** No additional constraints'')')
            WRITE(OP_STRING(3),'('' '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF !output

        ELSE IF(ACTN_IREGULARISE.EQ.2) THEN

C LKC 15-JUL-2002 New check
          CALL ASSERT(CALL_CALC_LAPL,' >> Evaluate Laplacian first',
     '      ERROR,*9999)

          IF(LIST_RESID.GE.2) THEN
            WRITE(OP_STRING(1),'('' '')')
            WRITE(OP_STRING(2),'('' *** Doing surface laplacian'')')
            WRITE(OP_STRING(3),'('' '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF !output

C*** LKC start surface lapl
C          TOT_LAPL=0.0d0
C          CALL SURFACELAPLACIAN(IBT,LIST_RESID,
C     '      NBH,NENP,NPLIST3,NPNE,nx,NXI,NYNP,
C     '      RESID,RESJAC,TOT_LAPL,XP,YP,
C     '      ERROR,*9999)


C*** Surfacelaplacian == L x tau
C*** Derivative == L^T L x tau

          nr=TRSF_NR_FIRST !heart surface
          TOT_LAPL=0.D0
          DO no_nynr=1,NYNR(0,0,1,nr) !same as heart nodes
            SUM2=0.D0
            SUM3=0.D0
            DO no_nynr1=1,NYNR(0,0,1,nr)
              ny1=NYNR(no_nynr1,0,1,nr)
              SUM2=SUM2+LAPL(no_nynr,no_nynr1)*YP(ny1,1)
              SUM3=SUM3+LAPLSQR(no_nynr,no_nynr1)*YP(ny1,1)
            ENDDO
            TOT_LAPL=TOT_LAPL+SUM2*SUM2
            RESID(NT_RES)=TOT_LAPL*ACTN_REG_PARAM_LAPLACE
            RESJAC(NT_RES,no_nynr)=SUM3*ACTN_REG_PARAM_LAPLACE*2
          ENDDO



        ENDIF !ACTN_IREGULARISE

C*** Output the residual values (have already been regularised)
        IF(LIST_RESID.GT.2) THEN
          WRITE(OP_STRING,'(/'' YP(noopti)'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(5D14.6)') (YP(NYNO(1,noopti,2,nr),1),
     '      noopti=1,NPLIST3(0))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          WRITE(OP_STRING(1),'('' '')')
          WRITE(OP_STRING(2),'('' RESID(nts)'')')
          WRITE(OP_STRING(3),'('' '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(5D14.6)') (RESID(nts),nts=1,NT_RES)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          DO nts=1,NT_RES
            WRITE(OP_STRING(1),'('' '')')
            WRITE(OP_STRING(2),'('' RESJAC('',I4,'',opt)'')') nts
            WRITE(OP_STRING(3),'('' '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(5D14.6)') (RESJAC(nts,noopti),
     '        noopti=1,NPLIST3(0))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          ENDDO !nts
        ENDIF ! IREGULARISE


C*** Output some statistics

        IF(LIST_RESID.GE.1) THEN
          WRITE(OP_STRING(1),'('' '')')
          WRITE(OP_STRING(2),'('' Activation Residual Summary'')')
          WRITE(OP_STRING(3),'('' '')')

C*** Ouput Phi Residuals

          WRITE(OP_STRING(4),
     '      '('' Total Phi Residual              : '',F11.4)') TOT_RES
          WRITE(OP_STRING(5),
     '      '('' Avg.  Phi Residual              : '',F11.4)')
     '      TOT_RES/DBLE(NT_RES)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          IF(ACTN_IREGULARISE.EQ.2) THEN

C*** Output Laplacian Residuals
            WRITE(OP_STRING(1),'('' '')')
            WRITE(OP_STRING(2),
     '        '('' Total Laplacian Residual        : '',F11.4)')
     '        TOT_LAPL
            WRITE(OP_STRING(3),
     '        '('' Avg. Laplacian Residual         : '',F11.4)')
     '        TOT_LAPL/DBLE(NTOPTI)
            WRITE(OP_STRING(4),
     '        '('' Reg. Avg. Laplacian Residual    : '',F11.4)')
     '        TOT_LAPL/DBLE(NTOPTI)*ACTN_REG_PARAM_LAPLACE
            WRITE(OP_STRING(5),
     '        '('' Reg. Tot. Laplacian Residual    : '',F11.4)')
     '        TOT_LAPL*ACTN_REG_PARAM_LAPLACE
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C*** Output ratios
C*** Output ratios
            WRITE(OP_STRING(1),'('' '')')
            IF(DABS(TOT_LAPL).GT.ZERO_TOL) THEN
              WRITE(OP_STRING(2),
     '          '('' Avg. Phi-Laplacian Ratio        : '',E11.4)')
     '          TOT_RES/DBLE(NT_RES)/TOT_LAPL
            ELSE
              WRITE(OP_STRING(2),
     '          '('' Avg. Phi-Laplacian Ratio        : No Laplacian'')')
            ENDIF

            IF(DABS(ACTN_REG_PARAM_LAPLACE*TOT_LAPL).GT.ZERO_TOL) THEN
              WRITE(OP_STRING(3),
     '          '('' Reg. Avg. Phi-Laplacian Ratio   : '',E11.4)')
     '          TOT_RES/DBLE(NT_RES)/
     '          TOT_LAPL*ACTN_REG_PARAM_LAPLACE
            ELSE
              WRITE(OP_STRING(3),
     '          '('' Regularised Phi-Laplacian Ratio : No Laplacian'')')
            ENDIF
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          ENDIF !ACTN_REGULARISE
        ENDIF !LIST_RESID

      ELSE IF((PARAMTYPE(1:6).EQ.'STRESS').OR.
     '  (PARAMTYPE(1:6).EQ.'STRAIN')) THEN
       CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '    nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)

        resid_index=1
        DO nolist=1,NELIST(0)
          ne=NELIST(nolist)
          nb=NBH(NH_LOC(1,nx),1,ne)
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
          CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),
     '      NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '      CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '      ERROR,*9999)
          CALL CPCG(1,NBH(1,1,ne),NPNE(1,1,ne),nr,nx,
     '      CE(1,ne),CG,CGE(1,1,ne),CP,PG,ERROR,*9999)

          CONNECTED=.FALSE.
          IF(KTYP29C.EQ.1) THEN        ! resids closest
            closest=RMAX
            NGLIST(0)=0
            DO nn=1,NNT(nb)
              np=NPNE(nn,nb,ne)
              DO no=1,NTOPTI
                ny=NYNO(1,no,1,nr)
                np2=NPNY(4,ny,0)
                IF(np2.EQ.np) THEN !find resids for no
                  CONNECTED=.TRUE.
                  NGLIST(0)=NGLIST(0)+1
                  DO ng=1,NGT(nb) ! find closest ng to to np
                    CALL CALC_NP_XI(IBT,INP,NBJ,ne,NENP,np,NPNE,
     '                NRLIST,XI,ERROR,*9999)
                    dist=DSQRT((XI(1)-XIG(1,ng,nb))**2+
     '                (XI(2)-XIG(2,ng,nb))**2+
     '                (XI(3)-XIG(3,ng,nb))**2)
                    IF(dist.LE.closest) THEN
                      closest=dist
                      NGLIST(NGLIST(0))=ng
                    ENDIF ! dist
                  ENDDO ! ng
                ENDIF ! .EQ.np
             ENDDO ! no
            ENDDO ! nn
          ELSE IF(KTYP29C.EQ.3) THEN   ! resids all connected weighted
            DO ig=1,NGLIST(0)
              WEIGHT(ig)=0.0d0
            ENDDO ! ig
            nu=1
            nk=1
            DO nn=1,NNT(nb)
              np=NPNE(nn,nb,ne)
              DO no=1,NTOPTI
                ny=NYNO(1,no,1,nr)
                np2=NPNY(4,ny,0)
                IF(np2.EQ.np) THEN !eval resid weights
                 CONNECTED=.TRUE.
                  DO ig=1,NGLIST(0)
                    ng=NGLIST(ig)
                    WEIGHT(ig)=WEIGHT(ig)+
     '                PSI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                INP(1,1,nb),nb,nn,nk,nu,XIG(1,ng,nb))
                  ENDDO ! ig
                ENDIF ! .EQ.np
              ENDDO ! no
            ENDDO ! nn
          ENDIF

          IF((KTYP27.LT.4).OR.(KTYP27.EQ.5)) THEN
            DO ig=1,NGLIST(0)
              ng=NGLIST(ig)
              CALL ZETX50(COORDS,'Cauchy',STRESSTYPE,
     '          IBT,IDO,INP,NAN,NBH(1,1,ne),NBJ(1,ne),ng,NHE(ne),
     '          NPNE(1,1,ne),nr,ne,nx,
     '          CE(1,ne),CG,CP,FEXT(1,ng,ne),PG,PHI,PST,RG(ng),RG2D,
     '          RGZ,RGZ2D,RM,TC,TG,TN,TNA,XE,XG,XIG(1,ng,nb),
     '          YG(1,ng,ne),ZE,ZG,ERROR,*9999)
              IF(CONNECTED) THEN
                IF(KTYP27.EQ.1) THEN
                  IF((KTYP29C.EQ.1).OR.(KTYP29C.EQ.2)) THEN ! resids not weighted
                    RESID(resid_index)=TG(1,1)
                    resid_index=resid_index+1
                  ELSE IF(KTYP29C.EQ.3) THEN   ! resids all connected weighted
                    RESID(resid_index)=TG(1,1)*WEIGHT(ig)
                    resid_index=resid_index+1
                  ENDIF
                ELSE IF (KTYP27.EQ.5) THEN
                  IF((KTYP29C.EQ.1).OR.(KTYP29C.EQ.2)) THEN ! resids not weighted
                    RESID(resid_index)=TG(2,2)
                    resid_index=resid_index+1
                  ELSE IF(KTYP29C.EQ.3) THEN   ! resids all connected weighted
                    RESID(resid_index)=TG(2,2)*WEIGHT(ig)
                    resid_index=resid_index+1
                  ENDIF
                ELSE
                  TGNG(ng)=TG(1,1)
                ENDIF
              ENDIF
            ENDDO !ig
          ELSE
            DO ig=1,NGLIST(0)
              ng=NGLIST(ig)
              CALL ZEEX50(COORDS,IBT,IDO,INP,NAN,NBH(1,1,ne),
     '          NBJ(1,ne),ng,NHE(ne),NPNE(1,1,ne),nr,nx,
     '          DZDX,CE(1,ne),CG,CP,EG,PG,PHI,PST,
     '          R,RG(ng),RI1,RI2,RI3,RM,U,
     '          XE,XG,XIG(1,ng,nb),ZE,ZG,ERROR,*9999)
              IF(CONNECTED) THEN
                IF(KTYP27.EQ.4) THEN
                  IF((KTYP29C.EQ.1).OR.(KTYP29C.EQ.2)) THEN ! resids not weighted
                    RESID(resid_index)=EG(2,3)
                    resid_index=resid_index+1
                  ELSE IF(KTYP29C.EQ.3) THEN   ! resids all connected weighted
                    RESID(resid_index)=EG(2,3)*WEIGHT(ig)+100
                    resid_index=resid_index+1
                  ENDIF
                ELSE
                  TGNG(ng)=EG(2,3)
                ENDIF
              ENDIF
            ENDDO !ig
          ENDIF


          IF((KTYP27.EQ.2).OR.(KTYP27.EQ.3)) THEN
            DO ig=1,NGLIST(0)
              ng=NGLIST(ig)
              IF(CONNECTED) THEN
                IF(ng.LE.NGAP(1,nb)*NGAP(2,nb)) THEN
                  NEIGHBOURS=1
                  NEIGHBOUR1=NGAP(1,nb)*NGAP(2,nb)+ng
                ELSE IF(ng.GT.(NGT(nb)-NGAP(1,nb)*NGAP(2,nb))) THEN
                  NEIGHBOURS=1
                  NEIGHBOUR1=NGT(nb)-2*(NGAP(1,nb)*NGAP(2,nb))
     '              +(ng-2*(NGAP(1,nb)*NGAP(2,nb)))
                ELSE
                  NEIGHBOURS=2
                  NEIGHBOUR1=ng-NGAP(1,nb)*NGAP(2,nb)
                  NEIGHBOUR2=ng+NGAP(1,nb)*NGAP(2,nb)
                ENDIF
                IF(NEIGHBOURS.EQ.1) THEN
                  DIFF=DABS(TGNG(ng)-TGNG(NEIGHBOUR1))
                ELSE
                  DIFF=0.5d0*DABS(TGNG(ng)-TGNG(NEIGHBOUR1))
     '                +0.5d0*DABS(TGNG(ng)-TGNG(NEIGHBOUR2))
                ENDIF
                IF(KTYP27.EQ.2) THEN
                  RESID(resid_index)=TGNG(ng)*WEIGHT(ig)+DIFF*WEIGHT(ig)
                ELSE IF (KTYP27.EQ.3) THEN
                  RESID(resid_index)=DIFF*WEIGHT(ig)
                ENDIF
                resid_index=resid_index+1
              ENDIF
            ENDDO !ig
          ENDIF
        ENDDO !nolist (ne)

      ELSE IF((PARAMTYPE(1:8).EQ.'PRESSURE')) THEN
        DO resid_index=1,NT_RES
          RESID(resid_index)=DABS(YP(resid_index,7))
     '      -DABS(YP(resid_index,8))
        ENDDO
        !residual calculated from vectors loaded into YP 
      
C news HS 17/4/04 changed index of YP from 8 to 9 since it interfered with the line search option in LSFUNC.f
C see also UPRESI.f line 143, had to be adjusted for residual calculation
      ELSE IF((PARAMTYPE(1:8).EQ.'REACTION')) THEN
        DO resid_index=1,NT_RES
          RESID(resid_index)=YP(resid_index,7)
     '      -YP(resid_index,9)
        ENDDO
C newe
      ENDIF


      CALL EXITS('RESFUN')
      RETURN
 9999 CALL ERRORS('RESFUN',ERROR)
      CALL EXITS('RESFUN')
      RETURN 1
      END


