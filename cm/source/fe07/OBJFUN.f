      SUBROUTINE OBJFUN(IBT,IDO,INP,LD,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '  NEELEM,NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,
     '  NLL,NLNO,NNB,NNF,NNL,NONL,NONY,NP_INTERFACE,NP1OPT,
     '  NPF,NPL,NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP,NVJE,NW,
     '  nx,NXI,NYNE,NYNO,NYNP,NYNR,PAOPTY,Z_CONT_LIST,CE,CG,CGE,CP,
     '  CURVCORRECT,DL,FEXT,FGRAD,PAOPTI,PG,RE1,RESIDM,
     '  RESJAC,RG,SE,WG,WU,XA,XE,XG,XID,XIG,XP,YG,YGF,YP,
     '  ZA,ZA1,Z_CONT,ZD,ZE,ZE1,ZG,ZP,ZP1,FIX,ERROR,*)

C#### Subroutine: OBJFUN
C###  Description:
C###    OBJFUN returns objective function for optimising material
C###    parameters etc.

C cpb 28/3/96 This comment needs to be updated for new iy scheme.
C**** For KTYP27=2 problems there are KTYP28 exptl measurements in fit
C****   and first set of coords & reactions are held in YP(ny,4) &
C****   YP(ny,5) and YP(ny,6/8/..) & YP(ny,7/9/..) hold
C****   remaining coords & reactions,respec.
C**** Note: Computed reactions from ZPRP are stored in YP(ny,1).
C**** For KTYP27=3 problems there are NTOPTI integrated fluxes which
C****   will all be zero when the saturated bulbs' surface has been
C****   located.
C****   For this case the solution for both regions is required.
C**** For KTYP26=2 AND KTYP27=5 The objective is for data fitting by
C**** optimisation using MINOS

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbst02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'opti00.cmn'
      INCLUDE 'stab00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(NDM),LGE(NHM*NSM,NRCM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NFF(6,NEM),NFFACE(0:NF_R_M,NRM),
     '  NGAP(NIM,NBM),NHE(NEM),NHP(NPM),
     '  NKB(2,2,2,NNM,NBFM),NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NLNO(NOPM),NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONL(NLM),NONY(0:NOYM,NYM,NRCM),NP1OPT(NOPM),
     '  NP_INTERFACE(0:NPM,0:3),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),
     '  nr,NRE(NEM),NSB(NKM,NNM,NBFM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),
     '  nx,NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM),PAOPTY(NOPM),Z_CONT_LIST(NDM,2,7)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),
     '  CP(NMM,NPM),
     '  CURVCORRECT(2,2,NNM,NEM),DL(3,NLM),
     '  FEXT(NIFEXTM,NGM,NEM),FGRAD(*),PAOPTI(*),
     '  PG(NSM,NUM,NGM,NBM),RE1(NSM,NHM),RESIDM(*),RESJAC(NREM,*),
     '  RG(NGM),SE(NSM,NBFM,NEM),WG(NGM,NBM),
     '  WU(0:NUM+1,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XIG(NIM,NGM,NBM),XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM),
     '  ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),Z_CONT(NDM,2,67),
     '  ZD(NJM,NDM),ZE(NSM,NHM),ZE1(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM)
      LOGICAL FIX(NYM,NIYFIXM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER GETNYR,IY,k,MAP(4,4),nb,nc,nd,
     '  ne,NE1,NE2,nj,njj2,nk,nl,nn,nn1,
     '  noelem,no_interface,nol,noopti,no_nynr,nores,
     '  np,NP1,NP2,NP3,ns,nv,ny,ny1,ny2,nyo
      REAL*8 ENERGY,PH3,PL2S1,PL2S3,PXI,RHO1,RHO2,
     '  SUM1(3),SUM2,SUM3,VELOC1,VELOC2,Z(3),DZ(3)
      LOGICAL CALCJAC,FOUND
      DATA MAP/3,4,0,0,0,0,3,4,1,0,2,0,0,1,0,2/

      CALL ENTERS('OBJFUN',*9999)

      nc=1 !temporary cpb 22/11/94

      IF(KTYP26.EQ.1.AND.KTYP27.EQ.1) THEN !Obj func is max princ stress difference
        CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '    nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)

        PRSTMAX=0.0d0
        PRSTMIN=1.0d6
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '      XA(1,1,ne),XE,XP,ERROR,*9999)
          CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '      CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '      ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
          IF(ITYP1(nr,nx).EQ.4) THEN !linear elasticity
            CALL CPCG(NW(ne,1),NBH(1,1,ne),NPNE(1,1,ne),nr,nx,
     '        CE(1,ne),CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
            nb=NBH(NH_LOC(1,nx),1,ne)
            CALL OPST40(nb,NBH(1,1,ne),NBJ(1,ne),ne,NHE(ne),
     '        NW(ne,1),nx,CG,ENERGY,PG,'BLANKS',WG,
     '        XE,XG,YG(1,1,ne),ZE,ZG,.TRUE.,'         ',ERROR,*9999)
          ELSE IF(ITYP1(nr,nx).EQ.5) THEN !finite elasticity
            nb=NBH(NH_LOC(1,nx),1,ne)
            CALL CPCG(1,nb,NPNE(1,1,ne),nr,nx,
     '        CE(1,ne),CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
            CALL OPC50(IBT,IDO,INP,NAN,NBH(1,1,ne),NBJ(1,ne),ne,
     '        NHE(ne),NPNE(1,1,ne),nr,nx,
     '        CE(1,ne),CG,CP,FEXT(1,1,ne),PG,RG,
     '        XE,XG,XIG(1,1,nb),YG(1,1,ne),ZE,ZG,ERROR,*9999)
          ENDIF
          FUNC=(PRSTMAX-PRSTMIN)**2
        ENDDO
        WRITE(OP_STRING,'('' Max principle stress  = '',D11.3)')
     '    PRSTMAX
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Min principle stress  = '',D11.3)')
     '    PRSTMIN
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Princ. stress diff.^2 = '',D11.3)')
     '    FUNC
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      ELSE IF(KTYP26.EQ.1.AND.KTYP27.EQ.2) THEN !Obj func is sum of squared residuals
        FUNC=0.0d0
        IF(KTYP28.EQ.0) THEN      !residuals calc.d during fit
          CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '      nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
          IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const Vol
C           Put reference state for cavity from YP(ny,10) into
C           ZA1,ZP1 for ZERE55
            CALL YPZP(10,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,
     '        NVHP,nx,NYNE,NYNP,YP,ZA1,ZP1,ERROR,*9999)
          ENDIF

          CALL ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,NFF,
     '      NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,NNF,
     '      NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP,NVJE,NW,nx,NXI,
     '      NYNE,NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,CP,
     '      CURVCORRECT,FEXT,
     '      PG,RE1,RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,
     '      %VAL(0),Z_CONT,ZE,ZE1,
     '      ZP,ZP1,%VAL(0),FIX,ERROR,*9999)

          DO no_nynr=1,NYNR(0,0,1) !loop over global variables
            ny1=NYNR(no_nynr,0,1) !is global variable number
            ny2=GETNYR(2,NPNY,nr,1,0,ny1,NYNE,NYNP)
            !ny2 is global row number
            IF(.NOT.FIX(NYNR(no_nynr,0,2),1)) THEN !not a displ b.c.
              FUNC=FUNC+(YP(ny2,4)-YP(ny1,1))**2
            ENDIF
          ENDDO
        ELSE IF(KTYP28.GT.0) THEN !use existing residuals
          DO k=1,KTYP28
            IF(k.EQ.1) THEN
              IY=1
            ELSE
C!!! This code needs checking for IY locations
              IY=2+2*k
            ENDIF
            CALL YPZP(IY,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '        nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
            IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const Vol
C             Put reference state for cavity from YP(ny,10) into
C             ZA1,ZP1 for ZERE55
              CALL YPZP(10,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,
     '          NVHP,nx,NYNE,NYNP,YP,ZA1,ZP1,ERROR,*9999)
            ENDIF

            CALL ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '        NEELEM,NFF,
     '        NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,NNF,
     '        NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP,NVJE,NW,
     '        nx,NXI,
     '        NYNE,NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,CP,
     '        CURVCORRECT,FEXT,
     '        PG,RE1,RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,
     '        %VAL(0),Z_CONT,ZE,ZE1,
     '        ZP,ZP1,%VAL(0),FIX,ERROR,*9999)

            DO no_nynr=1,NYNR(0,0,1) !loop over global variables
              ny1=NYNR(no_nynr,0,1) !is global variable number
              ny2=GETNYR(2,NPNY,nr,1,0,ny1,NYNE,NYNP)
              !ny2 is global row #
              IF(.NOT.FIX(NYNR(no_nynr,0,2),1)) THEN !not a displ b.c.
                FUNC=FUNC+(YP(ny2,4)-YP(ny1,IY+1))**2
              ENDIF
            ENDDO
          ENDDO
        ENDIF

      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.3) THEN
C ***   Objective function is difference between FE flux and BE flux
C ***   nr refers to the BE region.

        ERROR='>> Not implemented. See backup version'
        GOTO 9999

      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.4) THEN
C ***   Objective function is difference between calc.d height & bdry position
        CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '    nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
        FUNC=0.0D0
        DO noopti=1,NTOPTI
          NP1=NP1OPT(noopti)
          FUNC=FUNC+(ZP(1,1,1,NP1,nc)*40.0D0+XP(1,1,2,NP1))**2
        ENDDO

      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.5) THEN
C ***   Objective function is sum of squares of data errors for fitting

        nc=1
        G1SCALING=.FALSE.
C       Copy in the current estimates of XP and DL from PAOPTI
        DO noopti=1,NTOPTI
          IF(PAOPTY(noopti).EQ.1) THEN !Parameter is geometric dof
            DO nyo=1,NYNO(0,noopti,2)
              ny=NYNO(nyo,noopti,2)
              nk=NPNY(1,ny,0)
              nv=NPNY(2,ny,0)
              nj=NPNY(3,ny,0)
              np=NPNY(4,ny,0)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              XP(nk,nv,nj,np)=PAOPTI(noopti)
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
          ELSE
            ERROR=' >>Invalid PAOPTY for data fitting'
            GOTO 9999
          ENDIF
C Zero gradient vector
          FGRAD(noopti)=0.0d0
        ENDDO
        IF(G1SCALING) THEN
C         Update scale factor arrays
          DO nb=1,NBFT
            CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,DL,
     '        SE,ERROR,*9999)
          ENDDO
        ENDIF
C       Evaluate the objective
        FUNC=0.0d0
        DO nores=1,NDT
          nd=nores
          RESIDM(nores)=0.0d0
          IF(LD(nd).ne.0) THEN
            ne=LD(nd)
C ??? CPB 6/6/93 put xpxe call inside nj loop ?
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj2,nr)
              nb=NBJ(nj,ne)
              Z(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '          XID(1,nd),XE(1,nj))
              DZ(nj)=Z(nj)-ZD(nj,nd)
              RESIDM(nores)=RESIDM(nores)+DZ(nj)**2
            ENDDO
            RESIDM(nores)=DSQRT(RESIDM(nores))
            FUNC=FUNC+RESIDM(nores)
C           Evaluate the objective jacobian
C CPB 6/6/93 WARNING : Assuming same basis function in each geometric
C direction
            nb=NBJ(1,ne)
            IF(NIT(nb).EQ.1) THEN !1D Element
C             Calculate the jacobian wrt to element (nodal) parameters
              DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                nj=NJ_LOC(NJL_GEOM,njj2,nr)
                ns=0
                DO nn=1,NNT(nb)
                  DO nk=1,NKT(nn,nb)
                    ns=ns+1
                    ny=NYNP(nk,1,nj,NPNE(nn,nb,ne),0,nc,nr)
                    noopti=NONY(1,ny,2)
                    IF(noopti.GT.0) THEN !Geometric param. is optimised
                      SUM3=DZ(nj)*PH3(nn,nk,1,XID(1,nd))*SE(ns,nb,ne)
                      FGRAD(noopti)=FGRAD(noopti)+SUM3/RESIDM(nores)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
              IF(G1SCALING) THEN
C               Calculate the jacobian wrt to the line length
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
     '                    PH3(nn,2,1,XID(1,nd))*XP(2,1,nj,NPL(1+nn,1,
     '                    nl))
                      ENDDO
                      SUM1(nj)=SUM1(nj)*DZ(nj)
                      SUM2=SUM2+SUM1(nj)
                    ENDDO
                    FGRAD(noopti)=FGRAD(noopti)+2.0d0*SUM2
                  ENDIF
                ENDIF
              ENDIF
            ELSE IF(NIT(nb).EQ.2) THEN !2D Element
              IF(IBT(1,1,nb).EQ.2.AND.IBT(1,2,nb).EQ.2) THEN !BiHermite
C             Calculate the jacobian wrt to element (nodal) parameters
                DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj2,nr)
                  ns=0
                  DO nn=1,NNT(nb)
                    DO nk=1,NKT(nn,nb)
                      ns=ns+1
                      ny=NYNP(nk,1,nj,NPNE(nn,nb,ne),0,nc,nr)
                      noopti=NONY(1,ny,2)
                      IF(noopti.GT.0) THEN !Geometric param. is optimised
                        SUM3=DZ(nj)*
     '                    PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,XID(1,nd))*
     '                    PH3(INP(nn,2,nb),IDO(nk,nn,2,nb),1,XID(2,nd))*
     '                    SE(ns,nb,ne)
                        FGRAD(noopti)=FGRAD(noopti)+SUM3/RESIDM(nores)
                      ENDIF
                    ENDDO
                  ENDDO
                ENDDO
                IF(G1SCALING) THEN
C               Calculate the jacobian wrt to the line length
                  DO nol=1,4
                    nl=NLL(nol,ne)
                    IF(nl.NE.0) THEN
                      noopti=NONL(nl)
                      IF(noopti.GT.0) THEN !Line is optimised
                        SUM2=0.0d0
                        DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                          nj=NJ_LOC(NJL_GEOM,njj2,nr)
                          SUM1(nj)=0.0D0
                          DO nn1=1,2
                            np=NPL(1+nn1,1,nl)
                            nk=NPL(5,1,nl)
!                 Find local node nn of global node np on element ne
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
     '                        PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,
     '                        XID(1,nd))*
     '                        PH3(INP(nn,2,nb),IDO(nk,nn,2,nb),1,
     '                        XID(2,nd))*
     '                        XP(nk,1,nj,np)
                            nk=4
                            SUM1(nj)=SUM1(nj)+
     '                        PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,
     '                        XID(1,nd))*
     '                        PH3(INP(nn,2,nb),IDO(nk,nn,2,nb),1,
     '                        XID(2,nd))*
     '                        XP(nk,1,nj,np)*DL(1,NLL(MAP(nn,nol),ne))
                          ENDDO
                          SUM1(nj)=SUM1(nj)*DZ(nj)
                          SUM2=SUM2+SUM1(nj)
                        ENDDO
                        FGRAD(noopti)=FGRAD(noopti)+2.0D0*SUM2
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF
              ELSE IF(IBT(1,1,nb).EQ.3.AND.IBT(1,2,nb).EQ.3) THEN !Hermite-Simplex
C               Calculate the jacobian wrt to element (nodal) parameters
                DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj2,nr)
                  ns=0
                  DO nn=1,NNT(nb)
                    DO nk=1,NKT(nn,nb)
                      ns=ns+1
                      ny=NYNP(nk,1,nj,NPNE(nn,nb,ne),0,nc,nr)
                      noopti=NONY(1,ny,2)
                      IF(noopti.GT.0) THEN !Geometric param. is optimised
                        IF(NKT(1,nb).EQ.1) THEN !Apex at node 1
                          IF(ns.EQ.1) THEN
                            SUM3=DZ(nj)*
     '                        PL2S1(INP(nn,2,nb),IDO(nk,nn,2,nb),1,
     '                        XID(2,nd))*SE(ns,nb,ne)
                          ELSE
                            SUM3=DZ(nj)*
     '                        PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,
     '                        XID(1,nd))*
     '                        PL2S1(INP(nn,2,nb),IDO(nk,nn,2,nb),1,
     '                        XID(2,nd))*SE(ns,nb,ne)
                          ENDIF
                        ELSE IF(NKT(3,nb).EQ.1) THEN !Apex at node 3
                          IF(ns.EQ.9) THEN
                            SUM3=DZ(nj)*
     '                        PL2S3(INP(nn,2,nb),IDO(nk,nn,2,nb),1,
     '                        XID(2,nd))*SE(ns,nb,ne)
                          ELSE
                            SUM3=DZ(nj)*
     '                        PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,
     '                        XID(1,nd))*
     '                        PL2S3(INP(nn,2,nb),IDO(nk,nn,2,nb),1,
     '                        XID(2,nd))*SE(ns,nb,ne)
                          ENDIF
                        ELSE
                          ERROR='>>NKT incorrect for Hermite-'
     '                      //'Simplex element'
                          GOTO 9999
                        ENDIF
                        FGRAD(noopti)=FGRAD(noopti)+SUM3/RESIDM(nores)
                      ENDIF
                    ENDDO
                  ENDDO
                ENDDO
                IF(G1SCALING) THEN
C                 Calculate the jacobian wrt to the line length
                  DO nol=1,3
                    nl=NLL(nol,ne)
                    IF(nl.NE.0) THEN
                      noopti=NONL(nl)
                      IF(noopti.GT.0) THEN !Line is optimised
                        SUM2=0.0D0
                        DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                          nj=NJ_LOC(NJL_GEOM,njj2,nr)
                          SUM1(nj)=0.0D0
                          IF(NPL(1,nj,nl).EQ.4) THEN !Cubic Hermite line
                            IF(NKT(1,nb).EQ.1) THEN !Apex at node 1
                              DO nn=2,3
                                np=NPL(1+nn-1,1,nl)
                                nk=NPL(5,1,nl)
                                SUM1(nj)=SUM1(nj)+
     '                            PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,
     '                            XID(1,nd))*PL2S1(INP(nn,2,nb),
     '                            IDO(nk,nn,2,nb),1,XID(2,nd))*
     '                            XP(nk,1,nj,np)
                                nk=4
                                SUM1(nj)=SUM1(nj)+
     '                            PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,
     '                            XID(1,nd))*PL2S1(INP(nn,2,nb),
     '                            IDO(nk,nn,2,nb),1,XID(2,nd))*
     '                            XP(nk,1,nj,np)*DL(1,NLL(nn,ne))
                              ENDDO
                            ELSE IF(NKT(3,nb).EQ.1) THEN !Apex at node 3
                              DO nn=1,2
                                np=NPL(1+nn,1,nl)
                                nk=NPL(4,1,nl)
                                SUM1(nj)=SUM1(nj)+
     '                            PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,
     '                            XID(1,nd))*PL2S3(INP(nn,2,nb),
     '                            IDO(nk,nn,2,nb),1,XID(2,nd))*
     '                            XP(nk,1,nj,np)
                                nk=4
                                SUM1(nj)=SUM1(nj)+
     '                            PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,
     '                            XID(1,nd))*PL2S3(INP(nn,2,nb),
     '                            IDO(nk,nn,2,nb),1,XID(2,nd))*
     '                            XP(nk,1,nj,np)*DL(1,NLL(nn+1,ne))
                              ENDDO
                            ELSE
                              ERROR='>>NKT incorrect for Hermite-'
     '                          //'Simplex element'
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
     '                        PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,
     '                        XID(1,nd))*PL2S1(INP(nn,2,nb),
     '                        IDO(nk,nn,2,nb),1,XID(2,nd))*XP(nk,1,nj,
     '                        np)
                            nk=4
                            SUM1(nj)=SUM1(nj)+
     '                        PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,
     '                        XID(1,nd))*PL2S1(INP(nn,2,nb),
     '                        IDO(nk,nn,2,nb),1,XID(2,nd))*
     '                        XP(nk,1,nj,np)*DL(1,NLL(3,ne))
                          ELSE IF(NPL(1,nj,nl).EQ.7) THEN !Apex at node 3
                            np=NPL(2,1,nl)
                            nk=NPL(4,1,nl)
                            IF(NPNE(1,nb,ne).eq.np) THEN
                              nn=1
                            ELSE
                              nn=2
                            ENDIF
                            SUM1(nj)=SUM1(nj)+
     '                        PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,
     '                        XID(1,nd))*
     '                        PL2S3(INP(nn,2,nb),IDO(nk,nn,2,nb),1,
     '                        XID(2,nd))*XP(nk,1,nj,np)
                            nk=4
                            SUM1(nj)=SUM1(nj)+
     '                        PH3(INP(nn,1,nb),IDO(nk,nn,1,nb),1,
     '                        XID(1,nd))*
     '                        PL2S3(INP(nn,2,nb),IDO(nk,nn,2,nb),1,
     '                        XID(2,nd))*XP(nk,1,nj,np)*DL(1,NLL(1,ne))
                          ELSE
                            ERROR='>> Unknown line type'
                            GOTO 9999
                          ENDIF
                          SUM1(nj)=SUM1(nj)*DZ(nj)
                          SUM2=SUM2+SUM1(nj)
                        ENDDO
                        FGRAD(noopti)=FGRAD(noopti)+2.0D0*SUM2
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF
              ELSE
                ERROR='>>Element type not yet implemented'
                GOTO 9999
              ENDIF
            ELSE
              ERROR='>>Current value of NITB not yet implemented'
              GOTO 9999
            ENDIF
          ENDIF
        ENDDO
        IF(KTYP12.EQ.1) THEN !Sobolev Smoothing
          CALCJAC=.TRUE.
          CALL SOBOLEV(IBT,IDO,INP,NBJ,NEELEM,NKJE,NLL,NONL,
     '      NONY,NPF,NPNE,nr,NRE,NVJE,NYNO,NYNP,DL,RESJAC,
     '      FGRAD,PAOPTI,PG,RG,SE,RESIDM(NT_RES),WG,WU,XA,XE,
     '      XG,XIG,XP,CALCJAC,ERROR,*9999)
          FUNC=FUNC+RESIDM(NT_RES)
        ENDIF

      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.6) THEN
C ***   Objective function is fluid interface dynamic bdry condition
        DO nr=1,2  !!! Why does this go up to 2 only ?
          CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '      nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
        ENDDO
        NE1=NEELEM(1,1) !is first element of region 1
        NE2=NEELEM(1,2) !is first element of region 2
        RHO1=CE(1,NE1)  !is density for region 1
        RHO2=CE(1,NE2)  !is density for region 2
        VELOC1=CE(2,NE1)  !is velocity for region 1
        VELOC2=CE(2,NE2)  !is velocity for region 2
        WRITE(OP_STRING,'('' rho1='',e12.3,'' rho2='',e12.3,'
     '    //''' veloc1='',e12.3,'' veloc2='',e12.3)')
     '    rho1,rho2,veloc1,veloc2
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        FUNC=0.0D0
        DO no_interface=1,NP_INTERFACE(0,3)
          NP1=NP_INTERFACE(no_interface,1)
          NP2=NP_INTERFACE(no_interface,2)
          NP3=NP_INTERFACE(no_interface,3)
          FUNC=FUNC+(RHO1*((ZP(1,1,1,NP1,nc)-ZP(1,1,2,NP1,nc))/TINCR
     '            +VELOC1*ZP(2,1,1,NP1,nc)+G_ACCEL*ZP(1,1,1,NP3,nc))
     '            -RHO2*((ZP(1,1,1,NP2,nc)-ZP(1,1,2,NP2,nc))/TINCR
     '            +VELOC2*ZP(2,1,1,NP2,nc)+G_ACCEL*ZP(1,1,1,NP3,nc)))**2
          WRITE(OP_STRING,'('' interface node '',i2,'
     '      //''' func= '',e12.3)') no_interface,func
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('OBJFUN')
      RETURN
 9999 CALL ERRORS('OBJFUN',ERROR)
      CALL EXITS('OBJFUN')
      RETURN 1
      END


