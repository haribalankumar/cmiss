      SUBROUTINE D_RESFUN(PARAMTYPE,IBT,IDO,INP,LGE,NAN,NBH,NBJ,NBJF,
     '  NEELEM,NFF,NGAP,NHE,NHP,NKEF,NKH,NKHE,NKJE,NMNO,NNF,NPF,
     '  NPNE,NPNY,NPNODE,nr,NRE,NVHE,NVHP,NVJE,
     '  NW,nx,NXI,NYNE,NYNP,NYNR,CE,CG,CGE,CP,CURVCORRECT,D_RE,D_RI3,
     '  D_RP,D_TG,D_ZG,ES,FEXT,FIX,PG,RE1,RE2,RG,SE,WG,
     '  XA,XE,XG,XP,YG,YP,ZA,ZE,ZE1,ZG,ZG1,ZP,ERROR,*)

C#### Subroutine: D_RESFUN
C###  Description:
C###    D_RESFUN returns analytic derivatives of the residuals wrt the
C###    parameters in the optimising material parameter list.


C cpb 28/3/96 This comment needs updating wrt the new YP iy locations
C**** For KTYP27=2 problems there are KTYP28 exptl measurements
C****   in fit and first set of coords & reactions are held in
C****   YP(ny,1) & YP(ny,5) and YP(ny,6/8/..) & YP(ny,7/9/..)
C****   hold remaining coords & reactions, respec.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LGE(NHM*NSM,NRCM),NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NBJF(NJM,NFM),NEELEM(0:NE_R_M,0:NRM),NFF(6,NEM),
     '  NGAP(NIM,NBM),
     '  NHE(NEM),NHP(NPM),NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NMNO(1:2,0:NOPM),
     '  NNF(0:17,6,NBFM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNY(0:6,NYM,0:NRCM),NPNODE(0:NP_R_M,0:NRM),nr,NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNR(0:NY_R_M,0:NRCM,NCM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),CP(NMM,NPM),
     '  CURVCORRECT(2,2,NNM,NEM),D_RE(NSM,NHM,NOPM),
     '  D_RI3(NHM*NSM),D_RP(NYM,NYM),D_TG(3,3,NHM*NSM),
     '  D_ZG(NHM,NUM,NHM*NSM),ES(NHM*NSM,NHM*NSM),FEXT(NIFEXTM,NGM,NEM),
     '  PG(NSM,NUM,NGM,NBM),
     '  RE1(NSM,NHM),RE2(NSM,NHM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),YP(NYM,NIYM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZE1(NSM,NHM),ZG(NHM,NUM),
     '  ZG1(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER PARAMTYPE*(*),ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER iy,k

      CALL ENTERS('D_RESFUN',*9999)

      IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
        IF(KTYP27.EQ.2) THEN !Obj func is sum of squared residuals
          IF(KTYP28.EQ.0) THEN      !residuals calc.d during fit
            ERROR=' >>Not implemented'
            GO TO 9999
          ELSE IF(KTYP28.GT.0) THEN !use existing residuals

            DO k=1,KTYP28
              IF(k.EQ.1) THEN
                IY=1
              ELSE
                IY=4+2*k
              ENDIF
              CALL YPZP(IY,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '          nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
              CALL D_ZPRP(PARAMTYPE,IBT,IDO,INP,LGE,
     '          NAN,NBH,NBJ,NBJF,NEELEM,NFF,NGAP,NHE,
     '          NKEF,NKHE,NKJE,NMNO,NNF,NPF,NPNE,NPNY,
     '          nr,NRE,NVHE,NVJE,NW,nx,NXI,NYNE,NYNP,NYNR,
     '          CE,CG,CGE,CP,CURVCORRECT,D_RE,D_RI3,D_RP,D_TG,D_ZG,
     '          ES,FEXT,FIX,PG,RE1,RE2,RG,SE,WG,XA,XE,XG,XP,YG,
     '          ZA,ZE,ZE1,ZG,ZG1,ZP,ERROR,*9999)
            ENDDO
          ENDIF
        ENDIF

      ELSE IF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
        CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '    nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
        CALL D_ZPRP(PARAMTYPE,IBT,IDO,INP,LGE,
     '    NAN,NBH,NBJ,NBJF,NEELEM,NFF,NGAP,NHE,
     '    NKEF,NKHE,NKJE,NMNO,NNF,NPF,NPNE,NPNY,
     '    nr,NRE,NVHE,NVJE,NW,nx,NXI,NYNE,NYNP,NYNR,
     '    CE,CG,CGE,CP,CURVCORRECT,D_RE,
     '    D_RI3,D_RP,D_TG,D_ZG,ES,FEXT,FIX,
     '    PG,RE1,RE2,RG,SE,WG,XA,XE,XG,XP,YG,
     '    ZA,ZE,ZE1,ZG,ZG1,ZP,ERROR,*9999)
      ENDIF

      CALL EXITS('D_RESFUN')
      RETURN
 9999 CALL ERRORS('D_RESFUN',ERROR)
      CALL EXITS('D_RESFUN')
      RETURN 1
      END


