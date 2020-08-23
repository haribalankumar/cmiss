      SUBROUTINE MNBRAK(CONY,IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '  NEELEM,NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,NNF,
     '  NONY,NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP,NVJE,NW,nx,NXI,
     '  NYNE,NYNO,NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,CP,
     '  CURVCORRECT,FEXT,FIX,PG,RE,
     '  RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,ZE,ZE1,
     '  ZP,ZP1,AX,BX,
     '  CX,FA,FB,FC,F1,FTOL,LSSTAT,ERROR,*)

C#### Subroutine: MNBRAK
C###  Description:
C###    MNBRAK brackets a minimum in 1D. From Numerical Recipes.

C**** AX,FA represents the largest function value, and BX,FB represents
C**** the smallest. If FB is below our tolerance for relaxed line
C**** search, LSSTAT is returned 1. If the minimum is bracketed, LSSTAT
C**** is 0.  If a bracket cannot be found in a reasonable number of
C**** iterations, LSSTAT is 2, the first function value is returned in
C**** F1. This is generally for x=1.0

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LGE(NHM*NSM,NRCM),LSSTAT,NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NFF(6,NEM),NFFACE(0:NF_R_M,NRM),
     '  NGAP(NIM,NBM),NHE(NEM),NHP(NPM),
     '  NKB(2,2,2,NNM,NBFM),NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),
     '  nr,NRE(NEM),NSB(NKM,NNM,NBFM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),
     '  nx,NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNR(0:NY_R_M,0:NRCM,NCM),
     '  NYNO(0:NYOM,NOOPM,NRCM),NONY(0:NOYM,NYM,NRCM),
     '  Z_CONT_LIST(NDM,2,7)
      REAL*8 AX,BX,CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),
     '  CP(NMM,NPM),
     '  CURVCORRECT(2,2,NNM,NEM),CX,F1,FA,FB,FC,
     '  FEXT(NIFEXTM,NGM,NEM),FTOL,
     '  PG(NSM,NUM,NGM,NBM),RE(NSM,NHM),
     '  RG(NGM),SE(NSM,NBFM,NEM),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM),
     '  ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),
     '  Z_CONT(NDM,2,67),ZE(NSM,NHM),ZE1(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM),
     '  CONY(0:NOYM,NYM,NRCM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER iter
      REAL*8 DUM,FU,GLIMIT,GOLD,Q,R,TINY,U,ULIM
      PARAMETER (GLIMIT=20.0D0, GOLD=1.618034D0, TINY=1.0D-20)

      CALL ENTERS('MNBRAK',*9999)

      LSSTAT=0
      CALL LSFUNC(CONY,IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '  NEELEM,
     '  NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,NNF,
     '  NONY,NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP,NVJE,NW,
     '  nx,NXI,
     '  NYNE,NYNO,NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,CP,
     '  CURVCORRECT,FEXT,FIX,
     '  PG,RE,RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,
     '  ZE,ZE1,ZP,ZP1,
     '  BX,FB,ERROR,*9999)
      F1=FB
      IF(FB .LT. FTOL) THEN
        LSSTAT=1
        GOTO 9998
      ENDIF

      IF(FB.GT.FA) THEN
        DUM=AX
        AX=BX
        BX=DUM
        DUM=FB
        FB=FA
        FA=DUM
      ENDIF

      CX=BX+GOLD*(BX-AX)
      CALL LSFUNC(CONY,IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '  NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,NNF,NONY,
     '  NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP,NVJE,NW,nx,NXI,
     '  NYNE,NYNO,NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,CP,
     '  CURVCORRECT,FEXT,FIX,
     '  PG,RE,RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,
     '  ZE,ZE1,ZP,ZP1,
     '  CX,FC,ERROR,*9999)
      IF(FC .LT. FTOL) THEN
        BX=CX
        FB=FC
        LSSTAT=1
        GOTO 9998
      ENDIF
      DO 10 iter=1,8
        IF(DOP) THEN
          WRITE(OP_STRING,'(''Iteration '',I1)') iter
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''AX ='',D14.7,'' FA ='',D14.7)') AX,FA
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''BX ='',D14.7,'' FB ='',D14.7)') BX,FB
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''CX ='',D14.7,'' FC ='',D14.7)') CX,FC
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(FB.LT.FC) GOTO 9998
        R=(BX-AX)*(FB-FC)
        Q=(BX-CX)*(FB-FA)
        U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.0D0*DSIGN(DMAX1(DABS(Q-R),TINY),
     '    Q-R))
        ULIM=BX+GLIMIT*(CX-BX)
        IF((BX-U)*(U-CX).GT.0.0D0) THEN
          CALL LSFUNC(CONY,IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '      NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,
     '      NNF,NONY,NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP,
     '      NVJE,NW,nx,
     '      NXI,NYNE,NYNO,NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,
     '      CP,CURVCORRECT,
     '      FEXT,FIX,PG,RE,RG,SE,WG,XA,XG,XP,YG,YGF,YP,
     '      ZA,ZA1,Z_CONT,
     '      ZE,ZE1,ZP,ZP1,U,FU,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,'(''U  ='',D14.7,'' FU ='',D14.7)') U,FU
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          IF(FU .LT. FTOL) THEN
            BX=U
            FB=FU
            LSSTAT=1
            GOTO 9998
          ENDIF
          IF(FU.LT.FC) THEN
            AX=BX
            FA=FB
            BX=U
            FB=FU
            GOTO 9998
          ELSE IF(FU.GT.FB) THEN
            CX=U
            FC=FU
            GOTO 9998
          ENDIF
          U=CX+GOLD*(CX-BX)
          CALL LSFUNC(CONY,IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '      NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,
     '      NNF,NONY,NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP,NVJE,
     '      NW,nx,
     '      NXI,NYNE,NYNO,NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,
     '      CP,CURVCORRECT,
     '      FEXT,FIX,PG,RE,RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,
     '      ZA1,Z_CONT,
     '      ZE,ZE1,ZP,ZP1,U,FU,ERROR,*9999)
          IF(FU .LT. FTOL) THEN
            BX=U
            FB=FU
            LSSTAT=1
            GOTO 9998
          ENDIF
        ELSE IF((CX-U)*(U-ULIM).GT.0.0D0) THEN
          CALL LSFUNC(CONY,IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '      NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,
     '      NNF,NONY,NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP,NVJE,
     '      NW,nx,
     '      NXI,NYNE,NYNO,NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,
     '      CP,CURVCORRECT,
     '      FEXT,FIX,PG,RE,RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,
     '      ZA1,Z_CONT,
     '      ZE,ZE1,ZP,ZP1,U,FU,ERROR,*9999)
          IF(FU .LT. FTOL) THEN
            BX=U
            FB=FU
            LSSTAT=1
            GOTO 9998
          ENDIF
          IF(FU.LT.FC) THEN
            BX=CX
            CX=U
            U=CX+GOLD*(CX-BX)
            FB=FC
            FC=FU
            CALL LSFUNC(CONY,IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '        NEELEM,
     '        NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,
     '        NNF,NONY,NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP,NVJE,
     '        NW,nx,
     '        NXI,NYNE,NYNO,NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,
     '        CP,CURVCORRECT,
     '        FEXT,FIX,PG,RE,RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,
     '        ZA1,Z_CONT,
     '        ZE,ZE1,ZP,ZP1,U,FU,ERROR,*9999)
            IF(FU .LT. FTOL) THEN
              BX=U
              FB=FU
              LSSTAT=1
              GOTO 9998
            ENDIF
          ENDIF
        ELSE IF((U-ULIM)*(ULIM-CX).GE.0.0D0) THEN
          U=ULIM
          CALL LSFUNC(CONY,IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '      NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,
     '      NNF,NONY,NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP,NVJE,NW,
     '      nx,
     '      NXI,NYNE,NYNO,NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,
     '      CP,CURVCORRECT,
     '      FEXT,FIX,PG,RE,RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,
     '      ZA1,Z_CONT,
     '      ZE,ZE1,ZP,ZP1,U,FU,ERROR,*9999)
          IF(FU .LT. FTOL) THEN
            BX=U
            FB=FU
            LSSTAT=1
            GOTO 9998
          ENDIF
        ELSE
          U=CX+GOLD*(CX-BX)
          CALL LSFUNC(CONY,IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '      NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,
     '      NNF,NONY,NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP,NVJE,NW,
     '      nx,
     '      NXI,NYNE,NYNO,NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,
     '      CP,CURVCORRECT,
     '      FEXT,FIX,PG,RE,RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,
     '      ZA1,Z_CONT,
     '      ZE,ZE1,ZP,ZP1,U,FU,ERROR,*9999)
          IF(FU .LT. FTOL) THEN
            BX=U
            FB=FU
            LSSTAT=1
            GOTO 9998
          ENDIF
        ENDIF
        AX=BX
        BX=CX
        CX=U
        FA=FB
        FB=FC
        FC=FU
 10   CONTINUE
      LSSTAT=2

 9998 CALL EXITS('MNBRAK')
      RETURN
 9999 CALL ERRORS('MNBRAK',ERROR)
      CALL EXITS('MNBRAK')
      RETURN 1
      END


