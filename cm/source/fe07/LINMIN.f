      SUBROUTINE LINMIN(ALPHA,CE,CG,CGE,CONY,CP,CURVCORRECT,
     '  FEXT,FIX,GRR,IBT,IDO,INP,
     '  ITER1,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,NFF,NFFACE,
     '  NGAP,NHE,NHP,NKB,NKH,NKEF,NKHE,NNB,
     '  NNF,NKJE,NONY,NPF,NPNE,NPNODE,
     '  NPNY,NRLIST,nr,NRE,NSB,NVHE,NVHP,NVJE,NW,
     '  nx,NXI,NYNE,NYNO,NYNP,NYNR,PG,RE1,
     '  REF_ENERGY,RG,SE,WG,XA,XG,XP,YP,YG,YGF,ZA,ZA1,Z_CONT,
     '  Z_CONT_LIST,ZE,ZE1,ZP,ZP1,
     '  ERROR,*)

C#### Subroutine: LINMIN
C###  Description:
C###    Find alpha to minimize residual along
C###    search direction YP(...,5).  At the end of this procedure
C###    we will have a new residual vector (YP(...,4)), solution
C###    (YP(...,1)), and solution increment (YP(...,5))
      
      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
! Parameters
      REAL*8 ZE(NSM,NHM),ZE1(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM)
      INTEGER NKJE(NKM,NNM,NJM,NEM),NKH(NHM,NPM,NCM,0:NRM),
     '  NPNE(NNM,NBFM,NEM),NONY(0:NOYM,NYM,NRCM,0:NRM),
     '  NKEF(0:4,16,6,NBFM),NBHF(NHM,NCM,NFM),NPNY(0:6,NYM,0:NRCM),
     '  NKB(2,2,2,NNM,NBFM),NNF(0:17,6,NBFM),
     '  NBH(NHM,NCM,NEM),NBJF(NJM,NFM),
     '  NHE(NEM),NHP(NPM,0:NRM),NGAP(NIM,NBM),
     '  NFF(6,NEM),NFFACE(0:NF_R_M,NRM),NBJ(NJM,NEM),
     '  INP(NNM,NIM,NBFM),ITER1,
     '  NVJE(NNM,NBFM,NJM,NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NNB(4,4,4,NBFM),NSB(NKM,NNM,NBFM),NPF(9,NFM),NRE(NEM),
     '  IBT(3,NIM,NBFM),NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),nx,
     '  NEELEM(0:NE_R_M,0:NRM),
     '  LGE(NHM*NSM,NRCM),NAN(NIM,NAM,NBFM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NW(NEM,3),nr,NRLIST(0:NRM),
     '  IDO(NKM,NNM,0:NIM,NBFM),
     '  Z_CONT_LIST(NDM,2,7)
      REAL*8 ALPHA,CE(NMM,NEM),CG(NMM,NGM),
     '  CGE(NMM,NGM,NEM),CP(NMM,NPM),GRR(NOM),
     '  XA(NAM,NJM,NEM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM),
     '  CURVCORRECT(2,2,NNM,NEM),
     '  YP(NYM,NIYM),
     '  SE(NSM,NBFM,NEM),
     '  RE1(NSM,NHM),XG(NJM,NUM),REF_ENERGY,RG(NGM),
     '  PG(NSM,NUM,NGM,NBM),FEXT(NIFEXTM,NGM,NEM),WG(NGM,NBM),
     '  ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),Z_CONT(NDM,2,67)
      LOGICAL FIX(NYM,NIYFIXM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER counter,LSSTAT
      REAL*8 CGOLD,F0,FALPHA,FDUM,AX,BX,CX,FA,FB,FC,FTOL,ETA,F1
      PARAMETER (CGOLD=0.381966011250D0)

      CALL ENTERS('LINMIN',*9999)

C *** XSL 18Aug2010
C      ALPHA=0.0D0
C      CALL LSFUNC(CONY(0,1,1,nr),IBT,IDO,INP,LGE,NAN,NBH,NBHF,
C     '  NBJ,NBJF,NEELEM,NFF,NFFACE,NGAP,NHE,NHP(1,nr),NKB,
C     '  NKEF,NKH(1,1,1,nr),NKHE,NKJE,NNB,NNF,NONY(0,1,1,nr),
C     '  NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP(1,1,1,nr),
C     '  NVJE,NW,nx,NXI,NYNE,NYNO(0,1,1,nr),NYNP,NYNR(0,0,1,nr),
C     '  Z_CONT_LIST,
C     '  CE,CG,CGE,CP,CURVCORRECT,FEXT,FIX,PG,
C     '  RE1,RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,
C     '  Z_CONT,ZE,
C     '  ZE1,ZP,ZP1,ALPHA,F0,ERROR,*9999)
C      ETA=0.5D0                 !sets accuracy of linear minimisation
C      FTOL=ETA*F0
C      AX=0.0D0
C      FA=F0
C      BX=0.001D0
C      IF(IWRIT4(nr,nx).EQ.1) THEN
C        WRITE(OP_STRING,'(/'' Line Minimisation FTOL='',D14.7,/'
C     '    //''' Bracketing Minimum...'')')FTOL
C        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C      ENDIF
C      CALL MNBRAK(CONY(0,1,1,nr),IBT,IDO,INP,LGE,NAN,NBH,NBHF,
C     '  NBJ,NBJF,NEELEM,NFF,NFFACE,NGAP,NHE,NHP(1,nr),NKB,NKEF,
C     '  NKH(1,1,1,nr),NKHE,NKJE,NNB,NNF,NONY(0,1,1,nr),NPF,
C     '  NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP(1,1,1,nr),NVJE,
C     '  NW,nx,NXI,NYNE,NYNO(0,1,1,nr),NYNP,NYNR(0,0,1,nr),Z_CONT_LIST,
C     '  CE,CG,CGE,CP,CURVCORRECT,FEXT,FIX,PG,RE1,RG,SE,WG,
C     '  XA,XG,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,ZE,ZE1,ZP,ZP1,
C     '  AX,BX,CX,FA,FB,FC,F1,FTOL,LSSTAT,ERROR,*9999)
C      IF(LSSTAT.EQ.0) THEN
C        IF(IWRIT4(nr,nx).EQ.1) THEN
C          WRITE(OP_STRING,'(/''  Minimizing with Brents '
C     '      //'Method...'')')
C        ENDIF
C        CALL BRENT(CONY(0,1,1,nr),IBT,IDO,INP,LGE,NAN,NBH,NBHF,
C     '    NBJ,NBJF,NEELEM,NFF,NFFACE,NGAP,NHE,NHP(1,nr),NKB,NKEF,
C     '    NKH(1,1,1,nr),NKHE,NKJE,NNB,NNF,NONY(0,1,1,nr),NPF,
C     '    NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP(1,1,1,nr),NVJE,
C     '    NW,nx,NXI,NYNE,NYNO(0,1,1,nr),NYNP,NYNR(0,0,1,nr),
C     '    Z_CONT_LIST,
C     '    CE,CG,CGE,CP,CURVCORRECT,FEXT,FIX,PG,RE1,RG,
C     '    SE,WG,
C     '    XA,XG,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,ZE,ZE1,ZP,ZP1,
C     '    AX,BX,CX,FB,FTOL,ALPHA,FALPHA,LSSTAT,ERROR,*9999)
C        IF(LSSTAT.EQ.2) THEN
C          IF(FALPHA.GT.F1) THEN
C            ALPHA=1.0D0
C            CALL LSFUNC(CONY(0,1,1,nr),IBT,IDO,INP,LGE,NAN,NBH,
C     '        NBHF,NBJ,NBJF,NEELEM,NFF,NFFACE,NGAP,NHE,NHP(1,nr),NKB,
C     '        NKEF,NKH(1,1,1,nr),NKHE,NKJE,NNB,NNF,NONY(0,1,1,nr),NPF,
C     '        NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP(1,1,1,nr),
C     '        NVJE,NW,nx,NXI,NYNE,NYNO(0,1,1,nr),NYNP,
C     '        NYNR(0,0,1,nr),Z_CONT_LIST,CE,CG,CGE,CP,
C     '        CURVCORRECT,FEXT,FIX,PG,
C     '        RE1,RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,
C     '        Z_CONT,ZE,
C     '        ZE1,ZP,ZP1,ALPHA,FDUM,ERROR,*9999)
C          ENDIF
C        ENDIF
C      ELSE IF(LSSTAT.EQ.1) THEN
C        ALPHA=BX
C      ELSE IF(LSSTAT.EQ.2) THEN
C        ALPHA=1.0D0
C        CALL LSFUNC(CONY(0,1,1,nr),IBT,IDO,INP,LGE,NAN,NBH,NBHF,
C     '    NBJ,NBJF,NEELEM,NFF,NFFACE,NGAP,NHE,NHP(1,nr),NKB,NKEF,
C     '    NKH(1,1,1,nr),NKHE,NKJE,NNB,NNF,NONY(0,1,1,nr),NPF,NPNE,
C     '    NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP(1,1,1,nr),NVJE,
C     '    NW,nx,NXI,NYNE,NYNO(0,1,1,nr),NYNP,NYNR(0,0,1,nr),
C     '    Z_CONT_LIST,
C     '    CE,CG,CGE,CP,CURVCORRECT,FEXT,FIX,PG,RE1,RG,SE,WG,
C     '    XA,XG,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,ZE,ZE1,ZP,ZP1,
C     '    ALPHA,FDUM,ERROR,*9999)
C      ENDIF
C      IF(IWRIT4(nr,nx).EQ.1) THEN
C        WRITE(OP_STRING,'('' ALPHA = '',D12.3)') ALPHA
C        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C      ENDIF

C 14/07/08 XSL initialise counter to count the number of residual evaluations
C i.e. the number of times LSFUNC/ZPRP has been called
      counter=0
      AX=0.0D0
C 27/11/08 XSL NEWS When using energy norm as convergence check,
C use reference energy saved from nonlin to calculate tolerance
C otherwise FA=R.(delta*alpha), where alpha=0!
      IF(KTYP007.EQ.4) THEN !energy norm
        FA=REF_ENERGY
      ELSE
        CALL LSFUNC(CONY,IBT,IDO,
     &   INP,ITER1,LGE,NAN,NBH,NBHF,
     '    NBJ,NBJF,NEELEM,NFF,NFFACE,NGAP,NHE,NHP,NKB,
     '    NKEF,NKH,NKHE,NKJE,NNB,NNF,NONY,
     '    NPF,NPNE,NPNODE,NPNY,nr,NRE,NRLIST,NSB,NVHE,NVHP,
     '    NVJE,NW,nx,NXI,NYNE,NYNO,NYNP,NYNR,
     '    Z_CONT_LIST,
     '    CE,CG,CGE,CP,CURVCORRECT,FEXT,FIX,GRR,PG,
     '    RE1,RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,
     '    Z_CONT,ZE,
     '    ZE1,ZP,ZP1,AX,FA,ERROR,*9999)
        counter=counter+1
      ENDIF !KTYP007
C 27/11/08 NEWE
      ETA=0.5D0                 !sets accuracy of linear minimisation
      FTOL=ETA*FA
C      IF(IWRIT4(nr,nx).EQ.1) THEN
C        WRITE(OP_STRING,'(/'' Line Minimisation FTOL='',D14.7,/'
C     '    //''' Bracketing Minimum...'')')FTOL
C        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C      ENDIF

C 14/07/08 NEWS XSL
C Use Brent's method to find the min
C The starting bracket is defaulted to be [0,1]
C AX=0, CX=1, BX=1st golden section pt inside this interval
C 14/07/08 XSL Set starting alpha value to be 1, i.e. full Newton step
      CX=1.0D0
      CALL LSFUNC(CONY,IBT,IDO,
     & INP,ITER1,LGE,NAN,NBH,NBHF,
     & NBJ,NBJF,NEELEM,NFF,NFFACE,NGAP,NHE,NHP,NKB,
     & NKEF,NKH,NKHE,NKJE,NNB,NNF,NONY,
     '  NPF,NPNE,NPNODE,NPNY,nr,NRE,NRLIST,NSB,NVHE,
     '  NVHP,NVJE,NW,nx,NXI,
     '  NYNE,NYNO,NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,CP,
     '  CURVCORRECT,FEXT,FIX,GRR,
     '  PG,RE1,RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,
     '  ZE,ZE1,ZP,ZP1,
     '  CX,FC,ERROR,*9999)
      counter=counter+1
      IF(FC .LT. FTOL) THEN
        ALPHA=CX
        GOTO 9998
      ENDIF

C XSL Golden section method
      BX=AX+CGOLD*(CX-AX)
      CALL LSFUNC(CONY,IBT,IDO,
     & INP,ITER1,LGE,NAN,NBH,NBHF,
     & NBJ,NBJF,NEELEM,NFF,NFFACE,NGAP,NHE,NHP,NKB,
     & NKEF,NKH,NKHE,NKJE,NNB,NNF,NONY,
     '  NPF,NPNE,NPNODE,NPNY,nr,NRE,NRLIST,NSB,NVHE,
     '  NVHP,NVJE,NW,nx,NXI,
     '  NYNE,NYNO,NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,CP,
     '  CURVCORRECT,FEXT,FIX,GRR,
     '  PG,RE1,RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,
     '  ZE,ZE1,ZP,ZP1,
     '  BX,FB,ERROR,*9999)
      counter=counter+1
      IF(FB .LT. FTOL) THEN
        ALPHA=BX
        GOTO 9998
      ENDIF
C NEWE

C 14/07/08 XSL The outmost IF loop is not needed if MNBRAK is not called
C      IF(LSSTAT.EQ.0) THEN
C      IF(IWRIT4(nr,nx).EQ.1) THEN
C        WRITE(OP_STRING,'(/''  Minimizing with Brents '
C     '    //'Method...'')')
C      ENDIF
      CALL BRENT(CONY,counter,IBT,IDO,
     & INP,ITER1,LGE,NAN,NBH,
     '  NBHF,NBJ,NBJF,NEELEM,NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,
     '  NKH,NKHE,NKJE,NNB,NNF,NONY,NPF,
     '  NPNE,NPNODE,NPNY,nr,NRE,NRLIST,NSB,NVHE,NVHP,NVJE,
     '  NW,nx,NXI,NYNE,NYNO,NYNP,NYNR,
     '  Z_CONT_LIST,
     '  CE,CG,CGE,CP,CURVCORRECT,FEXT,FIX,GRR,PG,RE1,RG,
     '  SE,WG,
     '  XA,XG,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,ZE,ZE1,ZP,ZP1,
     '  AX,BX,CX,FB,FTOL,ALPHA,FALPHA,LSSTAT,ERROR,*9999)
      IF(LSSTAT.EQ.2) THEN
        IF(FALPHA.GT.F1) THEN
C XSL Can't find min after a certain number of iterations
C Take a very small Newton Step
          ALPHA=0.000001D0
          CALL LSFUNC(CONY,IBT,IDO,
     &     INP,ITER1,LGE,NAN,NBH,
     '      NBHF,NBJ,NBJF,NEELEM,NFF,NFFACE,NGAP,NHE,NHP,NKB,
     '      NKEF,NKH,NKHE,NKJE,NNB,NNF,NONY,NPF,
     '      NPNE,NPNODE,NPNY,nr,NRE,NRLIST,NSB,NVHE,NVHP,
     '      NVJE,NW,nx,NXI,NYNE,NYNO,NYNP,
     '      NYNR,Z_CONT_LIST,CE,CG,CGE,CP,
     '      CURVCORRECT,FEXT,FIX,GRR,PG,
     '      RE1,RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,
     '      Z_CONT,ZE,
     '      ZE1,ZP,ZP1,ALPHA,FDUM,ERROR,*9999)
          counter=counter+1
        ENDIF
      ENDIF !LSSTAT=2

C      IF(IWRIT4(nr,nx).EQ.1) THEN
C        WRITE(OP_STRING,'('' ALPHA = '',D12.3)') ALPHA
C        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C      ENDIF

 9998 CALL EXITS('LINMIN')
      RETURN
 9999 CALL ERRORS('LINMIN',ERROR)
      CALL EXITS('LINMIN')
      RETURN 1
      END
