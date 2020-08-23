C#### Module: CMISS_ARCHIVE0
C###  Description:
C###    Contains archived code from common include files, plus
C###    code which used BFGS  nonlinear solver option

CCOMM File       B09.CMN
C     File       B21.CMN
C     File       B22.CMN
C     File       CBDI01.CMN
C     File       CBSY01.CMN
C     File       cbzm00.cmn
C     File       CROS01.CMN
C     File       CSPI00.CMN
C     File       CURR00.CMN
C     File       DFT00.CMN
C     File       dial00.cmn
C     File       ditr.cmn
C     File       EDATA00.CMN
C     File       fgbez00.cmn
C     File       I03.CMN
C     File       ISO000.CMN
C     File       IMPPL01.CMN
C     File       LINE01.CMN
C     File       LOOP00.CMN
C     File       map3d.cmn
C     File       MATR00.CMN
C     File       ODE000.CMN
C     File       OPNOD00.CMN
C     File       OUTP00.CMN
C     File       OXS000.CMN
C     File       OXS002.CMN
C     File       PAPER00.CMN
C     File       pics00.cmn
C     File       pics01.cmn
C     File       REFI00.CMN
C     File       SOLV00.CMN
C     File       SOLV01.CMN
C     File       VIEW00.CMN
C     File       XPGK00.CMN
C     code which used BFGS  nonlinear solver option

File B09.CMN
================

!  b09.cmn
      REAL*8 PLTMIN(2),PLTMAX(2),PLTTOL,PENDIM(8),DIST
      COMMON /B09R/ PLTMIN,PLTMAX,PLTTOL,PENDIM,DIST

File B21.CMN
================

!  b21.cmn
      CHARACTER FOPTM*8,TOPTM*79
      COMMON /B21C/ FOPTM,TOPTM

File B22.CMN
================

!  b22.cmn
      CHARACTER FGEOM*8,TGEOM*79                            
      COMMON /B22C/ FGEOM,TGEOM                            

File CBDI01.CMN
================

!  cbdi01.cmn
      INTEGER NXCH
      COMMON /CBDI01I/ NXCH

File CBSY01.CMN
================

!  cbsy01.cmn
      REAL*8 MXIN,RXRE,RNRE,E
      COMMON /CBSY01R/ MXIN,RXRE,RNRE,E

File cbzm00.cmn
================

!  cbzm00.cmn
      INTEGER IZOOM
      COMMON /CBZM00I/ IZOOM(4)


File CROS01.CMN
================

!  cros01.cmn
      REAL*8 ZCROS(3)
      COMMON /CROS01R/ ZCROS


File CSPI00.CMN
================
!  cspi00.cmn
!      INTEGER MPLUN                 ! is MAP logical unit number
      LOGICAL MAPOPN
!      COMMON /CSPI0OI/ MPLUN
      COMMON /CSPI00L/ MAPOPN

File CURR00.CMN
===============
!  curr00.cmn
      INTEGER*2 COLOUR
      COMMON /CURR00I2/ COLOUR

File DFT00.CMN
===============

!  dft00.cmn
      REAL*8 ZDATA(500,2)
      COMMON /DFT00R/ ZDATA


File dial00.cmn
===============

!  dial00.cmn
      INTEGER OLD_HANDLER
      COMMON /DIAL00I/ OLD_HANDLER
      LOGICAL 
     '  FATAL_HANDLER,CHANGE_HANDLER,MENU,SMG_READ
      COMMON /DIAL00L/ 
     '  FATAL_HANDLER,CHANGE_HANDLER,MENU,SMG_READ

      File ditr00.cmn
===============

!  ditr00.cmn
      INTEGER ISIGNAL,ISIGLAST
      REAL*8 DATUM,FIDMARK(200),TTSTART,TTEND
      COMMON /DITR00I/ ISIGNAL,ISIGLAST
      COMMON /DITR00R/ DATUM,FIDMARK,TTSTART,TTEND


File EDATA00.CMN
================

!  edata00.cmn
C     note these need to be maintained identical with FE20
      INTEGER NAMX,NBMX,NEMX,NIMX,NJMX,NKMX,NNMX,NSMX
      PARAMETER (NAMX=1,NBMX=10,NEMX=60,NIMX=3,NJMX=6,NKMX=4,
     '  NNMX=16,NSMX=56)
C     Note: PJH 22-jan-92 changed ne to n_e
      INTEGER IBT2(2,NIMX,NBMX),IDO2(NKMX,0:NIMX,NBMX),
     '  INP2(NNMX,NIMX,NBMX),NBJ2(NJMX,NEMX),
     '  NAN2(NIMX,NAMX,NBMX),N_E,NITER
      COMMON /EDATA00I/ IBT2,IDO2,INP2,NBJ2,NAN2,N_E,NITER
      REAL*8 U1(NJMX),U2(NJMX),XENE(NSMX,NJMX)
      COMMON /EDATA00R/ U1,U2,XENE


File I03.CMN
=============

!  fgbez00.cmn
      INTEGER NSIMP
      REAL*8 XA(4),XB(4),XKSIA,SLEN,ACDIFF,YA(4),YB(4),XXA,YYA,XXB,YYB
      COMMON /FGBEZ00I/ NSIMP
      COMMON /FGBEZ00R/ XA,XB,XKSIA,SLEN,ACDIFF,YA,YB,XXA,YYA,XXB,YYB


File I03.CMN
=============

!  i03.cmn
      INTEGER NJTI
      COMMON /I03I/ NJTI


File IMPPL01.CMN
=================

!  imppl01.cmn
      REAL*8 ZVAL_MIN,ZVAL_MAX
      

File ISO000.CMN
================

!  iso000.cmn
      LOGICAL ISOAREA
      COMMON /ISO000L/ ISOAREA
      COMMON /IMPPL01R/ ZVAL_MIN,ZVAL_MAX


File LINE01.CMN
================

!  line01.cmn
      LOGICAL FIRST_DEF_LINE
      COMMON /LINE01L/ FIRST_DEF_LINE


File LOOP00.CMN
================

!  loop00.cmn
      LOGICAL LOOP
      COMMON /LOOP00L/ LOOP


File map3d.cmn
================

!  map3d.cmn
      INTEGER MAXNUMLEADS,MAXNUMMAPS
      PARAMETER(MAXNUMLEADS=256,MAXNUMMAPS=512)
      INTEGER DALBSPM
      PARAMETER (DALBSPM=0)
      INTEGER CHANNELS(MAXNUMLEADS),DATASET_OFFSET_O,
     '  NT_CHANNELS,NT_DATASETS,NT_LEADS,NT_SUB_DIV
      COMMON /MAP3DI/ CHANNELS,DATASET_OFFSET_O,
     '  NT_CHANNELS,NT_DATASETS,NT_LEADS,NT_SUB_DIV
C      CHARACTER DATAFILENAME*80,GEOMFILENAME*80,HEADERTEXT*256
C      COMMON /MAP3DC/ DATAFILENAME,GEOMFILENAME,HEADERTEXT
C      LOGICAL QPREFS
C      COMMON /MAP3DL/ QPREFS
C      LOGICAL*1 QXFIDS
C      COMMON /MAP3DL1/ QXFIDS


File MATR00.CMN
================

!  matr00.cmn
      REAL*8 
     '  B(2,2),C(2,2),E(2,2),F(2,2),P(2,2),R(2,2),T(2,2),U(2,2),V(2,2)
      COMMON /MATR00R/ 
     '  B,C,E,F,P,R,T,U,V


File ODE000.CMN
================

!  ode000.cmn
      LOGICAL ODE
      COMMON /ODE000L/ ODE

File OPNOD00.CMN
================

!  opnod00.cmn
      REAL*8 ZCOR(4)
      COMMON /OPNOD00R/ ZCOR

File OUTP00.CMN
===============

!  outp00.cmn
      REAL*8 d_Conc_N2_dVol_max
      COMMON /OUTP00R/ d_Conc_N2_dVol_max

File OXS000.CMN
================

!  oxs000.cmn
      LOGICAL FIRST_OXS
      COMMON /OXS000L/ FIRST_OXS

File OXS002.CMN
===============
!  oxs002.cmn
      REAL*8 F1
      COMMON /OXS002R/ F1

File PAPER00.CMN
================

!  paper00.cmn
      LOGICAL PORTRT
      COMMON /PAPER00L/ PORTRT
 

File pics00.CMN
================
!  pics00.cmn
C cpb 16/11/94 Moving I2P to cope with images
C      INTEGER*2 I2P(200,2048,4)
      INTEGER I2P(512,512,4)
      ! NUMSIGNALS must be <= 200, NUMSAMPLES <= 2048
      COMMON /PICS00I/ I2P

File pics01.CMN
================
!  pics01.cmn
      INTEGER IMGX,IMGY
      REAL*8 RINTENSITY_MAX
      COMMON /PICS01I/ IMGX,IMGY
      COMMON /PICS01R/ RINTENSITY_MAX


File REFI00.CMN
================

!  refi00.cmn
      LOGICAL SUBDIV
      COMMON /REFI00L/ SUBDIV


File SOLV00.CMN
================

!  solv00.cmn
      INTEGER INLET(0:20),OUTLET(0:20)
      REAL*8 APOLE,BPOLE,CPOLE,DPOLE
      COMMON /SOLV00I/ INLET,OUTLET
      COMMON /SOLV00R/ APOLE,BPOLE,CPOLE,DPOLE


File SOLV01.CMN
================

!  solv01.cmn
      COMPLEX*16 APOLE,BPOLE,CPOLE,DPOLE
      COMMON /SOLV01C16/ APOLE,BPOLE,CPOLE,DPOLE


File VIEW00.CMN
================

!  view00.cmn
      INTEGER ISVIEW
      COMMON /VIEW00I/ ISVIEW


File VSA00.CMN
================

!  vsa00.cmn
      INTEGER IPRI,IPRLEV,IPRESS,MSTOP,MSTART,MODIFY,
     '        IPRGOM,IPRNAD,IPRWAK,IPRCPV,IPRPPI,
     '        MODE,NPNMAX,NRBMAX,ITGSMX,IMERGE,NSUB,NSPMAX,NPCMAX,
     '        NWIT,NVPI,IBLTYP,NT,NHC,
     '        NORSET,NVORT,NPASUM,JETPAN,NBCHGE
      REAL*8 RSYM,RGPR,RNF,RFF,RCORE,SOLRES,TOL,
     '     ALDEG,YAWDEG,RMACH,VMOD,COMFAC,
     '     ALBAR,RFREQ,HX,HY,HZ,
     '     CBAR,SREF,SSPAN,RMPX,RMPY,RMPZ
      COMMON /VSA00I/ IPRI,IPRLEV,IPRESS,MSTOP,MSTART,MODIFY,
     '             IPRGOM,IPRNAD,IPRWAK,IPRCPV,IPRPPI,
     '             MODE,NPNMAX,NRBMAX,ITGSMX,IMERGE,NSUB,NSPMAX,NPCMAX,
     '             NWIT,NVPI,IBLTYP,NT,NHC,
     '             NORSET,NVORT,NPASUM,JETPAN,NBCHGE
      COMMON /VSA00R/ RSYM,RGPR,RNF,RFF,RCORE,SOLRES,TOL,
     '             ALDEG,YAWDEG,RMACH,VMOD,COMFAC,
     '             ALBAR,RFREQ,HX,HY,HZ,
     '             CBAR,SREF,SSPAN,RMPX,RMPY,RMPZ

File VSB00.CMN
================

!  vsb00.cmn
      INTEGER NTPAT,IDENT(9),MAKE(9),KOMP(9),KLASS(9),
     '        NTSECT(9),NTNODE(9),NODS(20,10,9),
     '        INMODE(9),NODES(10,9),NPS(9),             
     '        NODEC(9),NPC(9),INTC(9),MOVE(9)
      REAL*8 CTX(9),CTY(9),CTZ(9),SCAL(9),THET(9),
     '     CPX(9),CPY(9),CPZ(9),
     '     CHX(9),CHY(9),CHZ(9),
     '     STX(9),STY(9),STZ(9),SCALE(9),ALF(9),THETA(9)
      COMMON /VSB00I/ NTPAT,IDENT,MAKE,KOMP,KLASS,
     '             NTSECT,NTNODE,NODS,
     '             INMODE,NODES,NPS,             
     '             NODEC,NPC,INTC,MOVE
      COMMON /VSB00R/ CTX,CTY,CTZ,SCAL,THET,
     '             CPX,CPY,CPZ,CHX,CHY,CHZ,
     '             STX,STY,STZ,SCALE,ALF,THETA

File VSC00.CMN
================

!  vsc00.cmn
      INTEGER NODE,NPCW,INTCW,MARK,
     '        IDENTW(9),IFLEXW(9),IDEFW(9),
     '        KWPACH(9),KWSIDE(9),KWLINE(9),KWPAN1(9),KWPAN2(9),
     '        INPUT(9),NODEWS(2,9),NTSEPS(9),
     '        NODEWC(9),NPCP(9),INTCP(9)
      REAL*8 XMIN,XMAX,SWPX(20,9),SWPY(20,9),DELTAZ(20,9)
      COMMON /VSC00I/ NODE,NPCW,INTCW,MARK,
     '             IDENTW,IFLEXW,IDEFW,
     '             KWPACH,KWSIDE,KWLINE,KWPAN1,KWPAN2,
     '             INPUT,NODEWS,NTSEPS,
     '             NODEWC,NPCP,INTCP
      COMMON /VSC00R/ XMIN,XMAX,
     '             SWPX,SWPY,DELTAZ

File VSD00.CMN
================

!  vsd00.cmn
      INTEGER KP(100),NSTRM
      REAL*8 F(100)
      COMMON /VSD00I/ KP,NSTRM
      COMMON /VSD00R/ F

File VSE00.CMN
================

!  vse00.cmn
      REAL*8 RND,TRIPUP,TRIPOP,XPRINT,XSKIP
      COMMON /VSE00R/ RND,TRIPUP,TRIPOP,XPRINT,XSKIP

File VSF00.CMN
================

!  vsf00.cmn
      INTEGER 
     '  MOLD(100),MEET(100),NEAR(100),
     '  INCPRI(100),INCPRJ(100),INCPRK(100),
     '  NP(100),NP1(100),NP2(100),NP3(100),NBOX
      COMMON /VSF00I/ 
     '  MOLD,MEET,NEAR,
     '  INCPRI,INCPRJ,INCPRK,
     '  NP,NP1,NP2,NP3,NBOX
      REAL*8 
     '  X0(100),Y0(100),Z0(100),X1(100),Y1(100),Z1(100),
     '  X2(100),Y2(100),Z2(100),X3(100),Y3(100),Z3(100),
     '  AL1(10,100),AL2(10,100),AL3(10,100)
      COMMON /VSF00R/ 
     '  X0,Y0,Z0,X1,Y1,Z1,
     '  X2,Y2,Z2,X3,Y3,Z3,
     '  AL1,AL2,AL3

File VSG00.CMN
================

!  vsg00.cmn
      INTEGER NEAR2(100),NPOINT
      REAL*8 RSX(100),RSY(100),RSZ(100),SU(100),SD(100),DELS(100)
      COMMON /VSG00I/ NEAR2,NPOINT
      COMMON /VSG00R/ RSX,RSY,RSZ,SU,SD,DELS

File XPGK00.CMN
================

!  xpgk00.cmn
      INTEGER IEQUATION
      REAL*8 S
      COMMON /XPGK00I/ IEQUATION
      COMMON /XPGK00R/ S


C     code which used BFGS nonlinear solver option
==================================================
In FE20:
     '  QND(NOMX*NQNMX*USE_GRID+1),
     '  QNG(NOMX*NQNMX*USE_GRID+1),
     '  QNP(NQNMX*USE_GRID+1),
     '  QNA((NQNMX+1)*USE_GRID+1), ! QNA(0:nqn)
In FE21:Solve
        QNA(0:*),QND(NOM,*),QNG(NOM,*),QNP(*),
In FE21:Gensol
        QNA(0:*),QND(NOM,*),QNG(NOM,*),QNP(*),
In FE21:Nonlin
        QNA(0:*),QND(NOM,*),QNG(NOM,*),QNP(*),
          IF(KTYP9.EQ.3.AND.ITER2.GT.0) THEN !BFGS inverse
            DO NY=1,NYT(1,nr,nx)
              DO NOY=1,NONY(0,NY)
                NO=NONY(NOY,NY)
                QND(NO,ITER2) = YP(NY,1)
                QNG(NO,ITER2) = GR1(NY)-GR2(NY)
              ENDDO
            ENDDO
            SUM=0.0D0
            DO NO=1,NOT(nr,nx)
              SUM=SUM+QND(NO,ITER2)*QNG(NO,ITER2)
            ENDDO
            QNP(ITER2)=1.0D0/SUM
            DO NY=1,NYT(1,nr,nx)
              DO NOY=1,NONY(0,NY)
                NO=NONY(NOY,NY)
                YP(NO,1)=GR1(NY)
              ENDDO
            ENDDO
            DO NQN=ITER2-1,0,-1
              SUM=0.0D0
              DO NO=1,NOT(nr,nx)
                SUM=SUM+QND(NO,NQN+1)*YP(NO,1)
              ENDDO
              QNA(NQN)=QNP(NQN+1)*SUM
              DO NO=1,NOT(nr,nx)
                YP(NO,1)=YP(NO,1)-QNG(NO,NQN+1)*QNA(NQN)
              ENDDO
              IF(NQN.EQ.0) THEN
                DO NO=1,NOT(nr,nx)
                  GRR(NO)=YP(NO,1)
                ENDDO
              ENDIF
            ENDDO
          ELSE
          ENDIF
          IF(KTYP9.EQ.3.AND.ITER2.GT.0) THEN !BFGS inverse

            DO NQN=1,ITER2
              SUM=0.0D0
              DO NO=1,NOT(nr,nx)
                SUM=SUM+QNG(NO,NQN)*GRR(NO)
              ENDDO
              BETA=QNP(NQN)*SUM
              DO NO=1,NOT(nr,nx)
                GRR(NO)=GRR(NO)+QND(NO,NQN)*(QNA(NQN-1)-BETA)
              ENDDO
            ENDDO

          ENDIF
