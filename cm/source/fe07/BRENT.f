      SUBROUTINE BRENT(CONY,counter,IBT,IDO,
     & INP,ITER1,LGE,NAN,NBH,
     '  NBHF,NBJ,NBJF,NEELEM,NFF,
     '  NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,NNF,NONY,NPF,
     '  NPNE,NPNODE,NPNY,nr,NRE,NRLIST,
     '  NSB,NVHE,NVHP,NVJE,NW,nx,NXI,NYNE,NYNO,NYNP,
     '  NYNR,Z_CONT_LIST,CE,CG,CGE,CP,CURVCORRECT,
     '  FEXT,FIX,GRR,PG,RE,RG,SE,
     '  WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,ZE,ZE1,ZP,ZP1,AX,
     '  BX,CX,FB,
     '  FTOL,ALPHA,FALPHA,LSSTAT,ERROR,*)

C#### Subroutine: BRENT
C###  Description:
C###    BRENT zeros in on a minimum, given a bracketed minimum
C###    (e.g. output from mnbrak). From Numerical Recipes.  Modified
C###    to stop search if function value is below ftol, or if function
C###    value is within 1% of previous iteration.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER counter,IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ITER1,LGE(NHM*NSM,NRCM),LSSTAT,
     & NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NFF(6,NEM),NFFACE(0:NF_R_M,NRM),
     '  NGAP(NIM,NBM),NHE(NEM),NHP(NPM,0:NRM),NKB(2,2,2,NNM,NBFM),
     '  NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),
     '  nr,NRE(NEM),NRLIST(0:NRM),NSB(NKM,NNM,NBFM),
     & NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),
     '  nx,NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     & NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM),NONY(0:NOYM,NYM,NRCM,0:NRM),
     '  Z_CONT_LIST(NDM,2,7)
      REAL*8 ALPHA,AX,BX,CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),
     '  CP(NMM,NPM),CURVCORRECT(2,2,NNM,NEM),CX,
     '  FALPHA,FB,FEXT(NIFEXTM,NGM,NEM),
     '  FTOL,GRR(NOM),PG(NSM,NUM,NGM,NBM),
     '  RE(NSM,NHM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),
     '  YP(NYM,NIYM),ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),
     '  Z_CONT(NDM,2,67),ZE(NSM,NHM),ZE1(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER iter,ITMAX
      REAL*8 A,B,CGOLD,D,E,ETEMP,FU,FX,FV,FW,P,Q,R,RATIO,TOL,TOL1,TOL2,
     '  U,V,W,X,XM,ZEPS
      PARAMETER (ITMAX=50, CGOLD=0.381966011250D0, ZEPS=1.0D-10,
     &  TOL=1.0D-8)

      CALL ENTERS('BRENT',*9999)

      LSSTAT=0
      A=DMIN1(AX,CX)
      B=DMAX1(AX,CX)
      V=BX
      W=V
      X=V
      E=0.0D0
      FX=FB
      FV=FX
      FW=FX

      DO 11 iter=1,ITMAX
        IF(DOP) THEN
          WRITE(OP_STRING,'(''Iteration '',I2)') iter
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        XM=0.5D0*(A+B)
        TOL1=TOL*DABS(X)+ZEPS
        TOL2=2.0D0*TOL1
        IF(DABS(X-XM).LE.(TOL2-0.5D0*(B-A))) GOTO 3
        IF(DABS(E).GT.TOL1) THEN ! Is a parabolic possible
          R=(X-W)*(FX-FV) ! Yes, so fit parabola
          Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
          Q=2.0d0*(Q-R)
          IF(Q.GT.0.0d0) P=-P
          Q=DABS(Q)
          ETEMP=E
C Is the parabolic acceptable?
          IF(DABS(P).GE.DABS(0.5D0*Q*ETEMP).OR.P.LE.Q*(A-X).OR.
     '      P.GE.Q*(B-X)) GOTO 1
          D=P/Q ! Yes, parabolic interpolation step
          E=D
          U=X+D
          IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=DSIGN(TOL1,XM-X) !f must not be evaluated too close to ax or bx
          GOTO 2
        ENDIF
1       IF(X.GE.XM) THEN ! Not acceptable parabolic, must do a golden section step
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
2       IF(DABS(D).GE.TOL1) THEN ! The function must not be evaluated too close to X
          U=X+D
        ELSE
          U=X+DSIGN(TOL1,D)
        ENDIF

        CALL LSFUNC(CONY,IBT,IDO,
     &  INP,ITER1,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '    NEELEM,NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,
     '    NNF,NONY,NPF,NPNE,NPNODE,NPNY,nr,NRE,NRLIST,NSB,NVHE,NVHP,
     '    NVJE,NW,nx,NXI,NYNE,NYNO,NYNP,NYNR,
     '    Z_CONT_LIST,CE,CG,CGE,CP,
     '    CURVCORRECT,FEXT,FIX,GRR,PG,RE,RG,SE,WG,XA,XG,XP,
     '    YG,YGF,YP,
     '    ZA,ZA1,Z_CONT,ZE,ZE1,ZP,ZP1,U,FU,ERROR,*9999)
C XSL NEWS 18Aug2010
        counter=counter+1
C XSL NEWE
       
        IF(DOP) THEN
          WRITE(OP_STRING,'(''U  ='',D14.7,'' FU ='',D14.7)') U,FU
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(FU .LT. FTOL) THEN
          X=U
          FX=FU
          LSSTAT=1
          GOTO 3
        ENDIF

        RATIO=FU/FX
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
          IF(RATIO .GT. 0.95D0) THEN ! Check if the residual is within 1% of the previous value
            LSSTAT=1
            GOTO 3
          ENDIF
        ELSE
          IF(RATIO .LT. 1.05D0 .AND. X .ne. 0.0D0) THEN ! Check if the residual is within 1% of the previous value
            LSSTAT=1
            GOTO 3
          ENDIF
          IF(U.LT.X) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU.LE.FW .OR. W.EQ.X) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF
11    CONTINUE
      LSSTAT=2

3     ALPHA=X
      FALPHA=FX

      CALL EXITS('BRENT')
      RETURN
 9999 CALL ERRORS('BRENT',ERROR)
      CALL EXITS('BRENT')
      RETURN 1
      END


