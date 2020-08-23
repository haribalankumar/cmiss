      SUBROUTINE ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,NFF,
     '  NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,NNF,NPF,
     '  NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP,NVJE,NW,nx,NXI,NYNE,NYNP,
     '  NYNR,Z_CONT_LIST,CE,CG,CGE,
     '  CP,CURVCORRECT,FEXT,PG,RE,RG,SE,
     '  WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,ZAA,Z_CONT,ZE,ZE1,ZP,ZP1,ZPA,FIX,
     '  RET_ERROR,*)

C#### Subroutine: ZPRP
C###  Description:
C###    ZPRP calculates global residual vector YP(ny,4) at current
C###    solution.  RE(ns,nh) has been corrected with scaling factor
C###    SE(ns,nb,ne)

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'error0.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'mxch.inc'
      INCLUDE 'nonl00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List (LGE is not used)
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LGE(NHM*NSM,NRCM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NFF(6,NEM),NFFACE(0:NF_R_M,NRM),
     '  NGAP(NIM,NBM),NHE(NEM),NHP(NPM),NKB(2,2,2,NNM,NBFM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NKEF(0:4,16,6,NBFM),
     '  NKH(NHM,NPM,NCM),NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),
     '  nr,NRE(NEM),NSB(NKM,NNM,NBFM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),
     '  NW(NEM,3),NXI(-NIM:NIM,0:NEIM,0:NEM),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNR(0:NY_R_M,0:NRCM,NCM),
     '  Z_CONT_LIST(NDM,2,7)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM,NXM),
     '  CP(NMM,NPM),
     '  CURVCORRECT(2,2,NNM,NEM),FEXT(NIFEXTM,NGM,NEM),
     '  PG(NSM,NUM,NGM,NBM),RE(NSM,NHM),RG(NGM),
     '  SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM),
     '  ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),ZAA(NAM,NHM,NCM,NEM),
     '  Z_CONT(NDM,2,67),
     '  ZE(NSM,NHM),ZE1(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM),
     '  ZP1(NKM,NVM,NHM,NPM,NCM),ZPA(NKM,NVM,NHM,NPM,NCM)
      CHARACTER RET_ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER GETNYR,nc,ne,nf,nhx,nk,noelem,noface,nonode,
     '  no_nynr,np,nv,ny,ny1
c cpb 18/10/96 Adding dynamic allocation for parallel local arrays
      INTEGER*4 LGE_PTR,IDOXFT_PTR,MYMS_PTR,NSFE_PTR,NYNS_PTR,
     '  CG_PTR,RE_PTR,RDF_PTR,RG_PTR,SM_PTR,SN_PTR,
     '  XDF_PTR,XE_PTR,XG_PTR,ZDF_PTR,ZE_PTR,ZE1_PTR,ZG_PTR
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      CHARACTER ERROR*(ERRSTRLEN)
      LOGICAL ERROR_FLAG

      CALL ENTERS('ZPRP',*9999)

C CPB 18/10/96 Intialise pointers so the error free will work properly

      LGE_PTR=0
      CG_PTR=0
      RE_PTR=0
      RG_PTR=0
      XE_PTR=0
      XG_PTR=0
      ZE_PTR=0
      ZE1_PTR=0
      ZG_PTR=0

      IDOXFT_PTR=0
      MYMS_PTR=0
      NSFE_PTR=0
      NYNS_PTR=0
      RDF_PTR=0
      SM_PTR=0
      SN_PTR=0
      XDF_PTR=0
      ZDF_PTR=0

      nc=1 !LHS

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP   CRITICAL(ZPRP_1)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          WRITE(OP_STRING,'('' ZP(nk,nv,nh,'',I3,'',1): '',8D11.3,'
     '      //'/(16X,8D11.3))')
     '      np,(((ZP(nk,nv,NH_LOC(nhx,nx),np,1),
     '      nk=1,NKH(NH_LOC(nhx,nx),np,1)),
     '      nv=1,NVHP(NH_LOC(nhx,nx),np,1)),nhx=1,NHP(np))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$OMP   END CRITICAL(ZPRP_1)
      ENDIF

C     Initialise residuals
      DO no_nynr=1,NYNR(0,1,1) !loop over rows
        ny=NYNR(no_nynr,1,1) !is row number
!PJH 8Apr97 IF(NPNY(0,ny,0).EQ.1) THEN
!          np=NPNY(4,ny,0)
!          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
!        ENDIF
        YP(ny,4)=0.0d0
C 28/06/07 XSL Initialise YP(ny,6) - contact forces
        IF (KTYP5G(nr).GT.0) THEN !contact
          YP(ny,6)=0.0d0   
        ENDIF !KTYP5G(nr)   
      ENDDO !no_nynr

      CALL CPU_TIMER(CPU_USER,TIME_START)

      ERROR_FLAG=.FALSE.
C new MPN 1Feb2000: OMP parallel proc directive
C old
CC$DOACROSS local (ne,noelem,LGE_PTR,CG_PTR,RE_PTR,RG_PTR,XE_PTR,XG_PTR,
CC$&               ZE_PTR,ZE1_PTR,ZG_PTR,ERROR)
CC$&        share (nc,nr,nx,ERROR_FLAG)
C end old
C$OMP PARALLEL DO
C$OMP&  PRIVATE(ne,noelem,LGE_PTR,CG_PTR,RE_PTR,RG_PTR,XE_PTR,
C$OMP&          XG_PTR,ZE_PTR,ZE1_PTR,ZG_PTR,ERROR),
C$OMP&  SHARED(MEM_INIT,nc,NGM,nr,nx,ERROR_FLAG)
      DO noelem=1,NEELEM(0,nr) !is main element loop
        IF(.NOT.ERROR_FLAG) THEN
          ne=NEELEM(noelem,nr)

C CPB 18/10/96 Intialise pointers so they are zero in the parallel loop

          LGE_PTR=0
          CG_PTR=0
          RE_PTR=0
          RG_PTR=0
          XE_PTR=0
          XG_PTR=0
          ZE_PTR=0
          ZE1_PTR=0
          ZG_PTR=0

c cpb 18/10/96  Dynamically allocation parallel local arrays

          CALL ALLOCATE_MEMORY(NHM*NSM*NRCM,1,INTTYPE,LGE_PTR,
     '      MEM_INIT,ERROR,*100)
          CALL ALLOCATE_MEMORY(NMM*NGM,1,DPTYPE,CG_PTR,
     '      MEM_INIT,ERROR,*100)
          CALL ALLOCATE_MEMORY(NSM*NHM,1,DPTYPE,RE_PTR,
     '      MEM_INIT,ERROR,*100)
          CALL ALLOCATE_MEMORY(NGM,1,DPTYPE,RG_PTR,
     '      MEM_INIT,ERROR,*100)
          CALL ALLOCATE_MEMORY(NSM*NJM,1,DPTYPE,XE_PTR,
     '      MEM_INIT,ERROR,*100)
          CALL ALLOCATE_MEMORY(NJM*NUM,1,DPTYPE,XG_PTR,
     '      MEM_INIT,ERROR,*100)
          CALL ALLOCATE_MEMORY(NSM*NHM,1,DPTYPE,ZE_PTR,
     '      MEM_INIT,ERROR,*100)
          CALL ALLOCATE_MEMORY(NSM*NHM,1,DPTYPE,ZE1_PTR,
     '      MEM_INIT,ERROR,*100)
          CALL ALLOCATE_MEMORY(NHM*NUM,1,DPTYPE,ZG_PTR,
     '      MEM_INIT,ERROR,*100)

C cpb 18/10/96 Must create a dynam subroutine as dynamically allocated
C arrays are used within this subroutine level. This can be fixed on
C the move to F90.

          CALL ZPRP_DYNAM(IBT,IDO,INP,%VAL(LGE_PTR),NAN,NBH,NBJ,NBJF,
     '      nc,ne,NFF,NGAP,NHE,NKEF,NKHE,NKJE,NNF,NPF,NPNE,NPNY,nr,NRE,
     '      NVHE,NVJE,NW,nx,NXI,NYNE,NYNP,CE,%VAL(CG_PTR),CGE,CP,
     '      CURVCORRECT,FEXT,PG,%VAL(RE_PTR),%VAL(RG_PTR),SE,WG,XA,
     '      %VAL(XE_PTR),%VAL(XG_PTR),XP,YG,YP,ZA,ZA1,ZAA,%VAL(ZE_PTR),
     '      %VAL(ZE1_PTR),%VAL(ZG_PTR),ZP,ZP1,ZPA,ERROR,*100)

C cpb 18/10/96 Free dynamically allocated arrays

          CALL FREE_MEMORY(LGE_PTR,ERROR,*100)
          CALL FREE_MEMORY(CG_PTR,ERROR,*100)
          CALL FREE_MEMORY(RE_PTR,ERROR,*100)
          CALL FREE_MEMORY(RG_PTR,ERROR,*100)
          CALL FREE_MEMORY(XE_PTR,ERROR,*100)
          CALL FREE_MEMORY(XG_PTR,ERROR,*100)
          CALL FREE_MEMORY(ZE_PTR,ERROR,*100)
          CALL FREE_MEMORY(ZE1_PTR,ERROR,*100)
          CALL FREE_MEMORY(ZG_PTR,ERROR,*100)

          GO TO 102
C           This statement is designed to be skipped if no error
C           occur. However if a error occurs within a subroutine
C           the alternate return points to line 100 to set the flag
 100        CONTINUE
C$OMP       CRITICAL(ZPRP_2)
            ERROR_FLAG=.TRUE.
            WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
            CALL WRITES(IOER,OP_STRING,ERROR,*101)
            WRITE(OP_STRING,'(/'' >>An error occurred - '
     '        //'results may be unreliable!'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101        CONTINUE
C$OMP       END CRITICAL(ZPRP_2)
 102      CONTINUE
        ENDIF !.NOT.ERROR_FLAG
      ENDDO !noelem (ne)
C$OMP END PARALLEL DO

      CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '  //'element residual calculations',ERROR,*9999)

C      IF(IWRIT4(nr,nx).GE.1) THEN
C        WRITE(OP_STRING,
C     '    '(/'' CPU time of 1 thread for element residuals calcs: '','
C     '    //'D11.4,'' s'')') TIME2-TIME1
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      ENDIF
      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
      IF(IWRIT4(nr,nx).GE.1) THEN
        WRITE(OP_STRING,
     '    '(/'' CPU time of 1 thread for element residuals calcs: '','
     '    //'D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(ITYP5(nr,nx).EQ.5) THEN !wavefront path analysis
        IF(ITYP15(nr,nx).NE.0) THEN !upwinding

          CALL ASSERT(NJ_LOC(NJL_GEOM,0,nr).EQ.3,
     '      '>>Only 3D upwinding is implemented',ERROR,*9999)
          CALL CPU_TIMER(CPU_USER,TIME_START)

          ERROR_FLAG=.FALSE.
C new MPN 1Feb2000: OMP parallel proc directive
C old
CC$DOACROSS local (nf,noface,IDOXFT_PTR,MYMS_PTR,NSFE_PTR,
CC$&               NYNS_PTR,CG_PTR,RDF_PTR,SM_PTR,SN_PTR,
CC$&               XDF_PTR,ZDF_PTR,ZE_PTR,ERROR)
CC$&        share (nc,NPF,nr,nx,ERROR_FLAG)
C end old
C$OMP     PARALLEL DO
C$OMP&      PRIVATE(nf,noface,IDOXFT_PTR,MYMS_PTR,NSFE_PTR,
C$OMP&              NYNS_PTR,CG_PTR,RDF_PTR,SM_PTR,SN_PTR,
C$OMP&              XDF_PTR,ZDF_PTR,ZE_PTR,ERROR),
C$OMP&      SHARED(MEM_INIT,nc,NHM,NPF,nr,nx,ERROR_FLAG)
          DO noface=1,NFFACE(0,nr) !is main face loop
            nf=NFFACE(noface,nr)
            IF(.NOT.ERROR_FLAG) THEN
C     '        .AND.(ITYP15(nr,nx).EQ.3.EQV.NPF(5,nf).EQ.2)) THEN
CC             correct number of adjacent elements.

C             Dynamic allocation of parallel local arrays
C             Intialise pointers so they are zero in the parallel loop
              IDOXFT_PTR=0
              CALL ALLOCATE_MEMORY(NHM,1,INTTYPE,IDOXFT_PTR,
     '          MEM_INIT,ERROR,*200)
              MYMS_PTR=0
              CALL ALLOCATE_MEMORY((NSFM+1)*2*2*NHM,1,INTTYPE,MYMS_PTR,
     '          MEM_INIT,ERROR,*200)
              NSFE_PTR=0
              CALL ALLOCATE_MEMORY(NSM*2*2*NHM,1,INTTYPE,NSFE_PTR,
     '          MEM_INIT,ERROR,*200)
              NYNS_PTR=0
              CALL ALLOCATE_MEMORY((NSM+1)*2*2*NHM,1,INTTYPE,NYNS_PTR,
     '          MEM_INIT,ERROR,*200)
              CG_PTR=0
              CALL ALLOCATE_MEMORY(NMM*NGM,1,DPTYPE,CG_PTR,
     '          MEM_INIT,ERROR,*200)
              RDF_PTR=0
              CALL ALLOCATE_MEMORY(NSFM*2*NHM,1,DPTYPE,RDF_PTR,
     '          MEM_INIT,ERROR,*200)
              SM_PTR=0
              CALL ALLOCATE_MEMORY(NSFM*2*2*NHM,1,DPTYPE,SM_PTR,
     '          MEM_INIT,ERROR,*200)
              SN_PTR=0
              CALL ALLOCATE_MEMORY(NSM*2*2*NHM,1,DPTYPE,SN_PTR,
     '          MEM_INIT,ERROR,*200)
              XDF_PTR=0
              CALL ALLOCATE_MEMORY(NSFM*2*2*NJM,1,DPTYPE,XDF_PTR,
     '          MEM_INIT,ERROR,*200)
              ZDF_PTR=0
              CALL ALLOCATE_MEMORY(NSFM*2*2*NHM,1,DPTYPE,ZDF_PTR,
     '          MEM_INIT,ERROR,*200)

C             Must use a dynam subroutine as dynamically allocated
C             arrays are used within this subroutine level.
              CALL ZPRP_FACE(IBT,IDO,%VAL(IDOXFT_PTR),INP,
     '          %VAL(MYMS_PTR),NBH,NBHF(1,1,nf),NBJ,NBJF(1,nf),nc,nf,
     '          NHE,NKB,NKHE,NKJE,NNB,NNF,NPF(1,nf),NPNE,NPNY,nr,
     '          NSB,%VAL(NSFE_PTR),NVHE,NVJE,nx,
     '          NYNP(1,1,1,1,0,1,nr),%VAL(NYNS_PTR),
     '          CE,%VAL(CG_PTR),CP,PG,%VAL(RDF_PTR),SE,
     '          %VAL(SM_PTR),%VAL(SN_PTR),%VAL(XDF_PTR),
     '          XP,YGF(1,1,nf),YP,%VAL(ZDF_PTR),ERROR,*200)

C             Free dynamically allocated arrays
              CALL FREE_MEMORY(IDOXFT_PTR,ERROR,*200)
              CALL FREE_MEMORY(MYMS_PTR,ERROR,*200)
              CALL FREE_MEMORY(NSFE_PTR,ERROR,*200)
              CALL FREE_MEMORY(NYNS_PTR,ERROR,*200)
              CALL FREE_MEMORY(CG_PTR,ERROR,*200)
              CALL FREE_MEMORY(RDF_PTR,ERROR,*200)
              CALL FREE_MEMORY(SM_PTR,ERROR,*200)
              CALL FREE_MEMORY(SN_PTR,ERROR,*200)
              CALL FREE_MEMORY(XDF_PTR,ERROR,*200)
              CALL FREE_MEMORY(ZDF_PTR,ERROR,*200)

              GO TO 202
C               This statement is designed to be skipped if no error
C               occurs.  However if a error occurs within a subroutine
C               the alternate return points to line 200 to set the flag
 200            CONTINUE
C$OMP           CRITICAL(ZPRP_3)
                ERROR_FLAG=.TRUE.
                WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                CALL WRITES(IOER,OP_STRING,ERROR,*201)
 201            CONTINUE
C$OMP           END CRITICAL(ZPRP_3)
                IF(IDOXFT_PTR.NE.0)
     '            CALL FREE_MEMORY(IDOXFT_PTR,ERROR,*202)
                IF(MYMS_PTR.NE.0) CALL FREE_MEMORY(MYMS_PTR,ERROR,*202)
                IF(NSFE_PTR.NE.0) CALL FREE_MEMORY(NSFE_PTR,ERROR,*202)
                IF(NYNS_PTR.NE.0) CALL FREE_MEMORY(NYNS_PTR,ERROR,*202)
                IF(CG_PTR.NE.0) CALL FREE_MEMORY(CG_PTR,ERROR,*202)
                IF(RDF_PTR.NE.0) CALL FREE_MEMORY(RDF_PTR,ERROR,*202)
                IF(SM_PTR.NE.0) CALL FREE_MEMORY(SM_PTR,ERROR,*202)
                IF(SN_PTR.NE.0) CALL FREE_MEMORY(SN_PTR,ERROR,*202)
                IF(XDF_PTR.NE.0) CALL FREE_MEMORY(XDF_PTR,ERROR,*202)
                IF(ZDF_PTR.NE.0) CALL FREE_MEMORY(ZDF_PTR,ERROR,*202)
 202          CONTINUE
            ENDIF !.NOT.ERROR_FLAG
          ENDDO !noface (nf)
C$OMP     END PARALLEL DO

          CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '      //'face residual calculations',ERROR,*9999)

          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
          IF(IWRIT4(nr,nx).GE.1) THEN
            WRITE(OP_STRING,
     '        '('' CPU time of 1 thread for face residuals calcs: '','
     '        //'3X,D11.4,'' s'')') ELAPSED_TIME
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

        ENDIF
      ENDIF

C *** Add in global loads (including radiation/convection b.c. if set)
      DO no_nynr=1,NYNR(0,1,1) !loop over rows
        ny=NYNR(no_nynr,1,1) !row number
        IF(NPNY(0,ny,0).EQ.1) THEN
          np=NPNY(4,ny,0)
C GMH 8/1/97 Update cmgui link
          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
        ENDIF
        ny1=GETNYR(2,NPNY,nr,0,1,ny,NYNE,NYNP) !is RHS var #
C MPN/RY 15Nov06: add RHS loads only for force/moment balance residuals
        IF(NPNY(3,ny,0).LE.NH_LOC(0,nx)) THEN !ny is a force/moment eqn
           YP(ny,4)=YP(ny,4)-YP(ny1,1)
        ENDIF
      ENDDO !no_nynr

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP   CRITICAL(ZPRP_4)
        WRITE(OP_STRING,'(/'' Global Residuals:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO no_nynr=1,NYNR(0,1,1)
          ny=NYNR(no_nynr,1,1)
          WRITE(OP_STRING,'('' RP('',I6,'')='',D13.5)') ny,YP(ny,4)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$OMP   END CRITICAL(ZPRP_4)
      ENDIF
 
C 25/06/07 XSL Shifted from NONLIN 
C NEWS 26/06/05 JHC For Cont Mech prob, modifiy residual vector with contact contributions
      IF(KTYP5G(nr).GT.0) THEN ! contact

C *** Modify global residual vector YP(ny,4) with contact contribution
        IF(.NOT.ERROR_FLAG) THEN
          CALL CPU_TIMER(CPU_USER,TIME_START)
          CALL CONTACT_RESIDUAL(IBT,IDO,INP,NBH,NBJ,NBJF,NFF,NHE,
     '      NKEF,NKHE,NKJE,NNF,NPF,NPNE,nr,NRE,NVHE,NVJE,NW,nx,NYNP,
     '      Z_CONT_LIST,CURVCORRECT,SE,YP,XA,XP,Z_CONT,ZA,ZE,ZP,FIX,
     '      ERROR,*1112)
C*** 19/03/08 JHC Removed below GO TO statement. If error occurs during CONTACT_RESIDUAL.f then
C                 error gets displayed within CONTACT_RESIDUAL.f rather than displaying on here.
C          GO TO 301
CC ***     This statement is designed to be skipped if no error
CC ***     occurs.  However if a error occurs within a subroutine
CC ***     the alternate return points to line 300 to set the flag
C
C 300      CONTINUE
C          ERROR_FLAG=.TRUE.
C          WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
C          CALL WRITES(IOER,OP_STRING,ERROR,*301)
C 301      CONTINUE
        ENDIF !.NOT.ERROR_FLAG
        CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '    //'contact residual calculations',ERROR,*1112)
        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
        IF(IWRIT4(nr,nx).GE.1) THEN
          WRITE(OP_STRING,
     '    '('' CPU time of 1 thread for global residual contact '
     '     //'calcs: '',D11.4,'' s'')') ELAPSED_TIME
          CALL WRITES(IODI,OP_STRING,ERROR,*1112)
        ENDIF   
      ENDIF ! KTYP5G
C NEWE JHC

      CALL EXITS('ZPRP')
      RETURN

C cpb 18/10/96 Free dynamically allocated arrays
 9999 IF(LGE_PTR.NE.0) CALL FREE_MEMORY(LGE_PTR,ERROR,*1112)
      IF(CG_PTR.NE.0) CALL FREE_MEMORY(CG_PTR,ERROR,*1112)
      IF(RE_PTR.NE.0) CALL FREE_MEMORY(RE_PTR,ERROR,*1112)
      IF(RG_PTR.NE.0) CALL FREE_MEMORY(RG_PTR,ERROR,*1112)
      IF(XE_PTR.NE.0) CALL FREE_MEMORY(XE_PTR,ERROR,*1112)
      IF(XG_PTR.NE.0) CALL FREE_MEMORY(XG_PTR,ERROR,*1112)
      IF(ZE_PTR.NE.0) CALL FREE_MEMORY(ZE_PTR,ERROR,*1112)
      IF(ZE1_PTR.NE.0) CALL FREE_MEMORY(ZE1_PTR,ERROR,*1112)
      IF(ZG_PTR.NE.0) CALL FREE_MEMORY(ZG_PTR,ERROR,*1112)

 1112 CALL ERRORS('ZPRP',ERROR)
      RET_ERROR=ERROR
      CALL EXITS('ZPRP')
      RETURN 1
      END


