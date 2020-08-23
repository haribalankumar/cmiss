      SUBROUTINE MARCH1(IBT,IDO,INP,ISC_GD,ISC_GK,ISC_GKK,ISC_GM,
     '  ISC_GQ,ISR_GD,ISR_GK,ISR_GKK,ISR_GM,ISR_GQ,LGE,NBH,NBJ,NEELEM,
     '  NENP,NHE,NHP,NHQ,NKJE,NKH,NKHE,NONY,NORD,NPF,NP_INTERFACE,NPNE,
     '  NPNODE,NPNY,NQNY,nr,NRE,NRLIST,NRLIST2,NVHE,NVHP,NVJE,NW,nx,
     '  NYNE,NYNO,NYNP,NYNR,NYQNR,CE,CG,CGE,CONY,CP,CYNO,ED,EM,ER,ES,GD,
     '  GK,GKK,GM,GQ,GR,GRR,PG,RG,SE,WG,XA,XE,XG,XO,XP,YG,YP,YQ,YQS,ZA,
     '  ZE,ZG,ZP,FIX,ERROR,*)

C#### Subroutine: MARCH1
C###  Description:
C###    MARCH1 performs time integration of linear or nonlinear
C###    (in terms other than transient) probs (+ nonlinear b.c.s) by
C###    linear,quadratic or cubic algorithm.

C**** For 2nd order problems initial accel.ns are calculated from d.e.
C**** with  known initial velocities and displacements (cubic only).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'b08.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'marc00.cmn'
      INCLUDE 'time01.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISC_GD(NISC_GDM),ISC_GK(NISC_GKM),ISC_GKK(NISC_GKKM),
     '  ISC_GM(NISC_GMM),ISC_GQ(NISC_GQM),ISR_GD(NISR_GDM),
     '  ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),ISR_GM(NISR_GMM),
     '  ISR_GQ(NISR_GQM),LGE(NHM*NSM,NRCM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),
     '  NHE(NEM),NHP(NPM),NHQ(NRM),NKJE(NKM,NNM,NJM,NEM),
     '  NKH(NHM,NPM,NCM),NKHE(NKM,NNM,NHM,NEM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM),NORD(5,NE_R_M),NPF(9,NFM),
     '  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),NQNY(2,NYQM,0:NRCM),
     '  nr,NRE(NEM),NRLIST(0:NRM),NRLIST2(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),NYQNR(0:NYQM,0:NRCM,NCM,0:NRM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM),
     '  CP(NMM,NPM),CYNO(0:NYOM,NOOPM,NRCM,0:NRM),ED(NHM*NSM,NHM*NSM),
     '  EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),GD(NZ_GD_M),
     '  GK(NZ_GK_M),GKK(NZ_GKK_M),GM(NZ_GM_M),GQ(NZ_GQ_M),GR(NYROWM),
     '  GRR(NOM),PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XO(NOM),
     '  XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),YP(NYM,NIYM),
     '  YQ(NYQM,NIQM,NAM),YQS(NIQSM,NQM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),
     '  ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER HIST_file_counter,INCR,IP,k,KFAC,KOUNT,na,NIQLIST(0:1),
     '  NIQSLIST(0:1),NIYLIST(0:16),no_nynr,no_nynr2,np,NSTEP,
     '  NUMTIMEDATA,ny
      REAL*8 ERR,OLD_DT,SNORM,SUM,T,TOL,YPMAX(16),YPMIN(16)
      CHARACTER FILEFORMAT*6
      LOGICAL CONTINUE,CONVERGED,DYNAM1,DYNAM2,ENDFILE,FIRST_A,LINEAR,
     '  OUTPUT,UPDATE_MATRIX,YPDATA,YQDATA,YQSDATA

      SAVE NSTEP
      CALL ENTERS('MARCH1',*9999)

      IP=1
      INCR=0
C DMAL 27-JUNE-2002 Allows DT to be defined by "fem solve ... delta_t"
C      DT=TINCR
      IF(ITYP6(nr,nx).EQ.1) THEN      !linear equations
        LINEAR=.TRUE.
      ELSE IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear equations or bcs
        LINEAR=.FALSE.
      ENDIF !ityp6
      FIRST_A=.TRUE.
      UPDATE_MATRIX=.TRUE.

      IF(BINTIMEFILE.GT.0) THEN
        FILEFORMAT='BINARY'
      ELSE
        FILEFORMAT='ASCII'
      ENDIF !bintimefile

C    ..Determine whether or not problem needs the GD and GM matrices.
C    ..If the problem needs GD, DYNAM1 should be set the true, false
C    ..otherwise. If the problem needs GM DYNAM2 should be set to true
C    ..false otherwise.

C!!!cpb 9/12/94 Only Navier-Stokes covered at the moment.

      IF(ITYP2(nr,nx).EQ.3.OR.ITYP2(nr,nx).EQ.5) THEN !Advec-Diffusion
        DYNAM1=.TRUE.  !Use GD                        !or Navier-Stokes
        DYNAM2=.FALSE. !Don't use GM
      ELSE
        DYNAM1=.TRUE.  !Use GD
        DYNAM2=.TRUE.  !Use GM
      ENDIF !ityp2

      IF(RESTART.AND.T_RESTART(nx).NE.0.0d0) THEN
        ERROR='>>Not implemented yet'
        T=T_RESTART(nx)
        FIRST_A=.TRUE.
        UPDATE_MATRIX=.TRUE.

      ELSE IF((.NOT.RESTART).OR.(T_RESTART(nx).EQ.0.0d0)) THEN !perform initial tasks
        T=TSTART
        NSTEP=0
        T_RESTART(nx)=0.0d0
        YPDATA=.TRUE.
        YQDATA=.FALSE.
        YQSDATA=.FALSE.

C       ..CPB 29/3/96 Just output niy=1 for the moment
        NIYLIST(0)=1
        NIYLIST(1)=1
        NIQLIST(0)=0
        NIQLIST(1)=0
        NIQSLIST(0)=0
        NRLIST(0)=1
        NRLIST(1)=nr
        NRLIST2(0)=1
        NRLIST2(1)=nr
        na=1

        IF(HIST_file_intervals.GT.0) THEN !history file output
          CALL IOHIST(IOFILE1,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,
     '      NQNY,NRLIST,NRLIST2,NUMTIMEDATA,nx,NYNR,NYQNR,T,YP,YPMAX,
     '      YPMIN,YQ,YQS,'WRITE',FILEFORMAT,FILE02,'OPEN',ENDFILE,
     '      .TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)
          HIST_file_counter=0 !initialise counter
        ENDIF !history file output

C       ..put i.c.s into current soln YP(ny,1) & previous soln YP(ny,8)
        DO no_nynr=1,NYNR(0,0,1,nr)
          ny=NYNR(no_nynr,0,1,nr)
          IF(NPNY(0,ny,0).EQ.1) THEN
            np=NPNY(4,ny,0)
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          ENDIF !npny
          IF(.NOT.FIX(ny,1)) YP(ny,1)=YP(ny,3)
          YP(ny,8)=YP(ny,1)
          YP(ny,9)=0.d0
        ENDDO !no_nynr

C       ..write initial solution
        IF(HIST_file_intervals.GT.0) THEN !history file output
          CALL IOHIST(IOFILE1,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,
     '      NQNY,NRLIST,NRLIST2,NUMTIMEDATA,nx,NYNR,NYQNR,T,YP,YPMAX,
     '      YPMIN,YQ,YQS,'WRITE',FILEFORMAT,FILE02,'TIME_DATA',
     '      ENDFILE,.TRUE.,YPDATA,
     '      YQDATA,YQSDATA,ERROR,*9999)
        ENDIF !history file output

      ENDIF !.not.restart


C------------------------ Begin main time loop -------------------------

      CONTINUE=.TRUE.
      DO WHILE(CONTINUE) !loop start
        NSTEP=NSTEP+1
C       ..solve for the next time
        T=T+DT

C       ..check if current time is < finish time & if so repeat soln
        IF(T.LE.TFINISH) THEN
          T_RESTART(nx)=T
C         ..(re)calc stiffness matrices if needed
          IF(UPDATE_MATRIX.AND.NSTEP.EQ.1.AND.KTYP3C(nx).EQ.0)THEN !First time so assemble
            CALL ASSEMBLE3(IBT,IDO,INP,ISC_GD,ISC_GK,ISC_GM,ISR_GD,
     '        ISR_GK,ISR_GM,LGE,NBH,NBJ,NEELEM,NHE,NKHE,NKJE,NORD,NPF,
     '        NPNE,nr,NVHE,NVJE,NW,nx,NYNE,NYNP,NYNR(0,0,1,nr),CE,CG,
     '        CGE,CP,ED,EM,ER,ES,GD,GK,GM,GR,PG,RG,SE,WG,XA,XE,XG,XP,
     '        YG,ZA,ZE,ZG,ZP,DYNAM1,DYNAM2,UPDATE_MATRIX,ERROR,*9999)
          ELSEIF((KTYP3C(nx).EQ.1).OR.REASSEMBLE) THEN
            CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NVHP,nx,NYNE,
     '        NYNP,YP,ZA,ZP,ERROR,*9999)
            CALL ASSEMBLE3(IBT,IDO,INP,ISC_GD,ISC_GK,ISC_GM,ISR_GD,
     '        ISR_GK,ISR_GM,LGE,NBH,NBJ,NEELEM,NHE,NKHE,NKJE,NORD,NPF,
     '        NPNE,nr,NVHE,NVJE,NW,nx,NYNE,NYNP,NYNR(0,0,1,nr),CE,CG,
     '        CGE,CP,ED,EM,ER,ES,GD,GK,GM,GR,PG,RG,SE,WG,XA,XE,XG,XP,
     '        YG,ZA,ZE,ZG,ZP,DYNAM1,DYNAM2,UPDATE_MATRIX,ERROR,*9999)
          ENDIF !update_matrix & nstep=1
C dpn 06/07/98 - adding else for when output needs to switched off
          IF(IWRIT1(nr,nx).EQ.0) THEN
            OUTPUT=.FALSE.
          ELSE
            IF(MOD(NSTEP,IWRIT1(nr,nx)).EQ.0) THEN
              OUTPUT=.TRUE.
            ELSE
              OUTPUT=.FALSE.
            ENDIF
          ENDIF !iwrit1
C*** Adjust any incremental boundary conditions
          DO no_nynr=1,NYNR(0,0,1,nr) !loop over global variables of GK
            ny=NYNR(no_nynr,0,1,nr) !global variable #
            IF(FIX(ny,2)) THEN !incremental bc
              IF(NPNY(0,ny,0).EQ.1) THEN
                np=NPNY(4,ny,0)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              ENDIF !npny
              YP(ny,1)=YP(ny,1)+YP(ny,2)
            ENDIF !fix
C cpb 1/4/96 Don't understand what this part is doing.
C          IF(FIX(ny,1)) THEN
C            IF(KTYP22.EQ.1) THEN
C              YP(ny,1)=(YP(ny,1)-YP(ny,8))/DT
C            ELSE IF(KTYP22.EQ.2) THEN
C              YP(ny,1)=(YP(ny,1)-YP(ny,8))/(DT*DT)
C            ELSE IF(KTYP22.EQ.3) THEN
C              YP(ny,1)=(YP(ny,1)-YP(ny,8))/(DT*DT*DT)
C            ENDIF
C          ENDIF
          ENDDO !no_nynr

C*** Calculate any initial accelerations etc.
          IF(IP.EQ.1) THEN
            INIT=.TRUE.
            IF(KTYP22.EQ.3) THEN !cubic time-stepping algorithm
              DO no_nynr=1,NYNR(0,0,1,nr) !Loop over the glob vars of GK
                ny=NYNR(no_nynr,0,1,nr) !is global variable #
                IF(NPNY(0,ny,0).EQ.1) THEN
                  np=NPNY(4,ny,0)
C GMH 8/1/97 Update cmgui link
                  CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                ENDIF !npny
                YP(ny,10)=-YP(ny,8) !mean displacement = -previous disp
                YP(ny,12)=-YP(ny,11) !mean velocity = -previous velocity
              ENDDO !no_nynr (ny)
              CALL SOLVE3(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,
     '          LGE,NBH,NENP,NHE,NONY(0,1,1,nr),NP_INTERFACE,NPNE,
     '          NPNY,nr,NRE,NVHE,nx,NYNE,NYNO(0,1,1,nr),NYNP,
     '          NYNR(0,0,1,nr),0,CONY(0,1,1,nr),CYNO(0,1,1,nr),GD,GK,
     '          GKK,GM,GQ,GR,GRR,XO,YP,DYNAM1,DYNAM2,FIRST_A,FIX,
     '          UPDATE_MATRIX,ERROR,*9999)
              DO no_nynr=1,NYNR(0,0,1,nr) !Loop over the glob vars of GK
                ny=NYNR(no_nynr,0,1,nr) !is global variable #
                IF(NPNY(0,ny,0).EQ.1) THEN
                  np=NPNY(4,ny,0)
C GMH 8/1/97 Update cmgui link
                  CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                ENDIF !npny
                YP(ny,14)=YP(ny,1)
              ENDDO !no_nynr (ny)
            ENDIF !ktyp22
          ENDIF !IP
          A1=DT*THETA(1)
          IF(KTYP22.GE.2) A2=DT**2/2.0d0*THETA(2)
          IF(KTYP22.EQ.3) A3=DT**3/6.0d0*THETA(3)
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(MARCH1_1)
            WRITE(OP_STRING,'(/'' A1='',D10.3,'' A2='',D10.3,'' A3='','
     '        //'D10.3)') A1,A2,A3
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(MARCH1_1)
          ENDIF !dop
          DO no_nynr=1,NYNR(0,0,1,nr)
            ny=NYNR(no_nynr,0,1,nr)
            IF(NPNY(0,ny,0).EQ.1) THEN
              np=NPNY(4,ny,0)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            ENDIF !npny
            IF(KTYP22.EQ.1) THEN
              YP(ny,8)=YP(ny,1) !store the previous solution
            ELSE IF(KTYP22.EQ.2) THEN
              YP(ny,9)=YP(ny,8)+A1*YP(ny,11)
              YP(ny,12)=YP(ny,11)
            ELSE IF(KTYP22.EQ.3) THEN
              YP(ny,9)=YP(ny,8)+A1*YP(ny,11)+A2*YP(ny,14)
              YP(ny,12)=YP(ny,11)+A1*YP(ny,14)
              YP(ny,15)=YP(ny,14)
            ENDIF !ktyp22
          ENDDO !no_nynr (ny)

C*** Solve the problem for the current time, iterating on the
C*** solution process if necessary.
          CONVERGED=.FALSE.
          KOUNT=0
          IF(.NOT.LINEAR) THEN
            DO no_nynr=1,NYNR(0,0,1,nr) !Loop over global vars of GK
              ny=NYNR(no_nynr,0,1,nr) !global variable #
              IF(NPNY(0,ny,0).EQ.1) THEN
                np=NPNY(4,ny,0)
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              ENDIF !npny
              YP(ny,5)=0.0d0
            ENDDO !no_nynr
          ENDIF !.not.linear
          DO WHILE(.NOT.CONVERGED)
            KOUNT=KOUNT+1
            CALL SOLVE3(ISC_GK,ISC_GKK,ISC_GQ,ISR_GK,ISR_GKK,ISR_GQ,
     '        LGE,NBH,NENP,NHE,NONY(0,1,1,nr),NP_INTERFACE,NPNE,NPNY,
     '        nr,NRE,NVHE,nx,NYNE,NYNO(0,1,1,nr),NYNP,
     '        NYNR(0,0,1,nr),KTYP22,CONY(0,1,1,nr),CYNO(0,1,1,nr),
     '        GD,GK,GKK,GM,GQ,GR,GRR,XO,YP,DYNAM1,DYNAM2,FIRST_A,
     '        FIX,UPDATE_MATRIX,ERROR,*9999)
            FIRST_A=.FALSE.
            IF(LINEAR) THEN
              CONVERGED=.TRUE.
            ELSE IF(.NOT.LINEAR) THEN
              SUM=0.0d0
              IF(KOUNT.GT.1) THEN !check for convergence
                DO no_nynr=1,NYNR(0,0,1,nr) !Loop over global vars of GK
                  ny=NYNR(no_nynr,0,1,nr) !global variable #
                  SUM=SUM+(YP(ny,1)-YP(ny,5))**2
                ENDDO
                IF(SUM.LT.1.0d-6) CONVERGED=.TRUE.
              ENDIF !kount
              WRITE(OP_STRING,'('' Nonlinear iteration: Kount='',i3,'
     '          //''' sum='',D12.3)') KOUNT,SUM
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C Note: YP(ny,1) is, on entry to solve3, an estimate of new time step
C solution and on exit is time derivative
              IF(.NOT.CONVERGED) THEN
                DO no_nynr=1,NYNR(0,0,1,nr) !Loop over global vars of GK
                  ny=NYNR(no_nynr,0,1,nr) !global variable #
                  IF(NPNY(0,ny,0).EQ.1) THEN
                    np=NPNY(4,ny,0)
C GMH 8/1/97 Update cmgui link
                    CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                  ENDIF !npny
                  YP(ny,5)=YP(ny,1) !temporary storage of slope
                  YP(ny,1)=YP(ny,8)+DT*YP(ny,1) !is new time sol est.
                ENDDO !no_nynr
              ENDIF !not.converged
            ENDIF !linear
          ENDDO !not.converged
C*** Calculate the solution at the current time from the increment plus
C*** the previous time
          DO no_nynr=1,NYNR(0,0,1,nr)
            ny=NYNR(no_nynr,0,1,nr)
            IF(NPNY(0,ny,0).EQ.1) THEN
              np=NPNY(4,ny,0)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            ENDIF !npny
            IF(KTYP23.EQ.2) YP(ny,5)=YP(ny,1)
            IF(KTYP22.EQ.1) THEN
              IF(.NOT.FIX(ny,1))THEN
                YP(ny,9)=YP(ny,1) !store the incremental solution
                YP(ny,1)=YP(ny,8)+DT*YP(ny,1) !previous + change
              ENDIF
            ELSE IF(KTYP22.EQ.2) THEN
              IF(FIX(ny,1)) THEN
                YP(ny,1)=YP(ny,8)+DT*DT*YP(ny,1)
              ELSE
                YP(ny,10)=YP(ny,11)+DT*YP(ny,1)
                YP(ny,1)=YP(ny,8)+DT*YP(ny,11)+DT*DT/2.0d0*YP(ny,1)
              ENDIF !fix
C DMAL 26-JUNE-02 I think this should be "(KTYP22.EQ.3)"
            ELSE IF(KTYP23.EQ.3) THEN
              IF(FIX(ny,1)) THEN
                YP(ny,1)=YP(ny,8)+DT*DT*DT*YP(ny,1)
              ELSE
                YP(ny,10)=YP(ny,11)+DT*YP(ny,14)+DT*DT/2.0d0*YP(ny,1)
                YP(ny,13)=YP(ny,14)+DT*YP(ny,1)
                YP(ny,1)=YP(ny,8)+DT*YP(ny,11)+DT*DT/2.0d0*YP(ny,14)+
     '            DT*DT*DT/6.0d0*YP(ny,1)
              ENDIF !fix
            ENDIF !ktyp22
          ENDDO !no_nynr
C*** Adjust any time varying parameters.

C         ..write history file output
          IF(HIST_file_intervals.GT.0) THEN !history file output
            HIST_file_counter=HIST_file_counter+1
            IF(HIST_file_counter.EQ.HIST_file_intervals) THEN
              CALL IOHIST(IOFILE1,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY,NQNY,NRLIST,NRLIST2,NUMTIMEDATA,nx,NYNR,NYQNR,
     '          T,YP,YPMAX,YPMIN,YQ,YQS,'WRITE',FILEFORMAT,FILE02,
     '          'TIME_DATA',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '          ERROR,*9999)
              HIST_file_counter=0
            ENDIF !counter
          ENDIF !history file output


C         ..write output
          IF(OUTPUT) THEN
            WRITE(OP_STRING,'(/'' Solution at time T+DT='',D11.4,'
     '        //''' with DT='',D11.4)') T,DT
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,NVHP,nx,NYNE,
     '        NYNP,YP,ZA,ZP,ERROR,*9999)
            CALL ZPOP(4,NBH,1,NEELEM,NHP,NKH,NPNODE,nr,NVHP,nx,NYNE,
     '        NYNP,YP(1,4),ZA,ZP,FIX(1,4),ERROR,*9999)
          ENDIF !output

C         ..adjust time stepping parameters if required
          OLD_DT=DT
          IF(KTYP23.EQ.2) THEN
            KFAC=1
            DO K=1,KTYP22+1
              KFAC=KFAC*K
            ENDDO !k
            ERR=0.0d0
            DO no_nynr=1,NYNR(0,0,1,nr) !Loop over glob variables of GK
              ny=NYNR(no_nynr,0,1,nr) !global variable #
              IF(NPNY(0,ny,0).EQ.1) THEN
                np=NPNY(4,ny,0)
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              ENDIF !npny
              YP(ny,10)=(YP(ny,5)-YP(ny,10))*DT**KTYP22/KFAC
              ERR=ERR+YP(ny,10)**2
            ENDDO !no_nynr
            ERR=DSQRT(ERR)
            IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP         CRITICAL(MARCH1_2)
              WRITE(OP_STRING,'(/,'' Estimate of error = '',D11.4)') ERR
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP         END CRITICAL(MARCH1_2)
            ENDIF !dop
            SNORM=0.0d0
            DO no_nynr=1,NYNR(0,0,1,nr) !Loop over glob variables of GK
              ny=NYNR(no_nynr,0,1,nr) !global variable #
              SNORM=SNORM+YP(ny,1)**2
            ENDDO !no_nynr
            SNORM=DSQRT(SNORM)
            TOL=SNORM
            IF(ERR.GT.TOL) THEN
              DT=DT/2.0d0
              INCR=0
            ELSE IF(ERR.LT.TOL/2.0d0) THEN
              INCR=INCR+1
              IF(INCR.GT.1) THEN
                DT=DT*1.250d0
                INCR=0
              ENDIF !incr
            ENDIF !err
          ENDIF

C         ..transfer current soln to prev soln and increment time
          DO no_nynr2=1,NYNR(0,0,1,nr) !loop over global variables of GK
            ny=NYNR(no_nynr2,0,1,nr) !global variable #
            IF(NPNY(0,ny,0).EQ.1) THEN
              np=NPNY(4,ny,0)
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            ENDIF !npny
            YP(ny,8)=YP(ny,1) !is solution vector at time T
            IF(KTYP22.EQ.2) THEN
              YP(ny,11)=YP(ny,10) !is velocity vector at time T
            ELSE IF(KTYP22.EQ.3) THEN
              YP(ny,11)=YP(ny,10) !is velocity vector at time T
              YP(ny,14)=YP(ny,13) !is acceleration vector at time T
            ENDIF !ktyp22
          ENDDO !no_nynr2
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP       CRITICAL(MARCH1_3)
            WRITE(OP_STRING,'('' TIMES : '',/,'' Current = '',D11.4,'
     '        //'/,'' New Delta = '',D11.4,'', Old Delta = '',D11.4,'
     '        //'/,'' Initial = '',D11.4,'', Final = '',D11.4)')
     '        T,DT,OLD_DT,TSTART,TFINISH
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP       END CRITICAL(MARCH1_3)
          ENDIF !dop
        ELSE
          CONTINUE=.FALSE.
        ENDIF !T
      ENDDO !time loop

C-------------------------- End main time loop ------------------------

      IF(HIST_file_intervals.GT.0) THEN !close history file
        CALL IOHIST(IOFILE1,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,
     '    NQNY,NRLIST,NRLIST2,NUMTIMEDATA,nx,NYNR,NYQNR,T,YP,YPMAX,
     '    YPMIN,YQ,YQS,'CLOSE',FILEFORMAT,FILE02,' ',ENDFILE,.TRUE.,
     '    YPDATA,YQDATA,YQSDATA,ERROR,*9999)
      ENDIF !history file output

      CALL CLOSEF(IOFILE2,ERROR,*9999)
      CALL CLOSEF(IOFILE4,ERROR,*9999)

      CALL EXITS('MARCH1')
      RETURN

 9999 IF(HIST_file_intervals.GT.0) THEN !close history file
        CALL IOHIST(IOFILE1,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,
     '    NQNY,NRLIST,NRLIST2,NUMTIMEDATA,nx,NYNR,NYQNR,T,YP,YPMAX,
     '    YPMIN,YQ,YQS,'CLOSE',FILEFORMAT,FILE02,' ',ENDFILE,.TRUE.,
     '    YPDATA,YQDATA,YQSDATA,ERROR,*1111)
      ENDIF !history file output
 1111 CALL ERRORS('MARCH1',ERROR)
      CALL EXITS('MARCH1')
      RETURN 1
      END


