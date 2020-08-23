      SUBROUTINE IPSOLV(IBT,IDO,INP,NAN,NAQ,NBH,NBJ,NEELEM,NELIST,
     '  NENP,NHE,NKB,NKHE,NKJE,NLL,NNB,NNF,NNL,NONY,
     '  NP_INTERFACE,NPF,NPL,NPNE,NPNODE,NPNY,nr,NRE,NVHE,NVHP,NVJE,NW,
     '  NWP,NWQ,nx,NXI,NXQ,NYNE,NYNO,NYNP,NYNR,NYNY,NYQNR,AQ,CONY,
     '  CYNO,CYNY,SE,SP,XA,XE,XP,YP,FIX,LINE_NAME,SEND_GLOBAL,FILE_NAME,
     '  ERROR,*)

C#### Subroutine: IPSOLV
C###  Description:
C###    IPSOLV defines solution parameters for region nr.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'adam00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'cell00.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'integrator.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp60.cmn'
      INCLUDE 'ktyp70.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lsoda00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'ptr00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NAN(NIM,NAM,NBFM),NAQ(NQM,NAM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NHE(NEM),NKB(2,2,2,NNM,NBFM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM),NP_INTERFACE(0:NPM,0:3),NPF(9,NFM),
     '  NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM),nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),
     '  NWP(NPM,2),NWQ(8,0:NQM,NAM),nx,NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),NYNY(0:NYYM,NYM,NRM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM)
      REAL*8 AQ(NMAQM,NQM),CONY(0:NOYM,NYM,NRCM,0:NRM),
     '  CYNO(0:NYOM,NOOPM,NRCM,0:NRM),CYNY(0:NYYM,NYM,NRM),
     '  SE(NSM,NBFM,NEM),SP(NKM,NBFM,NPM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),FILE_NAME*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM),LINE_NAME,SEND_GLOBAL
!     Local Variables
      INTEGER ICHAR,INFO,maqdt,n,nb,nh,nkk,NOQUES,nonr,
     '  no_nrlist,npp,nq,nrr,nvv,IBEG2,IEND2
      CHARACTER CHAR11*20,STRING*255
      LOGICAL FILEIP,NOCROSS
C MHT 27Jun03 unreferenced
C     INTEGER m,ne,ne0,YLIST(0:99)
c     REAL*8 DIST
C     LOGICAL FOUND      

      CALL ENTERS('IPSOLV',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      IF(ITYP5(nr,nx).EQ.2) THEN      !time integration
        CALL IPSOLT(nr,nx,FILE_NAME,LINE_NAME,ERROR,*9999)
      ELSE IF(ITYP5(nr,nx).EQ.3) THEN !modal analysis
        CALL IPEIGE(NYNR(0,0,1,nr),FIX,ERROR,*9999)

c cpb 1/5/95 Replacing Fourier analysis with Quasi-static analysis
c      ELSE IF(ITYP5(nr,nx).EQ.4) THEN !Fourier analysis
c        CALL IPFOUR(ERROR,*9999)

      ELSE IF(ITYP5(nr,nx).EQ.4) THEN !Quasi-static analysis
        CALL IPQUAS(nr,nx,ERROR,*9999)
      ELSE IF(ITYP5(nr,nx).EQ.6) THEN !buckling analysis
      ENDIF

      IF(ITYP4(nr,nx).LE.3) THEN      !fem or bem
        IF(ITYP6(nr,nx).EQ.1) THEN      !linear analysis
          IWRIT3(nr,nx)=1
        ELSE IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear analysis
          CALL IPNONL(nr,nx,SEND_GLOBAL,ERROR,*9999)
        ENDIF
C SGM 26 Oct 2000 grid-based Finite element also
C MLT 29Nov02 grid finite volume also
      ELSE IF(ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '        ITYP4(nr,nx).EQ.7) THEN !collocation or grid-based FE/FV
        IF(ITYP5(nr,nx).EQ.2.AND.(ITYP2(nr,nx).EQ.8.
     '    OR.ITYP2(nr,nx).EQ.9)) THEN
        ELSE  !Not activation modelling
          FORMAT='('' Specify type of solution procedure [1]:'''//
     '      '/''   (1) Direct solve'''//
     '      '/''   (2) Gauss-Seidel iterations'''//
     '      '/''   (3) Gauss-Seidel with multigrid accel.'''//
     '      '/''   (4) Unused'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=ITYP9(nr,nx)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ITYP9(nr,nx)=IDATA(1)
        ENDIF
      ENDIF !ITYP4

      IF(ITYP9(nr,nx).EQ.3) THEN !multigrid acceleration
        IDEFLT(1)=1
        FORMAT='($,'' Enter the number of multigrid levels [1]: '','
     '    //'I2)'
        IF(IOTYPE.EQ.3) IDATA(1)=NMGT
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,11,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NMGT=IDATA(1)
        CALL ASSERT(NAM.GE.NMGT,'>>NAM too small',ERROR,*9999)

        CALL ConstructNXQ(NAQ,nr,NWQ,NXQ,ERROR,*9999)
      ENDIF !ITYP9=3

      IF(ITYP5(nr,nx).EQ.2.OR.ITYP5(nr,nx).EQ.3.OR.ITYP5(nr,nx).EQ.6)
     '  THEN
        IF(ITYP2(nr,nx).NE.9) THEN
          IF (ITYP4(nr,nx).NE.3) THEN
            FORMAT='($,'' Do you want mass lumping [N]? '',A)'
            IF(IOTYPE.EQ.3) THEN
              IF(LUMP) THEN
                ADATA(1)='Y'
              ELSE
                ADATA(1)='N'
              ENDIF
            ENDIF
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              IF(ADATA(1).EQ.'Y') THEN
                LUMP=.TRUE.
              ELSE IF(ADATA(1).EQ.'N') THEN
                LUMP=.FALSE.
              ENDIF
            ENDIF
          ENDIF !ityp4
        ENDIF
      ENDIF !ityp5

      IF(ITYP5(nr,nx).LE.5) THEN !static,time,modal,quasi,front analy.
        IF(.NOT.IS_COUPLED(nx).OR.nr.EQ.COUP_NRLIST(1,nx)) THEN
          IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear
            IF(IOTYPE.EQ.3) THEN
              IF(KTYP14.GT.0) THEN
                ADATA(1)='Y'
              ELSE IF(KTYP14.EQ.0) THEN
                ADATA(1)='N'
              ENDIF
            ENDIF
            FORMAT='($,'' Do you want to increment a material'
     '        //' parameter [N]? '',A)'
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,
     '        ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              IF(ADATA(1).EQ.'Y') THEN
                KTYP14=1
                ERROR='>>Not Implemented'
                GO TO 9999
              ELSE IF(ADATA(1).EQ.'N') THEN
                KTYP14=0
              ENDIF
            ENDIF
          ELSE
            KTYP14=0
          ENDIF
        ENDIF

CMLB 15-Oct-1998 IWRIT4 is set inside ipsolu
        IF(ITYP5(nr,nx).EQ.2.AND.(ITYP2(nr,nx).EQ.8
     '    .OR.ITYP2(nr,nx).EQ.9)) THEN
C         threshold cardiac activation
C          FORMAT='('' Specify option for output [0]: '''//
C     '      '/''   (0) No output'''//
C     '      '/''   (1) File o/p only'''//
C     '      '/''   (2) & terminal'''//
C     '      '/$,''    '',I1)'
C          IF(IOTYPE.EQ.3) IDATA(1)=IWRIT4(nr,nx)
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,2,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) IWRIT4(nr,nx)=IDATA(1)
C          IF(ITYP16(nr,nx).EQ.2) CALL IPSOLU(nr,nx,ERROR,*9999)
C          IF(ITYP16(nr,nx).EQ.3) CALL IPSOLU(nr,nx,ERROR,*9999)

C SGM 22Jan01 also call IPSOLU for fully explicit Grid-based FE
C MLT 29Nov02 added grid finite volume
          IF(THETA(1).GT.ZERO_TOL.OR.ITYP4(nr,nx).EQ.6
     '       .OR.ITYP4(nr,nx).EQ.7) THEN
            CALL IPSOLU(nr,nx,ERROR,*9999)
          ELSEIF(KTYP32.EQ.2) THEN
            CALL IPSOLU(nr,nx,ERROR,*9999)
          ENDIF
        ELSE
C News AJP 23-2-95
          IF(.NOT.IS_COUPLED(nx).OR.nr.EQ.COUP_NRLIST(1,nx)) THEN
C           Ask this question if not a coupled problem or for the
C           first regions only.  Assumes that the first region has
C           this question asked.
            CALL IPSOLU(nr,nx,ERROR,*9999)
          ENDIF !not coupled or first region
        ENDIF
      ENDIF !ityp5

      IF(ITYP5(nr,nx).EQ.2.OR.ITYP5(nr,nx).EQ.4.OR.ITYP6(nr,nx).EQ.2)
     '  THEN
        IF(ITYP2(nr,nx).NE.9) THEN
          IF(.NOT.IS_COUPLED(nx).OR.nr.EQ.COUP_NRLIST(1,nx)) THEN
            FORMAT='($,'' Enter the solution output frequency (0 for '
     '        //'no output) [1]: '',I5)'
            IF(IOTYPE.EQ.3) IDATA(1)=IWRIT1(nr,nx)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,
     &        10000,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) IWRIT1(nr,nx)=IDATA(1)
C cpb 2/2/96
C          FORMAT='($,'' Do you want computations to await prompt '
C     '      //'after every Nth step [N]? '',A)'
C          IF(IOTYPE.EQ.3) THEN
C            IF(PROMPT) THEN
C              ADATA(1)='Y'
C            ELSE
C              ADATA(1)='N'
C            ENDIF
C          ENDIF
C          CALL GINOUT(IOTYPE,1,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) THEN
C            IF(ADATA(1).EQ.'Y') THEN
C              PROMPT=.TRUE.
C            ELSE IF(ADATA(1).EQ.'N') THEN
C              PROMPT=.FALSE.
C            ENDIF
C          ENDIF
          ENDIF
          IF(ITYP2(nr,nx).EQ.11)THEN !pulmonary
            FIRST_SOLVE(nx)=.TRUE.
            IF(ITYP3(nr,nx).LE.2)THEN !inert gas or thermofluids
              IWRIT6=0
              FORMAT='($,'' Enter the timing output frequency (0 for '
     '          //'no output) [1]: '',I5)'
              IF(IOTYPE.EQ.3) IDATA(1)=IWRIT6
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,
     &          100,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) IWRIT6=IDATA(1)
              IF(IWRIT6.GT.0)THEN
                IDEFLT(1)=0
                IF(IOTYPE.EQ.3) IDATA(1)=0
                FORMAT=
     &            '($,'' Enter nodes (groups) for output (0 for '//
     &            'no output) [0]: '',I5)'
                CDATA(1)='NODES'
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,
     &            NP_R_M,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &            *9999)
                IF(IOTYPE.NE.3) THEN
                  NP_OP_LIST(0)=IDATA(0)
                  IF(ITYP7(nr,nx).GT.1)THEN !gas exchange
C Restrict size of output nodes that are written to file *.nodes
                    WRITE(STRING,'(''>>Max size of output nodes '//
     &                '(NP_OP_LIST) is 10 '')')
                    CALL ASSERT(NP_OP_LIST(0).LE.10,STRING,ERROR,*9999)
                  ENDIF
                  IF(IDATA(0).GT.0.AND.IDATA(1).EQ.0) NP_OP_LIST(0)=0
                  DO N=1,NP_OP_LIST(0)
c                    WRITE(STRING,'(''>>Increase NP_OP_LIST, try '',I6)')
c     &                NP_OP_LIST(0)
c                    CALL ASSERT(N.LE.,STRING,ERROR,*9999)
                    NP_OP_LIST(n)=IDATA(n)
                  ENDDO !N
                ENDIF !iotype.ne.3
              ENDIF
              IF(ITYP3(nr,nx).EQ.1)THEN !gas mixing
                IF(ITYP7(nr,nx).GT.1)THEN !gas exchange
                FORMAT='('' Specify option for output [0]:'''//
     '            '/''   (0) No output'''//
     '            '/''   (1) Alveolar and arterial PO2 (.history)'''//
     '            '/''   (2) & Mass error (.error)'''//
     '            '/''   (3) & Specified node history (.nodes)'''//
     '            '/$,''    '',I1)'
                ELSE
                FORMAT='('' Specify option for output [0]:'''//
     '            '/''   (0) No output'''//
     '            '/''   (1) Washout history (.history)'''//
     '            '/''   (2) & Normalised alveolar slope (.snhist)'''//
     '            '/''   (3) & Parameterisation for LPM (.acinus)'''//
     '            '/$,''    '',I1)'
                ENDIF
              ELSE IF(ITYP3(nr,nx).EQ.2)THEN !water-heat transfer
                FORMAT='('' Specify option for output [0]:'''//
     '            '/''   (0) No output'''//
     '            '/''   (1) Temperature history (.history)'''//
     '            '/''   (2) & Path profile end-phase (.end_phase)'''//
     '            '/''   (3) & Path profile for all (.all_end_phase)'''
     &            //'/$,''    '',I1)'
              ENDIF
              IDEFLT(1)=1
              IF(IOTYPE.EQ.3) IDATA(1)=LUNG_OP
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          3,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) LUNG_OP=IDATA(1)
            ENDIF
          ENDIF
        ELSEIF (ITYP2(nr,nx).EQ.9) THEN
          IDEFLT(1)=0
          FORMAT='('' Specify timing output for time integration'//
     '      ' [0]:'''//
     '      '/''   (0) No output  '''//
     '      '/''   (1) Simple     '''//
     '      '/''   (2) Verbose    '''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=IWRIT5(nr,nx)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,2,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) IWRIT5(nr,nx)=IDATA(1)
          IF(IWRIT5(nr,nx).EQ.2) THEN
          FORMAT='($,'' Enter the timing output frequency [1]: '',I5)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=IWRIT1(nr,nx)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        10000,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) IWRIT1(nr,nx)=IDATA(1)
          ENDIF
        ENDIF
      ENDIF

      IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.9) THEN !cellular modelling
         IF(KTYP23.EQ.1) THEN
            FORMAT='('' Specify type of integration procedure [1]:'''//
     '           '/''   (1) Euler'''//
     '           '/''   (2) Improved Euler'''//
     '           '/''   (3) Runge-Kutta (4th order)'''//
     '           '/''   (4) '''//
     '           '/''   (5) Adams-Moulton (variable order, adaptive '
     '           //'time step)'''//
     '           '/''   (6) LSODA (adaptive step, stiff switching)'''//
     '           '/$,''    '',I1)'
         ELSE IF(KTYP23.EQ.2) THEN
            FORMAT='('' Specify type of integration procedure [1]:'''//
     '           '/''   (1) Euler'''//
     '           '/''   (2) Improved Euler'''//
     '           '/''   (3) Runge-Kutta (4th order)'''//
     '           '/''   (4) '''//
     '           '/''   (5) Adams-Moulton (variable order, adaptive '
     '           //'time step)'''//
     '           '/''  *(6) LSODA (adaptive step, stiff switching)'''//
     '           '/$,''    '',I1)'
         ENDIF

        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP37
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,6,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP37=IDATA(1)
        IF((KTYP37.EQ.6).AND.(KTYP23.EQ.2)) THEN
           ERROR='Automatic stepping for LSODA not implemented'
           GOTO 9999
        ENDIF
        IF(KTYP37.EQ.5) THEN
C         Adams-Moulton. Get additional parameters and allocate work
C         arrays.
          FORMAT='($,'' Enter the maximum Adams polynomial order '
     '      //'[4]: '',I2)'
          IDEFLT(1)=4
          IF(IOTYPE.EQ.3) IDATA(1)=ADAMS_MAX_ORDER
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,12,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ADAMS_MAX_ORDER=IDATA(1)
          FORMAT='($,'' Enter the maximum Adams step size [0.100]:'
     '      //' '',D11.4)'
          RDEFLT(1)=0.1d0
          IF(IOTYPE.EQ.3) RDATA(1)=ADAMS_MAX_STEP
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ADAMS_MAX_STEP=RDATA(1)
          FORMAT='($,'' Enter the maximum number of Adams iterations '
     '      //'[100]: '',I3)'
          IDEFLT(1)=100
          IF(IOTYPE.EQ.3) IDATA(1)=ADAMS_MAX_ITERS
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,999,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ADAMS_MAX_ITERS=IDATA(1)
          FORMAT='('' Specify type of error control [1]:'''//
     '      '/''   (1) Pure absolute'''//
     '      '/''   (2) Relative to Y'''//
     '      '/''   (3) Relative to DY'''//
     '      '/''   (4) Mixed relative/absolute'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=4
          IF(IOTYPE.EQ.3) IDATA(1)=ADAMS_ERROR_TYPE
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,5,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ADAMS_ERROR_TYPE=IDATA(1)
C ***     DPN 18 May 1999 - you need to define this for all cases
c ***          IF(ADAMS_ERROR_TYPE.EQ.4) THEN
          FORMAT='($,'' Enter the absolute error component [0.05]:'
     '      //' '',D11.4)'
          RDEFLT(1)=0.05d0
          IF(IOTYPE.EQ.3) RDATA(1)=ADAMS_ABS_ERR
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ADAMS_ABS_ERR=RDATA(1)
          FORMAT='($,'' Enter the relative error component [0.05]:'
     '      //' '',D11.4)'
          RDEFLT(1)=0.05d0
          IF(IOTYPE.EQ.3) RDATA(1)=ADAMS_REL_ERR
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ADAMS_REL_ERR=RDATA(1)
c ***          ENDIF
          FORMAT='($,'' Use rounding control [N]? '',A)'
          IF(IOTYPE.EQ.3) THEN
            IF(ADAMS_USE_ROUND_CTRL) THEN
              ADATA(1)='Y'
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ADAMS_USE_ROUND_CTRL=ADATA(1).EQ.'Y'

C cpb 8/2/99 Should the memory allocation be done here or later?
C            ! Get YQ-Y mappings
C            ! DPN 17 May 1999 - only need this for MB's stuff ??
C            CALL GETYQYMAP(nr,nx,YLIST,ERROR,*9999)
C            ADAMS_LIWORK=7
C            ADAMS_LWORK=7+7*ADAMS_MAX_ORDER+YLIST(0)*(ADAMS_MAX_ORDER+6)
C            IF(ADAMS_USE_ROUND_CTRL) ADAMS_LWORK=ADAMS_LWORK+YLIST(0)*2
C            CALL ALLOCATE_MEMORY(ADAMS_LIWORK*NQM*NXM,1,INTTYPE,
C     '        ADAMS_IWORK_PTR,MEM_INIT,ERROR,*9999)
C            CALL ALLOCATE_MEMORY(ADAMS_LWORK*NQM*NXM,1,DPTYPE,
C     '        ADAMS_WORK_PTR,MEM_INIT,ERROR,*9999)


C            ! Get YQ-Y mappings
C            ! DPN 17 May 1999 - only need this for MB's stuff ??
C            CALL GETYQYMAP(nr,nx,YLIST,ERROR,*9999)
          INTEGRATOR_LIWORK=7
          INTEGRATOR_LWORK=7+7*ADAMS_MAX_ORDER+NIQST*
     '      (ADAMS_MAX_ORDER+6)
          IF(ADAMS_USE_ROUND_CTRL) INTEGRATOR_LWORK=INTEGRATOR_LWORK+
     '      NIQST*2

        ELSEIF(KTYP37.EQ.6) THEN !LSODA
C rgb 29th Nov 2000
C         LSODA Defaults
C         There should be some more questions
C         here to define these parameters
          LSODA_JACOBIAN_TYPE=2
          LSODA_ERROR_TYPE=1
          INTEGRATOR_LWORK=22+NIQST*MAX(16,NIQST+9)+
     '      LSODA_NUM_REAL_COMMON
          INTEGRATOR_LIWORK=20+NIQST+LSODA_NUM_INT_COMMON
          LSODA_ITASK=1
          LSODA_IOPT=1
C         Lsoda questions
          FORMAT='($,'' Enter the maximum LSODA step size [0.100]:'
     '      //' '',D11.4)'
          RDEFLT(1)=0.1d0
          IF(IOTYPE.EQ.3) RDATA(1)=LSODA_MAX_STEP
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) LSODA_MAX_STEP=RDATA(1)
          FORMAT='($,'' Enter the maximum number of LSODA iterations '
     '      //'[100]: '',I3)'
          IDEFLT(1)=100
          IF(IOTYPE.EQ.3) IDATA(1)=LSODA_MAX_ITERS
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,999,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) LSODA_MAX_ITERS=IDATA(1)
          FORMAT='('' Specify type of error control [1]:'''//
     '      '/''   (1) Pure absolute'''//
     '      '/''   (2) Relative to Y'''//
     '      '/''   (3) Mixed relative/absolute'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=LSODA_ERROR_CONTROL
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,5,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) LSODA_ERROR_CONTROL=IDATA(1)
          IF(LSODA_ERROR_CONTROL.EQ.1.OR.LSODA_ERROR_CONTROL.EQ.3) THEN
             FORMAT='($,'' Enter the absolute error component [0.005]:'
     '            //' '',D11.4)'
             RDEFLT(1)=0.005d0
             IF(IOTYPE.EQ.3) RDATA(1)=LSODA_ABS_ERR
             CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &         FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &         IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &         *9999)
             IF(IOTYPE.NE.3) LSODA_ABS_ERR=RDATA(1)
          ENDIF
          IF(LSODA_ERROR_CONTROL.EQ.2.OR.LSODA_ERROR_CONTROL.EQ.3) THEN
             FORMAT='($,'' Enter the relative error component [0.005]:'
     '            //' '',D11.4)'
             RDEFLT(1)=0.005d0
             IF(IOTYPE.EQ.3) RDATA(1)=LSODA_REL_ERR
             CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &         FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &         IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &         *9999)
             IF(IOTYPE.NE.3) LSODA_REL_ERR=RDATA(1)
          ENDIF
          IF(LSODA_ERROR_CONTROL.EQ.1) LSODA_REL_ERR = 0.0d0
          IF(LSODA_ERROR_CONTROL.EQ.2) LSODA_ABS_ERR = 0.0d0
        ELSE
C set up the minimal integration workspace sizes
          INTEGRATOR_LIWORK=1
          INTEGRATOR_LWORK=1
        ENDIF

C ***   DPN 2001-10-27 Only need one allocation block for the workspace
        IF(INTEGRATOR_IWORK_PTR.NE.0) THEN
          CALL FREE_MEMORY(INTEGRATOR_IWORK_PTR,ERROR,*9999)
          INTEGRATOR_IWORK_PTR=0
        ENDIF
C MLT 11May05 Trying to overcome a 64 bit limitation
C on memory allocation
C Old memory allocation
C        CALL ALLOCATE_MEMORY(INTEGRATOR_LIWORK*NQM,1,INTTYPE,
C     '    INTEGRATOR_IWORK_PTR,MEM_INIT,ERROR,*9999)
C New allocation that splits request into two integer parts
        CALL BIG_ALLOCATE_MEMORY(INTEGRATOR_LIWORK,NQM,1,INTTYPE,
     '    INTEGRATOR_IWORK_PTR,MEM_INIT,ERROR,*9999)
        IF(INTEGRATOR_WORK_PTR.NE.0) THEN
          CALL FREE_MEMORY(INTEGRATOR_WORK_PTR,ERROR,*9999)
          INTEGRATOR_WORK_PTR=0
        ENDIF
C MLT 11May05 Trying to overcome a 64 bit limitation
C on memory allocation
C Old memory allocation
C        CALL ALLOCATE_MEMORY(INTEGRATOR_LWORK*NQM,1,DPTYPE,
C     '    INTEGRATOR_WORK_PTR,MEM_INIT,ERROR,*9999)
C New allocation that splits request into two integer parts
        CALL BIG_ALLOCATE_MEMORY(INTEGRATOR_LWORK,NQM,1,DPTYPE,
     '    INTEGRATOR_WORK_PTR,MEM_INIT,ERROR,*9999)

        FORMAT='($,'' Use DTAR (Dynamic Tracking of Active Region)'
     '    //' [N]? '',A)'
        IF(IOTYPE.EQ.3) THEN
          IF(KTYP36.EQ.1) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(ADATA(1).EQ.'Y') THEN
            KTYP36=1
          ELSE
            KTYP36=0
          ENDIF
        ENDIF

        CALL ASSERT(NMAQM.GE.9,'>>Increase NMAQM (9)',ERROR,*9999)
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,maqdt,MAQ_DT,ERROR,*9999)
        IF(maqdt.EQ.0) THEN
          CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_TIME,maqdt,MAQ_DT,
     '      ERROR,*9999)
        ENDIF

        IF(KTYP23.EQ.2) THEN !automatic time stepping
          DO nq=1,NQT
            AQ(maqdt,nq)=TINCR/10.0d0
          ENDDO
        ENDIF

      ENDIF !ityp5 & ityp2

      IF(ITYP2(nr,nx).EQ.1.AND.ITYP6(nr,nx).EQ.2) THEN
        FORMAT='($,'' Are element pressures read from file [N]? '',A)'
        IF(IOTYPE.EQ.3) ADATA(1)='N'
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(ADATA(1).EQ.'Y') THEN
          KTYP71=1
          FORMAT='($,'' Enter file name [FILE]'
     '      //' (file extension is .PRESSURE): '',A)'
          CDEFLT(1)='FILE'
          IF(IOTYPE.EQ.3) CDATA(1)=FILE01(1:20)
          IF(IOTYPE.NE.3) FILE01(1:100)=CDATA(1)(1:100)
        ELSE
          KTYP71=0
        ENDIF
      ELSE
        KTYP71=0
      ENDIF !ityp2 & ityp6

      IF((ITYP4(nr,nx).LE.2.AND.ITYP5(nr,nx).NE.3).OR.
     '  (ITYP4(nr,nx).EQ.5)) THEN
C       FEM or BEM or FV solution and not modal
        IF(.NOT.IS_COUPLED(nx).OR.nr.EQ.COUP_NRLIST(1,nx)) THEN
          FORMAT='('' Specify option for solution matrices output '
     '      //'[0]: '''//
     '      '/''   (0) No output'''//
     '      '/''   (1) Solution matrix'''//
     '      '/''   (2) RHS vector'''//
     '      '/''   (3) Solution matrix & RHS vector'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=KTYP4
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) KTYP4=IDATA(1)
          IF(KTYP4.NE.0) THEN
            FORMAT='($,'' Are the matrices to be stored in a binary '
     '        //'file [N]? '',A)'
            IF(IOTYPE.EQ.3) THEN
              IF(KTYP4.LT.0) THEN
                ADATA(1)='Y'
              ELSE
                ADATA(1)='N'
              ENDIF
            ENDIF
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              IF(ADATA(1).EQ.'Y'.OR.ADATA(1).EQ.'y') THEN
                KTYP4=-1*KTYP4 !reverse sign of ktyp4 for binary
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF !ityp4


C SGM 26 Oct 2000 grid-based Finite element also
C MLT 29Nov02 grid finite volume also
      IF((ITYP2(nr,nx).EQ.3).AND.(ITYP3(nr,nx).EQ.1).AND.
     '  (ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6
     '     .OR.ITYP4(nr,nx).EQ.7).AND.
     '  (ITYP5(nr,nx).EQ.1).AND.
     '  (ITYP6(nr,nx).EQ.1).AND.(ITYP9(nr,nx).EQ.1)) THEN
         FORMAT='($,'' Do you wish to use flux boundary conditions '
     '    //' [N]? '',A)'
        IF(IOTYPE.EQ.3) THEN
          IF(SOLVE8_FLUXBC) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(ADATA(1).EQ.'Y'.OR.ADATA(1).EQ.'y') THEN
            SOLVE8_FLUXBC=.TRUE.
          ELSE
            SOLVE8_FLUXBC=.FALSE.
          ENDIF
        ENDIF
      ENDIF

C cpb 17/5/97 Adding Salu consistency from "Implementing a consistency
C criterion in numerical solution of the bioelectric forward problem",
C by Yehuda Salu, IEEETBME, 26(6):338-341, 1980.
C cpb 20/8/97 Modifying the Salu consistency criterion to allow for
C the fixing of a potential at an arbitary value. The formula for
C alpha now becomes
C
C alpha=(sum_{j=2}^{n}A_1j.phi_j^*-P_1+A_11.phi_1)/
C       (1-sum_{j=2}^{n}A_1j.phi_j^1)
C
      IF((ITYP5(nr,nx).EQ.1.OR.ITYP5(nr,nx).EQ.4).AND.
     '  (ITYP2(nr,nx).EQ.3.OR.
     '  (ITYP2(nr,nx).EQ.5.AND.ITYP3(nr,nx).EQ.2)).AND.
     '  (ITYP4(nr,nx).NE.4.AND.ITYP4(nr,nx).NE.6.AND.
     '   ITYP4(nr,nx).NE.7)) THEN 
C SGM 30 Oct 2000 ignore grid-based Finite element also
C MLT 10Dec02 ignore grid FV also
C       !Static or Quasi-static Laplace or Poisson (with special source)
        IF(.NOT.IS_COUPLED(nx).OR.nr.EQ.COUP_NRLIST(1,nx)) THEN
          FORMAT='($,'' Do you wish to use the Salu consistency '
     '      //'criterion [N]? '',A)'
          IF(IOTYPE.EQ.3) THEN
            IF(SALU_CONSISTENCY(nx)) THEN
              ADATA(1)='Y'
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(ADATA(1).EQ.'Y'.OR.ADATA(1).EQ.'y') THEN
              SALU_CONSISTENCY(nx)=.TRUE.
            ELSE
              SALU_CONSISTENCY(nx)=.FALSE.
            ENDIF
          ENDIF
          IF(SALU_CONSISTENCY(nx)) THEN
            FORMAT='($,'' Enter the region number of the fixed '
     '        //'potential [1]: '',I5)'
            IF(IOTYPE.EQ.3) IDATA(1)=NPNY(6,SALU_NY(nx),0)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NRM,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) nrr=IDATA(1)
            FORMAT='($,'' Enter the node number of the fixed potential '
     '        //'[1]: '',I5)'
C LKC 15-JAN-1999 - Set CDATA for grouped input
            CDATA(1)='NODES'
            IF(IOTYPE.EQ.3) IDATA(1)=NPNY(4,SALU_NY(nx),0)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NPM,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) npp=IDATA(1)
            FORMAT='($,'' Enter the version number of the fixed '
     '        //'potential [1]: '',I5)'
            IF(IOTYPE.EQ.3) IDATA(1)=NPNY(2,SALU_NY(nx),0)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NVM,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) nvv=IDATA(1)
            FORMAT='($,'' Enter the derivative number of the fixed '
     '        //'potential [1]: '',I5)'
            IF(IOTYPE.EQ.3) IDATA(1)=NPNY(1,SALU_NY(nx),0)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NKM,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              nkk=IDATA(1)
              nh=NH_LOC(1,nx)
              SALU_NY(nx)=NYNP(nkk,nvv,nh,npp,0,1,nrr)
              CALL ASSERT(FIX(SALU_NY(nx),1),
     '          '>>Fixed potential is not set as a boundary condition',
     '          ERROR,*9999)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      IF(ITYP4(nr,nx).EQ.2) THEN !BE solution
!News AJP 23-2-95
        nonr=1 !order of evaluation of the next test is uncertain
        IF(IS_COUPLED(nx)) THEN
          DO WHILE(ITYP4(COUP_NRLIST(nonr,nx),nx).NE.2) !not a be region
            nonr=nonr+1
            CALL ASSERT(nonr.LT.NRT,' Error in finding first BE region',
     '        ERROR,*9999)
          ENDDO
          !nonr is the first BE region
        ENDIF !is_coupled
        IF(.NOT.IS_COUPLED(nx).OR.nr.EQ.COUP_NRLIST(nonr,nx)) THEN
          !not coupled or nr is first be region
!newe
c cpb 10/5/95 This is now more general
C          FORMAT='($,'' Do you want the global matrices written out'
C     '      //' to a file [N]? '',A)'
C          IF(IOTYPE.EQ.3) ADATA(1)='N'
C          CALL GINOUT(IOTYPE,1,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF((ADATA(1).EQ.'Y').OR.(ADATA(1).EQ.'y'))THEN
C            KTYP4=1
C          ELSE
C            KTYP4=0
C          ENDIF
C          !New AJP 2-3-94
          CALL EQTYPE(IBT,NBH,NEELEM,nr,NW,nx,ERROR,*9999)
          !To see if hermite interpolation is used.
          IF(HERMITE)THEN
!           !at least one element containing Hermite interpolation
            FORMAT='($,'' Do you want to use the hypersingular '//
     '        'code [Y]? '',A)'
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(ADATA(1).EQ.'Y') THEN
C KAT 13Jan00: Cross derivatives may already be eliminated in ipbase.
              NOCROSS=.FALSE.
              DO nb=1,NBFT
                IF(NBCD(nb).EQ.1) NOCROSS=.TRUE.
              ENDDO
              HYP=.TRUE.
C              KTYP92=1 !Default value
              IF(KTYP93(1,nr).EQ.0.AND..NOT.NOCROSS) THEN
                IF(NJT.EQ.2)THEN !2d
                  FORMAT='(/'' Specify whether to use [1]:         '''//
     '              '/''   (1) BIE plus tangential derivative BIE  '''//
     '              '/''   (2) BIE plus normal derivative BIE      '''//
     '              '/''   (3) Normal and tangential derivative BIE'''//
     '              '/$,''    '',I1)'
                  IDEFLT(1)=1
                ELSE !IF (KTYP93(1,nr).EQ.0) THEN !3d and full interpolation
                  FORMAT='(/'' Specify whether to use [2]:         '''//
     '              '/'//
     '              '   (1) BIE plus s1,s2 and s1s2 derivative BIEs'''
     '              //'/''   (2) BIE plus s1,s2 and n derivative BIEs'''
     '              //'/''   (3) s1,s2,s1s2 and n derivative BIEs'''//
     '              '/$,''    '',I1)'
                  IDEFLT(1)=2
                ENDIF
C               IF(KTYP93(1,nr).EQ.0)THEN
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,1,3,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '            INFO,ERROR,*9999)
                KTYP92=IDATA(1)
              ELSEIF(NJT.EQ.2)THEN !2d
                KTYP92=1
              ELSE
                KTYP92=0 !Reduced set of derivative equations in 3d
              ENDIF
            ENDIF
          ENDIF
          FORMAT='($,'' Do you want to use the adaptive integration '
     '      //'[Y]? '',A)'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(ADATA(1).EQ.'Y') THEN
            ADAPINT=.TRUE.
          ELSE
            ADAPINT=.FALSE.
          ENDIF
C          FORMAT='($,'' Do you want the singular values of the '//
C     '    'solution matrix (no solution) [N]? '',A)'
C          CALL GINOUT(IOTYPE,1,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(ADATA(1).EQ.'Y') THEN
C            KTYP94=1
C          ELSE
C            KTYP94=0
C          ENDIF

C cpb 28/6/96 Adding node based outer loops
          FORMAT='('' Specify whether outer integration loop is '
     '      //'[1]: '''//
     '      '/''   (1) over elements'''//
     '      '/''   (2) over nodes'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=BEMLOOPTYPE
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) BEMLOOPTYPE=IDATA(1)

          IF((IS_COUPLED(nx).AND.KTYP92.EQ.2).OR.(ITYP5(nr,nx).EQ.1.AND.
     '      ITYP2(nr,nx).EQ.3).OR.(ITYP5(nr,nx).EQ.4.AND.
     '      (ITYP2(nr,nx).EQ.3.OR.
     '      (ITYP2(nr,nx).EQ.5.AND.ITYP3(NR,NX).EQ.2)))) THEN
            IF(NJT.EQ.2) THEN
C cpb adding logarithmic Gaussian quadrature
              FORMAT='($,'' Do you want to use logarithmic Gaussian '
     '          //'quadrature for GQ [Y]? '',A)'
              IF(IOTYPE.EQ.3) THEN
                IF(BEMLOGGAUSS) THEN
                  ADATA(1)='Y'
                ELSE
                  ADATA(1)='N'
                ENDIF
              ENDIF
              CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     &          ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                IF(ADATA(1).EQ.'Y') THEN
                  BEMLOGGAUSS=.TRUE.
                ELSE
                  BEMLOGGAUSS=.FALSE.
                ENDIF
              ENDIF
            ELSE IF(NJT.EQ.3) THEN
              FORMAT='($,'' Do you want to use the optimised code '
     '          //'[Y]? '',A)'
              IF(IOTYPE.EQ.3) THEN
                IF(OPTI3DLAPLACE) THEN
                  ADATA(1)='Y'
                ELSE
                  ADATA(1)='N'
                ENDIF
              ENDIF
              CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     &          ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                IF(ADATA(1).EQ.'Y') THEN
                  OPTI3DLAPLACE=.TRUE.
                ELSE
                  OPTI3DLAPLACE=.FALSE.
                ENDIF
              ENDIF
            ELSE
              BEMLOGGAUSS=.FALSE.
              OPTI3DLAPLACE=.FALSE.
            ENDIF
          ELSE
            BEMLOGGAUSS=.FALSE.
            OPTI3DLAPLACE=.FALSE.
          ENDIF
        ELSE
c cpb 28/2/95 Must set up IGREN for the other regions
          CALL EQTYPE(IBT,NBH,NEELEM,nr,NW,nx,ERROR,*9999)
        ENDIF
      ENDIF !End of BEM solution (ityp4)

      IF(ITYP4(nr,nx).EQ.5) THEN
C RGB 12/8/97 Finite Volume solver options
!Advection scheme
        IDEFLT(1)=0
        FORMAT='('' Specify the scheme for advection '
     '    //'[0]: '''//
     '    '/''   (0) Upwind differencing'''//
     '    '/''   (1) Central differencing'''//
     '    '/''   (2) Hybrid differencing'''//
     '    '/''   (3) Interpolated donor differencing'''//
     '    '/''   (4) QUICK differencing'''//
     '    '/''   (5) ULTRA-QUICK differencing'''//
     '    '/''   (6) QUICK-EST differencing'''//
     '    '/''   (7) ULTRA-QUICK-EST differencing'''//
     '    '/''  *(8) FULL-QUICK differencing'''//
     '    '/''  *(9) ULTRA-FULL-QUICK differencing'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP61
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,9,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP61=IDATA(1)

C       Courant and Reynolds numbers
        RDEFLT(1) = 0.3D0
        WRITE(CHAR11,'(E11.4)') RDEFLT(1)
        CALL STRING_TRIM(CHAR11,IBEG2,IEND2)
        FORMAT='($,'' The Courant number is '
     '    //'['//CHAR11(IBEG2:IEND2)//']: '',E11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=COURANT
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.ne.3) THEN
          COURANT=RDATA(1)
        ENDIF
        RDEFLT(1) = 1000.0D0
        WRITE(CHAR11,'(E11.4)') RDEFLT(1)
        CALL STRING_TRIM(CHAR11,IBEG2,IEND2)
        FORMAT='($,'' The Reynolds number is '
     '    //'['//CHAR11(IBEG2:IEND2)//']: '',E11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=REYNOLD
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.ne.3) THEN
          REYNOLD=RDATA(1)
        ENDIF
        INV_REYNOLD=1.0D0/REYNOLD

        FORMAT='($,'' Is the mesh fixed ? [N] '',A)'
        ADEFLT(1)='N'
        IF(IOTYPE.EQ.3) THEN
          IF(MESHFIXD) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF((ADATA(1).EQ.'Y').OR.(ADATA(1).EQ.'y')) THEN
            MESHFIXD=.TRUE.
          ELSE
            MESHFIXD=.FALSE.
          ENDIF
        ENDIF
        IF(.NOT.MESHFIXD) THEN
          FORMAT='($,'' Do you wish to use centroid based mesh '
     '      //'fluxes ? [N] '',A)'
          ADEFLT(1)='N'
          IF(IOTYPE.EQ.3) THEN
            IF(CENTMESH) THEN
              ADATA(1)='Y'
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF((ADATA(1).EQ.'Y').OR.(ADATA(1).EQ.'y')) THEN
              CENTMESH=.TRUE.
            ELSE
              CENTMESH=.FALSE.
            ENDIF
          ENDIF
          FORMAT='($,'' Do you wish to check mesh integrity and '
     '      //'quality at each iteration ? [N] '',A)'
          ADEFLT(1)='N'
          IF(IOTYPE.EQ.3) THEN
            IF(MESHCHECK) THEN
              ADATA(1)='Y'
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF((ADATA(1).EQ.'Y').OR.(ADATA(1).EQ.'y')) THEN
              MESHCHECK=.TRUE.
            ELSE
              MESHCHECK=.FALSE.
            ENDIF
          ENDIF
          FORMAT='($,'' Enter the number of node smoothing iterations'
     '      //' (0 for none) [0]: '',I5)'
          IF(IOTYPE.EQ.3) IDATA(1)=SMOOTHIT
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,10,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            SMOOTHIT=IDATA(1)
          ENDIF
        ENDIF
        RDEFLT(1)=1.d-6
        WRITE(CHAR11,'(E11.4)') RDEFLT(1)
        CALL STRING_TRIM(CHAR11,IBEG2,IEND2)
        FORMAT='($,'' The steady state convergence tolerance is '
     '    //'['//CHAR11(IBEG2:IEND2)//']: '',E11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=DVTOL
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.ne.3) THEN
          DVTOL=RDATA(1)
        ENDIF
      ENDIF !ityp4 - Finite volume solve options

      ASSEMBLE_SOLUTION(nr,nx)=.FALSE.

      IF(IS_COUPLED(nx).AND.nr.EQ.COUP_NRLIST(COUP_NRLIST(0,nx),nx))
     '  THEN
        ASSEMBLE_SOLUTION(0,nx)=.FALSE.
        IWRIT1(0,nx)=IWRIT1(COUP_NRLIST(1,nx),nx)
        IWRIT2(0,nx)=IWRIT2(COUP_NRLIST(1,nx),nx)
        IWRIT3(0,nx)=IWRIT3(COUP_NRLIST(1,nx),nx)
        IWRIT4(0,nx)=IWRIT4(COUP_NRLIST(1,nx),nx)
        DO no_nrlist=2,COUP_NRLIST(0,nx)
          nrr=COUP_NRLIST(no_nrlist,nx)
          IWRIT1(nrr,nx)=IWRIT1(COUP_NRLIST(1,nx),nx)
          IWRIT2(nrr,nx)=IWRIT2(COUP_NRLIST(1,nx),nx)
          IWRIT3(nrr,nx)=IWRIT3(COUP_NRLIST(1,nx),nx)
          IWRIT4(nrr,nx)=IWRIT4(COUP_NRLIST(1,nx),nx)
        ENDDO
      ENDIF

C SGM 30 Oct 2000 ignore grid-based Finite element also
C MLT 2Dec2002 ignore grid finite volume also
      IF(ITYP4(nr,nx).NE.4.AND.ITYP4(nr,nx).NE.6.AND.
     '   ITYP4(nr,nx).NE.7) THEN !not collocation
        CALL GLOBALH(IBT,IDO,INP,NAN,NBH,NBJ,NELIST,NENP,NHE,NKB,NKHE,
     '    NKJE,NLL,NNB,NNF,NNL,NONY,NP_INTERFACE,NPF,NPL,NPNE,NPNY,nr,
     '    NRE,NVHE,NVHP,NVJE,NWP,nx,NXI,NYNE,NYNO,NYNP,NYNR,NYNY,NYQNR,
     '    CONY,CYNO,CYNY,SE,SP,XA,XE,XP,YP,FIX,ERROR,*9999)

        IF(IS_COUPLED(nx).AND.nr.EQ.COUP_NRLIST(COUP_NRLIST(0,nx),nx))
     '    THEN
          !Calculate the coupled mapping arrays (only after all regional
          !mapping arrays have been set up).
          CALL GLOBALC(NAN,NBH,NHE,NONY,NP_INTERFACE,NPNODE,NPNY,NRE,NW,
     '      nx,NXI,NYNE,NYNO,NYNP,NYNR,CONY,CYNO,FIX,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('IPSOLV')
      RETURN
 9999 CALL ERRORS('IPSOLV',ERROR)
      CALL EXITS('IPSOLV')
      RETURN 1
      END


