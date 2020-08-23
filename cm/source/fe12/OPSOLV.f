      SUBROUTINE OPSOLV(NPNY,nr,nx,nxc,ERROR,*)

C#### Subroutine: OPSOLV
C###  Description:
C###    OPSOLV outputs solution parameters for region nr.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'adam00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'eige00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'host00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp60.cmn'
      INCLUDE 'ktyp70.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'quas00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'time01.cmn'
!     Parameter List
      INTEGER NPNY(0:6,NYM,0:NRCM),nr,nx,nxc
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,IBEG,IBEG_min,IEND,IEND_max,nhost
      CHARACTER GLOBSTAT*6,SOCKSTAT*6,TITLE2_1(3)*36,TITLE2_2(0:3)*39,
     '  TITL22(3)*9,TITL23(2)*18,
     '  QUASITITLE(2)*34,BEMLOOPTITLE(1:2)*13,
     '  ADVSCHME(0:9)*32,INTEGRATION_TYPE(6)*45,
     '  ADAMS_ERROR_CONTROL(4)*23

      DATA ADVSCHME /'Upwind differencing            ',
     '                'Central differencing           ',
     '                'Hybrid differencing            ',
     '                'Interpolated donor differencing',
     '                'QUICK differencing             ',
     '                'ULTRA-QUICK differencing       ',
     '                'QUICK-EST differencing         ',
     '                'ULTRA-QUICK-EST differencing   ',
     '                'FULL-QUICK differencing        ',
     '                'ULTRA-FULL-QUICK differencing  '/

      DATA QUASITITLE /'Changing sources                  ',
     '                  'Changing boundary condition values'/

      DATA BEMLOOPTITLE/'over elements',
     '                   'over nodes   '/

      DATA INTEGRATION_TYPE
     '                /'Euler                                        ',
     '                 'Improved Euler                               ',
     '                 'Runge-Kutta (4th Order)                      ',
     '                 '                                             ',
     '                 'Adams-Moulton (variable order, adaptive step)',
     '                 'LSODA (adaptive step, stiff switching'/

      DATA ADAMS_ERROR_CONTROL
     '                /'Pure absolute          ',
     '                 'Relative to Y          ',
     '                 'Relative to DY         ',
     '                 'Mixed relative/absolute'/

C      DATA TITLE6 /'constant with respect to time       ',
C     '             'defined in subroutine USER          ',
C     '             'read from file IPC at each time step'/
      DATA TITLE2_1 /'BIE plus tangential derivative BIE  ',
     '               'BIE plus normal derivative BIE      ',
     '               'Normal and tangential derivative BIE'/

      DATA TITLE2_2 /'BIE plus s1 and s2 derivative BIEs     ',
     '               'BIE plus s1,s2 and s1s2 derivative BIEs',
     '               'BIE plus s1,s2 and n derivative BIEs   ',
     '               's1,s2,s1s2 and n derivative BIEs       '/
C      DATA TITLE9_1/'Direct solve                      ',
C     '              'Gauss-Seidel iterations           ',
C     '              'Gauss-Seidel with multigrid accel.',
C     '              'Unused                            '/
C      DATA TITLE9_2/'Newton-Raphson  (with update at  every iteration)',
C     '              'Modified Newton (with update every 6th iteration)',
C     '              'BFGS inverse    (with update every 6th iteration)',
C     '              'Sequential quadratic programming (E04UPF)'/
C      DATA TITL10 /'no search       ',
C     '             'Linear search   ',
C     '             'Quadratic search'/
C      DATA TITL15 /'Impulse fn.',
C     '             'Step fn.   ',
C     '             'Sine wave  ',
C     '             'White noise',2*' '/
      DATA TITL22 /'linear   ','quadratic','cubic    '/
      DATA TITL23 /'fixed time step   ','automatic stepping'/

      CALL ENTERS('OPSOLV',*9999)

      WRITE(OP_STRING,'(/'' Class '',I1,'' (nx='',I1,'
     '  //''') Region '',I1,'':'')') nxc,nx,nr
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      IF(ITYP5(nr,nx).EQ.1) THEN !static
        IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear
          CALL STRING_TRIM(FILE01,IBEG,IEND)
          IF(IBEG.GT.IEND) THEN
            WRITE(OP_STRING,'(/''  Pressure output file name:'',A/)')
     '        FILE01(IBEG:IEND)//'.PRESSURE'
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF !ityp6

      ELSE IF(ITYP5(nr,nx).EQ.2) THEN !time integration
        WRITE(OP_STRING,'(''  #time step intervals for history file'
     '    //' o/p ='',I3)') HIST_file_intervals
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(HIST_file_intervals.GT.0) THEN
          CALL STRING_TRIM(FILE02,IBEG,IEND)
          WRITE(OP_STRING,'(''  Output history file name is '','
     '      //'A,''.iphist'')') FILE02(IBEG:IEND)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C         IF(NODE_HISTORY_OP) THEN
C           WRITE(OP_STRING,'(''  Time history o/p at nodes: '',8I6)')
C     '       (NODE_HISTORY(n),n=1,NODE_HISTORY(0))
C           CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C         ENDIF
          IF(BINTIMEFILE.GT.0) THEN
            WRITE(OP_STRING,'(''  The output file is stored as '
     '        //'a binary file'')')
          ELSE
            WRITE(OP_STRING,'(''  The output file is stored as an '
     '        //'ascii file'')')
          ENDIF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF !history file

        WRITE(OP_STRING,'(''  Time integration algorithm is '',A,'
     '    //''' with '',A)') TITL22(KTYP22),TITL23(KTYP23)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        WRITE(OP_STRING,'(''  Initial time in seconds    (T0)= '','
     '    //'D11.4,'
     '    //'/,''  Final time in seconds      (T1)= '',D11.4,'
     '    //'/,''  Initial time increment     (DT)= '',D11.4,'
     '    //'/,''  Restart time        (T_RESTART)= '',D11.4,'
     '    //'/,''  Time-integration params (THETA)='',3(1X,D11.4))')
     '    T0,T1,TINCR,T_RESTART(nx),(THETA(i),i=1,KTYP22)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(ITYP2(nr,nx).EQ.9) THEN !Cardiac activation
          WRITE(OP_STRING,'(''  Type of integration (KTYP37)= '',A)')
     '      INTEGRATION_TYPE(KTYP37)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(KTYP37.EQ.5) THEN
            WRITE(OP_STRING,'(''  Maximum Adams polynomial order = '','
     '        //'I3)') ADAMS_MAX_ORDER
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''  Maximum Adams step size = '','
     '        //'D11.4)') ADAMS_MAX_STEP
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''  Maximum Adams iterations = '','
     '        //'I3)') ADAMS_MAX_ITERS
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''  Adams error control is = '','
     '        //'A)') ADAMS_ERROR_CONTROL(ADAMS_ERROR_TYPE)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF(ADAMS_ERROR_TYPE.EQ.4) THEN
              WRITE(OP_STRING,'(''  Absolute error component = '','
     '          //'D11.4)') ADAMS_ABS_ERR
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''  Relative error component = '','
     '          //'D11.4)') ADAMS_ABS_ERR
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
            IF(ADAMS_USE_ROUND_CTRL) THEN
              WRITE(OP_STRING,'(''  Adams rounding control is used'')')
            ELSE
              WRITE(OP_STRING,'(''  Adams rounding control is not '
     '          //'used'')')
            ENDIF
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(KTYP36.EQ.1) THEN
            WRITE(OP_STRING,'(''  Using DTAR (Dynamic Tracking)'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

      ELSE IF(ITYP5(nr,nx).EQ.3) THEN !modal
        WRITE(OP_STRING,
     '      '(''  No. of eigenpairs required (KTYP17) = '',I5,'
     '    //'/''  No. of eigenvalues found            = '',I5,'
     '    //'/''  Minimum eigenvalue                  = '',D11.4,'
     '    //'/''  Maximum eigenvalue                  = '',D11.4)')
     '    KTYP17,NUMEIGEN,EIGEN_LIMITS(1),EIGEN_LIMITS(2)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

c cpb 1/5/95 Replacing Fourier analysis with Quasi-static analysis
c      ELSE IF(ITYP5(nr,nx).EQ.4) THEN !Fourier
c        WRITE(OP_STRING,'(''  Lower frequency limit   (F0)   = '','
c     '    //'E11.3,'
c     '    //'/''  Upper frequency limit   (F1)    = '',E11.3,'
c     '    //'/''  Log frequency increment (dFreq)= '',E11.3)')
c     '    F0,F1,dFreq
c        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c        WRITE(OP_STRING,
c     '    '(''  Damping mass      coeff = '',E12.3)') DAMPING_FACTOR1
c        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c        WRITE(OP_STRING,
c     '    '(''  Damping stiffness coeff = '',E12.3)') DAMPING_FACTOR2
c        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c        WRITE(OP_STRING,
c     '    '(''  Damping frequency       = '',E12.3)') DAMPING_FREQUENCY
c        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
c      ENDIF

      ELSE IF(ITYP5(nr,nx).EQ.4) THEN !Quasi-static analysis
        WRITE(OP_STRING,'(''  The type of quasi-static problem is '','
     '    //'A)') QUASITITLE(QUASIPROB)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(QUASIPROB.EQ.2) THEN
          CALL STRING_TRIM(QUASI_INHISTFILE,IBEG,IEND)
          WRITE(OP_STRING,'(''  Input history file name is '','
     '      //'A)') QUASI_INHISTFILE(IBEG:IEND)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        WRITE(OP_STRING,'(''  The type of time stepping is '',A)')
     '    TITL23(QUASITIMESTEP)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        CALL STRING_TRIM(QUASI_OUTHISTFILE,IBEG,IEND)
        WRITE(OP_STRING,'(''  Output history file name is '','
     '    //'A)') QUASI_OUTHISTFILE(IBEG:IEND)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(BINTIMEFILE.GT.0) THEN
          WRITE(OP_STRING,'(''  The output file is stored as a binary '
     '      //'file'')')
        ELSE
          WRITE(OP_STRING,'(''  The output file is stored as an ascii '
     '      //'file'')')
        ENDIF
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(OUTPUT_SOLTIMES) THEN
           WRITE(OP_STRING,'(''  Timing information for each time step '
     '      //' is given'')')
         ELSE
           WRITE(OP_STRING,'(''  Timing information for each time step '
     '       //' is not given'')')
        ENDIF
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C        IF(NODE_HISTORY_OP) THEN
C          WRITE(OP_STRING,'(''  Time history o/p at nodes: '',8I6)')
C     '      (NODE_HISTORY(n),n=1,NODE_HISTORY(0))
C          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C        ENDIF
        WRITE(OP_STRING,'('
     '    //'/ ''  Initial time in seconds            (T0)='',D10.3,'
     '    //'/ ''  Final time in seconds              (T1)='',D10.3,'
     '    //'/ ''  Initial time increment             (DT)='',D10.3)')
     '    T0,T1,TINCR
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF !ityp5

      IF(ITYP4(nr,nx).LE.3) THEN      !fem or bem
        IF(ITYP6(nr,nx).EQ.1) THEN      !linear analysis
          WRITE(OP_STRING,'('
     '      //'''  Restart from previous solution (RESTAR)=  '',L1/,'
     '      //'''  Diagnostic output              (DOP)   =  '',L1/,'
     '      //'''  Mass matrix lumping            (LUMP)  =  '',L1/,'
     '      //'''  Number of solution equations   (NOT1)  ='',I5/,'
     '      //'''  Number of solution variables   (NOT2)  ='',I5/,'
     '      //'''  Output printing frequency      (IWRIT1)='',I5/,'
     '      //'''  Solution diagnostic output     (IWRIT4)='',I5)')
     '      RESTAR,DOP,LUMP,NOT(1,1,nr,nx),NOT(2,1,nr,nx),IWRIT1(nr,nx),
     '      IWRIT4(nr,nx)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(ITYP5(nr,nx).NE.3) THEN
            CALL OPSOLU_SOLVER(nx,ERROR,*9999)
          ENDIF

          IF((ITYP5(nr,nx).EQ.1.OR.ITYP5(nr,nx).EQ.4).AND.
     '      ITYP2(nr,nx).EQ.3) THEN
            IF(SALU_CONSISTENCY(nx)) THEN
              WRITE(OP_STRING,'(''  Salu consistency criterion is '
     '          //'used'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''  Fixed potential set at nk='',I2,'
     '          //''', nv='',I2,'', np='',I5,'', nr='',I2)')
     '          NPNY(1,SALU_NY(nx),0),NPNY(2,SALU_NY(nx),0),
     '          NPNY(4,SALU_NY(nx),0),NPNY(6,SALU_NY(nx),0)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE
              WRITE(OP_STRING,'(''  Salu consistency criterion is not '
     '          //'used'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
          IF(ITYP4(nr,nx).EQ.2) THEN !bem
            WRITE(OP_STRING,'(''  Using adaptive integration'
     '        //'    (ADAPINT)=  '',L1)')
     '        ADAPINT
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF(HERMITE) THEN
              WRITE(OP_STRING,'(''  Using hypersingular BEM'
     '          //'           (HYP)=  '',L1)') HYP
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              IF(NJT.EQ.2)THEN
                WRITE(OP_STRING,'(''  Solution procedure is '',A)')
     '            TITLE2_1(KTYP92)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ELSE
                WRITE(OP_STRING,'(''  Solution procedure is '',A)')
     '            TITLE2_2(KTYP92)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
C cpb 28/6/96 Adding node based outer loops
            WRITE(OP_STRING,'(''  Outer integration loop is '',A)')
     '        BEMLOOPTITLE(BEMLOOPTYPE)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C cpb 5/11/96 Adding logarithmic Gaussian quadrature
            IF(BEMLOGGAUSS) THEN
              WRITE(OP_STRING,'(''  Using logarithmic Gaussian '
     '          //'quadrature for GQ'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
            IF(OPTI3DLAPLACE) THEN
              WRITE(OP_STRING,'(''  Using optimised code for 3D '
     '          //'Laplaces equation'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
          IF(KTYP4.EQ.0) THEN
            WRITE(OP_STRING,'(''  No solution matrices written to '
     '        //'a file'')')
          ELSE IF(IABS(KTYP4).EQ.1) THEN
            WRITE(OP_STRING,'(''  Solution matrix is written to a '
     '        //'file'')')
          ELSE IF(IABS(KTYP4).EQ.2) THEN
            WRITE(OP_STRING,'(''  Solution RHS vector is '
     '        //'written to a file'')')
          ELSE IF(IABS(KTYP4).EQ.3) THEN
            WRITE(OP_STRING,'(''  Solution matrix and RHS vector are '
     '        //'written to a file'')')
          ENDIF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(KTYP4.GT.0) THEN
            WRITE(OP_STRING,'(''  The matrix file is an ascii file'')')
          ELSE IF(KTYP4.LT.0) THEN
            WRITE(OP_STRING,'(''  The matrix file is a binary file'')')
          ENDIF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear analysis
          WRITE(OP_STRING,'('
     '      //'''  Restart from previous solution (RESTAR)=  '',L1/,'
     '      //'''  Diagnostic output                 (DOP)=  '',L1/,'
     '      //'''  Number of global d.o.f.s          (NOT)='',I5/,'
     '      //'''  Type of equil.m iter.n   (ITYP9(nr,nx))='',I3/,'
     '      //'''  Material parameter incremented (KTYP14)='',I3/,'
     '      //'''  Pressure loads read from file  (KTYP71)='',I3/,'
     '      //'''  Output printing frequency      (IWRIT1)='',I3/,'
     '      //'''  Equil solns /Interm solns also (IWRIT2)='',I3/,'
     '      //'''  Soln vector /Resid vector also (IWRIT3)='',I3/,'
     '      //'''  Solution diagnostic output     (IWRIT4)='',I3)')
     '      RESTAR,DOP,NOT(2,1,nr,nx),ITYP9(nr,nx),KTYP14,KTYP71,
     '      IWRIT1(nr,nx),IWRIT2(nr,nx),IWRIT3(nr,nx),IWRIT4(nr,nx)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          CALL OPSOLU_SOLVER(nx,ERROR,*9999)

          IF(KTYP4.EQ.0) THEN
            WRITE(OP_STRING,'(''  No solution matrices written to '
     '        //'a file'')')
          ELSE IF(IABS(KTYP4).EQ.1) THEN
            WRITE(OP_STRING,'(''  Solution matrix is written to a '
     '        //'file'')')
          ELSE IF(IABS(KTYP4).EQ.2) THEN
            WRITE(OP_STRING,'(''  Solution matrix and RHS vector are '
     '        //'written to a file'')')
          ENDIF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(KTYP4.GT.0) THEN
            WRITE(OP_STRING,'(''  The matrix file is an ascii file'')')
          ELSE IF(KTYP4.LT.0) THEN
            WRITE(OP_STRING,'(''  The matrix file is a binary file'')')
          ENDIF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          IF(ITYP9(nr,nx).EQ.4) THEN !Sequential quadratic prog. (E04UPF)
            WRITE(OP_STRING,'('
     '        //'''  Print      level for E04UPF    (KTYP10)='',I3/,'
     '        //'''  Derivative level for E04UPF    (KTYP1D)='',I3)')
     '      KTYP10,KTYP1D
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF(KTYP1D.EQ.0) THEN       !derivs to be calc'd using FD's
              IF(DIFF_INTERVAL.EQ.0.D0) THEN    !compute FD interval
                WRITE(OP_STRING,'('
     '            //'''  Finite difference interval (DIFF_INTERVAL)='
     '            //' computed'')')
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ELSE
                WRITE(OP_STRING,'('
     '            //'''  Finite difference interval (DIFF_INTERVAL)='','
     '            //'E12.3)') DIFF_INTERVAL
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ELSE                       !derivs calc'd analytically
              WRITE(OP_STRING,'('
     '          //'''  Verify     level for E04UPF    (KTYP15)='',I3)')
     '          KTYP15
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE  !ITYP9(nr,nx)<=3
            WRITE(OP_STRING,
     '        '(''  Type of line search            (KTYP10)='',I3)')
     '        KTYP10
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

            WRITE(OP_STRING,
     '        '(''  Type of convergence criteria  (KTYP007)='',I3)')
     '        KTYP007
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF(KTYP007.EQ.2) THEN
              WRITE(OP_STRING,
     '          '(''    Gauss point field         (KTYP007a)='',I3)')
     '          KTYP007a
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF

            IF(KTYP1D.EQ.1) THEN
              WRITE(OP_STRING,'(''  Derivatives calculated '
     '          //'algebraically (KTYP1D=1)'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE IF(KTYP1D.EQ.2) THEN
              WRITE(OP_STRING,'(''  Derivatives approximated by '
     '          //'one-sided differences (KTYP1D=2)'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE IF(KTYP1D.EQ.3) THEN
              WRITE(OP_STRING,'(''  Derivatives approximated by '
     '          //'central differences (KTYP1D=3)'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF

            IF(KTYP1A.EQ.1) THEN
              WRITE(OP_STRING,'(''  Element stiffness matrices '
     '          //'computed in series (KTYP1A=1)'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE IF(KTYP1A.EQ.2) THEN
              WRITE(OP_STRING,'(''  Element stiffness matrices '
     '          //'computed in parallel (KTYP1A=2)'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              IF(AUTO_SPAWN) THEN
                WRITE(OP_STRING,'(''  Remote processes will spawn '
     '            //'automatically and be used after '',I3,'
     '            //''' secs on'')') IDELAY_SOCKET
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ELSE
                WRITE(OP_STRING,'(''  Remote processes need to be '
     '            //'started manually on'')')
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
              WRITE(OP_STRING,'(''     slave host machines: '
     '          //'Host | Port | ConnID | Socket Status | '
     '          //'Global Arrays'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C             Find longest host name
              IBEG_min=100
              IEND_max=0
              DO nhost=1,NUMHOSTS
                CALL STRING_TRIM(HOSTS(nhost),IBEG,IEND)
                IF(IBEG_min.GT.IBEG) IBEG_min=IBEG
                IF(IEND_max.LT.IEND) IEND_max=IEND
              ENDDO
              DO nhost=1,NUMHOSTS
                IF(SOCKET_OPEN(nhost)) THEN
                  SOCKSTAT='open  '
                ELSE
                  SOCKSTAT='closed'
                ENDIF
                IF(GLOBAL_SENT(nhost)) THEN
                  GLOBSTAT='sent  '
                ELSE
                  GLOBSTAT='unsent'
                ENDIF
                WRITE(OP_STRING,'(25X,'''
     '            //HOSTS(nhost)(IBEG_min:IEND_max)//' |'',I5,'
     '            //''' |   '',I2,''   |     '',A6,''    |    '',A6)')
     '            IPORT(nhost),ICONNID(nhost),SOCKSTAT,GLOBSTAT
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF

          ENDIF
        ENDIF

C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volume also
      ELSE IF(ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '        ITYP4(nr,nx).EQ.7) THEN !collocation or grid-based FE
        WRITE(OP_STRING,'('
     '    //'''  Restart from previous solution (RESTAR)=  '',L1/,'
     '    //'''  Diagnostic output                 (DOP)=  '',L1/,'
     '    //'''  Mass matrix lumping              (LUMP)=  '',L1/,'
     '    //'''  Number of global d.o.f.s          (NOT)='',I5/,'
     '    //'''  Type of solver           (ITYP9(nr,nx))='',I3/,'
     '    //'''  Global matrices written to file (KTYP4)='',I3/,'
     '    //'''  Output printing frequency      (IWRIT1)='',I3/,'
     '    //'''  Equil solns only/Int sols also (IWRIT2)='',I3/,'
     '    //'''  Soln vector only/Resid v. also (IWRIT3)='',I3/,'
     '    //'''  Linear solver diagnostic output(IWRIT4)='',I3)')
     '    RESTAR,DOP,LUMP,NOT(2,1,nr,nx),ITYP9(nr,nx),KTYP4,
     '    IWRIT1(nr,nx),IWRIT2(nr,nx),IWRIT3(nr,nx),IWRIT4(nr,nx)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(ITYP9(nr,nx).EQ.3) THEN !multigrid acceleration
          WRITE(OP_STRING,'('
     '      //'/''  #multigrid relax.ns on descent (Relax1)='',I3/,'
     '      //' ''  #multigrid relax.ns on  ascent (Relax2)='',I3)')
     '      Relax1(nx),Relax2(nx)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF !ktyp9

      ELSE IF(ITYP4(nr,nx).EQ.5) THEN !Finite Volume
        WRITE(OP_STRING,'('
     '    //'''  Mesh fixed                     (MESHFIXD) = '',L1/,'
     '    //'''  Global matrices written to file   (KTYP4) ='',I2/,'
     '    //'''  Output printing frequency        (IWRIT1) ='',I2/,'
     '    //'''  Advection scheme                 (KTYP61) ='',I2,
     '    '' = '',A)')
     '    MESHFIXD,KTYP4,IWRIT1(nr,nx),KTYP61,ADVSCHME(KTYP61)

        CALL OPSOLU_SOLVER(nx,ERROR,*9999)

      ENDIF !ityp4

      CALL EXITS('OPSOLV')
      RETURN
 9999 CALL ERRORS('OPSOLV',ERROR)
      CALL EXITS('OPSOLV')
      RETURN 1
      END


