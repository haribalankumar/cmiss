      SUBROUTINE OPEQUA(NBH,NEELEM,NELIST,NHE,NHP,NKH,
     '  NPL,NPNODE,nr,NVHE,NVHP,NW,nx,nxc,FULL,ERROR,*)

C#### Subroutine: OPEQUA
C###  Description:
C###    OPEQUA outputs equation parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b10.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'b40.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp40.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'stab00.cmn'
      INCLUDE 'titl40.cmn'
      INCLUDE 'titl50.cmn'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NHE(NEM),NHP(NPM),NKH(NHM,NPM,NCM),NPL(5,0:3,NLM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM),NW(NEM,3),nx,nxc
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,ie,j,k,nb,NBP,nc,ne,nh,NHTT,nhx,nn,
     '  noelem,nonode,np,nvar
      LOGICAL FULL
      CHARACTER TITLE0(6)*21,TITLE1(9)*28,TITLE2(12)*28,TITLE3(9)*28,
     '  TITLE4(5)*28,TITLE5(5)*28,
     '  TITLE6(2)*28,TITLE7(0:3)*50,
     '  TITLE10(6)*26,TITLE10A(3)*41,
     '  TITLE1_02(6)*28,TITLE1_03(5)*28,TITLE1_04(2)*28,
     '  TITLE1_05(2)*28,
     '  TITLE1_08(3)*28,TITLE1_09(2)*28,TITLE2_09(7)*28,
     '  TITLE2_09_1(10)*28,TITLE2_09_2(5)*28,TITLE2_09_7(4)*28,
     '  TITLE2_05(4)*28,TITLE5_03(5)*28

      DATA TITLE0 /'Static analysis      ',
     '             'Time integration     ',
     '             'Modal analysis       ',
     '             'Quasi-static analysis',
     '             'Front path analysis  ',
     '             'Buckling analysis    '/

      DATA TITLE1 /'Linear elasticity           ', !Static analysis
     '             'Finite elasticity           ',
     '             'Laplace equation            ',
     '             'Helmholtz equation          ',
     '             'div(k.grad(u))=f            ',
     '             'Linear 2nd order elliptic   ',
     '             'Biharmonic equation         ',
     '             'Fluid mechanics             ',
     '             'Oxygen transport            '/

      DATA TITLE2 /'Linear elasticity           ', !Time integration
     '             'Finite elasticity           ',
     '             'Advection-diffusion         ',
     '             'Wave equation               ',
     '             'Navier-Stokes equations     ',
     '             'Bio-heat transfer           ',
     '             'Maxwell equations           ',
     '             'Huygens activation          ',
     '             'Cellular based modelling    ',
     '             'Oxygen transport            ',
     '             'Humidity transport in lung  ',
     '             '                            '/

      DATA TITLE3 /'Linear elasticity           ', !Modal analysis
     '             '                            ',
     '             'Laplace equation            ',
     '             'Helmholtz equation          ',
     '             'div(k.grad(u))=f            ',
     '             'Linear 2nd order elliptic   ',
     '             'Biharmonic equation         ',
     '             'Vocal tract equations       ',
     '             '                            '/
      DATA TITLE4 /'Linear elasticity           ', !Quasi-static analy
     '             '                            ',
     '             'Laplace equation            ',
     '             '                            ',
     '             'div(k.grad(u))=f            '/

      DATA TITLE5 /'                            ', !Wavefront path an.
     '             '                            ',
     '             'Elliptic 1st order eikonal  ',
     '             '                            ',
     '             '                            '/

      DATA TITLE6 /'Linear elasticity           ', !Buckling analysis
     '             'Finite elasticity           '/

      DATA TITLE7 /'Not included                                  ',
     '             'Included as fixed initial strain              ',
     '             'Included but constrained by displacement b.c.s',
     '             'Included with full coupling to mechanics      '/

      DATA TITLE10/'Galerkin technique        ',
     '             'Direct boundary elements  ',
C PM 02-JUL-01: Incorrect numerical method
C     '             'Indirect boundary elements',
     '             'Finite difference         ',
     '             'Collocation               ',
     '             'Finite volume technique   ',
     '             'Grid-based Finite Element '/

      DATA TITLE10A/'Petrov-Galerkin with material derivatives',
     '              'Petrov-Galerkin with xi derivatives      ',
     '              'Stabilized Galerkin                      '/

C      DATA TITLE11/'constant with respect to time       ',
C     '             'defined in subroutine USER          ',
C     '             'read from file IPC at each time step'/

C *** Static analysis (ITYP5(nr,nx)=1)
      DATA TITLE1_02 /'Plane stress                ', !(ITYP2=2,ITYP3=1)
     '                'Plane strain                ', !               2
     '                '3-dimensional               ', !               3
     '                'Membrane theory             ', !               4
     '                'String theory               ', !               5
     '                'Shell theory                '/ !               6
      DATA TITLE1_03 /'Standard Laplace            ', !(ITYP2=3,ITYP3=1)
     '                'Generalised Laplace         ', !               2
     '                'Potential flow (aerofoil)   ', !               3
     '                'Unused                      ', !               4
     '                'Unused                      '/ !               5
      DATA TITLE1_04 /'Standard Helmholtz          ', !(ITYP2=4,ITYP3=1)
     '                'Modified Helmholtz (Yukawa) '/ !               2
      DATA TITLE1_05 /'f=General                   ', !(ITYP2=5,ITYP3=1)
     '                'f=-div(s.grad(g))           '/ !               2
      DATA TITLE1_08 /'Prandtl boundary layer eqtns', !(ITYP2=8,ITYP3=1)
     '                'Fluid interface stability   ', !               2
     '                'Constant volume constraint  '/ !               3
      DATA TITLE1_09 /'Multi-field oxygen transport', !(ITYP2=9,ITYP3=1)
     '                'Glucose-Oxygen transport    '/ !               2
C *** Time integration (ITYP5=2)
      DATA TITLE2_05 /'Fluid in elastic tube       ', !(ITYP2=5,ITYP3=1)
     '                'Lung gas transport          ', !               2
     '                'General Navier-Stokes       ', !               3
     '                'Stokes flow (No advection)  '/ !               4
      DATA TITLE2_09 /'Electrical                  ', !(ITYP2=9,ITYP19=1
     '                'Mechanical                  ', !                2
     '                'Metabolism                  ', !                3
     '                'Signalling Pathways         ', !                4
     '                'Drug Interaction            ', !                5
     '                '                            ', !                6
     '                'Coupled                     '/ !                7
      DATA TITLE2_09_1 /'Cubic - no recovery         ', !(ITYP2=9,ITYP19=1,ITYP3=1)
     '                'FitzHugh-Nagumo             ', !               2
     '                'vanCapelle-Durrer           ', !               3
     '                'Beeler-Reuter               ', !               4
     '                'Jaffri-Rice-Winslow         ', !               5
     '                'Luo-Rudy                    ', !               6
     '                'diFrancesco-Noble           ', !               7
     '                'Noble98                     ', !               8
     '                'Hodgkin-Huxley              ', !               9
     '                'User defined                '/ !              10
      DATA TITLE2_09_2 /'Hunter-McCulloch-ter Keurs  ', !(ITYP2=9,ITYP19=2,ITYP3=1)
     '                'Fading Memory               ', !               2
     '                'Distribution Moment         ', !               3
     '                'Infarct                     ', !               4
     '                'User Defined                '/ !               5
      DATA TITLE2_09_7 /'LR - HMT             ', !(ITYP2=9,ITYP19=7,ITYP3=1)
     '                'Noble 98 - HMT      ', !               2
     '                'HMT - SPM           ', !               3
     '                'User Defined        '/ !               4
C *** Wavefront path analysis
      DATA TITLE5_03 /'an isotropic material       ', !(ITYP2=3,ITYP3=1)
     '                'an orthotropic monodomain   ', !               2
     '                'an orthotropic bidomain     ', !               3
     '                'an invalid material         ', !               4
     '                'an invalid material         '/ !               5

      CALL ENTERS('OPEQUA',*9999)
      WRITE(OP_STRING,'(/'' Class '',I1,'' (nx='',I1,'
     '  //''') Region '',I1,'':'')') nxc,nx,nr
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      IF(ITYP5(nr,nx).GT.0) THEN
        WRITE(OP_STRING,'('' Type of analysis: '',A)')
     '    TITLE0(ITYP5(nr,nx))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(ITYP5(nr,nx).EQ.1) THEN      !Static analysis
          WRITE(OP_STRING,'('' Type of equation: '',A)')
     '      TITLE1(ITYP2(nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(ITYP5(nr,nx).EQ.2) THEN !Time integration
          WRITE(OP_STRING,'('' Type of equation: '',A)')
     '      TITLE2(ITYP2(nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(ITYP5(nr,nx).EQ.3) THEN !Modal analysis
          WRITE(OP_STRING,'('' Type of equation: '',A)')
     '      TITLE3(ITYP2(nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

c cpb 1/5/95 Replacing Fourier analysis with Quasi-static analysis
c        ELSE IF(ITYP5(nr,nx).EQ.4) THEN !Fourier Transform
c          WRITE(OP_STRING,
c     '    '('' Type of equation: '',A)')  TITLE4(ITYP2(nr,nx))
c          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE IF(ITYP5(nr,nx).EQ.4) THEN !Quasi-static analysis
          WRITE(OP_STRING,'('' Type of equation: '',A)')
     '      TITLE4(ITYP2(nr,nx))
           CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE IF(ITYP5(nr,nx).EQ.5) THEN !Wavefront path analysis
          IF(ITYP2(nr,nx).LE.5) THEN
            WRITE(OP_STRING,'('' Type of equation: '',A)')
     '        TITLE5(ITYP2(nr,nx))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

        ELSE IF(ITYP5(nr,nx).EQ.6) THEN !Buckling analysis
          WRITE(OP_STRING,'('' Type of equation: '',A)')
     '      TITLE6(ITYP2(nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

      ENDIF

!  Linear elasticity
      IF(ITYP2(nr,nx).EQ.1) THEN
        DO ie=1,12
          IF(ETYP(ie)) THEN
            j=0
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              IF(NW(ne,1).EQ.ie) THEN
                j=j+1
                NELIST(j)=ne
              ENDIF
            ENDDO
            WRITE(OP_STRING,'(2X,A,16I3/,(31X,16I3))')
     '        TITL42(ie),(NELIST(k),k=1,j)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO
        WRITE(OP_STRING,'('' Thermal strain  : '',A)') TITLE7(KTYP43)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

!  Finite elasticity
      ELSE IF(ITYP2(nr,nx).EQ.2) THEN
        WRITE(OP_STRING,'('' Type of elements: '',A)')
     '    TITLE1_02(KTYP51(nr))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(JTYP10.GE.2) THEN
          IF(ITYP10(nr).EQ.2) THEN
            FORMAT='(/''  Note: radial interpolation is in r^2'')'
          ELSE IF(ITYP10(nr).EQ.3) THEN
            FORMAT='(/''  Note: radial interpolation is in r^3'')'
          ELSE IF(ITYP10(nr).EQ.4.AND.JTYP10.EQ.2) THEN
            FORMAT='(/''  Note: Lambda interpolation is in'','
     '        //''' Focus^2.Sinh^2(Lambda)'')'
          ELSE IF(ITYP10(nr).EQ.4.AND.JTYP10.EQ.3) THEN
            FORMAT='(/''  Note: Lambda interpolation is in'','
     '        //''' Focus^3.Cosh(Lambda).Sinh^2(Lambda)'')'
          ENDIF
          WRITE(OP_STRING,FORMAT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      IF(ITYP5(nr,nx).EQ.1) THEN !static analysis

!Static: Laplace
        IF(ITYP2(nr,nx).EQ.3) THEN
          WRITE(OP_STRING,'(19X,A)')  TITLE1_03(ITYP3(nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

!Static: Helmholtz
        ELSE IF(ITYP2(nr,nx).EQ.4) THEN
          WRITE(OP_STRING,'(19X,A)')  TITLE1_04(ITYP3(nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

!Static: div(k.grad(u))=f equation
        ELSE IF(ITYP2(nr,nx).EQ.5) THEN
          WRITE(OP_STRING,'(19X,A)')  TITLE1_05(ITYP3(nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

!Static: Linear 2nd order elliptic
        ELSE IF(ITYP2(nr,nx).EQ.6) THEN

!Static: Biharmonic equation
        ELSE IF(ITYP2(nr,nx).EQ.7) THEN

!Static: Fluid mechanics
        ELSE IF(ITYP2(nr,nx).EQ.8) THEN
          WRITE(OP_STRING,'(19X,A)')  TITLE1_08(ITYP3(nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

!Static: Oxygen transport
        ELSE IF(ITYP2(nr,nx).EQ.9) THEN
        ENDIF

      ELSE IF(ITYP5(nr,nx).EQ.2) THEN !time integration

!Time-int: Advection-diffusion
        IF(ITYP2(nr,nx).EQ.3) THEN

!Time-int: Wave equation
        ELSE IF(ITYP2(nr,nx).EQ.4) THEN

!Time-int: Navier-Stokes
        ELSE IF(ITYP2(nr,nx).EQ.5) THEN
          WRITE(OP_STRING,'(19X,A)')  TITLE2_05(ITYP3(nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C PM 02-JUL-01: Not relevant for fluid in elastic tubes
          IF(ITYP3(nr,nx).NE.1) THEN
            WRITE(OP_STRING,'(1X,''Courant Number:'',E11.4)')COURANT
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(1X,''Reynolds Number:'',E11.4)')REYNOLD
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
!Time-int: Bio-heat transfer
        ELSE IF(ITYP2(nr,nx).EQ.6) THEN

!Time-int: Maxwell equations
        ELSE IF(ITYP2(nr,nx).EQ.7) THEN

!Time-int: Huygen's Cardiac activation
        ELSE IF(ITYP2(nr,nx).EQ.8) THEN
          IF(ITYP3(nr,nx).EQ.1) THEN
            IF(KTYP31.EQ.1) THEN
              WRITE(OP_STRING,'('' The model is implemented'
     '          //' forward in time'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE IF(KTYP31.EQ.2) THEN
              WRITE(OP_STRING,
     '          '('' The model is implemented backward in time'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF

!Time-int: Cellular based modelling
        ELSE IF(ITYP2(nr,nx).EQ.9) THEN
          WRITE(OP_STRING,'(19X,A)')  TITLE2_09(ITYP19(nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(ITYP19(nr,nx).EQ.1) THEN !electrical
            WRITE(OP_STRING,'(19X,A)')  TITLE2_09_1(ITYP3(nr,nx))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSEIF(ITYP19(nr,nx).EQ.2) THEN !mechanical
            WRITE(OP_STRING,'(19X,A)')  TITLE2_09_2(ITYP3(nr,nx))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSEIF(ITYP19(nr,nx).EQ.7) THEN !coupled
            WRITE(OP_STRING,'(19X,A)')  TITLE2_09_7(ITYP3(nr,nx))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE
            CALL ASSERT(.FALSE.,'Invalid model',ERROR,*9999)
          ENDIF

          IF(ITYP19(nr,nx).EQ.1) THEN !electrical
            IF(KTYP32.EQ.1) THEN
              WRITE(OP_STRING,'(19X,'' Using a monodomain model. '')')
            ELSE
              WRITE(OP_STRING,'(19X,'' Using a bidomain model. '')')
            ENDIF
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF(ITYP3(nr,nx).EQ.1) THEN !Cubic (and polynomial)
              IF(KTYP33.EQ.1) THEN
                WRITE(OP_STRING,'(19X,'' Cubic '')')
              ELSE IF(KTYP33.EQ.2) THEN
                WRITE(OP_STRING,'(19X,'' Quintic '')')
              ELSE IF(KTYP33.EQ.3) THEN
                WRITE(OP_STRING,'(19X,'' Heptic (order 7) '')')
              ENDIF
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE IF(ITYP3(nr,nx).EQ.2) THEN !FHN
              IF(KTYP33.EQ.1) THEN
                WRITE(OP_STRING,'(19X,'' Standard FHN. '')')
              ELSE IF(KTYP33.EQ.2) THEN
                WRITE(OP_STRING,'(19X,'' Rogers-McCulloch modified'
     '            //' FHN. '')')
              ELSE IF(KTYP33.EQ.3) THEN
                WRITE(OP_STRING,'(19X,'' Panfilov modified'
     '            //' FHN. '')')
              ENDIF
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE IF(ITYP3(nr,nx).EQ.3) THEN !VCD
              IF(KTYP33.EQ.1) THEN
                WRITE(OP_STRING,'(19X,'' Standard VCD. '')')
              ELSE
                IF(KTYP34.EQ.1) THEN
                  WRITE(OP_STRING,
     '              '(19X,'' California mod VCD for normal'
     '              //' tissue. '')')
                ELSE
                  WRITE(OP_STRING,'(19X,'' California mod VCD for'
     '              //' ischemic tissue. '')')
                ENDIF
              ENDIF
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE IF(ITYP3(nr,nx).EQ.4) THEN !Beeler-Reuter
              IF(KTYP33.EQ.1) THEN
                WRITE(OP_STRING,'(19X,'' Standard B-R.'')')
              ELSEIF (KTYP33.EQ.2) THEN
                WRITE(OP_STRING,'(19X,'' Ebihara-Johnson.'')')
              ELSEIF (KTYP33.EQ.3) THEN
                WRITE(OP_STRING,'(19X,'' Drouhard-Roberge.'')')
              ELSEIF (KTYP33.EQ.4) THEN
                WRITE(OP_STRING,'(19X,'' Defibrillation BR.'')')
              ENDIF
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              IF(KTYP33.EQ.4) THEN !defib
                IF(KTYP39.EQ.1) THEN
                  WRITE(OP_STRING,'(21X,'' Electroporation is '
     '              //'included.'')')
                ELSE
                  WRITE(OP_STRING,'(21X,'' Electroporation is not '
     '              //'included.'')')
                ENDIF
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(21X,'' The scale factor for '
     '            //'reducing the APD is: '',F8.2)') KTYP39R
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ELSE IF(ITYP3(nr,nx).EQ.5) THEN !JRW
              IF(KTYP33.EQ.1) THEN
                WRITE(OP_STRING,'(19X,'' Standard JRW.'')')
              ELSE
                WRITE(OP_STRING,'(19X,'' Princeton JRW.'')')
              ENDIF
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE IF(ITYP3(nr,nx).EQ.10) THEN !User defined
              IF (KTYP33.NE.6) THEN
                WRITE(OP_STRING,'(19X,'' User defined model '',I1)')
     '            KTYP33
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ELSE
                WRITE(OP_STRING,'(19X,'' User defined cellML model'')')
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
          ENDIF !ITYP19.EQ.1
!Time-int: Oxygen transport
        ELSE IF(ITYP2(nr,nx).EQ.9) THEN
          WRITE(OP_STRING,'(19X,A)')  TITLE1_09(ITYP3(nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(ITYP3(nr,nx).EQ.1) THEN !Multi-field oxygen transport
            WRITE(OP_STRING,'(18X,'' Number of fields = '',I1)')
     '        KTYP15
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

      ELSE IF(ITYP5(nr,nx).EQ.3) THEN !Modal analysis
c cpb 1/5/95 Replacing Fourier analysis with Quasi-static analysi
C      ELSE IF(ITYP5(nr,nx).EQ.4) THEN !Fourier Transform
      ELSE IF(ITYP5(nr,nx).EQ.4) THEN !Quasi-static analysis

!Static: Laplace
        IF(ITYP2(nr,nx).EQ.3) THEN
          WRITE(OP_STRING,'(19X,A)')  TITLE1_03(ITYP3(nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

!Static: div(k.grad(u))=f equation
        ELSE IF(ITYP2(nr,nx).EQ.5) THEN
          WRITE(OP_STRING,'(19X,A)')  TITLE1_05(ITYP3(nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ENDIF

      ELSE IF(ITYP5(nr,nx).EQ.5) THEN !Wavefront path analysis
        IF(ITYP2(nr,nx).EQ.3) THEN !Elliptic eikonal
          WRITE(OP_STRING,'(19X,''for '',A)') TITLE5_03(ITYP3(nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(ITYP5(nr,nx).EQ.6) THEN !Buckling analysis
      ENDIF

! Numerical method
      IF(ITYP4(nr,nx).GT.0) THEN
        IF(ITYP5(nr,nx).EQ.5.AND.ITYP2(nr,nx).EQ.3.
     '    AND.ITYP4(nr,nx).EQ.1.AND.ITYP15(nr,nx).NE.0) THEN
C         FEM upwinding for eikonal equation
          WRITE(OP_STRING,
     '      '('' Numerical method: Finite elements'',A)')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '('' Upwinding method: '',A)') TITLE10A(ITYP15(nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,
     '      '('' Numerical method: '',A)') TITLE10(ITYP4(nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C PM 02-JUL-01: Added info for fluid in elastic tube problem
          IF(ITYP3(nr,nx).EQ.1) THEN
            IF(ITYP16(nr,nx).EQ.1) THEN
              WRITE(OP_STRING,'(18X,'' Fully implicit scheme '')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSEIF(ITYP16(nr,nx).EQ.2) THEN
              WRITE(OP_STRING,'(18X,'' Fully explicit scheme '')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSEIF(ITYP16(nr,nx).EQ.3) THEN
              WRITE(OP_STRING,'(18X,'' Crank_Nicholson method '')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSEIF(ITYP16(nr,nx).EQ.4) THEN
              WRITE(OP_STRING,'(18X,'' Lax Wendroff Scheme '')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDIF
      ENDIF

! Linear/Nonlinear
      IF(ITYP6(nr,nx).EQ.1) THEN
        WRITE(OP_STRING,'('' Problem is linear'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(ITYP6(nr,nx).EQ.2) THEN
        WRITE(OP_STRING,'('' Problem is nonlinear'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(ITYP7(nr,nx).EQ.2) THEN
          WRITE(OP_STRING,'('' Boundary conditions are nonlinear'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(ITYP8(nr,nx).EQ.1) THEN
            WRITE(OP_STRING,'('' Radiation bdy conditions used'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

! JWF 7-3-02 Contact Status

      IF(ITYP5(nr,nx).EQ.1) THEN !static analysis
        IF(ITYP2(nr,nx).EQ.2) THEN !finite elasticity
          IF(KTYP51(nr).EQ.3) THEN !3D
            IF(KTYP5G(nr).EQ.0) THEN ! No Contact
              WRITE(OP_STRING,'('' Contact Option: No Contact'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE IF(KTYP5G(nr).EQ.1) THEN ! Penalty Method
              WRITE(OP_STRING,'('' Contact Option: Penalty Method'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF !KTYP5G
          ENDIF !KTYP51
        ENDIF !ITYP2
      ENDIF !ITYP5

! Time-varying coefficients
      IF(ITYP5(nr,nx).EQ.2) THEN !Time integration problem
        IF(KTYP3_equa(nx).EQ.1) THEN
          WRITE(OP_STRING,'('' Equation coeffs are constant in time'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(KTYP3_equa(nx).EQ.2) THEN
          WRITE(OP_STRING,'('' Time-varying equation coeffs are defined'
     '      //' in USER_IPEQUA subroutine'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(KTYP3_equa(nx).EQ.3) THEN
          CALL STRING_TRIM(FILE03,IBEG,IEND)
          WRITE(OP_STRING,'('' Time-varying equation coeffs are defined'
     '      //' in '//FILE03(IBEG:IEND)//'.IPEQUA_time'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

! Additional output
      IF(ITYP5(nr,nx).EQ.1) THEN !static
        IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.2) THEN !Fluid interface
          WRITE(OP_STRING,'('' Time increment = '',E12.3)') TINCR
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Gravitational accel.n = '',E12.3)')
     '      G_ACCEL
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

! Output basis functions
      IF(ITYP1(nr,nx).EQ.4) THEN !FE40
        IF(ITYP6(nr,nx).EQ.1) THEN !linear
          DO ie=1,12
            IF(ETYP(ie)) THEN
              j=0
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(NW(ne,1).EQ.ie) THEN
                  j=j+1
                  NELIST(j)=ne
                ENDIF
              ENDDO
              WRITE(OP_STRING,'(/3X,A)') TITL42(ie)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              DO nvar=1,NVE(ie)
                IF(ie.EQ.2.AND.NHE(NELIST(1)).EQ.1) THEN
                  nhx=2
                ELSEIF(ie.EQ.6.AND.NHE(NELIST(1)).EQ.1) THEN
                  nhx=3
                ELSE
                  nhx=NHV(nvar,ie)
                ENDIF
                  nh=NH_LOC(nhx,nx)
                WRITE(OP_STRING,
     '            '(5X,A,'' basis function type no. = '',I1)')
     '            TITL43(nhx,ie),NBH(nh,1,NELIST(1))
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF
          ENDDO !ie
        ENDIF !ityp6

      ELSE IF(ITYP1(nr,nx).EQ.5) THEN !FE50
        WRITE(OP_STRING,'('' Material is: '',A)') TITL52(KTYP52(nr))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        NHTT=NHT50(KTYP52(nr),KTYP51(nr))
        DO nhx=1,NHTT
          nh=NH_LOC(nhx,nx)
          WRITE(OP_STRING,
     '      '('' Element basis function types for depend. var. '','
     '      //'I1,'' (NBH('',I1,'',1,ne)):'',(/10(1X,I3,'': '',I2)))')
     '      nh,nh,(NEELEM(noelem,nr),NBH(nh,1,NEELEM(noelem,nr)),
     '      noelem=1,NEELEM(0,nr))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !nhx (nh)
        NBP=NBH(NH_LOC(NH_LOC(0,nx),nx),1,NEELEM(1,nr)) !hyd press basis
        IF(NBP.GT.0) THEN
          IF(NNT(NBP).EQ.0) NHTT=NHTT-1
        ENDIF
        WRITE(OP_STRING,'('' '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nhx=1,NHTT
          nh=NH_LOC(nhx,nx)
          WRITE(OP_STRING,
     '      '('' Numbers of global node derivatives defined for'
     '      //' dependent var. '',I1,'' (NKH('',I1,'',np,1)):'','
     '      //'(/10(1X,I3,'': '',I1)))') nh,nh,(NPNODE(nonode,nr),
     '      NKH(nh,NPNODE(nonode,nr),1),nonode=1,NPNODE(0,nr))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !nhx (nh)

      ELSE !FE30 & FE90
        IF(ITYP4(nr,nx).EQ.1) THEN !finite elements
          DO nhx=1,NHE(NEELEM(1,nr))
            nh=NH_LOC(nhx,nx)
            WRITE(OP_STRING,'(/'' Element basis function types for'
     '        //' dependent variable'',I2,'' (NBH('',I1,'',1,ne)):''/,'
     '        //'(10(1X,I3,'': '',I1)))') nh,nh,(NEELEM(noelem,nr),
     '        NBH(nh,1,NEELEM(noelem,nr)),noelem=1,NEELEM(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !nh
          DO nhx=1,NHE(NEELEM(1,nr))
            nh=NH_LOC(nhx,nx)
            WRITE(OP_STRING,'(/'' Element basis function types for'
     '        //' equation number '',I2,'' (NBH('',I1,'',2,ne)):''/,'
     '        //'(10(1X,I3,'': '',I1)))') nh,nh,(NEELEM(noelem,nr),
     '        NBH(nh,2,NEELEM(noelem,nr)),noelem=1,NEELEM(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !nh

C PM 02-JUL-01: Incorrect numerical method
C       ELSE IF(ITYP4(nr,nx).EQ.2.OR.ITYP4(nr,nx).EQ.3) THEN !bdry elements
        ELSE IF(ITYP4(nr,nx).EQ.2) THEN !bdry elements
          DO nhx=1,NHE(NEELEM(1,nr))
            nh=NH_LOC(nhx,nx)
            WRITE(OP_STRING,'(/'' Element basis function types for'
     '        //' dependent variable'',I2,'' (NBH('',I1,'',1,ne)):''/,'
     '        //'(10(1X,I3,'': '',I1)))') nh,nh,(NEELEM(noelem,nr),
     '        NBH(nh,1,NEELEM(noelem,nr)),noelem=1,NEELEM(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,
     '        '('' Element basis function type for its '','
     '        //'''normal derivative (NBH('',I1,'',2,ne)):''/,'
     '        //'(10(1X,I3,'': '',I2)))') nh,(NEELEM(noelem,nr),
     '        NBH(nh,2,NEELEM(noelem,nr)),noelem=1,NEELEM(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !nh

C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volume also
        ELSE IF(ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '          ITYP4(nr,nx).EQ.7) THEN !collocation or grid-based FE/FV
        ENDIF !ityp4
      ENDIF !ityp2

      IF(ITYP5(nr,nx).EQ.1.AND   !static analysis
     '  .ITYP2(nr,nx).EQ.3.AND   !Laplace equation
     '  .ITYP3(nr,nx).EQ.3) THEN !Aerofoil analysis
        CALL OPAERO(NPL,ERROR,*9999)
      ENDIF !aerofoil problem

C cpb 20/11/96 Adding bem curvature corrections for cubic Hermite q's
      IF(BEMCURVATURECORRECTION) THEN
        WRITE(OP_STRING,'(/'' Using curvature corrections for '
     '    //'GQ'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(KTYP24.GT.0) THEN
        WRITE(OP_STRING,'(/'' Global matrices are stored as a '
     '    //'sparse matrix'')')
      ELSE
        WRITE(OP_STRING,'(/'' Global matrices are stored as a '
     '    //'fully populated matrix'')')
      ENDIF
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      IF(FULL) THEN
        WRITE(OP_STRING,'(/'' #Versions for nodes:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO nhx=1,NHP(np)
            nh=NH_LOC(nhx,nx)
            WRITE(OP_STRING,'('' NVHP(nh='',I1,'',np='',I5,'
     '        //''',nc=1,2):'',2I3)')
     '        nh,np,(NVHP(nh,np,nc),nc=1,2)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !nh
        ENDDO !nonode

        WRITE(OP_STRING,'(/'' Version#s for elements:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nhx=1,NHE(ne)
            nh=NH_LOC(nhx,nx)
            nb=NBH(nh,1,ne)
            WRITE(OP_STRING,'('' NVHE(nn=1..,nb='',I1,'
     '        //''',nh='',I1,'',ne='',I5,''):'',100I3)')
     '        nb,nh,ne,(NVHE(nn,nb,nh,ne),nn=1,NNT(nb))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !nh
        ENDDO !nonode
      ENDIF !full

      CALL EXITS('OPEQUA')
      RETURN
 9999 CALL ERRORS('OPEQUA',ERROR)
      CALL EXITS('OPEQUA')
      RETURN 1
      END


