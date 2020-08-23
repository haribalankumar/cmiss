      BLOCK DATA BLK30

C#### BlockData: BLK30
C###  Description:
C###    BLK30 sets up parameter titles in common block

      IMPLICIT NONE
      INCLUDE 'b32.cmn'
      INCLUDE 'iltot00.cmn'
      INCLUDE 'titl30.cmn'


C *** Static analysis ***   (ityp5=1)
C ***********************

! ityp2= 1      2      3      4      5      6      7      8       9
      DATA ILTOT1/ !(njt,ityp2,ityp3) are #s of parameters
     ' 0,0,0, 0,0,0, 0,0,0, 1,1,1, 2,3,4, 0,0,0, 1,1,1, 2,2,2, 15,15,15,
     ' 0,0,0, 0,0,0, 1,2,3, 1,1,1, 2,2,2, 0,0,0, 5,5,5, 2,2,2, 8,8,8,
     ' 0,0,0, 0,0,0, 2,2,2, 0,0,0, 0,0,0, 0,0,0, 5,5,5, 0,0,0, 0,0,0,
     ' 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     ' 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0/

      CHARACTER
     '  STATIC_01_1(20)*34,                   !Linear elasticity
     '  STATIC_02_1(20)*34,                   !Finite elasticity
     '  STATIC_03_1(20)*34,STATIC_03_2(20)*34,!Laplace,gen,
     '                     STATIC_03_3(20)*34,!aerofoil
     '  STATIC_04_1(20)*34,STATIC_04_2(20)*34,!Helmholtz, Yukawa eqtn
     '  STATIC_05_1(20)*34,STATIC_05_2(20)*34,!div(k.grad(u))=f
     '  STATIC_06_1(20)*34,                   !Linear 2nd order elliptic
     '  STATIC_07_1(20)*34,STATIC_07_2(20)*34,!Maxwells equations
     '                     STATIC_07_3(20)*34,
     '  STATIC_08_1(20)*34,STATIC_08_2(20)*34,!Fluid mechanics
     '  STATIC_09_1(20)*34,STATIC_09_2(20)*34,!Oxy trans,Multi/Glucose
     '  STATIC_10_1(20)*34,
     '  STATIC_11_1(20)*34,
     '  STATIC_12_1(20)*34

      EQUIVALENCE
     '  (TITL31(1, 1,1),STATIC_01_1),
     '  (TITL31(1, 2,1),STATIC_02_1),
     '  (TITL31(1, 3,1),STATIC_03_1),(TITL31(1, 3,2),STATIC_03_2),
     '                               (TITL31(1, 3,3),STATIC_03_3),
     '  (TITL31(1, 4,1),STATIC_04_1),(TITL31(1, 4,2),STATIC_04_2),
     '  (TITL31(1, 5,1),STATIC_05_1),(TITL31(1, 5,2),STATIC_05_2),
     '  (TITL31(1, 6,1),STATIC_06_1),
     '  (TITL31(1, 7,1),STATIC_07_1),(TITL31(1, 7,2),STATIC_07_2),
     '                               (TITL31(1, 7,3),STATIC_07_3),
     '  (TITL31(1, 8,1),STATIC_08_1),(TITL31(1, 8,2),STATIC_08_2),
     '  (TITL31(1, 9,1),STATIC_09_1),(TITL31(1, 9,2),STATIC_09_2),
     '  (TITL31(1,10,1),STATIC_10_1),
     '  (TITL31(1,11,1),STATIC_11_1),
     '  (TITL31(1,12,1),STATIC_12_1)

      REAL*8
     '  DSTATIC_01_1(20),
     '  DSTATIC_02_1(20),
     '  DSTATIC_03_1(20),DSTATIC_03_2(20),DSTATIC_03_3(20),
     '  DSTATIC_04_1(20),DSTATIC_04_2(20),
     '  DSTATIC_05_1(20),DSTATIC_05_2(20),
     '  DSTATIC_06_1(20),
     '  DSTATIC_07_1(20),DSTATIC_07_2(20),DSTATIC_07_3(20),
     '  DSTATIC_08_1(20),DSTATIC_08_2(20),
     '  DSTATIC_09_1(20),DSTATIC_09_2(20),
     '  DSTATIC_10_1(20),
     '  DSTATIC_11_1(20),
     '  DSTATIC_12_1(20)

      EQUIVALENCE
     '  (RDMATE(1,1, 1,1),DSTATIC_01_1(1)),
     '  (RDMATE(1,1, 2,1),DSTATIC_02_1(1)),
     '  (RDMATE(1,1, 3,1),DSTATIC_03_1(1)),
     '  (RDMATE(1,1, 3,2),DSTATIC_03_2(1)),
     '  (RDMATE(1,1, 3,3),DSTATIC_03_3(1)),
     '  (RDMATE(1,1, 4,1),DSTATIC_04_1(1)),
     '  (RDMATE(1,1, 4,2),DSTATIC_04_2(1)),
     '  (RDMATE(1,1, 5,1),DSTATIC_05_1(1)),
     '  (RDMATE(1,1, 5,2),DSTATIC_05_2(1)),
     '  (RDMATE(1,1, 6,1),DSTATIC_06_1(1)),
     '  (RDMATE(1,1, 7,1),DSTATIC_07_1(1)),
     '  (RDMATE(1,1, 7,2),DSTATIC_07_2(1)),
     '  (RDMATE(1,1, 7,3),DSTATIC_07_3(1)),
     '  (RDMATE(1,1, 8,1),DSTATIC_08_1(1)),
     '  (RDMATE(1,1, 8,2),DSTATIC_08_2(1)),
     '  (RDMATE(1,1, 9,1),DSTATIC_09_1(1)),
     '  (RDMATE(1,1, 9,2),DSTATIC_09_2(1)),
     '  (RDMATE(1,1,10,1),DSTATIC_10_1(1)),
     '  (RDMATE(1,1,11,1),DSTATIC_11_1(1)),
     '  (RDMATE(1,1,12,1),DSTATIC_12_1(1))

C *** Laplace's equation           (ITYP2=3,ITYP3=1,ITYP5=1)
      DATA STATIC_03_1/20*' '/
      DATA DSTATIC_03_1/20*0.0d0/
C *** Generalised Laplace        **(ITYP2=3,ITYP3=2,ITYP5=1)
      DATA STATIC_03_2/'conductivity s11',
     '                 'conductivity s22',
     '                 'conductivity s33',
     '                  17*' '/
      DATA DSTATIC_03_2/1.0d0,1.0d0,1.0d0,17*0.0d0/
C *** Aerofoil Laplace           **(ITYP2=3,ITYP3=3,ITYP5=1)
      DATA STATIC_03_3/'free stream velocity','fluid density',18*' '/
      DATA DSTATIC_03_3/1.0d0,1.0d0,18*0.0d0/
C *** Helmholtz equation         **(ITYP2=4,ITYP3=1,ITYP5=1)
      DATA STATIC_04_1/'Wave number k',19*' '/
      DATA DSTATIC_04_1/1.0d0,19*0.0d0/
C *** Yukawa equation            **(ITYP2=4,ITYP3=2,ITYP5=1)
      DATA STATIC_04_2/'s value',19*' '/
      DATA DSTATIC_04_2/1.0d0,19*0.0d0/
C *** div(kgrad(u))=f (f General)**(ITYP2=5,ITYP3=1,ITYP5=1)
      DATA STATIC_05_1/'domain source f',
     '                 'permeability k11',
     '                 'permeability k22',
     '                 'permeability k33',16*' '/
      DATA DSTATIC_05_1/0.0d0,1.0d0,1.0d0,1.0d0,16*0.0d0/
C *** div(kgrad(u))=f, f=-div(s.grad(g)) **(ITYP2=5,ITYP3=2,ITYP5=1)
      DATA STATIC_05_2/'permeability k',
     '                 'permeability s',18*' '/
      DATA DSTATIC_05_2/1.0d0,1.0d0,18*0.0d0/
C *** Linear 2nd order elliptic  **(ITYP2=6,ITYP3=1,ITYP5=1)
      DATA STATIC_06_1/20*' '/
      DATA DSTATIC_06_1/20*0.0d0/
C *** Electrostatic equation     **(ITYP2=7,ITYP3=1,ITYP5=1)
      DATA STATIC_07_1/'Source term',19*' '/
      DATA DSTATIC_07_1/20*0.0d0/
C *** Magnetostatic equation, no gauge  **(ITYP2=7,ITYP3=2,ITYP5=1)
      DATA STATIC_07_2/'Source term J1',
     &                 'Source term J2',
     &                 'Source term J3',
     &                 'Free space permeability',
     &                 'Relative permeability',15*' '/
      DATA DSTATIC_07_2/0.0d0,0.0d0,0.0d0,12.566371d-10,1.0d0,15*0.0d0/
C *** Magnetostatic equation, Coulomg gauge **(ITYP2=7,ITYP3=2,ITYP5=1)
      DATA STATIC_07_3/'Source term J1',
     &                 'Source term J2',
     &                 'Source term J3',
     &                 'Free space permeability',
     &                 'Relative permeability',15*' '/
      DATA DSTATIC_07_3/0.0d0,0.0d0,0.0d0,12.566371d-10,1.0d0,15*0.0d0/
C *** Boundary layer equations   **(ITYP2=8,ITYP3=1,ITYP5=1)
      DATA  STATIC_08_1/'Kinematic viscosity (m^2/s)',
     '                  'Fluid density      (kg/m^3)',18*' '/
      DATA DSTATIC_08_1/1.0d0,1.0d0,18*0.0d0/
C *** Fluid interface stability  **(ITYP2=8,ITYP3=2,ITYP5=1)
      DATA STATIC_08_2/'Fluid density ',
     '                 'Fluid velocity',18*' '/
      DATA DSTATIC_08_2/1.0d0,1.0d0,18*0.0d0/
C *** Multi-field oxy transport  **(ITYP2=9,ITYP3=1,ITYP5=1)
      DATA STATIC_09_1/'Oxy. m. solubility (nmol/kPa/mm^3)',
     '                 'Oxy. b. solubility (nmol/kPa/mm^3)',
     '                 'Oxy. m. diffusivity (mm^2/ks)     ',
     '                 'Oxy. b. diffusivity (mm^2/ks)     ',
     '                 'Myoglobin diffusivity (mm^2/ks)   ',
     '                 'Max myoglobin   conc. (nmol/mm^3) ',
     '                 'Max haemoglobin conc. (nmol/mm^3) ',
     '                 'Max mitochon. rate (nmol/ks/mm^3) ',
     '                 'Exchange coefficient (1/ks)       ',
     '                 'Arterial blood oxygen p.p. (kPa)  ',
     '                 'p_50 for myoglobin (kPa)          ',
     '                 'p_50 for haemoglobin (kPa)        ',
     '                 'n for haemoglobin                 ',
     '                 'p_50 for mitochondria (kPa)       ',
     '                 'Relative volume of vascular bed   ',5*' '/
      DATA DSTATIC_09_1/0.01125d0, !Oxy. m. sol. (nmol/kPa/mm^3)   ( 1)
     '                  0.0093d0,  !Oxy. b. sol. (nmol/kPa/mm^3)   ( 2)
     '                  2.3d0,     !Oxy. m. diffusivity (mm^2/ks)  ( 3)
     '                  1.5d0,     !Oxy. b. diffusivity (mm^2/ks)  ( 4)
     '                  0.07d0,    !Myoglobin diffusivity (mm^2/ks)( 5)
     '                  0.28d0,    !Max myoglobin conc. (nmol/mm^3)( 6)
     '                  8.8d0,     !Max haemoglob conc. (nmol/mm^3)( 7)
     '                  20.0d0,    !Max mitoch rate (nmol/ks/mm^3) ( 8)
     '                  0.1d0,     !Exchange coefficient (1/ks)    ( 9)
     '                  100.0d0,   !Arterial blood oxy p.p. (kPa)  (10)
     '                  0.6d0,     !p_50 for myoglobin (kPa)       (11)
     '                  3.0d0,     !p_50 for haemoglobin (kPa)     (12)
     '                  2.6d0,     !n for haemoglobin              (13)
     '                  0.1d0,     !p_50 for mitochondria (kPa)    (14)
     '                  0.0d0,     !Relative vol of vascular bed   (15)
     '                  5*0.0d0/
C *** Glucose-Oxygen transport   **(ITYP2=9,ITYP3=2,ITYP5=1)
      DATA STATIC_09_2/'Oxygen  solubility (nmol/kPa/mm^3)',
     '                 'Glucose solubility (nmol/kPa/mm^3)',
     '                 'Oxygen  diffusivity (mm^2/ks)     ',
     '                 'Glucose diffusivity (mm^2/ks)     ',
     '                 'Oxygen  Michaelis p.pressure (kPa)',
     '                 'Glucose Michaelis p.pressure (kPa)',
     '                 'Vmax (nmol/ks/mm^3)               ',
     '                 'Oxygen reaction weight (nu)       ',12*' '/
      DATA DSTATIC_09_2/20*0.0d0/


C *** Time integration ***   (ityp5=2)
C ************************

! ityp2= 1      2      3      4      5      6      7      8       9
      DATA ILTOT2/ !(njt,ityp2(nr),ityp3(nr)) are #s of parameters
     '  0,0,0, 0,0,0, 6,9,12, 1,2,3, 6,6,12, 5,5,5, 0,0,0, 5,6,7,
     '  8,8,8,  10,10,10, 2,2,2, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 2,2,2, 0,0,0, 0,0,0, 0,0,0,
     '  8,8,8,     8,8,8, 6,6,6, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 2,2,2, 0,0,0, 0,0,0, 0,0,0,
     '  8,8,8,     0,0,0, 3,3,3, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 3,4,5, 0,0,0, 0,0,0, 0,0,0,
     '  8,8,8,     0,0,0, 8,8,8, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 8,10,12, 0,0,0, 0,0,0, 0,0,0,
     '  8,8,8,     0,0,0, 0,0,0, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     '  8,8,8,        0,0,0, 3,3,3, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     '  8,8,8,     0,0,0, 7,7,7, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     '  8,8,8,     0,0,0, 0,0,0, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     '  8,8,8,        0,0,0, 0,0,0, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     '  8,8,8,        0,0,0, 0,0,0, 0,0,0/

      CHARACTER
     '  DYNAM_01_1(20)*34,                    !Linear elasticity
     '  DYNAM_02_1(20)*34,                    !Finite elasticity
     '  DYNAM_03_1(20)*34,                    !Advection-diff eqn
     '  DYNAM_04_1(20)*34,                    !Wave equation
     '  DYNAM_05_1(20)*34,                    !Fluid in elastic tube
     '                    DYNAM_05_2(20)*34,  !Lung gas transport
     '                    DYNAM_05_3(20)*34,  !General Navier-Stokes
     '                    DYNAM_05_4(20)*34,  !Stokes' Flow
     '                    DYNAM_05_5(20)*34,  !Bidomain Stokes' Flow
     '  DYNAM_06_1(20)*34,                    !Bio-heat transfer eqn
     '  DYNAM_07_1(20)*34,                    !Maxwell's equations
     '  DYNAM_08_1(20)*34,                    !Threshold activation
     '  DYNAM_09_1(20)*34,                    !Cubic w/o recovery
     '                    DYNAM_09_2(20)*34,  !FitzHugh-Nagumo eqns
     '                    DYNAM_09_3(20)*34,  !VanCapelle-Durrer
     '                    DYNAM_09_4(20)*34,  !Beeler-Reuter
     '                    DYNAM_09_5(20)*34,  !Jafri-Rice-Winslow
     '                    DYNAM_09_6(20)*34,  !Luo-Rudy
     '                    DYNAM_09_7(20)*34,  !diFrancesco-Noble
     '                    DYNAM_09_8(20)*34,  !Noble 98
     '                    DYNAM_09_9(20)*34,  !Hodgkin-Huxley
     '                    DYNAM_09_10(20)*34, !User defined
     '  DYNAM_10_1(20)*34,DYNAM_10_2(20)*34,  !Oxy transport: Multi/Glucose
     '  DYNAM_11_1(20)*34,                    !Lung gas mixing
     '                    DYNAM_11_2(20)*34,  !Lung temp/vapour transp
     '  DYNAM_11_3(20)*34,                    !Lung capillary blood
     '                    DYNAM_11_4(20)*34,  !Lung simple P-R-Flow
     '  DYNAM_11_6(20)*34,                    !Lung arteriole-capillary-venule blood
     '  DYNAM_11_7(20)*34,                    !Lung gas exchange with pre-defined transport
     '  DYNAM_12_1(20)*34

      EQUIVALENCE
     '  (TITL32(1, 1,1),DYNAM_01_1),
     '  (TITL32(1, 2,1),DYNAM_02_1),
     '  (TITL32(1, 3,1),DYNAM_03_1),
     '  (TITL32(1, 4,1),DYNAM_04_1),
     '  (TITL32(1, 5,1),DYNAM_05_1),(TITL32(1, 5,2),DYNAM_05_2),
     '                              (TITL32(1, 5,3),DYNAM_05_3),
     '                              (TITL32(1, 5,4),DYNAM_05_4),
     '                              (TITL32(1, 5,5),DYNAM_05_5),
     '  (TITL32(1, 6,1),DYNAM_06_1),
     '  (TITL32(1, 7,1),DYNAM_07_1),
     '  (TITL32(1, 8,1),DYNAM_08_1),
     '  (TITL32(1, 9,1),DYNAM_09_1),(TITL32(1, 9,2),DYNAM_09_2),
     '                              (TITL32(1, 9,3),DYNAM_09_3),
     '                              (TITL32(1, 9,4),DYNAM_09_4),
     '                              (TITL32(1, 9,5),DYNAM_09_5),
     '                              (TITL32(1, 9,6),DYNAM_09_6),
     '                              (TITL32(1, 9,7),DYNAM_09_7),
     '                              (TITL32(1, 9,8),DYNAM_09_8),
     '                              (TITL32(1, 9,9),DYNAM_09_9),
     '                              (TITL32(1, 9,10),DYNAM_09_10),
     '  (TITL32(1,10,1),DYNAM_10_1),(TITL32(1,10,2),DYNAM_10_2),
     '  (TITL32(1,11,1),DYNAM_11_1),(TITL32(1,11,2),DYNAM_11_2),
     '  (TITL32(1,11,3),DYNAM_11_3),(TITL32(1,11,4),DYNAM_11_4),
     '  (TITL32(1,11,6),DYNAM_11_6),
     '  (TITL32(1,11,7),DYNAM_11_7),
     '  (TITL32(1,12,1),DYNAM_12_1)

      REAL*8
     '  DDYNAM_01_1(20),
     '  DDYNAM_02_1(20),
     '  DDYNAM_03_1(20),
     '  DDYNAM_04_1(20),
     '  DDYNAM_05_1(20),DDYNAM_05_2(20),DDYNAM_05_3(20),DDYNAM_05_4(20),
     '  DDYNAM_05_5(20),
     '  DDYNAM_06_1(20),
     '  DDYNAM_07_1(20),
     '  DDYNAM_08_1(20),
     '  DDYNAM_09_1(20),DDYNAM_09_2(20),DDYNAM_09_3(20),DDYNAM_09_4(20),
     '  DDYNAM_09_5(20),DDYNAM_09_6(20),DDYNAM_09_7(20),DDYNAM_09_8(20),
     '  DDYNAM_09_9(20),DDYNAM_09_10(20),
     '  DDYNAM_10_1(20),DDYNAM_10_2(20),
     '  DDYNAM_11_1(20),DDYNAM_11_2(20),DDYNAM_11_3(20),DDYNAM_11_4(20),
     '  DDYNAM_11_6(20),DDYNAM_11_7(20),DDYNAM_12_1(20)

      EQUIVALENCE
     '  (RDMATE(1,2, 1,1),DDYNAM_01_1(1)),
     '  (RDMATE(1,2, 2,1),DDYNAM_02_1(1)),
     '  (RDMATE(1,2, 3,1),DDYNAM_03_1(1)),
     '  (RDMATE(1,2, 4,1),DDYNAM_04_1(1)),
     '  (RDMATE(1,2, 5,1),DDYNAM_05_1(1)),
     '  (RDMATE(1,2, 5,2),DDYNAM_05_2(1)),
     '  (RDMATE(1,2, 5,3),DDYNAM_05_3(1)),
     '  (RDMATE(1,2, 5,4),DDYNAM_05_4(1)),
     '  (RDMATE(1,2, 5,5),DDYNAM_05_5(1)),
     '  (RDMATE(1,2, 6,1),DDYNAM_06_1(1)),
     '  (RDMATE(1,2, 7,1),DDYNAM_07_1(1)),
     '  (RDMATE(1,2, 8,1),DDYNAM_08_1(1)),
     '  (RDMATE(1,2, 9,1),DDYNAM_09_1(1)),
     '  (RDMATE(1,2, 9,2),DDYNAM_09_2(1)),
     '  (RDMATE(1,2, 9,3),DDYNAM_09_3(1)),
     '  (RDMATE(1,2, 9,4),DDYNAM_09_4(1)),
     '  (RDMATE(1,2, 9,5),DDYNAM_09_5(1)),
     '  (RDMATE(1,2, 9,6),DDYNAM_09_6(1)),
     '  (RDMATE(1,2, 9,7),DDYNAM_09_7(1)),
     '  (RDMATE(1,2, 9,8),DDYNAM_09_8(1)),
     '  (RDMATE(1,2, 9,9),DDYNAM_09_9(1)),
     '  (RDMATE(1,2, 9,10),DDYNAM_09_10(1)),
     '  (RDMATE(1,2,10,1),DDYNAM_10_1(1)),
     '  (RDMATE(1,2,10,2),DDYNAM_10_2(1)),
     '  (RDMATE(1,2,11,1),DDYNAM_11_1(1)),
     '  (RDMATE(1,2,11,2),DDYNAM_11_2(1)),
     '  (RDMATE(1,2,11,3),DDYNAM_11_3(1)),
     '  (RDMATE(1,2,11,4),DDYNAM_11_4(1)),
     '  (RDMATE(1,2,11,6),DDYNAM_11_6(1)),
     '  (RDMATE(1,2,11,7),DDYNAM_11_7(1)),
     '  (RDMATE(1,2,12,1),DDYNAM_12_1(1))

C *** Advection-diffusion equation (ITYP2=3,ITYP3=1,ITYP5=2)
      DATA DYNAM_03_1 /'Source term (constant)           ',
     '                 'Source term (linear)             ',
     '                 'Storage coefficient              ',
     '                 'Diffusion coefficient (fibre_1)  ',
     '                 'Mass flow (fibre_1)              ',
     '                 'Grad Mass flow (fibre_1)         ',
     '                 'Diffusion coefficient (fibre_2)  ',
     '                 'Mass flow (fibre_2)              ',
     '                 'Grad Mass flow (fibre_2)         ',
     '                 'Diffusion coefficient (fibre_3)  ',
     '                 'Mass flow (fibre_3)              ',
     '                 'Grad Mass flow (fibre_3)',8*'    '/
      DATA DDYNAM_03_1/20*0.0d0/
C *** Wave equation (ITYP2=4,ITYP3=1,ITYP5=2)
      DATA DYNAM_04_1 /'Wave speed (x)',
     '                 'Wave speed (y)',
     '                 'Wave speed (z)',17*' '/
      DATA DDYNAM_04_1/20*0.0d0/
C *** Fluid in elastic tube (ITYP2=5,ITYP3=1,ITYP5=2)
      DATA DYNAM_05_1 /'Fluid density (kg/mm^3) ',
     '                 'Fluid viscosity  (mm^2/sec) ',
     '                 'Velocity profile factor Alpha   ',
     '                 'Unstrained radius (mm)',
     '                 'Elastic wall const Go_a (kPa)',
     '                 'Elastic wall const Beta_a ',
     '                 'Elastic wall const kappa ',
     '                 'Host element ',
     '                 'Elastic wall const beta_v ',
     '                 'Elastic wall const beta2_v ',
     '                 'Elastic wall const Go_v ',
     '                 'Elastic wall const beta2_a ',8*' '/
      DATA DDYNAM_05_1/0.00000105d0,3.2d0,1.1d0,1.0d0,21.20d0,2.0d0,
     '  2.0d0,1.0d0,0.0d0,2.0d0,21.20d0,2.0d0,8*0.0d0/
C *** Lung gas transport (ITYP2=5,ITYP3=2,ITYP5=2)
      DATA DYNAM_05_2 /'Gas diffusion coefficient (mm^2/s)',
     '                 'Oxygen loss coefficient           ',18*' '/
      DATA DDYNAM_05_2/74.3d0,0.d0,18*0.0d0/
C *** General Navier-Stokes (ITYP2=5,ITYP3=3,ITYP5=2)
      DATA DYNAM_05_3 /'Kinematic viscosity (m^2/s)',
     '                 'Fluid density (kg/m^3)',18*' '/
C *** Stokes' Flow (ITYP2=5,ITYP3=4,ITYP5=2)
      DATA DYNAM_05_4 /'Source term',
     '                 'Fluid density',
     '                 'Kinematic viscosity (fibre_1)',
     '                 'Kinematic viscosity (fibre_2)',
     '                 'Kinematic viscosity (fibre_3)',15*' '/
C *** Bidomain Stokes' Flow (ITYP2=5,ITYP3=5,ITYP5=2)
      DATA DYNAM_05_5 /'Ext. Source term',
     '                 'Int. Source term',
     '                 'Ext. Permeability',
     '                 'Int. Permeability',
     '                 'Ext. Fluid density',
     '                 'Int. Fluid density',
     '                 'Ext. Kinematic viscosity (fibre_1)',
     '                 'Int. Kinematic viscosity (fibre_1)',
     '                 'Ext. Kinematic viscosity (fibre_2)',
     '                 'Int. Kinematic viscosity (fibre_2)',
     '                 'Ext. Kinematic viscosity (fibre_3)',
     '                 'Int. Kinematic viscosity (fibre_3)',8*' '/
      DATA DDYNAM_05_3/1.d0,1.d0,18*0.0d0/
CC *** Bio-heat transfer equation (ITYP2=6,ITYP3=1,ITYP5=2)
      DATA DYNAM_06_1 /'Metabolic heat prod.n',
     '                 'Storage coefficient  ',
     '                 'Diffusion coefficient',
     '                 'Capillary temperature',
     '                 'Capillary heat flow  ',15*' '/
      DATA DDYNAM_06_1/20*0.0d0/
C *** Maxwell's equations (ITYP2=7,ITYP3=1,ITYP5=2)
      DATA DYNAM_07_1 /20*' '/
      DATA DDYNAM_07_1/20*0.0d0/
C *** Threshold activation (ITYP2=8,ITYP3=1,ITYP5=2)
      DATA DYNAM_08_1 /'Wavespeed -fibre direction (mm/ms)',
     '                 'Wavespeed -sheet direction (mm/ms)',
     '                 'Wavespeed -transverse dir. (mm/ms)',
     '                 'Speedup factor for subendo layers ',
     '                 '# of Gauss points searched in Xi1 ',
     '                 '# of Gauss points searched in Xi2 ',
     '                 '# of Gauss points searched in Xi3 ',13*' '/
      DATA DDYNAM_08_1/0.67d0,0.5d0,0.34d0,1.d0,1.d0,1.d0,1.d0,13*0.0d0/
C *** Cable equation (1D) (ITYP2=8,ITYP3=2,ITYP5=2)
C      DATA DYNAM_08_2 /'Capacitance coefficient (uF/mm)   ',
C     '                 'Diffusion coefficient (uF.m/s)    ',
C     '                 'Current scaling coefficient       ',
C     '                 '1st voltage root (mV)             ',
C     '                 '2nd voltage root (mV)             ',15*' '/
C      DATA DDYNAM_08_2/20*0.0d0/
C *** Noble-DiFrancesco equations (ITYP2=8,ITYP3=4,ITYP5=2)
C      DATA DYNAM_08_4 /'Axial resistance (Megohms/mm)     ',
C     '                 'Membrane capacitance (uF/mm)      ',
C     '                 'Sodium conductance (uMho/mm)      ',
C     '                 'Potassium conductance (uMho/mm)   ',16*' '/
C      DATA DDYNAM_08_4/20*0.0d0/
C *** Cubic w/o recovery (ITYP2=9,ITYP3=1,ITYP5=2)
C ***   default values rom Colli-Franzone in HPC
C *** DPN 14 March 2000 - fixing up units
C      DATA DYNAM_09_1 /'Membrane capacitance Cm (nf/mm^2) ',
C     '                 'Surface-to-volume ratio Am (1/mm) ',
C     '                 'Int. conductivity-fibre (Mho/mm)  ',
C     '                 'Int. conductivity-sheet (Mho/mm)  ',
C     '                 'Int. conductivity-cross (Mho/mm)  ',
C     '                 'Ext. conductivity-fibre (Mho/mm)  ',
C     '                 'Ext. conductivity-sheet (Mho/mm)  ',
C     '                 'Ext. conductivity-cross (Mho/mm)  ',
C     '                 'Resting potential (mV)            ',
C     '                 'Plateau potential (mV)            ',
C     '                 'Threshold potential (mV)          ',
C     '                 'Membrane cond (uMho/mm^2)         ',
C     '                 8*' '/
C      DATA DDYNAM_09_1/10.0d0,200.0d0,
C     '  0.25d-3,0.125d-3,0.125d-3,0.2d-3,0.0416d-3,0.0416d-3,
C     '  -85.0d0,15.0d0,-75.0d0,4.0d0,8*0.0d0/
      DATA DYNAM_09_1 /'Membrane capacitance Cm (uF/mm^2) ',
     '                 'Surface-to-volume ratio Am (1/mm) ',
     '                 'Int. conductivity-fibre (mS/mm)   ',
     '                 'Int. conductivity-sheet (mS/mm)   ',
     '                 'Int. conductivity-cross (mS/mm)   ',
     '                 'Ext. conductivity-fibre (mS/mm)   ',
     '                 'Ext. conductivity-sheet (mS/mm)   ',
     '                 'Ext. conductivity-cross (mS/mm)   ',
     '                 'Resting potential (mV)            ',
     '                 'Plateau potential (mV)            ',
     '                 'Threshold potential (mV)          ',
     '                 'Membrane cond (mS/mm^2)           ',
     '                 8*' '/
      DATA DDYNAM_09_1/10.0d-3,200.0d0,
     '  0.25d0,0.125d0,0.125d0,0.2d0,0.0416d0,0.0416d0,
     '  -85.0d0,15.0d0,-75.0d0,4.0d-3,8*0.0d0/
C *** FitzHugh-Nagumo w/ recovery (ITYP2=9,ITYP3=2,ITYP5=2)
C ***  default values: conductivities from Trayanova - HPC
C ***                  FHN consts from Rogers in submission
C ***  second constants and PFHN are for Panfilov FHN model only
C *** DPN 14 March 2000 - fixing up units
C      DATA DYNAM_09_2 /'Membrane capacitance Cm (nf/mm^2) ',      !1
C     '                 'Surface-to-volume ratio Am (1/mm) ',      !2
C     '                 'Int. conductivity-fibre (Mho/mm)  ',      !3
C     '                 'Int. conductivity-sheet (Mho/mm)  ',      !4
C     '                 'Int. conductivity-cross (Mho/mm)  ',      !5
C     '                 'Ext. conductivity-fibre (Mho/mm)  ',      !6
C     '                 'Ext. conductivity-sheet (Mho/mm)  ',      !7
C     '                 'Ext. conductivity-cross (Mho/mm)  ',      !8
C     '                 'Rest potential (Vr)      - (mV)   ',      !9
C     '                 'Plateau potential (Vp)   - (mV)   ',      !10
C     '                 'Threshold potential (Vt) - (mV)   ',      !11
C     '                 'Excitation rate (c1 /ms) (k)      ',      !12
C     '                 'Excitation decay (c2 /ms) (k2)    ',      !13
C     '                 'Recovery rate (b /ms) (eps1)      ',      !14
C     '                 'Recovery decay (d) (mu1)          ',      !15
C     '                 'Panfilov FHN (mu2)                ',      !16
C     '                 'Ca2+ Time Const (Panfilov only)   ',      !17
C     '                 'SAC conductance                   ',      !18
C     '                 2*' '/
      DATA DYNAM_09_2 /'Membrane capacitance Cm (uF/mm^2) ',      !1
     '                 'Surface-to-volume ratio Am (1/mm) ',      !2
     '                 'Int. conductivity-fibre (mS/mm)   ',      !3
     '                 'Int. conductivity-sheet (mS/mm)   ',      !4
     '                 'Int. conductivity-cross (mS/mm)   ',      !5
     '                 'Ext. conductivity-fibre (mS/mm)   ',      !6
     '                 'Ext. conductivity-sheet (mS/mm)   ',      !7
     '                 'Ext. conductivity-cross (mS/mm)   ',      !8
     '                 'Rest potential (Vr)      - (mV)   ',      !9
     '                 'Plateau potential (Vp)   - (mV)   ',      !10
     '                 'Threshold potential (Vt) - (mV)   ',      !11
     '                 'Excitation rate (c1 /ms) (k)      ',      !12
     '                 'Excitation decay (c2 /ms) (k2)    ',      !13
     '                 'Recovery rate (b /ms) (eps1)      ',      !14
     '                 'Recovery decay (d) (mu1)          ',      !15
     '                 'Panfilov FHN (mu2)                ',      !16
     '                 'Ca2+ Time Const (Panfilov only)   ',      !17
     '                 'SAC conductance                   ',      !18
     '                 2*' '/
c      DATA DDYNAM_09_2/10.0d0,200.0d0,
c     '  0.174d-3,0.019d-3,0.019d-3,0.625d-3,0.236d-3,0.236d-3,
c     '  -85.0d0,15.0d0,-75.0d0,0.26d0,0.01d0,0.013d0,1.0d0,1.0d0,
c     '  5.0d0,3*0.0d0/
C *** New defaults are for Panfilov FHN
C *** DPN 14 March 2000 - fixing up units
C      DATA DDYNAM_09_2/1.0d6,1.0d0,
C     '  1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,
C     '  -85.0d0,15.0d0,-75.0d0,0.8d1,0.1d1,0.002d0,0.2d0,0.3d0,
C     '  5.0d0,0.0d0,2*0.0d0/
      DATA DDYNAM_09_2/1.0d3,1.0d0,
     '  1.0d3,1.0d3,1.0d3,1.0d3,1.0d3,1.0d3,
     '  -85.0d0,15.0d0,-75.0d0,0.8d1,0.1d1,0.002d0,0.2d0,0.3d0,
     '  5.0d0,0.0d0,2*0.0d0/
C *** van Capelle-Durrer (ITYP2=9,ITYP3=3,ITYP5=2)
C *** VCD consts from Alan Garfinkel correspondence
C *** DPN 14 March 2000 - fixing up units
C      DATA DYNAM_09_3 /'Membrane capacitance Cm (nf/mm^2) ',
C     '                 'Surface-to-volume ratio Am (1/mm) ',
C     '                 'Int. conductivity-fibre (Mho/mm)  ',
C     '                 'Int. conductivity-sheet (Mho/mm)  ',
C     '                 'Int. conductivity-cross (Mho/mm)  ',
C     '                 'Ext. conductivity-fibre (Mho/mm)  ',
C     '                 'Ext. conductivity-sheet (Mho/mm)  ',
C     '                 'Ext. conductivity-cross (Mho/mm)  ',
C     '                 'Initial Resting potential (mV)    ',
C     '                 'Time Constant (Factor) T          ',
C     '                 'Initial Recovery State            ',
C     '                 'Ca2+ Time Const (Panfilov only)   ',
C     '                 'SAC conductance                   ',
C     '                 'Repolarisation time constant (ms) ',
C     '                 6*' '/
C      DATA DDYNAM_09_3/10.0d0,200.0d0,
C     '  0.174d-3,0.019d-3,0.019d-3,0.625d-3,0.236d-3,0.236d-3,
C     '  -75.0d0,1.00d0,0.08d0,5.0d0,0.0d0,0.066d0,6*0.0d0/
      DATA DYNAM_09_3 /'Membrane capacitance Cm (uF/mm^2) ',
     '                 'Surface-to-volume ratio Am (1/mm) ',
     '                 'Int. conductivity-fibre (mS/mm)   ',
     '                 'Int. conductivity-sheet (mS/mm)   ',
     '                 'Int. conductivity-cross (mS/mm)   ',
     '                 'Ext. conductivity-fibre (mS/mm)   ',
     '                 'Ext. conductivity-sheet (mS/mm)   ',
     '                 'Ext. conductivity-cross (mS/mm)   ',
     '                 'Initial Resting potential (mV)    ',
     '                 'Time Constant (Factor) T          ',
     '                 'Initial Recovery State            ',
     '                 'Ca2+ Time Const (Panfilov only)   ',
     '                 'SAC conductance                   ',
     '                 'Repolarisation time constant (ms) ',
     '                 6*' '/
      DATA DDYNAM_09_3/10.0d-3,200.0d0,
     '  0.174d0,0.019d0,0.019d0,0.625d0,0.236d0,0.236d0,
     '  -75.0d0,1.00d0,0.08d0,5.0d0,0.0d0,0.066d0,6*0.0d0/
C *** Beeler-Reuter (ITYP2=9,ITYP3=4,ITYP5=2)
C *** DPN 14 March 2000 - fixing up units
C      DATA DYNAM_09_4 /'Membrane capacitance Cm (nf/mm^2) ',
C     '                 'Surface-to-volume ratio Am (1/mm) ',
C     '                 'Int. conductivity-fibre (Mho/mm)  ',
C     '                 'Int. conductivity-sheet (Mho/mm)  ',
C     '                 'Int. conductivity-cross (Mho/mm)  ',
C     '                 'Ext. conductivity-fibre (Mho/mm)  ',
C     '                 'Ext. conductivity-sheet (Mho/mm)  ',
C     '                 'Ext. conductivity-cross (Mho/mm)  ',
C     '                 'Initial Resting potential (mV)    ',
C     '                 'Sodium Equil Pot (VNa)    (mV)    ',
C     '                 'Sodium Cond. (GNa)    (mMho/mm^2) ',
C     '                 's-s Sodium Cond (GNaC)(mMho/mm^2) ',
C     '                 'Slow Cond. (Gs)       (mMho/mm^2) ',
C     '                 7*' '/
C      DATA DDYNAM_09_4/1.0d0,200.0d0,
C     '  0.174d-3,0.019d-3,0.019d-3,0.625d-3,0.236d-3,0.236d-3,
C     '  -80.d0,50.d0,4.d0,0.003d0,0.09d0,7*0.0d0/
      DATA DYNAM_09_4 /'Membrane capacitance Cm (uF/mm^2) ',
     '                 'Surface-to-volume ratio Am (1/mm) ',
     '                 'Int. conductivity-fibre (mS/mm)   ',
     '                 'Int. conductivity-sheet (mS/mm)   ',
     '                 'Int. conductivity-cross (mS/mm)   ',
     '                 'Ext. conductivity-fibre (mS/mm)   ',
     '                 'Ext. conductivity-sheet (mS/mm)   ',
     '                 'Ext. conductivity-cross (mS/mm)   ',
     '                 'Initial Resting potential (mV)    ',
     '                 'Sodium Equil Pot (VNa)    (mV)    ',
     '                 'Sodium Cond. (GNa)    (mS/mm^2)   ',
     '                 's-s Sodium Cond (GNaC)(mS/mm^2)   ',
     '                 'Slow Cond. (Gs)       (mS/mm^2)   ',
     '                 7*' '/
      DATA DDYNAM_09_4/1.0d-3,200.0d0,
     '  0.174d0,0.019d0,0.019d0,0.625d0,0.236d0,0.236d0,
     '  -80.d0,50.d0,4.d0,0.003d0,0.09d0,7*0.0d0/
C *** Jafri-Rice-Winslow
C *** DPN 14 March 2000 - fixing up units
C      DATA DYNAM_09_5 /'Membrane capacitance Cm (nf/mm^2) ',
C     '                 'Surface-to-volume ratio Am (1/mm) ',
C     '                 'Int. conductivity-fibre (Mho/mm)  ',
C     '                 'Int. conductivity-sheet (Mho/mm)  ',
C     '                 'Int. conductivity-cross (Mho/mm)  ',
C     '                 'Ext. conductivity-fibre (Mho/mm)  ',
C     '                 'Ext. conductivity-sheet (Mho/mm)  ',
C     '                 'Ext. conductivity-cross (Mho/mm)  ',
C     '                 12*' '/
C      DATA DDYNAM_09_5/10.0d0,200.0d0,
C     '  0.174d-3,0.019d-3,0.019d-3,0.625d-3,0.236d-3,0.236d-3,
C     '  12*0.0d0/
      DATA DYNAM_09_5 /'Membrane capacitance Cm (uF/mm^2) ',
     '                 'Surface-to-volume ratio Am (1/mm) ',
     '                 'Int. conductivity-fibre (mS/mm)   ',
     '                 'Int. conductivity-sheet (mS/mm)   ',
     '                 'Int. conductivity-cross (mS/mm)   ',
     '                 'Ext. conductivity-fibre (mS/mm)   ',
     '                 'Ext. conductivity-sheet (mS/mm)   ',
     '                 'Ext. conductivity-cross (mS/mm)   ',
     '                 12*' '/
      DATA DDYNAM_09_5/10.0d-3,200.0d0,
     '  0.174d0,0.019d0,0.019d0,0.625d0,0.236d0,0.236d0,
     '  12*0.0d0/
C *** Luo-Rudy (ITYP2=9,ITYP3=6,ITYP5=2)
C *** DPN 14 March 2000 - fixing up units
C      DATA DYNAM_09_6 /'Membrane capacitance Cm (nf/mm^2) ',
C     '                 'Surface-to-volume ratio Am (1/mm) ',
C     '                 'Int. conductivity-fibre (Mho/mm)  ',
C     '                 'Int. conductivity-sheet (Mho/mm)  ',
C     '                 'Int. conductivity-cross (Mho/mm)  ',
C     '                 'Ext. conductivity-fibre (Mho/mm)  ',
C     '                 'Ext. conductivity-sheet (Mho/mm)  ',
C     '                 'Ext. conductivity-cross (Mho/mm)  ',
C     '                 12*' '/
CC *** Additionally, positions 17:20 contain the following computed
CC     parameters: GK, GK1, VK and VK1, respectively.
C      DATA DDYNAM_09_6/10.0d0,200.0d0,
C     '  0.174d-3,0.019d-3,0.019d-3,0.625d-3,0.236d-3,0.236d-3,
C     '  12*0.0d0/
      DATA DYNAM_09_6 /'Membrane capacitance Cm (uF/mm^2) ',
     '                 'Surface-to-volume ratio Am (1/mm) ',
     '                 'Int. conductivity-fibre (mS/mm)   ',
     '                 'Int. conductivity-sheet (mS/mm)   ',
     '                 'Int. conductivity-cross (mS/mm)   ',
     '                 'Ext. conductivity-fibre (mS/mm)   ',
     '                 'Ext. conductivity-sheet (mS/mm)   ',
     '                 'Ext. conductivity-cross (mS/mm)   ',
     '                 12*' '/
      DATA DDYNAM_09_6/10.0d-3,200.0d0,
     '  0.174d0,0.019d0,0.019d0,0.625d0,0.236d0,0.236d0,
     '  12*0.0d0/
C *** diFrancesco-Noble (ITYP2=9,ITYP3=7,ITYP5=2)
C *** DPN 14 March 2000 - fixing up units
C      DATA DYNAM_09_7 /'Membrane capacitance Cm (nf/mm^2) ',
C     '                 'Surface-to-volume ratio Am (1/mm) ',
C     '                 'Int. conductivity-fibre (Mho/mm)  ',
C     '                 'Int. conductivity-sheet (Mho/mm)  ',
C     '                 'Int. conductivity-cross (Mho/mm)  ',
C     '                 'Ext. conductivity-fibre (Mho/mm)  ',
C     '                 'Ext. conductivity-sheet (Mho/mm)  ',
C     '                 'Ext. conductivity-cross (Mho/mm)  ',
C     '                 'Initial Resting potential (mV)    ',
C     '                 '1/Membrane Res. (GL)  (uMho/mm^2) ',
C     '                 'Sodium Cond. (GNa)    (mMho/mm^2) ',
C     '                 'Sodium Equil Pot (VNa)    (mV)    ',
C     '                 'Initial value of m                ',
C     '                 'Initial value of h                ',
C     '                 6*' '/
C      DATA DDYNAM_09_7/10.0d0,200.0d0,
C     '  0.174d-3,0.019d-3,0.019d-3,0.625d-3,0.236d-3,0.236d-3,
C     '  -80.d0,0.5d0,0.28d0,33.45d0,0.d0,0.18d0,6*0.0d0/
      DATA DYNAM_09_7 /'Membrane capacitance Cm (uF/mm^2) ',
     '                 'Surface-to-volume ratio Am (1/mm) ',
     '                 'Int. conductivity-fibre (mS/mm)   ',
     '                 'Int. conductivity-sheet (mS/mm)   ',
     '                 'Int. conductivity-cross (mS/mm)   ',
     '                 'Ext. conductivity-fibre (mS/mm)   ',
     '                 'Ext. conductivity-sheet (mS/mm)   ',
     '                 'Ext. conductivity-cross (mS/mm)   ',
     '                 'Initial Resting potential (mV)    ',
     '                 '1/Membrane Res. (GL)  (uMho/mm^2) ',
     '                 'Sodium Cond. (GNa)    (mMho/mm^2) ',
     '                 'Sodium Equil Pot (VNa)    (mV)    ',
     '                 'Initial value of m                ',
     '                 'Initial value of h                ',
     '                 6*' '/
      DATA DDYNAM_09_7/10.0d-3,200.0d0,
     '  0.174d0,0.019d0,0.019d0,0.625d0,0.236d0,0.236d0,
     '  -80.d0,0.5d0,0.28d0,33.45d0,0.d0,0.18d0,6*0.0d0/
C *** Noble 98 (ITYP2=9,ITYP3=8,ITYP5=2)
C *** DPN 14 March 2000 - fixing up units
C      DATA DYNAM_09_8 /'Membrane capacitance Cm (nf/mm^2) ',
C     '                 'Surface-to-volume ratio Am (1/mm) ',
C     '                 'Int. conductivity-fibre (Mho/mm)  ',
C     '                 'Int. conductivity-sheet (Mho/mm)  ',
C     '                 'Int. conductivity-cross (Mho/mm)  ',
C     '                 'Ext. conductivity-fibre (Mho/mm)  ',
C     '                 'Ext. conductivity-sheet (Mho/mm)  ',
C     '                 'Ext. conductivity-cross (Mho/mm)  ',
C     '                 'Initial Resting potential (mV)    ',
C     '                 'Initial intracellular sodium      ',
C     '                 10*' '/
C      DATA DDYNAM_09_8 /10.0d0,200.0d0,0.174d-3,0.019d-3,0.019d-3,
C     '  0.625d-3,0.236d-3,0.236d-3,-93.0d0,5.6d0,10*0.0d0/
      DATA DYNAM_09_8 /'Membrane capacitance Cm (uF/mm^2) ',
     '                 'Surface-to-volume ratio Am (1/mm) ',
     '                 'Int. conductivity-fibre (mS/mm)   ',
     '                 'Int. conductivity-sheet (mS/mm)   ',
     '                 'Int. conductivity-cross (mS/mm)   ',
     '                 'Ext. conductivity-fibre (mS/mm)   ',
     '                 'Ext. conductivity-sheet (mS/mm)   ',
     '                 'Ext. conductivity-cross (mS/mm)   ',
     '                 'Initial Resting potential (mV)    ',
     '                 'Initial intracellular sodium      ',
     '                 10*' '/
      DATA DDYNAM_09_8 /10.0d-3,200.0d0,0.174d0,0.019d0,0.019d0,
     '  0.625d0,0.236d0,0.236d0,-93.0d0,5.6d0,10*0.0d0/
C *** Hodgkin-Huxley (ITYP2=9,ITYP3=9,ITYP5=2)
C *** DPN 14 March 2000 - fixing up units
C      DATA DYNAM_09_9 /'Membrane capacitance Cm (nf/mm^2) ',
C     '                 'Surface-to-volume ratio Am (1/mm) ',
C     '                 'Int. conductivity-fibre (Mho/mm)  ',
C     '                 'Int. conductivity-sheet (Mho/mm)  ',
C     '                 'Int. conductivity-cross (Mho/mm)  ',
C     '                 'Ext. conductivity-fibre (Mho/mm)  ',
C     '                 'Ext. conductivity-sheet (Mho/mm)  ',
C     '                 'Ext. conductivity-cross (Mho/mm)  ',
C     '                 12*' '/
C      DATA DDYNAM_09_9/10.0d0,200.0d0,
C     '  0.174d-3,0.019d-3,0.019d-3,0.625d-3,0.236d-3,0.236d-3,
C     '  12*0.0d0/
      DATA DYNAM_09_9 /'Membrane capacitance Cm (uF/mm^2) ',
     '                 'Surface-to-volume ratio Am (1/mm) ',
     '                 'Int. conductivity-fibre (mS/mm)  ',
     '                 'Int. conductivity-sheet (mS/mm)  ',
     '                 'Int. conductivity-cross (mS/mm)  ',
     '                 'Ext. conductivity-fibre (mS/mm)  ',
     '                 'Ext. conductivity-sheet (mS/mm)  ',
     '                 'Ext. conductivity-cross (mS/mm)  ',
     '                 12*' '/
      DATA DDYNAM_09_9/10.0d-3,200.0d0,
     '  0.174d0,0.019d0,0.019d0,0.625d0,0.236d0,0.236d0,
     '  12*0.0d0/
C *** User defined (ITYP2=9,ITYP3=10,ITYP5=2)
C *** DPN 14 March 2000 - fixing up units
C      DATA DYNAM_09_10 /'Membrane capacitance Cm (uf/mm^2) ',
C     '                  'Surface-to-volume ratio Am (1/mm) ',
C     '                  'Int. conductivity-fibre (mS/mm)  ',
C     '                  'Int. conductivity-sheet (mS/mm)  ',
C     '                  'Int. conductivity-cross (mS/mm)  ',
C     '                  'Ext. conductivity-fibre (mS/mm)  ',
C     '                  'Ext. conductivity-sheet (mS/mm)  ',
C     '                  'Ext. conductivity-cross (mS/mm)  ',
C     '                  12*' '/
C      DATA DDYNAM_09_10 /0.01d0,200.0d0,0.3d0,0.031525d0,0.06305d0,
C     '  0.2d0,0.13514d0,0.1756821d0,12*0.0d0/
      DATA DYNAM_09_10 /'Membrane capacitance Cm (uF/mm^2) ',
     '                  'Surface-to-volume ratio Am (1/mm) ',
     '                  'Int. conductivity-fibre (mS/mm)  ',
     '                  'Int. conductivity-sheet (mS/mm)  ',
     '                  'Int. conductivity-cross (mS/mm)  ',
     '                  'Ext. conductivity-fibre (mS/mm)  ',
     '                  'Ext. conductivity-sheet (mS/mm)  ',
     '                  'Ext. conductivity-cross (mS/mm)  ',
     '                  12*' '/
      DATA DDYNAM_09_10 /0.01d0,200.0d0,0.3d0,0.031525d0,0.06305d0,
     '  0.2d0,0.13514d0,0.1756821d0,12*0.0d0/
C *** Multi-field oxygen transport(ITYP2=10,ITYP3=1,ITYP5=2)
      DATA DYNAM_10_1 /'Oxy. m. solubility (nmol/kPa/mm^3)',
     '                 'Oxy. b. solubility (nmol/kPa/mm^3)',
     '                 'Oxy. m. diffusivity (mm^2/ks)     ',
     '                 'Oxy. b. diffusivity (mm^2/ks)     ',
     '                 'Myoglobin diffusivity (mm^2/ks)   ',
     '                 'Max myoglobin   conc. (nmol/mm^3) ',
     '                 'Max haemoglobin conc. (nmol/mm^3) ',
     '                 'Max mitochon. rate (nmol/ks/mm^3) ',
     '                 'Exchange coefficient (1/ks)       ',
     '                 'Arterial blood oxygen p.p. (kPa)  ',10*' '/
c    '                 'Relative volume of vascular bed   ',
      DATA DDYNAM_10_1/0.01125d0,
     '                 0.0093d0,
     '                 2.3d0,
     '                 1.5d0,
     '                 0.07d0,
     '                 0.28d0,
     '                 8.8d0,
     '                 20.0d0,
     '                 0.1d0,
     '                 100.0d0,
     '                 10*0.0d0/
C *** Glucose-Oxygen transport (ITYP2=10,ITYP3=2,ITYP5=2)
      DATA DYNAM_10_2 /'Oxygen  solubility (nmol/kPa/mm^3)',
     '                 'Glucose solubility (nmol/kPa/mm^3)',
     '                 'Oxygen  diffusivity (mm^2/ks)     ',
     '                 'Glucose diffusivity (mm^2/ks)     ',
     '                 'Oxygen  Michaelis p.pressure (kPa)',
     '                 'Glucose Michaelis p.pressure (kPa)',
     '                 'Vmax (nmol/ks/mm^3)               ',
     '                 'Oxygen reaction weight (nu)       ',12*' '/
      DATA DDYNAM_10_2/20*0.0d0/
C *** Lung gas transport (ITYP2=5,ITYP3=2,ITYP5=2)
      DATA DYNAM_11_1 /'Gas diffusion coefficient (mm^2/s)',
     '                 'Oxygen loss coefficient           ',18*' '/
      DATA DDYNAM_11_1/74.3d0,0.d0,18*0.0d0/
C *** Lung temperature/humidity transport (ITYP2=11,ITYP5=2)
      DATA DYNAM_11_2 /'Water vapour diffusivity (mm^2/s) ',
     '                 'Thermal diffusivity (mm^2/s)      ',
     '                 'Density of inspired gas (g/mm^3)  ',
     '                 'Specific heat of gas (kJ/g/K)     ',
     '                 'Thermal conductivity (kJ/s/mm/K)  ',
     '                 'Heat of evaporation (kJ/g)        ',14*' '/
      DATA DDYNAM_11_2/27.7d0,23.7d0,1.1d-6,1.08d-3,0.268d-7,2.411d0,
     '  14*0.0d0/
C *** Pulmonary Circulation   (ITYP2=11,ITYP3=3,ITYP5=2)
      DATA DYNAM_11_3 /'Flux cutoff parameter (r)         ',
     '                 'Preferential flux parameter (b)   ',
     '                 'Capillary wall stiffness (cmH20)  ',17*' '/
      DATA DDYNAM_11_3/0.05d0,1.15d0,22.5d0,17*0.0d0/
C *** Simple pulmonary pressure-resistance-flow model (ITYP2=11,ITYP3=4,ITYP5=2)
      DATA DYNAM_11_4 /'Density                               ',
     '                 'Viscosity                             ',
     '                 'elasticity parameter Go               ',
     '                 'elasticity parameter beta             ',
     '                 'elasticity parameter beta2            ',
     '                 'elasticity parameter Gov              ',
     '                 'elasticity parameter betav            ',
     '                 'Dummy                                 ',12*' '/
      DATA DDYNAM_11_4/1.05d-3,3.36d-3,6.7d0,1.2d0,0.d0,6.7d0,1.2d0,
     '    1.d0,12*0.0d0/
C *** Pulmonary Circulation (arteriole-cap-venule)   (ITYP2=11,ITYP3=6,ITYP5=2)
      DATA DYNAM_11_6 /'Flux cutoff parameter (r)         ',
     '                 'Preferential flux parameter (b)   ',
     '                 'Vessel compliance (cmH20)  ',17*' '/
      DATA DDYNAM_11_6/0.05d0,1.15d0,22.5d0,17*0.0d0/
C *** Pulmonary gas exchange with pre-defined transport(ITYP2=11,ITYP3=7)
      DATA DYNAM_11_7 /'atmospheric total pressure, mmHg  ',
     '                 'unit elastance E, mmHg/litre      ', 
     '                 'unit capillary blood hematocrit   ',
     '                 'unit capillary blood volume, litre',
     '                 'unit alveolar air volume, litre   ',
     '                 'unit air-blood surface area, m2   ',
     '                 'unit air-blood barrier thickness,m',
     '                 13*' '/
      DATA DDYNAM_11_7/760.0d0,
     '                 2.5d0,
     '                 0.4d0,
     '                 0.58333d-6,
     '                 0.187d-3,
     '                 0.11667d-2,
     '                 1.1d-6,13*0.0d0/

C *** Modal analysis ***   (ityp5=3)
C **********************

      DATA ILTOT3/ !(njt,ityp2,ityp3) are #s of parameters
     '  0,0,0, 0,0,0, 0,0,0, 1,1,1, 0,0,0, 1,1,1, 7,7,7, 0,0,0, 9,9,9,
     '  0,0,0, 0,0,0, 0,0,0, 1,1,1, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 8,8,8,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0/

      CHARACTER
     '  MODAL_01_1(20)*34,                    !Linear elasticity
     '  MODAL_02_1(20)*34,                    !Finite elasticity
     '  MODAL_03_1(20)*34,                    !Laplace's equation
     '  MODAL_04_1(20)*34,MODAL_04_2(20)*34,  !Helmholtz eqtn, Yukawa eqtn
     '  MODAL_05_1(20)*34,                    !div(kgrad(u))=f
     '  MODAL_06_1(20)*34,                    !Linear 2nd order elliptic
     '  MODAL_07_1(20)*34,                    !Biharmonic equation
     '  MODAL_08_1(20)*34,                    !Vocal tract equations
     '  MODAL_09_1(20)*34,
     '  MODAL_10_1(20)*34,
     '  MODAL_11_1(20)*34,
     '  MODAL_12_1(20)*34

      EQUIVALENCE
     '  (TITL33(1, 1,1),MODAL_01_1),
     '  (TITL33(1, 2,1),MODAL_02_1),
     '  (TITL33(1, 3,1),MODAL_03_1),
     '  (TITL33(1, 4,1),MODAL_04_1),(TITL33(1, 4,2),MODAL_04_2),
     '  (TITL33(1, 5,1),MODAL_05_1),
     '  (TITL33(1, 6,1),MODAL_06_1),
     '  (TITL33(1, 7,1),MODAL_07_1),
     '  (TITL33(1, 8,1),MODAL_08_1),
     '  (TITL33(1, 9,1),MODAL_09_1),
     '  (TITL33(1,10,1),MODAL_10_1),
     '  (TITL33(1,11,1),MODAL_11_1),
     '  (TITL33(1,12,1),MODAL_12_1)

      REAL*8
     '  DMODAL_01_1(20),
     '  DMODAL_02_1(20),
     '  DMODAL_03_1(20),
     '  DMODAL_04_1(20),DMODAL_04_2(20),
     '  DMODAL_05_1(20),
     '  DMODAL_06_1(20),
     '  DMODAL_07_1(20),
     '  DMODAL_08_1(20),
     '  DMODAL_09_1(20),
     '  DMODAL_10_1(20),
     '  DMODAL_11_1(20),
     '  DMODAL_12_1(20)

      EQUIVALENCE
     '  (RDMATE(1,3, 1,1),DMODAL_01_1(1)),
     '  (RDMATE(1,3, 2,1),DMODAL_02_1(1)),
     '  (RDMATE(1,3, 3,1),DMODAL_03_1(1)),
     '  (RDMATE(1,3, 4,1),DMODAL_04_1(1)),
     '  (RDMATE(1,3, 4,2),DMODAL_04_2(1)),
     '  (RDMATE(1,3, 5,1),DMODAL_05_1(1)),
     '  (RDMATE(1,3, 6,1),DMODAL_06_1(1)),
     '  (RDMATE(1,3, 7,1),DMODAL_07_1(1)),
     '  (RDMATE(1,3, 8,1),DMODAL_08_1(1)),
     '  (RDMATE(1,3, 9,1),DMODAL_09_1(1)),
     '  (RDMATE(1,3,10,1),DMODAL_10_1(1)),
     '  (RDMATE(1,3,11,1),DMODAL_11_1(1)),
     '  (RDMATE(1,3,12,1),DMODAL_12_1(1))

C *** Laplace's equation (ITYP2=3,ITYP3=1,ITYP5=3)
      DATA MODAL_03_1/20*' '/
      DATA DMODAL_03_1/20*0.0d0/
C *** Helmholtz equation (ITYP2=4,ITYP3=1,ITYP5=3)
      DATA MODAL_04_1/'Wave number k',19*' '/
      DATA DMODAL_04_1/1.0d0,19*0.0d0/
C *** Yukawa equation (ITYP2=4,ITYP3=2,ITYP5=3)
      DATA MODAL_04_2/'s value',19*' '/
      DATA DMODAL_04_2/1.0d0,19*0.0d0/
C *** div(kgrad(u))=f (ITYP2=5,ITYP3=1,ITYP5=3)
      DATA MODAL_05_1/'source term',19*' '/
      DATA DMODAL_05_1/20*0.0d0/
C *** Linear 2nd order elliptic (ITYP2=6,ITYP3=1,ITYP5=3)
      DATA MODAL_06_1/20*' '/
      DATA DMODAL_06_1/20*0.0d0/
C *** Biharmonic equation (ITYP2=7,ITYP3=1,ITYP5=3)
      DATA MODAL_07_1/'Source term',19*' '/
      DATA DMODAL_07_1/20*0.0d0/
C *** Vocal tract optimisation
      DATA MODAL_08_1/'Cross-section area',
     '                'Wave speed',18*' '/
      DATA DMODAL_08_1/20*0.0d0/


C *** Quasi-static analysis ***   (ityp5=4)
C *****************************

C gbs 9/9/95 copied iltot1
! ityp2= 1      2      3      4      5      6      7      8       9
      DATA ILTOT4/ !(njt,ityp2,ityp3) are #s of parameters
     ' 0,0,0, 0,0,0, 0,0,0, 1,1,1, 2,3,4, 0,0,0, 1,1,1, 2,2,2, 15,15,15,
     ' 0,0,0, 0,0,0, 1,2,3, 1,1,1, 2,2,2, 0,0,0, 0,0,0, 2,2,2, 8,8,8,
     ' 0,0,0, 0,0,0, 2,2,2, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     ' 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     ' 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0/
c      DATA ILTOT4/ !(njt,ityp2,ityp3) are #s of parameters
c     '  2,2,2, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
c     '  0,0,0, 0,0,0, 1,2,3, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
c     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
c     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
c     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0/

C KAT 2001-04-09: TITL31 is used not TITL34
C      CHARACTER
C     '  QUASIST_01_1(20)*34,
C     '  QUASIST_03_1(20)*34,QUASIST_03_2(20)*34,
C     '  QUASIST_05_1(20)*34,QUASIST_05_2(20)*34

C      EQUIVALENCE
C     '  (TITL34(1, 1,1),QUASIST_01_1),
C     '  (TITL34(1, 3,1),QUASIST_03_1),(TITL34(1, 3,2),QUASIST_03_2),
C     '  (TITL34(1, 5,1),QUASIST_05_1),(TITL34(1, 5,2),QUASIST_05_2)

      REAL*8
     '  DQUASIST_01_1(20),
     '  DQUASIST_02_1(20),
     '  DQUASIST_03_1(20),DQUASIST_03_2(20),
     '  DQUASIST_04_1(20),
     '  DQUASIST_05_1(20),DQUASIST_05_2(20),
     '  DQUASIST_06_1(20),
     '  DQUASIST_07_1(20),
     '  DQUASIST_08_1(20),
     '  DQUASIST_09_1(20),
     '  DQUASIST_10_1(20),
     '  DQUASIST_11_1(20),
     '  DQUASIST_12_1(20)

      EQUIVALENCE
     '  (RDMATE(1,4, 1,1),DQUASIST_01_1(1)),
     '  (RDMATE(1,4, 2,1),DQUASIST_02_1(1)),
     '  (RDMATE(1,4, 3,1),DQUASIST_03_1(1)),
     '  (RDMATE(1,4, 3,2),DQUASIST_03_2(1)),
     '  (RDMATE(1,4, 4,1),DQUASIST_04_1(1)),
     '  (RDMATE(1,4, 5,1),DQUASIST_05_1(1)),
     '  (RDMATE(1,4, 5,2),DQUASIST_05_2(1)),
     '  (RDMATE(1,4, 6,1),DQUASIST_06_1(1)),
     '  (RDMATE(1,4, 7,1),DQUASIST_07_1(1)),
     '  (RDMATE(1,4, 8,1),DQUASIST_08_1(1)),
     '  (RDMATE(1,4, 9,1),DQUASIST_09_1(1)),
     '  (RDMATE(1,4,10,1),DQUASIST_10_1(1)),
     '  (RDMATE(1,4,11,1),DQUASIST_11_1(1)),
     '  (RDMATE(1,4,12,1),DQUASIST_12_1(1))

C *** Laplace's equation           (ITYP2=3,ITYP3=1,ITYP5=4)
C      DATA QUASIST_03_1/20*' '/
      DATA DQUASIST_03_1/20*0.0d0/
C *** Generalised Laplace        **(ITYP2=3,ITYP3=2,ITYP5=4)
C      DATA QUASIST_03_2/'conductivity s11',
C     '                  'conductivity s22',
C     '                  'conductivity s33',
C     '                  17*' '/
      DATA DQUASIST_03_2/1.0d0,1.0d0,1.0d0,17*0.0d0/
C *** div(kgrad(u))=f (f General)**(ITYP2=5,ITYP3=1,ITYP5=1)
C      DATA QUASIST_05_1/'domain source f',
C     '                 'permeability k11',
C     '                 'permeability k22',
C     '                 'permeability k33',16*' '/
      DATA DQUASIST_05_1/0.0d0,1.0d0,1.0d0,1.0d0,16*0.0d0/
C *** div(kgrad(u))=f, f=-div(s.grad(g)) **(ITYP2=5,ITYP3=2,ITYP5=1)
C      DATA QUASIST_05_2/'permeability k',
C     '                 'permeability s',18*' '/
      DATA DQUASIST_05_2/1.0d0,1.0d0,18*0.0d0/


C *** Wavefront path analysis ***   (ityp5=5)
C *******************************

      DATA ILTOT5/ !(njt,ityp2,ityp3) are #s of parameters
     '  0,0,0, 0,0,0, 3,3,3, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     '  0,0,0, 0,0,0, 3,4,5, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     '  0,0,0, 0,0,0, 4,6,8, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0/

      CHARACTER FRONT_03_1(20)*34
      CHARACTER FRONT_03_2(20)*34
      CHARACTER FRONT_03_3(20)*34

      EQUIVALENCE (TITL35(1, 3,1),FRONT_03_1)
      EQUIVALENCE (TITL35(1, 3,2),FRONT_03_2)
      EQUIVALENCE (TITL35(1, 3,3),FRONT_03_3)

      REAL*8
     '  DFRONT_01_1(20),
     '  DFRONT_02_1(20),
     '  DFRONT_03_1(20),DFRONT_03_2(20),DFRONT_03_3(20),
     '  DFRONT_04_1(20),
     '  DFRONT_05_1(20),
     '  DFRONT_06_1(20),
     '  DFRONT_07_1(20),
     '  DFRONT_08_1(20),
     '  DFRONT_09_1(20),
     '  DFRONT_10_1(20),
     '  DFRONT_11_1(20),
     '  DFRONT_12_1(20)

      EQUIVALENCE
     '  (RDMATE(1,5, 1,1),DFRONT_01_1(1)),
     '  (RDMATE(1,5, 2,1),DFRONT_02_1(1)),
     '  (RDMATE(1,5, 3,1),DFRONT_03_1(1)),
     '    (RDMATE(1,5, 3,2),DFRONT_03_2(1)),
     '    (RDMATE(1,5, 3,3),DFRONT_03_3(1)),
     '  (RDMATE(1,5, 4,1),DFRONT_04_1(1)),
     '  (RDMATE(1,5, 5,1),DFRONT_05_1(1)),
     '  (RDMATE(1,5, 6,1),DFRONT_06_1(1)),
     '  (RDMATE(1,5, 7,1),DFRONT_07_1(1)),
     '  (RDMATE(1,5, 8,1),DFRONT_08_1(1)),
     '  (RDMATE(1,5, 9,1),DFRONT_09_1(1)),
     '  (RDMATE(1,5,10,1),DFRONT_10_1(1)),
     '  (RDMATE(1,5,11,1),DFRONT_11_1(1)),
     '  (RDMATE(1,5,12,1),DFRONT_12_1(1))

C     Homogeneous
      DATA FRONT_03_1/'Time constant source term',
     '                'Dimensionless coef. of advection  ',
     '                'Space constant                    ',
     '                17*' '/
      DATA DFRONT_03_1/0.9d0,1.8d0,0.13d0,17*0.0d0/

C     Monodomain
      DATA FRONT_03_2/'Time constant source term',
     '                'Dimensionless coef. of advection  ',
     '                'Space constant along fibres       ',
     '                'Sheet cross fibre space constant  ',
     '                'Cross sheet space constant        ',
     '                15*' '/
      DATA DFRONT_03_2/0.9d0,1.8d0,0.13d0,0.03d0,0.03d0,15*0.0d0/

C     Bidomain
      DATA FRONT_03_3/'Time constant source term',
     '                'Dimensionless coef. of advection  ',
     '                'Fibre dirn. intracellular coupling',
     '                'Fibre dirn. extracellular coupling',
     '                'Sheet cross fibre intra. coupling ',
     '                'Sheet cross fibre extra. coupling ',
     '                'Cross sheet intra. coupling       ',
     '                'Cross sheet extra. coupling       ',
     '                12*' '/
      DATA DFRONT_03_3/0.9d0,1.8d0,
     '  0.19d0,0.39d0,0.04d0,0.19d0,0.04d0,0.19d0,12*0.d0/


C *** Buckling analysis ***
C *************************

      DATA ILTOT6/ !(njt,ityp2,ityp3) are #s of parameters
     '  2,2,2, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
     '  0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0/

      CHARACTER BUCKLE_01_1(20)*34

      EQUIVALENCE (TITL36(1, 1,1),BUCKLE_01_1)

C      REAL*8
C     '  DBUCKLE_01_1(20)
C     '  DBUCKLE_02_1(20),
C     '  DBUCKLE_03_1(20),
C     '  DBUCKLE_04_1(20),
C     '  DBUCKLE_05_1(20),
C     '  DBUCKLE_06_1(20),
C     '  DBUCKLE_07_1(20),
C     '  DBUCKLE_08_1(20),
C     '  DBUCKLE_09_1(20),
C     '  DBUCKLE_10_1(20),
C     '  DBUCKLE_11_1(20),
C     '  DBUCKLE_12_1(20)

C      EQUIVALENCE
C     '  (RDMATE(1,6, 1,1),DBUCKLE_01_1(1)),
C     '  (RDMATE(1,6, 2,1),DBUCKLE_02_1(1)),
C     '  (RDMATE(1,6, 3,1),DBUCKLE_03_1(1)),
C     '  (RDMATE(1,6, 4,1),DBUCKLE_04_1(1)),
C     '  (RDMATE(1,6, 5,1),DBUCKLE_05_1(1)),
C     '  (RDMATE(1,6, 6,1),DBUCKLE_06_1(1)),
C     '  (RDMATE(1,6, 7,1),DBUCKLE_07_1(1)),
C     '  (RDMATE(1,6, 8,1),DBUCKLE_08_1(1)),
C     '  (RDMATE(1,6, 9,1),DBUCKLE_09_1(1)),
C     '  (RDMATE(1,6,10,1),DBUCKLE_10_1(1)),
C     '  (RDMATE(1,6,11,1),DBUCKLE_11_1(1)),
C     '  (RDMATE(1,6,12,1),DBUCKLE_12_1(1))

      DATA BUCKLE_01_1/20*'*** Unknown material constant *** '/

      END



