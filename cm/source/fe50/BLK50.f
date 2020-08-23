      BLOCK DATA BLK50

C#### BlockData: BLK50
C###  Description:
C###    BLK50 sets up parameters in common block

      IMPLICIT NONE
      INCLUDE 'b13.cmn'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'titl50.cmn'

      DATA TITL52 /'Compressible  ','Incompressible',
     '  'Incomp + fluid','Compr + fluid','Incomp + inext'/
      DATA TITL53 /'Reference (theta) coordinates            ',
     '             'Body (fibre/transverse) coordinates      ',
     '             'Body coordinates with active fibre stress'/
C      DATA TITL54 /'a Green strain energy function (hyperelastic)',
C     '             'a Cauchy-type dependence on deformation grads',
C     '             'a quasi-static Creep law                     '/
      DATA TITL55 /'the principal strain invariants',
     '             'the principal extension ratios ',
     '             'the fibre & transverse strains '/
      DATA TITL56 /
C **  Compressible(ktyp52(nr)=1), principle strain invars(ktyp55(nr)=1)

! Alex 19-Nov-02: replaced old Blatz-Ko function
C $$  W=0.5*u*a*((J1-3)+B*(J3**(-2/B)-1))+0.5*(1-a)*((J2-3)+B*(J3**(2/B)
C $$   -1))

     '             'Polynomial in the three strain invariants I1,I2,I3',
     '             'Blatz-Ko material                                 ',
     '             'Unused                                            ',
     '             'Unused                                            ',
     '             'Defined in subroutine USER51                      ',
C **  Incompress(ktyp52(nr)=2,3 or 5), princ strain invars(ktyp55(nr)=1)
     '             'Polynomial in I1,I2 (eg Neo-Hookian,Mooney-Rivlin)',
     '             'Exponential: a(exp(b(I1-3))-1) (biological tissue)',
     '             'Unused                                            ',
     '             'Unused                                            ',
     '             'Defined in subroutine USER51                      ',
C **  Compressible(ktyp52(nr)=1), principle extensions(ktyp55(nr)=2)
     '             'Polynomial in the three extension ratios L1,L2,L3 ',
     '             'Unused                                            ',
     '             'Unused                                            ',
     '             'Unused                                            ',
     '             'Defined in subroutine USER52                      ',
C **  Incompressible(ktyp52(nr)=2,3 or 5), principle extens(ktyp55(nr)=2)
     '             'Polynomial in the first two extension ratios L1,L2',
     '             'Ogden-type dependent on the extension ratios L1,L2',
     '             'Unused                                            ',
     '             'Unused                                            ',
     '             'Defined in subroutine USER52                      ',
C **  Compressible(ktyp52(nr)=1), fibre strains(ktyp55(nr)=3)
     '             'Polynomial  in fibre and transverse normal strains',
     '             'Exponential in fibre and transverse strains       ',
     '             'Linear stress-strain                              ',
     '             'Unused                                            ',
     '             'Defined in subroutine USER53                      ',
C **  Incompressible(ktyp52(nr)=2,3 or 5), fibre strains(ktyp55(nr)=3)
     '             'Polynomial  in fibre and transverse normal strains',
     '             'Exponential in fibre and transverse strains       ',
     '             'Pole-zero in fibre and transverse strains         ',
     '             'Fibre in fluid                                    ',
     '             'Defined in subroutine USER53                      '/
      DATA TITL57 /'I','L','E'/
C      DATA TITL58 /'Sarcomere length (SL) in the undeformed state (um)',
C     '             'Maximum peak isometric twitch tension (kPa)       ',
C     '             'Maximum peak intracellular Ca concentration (uM)  ',
C     '             'Exp. coeff. of SL-dependent Ca sensitivity (1/um) ',
C     '             'Sarcomere length for zero active tension (um)     ',
C     '             'Time to peak twitch tension (s)                   ',
C     '             'SL coefficient of relaxation time (s/um)          ',
C     '             'Constant coefficient of relaxation time (s)       '/

C             *************Default numbers of mat.l params*************
      DATA IMAT5/               !IMAT5(ktyp56(nr),ktyp55(nr),ktyp52(nr))
C **  Compressible  (ktyp52(nr)=1), princ strain invars   (ktyp55(nr)=1)
     '  1,1,1,1,1,1,
C **  Compressible  (ktyp52(nr)=1), principle extensions  (ktyp55(nr)=2)
     '  1,1,1,1,1,1,
C **  Compressible  (ktyp52(nr)=1), fibre strains         (ktyp55(nr)=3)
     '  1,1,1,1,1,1,
C **  Incompressible(ktyp52(nr)=2), princ strain invars   (ktyp55(nr)=1)
     '  1,1,1,1,1,1,
C **  Incompressible(ktyp52(nr)=2), principle extensions  (ktyp55(nr)=2)
     '  1,1,1,1,1,1,
C **  Incompressible(ktyp52(nr)=2), fibre strains         (ktyp55(nr)=3)
     '  1,1,16,1,1,1,
C **  Incomp + fluid(ktyp52(nr)=3), princ strain invars   (ktyp55(nr)=1)
     '  1,1,1,1,1,1,
C **  Incomp + fluid(ktyp52(nr)=3), principle extensions  (ktyp55(nr)=2)
     '  1,1,1,1,1,1,
C **  Incomp + fluid(ktyp52(nr)=3), fibre strains         (ktyp55(nr)=3)
     '  1,1,16,1,1,1,
C **  Compr  + fluid(ktyp52(nr)=4), princ strain invars   (ktyp55(nr)=1)
     '  1,1,1,1,1,1,
C **  Compr  + fluid(ktyp52(nr)=4), principle extensions  (ktyp55(nr)=2)
     '  1,1,1,1,1,1,
C **  Compr  + fluid(ktyp52(nr)=4), fibre strains         (ktyp55(nr)=3)
     '  1,1,16,1,1,1,
C **  Incomp + inext(ktyp52(nr)=5), princ strain invars   (ktyp55(nr)=1)
     '  1,1,1,1,1,1,
C **  Incomp + inext(ktyp52(nr)=5), principle extensions  (ktyp55(nr)=2)
     '  1,1,1,1,1,1,
C **  Incomp + inext(ktyp52(nr)=5), fibre strains         (ktyp55(nr)=3)
     '  1,1,16,1,1,1/

C             *************Material parameter names*************
      DATA CMAT5/            !CMAT5(nm,ktyp56(nr),ktyp55(nr),ktyp52(nr))
C **  Compressible  (ktyp52(nr)=1), princ strain invars   (ktyp55(nr)=1)
     '  35*' ',
     '  'mu','alpha','nu', ! Alex 21-Nov-02: 3 parameters (ktyp56(nr)=2)
     '  32*' ',
     '  140*' ',
C **  Compressible  (ktyp52(nr)=1), principle extensions  (ktyp55(nr)=2)
     '  210*' ',
C **  Compressible  (ktyp52(nr)=1), fibre strains         (ktyp55(nr)=3)
     '  35*' ',                           !polynomial     (ktyp56(nr)=1)
     '  35*' ',                           !exponential    (ktyp56(nr)=2)
     '  'Youngs modulus (kPa)',           !linear elastic (ktyp56(nr)=3)
     '  34*' ',
     '  35*' ',                           !unused         (ktyp56(nr)=4)
     '  35*' ',                           !USER53         (ktyp56(nr)=5)
     '  35*' ',                           !linear viscous (ktyp56(nr)=6)
C **  Incompressible(ktyp52(nr)=2), princ strain invars   (ktyp55(nr)=1)
     '  210*' ',
C **  Incompressible(ktyp52(nr)=2), principle extensions  (ktyp55(nr)=2)
     '  210*' ',
C **  Incompressible(ktyp52(nr)=2), fibre strains         (ktyp55(nr)=3)
     '  35*' ',                           !polynomial     (ktyp56(nr)=1)
     '  35*' ',                           !exponential    (ktyp56(nr)=2)
     '  'fibre axis coefficient (kPa)',   !pole-zero      (ktyp56(nr)=3)
     '  'fibre axis pole',
     '  'fibre axis curvature',
     '  'sheet axis coefficient (kPa)',
     '  'sheet axis pole',
     '  'sheet axis curvature',
     '  'sheet-normal axis coefficient (kPa)',
     '  'sheet-normal axis pole',
     '  'sheet-normal axis curvature',
     '  'fibre angle std dev (degrees)',
     '  'sheet angle std dev (degrees)',
     '  'sheet-normal std dev 1 (deg.s)',
     '  'sheet-normal std dev 2 (deg.s)',
     '  'initial fibre extension ratio',
     '  'initial sheet extension ratio',
     '  'initial sheet-normal extension ratio',
     '  19*' ',
     '  35*' ',                           !fibre in fluid (ktyp56(nr)=4)
     '  35*' ',                           !USER53         (ktyp56(nr)=5)
     '  'fibre elastic modulus',          !linear viscous (ktyp56(nr)=6)
     '  'fibre-direction viscosity',
     '  'transverse dirn viscosity',
     '  'velocity field',
     '  31*' ',
C **  Incomp + fluid(ktyp52(nr)=3), princ strain invars   (ktyp55(nr)=1)
     '  210*' ',
C **  Incomp + fluid(ktyp52(nr)=3), principle extensions  (ktyp55(nr)=2)
     '  210*' ',
C **  Incomp + fluid(ktyp52(nr)=3), fibre strains         (ktyp55(nr)=3)
     '  210*' ',
C **  Compr  + fluid(ktyp52(nr)=4), princ strain invars   (ktyp55(nr)=1)
     '  35*' ',
     '  35*' ',
     '  35*' ',
     '  35*' ',
     '  35*' ',
     '  35*' ',
C **  Compr  + fluid(ktyp52(nr)=4), principle extensions  (ktyp55(nr)=2)
     '  210*' ',
C **  Compr  + fluid(ktyp52(nr)=4), fibre strains         (ktyp55(nr)=3)
     '  210*' ',
C **  Incomp + inext(ktyp52(nr)=5), princ strain invars   (ktyp55(nr)=1)
     '  210*' ',
C **  Incomp + inext(ktyp52(nr)=5), principle extensions  (ktyp55(nr)=2)
     '  210*' ',
C **  Incomp + inext(ktyp52(nr)=5), fibre strains         (ktyp55(nr)=3)
     '  210*' '/

C             *************Default values of mat.l params*************
      DATA RMAT5/            !RMAT5(nm,ktyp56(nr),ktyp55(nr),ktyp52(nr))
C **  Compressible  (ktyp52(nr)=1), princ strain invars   (ktyp55(nr)=1)
     '  210*0.0d0,
C **  Compressible  (ktyp52(nr)=1), principle extensions  (ktyp55(nr)=2)
     '  210*0.0d0,
C **  Compressible  (ktyp52(nr)=1), fibre strains         (ktyp55(nr)=3)
     '  35*0.0d0, !                        polynomial     (ktyp56(nr)=1)
     '  35*0.0d0, !                        exponential    (ktyp56(nr)=2)
     '  10.0d0,   !Young's modulus         linear elastic (ktyp56(nr)=3)
     '  34*0.0d0,
     '  35*0.0d0, !                                       (ktyp56(nr)=4)
     '  35*0.0d0, !                                       (ktyp56(nr)=5)
     '  35*0.0d0, !                        linear viscous (ktyp56(nr)=6)
C **  Incompressible(ktyp52(nr)=2), princ strain invars   (ktyp55(nr)=1)
     '  210*0.0d0,
C **  Incompressible(ktyp52(nr)=2), principle extensions  (ktyp55(nr)=2)
     '  210*0.0d0,
C **  Incompressible(ktyp52(nr)=2), fibre strains         (ktyp55(nr)=3)
     '  35*0.0d0, !                           polynomial  (ktyp56(nr)=1)
     '  35*0.0d0, !                           exponential (ktyp56(nr)=2)
C CS 24 MAY 2001 set defaults to values from HMT papar
C     '  1.937d0,  !fibre axis coefficient     pole-zero   (ktyp56(nr)=3)
     '  0.2d0,  !fibre axis coefficient     pole-zero   (ktyp56(nr)=3)
C     '  0.523d0,  !fibre axis pole
     '  0.22d0,  !fibre axis pole
C     '  1.351d0,  !fibre axis curvature
     '  1.0d0,  !fibre axis curvature
     '  0.028d0,  !sheet axis coefficient
     '  0.681d0,  !sheet axis pole
     '  5.991d0,  !sheet axis curvature
     '  0.310d0,  !sheet-normal axis coefficient
     '  1.037d0,  !sheet-normal axis pole
     '  0.398d0,  !sheet-normal axis curvature
     '  5.0d0,    !fibre angle std dev
     '  10.0d0,   !sheet angle std dev
     '  20.0d0,   !sheet-normal std dev 1
     '  20.0d0,   !sheet-normal std dev 2
     '  1.0d0,    !initial fibre extension ratio
     '  1.0d0,    !initial sheet extension ratio
     '  1.0d0,    !initial sheet-normal extension ratio
     '  19*0.0d0,
     '  35*0.0d0, !                        fibre in fluid (ktyp56(nr)=4)
     '  35*0.0d0, !                        USER53         (ktyp56(nr)=5)
     '  35*0.0d0, !                        linear viscous (ktyp56(nr)=6)
C **  Incomp + fluid(ktyp52(nr)=3), princ strain invars   (ktyp55(nr)=1)
     '  210*0.0d0,
C **  Incomp + fluid(ktyp52(nr)=3), principle extensions  (ktyp55(nr)=2)
     '  210*0.0d0,
C **  Incomp + fluid(ktyp52(nr)=3), fibre strains         (ktyp55(nr)=3)
     '  210*0.0d0,
C **  Compr  + fluid(ktyp52(nr)=4), princ strain invars   (ktyp55(nr)=1)
     '  210*0.0d0,
C **  Compr  + fluid(ktyp52(nr)=4), principle extensions  (ktyp55(nr)=2)
     '  210*0.0d0,
C **  Compr  + fluid(ktyp52(nr)=4), fibre strains         (ktyp55(nr)=3)
     '  210*0.0d0,
C **  Incomp + inext(ktyp52(nr)=5), princ strain invars   (ktyp55(nr)=1)
     '  210*0.0d0,
C **  Incomp + inext(ktyp52(nr)=5), principle extensions  (ktyp55(nr)=2)
     '  210*0.0d0,
C **  Incomp + inext(ktyp52(nr)=5), fibre strains         (ktyp55(nr)=3)
     '  210*0.0d0/

      DATA CMAT6/'fibre axis coefficient (kPa)', !pole-zero  inc shear
     '  'fibre axis pole',
     '  'fibre axis curvature',
     '  'sheet axis coefficient (kPa)',
     '  'sheet axis pole',
     '  'sheet axis curvature',
     '  'sheet-normal axis coefficient (kPa)',
     '  'sheet-normal axis pole',
     '  'sheet-normal axis curvature',
     '  'fibre - sheet shear coefficient (kPa)',
     '  'fibre - sheet shear pole',
     '  'fibre - sheet shear curvature',
     '  'fibre - sheet-normal shear coefficient (kPa)',
     '  'fibre - sheet-normal shear pole',
     '  'fibre - sheet-normal shear curvature',
     '  'sheet - sheet-normal shear coefficient (kPa)',
     '  'sheet - sheet-normal shear pole',
     '  'sheet - sheet-normal shear curvature',
     '  'sheet - fibre shear coefficient (kPa)',
     '  'sheet - fibre shear pole',
     '  'sheet - fibre shear curvature',
     '  'sheet-normal - fibre shear coefficient (kPa)',
     '  'sheet-normal - fibre shear pole',
     '  'sheet-normal - fibre shear curvature',
     '  'sheet-normal - sheet shear coefficient (kPa)',
     '  'sheet-normal - sheet shear pole',
     '  'sheet-normal - sheet shear curvature',
     '  'initial fibre extension ratio',
     '  'initial sheet extension ratio',
     '  'initial sheet-normal extension ratio'/

       DATA RMAT6/0.2d0,   ! fibre axis coefficient (kPa)
     '           0.22d0,   ! fibre axis pole
     '           1.0d0,    ! fibre axis curvature
     '           0.028d0,  ! sheet axis coefficient (kPa)
     '           0.681d0,  ! sheet axis pole
     '           5.991d0,  ! sheet axis curvature
     '           0.310d0,  ! sheet-normal axis coefficient (kPa)
     '           1.037d0,  ! sheet-normal axis pole'
     '           0.398d0,  ! sheet-normal axis curvature
     '           1.0d0,    ! fibre - sheet shear coefficient (kPa)
     '           0.366d0,  ! fibre - sheet shear pole
     '           2.0d0,    ! fibre - sheet shear curvature
     '           1.0d0,    ! fibre - sheet-normal shear coefficient (kPa)
     '           0.366d0,  ! fibre - sheet-normal shear pole
     '           2.0d0,    ! fibre - sheet-normal shear curvature
     '           1.0d0,    ! sheet - sheet-normal shear coefficient (kPa)
     '           0.886d0,  ! sheet - sheet-normal shear pole
     '           2.0d0,    ! sheet - sheet-normal shear curvature
     '           1.0d0,    ! sheet - fibre shear coefficient (kPa)
     '           0.886d0,  ! sheet - fibre shear pole
     '           2.0d0,    ! sheet - fibre shear curvature
     '           1.0d0,    ! sheet-normal - fibre shear coefficient (kPa)
     '           1.183d0,  ! sheet-normal - fibre shear pole
     '           2.0d0,    ! sheet-normal - fibre shear curvature
     '           1.0d0,    ! sheet-normal - sheet shear coefficient (kPa)
     '           1.183d0,  ! sheet-normal - sheet shear pole
     '           2.0d0,    ! sheet-normal - sheet shear curvature
     '           1.0d0,    ! initial fibre extension ratio
     '           1.0d0,    ! initial sheet extension ratio
     '           1.0d0/    ! initial sheet-normal extension ratio


        DATA NHT50/2,2,3,4,3,2,2, ! compr/incomp etc for plane stress
     &             2,2,3,4,3,2,2, ! compr/incomp etc for plane strain
     &             3,4,4,4,5,3,4, ! compr/incomp etc for 3D
     &             3,3,4,4,4,3,3, ! compr/incomp etc for membrane
     &             3,3,4,4,4,3,3, ! compr/incomp etc for string
     &             3,3,4,4,4,3,3/ ! compr/incomp etc for shell

        END


