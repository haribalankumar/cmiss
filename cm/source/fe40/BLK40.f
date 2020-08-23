      BLOCK DATA BLK40

C#### BlockData: BLK40
C###  Description:
C###    BLK40 sets up parameter titles in common block

C***  Units: Forces are in kN
C***         Elastic constants are in GPa (or dimensionless)
C***         Density is entered in Kg/m^3 and later converted to Mg/m^3
C***         Ang. freq. is in (/ms)
C***         Displacements are in cm or um
C***         Equations are in kN, stiffness matrix is in GN/m

      IMPLICIT NONE
      INCLUDE 'b40.cmn'

      !           ie= 1 2 3 4 5 6 7 8 9 10 11 12
      DATA NMPP/1,1,2,0,2,2,2,2,2, 0, 2, 2, 12*5, 12*5, 12*9, 12*3/
      DATA NGP /1,2,3,4,1,1,1,1,0, 0, 0, 0/
      DATA NLP /0,2,4,0,1,1,3,0,0, 0, 2, 2/

      INCLUDE 'titl40.cmn'

      DATA TITL42 /'Truss or Cable elements      ',  ! (1)
     '             'Batten elements              ',  ! (2)
     '             'Beam elements (Kirchhoff)    ',  ! (3)
     '             'Link elements                ',  ! (4)
     '             'Membrane elements            ',  ! (5)
     '             'Plate elements (Kirchhoff)   ',  ! (6)
     '             'Shell elements               ',  ! (7)
     '             'Shell/fluid interface        ',  ! (8)
     '             '3D elasticity elements       ',  ! (9)
     '             'Shell tank bottom            ',  !(10)
     '             'Plane stress elements        ',  !(11)
     '             'Plane strain elements        '/  !(12)

C *** Dependent variables
C ***********************
      DATA TITL43/
C ** (1)Truss or Cable elements
     '  'Axial   displ.(u)','Transv. displ.(v)',
     '  'Transv. displ.(w)',3*' ',
C ** (2)Batten elements
     '  'Axial   displ.(u)','Transv. displ.(v)',
     '  'Transv. displ.(w)',3*' ',
C ** (3)Beam elements (Kirchhoff)
     '  'Axial   displ.(u)','Transv. displ.(v)',
     '  'Transv. displ.(w)',3*' ',
C ** (4)Link elements
     '  'X-axis displ. (u)','Y-axis displ. (v)',
     '  'Z-axis displ. (w)',3*' ',
C ** (5)Membrane elements
     '  'Transverse displ.','1st rotation     ',
     '  '2nd rotation     ',3*' ',
C ** (6)Plate elements (Kirchhoff)
     '  'In-plane displ(u)','In-plane displ(v)',
     '  'Transv. displ.(w)',3*' ',
C ** (7)Shell elements
     '  'In-plane displ(u)','In-plane displ(v)',
     '  'Transv. displ.(w)',3*' ',
C ** (8)Shell/Fluid Interface
     '  'In-plane displ(u)','In-plane displ(v)',
     '  'Transv. displ.(w)',3*' ',
C ** (9)3D elasticity elements
     '  '1st displ. cpt.  ','2nd displ. cpt.  ',
     '  '3rd displ. cpt.  ',3*' ',
C ** (10)Shell Tank Bottom
     '  'Transverse displ.',5*' ',
C ** (11)Plane stress elements
     '  '1st displ. cpt.  ','2nd displ. cpt.  ',4*' ',
C ** (12)Plane strain elements
     '  '1st displ. cpt.  ','2nd displ. cpt.  ',4*' '/

C *** Material anisotropy (IMT(ie))
C ***********************
      DATA TITL44 /'Isotropic',
     '             'Transversely isotropic wrt Xi_1',
     '             'Transversely isotropic wrt Xi_2',
     '             'Orthotropic',
     '             'Anisotropic (special case)'/

C *** Beam cross-section types
C ****************************
      DATA TITL45 /
     '  'Rectangular solid','Rectangular tube',
     '  'Ellipsoidal solid','Ellipsoidal tube'/

C *** Material elastic constants (NMPP(ie,IMT(ie)))
C ******************************
      DATA TITL46 /
     '  'Youngs mod.   (GPa)','Poissons ratio     ',7*' ', !isotropic

     '  'Youngs mod.1  (GPa)','Youngs mod.2 (GPa) ', !transv. isotropic 1
     '  'Shear mod.1   (GPa)','Poissons ratio 1   ',
     '  'Poissons ratio 2   ',4*' ',

     '  'Youngs mod.1  (GPa)','Youngs mod.2 (GPa) ', !transv. isotropic 2
     '  'Shear mod.1   (GPa)','Poissons ratio 1   ',
     '  'Poissons ratio 2   ',4*' ',

     '  'Youngs mod.1  (GPa)','Youngs mod.2 (GPa) ', !orthotropic
     '  'Shear mod.1-2 (GPa)','Poissons ratio 1-2 ',5*' ',

     '  'c11','c12','c44',6*' '/ !anisotropic (special case)

C *** Material thermal constants (NTPP(ie,IMT(ie)))
C ******************************
      DATA TITL41 /
     '  'Coeff. of thermal expansion','Temperature (above ambient)',
     '  2*' ',
     '  'Coeff. of thermal expansion','Temperature (above ambient)',
     '  2*' ',
     '  'Coeff. of thermal expansion','Temperature (above ambient)',
     '  2*' ',
     '  'Coeff. of thermal expansion','Temperature (above ambient)',
     '  2*' ',
     '  'Coeff. of thermal expansion','Temperature (above ambient)',
     '  2*' '/

C *** Geometric parameters NGP(ie)
C ************************
      DATA TITL47 /
C ** (1)Truss or Cable elements
     '  'Xsect.n area (cm^2)',5*' ',
C ** (2)Batten elements
     '  'Xsect.n area (cm^2)','2nd mom. Izz (cm^4)',4*' ',
C ** (3)Beam elements (Kirchhoff)
     '  'Depth (cm)','Width (cm)','Wall thickness (cm)',3*' ',
C ** (4)Link elements
     '  'X-dir. stiff.(kN/m)','Y-dir. stiff.(kN/m)',
     '  'Z-dir. stiff.(kN/m)',
     '  'X-Y stiff.  (kNm^2)',2*' ',
C ** (5)Membrane elements
     '  'Mem. thickness (mm)',5*' ',
C ** (6)Plate elements (Kirchhoff)
     '  'Plate thickness    ',5*' ',
C ** (7)Shell elements
     '  'Shell thickness    ',5*' ',
C ** (8)Shell/Fluid Interface
     '  'Shell thickness    ',5*' ',
C ** (9)3D elasticity
     '  6*' ',
C ** (10)Shell Tank Bottom
     '  6*' ',
C ** (11)Plane stress elements
     '  6*' ',
C ** (12)Plane strain elements
     '  6*' '/

C *** Element loads NLP(ie)
C *****************
      DATA TITL48 /
C ** (1)Truss or Cable elements
     '  6*' ',
C ** (2)Batten elements
     '  'Axial load     (kN)','Transverse load (y)',4*' ',
C ** (3)Beam elements (Kirchhoff)
     '  'Axial load (+ve t.)','Axial torque       ',
     '  'Transverse load (y)',
     '  'XY-plane moment    ','Transverse load (z)',
     '  'XZ-plane moment    ',
C ** (4)Link elements
     '  6*' ',
C ** (5)Membrane elements
     '  'Pressure load (kPa)','1st cpt of moment  ',
     '  '2nd cpt of moment  ',
     '  3*' ',
C ** (6)Plate elements (Kirchhoff)
     '  'Transverse loading ',5*' ',
C ** (7)Shell elements
     '  '1st in-plane load  ','2nd in-plane load  ',
     '  'Transverse loading ',3*' ',
C ** (8)Shell/Fluid Interface
     '  6*' ',
C ** (9)3D elasticity
     '  6*' ',
C ** (10)Shell Tank Bottom
     '  6*' ',
C ** (11)Plane stress elements
     '  '1st load component ','2nd load component ',4*' ',
C ** (12)Plane strain elements
     '  '1st load component ','2nd load component ',4*' '/

C *** Direction cosines
C*********************
      DATA TITL49 /
     '  'dir.n cosines of normal to 1st princ. bending axis',
     '  'dir.n cosines of normal to 2nd princ. bending axis'/
      END


