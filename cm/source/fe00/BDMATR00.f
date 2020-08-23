      BLOCK DATA BDMATR00

      IMPLICIT NONE
      INCLUDE 'matr00.cmn'

C***  MATRTYPE indicates whether the matrix is a matrix (=1) or a
C***  vector (=0).
      DATA MATRTYPE / 1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,0 /

C***  MATRNAME gives the name of the matrix.
      DATA MATRNAME /    'GK         ','GQ         ','GD         ',
     '     'GM         ','GR         ','GKK        ','GMM        ',
     '     'GRR        ','YP         ','T_BH       ','T_BH_INV   ',
     '     'PHI        ','ZCROSSING  ','PHI_H      ','YQ         ',
     '     'YQS        ','PHI_H_EXACT','MFI        ','LD_NP      '/
      END


