
! petg00.inc
C     Constants used in determining how much of the Petrov derivative
C     weight to use with the Galerkin weight.  They should depend on the
C     interpolation type but here they are only set up for cubic Hermites.
      REAL*8 PETRATIO,PETLIMIT
      PARAMETER (PETRATIO=0.02d0,PETLIMIT=2.0d0)
C     Constant used to smooth Petrov weights near singularities.
      REAL*8 PETSMOOTH
      PARAMETER (PETSMOOTH=0.25d0) !must be the in interval (0,1)
