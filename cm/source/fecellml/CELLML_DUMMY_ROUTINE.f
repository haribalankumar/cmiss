      SUBROUTINE CELLML_DUMMY_ROUTINE()

C#### Subroutine: CELLML_DUMMY_ROUTINE
C###  Description:
C###  A dummy routine to be passed into MARCH8 as the RHSROUTINE
C###  variable because for CellML we want to switch routine based on the
C###  variant of the grid points which is not done until INTEGRATOR is
C###  called.

      IMPLICIT NONE

      RETURN
      END

