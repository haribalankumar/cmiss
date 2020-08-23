      SUBROUTINE CELLML_INITIALISE(ERROR_CODE)

C#### Subroutine: CELLML_INITIALISE
C###  Description:
C###  CELLML_INITIALISE must be called before any CellML processing can
C###  occur and CELLML_TERMINATE must not be called until all CellML
C###  processing has been completed. ERROR_CODE is returned with a
C###  non-zero value if an error occurs during initialisation.

      IMPLICIT NONE

      !Parameter list
      INTEGER ERROR_CODE

      CALL CellMLProcessorInitialise(ERROR_CODE)

      RETURN
      END


