      SUBROUTINE CELLML_GET_NUM_FUNCTIONS(variant,NUM_FCNS)

C#### Subroutine: CELLML_GET_NUM_FUNCTIONS
C###  Description:
C###  Uses the math object that should have already been created for the
C###  given variant and sets the number of function variables in the
C###  model. Since this routine should never be called unless a model's
C###  math object has been successfully created the error checks have
C###  been left out.

      IMPLICIT NONE

      INCLUDE 'cellml.cmn'

      !Parameter list
      INTEGER variant,NUM_FCNS

      CALL CellMLProcessorGetNumberFunctions(CELLML_MODEL_MATH(variant),
     '  NUM_FCNS)

      RETURN
      END

