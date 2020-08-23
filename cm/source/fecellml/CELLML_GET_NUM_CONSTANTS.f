      SUBROUTINE CELLML_GET_NUM_CONSTANTS(variant,NUM_FCNS)

C#### Subroutine: CELLML_GET_NUM_CONSTANTS
C###  Description:
C###  Uses the math object that should have already been created for the
C###  given variant and sets the number of constant variables in the
C###  model. Since this routine should never be called unless a model's
C###  math object has been successfully created the error checks have
C###  been left out.

      IMPLICIT NONE

      INCLUDE 'cellml.cmn'

      !Parameter list
      INTEGER variant,NUM_FCNS

      CALL CellMLProcessorGetNumberConstants(CELLML_MODEL_MATH(variant),
     '  NUM_FCNS)

      RETURN
      END

