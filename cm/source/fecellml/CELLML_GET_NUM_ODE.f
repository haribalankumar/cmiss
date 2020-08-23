      SUBROUTINE CELLML_GET_NUM_ODE(variant,NUM_ODE)

C#### Subroutine: CELLML_GET_NUM_ODE
C###  Description:
C###  Uses the math object that should have already been created for the
C###  given variant and sets the number of differential equations in the
C###  model. Since this routine should never be called unless a model's
C###  math object has been successfully created the error checks have
C###  been left out.

      IMPLICIT NONE

      INCLUDE 'cellml.cmn'

      !Parameter list
      INTEGER variant,NUM_ODE

      CALL CellMLProcessorGetNumberODEs(CELLML_MODEL_MATH(variant),
     '  NUM_ODE)

      RETURN
      END

