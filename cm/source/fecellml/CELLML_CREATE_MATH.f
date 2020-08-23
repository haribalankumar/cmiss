      SUBROUTINE CELLML_CREATE_MATH(variant,ERROR_CODE)

C#### Subroutine: CELLML_CREATE_MATH
C###  Description:
C###  Takes the model object for the given variant and creates the
C###  corresponding math object which can then be used for querying the
C###  mathematics of a model and creating the RHS routine for use when
C###  solving the model.

      IMPLICIT NONE

      INCLUDE 'cellml.cmn'

      !Parameter list
      INTEGER variant,ERROR_CODE

      !Initialise error code
      ERROR_CODE = 0

      CALL CellMLProcessorCreateMath(CELLML_MODELS(variant),
     '  CELLML_MODEL_MATH(variant),ERROR_CODE)

      RETURN
      END


