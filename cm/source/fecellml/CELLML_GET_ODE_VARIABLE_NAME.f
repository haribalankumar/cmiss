      SUBROUTINE CELLML_GET_ODE_VARIABLE_NAME(variant,index,NAME)

C#### Subroutine: CELLML_GET_ODE_VARIABLE_NAME
C###  Description:
C###  Uses the math object that should have already been created for
C###  this variant and sets the index'th ODE variable's name.

      IMPLICIT NONE

      INCLUDE 'cellml.cmn'

      !Parameter list
      INTEGER variant,index,var
      CHARACTER NAME*(*)

      var = index-1 ! C array starts at 0
      CALL CellMLProcessorGetDifferentialVariableName(
     '  CELLML_MODEL_MATH(variant),var,NAME)

      RETURN
      END


