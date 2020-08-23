      SUBROUTINE CELLML_GET_CONST_VAR_INIT_VALUE(variant,
     '  index,INITIAL_VALUE)

C#### Subroutine: CELLML_GET_CONST_VAR_INIT_VALUE
C###  Description:
C###  Uses the math object that should have already been created for
C###  this variant and sets the index'th constant variable's initial
C###  value, if one has been provided in the CellML, otherwise a value
C###  of zero is used.

      IMPLICIT NONE

      INCLUDE 'cellml.cmn'

      !Parameter list
      INTEGER variant,index,var
      REAL*8 INITIAL_VALUE

      var = index-1 ! C array starts at 0
      CALL CellMLProcessorGetConstantVariableInitialValue(
     '  CELLML_MODEL_MATH(variant),var,INITIAL_VALUE)

      RETURN
      END


