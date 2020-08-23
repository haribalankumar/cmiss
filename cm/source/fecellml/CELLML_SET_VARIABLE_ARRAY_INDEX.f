      SUBROUTINE CELLML_SET_VARIABLE_ARRAY_INDEX(variant,NAME,INDEX)

C#### Subroutine: CELLML_SET_VARIABLE_ARRAY_INDEX
C###  Description:
C###  Uses the math object that should have already been created for
C###  this variant and sets the variable with the given name to be
C###  stored at the given index in its array.

      IMPLICIT NONE

      INCLUDE 'cellml.cmn'

      !Parameter list
      INTEGER variant,INDEX
      CHARACTER NAME*(*)
      !Local variables
      INTEGER IBEG,IEND

      CALL STRING_TRIM(NAME,IBEG,IEND)

      CALL CellMLProcessorSetVariableArrayIndex(
     '  CELLML_MODEL_MATH(variant),NAME(IBEG:IEND),INDEX)

      RETURN
      END

