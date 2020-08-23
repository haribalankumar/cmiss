      SUBROUTINE CELLML_SET_VARIABLE_ARRAY(variant,NAME,ARRAY)

C#### Subroutine: CELLML_SET_VARIABLE_ARRAY
C###  Description:
C###  Uses the math object that should have already been created for
C###  this variant and sets the variable with the given name to be
C###  stored in the given array.

      IMPLICIT NONE

      INCLUDE 'cellml.cmn'

      !Parameter list
      INTEGER variant
      CHARACTER NAME*(*),ARRAY*(*)
      !Local variables
      INTEGER IBEG,IEND

      CALL STRING_TRIM(NAME,IBEG,IEND)

      IF (ARRAY(1:1).EQ.'Y') THEN
        CALL CellMLProcessorSetVariableArrayY(
     '    CELLML_MODEL_MATH(variant),NAME(IBEG:IEND))
      ELSEIF (ARRAY(1:3).EQ.'FCN') THEN
        CALL CellMLProcessorSetVariableArrayFCN(
     '    CELLML_MODEL_MATH(variant),NAME(IBEG:IEND))
      ELSEIF (ARRAY(1:5).EQ.'PARAM') THEN
        CALL CellMLProcessorSetVariableArrayPARAM(
     '    CELLML_MODEL_MATH(variant),NAME(IBEG:IEND))
      ELSEIF (ARRAY(1:7).EQ.'DERIVED') THEN
        CALL CellMLProcessorSetVariableArrayDERIVED(
     '    CELLML_MODEL_MATH(variant),NAME(IBEG:IEND))
      ELSE
        CALL CellMLProcessorSetVariableArrayLOCAL(
     '    CELLML_MODEL_MATH(variant),NAME(IBEG:IEND))
      ENDIF

      RETURN
      END


