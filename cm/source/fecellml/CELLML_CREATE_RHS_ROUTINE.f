      SUBROUTINE CELLML_CREATE_RHS_ROUTINE(variant,DEBUG,SAVE_FILES,
     '  ERROR_CODE)

C#### Subroutine: CELLML_CREATE_RHS_ROUTINE
C###  Description:
C###  Uses the math object that should have already been created for
C###  this variant, and had all the variable array's defined for it, to
C###  create a RHS routine that can be called for the integration of the
C###  model. Will leave generated temporary files unless SAVE_FILES is
C###  set to false. Will generate debug code if DEBUG is true.

      IMPLICIT NONE

      INCLUDE 'cellml.cmn'

      !Parameter list
      INTEGER variant,ERROR_CODE
      LOGICAL DEBUG,SAVE_FILES
      !local variables
      INTEGER ISTIM,Vm,IBEG,IEND,SAVEF,DEB

      !don't know about passing logicals, so convert to integers
      IF (CELLML_CONTAINS_VM(variant)) THEN
        Vm = 1
      ELSE
        Vm = 0
      ENDIF
      IF (CELLML_CONTAINS_ISTIM(variant)) THEN
        ISTIM = 1
        CALL STRING_TRIM(CELLML_ISTIM_ARRAY_NAME(variant),IBEG,IEND)
      ELSE
        ISTIM = 0
        IBEG = 1
        IEND=1
      ENDIF
      IF (SAVE_FILES) THEN
        SAVEF = 1
      ELSE
        SAVEF = 0
      ENDIF
      IF (DEBUG) THEN
        DEB = 1
      ELSE
        DEB = 0
      ENDIF

      CALL CellMLProcessorCreateRHSRoutine(CELLML_MODEL_MATH(variant),
     '  CELLML_DSOS(variant),CELLML_ROUTINES(variant),Vm,ISTIM,
     '  CELLML_ISTIM_ARRAY_NAME(variant)(IBEG:IEND),
     '  CELLML_ISTIM_ARRAY_INDEX(variant),SAVEF,DEB,ERROR_CODE)

      !IF (ERROR_CODE.EQ.0) THEN
      !  CALL DUMMY_TEST(%val(CELLML_ROUTINES(variant)),1)
      !  CALL DUMMY_TEST(%val(CELLML_ROUTINES(variant)),2)
      !  CALL DUMMY_TEST(%val(CELLML_ROUTINES(variant)),3)
      !ENDIF

      RETURN
      END

