      SUBROUTINE EVAL_MODEL(ICQS,RCQS,YQS,ERROR,*)

C#### Subroutine: EVAL_MODEL
C###  Description:
C###    EVAL_MODEL evaluates a maths model

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'cell_reserved.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ptr00.cmn'
!     Parameter List
      INTEGER ICQS(NQIM)
      REAL*8 RCQS(NQRM),YQS(NIQSM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER VARIANT

      CALL ENTERS('EVAL_MODEL',*9999)
      
      VARIANT = ICQS(CELL_VARIANT_OFFSET)

C call the RHSROUTINE through EVAL_MODEL_DYNAM
      CALL EVAL_MODEL_DYNAM(ICQS,RCQS,
     &  %VAL(CELLML_ROUTINES(VARIANT)),
     &  VARIANT,YQS,ERROR,*9999)

      CALL EXITS('EVAL_MODEL')
      RETURN
 9999 CALL ERRORS('EVAL_MODEL',ERROR)
      CALL EXITS('EVAL_MODEL')
      RETURN 1
      END

