      SUBROUTINE EVALUATE_MATHS_SELF(MATHS_LABEL,XG,DXIX,ERROR,*)

C#### Subroutine: EVALUATE_MATHS_SELF
C###  Description:
C###    An indirect subroutine for EVALUATE_MATHS to call itself.
C###    It appears to be a Fortran77 quirk.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'maths00.cmn'
      
!     Parameter List
      REAL*8 DXIX(3,3),XG(NJM,NUM)
      CHARACTER MATHS_LABEL*(*),ERROR*(*)
      
      CALL ENTERS('EVALUATE_MATHS_SELF',*9999)
      
      CALL EVALUATE_MATHS(MATHS_LABEL,XG,DXIX,ERROR,*9999)
      
      CALL EXITS('EVALUATE_MATHS_SELF')
      RETURN
 9999 CALL ERRORS('EVALUATE_MATHS_SELF',ERROR)
      CALL EXITS('EVALUATE_MATHS_SELF')
      RETURN 1
      END

