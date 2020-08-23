      SUBROUTINE USER_CELL5(T,Y,DY,CONTROL,MODEL,SIZES,VARIANT,
     '  DERIVED,PARAM,PROTOCOL,AII,AIO,ARI,ARO,ERR_CODE)

C#### Subroutine: USER_CELL5
C###  Description:
C###    Calculates the RHS of the user defined cell 5 odes.

      IMPLICIT NONE

      INCLUDE 'cell_reserved.inc'

!     Parameters variables
      INTEGER SIZES(10)
      INTEGER CONTROL(*),MODEL(*),VARIANT,AII(*),
     '  AIO(*),ERR_CODE
      REAL*8 T,Y(*),DY(*),DERIVED(*),
     '  PARAM(*),PROTOCOL(*),ARI(*),
     '  ARO(*)

      ERR_CODE=1

      RETURN
      END
