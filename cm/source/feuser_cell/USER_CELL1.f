      SUBROUTINE USER_CELL1(TIME,Y,DY,CONTROL,MODEL,SIZES,VARIANT,
     & DERIVED,PARAM,PROTOCOL,AII,AIO,ARI,ARO,ERR_CODE)

C#### Subroutine: USER_CELL1
C###  Description:
C###    Calculates the RHS of the user defined cell 1 odes.

      IMPLICIT NONE
      INTEGER SIZES(11)
      INTEGER CONTROL(*),MODEL(*),VARIANT,AII(*),AIO(*),ERR_CODE 
      REAL*8 TIME(*),Y(*),DY(*),DERIVED(*),PARAM(*),PROTOCOL(*),ARI(*),
     & ARO(*)

      ERR_CODE=1
      
      RETURN
      END

