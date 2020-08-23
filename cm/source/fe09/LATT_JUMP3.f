      INTEGER FUNCTION LATT_JUMP3(i,i1,i2,j1,j2)

C#### Function: LATT_JUMP3 
C###  Type: INTEGER
C###  Description:
C###    LATT_JUMP3 return a value from a 4 part piecewise
C###    interpolation between jmax and jmin.      
C###    The smallest absolute gradient segment occurs
C###    immediately before or immediately after the lowest
C###    of the j1,j2 values      
      
      IMPLICIT NONE
!     Parameter List
      INTEGER i,i1,i2,j1,j2

      IF(j1.GE.j2) THEN
        IF(i.GE.i2) THEN
          LATT_JUMP3=j2
        ELSE
          LATT_JUMP3=MIN(j1,MAX(j2,i1-i+j1))
        ENDIF
      ELSE
        IF(i.LE.i1) THEN
          LATT_JUMP3=j1
        ELSE
          LATT_JUMP3=MIN(j2,MAX(j1,i-i2+j2))
        ENDIF
      ENDIF
      
      RETURN
      END
      
