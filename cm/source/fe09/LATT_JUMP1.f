      INTEGER FUNCTION LATT_JUMP1(i,i1,i2,j1,j2)

C#### Function: LATT_JUMP1
C###  Type: INTEGER
C###  Description:
C###    LATT_JUMP1 returns a value from a 4 part
C###    piecewise interpolation between jmax and jmin.
C###    This value is chosen so as to correctly position
C###    grid points between jmax and jmin.      
      
      IMPLICIT NONE
!     Parameter List
      INTEGER i,i1,i2,j1,j2

      IF(j1.GE.j2) THEN
        IF(i.LT.i2) THEN
          LATT_JUMP1=min(i2-i-1+max(j1-min(i2-i1-1,j1-j2),j2+1),j1)
        ELSE
          LATT_JUMP1=j2
        ENDIF
      ELSE
        IF(i.LE.i1) THEN
          LATT_JUMP1=j1
        ELSEIF(i.LE.max(i1,i2-1)) THEN
          LATT_JUMP1=min(i-i1+j1,j2)
        ELSE
          LATT_JUMP1=j2
        ENDIF
      ENDIF
      RETURN
      END
      
