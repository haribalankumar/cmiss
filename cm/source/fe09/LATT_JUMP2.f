      INTEGER FUNCTION LATT_JUMP2(i,i1,i2,j1,j2)

C#### Function: LATT_JUMP2 
C###  Type: INTEGER
C###  Description:
C###    LATT_JUMP2 return a value from a piecewise
C###    interpolation between j1 and j2.      
C###    This value is chosen so as to correctly position
C###    grid points between j1 and j2.
      
      IMPLICIT NONE
!     Parameter List
      INTEGER i,i1,i2,j1,j2
!     Local Variables
      INTEGER dj,jm,jstar
      
      IF(j1.GE.j2) THEN
        IF(i.LE.i1) THEN
          LATT_JUMP2=j1
        ELSEIF(i.LT.i2) THEN
          jstar=i1-i+j1
          dj=min(i2-i1-1,j1-j2)
          jm=max(j1-dj,j2+1)
          LATT_JUMP2=min(max(jstar,jm),j1)
        ELSE
          LATT_JUMP2=j2
        ENDIF
      ELSE
        IF(i.LE.i1) THEN
          LATT_JUMP2=j1
        ELSEIF(i.LE.max(i1,i2-1)) THEN
          jstar=i-i1+j1
          dj=min(i2-i1-1,j2-j1)
          jm=max(min(j1+dj,j2-1),j1+1)
          LATT_JUMP2=min(jstar,jm)
        ELSE
          LATT_JUMP2=j2
        ENDIF
      ENDIF
      
      RETURN
      END

