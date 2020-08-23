      REAL*8 FUNCTION DMEDIAN(A,B,C)

C#### Function: DMEDIAN
C###  Type: REAL*8
C###  Description: Double precision median

      REAL*8 A,B,C

      IF(A.GT.B) THEN
        IF(C.GT.A) THEN
          DMEDIAN=A
        ELSEIF(B.GT.C) THEN
          DMEDIAN=B
        ELSE
          DMEDIAN=C
        ENDIF
      ELSE
        IF(C.LT.A) THEN
          DMEDIAN=A
        ELSEIF(B.LT.C) THEN
          DMEDIAN=B
        ELSE
          DMEDIAN=C
        ENDIF
      ENDIF
      RETURN
      END

