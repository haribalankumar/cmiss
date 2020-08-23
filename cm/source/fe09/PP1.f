      REAL*8 FUNCTION PP1(NUM_TERMS,COEFFS,X,Y,Z)

C#### Function: PP1
C###  Type: REAL*8
C###  Description:
C###    PP1 evaluates a polynomial at coordinates (x,y,z).
C###
C****  The terms used in the polynomial expansion are
C****  NUM_TERMS=3 : P=[1,x,y]
C****  NUM_TERMS=4 : P=[1,x,y,xy]
C****  NUM_TERMS=6 : P=[1,x,y,xy,x^2,y^2]
C****  NUM_TERMS=8 : P=[1,x,y,xy,xz,yz,xyz]
C****  NUM_TERMS=9 : P=[1,x,y,xy,x^2,y^2,x^2y,xy^2,x^2y^2]

C**** Created by Carey Stevens 9 July 1997

      IMPLICIT NONE
!     Parameter List
      INTEGER NUM_TERMS
      REAL*8 COEFFS(*),X,Y,Z

      IF(NUM_TERMS.EQ.3) THEN
        PP1=COEFFS(1)+COEFFS(2)*X+COEFFS(3)*Y
      ELSE IF(NUM_TERMS.EQ.4) THEN
        PP1=COEFFS(1)+COEFFS(2)*X+COEFFS(3)*Y+COEFFS(4)*X*Y
      ELSE IF(NUM_TERMS.EQ.6) THEN
        PP1=COEFFS(1)+COEFFS(2)*X+COEFFS(3)*Y+COEFFS(4)*X*Y+
     '    COEFFS(5)*X*X+COEFFS(6)*Y*Y
      ELSE IF(NUM_TERMS.EQ.8) THEN
        PP1=COEFFS(1)+COEFFS(2)*X+COEFFS(3)*Y+COEFFS(4)*Z+
     '    COEFFS(5)*X*Y+COEFFS(6)*X*Z+COEFFS(7)*Y*Z+
     '    COEFFS(8)*X*Y*Z
      ELSE
        PP1=COEFFS(1)+COEFFS(2)*X+COEFFS(3)*Y+COEFFS(4)*X*Y+
     '    COEFFS(5)*X*X+COEFFS(6)*Y*Y+COEFFS(7)*X*X*Y+
     '    COEFFS(8)*X*Y*Y+COEFFS(9)*X*X*Y*Y
      ENDIF

      RETURN
      END


