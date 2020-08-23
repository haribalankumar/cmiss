      REAL*8 FUNCTION PHT3(I,J,K,XI)

C#### Function: PHT3
C###  Type: REAL*8
C###  Description:
C###    PHT3 evaluates 1D cubic Hermite basis function at XI.

C**** I is node position index
C**** J is nodal derivative index
C**** K is partial derivative index
C**** XI is Xi-coordinate

      IMPLICIT NONE
!     Parameter List
      INTEGER i,j,k
      REAL*8 XI

      GO TO (100,200,300,400),K
 100    GO TO (110,120),I
 110      GO TO (111,112),J
 111        PHT3=1.D0+XI*XI*(2.D0*XI-3.D0)
            RETURN
 112        PHT3=XI*(XI-1.D0)*(XI-1.D0)
            RETURN
 120      GO TO (121,122),J
 121        PHT3=XI*XI*(3.D0-2.D0*XI)
            RETURN
 122        PHT3=XI*XI*(XI-1.D0)
            RETURN
 200      GO TO (210,220),I
 210        GO TO (211,212),J
 211          PHT3=6.D0*XI*(XI-1.D0)
              RETURN
 212          PHT3=(XI-1.D0)*(3.D0*XI-1.D0)
              RETURN
 220        GO TO (221,222),J
 221          PHT3=6.D0*XI*(1.D0-XI)
              RETURN
 222          PHT3=XI*(3.D0*XI-2.D0)
              RETURN
 300        GO TO (310,320),I
 310          GO TO (311,312),J
 311            PHT3=12.D0*XI-6.D0
                RETURN
 312            PHT3=6.D0*XI-4.D0
                RETURN
 320          GO TO (321,322),J
 321            PHT3=6.D0-12.D0*XI
                RETURN
 322            PHT3=6.D0*XI-2.D0
                RETURN
 400          GO TO (410,420),I
 410            GO TO (411,412),J
 411              PHT3=12.D0
                  RETURN
 412              PHT3=6.D0
                  RETURN
 420            GO TO (421,422),J
 421              PHT3=-12.D0
                  RETURN
 422              PHT3=6.D0
                  RETURN
      END


