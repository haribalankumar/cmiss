      REAL*8 FUNCTION PH3(I,J,K,XI)

C#### Function: PH3
C###  Type: REAL*8
C###  Description:
C###    PH3(I,J,K,XI) evaluates 1D cubic Hermite basis function at XI
C###    where I is node position index, J is nodal derivative index
C###    and K is partial derivative index.

      IMPLICIT NONE
!     Parameter List
      INTEGER I,J,K
      REAL*8 XI

c cpb 28/6/95 Optimising mulitplication operations

      GO TO (100,200,300),K
 100    GO TO (110,120),I
 110      GO TO (111,112),J
 111        PH3=(2.0d0*XI-3.0d0)*XI*XI+1.0d0  ! 2xi^3-3xi^2+1
            RETURN
 112        PH3=((XI-2.0d0)*XI+1.0d0)*XI      ! xi^3-2xi^2+xi
            RETURN
 120      GO TO (121,122),J
 121        PH3=XI*XI*(3.0d0-2.0d0*XI)        ! -2xi^3+3xi^2
            RETURN
 122        PH3=XI*XI*(XI-1.0d0)              ! xi^3-xi^2
            RETURN
 200    GO TO (210,220),I
 210      GO TO (211,212),J
 211        PH3=6.0d0*XI*(XI-1.0d0)           ! 6xi^2-6xi
            RETURN
 212        PH3=(3.0d0*XI-4.0d0)*XI+1.0d0     ! 3xi^2-4xi+1
            RETURN
 220      GO TO (221,222),J
 221        PH3=6.0d0*XI*(1.0d0-XI)           ! -6xi^2+6xi
            RETURN
 222        PH3=XI*(3.0d0*XI-2.0d0)           ! 3xi^2-2xi
            RETURN
 300    GO TO (310,320),I
 310      GO TO (311,312),J
 311        PH3=12.0d0*XI-6.0d0               ! 12xi-6
            RETURN
 312        PH3=6.0d0*XI-4.0d0                ! 6xi-4
            RETURN
 320      GO TO (321,322),J
 321        PH3=6.0d0-12.0d0*XI               ! -12xi+6
            RETURN
 322        PH3=6.0d0*XI-2.0d0                ! 6xi-2
            RETURN
      END


