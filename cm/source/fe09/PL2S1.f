      REAL*8 FUNCTION PL2S1(I,J,K,XI)

C#### Function: PL2S1
C###  Type: REAL*8
C###  Description:
C###    PL2S1 evaluates special 1D quadratic basis functions at XI
C###    for use in special hermite simplex elements when apex is at
C###    node 1.
C**** I is node position index
C**** J is nodal derivative index
C**** K is partial derivative index
C**** XI is Xi-coordinate

      IMPLICIT NONE
!     Parameter List
      INTEGER I,J,K
      REAL*8 XI

c cpb 28/6/95 Optimising multiplication operations

      GO TO (100,200,300),K
 100    GO TO (110,120),I  !no node derivative at apex node (i.e. J=1)
C 110      PL2S1=1.0D0-2.0D0*XI+XI*XI
 110      PL2S1=(XI-2.0d0)*XI+1.0d0      ! xi^2-2xi+1
          RETURN
 120      GO TO (121,122),J
 121        PL2S1=XI*(2.0d0-XI)          ! -xi^2+2xi
            RETURN
 122        PL2S1=XI*(XI-1.0d0)          ! xi^2-xi
            RETURN
 200    GO TO (210,220),I
C 210      PL2S1=2.0D0*XI-2.0D0
 210      PL2S1=XI+XI-2.0d0              ! 2xi-2
          RETURN
 220      GO TO (221,222),J
C 221        PL2S1=2.0D0-2.0D0*XI
 221        PL2S1=2.0d0-XI-XI            ! -2xi+2
            RETURN
C 222        PL2S1=2.0D0*XI-1.0D0
 222        PL2S1=XI+XI-1.0d0            ! 2xi-1
            RETURN
 300    GO TO (310,320),I
 310      PL2S1=2.0d0                    ! 2
          RETURN
 320      GO TO (321,322),J
 321        PL2S1=-2.0d0                 ! -2
            RETURN
 322        PL2S1=2.0d0                  ! 2
            RETURN
      END


