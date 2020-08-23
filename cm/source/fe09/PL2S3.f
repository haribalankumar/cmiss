      REAL*8 FUNCTION PL2S3(I,J,K,XI)

C#### Function: PL2S3
C###  Type: REAL*8
C###  Description:
C###    PL2S3 evaluates special 1D quadratic basis functions at XI
C###    for use in special hermite simplex elements when apex is at
C###    node 3.
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
 100    GO TO (110,120),I
 110      GO TO (111,112),J
 111        PL2S3=1.0d0-XI*XI               ! -xi^2-1
            RETURN
 112        PL2S3=XI*(1.0d0-XI)             ! -xi^2+xi
            RETURN
 120      PL2S3=XI*XI !no node derivative at apex node (i.e. J=1)
          RETURN
 200    GO TO (210,220),I
 210      GO TO (211,212),J
 211        PL2S3=-2.0d0*XI                 ! -2xi
            RETURN
C 212        PL2S3=1.0d0-2.0d0*XI
 212        PL2S3=1.0d0-XI-XI               ! -2xi+1
            RETURN
 220      PL2S3=2.0d0*XI                    ! 2xi
          RETURN
 300    GO TO (310,320),I
 310      GO TO (311,312),J
 311        PL2S3=-2.0d0                    ! -2
            RETURN
 312        PL2S3=-2.0d0                    ! -2
            RETURN
 320      PL2S3=2.0d0                       ! 2
          RETURN
      END


