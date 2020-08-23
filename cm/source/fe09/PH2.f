      REAL*8 FUNCTION PH2(I,J,K,N,XI)

C#### Function: PH2
C###  Type: REAL*8
C###  Description:
C###    PH2 evaluates 1D quadratic Hermite basis function at XI.

C**** I is node position index
C**** J is nodal derivative index
C**** K is partial derivative index
C**** N is the local node number with no derivative term.
C**** XI is Xi-coordinate

      IMPLICIT NONE
!     Parameter List
      INTEGER I,J,K,N
      REAL*8 XI

      GOTO (1000,2000) N
 1000   GOTO (1100,1200,1300),K
 1100     GOTO (1110,1120),I
 1110       GOTO (1111,1112),J
 1111         PH2=(XI-2.0d0)*XI+1.0d0       ! xi^2-2xi+1
              RETURN
 1112         PH2=0.0d0                     ! 0
              RETURN
 1120       GOTO (1121,1122),J
 1121         PH2=(2.0d0-XI)*XI             ! -xi^2+2xi
              RETURN
 1122         PH2=(XI-1.0d0)*XI             ! xi^2-xi
              RETURN
 1200      GOTO (1210,1220),I
 1210        GOTO (1211,1212),J
 1211          PH2=2.0d0*XI-2.0d0           ! 2xi-2
               RETURN
 1212          PH2=0.0d0                    ! 0
               RETURN
 1220        GOTO (1221,1222),J
 1221          PH2=-2.0d0*XI+2.0d0          ! -2xi+2
               RETURN
 1222          PH2=2.0d0*XI-1.0d0           ! 2xi-1
               RETURN
 1300      GOTO (1310,1320),I
 1310        GOTO (1311,1312),J
 1311          PH2=2.0d0                    ! 2
               RETURN
 1312          PH2=0.0d0                    ! 0
               RETURN
 1320        GOTO (1321,1322),J
 1321          PH2=-2.0d0                   ! -2
               RETURN
 1322          PH2=2.0d0                    ! 2
               RETURN
 2000   GOTO (2100,2200,2300),K
 2100     GOTO (2110,2120),I
 2110       GOTO (2111,2112),J
 2111         PH2=1.0d0-XI*XI               ! -xi^2+1
              RETURN
 2112         PH2=XI*(1.0d0-XI)             ! -xi^2+xi
              RETURN
 2120       GOTO (2121,2122),J
 2121         PH2=XI*XI                     ! xi^2
              RETURN
 2122         PH2=0.0d0                     ! 0
              RETURN
 2200      GOTO (2210,2220),I
 2210        GOTO (2211,2212),J
 2211          PH2=-2.0d0*XI                ! -2xi
               RETURN
 2212          PH2=1.0d0-2.0d0*XI           ! -2xi+1
               RETURN
 2220        GOTO (2221,2222),J
 2221          PH2=2.0d0*XI                 ! 2xi
               RETURN
 2222          PH2=0.0d0                    ! 0
               RETURN
 2300      GOTO (2310,2320),I
 2310        GOTO (2311,2312),J
 2311          PH2=-2.0d0                   ! -2
               RETURN
 2312          PH2=-2.0d0                   ! -2
               RETURN
 2320        GOTO (2321,2322),J
 2321          PH2=2.0d0                    ! 2
               RETURN
 2322          PH2=0.0d0                    ! 0
               RETURN
      END


