      SUBROUTINE POLAR(id,F,R,U,ERROR,*)

C#### Subroutine: POLAR
C###  Description:
C###    POLAR calculates right polar decomposition of F for
C###    2D (ID=2) or 3D (ID=3).

C     Note: 3D case still needs to be tested !!    -dr 13/10/91

      IMPLICIT NONE
!     Parameter List
      INTEGER id
      REAL*8 F(3,3),R(3,3),U(3,3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ip,jp,kp
      REAL*8 C(3,3),C2D(2,2),DET,DUM1(3,3),E11,E22,E33,EVAL(3),
     '  EVAL2D(2),P(3,3),P2D(2,2),UI(3,3),XLAM(3,3)
      DATA XLAM/9*0.d0/

      CALL ENTERS('POLAR',*9999)
      IF(id.EQ.2) THEN            !2D case
        C2D(1,1)=F(1,1)*F(1,1)+F(2,1)*F(2,1)
        C2D(2,1)=F(1,2)*F(1,1)+F(2,2)*F(2,1)
        C2D(1,2)=C2D(2,1)
        C2D(2,2)=F(1,2)*F(1,2)+F(2,2)*F(2,2)

        CALL EIGEN1(2,2,C2D,EVAL2D,P2D,ERROR,*9999)
        E11=DSQRT(DABS(EVAL2D(1)))            !2 principal extensions
        E22=DSQRT(DABS(EVAL2D(2)))
        U(1,1)=E11*P2D(1,1)*P2D(1,1)+E22*P2D(1,2)*P2D(1,2)
        U(1,2)=E11*P2D(1,1)*P2D(2,1)+E22*P2D(1,2)*P2D(2,2)
        U(2,1)=U(1,2)
        U(2,2)=E11*P2D(2,1)*P2D(2,1)+E22*P2D(2,2)*P2D(2,2)
        DET=U(1,1)*U(2,2)-U(1,2)*U(2,1)
        UI(1,1)= U(2,2)/DET
        UI(1,2)=-U(1,2)/DET
        UI(2,1)=-U(2,1)/DET
        UI(2,2)= U(1,1)/DET
        R(1,1)=F(1,1)*UI(1,1)+F(1,2)*UI(2,1)
        R(1,2)=F(1,1)*UI(1,2)+F(1,2)*UI(2,2)
        R(2,1)=F(2,1)*UI(1,1)+F(2,2)*UI(2,1)
        R(2,2)=F(2,1)*UI(1,2)+F(2,2)*UI(2,2)

      ELSE IF(id.EQ.3) THEN        !3D case

        DO ip=1,3        !Establish C
          DO jp=1,3
            C(ip,jp)=0.d0
            DO kp=1,3
              C(ip,jp)=C(ip,jp)+F(kp,ip)*F(kp,jp)
            ENDDO
          ENDDO
        ENDDO

        CALL EIGEN1(3,3,C,EVAL,P,ERROR,*9999)
        E11=DSQRT(DABS(EVAL(1)))              !3 principal extensions
        E22=DSQRT(DABS(EVAL(2)))
        E33=DSQRT(DABS(EVAL(3)))

        XLAM(1,1)=E11  !components of U referred to principal axis sys.
        XLAM(2,2)=E22
        XLAM(3,3)=E33

        DO ip=1,3      !Mult. XLAM by P-transpose to get DUM1
          DO jp=1,3
            DUM1(ip,jp)=0.d0
            DO kp=1,3
              DUM1(ip,jp)=DUM1(ip,jp)+XLAM(ip,kp)*P(jp,kp)
            ENDDO
          ENDDO
        ENDDO
        DO ip=1,3    !Mult. P by DUM1 to get U (referred to global axes)
          DO jp=1,3
            U(ip,jp)=0.d0
            DO kp=1,3
              U(ip,jp)=U(ip,jp)+P(ip,kp)*DUM1(kp,jp)
            ENDDO
          ENDDO
        ENDDO
        CALL INVERT(3,U,UI,DET)         !Invert U to gt UI
        DO ip=1,3                       !Postmultiply F by UI to get R
          DO jp=1,3
            R(ip,jp)=0.d0
            DO kp=1,3
              R(ip,jp)=R(ip,jp)+F(ip,kp)*UI(kp,jp)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('POLAR')
      RETURN
 9999 CALL ERRORS('POLAR',ERROR)
      CALL EXITS('POLAR')
      RETURN 1
      END


