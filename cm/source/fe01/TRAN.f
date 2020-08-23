      SUBROUTINE TRAN(ip,A,B,C)

C#### Subroutine: TRAN
C###  Description:
C###    TRAN transforms 2nd order 3*3 tensor A into B with coordinate
C###    transformation matrix C,
C###              covariantly (C(t)AC) if IP=1,
C###           or contravariantly (CAC(t)) if IP=2,
C###           or mixed variantly (C(-1)AC -similarity trans.) if IP=3.

      IMPLICIT NONE
!     Parameter List
      INTEGER ip
      REAL*8  A(3,3),B(3,3),C(3,3)
!     Local Variables
      INTEGER i,j
      REAL*8 CC,D(3,3)

C     CALL ENTERS('TRAN',*9999)
      IF(ip.EQ.1) THEN
        DO i=1,3
          DO j=1,3
            B(i,j)=(A(1,1)*C(1,i)+A(2,1)*C(2,i)+A(3,1)*C(3,i))*C(1,j)
     '            +(A(1,2)*C(1,i)+A(2,2)*C(2,i)+A(3,2)*C(3,i))*C(2,j)
     '            +(A(1,3)*C(1,i)+A(2,3)*C(2,i)+A(3,3)*C(3,i))*C(3,j)
          ENDDO !j
        ENDDO !i

      ELSE IF(IP.EQ.2) THEN
        DO i=1,3
          DO j=1,3
            B(i,j)=(A(1,1)*C(i,1)+A(2,1)*C(i,2)+A(3,1)*C(i,3))*C(j,1)
     '            +(A(1,2)*C(i,1)+A(2,2)*C(i,2)+A(3,2)*C(i,3))*C(j,2)
     '            +(A(1,3)*C(i,1)+A(2,3)*C(i,2)+A(3,3)*C(i,3))*C(j,3)
          ENDDO !j
        ENDDO !i

      ELSE IF(IP.EQ.3) THEN
        CALL INVERT(3,C,D,CC)
        DO i=1,3
          DO j=1,3
            B(i,j)=(A(1,1)*D(i,1)+A(2,1)*D(i,2)+A(3,1)*D(i,3))*C(1,j)
     '            +(A(1,2)*D(i,1)+A(2,2)*D(i,2)+A(3,2)*D(i,3))*C(2,j)
     '            +(A(1,3)*D(i,1)+A(2,3)*D(i,2)+A(3,3)*D(i,3))*C(3,j)
          ENDDO !j
        ENDDO !i
      ENDIF

C     CALL EXITS('TRAN')
      RETURN
      END


