      SUBROUTINE XHERM(Q11,Q12,Q21,Q22,X1,X2,XGLOB,YGLOB,XX,ERROR,*)

C#### Subroutine: XHERM
C###  Description:
C###    XHERM returns cubic Hermite polyline XX(n,j),n=1,40 in global
C###    coordinates. Curve is defined over two line segments in local
C###    coordinates by X1(nk,j,nn) and X2(nk,j,nn),where line segments
C###    transform to global system with Q11,Q12,Q21,Q22 and XGLOB,YGLOB.

      IMPLICIT NONE
!     Parameter List
      REAL*8 Q11,Q12,Q21,Q22,XGLOB,XX(40,*),X1(2,2,*),X2(2,2,*),YGLOB
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER j,n
      REAL*8 X(2),XH3,XI

      CALL ENTERS('XHERM',*9999)
      DO n=1,21
        XI=DBLE(n-1)/20.D0
        DO j=1,2
          X(j)=XH3(j,X1,XI)
        ENDDO
        XX(n,1)=XGLOB+X(1)*Q11+X(2)*Q12
        XX(n,2)=YGLOB+X(1)*Q21+X(2)*Q22
      ENDDO
      DO n=21,40
        XI=DBLE(N-20)/20.D0
        DO j=1,2
          X(j)=XH3(j,X2,XI)
        ENDDO
        XX(n,1)=XGLOB+X(1)*Q11+X(2)*Q12
        XX(n,2)=YGLOB+X(1)*Q21+X(2)*Q22
      ENDDO

      CALL EXITS('XHERM')
      RETURN
 9999 CALL ERRORS('XHERM',ERROR)
      CALL EXITS('XHERM')
      RETURN 1
      END
