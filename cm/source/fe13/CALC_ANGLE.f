      REAL*8 FUNCTION CALC_ANGLE(U,V)

C#### Function: CALC_ANGLE
C###  Description:
C###    ANGLE calculates the angle between two vectors

      IMPLICIT NONE

!     Parameter List
      REAL*8 U(3),V(3)
!     Local variables
      INTEGER nj
      REAL*8 ANGLE,length_u,length_v,N_U(3),N_V(3)
      REAL*8 SCALAR

      length_u=0.d0
      length_v=0.d0
      DO nj=1,3
        N_U(nj)=U(nj)
        N_V(nj)=V(nj)
        length_u=length_u+N_U(nj)**2
        length_v=length_v+N_V(nj)**2
      ENDDO !nj
      length_u=DSQRT(length_u)
      length_v=DSQRT(length_v)
      DO nj=1,3
        N_U(nj)=N_U(nj)/length_u
        N_V(nj)=N_V(nj)/length_v
      ENDDO !nj

      ANGLE=SCALAR(3,N_U,N_V)
      ANGLE=MAX(-1.d0,ANGLE)
      ANGLE=MIN(1.d0,ANGLE)
      ANGLE=DACOS(ANGLE)

      CALC_ANGLE=ANGLE
      
      RETURN
      END


