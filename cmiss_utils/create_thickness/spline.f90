! -*-f90-*-
MODULE spline

  USE constants
  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------

  SUBROUTINE InterpLinear( x_in,x_out,n_in,n_out,m )

    ! On a linear spline constructed from the points x_in, 
    ! define a set of equally spaced points x_out.

    IMPLICIT NONE
    ! Subroutine arguments
    INTEGER(I4)::m,n_in,n_out
    REAL(DP)::x_in(n_in,m),x_out(n_out,m)
    ! Local variables
    INTEGER(I4)::i,j
    REAL(DP)::dx,x(n_in),xlen,x_pt

    ! Evaluate the distances of the points along the spline
    x(1)=0.0D0
    DO i=2,n_in
      dx=0.0D0
      DO j=1,m
        dx=dx+(x_in(i,j)-x_in(i-1,j))**2
      ENDDO
      dx=SQRT(dx)

      x(i)=x(i-1)+dx
    ENDDO
    xlen=x(n_in)

    ! Interpolate
    DO i=1,n_out
      IF(n_out.lt.2) THEN
        x_pt=0.5D0*xlen
      ELSE
        x_pt=xlen*DBLE(i-1)/DBLE(n_out-1)
      ENDIF
      CALL InterpPointLinear( x,x_in,n_in,m,x_pt,x_out(i,:) )
    ENDDO

    RETURN
  END SUBROUTINE InterpLinear

  !----------------------------------------------------------------

  SUBROUTINE InterpPointLinear( xa,ya,n,m,x,y )

    ! Evaluate a point on the linear spline
    IMPLICIT NONE
    ! Function arguments
    INTEGER(I4)::n,m
    REAL(DP)::x,y(m),xa(n),ya(n,m)
    ! Local variables
    INTEGER(I4)::k,khi,klo
    REAL(DP)::a,b,h

    klo=1
    khi=n
    DO WHILE (khi-klo.gt.1)
      k=(khi+klo)/2
      if(xa(k).gt.x)then
        khi=k
      else
        klo=k
      endif
    ENDDO

    h=xa(khi)-xa(klo)
    IF (h.eq.0.0D0) THEN
      WRITE(stderr,*) 'InterpPointLinear: error in input values',xa(khi),xa(klo)
      STOP
    ENDIF

    a=( xa(khi)-x )/h
    b=( x-xa(klo) )/h
    y(:)=a*ya(klo,:) + b*ya(khi,:)

    RETURN
  END SUBROUTINE InterpPointLinear

  !----------------------------------------------------------------

  SUBROUTINE InterpCubic( x_in,x_out,n_in,n_out,m )

    ! On a cubic spline constructed from the points x_in, 
    ! define a set of equally spaced points x_out.

    IMPLICIT NONE
    ! Subroutine arguments
    INTEGER(I4)::m,n_in,n_out
    REAL(DP)::x_in(n_in,m),x_out(n_out,m)
    ! Local variables
    INTEGER(I4)::i,j
    REAL(DP)::dx,x(n_in),y2(n_in,m),xlen,x_pt

    ! Evaluate the distances of the points along the spline
    x(1)=0.0D0
    DO i=1,n_in
      dx=0.0D0
      DO j=1,m
        dx=dx+(x_in(i,j)-x_in(i-1,j))**2
      ENDDO
      dx=SQRT(dx)

      x(i)=x(i-1)+dx
    ENDDO
    xlen=x(n_in)

    ! Evaluate the cubic spline coefficents
    CALL InterpCubicSpline( x,x_in,n_in,m,y2 )

    ! Interpolate
    DO i=1,n_out
      x_pt=xlen*DBLE(i-1)/DBLE(n_out-1)
      ! CALL InterpPointLinear( x,x_in,n_in,m,x_pt,x_out(i,:) )
      CALL InterpPointCubic( x,x_in,y2,n_in,m,x_pt,x_out(i,:) )
    ENDDO

    RETURN
  END SUBROUTINE InterpCubic

  !----------------------------------------------------------------

  SUBROUTINE InterpCubicSpline( x,y,n,m,y2 )

    ! Construct the equations for a cubic spline

    IMPLICIT none
    ! Function arguments
    INTEGER(I4)::n,m
    REAL(DP)::x(n),y(n,m),y2(n,m)
    ! Local parameters
    REAL(DP),PARAMETER::BIG_VAL=1.0D30
    ! Local variables
    INTEGER(I4)::i,j,k
    REAL(DP)::p,qn,sig,un,u(n)


    ! Evaluate for each axis
    DO j=1,m
      y2(1,j)=0.0
      u(1)=0.0

      DO i=2,n-1
        sig=( x(i)-x(i-1) )/( x(i+1)-x(i-1) )
        p=sig*y2(i-1,j) + 2.0
        y2(i,j)=( sig-1.0 )/p
        u(i)=( 6.0*( ( y(i+1,j)-y(i,j) )/( x(i+1)-x(i) )   &
          &        - ( y(i,j)-y(i-1,j) )/( x(i)-x(i-1) ) ) &
          &        / ( x(i+1)-x(i-1) ) &
          &  - sig*u(i-1) )/p
      ENDDO
      y2(n,j)=0.0

      DO k=n-1,1,-1
        y2(k,j)=y2(k,j)*y2(k+1,j) + u(k)
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE InterpCubicSpline

  !----------------------------------------------------------------

  SUBROUTINE InterpPointCubic( xa,ya,y2a,n,m,x,y )

    ! Interpolate to find a point in a cubic spline

    IMPLICIT NONE
    ! Function arguments
    INTEGER(I4)::n,m
    REAL(DP)::x,y(m),xa(n),y2a(n,m),ya(n,m)
    ! Local variables
    INTEGER(I4)::k,khi,klo
    REAL(DP)::a,b,h

    klo=1
    khi=n
    DO WHILE (khi-klo.gt.1)
      k=(khi+klo)/2
      if(xa(k).gt.x)then
        khi=k
      else
        klo=k
      endif
    ENDDO

    h=xa(khi)-xa(klo)
    IF (h.eq.0.0D0) THEN
      WRITE(stderr,*) 'InterpPointCubic: error in input values',xa(khi),xa(klo)
      STOP
    ENDIF

    a=( xa(khi)-x )/h
    b=( x-xa(klo) )/h
    y(:)=a*ya(klo,:) + b*ya(khi,:) &
      & + ( (a**3-a)*y2a(klo,:) + (b**3-b)*y2a(khi,:) )*(h**2)/6.0

    RETURN
  END SUBROUTINE InterpPointCubic

  !----------------------------------------------------------------

END MODULE spline
