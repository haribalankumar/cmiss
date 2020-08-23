C**** CMISS Module fzxblas.f: Dummy routines for the xblas library

      SUBROUTINE PDCOPY(N,X,INCX,Y,INCY)

      IMPLICIT none
      INTEGER N,INCX,INCY
      REAL*8 X(*),Y(*)

      CALL FLAG_ERROR(0,'Link with Xblas library: need sub PDCOPY')

      if (n.le.0 .or. incx.eq.0 .or. incy.eq.0) return      
      y(1) = 0d0/0d0

      RETURN
      END

      SUBROUTINE PDLOAD(n,a,x,incx)

      IMPLICIT none
      INTEGER n,incx
      DOUBLE PRECISION a,x(*)

      CALL FLAG_ERROR(0,'Link with Xblas library: need sub PDLOAD')

      if (n.le.0 .or. incx.eq.0) return
      x(1) = 0d0/0d0

      RETURN
      END

