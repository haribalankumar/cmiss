
c-------------------------------------------------------------------------------

      SUBROUTINE Lload( n,a,x,incx )

      IMPLICIT none
      INTEGER n,incx
      LOGICAL a,x(*)

c  Set a vector, x, to a scalar, a.

      INTEGER i,nx,nlow,nhigh

      if (n.le.0 .or. incx.eq.0) return

c     Step lengths = 1
      if (incx.eq.1) then
        do i = 1,n
          x(i) = a
        enddo

c     Other step lengths
      else
        nx = 1+( n-1 )*incx

        if (incx.gt.0) then
          nlow  = 1
          nhigh = nx
        else
          nlow  = nx
          nhigh = 1
        endif

        do i = nlow,nhigh,incx
          x(i) = a
        enddo

      endif

      return
      END
