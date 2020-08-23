
c-------------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION Damax( n,x,incx )

      IMPLICIT none
      INTEGER n,incx
      DOUBLE PRECISION x(*)

c  Return the value damax
c
c     damax = max( abs( x( i ) ) )

      INTEGER i,nx,nlow,nhigh
      DOUBLE PRECISION xmax

      if (n.lt.1 .or. incx.eq.0) then
        damax = 0.0d0
        return
      else if (n.eq.1) then
        damax = abs(x(1))
        return
      endif

c     Step lengths = 1
      if (incx.eq.1) then
        xmax = abs(x(1))

        do i = 1,n
          xmax = max(xmax,abs(x(i)))
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

        xmax = abs(x(nlow))
        do i = nlow,nhigh,incx
          xmax = max(xmax,abs(x(i)))
        enddo

      endif
      damax = xmax

      return
      END
