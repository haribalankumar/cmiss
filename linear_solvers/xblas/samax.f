
c-------------------------------------------------------------------------------

      REAL FUNCTION Samax( n,x,incx )

      IMPLICIT none
      INTEGER n,incx
      REAL x(*)

c  Return the value samax
c
c     samax = max( abs( x( i ) ) )

      INTEGER i,nx,nlow,nhigh
      REAL xmax

      if (n.lt.1 .or. incx.eq.0) then
        samax = 0.0
        return
      else if (n.eq.1) then
        samax = abs(x(1))
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
      samax = xmax

      return
      END
