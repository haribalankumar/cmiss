
c-------------------------------------------------------------------------------

      REAL FUNCTION Psamax( n,x,incx )

      IMPLICIT none
      INTEGER n,incx
      REAL x(*)

c  Return the value samax
c
c     samax = max( abs( x( i ) ) )

      INTEGER i,nx,nlow,nhigh
      REAL xmax

      if (n.lt.1 .or. incx.eq.0) then
        psamax = 0.0
        return
      else if (n.eq.1) then
        psamax = abs(x(1))
        return
      endif

c     Step lengths = 1
      if (incx.eq.1) then
        xmax = abs(x(1))

c$omp   parallel do private(i), reduction(max:xmax)
        do i = 1,n
          xmax = max(xmax,abs(x(i)))
        enddo
c$omp   end parallel do

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
c$omp   parallel do private(i), reduction(max:xmax)
        do i = nlow,nhigh,incx
          xmax = max(xmax,abs(x(i)))
        enddo
c$omp   end parallel do

      endif
      psamax = xmax

      return
      END
