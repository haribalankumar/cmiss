
c-------------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION Pzamax( n,x,incx )

      IMPLICIT none
      INTEGER n,incx
      DOUBLE COMPLEX x(*)

c  Return the value zamax
c
c     zamax = max( abs( x( i ) ) )

      INTEGER i,nx,nlow,nhigh
      DOUBLE PRECISION xmax

      if (n.lt.1 .or. incx.eq.0) then
        pzamax = 0.0d0
        return
      else if (n.eq.1) then
        pzamax = cdabs(x(1))
        return
      endif

c     Step lengths = 1
      if (incx.eq.1) then
        xmax = cdabs(x(1))

c$omp   parallel do private(i), reduction(max:xmax)
        do i = 1,n
          xmax = max(xmax,cdabs(x(i)))
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

        xmax = cdabs(x(nlow))
c$omp   parallel do private(i), reduction(max:xmax)
        do i = nlow,nhigh,incx
          xmax = max(xmax,cdabs(x(i)))
        enddo
c$omp   end parallel do

      endif
      pzamax = xmax

      return
      END
