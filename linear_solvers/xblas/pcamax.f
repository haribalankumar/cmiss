
c-------------------------------------------------------------------------------

      REAL FUNCTION Pcamax( n,x,incx )

      IMPLICIT none
      INTEGER n,incx
      COMPLEX x(*)

c  Return the value camax
c
c     camax = max( abs( x( i ) ) )

      INTEGER i,nx,nlow,nhigh
      REAL xmax

      if (n.lt.1 .or. incx.eq.0) then
        pcamax = 0.0
        return
      else if (n.eq.1) then
        pcamax = cabs(x(1))
        return
      endif

c     Step lengths = 1
      if (incx.eq.1) then
        xmax = cabs(x(1))

c$omp   parallel do private(i), reduction(max:xmax)
        do i = 1,n
          xmax = max(xmax,cabs(x(i)))
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

        xmax = cabs(x(nlow))
c$omp   parallel do private(i), reduction(max:xmax)
        do i = nlow,nhigh,incx
          xmax = max(xmax,cabs(x(i)))
        enddo
c$omp   end parallel do

      endif
      pcamax = xmax

      return
      END
