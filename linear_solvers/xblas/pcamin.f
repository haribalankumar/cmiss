
c-------------------------------------------------------------------------------

      REAL FUNCTION Pcamin( n,x,incx )

      IMPLICIT none
      INTEGER n,incx
      COMPLEX x(*)

c  Return the value camin
c
c     camin = min( abs( x( i ) ) ).

      INTEGER i,nx,nlow,nhigh
      REAL xmin

      if (n.lt.1 .or. incx.eq.0) then
        pcamin = 0.0
        return
      else if (n.eq.1) then
        pcamin = cabs(x(1))
        return
      endif

c     Step lengths = 1
      if (incx.eq.1) then
        xmin = cabs(x(1))

c$omp   parallel do private(i), reduction(min:xmin)
        do i = 1,n
          xmin = min(xmin,cabs(x(i)))
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

        xmin = cabs(x(nlow))
c$omp   parallel do private(i), reduction(min:xmin)
        do i = nlow,nhigh,incx
          xmin = min(xmin,cabs(x(i)))
        enddo
c$omp   end parallel do

      endif
      pcamin = xmin

      return
      END
