
c-------------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION Pzamin( n,x,incx )

      IMPLICIT none
      INTEGER n,incx
      DOUBLE COMPLEX x(*)

c  Return the value zamin
c
c     zamin = min( abs( x( i ) ) ).

      INTEGER i,nx,nlow,nhigh
      DOUBLE PRECISION xmin

      if (n.lt.1 .or. incx.eq.0) then
        pzamin = 0.0d0
        return
      else if (n.eq.1) then
        pzamin = cdabs(x(1))
        return
      endif

c     Step lengths = 1
      if (incx.eq.1) then
        xmin = cdabs(x(1))

c$omp   parallel do private(i), reduction(min:xmin)
        do i = 1,n
          xmin = min(xmin,cdabs(x(i)))
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

        xmin = cdabs(x(nlow))
c$omp   parallel do private(i), reduction(min:xmin)
        do i = nlow,nhigh,incx
          xmin = min(xmin,cdabs(x(i)))
        enddo
c$omp   end parallel do

      endif
      pzamin = xmin

      return
      END
