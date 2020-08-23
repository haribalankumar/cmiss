
c-------------------------------------------------------------------------------

      INTEGER FUNCTION Iasum( n,x,incx )

      IMPLICIT none
      INTEGER n,incx
      INTEGER x(*)

c  Takes the sum of the absolute values.

      INTEGER i,sum,nx,nlow,nhigh

      if (n.le.0 .or. incx.eq.0) then
        iasum = 0
        return
      endif

c     Step lengths = 1
      if (incx.eq.1) then
        sum = 0
        do i = 1,n
          sum = sum + x(i)
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

        sum = 0
        do i = nlow,nhigh,incx
          sum = sum + x(i)
        enddo

      endif
      iasum = sum

      return
      END
