
c-------------------------------------------------------------------------------

      REAL FUNCTION Camin( n,x,incx )

      IMPLICIT none
      INTEGER n,incx
      COMPLEX x(*)

c  Return the value camin
c
c     camin = min( abs( x( i ) ) ).

      INTEGER i,nx,nlow,nhigh
      REAL xmin

      if (n.lt.1 .or. incx.eq.0) then
        camin = 0.0
        return
      else if (n.eq.1) then
        camin = cabs(x(1))
        return
      endif

c     Step lengths = 1
      if (incx.eq.1) then
        xmin = cabs(x(1))

        do i = 1,n
          xmin = min(xmin,cabs(x(i)))
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

        xmin = cabs(x(nlow))
        do i = nlow,nhigh,incx
          xmin = min(xmin,cabs(x(i)))
        enddo

      endif
      camin = xmin

      return
      END
