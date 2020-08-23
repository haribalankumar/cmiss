
c-------------------------------------------------------------------------------

      INTEGER FUNCTION Iamin( n,x,incx )

      IMPLICIT none
      INTEGER n,incx
      INTEGER x(*)

c  Return the value iamin
c
c     iamin = min( abs( x( i ) ) ).

      INTEGER i,nx,nlow,nhigh
      INTEGER xmin

      if (n.lt.1 .or. incx.eq.0) then
        iamin = 0
        return
      else if (n.eq.1) then
        iamin = abs(x(1))
        return
      endif

c     Step lengths = 1
      if (incx.eq.1) then
        xmin = abs(x(1))

        do i = 1,n
          xmin = min(xmin,abs(x(i)))
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

        xmin = abs(x(nlow))
        do i = nlow,nhigh,incx
          xmin = min(xmin,abs(x(i)))
        enddo

      endif
      iamin = xmin

      return
      END
