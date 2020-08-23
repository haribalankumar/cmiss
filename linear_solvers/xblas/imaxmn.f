
c-------------------------------------------------------------------------------

      SUBROUTINE Imaxmn( n,x,incx,xmax,xmin )

      IMPLICIT none
      INTEGER n,incx
      INTEGER x(*),xmax,xmin

c  Return the values xmax and xmin
c
c     xmax = max( x( i ) )
c     xmin = min( x( i ) ).

      INTEGER i,nx,nlow,nhigh

      if (n.lt.1 .or. incx.eq.0) then
        xmax = 0
        xmin = 0
        return
      else if (n.eq.1) then
        xmax = x(1)
        xmin = x(1)
        return
      endif

c     Step lengths = 1
      if (incx.eq.1) then
        xmax = x(1)
        xmin = x(1)

        do i = 1,n
          xmax = max(xmax,x(i))
          xmin = min(xmin,x(i))
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

        xmax = x(nlow)
        xmin = x(nlow)
        do i = nlow,nhigh,incx
          xmax = max(xmax,x(i))
          xmin = min(xmin,x(i))
        enddo

      endif

      return
      END
