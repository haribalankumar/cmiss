
c-------------------------------------------------------------------------------

      INTEGER FUNCTION Iiamax( n,x,incx )

      IMPLICIT none
      INTEGER n,incx
      INTEGER x(*)

c  Find the index of element having max. absolute value.

      INTEGER xmax,nx,nlow,nhigh
      INTEGER i,imax

      if (n.lt.1 .or. incx.eq.0) then
        iiamax = 0
        return
      else if (n.eq.1) then
        iiamax = 1
        return
      endif

c     Step lengths = 1
      if (incx.eq.1) then
        xmax = x(1)
        imax = 1
        do i = 1,n
          if (abs(x(i)).gt.xmax) then
            imax = i
            xmax = abs(xmax)
          endif
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
        imax = nlow
        do i = nlow,nhigh,incx
          if (abs(x(i)).gt.xmax) then
            imax = i
            xmax = abs(xmax)
          endif
        enddo

      endif
      iiamax = imax

      return
      END
