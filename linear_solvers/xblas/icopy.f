
c-------------------------------------------------------------------------------

      SUBROUTINE Icopy( n,x,incx,y,incy )

      IMPLICIT none
      INTEGER n,incx,incy
      INTEGER x(*),y(*)

c  Copies a vector, x, to a vector, y.

      INTEGER i,j,k,nx,nlow,nhigh,ilow,jlow

      if (n.le.0 .or. incx.eq.0 .or. incy.eq.0) then
        return
      else if (n.eq.1) then
        y(1) = x(1)
        return
      endif

c     Step lengths = 1
      if (incx.eq.1 .and. incy.eq.1) then
        do i = 1,n
          y(i) = x(i)
        enddo

c     One step length = 1
      else if (incx.eq.1) then
        if (incy.gt.0) then
          jlow = 1
        else
          jlow = 1+( n-1 )*incy
        endif

        do i = 1,n
          j = jlow + (i-1)*incy
          y(j) = x(i)
        enddo

      else if (incy.eq.1) then
        if (incx.gt.0) then
          ilow = 1
        else
          ilow = 1+( n-1 )*incx
        endif

        do j = 1,n
          i = ilow + (j-1)*incx
          y(j) = x(i)
        enddo

c     Equal steps
      else if (incx.eq.incy) then
        nx = 1+( n-1 )*incx

        if (incx.gt.0) then
          nlow  = 1
          nhigh = nx
        else
          nlow  = nx
          nhigh = 1
        endif

        do i = nlow,nhigh,incx
          y(i) = x(i)
        enddo

c     Non-equal steps
      else
        if (incx.gt.0) then
          ilow = 1
        else
          ilow = 1+( n-1 )*incx
        endif

        if (incy.gt.0) then
          jlow = 1
        else
          jlow = 1+( n-1 )*incy
        endif

        do k = 0,n-1
          i = ilow + k*incx
          j = jlow + k*incy 
          y(j) = x(i)
        enddo

      endif

      return
      END
