
c-------------------------------------------------------------------------------

      SUBROUTINE Plcopy( n,x,incx,y,incy )

      IMPLICIT none
      INTEGER n,incx,incy
      LOGICAL x(*),y(*)

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
c$omp   parallel do private(i)
        do i = 1,n
          y(i) = x(i)
        enddo
c$omp   end parallel do

c     One step length = 1
      else if (incx.eq.1) then
        if (incy.gt.0) then
          jlow = 1
        else
          jlow = 1+( n-1 )*incy
        endif

c$omp   parallel do private(i,j)
        do i = 1,n
          j = jlow + (i-1)*incy
          y(j) = x(i)
        enddo
c$omp   end parallel do

      else if (incy.eq.1) then
        if (incx.gt.0) then
          ilow = 1
        else
          ilow = 1+( n-1 )*incx
        endif

c$omp   parallel do private(i,j)
        do j = 1,n
          i = ilow + (j-1)*incx
          y(j) = x(i)
        enddo
c$omp   end parallel do

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

c$omp   parallel do private(i)
        do i = nlow,nhigh,incx
          y(i) = x(i)
        enddo
c$omp   end parallel do

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

c$omp   parallel do private(i,j,k)
        do k = 0,n-1
          i = ilow + k*incx
          j = jlow + k*incy 
          y(j) = x(i)
        enddo
c$omp   end parallel do

      endif

      return
      END
