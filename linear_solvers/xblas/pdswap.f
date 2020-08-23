
      SUBROUTINE PDSWAP(N,X,INCX,Y,INCY)

c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)

      IMPLICIT none
!     Parameter List
      INTEGER N,INCX,INCY
      REAL*8 X(*),Y(*)
!     Local Variables
      INTEGER i,j,k,nx,nlow,nhigh,ilow,jlow
      REAL*8 swap

      if (n.le.0 .or. incx.eq.0 .or. incy.eq.0) then
        return
      else if (n.eq.1) then
        swap = y(1)
        y(1) = x(1)
        x(1) = swap
        return
      endif

c     Step lengths = 1
      if (incx.eq.1 .and. incy.eq.1) then
c$omp   parallel do private(i,swap)
        do i = 1,n
          swap = y(i)
          y(i) = x(i)
          x(i) = swap
        enddo
c$omp   end parallel do

c     One step length = 1
      else if (incx.eq.1) then
        if (incy.gt.0) then
          jlow = 1
        else
          jlow = 1+( n-1 )*incy
        endif

c$omp   parallel do private(i,j,swap)
        do i = 1,n
          j = jlow + (i-1)*incy
          swap = y(i)
          y(j) = x(i)
          x(i) = swap
        enddo
c$omp   end parallel do

      else if (incy.eq.1) then
        if (incx.gt.0) then
          ilow = 1
        else
          ilow = 1+( n-1 )*incx
        endif

c$omp   parallel do private(i,j,swap)
        do j = 1,n
          i = ilow + (j-1)*incx
          swap = y(i)
          y(j) = x(i)
          x(i) = swap
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

c$omp   parallel do private(i,swap)
        do i = nlow,nhigh,incx
          swap = y(i)
          y(i) = x(i)
          x(i) = swap
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

c$omp   parallel do private(i,j,k,swap)
        do k = 0,n-1
          i = ilow + k*incx
          j = jlow + k*incy 
          swap = y(i)
          y(j) = x(i)
          x(i) = swap
        enddo
c$omp   end parallel do

      endif

      RETURN
      END
