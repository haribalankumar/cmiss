
      REAL*8 FUNCTION PDDOT(N,X,INCX,Y,INCY)

c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)

      IMPLICIT none
!     Parameter List
      INTEGER N,INCX,INCY
      REAL*8 X(*),Y(*)
!     Local Variables
      INTEGER i,j,k,nx,nlow,nhigh,ilow,jlow
      REAL*8 sum

      if (n.le.0 .or. incx.eq.0 .or. incy.eq.0) then
        pddot = 0.0d0
        return
      else if (n.eq.1) then
        pddot = x(1)*y(1)
        return
      endif

c     Step lengths = 1
      if (incx.eq.1 .and. incy.eq.1) then
        sum = 0.0d0
c$omp   parallel do private(i), reduction(+:sum)
        do i = 1,n
          sum = sum + x(i)*y(i)
        enddo
c$omp   end parallel do
        pddot = sum

c     One step length = 1
      else if (incx.eq.1) then
        if (incy.gt.0) then
          jlow = 1
        else
          jlow = 1+( n-1 )*incy
        endif

        sum = 0.0d0
c$omp   parallel do private(i,j), reduction(+:sum)
        do i = 1,n
          j = jlow + (i-1)*incy
          sum = sum + x(i)*y(j)
        enddo
c$omp   end parallel do
        pddot = sum

      else if (incy.eq.1) then
        if (incx.gt.0) then
          ilow = 1
        else
          ilow = 1+( n-1 )*incx
        endif

        sum = 0.0d0
c$omp   parallel do private(i,j), reduction(+:sum)
        do j = 1,n
          i = ilow + (j-1)*incx
          sum = sum + x(i)*y(j)
        enddo
c$omp   end parallel do
        pddot = sum

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

        sum = 0.0d0
c$omp   parallel do private(i), reduction(+:sum)
        do i = nlow,nhigh,incx
          sum = sum + x(i)*y(i)
        enddo
c$omp   end parallel do
        pddot = sum

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

        sum = 0.0d0
c$omp   parallel do private(i,j,k), reduction(+:sum)
        do k = 0,n-1
          i = ilow + k*incx
          j = jlow + k*incy 
          sum = sum + x(i)*y(j)
        enddo
c$omp   end parallel do
        pddot = sum

      endif

      RETURN
      END
