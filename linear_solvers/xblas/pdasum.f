
      REAL*8 FUNCTION PDASUM(N,X,INCX)

c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c     converted to omp 6/6/02

      IMPLICIT none
!     Parameter List
      INTEGER N,INCX
      REAL*8 X(*)
!     Local Variables
      INTEGER i,nx,nlow,nhigh
      REAL*8 sum

      if (n.le.0 .or. incx.eq.0) then
        pdasum = 0.0d0

      else if (n.eq.1) then
        pdasum = abs(x(1))

c     Step lengths = 1
      else if (incx.eq.1) then
        sum = 0.0d0
c$omp   parallel do private(i), reduction(+:sum)
        do i = 1,n
          sum = sum + abs(x(i))
        enddo
c$omp   end parallel do
        pdasum = sum

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

        sum = 0.0d0
c$omp   parallel do private(i), reduction(+:sum)
        do i = nlow,nhigh,incx
          sum = sum + abs(x(i))
        enddo
c$omp   end parallel do
        pdasum = sum

      endif

      RETURN
      END
