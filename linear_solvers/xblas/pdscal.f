
      SUBROUTINE PDSCAL(N,ALPHA,X,INCX)

c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)

      IMPLICIT none
!     Parameter List
      INTEGER N,INCX
      REAL*8 ALPHA,X(*)
!     Local Variables
      INTEGER i,nx,nlow,nhigh

      if (n.le.0 .or. incx.eq.0) then
        return

      else if (n.eq.1) then
        x(1) = alpha*x(1)

c     Step lengths = 1
      else if (incx.eq.1) then
c$omp   parallel do private(i)
        do i = 1,n
          x(i) = alpha*x(i)
        enddo
c$omp   end parallel do

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

c$omp   parallel do private(i)
        do i = nlow,nhigh,incx
          x(i) = alpha*x(i)
        enddo
c$omp   end parallel do

      endif

      RETURN
      END
