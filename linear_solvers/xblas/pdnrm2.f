
      REAL*8 FUNCTION PDNRM2(N,X,INCX)

*  DNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DNRM2 := sqrt( x'*x )

      IMPLICIT none
!     Parameter List
      INTEGER N,INCX
      REAL*8 X(*)
!     Local Variables
      INTEGER i,nx,nlow,nhigh
      REAL*8 sum

      if (n.le.0 .or. incx.eq.0) then
        pdnrm2 = 0.0d0

      else if (n.eq.1) then
        pdnrm2 = abs(x(1))

c     Step lengths = 1
      else if (incx.eq.1) then
        sum = 0.0d0
c$omp   parallel do private(i), reduction(+:sum)
        do i = 1,n
          sum = sum + x(i)**2
        enddo
c$omp   end parallel do
        pdnrm2 = sqrt(sum)

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
          sum = sum + x(i)**2
        enddo
c$omp   end parallel do
        pdnrm2 = sqrt(sum)

      endif

      RETURN
      END
