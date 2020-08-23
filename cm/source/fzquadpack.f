      SUBROUTINE dqags(f,a,b,epsabs,epsrel,result,abserr,neval,ier,
     *   limit,lenw,last,iwork,work)

      double precision a,abserr,b,epsabs,epsrel,f,result,work
      integer ier,iwork,last,lenw,limit,neval
c
      dimension iwork(limit),work(lenw)
c
      external f

      call FLAG_ERROR(0,'Link with quadpack library'//
     '  ' (need subroutine dqags)')

      return
      END

