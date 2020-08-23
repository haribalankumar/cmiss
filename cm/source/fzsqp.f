      SUBROUTINE MinLSSQP( obj,x1,fv,fj,ldfj,m,n,iter,ieval,hess,ldh,
     &                     bl,bu,act,func,iuser,user,optlst,ierr )

      IMPLICIT none
      INTEGER m,n,ldfj,ldh,iter,ieval,optlst(*),ierr
      INTEGER act(n),iuser(*)
      REAL*8 obj,x1(n),fv(m),fj(ldfj,n),bl(n),bu(n),user(*)
      REAL*8 hess(ldh,n)
      EXTERNAL func

      call FLAG_ERROR(0,'Link with sqp library'//
     '  ' (need subroutine MinLSSQP)')
      obj = 0d0/0d0
      iter = 0
      ieval = 0
      ierr = -1

      return
      END

      SUBROUTINE MinLSSQPOption()
      call FLAG_ERROR(0,'Link with sqp library'//
     '  ' (need subroutine MinLSSQPOption)')
      return
      END

      SUBROUTINE MinLSSQPOptionR()
      call FLAG_ERROR(0,'Link with sqp library'//
     '  ' (need subroutine MinLSSQPOptionR)')
      return
      END

      SUBROUTINE MinLSSQPOptionI()
      call FLAG_ERROR(0,'Link with sqp library'//
     '  ' (need subroutine MinLSSQPOptionI)')
      return
      END

      SUBROUTINE MinLSSQPOptionL()
      call FLAG_ERROR(0,'Link with sqp library'//
     '  ' (need subroutine MinLSSQPOptionL)')
      return
      END

