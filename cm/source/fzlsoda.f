      subroutine lsoda(func,neq,y,t,tout,itol,rtol,atol,itask,
     1  istate,
     1  iopt,rwork,lrw,iwork,liw,jdum,jt,
     2  aii,aio,control,model,sizes,variant,
     3  ari,aro,derived,parameters,protocol,err)

      implicit none
      external func
      integer neq,itol,itask,istate,iopt,lrw,iwork(*),liw,jdum,jt,
     &  aii(*),aio(*),control(*),model(*),sizes(*),variant,err
      real*8 y(*),t,tout,rwork(*),ari(*),aro(*),derived(*),
     &  parameters(*),protocol(*),rtol,atol

      call FLAG_ERROR(0,'Link with lsoda library'//
     '  ' (need subroutine LSODA)')
      err=1
      end
