      subroutine lsoda_wrapper(func,neq,y,t,tout,itol,rtol,atol,itask,
     1  istate,iopt,rwork,lrw,iwork,liw,jdum,jt,aii,aio,control,model,
     2  sizes,variant,ari,aro,derived,parameters,protocol,err,
     3  rrwork,iiwork,error,*)

C#### Subroutine: LSODA_WRAPPER
C###  Description:
C###  Provides a CMISS wrapper to the LSODA ode solver
C###  Added by Richard Boyes 8th december 2000 upon
C###  request from physiome sciences. The lsoda package
C###  exists in an external library that is compiled into
C###  CMISS externally - modified 2001-10-28 by DPN to use dynamically
C###  allocated workspace and fixing multiprocessing.

C#### Comment: LSODA
C###  Description:
C###    <html>
C###      <p>
C###        The following is a brief description of the LSODA integrator as
C###        used in CMISS. For a more complete description of LSODA and all
C###        its configurable options see the LSODA source code.
C###      </p>
C###      <p>
C###        The version of the LSODA integrator used in CMISS is based on the
C###        1997-02-24 version of LSODA - <i>livermore solver for
C###        ordinary differential equations, with automatic method switching
C###        for stiff and nonstiff problems</i> - from netlib.org. It has been
C###        modified to enable its use in a multiprocessing environment - the
C###        main modification being the movement of common block variables to
C###        external storage space, enabling the existing common block
C###        variables to be separated for each integration point. The function
C###        calls used to evaluate the differentials has also been modified to
C###        match the existing CMISS standard - for example, see <a
C###        href="http://www.esc.auckland.ac.nz/Groups/Bioengineering/CMISS/scripts/cmiss-lookup.cgi?name=noble98_cell">NOBLE98_CELL</a>.
C###      </p>
C###      <p>
C###        For use in CMISS, the interface to the external LSODA library is via
C###        the wrapper routine <a
C###        href="http://www.esc.auckland.ac.nz/Groups/Bioengineering/CMISS/scripts/cmiss-lookup.cgi?name=lsoda_wrapper">LSODA_WRAPPER</a>.
C###        The following is a overview of some of the arguments to this wrapper
C###        routine.
C###        <ul>
C###          <li><b>neq</b> - the number of ODE variables to integrate
C###          </li>
C###          <li><b>y</b> - the array of state variable values. When istate
C###            is 1, y gives the intial values of the state variables, if
C###            istate is not 1 y is not used for input. On return, y contains
C###            the state variable values at time = tout
C###          </li>
C###          <li><b>t,tout</b> - on entry t is the current time and tout is
C###            the time you wish to integrate to. tout is unchanged on exit,
C###            and t is set to the time at which integration finished (usually
C###            tout, but if you check the lsoda documentation you'll see that
C###            it could take other values)
C###          </li>
C###          <li><b>itol,rtol,atol</b> - these variables control the tolerance
C###            of the integration, the estimated local error in y(i) will be
C###            controlled so as to be less than:
C###            <ul>
C###              <li>ewt(i) = rtol*abs(y(i)) + atol     if itol = 1, or
C###              </li>
C###              <li>ewt(i) = rtol*abs(y(i)) + atol(i)  if itol = 2
C###              </li>
C###            </ul>
C###            Currently in CMISS, itol has a fixed value of 1 - i.e. the same
C###            absolute tolerance is applied to all state variables. It should
C###            be trivial to modify CMISS to enable to specification of
C###            individual tolerances
C###          </li>
C###          <li><b>itask</b> - is set to 1 for normal computation of output
C###            values of y at t = tout
C###          </li>
C###          <li><b>istate</b> - integer flag used for input and output, gives
C###            the type of error on return
C###          </li>
C###          <li><b>iopt</b> - is set to 1 in CMISS to specify that some
C###            optional input parameters are set. The optional parameters used
C###            in CMISS are the intial and maximum time step sizes, and the
C###            maximum number of steps allowed during one call to the solver. These
C###            are all set inside <a  href="http://www.esc.auckland.ac.nz/Groups/Bioengineering/CMISS/scripts/cmiss-lookup.cgi?name=march8">MARCH8</a>
C###          </li>
C###          <li><b>rwork,lrw,iwork,liw</b> - workspace allocated for
C###            integrator use, also contains the optional input parameters and
C###            has some interesting output information (check the lsoda
C###            documentation)
C###          </li>
C###          <li><b>jdum,jt</b> - in CMISS these are used to specify that LSODA
C###            should use an internally generated (difference quotient) full
C###            jacobian (using neq extra calls to f per df/dy value)
C###          </li>
C###          <li><b>aii,aio,control,model,sizes,variant,ari,aro,derived,parameters,protocol,err</b>
C###            - CMISS arrays passed though to the cell RHS routine (cell
C###            model).
C###          </li>
C###          <li><b>rrwork,iiwork</b> - extra workspace used to store the
C###            variables which used to be stored in internal common blocks
C###            inside LSODA.
C###          </li>
C###        </ul>
C###      </p>
C###    </html>

      implicit none

      INCLUDE 'lsoda00.cmn'

      ! Parameter list
      external func
      integer neq,itol,itask,istate,iopt,lrw,iwork(*),liw,jdum,jt,
     1  aii(*),aio(*),control(*),model(*),sizes(*),variant,err,iiwork(*)
      real*8 y(*),t,tout,rwork(*),ari(*),aro(*),derived(*),
     1  parameters(*),protocol(*),rrwork(*),rtol,atol
      character error*(*)
      ! Local variables
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     4  mxstep, mxhnil, nhnil, ntrep, nslast, nyh,
     5  ialth, ipup, lmax, nqnyh, nslp,
     6  icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     7  maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu,
     1  insufr, insufi, ixpr, icount, irflag, jtyp, mused, mxordn,
     3  mxords,
     1  mesflg, lunit,
     1  conit, crate, el, elco,
     1  hold, rmax, tesco,
     2  ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     1  tsw, pdest, pdlast, ratio, cm1, cm2,
     1  pdnorm

      ! Real common block variables
      parameter(conit=1)
      parameter(crate=2)
      parameter(el=3) !(13)
      parameter(elco=16) !(13,12)
      parameter(hold=172)
      parameter(rmax=173)
      parameter(tesco=174) !(3,12)
      parameter(ccmax=210)
      parameter(el0=211)
      parameter(h=212)
      parameter(hmin=213)
      parameter(hmxi=214)
      parameter(hu=215)
      parameter(rc=216)
      parameter(tn=217)
      parameter(uround=218)
      parameter(tsw=219)
      parameter(pdest=220)
      parameter(pdlast=221)
      parameter(ratio=222)
      parameter(cm1=223) !(12)
      parameter(cm2=235) !(5)
      parameter(pdnorm=240)
      ! Integer common block variables
      parameter(illin=1)
      parameter(init=2)
      parameter(lyh=3)
      parameter(lewt=4)
      parameter(lacor=5)
      parameter(lsavf=6)
      parameter(lwm=7)
      parameter(liwm=8)
      parameter(mxstep=9)
      parameter(mxhnil=10)
      parameter(nhnil=11)
      parameter(ntrep=12)
      parameter(nslast=13)
      parameter(nyh=14)
      parameter(ialth=15)
      parameter(ipup=16)
      parameter(lmax=17)
      parameter(nqnyh=18)
      parameter(nslp=19)
      parameter(icf=20)
      parameter(ierpj=21)
      parameter(iersl=22)
      parameter(jcur=23)
      parameter(jstart=24)
      parameter(kflag=25)
      parameter(l=26)
      parameter(meth=27)
      parameter(miter=28)
      parameter(maxord=29)
      parameter(maxcor=30)
      parameter(msbp=31)
      parameter(mxncf=32)
      parameter(n=33)
      parameter(nq=34)
      parameter(nst=35)
      parameter(nfe=36)
      parameter(nje=37)
      parameter(nqu=38)
      parameter(insufr=39)
      parameter(insufi=40)
      parameter(ixpr=41)
      parameter(icount=42)
      parameter(irflag=43)
      parameter(jtyp=44)
      parameter(mused=45)
      parameter(mxordn=46)
      parameter(mxords=47)
      parameter(mesflg=48)
      parameter(lunit=49)

      iiwork(mesflg) = 1
      iiwork(lunit) = 6
      istate = 1

      call lsoda(func,neq,y,t,tout,itol,rtol,atol,itask,
     1  istate,
     1  iopt,rwork,lrw,iwork,liw,jdum,jt,
     2  aii,aio,control,model,sizes,variant,
     3  ari,aro,derived,parameters,protocol,err,
     1  rrwork(conit),rrwork(crate),rrwork(el),rrwork(elco),
     1  rrwork(hold),rrwork(rmax),rrwork(tesco),
     2  rrwork(ccmax),rrwork(el0),rrwork(h),rrwork(hmin),rrwork(hmxi),
     1  rrwork(hu),rrwork(rc),
     1  rrwork(tn),rrwork(uround),
     3  iiwork(illin),iiwork(init),iiwork(lyh),iiwork(lewt),
     1  iiwork(lacor),
     1  iiwork(lsavf),
     1  iiwork(lwm),iiwork(liwm),
     4  iiwork(mxstep),iiwork(mxhnil),iiwork(nhnil),iiwork(ntrep),
     1  iiwork(nslast),
     1  iiwork(nyh),
     5  iiwork(ialth),iiwork(ipup),iiwork(lmax),iiwork(nqnyh),
     1  iiwork(nslp),
     6  iiwork(icf),iiwork(ierpj),iiwork(iersl),iiwork(jcur),
     1  iiwork(jstart),
     1  iiwork(kflag),
     1  iiwork(l),iiwork(meth),iiwork(miter),
     7  iiwork(maxord),iiwork(maxcor),iiwork(msbp),iiwork(mxncf),
     1  iiwork(n),
     1  iiwork(nq),
     1  iiwork(nst),iiwork(nfe),iiwork(nje),iiwork(nqu),
     1  rrwork(tsw),rrwork(pdest),rrwork(pdlast),rrwork(ratio),
     1  rrwork(cm1),
     1  rrwork(cm2),
     1  rrwork(pdnorm),
     2  iiwork(insufr),iiwork(insufi),iiwork(ixpr),iiwork(icount),
     1  iiwork(irflag),
     1  iiwork(jtyp),iiwork(mused),iiwork(mxordn),
     3  iiwork(mxords),
     1  iiwork(mesflg),iiwork(lunit))

      IF(istate.LT.0) THEN
        IF(istate.EQ.-1) THEN
          ERROR='Too many steps taken'
          GOTO 9999
        ELSEIF(istate.EQ.-2) THEN
          ERROR='Precision too high'
          GOTO 9999
        ELSEIF(istate.EQ.-3) THEN
          WRITE(ERROR,'(A,E12.5)') 'Illegal input, tolsf =',rwork(14)
          GOTO 9999
        ELSEIF(istate.EQ.-4) THEN
          ERROR='Repeated error tests'
          GOTO 9999
        ELSEIF(istate.EQ.-5) THEN
          ERROR='Convergence test failures'
          GOTO 9999
        ELSEIF(istate.EQ.-6) THEN
          ERROR='EWT became zero'
          GOTO 9999
        ELSEIF(istate.EQ.-7) THEN
          ERROR='Length of RWORK and/or IWORK too small'
          GOTO 9999
        ENDIF
      ENDIF
      IF(err.NE.0) THEN
        ERROR='Error in RHSROUTINE inside lsoda'
        GOTO 9999
      ENDIF

      CALL EXITS('LSODA_WRAPPER')
      RETURN
 9999 CALL ERRORS('LSODA_WRAPPER',ERROR)
      CALL EXITS('LSODA_WRAPPER')
      RETURN 1
      END

!       SUBROUTINE LSODA_WRAPPER(AII,AIO,CONTROL,ERR,MODEL,
!      '  NEQ,SIZES,VARIANT,ARI,ARO,DERIVED,PARAMETERS,
!      '  PROTOCOL,T,TOUT,Y,F,ERROR,*)

! C#### Subroutine: LSODA_WRAPPER
! C###  Description:
! C###  Provides a CMISS wrapper to the lsoda ode solver
! C###  Added by Richard Boyes 8th december 2000 upon
! C###  request from physiome sciences. The lsoda package
! C###  exists in an external library that is compiled into
! C###  CMISS externally

!       IMPLICIT NONE
! !     Parameter list
!       INCLUDE 'cmiss$reference:cbdi02.cmn'
!       INCLUDE 'cmiss$reference:lsoda00.cmn'
!       INCLUDE 'cmiss$reference:integrator.cmn'
!       INTEGER AII(*),AIO(*),CONTROL(*),ERR,MODEL(*),NEQ,SIZES(*),
!      '  VARIANT
!       REAL*8 ARI(*),ARO(*),DERIVED(*),PARAMETERS(*),
!      '  PROTOCOL(*),T,TOUT,Y(*)
!       EXTERNAL F
!       CHARACTER ERROR*(*)
! !     Local variables
!       INTEGER MAX_NUM_VARS,LSODA_LIWORK_MAX,LSODA_LWORK_MAX
!       PARAMETER(MAX_NUM_VARS=100,LSODA_LIWORK_MAX=100,
!      '  LSODA_LWORK_MAX=5500)
!       INTEGER i,LSODA_IWORK(LSODA_LIWORK_MAX)
!       REAL*8 JAC,LSODA_WORK(LSODA_LWORK_MAX),TSTORE

!       CALL ENTERS('LSODA_WRAPPER',*9999)

!       IF(LSODA_LIWORK_MAX.LT.INTEGRATOR_LIWORK) THEN
!         ERROR='>>Increase LSODA_LIWORK_MAX parameter '
!      '    //'in LSODA_WRAPPER'
!         GOTO 9999
!       ENDIF
!       IF(LSODA_LWORK_MAX.LT.INTEGRATOR_LWORK) THEN
!         ERROR='>>Increase LSODA_LWORK_MAX parameter '
!      '    //'in LSODA_WRAPPER'
!         GOTO 9999
!       ENDIF

!       DO i=1,20 ! Only need to initialise optional inputs
!          LSODA_IWORK(i)=0
!       ENDDO
!       DO i=1,22 ! Only need to initialise optional inputs
!          LSODA_WORK(i)=0.d0
!       ENDDO

! !     Break up lsoda into an initialisation
! !     call and then a work call.
!       TSTORE=TOUT
!       TOUT=T
!       LSODA_ISTATE=1
!       LSODA_IOPT=0
!       CALL LSODA(F,NEQ,Y,T,TOUT,LSODA_ERROR_TYPE,
!      '  LSODA_REL_ERR,LSODA_ABS_ERR,LSODA_ITASK,
!      '  LSODA_ISTATE,LSODA_IOPT,LSODA_WORK,INTEGRATOR_LWORK,
!      '  LSODA_IWORK,INTEGRATOR_LIWORK,JAC,LSODA_JACOBIAN_TYPE,
!      '  AII,AIO,CONTROL,MODEL,SIZES,VARIANT,ARI,ARO,DERIVED,
!      '  PARAMETERS,PROTOCOL,ERR)
!       IF(LSODA_ISTATE.LT.0) THEN
!          IF(LSODA_ISTATE.EQ.-1) THEN
!             ERROR='Too many steps taken'
!             GOTO 9999
!          ELSEIF(LSODA_ISTATE.EQ.-2) THEN
!             ERROR='Precision too high'
!             GOTO 9999
!          ELSEIF(LSODA_ISTATE.EQ.-3) THEN
!             ERROR='Illegal input'
!             GOTO 9999
!          ELSEIF(LSODA_ISTATE.EQ.-4) THEN
!             ERROR='Repeated error tests'
!             GOTO 9999
!          ELSEIF(LSODA_ISTATE.EQ.-5) THEN
!             ERROR='Convergence test failures'
!             GOTO 9999
!          ELSEIF(LSODA_ISTATE.EQ.-6) THEN
!             ERROR='EWT became zero'
!             GOTO 9999
!          ELSEIF(LSODA_ISTATE.EQ.-7) THEN
!             ERROR='Length of RWORK and/or IWORK too small'
!             GOTO 9999
!          ENDIF
!       ENDIF
!       IF(ERR.NE.0) THEN
!         ERROR='Error in RHSROUTINE inside lsoda'
!         GOTO 9999
!       ENDIF
!       LSODA_IOPT=1
!       LSODA_WORK(5)=LSODA_MAX_STEP
!       LSODA_WORK(6)=LSODA_MAX_STEP
!       LSODA_WORK(7)=0.d0
!       LSODA_IWORK(6)=LSODA_MAX_ITERS
!       TOUT=TSTORE
!       CALL LSODA(F,NEQ,Y,T,TOUT,LSODA_ERROR_TYPE,LSODA_REL_ERR,
!      '  LSODA_ABS_ERR,LSODA_ITASK,LSODA_ISTATE,LSODA_IOPT,
!      '  LSODA_WORK,INTEGRATOR_LWORK,LSODA_IWORK,INTEGRATOR_LIWORK,
!      '  JAC,LSODA_JACOBIAN_TYPE,AII,AIO,CONTROL,MODEL,SIZES,VARIANT,
!      '  ARI,ARO,DERIVED,PARAMETERS,PROTOCOL,ERR)
!       IF(LSODA_ISTATE.LT.0) THEN
!         IF(LSODA_ISTATE.EQ.-1) THEN
!           ERROR='Too many steps taken'
!           GOTO 9999
!         ELSEIF(LSODA_ISTATE.EQ.-2) THEN
!           ERROR='Precision too high'
!           GOTO 9999
!         ELSEIF(LSODA_ISTATE.EQ.-3) THEN
!           ERROR='Illegal input'
!           GOTO 9999
!         ELSEIF(LSODA_ISTATE.EQ.-4) THEN
!           ERROR='Repeated error tests'
!           GOTO 9999
!         ELSEIF(LSODA_ISTATE.EQ.-5) THEN
!           ERROR='Convergence test failures'
!           GOTO 9999
!         ELSEIF(LSODA_ISTATE.EQ.-6) THEN
!           ERROR='EWT became zero'
!           GOTO 9999
!         ELSEIF(LSODA_ISTATE.EQ.-7) THEN
!           ERROR='Length of RWORK and/or IWORK too small'
!           GOTO 9999
!         ENDIF
!       ENDIF
!       IF(ERR.NE.0) THEN
!         ERROR='Error in RHSROUTINE inside lsoda'
!         GOTO 9999
!       ENDIF

!       CALL EXITS('LSODA_WRAPPER')
!       RETURN
!  9999 CALL ERRORS('LSODA_WRAPPER',ERROR)
!       CALL EXITS('LSODA_WRAPPER')
!       RETURN 1
!       END
