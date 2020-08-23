      subroutine aux1r6(a,ia,ja,u,f,ig,
     +                  nda,ndia,ndja,ndu,ndf,ndig,nnu,matrix,
     +                  eps,ifirst,iswtch,iout,iprint,
     +                  ierr)
c
c      ------------------------------------------------------------
c      ] interface to amg-module for solving linear systems l*u=f ]
c      ------------------------------------------------------------
c
c         --------------------------------------------------------------
c
c     assumptions on l:
c
c         the program requires:
c
c             - diagonal entries are always positive (on all grids);
c             - l is a square matrix which is either regular or singular
c               with rowsums=0.
c
c         for theoretical reasons the following should hold:
c
c             - l positive definite (or semi-definite with rowsum=0)
c             - l "essentially" positive type, i.e.,
c
c                  -- diagonal entries must be > 0 ;
c                  -- most of the off-diagonal entries <= 0 ;
c                  -- rowsums should be >= 0 .
c
c
c     the user has to provide the matrix l, the right hand side f and
c     certain pointer vectors ia and ja.
c
c         --------------------------------------------------------------
c
c     storage of l:
c
c         the non-zero entries of the matrix l are stored in
c         "compressed" sky-line fashion in a 1-d vector a, i.e., row
c         after row, each row starting with its diagonal element. the
c         other non-zero row entries follow their diagonal entry in any
c         order.
c
c         in order to identify each element in a, the user has to
c         provide two pointer arrays ia and ja. if nnu denotes the total
c         number of unknowns, the non-zero entries of any row i of l
c         (1.le.i.le.nnu) are stored in a(j) where the range of j
c         is given by
c
c                     ia(i) .le. j .le. ia(i+1)-1.
c
c         thus, ia(i) points to the position of the diagonal entry of
c         row i within the vector a. in particular,
c
c                     ia(1) = 1 ,  ia(nnu+1) = 1 + nna
c
c         where nna denotes the total number of matrix entries stored.
c         the pointer vector ja has to be defined such that
c         any entry a(j) corresponds to the unknown u(ja(j)), i.e.,
c         ja(j) points to the column index of a(j).
c         in particular, a(ia(i)) is the diagonal entry of row i
c         and corresponds to the unknown u(i): ja(ia(i))=i.
c
c         in this terminology, the i-th equation reads as follows
c         (for any i with  1.le.i.le.nnu):
c
c                  f(i) =        sum      a(j) * u(ja(j))
c                           j1.le.j.le.j2
c
c         where f(i) denotes the i-th component of the right hand
c         side and
c
c                     j1 = ia(i) ,  j2 = ia(i+1)-1.
c
c         notes: the entry ia(nnu+1) has to point to the first free
c                entry in vectors a and ja, respectively. otherwise,
c                amg cannot know the length of the last matrix row.
c
c                the input vectors a, ia and ja are changed by amg1r6.
c                so, after return from amg1r6, the package must not
c                be called a second time without having newly defined
c                the input vectors and using iswtch=4. otherwise, the
c                setup phase will fail.
c                  on the other hand, running amg a second time on the
c                same input data with iswtch=4 has no sense, because
c                the results of the first setup phase are still stored
c                and thus this phase can be skipped in a second call.
c                in order to do this, set iswtch to 1, 2 or 3.
c
c
c-----------------------------------------------------------------------
c
c         the form of the calling program has to be as follows:
c
c               program driver
c         c
c               real*8 a(#nda),u(#ndu),f(#ndf)
c               integer ia(#ndia),ja(#ndja),ig(#ndig)
c         c
c               nda  = #nda
c               ndu  = #ndu
c               ndf  = #ndf
c               ndia = #ndia
c               ndja = #ndja
c               ndig = #ndig
c         c
c         c     set up a, f, ia, ja and specify necessary parameters
c         c
c               ....
c               ....
c         c
c               call aux1r6(a,ia,ja,u,f,ig,
c        +                  nda,ndia,ndja,ndu,ndf,ndig,nnu,matrix,
c        +                  eps,ifirst,iswtch,iout,iprint,
c        +                  ierr)
c         c
c               ....
c               ....
c         c
c               stop
c               end
c
c-----------------------------------------------------------------------
c
c     input via arrays (see above):
c
c     a        -   matrix l
c
c     ia       -   pointer vector
c
c     ja       -   pointer vector
c
c     u        -   first approximation to solution
c
c     f        -   right hand side
c
c
c-----------------------------------------------------------------------
c
c
c     scalar input parameters of aux1r6:
c
c     the input parameters of aux1r6 in the list below are arranged
c     according to their importance to the general user. the parameters
c     preceeded by a * must be specified explicitely (definition of the
c     user-defined problem and dimensioning of vectors in the calling
c     program). the other parameters are set to standard values if zero
c     on input.
c
c
c  *  nda      -   dimensioning of vector a in calling program
c
c  *  ndia     -   dimensioning of vector ia in calling program
c
c  *  ndja     -   dimensioning of vector ja in calling program
c
c  *  ndu      -   dimensioning of vector u in calling program
c
c  *  ndf      -   dimensioning of vector f in calling program
c
c  *  ndig     -   dimensioning of vector ig in calling program
c
c  *  nnu      -   number of unknowns
c
c  *  matrix   -   integer value containing info about the matrix l.
c
c                  1st digit of matrix  --  isym:
c                    =1: l is symmetric;
c                    =2: l is not symmetric.
c
c                  2nd digit of matrix  --  irow0:
c                    =1: l has rowsum zero;
c                    =2: l does not have rowsum zero.
c
c
c  *  eps      -   convergence criterion for solution process. stop, if
c                  l2-norm of the residual of the user-defined problem
c                  is less than eps.
c
c     ifirst   -   parameter for first approximation.
c
c                  1st digit of ifirst: not used; has to be non-zero.
c
c                  2nd digit of ifirst  --  itypu:
c                    =0: no setting of first approximation,
c                    =1: first approximation constant to zero,
c                    =2: first approximation constant to one,
c                    =3: first approximation is random function with
c                        the concrete random sequence being determined
c                        by the follwing digits.
c
c                  rest of ifirst  --  rndu:
c                    determines the concrete random sequence used in
c                    the case itypu=3. (ifirst=13 is equivalent to
c                    ifirst=1372815)
c
c     iswtch   -   parameter controlling which modules of amg1r6 are to
c                  be used.
c                    =1:   call for -----, -----, -----, wrkcnt.
c                    =2:   call for -----, -----, solve, wrkcnt.
c                    =3:   call for -----, first, solve, wrkcnt.
c                    =4:   call for setup, first, solve, wrkcnt.
c                  setup defines the operators needed in the solution
c                         phase.
c                  first initializes the solution vector (see parameter
c                         ifirst).
c                  solve computes the solution by amg cycling (see
c                         amg1r6-parameter ncyc)
c                  wrkcnt provides the user with information about
c                         residuals, storage requirements and cp-times
c                         (see parameter iout).
c                  if aux1r6 is called the first time, iswtch has to
c                  be =4. independent of iswtch, single modules can be
c                  bypassed by a proper choice of the corresponding
c                  parameter.
c
c     iout     -   parameter controlling the amount of output during
c                  solution phase:
c
c                  1st digit: not used; has to be non-zero.
c
c                  2nd digit:
c                    =0: no output (except for messages)
c                    =1: residual before and after solution process
c                    =2: add.: statistics on cp-times and storage requi-
c                        rements
c                    =3: add.: residual after each amg-cycle
c
c     iprint   -   parameter specifying the fortran unit numbers for
c                  output:
c
c                  1st digit: not used; has to be non-zero
c
c                  2nd and 3rd digit  --  iup: unit number for results
c
c                  4th and 5th digit  --  ium: unit number for messages
c
c
c-----------------------------------------------------------------------
c
c     output:
c
c     u        -   contains the computed solution
c
c
c     ierr     -   error parameter:
c
c                    >0: fatal error (abnormal termination)
c                    <0: non-fatal error (execution continues)
c
c                  error codes in detail:
c
c                  1. dimensioning too small for vector
c                        a      (ierr = 1)
c                        ia     (ierr = 2)
c                        ja     (ierr = 3)
c                        u      (ierr = 4)
c                        f      (ierr = 5)
c                        ig     (ierr = 6)
c
c                     no yale-smp because of storage (nda too small):
c                               (ierr = -1)
c                     no yale-smp because of storage (ndja too small):
c                               (ierr = -3)
c                     no cg because of storage (ndu too small):
c                               (ierr = -4)
c                     no space for transpose of interpolation (nda or
c                                                     ndja too small):
c                               (ierr = -1)
c
c                  2. input data erroneous:
c
c                     a-entry missing, isym = 1:           (ierr = -11)
c                     parameter matrix may be erroneous:   (ierr = -12)
c                     diagonal element not stored first:   (ierr =  13)
c                     diagonal element not positiv:        (ierr =  14)
c                     pointer ia erroneous:                (ierr =  15)
c                     pointer ja erroneous:                (ierr =  16)
c                     parameter iswtch erroneous:          (ierr =  17)
c
c                  3. errors of the amg1r6-system (should not occur):
c
c                     transpose a-entry missing:           (ierr =  21)
c                     interpolation entry missing:         (ierr =  22)
c
c                  4. algorithmic errors:
c
c                     cg-correction not defined:           (ierr =  31)
c                     no yale-smp because of error in
c                     factorization:                       (ierr = -32)
c
c-----------------------------------------------------------------------
c
c     work space:
c
c     the integer vector ig has to be passed to aux1r6 as work space.
c
c-----------------------------------------------------------------------
c
c     dimensioning of input vectors and work space:
c
c     it's impossible to tell in advance the exact storage requirements
c     of amg. the following formulas thus give only reasonable guesses
c     for the required vector lengths as declared in the calling pro-
c     gram. in these formulas nna denotes the number of non-zero entries
c     in the input-matrix l and nnu is the number of unknowns.
c
c     vector         needed length (guess)
c       a               3*nna + 5*nnu
c       ja              3*nna + 5*nnu
c       ia              2.2*nnu
c       u               2.2*nnu
c       f               2.2*nnu
c       ig              5.4*nnu
c
c-----------------------------------------------------------------------
c
c
c     standard choices of parameters:
c
c          ifirst = 13
c          iswtch = 4
c          iout   = 12
c          iprint = 10606
c
c     if any one of these parameters is 0 on input, its corresponding
c     standard value is used by aux1r6.
c
c-----------------------------------------------------------------------
c
c     portability restrictions:
c
c     1. routine AMGTIME is machine dependent and has to be adapted to
c        your computer installation or replaced by a dummy routine.
c
c     2. the amg1r6 system uses integer parameters of up to five digits.
c        be sure that your computer can store five digits on an integer
c        variable.
c
c     3. apart from fortran intrinsic functions and service routines,
c        there is only one external reference to a program not contained
c        in the amg1r6 system, i.e. the linear system solver ndrv of
c        the yale sparse matrix package. if you havn't access to this
c        package, enter a dummy routine ndrv and replace the line
c             nsolco = 2
c        by
c             nsolco = 1
c        in aux1r6. then ndrv isn't called by aux1r6. in this case,
c        however, indefinite problems will not be solvable.
c          the yale sparse matrix package is freely available for non-
c        profit purposes. contact the department of computer science,
c        yale unitversity.
c
c-----------------------------------------------------------------------
c
c     authors:
c
c          john ruge, fort collins (usa),
c              institute for computational studies at csu;
c
c          klaus stueben, d-5205 st. augustin (w.-germany),
c              gesellschaft fuer mathematik und datenverarbeitung (gmd).
c
c          rolf hempel, d-5205 st. augustin (w.-germany),
c              gesellschaft fuer mathematik und datenverarbeitung (gmd).
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      real*8 a(*),u(*),f(*)
      integer ia(*),ja(*),ig(*)
c
c===> for zero-parameters standard values are used by amg1r6
c
      levelx = 0
      ncyc   = 10250
      madapt = 0
      nrd    = 0
      nsolco = 2
      nru    = 0
cc>>>
c     ecg1   = 0.d0
      ecg1   = 1.0d-2
cc<<<
      ecg2   = 0.25d0
      ewt2   = 0.35d0
      nwt    = 2
      ntr    = 0
c
      iout   = 13
c
      call amg1r6(a,ia,ja,u,f,ig,
     +            nda,ndia,ndja,ndu,ndf,ndig,nnu,matrix,
     +            iswtch,iout,iprint,
     +            levelx,ifirst,ifmg,
     +            ncyc,eps,madapt,nrd,nsolco,nru,
     +            ecg1,ecg2,ewt2,nwt,ntr,
     +            ierr)
      end
