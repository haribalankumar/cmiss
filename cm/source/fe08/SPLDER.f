      SUBROUTINE SPLDER(T,N,C,K,NU,X,Y,M,WRK,IER)
C#### Subroutine: SPLDER
C###  Description:
C###    SPLDER evaluates the n-th derivative of the B-spline defined by
C###    CURFIT. This routine and all subseqent calling routines were
C###    written by Paul Dierckx and are a subset of the FITPACK library.
C###    This routine is equivalent to the NAG E02BCF as it is written by
C###    the same author. These routines appear exactly as the author wrote
C###    them except the precision has been changed from REAL to REAL*8.
C###    Below contains the comments from the author or refer to NAG
C###    manuals.
c  subroutine splder evaluates in a number of points x(i),i=1,2,...,m
c  the derivative of order nu of a spline s(x) of degree k,given in
c  its b-spline representation.
c
c  calling sequence:
c     call splder(t,n,c,k,nu,x,y,m,wrk,ier)
c
c  input parameters:
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, which contains the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    nu   : integer, specifying the order of the derivative. 0<=nu<=k
c    x    : array,length m, which contains the points where the deriv-
c           ative of s(x) must be evaluated.
c    m    : integer, giving the number of points where the derivative
c           of s(x) must be evaluated
c    wrk  : real*8 array of dimension n. used as working space.
c
c  output parameters:
c    y    : array,length m, giving the value of the derivative of s(x)
c           at the different points.
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    0 <= nu <= k
c    m >= 1
c    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
c
c  other subroutines required: fpbspl
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
      IMPLICIT NONE
c  ..scalar arguments..
      INTEGER N,K,NU,M,IER
c  ..array arguments..
      REAL*8 T(N),C(N),X(M),Y(M),WRK(N)
c  ..local scalars..
      INTEGER I,J,KK,K1,K2,L,LL,L1,L2,NK1,NK2
      REAL*8 AK,ARG,FAC,SP,TB,TE
c  ..local arrays ..
      REAL*8 H(6)
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      IER = 10
      IF(NU.LT.0 .OR. NU.GT.K) GO TO 200
      IF(M-1) 200,30,10
  10  DO 20 I=2,M
        IF(X(I).LT.X(I-1)) GO TO 200
  20  CONTINUE
  30  IER = 0
c  fetch tb and te, the boundaries of the approximation interval.
      K1 = K+1
      NK1 = N-K1
      TB = T(K1)
      TE = T(NK1+1)
c  the derivative of order nu of a spline of degree k is a spline of
c  degree k-nu,the b-spline coefficients wrk(i) of which can be found
c  using the recurrence scheme of de boor.
      L = 1
      KK = K
      DO 40 I=1,NK1
         WRK(I) = C(I)
  40  CONTINUE
      IF(NU.EQ.0) GO TO 100
      NK2 = NK1
      DO 60 J=1,NU
         AK = KK
         NK2 = NK2-1
         L1 = L
         DO 50 I=1,NK2
            L1 = L1+1
            L2 = L1+KK
            FAC = T(L2)-T(L1)
            IF(FAC.LE.0.0d0) GO TO 50
            WRK(I) = AK*(WRK(I+1)-WRK(I))/FAC
  50     CONTINUE
         L = L+1
         KK = KK-1
  60  CONTINUE
      IF(KK.NE.0) GO TO 100
c  if nu=k the derivative is a piecewise constant function
      J = 1
      DO 90 I=1,M
         ARG = X(I)
  70     IF(ARG.LT.T(L+1) .OR. L.EQ.NK1) GO TO 80
         L = L+1
         J = J+1
         GO TO 70
  80     Y(I) = WRK(J)
  90  CONTINUE
      GO TO 200
 100  L = K1
      L1 = L+1
      K2 = K1-NU
c  main loop for the different points.
      DO 180 I=1,M
c  fetch a new x-value arg.
        ARG = X(I)
        IF(ARG.LT.TB) ARG = TB
        IF(ARG.GT.TE) ARG = TE
c  search for knot interval t(l) <= arg < t(l+1)
 140    IF(ARG.LT.T(L1) .OR. L.EQ.NK1) GO TO 150
        L = L1
        L1 = L+1
        GO TO 140
c  evaluate the non-zero b-splines of degree k-nu at arg.
 150    CALL FPBSPL(T,N,KK,ARG,L,H)
c  find the value of the derivative at x=arg.
        SP = 0.0d0
        LL = L-K1
        DO 160 J=1,K2
          LL = LL+1
          SP = SP+WRK(LL)*H(J)
 160    CONTINUE
        Y(I) = SP
 180  CONTINUE
 200  RETURN
      END

