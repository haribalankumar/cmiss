      SUBROUTINE FPCURF(IOPT,X,Y,W,M,XB,XE,K,S,NEST,TOL,MAXIT,K1,K2,
     ' N,T,C,FP,FPINT,Z,A,B,G,Q,NRDATA,IER)
C#### Subroutine: FPCURF
C###  Description:
C###    Auxililary routine for CURFIT.
c  ..
      IMPLICIT NONE
c  ..scalar arguments..
      REAL*8 XB,XE,S,TOL,FP
      INTEGER IOPT,M,K,NEST,MAXIT,K1,K2,N,IER
c  ..array arguments..
      REAL*8 X(M),Y(M),W(M),T(NEST),C(NEST),FPINT(NEST),
     ' Z(NEST),A(NEST,K1),B(NEST,K2),G(NEST,K2),Q(M,K1)
      INTEGER NRDATA(NEST)
c  ..local scalars..
      REAL*8 ACC,CON1,CON4,CON9,COS,HALF,FPART,FPMS,FPOLD,FP0,F1,F2,F3,
     ' ONE,P,PINV,PIV,P1,P2,P3,RN,SIN,STORE,TERM,WI,XI,YI
      INTEGER I,ICH1,ICH3,IT,ITER,I1,I2,I3,J,K3,L,L0,
     ' MK1,NEW,NK1,NMAX,NMIN,NPLUS,NPL1,NRINT,N8
c  ..local arrays..
      REAL*8 H(7)
c  ..function references
      REAL*8 ABS,FPRATI
      INTEGER MAX0,MIN0
c  ..subroutine references..
c    fpback,fpbspl,fpgivs,fpdisc,fpknot,fprota
c  ..
c  set constants
      ONE = 0.1D+01
      CON1 = 0.1d0
      CON9 = 0.9d0
      CON4 = 0.4D-01
      HALF = 0.5d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  part 1: determination of the number of knots and their position     c
c  **************************************************************      c
c  given a set of knots we compute the least-squares spline sinf(x),   c
c  and the corresponding sum of squared residuals fp=f(p=inf).         c
c  if iopt=-1 sinf(x) is the requested approximation.                  c
c  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
c    if fp <=s we will continue with the current set of knots.         c
c    if fp > s we will increase the number of knots and compute the    c
c       corresponding least-squares spline until finally fp<=s.        c
c    the initial choice of knots depends on the value of s and iopt.   c
c    if s=0 we have spline interpolation; in that case the number of   c
c    knots equals nmax = m+k+1.                                        c
c    if s > 0 and                                                      c
c      iopt=0 we first compute the least-squares polynomial of         c
c      degree k; n = nmin = 2*k+2                                      c
c      iopt=1 we start with the set of knots found at the last         c
c      call of the routine, except for the case that s > fp0; then     c
c      we compute directly the least-squares polynomial of degree k.   c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  determine nmin, the number of knots for polynomial approximation.
      NMIN = 2*K1
      IF(IOPT.LT.0) GO TO 60
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      ACC = TOL*S
c  determine nmax, the number of knots for spline interpolation.
      NMAX = M+K1
      IF(S.GT.0.0d0) GO TO 45
c  if s=0, s(x) is an interpolating spline.
c  test whether the required storage space exceeds the available one.
      N = NMAX
      IF(NMAX.GT.NEST) GO TO 420
c  find the position of the interior knots in case of interpolation.
  10  MK1 = M-K1
      IF(MK1.EQ.0) GO TO 60
      K3 = K/2
      I = K2
      J = K3+2
      IF(K3*2.EQ.K) GO TO 30
      DO 20 L=1,MK1
        T(I) = X(J)
        I = I+1
        J = J+1
  20  CONTINUE
      GO TO 60
  30  DO 40 L=1,MK1
        T(I) = (X(J)+X(J-1))*HALF
        I = I+1
        J = J+1
  40  CONTINUE
      GO TO 60
c  if s>0 our initial choice of knots depends on the value of iopt.
c  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
c  polynomial of degree k which is a spline without interior knots.
c  if iopt=1 and fp0>s we start computing the least squares spline
c  according to the set of knots found at the last call of the routine.
  45  IF(IOPT.EQ.0) GO TO 50
      IF(N.EQ.NMIN) GO TO 50
      FP0 = FPINT(N)
      FPOLD = FPINT(N-1)
      NPLUS = NRDATA(N)
      IF(FP0.GT.S) GO TO 60
  50  N = NMIN
      FPOLD = 0.0d0
      NPLUS = 0
      NRDATA(1) = M-2
c  main loop for the different sets of knots. m is a save upper bound
c  for the number of trials.
  60  DO 200 ITER = 1,M
        IF(N.EQ.NMIN) IER = -2
c  find nrint, tne number of knot intervals.
        NRINT = N-NMIN+1
c  find the position of the additional knots which are needed for
c  the b-spline representation of s(x).
        NK1 = N-K1
        I = N
        DO 70 J=1,K1
          T(J) = XB
          T(I) = XE
          I = I-1
  70    CONTINUE
c  compute the b-spline coefficients of the least-squares spline
c  sinf(x). the observation matrix a is built up row by row and
c  reduced to upper triangular form by givens transformations.
c  at the same time fp=f(p=inf) is computed.
        FP = 0.0d0
c  initialize the observation matrix a.
        DO 80 I=1,NK1
          Z(I) = 0.0d0
          DO 80 J=1,K1
            A(I,J) = 0.0d0
  80    CONTINUE
        L = K1
        DO 130 IT=1,M
c  fetch the current data point x(it),y(it).
          XI = X(IT)
          WI = W(IT)
          YI = Y(IT)*WI
c  search for knot interval t(l) <= xi < t(l+1).
  85      IF(XI.LT.T(L+1) .OR. L.EQ.NK1) GO TO 90
          L = L+1
          GO TO 85
c  evaluate the (k+1) non-zero b-splines at xi and store them in q.
  90      CALL FPBSPL(T,N,K,XI,L,H)
          DO 95 I=1,K1
            Q(IT,I) = H(I)
            H(I) = H(I)*WI
  95      CONTINUE
c  rotate the new row of the observation matrix into triangle.
          J = L-K1
          DO 110 I=1,K1
            J = J+1
            PIV = H(I)
            IF(PIV.EQ.0.0d0) GO TO 110
c  calculate the parameters of the givens transformation.
            CALL FPGIVS(PIV,A(J,1),COS,SIN)
c  transformations to right hand side.
            CALL FPROTA(COS,SIN,YI,Z(J))
            IF(I.EQ.K1) GO TO 120
            I2 = 1
            I3 = I+1
            DO 100 I1 = I3,K1
              I2 = I2+1
c  transformations to left hand side.
              CALL FPROTA(COS,SIN,H(I1),A(J,I2))
 100        CONTINUE
 110      CONTINUE
c  add contribution of this row to the sum of squares of residual
c  right hand sides.
 120      FP = FP+YI**2
 130    CONTINUE
        IF(IER.EQ.(-2)) FP0 = FP
        FPINT(N) = FP0
        FPINT(N-1) = FPOLD
        NRDATA(N) = NPLUS
c  backward substitution to obtain the b-spline coefficients.
        CALL FPBACK(A,Z,NK1,K1,C,NEST)
c  test whether the approximation sinf(x) is an acceptable solution.
        IF(IOPT.LT.0) GO TO 440
        FPMS = FP-S
        IF(ABS(FPMS).LT.ACC) GO TO 440
c  if f(p=inf) < s accept the choice of knots.
        IF(FPMS.LT.0.0d0) GO TO 250
c  if n = nmax, sinf(x) is an interpolating spline.
        IF(N.EQ.NMAX) GO TO 430
c  increase the number of knots.
c  if n=nest we cannot increase the number of knots because of
c  the storage capacity limitation.
        IF(N.EQ.NEST) GO TO 420
c  determine the number of knots nplus we are going to add.
        IF(IER.EQ.0) GO TO 140
        NPLUS = 1
        IER = 0
        GO TO 150
 140    NPL1 = NPLUS*2
        RN = NPLUS
        IF(FPOLD-FP.GT.ACC) NPL1 = INT(RN*FPMS/(FPOLD-FP))
        NPLUS = MIN0(NPLUS*2,MAX0(NPL1,NPLUS/2,1))
 150    FPOLD = FP
c  compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
c  t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        FPART = 0.0d0
        I = 1
        L = K2
        NEW = 0
        DO 180 IT=1,M
          IF(X(IT).LT.T(L) .OR. L.GT.NK1) GO TO 160
          NEW = 1
          L = L+1
 160      TERM = 0.0d0
          L0 = L-K2
          DO 170 J=1,K1
            L0 = L0+1
            TERM = TERM+C(L0)*Q(IT,J)
 170      CONTINUE
          TERM = (W(IT)*(TERM-Y(IT)))**2
          FPART = FPART+TERM
          IF(NEW.EQ.0) GO TO 180
          STORE = TERM*HALF
          FPINT(I) = FPART-STORE
          I = I+1
          FPART = STORE
          NEW = 0
 180    CONTINUE
        FPINT(NRINT) = FPART
        DO 190 L=1,NPLUS
c  add a new knot.
          CALL FPKNOT(X,M,T,N,FPINT,NRDATA,NRINT,NEST,1)
c  if n=nmax we locate the knots as for interpolation.
          IF(N.EQ.NMAX) GO TO 10
c  test whether we cannot further increase the number of knots.
          IF(N.EQ.NEST) GO TO 200
 190    CONTINUE
c  restart the computations with the new set of knots.
 200  CONTINUE
c  test whether the least-squares kth degree polynomial is a solution
c  of our approximation problem.
 250  IF(IER.EQ.(-2)) GO TO 440
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  part 2: determination of the smoothing spline sp(x).                c
c  ***************************************************                 c
c  we have determined the number of knots and their position.          c
c  we now compute the b-spline coefficients of the smoothing spline    c
c  sp(x). the observation matrix a is extended by the rows of matrix   c
c  b expressing that the kth derivative discontinuities of sp(x) at    c
c  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
c  ponding weights of these additional rows are set to 1/p.            c
c  iteratively we then have to determine the value of p such that      c
c  f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that    c
c  the least-squares kth degree polynomial corresponds to p=0, and     c
c  that the least-squares spline corresponds to p=infinity. the        c
c  iteration process which is proposed here, makes use of rational     c
c  interpolation. since f(p) is a convex and strictly decreasing       c
c  function of p, it can be approximated by a rational function        c
c  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
c  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
c  to calculate the new value of p such that r(p)=s. convergence is    c
c  guaranteed by taking f1>0 and f3<0.                                 c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  evaluate the discontinuity jump of the kth derivative of the
c  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
      CALL FPDISC(T,N,K2,B,NEST)
c  initial value for p.
      P1 = 0.0d0
      F1 = FP0-S
      P3 = -ONE
      F3 = FPMS
      P = 0.0d0
      DO 255 I=1,NK1
         P = P+A(I,1)
 255  CONTINUE
      RN = NK1
      P = RN/P
      ICH1 = 0
      ICH3 = 0
      N8 = N-NMIN
c  iteration process to find the root of f(p) = s.
      DO 360 ITER=1,MAXIT
c  the rows of matrix b with weight 1/p are rotated into the
c  triangularised observation matrix a which is stored in g.
        PINV = ONE/P
        DO 260 I=1,NK1
          C(I) = Z(I)
          G(I,K2) = 0.0d0
          DO 260 J=1,K1
            G(I,J) = A(I,J)
 260    CONTINUE
        DO 300 IT=1,N8
c  the row of matrix b is rotated into triangle by givens transformation
          DO 270 I=1,K2
            H(I) = B(IT,I)*PINV
 270      CONTINUE
          YI = 0.0d0
          DO 290 J=IT,NK1
            PIV = H(1)
c  calculate the parameters of the givens transformation.
            CALL FPGIVS(PIV,G(J,1),COS,SIN)
c  transformations to right hand side.
            CALL FPROTA(COS,SIN,YI,C(J))
            IF(J.EQ.NK1) GO TO 300
            I2 = K1
            IF(J.GT.N8) I2 = NK1-J
            DO 280 I=1,I2
c  transformations to left hand side.
              I1 = I+1
              CALL FPROTA(COS,SIN,H(I1),G(J,I1))
              H(I) = H(I1)
 280        CONTINUE
            H(I2+1) = 0.0d0
 290      CONTINUE
 300    CONTINUE
c  backward substitution to obtain the b-spline coefficients.
        CALL FPBACK(G,C,NK1,K2,C,NEST)
c  computation of f(p).
        FP = 0.0d0
        L = K2
        DO 330 IT=1,M
          IF(X(IT).LT.T(L) .OR. L.GT.NK1) GO TO 310
          L = L+1
 310      L0 = L-K2
          TERM = 0.0d0
          DO 320 J=1,K1
            L0 = L0+1
            TERM = TERM+C(L0)*Q(IT,J)
 320      CONTINUE
          FP = FP+(W(IT)*(TERM-Y(IT)))**2
 330    CONTINUE
c  test whether the approximation sp(x) is an acceptable solution.
        FPMS = FP-S
        IF(ABS(FPMS).LT.ACC) GO TO 440
c  test whether the maximal number of iterations is reached.
        IF(ITER.EQ.MAXIT) GO TO 400
c  carry out one more step of the iteration process.
        P2 = P
        F2 = FPMS
        IF(ICH3.NE.0) GO TO 340
        IF((F2-F3).GT.ACC) GO TO 335
c  our initial choice of p is too large.
        P3 = P2
        F3 = F2
        P = P*CON4
        IF(P.LE.P1) P=P1*CON9 + P2*CON1
        GO TO 360
 335    IF(F2.LT.0.0d0) ICH3=1
 340    IF(ICH1.NE.0) GO TO 350
        IF((F1-F2).GT.ACC) GO TO 345
c  our initial choice of p is too small
        P1 = P2
        F1 = F2
        P = P/CON4
        IF(P3.LT.0.0d0) GO TO 360
        IF(P.GE.P3) P = P2*CON1 + P3*CON9
        GO TO 360
 345    IF(F2.GT.0.0d0) ICH1=1
c  test whether the iteration process proceeds as theoretically
c  expected.
 350    IF(F2.GE.F1 .OR. F2.LE.F3) GO TO 410
c  find the new value for p.
        P = FPRATI(P1,F1,P2,F2,P3,F3)
 360  CONTINUE
c  error codes and messages.
 400  IER = 3
      GO TO 440
 410  IER = 2
      GO TO 440
 420  IER = 1
      GO TO 440
 430  IER = -1
 440  RETURN
      END

