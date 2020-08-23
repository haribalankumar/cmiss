      SUBROUTINE FPKNOT(X,M,T,N,FPINT,NRDATA,NRINT,NEST,ISTART)
C#### Subroutine: FPKNOT
C###  Description:
C###    Auxililary routine for CURFIT.
c  subroutine fpknot locates an additional knot for a spline of degree
c  k and adjusts the corresponding parameters,i.e.
c    t     : the position of the knots.
c    n     : the number of knots.
c    nrint : the number of knotintervals.
c    fpint : the sum of squares of residual right hand sides
c            for each knot interval.
c    nrdata: the number of data points inside each knot interval.
c  istart indicates that the smallest data point at which the new knot
c  may be added is x(istart+1)
c  ..
      IMPLICIT NONE
c  ..scalar arguments..
      INTEGER M,N,NRINT,NEST,ISTART
c  ..array arguments..
      REAL*8 X(M),T(NEST),FPINT(NEST)
      INTEGER NRDATA(NEST)
c  ..local scalars..
      REAL*8 AN,AM,FPMAX
      INTEGER IHALF,J,JBEGIN,JJ,JK,JPOINT,K,MAXBEG,MAXPT,
     ' NEXT,NRX,NUMBER
c  ..
      K = (N-NRINT-1)/2
c  search for knot interval t(number+k) <= x <= t(number+k+1) where
c  fpint(number) is maximal on the condition that nrdata(number)
c  not equals zero.
      FPMAX = 0.0d0
      JBEGIN = ISTART
      DO 20 J=1,NRINT
        JPOINT = NRDATA(J)
        IF(FPMAX.GE.FPINT(J) .OR. JPOINT.EQ.0) GO TO 10
        FPMAX = FPINT(J)
        NUMBER = J
        MAXPT = JPOINT
        MAXBEG = JBEGIN
  10    JBEGIN = JBEGIN+JPOINT+1
  20  CONTINUE
c  let coincide the new knot t(number+k+1) with a data point x(nrx)
c  inside the old knot interval t(number+k) <= x <= t(number+k+1).
      IHALF = MAXPT/2+1
      NRX = MAXBEG+IHALF
      NEXT = NUMBER+1
      IF(NEXT.GT.NRINT) GO TO 40
c  adjust the different parameters.
      DO 30 J=NEXT,NRINT
        JJ = NEXT+NRINT-J
        FPINT(JJ+1) = FPINT(JJ)
        NRDATA(JJ+1) = NRDATA(JJ)
        JK = JJ+K
        T(JK+1) = T(JK)
  30  CONTINUE
  40  NRDATA(NUMBER) = IHALF-1
      NRDATA(NEXT) = MAXPT-IHALF
      AM = MAXPT
      AN = NRDATA(NUMBER)
      FPINT(NUMBER) = FPMAX*AN/AM
      AN = NRDATA(NEXT)
      FPINT(NEXT) = FPMAX*AN/AM
      JK = NEXT+K
      T(JK) = X(NRX)
      N = N+1
      NRINT = NRINT+1
      RETURN
      END

