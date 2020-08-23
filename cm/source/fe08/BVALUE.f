      REAL*8 FUNCTION BVALUE(T,BCOEF,N,K,X,JDERIV)
C#### Function: BVALUE
C###  Type: REAL*8
C###  Description:
C###    BVALUE calculates the j-th derivative of a B-spline. This
C###    routine and all subsequent calling routines are a subset of the
C###    PPPACK library. These routines appear exactly as the author
C###    wrote them except the precision has been changed from REAL to
C###    REAL*8. Below contains the comments from the author.
C###  Reference:
C###    C. De. Boor, "A pritical quide to splines"
c
c  CALCULATES VALUE AT  X  OF  JDERIV-TH DERIVATIVE OF SPLINE FROM B-REPR.
c  the spline is taken to be continuous from the right, EXCEPT at the
c  rightmost knot, where it is taken to be continuous from the left.
c
C******  I N P U T ******
c  t, bcoef, n, k......forms the b-representation of the spline  f  to
c        be evaluated. specifically,
c  t.....knot sequence, of length  n+k, assumed nondecreasing.
c  bcoef.....b-coefficient sequence, of length  n .
c  n.....length of  bcoef  and dimension of spline(k,t),
c        a s s u m e d  positive .
c  k.....order of the spline .
c
c  w a r n i n g . . .   the restriction  k .le. kmax (=20)  is imposed
c        arbitrarily by the dimension statement for  aj, dl, dr  below,
c        but is  n o w h e r e  c h e c k e d  for.
c
c  x.....the point at which to evaluate .
c  jderiv.....integer giving the order of the derivative to be evaluated
c        a s s u m e d  to be zero or positive.
c
C******  O U T P U T  ******
c  bvalue.....the value of the (jderiv)-th derivative of  f  at  x .
c
C******  M E T H O D  ******
c     The nontrivial knot interval  (t(i),t(i+1))  containing  x  is lo-
c  cated with the aid of  interv . The  k  b-coeffs of  f  relevant for
c  this interval are then obtained from  bcoef (or taken to be zero if
c  not explicitly available) and are then differenced  jderiv  times to
c  obtain the b-coeffs of  (d**jderiv)f  relevant for that interval.
c  Precisely, with  j = jderiv, we have from x.(12) of the text that
c
c     (d**j)f  =  sum ( bcoef(.,j)*b(.,k-j,t) )
c
c  where
c                   / bcoef(.),                     ,  j .eq. 0
c                   /
c    bcoef(.,j)  =  / bcoef(.,j-1) - bcoef(.-1,j-1)
c                   / ----------------------------- ,  j .gt. 0
c                   /    (t(.+k-j) - t(.))/(k-j)
c
c     Then, we use repeatedly the fact that
c
c    sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) )
c  with
c                 (x - t(.))*a(.) + (t(.+m-1) - x)*a(.-1)
c    a(.,x)  =    ---------------------------------------
c                 (x - t(.))      + (t(.+m-1) - x)
c
c  to write  (d**j)f(x)  eventually as a linear combination of b-splines
c  of order  1 , and the coefficient for  b(i,1,t)(x)  must then be the
c  desired number  (d**j)f(x). (see x.(17)-(19) of text).
c
      IMPLICIT NONE
      INTEGER JDERIV,K,N,I,ILO,IMK,J,JC,JCMIN,JCMAX,JJ,KMAX,KMJ,KM1,
     '  MFLAG,NMI,JDRVP1
      PARAMETER (KMAX=20)
      REAL*8 BCOEF(N),T(N+K),X,AJ(KMAX),DDL(KMAX),DDR(KMAX),FKMJ

      BVALUE=0.0d0
      IF(JDERIV.GE.K) GO TO 99
c
c *** Find  i   s.t.   1 .le. i .lt. n+k   and   t(i) .lt. t(i+1)   and
c     t(i) .le. x .lt. t(i+1) . If no such i can be found,  x  lies
c     outside the support of  the spline  f , hence  bvalue = 0.
c     (The asymmetry in this choice of  i  makes  f  rightcontinuous, except
c     at  t(n+k) where it is leftcontinuous.)
      CALL INTERV(T,N+K,X,I,MFLAG)
      IF(MFLAG.NE.0) GO TO 99
c *** if k = 1 (and jderiv = 0), bvalue = bcoef(i).
      KM1=K-1
      IF(KM1.GT.0) GO TO 1
      BVALUE=BCOEF(I)
      GO TO 99
c
c *** store the k b-spline coefficients relevant for the knot interval
c    (t(i),t(i+1)) in aj(1),...,aj(k) and compute ddl(j) = x - t(i+1-j),
c    ddr(j) = t(i+j) - x, j=1,...,k-1 . set any of the aj not obtainable
c    from input to zero. set any t.s not obtainable equal to t(1) or
c    to t(n+k) appropriately.
    1 JCMIN=1
      IMK=I-K
      IF(IMK.GE.0) GO TO 8
      JCMIN=1-IMK
      DO 5 J=1,I
    5    DDL(J)=X-T(I+1-J)
      DO 6 J=I,KM1
         AJ(K-J)=0.0d0
    6    DDL(J)=DDL(I)
      GO TO 10
    8 DO 9 J=1,KM1
    9    DDL(J)=X-T(I+1-J)
c
   10 JCMAX=K
      NMI=N-I
      IF (NMI.GE.0) GO TO 18
      JCMAX=K+NMI
      DO 15 J=1,JCMAX
   15    DDR(J)=T(I+J)-X
      DO 16 J=JCMAX,KM1
         AJ(J+1)=0.0d0
   16    DDR(J)=DDR(JCMAX)
      GO TO 20
   18 DO 19 J=1,KM1
   19    DDR(J)=T(I+J)-X
c
   20 DO 21 JC=JCMIN,JCMAX
   21    AJ(JC)=BCOEF(IMK+JC)
c
c *** difference the coefficients  jderiv  times.
      IF (JDERIV.EQ.0) GO TO 30
      DO 23 J=1,JDERIV
         KMJ=K-J
         FKMJ=DBLE(KMJ)
         ILO=KMJ
         DO 23 JJ=1,KMJ
            AJ(JJ)=((AJ(JJ+1)-AJ(JJ))/(DDL(ILO)+DDR(JJ)))*FKMJ
   23       ILO=ILO-1
c
c *** compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative,
c     given its relevant b-spline coeffs in aj(1),...,aj(k-jderiv).
   30 IF (JDERIV.EQ.KM1) GO TO 39
      JDRVP1=JDERIV+1
      DO 33 J=JDRVP1,KM1
         KMJ=K-J
         ILO=KMJ
         DO 33 JJ=1,KMJ
            AJ(JJ)=(AJ(JJ+1)*DDL(ILO)+AJ(JJ)*DDR(JJ))/(DDL(ILO)+DDR(JJ))
   33       ILO=ILO-1
   39 BVALUE=AJ(1)
c
   99 RETURN
      END


