      SUBROUTINE INTERV(XT,LXT,X,LEFT,MFLAG)
C#### Subroutine: INTERV
C###  Description:
C###     Auxililary routine for BVALUE.

c COMPUTES  LEFT = MAX( I :  XT(I) .LT. XT(LXT) .AND.  XT(I) .LE. X )  .
c
C******  I N P U T  ******
c  xt.....a real sequence, of length  lxt , assumed to be nondecreasing
c  lxt.....number of terms in the sequence  xt .
c  x.....the point whose location with respect to the sequence  xt  is
c        to be determined.
c
C******  O U T P U T  ******
c  left, mflag.....both integers, whose value is
c
c   1     -1      if               x .lt.  xt(1)
c   i      0      if   xt(i)  .le. x .lt. xt(i+1)
c   i      0      if   xt(i)  .lt. x .eq. xt(i+1) .eq. xt(lxt)
c   i      1      if   xt(i)  .lt.        xt(i+1) .eq. xt(lxt) .lt. x
c
c        In particular,  mflag = 0  is the 'usual' case.  mflag .ne. 0
c        indicates that  x  lies outside the CLOSED interval
c        xt(1) .le. y .le. xt(lxt) . The asymmetric treatment of the
c        intervals is due to the decision to make all pp functions cont-
c        inuous from the right, but, by returning  mflag = 0  even if
C        X = XT(LXT), THERE IS THE OPTION OF HAVING THE COMPUTED PP FUNCTION
c        continuous from the left at  xt(lxt) .
c
C******  M E T H O D  ******
c  The program is designed to be efficient in the common situation that
c  it is called repeatedly, with  x  taken from an increasing or decrea-
c  sing sequence. This will happen, e.g., when a pp function is to be
c  graphed. The first guess for  left  is therefore taken to be the val-
c  ue returned at the previous call and stored in the  l o c a l  varia-
c  ble  ilo . A first check ascertains that  ilo .lt. lxt (this is nec-
c  essary since the present call may have nothing to do with the previ-
c  ous call). Then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
c  ilo  and are done after just three comparisons.
c     Otherwise, we repeatedly double the difference  istep = ihi - ilo
c  while also moving  ilo  and  ihi  in the direction of  x , until
c                      xt(ilo) .le. x .lt. xt(ihi) ,
c  after which we use bisection to get, in addition, ilo+1 = ihi .
c  left = ilo  is then returned.
c
      IMPLICIT NONE
      INTEGER LEFT,LXT,MFLAG,IHI,ILO,ISTEP,MIDDLE
      REAL*8 X,XT(LXT)
      DATA ILO /1/
      SAVE ILO

      IHI=ILO+1
      IF(IHI.LT.LXT) GO TO 20
         IF(X.GE.XT(LXT)) GO TO 110
         IF(LXT.LE.1) GO TO 90
         ILO=LXT-1
         IHI=LXT
c
   20 IF(X.GE.XT(IHI)) GO TO 40
      IF(X.GE.XT(ILO)) GO TO 100
c
c **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
      ISTEP=1
   31    IHI=ILO
         ILO=IHI-ISTEP
         IF(ILO.LE.1) GO TO 35
         IF(X.GE.XT(ILO)) GO TO 50
         ISTEP=ISTEP*2
         GO TO 31
   35 ILO=1
      IF(X.LT.XT(1)) GO TO 90
      GO TO 50
c **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
   40 ISTEP=1
   41    ILO=IHI
         IHI=ILO+ISTEP
         IF(IHI.GE.LXT) GO TO 45
         IF(X.LT.XT(IHI)) GO TO 50
         ISTEP=ISTEP*2
         GO TO 41
   45 IF(X.GE.XT(LXT)) GO TO 110
      IHI=LXT
c
c **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
   50 MIDDLE=(ILO+IHI)/2
      IF(MIDDLE.EQ.ILO) GO TO 100
c     note. it is assumed that middle = ilo in case ihi = ilo+1 .
      IF(X.LT.XT(MIDDLE)) GO TO 53
         ILO=MIDDLE
         GO TO 50
   53    IHI = MIDDLE
         GO TO 50
C**** SET OUTPUT AND RETURN.
   90 MFLAG=-1
      LEFT=1
      RETURN
  100 MFLAG=0
      LEFT=ILO
      RETURN
  110 MFLAG=1
      IF(X.EQ.XT(LXT)) MFLAG=0
      LEFT=LXT
  111 IF(LEFT.EQ.1) RETURN
      LEFT=LEFT-1
      IF(XT(LEFT).LT.XT(LXT)) RETURN
      GO TO 111

C LKC 8-DEC-2000 FTNCHK says there is no path to this statement
C      RETURN

      END

