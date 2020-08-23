      SUBROUTINE USER3_CORONARY1(DENOM,DERIV,
     '  DYDX,NUM,VALUE,X,Y,ERROR,*)

C#### Subroutine: USER3_CORONARY1
C###  Description:
C###    USER3_CORONARY1 returns the vlaue and derivatives of the
C###    function defining the lumped paramter models
C###    of coronary arteriol and venule resistance and capcitance

      IMPLICIT NONE
!     Parameter List
      REAL*8 DENOM(6),DERIV,DFDX,DGDX,DYDX,F,G,NUM(6),VALUE,X,Y
      CHARACTER ERROR*(*)

      CALL ENTERS('USER3_CORONARY1',*9999)

      F=NUM(1)+(NUM(2)*X)+(NUM(3)*Y)+(NUM(4)*X*Y)+
     '  (NUM(5)*(X**2.d0))+(NUM(6)*(Y**2.d0))

      G=DENOM(1)+(DENOM(2)*X)+(DENOM(3)*Y)+(DENOM(4)*X*Y)+
     '  (DENOM(5)*(X**2.d0))+(DENOM(6)*(Y**2.d0))

      DFDX=(NUM(2)+(NUM(4)*Y)+(2.0D0*NUM(5)*X))+(DYDX*
     '  (NUM(3)+(NUM(4)*X)+(2.0D0*NUM(6)*Y)))

      DGDX=(DENOM(2)+(DENOM(4)*Y)+(2.0D0*DENOM(5)*X))+(DYDX*
     '  (DENOM(3)+(DENOM(4)*X)+(2.0D0*DENOM(6)*Y)))

      VALUE=F/G

      DERIV=((DFDX*G)-(DGDX*F))/(G**2.0D0)

      CALL EXITS('USER3_CORONARY1')
      RETURN
 9999 CALL ERRORS('USER3_CORONARY1',ERROR)
      CALL EXITS('USER3_CORONARY1')
      RETURN 1
      END



