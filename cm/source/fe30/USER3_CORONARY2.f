      SUBROUTINE USER3_CORONARY2(DENOM,
     '  NUM,PDERIV1,PDERIV2,VALUE,X,Y,ERROR,*)

C#### Subroutine: USER3_CORONARY2
C###  Description:
C###   USER3_CORONARY2  returns the value and derivatives of the
C###    function defining the lumped paramter models
C###    of coronary capillary resistance

      IMPLICIT NONE
!     Parameter List
      REAL*8 DENOM(6),DFDX,DGDX,DFDY,DGDY
      REAL*8 F,G,NUM(6),PDERIV1,PDERIV2,VALUE
      REAL*8 X,Y
      CHARACTER ERROR*(*)

      CALL ENTERS('USER3_CORONARY2',*9999)

      F=NUM(1)+(NUM(2)*X)+(NUM(3)*Y)+(NUM(4)*X*Y)+
     '  (NUM(5)*(X**2.d0))+(NUM(6)*(Y**2.d0))

      G=DENOM(1)+(DENOM(2)*X)+(DENOM(3)*Y)+(DENOM(4)*X*Y)+
     '  (DENOM(5)*(X**2.d0))+(DENOM(6)*(Y**2.d0))

      DFDX=(NUM(2)+(NUM(4)*Y)+(2.0D0*NUM(5)*X))

      DGDX=(DENOM(2)+(DENOM(4)*Y)+(2.0D0*DENOM(5)*X))

      VALUE=F/G !RC

      PDERIV1=((DFDX*G)-(DGDX*F))/(G**2.0D0) ! dRc/dP1

      DFDY=(NUM(3)+(NUM(4)*X)+(2.0D0*NUM(6)*Y))

      DGDY=(DENOM(3)+(DENOM(4)*X)+(2.0D0*DENOM(6)*Y))

      PDERIV2=((DFDY*G)-(DGDY*F))/(G**2.0D0) ! dRc/dP2

      CALL EXITS('USER3_CORONARY2')
      RETURN
 9999 CALL ERRORS('USER3_CORONARY2',ERROR)
      CALL EXITS('USER3_CORONARY2')
      RETURN 1
      END



