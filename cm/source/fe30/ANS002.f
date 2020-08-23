      SUBROUTINE ANS002(DSUMDX,SUM,TIME,X,ERROR,*)

C#### Subroutine: ANS002
C###  Description:
C###    ANS002 calculates the analytical solution for diffusion of
C###    temperature or concentations of some particle in a
C###    rod where the boundary conditions at the two end
C###    are fixed at T1 and T2. The initial conditions is a half sine
C###    over over the length, L of the rod.
C###    This analytical solution does not solve for advection.

      IMPLICIT NONE
      INCLUDE 'anal00.cmn'
      INCLUDE 'b00.cmn'
C      INCLUDE 'cmiss$reference:cbdi02.cmn'
C      INCLUDE 'cmiss$reference:cspi00.cmn'
C      INCLUDE 'cmiss$reference:geom00.cmn'
C      INCLUDE 'cmiss$reference:ityp00.cmn'

      !Passed variables
      REAL*8 DSUMDX,SUM,TIME,X
      CHARACTER ERROR*(*)
C      LOGICAL BRANCH
      !Local variables
      INTEGER i,N
      REAL*8 AN,DSSDX,STEADY_STATE,TN
C      REAL*8 CHANGE,SUM_OLD,TOL


      CALL ENTERS('ANS002',*9999)

      SUM=0.0d0 !value
      DSUMDX=0.0d0 !derivative wrt x

      !loop over number of terms in solution
C      TOL=1.0d-12
C      CHANGE=1.0d0
      N=2000
C      DO WHILE(CHANGE>TOL)
      DO i=1,N
        AN=-(ANAL_A*(ANAL_C-ANAL_B)*SIN(i*PI))/((i*PI)**2) +
     '    (ANAL_A*(ANAL_C-ANAL_B)*COS(i*PI))/(i*PI) +
     '    (ANAL_A*ANAL_B*COS(i*PI))/(i*PI) -
     '    (ANAL_A*ANAL_B)/(i*PI)-COS(i*PI)*ANAL_A/(i*PI) +
     '    (ANAL_A/(i*PI))
        TN=-(i**2*PI**2*ANAL_D**2)/(ANAL_A**2)
C        SUM_OLD=SUM
        SUM=SUM+AN*SIN(i*PI*X/ANAL_A)* EXP(TN*TIME)
C        CHANGE=ABS(SUM_OLD-SUM)
        DSUMDX=DSUMDX+((AN*ANAL_A)/(i*PI))*SIN(i*PI*X/ANAL_A)*
     '    EXP(TN*TIME)
C WRITE(*,'('' i='',I6,'' CHANGE='',D12.6)')
C     '                  i,CHANGE
      ENDDO

      SUM=(2.0d0/ANAL_A)*SUM
      DSUMDX=(2.0d0/ANAL_A)*DSUMDX
      STEADY_STATE=(1.0d0/ANAL_A)*(ANAL_C-ANAL_B)*X+ANAL_B
      DSSDX=(1.0d0/ANAL_A)*(ANAL_C-ANAL_B)
      SUM=SUM+STEADY_STATE
      DSUMDX=DSUMDX+DSSDX


      CALL EXITS('ANS002')
      RETURN
 9999 CALL ERRORS('ANS002',ERROR)
      CALL EXITS('ANS002')
      RETURN 1
      END


