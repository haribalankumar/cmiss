      SUBROUTINE ANS001(D2SUMDXDY,DSUMDX,DSUMDY,SUM,TIME,X,Y,ERROR,*)

C#### Subroutine: ANS001
C###  Description:
C###    ANS001 calculates the analytical solution for diffusion of
C###    temperature or concentations of some particle over a
C###    rectangular plate where the boundary conditions on three sides
C###    are zero and the remaining side has a half sine wave profile.
C###    This analytical solution does not solve for advection.

      IMPLICIT NONE
      INCLUDE 'anal00.cmn'
      INCLUDE 'b00.cmn'
C      INCLUDE 'cmiss$reference:cbdi02.cmn'
C      INCLUDE 'cmiss$reference:cspi00.cmn'
C      INCLUDE 'cmiss$reference:geom00.cmn'
C      INCLUDE 'cmiss$reference:ityp00.cmn'

      !Passed variables
C      INTEGER
      REAL*8 D2SUMDXDY,DSUMDX,DSUMDY,SUM,TIME,X,Y
      CHARACTER ERROR*(*)
C      LOGICAL BRANCH
      !Local variables
      INTEGER i,N
      REAL*8 AN,D2SSDXDY,DSSDX,DSSDY,STEADY_STATE,TN
C      LOGICAL


      CALL ENTERS('ANS001',*9999)

      SUM=0.0d0 !value
      DSUMDX=0.0d0 !derivative wrt x
      DSUMDY=0.0d0 !derivative wrt y
      D2SUMDXDY=0.0d0 !cross derivative

      !loop over number of terms in solution
      N=100
      DO i=1,N
        AN=(2.0d0*i*ANAL_B*ANAL_B*(-1.0d0)**i)/
     '   (PI*(ANAL_A*ANAL_A+i*i*ANAL_B*ANAL_B))
        TN=-ANAL_D*PI*PI*(i*i/(ANAL_A*ANAL_A)+1.0d0
     '   /(ANAL_B*ANAL_B))*TIME
        SUM=SUM+AN*sin(i*PI*X/ANAL_A)*sin(PI*Y/ANAL_B)*EXP(TN)
        DSUMDX=DSUMDX+AN*(i*PI/ANAL_A)*cos(i*PI*X/ANAL_A)
     '    *sin(PI*Y/ANAL_B)*EXP(TN)
        DSUMDY=DSUMDY+AN*sin(i*PI*X/ANAL_A)*(PI/ANAL_B)
     '    *cos(PI*Y/ANAL_B)*EXP(TN)
        D2SUMDXDY=D2SUMDXDY+AN*(i*PI/ANAL_A)*cos(i*PI*X/ANAL_A)
     '    *(PI/ANAL_B)*cos(PI*Y/ANAL_B)*EXP(TN)
      ENDDO

      STEADY_STATE=SINH(PI*X/ANAL_B)*SIN(PI*Y/ANAL_B)/
     ' SINH(PI*ANAL_A/ANAL_B)
      DSSDX=(PI/ANAL_B)*COSH(PI*X/ANAL_B)*SIN(PI*Y/ANAL_B)/
     ' SINH(PI*ANAL_A/ANAL_B)
      DSSDY=SINH(PI*X/ANAL_B)*(PI/ANAL_B)*COS(PI*Y/ANAL_B)/
     ' SINH(PI*ANAL_A/ANAL_B)
      D2SSDXDY=(PI/ANAL_B)*COSH(PI*X/ANAL_B)*(PI/ANAL_B)
     ' *COS(PI*Y/ANAL_B)/SINH(PI*ANAL_A/ANAL_B)

      SUM=SUM+STEADY_STATE
      DSUMDX=DSUMDX+DSSDX
      DSUMDY=DSUMDY+DSSDY
      D2SUMDXDY=D2SUMDXDY+D2SSDXDY

      CALL EXITS('ANS001')
      RETURN
 9999 CALL ERRORS('ANS001',ERROR)
      CALL EXITS('ANS001')
      RETURN 1
      END

C time integration advection diffusion analytical solution
C DMAL 24-MAY-2002
