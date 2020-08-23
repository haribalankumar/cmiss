      SUBROUTINE LINREGRESS(N,R_SQUARED,slope,X,Y,ERROR,*)

C#### Subroutine: LINREGRESS
C###  Description:
C###    Calculates linear regression equation and r-squared 
C###    correlation coefficient for a set of data.

C*** Created by Kelly Burrowes, March 2003.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter list
      INTEGER N
      REAL*8 R_SQUARED,slope,X(*),Y(*)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER i
      REAL*8 AX,AY,intercept,R,SXX,SXY,SYY,XSUM,XT,XXSUM,XYSUM,YSUM,YT

      CALL ENTERS('LINREGRESS',*9999)
      
      YSUM=0.d0
      XSUM=0.d0
      XXSUM=0.d0
      XYSUM=0.d0
      DO i=1,N
        YSUM=YSUM+Y(i)
        XSUM=XSUM+X(i)
        XYSUM=XYSUM+X(i)*Y(i)
        XXSUM=XXSUM+X(i)*X(i)
      ENDDO !N
C... calculate least squares estimate of straight line thru solution
      slope=(XYSUM-XSUM*YSUM/N)/(XXSUM-XSUM*XSUM/N) 
      intercept=YSUM/N-slope*XSUM/N 
C... calculate r-squared correlation coefficient
c... see Numerical Recipes, Fortran 77, 2nd edition, page 632.      
      AX=XSUM/N !mean of X
      AY=YSUM/N !mean of Y
      SXX=0.d0
      SYY=0.d0
      SXY=0.d0
      DO i=1,N
        XT=X(i)-AX
        YT=Y(i)-AY
        SXX=SXX+XT**2.d0
        SYY=SYY+YT**2.d0
        SXY=SXY+XT*YT
      ENDDO
      R=SXY/DSQRT(SXX*SYY)
      R_SQUARED=R**2
      IF(DOP)THEN
        WRITE(OP_STRING,'('' Gradient: '',F8.3,'' Intercept: '','
     '    //' F8.3,''R-squared correlation'',F8.3)')
     '    slope,intercept,R_SQUARED
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      
      
      CALL EXITS('LINREGRESS')
      RETURN
 9999 CALL ERRORS('LINREGRESS',ERROR)
      CALL EXITS('LINREGRESS')
      RETURN 1
      END

      
