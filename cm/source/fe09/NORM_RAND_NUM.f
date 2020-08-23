      REAL*8 FUNCTION NORM_RAND_NUM(ISEED,MEAN,SD)

C#### Function: NORM_RAND_NUM
C###  Description:
C###    This function returns a normally distributed random
C###    number with mean MEAN and std dev SD.

      INTEGER i,ISEED
      REAL*8 MEAN,CM_RANDOM_NUMBER,SD,SUM,SUM2,SUMTEMP

      SUM=-6.0d0
      DO i=1,12
        SUM=SUM+CM_RANDOM_NUMBER(ISEED)
      ENDDO
      SUM2=SUM*SD
      SUMTEMP=MEAN+SUM2
      NORM_RAND_NUM=SUMTEMP

      RETURN
      END


C FE04 Functions
C ==============

