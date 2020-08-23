      REAL*8 FUNCTION RANDSEEDS(i,j,k)
C#### Function: RANDSEEDS
C###  Type: REAL*8
C###  Description:
C###    Performs some numerical manipulation of 3 seeds to produce a
C###    random number.
C###
C###    NOTE : Input Parameters i,j,k refer to the 3 seeds used
C###           generate the random number

      IMPLICIT NONE
!     Paramter List
      INTEGER i,j,k
!     Local Variables
      INTEGER i_lower
      REAL*8  temp

      i=171*(MOD(i,177))-2*(i/177)
      j=172*(MOD(j,176))-35*(j/176)
      k=170*(MOD(k,178))-63*(k/178)

      IF(i.LT.0.D0) i=i+30269
      IF(j.LT.0.D0) j=j+30307
      IF(k.LT.0.D0) k=k+30323

      temp=DBLE(i)/30269.0D0+ DBLE(j)/30307.0D0+ DBLE(k)/30323.0D0

C***  Below follows a sequence to imitate the FLOOR function
      IF(temp.GT.0.D0) THEN
        i_lower=INT(temp)
      ELSE
        i_lower=INT(temp)-1
      ENDIF
      RANDSEEDS=temp-DBLE(i_lower)
      END


