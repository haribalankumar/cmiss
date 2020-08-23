      LOGICAL FUNCTION LFROMC(CDATA)

C#### Function: LFROMC
C###  Type: LOGICAL
C###  Description:
C###    LRFOMC converts character string CDATA to a logical variable.

      IMPLICIT NONE
!     Parameter List
      CHARACTER CDATA*(*)
!     Local Variables
      LOGICAL LDATA

C     CALL ENTERS('LFROMC',*9999)
      READ(UNIT=CDATA,FMT=*) LDATA
      LFROMC=LDATA

C     CALL EXITS('LFROMC')
      RETURN
      END


