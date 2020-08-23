      REAL*8 FUNCTION RFROMC(CDATA)

C#### Function: RFROMC
C###  Type: REAL*8
C###  Description:
C###    RFROMC converts character string CDATA to a REAL*8 variable.

      IMPLICIT NONE
!     Parameter List
      CHARACTER CDATA*(*)
!     Local Variables
      REAL*8 RDATA

C     CALL ENTERS('RFROMC',*9999)
C!!! KAT 14/3/00: Should we check that CDATA is not blank?
      READ(UNIT=CDATA,FMT=*) RDATA
      RFROMC=RDATA

C     CALL EXITS('RFROMC')
      RETURN
      END

