      INTEGER FUNCTION IFROMC(CDATA)

C#### Function: IFROMC
C###  Type: INTEGER
C###  Description:
C###    Converts character string CDATA to an integer.
C###    If CDATA does not look like an integer then a warning is
C###    produced and a semi-reasonable guess is returned.  It is
C###    probably preferrable to use INTFROMCHAR instead and check the
C###    error code.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER (ROUTINENAME='IFROMC')
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      CHARACTER CDATA*(*)
!     Local Variables
      INTEGER ERR,IRESULT

C     CALL ENTERS('IFROMC',*9999)

      CALL INTFROMCHAR(IRESULT,CDATA,ERR)

      IF(ERR.NE.0) THEN
        CALL ERRORIN(ROUTINENAME)
C       Flush the call stack as this error is being ignored
        CALL ERRORIN(' ')
      ENDIF

      IFROMC=IRESULT
C     CALL EXITS('IFROMC')
      RETURN
      END


