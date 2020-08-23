      CHARACTER*(*) FUNCTION CFROMI(IDATA,FORMAT)

C#### Function: CFROMI
C###  Type: CHARACTER
C###  Description:
C###    Converts integer variable IDATA to character string CFROMI
C###    as specified by the character string FORMAT.

      IMPLICIT NONE
!     Parameter List
      INTEGER IDATA
      CHARACTER FORMAT*(*)
!     Local Variables
      CHARACTER CDATA*500

      WRITE(UNIT=CDATA,FMT=FORMAT) IDATA
      CFROMI=CDATA

      RETURN
      END


