      SUBROUTINE CMDESTROY(ERR)

C#### Subroutine: CMDESTROY
C###  Description:
C###    Destroys the initial global variables for CMISS.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cmgui00.cmn'
      INCLUDE 'gtstr00.cmn'
      INCLUDE 'mxch.inc'
!     Parameter List
      INTEGER ERR
!     Local Variables
      CHARACTER ERROR*(MXCH)

C     Destroy the cm-cmgui link
      IF(CMGUI_LINK) CALL CMGUI_LINK_DESTROY(ERROR,*9999)
C KAT 25/2/00: closing open comfile
      IF(IREC_COMFILE(0).NE.-1) CALL CLOSEF(OPEN_COM_UNIT,ERROR,*9999)

      ERR=0
      RETURN
 9999 CALL ERRORS('CMDESTROY',ERROR)
      CALL ERRORIN(' ') ! end of fortran cm call stack
      ERR=1
      RETURN
      END


