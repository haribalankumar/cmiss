      SUBROUTINE CABASE(STRING,ERROR,*)

C#### Subroutine: CABASE
C###  Description:
C###    CABASE cancels basis functions.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,nb

      CALL ENTERS('CABASE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel bases
C###  Description:
C###    Cancel all bases defined.

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CABASE',ERROR,*9999)
      ELSE
        DO nb=1,NBFT
          NBASEF(nb,0)=0
        ENDDO
        NBT=0
        NBFT=0
      ENDIF

      CALL EXITS('CABASE')
      RETURN
 9999 CALL ERRORS('CABASE',ERROR)
      CALL EXITS('CABASE')
      RETURN 1
      END


