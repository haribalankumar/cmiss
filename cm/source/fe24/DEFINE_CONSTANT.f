      SUBROUTINE DEFINE_CONSTANT(STRING,ERROR,*)

C#### Subroutine: DEFINE_CONSTANT
C###  Description:
C###    Defines a constant and it's value.
C###    If the constant or its component does not exist then they
C###    would be created. If the constant or its component exists
C###    then the value would be updated.

C Author: Duane Malcolm
C Created: 11 March 2004

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'constant00.cmn'
      
!     Parameter List
      CHARACTER STRING*(MXCH),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,INDEX,N3CO
      REAL*8 RFROMC,VALUE
      CHARACTER COMPONENT_LABEL*(CONST_COMPNT_LABELS_LEN),
     &  CONSTANT_LABEL*(CONSTANT_LABELS_LEN)
      LOGICAL CBBREV

      CALL ENTERS('DEFINE_CONSTANT',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM define CONSTANT
C###  Parameter:      NAME
C###    Specifies the label for the constant
C###  Parameter:      <component NAME [default]>
C###    Specifies the label for the constant component
C###  Parameter:      <value VALUE [0.0]>
C###    Specifies the value of the constant component
C###  Description:
C###    This command defines constant and its value.
C###    If the constant or its component does not exist then they
C###    would be created. If the constant or its component exists
C###    then the value would be updated.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'NAME'
        OP_STRING(3)=BLANK(1:15)//'<component NAME [default]>'
        OP_STRING(4)=BLANK(1:15)//'<value VALUE [0.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','DEFINE_CONSTANT',ERROR,*9999)
      ELSE
        
        CONSTANT_LABEL=CO(noco+1)
        IF(CBBREV(CO,'COMPONENT',4,noco+1,NTCO,N3CO)) THEN
          COMPONENT_LABEL=CO(N3CO+1)
        ELSE
          COMPONENT_LABEL='default'
        ENDIF
          
        IF(CBBREV(CO,'VALUE',3,noco+1,NTCO,N3CO)) THEN
          VALUE=RFROMC(CO(N3CO+1))
        ELSE
          VALUE=0.0D0
        ENDIF
          
        CALL GET_CONSTANT_INDEX(CONSTANT_LABEL,INDEX,ERROR,*9999)
        IF(INDEX.EQ.0)THEN ! Constant doesn't exist
          CALL ADD_CONSTANT(CONSTANT_LABEL,ERROR,*9999)
        ENDIF

        CALL GET_CONSTANT_COMPONENT_INDEX(CONSTANT_LABEL,
     &    COMPONENT_LABEL,INDEX,ERROR,*9999)
        IF(INDEX.EQ.0)THEN ! Constant component doesn't exist
          CALL ADD_CONSTANT_COMPONENT(CONSTANT_LABEL,COMPONENT_LABEL,
     &      ERROR,*9999)
        ENDIF
        
        CALL SET_CONSTANT_COMPONENT_VALUE(CONSTANT_LABEL,
     &    COMPONENT_LABEL,VALUE,ERROR,*9999)
     
      ENDIF

      CALL EXITS('DEFINE_CONSTANT')
      RETURN
 9999 CALL ERRORS('DEFINE_CONSTANT',ERROR)
      CALL EXITS('DEFINE_CONSTANT')
      RETURN 1
      END


