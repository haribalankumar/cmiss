      SUBROUTINE LIST_LABELS(STRING,ERROR,*)

C#### Subroutine: LIST_LABELS
C###  Description:
C###    Lists data objects.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'

!     Parameter List
c      INTEGER 
c      REAL*8 
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER IBEG,IEND,N3CO
      CHARACTER OBJECT*16
      LOGICAL CBBREV

      CALL ENTERS('LIST_LABELS',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list labels 
C###  Description:
C###    Lists the labels and parameter labels for equations, maths,
C###    models and fields.
C###  Parameter:    <(all/equation/maths/model/fields)[all]>
C###    Specifies which object labels to list.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<(all/equations/maths/models/'
     &    //'fields)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIST_LABELS',ERROR,*9999)
      ELSE
        
        IF(CBBREV(CO,'ALL',3,noco+1,NTCO,N3CO)) THEN
          OBJECT='ALL'
        ELSEIF(CBBREV(CO,'CONSTANTS',4,noco+1,NTCO,N3CO)) THEN
          OBJECT='CONSTANTS'
        ELSEIF(CBBREV(CO,'EQUATIONS',3,noco+1,NTCO,N3CO)) THEN
          OBJECT='EQUATIONS'
        ELSEIF(CBBREV(CO,'MATHS',4,noco+1,NTCO,N3CO)) THEN
          OBJECT='MATHS'
        ELSEIF(CBBREV(CO,'MODELS',5,noco+1,NTCO,N3CO)) THEN
          OBJECT='MODELS'
        ELSEIF(CBBREV(CO,'FIELDS',5,noco+1,NTCO,N3CO)) THEN
          OBJECT='FIELDS'
        ELSE
          OBJECT='ALL'
        ENDIF

        IF((OBJECT.EQ.'ALL').OR.(OBJECT.EQ.'CONSTANTS'))THEN
          CALL  PRINT_CONSTANTS(ERROR,*9999)       
        ENDIF
           
        IF((OBJECT.EQ.'ALL').OR.(OBJECT.EQ.'FIELDS'))THEN
          CALL PRINT_FIELDS(ERROR,*9999)          
        ENDIF
           
        IF((OBJECT.EQ.'ALL').OR.(OBJECT.EQ.'MATHS'))THEN
          CALL PRINT_MATHS(ERROR,*9999)
        ENDIF

        IF((OBJECT.EQ.'ALL').OR.(OBJECT.EQ.'EQUATIONS'))THEN
          CALL PRINT_EQUATIONS(ERROR,*9999)
        ENDIF
           
      ENDIF

      CALL EXITS('LIST_LABELS')
      RETURN
 9999 CALL ERRORS('LIST_LABELS',ERROR)
      CALL EXITS('LIST_LABELS')
      RETURN 1
      END


