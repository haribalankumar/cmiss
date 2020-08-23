      SUBROUTINE DEFINE_MATHS(CELL_RCQS_NAMES,CELL_RCQS_VALUE,
     &  CELL_YQS_NAMES,STRING,ERROR,*)

C#### Subroutine: DEFINE_MATHS
C###  Description:
C###    DEFINE_MATHS reads in a maths file, like a cellml, and generates
C###    the code required for the maths to be used in cmiss

C Author: Duane Malcolm
C Created: 11 March 2004

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
!      INTEGER 
      REAL*8 CELL_RCQS_VALUE(NQRM,NQVM)
      CHARACTER CELL_RCQS_NAMES(NQRM,NQVM)*(CELL_NAME_LENGTH),
     &  CELL_YQS_NAMES(NIQSM,NQVM)*(CELL_NAME_LENGTH),
     &  ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER ERROR_CODE,IBEG,IEND,IFROMC,INDEX1,N3CO,NCONSTS,
     '  NDERIVED,NFCNS,NODE,variable,variant,VERBOSE
      REAL*8 INITIAL_VALUE
      CHARACTER URI*132,VARNAME*255
      LOGICAL CBBREV,DEBUG,INITIALISED,SAVE_FILES,NIQSMEXCEEDED,
     '  TERMINATE_CALLED


      CALL ENTERS('DEFINE_MATHS',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM define maths
C###  Parameter:      NAME
C###    Name to label the maths
C###  Parameter:      uri PATH
C###    Path to the maths file
C###  Parameter:      <verbose verbosity#>
C###    Specifies the verbosity level
C###  Parameter:      <DEBUG>
C###    Turns on debugging
C###  Parameter:      <KEEPFILES>
C###    Keeps the generated fortran and object files
C###  Description:
C###    This command reads in a maths file, like a cellml, and
C###    generates the code required for the maths to be used in cmiss

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'NAME'
        OP_STRING(3)=BLANK(1:15)//'uri PATH'
        OP_STRING(4)=BLANK(1:15)//'<verbose verbosity#>'
        OP_STRING(5)=BLANK(1:15)//'<debug>'
        OP_STRING(6)=BLANK(1:15)//'<keepfiles>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','DEFINE_MATHS',ERROR,*9999)
      ELSE

        IF(CBBREV(CO,'URI',3,noco+1,NTCO,N3CO)) THEN
          URI=CO(N3CO+1)(1:132)
        ELSE
          ERROR='>> No URI defined'
          GOTO 9999
        ENDIF

        IF(CBBREV(CO,'VERBOSE',3,noco+1,NTCO,N3CO)) THEN
          VERBOSE=IFROMC(CO(N3CO+1))
        ELSE
          VERBOSE=0
        ENDIF
        
        IF(CBBREV(CO,'DEBUG',3,noco+1,NTCO,N3CO)) THEN
          DEBUG=.TRUE.
        ELSE
          DEBUG=.FALSE.
        ENDIF

        IF(CBBREV(CO,'KEEPFILES',4,noco+1,NTCO,N3CO)) THEN
          SAVE_FILES=.TRUE.
        ELSE
          SAVE_FILES=.FALSE.
        ENDIF        
        
        NDERIVED = 0
        NIQSMEXCEEDED = .FALSE.
        
C *** Initialise the CellML processor
        TERMINATE_CALLED = .FALSE.
        INITIALISED = .FALSE.
        CALL CELLML_INITIALISE(ERROR_CODE)
        CALL ASSERT(ERROR_CODE.EQ.0,
     '    'Failed to initialise CellML processor',ERROR,*9999)
        INITIALISED = .TRUE.

        CELL_NUM_VARIANTS=CELL_NUM_VARIANTS+1
        variant=CELL_NUM_VARIANTS

        CALL ADD_MATHS(CO(noco+1),ERROR,*9999)
        CALL SET_MATHS_VARIANT(CO(noco+1),variant,ERROR,*9999)
        CALL SET_MATHS_VERBOSE(CO(noco+1),VERBOSE,ERROR,*9999)
        
        CALL ASSERT(CELL_NUM_VARIANTS.LE.CELLML_MAX_MODELS,
     '    'Increase CELLML_MAX_MODELS in cellml.cmn',ERROR,*9999)

        CALL ASSERT(CELL_NUM_VARIANTS.LE.NQVM,'>>Increase NQVM',ERROR,
     '    *9999)
      
        CELLML_URIS(variant) = URI

C ***   Parse the cellML file - if we're gonna have multiple models in
C       memory at the same time, we're most likely gonna need some kind
C       of handle on each one...which is the variant number ??
        CALL CELLML_PARSE(variant,ERROR_CODE)
        CALL ASSERT(ERROR_CODE.EQ.0,'The file failed to parse',ERROR,
     '    *9999)

C ***   And create the math object for the model
        CALL CELLML_CREATE_MATH(variant,ERROR_CODE)
        CALL ASSERT(ERROR_CODE.EQ.0,'Error creating math object',ERROR,
     '    *9999)

C *** *** Set up the ODE variables - these need to be put
C       first into the YQS array in order for them to be safely
C       integrated.
        CALL CELLML_GET_NUM_ODE(variant,NODE)
        CELL_STATE_OFFSET(VARIANT)=1
        CELL_DERIVED_OFFSET(VARIANT)=1
        
        CELL_CONTROL_OFFSET(VARIANT)=1
        CELL_MODEL_OFFSET(VARIANT)=1
        CELL_PARAMETERS_OFFSET(VARIANT)=1
        CELL_PROTOCOL_OFFSET(VARIANT)=1
        CELL_AII_OFFSET(VARIANT)=1
        CELL_AIO_OFFSET(VARIANT)=1
        CELL_ARI_OFFSET(VARIANT)=1
        CELL_ARO_OFFSET(VARIANT)=1
      
        IF(NODE.GT.0)THEN
          WRITE(OP_STRING,'('' >>WARNING: ODES are not implemented'')')
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
        ENDIF
        
C ***   Set up the function variables - if the user wants to
C       save them they need to go into the derived array, otherwise they
C       are local variables in the cellular RHS routine
        CALL CELLML_GET_NUM_FUNCTIONS(variant,NFCNS)

        DO variable=1,NFCNS
          ! Get the variable name
          CALL CELLML_GET_FUNCTION_VAR_NAME(variant,variable,
     '      VARNAME)
          CALL STRING_TRIM(VARNAME,IBEG,IEND)
          
          NDERIVED = NDERIVED+1
          
          IF(NDERIVED.GT.NIQSM) NIQSMEXCEEDED=.TRUE.
          
          IF(.NOT.NIQSMEXCEEDED)THEN
            
            CELL_YQS_NAMES(NDERIVED,variant) = VARNAME(IBEG:IEND)
            
            CALL CELLML_SET_VARIABLE_ARRAY(variant,VARNAME(IBEG:IEND),
     '        'DERIVED')
          
            CALL CELLML_SET_VARIABLE_ARRAY_INDEX(variant,
     '        VARNAME(IBEG:IEND),NDERIVED)
            
            CALL ADD_MATHS_OUTPUT(CO(noco+1),VARNAME(IBEG:IEND),
     &        ERROR,*9999)
            CALL SET_MATHS_OUTPUT_YQS_INDEX(CO(noco+1),
     &        VARNAME(IBEG:IEND),NDERIVED,ERROR,*9999)
          
          ENDIF

        ENDDO

C ***   Set up the constant variables - if the user wants to save them,
C       they need to be in the DERIVED array; otherwise they go in the
C       PARAM array. - Oops, if you put a parameter in the DERIVED array
C       you cannot specify spatial variation...
        CALL CELLML_GET_NUM_CONSTANTS(variant,NCONSTS)

        DO variable=1,NCONSTS
          INDEX1=variable
          ! Get the variable name
          CALL CELLML_GET_CONST_VAR_NAME(variant,variable,
     '      VARNAME)
          CALL STRING_TRIM(VARNAME,IBEG,IEND)
          CELL_RCQS_NAMES(INDEX1,variant) = VARNAME(IBEG:IEND)
          CALL CELLML_GET_CONST_VAR_INIT_VALUE(variant,
     '      variable,INITIAL_VALUE)
          CELL_RCQS_VALUE(INDEX1,variant) =INITIAL_VALUE
          CALL CELLML_SET_VARIABLE_ARRAY(variant,VARNAME(IBEG:IEND),
     '      'PARAM')
          CALL CELLML_SET_VARIABLE_ARRAY_INDEX(variant,
     '      VARNAME(IBEG:IEND),INDEX1)
          CALL ADD_MATHS_INPUT(CO(noco+1),VARNAME(IBEG:IEND),
     &        ERROR,*9999)
          CALL SET_MATHS_INPUT_INITIAL_VALUE(CO(noco+1),
     &      VARNAME(IBEG:IEND),INITIAL_VALUE,ERROR,*9999)
          CALL SET_MATHS_INPUT_RCQS_INDEX(CO(noco+1),VARNAME(IBEG:IEND),
     &      INDEX1,ERROR,*9999)
          CALL SET_MATHS_INPUT_MAP(CO(noco+1),VARNAME(IBEG:IEND),
     &      'INITIAL_VALUE',' ',' ',ERROR,*9999)
          
          
        ENDDO
        
C ***   Now have enough information to create the RHS routine
        CALL CELLML_CREATE_RHS_ROUTINE(variant,DEBUG,SAVE_FILES,
     '    ERROR_CODE)
        CALL ASSERT(ERROR_CODE.EQ.0,
     '    'Failed to create RHS routine',ERROR,*9999)

C *** Terminate the CellML processor
        CALL CELLML_TERMINATE(ERROR_CODE)
        TERMINATE_CALLED=.TRUE.
        CALL ASSERT(ERROR_CODE.EQ.0,
     '    'Failed to terminate CellML processor',ERROR,*9999)

        IF(NDERIVED.GT.NIQSM) THEN
          WRITE(ERROR,'(''>>Increase NIQSM to >= '',I12)') NDERIVED
          GOTO 9999
        ENDIF

      ENDIF

      CALL EXITS('DEFINE_MATHS')
      RETURN
 9999 CALL ERRORS('DEFINE_MATHS',ERROR)

C *** Terminate the CellML processor if it was initialised and not
C     terminated.
      IF (INITIALISED.AND..NOT.TERMINATE_CALLED) THEN
        CALL CELLML_TERMINATE(ERROR_CODE)
        TERMINATE_CALLED=.TRUE.
        CALL ASSERT(ERROR_CODE.EQ.0,
     '    'Failed to terminate CellML processor',ERROR,*9999)
      ENDIF

      CALL EXITS('DEFINE_MATHS')
      RETURN 1
      END


