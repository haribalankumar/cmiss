      SUBROUTINE EVALUATE_MATHS(MATHS_LABEL,XG,DXIX,ERROR,*)

C#### Subroutine: EVALUATE_MATHS
C###  Description:
C###    Evaluates a maths object and updates the 
C###    maths output variable values.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'constant00.cmn'
      INCLUDE 'field00.cmn'
      INCLUDE 'equation00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'maths00.cmn'
      
!     Parameter List
      REAL*8 DXIX(3,3),XG(NJM,NUM)
      CHARACTER MATHS_LABEL*(*),ERROR*(*)
!     Local Variables
      INTEGER IB1,IB2,IE1,IE2,INPUT_TOTAL,index,nj,nu,OUTPUT_TOTAL,
     &  RCQS_INDEX,SOLVED,VARIANT,VERBOSE,YQS_INDEX
      INTEGER ICQS(NQIM)
      REAL*8 RCQS(NQRM),YQS(NIQSM),VALUE
      CHARACTER MAP_TYPE*(MATHS_IN_MAP_TYPE_LEN),
     &  MAP_PARENT*(MATHS_IN_MAP_PARENT_LEN),
     &  MAP_CHILD*(MATHS_IN_MAP_CHILD_LEN),
     &  MATHS_OUT_LABEL*(MATHS_OUT_LABELS_LEN),
     &  MATHS_IN_LABEL*(MATHS_IN_LABELS_LEN)
      
      CALL ENTERS('EVALUATE_MATHS',*9999)
      
      CALL GET_MATHS_VARIANT(MATHS_LABEL,VARIANT,ERROR,*9999)
      
      CALL GET_MATHS_INPUT_TOTAL(MATHS_LABEL,INPUT_TOTAL,
     &  ERROR,*9999)

      DO index=1,INPUT_TOTAL
        CALL GET_MATHS_INPUT_LABEL(MATHS_LABEL,index,
     &    MATHS_IN_LABEL,ERROR,*9999)
        CALL GET_MATHS_INPUT_MAP(MATHS_LABEL,
     &      MATHS_IN_LABEL,MAP_TYPE,MAP_PARENT,MAP_CHILD,
     &      ERROR,*9999)
     
        IF(MAP_TYPE.EQ.'INITIAL_VALUE')THEN
          CALL GET_MATHS_INPUT_INITIAL_VALUE(MATHS_LABEL,
     &      MATHS_IN_LABEL,VALUE,ERROR,*9999)
        ELSEIF(MAP_TYPE.EQ.'CONSTANT')THEN
          CALL GET_CONSTANT_COMPONENT_VALUE(MAP_PARENT,MAP_CHILD,
     &      VALUE,ERROR,*9999)
        ELSEIF(MAP_TYPE.EQ.'FIELD')THEN
          CALL GET_FIELD_NJ(MAP_PARENT,nj,ERROR,*9999)
          CALL GET_FIELD_COMPONENT_NU(MAP_PARENT,MAP_CHILD,
     &      nu,ERROR,*9999)

C DMAL 13 APRIL 2004: This code is dodgy and hasn't been checked
C for non-unit scale factors
          IF(nu.EQ.1)THEN
            VALUE=XG(nj,nu)
          ELSEIF(nu.EQ.2)THEN
            VALUE=XG(nj,nu)*DXIX(1,1)
          ELSEIF(nu.EQ.3)THEN
            VALUE=XG(nj,nu)*DXIX(1,1)*DXIX(1,1)
          ENDIF
          CALL STRING_TRIM(MAP_PARENT,IB1,IE1)
          CALL STRING_TRIM(MAP_CHILD,IB2,IE2)
C DMAL 13 APRIL 2004

        ELSEIF(MAP_TYPE.EQ.'MATHS')THEN
          ! Fortran 77 Quirk: a routine can't call itself, so I created 
          ! another routine called EVALUATE_MATHS_SELF to call EVALUATE_MATHS.
          CALL GET_MATHS_SOLVED(MAP_PARENT,SOLVED,ERROR,*9999)
          IF(SOLVED.EQ.0)THEN
            CALL EVALUATE_MATHS_SELF(MAP_PARENT,XG,DXIX,ERROR,*9999)
          ENDIF
            
          CALL GET_MATHS_OUTPUT_VALUE(MAP_PARENT,MAP_CHILD,VALUE,
     &      ERROR,*9999)
        ELSEIF(MAP_TYPE.EQ.'EQUATION')THEN
          CALL ASSERT(.FALSE.,'Not Implemented. '//
     &      'Due to the way cmiss is I can not get to solution '//
     &      'fields (YG) without a lot of work',ERROR,*9999)
        ELSE
          CALL ASSERT(.FALSE.,'Error: Unidentified MAP_TYPE',
     &      ERROR,*9999)
        ENDIF
              
        CALL GET_MATHS_INPUT_RCQS_INDEX(MATHS_LABEL,
     &    MATHS_IN_LABEL,RCQS_INDEX,ERROR,*9999)
        CALL SET_MATHS_INPUT_VALUE(MATHS_LABEL,
     &    MATHS_IN_LABEL,VALUE,ERROR,*9999)
        RCQS(RCQS_INDEX)=VALUE
        
      ENDDO
      

C call the RHSROUTINE through EVAL_MODEL_DYNAM
      CALL EVAL_MODEL_DYNAM(ICQS,RCQS,
     &  %VAL(CELLML_ROUTINES(VARIANT)),
     &  VARIANT,YQS,ERROR,*9999)

      CALL GET_MATHS_OUTPUT_TOTAL(MATHS_LABEL,OUTPUT_TOTAL,
     &  ERROR,*9999)

      DO index=1,OUTPUT_TOTAL
        CALL GET_MATHS_OUTPUT_LABEL(MATHS_LABEL,index,
     &    MATHS_OUT_LABEL,ERROR,*9999)
        CALL GET_MATHS_OUTPUT_YQS_INDEX(MATHS_LABEL,
     &    MATHS_OUT_LABEL,YQS_INDEX,ERROR,*9999)
        CALL SET_MATHS_OUTPUT_VALUE(MATHS_LABEL,
     &    MATHS_OUT_LABEL,YQS(YQS_INDEX),ERROR,*9999)
      ENDDO
      
      CALL SET_MATHS_SOLVED(MATHS_LABEL,1,ERROR,*9999)
      
C DMAL 26 MARCH 2004: Personally I think this block of code has no business here
C since this subroutine should do only one job, update the values in CG. When
C this is removed, remove the COMMENT variable above.        
      CALL GET_MATHS_VERBOSE(MATHS_LABEL,VERBOSE,ERROR,*9999)
      IF(VERBOSE.GT.0) THEN
        WRITE(OP_STRING,'()')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        CALL STRING_TRIM(MATHS_LABEL,IB1,IE1)
        WRITE(OP_STRING,'(''Maths '',A,'' Inputs:'')') 
     &    MATHS_LABEL(IB1:IE1)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        CALL PRINT_MATHS_INPUT_VALUES(MATHS_LABEL,ERROR,*9999)
        WRITE(OP_STRING,'(''Maths '',A,'' Outputs:'')') 
     &    MATHS_LABEL(IB1:IE1)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        CALL PRINT_MATHS_OUTPUT_VALUES(MATHS_LABEL,ERROR,*9999)
      ENDIF
      

      CALL EXITS('EVALUATE_MATHS')
      RETURN
 9999 CALL ERRORS('EVALUATE_MATHS',ERROR)
      CALL EXITS('EVALUATE_MATHS')
      RETURN 1
      END

