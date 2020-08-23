      SUBROUTINE UPDATE_EQUATION_INPUTS(nx,CG,XG,DXIX,COMMENT,ERROR,*)

C#### Subroutine: UPDATE_EQUATION_INPUTS
C###  Description:
C###    Updates the values in CG for an equation using the 
C###    mapping.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'equation00.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER nx
      REAL*8 CG(NMM),DXIX(3,3),XG(NJM,NUM)
      CHARACTER COMMENT*64,ERROR*(*)

!     Local Parameters
      INTEGER EQN_IN_TOTAL,EQUATION_TOTAL,IB1,IE1,index,
     &  lnx,nj,nm,nu,SOLVED,VERBOSE
      REAL*8 VALUE
      CHARACTER EQUATION_LABEL*(EQUATION_LABELS_LEN),
     &  EQN_IN_LABEL*(EQN_IN_LABELS_LEN),
     &  MAP_TYPE*(EQN_IN_MAP_TYPE_LEN),
     &  MAP_PARENT*(EQN_IN_MAP_PARENT_LEN),
     &  MAP_CHILD*(EQN_IN_MAP_CHILD_LEN)
     
      LOGICAL FOUND

      CALL ENTERS('UPDATE_EQUATION_INPUTS',*9999)

      CALL GET_EQUATION_TOTAL(EQUATION_TOTAL,ERROR,*9999)
            
      FOUND=.FALSE.
      index=0
      DO WHILE ((.NOT.FOUND).AND.(index.LE.EQUATION_TOTAL))
        index=index+1
        CALL GET_EQUATION_LABEL(index,EQUATION_LABEL,ERROR,*9999)
        CALL GET_EQUATION_NX(EQUATION_LABEL,lnx,ERROR,*9999)
        IF (LNX.EQ.nx) FOUND=.TRUE.
      ENDDO
      
      IF(FOUND)THEN
        
        
        CALL GET_EQUATION_INPUT_TOTAL(EQUATION_LABEL,
     &    EQN_IN_TOTAL,ERROR,*9999)
 
        DO index=1,EQN_IN_TOTAL
          
          CALL GET_EQUATION_INPUT_LABEL(EQUATION_LABEL,index,
     &      EQN_IN_LABEL,ERROR,*9999)
          CALL GET_EQUATION_INPUT_NM(EQUATION_LABEL,
     &      EQN_IN_LABEL,nm,ERROR,*9999)
          CALL GET_EQUATION_INPUT_MAP(EQUATION_LABEL,
     &      EQN_IN_LABEL,MAP_TYPE,MAP_PARENT,MAP_CHILD,
     &      ERROR,*9999)
          
          CALL STRING_TRIM(EQN_IN_LABEL,IB1,IE1)
          
          IF(MAP_TYPE.EQ.'INITIAL_VALUE')THEN
            
            CALL GET_EQUATION_INPUT_INIT_VALUE(EQUATION_LABEL,
     &        EQN_IN_LABEL,VALUE,ERROR,*9999)
            CG(nm)=VALUE

            CALL SET_EQUATION_INPUT_VALUE(EQUATION_LABEL,
     &        EQN_IN_LABEL,CG(nm),ERROR,*9999)
            
          ELSEIF(MAP_TYPE.EQ.'CONSTANT')THEN
            
            CALL GET_CONSTANT_COMPONENT_VALUE(MAP_PARENT,MAP_CHILD,
     &        VALUE,ERROR,*9999)
            CG(nm)=VALUE

            CALL SET_EQUATION_INPUT_VALUE(EQUATION_LABEL,
     &        EQN_IN_LABEL,CG(nm),ERROR,*9999)
            
          ELSEIF(MAP_TYPE.EQ.'FIELD')THEN
            
            CALL GET_FIELD_NJ(MAP_PARENT,nj,ERROR,*9999)
            CALL GET_FIELD_COMPONENT_NU(MAP_PARENT,MAP_CHILD,nu,
     &        ERROR,*9999)

C DMAL 13 APRIL 2004: This code is dodgy and hasn't been checked
C for non-unit scale factors
            IF(nu.EQ.1)THEN
              CG(nm)=XG(nj,nu)
            ELSEIF(nu.EQ.2)THEN
              CG(nm)=XG(nj,nu)/DXIX(1,1)
            ELSEIF(nu.EQ.3)THEN
              CG(nm)=XG(nj,nu)/DXIX(1,nj)/DXIX(1,nj)
            ENDIF
C DMAL 13 APRIL 2004

            CALL SET_EQUATION_INPUT_VALUE(EQUATION_LABEL,
     &        EQN_IN_LABEL,CG(nm),ERROR,*9999)
            
          ELSEIF(MAP_TYPE.EQ.'MATHS')THEN
            
            CALL GET_MATHS_SOLVED(MAP_PARENT,SOLVED,ERROR,*9999)
            IF(SOLVED.EQ.0)THEN
              CALL EVALUATE_MATHS(MAP_PARENT,XG,DXIX,ERROR,*9999)
            ENDIF
            
            CALL GET_MATHS_OUTPUT_VALUE(MAP_PARENT,MAP_CHILD,VALUE,
     &        ERROR,*9999)
            CG(nm)=VALUE

            CALL SET_EQUATION_INPUT_VALUE(EQUATION_LABEL,
     &        EQN_IN_LABEL,CG(nm),ERROR,*9999)
            
          ELSEIF(MAP_TYPE.EQ.'EQUATION')THEN
            CALL ASSERT(.FALSE.,'Not Implemented. '//
     &      'Due to the way cmiss is I can not get to solution '//
     &      'fields (YG) without a lot of work',ERROR,*9999)
          ELSE
            CALL ASSERT(.FALSE.,'Error: Unidentified MAP_TYPE',
     &        ERROR,*9999)
          ENDIF
        ENDDO
        
        
C DMAL 26 MARCH 2004: Personally I think this block of code has no business here
C since this subroutine should do only one job, update the values in CG. When
C this is removed, remove the COMMENT variable above.        
        CALL GET_EQUATION_VERBOSE(EQUATION_LABEL,VERBOSE,ERROR,*9999)
        IF(VERBOSE.GT.0) THEN
          WRITE(OP_STRING,'()')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          CALL STRING_TRIM(COMMENT,IB1,IE1)
          WRITE(OP_STRING,'(A)') COMMENT(IB1:IE1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          CALL STRING_TRIM(EQUATION_LABEL,IB1,IE1)
          WRITE(OP_STRING,'(''Equation '',A,'' Inputs:'')') 
     &      EQUATION_LABEL(IB1:IE1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          CALL PRINT_EQUATION_INPUT_VALUES(EQUATION_LABEL,ERROR,*9999)
        ENDIF
      
      ENDIF

      CALL EXITS('UPDATE_EQUATION_INPUTS')
      RETURN
 9999 CALL ERRORS('UPDATE_EQUATION_INPUTS',ERROR)
      CALL EXITS('UPDATE_EQUATION_INPUTS')
      RETURN 1
      END


