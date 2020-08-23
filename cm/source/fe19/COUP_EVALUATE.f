      SUBROUTINE COUP_EVALUATE(COUP_TYPE,COUP_SUBTYPE,RVALUE,ERROR,*)


C#### Subroutine: COUP_QUES
C###  Description:
C###    COUP_QUES queries the model for questions to ask
C**** Written by Duane Malcolm, 26 August 2002

      IMPLICIT NONE

      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'

!     Parameter List
      INTEGER COUP_TYPE,COUP_SUBTYPE
      REAL*8 RVALUE(99)
      CHARACTER ERROR*(*)
!      LOGICAL
!     Local Variables
!      INTEGER
!      REAL*8

      CALL ENTERS('COUP_EVALUATE',*9999)

      IF(COUP_TYPE.EQ.2) THEN ! membrane channel models
        IF(COUP_SUBTYPE.EQ.1) THEN ! INa channel
          CALL CH_INA('EVALU',%VAL(0),%VAL(0),RVALUE,ERROR,*9999)
        ELSEIF(COUP_SUBTYPE.EQ.2) THEN ! IKNa channel
          CALL CH_IK('EVALU',%VAL(0),%VAL(0),RVALUE,ERROR,*9999)
        ELSEIF(COUP_SUBTYPE.EQ.3) THEN ! IKNa channel
          CALL CH_ICL('EVALU',%VAL(0),%VAL(0),RVALUE,ERROR,*9999)
        ENDIF
      ELSEIF(COUP_TYPE.EQ.3) THEN ! membrane pump models
        IF(COUP_SUBTYPE.EQ.1) THEN ! IKNa-ATPase pump
          CALL PUMP_IKNA('EVALU',%VAL(0),%VAL(0),RVALUE,ERROR,*9999)
        ENDIF
      ELSEIF(COUP_TYPE.EQ.6) THEN ! physical process models
        IF(COUP_SUBTYPE.EQ.1) THEN ! Osmotic pressure
          CALL PP_OSMOTIC_PRESSURE('EVALU',%VAL(0),%VAL(0),RVALUE,
     '      ERROR,*9999)
        ENDIF
      ELSEIF(COUP_TYPE.EQ.7) THEN ! polynomial functions
        IF(COUP_SUBTYPE.EQ.1) THEN ! linear
          CALL PF_LINEAR('EVALU',%VAL(0),%VAL(0),RVALUE,
     '      ERROR,*9999)
        ELSEIF(COUP_SUBTYPE.EQ.2) THEN ! quadratic
          CALL PF_QUADRATIC('EVALU',%VAL(0),%VAL(0),RVALUE,
     '      ERROR,*9999)
        ELSEIF(COUP_SUBTYPE.EQ.3) THEN ! cubic
          CALL PF_CUBIC('EVALU',%VAL(0),%VAL(0),RVALUE,
     '      ERROR,*9999)
        ENDIF
      ENDIF


      CALL EXITS('COUP_EVALUATE')
      RETURN
 9999 CALL ERRORS('COUP_EVALUATE',ERROR)
      CALL EXITS('COUP_EVALUATE')
      RETURN 1
      END


