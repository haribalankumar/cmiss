      SUBROUTINE COUP_QUES(COUP_TYPE,COUP_SUBTYPE,QTYPE,
     '  QUESTIONS,ERROR,*)


C#### Subroutine: COUP_QUES
C###  Description:
C###    COUP_QUES queries the model for questions to ask
C**** Written by Duane Malcolm, 26 August 2002

      IMPLICIT NONE

      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'

!     Parameter List
      INTEGER COUP_TYPE,COUP_SUBTYPE,QTYPE(0:10,2)
!      REAL*8
      CHARACTER QUESTIONS(10)*64,ERROR*(*)
!      LOGICAL
!     Local Variables
      INTEGER i
!      REAL*8

      CALL ENTERS('COUP_QUES',*9999)

      QTYPE(0,1)=0
      DO i=1,10
        QTYPE(i,1)=0
        QTYPE(i,2)=0
        WRITE(QUESTIONS(i),'('''')')
      ENDDO

      IF(COUP_TYPE.EQ.2) THEN ! membrane channel models
        IF(COUP_SUBTYPE.EQ.1) THEN ! INa channel
          CALL CH_INA('IPMATE',QTYPE,QUESTIONS,
     '      %VAL(0),ERROR,*9999)
        ELSEIF(COUP_SUBTYPE.EQ.2) THEN ! IK channel
          CALL CH_IK('IPMATE',QTYPE,QUESTIONS,
     '      %VAL(0),ERROR,*9999)
        ELSEIF(COUP_SUBTYPE.EQ.3) THEN ! ICl channel
          CALL CH_ICL('IPMATE',QTYPE,QUESTIONS,
     '      %VAL(0),ERROR,*9999)
        ENDIF
      ELSEIF(COUP_TYPE.EQ.3) THEN ! membrane pump models
        IF(COUP_SUBTYPE.EQ.1) THEN ! IKNa-ATPase pump
          CALL PUMP_IKNA('IPMATE',QTYPE,QUESTIONS,
     '      %VAL(0),ERROR,*9999)
        ENDIF
      ELSEIF(COUP_TYPE.EQ.6) THEN ! physical process models
        IF(COUP_SUBTYPE.EQ.1) THEN ! Osmotic pressure
          CALL PP_OSMOTIC_PRESSURE('IPMATE',QTYPE,QUESTIONS,
     '      %VAL(0),ERROR,*9999)
        ENDIF
      ELSEIF(COUP_TYPE.EQ.7) THEN ! polynomial functions
        IF(COUP_SUBTYPE.EQ.1) THEN ! linear
          CALL PF_LINEAR('IPMATE',QTYPE,QUESTIONS,
     '      %VAL(0),ERROR,*9999)
        ELSEIF(COUP_SUBTYPE.EQ.2) THEN ! quadratic
          CALL PF_QUADRATIC('IPMATE',QTYPE,QUESTIONS,
     '      %VAL(0),ERROR,*9999)
        ELSEIF(COUP_SUBTYPE.EQ.3) THEN ! cubic
          CALL PF_CUBIC('IPMATE',QTYPE,QUESTIONS,
     '      %VAL(0),ERROR,*9999)
        ENDIF
      ENDIF


      CALL EXITS('COUP_QUES')
      RETURN
 9999 CALL ERRORS('COUP_QUES',ERROR)
      CALL EXITS('COUP_QUES')
      RETURN 1
      END


