      SUBROUTINE CH_ICL(CALL_FLAG,QTYPE,QUESTIONS,RVALUE,ERROR,*)

C#### Subroutine: CH_CL
C###  Description:
C###    CH_CL models the kinetics of the ICl (chloride) channel.
C**** Written by Duane Malcolm, 28 September 2002

      IMPLICIT NONE

      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER QTYPE(0:10,2)
      REAL*8 RVALUE(99)
      CHARACTER CALL_FLAG*(*),QUESTIONS(10)*64,ERROR*(*)
!     Local Variables
!      INTEGER
!      REAL*8
!      LOGICAL

      CALL ENTERS('CH_ICL',*9999)

      IF(CALL_FLAG(1:6).EQ.'IPMATE') THEN
        QTYPE(0,1)=5 ! number of questions
        QUESTIONS(1)='chloride current, ICl'
        QTYPE(1,1)=0 ! I/O flag, 0=output, 1=input
        QTYPE(1,2)=1 ! Multi ask flag, 1=ask once, 2=ask until exit
        QUESTIONS(2)='chloride concentration, [Cl-]o (outside)'
        QTYPE(2,1)=1
        QTYPE(2,2)=1
        QUESTIONS(3)='chloride concentration, [Cl-]i (inside)'
        QTYPE(3,1)=1
        QTYPE(3,2)=1
        QUESTIONS(4)='chloride channel conductance, gcl'
        QTYPE(4,1)=1
        QTYPE(4,2)=1
        QUESTIONS(5)='density of the chloride channels'
        QTYPE(5,1)=1
        QTYPE(5,2)=1
      ELSEIF(CALL_FLAG(1:5).EQ.'EVALU') THEN
        RVALUE(1)=RVALUE(5)*RVALUE(4)*(RVALUE(2)-RVALUE(3))
      ENDIF

      CALL EXITS('CH_ICL')
      RETURN
 9999 CALL ERRORS('CH_ICL',ERROR)
      CALL EXITS('CH_ICL')
      RETURN 1
      END


