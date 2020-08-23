      SUBROUTINE PUMP_IKNA(CALL_FLAG,QTYPE,QUESTIONS,RVALUE,ERROR,*)

C#### Subroutine: PUMP_IKNA
C###  Description:
C###    PUMP_IKNA models the kinetics of the IKNa-ATPase
C###    (sodium/ptassium ATPase) pump.
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

      CALL ENTERS('PUMP_IKNA',*9999)

      IF(CALL_FLAG(1:6).EQ.'IPMATE') THEN
        QTYPE(0,1)=7 ! number of questions
        QUESTIONS(1)='sodium-potassium current, INaK'
        QTYPE(1,1)=0 ! I/O flag, 0=output, 1=input
        QTYPE(1,2)=1 ! Multi ask flag, 1=ask once, 2=ask until exit
        QUESTIONS(2)='potassium concentration, [K+]o (outside)'
        QTYPE(2,1)=1
        QTYPE(2,2)=1
        QUESTIONS(3)='sodium concentration, [Na+]i (inside)'
        QTYPE(3,1)=1
        QTYPE(3,2)=1
        QUESTIONS(4)='potassium half-maximal concentration, KMK'
        QTYPE(4,1)=1
        QTYPE(4,2)=1
        QUESTIONS(5)='sodium half-maximal concentration, KMNa'
        QTYPE(5,1)=1
        QTYPE(5,2)=1
        QUESTIONS(6)='maximum IKNa-ATPase pump current, IKNaMax'
        QTYPE(6,1)=1
        QTYPE(6,2)=1
        QUESTIONS(7)='density of the IKNA-ATPase pumps'
        QTYPE(7,1)=1
        QTYPE(7,2)=1
      ELSEIF(CALL_FLAG(1:5).EQ.'EVALU') THEN
        RVALUE(1)=RVALUE(7)*RVALUE(6)*(RVALUE(2)/(RVALUE(2)+RVALUE(4)))*
     '            (RVALUE(3)/(RVALUE(3)+RVALUE(5)))
      ENDIF

      CALL EXITS('PUMP_IKNA')
      RETURN
 9999 CALL ERRORS('PUMP_IKNA',ERROR)
      CALL EXITS('PUMP_IKNA')
      RETURN 1
      END


