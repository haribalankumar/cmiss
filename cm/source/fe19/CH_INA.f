      SUBROUTINE CH_INA(CALL_FLAG,QTYPE,QUESTIONS,RVALUE,ERROR,*)

C#### Subroutine: CH_INA
C###  Description:
C###    CH_INA models the kinetics of the INa (sodium) channel.
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
      REAL*8 ENa,ROF
!      LOGICAL
!     Parameter list
      PARAMETER(ROF=8.617342D-2) !R/F

      CALL ENTERS('CH_INA',*9999)

      IF(CALL_FLAG(1:6).EQ.'IPMATE') THEN
        QTYPE(0,1)=6 ! number of questions
        QUESTIONS(1)='sodium current, INa'
        QTYPE(1,1)=0 ! I/O flag, 0=output, 1=input
        QTYPE(1,2)=1 ! Multi ask flag, 1=ask once, 2=ask until exit
        QUESTIONS(2)='transmembrane voltage, Vm'
        QTYPE(2,1)=1
        QTYPE(2,2)=1
        QUESTIONS(3)='sodium concentration, [Na+]o (outside)'
        QTYPE(3,1)=1
        QTYPE(3,2)=1
        QUESTIONS(4)='sodium concentration, [Na+]i (inside)'
        QTYPE(4,1)=1
        QTYPE(4,2)=1
        QUESTIONS(5)='sodium channel conductance, gna'
        QTYPE(5,1)=1
        QTYPE(5,2)=1
        QUESTIONS(6)='temperature (K)'
        QTYPE(6,1)=1
        QTYPE(6,2)=1
      ELSEIF(CALL_FLAG(1:5).EQ.'EVALU') THEN
        ENa=(ROF*RVALUE(6))*DLOG(RVALUE(3)/RVALUE(4))
        RVALUE(1)=RVALUE(5)*(RVALUE(2)-ENa)
      ENDIF

      CALL EXITS('CH_INA')
      RETURN
 9999 CALL ERRORS('CH_INA',ERROR)
      CALL EXITS('CH_INA')
      RETURN 1
      END


