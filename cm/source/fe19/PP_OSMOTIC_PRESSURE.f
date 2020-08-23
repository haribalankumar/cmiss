      SUBROUTINE PP_OSMOTIC_PRESSURE(CALL_FLAG,QTYPE,QUESTIONS,RVALUE,
     '  ERROR,*)

C#### Subroutine: PP_OSMOTIC_PRESSURE
C###  Description:
C###    PP_OSMOTIC_PRESSURE calculates the osmotic pressure between two
C###      domains.
C**** Written by Duane Malcolm, 1 November 2002

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
      REAL*8 Faraday,Rgas
!      LOGICAL
      PARAMETER(Faraday=96.5d0)
      PARAMETER(Rgas=8.314d0)


      CALL ENTERS('PP_OSMOTIC_PRESSURE',*9999)

      IF(CALL_FLAG(1:6).EQ.'IPMATE') THEN
        QTYPE(0,1)=4 ! number of questions
        QUESTIONS(1)='osmotic pressure'
        QTYPE(1,1)=0 ! I/O flag, 0=output, 1=input
        QTYPE(1,2)=1 ! Multi ask flag, 1=ask once, 2=ask until exit
        QUESTIONS(2)='total concentration, [C]o (outside)'
        QTYPE(2,1)=1
        QTYPE(2,2)=1
        QUESTIONS(3)='total concentration, [C]i (inside)'
        QTYPE(3,1)=1
        QTYPE(3,2)=1
        QUESTIONS(4)='temperature, T'
        QTYPE(4,1)=1
        QTYPE(4,2)=1
      ELSEIF(CALL_FLAG(1:5).EQ.'EVALU') THEN
        RVALUE(1)=Rgas*RVALUE(4)*(RVALUE(2)-RVALUE(3))
      ENDIF

      CALL EXITS('PP_OSMOTIC_PRESSURE')
      RETURN
 9999 CALL ERRORS('PP_OSMOTIC_PRESSURE',ERROR)
      CALL EXITS('PP_OSMOTIC_PRESSURE')
      RETURN 1
      END


