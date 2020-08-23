      SUBROUTINE PF_LINEAR(CALL_FLAG,QTYPE,QUESTIONS,RVALUE,ERROR,*)

C#### Subroutine: PF_LINEAR
C###  Description:
C###    PF_LINEAR is a linear polynomial function y=ax+b.
C**** Written by Duane Malcolm, 17 January 2003

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

      CALL ENTERS('PF_LINEAR',*9999)

      IF(CALL_FLAG(1:6).EQ.'IPMATE') THEN
        QTYPE(0,1)=4 ! number of questions
        QUESTIONS(1)='y'
        QTYPE(1,1)=0 ! I/O flag, 0=output, 1=input
        QTYPE(1,2)=1 ! Multi ask flag, 1=ask once, 2=ask until exit
        QUESTIONS(2)='x'
        QTYPE(2,1)=1
        QTYPE(2,2)=1
        QUESTIONS(3)='a'
        QTYPE(3,1)=1
        QTYPE(3,2)=1
        QUESTIONS(4)='b'
        QTYPE(4,1)=1
        QTYPE(4,2)=1
      ELSEIF(CALL_FLAG(1:5).EQ.'EVALU') THEN
        RVALUE(1)=RVALUE(3)*RVALUE(2)+RVALUE(4)
      ENDIF

      CALL EXITS('PF_LINEAR')
      RETURN
 9999 CALL ERRORS('PF_LINEAR',ERROR)
      CALL EXITS('PF_LINEAR')
      RETURN 1
      END


