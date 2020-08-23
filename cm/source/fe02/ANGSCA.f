      SUBROUTINE ANGSCA(NPL,DL,XP,ERROR,*)

C#### Subroutine: ANGSCA
C###  Description:
C###    ANGSCA calculates global scale factors from angle change or
C###    change in radial coordinate if Xi(3) direction.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER NPL(5,0:3)
      REAL*8 DL(3),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER N,ni,NP1,NP2

      CALL ENTERS('ANGSCA',*9999)

      ni=NPL(1,0)
      NP1=NPL(2,1)
      NP2=NPL(3,1)
      IF(ITYP10(1).EQ.1) THEN       !rectangular cartesian
        ERROR=' Scale factors cannot be calc.d from angle change'
     '    //' when coords are rect. cart.'
        GO TO 9999
      ELSE IF(ITYP10(1).EQ.2) THEN  !cylindrical polar
        IF(ni.EQ.1) THEN
          DL(1)=XP(1,1,2,NP2)-XP(1,1,2,NP1)
          IF(DL(1).LT.0.0d0) DL(1)=DL(1)+2.0d0*PI
        ELSE IF(ni.EQ.2) THEN
          IF(NJT.EQ.2) THEN
            DL(1)=XP(1,1,1,NP2)-XP(1,1,1,NP1)
          ELSE IF(NJT.EQ.3) THEN
            DL(1)=XP(1,1,3,NP2)-XP(1,1,3,NP1)
          ENDIF
        ENDIF
      ELSE IF(ITYP10(1).EQ.3) THEN  !spherical polar
        ERROR=' Not implemented'
        GO TO 9999
      ELSE IF(ITYP10(1).EQ.4) THEN  !prolate spheroidal
        IF(ni.EQ.1) THEN
          DL(1)=XP(1,1,3,NP1)-XP(1,1,3,NP2)
          IF(DL(1).LT.0.0d0) DL(1)=DL(1)+2.0d0*PI
        ELSE IF(ni.EQ.2) THEN
          DL(1)=XP(1,1,2,NP2)-XP(1,1,2,NP1)
          IF(DL(1).LT.0.0d0) DL(1)=DL(1)+2.0d0*PI
        ELSE IF(ni.EQ.3) THEN
          DL(1)=XP(1,1,1,NP2)-XP(1,1,1,NP1)
        ENDIF
      ELSE IF(ITYP10(1).EQ.5) THEN  !oblate spheroidal
        ERROR=' Not implemented'
        GO TO 9999
      ENDIF
      DL(2)=DL(1)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' DL(n):'',3E13.5)')
     '    (DL(N),N=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ANGSCA')
      RETURN
 9999 CALL ERRORS('ANGSCA',ERROR)
      CALL EXITS('ANGSCA')
      RETURN 1
      END


