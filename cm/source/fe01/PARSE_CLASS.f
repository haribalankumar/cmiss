      SUBROUTINE PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*)

C#### Subroutine: PARSE_CLASS
C###  Description:
C###    PARSE_CLASS parses command CO for keyword 'class' and
C###    returns list of classes specified in
C###    NXLIST(nolist),nolist=1,NXLIST(0).
C***  Default is 1. Max number is 9.
C***  Now handles lists MLB 1 September 1997

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER NXLIST(0:NXM),noco,NTCO
      CHARACTER CO(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER N3CO,nxc
      LOGICAL CBBREV

      CALL ENTERS('PARSE_CLASS',*9999)

      IF(CBBREV(CO,'CLASS',1,noco+1,NTCO,N3CO)) THEN
        CALL PARSIL(CO(N3CO+1),NXM,NXLIST(0),NXLIST(1),ERROR,*9999)
      ELSE
        NXLIST(0)=1
        NXLIST(1)=1
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(''NXLIST: '',9I3)')
     '    (NXLIST(nxc),nxc=1,NXLIST(0))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('PARSE_CLASS')
      RETURN
 9999 CALL ERRORS('PARSE_CLASS',ERROR)
      CALL EXITS('PARSE_CLASS')
      RETURN 1
      END


