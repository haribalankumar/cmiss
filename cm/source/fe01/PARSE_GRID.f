      SUBROUTINE PARSE_GRID(NQLIST,noco,NTCO,CO,ERROR,*)

C#### Subroutine: PARSE_GRID
C###  Description:
C###    PARSE_GRID parses command CO for keyword 'grid' and
C###    returns list of grid points in
C###    NQLIST(nolist),nolist=1,NQLIST(0).
C**** Keyword 'grid' can be followed by a list or by a group name.
C**** Default is all grid points. Max number in list is NQM.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER NQLIST(0:NQM),noco,NTCO
      CHARACTER CO(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER N3CO,nolist,nq
      LOGICAL CBBREV

      CALL ENTERS('PARSE_GRID',*9999)
      IF(CBBREV(CO,'GRIDS',1,noco+1,NTCO,N3CO)) THEN
        CDATA(1)='GRIDS'
        CALL PARSILG(NQLIST,NQM,CDATA(1),CO(N3CO+1),ERROR,*9999)

C 21/2/97 LC section to archive : C old MPN 4-Apr-96: now uses PARSILG

      ELSE
C       Return all grid points
        DO nq=1,NQT
          NQLIST(nq)=nq
        ENDDO
        NQLIST(0)=NQT
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,*)' nqlist(1..NQLIST(0))=',
     '    (NQLIST(nolist),nolist=1,NQLIST(0))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('PARSE_GRID')
      RETURN
 9999 CALL ERRORS('PARSE_GRID',ERROR)
      CALL EXITS('PARSE_GRID')
      RETURN 1
      END


