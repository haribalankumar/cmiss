      SUBROUTINE PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*)

C#### Subroutine: PARSE_REGIONS
C###  Description:
C###    PARSE_REGIONS parses command CO for keyword 'region' and
C###    returns list of regions in NRLIST(nolist),nolist=1,NRLIST(0).
C**** Keyword 'region' can be followed by a list or by 'all'.
C**** Default is 1. Max number in list is 9.


C*** LKC 22-MAR-1998 Modified "ALL_REGIONS" so that it refers
C    to more than one region, ie not necessarily all. Now if
C    we define a group of regions eg reg 1..3 it will look for an
C    *ir* file. This needs to be tidied further !

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER noco,NTCO,NRLIST(0:NRM)
      CHARACTER CO(*)*(*),ERROR*(*)
      LOGICAL ALL_REGIONS
!     Local Variables
      INTEGER N3CO,N4CO,nolist,nr
      LOGICAL CBBREV

      CALL ENTERS('PARSE_REGIONS',*9999)
      IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
        IF(CBBREV(CO,'ALL',1,N3CO+1,N3CO+2,N4CO)) THEN

C LKC 2-JUN-1998 if 'ALL' yet there really is only 1 region look for ip* file
          IF(NRT.GT.1) THEN
            ALL_REGIONS=.TRUE.
          ELSE
            ALL_REGIONS=.FALSE.
          ENDIF

          NRLIST(0)=NRT
          DO nolist=1,NRLIST(0)
            NRLIST(nolist)=nolist
          ENDDO
        ELSE
C LKC 22-MAR-1998          ALL_REGIONS=.FALSE.
          CALL PARSIL(CO(N3CO+1),NRM,NRLIST(0),NRLIST(1),ERROR,*9999)
          DO nolist=1,NRLIST(0)
            nr=NRLIST(nolist)
            CALL ASSERT(nrt.ge.nr,' >>Region not defined',
     '        ERROR,*9999)
C LKC 10-MAY-2002 new check for regions
            CALL ASSERT(nr.GT.0,' >>Region number 0 specified',
     '        ERROR,*9999)
          ENDDO

C LKC start  22-MAR-1998
C changed so ALL_REGIONS really refers to multiple regions
          IF(NRLIST(0).GT.1)THEN
            ALL_REGIONS=.TRUE.
          ELSE
            ALL_REGIONS=.FALSE.
          ENDIF
C LKC end
        ENDIF
      ELSE
        ALL_REGIONS=.FALSE.
        NRLIST(0)=1 !Default number
        NRLIST(1)=1
      ENDIF

      CALL EXITS('PARSE_REGIONS')
      RETURN
 9999 CALL ERRORS('PARSE_REGIONS',ERROR)
      CALL EXITS('PARSE_REGIONS')
      RETURN 1
      END


