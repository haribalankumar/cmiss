      SUBROUTINE PARSE_WINDOWS(IWK,MAXNTIW,noco,NTCO,NTIW,CO,
     '  ALL_WINDOWS,ERROR,*)

C#### Subroutine: PARSE_WINDOWS
C###  Description:
C###    PARSE_WINDOES parses command CO for keyword 'window' and
C###    returns list of windows in IWK(nolist),nolist=1,NTIW.
C**** Keyword 'window' can be followed by a list or by 'all'.
C**** Default is ALL. Max number in list is MAXNTIW.

      IMPLICIT NONE
      INCLUDE 'cbwk01.cmn'

!     Parameter List
      INTEGER IWK(*),MAXNTIW,noco,NTCO,NTIW
      CHARACTER CO(*)*(*),ERROR*(*)
      LOGICAL ALL_WINDOWS
!     Local Variables
      INTEGER N3CO,N4CO,nolist,iw
      LOGICAL CBBREV

      CALL ENTERS('PARSE_WINDOWS',*9999)

      IF(CBBREV(CO,'WINDOW',1,noco+1,NTCO,N3CO)) THEN
        IF(CBBREV(CO,'ALL',1,N3CO+1,N3CO+2,N4CO)) THEN
          ALL_WINDOWS=.TRUE.
        ELSE
          ALL_WINDOWS=.FALSE.
          CALL PARSIL(CO(N3CO+1),MAXNTIW,NTIW,IWK,ERROR,*9999)
        ENDIF
      ELSE
        ALL_WINDOWS=.TRUE.
      ENDIF

      IF(ALL_WINDOWS) THEN
        NTIW=0
        DO nolist=1,IWKDEF(0)
          iw=IWKDEF(nolist)
          IF(IWKG(iw).EQ.1) THEN !graphics ouput window
            NTIW=NTIW+1
            IF(NTIW.LE.MAXNTIW) IWK(NTIW)=iw
          ENDIF
        ENDDO !nolist
        CALL ASSERT(NTIW.LE.MAXNTIW,'>>Too many values to store in '
     '      //'the list',ERROR,*9999)
      ENDIF

      CALL EXITS('PARSE_WINDOWS')
      RETURN
 9999 CALL ERRORS('PARSE_WINDOWS',ERROR)
      CALL EXITS('PARSE_WINDOWS')
      RETURN 1
      END


