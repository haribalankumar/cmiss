      SUBROUTINE PARSE_LINES(NLLINE,NLLIST,noco,NRLIST,NTCO,CO,
     '  ERROR,*)

C#### Subroutine: PARSE_LINES
C###  Description:
C###    PARSE_LINES parses command CO for keyword 'lines'
C###    and returns list of lines in
C###    NLLIST(nolist),nolist=1,NLLIST(0).
C**** Keyword 'lines' can be followed by a list or by a group name.
C**** Default is all lines. Max number in list is NL_R_M.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER NLLINE(0:NL_R_M,0:NRM),NLLIST(0:NLM),noco,
     '  NRLIST(0:NRM),NTCO
      CHARACTER CO(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER N3CO,
     '  nl,NL1,noline,nolist,no_nrlist,nr
      LOGICAL CBBREV,DEFINED
      CHARACTER CHAR1*4

      CALL ENTERS('PARSE_LINES',*9999)
      IF(CBBREV(CO,'LINES',1,noco+1,NTCO,N3CO)) THEN
        CDATA(1)='LINES'
        CALL PARSILG(NLLIST,NLM,CDATA(1),CO(N3CO+1),ERROR,*9999)

C 21/2/97 LC section archived :   old MPN 4-Apr-96: now uses PARSILG

C ***   Check that lines in list are defined in region list
        DO nolist=1,NLLIST(0)
          nL=NLLIST(nolist)
          DEFINED=.FALSE.
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO noline=1,NLLINE(0,nr)
              NL1=NLLINE(noline,nr)
              IF(nl.EQ.NL1) THEN
                DEFINED=.TRUE.
                GO TO 5
              ENDIF
            ENDDO
          ENDDO
  5       IF(.NOT.DEFINED) THEN
            WRITE(CHAR1,'(I4)') nl
            ERROR=' >>Line '//CHAR1(1:4)//' is not defined'
            GO TO 9999
          ENDIF
        ENDDO

      ELSE
C       Return all lines in region list
        NLLIST(0)=0
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          DO noline=1,NLLINE(0,nr)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' nr='',I2,'' noline='',I5,'
     '          //''' NLLINE(noline,nr)='',I5)')
     '          nr,noline,NLLINE(noline,nr)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            NLLIST(NLLIST(0)+noline)=NLLINE(noline,nr)
          ENDDO
          NLLIST(0)=NLLIST(0)+NLLINE(0,nr)
        ENDDO
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,*)' nllist(1..NLLIST(0))=',
     '    (NLLIST(nolist),nolist=1,NLLIST(0))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('PARSE_LINES')
      RETURN
 9999 CALL ERRORS('PARSE_LINES',ERROR)
      CALL EXITS('PARSE_LINES')
      RETURN 1
      END


