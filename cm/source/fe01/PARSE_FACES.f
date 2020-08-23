      SUBROUTINE PARSE_FACES(NFFACE,NFLIST,noco,NRLIST,NTCO,CO,
     '  ERROR,*)

C#### Subroutine: PARSE_FACES
C###  Description:
C###    PARSE_FACES parses command CO for keyword 'faces'
C###    and returns list of faces in
C###    NFLIST(nolist),nolist=1,NFLIST(0).
C**** Keyword 'faces' can be followed by a list or by a group name.
C**** Default is all faces. Max number in list is NF_R_M.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER NFFACE(0:NF_R_M,NRM),NFLIST(0:NFM),noco,
     '  NRLIST(0:NRM),NTCO
      CHARACTER CO(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER N3CO,
     '  nf,NF1,noface,nolist,no_nrlist,nr
      LOGICAL CBBREV,DEFINED
      CHARACTER CHAR1*4

      CALL ENTERS('PARSE_FACES',*9999)
      IF(CBBREV(CO,'FACES',1,noco+1,NTCO,N3CO)) THEN
        CDATA(1)='FACES'
        CALL PARSILG(NFLIST,NFM,CDATA(1),CO(N3CO+1),ERROR,*9999)

C 21/2/97 LC section archived :   old MPN 4-Apr-96: now uses PARSILG

C ***   Check that faces in list are defined in region list
        DO nolist=1,NFLIST(0)
          nf=NFLIST(nolist)
          DEFINED=.FALSE.
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO noface=1,NFFACE(0,nr)
              NF1=NFFACE(noface,nr)
              IF(nf.EQ.NF1) THEN
                DEFINED=.TRUE.
                GO TO 5
              ENDIF
            ENDDO
          ENDDO
  5       IF(.NOT.DEFINED) THEN
            WRITE(CHAR1,'(I4)') nf
            ERROR=' >>Face '//CHAR1(1:4)//' is not defined'
            GO TO 9999
          ENDIF
        ENDDO

      ELSE
C       Return all faces in region list
        NFLIST(0)=0
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          DO noface=1,NFFACE(0,nr)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' nr='',I2,'' noface='',I5,'
     '          //''' NFFACE(noface,nr)='',I5)')
     '          nr,noface,NFFACE(noface,nr)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            NFLIST(NFLIST(0)+noface)=NFFACE(noface,nr)
          ENDDO
          NFLIST(0)=NFLIST(0)+NFFACE(0,nr)
        ENDDO
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,*)' nflist(1..NFLIST(0))=',
     '    (NFLIST(nolist),nolist=1,NFLIST(0))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('PARSE_FACES')
      RETURN
 9999 CALL ERRORS('PARSE_FACES',ERROR)
      CALL EXITS('PARSE_FACES')
      RETURN 1
      END


