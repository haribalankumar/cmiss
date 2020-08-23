      SUBROUTINE PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,ERROR,*)

C#### Subroutine: PARSE_NODES
C###  Description:
C###    PARSE_NODES parses command CO for keyword 'nodes' and
C###    returns list of nodes in NPLIST(nolist),nolist=1,NPLIST(0).
C**** Keyword 'nodes' can be followed by a list or by a group name.
C**** Default is all nodes. Max number in list is NP_R_M.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER NPNODE(0:NP_R_M,0:NRM),NPLIST(0:NPM),noco,
     '  NRLIST(0:NRM),NTCO
      CHARACTER CO(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER
     '  N3CO,N1LIST,
     '  nolist,nonode,no_nrlist,np,NP1,nr
      LOGICAL CBBREV,DEFINED,INLIST
      CHARACTER CHAR1*4

      CALL ENTERS('PARSE_NODES',*9999)
      IF(CBBREV(CO,'NODES',1,noco+1,NTCO,N3CO)) THEN
        CDATA(1)='NODES'
        CALL PARSILG(NPLIST,NPM,CDATA(1),CO(N3CO+1),ERROR,*9999)

C 21/2/97 LC archived section : old MPN 4-Apr-96: now uses PARSILG

C       Check that nodes in list are defined in region list
        DO nolist=1,NPLIST(0)
          np=NPLIST(nolist)
          DEFINED=.FALSE.
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO nonode=1,NPNODE(0,nr)
              NP1=NPNODE(nonode,nr)
              IF(np.EQ.NP1) THEN
                DEFINED=.TRUE.
                GO TO 5
              ENDIF
            ENDDO
          ENDDO
  5       IF(.NOT.DEFINED) THEN
            WRITE(CHAR1,'(I4)') np
            ERROR=' >>Node '//CHAR1(1:4)//' is not defined'
            GO TO 9999
          ENDIF
        ENDDO

      ELSE
C       Return all nodes in region list
        NPLIST(0)=0
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          DO nonode=1,NPNODE(0,nr)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' nr='',I2,'' nonode='',I5,'
     '          //''' NPNODE(nonode,nr)='',I5)')
     '          nr,nonode,NPNODE(nonode,nr)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            IF(.NOT.INLIST(NPNODE(nonode,nr),NPLIST(1),NPLIST(0),
     '        N1LIST)) THEN
              NPLIST(0)=NPLIST(0)+1
              IF(NPLIST(0).LE.NP_R_M) THEN
                NPLIST(NPLIST(0))=NPNODE(nonode,nr)
              ENDIF
            ENDIF
          ENDDO !nonode (np)
        ENDDO !nr
        CALL ASSERT(NPLIST(0).LE.NP_R_M,'>>Increase NP_R_M',ERROR,
     '    *9999)
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,*)' nplist(1..NPLIST(0))=',
     '    (NPLIST(nolist),nolist=1,NPLIST(0))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('PARSE_NODES')
      RETURN
 9999 CALL ERRORS('PARSE_NODES',ERROR)
      CALL EXITS('PARSE_NODES')
      RETURN 1
      END


