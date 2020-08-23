      SUBROUTINE PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '  ERROR,*)

C#### Subroutine: PARSE_ELEMENTS
C###  Description:
C###    PARSE_ELEMENTS parses command CO for keyword 'elements'
C###    and returns list of elements in
C###    NELIST(nolist),nolist=1,NELIST(0).
C**** Keyword 'elements' can be followed by a list or by a group name.
C**** Default is all elements. Max number in list is NE_R_M.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),noco,
     '  NRLIST(0:NRM),NTCO
      CHARACTER CO(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER N3CO,
     '  ne,NE1,noelem,nolist,no_nrlist,nr
      LOGICAL CBBREV,DEFINED
      CHARACTER CHAR1*4

      CALL ENTERS('PARSE_ELEMENTS',*9999)
      IF(CBBREV(CO,'ELEMENTS',1,noco+1,NTCO,N3CO)) THEN
        CDATA(1)='ELEMENTS'
        CALL PARSILG(NELIST,NEM,CDATA(1),CO(N3CO+1),ERROR,*9999)

C 21/2/97 LC section archived :   old MPN 4-Apr-96: now uses PARSILG

C ***   Check that elements in list are defined in region list
        DO nolist=1,NELIST(0)
          ne=NELIST(nolist)
          DEFINED=.FALSE.
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO noelem=1,NEELEM(0,nr)
              NE1=NEELEM(noelem,nr)
              IF(ne.EQ.NE1) THEN
                DEFINED=.TRUE.
                GO TO 5
              ENDIF
            ENDDO
          ENDDO
  5       IF(.NOT.DEFINED) THEN
            WRITE(CHAR1,'(I4)') ne
            ERROR=' >>Element '//CHAR1(1:4)//' is not defined'
            GO TO 9999
          ENDIF
        ENDDO

      ELSE
C       Return all elements in region list
        NELIST(0)=0
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          DO noelem=1,NEELEM(0,nr)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' nr='',I2,'' noelem='',I5,'
     '          //''' NEELEM(noelem,nr)='',I5)')
     '          nr,noelem,NEELEM(noelem,nr)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            NELIST(NELIST(0)+noelem)=NEELEM(noelem,nr)
          ENDDO
          NELIST(0)=NELIST(0)+NEELEM(0,nr)
        ENDDO
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,*)' nelist(1..NELIST(0))=',
     '    (NELIST(nolist),nolist=1,NELIST(0))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('PARSE_ELEMENTS')
      RETURN
 9999 CALL ERRORS('PARSE_ELEMENTS',ERROR)
      CALL EXITS('PARSE_ELEMENTS')
      RETURN 1
      END


