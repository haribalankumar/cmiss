      SUBROUTINE CMGUI_LINK_INITIALISE(ERROR,*)

C#### Subroutine: CMGUI_LINK_INITIALISE
C###  Description:
C###    Initialise any variables used by the CM-CMGUI
C###    link.  This should be called before any other
C###    calls to link functions.

C###  Note  We expect at least a 4 digit number in the
C###    argument.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cmgui00.cmn'

!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,NUMBER
      REAL*8 WORMHOLE_PROMPT_TIMEOUT,WORMHOLE_DATA_TIMEOUT,
     '  WORMHOLE_COMMAND_TIMEOUT
      PARAMETER(WORMHOLE_PROMPT_TIMEOUT=300.0d0,
     '  WORMHOLE_DATA_TIMEOUT=300.0d0,
     '  WORMHOLE_COMMAND_TIMEOUT=300.0d0)
      CHARACTER WORMHOLE_ARG*50

      CALL ENTERS('CMGUI_LINK_INITIALISE',*9999)

C     Pointers are initialized at startup in BLKCMGUI
      CMGUI_DATA_LEN=0
      CMGUI_NODE_LEN=0
      CMGUI_REGION_LEN=0
C     Initialise any status variables.
      IDLE_SENT=.FALSE.
      CMGUI_MESSAGE_TIME=0
      DATA_CHANGED=.FALSE.
      DATA_DELETED=.FALSE.
      NODES_CHANGED=.FALSE.
      NODES_DELETED=.FALSE.
      REGIONS_CHANGED=.FALSE.
      REGIONS_DELETED=.FALSE. !not sure if used / useable
C     Will probably have to start up the link here (command option)
      CMGUI_DATA_INPUT_STATE=1 !type
      IF(CMGUI_LINK) THEN
C       Check the arguments to create the wormholes
        CALL STRING_TRIM(CMGUI_LINK_ARG,IBEG,IEND)
        IF(IEND-IBEG.GE.3) THEN
          NUMBER=IFROMC(CMGUI_LINK_ARG(IEND-3:IEND))
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' CMGUI link number: '',I4)')
     '        NUMBER
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          IF(CMGUI_LINK_ARG(1:4).EQ.'file') THEN
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' CMGUI link is file based'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            WRITE(WORMHOLE_ARG,'('''//CMGUI_LINK_ARG(IBEG:IEND-4)
     '        //''',I4)') NUMBER+1
            CALL WH_INPUT_F_CREATE(WORMHOLE_ARG,' ',
     '        CMGUI_COMMAND_I,WORMHOLE_COMMAND_TIMEOUT)
            WRITE(WORMHOLE_ARG,'('''//CMGUI_LINK_ARG(IBEG:IEND-4)
     '        //''',I4)') NUMBER
            CALL WH_OUTPUT_F_CREATE(WORMHOLE_ARG,' ',
     '        CMGUI_COMMAND_O)
            WRITE(WORMHOLE_ARG,'('''//CMGUI_LINK_ARG(IBEG:IEND-4)
     '        //''',I4)') NUMBER+3
            CALL WH_INPUT_F_CREATE(WORMHOLE_ARG,' ',
     '        CMGUI_PROMPT_I,WORMHOLE_PROMPT_TIMEOUT)
            WRITE(WORMHOLE_ARG,'('''//CMGUI_LINK_ARG(IBEG:IEND-4)
     '        //''',I4)') NUMBER+2
            CALL WH_OUTPUT_F_CREATE(WORMHOLE_ARG,' ',
     '        CMGUI_PROMPT_O)
            WRITE(WORMHOLE_ARG,'('''//CMGUI_LINK_ARG(IBEG:IEND-4)
     '        //''',I4)') NUMBER+5
            CALL WH_INPUT_F_CREATE(WORMHOLE_ARG,' ',
     '        CMGUI_DATA_I,WORMHOLE_DATA_TIMEOUT)
            WRITE(WORMHOLE_ARG,'('''//CMGUI_LINK_ARG(IBEG:IEND-4)
     '        //''',I4)') NUMBER+4
            CALL WH_OUTPUT_F_CREATE(WORMHOLE_ARG,' ',
     '        CMGUI_DATA_O)
          ELSEIF(CMGUI_LINK_ARG(1:4).EQ.'sock') THEN
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' CMGUI link is socket based'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            WRITE(WORMHOLE_ARG,'('''//CMGUI_LINK_ARG(IBEG:IEND-4)
     '        //''',I4)') NUMBER+1
            CALL WH_INPUT_F_CREATE(WORMHOLE_ARG,' ',
     '        CMGUI_COMMAND_I,5.0d0)
            WRITE(WORMHOLE_ARG,'('''//CMGUI_LINK_ARG(IBEG:IEND-4)
     '        //''',I4)') NUMBER
            CALL WH_OUTPUT_F_CREATE(WORMHOLE_ARG,' ',
     '        CMGUI_COMMAND_O)
            WRITE(WORMHOLE_ARG,'('''//CMGUI_LINK_ARG(IBEG:IEND-4)
     '        //''',I4)') NUMBER+3
            CALL WH_INPUT_F_CREATE(WORMHOLE_ARG,' ',
     '        CMGUI_PROMPT_I,5.0d0)
            WRITE(WORMHOLE_ARG,'('''//CMGUI_LINK_ARG(IBEG:IEND-4)
     '        //''',I4)') NUMBER+2
            CALL WH_OUTPUT_F_CREATE(WORMHOLE_ARG,' ',
     '        CMGUI_PROMPT_O)
            WRITE(WORMHOLE_ARG,'('''//CMGUI_LINK_ARG(IBEG:IEND-4)
     '        //''',I4)') NUMBER+5
            CALL WH_INPUT_F_CREATE(WORMHOLE_ARG,' ',
     '        CMGUI_DATA_I,5.0d0)
            WRITE(WORMHOLE_ARG,'('''//CMGUI_LINK_ARG(IBEG:IEND-4)
     '        //''',I4)') NUMBER+4
            CALL WH_OUTPUT_F_CREATE(WORMHOLE_ARG,' ',
     '        CMGUI_DATA_O)
          ELSE
            ERROR='Invalid CMGUI link argument'
            GOTO 9999
          ENDIF
          IF((CMGUI_COMMAND_I.EQ.0).OR.(CMGUI_COMMAND_O.EQ.0).OR.
     '      (CMGUI_PROMPT_I.EQ.0).OR.(CMGUI_PROMPT_O.EQ.0).OR.
     '      (CMGUI_DATA_I.EQ.0).OR.(CMGUI_DATA_O.EQ.0)) THEN
            ERROR='Could not create CMGUI wormholes'
            GOTO 9999
          ENDIF
        ELSE
          ERROR='Invalid CMGUI number'
          GOTO 9999
        ENDIF
      ENDIF !CMGUI_LINK

      CALL EXITS('CMGUI_LINK_INITIALISE')
      RETURN
 9999 CALL ERRORS('CMGUI_LINK_INITIALISE',ERROR)
      CALL EXITS('CMGUI_LINK_INITIALISE')
      RETURN 1
      END


