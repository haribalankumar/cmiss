      SUBROUTINE OPIMPO(ERROR,*)

C#### Subroutine: OPIMPO
C###  Description:
C###    OPIMPO outputs import parameters.

      IMPLICIT NONE
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'emap00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'impo00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj

      CALL ENTERS('OPIMPO',*9999)

      CALL ASSERT(CALL_IMPO,
     '  '>>Import parameters not defined',ERROR,*9999)

      IF(DEFIMPO_TYPE.EQ.1) THEN !Signal
        IF(SIGIMPO_TYPE.EQ.1) THEN !UNEMAP
          WRITE(OP_STRING,'(/'' UNEMAP import signal paramters'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(/'' Sock origin: '',3D12.4)')
     '      (EMAP_SOCKORIGIN(nj),nj=1,NJT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Sock apex: '',3D12.4)')
     '      (EMAP_SOCKAPEX(nj),nj=1,NJT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Sock theta zero: '',3D12.4)')
     '      (EMAP_SOCKTHETAZERO(nj),nj=1,NJT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Sock focus: '',D12.4)') EMAP_SOCKFOCUS(0)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('OPIMPO')
      RETURN
 9999 CALL ERRORS('OPIMPO',ERROR)
      CALL EXITS('OPIMPO')
      RETURN 1
      END


