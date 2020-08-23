      SUBROUTINE WS_LIST(IWK,IW_DEFAULT,NTIW,noco,NTCO,CO,ERROR,*)

C#### Subroutine: WS_LIST
C###  Description:
C###    WS_LIST returns list of workstations for segment creation.
C**** List is limited to IW_DEFAULT if this is non-zero, else all open
C**** windows.

      IMPLICIT NONE
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'post00.cmn'
!     Parameter List
      INTEGER IWK(6),IW_DEFAULT,noco,NTCO,NTIW
      CHARACTER CO(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER iw,noiw,N3CO
      CHARACTER CHAR2*2
      LOGICAL ABBREV,CBBREV

      CALL ENTERS('WS_LIST',*9999)
      CALL ASSERT(IWKDEF(0).GT.0,'  >>no workstations open',ERROR,*9999)
      IF(CBBREV(CO,'ON',1,noco+1,NTCO,N3CO)) THEN
        IF(ABBREV(CO(N3CO+1),'ENCAPSULATED',1).OR.
     1     ABBREV(CO(N3CO+1),'POSTSCRIPT',1)) THEN
          NTIW=1
C CPB 10/9/92 Changed default postscript workstation to be njt dependent
          IF(NJT.LE.2) THEN
            IWK(1)=15
          ELSE
            IWK(1)=16
          ENDIF
        ELSE !parse list of workstations in command
          CALL PARSIL(CO(N3CO+1),6,NTIW,IWK,ERROR,*9999)
        ENDIF
        DO noiw=1,NTIW !check specified workstations are defined
          iw=IWK(noiw)
          WRITE(CHAR2,'(I2)') iw
          CALL ASSERT(IWKS(iw).GE.1.AND.IWKS(iw).LE.4,
     '      ' >>Workstation '//CHAR2//' is not defined ',ERROR,*9999)
        ENDDO

      ELSE IF(CBBREV(CO,'POST',1,noco+1,NTCO,N3CO)) THEN
        POSTSCRIPT=.TRUE.
        !parse list of workstations in command
        CALL PARSIL(CO(N3CO+1),6,NTIW,IWK,ERROR,*9999)
        DO noiw=1,NTIW !check specified workstations are defined
          iw=IWK(noiw)
          WRITE(CHAR2,'(I2)') iw
          CALL ASSERT(IWKS(iw).GE.1.AND.IWKS(iw).LE.3,
     '      ' >>Workstation '//CHAR2//' is not defined ',ERROR,*9999)
        ENDDO

      ELSE !use list of defined workstations if IW_DEFAULT = 0
        IF(IW_DEFAULT.EQ.0) THEN
          NTIW=0
          DO noiw=1,IWKDEF(0)
            iw=IWKDEF(noiw)
C cpb 25/2/97 Restrict list to the first 3 workstations
            IF(IWKG(iw).EQ.1.AND.iw.GE.1.AND.iw.LE.3) THEN
              NTIW=NTIW+1
              IWK(NTIW)=iw
            ENDIF
          ENDDO
        ELSE                                !or default number
          NTIW=1
          IWK(1)=IW_DEFAULT
          WRITE(CHAR2,'(I2)') IWK(1)
          CALL ASSERT(IWKS(IWK(1)).GE.1.AND.IWKS(IWK(1)).LE.3,
     '      ' >>Workstation '//CHAR2//' is not defined ',ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('WS_LIST')
      RETURN
 9999 CALL ERRORS('WS_LIST',ERROR)
      CALL EXITS('WS_LIST')
      RETURN 1
      END




