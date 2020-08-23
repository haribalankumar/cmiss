      SUBROUTINE OPEN_FILE(ALL_REGIONS,IPFILE,FILE,IPEXTEN,STATUS,
     '  ERR,ERROR,*)

C#### Subroutine: OPEN_FILE
C###  Description:
C###    OPEN_FILE opens FILE.ip**** (if .not.ALL) or FILE.ir****
C###    file (if ALL) and reads/writes version number and heading.
C**** IPFILE is current version #.
C**** If writing file, this version # is printed to file.
C**** If reading file, IP_FILE is read and compared to IPFILE.
C   CS 30/1/2000 new
C**** ERR is a return code indicating success or failure types
C**** where a negative number is a failure type and a positive
C**** is a success type.
C**** The only error code implemented so far is -1 to indicate
C**** a mismatch in file versions


      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'head00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'mxch.inc'
!     Parameter List
      INTEGER ERR,IPFILE
      CHARACTER FILE*(*),IPEXTEN*(*),STATUS*3,ERROR*(*)
      LOGICAL ALL_REGIONS
!     Local Variables
      INTEGER IBEG,IBEG2,IBEG4,ICHAR,IEND,IEND2,IEND4,INFO,
     '  IOSTAT,IOTYPTEMP,IP_FILE,NOQUES
      CHARACTER CHAR*3,FILENAME*(MXCH),FMT_LOCAL*(MXCH),ISTAT*3
      LOGICAL FILEIP

      CALL ENTERS('OPEN_FILE',*9999)
      ERR=0
      NOQUES=0
      FILEIP=.FALSE.
      ICHAR=999
      CALL STRING_TRIM(FILE,IBEG,IEND)
      CALL STRING_TRIM(IPEXTEN,IBEG2,IEND2)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        FILENAME=FILE(IBEG:IEND)//'.ip'//IPEXTEN(IBEG2:IEND2)
        CALL STRING_TRIM(FILENAME,IBEG,IEND)
        WRITE(OP_STRING,'('' IOTYPE='',I2,'' IFILE='',I3,'
     '    //''' IPFILE='',I3,'' FILE='',A)')
     '    IOTYPE,IFILE,IPFILE,FILENAME(IBEG:IEND)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      IF(ALL_REGIONS) THEN !Write/read for every region
        IF(IOTYPE.NE.5) THEN !open file when not default option
          CALL STRING_TRIM(FILE,IBEG,IEND)
          FILENAME=FILE(IBEG:IEND)//'.ir'//IPEXTEN(IBEG2:IEND2)
          CALL STRING_TRIM(FILENAME,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILENAME(IBEG:IEND),
     '      STATUS,'DIRECT','FORMATTED',132,ERROR,*9999)
        ENDIF

        IF(IOTYPE.EQ.1.OR.IOTYPE.EQ.3) THEN
          FMT_LOCAL='('' CMISS Version '//CMISS//
     '      ' ir'//IPEXTEN(IBEG2:IEND2)//' File Version '',I1)'
          WRITE(UNIT=IFILE,REC=1,FMT=FMT_LOCAL,IOSTAT=IOSTAT,ERR=9998)
     '      IPFILE
          CALL STRING_TRIM(HEADING,IBEG4,IEND4)
          FMT_LOCAL='('' Heading: '//HEADING(IBEG4:IEND4)//''')'
          WRITE(UNIT=IFILE,REC=2,FMT=FMT_LOCAL,IOSTAT=IOSTAT,ERR=9998)
          FMT_LOCAL='(X)'
          WRITE(UNIT=IFILE,REC=3,FMT=FMT_LOCAL,IOSTAT=IOSTAT,ERR=9998)
        ELSE IF(IOTYPE.EQ.2.OR.IOTYPE.EQ.4) THEN
          WRITE(CHAR,'(I3)') 35+(IEND2-IBEG2)
C          CHAR=CFROMI(35+(IEND2-IBEG2),'(I3)')
          READ(UNIT=IFILE,REC=1,FMT='(1X,'//CHAR//'X,I2)',IOSTAT=IOSTAT,
     '      ERR=9998) IP_FILE
          IF(IP_FILE.NE.IPFILE) THEN
            ERR=-1
            WRITE(OP_STRING,'('' Warning!!! May need to redo file '
     '        //'input:'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''  Old file version is '',I3)') IP_FILE
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''  New file version is '',I3)') IPFILE
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          READ(UNIT=IFILE,REC=2,FMT='(1X,9X,A)',IOSTAT=IOSTAT,ERR=9998)
     '      HEADING
          WRITE(OP_STRING,'('' File heading: '',A)') HEADING
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          READ(UNIT=IFILE,REC=3,FMT='(1X)',IOSTAT=IOSTAT,ERR=9998)
        ENDIF

        WRITE(CHAR,'(I3)') NRT
C        CHAR=CFROMI(NRT,'(I3)')
        IDATA(1)=NRT
        IOTYPTEMP=3
        IF(IOTYPE.EQ.1) THEN !prompted input
          OP_STRING(1)=' '
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          OP_STRING(1)=' Total number of regions :'//CHAR(1:3)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          FORMAT='($,'' Total number of regions :'',I3)'
          CALL GINOUT(IOTYPTEMP,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NRM,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        ELSE
          FORMAT='($,'' Total number of regions : '',I3)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NRM,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        ENDIF

      ELSE IF(.NOT.ALL_REGIONS) THEN !Write/read for one region
        IF(IOTYPE.NE.5) THEN !open file when not default option
          CALL STRING_TRIM(FILE,IBEG,IEND)
          FILENAME=FILE(IBEG:IEND)//'.ip'//IPEXTEN(IBEG2:IEND2)
          CALL STRING_TRIM(FILENAME,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILENAME(IBEG:IEND),STATUS,
     '      'DIRECT','FORMATTED',132,ERROR,*9999)
        ENDIF
        IF(IOTYPE.EQ.1.OR.IOTYPE.EQ.3) THEN
          FMT_LOCAL='('' CMISS Version '//CMISS//
     '      ' ip'//IPEXTEN(IBEG2:IEND2)//' File Version '',I1)'
          WRITE(UNIT=IFILE,REC=1,FMT=FMT_LOCAL,
     '      IOSTAT=IOSTAT,ERR=9998) IPFILE
          CALL STRING_TRIM(HEADING,IBEG4,IEND4)
          WRITE(UNIT=IFILE,REC=2,FMT='('' Heading: '
     '      //HEADING(IBEG4:IEND4)//''')',IOSTAT=IOSTAT,ERR=9998)
          WRITE(UNIT=IFILE,REC=3,FMT='(X)',IOSTAT=IOSTAT,ERR=9998)
        ELSE IF(IOTYPE.EQ.2.OR.IOTYPE.EQ.4) THEN
          WRITE(CHAR,'(I3)') 35+(IEND2-IBEG2)
C          CHAR=CFROMI(35+(IEND2-IBEG2),'(I3)')
          READ(UNIT=IFILE,REC=1,FMT='(1X,'//CHAR//'X,I2)',IOSTAT=IOSTAT,
     '      ERR=9998) IP_FILE
          IF(IP_FILE.NE.IPFILE) THEN
            ERR=-1
            WRITE(OP_STRING,'('' Warning!!! May need to redo file '
     '        //'input:'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''  Old file version is '',I3)') IP_FILE
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''  New file version is '',I3)') IPFILE
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          READ(UNIT=IFILE,REC=2,FMT='(1X,9X,A)',IOSTAT=IOSTAT,ERR=9998)
     '      HEADING
          WRITE(OP_STRING,'('' File heading: '',A)') HEADING
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          READ(UNIT=IFILE,REC=3,FMT='(1X)',IOSTAT=IOSTAT,ERR=9998)
        ENDIF
      ENDIF

      CALL EXITS('OPEN_FILE')
      RETURN
 9998 WRITE(ISTAT,'(I3)') IOSTAT
C 9998 ISTAT=CFROMI(IOSTAT,'(I3)')
      ERROR='>> Error : Iostat = '//ISTAT
 9999 CALL ERRORS('OPEN_FILE',ERROR)
      CALL EXITS('OPEN_FILE')
      RETURN 1
      END


