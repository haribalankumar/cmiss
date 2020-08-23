      SUBROUTINE IMSIGN(LD,NBJ,NDDATA,NRLIST,NRLIST2,WD,XID,ZD,ZD2,
     '  STRING,ERROR,*)


C#### Subroutine: IMSIGN
C###  Description:
C###    IMSIGN imports signals into the FE data base.

C MPN 17JUL2002 NOTE: Generalised region handling, so changed nr-type
C                     indices in many places in the routine.

C!!!! WARNING!!!!   MPN 23 JULY 2002
C     THIS ROUTINE DOES NOT WORK CORRECTLY FOR OPTIMISED VERSIONS OF CM

C LKC 2-JUN-1998 NOTE: Changed several sections due to incorrect
C                      reading in of information


      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b14.cmn'
      INCLUDE 'binf00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'emap00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'impo00.cmn'
      INCLUDE 'mach00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'sign00.cmn'
!     Parameter List
      INTEGER LD(NDM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),NRLIST(0:NRM),
     '  NRLIST2(0:NRM)
      REAL*8 WD(NJM,NDM),XID(NIM,NDM),ZD(NJM,NDM),ZD2(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER CERROR(50),CERRLEN,device,DUMMY_DEVICE_TYPE(1),
     '  DUMMY_NUM_DEVICES(9),DUMMY_RIG_TYPE(1),ERR,IBEG,IBEG1,
     '  IEND,IEND1,INPUTTYPE,N3CO,naccept,ne,nd,nj,nodata,nonr_signal,
     '  nr,nreject,nr_signal,NUMSKIPBYTES,NUMSKIPBYTES1,
     '  NUMTIMEDATA,NUM_PAGES(2),NUM_REGIONS(2),
     '  NUM_RIGNAMEBYTES(2),NUM_SAMPLES(2),NUM_SIGNALS(2),
     '  page,POSITION,sample,sigdevice
      REAL*8 DATAN_MOD,PHI,PHIVECT(3),SIGNALMAX(9),SIGNALMIN(9),SUM1,
     '  SUM2,TIME,THETA,THETAVECT(3),TRANSMAT(3,4),X(3),Y(3),
     '  XVECT(3),YVECT(3),ZVECT(3)
      CHARACTER EMAPFNAME*100,EXTENSION*20,
     '  INSIGFNAME*100,OUTSIGFNAME*100
      LOGICAL ABBREV,ACTIVATION,ALL_REGIONS,ALL_SIG_REGIONS,CBBREV,
     '  EMAP,INLIST,ISBINFILEOPEN,ISENDBINFILE,ENDFILE,UPDATE_LOCATIONS,
     '  UPDATE_PROJECTIONS

      CALL ENTERS('IMSIGN',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM import signal<;FILENAME[default]>
C###  Parameter:      <signal_output FILENAME[file]>
C###    Specify the filename for the output cm signal file
C###  Parameter:      <extension EXTENSION[signal]>
C###    Specify the filename extension on the input signal file
C###  Parameter:      <update_locations>
C###    Specify that the electrode locations in the input signal file
C###    will be updated to be the same as the current electrodes
C###    locations
C###  Parameter:      <update_projections>
C###    Specify that the electrode projections in the input signal file
C###    will be updated to be the same as the current electrode
C###    projections
C###  Parameter:      <signal_region (#s/all)[all]>
C###    Specify the signal regions from which to import.
C###  Parameter:      <region (#s/all)[all]>
C###    Specify the regions to which imported signals should be stored.
C###  Description:
C###    Imports signals into the FE data base i.e., convert an input
C###    signal file into a cm signal file.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME['//FILE00(IBEG1:IEND1)
     '    //']>'
        OP_STRING(2)=BLANK(1:15)//'<signal_output FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(3)=BLANK(1:15)//'<extension EXTENSION[signal]>'
        OP_STRING(4)=BLANK(1:15)//'<update_locations>'
        OP_STRING(5)=BLANK(1:15)//'<update_projections>'
        OP_STRING(6)=BLANK(1:15)//'<signal_region (#s/all)[1]>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM import signal<;FILENAME[default]> activation
C###  Parameter:      <signal_output FILENAME[file]>
C###    Specify the filename for the output cm signal file
C###  Parameter:      <extension EXTENSION[signal]>
C###    Specify the filename extension on the input signal file
C###  Description:
C###    Imports signal events (activation times) into the FE data base
C###    i.e., converts an input signal events into a cm signal file.

      OP_STRING(1)=STRING(1:IEND)//'<;FILENAME['//FILE00(IBEG1:IEND1)
     '    //'] activation>'
      OP_STRING(2)=BLANK(1:15)//'<signal_output FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
      OP_STRING(3)=BLANK(1:15)//'<extension EXTENSION[signal]>'
      CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','IMSIGN',ERROR,*9999)
      ELSE

        CALL ASSERT(CALL_IMPO,'>>Define import first',ERROR,*9999)

C LKC 10-JUN-2002 Check memory is present for reading in signal
        CALL ASSERT(NJT+1.LE.NJM,'>>Increase NJM to NJT+1',ERROR,*9999)

C MPN 22JUL2002: Generalise region handling
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        IF(CBBREV(CO,'SIGNAL_REGION',8,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'ALL',1)) THEN
            ALL_SIG_REGIONS=.TRUE.
            NRLIST2(0)=0
          ELSE
            ALL_SIG_REGIONS=.FALSE.
            CALL PARSIL(CO(N3CO+1),NRM,NRLIST2(0),NRLIST2(1),
     '        ERROR,*9999)
          ENDIF
        ELSE
C LKC 11-JAN-2011
C!!!! This part desn't look right ... should default not be 1 region??
          
          ALL_SIG_REGIONS=.TRUE.
          NRLIST2(0)=0
        ENDIF
C OLD
CC cpb 10/2/00 Do it this way instead of using PARSE_REGIONS as we
CC are talking about signal regions and not cm regions
C        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
C          IF(ABBREV(CO(N3CO+1),'ALL',1)) THEN
C            ALL_SIG_REGIONS=.TRUE.
C            NRLIST(0)=0
C          ELSE
C            ALL_SIG_REGIONS=.FALSE.
C            CALL PARSIL(CO(N3CO+1),NRM,NRLIST(0),NRLIST(1),ERROR,*9999)
C          ENDIF
C        ELSE
C          ALL_SIG_REGIONS=.TRUE.
C          NRLIST(0)=0
C        ENDIF

        IF(CBBREV(CO,'SIGNAL_OUTPUT',3,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          OUTSIGFNAME=CO(N3CO+1)(IBEG:IEND)
        ELSE
          CALL STRING_TRIM(FILE00,IBEG,IEND)
          OUTSIGFNAME=FILE00(IBEG:IEND)
        ENDIF

        CALL CHECKF(1,noco,NTCOQU,CO,COQU,INSIGFNAME,STRING,*1)

        IF(SIGIMPO_TYPE.EQ.1) THEN
          EMAP=.TRUE.
        ELSE
          ERROR='>>Unknown signal type'
          GOTO 9999
        ENDIF

C LKC 18-AUG-1998 Can specify UNEMAP file extension
        IF(CBBREV(CO,'EXTENSION',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
          EXTENSION=CO(N3CO+1)(IBEG1:IEND1)
        ELSE
          EXTENSION='signal'
        ENDIF

C cpb 10/02/00 Adding update locations
        IF(CBBREV(CO,'UPDATE_LOCATIONS',8,noco+1,NTCO,N3CO)) THEN
          UPDATE_LOCATIONS=.TRUE.
        ELSE
          UPDATE_LOCATIONS=.FALSE.
        ENDIF
C cpb 10/02/00 Adding update projections
        IF(CBBREV(CO,'UPDATE_PROJECTIONS',8,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(CALC_XI,'>>Calculate Xi first',ERROR,*9999)
          UPDATE_PROJECTIONS=.TRUE.
        ELSE
          UPDATE_PROJECTIONS=.FALSE.
        ENDIF

C JMB 15-AUG-2000 Adding activation events
        IF(CBBREV(CO,'ACTIVATION',2,noco+1,NTCO,N3CO)) THEN
          ACTIVATION=.TRUE.
        ELSE
          ACTIVATION=.FALSE.
        ENDIF

C LKC 4-APR-1998 Initialise time
        TIME=0.D0

        IF(EMAP) THEN

c cpb 28/10/95 An emap signal file consists of three main types of
C data, namely configuration data, signal data, and then analysis data.
C The configuration data contains information on the recording rig,
C the regions involved and the devices within the regions. The signal
C data contains the actual data and information from recording the data
C (e.g. recording frequency). The analysis data conatins post recording
C information about the signals (e.g. bad leads).
C
C The bad lead information in a cmiss signal file is at the beginning
C of the file. Hence the procedure is read the configuration
C information then skip the signal data and read the analysis
C information (if any) write the first two tags of the cmiss signal
C file and then back up to read the signal data and write it out to
C the cmiss signal file instant by instant. As emap files contain no
C information about the size of the quantities stored in the file it
C is assumed that they are the same size as the native sizes.

          FILEMACHTYPE(IOFILE1)=MACHTYPE
          FILEOSTYPE(IOFILE1)=OSTYPE
          FILEENDIANTYPE(IOFILE1)=CHAR(MACH_BIGENDIAN)
          FILESPFORMTYPE(IOFILE1)=SPFORMTYPE
          FILEDPFORMTYPE(IOFILE1)=DPFORMTYPE
          FILECHARSIZE(IOFILE1)=CHAR(CHARSIZE)
          FILEINTSIZE(IOFILE1)=CHAR(INTSIZE)
          FILESINTSIZE(IOFILE1)=CHAR(SINTSIZE)
          FILESPSIZE(IOFILE1)=CHAR(SPSIZE)
          FILEDPSIZE(IOFILE1)=CHAR(DPSIZE)
          FILELOGSIZE(IOFILE1)=CHAR(LOGSIZE)

C***      Open up the emap signal file

          CALL STRING_TRIM(INSIGFNAME,IBEG,IEND)
          CALL STRING_TRIM(EXTENSION,IBEG1,IEND1)
          EMAPFNAME=INSIGFNAME(IBEG:IEND)//'.'//EXTENSION(IBEG1:IEND1)
          CALL BINOPENFILE(IOFILE1,'READ',EMAPFNAME,ERROR,*9999)

C***      Read in the configuration

C***      Read in the global rig type
          CALL BINREADFILE(IOFILE1,INTTYPE,1,EMAP_RIGTYPE(0),REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

C***      Read in the rig name
C LKC 2-JUN-1998 swap around
C          NUM_RIGNAMEBYTES(1)=EMAP_NUMRIGNAMEBYTES
C          CALL BINREADFILE(IOFILE1,INTTYPE,1,NUM_RIGNAMEBYTES,
C     '      REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

          CALL BINREADFILE(IOFILE1,INTTYPE,1,NUM_RIGNAMEBYTES,
     '      REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          EMAP_NUMRIGNAMEBYTES=NUM_RIGNAMEBYTES(1)
          CALL BINREADFILE(IOFILE1,CHARTYPE,EMAP_NUMRIGNAMEBYTES,
     '      INTDATA,REAL4DATA,REAL8DATA,EMAP_RIGNAME,LOGDATA,SINTDATA,
     '      ERROR,*9999)
          NUMSKIPBYTES=INTSIZE+INTSIZE+EMAP_NUMRIGNAMEBYTES*CHARSIZE
          IF(EMAP_RIGTYPE(0).EQ.EMAP_SOCK) THEN
C***        Read in the sock focus
            CALL BINREADFILE(IOFILE1,SPTYPE,1,INTDATA,REAL4DATA,
     '        REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
            EMAP_SOCKFOCUS(1)=DBLE(REAL4DATA(1))
            NUMSKIPBYTES=NUMSKIPBYTES+SPSIZE
          ENDIF

C***      Read in the number of regions
C LKC 2-JUN-1998 swap around
C          NUM_REGIONS(1)=EMAP_NUMREGIONS
C          CALL BINREADFILE(IOFILE1,INTTYPE,1,NUM_REGIONS,
C     '      REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C          NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE
          CALL BINREADFILE(IOFILE1,INTTYPE,1,NUM_REGIONS,
     '      REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE
          EMAP_NUMREGIONS=NUM_REGIONS(1)

          IF(EMAP_NUMREGIONS.LE.0) THEN
            WRITE(OP_STRING,'(''>>WARNING: EMAP_NUMREGIONS <= 0 '')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
CC JMBs 15-AUG-2000 Check that only one region exist for the activation
CC times.
            IF(ACTIVATION) THEN
              CALL ASSERT(EMAP_NUMREGIONS.EQ.1,
     '          '>>Only one region valid for activation times',ERROR,
     '          *9999)
            ENDIF
CC JMBe
C LKC 12-JUN-1998
          WRITE(OP_STRING(1),'('' '')')
          WRITE(OP_STRING(2),'('' Importing from  '',A50)') INSIGFNAME
          WRITE(OP_STRING(3),'('' Total Unemap regions    '',I6)')
     '      EMAP_NUMREGIONS
          WRITE(OP_STRING(4),'('' Importing:'')')
          CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C MPN 22JUL2002: Generalise region handling
          IF(ALL_REGIONS.AND.ALL_SIG_REGIONS) THEN
C LKC 7-JUL-2004
C  There is no reason why this needs to be enforced. Should be able to
C  import a multiregion signal without any regions having being defined first
C            
C            CALL ASSERT(EMAP_NUMREGIONS.LE.NRT,
C     '        '>>Too many signal regions (> num cm regions)',
C     '        ERROR,*9999)

            CALL ASSERT(EMAP_NUMREGIONS.LE.NRM,'>>Increase NRM',
     &        ERROR,*9999)
            
            NRLIST(0)=EMAP_NUMREGIONS
            NRLIST2(0)=EMAP_NUMREGIONS
            DO nr_signal=1,EMAP_NUMREGIONS
              NRLIST2(nr_signal)=nr_signal
            ENDDO !nr_signal
          ELSE IF(.NOT.ALL_REGIONS.AND.ALL_SIG_REGIONS) THEN
C LKC 7-JUL-2004
C  There is no reason why this needs to be enforced. Should be able to
C  import a multiregion signal without any regions having being defined first
C                        
C            CALL ASSERT(EMAP_NUMREGIONS.LE.NRT,
C     '        '>>Too many signal regions (> num cm regions)',
C     '        ERROR,*9999)
            NRLIST2(0)=NRLIST(0)
            DO nonr_signal=1,NRLIST(0)
              nr=NRLIST(nonr_signal)
              NRLIST2(nonr_signal)=nr
              CALL ASSERT(nr.LE.EMAP_NUMREGIONS,
     '          'specified cm region is not in the signal file',
     '          ERROR,*9999)
            ENDDO
          ELSE IF(.NOT.ALL_SIG_REGIONS) THEN
            IF(ALL_REGIONS) NRLIST(0)=NRLIST2(0)
            DO nonr_signal=1,NRLIST2(0)
              nr_signal=NRLIST2(nonr_signal)
              IF(ALL_REGIONS) NRLIST(nonr_signal)=nr_signal
              CALL ASSERT(nr_signal.GE.1.AND.
     '          nr_signal.LE.EMAP_NUMREGIONS,
     '          'specified signal region is not in the signal file',
     '          ERROR,*9999)

C LKC 7-JUL-2004
C  There is no reason why this needs to be enforced. Should be able to
C  import a multiregion signal without any regions having being defined first
C
C              CALL ASSERT(nr_signal.LE.NRT,
C     '          '>>specified signal region > number of cm regions',
C     '          ERROR,*9999)
C            ENDDO !nonr_signal            
C            CALL ASSERT(NRLIST(0).LE.NRT,
C     '        '>>number of signal regions > number of cm regions',
C     '        ERROR,*9999)
            ENDDO !nonr_signal            
          ENDIF
          
          CALL ASSERT(NRLIST(0).EQ.NRLIST2(0),
     '      '>>Incompatible signal/cm region lists',ERROR,*9999)
C OLD
c cpb 10/02/00 Adding read of specific signal regions
c          IF(.NOT.ALL_SIG_REGIONS) THEN
c            DO nonr=1,NRLIST(0)
c              nr=NRLIST(nonr)
c              CALL ASSERT(nr.GE.1.AND.nr.LE.EMAP_NUMREGIONS,
c     '          'Specified region is not in the signal file',
c     '          ERROR,*9999)
c            ENDDO !nornr
c          ENDIF

C cpb 10/02/00 Adding update_locations
          IF(UPDATE_LOCATIONS) THEN
            CALL ASSERT(NDT.GT.0,'No electrodes currently defined',
     '        ERROR,*9999)
C***        Copy current locations to ZD2 so ZD can be overwritten
            DO nd=1,NDT
              DO nj=1,NJT
                ZD2(nj,nd)=ZD(nj,nd)
              ENDDO !nj
            ENDDO !nd
          ENDIF

          nd=0
          DO nr_signal=1,EMAP_NUMREGIONS

C LKC 19-SEP-2000 fixing up this section - EMAP_RIGTYPE(nr)
C   is not initialised.
C
C cpb 10/02/00 Adding read of specific regions
C            IF(INLIST(nr_signal,NRLIST2(1),NRLIST2(0),POSITION)) THEN
C              nr=NRLIST(POSITION)
C              nodata=0
C              IF(EMAP_RIGTYPE(0).EQ.EMAP_MIXED) THEN
CC***            Read in the region rig type
C                CALL BINREADFILE(IOFILE1,INTTYPE,1,EMAP_RIGTYPE(nr_signal),
C     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
C     '            ERROR,*9999)
C                NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE
C              ELSE IF(EMAP_RIGTYPE(0).EQ.EMAP_SOCK) THEN
C                EMAP_SOCKFOCUS(nr_signal)=EMAP_SOCKFOCUS(1)
C              ENDIF

            IF(INLIST(nr_signal,NRLIST2(1),NRLIST2(0),POSITION)) THEN
              nr=NRLIST(POSITION)
              nodata=0
              IF(EMAP_RIGTYPE(0).EQ.EMAP_MIXED) THEN
C***            Read in the region rig type
                CALL BINREADFILE(IOFILE1,INTTYPE,1,
     '            EMAP_RIGTYPE(nr_signal),REAL4DATA,REAL8DATA,CHARDATA,
     '            LOGDATA,SINTDATA,ERROR,*9999)
                NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE
              ELSE IF(EMAP_RIGTYPE(0).EQ.EMAP_SOCK) THEN
                EMAP_SOCKFOCUS(nr_signal)=EMAP_SOCKFOCUS(1)
                EMAP_RIGTYPE(nr_signal)=EMAP_SOCK
              ELSE IF(EMAP_RIGTYPE(0).EQ.EMAP_TORSO) THEN
                EMAP_RIGTYPE(nr_signal)=EMAP_TORSO

C LKC 24-FEB-2003 EMAP_PATCH not initialised till now
              ELSE IF(EMAP_RIGTYPE(0).EQ.EMAP_PATCH) THEN
                EMAP_RIGTYPE(nr_signal)=EMAP_PATCH                
              ENDIF

C***          Read in the region name
              CALL BINREADFILE(IOFILE1,INTTYPE,1,INTDATA,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              CALL BINREADFILE(IOFILE1,CHARTYPE,INTDATA(1),INTDATA,
     '          REAL4DATA,REAL8DATA,EMAP_REGIONNAME(nr_signal),LOGDATA,
     '          SINTDATA,ERROR,*9999)
              NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE+INTDATA(1)*CHARSIZE
              IF(EMAP_RIGTYPE(0).EQ.EMAP_MIXED.AND.
     '          EMAP_RIGTYPE(nr_signal).EQ.EMAP_SOCK) THEN
C***            Read in the sock focus
                CALL BINREADFILE(IOFILE1,SPTYPE,1,INTDATA,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                EMAP_SOCKFOCUS(nr_signal)=DBLE(REAL4DATA(1))
                NUMSKIPBYTES=NUMSKIPBYTES+SPSIZE
              ENDIF
C***          Read in the number of devices
              CALL BINREADFILE(IOFILE1,INTTYPE,1,
     '          EMAP_NUMDEVICES(nr_signal),REAL4DATA,REAL8DATA,
     '          CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              CALL ASSERT(EMAP_NUMDEVICES(nr_signal).LE.NEMAPDEVICESMX,
     '          '>>Increase NEMAPDEVICESMX in emap00.cmn',ERROR,*9999)
              NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE

C LKC 12-JUN-1998
C            WRITE(OP_STRING(1),'(''   Region        '',I12)') nr_signal
C            WRITE(OP_STRING(2),'(''     Total Channels'',I12)')
C     '        EMAP_NUMDEVICES(nr_signal)
C            CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

              DO device=1,EMAP_NUMDEVICES(nr_signal)
C***            Read the device type
                CALL BINREADFILE(IOFILE1,INTTYPE,1,
     '            EMAP_DEVICETYPE(device,nr_signal),REAL4DATA,REAL8DATA,
     '            CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C***            Read the device number
                CALL BINREADFILE(IOFILE1,INTTYPE,1,INTDATA,
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)
C***            Read the device name
                CALL BINREADFILE(IOFILE1,INTTYPE,1,INTDATA,
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)
                CALL BINREADFILE(IOFILE1,CHARTYPE,INTDATA(1),INTDATA,
     '            REAL4DATA,REAL8DATA,EMAP_DEVICENAME(device,nr_signal),
     '            LOGDATA,SINTDATA,ERROR,*9999)
                NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE+INTSIZE+INTSIZE+
     '            INTDATA(1)*CHARSIZE
C***            Read the channel number (not needed)
                CALL BINREADFILE(IOFILE1,INTTYPE,1,
     '            EMAP_CHANNELNUM(device,nr_signal),
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)
                NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE
C***            Read device dependent properties
                IF(EMAP_DEVICETYPE(device,nr_signal).EQ.EMAP_ELECTRODE)
     '            THEN
C***              Is an electrode so add the device to the list of
C***              CMISS electrodes
                  nd=nd+1
                  IF(nd.LE.NDM) THEN
                    nodata=nodata+1
                    NDDATA(0,nr)=nodata
                    NDDATA(nodata,nr)=nd

C LKC 7-SEP-1999 Changing NDT
C                  IF(nd.GT.NDT) NDT=nd
                    IF(nd.GT.NDT) NDT=nodata

                  ENDIF
                  
C***              Read position
                  IF(EMAP_RIGTYPE(nr_signal).EQ.EMAP_SOCK.OR.
     '              EMAP_RIGTYPE(nr_signal).EQ.EMAP_TORSO) THEN
                    CALL ASSERT(NJT.EQ.3,
     '                '>>NJT must be 3 for sock or torso electrodes',
     '                ERROR,*9999)
                    CALL BINREADFILE(IOFILE1,SPTYPE,3,INTDATA,
     '                REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '                ERROR,*9999)
                    IF(nd.LE.NDM) THEN
                      IF(UPDATE_LOCATIONS) THEN
                        ZD(1,nd)=ZD2(1,nd)
                        ZD(2,nd)=ZD2(2,nd)
                        ZD(3,nd)=ZD2(3,nd)
                      ELSE
                        ZD(1,nd)=DBLE(REAL4DATA(1))
                        ZD(2,nd)=DBLE(REAL4DATA(2))
                        ZD(3,nd)=DBLE(REAL4DATA(3))
                      ENDIF
                      WD(1,nd)=1.0d0
                      WD(2,nd)=1.0d0
                      WD(3,nd)=1.0d0
                    ENDIF
                    NUMSKIPBYTES=NUMSKIPBYTES+3*SPSIZE
                  ELSE IF(EMAP_RIGTYPE(nr_signal).EQ.EMAP_PATCH) THEN
                    CALL ASSERT(NJT.EQ.2,
     '                '>>NJT must be 2 for patch electrodes',
     '                ERROR,*9999)
                    CALL BINREADFILE(IOFILE1,SPTYPE,2,INTDATA,
     '                REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '                ERROR,*9999)
                    IF(nd.LE.NDM) THEN
                      IF(UPDATE_LOCATIONS) THEN
                        ZD(1,nd)=ZD2(1,nd)
                        ZD(2,nd)=ZD2(2,nd)
                      ELSE
                        ZD(1,nd)=DBLE(REAL4DATA(1))
                        ZD(2,nd)=DBLE(REAL4DATA(2))
                      ENDIF
                      WD(1,nd)=1.0d0
                      WD(2,nd)=1.0d0
                    ENDIF
                    NUMSKIPBYTES=NUMSKIPBYTES+2*SPSIZE
                  ENDIF
                ENDIF
              ENDDO !ndev
            ELSE
              IF(EMAP_RIGTYPE(0).EQ.EMAP_MIXED) THEN
C***            Skip the the region rig type
                CALL BINREADFILE(IOFILE1,INTTYPE,1,DUMMY_RIG_TYPE(1),
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)
                NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE
              ENDIF
C***          Skip the region name
              CALL BINREADFILE(IOFILE1,INTTYPE,1,INTDATA,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              NUMSKIPBYTES1=INTDATA(1)*CHARSIZE
              CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
              NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE+INTDATA(1)*CHARSIZE
              IF(EMAP_RIGTYPE(0).EQ.EMAP_MIXED.AND.
     '          DUMMY_RIG_TYPE(1).EQ.EMAP_SOCK) THEN
C***            Skip the sock focus
                NUMSKIPBYTES1=SPSIZE
                CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
                NUMSKIPBYTES=NUMSKIPBYTES+SPSIZE
              ENDIF
C***          Skip the number of devices
              CALL BINREADFILE(IOFILE1,INTTYPE,1,
     '          DUMMY_NUM_DEVICES(nr_signal),REAL4DATA,REAL8DATA,
     '          CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE

              DO device=1,DUMMY_NUM_DEVICES(nr_signal)
C***            Skip the device type
                CALL BINREADFILE(IOFILE1,INTTYPE,1,DUMMY_DEVICE_TYPE(1),
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)
C***            Skip the device number
                NUMSKIPBYTES1=INTSIZE
                CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
C***            Skip the device name
                CALL BINREADFILE(IOFILE1,INTTYPE,1,INTDATA,
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)
                NUMSKIPBYTES1=INTDATA(1)*CHARSIZE
                CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
                NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE+INTSIZE+INTSIZE+
     '            INTDATA(1)*CHARSIZE
C***            Skip the channel number (not needed)
                NUMSKIPBYTES1=INTSIZE
                CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
                NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE
C***            Read device dependent properties
                IF(DUMMY_DEVICE_TYPE(1).EQ.EMAP_ELECTRODE) THEN
C***              Skip position
                  IF(DUMMY_RIG_TYPE(1).EQ.EMAP_SOCK.OR.
     '              DUMMY_RIG_TYPE(1).EQ.EMAP_TORSO) THEN
                    NUMSKIPBYTES1=3*SPSIZE
                    CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
                    NUMSKIPBYTES=NUMSKIPBYTES+3*SPSIZE
                  ELSE IF(DUMMY_RIG_TYPE(1).EQ.EMAP_PATCH) THEN
                    NUMSKIPBYTES1=2*SPSIZE
                    CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
                    NUMSKIPBYTES=NUMSKIPBYTES+2*SPSIZE
                  ENDIF
                ENDIF
              ENDDO !ndev
            ENDIF
          ENDDO !nr_signal

C***      Read the number of pages
C LKC 2-JUN-1998 swap around
C          NUM_PAGES(1)=EMAP_NUMPAGES
          CALL BINREADFILE(IOFILE1,INTTYPE,1,NUM_PAGES(1),
     '      REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '      *9999)
           EMAP_NUMPAGES=NUM_PAGES(1)

          NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE
C cpb 29/10/95 Pages are just display info so just skip over the data
          DO page=1,EMAP_NUMPAGES
C***        Read the page name
            CALL BINREADFILE(IOFILE1,INTTYPE,1,INTDATA,
     '        REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '        *9999)
            NUMSKIPBYTES1=INTDATA(1)*CHARSIZE
            CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
            NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE+NUMSKIPBYTES1
C***        Read the number of devices in the device list
            CALL BINREADFILE(IOFILE1,INTTYPE,1,INTDATA,
     '        REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '        *9999)
C***        Read the device numbers
            NUMSKIPBYTES1=INTDATA(1)*INTSIZE
            CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
            NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE+NUMSKIPBYTES1
          ENDDO !page

          NDT=nd
          CALL ASSERT(NDT.LE.NDM,'>>Increase NDM',ERROR,*9999)

C***      Read the signal data

C***      Read the number of signals
C LKC 2-JUN-1998 swap around
C          NUM_SIGNALS(1)=EMAP_NUMSIGNALS
C          CALL BINREADFILE(IOFILE1,INTTYPE,1,NUM_SIGNALS,
C     '      REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
C     '      *9999)
          CALL BINREADFILE(IOFILE1,INTTYPE,1,NUM_SIGNALS,
     '      REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '      *9999)
          EMAP_NUMSIGNALS=NUM_SIGNALS(1)

C??? LKC do we need an absolute qualifier?
          CALL ASSERT(EMAP_NUMSIGNALS.LE.NEMAPSIGNALSMX,
     '      '>>Increase NEMAPSIGNALSMX in emap00.cmn',ERROR,*9999)
          IF(EMAP_NUMSIGNALS.LT.0) THEN
            INPUTTYPE=EMAP_FLOAT_VALUE
            EMAP_NUMSIGNALS=-EMAP_NUMSIGNALS
          ELSE
            INPUTTYPE=EMAP_SHORT_INT_VALUE
          ENDIF
C***      Read in the number of samples
C          NUM_SAMPLES(1)=EMAP_NUMSAMPLES
C          CALL BINREADFILE(IOFILE1,INTTYPE,1,NUM_SAMPLES,
C     '      REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
C     '      *9999)
          CALL BINREADFILE(IOFILE1,INTTYPE,1,NUM_SAMPLES,
     '      REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '      *9999)
          EMAP_NUMSAMPLES=NUM_SAMPLES(1)

C***      Read in the sampling frequency
          CALL BINREADFILE(IOFILE1,SPTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          EMAP_FREQUENCY=DBLE(REAL4DATA(1))
          CALL ASSERT(EMAP_FREQUENCY.GT.0.0d0,
     '      '>>Invalid frequency (< 0)',ERROR,*9999)

C***      Read in the sample times
C cpb 29/10/95 These aren't needed at the moment so just skip them
          NUMSKIPBYTES1=EMAP_NUMSAMPLES*INTSIZE
          CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
          NUMSKIPBYTES=NUMSKIPBYTES+INTSIZE+INTSIZE+SPSIZE+
     '      NUMSKIPBYTES1

C***      Skip over the signal data for now (will read it later)

          IF(INPUTTYPE.EQ.EMAP_FLOAT_VALUE) THEN
            NUMSKIPBYTES1=EMAP_NUMSIGNALS*EMAP_NUMSAMPLES*SPSIZE
          ELSE IF(INPUTTYPE.EQ.EMAP_SHORT_INT_VALUE) THEN
            NUMSKIPBYTES1=EMAP_NUMSIGNALS*EMAP_NUMSAMPLES*SINTSIZE
          ENDIF
          CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)

          DO nr_signal=1,EMAP_NUMREGIONS
C cpb 10/02/00 Adding read of specific regions
            IF(INLIST(nr_signal,NRLIST2(1),NRLIST2(0),POSITION)) THEN
              DO device=1,EMAP_NUMDEVICES(nr_signal)
C***            Read in the device index
                CALL BINREADFILE(IOFILE1,INTTYPE,1,
     '            EMAP_SIGNALDEVICE(device,nr_signal),
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '            *9999)
C***            Read in the channel offset
                CALL BINREADFILE(IOFILE1,SPTYPE,1,INTDATA,
     '            EMAP_CHANNELOFFSET(device,nr_signal),REAL8DATA,
     '            CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C***            Read in the channel gain
                CALL BINREADFILE(IOFILE1,SPTYPE,1,INTDATA,
     '            EMAP_CHANNELGAIN(device,nr_signal),REAL8DATA,CHARDATA,
     '            LOGDATA,SINTDATA,ERROR,*9999)

C LKC uncomment EMAP_CHANNELGAIN 4-JUN-1998
CC cpb 12/11/95 Hard coding gain for the moment to be +/- 300mv fsd
C              EMAP_CHANNELGAIN(device,nr_signal)=0.3/REAL(32767)
              ENDDO !device
            ELSE
              DO device=1,DUMMY_NUM_DEVICES(nr_signal)
C***            Skip the device index, channel offset and gain
                NUMSKIPBYTES1=INTSIZE+2*SPSIZE
                CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
              ENDDO !device
            ENDIF
          ENDDO !nr_signal

C***      Check if there is any analysis data (ie. if not at eof)
          IF(.NOT.ISENDBINFILE(IOFILE1)) THEN
            IF(.NOT.ACTIVATION) THEN
C cpb 29/10/95 Only interested in the signal status at the moment so
C skip over analysis data.
              NUMSKIPBYTES1=12*INTSIZE+CHARSIZE
              CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
              DO nr_signal=1,EMAP_NUMREGIONS
C cpb 10/02/00 Adding read of specific regions
                IF(INLIST(nr_signal,NRLIST2(1),NRLIST2(0),POSITION))THEN
                  nr=NRLIST(POSITION)
                  nodata=0

C*** LKC 12-JUN-1998 Calculates number of accepted and rejected
                  nreject=0
                  naccept=0

                  DO device=1,EMAP_NUMDEVICES(nr_signal)
C***                Read the signal status
                    CALL BINREADFILE(IOFILE1,INTTYPE,1,
     '                EMAP_SIGNALTYPE(device,nr_signal),REAL4DATA,
     '                REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                    IF(EMAP_DEVICETYPE(device,nr_signal).EQ.
     '                EMAP_ELECTRODE) THEN
                      nodata=nodata+1
                      nd=NDDATA(nodata,nr)
                      IF(EMAP_SIGNALTYPE(device,nr_signal).EQ.
     '                  EMAP_REJECTED) THEN
                        WD(NJT+1,nd)=0.0d0
                        nreject=nreject+1
                      ELSE
                        WD(NJT+1,nd)=1.0d0
                        naccept=naccept+1
                      ENDIF
                    ENDIF
C***                Skip the signal min and max
                    NUMSKIPBYTES1=2*SPSIZE
                    CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
C***                Read the number of events
                    CALL BINREADFILE(IOFILE1,INTTYPE,1,INTDATA,
     '                REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '                ERROR,*9999)
C***                Skip the event data
                    NUMSKIPBYTES1=3*INTDATA(1)*INTSIZE
                    CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
                  ENDDO !device

C LKC 12-JUN-1998

                  WRITE(OP_STRING(1),'(''   - signal region '',I12)')
     '              nr_signal
                  WRITE(OP_STRING(2),'(''     into cm region'',I12)')
     '              nr
                  WRITE(OP_STRING(3),'(''     Number of Channels  '',
     '              I6)') EMAP_NUMDEVICES(nr_signal)
                  WRITE(OP_STRING(4),'(''     Accepted Channels   '',
     '              I6)') naccept
                  WRITE(OP_STRING(5),'(''     Rejected Channels   '',
     '              I6)') nreject
                  WRITE(OP_STRING(6),'(''     Frequency     '',F12.1)')
     '              EMAP_FREQUENCY
                  CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
                ELSE
                  DO device=1,DUMMY_NUM_DEVICES(nr_signal)
C***                Skip the signal status and signal min and max
                    NUMSKIPBYTES1=INTSIZE+2*SPSIZE
                    CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
C***                Skip the number of events
                    CALL BINREADFILE(IOFILE1,INTTYPE,1,INTDATA,
     '                REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '                ERROR,*9999)
C***                Skip the event data
                    NUMSKIPBYTES1=3*INTDATA(1)*INTSIZE
                    CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
                  ENDDO !device
                ENDIF
              ENDDO !nr_signal
            ELSE
CC JMB 15-AUG-2000 The event data contains the activation times from
CC maximum -dV/dt. Only one region must be present in the signal file.
              NUMSKIPBYTES1=12*INTSIZE+CHARSIZE
              CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
              nr=1
              nodata=0
              nr_signal=1 ! Default only one region
              DO device=1,EMAP_NUMDEVICES(nr_signal)
C***            Read the signal status
                CALL BINREADFILE(IOFILE1,INTTYPE,1,
     '            EMAP_SIGNALTYPE(device,nr_signal),REAL4DATA,REAL8DATA,
     '            CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                IF(EMAP_DEVICETYPE(device,nr_signal).EQ.EMAP_ELECTRODE)
     '            THEN
                  nodata=nodata+1
                  nd=NDDATA(nodata,nr)
                  IF(EMAP_SIGNALTYPE(device,nr_signal).EQ.EMAP_REJECTED)
     '              THEN
                    WD(NJT+1,nd)=0.0d0
                  ELSE
                    WD(NJT+1,nd)=1.0d0
                  ENDIF
                ENDIF
C***            Skip the signal min and max
                NUMSKIPBYTES1=2*SPSIZE
                CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
C***            Read the number of events
                CALL BINREADFILE(IOFILE1,INTTYPE,1,INTDATA,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                CALL ASSERT(INTDATA(1).EQ.1,
     '            '>>Only 1 event valid for activation times',
     '            ERROR,*9999)
C***            Read the events data
                CALL BINREADFILE(IOFILE1,INTTYPE,1,INTDATA,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                EMAP_FSIGBUFF(device)=FLOAT(INTDATA(1))
C***            Skip the remaining event data
                NUMSKIPBYTES1=2*INTSIZE
                CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES1,ERROR,*9999)
              ENDDO !device
            ENDIF
          ELSE
C***        No analysis data so assume all electrodes are good.
            DO nd=1,NDT
              WD(NJT+1,nd)=1.0d0
            ENDDO !nd
          ENDIF

C***      Write the CMISS signal file header

          SIGNAL_HOWGENERATED(IOFILE2)=1 !Generated by EMAP
          SIGNAL_HEADER(IOFILE2)='From emap file : '
     '      //INSIGFNAME(IBEG:IEND)//'.signal'
          IF(ALL_SIG_REGIONS) THEN
            SIGNAL_NUMREGIONS(IOFILE2)=EMAP_NUMREGIONS
          ELSE
            SIGNAL_NUMREGIONS(IOFILE2)=NRLIST(0)
          ENDIF
          DO nonr_signal=1,SIGNAL_NUMREGIONS(IOFILE2)
            nr=NRLIST(nonr_signal)
            nr_signal=NRLIST2(nonr_signal)
            SIGNAL_NUMELEC(nr,IOFILE2)=NDDATA(0,nr)
            SIGNAL_REGNAME(nr,IOFILE2)=EMAP_REGIONNAME(nr_signal)
            SIGNAL_REGTYPE(nr,IOFILE2)=0 !irregular rect. cartesian
          ENDDO !nonr_signal
          IF(UPDATE_PROJECTIONS) THEN
            SIGNAL_ELEMLOC(IOFILE2)=1
            DO nonr_signal=1,SIGNAL_NUMREGIONS(IOFILE2)
              nr=NRLIST(nonr_signal)
              ne=0
              nodata=1
              DO WHILE(ne.EQ.0.AND.nodata.LE.NDDATA(0,nr))
                nd=NDDATA(nodata,nr)
                ne=LD(nd)
                nodata=nodata+1
              ENDDO
              CALL ASSERT(ne.NE.0,
     '          '>>Could not find projection element',ERROR,*9999)
              SIGNAL_NUMXI(nr,IOFILE2)=NIT(NBJ(1,ne))
            ENDDO !nonr_signal
          ELSE
            SIGNAL_ELEMLOC(IOFILE2)=0
          ENDIF

          CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'WRITE','BINARY',OUTSIGFNAME,
     '      'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

C***      Transform the electrode positions into rectangular cartesian
C***      coordinates.
          DO nr_signal=1,EMAP_NUMREGIONS
C cpb 10/02/00 Adding read of specific regions
            IF(INLIST(nr_signal,NRLIST2(1),NRLIST2(0),POSITION)) THEN
              nr=NRLIST(POSITION)
              IF(EMAP_RIGTYPE(nr_signal).EQ.EMAP_SOCK) THEN
C***            Determine the sock transformation parameters

C***            First adjust the focus of the data
                nodata=0
                DO device=1,EMAP_NUMDEVICES(nr_signal)
                  IF(EMAP_DEVICETYPE(device,nr_signal).EQ.
     '              EMAP_ELECTRODE) THEN
                    nodata=nodata+1
                    nd=NDDATA(nodata,nr)
                    FOCUS=EMAP_SOCKFOCUS(nr_signal)
                    CALL ZX(4,ZD(1,nd),X)
                    FOCUS=EMAP_SOCKFOCUS(0)
                    CALL XZ(4,X,Y)
C cpb 9/11/95 Convert to body coordinates
                    ZD(1,nd)=-Y(2)
                    ZD(2,nd)=Y(3)
                    ZD(3,nd)=-Y(1)
                  ENDIF
                ENDDO !device

c cpb 24/7/02 Don't do this if you are updating positions
                IF(.NOT.UPDATE_LOCATIONS) THEN
C***              Rotate and translate the axes of the prolate
C***              Determine the theta, and phi rotation
                  THETAVECT(1)=EMAP_SOCKORIGIN(1)-EMAP_SOCKTHETAZERO(1)
                  THETAVECT(2)=EMAP_SOCKORIGIN(2)-EMAP_SOCKTHETAZERO(2)
                  THETAVECT(3)=EMAP_SOCKORIGIN(3)-EMAP_SOCKTHETAZERO(3)
                  PHIVECT(1)=EMAP_SOCKORIGIN(1)-EMAP_SOCKAPEX(1)
                  PHIVECT(2)=EMAP_SOCKORIGIN(2)-EMAP_SOCKAPEX(2)
                  PHIVECT(3)=EMAP_SOCKORIGIN(3)-EMAP_SOCKAPEX(3)
C***              Normalise vectors
                  SUM1=0.0d0
                  SUM2=0.0d0
                  SUM1=THETAVECT(1)*THETAVECT(2)+THETAVECT(2)*
     '              THETAVECT(2)+THETAVECT(3)*THETAVECT(3)
                  SUM2=PHIVECT(1)*PHIVECT(1)+PHIVECT(2)*PHIVECT(2)+
     '              PHIVECT(3)*PHIVECT(3)
                  SUM1=DSQRT(SUM1)
                  SUM2=DSQRT(SUM2)
                  IF(DABS(SUM1).GT.1.0d-6) THEN
                    THETAVECT(1)=THETAVECT(1)/SUM1
                    THETAVECT(2)=THETAVECT(2)/SUM1
                    THETAVECT(3)=THETAVECT(3)/SUM1
                  ENDIF
                  IF(DABS(SUM2).GT.1.0d-6) THEN
                    PHIVECT(1)=PHIVECT(1)/SUM2
                    PHIVECT(2)=PHIVECT(2)/SUM2
                    PHIVECT(3)=PHIVECT(3)/SUM2
                  ENDIF
C***              Make sure THETAVECT is orthogonal to PHIVECT
                  SUM1=THETAVECT(1)*PHIVECT(1)+THETAVECT(2)*PHIVECT(2)
                  THETAVECT(3)=-SUM1/PHIVECT(3)
C***              Renormalise THETAVECT
                  SUM1=THETAVECT(1)*THETAVECT(1)+THETAVECT(2)*
     '              THETAVECT(2)+THETAVECT(3)*THETAVECT(3)
                  SUM1=DSQRT(SUM1)
                  IF(DABS(SUM1).GT.1.0d-6) THEN
                    THETAVECT(1)=THETAVECT(1)/SUM1
                    THETAVECT(2)=THETAVECT(2)/SUM1
                    THETAVECT(3)=THETAVECT(3)/SUM1
                  ENDIF
C***              Find the rotation angles
                  XVECT(1)=THETAVECT(1)
                  XVECT(2)=THETAVECT(2)
                  XVECT(3)=THETAVECT(3)
                  ZVECT(1)=PHIVECT(1)
                  ZVECT(2)=PHIVECT(2)
                  ZVECT(3)=PHIVECT(3)
                  YVECT(1)=XVECT(2)*ZVECT(3)-XVECT(3)*ZVECT(2)
                  YVECT(2)=XVECT(3)*ZVECT(1)-XVECT(1)*ZVECT(3)
                  YVECT(3)=XVECT(1)*ZVECT(2)-XVECT(2)*ZVECT(1)
                  THETA=DATAN_MOD(THETAVECT(1),THETAVECT(2))
                  PHI=DATAN_MOD(PHIVECT(3),PHIVECT(1))
                  CALL SETTRANSMAT(THETA,PHI,EMAP_SOCKORIGIN,TRANSMAT,
     '              ERROR,*9999)
                  TRANSMAT(1,1)=XVECT(1)
                  TRANSMAT(2,1)=XVECT(2)
                  TRANSMAT(3,1)=XVECT(3)
                  TRANSMAT(1,2)=YVECT(1)
                  TRANSMAT(2,2)=YVECT(2)
                  TRANSMAT(3,2)=YVECT(3)
                  TRANSMAT(1,3)=ZVECT(1)
                  TRANSMAT(2,3)=ZVECT(2)
                  TRANSMAT(3,3)=ZVECT(3)
C***              Transform the geometric positions of the data points
                  nodata=0
                  DO device=1,EMAP_NUMDEVICES(nr_signal)
                    IF(EMAP_DEVICETYPE(device,nr_signal).EQ.
     '                EMAP_ELECTRODE) THEN
                      nodata=nodata+1
                      nd=NDDATA(nodata,nr)
                      CALL TRANSVECT(TRANSMAT,ZD(1,nd),ERROR,*9999)
                    ENDIF
                  ENDDO !device
                ENDIF
              ENDIF
            ENDIF
          ENDDO !nr_signal

C***      Write the CMISS electrode data

          CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'WRITE','BINARY',OUTSIGFNAME,
     '      'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C***      Write the signal data

          IF(.NOT.ACTIVATION) THEN
CC JMB 15-AUG-2000 Write out only the signal data

C***        Rewind and skip to the beginning of the signal data
            CALL BINSETFILE(IOFILE1,0,ERROR,*9999)
            CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES,ERROR,*9999)

            DO sample=1,EMAP_NUMSAMPLES
C***          Read the signal buffer
              IF(INPUTTYPE.EQ.EMAP_FLOAT_VALUE) THEN
                CALL BINREADFILE(IOFILE1,SPTYPE,EMAP_NUMSIGNALS,INTDATA,
     '            EMAP_FSIGBUFF,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)
                DO nr_signal=1,EMAP_NUMREGIONS
C cpb 10/02/00 Adding read of specific regions
                  IF(INLIST(nr_signal,NRLIST2(1),NRLIST2(0),POSITION))
     '              THEN
                    nr=NRLIST(POSITION)
                    nodata=0
                    DO device=1,EMAP_NUMDEVICES(nr_signal)
                      IF(EMAP_DEVICETYPE(device,nr_signal).EQ.
     '                  EMAP_ELECTRODE) THEN
                        sigdevice=EMAP_SIGNALDEVICE(device,nr_signal)+1
                        nodata=nodata+1
                        nd=NDDATA(nodata,nr)
C cpb 28/10/99 Unemap now uses offsets and gains for floats as well.
C                    ZD(NJT+1,nd)=DBLE(EMAP_FSIGBUFF(sigdevice))
                        ZD(NJT+1,nd)=(DBLE(EMAP_FSIGBUFF(sigdevice))-
     '                    DBLE(EMAP_CHANNELOFFSET(device,nr_signal)))*
     '                    DBLE(EMAP_CHANNELGAIN(device,nr_signal))
                      ENDIF
                    ENDDO !device
                  ENDIF
                ENDDO !nr_signal
              ELSE
                CALL BINREADFILE(IOFILE1,SINTTYPE,EMAP_NUMSIGNALS,
     '            INTDATA,REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,
     '            EMAP_SISIGBUFF,ERROR,*9999)
                DO nr_signal=1,EMAP_NUMREGIONS
C cpb 10/02/00 Adding read of specific regions
                  IF(INLIST(nr_signal,NRLIST2(1),NRLIST2(0),POSITION))
     '              THEN
                    nr=NRLIST(POSITION)
                    nodata=0
                    DO device=1,EMAP_NUMDEVICES(nr_signal)
                      IF(EMAP_DEVICETYPE(device,nr_signal).EQ.
     '                  EMAP_ELECTRODE) THEN
                        sigdevice=EMAP_SIGNALDEVICE(device,nr_signal)+1
                        nodata=nodata+1
                        nd=NDDATA(nodata,nr)
                        ZD(NJT+1,nd)=(DBLE(EMAP_SISIGBUFF(sigdevice))-
     '                    DBLE(EMAP_CHANNELOFFSET(device,nr_signal)))*
     '                    DBLE(EMAP_CHANNELGAIN(device,nr_signal))
                      ENDIF
                    ENDDO !device
                  ENDIF
                ENDDO !nr_signal
              ENDIF

C***          Write out the signal to the CMISS signal file
              TIME=DBLE(sample-1)/EMAP_FREQUENCY
              CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '          SIGNALMIN,TIME,WD,XID,ZD,'WRITE','BINARY',OUTSIGFNAME,
     '          'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

            ENDDO !sample


            WRITE(OP_STRING(1),'(''   Num samples      '',I11)')
     '        EMAP_NUMSAMPLES
            WRITE(OP_STRING(2),'(''   Final time       '',F11.3)')
     '        TIME
            CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
          ELSE
CC JMB 15-AUG-2000 Write out the activation times stored in ZD at one
CC timestep. The epicardial potential distribution is not required
CC within CMISS once the activation times have been obtained.
            nr_signal=1
            nr=1
            nodata=0
            DO device=1,EMAP_NUMDEVICES(nr_signal)
              IF(EMAP_DEVICETYPE(device,nr_signal).EQ.EMAP_ELECTRODE)
     '          THEN
                nodata=nodata+1
                nd=NDDATA(nodata,nr)
                ZD(NJT+1,nd)=DBLE(EMAP_FSIGBUFF(device))
              ENDIF
            ENDDO ! device

            TIME=1.0d0/EMAP_FREQUENCY
            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'WRITE','BINARY',OUTSIGFNAME,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
          ENDIF

C***      Close emap signal file
          CALL BINARYCLOSEFILE(IOFILE1,ERR,CERROR)
          IF(ERR.NE.0) THEN
            CALL CSTRINGLEN(CERRLEN,CERROR)
            CALL C2FSTRING(CERROR,CERRLEN,ERROR)
            GOTO 9999
          ENDIF
        ENDIF
C***    Close CMISS signal file
        CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE','BINARY',OUTSIGFNAME,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)

      ENDIF

      CALL EXITS('IMSIGN')
      RETURN
 9999 CALL ERRORS('IMSIGN',ERROR)
      CALL EXITS('IMSIGN')
      IF(ISBINFILEOPEN(IOFILE1))
     '  CALL BINARYCLOSEFILE(IOFILE1,ERR,CERROR)
      IF(ISBINFILEOPEN(IOFILE2))
     '  CALL BINARYCLOSEFILE(IOFILE2,ERR,CERROR)
      RETURN 1
      END


