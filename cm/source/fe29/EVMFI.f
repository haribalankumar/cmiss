      SUBROUTINE EVMFI(ISIZE_MFI,LD,NBJ,NDDATA,MFI,
     '  WD,XID,ZD,STRING,ERROR,*)

C#### Subroutine: EVMFI
C###  Description:
C###    <HTML>
C###    EVMFI evaluates the MFI array from a SQUID signal file.
C###    The signal file is assumed to either:
C###    <OL>
       
C###    <LI> have at least 37 channels - 8 noise channels, 5 x component,
C###       5 y component and 19 in the z component.
C###      <OL>
C###      <LI> X -- 10,13,20,21,34
C###      <LI> Y -- 4,11,24,27,37
C###      <LI> Z -- 3,5,7,9,12,14,15,16,17,18,19,22,23,25,26,30,33,35,36
C###      <LI> Noise -- 1,2,6,8,28,29,31,32
C###      </OL>
C###    <LI> have 19 channels all representing z channel recordings
C###    </OL>
C###    Note: use with caution - will need to be updated for the new
C###    2011 SQUID.
C###    </HTML>  

C*** Created by Leo Cheng 10 February 2003
      
      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'sign00.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'trsf00.cmn'

!     Parameter List
      INTEGER ISIZE_MFI(3,NSSM),LD(NDM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM)
      REAL*8 MFI(NDM,NTSM,3,NSSM),WD(NJM,NDM),XID(NIM,NDM),ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,N3CO,nd,nj,notime,nss,
     '  NUMTIMEDATA
      REAL*8 SIGNALMAX(9),SIGNALMIN(9),TIME
      CHARACTER SIGNALFNAME*(MXCH),FILEFORMAT*6
      LOGICAL CBBREV,COMBINE,ENDFILE,NOISE,ZONLY

      INTEGER CHANNELTYPES,SQUIDCHANNELS,CMISSCHANNELS
      PARAMETER(CHANNELTYPES=4,CMISSCHANNELS=20,SQUIDCHANNELS=37)
      INTEGER SQUIDCOMB(0:SQUIDCHANNELS),
     '  SQUIDELEC(CHANNELTYPES,0:19)
      REAL*8 SQUIDCOMPONENTS(3,SQUIDCHANNELS)

!     Functions
      INTEGER IFROMC

      DATA SQUIDELEC /
     '  5,  5,19, 8,
     '  10, 4, 3, 1,
     '  13,11, 5, 2,
     '  20,24, 7, 6,
     '  21,27, 9, 8,
     '  34,37,12,28,
     '   0, 0,14,29,
     '   0, 0,15,31,
     '   0, 0,16,32,
     '   0, 0,17, 0,
     '   0, 0,18, 0,
     '   0, 0,19, 0,
     '   0, 0,22, 0,
     '   0, 0,23, 0,
     '   0, 0,25, 0,
     '   0, 0,26, 0,
     '   0, 0,30, 0,
     '   0, 0,33, 0,
     '   0, 0,35, 0,
     '   0, 0,36, 0
     '  /

C Global components of each squid channel from gut2posor.m from
C     Vanderbilt
C NOTE: these need to be rotated as the gut2posor file is WRONG

      DATA SQUIDCOMPONENTS /
     '  -.5,      .5, -.7071, !sensor 1
     '  0.,       0.,    -1.,
     '  0.,       0.,    -1.,
     '  -.7071,.7071,     0.,
     '  0.,       0.,    -1., !sensor 5
     '  -.7071,.7071,     0.,
     '  0.,       0.,    -1.,
     '  0.,       1.,     0.,
     '  0.,       0.,    -1.,
     '  .7071, .7071,     0., !sensor 10
     '  -.7071,.7071,     0.,
     '  0,        0.,    -1.,
     '  .7071, .7071,     0., 
     '  0.,       0.,    -1.,
     '  0.,       0.,    -1., !sensor 15
     '  0.,       0.,    -1.,
     '  0.,       0.,    -1.,
     '  0.,       0.,    -1.,
     '  0.,       0.,    -1.,
     '  .7071, .7071,     0., !sensor 20
     '  .7071, .7071,     0.,
     '  0.,       0.,    -1.,
     '  0.,       0.,    -1.,
     '  -.7071,.7071,     0.,
     '  0.,       0.,    -1., !sensor 25
     '  0.,       0.,    -1.,
     '  -.7071,.7071,     0.,
     '  -.7071,.7071,     0.,
     '  0.,       0.,    -1.,
     '  0.,       0.,    -1., !sensor 30
     '  .7071, .7071,     0.,
     '  .5,       .5,  .7071,
     '  0.,       0.,    -1.,
     '  .7071, .7071,     0.,
     '  0.,       0.,    -1., !sensor 35
     '  0.,       0.,    -1.,
     '  -.7071,.7071,     0.  !sensor 37
     '  /
      
C Mapping for combining the 37 squid channels locations
C     into the 20 cmiss channel locations      
      DATA SQUIDCOMB / 20,
     '  20,20, 2,11,11,
     '  20,12,20, 4,11,
     '  14,15,14,16, 5,
     '   8,17,18, 7, 8,
     '  17,18, 7, 1, 6,
     '   8, 8,20,20, 9,
     '  20,20, 1, 1, 3,
     '  10,17
     '  /
      
      CALL ENTERS('EVMFI',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM evaluate MFI<;FILENAME>
C###  Parameter:        <tstart (#/beginning)[beginning]>
C###    Specify start time
C###  Parameter:        <tend (#/end)[end]>
C###    Specify end time
C###  Parameter:        <(ascii/binary)[ascii]>
C###    Specify whether the file is stored as a binary or
C###    ascii file.
C###  Parameter:        <nss #[1]>
C###    Specify the data set number for the MFI array.
C###  Parameter:        <region #[1]>
C###    Specify the region where the sensor positions are stored.
C###  Parameter:        <combine>
C###    Combined the 37 SQUID Channels into the 20 CMISS Channels
C###  Parameter:        <noise>
C###    Reject the 8 noise channels
C###  Parameter:        <zonly>
C###    Manipulate the 19 z compent SQUID Channels
C###  Description:
C###    Evaluates the MFI matrix from a SQUID signal file.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<signal FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(3)=BLANK(1:15)//'<tstart (#/beginning)[beginning]>'
        OP_STRING(4)=BLANK(1:15)//'<tend (#/end)[end]>'
        OP_STRING(5)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(6)=BLANK(1:15)//'<nss [1]'
        OP_STRING(7)=BLANK(1:15)//'<region #[1]'
        OP_STRING(8)=BLANK(1:15)//'<combine>'
        OP_STRING(9)=BLANK(1:15)//'<noise>'
        OP_STRING(10)=BLANK(1:15)//'<zonly>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
 
C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVPHI',ERROR,*9999)
      ELSE

        CALL ASSERT(USE_MAGNETIC.EQ.1,
     '    '>>Set USE_MAGNETIC to 1 in parameter file',ERROR,*9999)

        IF(CBBREV(CO,'SIGNAL',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          SIGNALFNAME=CO(N3CO+1)(IBEG:IEND)
        ELSE
          CALL STRING_TRIM(FILE00,IBEG,IEND)
          SIGNALFNAME=FILE00(IBEG:IEND)
        ENDIF

        IF(CBBREV(CO,'BINARY',3,noco+1,NTCO,N3CO)) THEN
          FILEFORMAT='BINARY'
        ELSE
          FILEFORMAT='ASCII'
          ERROR='Only implemented for Binary files'
          GOTO 9999
        ENDIF

        IF(CBBREV(CO,'NSS',3,noco+1,NTCO,N3CO)) THEN
          nss=IFROMC(CO(N3CO+1))
          IF(nss.LT.NSSM) THEN
            IEND=0
            CALL APPENDC(IEND,'>> Increase NSSM to at least ',ERROR)
            CALL APPENDI(IEND,nss,ERROR)
            GOTO 9999
          ENDIF
        ELSE
          nss=1
        ENDIF
        
        IF(CBBREV(CO,'NOISE',3,noco+1,NTCO,N3CO)) THEN
          NOISE=.TRUE.
        ELSE
          NOISE=.FALSE.
        ENDIF

        IF(CBBREV(CO,'COMBINE',3,noco+1,NTCO,N3CO)) THEN
          COMBINE=.TRUE.
        ELSE
          COMBINE=.FALSE.
        ENDIF

        IF(CBBREV(CO,'ZONLY',3,noco+1,NTCO,N3CO)) THEN
          ZONLY=.TRUE.
        ELSE
          ZONLY=.FALSE.
        ENDIF

        IF(.NOT.ZONLY) THEN
          OP_STRING(1)=' '
          OP_STRING(2)='Evaluating SQUID Channels'
          OP_STRING(3)='-------------------------'
          CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
          
          CALL WRITE_LONG(INTTYPE,1,CHANNELTYPES,IOOP,
     '      SQUIDELEC(3,0)*CHANNELTYPES,5,5,SQUIDELEC(3,1),%VAL(0),
     '      '('' z-component Channels:'',5I5)','(22X,5I5)',ERROR,*9999)
          CALL WRITE_LONG(INTTYPE,1,CHANNELTYPES,IOOP,
     '      SQUIDELEC(1,0)*CHANNELTYPES,5,5,SQUIDELEC(1,1),%VAL(0),
     '      '('' x-component Channels:'',5I5)','(22X,5I5)',ERROR,*9999)

          CALL WRITE_LONG(INTTYPE,1,CHANNELTYPES,IOOP,
     '      SQUIDELEC(2,0)*CHANNELTYPES,5,5,SQUIDELEC(2,1),%VAL(0),
     '      '('' y-component Channels:'',5I5)','(22X,5I5)',ERROR,*9999)
        
          CALL WRITE_LONG(INTTYPE,1,CHANNELTYPES,IOOP,
     '      SQUIDELEC(3,0)*CHANNELTYPES,5,5,SQUIDELEC(3,1),%VAL(0),
     '      '('' z-component Channels:'',5I5)','(22X,5I5)',ERROR,*9999)
          
          CALL WRITE_LONG(INTTYPE,1,CHANNELTYPES,IOOP,
     '      SQUIDELEC(4,0)*CHANNELTYPES,5,5,SQUIDELEC(4,1),%VAL(0),
     '      '('' Noise Channels:      '',5I5)','(22X,5I5)',ERROR,*9999)
        ENDIF
        
        TIME=0.d0 !initialise TIME.
        
C***    Open up the signal file
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,SIGNALFNAME,
     '    'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

C***    Read electrode geometry etc.
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,SIGNALFNAME,
     '    'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)
        
        CALL ASSERT(SIGNAL_NUMREGIONS(IOFILE1).EQ.1,
     '    '>> Not implemented for multiple regions',ERROR,*9999)

C LKC 18-MAR-2003  Can't use NJT because the signals are stored in
C     a 2D coordinate system (patch rig)
C       DO nj=1,NJT
        
C*** Check the MFI array is large enough        
        IF(NUMTIMEDATA.GT.NTSM) THEN
          IEND=0
          CALL APPENDC(IEND,'>> Increase NTSM to at least ',ERROR)
          CALL APPENDI(IEND,NUMTIMEDATA,ERROR)
          GOTO 9999
        ENDIF

        IF(SQUIDCHANNELS.GT.NDM) THEN
          IEND=0
          CALL APPENDC(IEND,'>> Increase NDM to at least ',ERROR)
          CALL APPENDI(IEND,SQUIDCHANNELS,ERROR)
          GOTO 9999
        ENDIF

C***    Initialise the appropriate parts of MFI
        DO nj=1,3
          DO notime=1,NUMTIMEDATA !Loop over the times
            DO nd=1,SQUIDCHANNELS
              MFI(nd,notime,nj,nss)=0.d0
            ENDDO
          ENDDO
        ENDDO


        IF(COMBINE) THEN
          IF(NOISE) THEN
            ISIZE_MFI(1,nss)=CMISSCHANNELS
          ELSE
            ISIZE_MFI(1,nss)=CMISSCHANNELS-1
          ENDIF
        ELSE
          IF(ZONLY) THEN
            ISIZE_MFI(1,nss)=19            
          ELSE
            ISIZE_MFI(1,nss)=SQUIDCHANNELS
          ENDIF
        ENDIF
        ISIZE_MFI(2,nss)=NUMTIMEDATA
        ISIZE_MFI(3,nss)=3

        CALL ASSERT(NUMTIMEDATA.LE.NTSM,'>> Increase NTSM',ERROR,*9999)
        IF(NOISE) THEN
          CALL ASSERT(ISIZE_MFI(1,nss).LE.NDM,
     '      '>> Increase NDM',ERROR,*9999)
        ELSE
C +1 incase there is a repiration channel - which we currently ignore.
          CALL ASSERT(ISIZE_MFI(1,nss)+1.LE.NDM,
     '      '>> Increase NDM',ERROR,*9999)
        ENDIF
        
        DO notime=1,NUMTIMEDATA !Loop over the times

C***      Read the signal data in
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,SIGNALFNAME,
     '      'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C*** Note Currently implemented in Vanderbilt coordinates
C***          (y +ve == patients head, z +ve to ceiling
C*** x channels are 2 oclock, y channels are 10 oclock

          IF(COMBINE) THEN !output combined vectors

            DO nd=1,SQUIDCHANNELS !setup the combined channels
C LKC 18-MAR-2003  Can't use NJT because the signals are stored in
C     a 2D coordinate system (patch rig)
C              DO nj=1,NJT
              DO nj=1,3
                MFI(SQUIDCOMB(nd),notime,nj,nss)=
     '            MFI(SQUIDCOMB(nd),notime,nj,nss)+
     '            SQUIDCOMPONENTS(nj,nd)*ZD(NJT+1,nd)
              ENDDO
            ENDDO !nd
            
          ELSEIF(ZONLY) THEN !only maniputate the 19 z channels

C Reorder from VU to CMISS numbering            
C            nj=2 ! only the "z" compnents            
C            MFI(2,notime,nj,nss)=ZD(NJT+1,1)
C            MFI(11,notime,nj,nss)=ZD(NJT+1,2)
C            MFI(12,notime,nj,nss)=ZD(NJT+1,3)
C            MFI(4,notime,nj,nss)=ZD(NJT+1,4)
C            MFI(15,notime,nj,nss)=ZD(NJT+1,5)
C            MFI(16,notime,nj,nss)=ZD(NJT+1,6)
C            MFI(5,notime,nj,nss)=ZD(NJT+1,7)
C            MFI(13,notime,nj,nss)=ZD(NJT+1,8)
C            MFI(14,notime,nj,nss)=ZD(NJT+1,9)
C            MFI(17,notime,nj,nss)=ZD(NJT+1,10)
C            MFI(19,notime,nj,nss)=ZD(NJT+1,11)
C            MFI(18,notime,nj,nss)=ZD(NJT+1,12)
C            MFI(7,notime,nj,nss)=ZD(NJT+1,13)
C            MFI(6,notime,nj,nss)=ZD(NJT+1,14)
C            MFI(8,notime,nj,nss)=ZD(NJT+1,15)
C            MFI(9,notime,nj,nss)=ZD(NJT+1,16)
C            MFI(1,notime,nj,nss)=ZD(NJT+1,17)
C            MFI(3,notime,nj,nss)=ZD(NJT+1,18)
C            MFI(10,notime,nj,nss)=ZD(NJT+1,19)
C
C!!! 7-FEB-2007, unsure about this stuff. Use with caution!
C Now assume the signal numbers have already been reordered
C
            nj=2
            DO nd=1,ISIZE_MFI(1,nss) !19
              MFI(nd,notime,nj,nss)=ZD(NJT+1,nd)
            ENDDO
            
          ELSE !all individual components -- should be 37 channels
            DO nd=1,SQUIDCHANNELS
              DO nj=1,3
                MFI(nd,notime,nj,nss)=SQUIDCOMPONENTS(nj,nd)*
     '            ZD(NJT+1,nd)
              ENDDO
            ENDDO !nd
          ENDIF
        ENDDO !notime

C*** Close signal file
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,SIGNALFNAME,
     '    ' ',ENDFILE,.TRUE.,ERROR,*9999)
        
      ENDIF

      CALL EXITS('EVMFI')
      RETURN      
 9999 CALL ERRORS('EVMFI',ERROR)
      CALL EXITS('EVMFI')
      RETURN
      END

      
