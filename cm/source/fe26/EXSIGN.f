      SUBROUTINE EXSIGN(ISIZE_MFI,ISIZE_PHI,ISIZE_PHIH,ISIZE_TBH,
     '  LD,MFI,NBJ,NDDATA,NDP,NPLIST3,PHI,PHI_H,WD,XID,
     '  ZCROSSING,ZD,CELL_YQS_NAMES,STRING,ERROR,*)

C#### Subroutine: EXSIGN
C###  Description:
C###    <HTML>
C###       <p>EXSIGN exports signals from the FE data base. Signals
C###       can either be exported to UNEMAP via a .signal/.sig file
C###       or individual leads exported to CMGUI via a graphics
C###       object.
C###       </P><p>
C###       When exporting to CMGUI the lead is denoted by
C###       EMAP_STARTELEC(1)
C###       <OL>
C###        <LI> IOFILE1 - input signal file
C###        <LI> IOFILE2 - output for the signal trace
C###        <LI> IOFILE3 - output for the pointer
C###       </OL></p>
C###       <p>
C###       Note that when exporting magnetic fields,
C###       by default CMISS calculates H [uA/mm]. To export as
C###       Tesla (using B=uH) requires a conversion factor of   
C###       4.PI.1E-10 as the SI units for H is [A/m].
C###       </p>      
C###    </HTML>


      
      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'binf00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'cmgui01.cmn'
      INCLUDE 'emap00.cmn'
      INCLUDE 'expo00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'lead00.cmn'
      INCLUDE 'mach00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'sign00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER ISIZE_MFI(3,NSSM),ISIZE_PHI(2),ISIZE_PHIH(2),
     &  ISIZE_TBH(2),LD(NDM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),
     &  NDP(NDM),NPLIST3(0:NPM)
      REAL*8 MFI(NDM,NTSM,3,NSSM),PHI(NY_TRANSFER_M,NTSM),
     &  PHI_H(NY_TRANSFER_M,NTSM),WD(NJM,NDM),XID(NIM,NDM),
     '  ZCROSSING(NY_TRANSFER_M,NTSM),
     '  ZD(NJM,NDM)
      CHARACTER ERROR*(*),CELL_YQS_NAMES(NIQSM,NQVM)*(*),STRING*(MXCH)

!     Local Variables
      INTEGER CERROR(50),COMPONENT,DEVICE,ERR,EMAP_NUMREGIONS_TOEXPORT,
     '  IBEG,IBEG1,IDUMMY(1),
     '  IEND,IEND1,N3CO,ncol,nd,ndevice,nelec,noelec,nonr,
     '  nj,nlead,nleadelec,nr,nrdata,
     '  nsample,nss,NUM_SAMPLES,NUM_SIGNALS,NUMTIMEDATA,NUMTIMEDATA1,
     '  UNEMAP_TSTART,NUMRIGNAMEBYTES(2),NUMREGIONS(2),sample,
     '  ZCROSSING_INDEX(NY_TRANSFER_M)
      REAL*8 DELTATIME,RDUMMY(1),FREQUENCY,POTENTIAL,PREVTIME,
     '  SIGNALMAX(9),SIGNALMIN(9),RMSSUM,
     '  ZCROSSING_RANGE(NY_TRANSFER_M,3),TIME,TSTART
      CHARACTER CHAR1*4,COMPFNAME*100,
     '  EMAPFNAME*100,FILE*100,SIGNALFNAME*100,
     '  FILENAME*100,DEVICE_NAMES*10
      LOGICAL ABBREV,CBBREV,COMPSIGN,ELECTRODES,EMAP,EXPORT_MFI,CMGUI,
     &  DATAFILE,ENDFILE,ISBINFILEOPEN,LEADS,
     '  PHI_OR_PHIH,PHIH,TESLA,RMSLEAD,WARNING,ZERO_CROSSING

!     Functions
      INTEGER IFROMC
      
      CALL ENTERS('EXSIGN',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM export signal<;FILENAME[default]>
C###  Parameter:      <electrodes/leads[both]>
C###    Specify whether to export electrode or lead signals or both.
C###    Electrodes are point values and leads are combinations of
C###    electrodes.
C###  Parameter:      <signal FILENAME[default]>
C###    Specify the binary signal (binsig)
C###    file to export from. Not implemented for ascii signals.
C###  Parameter:      <compare_signal FILENAME[default]>
C###    An addition signal comparison file if exporting
C###    to a DATA file
C###  Parameter:      <rmslead>
C###    Export an additional lead which is calculated from the rms
C###    value of the specified electrodes at each point in time
C###  Parameter:      <names (numbers/cell/data) [numbers]>
C###    Specify whether to use sequential numbers, data point labels
C###        or cell parameter names as the device names in UnEMAP.
C###  Description:
C###    Exports CMISS signals from the FE data base.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME['//FILE00(IBEG1:IEND1)
     '    //']> <(electrodes/leads)[both]>'
        OP_STRING(2)=BLANK(1:15)//'<signal FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(3)=BLANK(1:15)//'<compare_signal FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(4)=BLANK(1:15)//'<rmslead>'
        OP_STRING(5)=BLANK(1:15)//'<names (numbers/cell/data)[numbers]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM export signal<;FILENAME[default]>
C###  Description:
C###    Exports CMISS signals from the (zero_crossing/phi/phi_h/mfi) array.
C###  Parameter:      <tesla>
C###    Convert MFI (magnetic field intensity, H) solutions to magnetic
C###    flux intensity (B with units Tesla) via the equation B=uH.
C###  Parameter:      <component #[2]>
C###    Determine which component magnetic field to export.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME['//FILE00(IBEG1:IEND1)
     '    //']>'
        OP_STRING(2)=BLANK(1:15)//'<(zero_crossing/phi/phi_h/mfi)>'
        OP_STRING(3)=BLANK(1:15)//'<tesla>'
        OP_STRING(4)=BLANK(1:15)//'<component #[2]>' 
       CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','EXSIGN',ERROR,*9999)
      ELSE


        CALL ASSERT(CALL_EXPO,'>>Define export first',ERROR,*9999)
        CALL CHECKF(1,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

C LKC 31-AUG-1999 This is a local variable to that it does not
C  get overwritten when the signal to export is read in.
C  This allows a different number of signals to be exported.
C  Needs to be generalised to allow a start and end region to export.
        EMAP_NUMREGIONS_TOEXPORT=EMAP_NUMREGIONS



        IF(CBBREV(CO,'COMPARE_SIGNAL',3,noco+1,NTCO,N3CO)) THEN
          COMPSIGN=.TRUE.
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          COMPFNAME=CO(N3CO+1)(IBEG:IEND)
        ELSE
          COMPSIGN=.FALSE.
        ENDIF

        IF(CBBREV(CO,'ZERO_CROSSING',3,noco+1,NTCO,N3CO)) THEN
          ZERO_CROSSING=.TRUE.
        ELSE
          ZERO_CROSSING=.FALSE.
        ENDIF

C GBS 8-Aug-2000
C PHI_OR_PHIH is true for both "ex sig phi" and "ex sig phi_h"
C PHIH is true only for the secong of those
        PHI_OR_PHIH=.FALSE.

        IF(CBBREV(CO,'PHI',3,noco+1,NTCO,N3CO)) THEN
          PHI_OR_PHIH=.TRUE.
          PHIH=.FALSE.
        ENDIF
        IF(CBBREV(CO,'PHI_H',5,noco+1,NTCO,N3CO)) THEN
          PHI_OR_PHIH=.TRUE.
          PHIH=.TRUE.
        ENDIF

        TESLA=.FALSE.
        IF(CBBREV(CO,'MFI',3,noco+1,NTCO,N3CO)) THEN
          EXPORT_MFI=.TRUE.
          LEADS=.FALSE.

          IF(CBBREV(CO,'TESLA',3,noco+1,NTCO,N3CO)) THEN
            TESLA=.TRUE.
          ELSE
            TESLA=.FALSE.
          ENDIF

          IF(CBBREV(CO,'COMPONENT',4,noco+1,NTCO,N3CO)) THEN
            COMPONENT=IFROMC(CO(N3CO+1))
          ELSE
            COMPONENT=2
          ENDIF
        ELSE
          EXPORT_MFI=.FALSE.
        ENDIF
        nss=1 !tmp
        WRITE(OP_STRING,
     &    '(''>>WARNING: Only implimented to export nss 1 ''
     &    )')

C *** DPN 19 February 2000 - Adding ability to use names rather than
C *** data point numbers for the device names in UnEMAP signal files.
C *** Start with the cellular name arrays, but should be able to extend
C *** for any other set of strings.
        IF(CBBREV(CO,'NAMES',4,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'NUMBERS',4)) THEN
            DEVICE_NAMES='NUMBERS'
          ELSEIF(ABBREV(CO(N3CO+1),'CELL',4)) THEN
            DEVICE_NAMES='CELL'
          ELSEIF(ABBREV(CO(N3CO+1),'DATA',4)) THEN
            DEVICE_NAMES='DATA'
          ELSE
            ERROR='Invalid name type'
            GOTO 9999
          ENDIF
        ELSE
          DEVICE_NAMES='NUMBERS'
        ENDIF

        IF(.NOT.ZERO_CROSSING.AND..NOT.PHI_OR_PHIH.AND..NOT.EXPORT_MFI)
     &    THEN
          IF(CBBREV(CO,'ELECTRODES',2,noco+1,NTCO,N3CO)) THEN
            ELECTRODES=.TRUE.
            LEADS=.FALSE.
          ELSE IF(CBBREV(CO,'LEADS',2,noco+1,NTCO,N3CO)) THEN
            ELECTRODES=.FALSE.
            LEADS=.TRUE.
          ELSE
            ELECTRODES=.TRUE.
            LEADS=.TRUE.
          ENDIF

          IF(LEADS) THEN
            CALL ASSERT(NUMLEADS.GT.0,'>>No leads defined',ERROR,*9999)
          ENDIF

          IF(CBBREV(CO,'SIGNAL',2,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
            SIGNALFNAME=CO(N3CO+1)(IBEG:IEND)
          ELSE
            CALL STRING_TRIM(FILE00,IBEG,IEND)
            SIGNALFNAME=FILE00(IBEG:IEND)
          ENDIF

c cpb 5/9/00 Adding rmslead
          IF(CBBREV(CO,'RMSLEAD',1,noco+1,NTCO,N3CO)) THEN
            RMSLEAD=.TRUE.
          ELSE
            RMSLEAD=.FALSE.
          ENDIF

          IF(RMSLEAD.AND.LEADS) THEN
            CALL ASSERT(NUMLEADS+1.LT.NUMLEADMX,
     '        '>>Increase NUMLEADMX in lead00.cmn',ERROR,*9999)
          ENDIF

        ELSE
          ELECTRODES=.TRUE.
          LEADS=.FALSE.
          RMSLEAD=.FALSE.
        ENDIF ! exporting normal signals


        EMAP=.FALSE.
        CMGUI=.FALSE.
        DATAFILE=.FALSE.
        IF(.NOT.ZERO_CROSSING.AND..NOT.PHI_OR_PHIH) THEN
          IF(SIGEXPO_TYPE.EQ.1) THEN
            EMAP=.TRUE.
          ELSEIF(SIGEXPO_TYPE.EQ.2) THEN
            CMGUI=.TRUE.
          ELSEIF(SIGEXPO_TYPE.EQ.3) THEN
            DATAFILE=.TRUE.
          ELSE
            ERROR='>>Unknown signal export type'
            GOTO 9999
          ENDIF
C GBS  9-Aug-2000 Only Unemap output supported
        ELSEIF(PHI_OR_PHIH) THEN
          IF(SIGEXPO_TYPE.EQ.1) THEN
            EMAP=.TRUE.
          ELSEIF(SIGEXPO_TYPE.EQ.2) THEN
          ELSEIF(SIGEXPO_TYPE.EQ.3) THEN
          ELSE
            ERROR='>>Unknown signal export type'
            GOTO 9999
          ENDIF
        ELSEIF(ZERO_CROSSING) THEN
          IF(SIGEXPO_TYPE.EQ.1) THEN
            EMAP=.TRUE.
          ELSEIF(SIGEXPO_TYPE.EQ.2) THEN
            DATAFILE=.TRUE.
          ELSEIF(SIGEXPO_TYPE.EQ.3) THEN
            DATAFILE=.TRUE.
          ELSE
            ERROR='>>Unknown signal export type'
            GOTO 9999
          ENDIF
        ELSE
          ERROR='Unknown export type'
          GOTO 9999
        ENDIF


        
C*** Determining the number of samples

        TIME=0.0d0
        IDUMMY(1)=0
        RDUMMY(1)=0.0d0        
        IF(.NOT.ZERO_CROSSING.AND..NOT.PHI_OR_PHIH.AND..NOT.EXPORT_MFI)
     &    THEN
C***      Open signal file.
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ','BINARY',SIGNALFNAME,
     '      'OPEN',ENDFILE,.TRUE.,ERROR,*9999)
          NUM_SAMPLES=NUMTIMEDATA

C LKC new warning if no data in signal file
          IF(NUM_SAMPLES.LE.0) THEN
            WRITE(OP_STRING,'(''>>WARNING: Num. Time Samples <= 0 '')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

C***      Read electrode geometry etc.
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ','BINARY',SIGNALFNAME,
     '      'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

          CALL ASSERT(EMAP_NUMREGIONS_TOEXPORT.LE.EMAP_NUMREGIONS,
     '      '>> Trying to export too many regions',ERROR,*9999)

          IF(COMPSIGN) THEN
C***        Open signal file.
            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ','BINARY',COMPFNAME,
     '        'OPEN',ENDFILE,.TRUE.,ERROR,*9999)
            NUM_SAMPLES=NUMTIMEDATA

C***        Read electrode geometry etc.
            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA1,1,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ','BINARY',COMPFNAME,
     '        'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

            CALL ASSERT(NUMTIMEDATA.EQ.NUMTIMEDATA1,
     '        '>>Number of time samples must be the same',ERROR,*9999)
          ENDIF

C LKC 7-SEP-1999 adding a simple check for upper bound
C check that not too many electrodes are being exported
C Assumes electrodes exported in order.

          DO nr=1,EMAP_NUMREGIONS_TOEXPORT
            nonr=EXPORT_NRLIST(nr)

C!!! LKC 9-SEP-1999 Note this assumes the data is in sequential order
C            IF(EMAP_STOPELEC(nr).GT.NDDATA(0,nr)) THEN
            IF (EMAP_STOPELEC(nr)-EMAP_STARTELEC(nr)+1
     '        .GT.        ! Last - First
     '        NDDATA(NDDATA(0,nonr),nonr)-NDDATA(1,nonr)+1) THEN

              IEND=0
              CALL APPENDC(IEND,'>>   Exporting ',OP_STRING(1))
              CALL APPENDI(IEND,
     '          EMAP_STOPELEC(nr)-EMAP_STARTELEC(nr)+1,
     '          OP_STRING(1))
              CALL APPENDC(IEND,' electrodes in region ',OP_STRING(1))
              CALL APPENDI(IEND,nr,OP_STRING(1))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

              IEND=0
              CALL APPENDC(IEND,'>>   Signal file contains ',
     '          OP_STRING(1))
              CALL APPENDI(IEND,
     '          NDDATA(NDDATA(0,nonr),nonr)-NDDATA(1,nonr)+1,
     '          OP_STRING(1))
              CALL APPENDC(IEND,' electrodes in region ',OP_STRING(1))
              CALL APPENDI(IEND,nonr,OP_STRING(1))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

              ERROR=' Too many electrodes being exported?'
              GOTO 9999
            ENDIF
          ENDDO !nr

          
C LKC 28-NOV-2002 Explicitly set for PHI and PHI_H
C        ELSE
C         NUM_SAMPLES=NTST
        ELSEIF(PHI_OR_PHIH) THEN
          IF(PHIH) THEN
            NUM_SAMPLES=ISIZE_PHIH(2)
          ELSE
            NUM_SAMPLES=ISIZE_PHI(2)
          ENDIF
        ELSEIF(EXPORT_MFI) THEN
          NUM_SAMPLES=ISIZE_MFI(2,nss)          
        ELSE
          NUM_SAMPLES=NTST
        ENDIF
        
        CALL ASSERT(NUM_SAMPLES.GT.0,
     '      '>> No time steps to export for PHI or PHI_H',ERROR,*9999)


        IF(EMAP) THEN

          TSTART=0.d0
          IF(SIGNAL_HOWGENERATED(IOFILE1).NE.1.AND.
     '      .NOT.ZERO_CROSSING.AND..NOT.PHI_OR_PHIH
     &      .AND..NOT.EXPORT_MFI) THEN

C LKC 21-NOV-2003 determining what the initial time step
C           in the cm signal file is and using that
C           as the start time for the unemap signal file latter on.
            
C*** Want to determine what the first times step is ... 
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,
     '        SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'READ','BINARY',
     '        SIGNALFNAME,'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
            TSTART=TIME

C***        Reset file
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ','BINARY',SIGNALFNAME,
     '        'RESET',ENDFILE,.TRUE.,ERROR,*9999)

C LKC 13-SEP-1999 Adding frequency override
            IF(.NOT.SET_FREQUENCY) THEN
C***          Determine the frequency of the signal
              PREVTIME=0.0d0
              DELTATIME=0.0d0
              WARNING=.FALSE.
              DO sample=1,NUM_SAMPLES
                CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,
     '            SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'READ','BINARY',
     '            SIGNALFNAME,'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

                IF(PREVTIME.NE.0.0d0) THEN
                  IF(DELTATIME.EQ.0.0d0) THEN
                    DELTATIME=TIME-PREVTIME
                  ELSE
                    IF(DABS((TIME-PREVTIME)-DELTATIME).GT.1.0d-6)
     '                WARNING=.TRUE.
                  ENDIF
                ENDIF
                PREVTIME=TIME
              ENDDO !tag
              IF(WARNING) THEN
                WRITE(OP_STRING,'('' >>Warning: Time steps are not '
     '            //'consistent'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
              IF(DELTATIME.NE.0.0d0) THEN
                FREQUENCY=1.0d0/DELTATIME
              ELSE
                ERROR='>>Zero delta time'
                GOTO 9999
              ENDIF
C***          Reset file
              CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '          SIGNALMIN,TIME,WD,XID,ZD,'READ','BINARY',SIGNALFNAME,
     '          'RESET',ENDFILE,.TRUE.,ERROR,*9999)
            ELSE

C*** Setting the frequency
              FREQUENCY=CM_FREQUENCY
            ENDIF ! SET_FREQUENCY


          ELSEIF(ZERO_CROSSING.OR.PHI_OR_PHIH.OR.EXPORT_MFI) THEN

C LKC 1-MAR-1999 input frequency from def export zeroxing
C            FREQUENCY=1 !?
            IF(CM_FREQUENCY.LT.DABS(ZERO_TOL)) THEN
              FREQUENCY=1
            ELSE
              FREQUENCY=CM_FREQUENCY              
            ENDIF
          ELSE
            FREQUENCY=EMAP_FREQUENCY
          ENDIF
          
          CALL ASSERT(FREQUENCY.GT.DABS(ZERO_TOL),
     '      '>> Frequency is zero',ERROR,*9999)
          CALL ASSERT(FREQUENCY.LT.10D3,
     '      '>> Frequency is GT 10,000',ERROR,*9999)


C*** Write EMAP signal file
          CALL STRING_TRIM(FILE,IBEG,IEND)
          EMAPFNAME=FILE(IBEG:IEND)//'.signal'

          CALL BINOPENFILE(IOFILE2,'WRITE',EMAPFNAME,ERROR,*9999)

C*** Write out the configurations

C*** Write out the rig type
          CALL BINWRITEFILE(IOFILE2,INTTYPE,1,EMAP_RIGTYPE(0),
     '      REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

          
C*** Write out the rig name
          CALL STRING_TRIM(EMAP_RIGNAME,IBEG,IEND)
          IF(IBEG.GT.IEND) THEN
            EMAP_NUMRIGNAMEBYTES=0
          ELSE IF(IBEG.EQ.IEND.AND.EMAP_RIGNAME(IBEG:IEND).EQ.' ') THEN
            EMAP_NUMRIGNAMEBYTES=0
          ELSE
            EMAP_NUMRIGNAMEBYTES=IEND-IBEG+1
          ENDIF

C          WRITE(*,*) 'Rig type',EMAP_RIGTYPE(0)
C          WRITE(*,*) 'Rig name',EMAP_RIGNAME(IBEG:IEND)

          
          NUMRIGNAMEBYTES(1)=EMAP_NUMRIGNAMEBYTES
          CALL BINWRITEFILE(IOFILE2,INTTYPE,1,NUMRIGNAMEBYTES,
     '      REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          CALL BINWRITEFILE(IOFILE2,CHARTYPE,EMAP_NUMRIGNAMEBYTES,
     '      INTDATA,REAL4DATA,REAL8DATA,EMAP_RIGNAME,LOGDATA,SINTDATA,
     '      ERROR,*9999)
          IF(EMAP_RIGTYPE(0).EQ.EMAP_SOCK) THEN
C*** Write out the sock focus if necessary
            REAL4DATA(1)=REAL(EMAP_SOCKFOCUS(1))
            CALL BINWRITEFILE(IOFILE2,SPTYPE,1,INTDATA,
     '        REAL4DATA,REAL8DATA,CHARDATA,
     '        LOGDATA,SINTDATA,ERROR,*9999)
C            CALL BINWRITEFILE(IOFILE2,SPTYPE,1,INTDATA,
C     '        REAL(EMAP_SOCKFOCUS(1)),REAL8DATA,CHARDATA,
C     '        LOGDATA,SINTDATA,ERROR,*9999)
          ENDIF
C*** Write out the number of regions
          NUMREGIONS(1)=EMAP_NUMREGIONS_TOEXPORT
          CALL BINWRITEFILE(IOFILE2,INTTYPE,1,NUMREGIONS,
     '      REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C cpb 3/2/97 Fixing export signal
          NUM_SIGNALS=0



C*** Write out the region information

C LKC 7-SEP-1999 Modification to export specific regions
C          DO nr=1,EMAP_NUMREGIONS_TOEXPORT

          DO nr=1,EMAP_NUMREGIONS_TOEXPORT
            nonr=EXPORT_NRLIST(nr)

            IF(EMAP_RIGTYPE(0).EQ.EMAP_MIXED) THEN
C*** Write out the rig region type
              CALL BINWRITEFILE(IOFILE2,INTTYPE,1,EMAP_RIGTYPE(nr),
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
            ENDIF
C*** Write out the region name
            CALL STRING_TRIM(EMAP_REGIONNAME(nr),IBEG,IEND)
            IF(IBEG.GT.IEND) THEN
              INTDATA(1)=0
            ELSE IF(IBEG.EQ.IEND.AND.EMAP_REGIONNAME(nr)(IBEG:IEND).EQ.
     '        ' ') THEN
              INTDATA(1)=0
            ELSE
              INTDATA(1)=IEND-IBEG+1
            ENDIF
            CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '        REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
            CALL BINWRITEFILE(IOFILE2,CHARTYPE,INTDATA(1),INTDATA,
     '        REAL4DATA,REAL8DATA,EMAP_REGIONNAME(nr),LOGDATA,SINTDATA,
     '        ERROR,*9999)
            IF(EMAP_RIGTYPE(0).EQ.EMAP_MIXED.AND.EMAP_RIGTYPE(nr).EQ.
     '        EMAP_SOCK) THEN
C***          Write out the sock focus if necessary
              REAL4DATA(1)=REAL(EMAP_SOCKFOCUS(nr))
              CALL BINWRITEFILE(IOFILE2,SPTYPE,1,INTDATA,
     '          REAL4DATA,REAL8DATA,CHARDATA,
     '          LOGDATA,SINTDATA,ERROR,*9999)
C              CALL BINWRITEFILE(IOFILE2,SPTYPE,1,INTDATA,
C     '          REAL(EMAP_SOCKFOCUS(nr)),REAL8DATA,CHARDATA,
C     '          LOGDATA,SINTDATA,ERROR,*9999)
            ENDIF
            
            IF(.NOT.ZERO_CROSSING.AND..NOT.PHI_OR_PHIH
     &        .AND..NOT.EXPORT_MFI) THEN
              IF(ELECTRODES) THEN
                EMAP_NUMELEC(nr)=EMAP_STOPELEC(nr)-EMAP_STARTELEC(nr)+1
              ELSE
                EMAP_NUMELEC(nr)=0
              ENDIF
              IF(LEADS) THEN
                EMAP_NUMLEADS(nr)=NUMLEADS
              ELSE
                EMAP_NUMLEADS(nr)=0
              ENDIF
              IF(RMSLEAD) THEN
                EMAP_NUMLEADS(nr)=EMAP_NUMLEADS(nr)+1
              ENDIF
            ELSEIF(ZERO_CROSSING) THEN

              CALL ASSERT(ISIZE_TBH(2).GT.0,
     '          'No electrodes to export',ERROR,*9999)
              EMAP_NUMELEC(nr)=ISIZE_TBH(2)
              EMAP_NUMLEADS(nr)=0
              EMAP_STARTELEC(nr)=1
              EMAP_STOPELEC(nr)=ISIZE_TBH(2)

              IF (NDT.EQ.EMAP_NUMELEC(nr)) THEN
                WRITE(OP_STRING,'(/'' >>WARNING: Wrong number of '
     '            //'data points for electrodes!'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF

            ELSEIF(PHI_OR_PHIH) THEN

              EMAP_NUMLEADS(nr)=0
              EMAP_STARTELEC(nr)=1
              IF(PHIH) THEN
                EMAP_NUMELEC(nr)=ISIZE_PHIH(1)
                EMAP_STOPELEC(nr)=ISIZE_PHIH(1)
              ELSE
                EMAP_NUMELEC(nr)=ISIZE_PHI(1)
                EMAP_STOPELEC(nr)=ISIZE_PHI(1)
              ENDIF

              IF(EVALUATE_PHI_NEAREST) THEN
                nrdata=0
              ELSE
                nrdata=1
              ENDIF
              
            ELSEIF(EXPORT_MFI) THEN
              EMAP_NUMLEADS(nr)=0
              EMAP_NUMELEC(nr)=ISIZE_MFI(1,nss)
              EMAP_STARTELEC(nr)=1
              EMAP_STOPELEC(nr)=ISIZE_MFI(1,nss)
            ENDIF

            EMAP_NUMDEVICES(nr)=EMAP_NUMELEC(nr)+EMAP_NUMLEADS(nr)

            WRITE(*,*) 'Num devices: ',EMAP_NUMDEVICES(nr)
            
            CALL BINWRITEFILE(IOFILE2,INTTYPE,1,EMAP_NUMDEVICES(nr),
     '        REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '        *9999)


            NUM_SIGNALS=NUM_SIGNALS+EMAP_NUMDEVICES(nr)
            DEVICE=0
            IF(ELECTRODES) THEN
              DO nelec=EMAP_STARTELEC(nr),EMAP_STOPELEC(nr)
                DEVICE=DEVICE+1
C*** Write out the device type
                INTDATA(1)=EMAP_ELECTRODE
                CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C*** Write the device number
                INTDATA(1)=DEVICE
                CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

C*** Write the device name
                IF(DEVICE_NAMES(1:7).EQ.'NUMBERS') THEN
                  !use the electrode number as the device name
                  WRITE(CHARDATA,'(I9)') nelec
                ELSEIF(DEVICE_NAMES(1:4).EQ.'CELL') THEN
                  !if possible, use the cell name
                  IF(nelec.LE.CELL_NUM_STATE(1)+CELL_NUM_DERIVED(1))
     '              THEN
                    !use variant 1 ???
                    WRITE(CHARDATA,'(A)') CELL_YQS_NAMES(nelec,1)
                    CALL STRING_TRIM(CHARDATA,IBEG,IEND)

C LKC cannot have zero-length strings                    
C                    IF (CHARDATA(IBEG:IEND).EQ.'') THEN
C                      WRITE(CHARDATA,'(I9)') nelec
C                    ENDIF
                    IF(IBEG.EQ.IEND) THEN
                      WRITE(CHARDATA,'(I9)') nelec
                    ENDIF
                    
                  ELSE
                    WRITE(CHARDATA,'(I9)') nelec
                  ENDIF
                ELSEIF(DEVICE_NAMES(1:4).EQ.'DATA') THEN
C LKC 29-NOV-2002 use the data labels rather than the
C                 sequential nd points
                  WRITE(CHARDATA,'(I9)') NDP(nelec)
                ELSE
                  ERROR='Invalid device name type'
                  GOTO 9999
                ENDIF

                CALL STRING_TRIM(CHARDATA,IBEG,IEND)
                INTDATA(1)=IEND-IBEG+1
                CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                CALL BINWRITEFILE(IOFILE2,CHARTYPE,INTDATA(1),INTDATA,
     '            REAL4DATA,REAL8DATA,CHARDATA(IBEG:IEND),LOGDATA,
     '            SINTDATA,ERROR,*9999)
C*** Write the device channel
C cpb 1/2/00 Use the device number as the channel number
C                INTDATA(1)=1
                INTDATA(1)=DEVICE
                CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

C*** Write the device dependent properties
C!!! cpb 17/3/97 Just use nr=1 to access NDDATA until data points are
C!!! truely region dependent.
C                nd=NDDATA(nelec,nr)
                IF(.NOT.ZERO_CROSSING.AND..NOT.PHI_OR_PHIH) THEN

C!!! 7-SEP-1999 LKC adding region dependence
C!!!  Needs to be changed when we have multiple region data
C!!!  A little trick to get multiple regions being exported
C                  nd=NDDATA(nelec,1)
C                  IF(nd.GT.0.AND.nd.LE.NDT) THEN
                  nd=nelec+NDDATA(1,nonr)-1


C LKC 7-SEP-1999 doesn't work if data points don't start at 1.
C  ie. multiregion signal files.
C  Can't really find an upper bound for the data
C
C                  IF(nd.GT.0.AND.nd.LE.NDDATA(0,nr)) THEN
                  IF(nd.GT.0) THEN
                    DO nj=1,NJT
                      REAL4DATA(nj)=REAL(ZD(nj,nd))
                    ENDDO !nj
                  ELSE
                    ERROR='>>Invalid electrode number'
                    GOTO 9999
                  ENDIF
                ELSEIF(ZERO_CROSSING) THEN

C*** Need the geometric coordinate of current ny of ZCROSSING matrix
C LKC 10-MAR-99 use data points from ZD as the electrode positions

                  DO nj=1,NJT
                    REAL4DATA(nj)=SNGL(ZD(nj,nelec))
                  ENDDO

                ELSE !PHI_OR_PHIH

C Not exporting PHI positions correctly at the moment
C GBS 8-Aug-2000 Might be doing so now

C LKC 28-NOV-2002 Better put some checks in!
                  CALL ASSERT(NDDATA(0,nrdata).GT.0,
     '              '>> No electrode positions to export',ERROR,*9999)


C LKC 24-JUN-2008 This array is not really reliable. It may have
C different information in it depending on what update was done
C last via UPPHI, UPDATA etc. Should be okay just to hard code some
C channel numbers in ....                  
C                  nd=NDDATA(nelec,nrdata)
                  nd=nelec
                  
                  DO nj=1,NJT
                    REAL4DATA(nj)=SNGL(ZD(nj,nd))
                  ENDDO
                ENDIF ! ZCROSSING or PHI_OR_PHIH

                CALL BINWRITEFILE(IOFILE2,SPTYPE,NJT,INTDATA,
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)
              ENDDO !nelec              
            ENDIF !electrodes

            
            IF(LEADS) THEN
              DO nlead=1,NUMLEADS
                DEVICE=DEVICE+1
C*** Write out the device type
                INTDATA(1)=EMAP_AUXILLIARY
                CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C*** Write the device number
                INTDATA(1)=DEVICE
                CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C*** Write the device name
                CHARDATA=LEADTITLE(nlead)
                CALL STRING_TRIM(CHARDATA,IBEG,IEND)
                INTDATA(1)=IEND-IBEG+1
                CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                CALL BINWRITEFILE(IOFILE2,CHARTYPE,INTDATA(1),INTDATA,
     '            REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '            ERROR,*9999)

C*** Write the device channel

C LKC 7-SEP-1990 - doesn't really matter but not really correct
C                INTDATA(1)=1
                INTDATA(1)=DEVICE
                CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C*** Write the device dependent properties
C             Nothing for auxillaries
              ENDDO !nlead
            ENDIF !leads
            IF(RMSLEAD) THEN
              DEVICE=DEVICE+1
C*** Write out the device type
              INTDATA(1)=EMAP_AUXILLIARY
              CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C*** Write the device number
              INTDATA(1)=DEVICE
              CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C*** Write the device name
              CHARDATA='RMS'
              CALL STRING_TRIM(CHARDATA,IBEG,IEND)
              INTDATA(1)=IEND-IBEG+1
              CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              CALL BINWRITEFILE(IOFILE2,CHARTYPE,INTDATA(1),INTDATA,
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
C*** Write the device channel
              INTDATA(1)=DEVICE
              CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C*** Write the device dependent properties
C             Nothing for auxillaries
            ENDIF
          ENDDO !nrr


          
C*** Write out the number of pages
          INTDATA(1)=0
          CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)

C*** Write out the signals

C*** Write out the number of signals
          INTDATA(1)=-NUM_SIGNALS
          CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C*** Write out the number of samples
          INTDATA(1)=NUM_SAMPLES
          CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C*** Write out the frequency
          REAL4DATA(1)=REAL(FREQUENCY)
          CALL BINWRITEFILE(IOFILE2,SPTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)


C LKC - what if the times don't start at 1? We need to assume that
C the original signal had times starting from 0.0          
C           DO nsample=1,NUM_SAMPLES
C             INTDATA(1)=nsample
C             CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
C      '        REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C         ENDDO !nsample
          
C*** Write out the times - offset by TSTART
          UNEMAP_TSTART=INT(TSTART*FREQUENCY)
C          UNEMAP_TSTART=1
          
          DO nsample=UNEMAP_TSTART+1,NUM_SAMPLES+UNEMAP_TSTART
            INTDATA(1)=nsample
            CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '        REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          ENDDO !nsample

          
C*** Write out the samples
          IF(.NOT.ZERO_CROSSING.AND..NOT.PHI_OR_PHIH
     &      .AND..NOT.EXPORT_MFI) THEN
            DO nsample=1,NUM_SAMPLES

C*** Read in the electrode data
              CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '          SIGNALMIN,TIME,WD,XID,ZD,'READ','BINARY',SIGNALFNAME,
     '          'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C*** Write out the electrode samples if necessary
              IF(ELECTRODES) THEN
                IF(RMSLEAD) RMSSUM=0.0d0
                DO nr=1,EMAP_NUMREGIONS_TOEXPORT

C LKC 7-SEP-1999 adding region export dependence
                  nonr=EXPORT_NRLIST(nr)

                  DO noelec=EMAP_STARTELEC(nr),EMAP_STOPELEC(nr)
C!!! cpb 17/3/97 Just use nr=1 to access NDDATA until data points are
C!!! truely region dependent.
C                  nelec=NDDATA(noelec,nr)

C LKC 7-SEP-1999 fixing up for multiple regions
C Assumes sequential data points
C                    nelec=NDDATA(noelec,1)

                    IF(nonr.EQ.1) THEN
                      nelec=noelec
                    ELSE
                      nelec=noelec+NDDATA(0,nonr-1)
                    ENDIF

                    REAL4DATA(1)=REAL(ZD(NJT+1,nelec))
                    IF(RMSLEAD) RMSSUM=RMSSUM+WD(NJT+1,nelec)*
     '                ZD(NJT+1,nelec)**2
                    CALL BINWRITEFILE(IOFILE2,SPTYPE,1,INTDATA,
     '                REAL4DATA,
     '                REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                  ENDDO !noelec
                ENDDO !nr
                IF(RMSLEAD) THEN
                  IF(DABS(RMSSUM).GT.ZERO_TOL) THEN
                    RMSSUM=DSQRT(RMSSUM)
                  ENDIF
                ENDIF
              ENDIF !electrodes
C*** Write out the lead samples if necessary
              IF(LEADS) THEN
                DO nlead=1,NUMLEADS
                  POTENTIAL=0.0d0
                  DO nleadelec=1,LEADELECS(0,nlead)
                    nelec=LEADELECS(nleadelec,nlead)
C                    CALL ASSERT(nelec.GT.0.AND.nelec.LE.NDT,
C     '                '>>Lead electrode not found',ERROR,*9999)
C CPB 15/02/00 Bug fix
C                    CALL ASSERT(nelec.GT.0.AND.
C     '                nelec.LT.NDDATA(NDDATA(0,nonr),nonr),
C     '                '>>Lead electrode not found',ERROR,*9999)
                    CALL ASSERT(nelec.GT.0.AND.
     '                nelec.LE.NDDATA(NDDATA(0,nonr),nonr),
     '                '>>Lead electrode not found',ERROR,*9999)
                    POTENTIAL=POTENTIAL+LEADCOUP(nleadelec,nlead)*
     '                ZD(NJT+1,nelec)
                  ENDDO !nleadelec
                  POTENTIAL=POTENTIAL+LEADCOUP(0,nlead)
                  REAL4DATA(1)=REAL(POTENTIAL)
                  CALL BINWRITEFILE(IOFILE2,SPTYPE,1,INTDATA,REAL4DATA,
     '              REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                ENDDO !nlead
              ENDIF !leads
              IF(RMSLEAD) THEN
C LKC 2-NOV-2000 RMSSUM needs casting
C                REAL4DATA(1)=RMSSUM
                REAL4DATA(1)=SNGL(RMSSUM)
                CALL BINWRITEFILE(IOFILE2,SPTYPE,1,INTDATA,REAL4DATA,
     '            REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              ENDIF
            ENDDO !nsample
          ELSEIF(ZERO_CROSSING) THEN !ZERO_CROSSING
            DO nsample=1,NUM_SAMPLES
C*** Write out the electrode samples if necessary
              IF(ELECTRODES) THEN
                DO nr=1,EMAP_NUMREGIONS_TOEXPORT
                  DO noelec=EMAP_STARTELEC(nr),EMAP_STOPELEC(nr)
                    REAL4DATA(1)=REAL(ZCROSSING(noelec,nsample))
                    CALL BINWRITEFILE(IOFILE2,SPTYPE,1,INTDATA,
     '                REAL4DATA,
     '                REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                  ENDDO !noelec
                ENDDO !nr
              ENDIF !electrodes
            ENDDO !nsample

C LKC 22-APR-1999
C*** Output the lead with the highest range
C    NOTE: ZCROSSING_RANGE is re-ordered
C*** Update the ZCROSSING_RANGE array
            DO noelec=1,ISIZE_TBH(2)
              DO ncol=1,3
                ZCROSSING_RANGE(noelec,ncol)=0.0d0
              ENDDO
            ENDDO

            nr=1
            DO noelec=1,ISIZE_TBH(2)
              DO nsample=1,NUM_SAMPLES !loop over time
                IF(ZCROSSING(noelec,nsample).GT.
     '            ZCROSSING_RANGE(noelec,1)) THEN
                  ZCROSSING_RANGE(noelec,1)=ZCROSSING(noelec,nsample)
                ENDIF !Calculate maximum value

                IF(ZCROSSING(noelec,nsample).LT.
     '            ZCROSSING_RANGE(noelec,2)) THEN
                  ZCROSSING_RANGE(noelec,2)=ZCROSSING(noelec,nsample)
                ENDIF !Calculate minimum value

              ENDDO !nts
            ENDDO

            DO noelec=1,ISIZE_TBH(2)
              ZCROSSING_RANGE(noelec,3)=
     '          ZCROSSING_RANGE(noelec,1)-ZCROSSING_RANGE(noelec,2)
            ENDDO !noelec

            DO noelec=1,ISIZE_TBH(2)
              ZCROSSING_INDEX(noelec)=noelec
            ENDDO
            noelec=ISIZE_TBH(2)
            CALL RSORT_REV(noelec,
     '        ZCROSSING_RANGE(1,3),ZCROSSING_INDEX)


            IF((ZCROSSING_INDEX(1).EQ.1).AND.
     '        (ZCROSSING_INDEX(2).EQ.2).AND.
     '        (ZCROSSING_INDEX(3).EQ.3)) THEN
              WRITE(OP_STRING,
     '          '(/'' >>WARNING: ZCROSSING already sorted?'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,
     '          '('' >>  ZCROSSING arrays need re-updating'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF


            noelec=MIN(20,EMAP_NUMELEC(nr))
            WRITE(OP_STRING,
     '        '(/'' Ordered Zero-Crossing electrodes: '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            CALL WRITE_LONG(INTTYPE,1,1,IOFI,noelec,5,5,ZCROSSING_INDEX,
     '        RDUMMY,'(1X,5I14)','(1X,5I14)',ERROR,*9999)
            WRITE(OP_STRING,'('' '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

            WRITE(OP_STRING,
     '        '(/'' Corresponding Ranges of: '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            CALL WRITE_LONG(DPTYPE,1,1,IOFI,noelec,5,5,IDUMMY,
     '        ZCROSSING_RANGE(1,3),
     '        '(1X,5D14.4)','(1X,5D14.4)',ERROR,*9999)


          ELSEIF(PHI_OR_PHIH) THEN !PHI or PHI_H

            DO nsample=1,NUM_SAMPLES
C*** Write out the electrode samples if necessary
              IF(ELECTRODES) THEN
                DO nr=1,EMAP_NUMREGIONS_TOEXPORT
                  DO noelec=EMAP_STARTELEC(nr),EMAP_STOPELEC(nr)
                    IF(PHIH) THEN
                      REAL4DATA(1)=REAL(PHI_H(noelec,nsample))
                    ELSE
                      REAL4DATA(1)=REAL(PHI(noelec,nsample))
                    ENDIF
                    CALL BINWRITEFILE(IOFILE2,SPTYPE,1,INTDATA,
     '                REAL4DATA,
     '                REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
                  ENDDO !noelec
                ENDDO !nr
              ENDIF !electrodes
            ENDDO !nsample

            
          ELSEIF(EXPORT_MFI) THEN
            WRITE(OP_STRING,
     &        '(''>>Note: Exporting MFI component '',I1)') COMPONENT
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)

C LKC 6-DEC-2006. Maybe this information should be added to
C the "gains" rather than scaling the actual values, esp since the
C values are stored single precision. I'm not sure.

            IF(TESLA) THEN
              
              DO nsample=1,NUM_SAMPLES
                DO nr=1,EMAP_NUMREGIONS_TOEXPORT
                  DO noelec=EMAP_STARTELEC(nr),EMAP_STOPELEC(nr)
                    REAL4DATA(1)=
     &                REAL(MFI(noelec,nsample,COMPONENT,nss)
     &                *4.D0*PI*1D-10)
                    CALL BINWRITEFILE(IOFILE2,SPTYPE,1,INTDATA,
     &                REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     &                ERROR,*9999)
                  ENDDO !noelec
                ENDDO !nr
              ENDDO !nsample
            ELSE
              
              DO nsample=1,NUM_SAMPLES
                DO nr=1,EMAP_NUMREGIONS_TOEXPORT
                  DO noelec=EMAP_STARTELEC(nr),EMAP_STOPELEC(nr)
                    REAL4DATA(1)=
     &                REAL(MFI(noelec,nsample,COMPONENT,nss))
                    CALL BINWRITEFILE(IOFILE2,SPTYPE,1,INTDATA,
     &                REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     &                ERROR,*9999)
                  ENDDO !noelec
                ENDDO !nr
              ENDDO !nsample

            ENDIF !TESLA
          ENDIF !EXPORT_MFI

C*** Write out the channel indicies, gains and offsets
          REAL4DATA(1)=0.0 !offset
          REAL4DATA(2)=1.0 !gain
          DO ndevice=0,NUM_SIGNALS-1
            INTDATA(1)=ndevice
            CALL BINWRITEFILE(IOFILE2,INTTYPE,1,INTDATA,REAL4DATA,
     '        REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
            CALL BINWRITEFILE(IOFILE2,SPTYPE,2,INTDATA,REAL4DATA,
     '        REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          ENDDO !ndevice
C*** Close the input signal file and binary output
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'CLOSE','BINARY',SIGNALFNAME,
     '      ' ',ENDFILE,.TRUE.,ERROR,*9999)
          CALL BINCLOSEFILE(IOFILE2,ERROR,*9999)



C*** Export signal to Cmgui (.exobj)
        ELSE IF(CMGUI) THEN !export to CMGUi

C*** Open files for signal trace and pointer
          CALL STRING_TRIM(FILE00,IBEG1,IEND1)
          CALL OPENF(IOFILE2,'DISK',FILE00(IBEG1:IEND1)//'.exgobj',
     '      'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)

C*** Create header for exgobj file
          OP_STRING(1)=NAME_EXGOBJ(1)      !name of gfxobj
          OP_STRING(2)='POLYLINE'          !type of object
          CALL WRITES(IOFILE2,OP_STRING,ERROR,*9999)
          WRITE(IOFILE2,'(I1)') 1          !don't know
          WRITE(IOFILE2,'(I1)') 1          !don't know
          OP_STRING(1)=MATERIAL_EXGOBJ(1)  !material
          CALL WRITES(IOFILE2,OP_STRING,ERROR,*9999)

C*** Transformation Matrix
C          DO nrow=1,4
C            WRITE(IOFILE2,'(4E12.5)')
C     '        (TRANSFORM_EXGOBJ(nrow,ncol),ncol=1,4)
C          ENDDO
          OP_STRING(1)='1 0 0 0'
          OP_STRING(2)='0 1 0 0'
          OP_STRING(3)='0 0 1 0'
C          OP_STRING(4)='0 0 0 1'
          CALL WRITES(IOFILE2,OP_STRING,ERROR,*9999)
C*** Displacement
          WRITE(IOFILE2,'(4E12.5)')
     '      (TRANSFORM_EXGOBJ(4,ncol),ncol=1,4)

C*** Line type
          OP_STRING(1)='PLAIN'
          CALL WRITES(IOFILE2,OP_STRING,ERROR,*9999)

          CALL ASSERT(EMAP_STARTELEC(1).LE.NDT,
     '      '>> Invalid electrode number',ERROR,*9999)

C*** Output the actual signal now
          WRITE(IOFILE2,'(I5)') NUM_SAMPLES !Num line segments
          DO nsample=1,NUM_SAMPLES

            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ','BINARY',SIGNALFNAME,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

            IF(LEADS) THEN
              nlead=ELEC_EXGOBJ
              POTENTIAL=0.0d0
              DO nleadelec=1,LEADELECS(0,nlead)
                nelec=LEADELECS(nleadelec,nlead)
                CALL ASSERT(nelec.GT.0.AND.nelec.LE.NDT,
     '            '>>Lead electrode not found',ERROR,*9999)
                POTENTIAL=POTENTIAL+LEADCOUP(nleadelec,nlead)*
     '            ZD(NJT+1,nelec)
              ENDDO !nleadelec
              POTENTIAL=POTENTIAL+LEADCOUP(0,nlead)
              IF(NJT.EQ.2) THEN
                WRITE(IOFILE2,'(E12.5,''  '',E12.5,''  '',E12.5)')
     '            nsample*TRANSFORM_EXGOBJ(1,1),
     '            POTENTIAL*TRANSFORM_EXGOBJ(3,3), 0.0
              ELSE
                WRITE(IOFILE2,'(E12.5,''  '',E12.5,''  '',E12.5)')
     '            nsample*TRANSFORM_EXGOBJ(1,1), 0.0,
     '            POTENTIAL*TRANSFORM_EXGOBJ(3,3)
              ENDIF
            ELSE
              IF(NJT.EQ.2) THEN
                WRITE(IOFILE2,'(E12.5,''  '',E12.5,''  '',E12.5)')
     '            nsample*TRANSFORM_EXGOBJ(1,1),
     '            ZD(NJT+1,ELEC_EXGOBJ)*TRANSFORM_EXGOBJ(3,3), 0.0
              ELSE
                WRITE(IOFILE2,'(E12.5,''  '',E12.5,''  '',E12.5)')
     '            nsample*TRANSFORM_EXGOBJ(1,1), 0.0,
     '            ZD(NJT+1,ELEC_EXGOBJ)*TRANSFORM_EXGOBJ(3,3)
              ENDIF
            ENDIF !leads

C*** Create and output the pointer over time
            WRITE(CHAR1,'(I4)') nsample+1000
            CALL STRING_TRIM(NAME_EXGOBJ(2),IBEG1,IEND1)
            CALL STRING_TRIM(CHAR1,IBEG,IEND)
            FILENAME=NAME_EXGOBJ(2)(IBEG1:IEND1)//'_'//CHAR1(IBEG:IEND)
     '        //'.exgobj'
            CALL STRING_TRIM(CHAR1,IBEG,IEND)
            CALL OPENF(IOFILE3,'DISK',FILENAME,
     '        'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
            OP_STRING(1)=NAME_EXGOBJ(2)          !name of object
            OP_STRING(2)='POLYLINE'              !type of object
            CALL WRITES(IOFILE3,OP_STRING,ERROR,*9999)
            WRITE(IOFILE3,'(I1)') 1              !don't know
            WRITE(IOFILE3,'(I1)') 1              !don't know
            OP_STRING(1)=MATERIAL_EXGOBJ(2)      !material
            CALL WRITES(IOFILE3,OP_STRING,ERROR,*9999)

C*** Transformation Matrix - incorporate scalings latter?
C            DO nrow=1,4
C              WRITE(IOFILE3,'(4E12.5)')
C     '          (TRANSFORM_EXGOBJ(nrow,ncol),ncol=1,4)
C            ENDDO
            OP_STRING(1)='1 0 0 0'
            OP_STRING(2)='0 1 0 0'
            OP_STRING(3)='0 0 1 0'
C            OP_STRING(4)='0 0 0 1'
            CALL WRITES(IOFILE3,OP_STRING,ERROR,*9999)
C*** Displacement
            WRITE(IOFILE3,'(4E12.5)')
     '        (TRANSFORM_EXGOBJ(4,ncol),ncol=1,4)

            OP_STRING(1)='PLAIN'          !Line type
            CALL WRITES(IOFILE3,OP_STRING,ERROR,*9999)
            WRITE(IOFILE3,'(I5)') 2       !Num line segments
            IF(NJT.EQ.2) THEN
              WRITE(IOFILE3,'(E12.5,''  '',E12.5,''  '',E12.5)')
     '          nsample*TRANSFORM_EXGOBJ(1,1), -30.0, 0.0
              WRITE(IOFILE3,'(E12.5,''  '',E12.5,''  '',E12.5)')
     '          nsample*TRANSFORM_EXGOBJ(1,1),  30.0, 0.0
            ELSE
              WRITE(IOFILE3,'(E12.5,''  '',E12.5,''  '',E12.5)')
     '          nsample*TRANSFORM_EXGOBJ(1,1), 0.0, -30.0
              WRITE(IOFILE3,'(E12.5,''  '',E12.5,''  '',E12.5)')
     '          nsample*TRANSFORM_EXGOBJ(1,1), 0.0,  30.0
            ENDIF

            CALL CLOSEF(IOFILE3,ERROR,*9999)
          ENDDO

C*** Close the input signal file and output file
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'CLOSE','BINARY',SIGNALFNAME,
     '      ' ',ENDFILE,.TRUE.,ERROR,*9999)
          CALL CLOSEF(IOFILE2,ERROR,*9999)


C*** Export to datafile
        ELSEIF(DATAFILE) THEN !export to datafile

          CALL STRING_TRIM(FILE00,IBEG1,IEND1)
          CALL OPENF(IOFILE1,'DISK',FILE00(IBEG1:IEND1)//'.dat',
     '      'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)

          nr=1

C LKC 31-JUL-2000
C moving this section to be be inside the if loop
C
CC MPN 13Mar2000: handled in loops below
CC          nelec=1
CC          CALL ASSERT(nelec.GE.0,'>> No electrodes to export ',
CC     '      ERROR,*9999)
C
C          IF(COMPSIGN) THEN
C            CALL STRING_TRIM(COMPFNAME,IBEG,IEND)
C            CALL STRING_TRIM(SIGNALFNAME,IBEG1,IEND1)
C
C            OP_STRING(1)=' Sample    Time        '
C     '        //COMPFNAME(IBEG:IEND)//'    '//SIGNALFNAME(IBEG1:IEND1)
C            CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)
C          ELSE
CC MPN 13Mar2000: added implied loop for multiple leads
C            WRITE(OP_STRING,'('' Sample    Time         '',
C     '        10(1X,A15),:/(1X,10(1X,A15)))')
C     '        (LEADTITLE(noelec),noelec=EMAP_STARTELEC(nr),
C     '        EMAP_STOPELEC(nr))
CC old       WRITE(OP_STRING,'('' Sample    Time         Electrode'')')
C            CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)
C          ENDIF

C*** Export signal to datafile



C LKC 20-SEP-2000 This is NOT handled for the previously existing
C  re-setting the variable
C
C MPN 13Mar2000: handled in loops below
C OLD          noelec=EMAP_STARTELEC(nr)
          noelec=EMAP_STARTELEC(nr)

          IF(.NOT.ZERO_CROSSING.AND..NOT.PHI_OR_PHIH) THEN
            IF(COMPSIGN) THEN
              CALL STRING_TRIM(COMPFNAME,IBEG,IEND)
              CALL STRING_TRIM(SIGNALFNAME,IBEG1,IEND1)

C LKC 30-SEP-2000 The filename were the wrong way around
C              OP_STRING(1)=' Sample    Time        '
C     '          //COMPFNAME(IBEG:IEND)//'    '//SIGNALFNAME(IBEG1:IEND1)
              OP_STRING(1)=' Sample    Time        '
     '          //SIGNALFNAME(IBEG1:IEND1)//'    '//COMPFNAME(IBEG:IEND)
              CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)
            ELSE
              IF(LEADS) THEN
                WRITE(OP_STRING,'('' Sample    Time         '','
     '            //'10(X,A15),:/(X,10(X,A15)))')
     '            (LEADTITLE(noelec),noelec=EMAP_STARTELEC(nr),
     '            EMAP_STOPELEC(nr))
              ELSE
                 WRITE(OP_STRING,'('' Sample    Time         '')')
              ENDIF
              CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)
            ENDIF

            DO nsample=1,NUM_SAMPLES
              CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ','BINARY',SIGNALFNAME,
     '          'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
              IF(COMPSIGN) THEN
                CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,1,
     '            SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'READ','BINARY',
     '            COMPFNAME,'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
                WRITE(OP_STRING,'(I5,D14.4,5D14.4)')
     '            nsample,time,(ZD(NJT+nj,noelec),nj=1,2)
                CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)
              ELSE
C MPN 13Mar2000: added implied loop for multiple leads
C                (assumes 1 elecrode per lead)
c                WRITE(OP_STRING,'(I5,D16.6,:10D16.6)')
c     '            nsample,time,(ZD(NJT+1,LEADELECS(1,noelec)),
c     '            noelec=EMAP_STARTELEC(nr),EMAP_STOPELEC(nr))
c                CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)
C MLB 9Aug2000: Put back old code as new code doesn't work
C               when you don't have any leads.
c                noelec=EMAP_STARTELEC(nr)
c                WRITE(OP_STRING,'(I5,D14.4,5D14.4)')
c     '            nsample,time,ZD(NJT+1,noelec)
c                CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)
C *** DPN 28 August 2000 - I need MPN's version, but I'll fix it for
C *** both of you!!
                IF(LEADS) THEN
                  ! Use MPN's method...
                  WRITE(OP_STRING,'(I5,D16.6,:10D16.6)')
     '              nsample,time,(ZD(NJT+1,LEADELECS(1,noelec)),
     '              noelec=EMAP_STARTELEC(nr),EMAP_STOPELEC(nr))
                ELSE
                  ! Use MLB's method...
                  noelec=EMAP_STARTELEC(nr)
                  ! DPN 19 May 2001 - increasing output precision
                  WRITE(OP_STRING,'(I5,E16.6,5E16.6)')
     '              nsample,time,ZD(NJT+1,noelec)
                ENDIF
                CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO !nsample

C*** Export zcrossing signals to datafile
C*** LKC 31-JUL-2000

          ELSEIF(ZERO_CROSSING) THEN

C LKC 2-AUG-2000 Adding the ability to export the zcrossing ranges
C  large ranges should correspond to critical points

            IF(SIGEXPO_TYPE.EQ.2) THEN

C*** Open .dat file
              CALL STRING_TRIM(FILE00,IBEG1,IEND1)
              CALL OPENF(IOFILE1,'DISK',FILE00(IBEG1:IEND1)//'.dat',
     '          'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)

              WRITE(OP_STRING,'(''#'')')
              CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,
     '          '(''# PointNum    NodeNum     Range ..'')')
              CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''#'')')
              CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)

              DO noelec=1,NPLIST3(0) ! heart nodes
                WRITE(OP_STRING,'(I5,I12,:10E16.5)')
     '            noelec,ZCROSSING_INDEX(noelec),
     '            ZCROSSING_RANGE(noelec,3)
                CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)
              ENDDO

C*** Close .dat file
              CALL CLOSEF(IOFILE1,ERROR,*9999)

            ELSEIF (SIGEXPO_TYPE.EQ.3) THEN
              CALL ASSERT(EMAP_STOPELEC(nr).LE.ISIZE_TBH(2),
     '          '>> Invalid number of electrodes for zcross',
     '          ERROR,*9999)
              CALL ASSERT(NPLIST3(0).EQ.ISIZE_TBH(2),
     '          '>> Invalid setup of NPLIST3 and TBH',ERROR,*9999)

C*** Note the output is currently limited to 10 - this is probably limited to
C***  132 characters in the file output.
              CALL ASSERT(EMAP_STOPELEC(nr)-EMAP_STARTELEC(nr)+1.LE.10,
     '          '>> Change WRITE formating for exporting > 10 elect'
     '          ,ERROR,*9999)


C*** Open .dat file
              CALL STRING_TRIM(FILE00,IBEG1,IEND1)
              CALL OPENF(IOFILE1,'DISK',FILE00(IBEG1:IEND1)//'.dat',
     '          'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)

              DO nonr=1,EMAP_NUMREGIONS_TOEXPORT !should always be 1 for now
                nr=EXPORT_NRLIST(nonr)

                WRITE(OP_STRING,'(''#'')')
                CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,
     '            '(''# Num     Time       Zcrossing ..'')')
                CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(''# Node        '',(10(1X,I11)))')
     '            (NPLIST3(noelec),noelec=EMAP_STARTELEC(nr),
     '            EMAP_STOPELEC(nr))
                CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(''#'')')
                CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)

                TIME=0.D0
                DO nsample=1,ISIZE_PHI(2) !time
                  WRITE(OP_STRING,'(I5,D14.4,:10E12.3)')
     '              nsample,time,(ZCROSSING(noelec,nsample),
     '              noelec=EMAP_STARTELEC(nr),EMAP_STOPELEC(nr))
                  CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)
                  TIME=TIME+1/CM_FREQUENCY
                ENDDO ! nsample
              ENDDO ! nonr

C*** Close .dat file
              CALL CLOSEF(IOFILE1,ERROR,*9999)
            ELSE
              ERROR='>> Unknown signal export type'
              GOTO 9999
            ENDIF ! SIGEXPO_TYPE
          ELSE
            ERROR='>> Unknown datafile export type'
            GOTO 9999
          ENDIF


C*** Close the input signal file and output file
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'CLOSE','BINARY',SIGNALFNAME,
     '      ' ',ENDFILE,.TRUE.,ERROR,*9999)
          CALL CLOSEF(IOFILE1,ERROR,*9999)

          IF(COMPSIGN) THEN
            CALL IOSIGN(IOFILE2,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'CLOSE','BINARY',COMPFNAME,
     '        ' ',ENDFILE,.TRUE.,ERROR,*9999)
          ENDIF

        ENDIF !emap

      ENDIF

      CALL EXITS('EXSIGN')
      RETURN
 9999 CALL ERRORS('EXSIGN',ERROR)

      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     &  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE','BINARY',SIGNALFNAME,
     &  ' ',ENDFILE,.TRUE.,ERROR,*9998)
      IF(ISBINFILEOPEN(IOFILE2))
     &  CALL BINARYCLOSEFILE(IOFILE2,ERR,CERROR)

 9998 CALL EXITS('EXSIGN')
      RETURN 1
      END


