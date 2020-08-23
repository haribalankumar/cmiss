      SUBROUTINE DILEAD(ISEG,ISLEAD,LD,LEAD_LIST,NBJ,NDDATA,NRLIST,
     '  WD,XID,ZD,CSEG,STRING,ERROR,*)

C#### Subroutine: DILEAD
C###  Description:
C###    DILEAD displays time history of electrocardiographic leads.

C**** Note: ISLEAD(0) is segment containing axes and labels.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lead00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER ISEG(*),ISLEAD(0:NUMLEADMX),LD(NDM),LEAD_LIST(0:NLISTM),
     '  NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),NRLIST(0:NRM)
      REAL*8 WD(NJM,NDM),XID(NIM,NDM),ZD(NJM,NDM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER AXES_INDEX,IBEG,IEND,IFROMC,
     '  INDEX_POLYLINE,lead,LINE_INDEX,N1LIST,N3CO,nelec,nlead,
     '  nleadelec,nolead,nosig,nosign,nr,NSIGNALS,NTSIG,NUMTIMEDATA,
     '  PLOTDEPTH,PLOTHEIGHT,PLOTORIENT,PLOTWIDTH,PS_TYPE,tag
      REAL*8 FACTOR,POTENTIAL,POTMIN,POTMAX,SIGNALMAX(9),
     '  SIGNALMIN(9),RFROMC,TEND,TIME,TMIN,TMAX,TSTART,VMAX,VMIN,
     '  XSIG(8192),YSIG(8192)
      CHARACTER AXES_COLOUR*10,AXES_STYLE*8,AXES_WIDTH*8,
     '  ERROR1*255,FILE*100,FILEFORMAT*6,
     '  LINE_COLOUR*10,LINE_STYLE*8,LINE_WIDTH*8,PRINTEXT*4,
     '  PRINTFILENAME*80,PRINTFILETYPE*10,PTITLE*100,SIG_LIST(5)*32
      LOGICAL ABBREV,ALL_REGIONS,AXESCOLOURSET,AXESSTYLESET,
     '  AXESWIDTHSET,CALCPOTMAX,CALCPOTMIN,
     '  CBBREV,ENDFILE,INLIST,LINECOLOURSET,LINESTYLESET,
     '  LINEWIDTHSET,PLOTAXESSCALE,PLOTLABELLEAD,PLOTTITLE,
     '  PRINTFILE

      CALL ENTERS('DILEAD',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C MPN 21Mar2002: added lots more options.

C#### Command: FEM display leads
C###  Parameter:       <leads (#s/all)[all]>
C###  Parameter:       <signal FILE(s)[$current]>
C###  Parameter:       <(ascii/binary)[ascii]>
C###  Parameter:       <tstart #[beginning]>
C###  Parameter:       <tend #[end]>
C###  Parameter:       <vmin #[signalmin]>
C###  Parameter:       <vmax #[signalmax]>
C###  Parameter:       <scale #[1.0]>
C###  Parameter:       <(title NAME/notitle)[SignalPlot]>
C###  Parameter:       <(labellead/nolabellead)[labellead]>
C###  Parameter:       <rgb=RGB[red]>
C###  Parameter:       <linestyle=STYLE[solid]>
C###  Parameter:       <linewidth=WIDTH[width1]>
C###  Parameter:       <axesrgb=RGB[black]>
C###  Parameter:       <axesstyle=STYLE[solid]>
C###  Parameter:       <axeswidth=WIDTH[width1]>
C###  Parameter:       <(axesscale/noaxesscale)[axesscale]>
C###  Parameter:       <printfile <portable/postscript[portable]><plotwidth #[400]><plotheight #[400]>>
C###  Parameter:       <region #[1]>
C###  Description:
C###    <HTML>
C###      Displays a set of predefined leads (DEFINE LEAD) with GX.
C###      <UL>
C###      <LI> SIGNALS specifies the BINARY/ASCII signal file(s) from
C###            which to read the specified LEADS.
C###      <LI> TSTART/TEND specify the minimum/maximum values for the
C###            time axis and defaults to the min/max time values for
C###            the signal
C###      <LI> VMIN/VMAX specify the minimum/maximum values for the
C###            potential axis and defaults to the min/max potential
C###            values for the signal
C###      <LI> SCALE scales the signals by a specified magnitude.
C###      <LI> TITLE/NOTITLE species whether a plot title NAME will
C###            be displayed
C###      <LI> NOLABELLEAD does not plot the name of the lead
C###      <LI> The lead from each file will be drawn using a
C###            different RGB colour from that specified (red) and
C###            using the specified LINESTYLE and LINEWIDTH.
C###      <LI> The axes will be drawn using the specified AXESRGB
C###           colour, AXESSTYLE and AXESWIDTH, and NOAXESSCALE can
C###           be used to omit the numbers from the axes.
C###      <LI> PRINTFILE will print the plot to a appropriately named
C###           PORTABLE pixmap or POSTscript file with specified
C###           PLOTWIDTH and/or PLOTHEIGHT.
C###      </UL>
C###    </HTML>

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<leads (#s/ALL)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<signal FILE(s)[$current]>'
        OP_STRING(4)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(5)=BLANK(1:15)//'<tstart #[beginning]>'
        OP_STRING(6)=BLANK(1:15)//'<tend #[end]>'
        OP_STRING(7)=BLANK(1:15)//'<vmin #[signalmin]>'
        OP_STRING(8)=BLANK(1:15)//'<vmax #[signalmax]>'
        OP_STRING(9)=BLANK(1:15)//'<scale #[1.0]>'
        OP_STRING(10)=BLANK(1:15)//'<(title NAME/notitle)[SignalPlot]>'
        OP_STRING(11)=BLANK(1:15)//'<(labellead/nolabellead)'
     '    //'[labellead]>'
        OP_STRING(12)=BLANK(1:15)//'<rgb=RGB[red]>'
        OP_STRING(13)=BLANK(1:15)//'<linestyle=STYLE[solid]>'
        OP_STRING(14)=BLANK(1:15)//'<linewidth=WIDTH[width1]>'
        OP_STRING(15)=BLANK(1:15)//'<axesrgb=RGB[black]>'
        OP_STRING(16)=BLANK(1:15)//'<axesstyle=STYLE[solid]>'
        OP_STRING(17)=BLANK(1:15)//'<axeswidth=WIDTH[width1]>'
        OP_STRING(18)=BLANK(1:15)//'<(axesscale/noaxesscale)'
     '    //'[axesscale]>'
        OP_STRING(19)=BLANK(1:15)//'<printfile '
     '    //'<portable/postscript[portable]>'
     '    //'<plotwidth #[400]><plotheight #[400]>>'
        OP_STRING(20)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','DILEAD',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL ASSERT(CALL_LEAD,'>> Define Leads first',ERROR,*9999)
        CALL ASSERT(NJT+1.LE.NJM,'>>Increase NJM',ERROR,*9999)

C LKC Note: 23-JUL-98 no true nr dependency yet
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '    ERROR,*9999)
        nr=NRLIST(1)

C*** List of Signal files
        IF(CBBREV(CO,'SIGNAL',3,noco+1,NTCO,N3CO)) THEN
          CALL PARSSL(CO(N3CO+1),32,NSIGNALS,SIG_LIST,ERROR,*9999)
        ELSE
C need to set the default
          CALL STRING_TRIM(FILE00,IBEG,IEND)
          SIG_LIST(1)=FILE00(IBEG:IEND)
          NSIGNALS=1
        ENDIF

        IF(CBBREV(CO,'LEADS',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'ALL',1)) THEN
            LEAD_LIST(0)=NUMLEADS
            DO nlead=1,NUMLEADS
              LEAD_LIST(nlead)=nlead
            ENDDO !nlead
          ELSE
            CALL PARSIL(CO(N3CO+1),NLISTM,LEAD_LIST(0),
     '        LEAD_LIST(1),ERROR,*9999)
            DO nlead=1,LEAD_LIST(0)
              lead=LEAD_LIST(nlead)
              CALL ASSERT(lead.GT.0.AND.lead.LE.NUMLEADS,
     '          '>>Lead number not defined',ERROR,*9999)
            ENDDO !nlead
          ENDIF
        ELSE
          LEAD_LIST(0)=NUMLEADS
          DO nlead=1,NUMLEADS
            LEAD_LIST(nlead)=nlead
          ENDDO !nlead
        ENDIF

C LKC 23-JUL-98 removed
C        IF(NTCOQU(noco).EQ.1) THEN
C          FILE=COQU(noco,1)
C        ELSE
C          FILE=SIGFNAME
C        ENDIF
C        CALL STRING_TRIM(FILE,IBEG,IEND)

        IF(CBBREV(CO,'BINARY',1,noco+1,NTCO,N3CO)) THEN
          FILEFORMAT='BINARY'
        ELSE
          FILEFORMAT='ASCII'
        ENDIF

        IF(CBBREV(CO,'SCALE',2,noco+1,NTCO,N3CO)) THEN
          FACTOR=RFROMC(CO(N3CO+1))
        ELSE
          FACTOR=1.0d0
        ENDIF

        IF(CBBREV(CO,'TSTART',2,noco+1,NTCO,N3CO)) THEN
          TSTART=RFROMC(CO(N3CO+1))
        ELSE
          TSTART=-RMAX
        ENDIF

        IF(CBBREV(CO,'TEND',2,noco+1,NTCO,N3CO)) THEN
          TEND=RFROMC(CO(N3CO+1))
        ELSE
          TEND=RMAX
        ENDIF

C MPN 21Mar2002: added lots more options.
        IF(CBBREV(CO,'VMAX',3,noco+1,NTCO,N3CO)) THEN
          VMAX=RFROMC(CO(N3CO+1))
          CALCPOTMAX=.FALSE.
        ELSE
          VMAX=-RMAX
          CALCPOTMAX=.TRUE.
        ENDIF

        IF(CBBREV(CO,'VMIN',3,noco+1,NTCO,N3CO)) THEN
          VMIN=RFROMC(CO(N3CO+1))
          CALCPOTMIN=.FALSE.
        ELSE
          VMIN=RMAX
          CALCPOTMIN=.TRUE.
        ENDIF

        IF(CBBREV(CO,'TITLE',3,noco+1,NTCO,N3CO)) THEN
          PLOTTITLE=.TRUE.
          PTITLE=CO(N3CO+1)(1:100)
        ELSE IF(CBBREV(CO,'NOTITLE',3,noco+1,NTCO,N3CO)) THEN
          PLOTTITLE=.FALSE.
          PTITLE=' '
        ELSE
          PLOTTITLE=.TRUE.
          PTITLE='SignalPlot'
        ENDIF

        IF(CBBREV(CO,'NOLABELLEAD',3,noco+1,NTCO,N3CO)) THEN
          PLOTLABELLEAD=.FALSE.
        ELSE
          PLOTLABELLEAD=.TRUE.
        ENDIF

        IF(CBBREV(CO,'NOAXESSCALE',3,noco+1,NTCO,N3CO)) THEN
          PLOTAXESSCALE=.FALSE.
        ELSE
          PLOTAXESSCALE=.TRUE.
        ENDIF

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          LINE_COLOUR=CO(N3CO+1)(IBEG:IEND)
          LINECOLOURSET=.TRUE.
        ELSE
          LINE_COLOUR='RED'
          LINECOLOURSET=.FALSE.
        ENDIF
        IF(CBBREV(CO,'LINESTYLE',5,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          LINE_STYLE=CO(N3CO+1)(IBEG:IEND)
          LINESTYLESET=.TRUE.
        ELSE
          LINE_STYLE='SOLID'
          LINESTYLESET=.FALSE.
        ENDIF
        IF(CBBREV(CO,'LINEWIDTH',5,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          LINE_WIDTH=CO(N3CO+1)(IBEG:IEND)
          LINEWIDTHSET=.TRUE.
        ELSE
          LINE_WIDTH='WIDTH1'
          LINEWIDTHSET=.FALSE.
        ENDIF
        LINE_INDEX=INDEX_POLYLINE(0,LINE_STYLE,LINE_WIDTH,LINE_COLOUR)

        IF(CBBREV(CO,'AXESRGB',5,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          AXES_COLOUR=CO(N3CO+1)(IBEG:IEND)
          AXESCOLOURSET=.TRUE.
        ELSE
          AXES_COLOUR='BLACK'
          AXESCOLOURSET=.FALSE.
        ENDIF
        IF(CBBREV(CO,'AXESSTYLE',5,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          AXES_STYLE=CO(N3CO+1)(IBEG:IEND)
          AXESSTYLESET=.TRUE.
        ELSE
          AXES_STYLE='SOLID'
          AXESSTYLESET=.FALSE.
        ENDIF
        IF(CBBREV(CO,'AXESWIDTH',5,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          AXES_WIDTH=CO(N3CO+1)(IBEG:IEND)
          AXESWIDTHSET=.TRUE.
        ELSE
          AXES_WIDTH='WIDTH1'
          AXESWIDTHSET=.FALSE.
        ENDIF
        AXES_INDEX=INDEX_POLYLINE(0,AXES_STYLE,AXES_WIDTH,AXES_COLOUR)

        IF(CBBREV(CO,'PRINTFILE',5,noco+1,NTCO,N3CO)) THEN
          PRINTFILE=.TRUE.
          IF(CBBREV(CO,'POSTSCRIPT',4,noco+1,NTCO,N3CO)) THEN
            PRINTFILETYPE='POSTSCRIPT'
            PS_TYPE=2 !1=ps; 2=eps??; BUT unused in GX?
            PRINTEXT='.eps'
          ELSE
            PRINTFILETYPE='PORTABLE'
            PRINTEXT='.ppm'
          ENDIF
          IF(CBBREV(CO,'PLOTWIDTH',5,noco+1,NTCO,N3CO)) THEN
            PLOTWIDTH=IFROMC(CO(N3CO+1))
          ELSE
            PLOTWIDTH=400
          ENDIF
          IF(CBBREV(CO,'PLOTHEIGHT',5,noco+1,NTCO,N3CO)) THEN
            PLOTHEIGHT=IFROMC(CO(N3CO+1))
          ELSE
            PLOTHEIGHT=400
          ENDIF
          PLOTDEPTH=1  !1=colour 2=monochrome
          PLOTORIENT=1 !1=portrait 2=landscape
        ELSE
          PRINTFILE=.FALSE.
        ENDIF

C MPN 21Mar2002: end

        POTMAX=VMAX
        POTMIN=VMIN
        TMAX=-RMAX
        TMIN=RMAX
        TIME=0.D0

C***    Find out max and min potentials for the scale.
        DO nosign=1,NSIGNALS
          FILE=SIG_LIST(nosign)
C          CALL STRING_TRIM(FILE,IBEG,IEND)

C***      Open signal file and read in electrode geometry etc.
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,FILE,
     '      'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,FILE,
     '      'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

C***      Calculate the number of time data points
          IF(FILEFORMAT(1:5).EQ.'ASCII') THEN
            NUMTIMEDATA=0
            ENDFILE=.FALSE.
            DO WHILE(.NOT.ENDFILE)
              CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '          SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,FILE,
     '          'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
              IF(.NOT.ENDFILE) NUMTIMEDATA=NUMTIMEDATA+1
            ENDDO
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,FILE,
     '        'RESET',ENDFILE,.TRUE.,ERROR,*9999)
          ENDIF

          DO nolead=1,LEAD_LIST(0)
            nlead=LEAD_LIST(nolead)
            DO nleadelec=1,LEADELECS(0,nlead)
              nelec=LEADELECS(nleadelec,nlead)
              CALL ASSERT(INLIST(nelec,NDDATA(1,nr),NDDATA(0,nr),
     '          N1LIST),
     '          '>>Lead electrode not found in region electrode list',
     '          ERROR,*9999)
            ENDDO !nleadelec
          ENDDO !nlead

          DO tag=1,NUMTIMEDATA
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,FILE,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

            IF(TIME.GE.TSTART.AND.TIME.LE.TEND) THEN
              DO nolead=1,LEAD_LIST(0)
                nlead=LEAD_LIST(nolead)
                POTENTIAL=0.0d0
                DO nleadelec=1,LEADELECS(0,nlead)
                  nelec=LEADELECS(nleadelec,nlead)
                  POTENTIAL=POTENTIAL+LEADCOUP(nleadelec,nlead)*
     '              ZD(NJT+1,nelec)
                ENDDO !nleadelec
                POTENTIAL=POTENTIAL+LEADCOUP(0,nlead)
C MPN 21Mar2002: added lots more options.
                IF(CALCPOTMAX.AND.POTENTIAL.GT.POTMAX) POTMAX=POTENTIAL
                IF(CALCPOTMIN.AND.POTENTIAL.LT.POTMIN) POTMIN=POTENTIAL
C old            IF(POTENTIAL.GT.POTMAX) POTMAX=POTENTIAL
C old            IF(POTENTIAL.LT.POTMIN) POTMIN=POTENTIAL
              ENDDO !nlead
              IF(TIME.GT.TMAX) TMAX=TIME
              IF(TIME.LT.TMIN) TMIN=TIME
            ENDIF
          ENDDO !tag
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,FILE,
     '      ' ',ENDFILE,.TRUE.,ERROR,*9999)

        ENDDO !nosign

C***    Open window
        CALL ACWK(9,0,ERROR,*9999)

C***    Calculate and draw signals
        DO nosign=1,NSIGNALS
          FILE=SIG_LIST(nosign)
C          CALL STRING_TRIM(FILE,IBEG,IEND)

C***      Open signal file and read in electrode geometry etc.
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,FILE,
     '      'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,FILE,
     '      'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

          DO nolead=1,LEAD_LIST(0)
            nlead=LEAD_LIST(nolead)

            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,FILE,
     '        'RESET',ENDFILE,.TRUE.,ERROR,*9999)

            nosig=0
            DO tag=1,NUMTIMEDATA

              CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '          SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,FILE,
     '          'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

              IF(TIME.GE.TSTART.AND.TIME.LE.TEND) THEN
                nosig=nosig+1
                CALL ASSERT(nosig.LE.8192,'>>Too many time points',
     '            ERROR,*9999) !1 kHz sampling freq. for 8 seconds
                XSIG(nosig)=(TIME-DABS(TMIN))/(DABS(TMAX)+DABS(TMIN))
                POTENTIAL=0.0d0
                DO nleadelec=1,LEADELECS(0,nlead)
                  nelec=LEADELECS(nleadelec,nlead)
                  POTENTIAL=POTENTIAL+LEADCOUP(nleadelec,nlead)*
     '              ZD(NJT+1,nelec)
                ENDDO !nleadelec
                POTENTIAL=POTENTIAL+LEADCOUP(0,nlead)
                YSIG(nosig)=POTENTIAL
              ENDIF
            ENDDO !tag
            NTSIG=nosig

C LKC change for plotting multiple signal files
C            CALL SGLEAD(LINE_INDEX,ISEG,ISLEAD(nolead),nlead,nolead,
C     '        NTSIG,FACTOR,XSIG,YSIG,TMAX,TMIN,POTMAX,POTMIN,CSEG,
C     '        PTITLE,PLOTAXESSCALE,PLOTLABELLEAD,PLOTTITLE,ERROR,*9999)
            CALL SGLEAD(LINE_INDEX,ISEG,ISLEAD(nolead+nosign-1),
     '        nlead,nolead,
     '        NTSIG,FACTOR,XSIG,YSIG,TMAX,TMIN,POTMAX,POTMIN,CSEG,
     '        PTITLE,PLOTAXESSCALE,PLOTLABELLEAD,PLOTTITLE,ERROR,*9999)

C LKC 24-JUL-98 what does this do ? Redundant
C            CALL BINSETFILE(IOFILE1,0,ERROR,*9999)
C            CALL BINSKIPFILE(IOFILE1,NUMSKIPBYTES,ERROR,*9999)

          ENDDO !nlead
          LINE_INDEX=LINE_INDEX+4 !adjust the line colour
                                  !(sets of 3 refer POLYLINE)

          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,FILE,
     '      ' ',ENDFILE,.TRUE.,ERROR,*9999)

        ENDDO !nosign

C*** Draw axes and labels (if required)
        CALL SGLEAD(AXES_INDEX,ISEG,ISLEAD(0),0,0,NTSIG,FACTOR,
     '    XSIG,YSIG,TMAX,TMIN,POTMAX,POTMIN,CSEG,
     '    PTITLE,PLOTAXESSCALE,PLOTLABELLEAD,PLOTTITLE,ERROR,*9999)

C MPN 21Mar2002: added lots more options.
        IF(PRINTFILE) THEN
          CALL STRING_TRIM(FILE,IBEG,IEND)
          PRINTFILENAME=FILE(IBEG:IEND)//PRINTEXT(1:4)
          CALL PRINT_SCREEN_TO_FILE(PLOTDEPTH,PLOTHEIGHT,PLOTORIENT,
     '      PS_TYPE,PLOTWIDTH,PRINTFILENAME,PRINTFILETYPE,ERROR,*9999)
        ENDIF

        CALL DAWK(9,0,ERROR,*9999)

      ENDIF

      CALL EXITS('DILEAD')
      RETURN
 9999 CALL ERRORS('DILEAD',ERROR)
      CLOSE(UNIT=9)
C LKC 24-JUL-98 Close file properly
C      IF(ISBINFILEOPEN(IOFILE1))
C     '  CALL BINARYCLOSEFILE(IOFILE1,ERR,CERROR)
C      RETURN 1
      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,FILE,
     '  ' ',ENDFILE,.TRUE.,ERROR1,*9998)
 9998 CALL EXITS('DILEAD')
      RETURN 1
      END


