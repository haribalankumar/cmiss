      SUBROUTINE SETUP(iw,ERROR,*)

C#### Subroutine: SETUP
C###  Description:
C###    SETUP performs setup operations for workstation iw.

C**** Note that default workstation type is used (GWSDEF=41 on
C**** microvax). If using microvax from another workstation use
C**** DCL command Define GKS$WSTYPE to set appropriate value.
C****
C**** IWKDEF(0) is number of defined workstations
C**** IWKDEF(noiw),noiw=1,IWKDEF(0) is list of defined workstations
C**** IWKG(iw) is 0 for nongraphics window (eg menu)
C****    "        1 for graphics output window
C**** IWKS(iw) is 0 when workstation iw is not defined
C****    "        1  "      "        iw is defined but not active
C****    "        2  "      "        iw is defined and active
C**** IWKT(iw) is 1 for GX or GKS workstation
C****    "        2 for PHIGS workstation
C****    "        3 for Frame-grabber display screen
C****
C**** Workstation_id=1 for x,y viewport in 2D or 3D
C****      "         2  "  y,z    "      "  "
C****      "         3  "  x,z    "      "  "
C****      "         4  "  map viewport (3D only)
C****      "         5  "  DT2651 frame-buffer 1
C****      "         6  "  DT2651 frame-buffer 2
C****      "         7  "  choice viewport
C****      "         8  "  help documentation viewport
C****      "         9  "  PHIGS 3D plot viewport (now elect. lead)
C****      "        10  "  nodal time history viewport
C****      "        11  "  sections viewport
C****      "        12  "  fibre angle distribution (3D only)
C****      "        13  "  sheet angle distribution (3D only)
C****      "        15  "  GKS postscript file
C****      "        16  "  Phigs postscript file & .POST file
C****      "        17  "  Metafile
C****      "        18  "  x,y viewport in 3D for chords
C****      "        19  "  z,y viewport in 3D for chords
C****      "        20  "  record output
C****      "        21  "  ode output
C****      "        22  "  ode output
C****      "        23  "  ode output
C****      "        24  "  ode output
C****      "        31  "  stress/strain profile viewport
C****      "        32  "  stress/strain profile viewport
C****      "        33  "  fibre angle profile viewport
C****      "        34  "  signals
C****      "        35  "  signal trace
C****      "        36  "  trace choice viewport
C****      "        40  "  s-plane plot
C****      "        41  "  Bode plot
C****      "        42  "  Nyquist plot
C****      "        43  "  transient response plot
C****      "        44  "  ionic current model - voltage
C****      "        45  "  ionic current model - gates
C****      "        46  "  ionic current model - currents
C****      "        47  "  ionic current model - concentrations
C****      "        50  "  Matrix decomposition window
C****      "        51  "  Homogeneous strain window
C****      "        55  "  String input window
C****      "        56  "  String input window
C****      "        60  "  Spreadsheet window
C****      "        61  "  Spreadsheet GKS-plot window
C****      "        62  "  Spreadsheet GKS-plot window
C****      "        63  "  Spreadsheet GKS-plot window
C****      "        64  "  Spreadsheet GKS-plot window
C****      "        65  "  Spreadsheet GKS-plot window
C****      "        66  "  Spreadsheet GKS-plot window
C****      "        67  "  Spreadsheet PHIGS-plot window
C****      "        68  "  Spreadsheet PHIGS-plot window
C****      "        69  "  Spreadsheet GKS-plot window
C****      "        70  "  choice viewport (main spreadsheet menu)
C****      "        71  "  choice viewport (3D rotation etc)
C****      "        72  "  choice viewport (page scroll/col menu)
C****      "        73  "  choice viewport (parameter display)
C****      "        74  "  choice viewport (various)
C****      "        75  "  choice viewport (ODE choice menu)
C****      "        76  "  choice viewport (s-plane choice menu)
C****      "        77  "  choice viewport
C****      "        78  "  choice viewport
C****      "        79  "  choice viewport
C****      "        81  "  valuator
C****      "        82  "  valuator
C****      "        83  "  valuator for 1st window
C****      "        84  "  valuator  "   "     "
C****      "        85  "  valuator  "  2nd    "
C****      "        86  "  valuator  "   "     "
C****      "        87  "  valuator  "  3rd    "
C****      "        88  "  valuator  "   "     "
C****      "        89  "  valuator  "   "     "
C****      "        91  "  choice viewport (window 1 menu; row bar plot)
C****      "        92  "  choice viewport (window 2 menu; histogram)
C****      "        93  "  choice viewport (window 3 menu; line/area plot)
C****      "        94  "  choice viewport (window 4 menu; phase plane)
C****      "        95  "  choice viewport (window 5 menu; point/bar plot)
C****      "        96  "  choice viewport (2D scatter plot)
C****      "        97  "  choice viewport (GKS_DRAW)
C****      "        98  "  choice viewport (GKS_DRAW)
C****      "        99  "  choice viewport (GKS_DRAW)

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'disp00.cmn'
      INCLUDE 'gks000.cmn'
      INCLUDE 'graf00.cmn'
!     Parameter List
      INTEGER iw
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER iiw,noiw
      REAL WINRECT(4),WLDRECT(4)
      INTEGER ERR

      CALL ENTERS('SETUP',*9999)

      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' iw='',I3,'' IWKS(iw)='',I2,'' IWKT(iw)='',I2,'
     '    //''' IWKG(iw)='',I2)') iw,IWKS(iw),IWKT(iw),IWKG(iw)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(iw.EQ.16.OR.iw.EQ.67.OR.iw.EQ.68) THEN !DPB. 18feb95 got
                                                !rid of iw.eq.3
        IWKT(iw)=2 !indicates iw is PHIGS workstation
      ELSE IF(iw.EQ.5.OR.iw.EQ.6) THEN
        IWKT(iw)=3 !indicates iw is Frame-grabber screen
      ELSE
        IWKT(iw)=1 !indicates iw is GKS   workstation
      ENDIF

      IF(IWKT(iw).EQ.1) THEN      !GKS
        CALL CHECK_GX_OPEN(ERROR,*9999)
      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
      ENDIF


      IF(IW.EQ.1.AND.IWKS(1).EQ.0) THEN ! XY
        IWKS(1)=1 !indicates window open
        IWKG(1)=1 !indicates graphics output window
        IWKT(1)=1 !indicates GKS workstation
        WINRECT(1)=0.78
        WINRECT(2)=0.25
        WINRECT(3)=1.18
        WINRECT(4)=0.65
        WLDRECT(1)=XMIN
        WLDRECT(2)=YMIN
        WLDRECT(3)=XMAX
        WLDRECT(4)=YMAX
      ELSE IF(IW.EQ.2.AND.IWKS(2).EQ.0) THEN ! YZ
        IWKS(2)=1 !indicates window open
        IWKG(2)=1 !indicates graphics output window
        IWKT(2)=1 !indicates GKS workstation
        WINRECT(1)=0.78
        WINRECT(2)=0.25
        WINRECT(3)=1.18
        WINRECT(4)=0.65
        WLDRECT(1)=YMIN
        WLDRECT(2)=ZMIN
        WLDRECT(3)=YMAX
        WLDRECT(4)=ZMAX
      ELSE IF(IW.EQ.3.AND.IWKS(3).EQ.0) THEN ! XZ
        IWKS(3)=1 !indicates window open
        IWKG(3)=1 !indicates graphics output window
        IWKT(3)=1 !indicates GKS workstation
        WINRECT(1)=0.78
        WINRECT(2)=0.25
        WINRECT(3)=1.18
        WINRECT(4)=0.65
        WLDRECT(1)=XMIN
        WLDRECT(2)=ZMIN
        WLDRECT(3)=XMAX
        WLDRECT(4)=ZMAX
      ELSE IF(IW.EQ.7.AND.IWKS(7).EQ.0) THEN
        IWKS(7)=1 !indicates window open
        IWKG(7)=1 !indicates graphics output window
        IWKT(7)=1 !indicates GKS workstation
        WINRECT(1)=0.78
        WINRECT(2)=0.25
        WINRECT(3)=1.18
        WINRECT(4)=0.65
        WLDRECT(1)=0.0
        WLDRECT(2)=0.0
        WLDRECT(3)=1.0
        WLDRECT(4)=1.0
      ELSE IF(IW.EQ.9.AND.IWKS(9).EQ.0) THEN
        IWKS(9)=1 !indicates window open
        IWKG(9)=1 !indicates graphics output window
        IWKT(9)=1 !indicates GKS workstation
        WINRECT(1)=0.38
        WINRECT(2)=0.25
        WINRECT(3)=1.18
        WINRECT(4)=0.65
        WLDRECT(1)=-0.14
        WLDRECT(2)=-0.15
        WLDRECT(3)=1.1
        WLDRECT(4)=1.1
      ELSE IF(IW.EQ.10.AND.IWKS(10).EQ.0) THEN
        IWKS(10)=1 !indicates window open
        IWKG(10)=1 !indicates graphics output window
        IWKT(10)=1 !indicates GKS workstation
        WINRECT(1)=0.38
        WINRECT(2)=0.25
        WINRECT(3)=1.18
        WINRECT(4)=0.65
        WLDRECT(1)=-0.15
        WLDRECT(2)=-0.3
        WLDRECT(3)=1.1
        WLDRECT(4)=1.1
      ELSE IF(IW.EQ.11.AND.IWKS(11).EQ.0) THEN
        IWKS(11)=1 !indicates window open
        IWKG(11)=1 !indicates graphics output window
        IWKT(11)=1 !indicates GKS workstation
        WINRECT(1)=0.78
        WINRECT(2)=0.25
        WINRECT(3)=1.18
        WINRECT(4)=0.65
        WLDRECT(1)=-0.14
        WLDRECT(2)=1.1
        WLDRECT(3)=-2.3
        WLDRECT(4)=2.3
      ELSE IF(IW.EQ.12.AND.IWKS(12).EQ.0) THEN
        IWKS(12)=1 !indicates window open
        IWKG(12)=1 !indicates graphics output window
        IWKT(12)=1 !indicates GKS workstation
        WINRECT(1)=0.78
        WINRECT(2)=0.25
        WINRECT(3)=1.18
        WINRECT(4)=0.65
        WLDRECT(1)=-0.20
        WLDRECT(2)=-1.1*REAL(PI)/2.0
        WLDRECT(3)=1.1
        WLDRECT(4)=1.1*REAL(PI)/2.0
      ELSE IF(IW.EQ.13.AND.IWKS(13).EQ.0) THEN
        IWKS(13)=1 !indicates window open
        IWKG(13)=1 !indicates graphics output window
        IWKT(13)=1 !indicates GKS workstation
        WINRECT(1)=0.78
        WINRECT(2)=0.25
        WINRECT(3)=1.18
        WINRECT(4)=0.65
        WLDRECT(1)=XMIN
        WLDRECT(2)=-0.2*(XMAX-XMIN)
        WLDRECT(3)=XMAX
        WLDRECT(4)=0.8*(XMAX-XMIN)
      ENDIF

      IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
        CALL OPWIN(IW,WINRECT,ERR)
        IF(err.ne.0) THEN
          ERROR='>>Error from GX'
          GOTO 9999
        ENDIF
        CALL SETWIN(IW,ERR)
        IF(err.ne.0) THEN
          ERROR='>>Error from GX'
          GOTO 9999
        ENDIF
        CALL SETVIEW(WLDRECT,ERR)
        IF(err.ne.0) THEN
          ERROR='>>Error from GX'
          GOTO 9999
        ENDIF
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,'('' Reset to:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' IW='',I3,'' IWKS(iw)='',I2,'
     '    //''' IWKT(iw)='',I2,'' IWKG(iw)='',I2)')
     '    IW,IWKS(iw),IWKT(iw),IWKG(iw)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      noiw=0
      DO iiw=1,99      !to record defined workstations in IWKDEF
        IF(IWKS(iiw).GT.0) THEN !iiw is defined
          noiw=noiw+1
          IWKDEF(noiw)=iiw
        ENDIF
      ENDDO
      IWKDEF(0)=noiw

      CALL EXITS('SETUP')
      RETURN
 9999 CALL ERRORS('SETUP',ERROR)
      CALL EXITS('SETUP')
      RETURN 1
      END


