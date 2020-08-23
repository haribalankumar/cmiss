C#### Module: CMISS_ARCHIVE1
C###  Description:
C###    Contains archived code from modules FE00 -> FE09

CFE00 Subroutine AINOUT   handle input/output of answers (yes/no)
C     Subroutine AINPUT   read answer values (yes/no)
C     Subroutine CINOUT   handle input/output of character strings
C     Subroutine CINPUT   read character strings
C     Subroutine IINOUT   handle i/o of integer values
C     Subroutine IINOUT1  handle i/o of integer values (backup)
C     Subroutine IINPUT   read integer values
C     Subroutine INOUT    handle prompting input/output
C     Subroutine IOPOTN   handle input/output of potential data  
C     Subroutine IPGEOM     read previous version of iod files
C     Subroutine LINOUT   handle i/o of logical values
C     Subroutine LINPUT   read logical values
C     Subroutine RINOUT   handle input/output of real values
C     Subroutine RINOUT1  handle input/output of real values (backup)
C     Subroutine RINPUT   read real values
C     Subroutine RINPUT1  read real values (with file i/p provision)
C###  Routine: SOCKET_READ_SIGNAL read fid. marker times from socket
C     Subroutine WRITE_POLYLINE writes polyline data to socket 
C###  Routine: WRITE_SIGNAL writes signal data to socket

CFE01 Function AFFIRM      true if string argument is 'true' or 'yes'
C     Function BINSEARCH   perform a binary search on an ordered list 
C     Function CFROML      character string from logical variable
C     Function CI1         number of characters for formatting         
C     Function EXIST       true if file name in argument exists 
C     Function MAXD        position of abs max of real*8 array
C     Function MAXIM       position of abs maximum of REAL*8 array
C     Function MAXIM2      position of maximum of REAL*8 array
C     Function MINIM       position of abs minimum of REAL*8 array
C     Function MINIM2      position of minimum of REAL*8 array
C     Function NEGATE      true if argument is false
C     Function NKNU        deriv index nk given partial deriv index nu
C     Function QUAD        performs quadrature  
C###  Routine: CALC_FID   calculate the position of fiducial markers
C     Subroutine CMS       replaces CMS routine on VAX         
C     Subroutine CONSTANTS displays values for physical constants
C     Subroutine FREELU    returns first available logical unit number 
C     Subroutine FREETIMER   frees a timer handle
C     Subroutine GETTIMER    returns a handle for use with timer
C     Subroutine HELP      display help information on screen
C###  Routine: IMAGE_READ  reads  unformatted image file 512x512
C###  Routine: IMAGE_WRITE  writes unformatted image file 512x512
C     Subroutine MAP4      map coords onto vp 4 (cartography projection)
C     Subroutine MULT      calc inner product of two matrices       
C     Subroutine PARSIA    put nested item list into integer array
C     Subroutine PARSLO    put item in string into logical variable
C     Subroutine POLY      eval complex poly using Horners algorithm
C     Subroutine POST      convert linefeeds to returns in ps files  
C     Subroutine PUTRAN    write out transformation matrix to file     
C     Subroutine RESETTIMERS resets all timer handles for later use
C     Subroutine ROOTSC    finds complex roots of polynomial
C     Subroutine SCALE0    update transformation matrix by scale factor
C###  Routine: SEGPAC     pack segments down to elim cancelled segments
C###  Routine: SIGNAL_READ reads  unformatted signal file 128*2048
C     Subroutine STACK     Controls operations defined in command stack      
C     Subroutine UNITS     gives unit conversions
C     Subroutine ZERO_SPARSE_MATRIX  Zeroes a row/entire matrix
      
CFE02
C###  Routine: CALC_IMAGE_DATA define reduced data set ZD for image
C     Subroutine CALC_STRIPE_PTS calc SPAMM stripe data points
C     Subroutine CALC_NVNEP calc the NVNE & NVNP soln mapping array
C     Subroutine CPCP2    subdivide control point array  
C###  Routine: FEMREINI reinitialise all variables & arrays
C     Subroutine GAUSS3   Gauss coords & wgts for B-spline elements
C     Subroutine GAUSSPWB Gauss coords & wgts using NAG 
C     Subroutine GLOBALH    calc various global arrays from elem arrays
C     Subroutine SIMPS    Simpson's rule calc of arc length (D Rahdert)
C     Subroutine SPLINE   calc f.e. B-spline polynomial coefficients   



CFE03 Function   NYPJK      global dof assoc with deriv,opt.var & node
C     Function   NYF        global dof for Fourier basis
C     Function   NYPHK      global dof assoc with deriv,opt.var & node
C     Subroutine CLOS12     find closest Xi pt for 1D cyl.polar elements
C     Subroutine CLOS13     find closest Xi pt for 1D sph.polar elements
C     Subroutine CLOS14     find closest Xi pt for 1D prolate elements
C     Subroutine CLOS22     find closest Xi pt for 2D cyl.polar elements
C     Subroutine CLOS23     find closest Xi pt for 2D sph.polar elements 
C     Subroutine CONT_PLANE_X finds contour intersection xi points
C     Subroutine FITGAU     fit field variables to Gauss point data
C     Subroutine FITGEO     fit geometric variables to data points
C     Subroutine FITFOU     fit fourier parameters to time data
C     Subroutine IOHESS     I/O of Hessian matrix for nonlinear opt.n
C     Subroutine IPDATA     read data points from keyboard
C     Subroutine LNLXI      find elements surrounding each element
C     Subroutine NEWXID     calculate Xi coords of data point projection
C     Subroutine NY_NO      calc optimisation dofs no from mesh dofs ny
C     Subroutine PROJ21     find closest Xi proj pt for 2D r.c. elements
C     Subroutine PROJ_FUNCT function evaluation for proj21
C     Subroutine PROJ_HESS  hessian evaluation for proj21
C     Subroutine PROJ_MONIT monitor function evaluation for proj21
C     Subroutine XPES       eval elem stiffness matrix in least sqs fit
C     Subroutine ZDERF      eval elem load vector for l.s. fourier fit
C     Subroutine ZDESF      ditto for least sqs fourier fit

CFE04 
C     Subroutine DIVH1    Subdivide Lagrange/Hermite tensor product element
C     Subroutine BASIS6   Singular Boundary Element basis function routine
C     Subroutine BASIS7   Inputs param for Lagrange/Hermite tensor prod 
C                basis funct & Gauss-Legendre or Gauss-Lobatto quad.

CFE05
C     Function PS1        eval 1D linear B-spline polynomial coeffs
C     Function PS2        eval 1D quadratic B-spline polynomial coeffs
C     Function PS3        eval 1D cubic B-spline polynomial coeffs
C     Function PSE        eval tensor product B-spline basis functions 
C###  Routine: PSI3     (fn) eval tensor product B-spline basis functions 
C###  Routine: PSL      (fn) eval lin,quadratic or cubic spline poly terms
    
CFE06
C     Subroutine SHEET3
       
                
CFE07 Subroutine ASSEMBLE3  Assem matrices for Time dependent FEM probs
C     Subroutine ASSEMBLE8  Assem matrices for Laplace transform probs
C     Subroutine FUNCT1     define fn to minimize (called by NAG E04JBF)
C     Subroutine FUNCT3     define fn to zero (called by NAG C05NCF)
C     Subroutine GEOMIN     nonlinear geometric fitting
C     Subroutine MARCH1     solve parabolic  eqns by time integration
C     Subroutine MARCH4     cardiac activ.n problems (variable step)
C     Subroutine MARCH19    solver for cellular models
C     Subroutine MODAL      finds eigenvalues & eigenvectors
C     Subroutine MONIT      monitor minimisation (called by NAG E04JBF) 
C     Subroutine SOLVE1     (old) solve symm, pos. def. system of eqtns
C     Subroutine SOLVE1C    solve symm, pos. def. complex eqtns
C     Subroutine SOLVE2     solve unsymm sparse linear equations
C     Subroutine SOLVE3     solve unsymm sparse linear equations
C     Subroutine SOLVE5     solve symm sparse linear equations
C     Subroutine SOLVE6     solve complex eqtns for freq domain analysis

CFE09
C###  Routine: PYTHAG   (fn) calc dist to pt without under/overflow
C     Function TIMER       calls a 'c' routine to get system times

Module FE00
=========== 

      SUBROUTINE AINOUT(IO_TYPE,IFREE,IFIXED,A_FORMAT,NDATA,A_DATA,
     '  A_DEFLT,INFO,ERROR,*)

C#### Subroutine: AINOUT
C###  Description:
C###    AINOUT handle input/output of answers (yes/no).
C**** If IO_TYPE=0 prompts are written to IFREE
C****  "    "    1   "      "  written to IFREE & (with data) to IFIXED
C****  "    "    2   "      "  read with data from IFIXED
C****  "    "    3   "      "  written with data to IFIXED
C****  "    "    4   "      "  with data from IFIXED & to IFREE
C****  "    "    5   default values are used
C**** If IFIXED is not open no IO is performed on that file.
C**** Answers separated by commas or blanks are read into
C**** A_DATA(NDAT) NDAT=1,NDATA.
C**** If an end of file or carriage return is detected while attempting
C**** to read A_DATA(NDAT) it is assigned the value A_DEFLT(NDAT).
C**** If an error is detected a diagnostic message is returned in ERROR
C**** and control is returned to the statement number of the asterisk.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
!     Parameter List
      INTEGER IFIXED,IFREE,INFO,IO_TYPE,NDATA
      CHARACTER A_FORMAT*(*),A_DATA(*),A_DEFLT(*),ERROR*(*)
!     Local Variables
      INTEGER CLOCAT,IBEG,IEND,IOSTAT,IPOS,irec,JPOS,m,MAXTRY,n,nd
      CHARACTER CHAR*1,CIOSTA*5,CIREC*5,CIUNIT*10,LINE*132,SUBSTR*1
      LOGICAL REVERT,OPENED 
      DATA MAXTRY/8/

      CALL ENTERS('AINOUT',*9999)
      IF(NDATA.GT.IOAM) THEN
        ERROR=' NDATA exceeds the length of the array A_DATA'
        GOTO 9999
      ENDIF
      REVERT=.FALSE.
 10   IF(IO_TYPE.LE.1.OR.REVERT) THEN
        DO m=1,MAXTRY
          CALL WRITE_LINE(IOOP,A_FORMAT,ERROR,*9999)
          CALL AINPUT(IFREE,NDATA,A_DATA,A_DEFLT,INFO,ERROR,*9999)
          IF(ERROR.EQ.' ') GOTO 2
          IF(ERROR(1:2).NE.' ?') THEN
            WRITE(OP_STRING,*) ERROR
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO
        ERROR=' More than 8 invalid input attempts made'
        GOTO 9999
    2   IF(IO_TYPE.EQ.1) THEN
          INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
          IF(OPENED.AND.IFIXED.NE.0) THEN
            WRITE(UNIT=IFIXED,FMT=A_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '        (A_DATA(n),n=1,NDATA)
            IF(IOSTAT.EQ.0) THEN
              ERROR=' '
            ELSE
              ERROR=' Character data write error'
              GOTO 9999
            ENDIF
          ELSE
            ERROR=' '
          ENDIF
        ENDIF

      ELSE IF(IO_TYPE.EQ.2.OR.IO_TYPE.EQ.4) THEN
        IF(DOP) THEN
          WRITE(OP_STRING,'(A)') A_FORMAT(1:80)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(OPENED) THEN
 3        READ(UNIT=IFIXED,FMT='(A)',REC=irec,IOSTAT=IOSTAT) LINE
          IF(LINE(1:1).EQ.'!') THEN   !This line is a comment
            IF(LINE(2:2).EQ.'!') THEN !Do not print
            ELSE                      !Print comment line
              WRITE(OP_STRING,'(A)') LINE
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            irec=irec+1
            GO TO 3
          ENDIF
          SUBSTR=''''
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' SUBSTR=',SUBSTR
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IPOS=CLOCAT(SUBSTR,A_FORMAT)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' IPOS=',IPOS
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CHAR=' '
          IF(IPOS.GT.0) THEN
            SUBSTR='/'
            JPOS=CLOCAT(SUBSTR,A_FORMAT(1:4))
            IF(JPOS.GT.0) THEN
              READ(UNIT=IFIXED,FMT='(/A,'''//A_FORMAT(IPOS+2:),
     '          REC=irec,IOSTAT=IOSTAT) CHAR,(A_DATA(n),n=1,NDATA)
            ELSE
              READ(UNIT=IFIXED,FMT='(A,'''//A_FORMAT(IPOS+2:),
     '          REC=irec,IOSTAT=IOSTAT) CHAR,(A_DATA(n),n=1,NDATA)
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' CHAR='',A)') CHAR
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE
            READ(UNIT=IFIXED,FMT=A_FORMAT,
     '        REC=irec,IOSTAT=IOSTAT) (A_DATA(n),n=1,NDATA)
          ENDIF
          IF(IOSTAT.EQ.0) THEN
            ERROR=' '
          ELSE IF(IOSTAT.EQ.36) THEN
            ERROR='End of file'
            GO TO 9999
          ELSE
            ERROR=' Character data read error'
            GOTO 9999
          ENDIF
          IF(CHAR.EQ.'*') THEN
            REVERT=.TRUE.
            GO TO 10
          ENDIF
        ELSE
          WRITE(UNIT=CIUNIT,FMT='(I10)') IFIXED
          CALL TRIM(CIUNIT,IBEG,IEND)
          ERROR=' Unit '//CIUNIT(IBEG:IEND)//' is not open'
          GOTO 9999
        ENDIF
        IF(IO_TYPE.EQ.4) THEN
          WRITE(OP_STRING,A_FORMAT) (A_DATA(n),n=1,NDATA)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(IO_TYPE.EQ.3) THEN !write data to file
        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(OPENED.AND.IFIXED.NE.0) THEN
          WRITE(UNIT=IFIXED,FMT=A_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '      (A_DATA(n),n=1,NDATA)
          IF(IOSTAT.EQ.0) THEN
            ERROR=' '
          ELSE
            ERROR=' Character data write error'
            GOTO 9999
          ENDIF
        ELSE
          ERROR=' '
        ENDIF

      ELSE IF(IO_TYPE.EQ.5) THEN !set to defaults
        DO nd=1,NDATA
          A_DATA(nd)=A_DEFLT(nd)
        ENDDO
      ENDIF

      CALL EXITS('AINOUT')
      RETURN

 9999 CALL TRIM(ERROR,IBEG,IEND)
      IF(ERROR(1:IEND).EQ.'End of file') THEN
        IO_TYPE=1
        WRITE(OP_STRING,'('' End of file: revert to interactive '','
     '    //'''input'')')
      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GOTO 10
      ELSE IF(ERROR(1:4).EQ.'Edit') THEN
        ERROR=ERROR(1:IEND)//'>AINOUT'
        CLOSE(UNIT=IFIXED)
      ELSE IF(ERROR(1:7).EQ.'Restart') THEN
        ERROR=ERROR(1:IEND)//'>AINOUT'
        CLOSE(UNIT=IFIXED,STATUS='DELETE')
      ELSE IF(ERROR(1:14).EQ.'File not found') THEN
        ERROR='Help file not found>AINOUT'
      ELSE
        WRITE(UNIT=CIOSTA,FMT='(I5)') IOSTAT
        WRITE(UNIT=CIREC ,FMT='(I5)') irec
        ERROR=ERROR(1:IEND)//': IOSTAT='//CIOSTA(1:5)//' irec='//
     '    CIREC(1:5)//' A_FORMAT='//A_FORMAT(1:10)//'>AINOUT'
        CLOSE(UNIT=IFIXED)
      ENDIF
      CALL EXITS('AINOUT')
      RETURN 1
      END


      SUBROUTINE AINPUT(IUNIT,NDATA,ADATA,ADEFLT,INFO,ERROR,*)

C#### Subroutine: AINPUT
C###  Description:
C###    AINPUT reads answer values (yes/no) from IUNIT.
C**** IUNIT has a logical record length of IRECL.
C**** Answers separated by blanks or commas are read into
C**** ADATA(NDAT) NDAT=1,NDATA.
C**** If an end of file or carriage return is detected while attempting
C**** to read ADATA(NDAT) it is assigned the value ADEFLT(NDAT).
C**** No range check is made on the default values.
C**** If help request is made (?) a help file (based on INFO) is opened.
C**** If stop request is made (STOP/stop/S/s) the program is stopped.
C**** If an edit request is made (EDIT/edit/E/e) control returns to
C**** the command line and the file is closed.
C**** If restart request is made (RESTART/restart/R/r) control returns 
C**** to the command line and the file is deleted.
C**** If an error is detected or a help request is made ERROR is
C**** returned with a diagnostic message.
C**** NOTE: The Q edit descriptor in the read format is nonstandard!

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER INFO,IUNIT,NDATA
      CHARACTER ADATA(*),ADEFLT(*),ERROR*(*)
!     Local Variables
      INTEGER ic,ICOL,IEND,IOSTAT,irec,ISTART,NCHAR,nd,NDAT
      CHARACTER CDATA*500,CHAR*5,CLINE*500,FILE*6

      CALL ENTERS('AINPUT',*9999)
      NDAT=0
      DO 4 irec=1,NDATA
        CALL READ_LINE(IUNIT,IOSTAT,NCHAR,CLINE,ERROR,*9999)
        IF(IOSTAT.EQ.0.AND.NCHAR.GT.0) THEN
          ICOL=1
          DO 3 nd=NDAT+1,NDATA
            ISTART=0
            DO 1 ic=ICOL,NCHAR
              IF((CLINE(ic:ic).EQ.' ').OR.(CLINE(ic:ic).EQ.',')) THEN
                IF(ISTART.GT.0) THEN
                  IEND=ic-1
                  GOTO 2
                ENDIF
              ELSE
                IF(ISTART.EQ.0) THEN
                  ISTART=ic
                ENDIF
              ENDIF
    1       CONTINUE
            IEND=NCHAR
    2       IF(ISTART.GT.0) THEN
              ICOL=IEND+2
              CDATA=CLINE(ISTART:IEND)
              IF((CDATA.EQ.'EDIT').OR.
     '           (CDATA.EQ.'edit').OR.
     '           (CDATA.EQ.'E').OR.
     '           (CDATA.EQ.'e'))THEN
                ERROR='Edit'
                GO TO 9999
              ELSE IF((CDATA.EQ.'RESTART').OR.
     '                (CDATA.EQ.'restart').OR.
     '                (CDATA.EQ.'R').OR.
     '                (CDATA.EQ.'r'))THEN
                ERROR='Restart'
                GO TO 9999
              ELSE IF((CDATA.EQ.'STOP').OR.
     '                (CDATA.EQ.'stop').OR.
     '                (CDATA.EQ.'S').OR.
     '                (CDATA.EQ.'s'))THEN
                CLOSE(UNIT=IUNIT)
                STOP
              ELSE IF(CDATA.EQ.'?') THEN
                WRITE(UNIT=CHAR,FMT='(I5)') INFO
                FILE='FE'//CHAR(1:4)
                CALL DOCUM(FILE,'doc','Parameter number '//CHAR(5:5),
     '            ERROR,*9997)
 9997           IF(ERROR.EQ.' ') THEN
                  ERROR=' ?'
                  GO TO 9998
                ELSE IF(ERROR(1:14).EQ.'File not found') THEN
                  WRITE(OP_STRING,*) '>>Help not available'
      		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ERROR=' ?'
                  GO TO 9998
                ELSE
                  GO TO 9999
                ENDIF
              ELSE
                IF((CDATA.EQ.'YES').OR.
     '             (CDATA.EQ.'yes').OR.
     '             (CDATA.EQ.'Y').OR.
     '             (CDATA.EQ.'y')) THEN
                  NDAT=NDAT+1
                  ADATA(NDAT)='Y'
                  IF(NDAT.EQ.NDATA) THEN
                    ERROR=' '
                    GO TO 9998
                  ENDIF
                ELSE IF((CDATA.EQ.'NO').OR.
     '                  (CDATA.EQ.'no').OR.
     '                  (CDATA.EQ.'N').OR.
     '                  (CDATA.EQ.'n')) THEN
                  NDAT=NDAT+1
                  ADATA(NDAT)='N'
                  IF(NDAT.EQ.NDATA) THEN
                    ERROR=' '
                    GO TO 9998
                  ENDIF
                ELSE
                  ERROR=' Answer data range error'
                  GO TO 9998
                ENDIF
              ENDIF
            ELSE
              GOTO 4
            ENDIF
    3     CONTINUE
        ELSE IF(IOSTAT.EQ.0.AND.NCHAR.EQ.0) THEN
          NDAT=NDAT+1
          DO nd=1,NDATA
            ADATA(nd)=ADEFLT(nd)
          ENDDO
          ERROR=' '
          GO TO 9998
        ELSE
          REWIND(UNIT=IUNIT)
          ERROR=' Character data read error'
          GO TO 9998
        ENDIF
    4 CONTINUE

      ERROR=' '
 9998 CALL EXITS('AINPUT')
      RETURN
 9999 CALL ERRORS('AINPUT',ERROR)
      CALL EXITS('AINPUT')
      RETURN 1
      END


      SUBROUTINE CINOUT(IO_TYPE,IFREE,IFIXED,C_FORMAT,NDATA,C_DATA,
     '  C_DEFLT,LENGTH,INFO,ERROR,*)

C#### Subroutine: CINOUT
C###  Description:
C###    CINOUT handle input/output of character strings.
C**** If IO_TYPE=0 prompts are written to IFREE
C****  "    "    1   "      "  written to IFREE and (with data) to IFIXED
C****  "    "    2   "      "  read with data from IFIXED
C****  "    "    3   "      "  written with data to IFIXED
C****  "    "    4   "      "  read with data from IFIXED and written to IFREE
C****  "    "    5   default values are used
C**** If IFIXED is not open no IO is performed on that file.
C**** Strings separated by commas are read into C_DATA(NDAT) NDAT=1,NDATA
C**** Strings are truncated after LENGTH characters.
C**** If an end of file or carriage return is detected while attempting
C**** to read C_DATA(NDAT) it is assigned the value C_DEFLT(NDAT).
C**** If an error is detected a diagnostic message is returned in ERROR
C**** and control is returned to the statment number of the asterisk.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
!     Parameter List
      INTEGER IFIXED,IFREE,INFO,IO_TYPE,LENGTH,NDATA
      CHARACTER C_FORMAT*(*),C_DATA(*)*(*),C_DEFLT(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER CLOCAT,IBEG,IEND,IOSTAT,IPOS,irec,JPOS,m,MAXTRY,n,
     '  NCHAR,nd
      CHARACTER CHAR*1,CIOSTA*5,CIREC*5,CIUNIT*10,LINE*132,SUBSTR*1
      LOGICAL OPENED,REVERT
      DATA MAXTRY/8/

      CALL ENTERS('CINOUT',*9999)
      IF(NDATA.GT.IOCM) THEN
        ERROR=' NDATA exceeds the length of the array C_DATA'
        GOTO 9999
      ENDIF
      REVERT=.FALSE.
 10   IF(IO_TYPE.LE.1.OR.REVERT) THEN
        DO 1 m=1,MAXTRY
          CALL WRITE_LINE(IOOP,C_FORMAT,ERROR,*9999)
C ajp      WRITE(OP_STRING,C_FORMAT)
C 17-5-93 CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL CINPUT(IFREE,NCHAR,NDATA,C_DATA,C_DEFLT,INFO,ERROR,*9999)
          IF(ERROR.EQ.' ') GOTO 2
          IF(ERROR(1:2).NE.' ?') THEN
            WRITE(OP_STRING,*) ERROR
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
    1   CONTINUE
        ERROR=' More than 8 invalid input attempts made'
        GOTO 9999
    2   IF(IO_TYPE.EQ.1) THEN
          INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
          IF(OPENED.AND.IFIXED.NE.0) THEN
            IF(NCHAR.EQ.0) NCHAR=LENGTH
            WRITE(UNIT=IFIXED,FMT=C_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '        (C_DATA(n)(1:NCHAR),n=1,NDATA)
            IF(IOSTAT.EQ.0) THEN
              ERROR=' '
            ELSE
              ERROR=' Character data write error'
              GOTO 9999
            ENDIF
          ELSE
            ERROR=' '
          ENDIF
        ENDIF

      ELSE IF(IO_TYPE.EQ.2.OR.IO_TYPE.EQ.4) THEN
        IF(DOP) THEN
          WRITE(OP_STRING,'(A)') C_FORMAT(1:80)
      	  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(OPENED) THEN
 3        READ(UNIT=IFIXED,FMT='(A)',REC=irec,IOSTAT=IOSTAT) LINE
          IF(LINE(1:1).EQ.'!') THEN   !This line is a comment
            IF(LINE(2:2).EQ.'!') THEN !Do not print
            ELSE                      !Print comment line
              WRITE(OP_STRING,'(A)') LINE
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            irec=irec+1
            GO TO 3
          ENDIF
          SUBSTR=''''
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' SUBSTR=',SUBSTR
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IPOS=CLOCAT(SUBSTR,C_FORMAT)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' IPOS=',IPOS
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CHAR=' '
          IF(IPOS.GT.0) THEN
            SUBSTR='/'
            JPOS=CLOCAT(SUBSTR,C_FORMAT(1:4))
            IF(JPOS.GT.0) THEN
              READ(UNIT=IFIXED,FMT='(/A,'''//C_FORMAT(IPOS+2:),
     '          REC=irec,IOSTAT=IOSTAT) CHAR,(C_DATA(n)(1:LENGTH),
     '          n=1,NDATA)
            ELSE
              READ(UNIT=IFIXED,FMT='(A,'''//C_FORMAT(IPOS+2:),
     '          REC=irec,IOSTAT=IOSTAT) CHAR,(C_DATA(n)(1:LENGTH),
     '          n=1,NDATA)
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' CHAR='',A)') CHAR
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE
            READ(UNIT=IFIXED,FMT=C_FORMAT,
     '        REC=irec,IOSTAT=IOSTAT) (C_DATA(n)(1:LENGTH),n=1,NDATA)
          ENDIF
          IF(IOSTAT.EQ.0) THEN
            ERROR=' '
          ELSE IF(IOSTAT.EQ.36) THEN
            ERROR='End of file'
            GO TO 9999
          ELSE
            ERROR=' Character data read error'
            GOTO 9999
          ENDIF
          IF(CHAR.EQ.'*') THEN
            REVERT=.TRUE.
            GO TO 10
          ENDIF
        ELSE
          WRITE(UNIT=CIUNIT,FMT='(I10)') IFIXED
          CALL TRIM(CIUNIT,IBEG,IEND)
          ERROR=' Unit '//CIUNIT(IBEG:IEND)//' is not open'
          GOTO 9999
        ENDIF
        IF(IO_TYPE.EQ.4) THEN
          WRITE(OP_STRING,C_FORMAT) (C_DATA(n)(1:LENGTH),n=1,NDATA)
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(IO_TYPE.EQ.3) THEN !write data to file
        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(OPENED.AND.IFIXED.NE.0) THEN
          NCHAR=LENGTH
          WRITE(UNIT=IFIXED,FMT=C_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '         (C_DATA(n)(1:NCHAR),n=1,NDATA)
          IF(IOSTAT.EQ.0) THEN
            ERROR=' '
          ELSE
            ERROR=' Character data write error'
            GOTO 9999
          ENDIF
        ELSE
          ERROR=' '
        ENDIF

      ELSE IF(IO_TYPE.EQ.5) THEN !set to defaults
        DO nd=1,NDATA
          C_DATA(nd)=C_DEFLT(nd)
        ENDDO
      ENDIF

      CALL EXITS('CINOUT')
      RETURN

 9999 CALL TRIM(ERROR,IBEG,IEND)
      IF(ERROR(1:IEND).EQ.'End of file') THEN
        IO_TYPE=1
        WRITE(OP_STRING,'('' End of file: revert to interactive '','
     '    //'''input'')')
      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GOTO 10
      ELSE IF(ERROR(1:4).EQ.'Edit') THEN
        ERROR=ERROR(1:IEND)//'>CINOUT'
        CLOSE(UNIT=IFIXED)
      ELSE IF(ERROR(1:7).EQ.'Restart') THEN
        ERROR=ERROR(1:IEND)//'>CINOUT'
        CLOSE(UNIT=IFIXED,STATUS='DELETE')
      ELSE IF(ERROR(1:14).EQ.'File not found') THEN
        ERROR='Help file not found>CINOUT'
      ELSE
        WRITE(UNIT=CIOSTA,FMT='(I5)') IOSTAT
        WRITE(UNIT=CIREC ,FMT='(I5)') irec
        ERROR=ERROR(1:IEND)//': IOSTAT='//CIOSTA(1:5)//' irec='//
     '    CIREC(1:5)//' C_FORMAT='//C_FORMAT(1:10)//'>CINOUT'
        CLOSE(UNIT=IFIXED)
      ENDIF
      CALL EXITS('CINOUT')
      RETURN 1
      END


      SUBROUTINE CINPUT(IUNIT,NCHAR,NDATA,CDATA,CDEFLT,INFO,ERROR,*)

C#### Subroutine: CINPUT
C###  Description:
C###    CINPUT reads character strings from IUNIT.
C**** IUNIT has a logical record length of IRECL.
C**** Strings separated by commas are read into CDATA(NDAT) NDAT=1,NDATA
C**** If an end of file or carriage return is detected while attempting
C**** to read CDATA(NDAT) it is assigned the value CDEFLT(NDAT).
C**** If a help request is made (?) a help file (based on INFO) is opened.
C**** If a stop request is made (STOP/stop/S/s) the program is stopped.
C**** If an edit request is made (EDIT/edit/E/e) control returns to
C**** the command line and the file is closed.
C**** If a restart request is made (RESTART/restart/R/r) control returns to
C**** the command line and the file is deleted.
C**** If an error is detected or a help request is made ERROR is
C**** returned with a diagnostic message.
C**** NOTE: The Q edit descriptor in the read format is nonstandard!

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER INFO,IUNIT,NCHAR,NDATA
      CHARACTER CDATA(*)*(*),CDEFLT(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER ic,ICOL,IEND,irec,IOSTAT,ISTART,nd,NDAT
      CHARACTER CHAR*5,CLINE*500,FILE*6

      CALL ENTERS('CINPUT',*9999)
      NDAT=0
      DO 4 irec=1,NDATA
        CALL READ_LINE(IUNIT,IOSTAT,NCHAR,CLINE,ERROR,*9999)
        IF(IOSTAT.EQ.0.AND.NCHAR.GT.0) THEN
          ICOL=1
          DO 3 nd=NDAT+1,NDATA
            ISTART=0
            DO 1 ic=ICOL,NCHAR
              IF(CLINE(ic:ic).EQ.',') THEN
                IF(ISTART.GT.0) THEN
                  IEND=ic-1
                  GOTO 2
                ENDIF
              ELSE
                IF(ISTART.EQ.0) THEN
                  ISTART=ic
                ENDIF
              ENDIF
    1       CONTINUE
            IEND=NCHAR
    2       IF(ISTART.GT.0) THEN
              ICOL=IEND+2
              NDAT=NDAT+1
              CDATA(NDAT)=CLINE(ISTART:IEND)
              IF((CDATA(NDAT).EQ.'EDIT').OR.
     '           (CDATA(NDAT).EQ.'edit').OR.
     '           (CDATA(NDAT).EQ.'E').OR.
     '           (CDATA(NDAT).EQ.'e'))THEN
                ERROR='Edit'
                GO TO 9999
              ELSE IF((CDATA(NDAT).EQ.'RESTART').OR.
     '                (CDATA(NDAT).EQ.'restart').OR.
     '                (CDATA(NDAT).EQ.'R').OR.
     '                (CDATA(NDAT).EQ.'r'))THEN
                ERROR='Restart'
                GO TO 9999
              ELSE IF((CDATA(NDAT).EQ.'STOP').OR.
     '                (CDATA(NDAT).EQ.'stop').OR.
     '                (CDATA(NDAT).EQ.'S').OR.
     '                (CDATA(NDAT).EQ.'s'))THEN
                CLOSE(UNIT=IUNIT)
                STOP
              ELSE IF(CDATA(NDAT).EQ.'?') THEN
                WRITE(UNIT=CHAR,FMT='(I5)') INFO
                FILE='FE'//CHAR(1:4)
                CALL DOCUM(FILE,'doc','Parameter number '//CHAR(5:5),
     '            ERROR,*9997)
 9997           IF(ERROR.EQ.' ') THEN
                  ERROR=' ?'
                  GO TO 9998
                ELSE IF(ERROR(1:14).EQ.'File not found') THEN
                  WRITE(OP_STRING,*) '>>Help not available'
      		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ERROR=' ?'
                  GO TO 9998
                ELSE
                  GO TO 9999
                ENDIF
              ELSE
                IF(NDAT.EQ.NDATA) THEN
                  ERROR=' '
                  GO TO 9998
                ENDIF
              ENDIF
            ELSE
              GOTO 4
            ENDIF
    3     CONTINUE
        ELSE IF(IOSTAT.EQ.0.AND.NCHAR.EQ.0) THEN
C         REWIND(UNIT=IUNIT)
          NDAT=NDAT+1
C         CDATA(NDAT)=CDEFLT(NDAT)
C         IF(NDAT.EQ.NDATA) THEN
          DO nd=1,NDATA
            CDATA(nd)=CDEFLT(nd)
          ENDDO
          ERROR=' '
          GO TO 9998
        ELSE
          REWIND(UNIT=IUNIT)
          ERROR=' Character data read error'
          GO TO 9998
        ENDIF
    4 CONTINUE

      ERROR=' '
 9998 CALL EXITS('CINPUT')
      RETURN
 9999 CALL ERRORS('CINPUT',ERROR)
      CALL EXITS('CINPUT')
      RETURN 1
      END


      SUBROUTINE IINOUT(IO_TYPE,IFREE,IFIXED,I_FORMAT,NDATA,I_DATA,
     '  I_DEFLT,IMIN,IMAX,INFO,ERROR,*)

C#### Subroutine: IINOUT
C###  Description:
C###    IINOUT handles input/output of integer values.
C**** If IO_TYPE=0 prompts written to IFREE
C****  "    "    1   "     written to IFREE and (with data) to IFIXED
C****  "    "    2   "     read with data from IFIXED
C****  "    "    3   "     written with data to IFIXED
C****  "    "    4   "     read with data from IFIXED & written to IFREE
C****  "    "    5   default values are used
C**** If IFIXED is not open no IO is performed on that file.
C**** Integer values separated by commas or blanks are read into
C**** I_DATA(NDAT) NDAT=1,NDATA.
C**** The limiting values for the read are defined by IMIN and IMAX.
C**** If an end of file or carriage return is detected while attempting
C**** to read I_DATA(NDAT) it is assigned the value I_DEFLT(NDAT).
C**** If an error is detected a diagnostic message is returned in ERROR
C**** and control is returned to the statement number of the asterisk.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
!     Parameter List
      INTEGER IFIXED,IFREE,IMAX,IMIN,INFO,IO_TYPE,I_DATA(*),
     '  I_DEFLT(*),NDATA
      CHARACTER I_FORMAT*(*),ERROR*(*)
!     Local Variables
      INTEGER CLOCAT,IBEG,IEND,IOSTAT,IPOS,irec,IUNIT,JPOS,m,MAXTRY,n,nd
      CHARACTER CHAR*1,CIOSTA*5,CIREC*5,CIUNIT*10,LINE*132,SUBSTR*1
      LOGICAL OPENED,REVERT
      DATA MAXTRY/8/

      CALL ENTERS('IINOUT',*9999)
      IF(NDATA.GT.IOIM) THEN
        ERROR=' NDATA exceeds the length of the array I_DATA'
        GOTO 9999
      ENDIF
      REVERT=.FALSE.
 10   IF(IO_TYPE.LE.1.OR.REVERT) THEN
        DO 1 m=1,MAXTRY
          CALL WRITE_LINE(IOOP,I_FORMAT,ERROR,*9999)
          CALL IINPUT(IFREE,NDATA,I_DATA,I_DEFLT,IMIN,IMAX,
     '      INFO,ERROR,*9999)
          IF(ERROR.EQ.' ') GOTO 2
          IF(ERROR(1:2).NE.' ?') THEN
            WRITE(OP_STRING,*) ERROR
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
    1   CONTINUE
        ERROR=' More than 8 invalid input attempts made'
        GOTO 9999
    2   IF(IO_TYPE.EQ.1) THEN
          INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
          IF(OPENED.AND.IFIXED.NE.0) THEN
            WRITE(UNIT=IFIXED,FMT=I_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '        (I_DATA(n),n=1,NDATA)
            IF(IOSTAT.EQ.0) THEN
              ERROR=' '
            ELSE
              ERROR=' Integer data write error'
              GOTO 9999
            ENDIF
          ELSE
            ERROR=' '
          ENDIF
        ENDIF

      ELSE IF(IO_TYPE.EQ.2.OR.IO_TYPE.EQ.4) THEN
        IF(DOP) THEN
          WRITE(OP_STRING,'(A)') I_FORMAT(1:80)
      	  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(OPENED) THEN
 3        READ(UNIT=IFIXED,FMT='(A)',REC=irec,IOSTAT=IOSTAT) LINE
          IF(LINE(1:1).EQ.'!') THEN   !This line is a comment
            IF(LINE(2:2).EQ.'!') THEN !Do not print
            ELSE                      !Print comment line
              WRITE(OP_STRING,'(A)') LINE
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            irec=irec+1
            GO TO 3
          ENDIF
          SUBSTR=''''
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' SUBSTR=',SUBSTR
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IPOS=CLOCAT(SUBSTR,I_FORMAT)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' IPOS=',IPOS
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CHAR=' '
          IF(IPOS.GT.0) THEN
            SUBSTR='/'
            JPOS=CLOCAT(SUBSTR,I_FORMAT(1:4))
            IF(JPOS.GT.0) THEN
              READ(UNIT=IFIXED,FMT='(/A,'''//I_FORMAT(IPOS+2:),
     '          REC=irec,IOSTAT=IOSTAT) CHAR,(I_DATA(n),n=1,NDATA)
            ELSE
              READ(UNIT=IFIXED,FMT='(A,'''//I_FORMAT(IPOS+2:),
     '          REC=irec,IOSTAT=IOSTAT) CHAR,(I_DATA(n),n=1,NDATA)
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' CHAR='',A)') CHAR
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE
            READ(UNIT=IFIXED,FMT=I_FORMAT,
     '        REC=irec,IOSTAT=IOSTAT) (I_DATA(n),n=1,NDATA)
          ENDIF
          IF(IOSTAT.EQ.0) THEN
            ERROR=' '
          ELSE IF(IOSTAT.EQ.36) THEN
            ERROR='End of file'
            GO TO 9999
          ELSE
            ERROR=' Integer data read error'
            GOTO 9999
          ENDIF
          IF(CHAR.EQ.'*') THEN
            REVERT=.TRUE.
            GO TO 10
          ENDIF
        ELSE
          WRITE(UNIT=CIUNIT,FMT='(I10)') IUNIT
          CALL TRIM(CIUNIT,IBEG,IEND)
          ERROR=' Unit '//CIUNIT(IBEG:IEND)//' is not open'
          GOTO 9999
        ENDIF
        IF(IO_TYPE.EQ.4) THEN
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' IO_TYPE=',IO_TYPE,' IFREE=',IFREE
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*) I_FORMAT
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*) ' NDATA=',NDATA
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*) (I_DATA(n),n=1,NDATA)
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(OP_STRING,I_FORMAT) (I_DATA(n),n=1,NDATA)
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(IO_TYPE.EQ.3) THEN !write data to file
        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(OPENED.AND.IFIXED.NE.0) THEN
          WRITE(UNIT=IFIXED,FMT=I_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '      (I_DATA(N),N=1,NDATA)
          IF(IOSTAT.EQ.0) THEN
            ERROR=' '
          ELSE
            ERROR=' Integer data write error'
            GOTO 9999
          ENDIF
        ELSE
          ERROR=' '
        ENDIF

      ELSE IF(IO_TYPE.EQ.5) THEN !set to defaults
        DO nd=1,NDATA
          I_DATA(nd)=I_DEFLT(nd)
        ENDDO
      ENDIF

      CALL EXITS('IINOUT')
      RETURN

 9999 CALL TRIM(ERROR,IBEG,IEND)
      IF(ERROR(1:IEND).EQ.'End of file') THEN
        IO_TYPE=1
        WRITE(OP_STRING,'('' End of file: revert to interactive '','
     '    //'''input'')')
      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GOTO 10
      ELSE IF(ERROR(1:4).EQ.'Edit') THEN
        ERROR=ERROR(1:IEND)//'>IINOUT'
        CLOSE(UNIT=IFIXED)
      ELSE IF(ERROR(1:7).EQ.'Restart') THEN
        ERROR=ERROR(1:IEND)//'>IINOUT'
        CLOSE(UNIT=IFIXED,STATUS='DELETE')
      ELSE IF(ERROR(1:14).EQ.'File not found') THEN
        ERROR='Help file not found>IINOUT'
      ELSE
        WRITE(UNIT=CIOSTA,FMT='(I5)') IOSTAT
        WRITE(UNIT=CIREC ,FMT='(I5)') irec
        ERROR=ERROR(1:IEND)//': IOSTAT='//CIOSTA(1:5)//' irec='//
     '    CIREC(1:5)//' I_FORMAT='//I_FORMAT(1:10)//'>IINOUT'
        CLOSE(UNIT=IFIXED)
      ENDIF
      CALL EXITS('IINOUT')
      RETURN 1
      END


      SUBROUTINE IINPUT(IUNIT,NDATA,IDATA,IDEFLT,IMIN,IMAX,
     '  INFO,ERROR,*)

C#### Subroutine: IINPUT
C###  Description:
C###    IINPUT reads integer values from IUNIT.
C**** Integers separated by blanks or commas are read into
C**** IDATA(NDAT) NDAT=1,NDATA.
C**** The limiting values for the read are defined by IMIN and IMAX.
C**** If an end of file or carriage return is detected while attempting
C**** to read IDATA(NDAT) it is assigned the value IDEFLT(NDAT).
C**** No range check is made on the default integers.
C**** If help request is made (?) a help file (based on INFO) is opened.
C**** If stop request is made (STOP/stop/S/s) the program is stopped.
C**** If an edit request is made (EDIT/edit/E/e) control returns to
C**** the command line and the file is closed.
C**** If restart request is made (RESTART/restart/R/r) control returns 
C**** to the command line and the file is deleted.
C**** If an error is detected or a help request is made ERROR is
C**** returned with a diagnostic message.
C**** NOTE: The Q edit descriptor in the read format is nonstandard!

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER IDATA(*),IDEFLT(*),IMAX,IMIN,INFO,IUNIT,NDATA
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ic,ICOL,IEND,IOSTAT,irec,ISTART,LENGTH,NCHAR,nd,NDAT
      CHARACTER CDATA*16,CHAR*5,CLINE*500,FILE*6
      LOGICAL IVALID

      CALL ENTERS('IINPUT',*9999)
      NDAT=0
      DO 4 irec=1,NDATA
        CALL READ_LINE(IUNIT,IOSTAT,NCHAR,CLINE,ERROR,*9999)
        IF(IOSTAT.EQ.0.AND.NCHAR.GT.0) THEN
          ICOL=1
          DO 3 nd=NDAT+1,NDATA
            ISTART=0
            DO 1 ic=ICOL,NCHAR
              IF((CLINE(ic:ic).EQ.' ').OR.(CLINE(ic:ic).EQ.',')) THEN
                IF(ISTART.GT.0) THEN
                  IEND=ic-1
                  GOTO 2
                ENDIF
              ELSE
                IF(ISTART.EQ.0) THEN
                  ISTART=ic
                ENDIF
              ENDIF
    1       CONTINUE
            IEND=NCHAR
    2       IF(ISTART.GT.0) THEN
              ICOL=IEND+2
              CDATA=' '
              LENGTH=IEND-ISTART+1
              CDATA(17-LENGTH:16)=CLINE(ISTART:IEND)
              IF((CDATA.EQ.'            EDIT').OR.
     '           (CDATA.EQ.'            edit').OR.
     '           (CDATA.EQ.'               E').OR.
     '           (CDATA.EQ.'               e'))THEN
                ERROR='Edit'
                GO TO 9999
              ELSE IF((CDATA.EQ.'         RESTART').OR.
     '                (CDATA.EQ.'         restart').OR.
     '                (CDATA.EQ.'               R').OR.
     '                (CDATA.EQ.'               r'))THEN
                ERROR='Restart'
                GO TO 9999
              ELSE IF((CDATA.EQ.'            STOP').OR.
     '                (CDATA.EQ.'            stop').OR.
     '                (CDATA.EQ.'               S').OR.
     '                (CDATA.EQ.'               s'))THEN
                CLOSE(UNIT=IUNIT)
                STOP
              ELSE IF(CDATA.EQ.'               ?') THEN
                WRITE(UNIT=CHAR,FMT='(I5)') INFO
                FILE='FE'//CHAR(1:4)
                CALL DOCUM(FILE,'doc','Parameter number '//CHAR(5:5),
     '            ERROR,*9997)
 9997           IF(ERROR.EQ.' ') THEN
                  ERROR=' ?'
                  GO TO 9998
                ELSE IF(ERROR(1:14).EQ.'File not found') THEN
                  WRITE(OP_STRING,*) '>>Help not available'
      		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ERROR=' ?'
                  GO TO 9998
                ELSE
                  GO TO 9999
                ENDIF
              ELSE
                IF(IVALID(CDATA)) THEN
                  NDAT=NDAT+1
                  READ(UNIT=CDATA,FMT='(I16)',IOSTAT=IOSTAT) IDATA(NDAT)
                  IF(IOSTAT.EQ.0) THEN
                    IF((IDATA(NDAT).GE.IMIN).AND.
     '                 (IDATA(NDAT).LE.IMAX)) THEN
                      IF(NDAT.EQ.NDATA) THEN
                        ERROR=' '
                        GO TO 9998
                      ENDIF
                    ELSE
                      ERROR=' Integer data range error'
                      GO TO 9998
                    ENDIF
                  ELSE
                    ERROR=' Integer data read error'
                    GO TO 9998
                  ENDIF
                ELSE
                  ERROR=' Integer data type error'
                  GO TO 9998
                ENDIF
              ENDIF
            ELSE
              GOTO 4
            ENDIF
    3     CONTINUE
        ELSE IF(IOSTAT.EQ.0.AND.NCHAR.EQ.0) THEN
          NDAT=NDAT+1
          DO nd=1,NDATA
            IDATA(nd)=IDEFLT(nd)
          ENDDO
          ERROR=' '
          GO TO 9998
        ELSE
          REWIND(UNIT=IUNIT)
          ERROR=' Character data read error'
          GO TO 9998
        ENDIF
    4 CONTINUE

      ERROR=' '
 9998 CALL EXITS('IINPUT')
      RETURN
 9999 CALL ERRORS('IINPUT',ERROR)
      CALL EXITS('IINPUT')
      RETURN 1
      END


      SUBROUTINE INOUT(IO_TYPE,IFREE,IFIXED,FORMAT,ERROR,*)

C#### Subroutine: INOUT
C###  Description:
C###    INOUT handles prompting input/output.
C**** If IO_TYPE=0 prompts written to IFREE
C****  "    "    1   "     written to IFREE and (with data) to IFIXED
C****  "    "    2   "     read with data from IFIXED
C****  "    "    3   "     written with data to IFIXED
C****  "    "    4   "     read with data from IFIXED & written to IFREE
C****  "    "    5   nothing is done
C**** If IFIXED is not open no IO is performed on that file.
C**** If an error is detected a diagnostic message is returned in ERROR
C**** and control is returned to the statement number of the first
C**** asterisk.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER IFIXED,IFREE,IO_TYPE
      CHARACTER FORMAT*(*),ERROR*(*)
!     Local Variables
      INTEGER CLOCAT,IBEG,IEND,IOSTAT,IPOS,irec,IUNIT,JPOS
      CHARACTER CHAR*1,CIOSTA*5,CIREC*5,CIUNIT*10,LINE*132,SUBSTR*1
      LOGICAL OPENED,REVERT

      CALL ENTERS('INOUT',*9999)
      REVERT=.FALSE.
 10   IF(IO_TYPE.LE.1.OR.REVERT) THEN
      	CALL WRITE_LINE(IOOP,FORMAT,ERROR,*9999)
        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(IO_TYPE.EQ.1) THEN
          IF(OPENED) THEN
            WRITE(UNIT=IFIXED,FMT=FORMAT,REC=irec,IOSTAT=IOSTAT)
            IF(IOSTAT.EQ.0) THEN
              ERROR=' '
            ELSE
              ERROR=' Write error'
              GOTO 9999
            ENDIF
          ELSE
            ERROR=' '
          ENDIF
        ENDIF

      ELSE IF(IO_TYPE.EQ.2.OR.IO_TYPE.EQ.4) THEN
        IF(DOP) THEN
          WRITE(OP_STRING,'(A)') FORMAT(1:80)
      	  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(OPENED) THEN
 3        READ(UNIT=IFIXED,FMT='(A)',REC=irec,IOSTAT=IOSTAT) LINE
          IF(LINE(1:1).EQ.'!') THEN   !This line is a comment
            IF(LINE(2:2).EQ.'!') THEN !Do not print
            ELSE                      !Print comment line
              WRITE(OP_STRING,'(A)') LINE
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            irec=irec+1
            GO TO 3
          ENDIF
          SUBSTR=''''
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' SUBSTR=',SUBSTR
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IPOS=CLOCAT(SUBSTR,FORMAT)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' IPOS=',IPOS
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CHAR=' '
          IF(IPOS.GT.0) THEN
            SUBSTR='/'
            JPOS=CLOCAT(SUBSTR,FORMAT(1:4))
            IF(JPOS.GT.0) THEN
              READ(UNIT=IFIXED,FMT='(/A,'''//FORMAT(IPOS+2:),
     '          REC=irec,IOSTAT=IOSTAT) CHAR
            ELSE
              READ(UNIT=IFIXED,FMT='(A,'''//FORMAT(IPOS+2:),
     '          REC=irec,IOSTAT=IOSTAT) CHAR
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' CHAR='',A)') CHAR
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE
            READ(UNIT=IFIXED,FMT=FORMAT,
     '        REC=irec,IOSTAT=IOSTAT)
          ENDIF
          IF(IOSTAT.EQ.0) THEN
            ERROR=' '
          ELSE IF(IOSTAT.EQ.36) THEN
            ERROR=' End of file'
            GOTO 9999
          ELSE
            ERROR=' Read error'
            GOTO 9999
          ENDIF
          IF(CHAR.EQ.'*') THEN
            REVERT=.TRUE.
            GO TO 10
          ENDIF
        ELSE
          WRITE(UNIT=CIUNIT,FMT='(I10)') IFIXED
          CALL TRIM(CIUNIT,IBEG,IEND)
          ERROR=' Unit '//CIUNIT(IBEG:IEND)//' is not open'
          GOTO 9999
        ENDIF
        IF(IO_TYPE.EQ.4) THEN
          WRITE(OP_STRING,FORMAT)
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(IO_TYPE.EQ.3) THEN !write data to file
        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(OPENED) THEN
          WRITE(UNIT=IFIXED,FMT=FORMAT,REC=irec,IOSTAT=IOSTAT)
          IF(IOSTAT.EQ.0) THEN
            ERROR=' '
          ELSE
            ERROR=' Write error'
            GOTO 9999
          ENDIF
        ELSE
          ERROR=' '
        ENDIF

      ELSE IF(IO_TYPE.EQ.5) THEN !set to defaults (do nothing here)
      ENDIF

      CALL EXITS('INOUT')
      RETURN

 9999 CALL TRIM(ERROR,IBEG,IEND)
      IF(ERROR(1:IEND).EQ.'End of file') THEN
        IO_TYPE=1
        WRITE(OP_STRING,'('' End of file: revert to interactive '','
     '    //'''input'')')
      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GOTO 10
      ENDIF    
      WRITE(UNIT=CIOSTA,FMT='(I5)') IOSTAT
      WRITE(UNIT=CIREC ,FMT='(I5)') irec
      ERROR=ERROR(1:IEND)//': IOSTAT='//CIOSTA(1:5)//' irec='//
     '  CIREC(1:5)//' FORMAT='//FORMAT(1:10)//'>INOUT'
      CLOSE(UNIT=IFIXED)
      CALL EXITS('INOUT')
      RETURN 1
      END

      SUBROUTINE IOPOTN(COMAND,DATAFORM,IUNIT,DATASET,TITLE,
     '  PD,WD,ERROR,*)

C#### Subroutine: IOPOTN
C###  Description:
C###    IOPOTN handles I/O of measured potential data for processing.
C**** COMAND determines whether data is to be read 'READ' from or
C**** written 'WRITE' to IUNIT.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER DATASET,IUNIT
      REAL*8 PD(*),WD(NJM,NDM)
      CHARACTER COMAND*(*),DATAFORM*(*),ERROR*(*),TITLE*80
!     Local Variables
      INTEGER IBEG,IEND,IOSTAT,nd,NDD,NDPT,nr
      CHARACTER ACCESS*10,CFROMI*5,CI*5,FMTI*8,FORM*11,IOS*5
      LOGICAL OPENED
      CALL ENTERS('IOPOTN',*9999)
      nr=1 !temporary
      IF(DOP) THEN
        WRITE(OP_STRING,'('' COMMAND='',A,'', DATAFORM='',A,'
     '    //''', IUNIT='',I3,'', DATASET='',I4)') COMAND,DATAFORM,
     '    IUNIT,DATASET
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(DATAFORM(1:4).EQ.'EMAP') THEN
      ELSE IF(DATAFORM(1:6).EQ.'IPFILE') THEN
        INQUIRE(UNIT=IUNIT,OPENED=OPENED)
        IF(OPENED) THEN
          INQUIRE(UNIT=IUNIT,FORM=FORM)
          IF(FORM.EQ.'FORMATTED') THEN
            INQUIRE(UNIT=IUNIT,ACCESS=ACCESS)
            IF(ACCESS.EQ.'SEQUENTIAL') THEN
              IF(COMAND.EQ.'READ') THEN          
               READ(UNIT=IUNIT,FMT='(A80)') TITLE(1:80)
                NDPT=0
                DO nd=1,NDT
                  READ(UNIT=IUNIT,FMT=*,IOSTAT=IOSTAT,END=9998)
     '              ndd,PD(nd),WD(NJT+1,nd)
                  IF(IOSTAT.GT.0) THEN
                    WRITE(UNIT=IOS,FMT='(I5)') IOSTAT
                    ERROR=' Read error: IOSTAT='//IOS
                    GOTO 9999
                  ENDIF
                  NDPT=NDPT+1
                ENDDO
                IF(NDPT.NE.NDT) THEN
                  WRITE(OP_STRING,'('' >>WARNING: number of '','
     '              //'''potential values read does not equal NDT'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
              ELSE IF(COMAND.EQ.'WRITE') THEN
                WRITE(UNIT=IUNIT,FMT='('' Potential file:'')')
                DO nd=1,NDT
                  WRITE(UNIT=IUNIT,FMT=*,IOSTAT=IOSTAT) nd,PD(nd),
     '              WD(NJT+1,nd)
                  IF(IOSTAT.NE.0) THEN
                    WRITE(UNIT=IOS,FMT='(I5)') IOSTAT
                    ERROR=' Write error: IOSTAT='//IOS
                    GOTO 9999
                  ENDIF
                ENDDO
              ELSE
                ERROR=' Command error: COMAND='//COMAND
                GOTO 9999
              ENDIF
            ELSE
              CI=CFROMI(IUNIT,FMTI)
              CALL TRIM(CI,IBEG,IEND)
              ERROR=' File '//CI(IBEG:IEND)//
     '          ' is not connected for sequential access'
                GOTO 9999
            ENDIF
          ELSE            
                CI=CFROMI(IUNIT,FMTI)
            CALL TRIM(CI,IBEG,IEND)
            ERROR=' File '//
     '        CI(IBEG:IEND)//' is not connected for formatted i/o'
            GOTO 9999
          ENDIF
        ELSE
          CI=CFROMI(IUNIT,FMTI)
          CALL TRIM(CI,IBEG,IEND)
          ERROR=' File '//CI(IBEG:IEND)//' is not opened'
          GOTO 9999
        ENDIF
      ELSE IF(DATAFORM(1:5).EQ.'MAP3D') THEN
      ELSE
        ERROR='Data format error : DATAFORM='//DATAFORM
        GOTO 9999
      ENDIF 
 9998 ERROR=' '
      CALL EXITS('IOPOTN')
      RETURN
 9999 CALL ERRORS('IOPOTN',ERROR)
      CALL EXITS('IOPOTN')
      RETURN 1
      END

                  
      SUBROUTINE IPGEOM(COMAND,IUNIT,IBT,IDO,INP,ITHRES,NAN,NBH,NBJ,
     '  NCO,NCNP,NEELEM,NFF,NGAP,NHE,NHP,NJE,NJP,NKE,NKH,NKJ,
     '  NLL,NNF,NNL,NPE,NPF,NPL,NPNODE,NQE,NW,NYNE,NYNP,
     '  CE,CP,DL,FEXT,PE,PF,PG,SE,THRES,VE,WG,XA,XIG,XP,YP,
     '  FIX,FIXP,ERROR,*)

C**** Reads previous version of IOD files.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b10.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:b13.cmn'
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:b22.cmn'
      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:four00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'            
      INCLUDE 'cmiss$reference:ktyp30.cmn'
      INCLUDE 'cmiss$reference:ktyp40.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:ktyp60.cmn'
      INCLUDE 'cmiss$reference:ktyp70.cmn'
!     Parameter List
      INTEGER IBT(2,NIM,*),IDO(NKM,0:NIM,*),INP(NNM,NIM,*),
     '  ITHRES(3,NGM,*),NAN(NIM,NAM,*),NBH(NCM,NHM,*),NBJ(NCM,NJM,*),
     '  NCNP(0:NCM,NJM,NPM,*),NCO(*),NEELEM(0:NEM,0:*),NFF(6,*),
     '  NGAP(NIM,*),NHE(*),NHP(*),NJE(*),NJP(*),
     '  NKE(NKM,NNM,NBM,*),NKH(NCM,NHM,*),NKJ(NCM,NJM,*),
     '  NLL(12,*),NNF(0:14,6,*),NNL(4,12,*),
     '  NPE(NNM,NBM,*),NPF(12,*),NPL(20,*),
     '  NPNODE(0:NPM,0:*),NQE(NSM,NBM,*),NW(NEM,*),
     '  NYNE(NAM,NCM,NHM,*),NYNP(NKM,NCM,NHM,0:NCRM,NPM,0:*)
      REAL*8 CE(NMM,*),CP(NMM,*),DL(3,*),FEXT(8,NGM,*),
     '  PE(2,*),PF(2,*),PG(NSM,NUM,NGM,*),SE(NSM,NBM,*),THRES(3,NGM,*),
     '  VE(NSM,NKM,*),WG(NGM,*),XA(NAM,NJM,*),XIG(NIM,NGM,*),
     '  XP(NKM,NCM,NJM,*),YP(NYM,NIYM,0:*)
      CHARACTER COMAND*(*),ERROR*(*)
      LOGICAL FIX(NYM,5,0:*),FIXP(2,*)
!     Local Variables
      INTEGER I,ID,ID1,ID2,IE,IEE,IE_TOT,IFACE,IL,IUNIT,IY,
     '  J,MG,NA,NB,NC,NCC,ncr,NE,NEE,NF,NG,NH,NHH,
     '  NI,NJ,NK,NL,NM,NN,NOELEM,NONODE,NP,NPP,NQ,NR,NRR,NS,NV,nx,
     '  NY,NYY
      CHARACTER FMT*500
      LOGICAL DEPENDENT

c      CALL ENTERS('IPGEOM',*9999)
cC AJP 3-9-92      NC=1 !Temporary AJP 17-12-91
c      ncr=1 ! Temporary cpb 13/11/94
c      nx=1 ! temporary
c      IF(COMAND.EQ.'READ') THEN
c        FMT='('' NBT='',I6,'' NDT='',I6,'' NFT='',I6,'
c     '     //''' NJT='',I6,'' NLT='',I6,'' NOT='',I6,/'
c     '     //''' NQT='',I6,'' NRT='',I6,'' NZT='',I6,'
c     '     //''' NBFT='',I5)'
c        READ(IUNIT,FMT) NBT,NDT,NFT,NJT,NLT,NOT(1,1),NQT,NRT,NZT,NBFT
c
c        FMT='(/'' JTYP1 ='',I2,'' JTYP2 ='',I2,'' JTYP3 ='',I2,'
c     '      //''' JTYP4 ='',I2,'' JTYP5 ='',I2,'' JTYP6 ='',I2,/'
c     '      //''' JTYP7 ='',I2,'' JTYP8 ='',I2,'' JTYP9 ='',I2,'
c     '      //''' JTYP10='',I2,'' JTYP11='',I2,'' JTYP12='',I2,/'
c     '      //''' JTYP13='',I2,'' JTYP14='',I2,'' JTYP15='',I2)'
c        READ(IUNIT,FMT) JTYP1,JTYP2,JTYP3,JTYP4,JTYP5,JTYP6,
c     '                  JTYP7,JTYP8,JTYP9,JTYP10,JTYP11,JTYP12,
c     '                  JTYP13,JTYP14,JTYP15
c
c        FMT='(/'' NET(0)='',I6,'' NPT(0)='',I6,'' NYT(1,0)='',I6)'
c        READ(IUNIT,FMT) NET(0),NPT(0),NYT(1,0,1)
c
c        DO NR=1,NRT
c          FMT='(/'' Region '',I1,'': ITYP1 ='',I2,'' ITYP2 ='',I2,'
c     '                         //''' ITYP3 ='',I2,'' ITYP4 ='',I2,/'
c     '               //'''           ITYP5 ='',I2,'' ITYP6 ='',I2,'
c     '                         //''' ITYP7 ='',I2,'' ITYP8 ='',I2,)'
c          READ(IUNIT,FMT)   NRR,ITYP1(nr),ITYP2(nr),ITYP3(nr),
c     '      ITYP4(nr),ITYP5(nr),ITYP6(nr),ITYP7(nr),ITYP8(nr)
c          FMT='('' NET('',I2,'')='',I6,'' NPT('',I2,'')='',I6,'
c     '      //' '' NYT(1,'',I2,'')='',I6)'
c          READ(IUNIT,FMT) NRR,NET(nr),NRR,NPT(nr),NRR,NYT(1,NR,1)
c          FMT='('' NCT('',I2,'')='',I2)'
c          READ(IUNIT,FMT) NRR,NCT(nr) 
c        ENDDO
c
c        DO NB=1,NBT
c          FMT='(/'' Basis type nb='',I1)'
c          READ(IUNIT,FMT) ID
c          FMT='('' NAT(nb)='',I3,''  NBI(nb)='',I3,''  NFE(nb)='',I3,'
c     '      //'''  NGT(nb)='',I3,''  NIT(nb)='',I3,''  NKT(nb)='',I3,/'
c     '       //''' NLE(nb)='',I3,''  NNT(nb)='',I3,''  NST(nb)='',I3,'
c     '      //'''  NUT(nb)='',I3,''  NBC(nb)='',I3,''  NBF(nb)='',I3,/'
c     '       //''' NABTYP(nb)='',I3)'
c          READ(IUNIT,FMT) NAT(NB),NBI(NB),NFE(NB),NGT(NB),NIT(NB),
c     '      NKT(0,NB),NLE(NB),NNT(NB),NST(NB),NUT(NB),NBC(NB),NBF(NB),
c     '      NABTYP(NB)
c          FMT='('' IBT( i=1..,ni=1..,nb):'',40I2)'
c          READ(IUNIT,FMT) ((IBT(I,NI,NB),I=1,2),NI=1,NIT(NB))
c          FMT='('' IDO(nk=1..,ni=1..,nb):'',40I2)'
c          READ(IUNIT,FMT) ((IDO(NK,NI,NB),NK=1,NKT(0,NB)),NI=0,NIT(NB))
c          FMT='('' INP(nn=1..,ni=1..,nb):'',40I2)'
c          READ(IUNIT,FMT) ((INP(NN,NI,NB),NN=1,NNT(NB)),NI=1,NIT(NB))
c          FMT='('' NAN(ni=1..,na=1..,nb):'',40I2)'
c          READ(IUNIT,FMT) ((NAN(NI,NA,NB),NI=1,NIT(NB)),NA=1,NAT(NB))
c          FMT='('' NGAP(ni=1..,nb):'',20I4)'
c          READ(IUNIT,FMT) (NGAP(NI,NB),NI=1,NIT(NB))
c          FMT='('' NNL(i=1..,j=1..,nb):'',/(20I4))'
c          READ(IUNIT,FMT) ((NNL(I,J,NB),I=1,4),J=1,12)
c          FMT='('' NNF(i=1..,j=1..,nb):'',/(20I4))'
c          READ(IUNIT,FMT) ((NNF(I,J,NB),I=0,14),J=1,6)
c          IF(NBC(NB).EQ.4) THEN
c            FMT='('' OMEGA='',E25.17)'
c            READ(IUNIT,FMT) OMEGA
c          ENDIF
c        ENDDO
c
c        FMT='(/'' NPNODE(0,0):'',I5,''  NEELEM(0,0):'',I5)'
c        READ(IUNIT,FMT) NPNODE(0,0),NEELEM(0,0)
c        DO NR=1,NRT
c          FMT='(/'' Region nr='',I2)'
c          READ(IUNIT,FMT) NRR
c          FMT='(/'' NPNODE(0,nr):'',I5)'
c          READ(IUNIT,FMT) NPNODE(0,NR)
c          FMT='('' NPNODE(nonode=1..,'',I1,''):''/,(20I5))'
c          READ(IUNIT,FMT) NRR,(NPNODE(NONODE,NR),NONODE=1,NPNODE(0,NR))
c          FMT='(/'' NEELEM(0,nr):'',I5)'
c          READ(IUNIT,FMT) NEELEM(0,NR)
c          FMT='('' NEELEM(noelem=1..,'',I1,''):''/,(20I5))'
c          READ(IUNIT,FMT) NRR,(NEELEM(NOELEM,NR),NOELEM=1,NEELEM(0,NR))
c        ENDDO
c
c        FMT='(/'' Focus:'',E25.17)'
c        READ(IUNIT,FMT) FOCUS
c
c        DO NR=1,NRT
c          DO NONODE=1,NPNODE(0,NR)
c            NP=NPNODE(NONODE,NR)
c            FMT='(/'' Node np='',I4)'
c            READ(IUNIT,FMT) ID
c            FMT='('' NJP(np)='',I2)'
c            READ(IUNIT,FMT) NJP(NP)
c            FMT='('' NKJ(1,nj=1..,np):'',10I3)'
c            READ(IUNIT,FMT) (NKJ(1,NJ,NP),NJ=1,NJP(NP)+JTYP9+JTYP11)
c            DO NJ=1,NJP(NP)+JTYP9+JTYP11
c              FMT='('' XP(nk=1..,1,nj='',I1,'',np):'',3E25.17,'
c     '          //':/3E25.17)'
c              READ(IUNIT,FMT) ID,(XP(NK,1,NJ,NP),NK=1,NKJ(1,NJ,NP))
c            ENDDO
c          ENDDO
c        ENDDO
c
cC       NEFMAX=MAX0(NET(1),NFT)
c        DO NR=1,NRT
c          DO NOELEM=1,NEELEM(0,NR)
c            NE=NEELEM(NOELEM,NR)
c            FMT='(/'' Element ne='',I4)'
c            READ(IUNIT,FMT) ID
c            FMT='('' NJE(ne)='',I2)'
c            READ(IUNIT,FMT) NJE(NE)
c            FMT='('' NBJ(1,nj=1..,ne):'',10I2)'
c            READ(IUNIT,FMT) (NBJ(1,NJ,NE),NJ=1,NJE(NE)+JTYP9+JTYP11)
c            FMT='('' NCO(ne)='',I2)'
c            READ(IUNIT,FMT) NCO(NE)
c            FMT='('' NLL(i=1..,ne):'',20I4)'
c            READ(IUNIT,FMT) (NLL(I,NE),I=1,12)
c            FMT='('' NFF(i=1..,ne):'',20I4)'
c            READ(IUNIT,FMT) (NFF(I,NE),I=1,6)
c            DO NB=1,NBT
c              FMT='('' NPE(nn=1..,nb='',I1,'',ne):'',20I4)'
c              READ(IUNIT,FMT) ID,(NPE(NN,NB,NE),NN=1,NNT(NB))
c              DO NN=1,NNT(NB)
c                FMT='('' NKE(nk=1..,nn='',I2,'',nb='',I1,'',ne):'','
c     '            //'(20I4))'
c                READ(IUNIT,FMT) ID1,ID2,
c     '             (NKE(NK,NN,NB,NE),NK=1,NKT(0,NB))
c              ENDDO
c              FMT='('' NQE(ns=1..,nb='',I1,'',ne):'',/(20I4))'
c              READ(IUNIT,FMT) ID,(NQE(NS,NB,NE),NS=1,NST(NB))
c              FMT='('' SE(ns=1..,nb='',I1,'',ne):'',/(3E25.17))'
c              READ(IUNIT,FMT) ID,(SE(NS,NB,NE),NS=1,NST(NB))
c              IF(IBT(1,1,NB).EQ.5) THEN
c                FMT='('' VE(ns=1..,nk=1..,ne):'',/(3E25.17))'
c                READ(IUNIT,FMT) ((VE(NS,NK,NE),NS=1,NST(NB)),NK=1,
c     '            NKT(0,NB))
c              ENDIF
c            ENDDO
c          ENDDO
c        ENDDO
c
c        FMT='(/'' XA(1,nj=1..,nq=1..):'',/(3E25.17))'
c        READ(IUNIT,FMT) ((XA(1,NJ,NQ),NJ=1,NJM),NQ=1,NQT)
c
c        DO NF=1,NFT
c          FMT='(/'' Face nf='',I4)'
c          READ(IUNIT,FMT) ID
c          FMT='('' NPF(i=1..,nf):'',10I4)'
c          READ(IUNIT,FMT) (NPF(I,NF),I=1,10)
c        ENDDO
c
c        DO NL=1,NLT
c          FMT='(/'' Line nl='',I4)'
c          READ(IUNIT,FMT) ID
c          FMT='('' NPL(i=1..,nl):'',18I4)'
c          READ(IUNIT,FMT) (NPL(I,NL),I=1,18)
c          FMT='('' DL(i=1..,nl):'',3E25.17)'
c          READ(IUNIT,FMT) (DL(I,NL),I=1,3)
c        ENDDO
c
c        FMT='(/'' CALL_BASE='',L1,'' CALL_ELEM='',L1,'
c     '    //'  '' CALL_LINE='',L1,'' CALL_NODE='',L1)'
c        READ(IUNIT,FMT) CALL_BASE,CALL_ELEM,CALL_LINE,CALL_NODE
c
c        DEPENDENT=.FALSE.
c        DO NR=1,NRT
c          IF(ITYP1(nr).GT.1) DEPENDENT=.TRUE.
c        ENDDO
c
c        IF(DEPENDENT) THEN
c          FMT='(/'' KTYP1 ='',I2,'' KTYP2 ='',I2,'' KTYP3 ='',I2,'
c     '        //''' KTYP4 ='',I2,'' KTYP5 ='',I2,'' KTYP6 ='',I2,/'
c     '        //''' KTYP7 ='',I2,'' KTYP8 ='',I2,'' KTYP9 ='',I2,'
c     '        //''' KTYP10='',I2,'' KTYP11='',I2,'' KTYP12='',I2,/'
c     '        //''' KTYP13='',I2,'' KTYP14='',I2)'
c          READ(IUNIT,FMT) KTYP1,KTYP2,KTYP3,KTYP4,KTYP5,KTYP6,
c     '      KTYP7,KTYP8,KTYP9,KTYP10,KTYP11,KTYP12,KTYP13,KTYP14
c          FMT='('' KTYP15='',I2,'' KTYP16='',I2,'' KTYP17='',I2,'
c     '       //''' KTYP18='',I2,'' KTYP19='',I2,'' KTYP20='',I2,/'
c     '       //''' KTYP21='',I2,'' KTYP22='',I2,'' KTYP23='',I2,'
c     '       //''' KTYP24='',I2,'' KTYP25='',I2,'' KTYP26='',I2,/'
c     '       //''' KTYP27='',I2,'' KTYP28='',I2,'' KTYP29='',I2)'
c          READ(IUNIT,FMT) KTYP15,KTYP16,KTYP17,KTYP18,KTYP19,
c     '      KTYP20,KTYP21,KTYP22,KTYP23,KTYP24,KTYP25,KTYP26,
c     '      KTYP27,KTYP28,KTYP29
c          FMT='('' KTYP30='',I2,'' KTYP31='',I2,'' KTYP32='',I2,'
c     '       //''' KTYP33='',I2,'' KTYP34='',I2,'' KTYP35='',I2)'
c          READ(IUNIT,FMT) KTYP30,KTYP31,KTYP32,KTYP33,KTYP34,
c     '      KTYP35
c          FMT='('' KTYP40='',I2,'' KTYP41='',I2,'' KTYP42='',I2,'
c     '       //''' KTYP43='',I2,'' KTYP44='',I2,'' KTYP45='',I2)'
c          READ(IUNIT,FMT) KTYP40,KTYP41,KTYP42,KTYP43,KTYP44,
c     '      KTYP45
c          FMT='('' KTYP50='',I2,'' KTYP51='',I2,'' KTYP52='',I2,'
c     '       //''' KTYP53='',I2,'' KTYP54='',I2,'' KTYP55='',I2,'
c     '      //'/'' KTYP56='',I2,'' KTYP57='',I2,'' KTYP58='',I2,'
c     '       //''' KTYP59='',I2,'' KTYP5A='',I2,'' KTYP5B='',I2,'
c     '      //'/'' KTYP5C='',I2,'' KTYP5D='',I2)'
c          READ(IUNIT,FMT) KTYP50,KTYP51,KTYP52,KTYP53,KTYP54,KTYP55,
c     '      KTYP56,KTYP57,KTYP58,KTYP59,KTYP5A,KTYP5B,KTYP5C,KTYP5D
c          FMT='('' KTYP60='',I2,'' KTYP61='',I2,'' KTYP62='',I2,'
c     '       //''' KTYP63='',I2,'' KTYP64='',I2,'' KTYP65='',I2)'
c          READ(IUNIT,FMT) KTYP60,KTYP61,KTYP62,KTYP63,KTYP64,
c     '      KTYP65
c          FMT='('' KTYP70='',I2,'' KTYP71='',I2,'' KTYP72='',I2,'
c     '       //''' KTYP73='',I2,'' KTYP74='',I2,'' KTYP75='',I2)'
c          READ(IUNIT,FMT) KTYP70,KTYP71,KTYP72,KTYP73,KTYP74,
c     '      KTYP75
c
c          FMT='(/'' IWRIT1='',I1,'' IWRIT2='',I1,'' IWRIT3='',I1,'
c     '      //''' IWRIT4='',I1)'
c          READ(IUNIT,FMT) IWRIT1,IWRIT2,IWRIT3,IWRIT4
c
c          FMT='(/'' DT='',E25.17,'' TINCR='',E25.17,'' T0='',E25.17,'
c     '      //''' T1='',E25.17,/'' THETA(1)='',E25.17,'' THETA(2)='','
c     '      //'E25.17,'//''' THETA(3)='',E25.17)'
c          READ(IUNIT,FMT) DT,TINCR,T0,T1,THETA(1),THETA(2),THETA(3)
c
c          FMT='(/'' INIT='',L1,'' LUMP='',L1,'' PROMPT='',L1,'
c     '      //''' LRESID='',L1,'' LSTEP='',L1,'' RESTAR='',L1)'
c          READ(IUNIT,FMT) INIT,LUMP,PROMPT,LRESID,LSTEP,RESTAR
c
c          IF(ITYP1(1).EQ.3) THEN !pdes
c            IE_TOT=1
c          ELSE IF(ITYP1(1).EQ.4) THEN !linear elasticity
c            IE_TOT=12
c          ELSE IF(ITYP1(1).EQ.5) THEN !finite elasticity
c            IE_TOT=1
c          ENDIF
c          DO IE=1,IE_TOT
c            FMT='(/'' Element type ie='',I3)'
c            READ(IUNIT,FMT) IEE
c            FMT='('' ETYP(ie)='',L1,'' ILT(ie)='',I3,'' IMT(ie)='','
c     '        //'I3,'//''' NGP(ie)='',I3,'' NLP(ie)='',I3,'
c     '        //''' NMP(ie)='',I3,'//''' NVE(ie)='',I3,'
c     '        //''' NHV(nv=1..6,ie):'',6I2)'
c            READ(IUNIT,FMT) ETYP(IE),
c     '        ILT(IE),IMT(IE),NGP(IE),NLP(IE),NMP(IE),
c     '        NVE(IE),(NHV(NV,IE),NV=1,6)
c            FMT='('' ILP(il=1..,ie):'',20I2)'
c            READ(IUNIT,FMT) (ILP(IL,IE),IL=1,ILT(IE))
c            FMT='('' NMB(il=1..,ie):'',20I2)'
c            READ(IUNIT,FMT) (NMB(IL,IE),IL=1,ILT(IE))
c          ENDDO
c          FMT='('' IT(1..3):'',3I3)'
c          READ(IUNIT,FMT) IT(1),IT(2),IT(3)
c
c          DO NR=1,NRT
c            DO NONODE=1,NPNODE(0,NR)
c              NP=NPNODE(NONODE,NR)
c              FMT='(/'' Region nr='',I2,''   Node np='',I4)'
c              READ(IUNIT,FMT) NRR,NPP
c              FMT='('' NHP(np)='',I2)'
c              READ(IUNIT,FMT) NHP(NP)
cC OLD AJP 3-9-92  FMT='('' NCNP(nh=1..,np,nr):'',6I3)'
cC                 READ(IUNIT,FMT) (NCNP(NH,NP,NR),NH=1,NHP(NP))
c              FMT='('' NCNP(0,nh=1..,np,nr):'',6I3)'
c              READ(IUNIT,FMT) (NCNP(0,NH,NP,NR),NH=1,NHP(NP))
c              DO NC=1,NCNP(0,NH,NP,NR)
c                FMT='('' NCNP(nc,nh=1..,np,nr):'',6I3)'
c                READ(IUNIT,FMT) (NCNP(NC,NH,NP,NR),NH=1,NHP(NP))
c              ENDDO
c              DO NC=1,NCT(nr)
c                FMT='('' NKH(nc='',I2,'',nh=1..,np):'',10I3)'
c                READ(IUNIT,FMT) NCC,(NKH(NC,NH,NP),NH=1,NHP(NP))
c              ENDDO
c              DO NC=1,NCT(nr)
c                DO NH=1,NHP(NP)
c                  FMT='('' NYNP(nk=1..,nh='',I2,'',nc='',I2,'
c     '              //''',np,nr):'',10I3)'
c                  READ(IUNIT,FMT) NHH,NCC,
c     '              (NYNP(NK,nc,NH,ncr,NP,NR),NK=1,NKH(NC,NH,NP))
c                ENDDO
c              ENDDO
c              FMT='('' CP(nm=1..,np):'',/(3E25.17))'
c              READ(IUNIT,FMT) (CP(NM,NP),NM=1,NMM)
c            ENDDO
c          ENDDO
c
c          DO NR=1,NRT
c            DO NOELEM=1,NEELEM(0,NR)
c              NE=NEELEM(NOELEM,NR)
c              FMT='(/'' Region nr='',I2,''   Element ne='',I4)'
c              READ(IUNIT,FMT) NRR,NEE
c              FMT='('' NHE(ne)='',I2)'
c              READ(IUNIT,FMT) NHE(NE)
c              DO NC=1,NCT(nr)
c                FMT='('' NBH(nc='',I2,'',nh=1..,ne):'',10I3)'
c                READ(IUNIT,FMT) NCC,(NBH(NC,NH,NE),NH=1,NHE(NE))
c              ENDDO
c              DO NC=1,NCT(nr)
c                DO NH=1,NHE(NE)
c                  NB=NBH(NC,NH,NE)
c                  FMT='('' NYNE(na=1..,nc='',I2,'',nh='',I2,'
c     '              //''',ne):'',10I3)'
c                  READ(IUNIT,FMT) NCC,NHH,
c     '              (NYNE(NA,NC,NH,NE),NA=1,NAT(NB))
c                ENDDO
c              ENDDO
c              FMT='('' NW(ne,1)='',I3)'
c              READ(IUNIT,FMT) NW(NE,1)
c              FMT='('' CE(nm=1..,ne):'',/(3E25.17))'
c              READ(IUNIT,FMT) (CE(NM,NE),NM=1,NMM)
c            ENDDO
c          ENDDO
c
c          DO NR=1,NRT
cC           IF(KTYP8.EQ.5) NYTOT=NJT*NPT(nr)*NKT(0,NR)
c            FMT='(X)'
c            READ(IUNIT,FMT)
c            DO NY=1,NYT(1,NR,1)
c              FMT='('' FIX('',I5,'',iy=1..,nr='',I1,''): '',5L1)'
c              READ(IUNIT,FMT) NYY,NRR,(FIX(NY,IY,NR),IY=1,5)
c              FMT='('' YP('',I5,'',iy=1..,nr='',I1,''): '',/,(3E25.17))'
c              READ(IUNIT,FMT) NYY,NRR,(YP(NY,IY,NR),IY=1,NIYM)
c            ENDDO
c          ENDDO
c
c          IF(KTYP57.EQ.2) THEN     !Boundary pressure increments entered
c            DO NR=1,NRT
c              FMT='(/'' Region '',I1,'':Pressure Boundary Conditions'')'
c              READ(IUNIT,FMT) NRR
c              DO NOELEM=1,NEELEM(0,NR)
c                NE=NEELEM(NOELEM,NR)
c                FMT='('' FIXP(iface=1..2,ne='',I3,''): '',2L1)'
c                READ(IUNIT,FMT) NEE,(FIXP(IFACE,NE),IFACE=1,2)
c                FMT='('' PE(iface=1..2,ne='',I3,''): '',2E25.17)'
c                READ(IUNIT,FMT) NEE,(PE(IFACE,NE),IFACE=1,2)
c              ENDDO
c              DO NOELEM=1,NEELEM(0,NR)
c                NE=NEELEM(NOELEM,NR)
c                FMT='(/'' PF(iface=1..2,ne='',I3,''): '',2E25.17)'
c                READ(IUNIT,FMT) NEE,(PF(IFACE,NE),IFACE=1,2)
c              ENDDO
c            ENDDO
c          ENDIF
c
c          FMT='(/'' FEXT(i=1..7,ng,ne):'')'
c          READ(IUNIT,FMT)
c          DO NR=1,NRT
c            DO NOELEM=1,NEELEM(0,NR)
c              NE=NEELEM(NOELEM,NR)
c              IF(NBH(1,1,NE).GT.0) THEN
c                DO NG=1,NGT(NBH(1,1,NE))
c                  FMT='('' FEXT(i,'',I3,'','',I4,''): '',3E25.17,'
c     '              //'/(4E25.17))'
c                  READ(IUNIT,FMT) MG,NE,(FEXT(I,NG,NE),I=1,7)
c                ENDDO
c              ENDIF
c            ENDDO
c          ENDDO
c
c          IF(KTYP53.EQ.3)THEN      !Active stress included
c            FMT='(/'' CALL_ACTI='',L1)'
c            READ(IUNIT,FMT) CALL_ACTI
c          ENDIF
c
c          IF(KTYP1.EQ.8.AND.KTYP3.EQ.2) THEN !activation variables
c            DO NR=1,NRT
c              DO NOELEM=1,NEELEM(0,NR)
c                NE=NEELEM(NOELEM,NR)
c                FMT='(/'' Element ne='',I4)'
c                READ(IUNIT,FMT) NE
c                NB=NBH(1,1,NE)
c                IF(NB.GT.0) THEN
c                  FMT='(/'' ITHRES(1,ng,ne):'',/(1X,75I1))'
c                  READ(IUNIT,FMT) (ITHRES(1,NG,NE),NG=1,NGT(NB))
c                  FMT='(/'' ITHRES(2,ng,ne):'',/(1X,75I1))'
c                  READ(IUNIT,FMT) (ITHRES(2,NG,NE),NG=1,NGT(NB))
c                  FMT='(/'' THRES(1,ng,ne):'',/(1X,10E12.4))'
c                  READ(IUNIT,FMT) (THRES(1,NG,NE),NG=1,NGT(NB))
c                  FMT='(/'' THRES(2,ng,ne):'',/(1X,10E12.4))'
c                  READ(IUNIT,FMT) (THRES(2,NG,NE),NG=1,NGT(NB))
c                  FMT='(/'' THRES(3,ng,ne):'',/(1X,10E12.4))'
c                  READ(IUNIT,FMT) (THRES(3,NG,NE),NG=1,NGT(NB))
c                ENDIF
c              ENDDO
c            ENDDO
c          ENDIF
c
c          FMT='(/'' CALL_EQUA='',L1,'' CALL_INIT='',L1,'
c     '      //'  '' CALL_MATE='',L1,'' CALL_SOLV='',L1)'
c          READ(IUNIT,FMT) CALL_EQUA,CALL_INIT,CALL_MATE,CALL_SOLV
c          FMT='(/'' CALL_DATA='',L1,'' CALL_FIT ='',L1,'
c     '      //'  '' CALL_GROW='',L1,'' CALL_MOTI='',L1)'
c          READ(IUNIT,FMT) CALL_DATA,CALL_FIT ,CALL_GROW,CALL_MOTI
c
c          IF(KTYP1.EQ.1) THEN !linear elasticity
c            IF(ETYP(7).OR.ETYP(8)) THEN !shell elements
c              NB=NBJ(1,1,NEELEM(1,1)) !must be cubic Hermite
c              CALL ASSERT(NKT(0,NB).EQ.4,
c     '          '>>Geometry must be cubic Hermite',ERROR,*9999)
cc D3PG is not passed through - PJH 5-Jan-1991
cc             CALL GAUS20(IBT(1,1,NB),IDO(1,0,NB),INP(1,1,NB),
cc    '          NB,NGAP(1,NB),D3PG(1,1,1,NB))
c            ENDIF
c          ENDIF
c
c        ENDIF !dependent
c
c        CALL DIMCHK(IOOP,nx,ERROR,*9999)
c
c        DO NB=1,NBT
c          IF(NBC(NB).EQ.1) THEN      !Lagrange/Hermite tensor prod basis
c            CALL GAUSS1(IBT(1,1,NB),IDO(1,0,NB),INP(1,1,NB),NB,
c     '        NGAP(1,NB),PG(1,1,1,NB),WG(1,NB),XIG(1,1,NB),ERROR,*9999)
c          ELSE IF(NBC(NB).EQ.2) THEN !Simplex/Serendipity/Lagrange basis
c            CALL GAUSS2(IBT,INP,NB,NGAP(1,NB),PG,WG,XIG,ERROR,*9999)
c          ELSE IF(NBC(NB).EQ.3) THEN !B-spline tensor product basis
c            WRITE(OP_STRING,'('' Warning: No call to Gauss3 from '','
c     '        //'''IOGEOM!!'')')
c      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c          ELSE IF(NBC(NB).EQ.4) THEN !Fourier Series tensor prod basis
c            CALL GAUSS4(IBT(1,1,NB),IDO(1,0,NB),INP(1,1,NB),NB,
c     '        NGAP(1,NB),PG(1,1,1,NB),WG(1,NB),XIG(1,1,NB),ERROR,*9999)
c          ELSE IF(NBC(NB).EQ.5) THEN !Boundary element basis
c            WRITE(OP_STRING,'('' Warning: No call to Gauss5 etc '','
c     '        //'''from IOGEOM'')')
c      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c          ELSE IF(NBC(NB).EQ.6) THEN !Singular basis function
c            WRITE(OP_STRING,'('' Warning: No call to Gauss6 etc '','
c     '        //'''from IOGEOM'')')
c      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c          ELSE IF(NBC(NB).EQ.7) THEN !Extended Lagrange basis
c            CALL GAUSS7(IBT(1,1,NB),IDO(1,0,NB),INP(1,1,NB),NB,
c     '        NGAP(1,NB),PG(1,1,1,NB),WG(1,NB),XIG(1,1,NB),ERROR,*9999)
c          ENDIF
c          IF(NBC(NB).LE.2) THEN !Auxillary basis
c            IF(NAT(NB).GT.0) THEN
c              CALL GAUSS8(NAN(1,1,NB),NB,NGAP(1,NB),PG(1,1,1,NB),
c     '          WG(1,NB),XIG(1,1,NB),ERROR,*9999)
c            ENDIF
c          ENDIF
c        ENDDO
c      ENDIF

      CALL EXITS('IPGEOM')
      RETURN
 9999 CALL ERRORS('IPGEOM',ERROR)
      CALL EXITS('IPGEOM')
      RETURN 1
      END


      SUBROUTINE LINOUT(IO_TYPE,IFREE,IFIXED,L_FORMAT,NDATA,L_DATA,
     '  L_DEFLT,INFO,ERROR,*)

C#### Subroutine: LINOUT
C###  Description:
C###    LINOUT handles input/output of logical values (.TRUE./.FALSE.)
C**** If IO_TYPE=0 prompts are written to IFREE
C****  "    "    1   "      "  written to IFREE and (with data) to IFIXED
C****  "    "    2   "      "  read with data from IFIXED
C****  "    "    3   "      "  written with data to IFIXED
C****  "    "    4   "      "  read with data from IFIXED and written to IFREE
C****  "    "    5   default values are used
C**** If IFIXED is not open no IO is performed on that file.
C**** Logical values separated by commas or blanks are read into
C**** L_DATA(NDAT) NDAT=1,NDATA.
C**** If an end of file or carriage return is detected while attempting
C**** to read L_DATA(NDAT) it is assigned the value L_DEFLT(NDAT).
C**** If an error is detected a diagnostic message is returned in ERROR
C**** and control is returned to the statment number of the asterisk.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
!     Parameter List
      INTEGER IFIXED,IFREE,INFO,IO_TYPE,NDATA
      CHARACTER ERROR*(*),L_FORMAT*(*)
      LOGICAL L_DATA(*),L_DEFLT(*)
!     Local Variables
      INTEGER CLOCAT,IBEG,IEND,IOSTAT,IPOS,irec,IUNIT,JPOS,m,MAXTRY,n,nd
      CHARACTER CHAR*1,CIOSTA*5,CIREC*5,CIUNIT*500,LINE*132,SUBSTR*1
      LOGICAL OPENED,REVERT
      DATA MAXTRY/8/

      CALL ENTERS('LINOUT')
      IF(NDATA.GT.IOLM) THEN
        ERROR=' NDATA exceeds the length of the array L_DATA'
        GOTO 9999
      ENDIF
      REVERT=.FALSE.
 10   IF(IO_TYPE.LE.1.OR.REVERT) THEN
        DO 1 m=1,MAXTRY
          CALL WRITE_LINE(IOOP,L_FORMAT,ERROR,*9999)
C ajp     WRITE(OP_STRING,L_FORMAT)
C 17-5-93 CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL LINPUT(IFREE,NDATA,L_DATA,L_DEFLT,INFO,ERROR,*9999)
          IF(ERROR.EQ.' ') GOTO 2
          IF(ERROR(1:2).NE.' ?') THEN
            WRITE(OP_STRING,*) ERROR
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
    1   CONTINUE
        ERROR=' More than 8 invalid input attempts made'
        GOTO 9999
    2   IF(IO_TYPE.EQ.1) THEN
          INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
          IF(OPENED.AND.IFIXED.NE.0) THEN
            WRITE(UNIT=IFIXED,FMT=L_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '        (L_DATA(n),n=1,NDATA)
            IF(IOSTAT.EQ.0) THEN
              ERROR=' '
            ELSE
              ERROR=' Logical data write error'
              GOTO 9999
            ENDIF
          ELSE
            ERROR=' '
          ENDIF
        ENDIF

      ELSE IF(IO_TYPE.EQ.2.OR.IO_TYPE.EQ.4) THEN
        IF(DOP) THEN
          WRITE(OP_STRING,'(A)') L_FORMAT(1:80)
      	  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(OPENED) THEN
 3        READ(UNIT=IFIXED,FMT='(A)',REC=irec,IOSTAT=IOSTAT) LINE
          IF(LINE(1:1).EQ.'!') THEN   !This line is a comment
            IF(LINE(2:2).EQ.'!') THEN !Do not print
            ELSE                      !Print comment line
              WRITE(OP_STRING,'(A)') LINE
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            irec=irec+1
            GO TO 3
          ENDIF
          SUBSTR=''''
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' SUBSTR=',SUBSTR
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IPOS=CLOCAT(SUBSTR,L_FORMAT)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' IPOS=',IPOS
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CHAR=' '
          IF(IPOS.GT.0) THEN
            SUBSTR='/'
            JPOS=CLOCAT(SUBSTR,L_FORMAT(1:4))
            IF(JPOS.GT.0) THEN
              READ(UNIT=IFIXED,FMT='(/A,'''//L_FORMAT(IPOS+2:),
     '          REC=irec,IOSTAT=IOSTAT) CHAR,(L_DATA(n),n=1,NDATA)
            ELSE
              READ(UNIT=IFIXED,FMT='(A,'''//L_FORMAT(IPOS+2:),
     '          REC=irec,IOSTAT=IOSTAT) CHAR,(L_DATA(n),n=1,NDATA)
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' CHAR='',A)') CHAR
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE
            READ(UNIT=IFIXED,FMT=L_FORMAT,
     '        REC=irec,IOSTAT=IOSTAT) (L_DATA(n),n=1,NDATA)
          ENDIF
          IF(IOSTAT.EQ.0) THEN
            ERROR=' '
          ELSE IF(IOSTAT.EQ.36) THEN
            ERROR='End of file'
            GO TO 9999
          ELSE
            ERROR=' Logical data read error'
            GOTO 9999
          ENDIF
          IF(CHAR.EQ.'*') THEN
            REVERT=.TRUE.
            GO TO 10
          ENDIF
        ELSE
          WRITE(UNIT=CIUNIT,FMT='(I10)') IUNIT
          CALL TRIM(CIUNIT,IBEG,IEND)
          ERROR=' Unit '//CIUNIT(IBEG:IEND)//' is not open'
          GOTO 9999
        ENDIF
        IF(IO_TYPE.EQ.4) THEN
          WRITE(OP_STRING,L_FORMAT) (L_DATA(n),n=1,NDATA)
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(IO_TYPE.EQ.3) THEN !write data to file
        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(OPENED.AND.IFIXED.NE.0) THEN
          WRITE(UNIT=IFIXED,FMT=L_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '      (L_DATA(n),n=1,NDATA)
          IF(IOSTAT.EQ.0) THEN
            ERROR=' '
          ELSE
            ERROR=' Logical data write error'
            GOTO 9999
          ENDIF
        ELSE
          ERROR=' '
        ENDIF

      ELSE IF(IO_TYPE.EQ.5) THEN !set to defaults
        DO nd=1,NDATA
          L_DATA(nd)=L_DEFLT(nd)
        ENDDO
      ENDIF

      CALL EXITS('LINOUT')
      RETURN

 9999 CALL TRIM(ERROR,IBEG,IEND)
      IF(ERROR(1:IEND).EQ.'End of file') THEN
        IO_TYPE=1
        WRITE(OP_STRING,'('' End of file: revert to interactive '','
     '    //'''input'')')
      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GOTO 10
      ELSE IF(ERROR(1:4).EQ.'Edit') THEN
        ERROR=ERROR(1:IEND)//'>LINOUT'
        CLOSE(UNIT=IFIXED)
      ELSE IF(ERROR(1:7).EQ.'Restart') THEN
        ERROR=ERROR(1:IEND)//'>LINOUT'
        CLOSE(UNIT=IFIXED,STATUS='DELETE')
      ELSE IF(ERROR(1:14).EQ.'File not found') THEN
        ERROR='Help file not found>LINOUT'
      ELSE
        WRITE(UNIT=CIOSTA,FMT='(I5)') IOSTAT
        WRITE(UNIT=CIREC ,FMT='(I5)') irec
        ERROR=ERROR(1:IEND)//': IOSTAT='//CIOSTA(1:5)//' irec='//
     '    CIREC(1:5)//' L_FORMAT='//L_FORMAT(1:10)//'>LINOUT'
        CLOSE(UNIT=IFIXED)
      ENDIF
      CALL EXITS('LINOUT')
      RETURN 1
      END


      SUBROUTINE LINPUT(IUNIT,NDATA,LDATA,LDEFLT,INFO,ERROR,*)

C#### Subroutine: LINPUT
C###  Description:
C###    LINPUT reads logical values (true/false) from IUNIT.
C**** Logical values,separated by blanks or commas are read into
C**** LDATA(NDAT) NDAT=1,NDATA.
C**** If an end of file or carriage return is detected while attempting
C**** to read LDATA(NDAT) it is assigned the value ADEFLT(NDAT).
C**** If a help request is made (?) a help file (based on INFO) is opened.
C**** If a stop request is made (STOP/stop/S/s) the program is stopped.
C**** If an edit request is made (EDIT/edit/E/e) control returns to
C**** the command line and the file is closed.
C**** If a restart request is made (RESTART/restart/R/r) control returns to
C**** the command line and the file is deleted.
C**** If an error is detected or a help request is made ERROR is
C**** returned with a diagnostic message.
C**** NOTE: The Q edit descriptor in the read format is nonstandard!

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER INFO,IUNIT,NDATA
      CHARACTER ERROR*(*)
      LOGICAL LDATA(*),LDEFLT(*)
!     Local Variables
      INTEGER ic,ICOL,IEND,IOSTAT,irec,ISTART,NCHAR,nd,NDAT
      CHARACTER CDATA*500,CHAR*5,CLINE*500,FILE*6

      CALL ENTERS('LINPUT')
      NDAT=0
      DO 4 irec=1,NDATA
        CALL READ_LINE(IUNIT,IOSTAT,NCHAR,CLINE,ERROR,*9999)
        IF(IOSTAT.EQ.0.AND.NCHAR.GT.0) THEN
          ICOL=1
          DO 3 nd=NDAT+1,NDATA
            ISTART=0
            DO 1 ic=ICOL,NCHAR
              IF((CLINE(ic:ic).EQ.' ').OR.(CLINE(ic:ic).EQ.',')) THEN
                IF(ISTART.GT.0) THEN
                  IEND=ic-1
                  GOTO 2
                ENDIF
              ELSE
                IF(ISTART.EQ.0) THEN
                  ISTART=ic
                ENDIF
              ENDIF
    1       CONTINUE
            IEND=NCHAR
    2       IF(ISTART.GT.0) THEN
              ICOL=IEND+2
              CDATA=CLINE(ISTART:IEND)
              IF((CDATA.EQ.'EDIT').OR.
     '           (CDATA.EQ.'edit').OR.
     '           (CDATA.EQ.'E').OR.
     '           (CDATA.EQ.'e'))THEN
                ERROR='Edit'
                GO TO 9999
              ELSE IF((CDATA.EQ.'RESTART').OR.
     '                (CDATA.EQ.'restart').OR.
     '                (CDATA.EQ.'R').OR.
     '                (CDATA.EQ.'r'))THEN
                ERROR='Restart'
                GO TO 9999
              ELSE IF((CDATA.EQ.'STOP').OR.
     '                (CDATA.EQ.'stop').OR.
     '                (CDATA.EQ.'S').OR.
     '                (CDATA.EQ.'s'))THEN
                CLOSE(UNIT=IUNIT)
                STOP
              ELSE IF(CDATA.EQ.'   ?') THEN
                WRITE(UNIT=CHAR,FMT='(I5)') INFO
                FILE='FE'//CHAR(1:4)
                CALL DOCUM(FILE,'doc','Parameter number '//CHAR(5:5),
     '            ERROR,*9997)
 9997           IF(ERROR.EQ.' ') THEN
                  ERROR=' ?'
                  GO TO 9998
                ELSE IF(ERROR(1:14).EQ.'File not found') THEN
                  WRITE(OP_STRING,*) '>>Help not available'
      		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ERROR=' ?'
                  GO TO 9998
                ELSE
                  GO TO 9999
                ENDIF
              ELSE
                IF((CDATA.EQ.'TRUE').OR.
     '             (CDATA.EQ.'true').OR.
     '             (CDATA.EQ.'T').OR.
     '             (CDATA.EQ.'t')) THEN
                  NDAT=NDAT+1
                  LDATA(NDAT)=.TRUE.
                  IF(NDAT.EQ.NDATA) THEN

                    ERROR=' '
                    GO TO 9998
                  ENDIF
                ELSE IF((CDATA.EQ.'FALSE').OR.
     '                  (CDATA.EQ.'false').OR.
     '                  (CDATA.EQ.'F').OR.
     '                  (CDATA.EQ.'f')) THEN
                  NDAT=NDAT+1
                  LDATA(NDAT)=.FALSE.
                  IF(NDAT.EQ.NDATA) THEN
                    ERROR=' '
                    GO TO 9998
                  ENDIF
                ELSE
                  ERROR=' Answer data range error'
                  GO TO 9998
                ENDIF
              ENDIF
            ELSE
              GOTO 4
            ENDIF
    3     CONTINUE
        ELSE IF(IOSTAT.EQ.0.AND.NCHAR.EQ.0) THEN
C         REWIND(UNIT=IUNIT)
          NDAT=NDAT+1
C         LDATA(NDAT)=LDEFLT(NDAT)
C         IF(NDAT.EQ.NDATA) THEN
          DO nd=1,NDATA
            LDATA(nd)=LDEFLT(nd)
          ENDDO
          ERROR=' '
          GO TO 9998
C         ENDIF
        ELSE
          REWIND(UNIT=IUNIT)
          ERROR=' Character data read error'
          GO TO 9998
        ENDIF
    4 CONTINUE

      ERROR=' '
 9998 CALL EXITS('LINPUT')
      RETURN
 9999 CALL ERRORS('LINPUT',ERROR)
      CALL EXITS('LINPUT')
      RETURN 1
      END


      SUBROUTINE RINOUT(IO_TYPE,IFREE,IFIXED,R_FORMAT,NDATA,R_DATA,
     '  R_DEFLT,RMIN,RMAX,INFO,ERROR,*)

C#### Subroutine: RINOUT
C###  Description:
C###    RINOUT handles input/output of real*8 values.
C**** If IO_TYPE=0 prompts written to IFREE
C****  "    "    1   "     written to IFREE & (with data) to IFIXED
C****  "    "    2   "     read with data from IFIXED
C****  "    "    3   "     written with data to IFIXED
C****  "    "    4   "     read with data from IFIXED & written to IFREE
C****  "    "    5   default values are used
C**** If IFIXED is not open no IO is performed on that file.
C**** Real values separated by commas or blanks are read into
C**** R_DATA(NDAT) NDAT=1,NDATA.
C**** The limiting values for the read are defined by RMIN and RMAX.
C**** If an end of file or carriage return is detected while attempting
C**** to read R_DATA(NDAT) it is assigned the value R_DEFLT(NDAT).
C**** If an error is detected a diagnostic message is returned in ERROR
C**** and control is returned to the statment number of the asterisk.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
!     Parameter List
      INTEGER IFIXED,IFREE,IO_TYPE,INFO,NDATA
      REAL*8 RMAX,RMIN,R_DATA(*),R_DEFLT(*)
      CHARACTER ERROR*(*),R_FORMAT*(*)
!     Local Variables
      INTEGER CLOCAT,IBEG,IEND,IOSTAT,IPOS,irec,JPOS,m,MAXTRY,n,nd
      CHARACTER CHAR*1,CIOSTA*5,CIREC*5,CIUNIT*10,LINE*132,SUBSTR*1
      LOGICAL OPENED,REVERT
      DATA MAXTRY/8/

      CALL ENTERS('RINOUT')
      IF(NDATA.GT.IORM) THEN
        ERROR=' NDATA exceeds the length of the array R_DATA'
        GOTO 9999
      ENDIF
      REVERT=.FALSE.
 10   IF(IO_TYPE.LE.1.OR.REVERT) THEN
        DO 1 m=1,MAXTRY
          CALL WRITE_LINE(IOOP,R_FORMAT,ERROR,*9999)
C ajp     WRITE(OP_STRING,R_FORMAT)
C 17-5-93 CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL RINPUT(IFREE,NDATA,R_DATA,R_DEFLT,RMIN,RMAX,
     '      INFO,ERROR,*9999)
          IF(ERROR.EQ.' ') GOTO 2
          IF(ERROR(1:2).NE.' ?') THEN
            WRITE(OP_STRING,*) ERROR
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
    1   CONTINUE
        ERROR=' More than 8 invalid input attempts made'
        GOTO 9999
    2   IF(IO_TYPE.EQ.1) THEN
          INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
          IF(OPENED.AND.IFIXED.NE.0) THEN
            WRITE(UNIT=IFIXED,FMT=R_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '        (R_DATA(n),n=1,NDATA)
            IF(IOSTAT.EQ.0) THEN
              ERROR=' '
            ELSE
              ERROR=' Real data write error'
              GOTO 9999
            ENDIF
          ELSE
            ERROR=' '
          ENDIF
        ENDIF

      ELSE IF(IO_TYPE.EQ.2.OR.IO_TYPE.EQ.4) THEN
        IF(DOP) THEN
          WRITE(OP_STRING,'(A)') R_FORMAT(1:80)
      	  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(OPENED) THEN
 3        READ(UNIT=IFIXED,FMT='(A)',REC=irec,IOSTAT=IOSTAT) LINE
          IF(LINE(1:1).EQ.'!') THEN   !This line is a comment
            IF(LINE(2:2).EQ.'!') THEN !Do not print
            ELSE                      !Print comment line
              WRITE(OP_STRING,'(A)') LINE
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            irec=irec+1
            GO TO 3
          ENDIF
          SUBSTR=''''
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' SUBSTR=',SUBSTR
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IPOS=CLOCAT(SUBSTR,R_FORMAT)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' IPOS=',IPOS
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CHAR=' '
          IF(IPOS.GT.0) THEN
            SUBSTR='/'
            JPOS=CLOCAT(SUBSTR,R_FORMAT(1:4))
            IF(JPOS.GT.0) THEN
              READ(UNIT=IFIXED,FMT='(/A,'''//R_FORMAT(IPOS+2:),
     '          REC=irec,IOSTAT=IOSTAT) CHAR,(R_DATA(n),n=1,NDATA)
            ELSE
              READ(UNIT=IFIXED,FMT='(A,'''//R_FORMAT(IPOS+2:),
     '          REC=irec,IOSTAT=IOSTAT) CHAR,(R_DATA(n),n=1,NDATA)
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' CHAR='',A)') CHAR
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE
            READ(UNIT=IFIXED,FMT=R_FORMAT,
     '        REC=irec,IOSTAT=IOSTAT) (R_DATA(n),n=1,NDATA)
          ENDIF
          IF(IOSTAT.EQ.0) THEN
            ERROR=' '
          ELSE IF(IOSTAT.EQ.36) THEN
            ERROR='End of file'
            GO TO 9999
          ELSE
            ERROR=' Real data read error'
            GOTO 9999
          ENDIF
          IF(CHAR.EQ.'*') THEN
            REVERT=.TRUE.
            GO TO 10
          ENDIF
        ELSE
          WRITE(UNIT=CIUNIT,FMT='(I10)') IFIXED
          CALL TRIM(CIUNIT,IBEG,IEND)
          ERROR=' Unit '//CIUNIT(IBEG:IEND)//' is not open'
          GOTO 9999
        ENDIF
        IF(IO_TYPE.EQ.4) THEN
          WRITE(OP_STRING,R_FORMAT) (R_DATA(n),n=1,NDATA)
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(IO_TYPE.EQ.3) THEN !write data to file
        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(OPENED.AND.IFIXED.NE.0) THEN
          WRITE(UNIT=IFIXED,FMT=R_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '      (R_DATA(n),n=1,NDATA)
          IF(IOSTAT.EQ.0) THEN
            ERROR=' '
          ELSE
            ERROR=' Real data write error'
            GOTO 9999
          ENDIF
        ELSE
          ERROR=' '
        ENDIF

      ELSE IF(IO_TYPE.EQ.5) THEN !set to defaults
        DO nd=1,NDATA
          R_DATA(nd)=R_DEFLT(nd)
        ENDDO
      ENDIF

      CALL EXITS('RINOUT')
      RETURN

 9999 CALL TRIM(ERROR,IBEG,IEND)
      IF(ERROR(1:IEND).EQ.'End of file') THEN
        IO_TYPE=1
        WRITE(OP_STRING,'('' End of file: revert to interactive '','
     '    //'''input'')')
      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GOTO 10
      ELSE IF(ERROR(1:4).EQ.'Edit') THEN
        ERROR=ERROR(1:IEND)//'>RINOUT'
        CLOSE(UNIT=IFIXED)
      ELSE IF(ERROR(1:7).EQ.'Restart') THEN
        ERROR=ERROR(1:IEND)//'>RINOUT'
        CLOSE(UNIT=IFIXED,STATUS='DELETE')
      ELSE IF(ERROR(1:14).EQ.'File not found') THEN
        ERROR='Help file not found>RINOUT'
      ELSE
        WRITE(UNIT=CIOSTA,FMT='(I5)') IOSTAT
        WRITE(UNIT=CIREC ,FMT='(I5)') irec
        ERROR=ERROR(1:IEND)//': IOSTAT='//CIOSTA(1:5)//' irec='//
     '    CIREC(1:5)//' R_FORMAT='//R_FORMAT(1:10)//'>RINOUT'
        CLOSE(UNIT=IFIXED)
      ENDIF
      CALL EXITS('RINOUT')
      RETURN 1
      END


      SUBROUTINE RINOUT1(IFREE,IFIXED,N1char,N2char,R_FORMAT,
     '  NDATA,R_DATA,R_DEFLT,RMIN,RMAX,FILEIP,INFO,ERROR,*)

C#### Subroutine: RINOUT1
C###  Description:
C###    RINOUT handles input/output of real*8 values.
C**** If IO_TYPE=0 prompts are written to IFREE
C****  "    "    1   "      "  written to IFREE and (with data) to IFIXED
C****  "    "    2   "      "  read with data from IFIXED
C****  "    "    3   "      "  written with data to IFIXED
C****  "    "    4   "      "  read with data from IFIXED and written to IFREE
C**** If IFIXED is not open no IO is performed on that file.
C**** If NTCHAR is >0 it indicates the position in the R_FORMAT string at the
C**** end of the prompt.
C**** Real values separated by commas or blanks are read into
C**** R_DATA(NDAT) NDAT=1,NDATA.
C**** The limiting values for the read are defined by RMIN and RMAX.
C**** If an end of file or carriage return is detected while attempting
C**** to read R_DATA(NDAT) it is assigned the value R_DEFLT(NDAT).
C**** FILEIP is .true. if response indicates file input.
C**** If an error is detected a diagnostic message is returned in ERROR
C**** and control is returned to the statment number of the asterisk.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
!     Parameter List
      INTEGER IFIXED,IFREE,INFO,N1char,N2char,NDATA
      REAL*8 RMAX,RMIN,R_DATA(*),R_DEFLT(*)
      CHARACTER ERROR*(*),R_FORMAT*(*)
      LOGICAL FILEIP
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,IOSTAT,irec,m,MAXTRY,n,N3CO
      REAL*8 RFROMC
      CHARACTER CHAR*1,CHAR1*20,CIOSTA*5,CIREC*5,CIUNIT*10,LINE*132
      LOGICAL CBBREV,FILEIP_FIRST,OPENED,REVERT
      DATA MAXTRY/8/

      CALL ENTERS('RINOUT1')
      IF(NDATA.GT.IORM) THEN
        ERROR=' NDATA exceeds the length of the array R_DATA'
        GOTO 9999
      ENDIF
      REVERT=.FALSE.
 10   IF(IOTYPE.LE.1.OR.REVERT) THEN
        DO 1 m=1,MAXTRY
          CALL WRITE_LINE(IOOP,R_FORMAT,ERROR,*9999)
C ajp     WRITE(OP_STRING,R_FORMAT)
C 17-5-93 CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL RINPUT1(IFREE,NDATA,R_DATA,R_DEFLT,RMIN,RMAX,
     '      FILEIP,FILEIP_FIRST,INFO,ERROR,*9999)
          IF(ERROR.EQ.' ') GOTO 2
          IF(ERROR(1:2).NE.' ?') THEN
            WRITE(OP_STRING,*) ERROR
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
    1   CONTINUE
        ERROR=' More than 8 invalid input attempts made'
        GOTO 9999
    2   IF(IOTYPE.EQ.1) THEN
          INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
          IF(OPENED.AND.IFIXED.NE.0) THEN
            IF(FILEIP) THEN
              R_FORMAT=R_FORMAT(1:N2char)//'FILE'
              WRITE(UNIT=IFIXED,FMT='(A)',REC=irec,IOSTAT=IOSTAT)
     '          R_FORMAT(N1char:N2char+4)
              IF(FILEIP_FIRST) THEN
                WRITE(UNIT=IFIXED,FMT='('' Enter FILENAME: '',A)',
     '            REC=irec+1,IOSTAT=IOSTAT) FILE07
              ENDIF
            ELSE IF(.NOT.FILEIP) THEN
              WRITE(UNIT=IFIXED,FMT=R_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '          (R_DATA(n),n=1,NDATA)
            ENDIF
            IF(IOSTAT.EQ.0) THEN
              ERROR=' '
            ELSE
              ERROR=' Real data write error'
              GOTO 9999
            ENDIF
          ELSE
            ERROR=' '
          ENDIF
        ENDIF

      ELSE IF(IOTYPE.EQ.2.OR.IOTYPE.EQ.4) THEN
        IF(DOP) THEN
          WRITE(OP_STRING,'(A)') R_FORMAT(1:80)
      	  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(OPENED) THEN
 3        READ(UNIT=IFIXED,FMT='(A)',REC=irec,IOSTAT=IOSTAT) LINE
          IF(LINE(1:1).EQ.'!') THEN   !This line is a comment
            IF(LINE(2:2).EQ.'!') THEN !Do not print
            ELSE                      !Print comment line
              WRITE(OP_STRING,'(A)') LINE
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            irec=irec+1
            GO TO 3
          ENDIF
c         SUBSTR=''''
c         IF(DOP) THEN
c           WRITE(OP_STRING,*) ' SUBSTR=',SUBSTR
c     	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c         ENDIF
c         IPOS=CLOCAT(SUBSTR,R_FORMAT)
c         IF(DOP) THEN
c	    WRITE(OP_STRING,*) ' IPOS=',IPOS
c      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c         ENDIF
          CHAR=' '
c         IF(IPOS.GT.0) THEN
c           SUBSTR='/'
c           JPOS=CLOCAT(SUBSTR,R_FORMAT(1:4))
c           IF(JPOS.GT.0) THEN
c             READ(UNIT=IFIXED,FMT='(/A,'''//R_FORMAT(IPOS+2:),
c    '          REC=irec,IOSTAT=IOSTAT) CHAR,(R_DATA(n),n=1,NDATA)
c           ELSE
c             READ(UNIT=IFIXED,FMT='(A,'''//R_FORMAT(IPOS+2:),
c    '          REC=irec,IOSTAT=IOSTAT) CHAR,(R_DATA(n),n=1,NDATA)
c           ENDIF
c           IF(DOP) THEN
c	      WRITE(OP_STRING,'('' CHAR='',A)') CHAR
c             CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c           ENDIF
c         ELSE
            IF(N1char.EQ.0) THEN
              READ(UNIT=IFIXED,FMT=R_FORMAT,
     '          REC=irec,IOSTAT=IOSTAT) (R_DATA(n),n=1,NDATA)
            ELSE IF(N1char.GT.0) THEN
!old          IF(DOP) THEN
!               WRITE(OP_STRING,*) ' FMT=',R_FORMAT(1:N2char)//''',A)'
!               CALL WRITES(IODI,OP_STRING,ERROR,*9999)
!             ENDIF
              READ(UNIT=IFIXED,FMT=R_FORMAT(1:N2char)//''',A)',
     '          REC=irec,IOSTAT=IOSTAT) CHAR1(1:20)
              CALL TRIM(CHAR1,IBEG1,IEND1)
              IF(DOP) THEN
                WRITE(OP_STRING,*) ' CHAR1=',CHAR1(IBEG1:IEND1)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(CBBREV(CHAR1(IBEG1:IEND1),'FILE',1,1,20,N3CO)) THEN
                IF(.NOT.FILEIP) THEN
                  FILEIP=.TRUE.
                  READ(UNIT=IFIXED,FMT='(17X,A)',REC=irec+1,
     '              IOSTAT=IOSTAT)FILE07
                  IF(DOP) THEN
                    WRITE(OP_STRING,*) ' FILE07=',FILE07
      		    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  CALL TRIM(FILE07,IBEG,IEND)
                  CALL OPENF(7,'DISK',FILE07(IBEG:IEND)//'.data','OLD',
     '              'SEQUEN','FORMATTED',132,ERROR,*9999)
                ENDIF
                READ(7,*) (R_DATA(n),n=1,NDATA)
              ELSE
                R_DATA(1)=RFROMC(CHAR1(IBEG1:IEND1))
              ENDIF
            ENDIF
c         ENDIF
          IF(IOSTAT.EQ.0) THEN
            ERROR=' '
          ELSE IF(IOSTAT.EQ.36) THEN
            ERROR='End of file'
            GO TO 9999
          ELSE
            ERROR=' Real data read error'
            GOTO 9999
          ENDIF
          IF(CHAR.EQ.'*') THEN
            REVERT=.TRUE.
            GO TO 10
          ENDIF
        ELSE
          WRITE(UNIT=CIUNIT,FMT='(I10)') IFIXED
          CALL TRIM(CIUNIT,IBEG,IEND)
          ERROR=' Unit '//CIUNIT(IBEG:IEND)//' is not open'
          GOTO 9999
        ENDIF
        IF(IOTYPE.EQ.4) THEN
          WRITE(OP_STRING,R_FORMAT) (R_DATA(n),n=1,NDATA)
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(IOTYPE.EQ.3) THEN
        INQUIRE(UNIT=IFIXED,OPENED=OPENED,NEXTREC=irec)
        IF(OPENED.AND.IFIXED.NE.0) THEN
          WRITE(UNIT=IFIXED,FMT=R_FORMAT,REC=irec,IOSTAT=IOSTAT)
     '      (R_DATA(n),n=1,NDATA)
          IF(IOSTAT.EQ.0) THEN
            ERROR=' '
          ELSE
            ERROR=' Real data write error'
            GOTO 9999
          ENDIF
        ELSE
          ERROR=' '
        ENDIF
      ENDIF

      CALL EXITS('RINOUT1')
      RETURN

 9999 CALL TRIM(ERROR,IBEG,IEND)
      IF(ERROR(1:IEND).EQ.'End of file') THEN
        IOTYPE=1
        WRITE(OP_STRING,'('' End of file: revert to interactive '','
     '    //'''input'')')
      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GOTO 10
      ELSE IF(ERROR(1:4).EQ.'Edit') THEN
        ERROR=ERROR(1:IEND)//'>RINOUT1'
        CLOSE(UNIT=IFIXED)
      ELSE IF(ERROR(1:7).EQ.'Restart') THEN
        ERROR=ERROR(1:IEND)//'>RINOUT1'        
        CLOSE(UNIT=IFIXED,STATUS='DELETE')
      ELSE IF(ERROR(1:14).EQ.'File not found') THEN
        ERROR='Help file not found>RINOUT1'
      ELSE
        WRITE(UNIT=CIOSTA,FMT='(I5)') IOSTAT
        WRITE(UNIT=CIREC ,FMT='(I5)') irec
        ERROR=ERROR(1:IEND)//': IOSTAT='//CIOSTA(1:5)//' irec='//
     '    CIREC(1:5)//' R_FORMAT='//R_FORMAT(1:10)//'>RINOUT1'
        CLOSE(UNIT=IFIXED)
      ENDIF
      CALL EXITS('RINOUT1')
      RETURN 1
      END


      SUBROUTINE RINPUT(IUNIT,NDATA,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*)

C#### Subroutine: RINPUT
C###  Description:
C###    RINPUT reads real*8 values from IUNIT.
C**** Reals separated by commas or blanks are read into
C**** RDATA(NDAT) NDAT=1,NDATA.
C**** If an end of file or carriage return is detected while attempting
C**** to read RDATA(NDAT) it is assigned the value RDEFLT(NDAT).
C**** The limiting values for the read are defined by RMIN and RMAX.
C**** If help request is made (?) a help file (based on INFO) is opened.
C**** If stop request is made (STOP/stop/S/s) the program is stopped.
C**** If an edit request is made (EDIT/edit/E/e) control returns to
C**** the command line and the file is closed.
C**** If restart request is made (RESTART/restart/R/r) control returns 
C**** to the command line and the file is deleted.
C**** If a file request is made (FILE/file/F/f) the program prompts for
C**** the file name (FILE07) and data is read from FILE07.iphist during
C**** a time-dependent iteration in subroutine MARCH1.
C**** If a file request is made (GET/get/G/g) the program prompts for
C**** the file name (FILE.ext) and the data is read from this file
C**** immediately.
C**** If an error is detected or a help request is made ERROR is
C**** returned with a diagnostic message.
C**** NOTE: The Q edit descriptor in the read format is nonstandard!

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:file01.cmn'
!     Parameter List
      INTEGER INFO,IUNIT,NDATA
      REAL*8 RDATA(*),RDEFLT(*),RMAX,RMIN
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,ic,ICOL,IEND,IOSTAT,irec,ISTART,LENGTH,NCHAR,nd,NDAT
      CHARACTER CDATA*24,CHAR*5,CLINE*500,FILE*6
      LOGICAL RVALID

      CALL ENTERS('RINPUT')
      FILEIP=.FALSE.
      NDAT=0
      DO 4 irec=1,NDATA
        CALL READ_LINE(IUNIT,IOSTAT,NCHAR,CLINE,ERROR,*9999)
        IF(IOSTAT.EQ.0.AND.NCHAR.GT.0) THEN
          ICOL=1
          DO 3 nd=NDAT+1,NDATA
            ISTART=0
            DO 1 ic=ICOL,NCHAR
              IF((CLINE(ic:ic).EQ.' ').OR.(CLINE(ic:ic).EQ.',')) THEN
                IF(ISTART.GT.0) THEN
                  IEND=ic-1
                  GOTO 2
                ENDIF
              ELSE
                IF(ISTART.EQ.0) THEN
                  ISTART=ic
                ENDIF
              ENDIF
    1       CONTINUE
            IEND=NCHAR
    2       IF(ISTART.GT.0) THEN
              ICOL=IEND+2
              CDATA=' '
              LENGTH=IEND-ISTART+1
              CDATA(25-LENGTH:24)=CLINE(ISTART:IEND)
              IF((CDATA.EQ.'                    EDIT').OR.
     '           (CDATA.EQ.'                    edit').OR.
     '           (CDATA.EQ.'                       E').OR.
     '           (CDATA.EQ.'                       e'))THEN
                ERROR='Edit'
                GO TO 9999
              ELSE IF((CDATA.EQ.'                 RESTART').OR.
     '                (CDATA.EQ.'                 restart').OR.
     '                (CDATA.EQ.'                       R').OR.
     '                (CDATA.EQ.'                       r'))THEN
                ERROR='Restart'
                GO TO 9999
              ELSE IF((CDATA.EQ.'                    STOP').OR.
     '                (CDATA.EQ.'                    stop').OR.
     '                (CDATA.EQ.'                       S').OR.
     '                (CDATA.EQ.'                       s'))THEN
                CLOSE(UNIT=IUNIT)
                STOP
              ELSE IF((CDATA.EQ.'                    FILE').OR.
     '                (CDATA.EQ.'                    file').OR.
     '                (CDATA.EQ.'                       F').OR.
     '                (CDATA.EQ.'                       f'))THEN
                WRITE(IUNIT,'($,'' Enter FILENAME(.iphist): '')')
                READ(IUNIT,'(A)') FILE07
                FILEIP=.TRUE.
              ELSE IF((CDATA.EQ.'                     GET').OR.
     '                (CDATA.EQ.'                     get').OR.
     '                (CDATA.EQ.'                       G').OR.
     '                (CDATA.EQ.'                       g'))THEN
                WRITE(IUNIT,'($,'' Enter FILENAME(.data): '')')
                READ(IUNIT,'(A)') FILE07
                CALL TRIM(FILE07,IBEG,IEND)
                CALL OPENF(IOFILE2,'DISK',FILE07(IBEG:IEND)//'.data',
     '            'OLD','SEQUEN','FORMATTED',132,ERROR,*9999)
                NDAT=NDAT+1
                READ(IOFILE2,*) RDATA(NDAT)
                CALL CLOSEF(IOFILE2,ERROR,*9999)
                FILEIP=.TRUE.
                ERROR=' '
                GO TO 9998
              ELSE IF(CDATA.EQ.'                       ?') THEN
                WRITE(UNIT=CHAR,FMT='(I5)') INFO
                FILE='FE'//CHAR(1:4)
                CALL DOCUM(FILE,'doc','Parameter number '//CHAR(5:5),
     '            ERROR,*9997)
 9997           IF(ERROR.EQ.' ') THEN
                  ERROR=' ?'
                  GO TO 9998
                ELSE IF(ERROR(1:14).EQ.'File not found') THEN
                  WRITE(OP_STRING,*) '>>Help not available'
      		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ERROR=' ?'
                  GO TO 9998
                ELSE
                  GO TO 9999
                ENDIF
              ELSE
                IF(RVALID(CDATA)) THEN
                  NDAT=NDAT+1
                  READ(UNIT=CDATA,FMT='(G24.0)',IOSTAT=IOSTAT)
     '              RDATA(NDAT)
                  IF(IOSTAT.EQ.0) THEN
                    IF((RDATA(NDAT).GE.RMIN).AND.
     '                 (RDATA(NDAT).LE.RMAX)) THEN
                      IF(NDAT.EQ.NDATA) THEN
                        ERROR=' '
                        GO TO 9998
                      ENDIF
                    ELSE
                      ERROR=' Real data range error'
                      GO TO 9998
                    ENDIF
                  ELSE
                    ERROR=' Real data read error'
                    GO TO 9998
                  ENDIF
                ELSE
                  ERROR=' Real data type error'
                  GO TO 9998
                ENDIF
              ENDIF
            ELSE
              GOTO 4
            ENDIF
    3     CONTINUE
        ELSE IF(IOSTAT.EQ.0.AND.NCHAR.EQ.0) THEN
C         REWIND(UNIT=IUNIT)
          NDAT=NDAT+1
C         RDATA(NDAT)=RDEFLT(NDAT)
C         IF(NDAT.EQ.NDATA) THEN
          DO nd=1,NDATA
            RDATA(nd)=RDEFLT(nd)
          ENDDO
          ERROR=' '
          GO TO 9998
C         ENDIF
        ELSE
          REWIND(UNIT=IUNIT)
          ERROR=' Character data read error'
          GO TO 9998
        ENDIF
    4 CONTINUE

      ERROR=' '
 9998 CALL EXITS('RINPUT')
      RETURN
 9999 CALL ERRORS('RINPUT',ERROR)
      CALL EXITS('RINPUT')
      RETURN 1
      END


      SUBROUTINE RINPUT1(IUNIT,NDATA,RDATA,RDEFLT,RMIN,RMAX,
     '  FILEIP,FILEIP_FIRST,INFO,ERROR,*)

C#### Subroutine: RINPUT1
C###  Description:
C###    RINPUT reads real*8 values from IUNIT.
C**** Reals separated by commas or blanks are read into
C**** RDATA(NDAT) NDAT=1,NDATA.
C**** If an end of file or carriage return is detected while attempting
C**** to read RDATA(NDAT) it is assigned the value RDEFLT(NDAT).
C**** The limiting values for the read are defined by RMIN and RMAX.
C**** If a help request is made (?) a help file (based on INFO) is opened.
C**** If a stop request is made (STOP/stop/S/s) the program is stopped.
C**** If an edit request is made (EDIT/edit/E/e) control returns to
C**** the command line and the file is closed.
C**** If a restart request is made (RESTART/restart/R/r) control returns to
C**** the command line and the file is deleted.
C**** If a file request is made (FILE/file/F/f) the program prompts for
C**** the file name (FILE07) and FILEIP is set .true.
C**** If a file request is made (GET/get/G/g) the program prompts for
C**** the file name (FILE07) (first time only) and the data is read from
C**** this file immediately.
C**** If an error is detected or a help request is made ERROR is
C**** returned with a diagnostic message.
C**** NOTE: The Q edit descriptor in the read format is nonstandard!

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
!     Parameter List
      INTEGER INFO,IUNIT,NDATA
      REAL*8 RDATA(*),RDEFLT(*),RMAX,RMIN
      CHARACTER ERROR*(*)
      LOGICAL FILEIP,FILEIP_FIRST
!     Local Variables
      INTEGER IBEG,ic,ICOL,IEND,IOSTAT,irec,ISTART,LENGTH,NCHAR,nd,NDAT
      CHARACTER CDATA*24,CHAR*5,CLINE*500,FILE*6
      LOGICAL RVALID

      CALL ENTERS('RINPUT1')
      NDAT=0
      DO 4 irec=1,NDATA
        CALL READ_LINE(IUNIT,IOSTAT,NCHAR,CLINE,ERROR,*9999)
        IF(IOSTAT.EQ.0.AND.NCHAR.GT.0) THEN
          ICOL=1
          DO 3 nd=NDAT+1,NDATA
            ISTART=0
            DO 1 ic=ICOL,NCHAR
              IF((CLINE(ic:ic).EQ.' ').OR.(CLINE(ic:ic).EQ.',')) THEN
                IF(ISTART.GT.0) THEN
                  IEND=ic-1
                  GOTO 2
                ENDIF
              ELSE
                IF(ISTART.EQ.0) THEN
                  ISTART=ic
                ENDIF
              ENDIF
    1       CONTINUE
            IEND=NCHAR
    2       IF(ISTART.GT.0) THEN
              ICOL=IEND+2
              CDATA=' '
              LENGTH=IEND-ISTART+1
              CDATA(25-LENGTH:24)=CLINE(ISTART:IEND)
              IF((CDATA.EQ.'                    EDIT').OR.
     '           (CDATA.EQ.'                    edit').OR.
     '           (CDATA.EQ.'                       E').OR.
     '           (CDATA.EQ.'                       e'))THEN
                ERROR='Edit'
                GO TO 9999
              ELSE IF((CDATA.EQ.'                 RESTART').OR.
     '                (CDATA.EQ.'                 restart').OR.
     '                (CDATA.EQ.'                       R').OR.
     '                (CDATA.EQ.'                       r'))THEN
                ERROR='Restart'
                GO TO 9999
              ELSE IF((CDATA.EQ.'                    STOP').OR.
     '                (CDATA.EQ.'                    stop').OR.
     '                (CDATA.EQ.'                       S').OR.
     '                (CDATA.EQ.'                       s'))THEN
                CLOSE(UNIT=IUNIT)
                STOP
              ELSE IF((CDATA.EQ.'                    FILE').OR.
     '                (CDATA.EQ.'                    file').OR.
     '                (CDATA.EQ.'                       F').OR.
     '                (CDATA.EQ.'                       f'))THEN
                FILEIP=.TRUE.
                WRITE(IUNIT,'($,'' Enter FILENAME: '')')
                READ(IUNIT,'(A)') FILE07
              ELSE IF((CDATA.EQ.'                     GET').OR.
     '                (CDATA.EQ.'                     get').OR.
     '                (CDATA.EQ.'                       G').OR.
     '                (CDATA.EQ.'                       g'))THEN
                IF(.NOT.FILEIP) THEN
                  FILEIP=.TRUE.
                  FILEIP_FIRST=.TRUE.
                  WRITE(IUNIT,'($,'' Enter FILENAME(.data): '')')
                  READ(IUNIT,'(A)') FILE07
                  CALL TRIM(FILE07,IBEG,IEND)
                  CALL OPENF(7,'DISK',FILE07(IBEG:IEND)//'.data','OLD',
     '              'SEQUEN','FORMATTED',132,ERROR,*9999)
                ELSE
                  FILEIP_FIRST=.FALSE.
                ENDIF
                NDAT=NDAT+1
                READ(7,*) RDATA(NDAT)
                ERROR=' '
                GO TO 9998
              ELSE IF(CDATA.EQ.'                       ?') THEN
                WRITE(UNIT=CHAR,FMT='(I5)') INFO
                FILE='FE'//CHAR(1:4)
                CALL DOCUM(FILE,'doc','Parameter number '//CHAR(5:5),
     '            ERROR,*9997)
 9997           IF(ERROR.EQ.' ') THEN
                  ERROR=' ?'
                  GO TO 9998
                ELSE IF(ERROR(1:14).EQ.'File not found') THEN
                  WRITE(OP_STRING,*) '>>Help not available'
      		  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ERROR=' ?'
                  GO TO 9998
                ELSE
                  GO TO 9999
                ENDIF
              ELSE
                IF(RVALID(CDATA)) THEN
                  NDAT=NDAT+1
                  READ(UNIT=CDATA,FMT='(G24.0)',IOSTAT=IOSTAT)
     '              RDATA(NDAT)
                  IF(IOSTAT.EQ.0) THEN
                    IF((RDATA(NDAT).GE.RMIN).AND.
     '                 (RDATA(NDAT).LE.RMAX)) THEN
                      IF(NDAT.EQ.NDATA) THEN
                        ERROR=' '
                        GO TO 9998
                      ENDIF
                    ELSE
                      ERROR=' Real data range error'
                      GO TO 9998
                    ENDIF
                  ELSE
                    ERROR=' Real data read error'
                    GO TO 9998
                  ENDIF
                ELSE
                  ERROR=' Real data type error'
                  GO TO 9998
                ENDIF
              ENDIF
            ELSE
              GOTO 4
            ENDIF
    3     CONTINUE
        ELSE IF(IOSTAT.EQ.0.AND.NCHAR.EQ.0) THEN
C         REWIND(UNIT=IUNIT)
          NDAT=NDAT+1
C         RDATA(NDAT)=RDEFLT(NDAT)
C         IF(NDAT.EQ.NDATA) THEN
          DO nd=1,NDATA
            RDATA(nd)=RDEFLT(nd)
          ENDDO
          ERROR=' '
          GO TO 9998
C         ENDIF
        ELSE
          REWIND(UNIT=IUNIT)
          ERROR=' Character data read error'
          GO TO 9998
        ENDIF
    4 CONTINUE

      ERROR=' '
 9998 CALL EXITS('RINPUT1')
      RETURN
 9999 CALL ERRORS('RINPUT1',ERROR)
      CALL EXITS('RINPUT1')
      RETURN 1
      END

                     
      SUBROUTINE SOCKET_READ_SIGNAL(ERROR,*)

C#### Subroutine: SOCKET_READ_SIGNAL
C###  Description:

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:ditr00.cmn'
      INCLUDE 'cmiss$reference:fsklib.inc'
      INCLUDE 'cmiss$reference:read00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER DATA_TYPE,ACCEPT_INT,NUMSIGNALS,i
      REAL TIME
      PARAMETER(DATA_TYPE = 3)  !Fiducial Marker information

      CALL ENTERS('SOCKET_READ_SIGNAL',*9999)

      ! Tell the front end to send data if it has it
      IF(FSKWRITE(DATA_TYPE,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999

      IF(FSKREAD(NUMSIGNALS, SK_LONG_INT, 1, CONNID2).EQ.-1) GOTO 9999

      DO i=1,NUMSIGNALS
      	IF(FSKREAD(TIME,SK_SINGLE_FLOAT,1,CONNID2).EQ.-1) GOTO 9999
        FIDMARK(i)=DBLE(TIME)
      
      	IF(FSKREAD(ACCEPT_INT,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999
        IF(ACCEPT_INT.EQ.1) THEN
          ACCEPT(i)=.TRUE.
        ELSE
          ACCEPT(i)=.FALSE.
        ENDIF
      ENDDO  
      
      CALL EXITS('SOCKET_READ_SIGNAL')
      RETURN
 9999 CALL ERRORS('SOCKET_READ_SIGNAL',ERROR)
      CALL EXITS('SOCKET_READ_SIGNAL')
      RETURN 1
      END


      SUBROUTINE WRITE_POLYLINE(IBUNDLE,IW,NPOINTS,X,Y,ERROR,*)

C#### Subroutine: WRITE_POLYLINE
C###  Description:

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:fsklib.inc'
!     Parameter List
      INTEGER IBUNDLE,IW,NPOINTS
      REAL X(NPOINTS),Y(NPOINTS)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER DATA_TYPE
      PARAMETER(DATA_TYPE = 4)  !Polyline Data

      CALL ENTERS('WRITE_POLYLINE',*9999)

      IF(FSKWRITE(DATA_TYPE,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999

      IF(FSKWRITE(IW, SK_LONG_INT, 1, CONNID2).EQ.-1) GOTO 9999
      IF(FSKWRITE(IBUNDLE, SK_LONG_INT, 1, CONNID2).EQ.-1) GOTO 9999
      IF(FSKWRITE(NPOINTS, SK_LONG_INT, 1, CONNID2).EQ.-1) GOTO 9999

      IF(FSKWRITE(X, SK_SINGLE_FLOAT, NPOINTS, CONNID2).EQ.-1) GOTO 9999
      IF(FSKWRITE(Y, SK_SINGLE_FLOAT, NPOINTS, CONNID2).EQ.-1) GOTO 9999

      CALL EXITS('WRITE_POLYLINE')
      RETURN
 9999 CALL ERRORS('WRITE_POLYLINE',ERROR)
      CALL EXITS('WRITE_POLYLINE')
      RETURN 1
      END
      	          

C LKC 6-OCT-1999 archived with fiducial

      SUBROUTINE WRITE_SIGNAL(SIGNALDATA,NUMSIGNALS,NUMSAMPLES,
     '    NUMSAMPLESPERSECOND,ERROR,*)

C#### Subroutine: WRITE_SIGNAL
C###  Description:

C GMH 20-5-96 Changing SIGNALDATA from int*2 to int

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:fsklib.inc'
      INCLUDE 'cmiss$reference:read00.cmn'
!     Parameter List
      INTEGER NUMSIGNALS,NUMSAMPLES,NUMSAMPLESPERSECOND,
     '  SIGNALDATA(NUMSIGNALS,NUMSAMPLES)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CLEN,INTSTR(1024),i,j,
     '  SIGNAL(2048)
      INTEGER DATA_TYPE,ACCEPT_INT
      CHARACTER HEADER*5
      PARAMETER(DATA_TYPE = 2)  !Signal Data

      CALL ENTERS('WRITE_SIGNAL',*9999)

      IF(FSKWRITE(DATA_TYPE,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999

      IF(FSKWRITE(NUMSIGNALS, SK_LONG_INT, 1, CONNID2).EQ.-1) GOTO 9999
      IF(FSKWRITE(NUMSAMPLES, SK_LONG_INT, 1, CONNID2).EQ.-1) GOTO 9999
      IF(FSKWRITE(NUMSAMPLESPERSECOND, SK_LONG_INT, 1,
     '   CONNID2).EQ.-1) GOTO 9999

      DO i=1,NUMSIGNALS
        DO j=1,NUMSAMPLES
          SIGNAL(j)=SIGNALDATA(i,j)
        ENDDO
        IF(FSKWRITE(SIGNAL, SK_LONG_INT, NUMSAMPLES, CONNID2)
     '    .EQ.-1) GOTO 9999
        HEADER=SIGHEADER(1,i)
        CLEN=FSKLEN(HEADER)
        CALL FSKF2C(HEADER,CLEN,INTSTR)
        IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999
        IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1) GOTO 9999

        IF(ACCEPT(i)) THEN
          ACCEPT_INT=1
        ELSE
          ACCEPT_INT=0
        ENDIF
        IF(FSKWRITE(ACCEPT_INT,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999
      ENDDO

      CALL EXITS('WRITE_SIGNAL')
      RETURN
 9999 CALL ERRORS('WRITE_SIGNAL',ERROR)
      CALL EXITS('WRITE_SIGNAL')
      RETURN 1
      END


Module FE01
=========== 
        
C LC removed section from routine:
C#### Subroutine: PARSE_ELEMENTS
C###  Description:
C###    PARSE_ELEMENTS parses command CO for keyword 'elements' 
C###    and returns list of elements in 
C###    NELIST(nolist),nolist=1,NELIST(0).
C**** Keyword 'elements' can be followed by a list or by a group name.
C**** Default is all elements. Max number in list is NE_R_M.
C
C old MPN 4-Apr-96: now uses PARSILG
C        CALL TRIM(CO(N3CO+1),IBEG1,IEND1)
C        IPOS1=IBEG1
C        NELIST(0)=0 !current total number of nodes in list
C        DO WHILE(IPOS1.LE.IEND1)
C          IF(CO(N3CO+1)(IPOS1:IPOS1).GE.'0'.AND.
C     '      CO(N3CO+1)(IPOS1:IPOS1).LE.'9') THEN !list of #'s follows
C            IPOS2=IPOS1
C            CHARAC=CO(N3CO+1)(IPOS2:IPOS2)
C            DO WHILE(((CHARAC.GE.'0'.AND.CHARAC.LE.'9').OR. !a digit
C     '        CHARAC.EQ.'.'.OR.CHARAC.EQ.',').AND. !a period or comma
C     '        IPOS2.LT.IEND1)
C              IPOS2=IPOS2+1
C              CHARAC=CO(N3CO+1)(IPOS2:IPOS2)
C            ENDDO
C            IF(IPOS2.NE.IEND1) IPOS2=IPOS2-2
CC           Parse and add current list of nodes to list
C            CALL PARSIL(CO(N3CO+1)(IPOS1:IPOS2),NE_R_M-NELIST(0),
C     '        NUMINLIST,NELIST(NELIST(0)+1),ERROR,*9999)
C            NELIST(0)=NELIST(0)+NUMINLIST
C          ELSE !possibly a group name entered
C            IPOS2=IPOS1+CLOCAT(',',CO(N3CO+1)(IPOS1:IEND1))
C            IF(IPOS2.EQ.IPOS1) THEN
C              IPOS2=IEND1
C            ELSE
C              IPOS2=IPOS2-2
C            ENDIF
C            CHAR=CUPPER(CO(N3CO+1)(IPOS1:IPOS2)) !group name
C            CALL TRIM(CHAR,IBEG2,IEND2)
C            N1GREL=0
C            DO nogrel=1,NTGREL
C              LABEL=CUPPER(LAGREL(nogrel))
C              CALL TRIM(LABEL,IBEG3,IEND3)
C              IF(CHAR(IBEG2:IEND2).EQ.LABEL(IBEG3:IEND3)) THEN
C                N1GREL=nogrel
C                GO TO 100
C              ENDIF
C            ENDDO
C 100        IF(N1GREL.GT.0) THEN
C              DO nolist=1,LIGREL(0,N1GREL)
C                NELIST(NELIST(0)+nolist)=LIGREL(nolist,N1GREL)
C              ENDDO !nolist
C              NELIST(0)=NELIST(0)+LIGREL(0,N1GREL)
C            ELSE
C              ERROR=' >>Group name '//CHAR(IBEG2:IEND2)//' not defined'
C              GO TO 9999
C            ENDIF
C          ENDIF
C          IPOS1=IPOS2+2
C        ENDDO !while(IPOS1.LE.IEND)



C LC removed section from routine :
C#### Subroutine: PARSE_GRID
C###  Description:
C###    PARSE_GRID parses command CO for keyword 'grid' and 
C###    returns list of grid points in 
C###    NQLIST(nolist),nolist=1,NQLIST(0).
C**** Keyword 'grid' can be followed by a list or by a group name.
C**** Default is all grid points. Max number in list is NQM.
C
C old MPN 4-Apr-96: now uses PARSILG
C        CALL TRIM(CO(N3CO+1),IBEG1,IEND1)
C        IPOS1=IBEG1
C        NQLIST(0)=0 !current total number of grid points in list
C        DO WHILE(IPOS1.LE.IEND1)
C          IF(CO(N3CO+1)(IPOS1:IPOS1).GE.'0'.AND.
C     '      CO(N3CO+1)(IPOS1:IPOS1).LE.'9') THEN !list of #'s follows
C            IPOS2=IPOS1
C            CHARAC=CO(N3CO+1)(IPOS2:IPOS2)
C            DO WHILE(((CHARAC.GE.'0'.AND.CHARAC.LE.'9').OR. !a digit
C     '        CHARAC.EQ.'.'.OR.CHARAC.EQ.',').AND. !a period or comma
C     '        IPOS2.LT.IEND1)
C              IPOS2=IPOS2+1
C              CHARAC=CO(N3CO+1)(IPOS2:IPOS2)
C            ENDDO
C            IF(IPOS2.NE.IEND1) IPOS2=IPOS2-2
CC           Parse and add current list of grid points to list
C            CALL PARSIL(CO(N3CO+1)(IPOS1:IPOS2),NQM-NQLIST(0),
C     '        NUMINLIST,NQLIST(NQLIST(0)+1),ERROR,*9999)
C            NQLIST(0)=NQLIST(0)+NUMINLIST
C          ELSE !possibly a group name entered
C            IPOS2=IPOS1+CLOCAT(',',CO(N3CO+1)(IPOS1:IEND1))
C            IF(IPOS2.EQ.IPOS1) THEN
C              IPOS2=IEND1
C            ELSE
C              IPOS2=IPOS2-2
C            ENDIF
C            CHAR=CUPPER(CO(N3CO+1)(IPOS1:IPOS2)) !group name
C            CALL TRIM(CHAR,IBEG2,IEND2)
C            N1GRGR=0
C            DO nogrgr=1,NTGRGR
C              LABEL=CUPPER(LAGRGR(nogrgr))
C              CALL TRIM(LABEL,IBEG3,IEND3)
C              IF(CHAR(IBEG2:IEND2).EQ.LABEL(IBEG3:IEND3)) THEN
C                N1GRGR=nogrgr
C                GO TO 100
C              ENDIF
C            ENDDO
C 100        IF(N1GRGR.GT.0) THEN
C              DO nolist=1,LIGRGR(0,N1GRGR)
C                NQLIST(NQLIST(0)+nolist)=LIGRGR(nolist,N1GRGR)
C              ENDDO !nolist
C              NQLIST(0)=NQLIST(0)+LIGRGR(0,N1GRGR)
C            ELSE
C              ERROR=' >>Group name '//CHAR(IBEG2:IEND2)//' not defined'
C              GO TO 9999
C            ENDIF
C          ENDIF
C          IPOS1=IPOS2+2
C        ENDDO !while(IPOS1.LE.IEND)



C LC removed section from routine :
C
C#### Subroutine: PARSE_NODES
C###  Description:
C###    PARSE_NODES parses command CO for keyword 'nodes' and 
C###    returns list of nodes in NPLIST(nolist),nolist=1,NPLIST(0).
C**** Keyword 'nodes' can be followed by a list or by a group name.
C**** Default is all nodes. Max number in list is NP_R_M.

C old MPN 4-Apr-96: now uses PARSILG
C        CALL TRIM(CO(N3CO+1),IBEG1,IEND1)
C        IPOS1=IBEG1
C        NPLIST(0)=0 !current total number of nodes in list
C        DO WHILE(IPOS1.LE.IEND1)
C          IF(CO(N3CO+1)(IPOS1:IPOS1).GE.'0'.AND.
C     '      CO(N3CO+1)(IPOS1:IPOS1).LE.'9') THEN !list of #'s follows
C            IPOS2=IPOS1
C            CHARAC=CO(N3CO+1)(IPOS2:IPOS2)
C            DO WHILE(((CHARAC.GE.'0'.AND.CHARAC.LE.'9').OR. !a digit
C     '        CHARAC.EQ.'.'.OR.CHARAC.EQ.',').AND. !a period or comma
C     '        IPOS2.LT.IEND1)
C              IPOS2=IPOS2+1
C              CHARAC=CO(N3CO+1)(IPOS2:IPOS2)
C            ENDDO
C            IF(IPOS2.NE.IEND1) IPOS2=IPOS2-2
CC           Parse and add current list of nodes to list
C            CALL PARSIL(CO(N3CO+1)(IPOS1:IPOS2),NP_R_M-NPLIST(0),
C     '        NUMINLIST,NPLIST(NPLIST(0)+1),ERROR,*9999)
C            NPLIST(0)=NPLIST(0)+NUMINLIST
C          ELSE !possibly a group name entered
C            IPOS2=IPOS1+CLOCAT(',',CO(N3CO+1)(IPOS1:IEND1))
C            IF(IPOS2.EQ.IPOS1) THEN
C              IPOS2=IEND1
C            ELSE
C              IPOS2=IPOS2-2
C            ENDIF
C            CHAR=CUPPER(CO(N3CO+1)(IPOS1:IPOS2)) !group name
C            CALL TRIM(CHAR,IBEG2,IEND2)
C            N1GRNO=0
C            DO nogrno=1,NTGRNO
C              LABEL=CUPPER(LAGRNO(nogrno))
C              CALL TRIM(LABEL,IBEG3,IEND3)
C              IF(CHAR(IBEG2:IEND2).EQ.LABEL(IBEG3:IEND3)) THEN
C                N1GRNO=nogrno
C                GO TO 100
C              ENDIF
C            ENDDO
C 100        IF(N1GRNO.GT.0) THEN
C              DO nolist=1,LIGRNO(0,N1GRNO)
C                NPLIST(NPLIST(0)+nolist)=LIGRNO(nolist,N1GRNO)
C              ENDDO !nolist
C              NPLIST(0)=NPLIST(0)+LIGRNO(0,N1GRNO)
C            ELSE
C              ERROR=' >>Group name '//CHAR(IBEG2:IEND2)//' not defined'
C              GO TO 9999
C            ENDIF
C          ENDIF
C          IPOS1=IPOS2+2
C        ENDDO !while(IPOS1.LE.IEND)



C 24/2/97 LC removed section from routine :
C#### Subroutine: PARSILG
C###  Description:
C###    PARSILG handles input of a sequence of integer data including 
C###    group names in the character string CLINE.   
C###    Set CDATA to 'NODES','ELEMENTS','FACES','GRIDS' or 'LINES' 
C###    to pick up appropriate group name. 
C###    Values are output in ILIST(n),n=1,ILIST(0).
C***    Rewritten MPN 4-4-96
C
C old : rewritten MPN 4-4-96
C      CALL TRIM(CLINE,IBEG,NCHAR)
CC     Check for null string
C      IF(NCHAR.EQ.0) THEN !no data
C        ILIST(0)=0
C        ILIST(1)=0
C        GOTO 4
C      ENDIF
CC     Deal with each data value in turn
C      nValues=0 !keeps track of #values read in
C      ICOL=1
CC AJP 3/11/95  DO nd=nValues+1,99 !deal with remaining data values
C      DO nd=nValues+1,999 !deal with remaining data values
C        IF(ICOL.GT.NCHAR) GOTO 4
C        ISTART=0
C        DO ic=ICOL,NCHAR !examine characters left in CLINE
C          IF((CLINE(ic:ic).EQ.' ').OR.(CLINE(ic:ic).EQ.',')) THEN
C            IF(ISTART.GT.0) THEN
C              IEND=ic-1
C              GOTO 2
C            ENDIF
C          ELSE !character is not blank or ',' 
C            IF(ISTART.EQ.0) THEN
C              ISTART=ic
C            ENDIF
C          ENDIF
C        ENDDO !ic
C        IEND=NCHAR
C    2   IF(ISTART.GT.0) THEN
C          IF(DOP) THEN
C            WRITE(OP_STRING,'('' istart='',I4,'' iend='',I4)')
C     '        ISTART,IEND
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
C          ICOL=IEND+2
C          CHDATA=' '
C          LENGTH=IEND-ISTART+1
C          ISTART1=MAX(17-LENGTH,1)
C          IEND1=MIN(IEND,ISTART+15)
C          CHDATA(ISTART1:16)=CLINE(ISTART:IEND1) !I16 format
C          IF(IVALID(CHDATA)) THEN !valid integer data
C            nValues=nValues+1
C            READ(UNIT=CHDATA,FMT='(I16)',IOSTAT=IOSTAT) ILIST(nValues)
C            IF(IOSTAT.EQ.0) THEN
C              ILIST(0)=nValues !is #integers read in
C              IF(DOP) THEN
C                WRITE(OP_STRING,'('' ILIST(0)='',I3,'
C     '            //''' ILIST('',I2,'')='',I5)') 
C     '            ILIST(0),nValues,ILIST(nValues)
C                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C              ENDIF
C            ELSE
C              ERROR=' Integer data read error'
C              GO TO 9999
C            ENDIF
C          
CC         check for list in form FIRST..LAST:INCR 
C          ELSE IF(CLOCAT('..',CHDATA(1:16)).GT.0) THEN
C            IF(DOP) THEN
C              WRITE(OP_STRING,'('' CHDATA(1:16) is'',A)') CHDATA(1:16)
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C            ENDIF !dop
CC!!! Passing 999-nValues is very dangerous: should pass max dimension
CC!!! of ILIST (ILIST_MAX) from calling routing and use ILIST_MAX-nValues
C            CALL PARSIL(CHDATA(1:16),999-nValues,NTIL,ILIST(nValues+1),
C     '        ERROR,*9999)  
C            nValues=nValues+NTIL
C            ILIST(0)=nValues
C            IF(DOP) THEN
C              WRITE(OP_STRING,'('' NTIL='',I4,'' nValues='',I4)') 
C     '          NTIL,nValues
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C              WRITE(OP_STRING,'('' ILIST:'',10I5,/:(10I5))') 
C     '          (ILIST(n),n=1,ILIST(0))
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C            ENDIF !dop
C          
CC         check for group name
C          ELSE 
C            IF(DOP) THEN
C              WRITE(OP_STRING,'('' CLINE is '',A)') CLINE(istart:iend)
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C              WRITE(OP_STRING,'('' CDATA='',A)') CDATA(1:10)
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C            ENDIF !dop
C            FOUND_GROUP=.FALSE.
C            IF(CDATA(1:5).EQ.'NODES') THEN
CC             put group nodes into ILIST(n)
C              DO nogrno=1,NTGRNO
C                IF(CLINE(istart:iend).EQ.
C     '            LAGRNO(nogrno)(1:iend-istart+1)) THEN
C                  IF(DOP) THEN
C                    WRITE(OP_STRING,'('' Group name is '',A)')
C     '                LAGRNO(nogrno)
C                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                    WRITE(OP_STRING,'('' Group nodes are: '',10I3)') 
C     '                (LIGRNO(n,nogrno),n=1,LIGRNO(0,nogrno))
C                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                  ENDIF !dop
C                  DO n=1,LIGRNO(0,nogrno)
C                    ILIST(nValues+n)=LIGRNO(n,nogrno)
CC                   Flag warning if nodes are repeated in the list
C                    IF(INLIST(ILIST(nValues+n),ILIST(1),
C     '                nValues+n-1,N1)) THEN
C                      WRITE(OP_STRING,'('' >>WARNING: Node '',I5,'
C     '                  //''' is repeated in the list'')') 
C     '                  ILIST(nValues+n)
C                      CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C                    ENDIF
C                  ENDDO !n
C                  nValues=nValues+LIGRNO(0,nogrno)
C                  ILIST(0)=nValues
C                  FOUND_GROUP=.TRUE.
C                ENDIF !cline
C              ENDDO !nogrno
C
C            ELSE IF(CDATA(1:8).EQ.'ELEMENTS') THEN
CC             put group elements into ILIST(n)
C              DO nogrel=1,NTGREL
C                IF(CLINE(istart:iend).EQ.
C     '            LAGREL(nogrel)(1:iend-istart+1)) THEN
C                  IF(DOP) THEN
C                    WRITE(OP_STRING,'('' Group name is '',A)')
C     '                LAGREL(nogrel)
C                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                    WRITE(OP_STRING,'('' Group elements are: '',10I3)') 
C     '                (LIGREL(n,nogrel),n=1,LIGREL(0,nogrel))
C                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                  ENDIF !dop
C                  DO n=1,LIGREL(0,nogrel)
C                    ILIST(nValues+n)=LIGREL(n,nogrel)
CC                   Flag warning if elementss are repeated in the list
C                    IF(INLIST(ILIST(nValues+n),ILIST(1),
C     '                nValues+n-1,N1)) THEN
C                      WRITE(OP_STRING,'('' >>WARNING: Element '',I5,'
C     '                  //''' is repeated in the list'')') 
C     '                  ILIST(nValues+n)
C                      CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C                    ENDIF
C                  ENDDO
C                  nValues=nValues+LIGREL(0,nogrel)
C                  ILIST(0)=nValues
C                  FOUND_GROUP=.TRUE.
C                ENDIF !cline
C              ENDDO !nogrel
C
C            ELSE IF(CDATA(1:5).EQ.'FACES') THEN
CC             put group faces into ILIST(n)
C              DO nogrfa=1,NTGRFA
C                IF(CLINE(istart:iend).EQ.
C     '            LAGRFA(nogrfa)(1:iend-istart+1)) THEN
C                  IF(DOP) THEN
C                    WRITE(OP_STRING,'('' Group name is '',A)')
C     '                LAGRFA(nogrfa)
C                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                    WRITE(OP_STRING,'('' Group faces are: '',10I3)') 
C     '                (LIGRFA(n,nogrfa),n=1,LIGRFA(0,nogrfa))
C                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                  ENDIF !dop
C                  DO n=1,LIGRFA(0,nogrfa)
C                    ILIST(nValues+n)=LIGRFA(n,nogrfa)
CC                   Flag warning if faces are repeated in the list
C                    IF(INLIST(ILIST(nValues+n),ILIST(1),
C     '                nValues+n-1,N1)) THEN
C                      WRITE(OP_STRING,'('' >>WARNING: Face '',I5,'
C     '                  //''' is repeated in the list'')') 
C     '                  ILIST(nValues+n)
C                      CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C                    ENDIF
C                  ENDDO
C                  nValues=nValues+LIGRFA(0,nogrfa)
C                  ILIST(0)=nValues
C                  FOUND_GROUP=.TRUE.
C                ENDIF !cline
C              ENDDO !nogrfa
C
C            ELSE IF(CDATA(1:5).EQ.'GRIDS') THEN
CC             put group grid points into ILIST(n)
C              DO nogrgr=1,NTGRGR
C                IF(CLINE(istart:iend).EQ.
C     '            LAGRGR(nogrgr)(1:iend-istart+1)) THEN
C                  IF(DOP) THEN
C                    WRITE(OP_STRING,'('' Group name is '',A)')
C     '                LAGRGR(nogrgr)
C                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                    WRITE(OP_STRING,'('' Group grid pts are: '',10I3)') 
C     '                (LIGRGR(n,nogrgr),n=1,LIGRGR(0,nogrgr))
C                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                  ENDIF !dop
C                  DO n=1,LIGRGR(0,nogrgr)
C                    ILIST(nValues+n)=LIGRGR(n,nogrgr)
CC                   Flag warning if grid pts are repeated in the list
C                    IF(INLIST(ILIST(nValues+n),ILIST(1),
C     '                nValues+n-1,N1)) THEN
C                      WRITE(OP_STRING,'('' >>WARNING: Grid pt '',I5,'
C     '                  //''' is repeated in the list'')') 
C     '                  ILIST(nValues+n)
C                      CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C                    ENDIF
C                  ENDDO
C                  nValues=nValues+LIGRGR(0,nogrgr)
C                  ILIST(0)=nValues
C                  FOUND_GROUP=.TRUE.
C                ENDIF !cline
C              ENDDO !nogrgr
C
C            ELSE IF(CDATA(1:5).EQ.'LINES') THEN
CC             put group lines into ILIST(n)
C              DO nogrli=1,NTGRLI
C                IF(CLINE(istart:iend).EQ.
C     '            LAGRLI(nogrli)(1:iend-istart+1)) THEN
C                  IF(DOP) THEN
C                    WRITE(OP_STRING,'('' Group name is '',A)')
C     '                LAGRLI(nogrli)
C                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                    WRITE(OP_STRING,'('' Group lines are: '',10I3)') 
C     '                (LIGRLI(n,nogrli),n=1,LIGRLI(0,nogrli))
C                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                  ENDIF !dop
C                  DO n=1,LIGRLI(0,nogrli)
C                    ILIST(nValues+n)=LIGRLI(n,nogrli)
CC                   Flag warning if lines are repeated in the list
C                    IF(INLIST(ILIST(nValues+n),ILIST(1),
C     '                nValues+n-1,N1)) THEN
C                      WRITE(OP_STRING,'('' >>WARNING: Line '',I5,'
C     '                  //''' is repeated in the list'')') 
C     '                  ILIST(nValues+n)
C                      CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C                    ENDIF
C                  ENDDO
C                  nValues=nValues+LIGRLI(0,nogrli)
C                  ILIST(0)=nValues
C                  FOUND_GROUP=.TRUE.
C                ENDIF !cline
C              ENDDO !nogrli
C
C            ENDIF !cdata
C
C            IF(.NOT.FOUND_GROUP) THEN
C              ERROR=' Group "'//CLINE(istart:iend)//'" not defined'
C              GO TO 9999
C            ENDIF
C          ENDIF !valid integer data or group name
C        ELSE  !istart=0
C          GOTO 4
C        ENDIF !istart>0
C      ENDDO !nd (i.e. index for list of data values in CLINE)

      

      LOGICAL FUNCTION AFFIRM(STRING)

C#### Function: AFFIRM
C###  Type: LOGICAL
C###  Description:
C###    AFFIRM returns .TRUE. if the string STRING is an 
C###    affirmative answer.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      CHARACTER STRING*(*)
!     Local Variables
      CHARACTER CUPPER*(MXCH)
      LOGICAL ABBREV

      STRING=CUPPER(STRING)
      IF((ABBREV(STRING,'YES',1))
     '  .OR.(ABBREV(STRING,'TRUE',1))) THEN
        AFFIRM=.TRUE.
      ELSE
        AFFIRM=.FALSE.
      ENDIF

      RETURN
      END

      INTEGER FUNCTION BINSEARCH(ITEM,LIST,N,FOUND)

C#### Function: BINSEARCH
C###  Type: INTEGER
C###  Description: 
C###    BINSEARCH perfoms a binary search on an ordered list
C###    of integers of size N. The function returns the position of the
C###    desired item in the list and sets FOUND to .TRUE. if it is 
C###    found otherwise the function returns the position before the
C###    location in the ordered list where the iterm should go and
C###    sets FOUND to .FALSE..

      IMPLICIT NONE
!     Parameter List
      INTEGER ITEM,LIST(*),N
      LOGICAL FOUND
!     Local Variables
      INTEGER LOWLIMIT,MIDPOINT,UPLIMIT


      LOWLIMIT=0
      IF(N.GT.0) THEN
        UPLIMIT=N+1
        DO WHILE((UPLIMIT-LOWLIMIT).GT.1)
          MIDPOINT=(UPLIMIT+LOWLIMIT)/2
          IF(LIST(MIDPOINT).GT.ITEM) THEN
            UPLIMIT=MIDPOINT
          ELSE
            LOWLIMIT=MIDPOINT       
          ENDIF
        ENDDO
        IF(LIST(LOWLIMIT).EQ.ITEM) THEN
          FOUND=.TRUE.
        ELSE
          FOUND=.FALSE.
        ENDIF
      ELSE
        FOUND=.FALSE.      
      ENDIF

      BINSEARCH=LOWLIMIT

      RETURN
      END
                     
      CHARACTER*(*) FUNCTION CFROML(LDATA,FORMAT)

C#### Function: CFROML
C###  Type: CHARACTER
C###  Description:
C###    Converts logical variable LDATA to character string CFROML
C###    as specified by the character string FORMAT.

      IMPLICIT NONE
!     Parameter List
      CHARACTER FORMAT*(*)
      LOGICAL LDATA
!     Local Variables
      CHARACTER CDATA*500

      WRITE(UNIT=CDATA,FMT=FORMAT) LDATA
      CFROML=CDATA

      RETURN
      END
           
                   
!news AAY 21 Feb 92 This function needed in GETSTR due to f77 compiler
      CHARACTER*(*) FUNCTION CONCAT(STR1,IBEG1,IEND1,STR2,IBEG2,IEND2)

C**** concatenation via function
C     Workaround for f77 compiler bug which does not allow concatenation
C     of dynamically dimensioned strings in function calls

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER IBEG1,IEND1,IBEG2,IEND2
      CHARACTER*(*) STR1,STR2
!     Local declarations
      INTEGER LEN
      CHARACTER ERROR*10
                                         
      IF(IEND1-IBEG1+1+IEND2-IBEG2+1.GT.LEN(CONCAT))THEN
        WRITE(OP_STRING,*)' Error from CONCAT- resultant too large'
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ELSE
        CONCAT=STR1(IBEG1:IEND1)//STR2(IBEG2:IEND2)
      ENDIF

 9999 RETURN
      END


      CHARACTER*1 FUNCTION CI1(i)

C#### Function: CI1
C###  Type: CHARACTER
C###  Description: 
C###    CI1 returns number of decimal places required for I for 
C###    formatting purposes.

      IMPLICIT NONE
!     Parameter List
      INTEGER i
!     Local Variables
      INTEGER ib,ie
      CHARACTER CFROMI*9,CI9*9

      WRITE(UNIT=CI9,FMT='(I9)') i !turns i into character string CI9
      CALL TRIM(CI9,ib,ie)         !returns beginning and end of CI9
      WRITE(UNIT=CI1,FMT='(I1)') ie-ib+1 !#decimals required for i

      RETURN
      END


      LOGICAL FUNCTION EXIST(FILE,ERROR)

C#### Function: EXIST
C###  Type: LOGICAL
C###  Description:
C###    EXIST returns .TRUE. if a file corresponding to the character 
C###    string FILE exists. FILE is defined by Filename Filetype 
C###    [Filemode] where the default Filemode A.

      IMPLICIT NONE
!     Parameter List
      CHARACTER FILE*100,ERROR*(*)
!     Local Variables
      INTEGER IERROR

C     CALL ENTERS('EXIST',*9999)
      CALL CMS('STATE '//FILE,IERROR)
      IF(IERROR.EQ.0) THEN
        EXIST=.TRUE.
        ERROR=' '
      ELSE IF(IERROR.EQ.28) THEN
        EXIST=.FALSE.
        ERROR=' '
      ELSE
        ERROR=' Error occurred in CMS(STATE('//FILE//')'
      ENDIF
C     CALL EXITS('EXIST')
      RETURN
      END


      INTEGER FUNCTION MAXD(NTLIST,A)

C#### Function: MAXD
C###  Type: INTEGER
C###  Description:
C###    MAXD returns integer position of absolute maximum of 
C###    A(i),i=1,NTLIST.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER NTLIST
      REAL*8 A(NTLIST)
!     Local Variables
      INTEGER nolist
      REAL*8 AMAX
      CHARACTER ERROR*10

      AMAX=DABS(A(1))
      MAXD=1
      DO nolist=2,NTLIST
        IF(DABS(A(nolist)).GT.AMAX) THEN
          AMAX=DABS(A(nolist))
          MAXD=nolist
        ENDIF
      ENDDO
      WRITE(OP_STRING,'('' ntlist='',I4,'' maxd='',I4,'' amax='','
     '  //'E10.4)') NTLIST,MAXD,A(MAXD)
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

 9999 RETURN
      END
          
      INTEGER FUNCTION MAXIM(NTLIST,A)

C#### Function: MAXIM
C###  Type: INTEGER
C###  Description:
C###    MAXIM returns integer position of absolute maximum of 
C###    A(i),i=1,NTLIST.

      IMPLICIT NONE
!     Parameter List
      INTEGER NTLIST
      REAL*8 A(NTLIST)
!     Local Variables
      INTEGER nolist
      REAL*8 AMAX

      AMAX=DABS(A(1))
      MAXIM=1
      DO nolist=2,NTLIST
        IF(DABS(A(nolist)).GT.AMAX) THEN
          AMAX=DABS(A(nolist))
          MAXIM=nolist
        ENDIF
      ENDDO

      RETURN
      END
      

                     
      INTEGER FUNCTION MAXIM2(NTLIST,A)

C#### Function: MAXIM2
C###  Type: INTEGER
C###  Description: 
C###    MAXIM2 returns integer position of maximum of A(i),i=1,NTLIST.

      IMPLICIT NONE
!     Parameter List
      INTEGER NTLIST
      REAL*8 A(NTLIST)
!     Local Variables
      INTEGER nolist
      REAL*8 AMAX

      AMAX=A(1)
      MAXIM2=1
      DO nolist=2,NTLIST
        IF(A(nolist).GT.AMAX) THEN
          AMAX=A(nolist)
          MAXIM2=nolist
        ENDIF
      ENDDO

      RETURN
      END



      
      INTEGER FUNCTION MINIM(NTLIST,A)

C#### Function: MINIM
C###  Type: INTEGER
C###  Description: 
C###    MINIM returns integer position of absolute minimum of 
C###    A(i),i=1,NTLIST.

      IMPLICIT NONE
!     Parameter List
      INTEGER NTLIST
      REAL*8 A(NTLIST)
!     Local Variables
      INTEGER nolist
      REAL*8 AMIN          

      AMIN=DABS(A(1))
      MINIM=1
      DO nolist=2,NTLIST
        IF(DABS(A(nolist)).LT.AMIN) THEN
          AMIN=DABS(A(nolist))
          MINIM=nolist
        ENDIF
      ENDDO

      RETURN
      END




      INTEGER FUNCTION MINIM2(NTLIST,A)

C#### Function: MINIM2
C###  Type: INTEGER
C###  Description: 
C###    MINIM2 returns integer position of minimum of A(i),i=1,NTLIST.

      IMPLICIT NONE
!     Parameter List
      INTEGER NTLIST
      REAL*8 A(NTLIST)
!     Local Variables
      INTEGER nolist
      REAL*8 AMIN

      AMIN=A(1)
      MINIM2=1
      DO nolist=2,NTLIST
        IF(A(nolist).LT.AMIN) THEN
          AMIN=A(nolist)
          MINIM2=nolist
        ENDIF
      ENDDO

      RETURN
      END


       
      LOGICAL FUNCTION NEGATE(STRING)

C#### Function: NEGATE
C###  Type: LOGICAL
C###  Description: 
C###    NEGATE returns .TRUE. if the string STRING is a negative answer.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      CHARACTER STRING*(*)
!     Local Variables
      CHARACTER CUPPER*(MXCH)
      LOGICAL ABBREV

C     CALL ENTERS('NEGATE',*9999)
      STRING=CUPPER(STRING)
      IF((ABBREV(STRING,'no',1))
     '  .OR.(ABBREV(STRING,'FALSE',1))) THEN
        NEGATE=.TRUE.
      ELSE
        NEGATE=.FALSE.
      ENDIF

C     CALL EXITS('NEGATE')
      RETURN
      END


      
     
      INTEGER FUNCTION NKNU(IDO,nb,nu)

C#### Function: NKNU
C###  Type: INTEGER
C###  Description: 
C###    NKNU returns derivative index nk given partial derivative 
C###    index nu.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IDO(NKM),nb,nu
!     Local Variables
      CHARACTER ERROR*10
      INTEGER nk

      DO nk=1,NKT(0,nb)
        IF(IDO(nk).eq.nu) THEN
          NKNU=nk
          RETURN
        ENDIF
      ENDDO
      WRITE(OP_STRING,*) ' >>NKNU=0'
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      NKNU=0

 9999 RETURN
      END
           

      REAL*8 FUNCTION QUAD(TINCR,A)

C#### Function: QUAD
C###  Type: REAL*8
C###  Description: 
C###    QUAD performs quadrature.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      REAL*8 A(*),TINCR
!     Local Variables
      INTEGER nd

      QUAD=A(1)/2.0D0
      DO nd=2,NDT-1
        QUAD=QUAD+A(nd)
      ENDDO
      QUAD=QUAD+A(NDT)/2.0D0
      QUAD=QUAD*TINCR

      RETURN
      END


C LKC 6-OCT-1999 Archived

      SUBROUTINE CALC_FID(ISAMPLES,ISIGNAL,NAV,FIDMARK,TTSTART,TTEND,
     '  ACCEPT,ERROR,*)

C#### Subroutine: CALC_FID
C###  Description:
C###    CALC_FID calculates position of fiducial markers.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:pics00.cmn'
!     Parameter List
      INTEGER ISAMPLES,ISIGNAL,NAV
      REAL*8 FIDMARK,TTEND,TTSTART
      CHARACTER ERROR*(*)
      LOGICAL ACCEPT
!     Local Variables
      INTEGER i,IEND,ISTART,j,MARK
      REAL*8 AV1,AV2,SAMPLES,SLOPE,STEEPEST

      CALL ENTERS('CALC_FID',*9999)

      SAMPLES=DBLE(ISAMPLES)
      ISTART=INT(TTSTART*SAMPLES)
      IEND=INT(TTEND*SAMPLES)
      STEEPEST=0.0D0
      MARK=0

      DO i=ISTART+NAV+1,IEND-NAV-1
        AV1=0.0D0
        AV2=0.0D0
        DO j=i-NAV,i-1
          AV1=AV1+I2P(ISIGNAL,j,1)
        ENDDO
        DO j=i+1,i+NAV
          AV2=AV2+I2P(ISIGNAL,j,1)
        ENDDO
        SLOPE=DABS(AV2-AV1)
        IF(SLOPE.GT.STEEPEST) THEN
          MARK=i
          STEEPEST=SLOPE
        ENDIF
      ENDDO
      IF(MARK.EQ.0) THEN
        ACCEPT=.FALSE.
      ELSE
        FIDMARK=DBLE(MARK)/SAMPLES
      ENDIF

      CALL EXITS('CALC_FID')
      RETURN
 9999 CALL ERRORS('CALC_FID',ERROR)
      CALL EXITS('CALC_FID')
      RETURN 1
      END

      
      SUBROUTINE CMS(STRING,IRET)

C#### Subroutine: CMS
C###  Description:
C###    CMS replaces CMS routine when using VAX machines.

      IMPLICIT NONE
!     Parameter List
      INTEGER IRET
      CHARACTER STRING*(*)

      IRET=0
      IF(STRING(1:7).EQ.'CLRSCRN') THEN
        RETURN
      ELSE IF(STRING(1:7).EQ.'FILEDEF') THEN
        RETURN
      ELSE IF(STRING(1:5).EQ.'STATE') THEN
        RETURN
      ELSE IF(STRING(1:5).EQ.'XEDIT') THEN
        RETURN
      ENDIF
      END
               

      SUBROUTINE CONSTANTS(CO,noco,STRING,ERROR,*)

C#### Subroutine: CONSTANTS
C###  Description:
C###    Displays physical constants.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER noco
      CHARACTER CO(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER INSTAT,noch
      CHARACTER  OPTION(20)*60
      LOGICAL ABBREV

      CALL ENTERS('CONSTANTS',*9999)
      IF(ABBREV(CO(noco),'?',1)) THEN
        OPTION( 1)='Absolute zero of temperature = -273.2 degC'
        OPTION( 2)='Acceleration due to gravity = 9.807 m/s^2'
        OPTION( 3)='Avagadro''s number = 6.022e+23 /mol'
        OPTION( 4)='Base of natural logarithms = 2.718'
        OPTION( 5)='Boltzmann''s constant = 1.381e-23'
        OPTION( 6)='Faraday''s constant = 9.648e+4 C/mol'
        OPTION( 7)='Gas constant = 8.314 J/mol/degK'
        OPTION( 8)='Permeability of vacuum = 1.257e-6 H/m'
        OPTION( 9)='Permittivity of vacuum = 8.854e-12 F/m'
        OPTION(10)='pi = 3.14159'
        OPTION(11)='Planck''s constant = 6.626e-34 J/s'
        OPTION(12)='Velocity of light in vacuum = 2.998e+8 m/s'
        OPTION(13)='Volume of perfect gas at STP = 22.41e-3 m^3/mol'
        OPTION(14)='Return'
        CALL CHOICE('CONSTANTS',1,1,INSTAT,7,'REQUEST',14,14,noch,noco,
     '    5,CO,OPTION,STRING,0.2,0.,ERROR,*9999)

      ENDIF

      CALL EXITS('CONSTANTS')
      RETURN
 9999 CALL ERRORS('CONSTANTS',ERROR)
      CALL EXITS('CONSTANTS')
      RETURN 1
      END

      SUBROUTINE FREELU(LUMIN,LUMAX,LU,ERROR,*)

C#### Subroutine: FREELU
C###  Description:
C###    FREELU returns the first available free logical unit number, LU,
C###    from LUMIN to LUMAX.
C**** If no free unit exists control is returned via the
C**** alternate return argument.

      IMPLICIT NONE
!     Parameter List
      INTEGER LU,LUMAX,LUMIN
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IOSTAT,lucur
      LOGICAL EXIST,OPENED

      CALL ENTERS('FREELU',*9999)
      DO 1 lucur=LUMIN,LUMAX
        INQUIRE(UNIT=lucur,EXIST=EXIST,OPENED=OPENED,IOSTAT=IOSTAT)
        IF(IOSTAT.EQ.0) THEN
          IF(EXIST) THEN
            IF(.NOT.OPENED) THEN
              LU=lucur
              GOTO 2
            ENDIF
          ENDIF
        ELSE
          ERROR=' Error occurred when searching for a free logical unit'
          GOTO 9999
        ENDIF
 1    CONTINUE
      ERROR=' no free logical units available'
      GOTO 9999

 2    CALL EXITS('FREELU')
      RETURN
 9999 CALL ERRORS('FREELU',ERROR)          
      CALL EXITS('FREELU')
      RETURN 1
      END
     

      SUBROUTINE FREETIMER(IHANDLE,ERROR,*)

C#### Subroutine: FREETIMER
C###  Description:
C###    Frees a timer allocated via GETTIMER
C###  See-Also: GETTIMER

      IMPLICIT NONE
!     Parameter List
      INTEGER IHANDLE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CERROR(50),CERRLEN,ERR

      CALL ENTERS('FREETIMER',*9999)

      ERR=0
      CALL CFREETIMER(IHANDLE,ERR,CERROR)
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,ERROR)
        ERROR(CERRLEN:)=' '
        GOTO 9999
      ENDIF

      CALL EXITS('FREETIMER')
      RETURN
 9999 CALL ERRORS('FREETIMER',ERROR)
      CALL EXITS('FREETIMER')
      RETURN 1
      END


      SUBROUTINE GETTIMER(IHANDLE,ERROR,*)

C#### Subroutine: GETTIMER
C###  Description:
C###    Returns a new timer handle for use with timer

      IMPLICIT NONE
!     Parameter List
      INTEGER IHANDLE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CERROR(50),CERRLEN,ERR

      CALL ENTERS('GETTIMER',*9999)

      ERR=0
      CALL CGETTIMER(IHANDLE,ERR,CERROR)
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,ERROR)
        ERROR(CERRLEN:)=' '
        GOTO 9999
      ENDIF

      CALL EXITS('GETTIMER')
      RETURN
 9999 CALL ERRORS('GETTIMER',ERROR)
      CALL EXITS('GETTIMER')
      RETURN 1
      END


      SUBROUTINE HELP(STRING,noco,NTCO,CO,ERROR,*)

C#### Subroutine: HELP
C###  Description:
C###    Provides an on-line help facility.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER noco,NTCO
      CHARACTER CO(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,INSTAT,noch,NTCH
      CHARACTER COMMAND*50,CUPPER*50,OPTION(100)*80
      LOGICAL ABBREV,CELEM

      CALL ENTERS('HELP',*9999)
      noco=noco+1
      IF(ABBREV(CO(noco),'?',1)) THEN
        OPTION( 1)='Commands'
        OPTION( 2)='Common blocks'
        OPTION( 3)='Examples'
        OPTION( 4)='Files'
        OPTION( 5)='Image commands'
        OPTION( 6)='Image trans'
        OPTION( 7)='Include files'
        OPTION( 8)='Introduction'
        OPTION( 9)='Modules'
        OPTION(10)='New features'
        OPTION(11)='Subroutines'
        OPTION(12)='Variables'
        OPTION(13)='Exit'
        NTCH=13
        CALL CHOICE('HELP',1,1,INSTAT,7,'REQUEST',NTCH,NTCH,noch,noco,0,
     '    CO,OPTION,STRING,0.,0.,ERROR,*9999)
      ENDIF
      IF(ABBREV(CO(noco),'COMMANDS',5)) THEN
        IF(NTCO.EQ.noco) THEN !display whole file
          CALL DOCUM('commands','doc','COMMANDS',ERROR,*9999)
        ELSE                  !display specified command
          noco=noco+1
          CALL TRIM(CO(noco),IBEG,IEND)
          COMMAND=CUPPER(CO(noco)(IBEG:IEND))
          CALL DOCUM('commands','doc','Command '//COMMAND(IBEG:IEND),
     '      ERROR,*9999)
        ENDIF
      ELSE IF(ABBREV(CO(noco),'COMMON',5)) THEN
        IF(NTCO.EQ.noco) THEN !display whole file
          CALL DOCUM('common_blocks','doc','COMMON_BLOCKS',ERROR,*9999)
        ELSE                  !display specified common_block or variable
          noco=noco+1
          CALL TRIM(CO(noco),IBEG,IEND)
          COMMAND=CUPPER(CO(noco)(IBEG:IEND))
          IF(CELEM('/',COMMAND(IBEG:IEND))) THEN !specified common_block
            CALL DOCUM('common_blocks','doc','COMMON '
     '        //COMMAND(IBEG:IEND),ERROR,*9999)
          ELSE                                    !specified variable
            CALL DOCUM('common_blocks','doc','VARIABLE '
     '        //COMMAND(IBEG:IEND),ERROR,*9999)
          ENDIF
        ENDIF
      ELSE IF(ABBREV(CO(noco),'EXAMPLES',1)) THEN
        CALL DOCUM('examples','doc','EXAMPLES',ERROR,*9999)
      ELSE IF(ABBREV(CO(noco),'FILES',1)) THEN
        CALL DOCUM('help','doc','FILES',ERROR,*9999)
      ELSE IF(ABBREV(CO(noco),'IMAGE COMMANDS',7)) THEN
        CALL DOCUM('help','doc','IMAGE_COMMANDS',ERROR,*9999)
      ELSE IF(ABBREV(CO(noco),'IMAGE TRANS',5)) THEN
        CALL DOCUM('help','doc','IMAGE_TRANS',ERROR,*9999)
      ELSE IF(ABBREV(CO(noco),'INCLUDE FILES',3)) THEN
        IF(NTCO.EQ.noco) THEN !display whole file
          CALL DOCUM('common_blocks','doc','COMMON_BLOCKS',ERROR,*9999)
        ELSE                  !display specified common_block
          noco=noco+1
          CALL TRIM(CO(noco),IBEG,IEND)
          COMMAND=CUPPER(CO(noco)(IBEG:IEND))
          CALL DOCUM('common_blocks','doc','Include '
     '      //COMMAND(IBEG:IEND),ERROR,*9999)
        ENDIF
      ELSE IF(ABBREV(CO(noco),'INTRODUCTION',4)) THEN
        CALL DOCUM('introduction','doc','INTRODUCTION',ERROR,*9999)
      ELSE IF(ABBREV(CO(noco),'MODULES',1)) THEN
        IF(NTCO.EQ.noco) THEN !display whole file
          CALL DOCUM('modules','doc','MODULES',ERROR,*9999)
        ELSE                  !display specified module
          noco=noco+1
          CALL TRIM(CO(noco),IBEG,IEND)
          COMMAND=CUPPER(CO(noco)(IBEG:IEND))
          CALL DOCUM('modules','doc','MODULE '//COMMAND(IBEG:IEND),
     '      ERROR,*9999)
        ENDIF
      ELSE IF(ABBREV(CO(noco),'NEW FEATURES',1)) THEN
        CALL DOCUM('help','doc','NEW_FEATURES',ERROR,*9999)
      ELSE IF(ABBREV(CO(noco),'SUBROUTINES',1)) THEN
        IF(NTCO.EQ.noco) THEN !display whole file
          CALL DOCUM('subroutines','doc','SUBROUTINES',ERROR,*9999)
        ELSE                  !display specified subroutine
          noco=noco+1
          CALL TRIM(CO(noco),IBEG,IEND)
          COMMAND=CUPPER(CO(noco)(IBEG:IEND))
          CALL DOCUM('subroutines','doc','Subroutine '
     '      //COMMAND(IBEG:IEND),ERROR,*9999)
        ENDIF
      ELSE IF(ABBREV(CO(noco),'VARIABLES',1)) THEN
        IF(NTCO.EQ.noco) THEN !display whole file
          CALL DOCUM('variables','doc','VARIABLES',ERROR,*9999)
        ELSE                  !display specified variable
          noco=noco+1
          CALL TRIM(CO(noco),IBEG,IEND)
          COMMAND=CUPPER(CO(noco)(IBEG:IEND))
          CALL DOCUM('variables','doc','Variable '//COMMAND(IBEG:IEND),
     '      ERROR,*9999)
        ENDIF
      ELSE IF(ABBREV(CO(noco),'EXIT',3)) THEN
      ELSE
        CALL TRIM(CO(noco),IBEG,IEND)
        WRITE(OP_STRING,'('' >>no help available on '',A)')
     '    CO(noco)(IBEG:IEND)
      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('HELP')
      RETURN
 9999 CALL ERRORS('HELP',ERROR)
      CALL EXITS('HELP')
      RETURN 1
      END


C KAT 2001-12-13
      SUBROUTINE IMAGE_READ(IMP,IUNIT)

C#### Subroutine: IMAGE_READ
C###  Description:
C###    IMAGE_READ reads an unformatted image array from file IUNIT.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:pics01.cmn'
!     Parameter List
      INTEGER IMP(512,512),IUNIT
!     Local Variables
      INTEGER i,j
      CHARACTER C2P(512,512)

      READ(IUNIT) ((C2P(i,j),i=1,IMGX),j=1,IMGY)
      DO i=1,IMGX
        DO j=1,IMGY
          IMP(i,j)=ICHAR(C2P(i,j))
        ENDDO
      ENDDO

      RETURN
      END


C KAT 2001-12-13
      SUBROUTINE IMAGE_WRITE(IMP,IUNIT)

C#### Subroutine: IMAGE_WRITE
C###  Description:
C###    IMAGE_WRITE reads an unformatted image array from file IUNIT.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:pics01.cmn'
!     Parameter List
      INTEGER IMP(512,512),IUNIT
!     Local Variables
      INTEGER i,j

      WRITE(IUNIT) ((CHAR(IMP(i,j)),i=1,IMGX),j=1,IMGY)

      RETURN
      END


      SUBROUTINE MAP4(NTPTS,POINTS,PTS,ERROR,*)

C#### Subroutine: MAP4
C###  Description:
C###    MAP4 maps POINTS(nj,nopts),nopts=1,NTPTS in ityp10(1)
C###    coordinates to map coordinates PTS(nj,nopts),nopts=1,NTPTS
C###    by Cylindrical, Hammer, Polar, Rectangular or Xi projection.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:map000.cmn'
!     Parameter List
      INTEGER NTPTS
      REAL*8 POINTS(3,NTPTS),PTS(3,NTPTS)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nopts                    
      REAL*8 COSM,RK,RLAMDA,RMUU,THETA

      CALL ENTERS('MAP4',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' Projec = '',A)') PROJEC
      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      DO nopts=1,NTPTS
        IF(ITYP10(1).EQ.1) THEN
          IF(PROJEC(1:11).EQ.'RECTANGULAR') THEN
            PTS(1,nopts)=XI1MAP*(POINTS(1,nopts)-1.0D0)
     '                   +XI1MAP*POINTS(1,nopts)
            PTS(2,nopts)=XI2MAP*(POINTS(2,nopts)-1.0D0)
     '                   +XI2MAP*POINTS(2,nopts)
          ENDIF
            
        ELSE IF(ITYP10(1).EQ.2) THEN       
        ELSE IF(ITYP10(1).EQ.3) THEN
        ELSE IF(ITYP10(1).EQ.4) THEN
          IF(PROJEC(1:11).EQ.'CYLINDRICAL') THEN
            THETA=POINTS(2,nopts)+RMAP*PI/180.0D0
            IF(THETA.LT.-1.0D-2) THEN
              THETA=THETA+2.0D0*PI
            ELSE IF(THETA.GT.2.0D0*PI+1.0D-2) THEN
              THETA=THETA-2.0D0*PI
            ENDIF
            RMUU=POINTS(1,nopts)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' theta='',E12.3,'' mu='',E12.3)') 
     '          THETA,RMUU
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            PTS(1,nopts)=THETA/PI-1.0D0
            PTS(2,nopts)=2.0D0*(RMUU/PI)-1.0D0

          ELSE IF(PROJEC(1:6).EQ.'HAMMER') THEN
            THETA=POINTS(3,nopts)+RMAP*PI/180.0D0
            IF(THETA.LT.-1.0D-2) THEN
              THETA=THETA+2.0D0*PI
            ELSE IF(THETA.GT.2.0D0*PI+1.0D-2) THEN
              THETA=THETA-2.0D0*PI
            ENDIF
            RMUU=POINTS(2,nopts)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' theta='',E12.3,'' mu='',E12.3)') 
     '          THETA,RMUU
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            COSM=DCOS(RMUU-PI/2.d0)
            RK=1.0D0/DSQRT(1.0D0+COSM*DCOS((THETA-PI)/2.d0))
            PTS(1,nopts)=-RK*COSM*DSIN((THETA-PI)/2.d0)
            PTS(2,nopts)= RK*DSIN(RMUU-PI/2.d0)
                          
                                
          ELSE IF(PROJEC(1:5).EQ.'POLAR') THEN
            RLAMDA=POINTS(1,nopts)
            RMUU  =POINTS(2,nopts)
            THETA =POINTS(3,nopts)
            IF(THETA.LT.-1.0D-2) THEN
              THETA=THETA+2.0D0*PI
            ELSE IF(THETA.GT.2.0D0*PI+1.0D-2) THEN
              THETA=THETA-2.0D0*PI
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' lamda='',E12.3,'' mu='',E12.3,'
     '          //''' theta='',E12.3)') RLAMDA,RMUU,THETA
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            PTS(1,nopts)=FOCUS*DSINH(RLAMDA)*DSIN(RMUU)*DCOS(THETA)
     '                       /RMAP
            PTS(2,nopts)=FOCUS*DSINH(RLAMDA)*DSIN(RMUU)*DSIN(THETA)
     '                       /RMAP

          ELSE IF(PROJEC(1:11).EQ.'RECTANGULAR') THEN
            PTS(1,nopts)=XI1MAP*(POINTS(1,nopts)-1.0D0)
     '                   +XI1MAP*POINTS(1,nopts)
            PTS(2,nopts)=XI2MAP*(POINTS(2,nopts)-1.0D0)
     '                   +XI2MAP*POINTS(2,nopts)

          ELSE IF(PROJEC(1:2).EQ.'XI') THEN
          ENDIF

          IF(DOP) THEN
            WRITE(OP_STRING,'('' x='',E12.3,'' y='',E12.3)') 
     '        PTS(1,nopts),PTS(2,nopts)
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE IF(ITYP10(1).EQ.5) THEN
        ENDIF
      ENDDO
         CALL EXITS('MAP4')
      RETURN
 9999 CALL ERRORS('MAP4',ERROR)
      CALL EXITS('MAP4')
      RETURN 1
      END


      
      SUBROUTINE MULT(L,M,N,A,B,C)

C#### Subroutine: MULT  
C###  Description:
C###    MULT returns matrix C(l*n) as the inner product of 
C###    A(l*m) and B(m*n).

      IMPLICIT NONE
!     Parameter List
      INTEGER L,M,N
      REAL*8 A(L,M),B(M,N),C(L,N)
!     Local Variables
      INTEGER i,j,k
      REAL*8 SUM

C     CALL ENTERS('MULT',*9999)
      DO i=1,L
        DO j=1,N
          SUM=0.0D0
          DO k=1,M
            SUM=SUM+A(i,k)*B(k,j)
          ENDDO
          C(i,j)=SUM
        ENDDO
      ENDDO

C     CALL EXITS('MULT')
      RETURN
      END


      SUBROUTINE PARSIA(STRING,NXLV,NTLV,NTITLV,NXIA,IA,ERROR,*)

C#### Subroutine: PARSIA
C###  Description:
C###    PARSIA breaks a nested list of items found in STRING
C###    into separate integer numbers stored in the integer array IA.
C**** The maximum level of nesting is returned as NTLV.
C**** The size of each dimension is carried in the first NTLV
C**** components of the integer array NTITLV.
C**** Each nesting is delimited by left and right brackets
C**** while each item is separated by commas.

      IMPLICIT NONE
!     Parameter List
      INTEGER NTLV,NXIA,NXLV,IA(NXIA),NTITLV(NXLV)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER MXLV
      PARAMETER (MXLV=16)
      INTEGER IBEG,IEND,n1ch,n2ch,NOIA,NOITLV(MXLV),nolv,NTIA
      CHARACTER CHAR

      CALL ENTERS('PARSIA',*9999)
      CALL TRIM(STRING,IBEG,IEND)
      DO nolv=1,NXLV
        NTITLV(nolv)=-1
      ENDDO
      NTLV=-1
      nolv=0
      NTIA=0
      DO n2ch=IBEG,IEND+1    

            CHAR=STRING(n2ch:n2ch)
        IF(CHAR.EQ.'[') THEN
          nolv=nolv+1
          IF((nolv.GT.NXLV).OR.(nolv.GT.MXLV)) THEN
            ERROR=' Maximum array nesting exceeded'
            GOTO 9999
          ENDIF
          NOITLV(nolv)=0
          n1ch=n2ch
        ELSE IF(CHAR.EQ.']') THEN
          IF(NTLV.LT.0) THEN
            NTLV=nolv
          ENDIF
          IF(nolv.EQ.NTLV) THEN
            CALL PARSIL(STRING(n1ch+1:n2ch-1),NXIA-NTIA,NOIA,IA(NTIA+1),
     '        ERROR,*9999)
            NOITLV(nolv)=NOIA
            NTIA=NTIA+NOIA
          ENDIF
          IF(NTITLV(nolv).LT.0) THEN
            NTITLV(nolv)=NOITLV(nolv)
          ELSE IF(NTITLV(nolv).NE.NOITLV(nolv)) THEN
            ERROR=' Inconsistent number of items in array'
            GOTO 9999
          ENDIF
          nolv=nolv-1
          IF(nolv.LT.0) THEN
            ERROR=' Inconsistent nesting of array'
            GOTO 9999
          ENDIF
          NOITLV(nolv)=NOITLV(nolv)+1
        ENDIF
      ENDDO !n2ch

      IF(nolv.EQ.0) THEN 

               CHAR=STRING(n2ch:n2ch)
        IF(CHAR.EQ.'[') THEN
          nolv=nolv+1
          IF((nolv.GT.NXLV).OR.(nolv.GT.MXLV)) THEN
            ERROR=' Maximum array nesting exceeded'
            GOTO 9999
          ENDIF
          NOITLV(nolv)=0
          n1ch=n2ch
        ELSE IF(CHAR.EQ.']') THEN
          IF(NTLV.LT.0) THEN
            NTLV=nolv
          ENDIF
          IF(nolv.EQ.NTLV) THEN
            CALL PARSIL(STRING(n1ch+1:n2ch-1),NXIA-NTIA,NOIA,IA(NTIA+1),
     '        ERROR,*9999)
            NOITLV(nolv)=NOIA
            NTIA=NTIA+NOIA
          ENDIF
          IF(NTITLV(nolv).LT.0) THEN
            NTITLV(nolv)=NOITLV(nolv)
          ELSE IF(NTITLV(nolv).NE.NOITLV(nolv)) THEN
            ERROR=' Inconsistent number of items in array'
            GOTO 9999
          ENDIF
          nolv=nolv-1
          IF(nolv.LT.0) THEN
            ERROR=' Inconsistent nesting of array'
            GOTO 9999
          ENDIF
          NOITLV(nolv)=NOITLV(nolv)+1
        ENDIF
      ENDDO !n2ch

      IF(nolv.EQ.0) THEN
         CHAR=STRING(n2ch:n2ch)
        IF(CHAR.EQ.'[') THEN
          nolv=nolv+1
          IF((nolv.GT.NXLV).OR.(nolv.GT.MXLV)) THEN
            ERROR=' Maximum array nesting exceeded'
            GOTO 9999
          ENDIF
          NOITLV(nolv)=0
          n1ch=n2ch
        ELSE IF(CHAR.EQ.']') THEN
          IF(NTLV.LT.0) THEN
            NTLV=nolv
          ENDIF
          IF(nolv.EQ.NTLV) THEN
            CALL PARSIL(STRING(n1ch+1:n2ch-1),NXIA-NTIA,NOIA,IA(NTIA+1),
     '        ERROR,*9999)
            NOITLV(nolv)=NOIA
            NTIA=NTIA+NOIA
          ENDIF
          IF(NTITLV(nolv).LT.0) THEN
            NTITLV(nolv)=NOITLV(nolv)
          ELSE IF(NTITLV(nolv).NE.NOITLV(nolv)) THEN
            ERROR=' Inconsistent number of items in array'
            GOTO 9999
          ENDIF
          nolv=nolv-1
          IF(nolv.LT.0) THEN
            ERROR=' Inconsistent nesting of array'
            GOTO 9999
          ENDIF
          NOITLV(nolv)=NOITLV(nolv)+1
        ENDIF
      ENDDO !n2ch

      IF(nolv.EQ.0) THEN 
        NTITLV(1)=NOITLV(1)
      ELSE
        ERROR=' Inconsistent nesting of array'
        GOTO 9999
      ENDIF

      CALL EXITS('PARSIA')
      RETURN
 9999 CALL ERRORS('PARSIA',ERROR)
      CALL EXITS('PARSIA')
      RETURN 1
      END 


      SUBROUTINE PARSLO(STRING,LO,ERROR,*)

C#### Subroutine: PARSLO
C###  Description:
C###    PARSLO converts the character string STRING into a logical 
C###    variable LO.

      IMPLICIT NONE
!     Parameter List
      CHARACTER ERROR*(*),STRING*(*)
      LOGICAL LO
!     Local Variables
      INTEGER IOSTAT

      CALL ENTERS('PARSLO',*9999)
      READ(UNIT=STRING,FMT=*,IOSTAT=IOSTAT) LO
      IF(IOSTAT.LT.0) THEN
        ERROR=' no items were found in STRING'
        GOTO 9999
      ELSE IF(IOSTAT.GT.0) THEN
        ERROR=' An error occurred while parsing STRING'
        GOTO 9999
      ENDIF

      CALL EXITS('PARSLO')
      RETURN
 9999 CALL ERRORS('PARSLO',ERROR)
      CALL EXITS('PARSLO')
      RETURN 1
      END


       SUBROUTINE POLY(X,na,A,nd,PX)

C#### Subroutine: POLY
C###  Description:
C###    POLY evaluates the complex polynomial 
C###    A(0) + A(1)*X + .. A(na)*(X**na)
C###    and its first nd derivatives using a nested Horner's scheme.
C***  On exit the nd'th derivative is stored in PX(nd).
C**** NOTE: The dimension of PX must be at least MAX(na,nd).

      IMPLICIT NONE
!     Parameter List
      INTEGER na,nd
      COMPLEX*16 A(0:*),PX(0:*),X
!     Local Variables
      INTEGER IFACT,n1,n2,n3,n4,n5

C     CALL ENTERS('POLY',*9999)
      DO n1=0,na
        PX(n1)=A(n1)
      ENDDO
      DO n3=0,MIN((na-1),nd)
        DO n2=(na-1),n3,-1
          PX(n2)=PX(n2)+PX(n2+1)*X
        ENDDO
      ENDDO
      IFACT=1
      DO n4=2,MIN(na,nd)
        IFACT=IFACT*n4
        PX(n4)=PX(n4)*IFACT
      ENDDO
      DO n5=MIN(na,nd)+1,MAX(na,nd)
        PX(n5)=(0.0D0,0.0D0)
      ENDDO

C     CALL EXITS('POLY')
      RETURN
      END
      

      SUBROUTINE POST(iw,FILE,ERROR,*)

C#### Subroutine: POST
C###  Description:
C###    POST converts linefeeds to returns in postscript files.
C**** If IW=3 then PHIGS files, else GKS files. !!!OLD!
C**** Note change to 256 char/line files for new GKS & PHIGS (8-Feb-90)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:post00.cmn'
!     Parameter List
      INTEGER iw
      CHARACTER ERROR*(*),FILE*(*)
!     Local Variables
      INTEGER i,IBEG,IEND,IS,ISLAST
      CHARACTER FILENAME*100,LAST*256,LINE*256
      LOGICAL MORE

      CALL ENTERS('POST',*9999)
      CALL TRIM(FILE,IBEG,IEND)
      FILENAME=FILE(IBEG:IEND)
      IF(IWKT(iw).EQ.1) THEN      !GKS
        CALL OPENF(15,'DISK',FILENAME(IBEG:IEND)//'.GKS','OLD',
     '    'SEQUEN','FORMATTED',256,ERROR,*9999)
      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
        CALL OPENF(15,'DISK',FILENAME(IBEG:IEND)//'.PHIGS','OLD',
     '    'SEQUEN','FORMATTED',256,ERROR,*9999)
      ENDIF
      IF(EPS) THEN
     	CALL OPENF(16,'DISK',FILENAME(IBEG:IEND)//'.EPS','NEW',
     '	  'SEQUEN','FORMATTED',256,ERROR,*9999)
      ELSE
      	CALL OPENF(16,'DISK',FILENAME(IBEG:IEND)//'.PS','NEW',
     '	  'SEQUEN','FORMATTED',256,ERROR,*9999)
      ENDIF
      MORE=.FALSE.
 10  READ(15,'(A)',END=9997) LINE
        CALL TRIM(LINE,IBEG,IEND)
        IS=1
        i=0
        DO WHILE(I.LT.IEND)
          i=i+1
          IF(ICHAR(LINE(i:i)).EQ.10) THEN
            IF(MORE) THEN
              MORE=.FALSE.
              IF(IS.EQ.i)THEN
                WRITE(16,'(A)') ' '//LAST(ISLAST:)
              ELSE
                WRITE(16,'(A)') ' '//LAST(ISLAST:)//LINE(IS:i-1)
              ENDIF
            ELSE IF(IS.LT.i)THEN
              WRITE(16,'(A)') ' '//LINE(IS:i-1)
            ENDIF
            IS=i+1
          ENDIF
        ENDDO
        IF(IS.GT.IEND) THEN
          MORE=.FALSE.
        ELSE
          MORE=.TRUE.
          ISLAST=IS
          LAST(IS:)=LINE(IS:)
        ENDIF
        GOTO 10
 9997 CONTINUE
      CALL CLOSEF(15,ERROR,*9999)
      CALL CLOSEF(16,ERROR,*9999)

 9998 CALL EXITS('POST')
      RETURN
 9999 CALL ERRORS('POST',ERROR)
      CALL EXITS('POST')
      RETURN 1
      END


  
      SUBROUTINE RESETTIMERS(ERROR,*)

C#### Subroutine: RESETTIMERS
C###  Description:
C###    Resets all timer handles for later use with timer

      IMPLICIT NONE
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CERROR(50),CERRLEN,ERR

      CALL ENTERS('RESETTIMERS',*9999)

      ERR=0
      CALL CRESETTIMERS(ERR,CERROR)
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,ERROR)
        ERROR(CERRLEN:)=' '
        GOTO 9999
      ENDIF

      CALL EXITS('RESETTIMERS')
      RETURN
 9999 CALL ERRORS('RESETTIMERS',ERROR)
      CALL EXITS('RESETTIMERS')
      RETURN 1
      END


      SUBROUTINE PUTRAN(IUNIT,TRANS)

C#### Subroutine: PUTRAN
C###  Description:
C###    PUTRAN writes out the transformation matrix TRANS to IUNIT.

      IMPLICIT NONE
!     Parameter List
      INTEGER IUNIT
      REAL*8 TRANS(3,4)
!     Local Variables
      INTEGER i,j

C     CALL ENTERS('PUTRAN',*9999)
      REWIND( IUNIT )
      WRITE(IUNIT,'(/)')
      DO 1 i=1,3
        WRITE(IUNIT,'('' TRANS( '',I1,'',..) ='',4F16.8)')
     '       i,(TRANS(i,j),j=1,4)
 1    CONTINUE

C     CALL EXITS('PUTRAN')
      RETURN
      END
                        
       
      SUBROUTINE ROOTSC(NA,A,NX,X,TOL,ITMAX,IFAIL)

C**** Finds NX complex roots of the polynomial
C**** A(0) + A(1)*X + A(2)*(X**2) + ..+ A(NA)*(X**NA)
C**** using the method of Laguerre.
C**** On exit the roots are stored in X(1)..X(NX).
C**** IFAIL=1 if convergence is not reached, otherwise IFAIL=0.
C**** NOTE: NA must be less than or equal to 15.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER NA,NX,IFAIL,ITMAX
      COMPLEX*16 A(0:*),X(*)
      REAL*8 TOL
!     Local Variables
      CHARACTER ERROR*10 
      COMPLEX*16 D,DX,DXOLD,DP,DN,PX(0:15)
      INTEGER I1,I1X,I2,IT,N1,N2

C     CALL ENTERS('ROOTSC',*9999)
      NX=NA
      DO I1=NA,0,-1
        IF(CDABS(A(I1)).GT.1.0D-15) GOTO 20
        NX=NX-1
        X(I1)=(0.0D0,0.0D0)
      ENDDO
 20   DO I1=1,NX
        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' NROOT='',I4)') I1
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        N1=NX-I1+1
        N2=I1-1
        X(I1)=(0.0D0,0.0D0)
        DX   =(0.0D0,0.0D0)
        DO IT=1,ITMAX
          CALL POLY(X(I1),N1,A,2,PX)
          DXOLD=DX
          D=CDSQRT((N1-1)*((N1-1)*(PX(1)*PX(1))-DBLE(N1)*PX(0)*
     '      PX(2)))
          DP=PX(1)+D
          DN=PX(1)-D
          IF(CDABS(DP).GT.CDABS(DN)) THEN
            IF(CDABS(DP).GT.1.0D-15) THEN
              DX=N1*PX(0)/DP
            ELSE
              DX=(0.0D0,0.0D0)
            ENDIF
          ELSE
            IF(CDABS(DN).GT.1.0D-15) THEN
              DX=N1*PX(0)/DN
            ELSE
              DX=(0.0D0,0.0D0)
            ENDIF
          ENDIF
          X(I1)=X(I1)-DX
          IF(DOP) THEN
            WRITE(OP_STRING,'(''    IT='',I4,3X,
     '        ''   X= ('',D12.6,'','',D12.6,'')'',
     '        ''  DX= ('',D12.6,'','',D12.6,'')'',
     '        ''  PX= ('',D12.6,'','',D12.6,'')'',
     '        /14X, '' DPX= ('',D12.6,'','',D12.6,'')'',
     '        ''D2PX= ('',D12.6,'','',D12.6,'')'')')
     '        IT,X(I1),DX,PX(0),PX(1),PX(2)
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)      
          ENDIF
          IF((CDABS(DX).LE.TOL).AND.(CDABS(PX(0)).LE.TOL)) GO TO 40
        ENDDO
 40     DO I2=(NA-1),N2,-1
          A(I2)=A(I2+1)*X(I1)+A(I2)
        ENDDO
        IF(DOP) THEN
          WRITE(OP_STRING,'(12X,''   Root= ('',D12.6,'','',D12.6,
     '      '')'')') (X(I1X),I1X=1,NX)
      	  CALL WRITES(IODI,OP_STRING,ERROR,*9999)      
        ENDIF
      ENDDO

C     CALL EXITS('ROOTSC')
 9999 RETURN
      END

        SUBROUTINE SCALE0(AXIS,AMOUNT,TRANS)

C#### Subroutine: SCALE0
C###  Description:
C###    SCALE0 updates the transformation matrix TRANS with a scaling
C###    by AMOUNT along the axis defined by AXIS.
C**** If AXIS is zero a scaling by AMOUNT is performed in all directions

      IMPLICIT NONE
!     Parameter List
      REAL*8 AMOUNT,AXIS(3),TRANS(3,4)
!     Local Variables
      REAL*8 AMOUN0

C     CALL ENTERS('SCALE0',*9999)
      IF(AXIS(1).NE.0.0D0) THEN
        IF(AXIS(2).NE.0.0D0) THEN
          IF(AXIS(3).NE.0.0D0) THEN
            AMOUN0=AMOUNT/DSQRT(AXIS(1)**2+AXIS(2)**2+AXIS(3)**2)
            CALL SCALE1(AMOUN0*AXIS(1),TRANS)
            CALL SCALE2(AMOUN0*AXIS(2),TRANS)
            CALL SCALE3(AMOUN0*AXIS(3),TRANS)
          ELSE
            AMOUN0=AMOUNT/DSQRT(AXIS(1)**2+AXIS(2)**2)
            CALL SCALE1(AMOUN0*AXIS(1),TRANS)
            CALL SCALE2(AMOUN0*AXIS(2),TRANS)
          ENDIF      
        ELSE
          IF(AXIS(3).NE.0.0D0) THEN
            AMOUN0=AMOUNT/DSQRT(AXIS(1)**2+AXIS(3)**2)
            CALL SCALE1(AMOUN0*AXIS(1),TRANS)            
            CALL SCALE3(AMOUN0*AXIS(3),TRANS)
          ELSE
            CALL SCALE1( AMOUNT,TRANS)
          ENDIF
        ENDIF
      ELSE
        IF(AXIS(2).NE.0.0D0) THEN
          IF(AXIS(3).NE.0.0D0) THEN
            AMOUN0=AMOUNT/DSQRT(AXIS(2)**2+AXIS(3)**2)
            CALL SCALE2(AMOUN0*AXIS(2),TRANS)
            CALL SCALE3(AMOUN0*AXIS(3),TRANS)
          ELSE
            CALL SCALE2(AMOUNT,TRANS)
          ENDIF
        ELSE
          IF(AXIS(3).NE.0.0D0) THEN
            CALL SCALE3(AMOUNT,TRANS)
          ELSE
            CALL SCALE1(AMOUNT,TRANS)
            CALL SCALE2(AMOUNT,TRANS)
            CALL SCALE3(AMOUNT,TRANS)
          ENDIF
        ENDIF
      ENDIF

C     CALL EXITS('SCALE0')
      RETURN
      END


C KAT 2001-12-14
      SUBROUTINE SEGPAC(ISALIG,ISAXES,ISBASE,ISCLOC,ISCONO,ISCONT,
     '  ISCROS,ISDANO,ISDAPR,ISDATA,ISDATR,ISEG,ISELEC,ISELNO,ISFACE,
     '  ISFANO,ISFIBR,ISFIEL,ISGAUS,ISGRAD,ISGRID,ISHIST,ISIMAG,ISINCR,
     '  ISISOC,ISLINE,ISLINO,ISL2BE,ISL3BE,ISMAP,ISMATE,ISNONO,ISN2BE,
     '  ISN3BE,ISOBJE,ISPLIN,ISPLOT,ISPROF,ISREAC,ISRULE,ISSCAL,ISSECT,
     '  ISSHEE,ISSIGN,ISSTRA,ISSTRE,ISSTRM,ISSURF,ISTEXT,ISTRAC,ISVELO,
     '  NEELEM,NPNODE,CSEG,ERROR,*)

C#### Subroutine: SEGPAC
C###  Description:
C###    SEGPAC packs segments down to eliminate cancelled segments.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:map000.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'cmiss$reference:obje00.cmn'
!     Parameter List
      INTEGER ISALIG(NWM),ISAXES(NWM),ISBASE(99),
     '  ISCLOC(NWM),ISCONO(NHM,NEM),ISCONT(NHM,NEM,NGRSEGM),
     '  ISCROS(NWM,NGRSEGM),
     '  ISDANO(NWM,NEM),ISDAPR(NWM,NEM),ISDATA(NWM,NGRSEGM),
     '  ISDATR(NWM,NEM),
     '  ISEG(*),ISELEC(128),ISELNO(NWM,NEM),ISFACE(NWM,NFM),
     '  ISFANO(NWM,NFM),ISFIBR(NWM,NEM,NGRSEGM),ISFIEL(NWM,NEM),
     '  ISGAUS(NWM,NGM,NEM),ISGRAD(NEM,NGRSEGM),ISGRID(NWM),
     '  ISHIST(0:NPM),ISIMAG(NWM),ISINCR(NWM),ISISOC(NWM),
     '  ISLINE(NWM,2*NGRSEGM),ISLINO(NWM),
     '  ISL2BE(NLM),ISL3BE(NLM),ISMAP(NGRSEGM),
     '  ISMATE(NWM,NEM),ISNONO(NWM,NPM),ISN2BE(NLM),ISN3BE(NLM),
     '  ISOBJE(NWM,NGRSEGM,NGRSEGM),ISPLIN(NWM,NGRSEGM),
     '  ISPLOT(NHM,0:NEM,NGRSEGM),
     '  ISPROF,ISREAC(NWM),ISRULE(NWM),ISSCAL(NWM,NGRSEGM),
     '  ISSECT(NGRSEGM),
     '  ISSHEE(NWM,NEM,NGRSEGM),ISSIGN(2,128),ISSTRA(NEM,NGRSEGM),
     '  ISSTRE(NEM,NGRSEGM),ISSTRM(NEM,NGRSEGM),ISSURF(NWM,NGRSEGM),
     '  ISTEXT(NWM),ISTRAC(10),ISVELO(NEM,NGRSEGM),
     '  NEELEM(0:NE_R_M,0:NRM),NPNODE(0:NP_R_M,0:NRM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER i,iw,n1sg,nb,ne,nf,ng,nh,nl,np,nocont,nocros,
     '  nodata,noelec,noelem,nofibr,nograd,noline,nomap,nonode,noobje,
     '  nopart,noplin,noplot,noscal,nosect,nosg,noshee,nostra,nostre,
     '  nostrm,nosurf,novelo,nr,NTSG0
      LOGICAL CHANGE

      CALL ENTERS('SEGPAC',*9999)
 100  CHANGE=.FALSE.
      NTSG0=NTSG
      DO nosg=1,NTSG0
        IF(ISEG(nosg).EQ.0) THEN !shift all segs above NOSG down by 1
          CHANGE=.TRUE.
          NTSG=NTSG-1
          DO n1sg=nosg,NTSG
            ISEG(n1sg)=ISEG(n1sg+1)
            CSEG(n1sg)=CSEG(n1sg+1)
          ENDDO
          IF(ISPROF.GT.nosg) ISPROF=ISPROF-1
          DO iw=1,NWM
            IF(ISALIG(iw).GT.nosg) ISALIG(iw)=ISALIG(iw)-1
            IF(ISAXES(iw).GT.nosg) ISAXES(iw)=ISAXES(iw)-1
            IF(ISCLOC(iw).GT.nosg) ISCLOC(iw)=ISCLOC(iw)-1
            IF(ISRULE(iw).GT.nosg) ISRULE(iw)=ISRULE(iw)-1
            IF(ISIMAG(iw).GT.nosg) ISIMAG(iw)=ISIMAG(iw)-1
            IF(ISINCR(iw).GT.nosg) ISINCR(iw)=ISINCR(iw)-1
            IF(ISISOC(iw).GT.nosg) ISISOC(iw)=ISISOC(iw)-1
            IF(ISLINO(iw).GT.nosg) ISLINO(iw)=ISLINO(iw)-1
            IF(ISREAC(iw).GT.nosg) ISREAC(iw)=ISREAC(iw)-1
            IF(ISTEXT(iw).GT.nosg) ISTEXT(iw)=ISTEXT(iw)-1
            DO nocros=1,NTCROS
              IF(ISCROS(iw,nocros).GT.nosg)
     '          ISCROS(iw,nocros)=ISCROS(iw,nocros)-1
            ENDDO
            DO nodata=1,NTDATA
              IF(ISDATA(iw,nodata).GT.nosg)
     '          ISDATA(iw,nodata)=ISDATA(iw,nodata)-1
            ENDDO
            DO nf=1,NFT
              IF(ISFACE(iw,nf).GT.nosg) ISFACE(iw,nf)=ISFACE(iw,nf)-1
              IF(ISFANO(iw,nf).GT.nosg) ISFANO(iw,nf)=ISFANO(iw,nf)-1
            ENDDO
            IF(ISGRID(iw).GT.nosg) ISGRID(iw)=ISGRID(iw)-1
            DO noline=1,NTLINE
              IF(ISLINE(iw,noline).GT.nosg)
     '          ISLINE(iw,noline)=ISLINE(iw,noline)-1
            ENDDO
            DO noplin=1,NTPLIN
              IF(ISPLIN(iw,noplin).GT.nosg)
     '          ISPLIN(iw,noplin)=ISPLIN(iw,noplin)-1
            ENDDO
            DO noscal=1,NTSCAL
              IF(ISSCAL(iw,noscal).GT.nosg)
     '          ISSCAL(iw,noscal)=ISSCAL(iw,noscal)-1
            ENDDO
            DO nosurf=1,NTSURF
              IF(ISSURF(iw,nosurf).GT.nosg)
     '          ISSURF(iw,nosurf)=ISSURF(iw,nosurf)-1
            ENDDO
            DO nr=1,NRT
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(ISELNO(iw,ne).GT.nosg) ISELNO(iw,ne)=ISELNO(iw,ne)-1
                IF(ISDANO(iw,ne).GT.nosg) ISDANO(iw,ne)=ISDANO(iw,ne)-1
                IF(ISDAPR(iw,ne).GT.nosg) ISDAPR(iw,ne)=ISDAPR(iw,ne)-1
                IF(ISDATR(iw,ne).GT.nosg) ISDATR(iw,ne)=ISDATR(iw,ne)-1
                IF(ISFIEL(iw,ne).GT.nosg) ISFIEL(iw,ne)=ISFIEL(iw,ne)-1
                IF(ISMATE(iw,ne).GT.nosg) ISMATE(iw,ne)=ISMATE(iw,ne)-1
                DO nofibr=1,NTFIBR
                  IF(ISFIBR(iw,ne,nofibr).GT.nosg)
     '              ISFIBR(iw,ne,nofibr)=ISFIBR(iw,ne,nofibr)-1
                ENDDO
                DO noshee=1,NTSHEE
                  IF(ISSHEE(iw,ne,noshee).GT.nosg)
     '              ISSHEE(iw,ne,noshee)=ISSHEE(iw,ne,noshee)-1
                ENDDO
                DO ng=1,NGM
                  IF(ISGAUS(iw,ng,ne).GT.nosg)
     '              ISGAUS(iw,ng,ne)=ISGAUS(iw,ng,ne)-1
                ENDDO
              ENDDO
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                IF(ISNONO(iw,np).GT.nosg) ISNONO(iw,np)=ISNONO(iw,np)-1
              ENDDO
            ENDDO
            DO nopart=1,NRM
              DO noobje=1,NTOBJE
                IF(ISOBJE(iw,noobje,nopart).GT.nosg)
     '            ISOBJE(iw,noobje,nopart)=ISOBJE(iw,noobje,nopart)-1
              ENDDO
            ENDDO
          ENDDO !iw
          DO nh=1,NHM
            DO noplot=1,NTPLOT
              IF(ISPLOT(nh,0,noplot).GT.nosg)
     '          ISPLOT(nh,0,noplot)=ISPLOT(nh,0,noplot)-1
            ENDDO
          ENDDO
          DO nr=1,NRT
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              DO nh=1,NHM
                IF(ISCONO(nh,ne).GT.nosg) ISCONO(nh,ne)=ISCONO(nh,ne)-1
                DO nocont=1,NTCONT
                  IF(ISCONT(nh,ne,nocont).GT.nosg)
     '              ISCONT(nh,ne,nocont)=ISCONT(nh,ne,nocont)-1
                ENDDO
                DO noplot=1,NTPLOT
                  IF(ISPLOT(nh,ne,noplot).GT.nosg)
     '              ISPLOT(nh,ne,noplot)=ISPLOT(nh,ne,noplot)-1
                ENDDO
              ENDDO
              DO nograd=1,NTGRAD
                IF(ISGRAD(ne,nograd).GT.nosg)
     '            ISGRAD(ne,nograd)=ISGRAD(ne,nograd)-1
              ENDDO
              DO nostra=1,NTSTRA
                IF(ISSTRA(ne,nostra).GT.nosg)
     '            ISSTRA(ne,nostra)=ISSTRA(ne,nostra)-1
              ENDDO
              DO nostre=1,NTSTRE
                IF(ISSTRE(ne,nostre).GT.nosg)
     '            ISSTRE(ne,nostre)=ISSTRE(ne,nostre)-1
              ENDDO
              DO nostrm=1,NTSTRM
                IF(ISSTRM(ne,nostrm).GT.nosg)
     '            ISSTRM(ne,nostrm)=ISSTRM(ne,nostrm)-1
              ENDDO
              DO novelo=1,NTVELO
                IF(ISVELO(ne,novelo).GT.nosg)
     '            ISVELO(ne,novelo)=ISVELO(ne,novelo)-1
              ENDDO
            ENDDO !noelem
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              IF(ISHIST(np).GT.nosg) ISHIST(np)=ISHIST(np)-1
            ENDDO !nonode
          ENDDO !nr
         
          DO nb=1,NBFT
            IF(ISBASE(nb).GT.nosg) ISBASE(nb)=ISBASE(nb)-1
          ENDDO
          DO noelec=1,128
            IF(ISELEC(noelec).GT.nosg) 
     '        ISELEC(noelec)=ISELEC(noelec)-1
            IF(ISSIGN(1,noelec).GT.nosg) 
     '        ISSIGN(1,noelec)=ISSIGN(1,noelec)-1
            IF(ISSIGN(2,noelec).GT.nosg) 
     '        ISSIGN(2,noelec)=ISSIGN(2,noelec)-1
          ENDDO
          DO i=1,10
            IF(ISTRAC(i).GT.nosg) ISTRAC(i)=ISTRAC(i)-1
          ENDDO
          DO np=0,NPT(nr)
            IF(ISHIST(np).GT.nosg) ISHIST(np)=ISHIST(np)-1
          ENDDO
          DO nl=1,NLM
            IF(ISL2BE(nl).GT.nosg) ISL2BE(nl)=ISL2BE(nl)-1
            IF(ISL3BE(nl).GT.nosg) ISL3BE(nl)=ISL3BE(nl)-1
            IF(ISN2BE(nl).GT.nosg) ISN2BE(nl)=ISN2BE(nl)-1
            IF(ISN3BE(nl).GT.nosg) ISN3BE(nl)=ISN3BE(nl)-1
          ENDDO
          DO nomap=1,NTMAP
            IF(ISMAP(nomap).GT.nosg) ISMAP(nomap)=ISMAP(nomap)-1
          ENDDO
          DO nosect=1,NRM
            IF(ISSECT(nosect).GT.nosg) ISSECT(nosect)=ISSECT(nosect)-1
          ENDDO
        ENDIF
      ENDDO
      IF(CHANGE) GO TO 100

      CALL EXITS('SEGPAC')
      RETURN
 9999 CALL ERRORS('SEGPAC',ERROR)
      CALL EXITS('SEGPAC')
      RETURN 1
      END


C KAT 2001-12-13
      SUBROUTINE SIGNAL_READ(ISIG,IUNIT,NUMSIGNALS,NUMSAMPLES)

C#### Subroutine: SIGNAL_READ
C###  Description:
C###    SIGNAL_READ reads an unformatted signal array from file IUNIT.

      IMPLICIT NONE
!     Parameter List
      INTEGER NUMSIGNALS,NUMSAMPLES,IUNIT
      INTEGER ISIG(NUMSIGNALS,NUMSAMPLES)
      
      READ(IUNIT) ISIG

      RETURN
      END


      SUBROUTINE STACK(COMAND,NVALUE,VALUE,FAIL)

C#### Subroutine: STACK
C###  Description:
C###    STACK controls 'PUSH', 'PULL' and 'CLEAR' operations 
C###    defined in COMAND on the stack.
C**** The input values provided for the command 'PUSH' and
C**** the output values for 'PULL' are carried in VALUE(1..NVALUE).
C**** If an attempt is made to push a value on the stack when the
C**** stack is full FAIL will return 'STACK OVERFLOW'.
C**** If an attempt is made to pull a value off the stack when the
C**** stack is empty FAIL will return 'STACK UNDERFLOW'.
C**** If a command other than 'PUSH' or 'PULL' is given to STACK
C**** via COMAND FAIL will return 'COMMAND ERROR'.
C**** With a successful invocation of STACK FAIL returns ''.

      IMPLICIT NONE
!     Parameter List
      INTEGER NVALUE
      REAL*8 VALUE(*)
      CHARACTER COMAND*(*),FAIL*(*)
!     Local Variables
      INTEGER IDXSTK,MAXSTK,nval
      PARAMETER (MAXSTK=100)
      REAL*8 STK(MAXSTK)
      DATA IDXSTK/0/

C     CALL ENTERS('STACK',*9999)
      IF(COMAND.EQ.'PUSH') THEN
        IF((IDXSTK+NVALUE).GT.MAXSTK) THEN
          FAIL='STACK OVERFLOW'
        ELSE
          FAIL=' '
          DO nval=1,NVALUE
            IDXSTK=IDXSTK+1
            STK(IDXSTK)=VALUE(nval)
          ENDDO
        ENDIF  
      ELSE IF(COMAND.EQ.'PULL') THEN
        IF((IDXSTK-NVALUE).LT.0) THEN
          FAIL='STACK UNDERFLOW'
        ELSE
          FAIL=' '
          DO nval=NVALUE,1,-1
            VALUE(nval)=STK(IDXSTK)
            IDXSTK=IDXSTK-1
          ENDDO
        ENDIF
      ELSE IF(COMAND.EQ.'CLEAR') THEN
        IDXSTK=0
      ELSE
        FAIL='COMMAND ERROR'
      ENDIF

C     CALL EXITS('STACK')
      RETURN
      END       

      SUBROUTINE SIMPS(XA_LOCAL,YA,NSIMP,XKSIA,XKSIB,ACDIFF)

C#### Subroutine: SIMPS
C###  Description:
C###    SIMPS does Simpson's rule calculation of arc length
C###    (DAVE RAHDERT).

      IMPLICIT NONE
!     Parameter List
      INTEGER NSIMP
      REAL*8 XA_LOCAL(4),XKSIA,YA(4),XKSIB
!     Local Variables
      INTEGER kk,KK1
      REAL*8 ACDIFF,ACN,BCN,W1,XFA

      ACDIFF=0.0d0
      DO kk=1,NSIMP+1      !SIMPSON'S RULE
        XFA=XKSIA+(XKSIB-XKSIA)*(DBLE(kk)-1.D0)/DBLE(NSIMP)
        ACN=XA_LOCAL(1)*(-3.0D0*(1.0D0-XFA)**2)
        ACN=ACN+XA_LOCAL(2)*(3.0D0*(1.0D0-XFA)*(1.0D0-3.0D0*XFA))
        ACN=ACN+XA_LOCAL(3)*(3.0D0*XFA*(2.0D0-3.0D0*XFA))
        ACN=ACN+XA_LOCAL(4)*3.0D0*XFA*XFA

        BCN=YA(1)*(-3.0D0*(1.0D0-XFA)**2)
        BCN=BCN+YA(2)*(3.0D0*(1.0D0-XFA)*(1.0D0-3.0D0*XFA))
        BCN=BCN+YA(3)*(3.0D0*XFA*(2.0D0-3.0D0*XFA))
        BCN=BCN+YA(4)*3.0D0*XFA*XFA

        W1=2.0D0                  !PRESUME KK IS ODD
        KK1=kk/2
        KK1=KK1*2
        IF(KK1.EQ.kk) W1=4.0D0   !KK IS EVEN
        IF(kk.EQ.1.OR.kk.EQ.NSIMP+1) W1=1.0D0    !KK IS ENDPT.
        ACDIFF=ACDIFF+((ACN*ACN+BCN*BCN)**0.5D0)*(XKSIB-XKSIA)
     '         *W1/DBLE(NSIMP)/3.0D0
      ENDDO   !end KK

      RETURN
      END

      SUBROUTINE SPLINE(IBT,nb,NEELEM,VE,ERROR,*)

C#### Subroutine: SPLINE
C###  Description:
C###    SPLINE calculates B-spline polynomial coefficients and stores
C###    them in VE(ns,nk,ne),  where nk=1,NKT(0,nb) are polynomial 
C###    coefficients for each dof ns=1,NST(nb) for each element ne.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),nb,NEELEM(0:NE_R_M,0:NRM)
      REAL*8 VE(NSM,NKM,NEFM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER k1,k2,k3,n1,n2,n3,ne,nk,noelem,nr,ns
      REAL*8 PSE

      CALL ENTERS('SPLINE',*9999)
      DO nr=1,NRT
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          ns=0
          DO 48 n3=1,IBT(2,3)+1
          DO 48 n2=1,IBT(2,2)+1
          DO 48 n1=1,IBT(2,1)+1
            ns=ns+1
            nk=0
            DO 45 k3=1,IBT(2,3)+1
            DO 45 k2=1,IBT(2,2)+1
            DO 45 k1=1,IBT(2,1)+1
              nk=nk+1
              VE(ns,nk,ne)=PSE(k1,k2,k3,n1,n2,n3,IBT,NIT(nb))
 45         CONTINUE
 48       CONTINUE
        ENDDO
      ENDDO

      CALL EXITS('SPLINE')
      RETURN
 9999 CALL ERRORS('SPLINE',ERROR)
      CALL EXITS('SPLINE')
      RETURN 1
      END

                               
      SUBROUTINE UNITS(CO,noco,STRING,ERROR,*)

C#### Subroutine: UNITS
C###  Description:
C###    Displays unit conversions.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER noco
      CHARACTER CO(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IDATA(1),ID_DEVICE,ID_STATUS,ID_WS,INPUT_CHOICE,INSTAT,
     '  noch
      REAL R4DATA
      CHARACTER CLASS*8,OPTION(20)*30,OPTION1(20)*30,OPTION2(20)*30,
     '  SDATA*10
      LOGICAL ABBREV,CONTINUE

      CALL ENTERS('UNITS',*9999)
      IF(ABBREV(CO(noco),'?',1)) THEN
        OPTION( 1)='Energy'
        OPTION( 2)='Force'
        OPTION( 3)='Length'
        OPTION( 4)='Mass'
        OPTION( 5)='Power'
        OPTION( 6)='Specific heat'
        OPTION( 7)='Stress'
        OPTION( 8)='Viscosity'
        OPTION( 9)='Return'
        CALL CHOICE('UNITS',1,1,INSTAT,7,'EVENT',9,9,noch,noco,5,
     '    CO,OPTION,STRING,0.2,0.,ERROR,*9999)
        CONTINUE=.TRUE.
        DO WHILE (CONTINUE)
          CALL EVENT(ID_WS,ID_DEVICE,ID_STATUS,CLASS,IDATA,R4DATA,SDATA,
     '      ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,*) ' ID_WS=',ID_WS
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(CLASS(1:6).EQ.'CHOICE') THEN
            INPUT_CHOICE=IDATA(1)
            IF(DOP) THEN
              WRITE(OP_STRING,*) ' INPUT_CHOICE=',INPUT_CHOICE
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF

            IF(ID_WS.EQ.7) THEN
              CALL INPUT_MODE(72,1,'CHOICE','REQUEST',ERROR,*9999)
              IF(ABBREV(OPTION(Input_Choice),'ENERGY',2)) THEN
                OPTION1( 1)='1 Btu'
                OPTION1( 2)='1 cal'
                OPTION1( 3)='1 erg'
                OPTION1( 4)='1 eV'
                OPTION1( 5)='1 ft lbf'
                OPTION1( 6)='1 J (N.m)'
                CALL CHOICE('UNITS',1,1,INSTAT,71,'EVENT',6,6,noch,
     '            noco,5,CO,OPTION1,STRING,0.25,0.,ERROR,*9999)
              ELSE IF(ABBREV(OPTION(Input_Choice),'FORCE',2)) THEN
                OPTION1( 1)='1 dyne'
                OPTION1( 1)='1 kgf'
                OPTION1( 2)='1 lbf'
                OPTION1( 3)='1 N'
                CALL CHOICE('UNITS',1,1,INSTAT,71,'EVENT',4,4,noch,
     '            noco,5,CO,OPTION1,STRING,0.25,0.,ERROR,*9999)
              ELSE IF(ABBREV(OPTION(Input_Choice),'LENGTH',2)) THEN
                OPTION1( 1)='1 ft'
                OPTION1( 2)='1 in'
                OPTION1( 3)='1 m'
                OPTION1( 4)='1 mile'
                CALL CHOICE('UNITS',1,1,INSTAT,71,'EVENT',4,4,noch,
     '            noco,5,CO,OPTION1,STRING,0.25,0.,ERROR,*9999)
              ELSE IF(ABBREV(OPTION(Input_Choice),'MASS',2)) THEN
                OPTION1( 1)='1 kg mass'
                OPTION1( 2)='1 lb mass'
                OPTION1( 3)='1 long ton'
                OPTION1( 4)='1 short ton'
                OPTION1( 5)='1 tonne'
                CALL CHOICE('UNITS',1,1,INSTAT,71,'EVENT',5,5,noch,
     '            noco,5,CO,OPTION1,STRING,0.25,0.,ERROR,*9999)
              ELSE IF(ABBREV(OPTION(Input_Choice),'POWER',2)) THEN
                OPTION1( 1)='1 erg/s'
                OPTION1( 2)='1 ft lbf/s'
                OPTION1( 3)='1 hp'
                OPTION1( 4)='1 kW (kJ/s)'
                CALL CHOICE('UNITS',1,1,INSTAT,71,'EVENT',4,4,noch,
     '            noco,5,CO,OPTION1,STRING,0.25,0.,ERROR,*9999)
              ELSE IF(ABBREV(OPTION(Input_Choice),'SPECIFIC HEAT',2)) 
     '          THEN
                OPTION1( 1)='1 Btu/lb/degF'
                OPTION1( 2)='1 cal/g/degC'
                CALL CHOICE('UNITS',1,1,INSTAT,71,'EVENT',2,2,noch,
     '            noco,5,
     '            CO,OPTION1,STRING,0.25,0.,ERROR,*9999)
              ELSE IF(ABBREV(OPTION(Input_Choice),'STRESS',2)) THEN
                OPTION1( 1)='1 atm (bar)'
                OPTION1( 2)='1 dyn/cm^2'
                OPTION1( 3)='1 kgf/mm^2'
                OPTION1( 4)='1 kPa'
                OPTION1( 5)='1 mm Hg (torr)'
                OPTION1( 6)='1 psi'
                CALL CHOICE('UNITS',1,1,INSTAT,71,'EVENT',6,6,noch,
     '            noco,5,CO,OPTION1,STRING,0.25,0.,ERROR,*9999)
              ELSE IF(ABBREV(OPTION(Input_Choice),'VISCOSITY',2)) THEN
                OPTION1( 1)='1 lb ft.s'
                OPTION1( 2)='1 poise'
                OPTION1( 3)='1 N.s/m^2'
                CALL CHOICE('UNITS',1,1,INSTAT,71,'EVENT',3,3,noch,
     '            noco,5,CO,OPTION1,STRING,0.25,0.,ERROR,*9999)
              ELSE IF(ABBREV(OPTION(Input_Choice),'RETURN',2)) THEN
                CALL INPUT_MODE(72,1,'CHOICE','REQUEST',ERROR,*9999)
                CALL INPUT_MODE(71,1,'CHOICE','REQUEST',ERROR,*9999)
                CALL INPUT_MODE(7,1,'CHOICE','REQUEST',ERROR,*9999)
                CONTINUE=.FALSE.
              ENDIF

            ELSE IF(ID_WS.EQ.71) THEN
              IF(OPTION1(Input_Choice)(1:5).EQ.'1 Btu') THEN
                OPTION2( 1)='= 2.52e+2 cal'
                OPTION2( 2)='= 1.06e+10 erg'
                OPTION2( 3)='= 6.59e+21 eV'
                OPTION2( 4)='= 7.78e+2 ft lbf'
                OPTION2( 5)='= 1.06 J'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',5,5,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:5).EQ.'1 cal') THEN
                OPTION2( 1)='= 3.97e-3 Btu'
                OPTION2( 2)='= 4.19e+7 erg'
                OPTION2( 3)='= 2.61e+19 eV'
                OPTION2( 4)='= 3.09 ft lbf'
                OPTION2( 5)='= 4.19 J'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',5,5,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:5).EQ.'1 erg') THEN
                OPTION2( 1)='= 9.48e-11 Btu'
                OPTION2( 2)='= 2.39e-8 cal'
                OPTION2( 3)='= 6.24e+11 eV'
                OPTION2( 4)='= 7.38e-8 ft lbf'
                OPTION2( 5)='= 1e-7 J'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',5,5,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:4).EQ.'1 eV') THEN
                OPTION2( 1)='= 1.52e-22 Btu'
                OPTION2( 2)='= 3.83e-20 cal'
                OPTION2( 3)='= 1.60e-12 erg'
                OPTION2( 4)='= 1.18e-19 ft lbf'
                OPTION2( 5)='= 1.60e-19 J'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',5,5,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:8).EQ.'1 ft lbf') THEN
                OPTION2( 1)='= 1.29e-3 Btu'
                OPTION2( 2)='= 0.324 cal'
                OPTION2( 3)='= 1.36e+7 erg'
                OPTION2( 4)='= 8.46e+18 eV'
                OPTION2( 5)='= 1.36 J'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',5,5,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:3).EQ.'1 J') THEN
                OPTION2( 1)='= 9.48e-4 Btu'
                OPTION2( 2)='= 0.239 cal'
                OPTION2( 3)='= e+7 erg'
                OPTION2( 4)='= 6.24e+18 eV'
                OPTION2( 5)='= 0.738 ft lbf'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',5,5,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)

              ELSE IF(OPTION1(Input_Choice)(1:6).EQ.'1 dyne') THEN
                OPTION2( 1)='=  kgf'
                OPTION2( 2)='=  lbf'
                OPTION2( 3)='= 1e-5 N'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',3,3,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:5).EQ.'1 kgf') THEN
                OPTION2( 1)='=  dyne'
                OPTION2( 2)='= 2.205 lbf'
                OPTION2( 3)='= 9.807 N'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',3,3,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:5).EQ.'1 lbf') THEN
                OPTION2( 1)='=  dyne'
                OPTION2( 2)='= 0.4535 kgf'
                OPTION2( 3)='= 4.448 N'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',3,3,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:3).EQ.'1 N') THEN
                OPTION2( 1)='= 1e+5 dyne'
                OPTION2( 2)='= 0.1020 kgf'
                OPTION2( 3)='= 0.2248 lbf'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',3,3,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)

              ELSE IF(OPTION1(Input_Choice)(1:4).EQ.'1 ft') THEN
                OPTION2( 1)='= 12 in'
                OPTION2( 2)='= 0.3048 m'
                OPTION2( 3)='= mile'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',3,3,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:4).EQ.'1 in') THEN
                OPTION2( 1)='= 0.0833 ft'
                OPTION2( 2)='= 0.0254 m'
                OPTION2( 3)='= mile'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',3,3,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:3).EQ.'1 m') THEN
                OPTION2( 1)='= 3.281 ft'
                OPTION2( 2)='= 39.37 in'
                OPTION2( 3)='= mile'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',3,3,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:6).EQ.'1 mile') THEN
                OPTION2( 1)='= 5280 ft'
                OPTION2( 2)='= 0.6336e+5 in'
                OPTION2( 3)='= m'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',3,3,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)

              ELSE IF(OPTION1(Input_Choice)(1:9).EQ.'1 kg mass') THEN
                OPTION2( 1)='= 2.2046 kg'
                OPTION2( 2)='= long ton'
                OPTION2( 3)='= short ton'
                OPTION2( 4)='= tonne'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',4,4,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)

              ELSE IF(OPTION1(Input_Choice)(1:9).EQ.'1 lb mass') THEN
                OPTION2( 1)='= 0.4536 kg'
                OPTION2( 2)='= long ton'
                OPTION2( 3)='= short ton'
                OPTION2( 4)='= tonne'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',4,4,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:10).EQ.'1 long ton') THEN
                OPTION2( 1)='= kg'
                OPTION2( 2)='= lb mass'
                OPTION2( 3)='= short ton'
                OPTION2( 4)='= tonne'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',4,4,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:11).EQ.'1 short ton') THEN
                OPTION2( 1)='= kg'
                OPTION2( 2)='= lb mass'
                OPTION2( 3)='= long ton'
                OPTION2( 4)='= tonne'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',4,4,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:7).EQ.'1 tonne') THEN
                OPTION2( 1)='= kg'
                OPTION2( 2)='= lb mass'
                OPTION2( 3)='= long ton'
                OPTION2( 4)='= short ton'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',4,4,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)

              ELSE IF(OPTION1(Input_Choice)(1:7).EQ.'1 erg/s') THEN
                OPTION2( 1)='= 7.38e-8 ft lbf/s'
                OPTION2( 2)='= 1.34e-10 hp'
                OPTION2( 3)='= 1.e-10 kW (kJ/s)'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',3,3,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:10).EQ.'1 ft lbf/s') THEN
                OPTION2( 1)='= 1.36e+7 erg/s'
                OPTION2( 2)='= 1.82e-3 hp'
                OPTION2( 3)='= 1.36e-3 kW (kJ/s)'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',3,3,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:4).EQ.'1 hp') THEN
                OPTION2( 1)='= 7.46e+9 erg/s'
                OPTION2( 2)='= 5.50e+2 ft lbf/s'
                OPTION2( 3)='= 0.7457 kW (kJ/s)'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',3,3,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:4).EQ.'1 kW') THEN
                OPTION2( 1)='= 1.e+10 erg/s'
                OPTION2( 2)='= 7.38e+2 ft lbf/s'
                OPTION2( 3)='= 1.341 hp'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',3,3,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)

              ELSE IF(OPTION1(Input_Choice)(1:13).EQ.'1 Btu/lb/degF') 
     '          THEN
                OPTION2( 1)='=  cal/g/degC'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',1,1,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:12).EQ.'1 cal/g/degC') 
     '          THEN
                OPTION2( 1)='=  Btu/lb/degF'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',1,1,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:13).EQ.'1 Btu/lb/degF') 
     '          THEN

              ELSE IF(OPTION1(Input_Choice)(1:11).EQ.'1 atm (bar)') 
     '          THEN
                OPTION2( 1)='= 1.e+6 dyn/cm^2'
                OPTION2( 2)='= 1.02e-2 kgf/mm^2'
                OPTION2( 3)='= 101.325 kPa'
                OPTION2( 4)='=  mm Hg (torr)'
                OPTION2( 5)='= 14.48 psi'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',5,5,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:10).EQ.'1 dyn/cm^2') THEN
                OPTION2( 1)='= 1.e-6 atm (bar)'
                OPTION2( 2)='= 1.02e-8 kgf/mm^2'
                OPTION2( 3)='=  kPa'
                OPTION2( 4)='=  mm Hg (torr)'
                OPTION2( 5)='= 1.45e-5 psi'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',5,5,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:10).EQ.'1 kgf/mm^2') THEN
                OPTION2( 1)='= 98.1 atm (bar)'
                OPTION2( 2)='=  dyn/cm^2'
                OPTION2( 3)='=  kPa'
                OPTION2( 4)='=  mm Hg (torr)'
                OPTION2( 5)='= 1.42e+3 psi'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',5,5,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,99999)
              ELSE IF(OPTION1(Input_Choice)(1:5).EQ.'1 kPa') THEN
                OPTION2( 1)='= 0.9869E-02 atm (bar)'
                OPTION2( 2)='=  dyn/cm^2'
                OPTION2( 3)='=  kgf/mm^2'
                OPTION2( 4)='=  mm Hg (torr)'
                OPTION2( 5)='= 0.1450 psi'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',5,5,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:14).EQ.'1 mm Hg (torr)')
     '          THEN
                OPTION2( 1)='=  atm (bar)'
                OPTION2( 2)='=  dyn/cm^2'
                OPTION2( 3)='=  kgf/mm^2'
                OPTION2( 4)='=  kPa'
                OPTION2( 5)='=  psi'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',5,5,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:5).EQ.'1 psi') THEN
                OPTION2( 1)='= 6.89e-2 atm (bar)'
                OPTION2( 2)='=  dyn/cm^2'
                OPTION2( 3)='= 7.03e-4 kgf/mm^2'
                OPTION2( 4)='= 6.895 kPa'
                OPTION2( 5)='=  mm Hg (torr)'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',5,5,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)

              ELSE IF(OPTION1(Input_Choice)(1:9).EQ.'1 lb ft.s') THEN
                OPTION2( 1)='=  poise'
                OPTION2( 2)='= 0.1517 N.s/m^2'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',2,2,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)
              ELSE IF(OPTION1(Input_Choice)(1:7).EQ.'1 poise') THEN
                OPTION2( 1)='=  lb ft.s'
                OPTION2( 2)='= 0.1 N.s/m^2'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',2,2,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,99999)
              ELSE IF(OPTION1(Input_Choice)(1:9).EQ.'1 N.s/m^2') THEN
                OPTION2( 1)='= 6.592 lb ft.s'
                OPTION2( 2)='= 10 poise'
                CALL CHOICE('UNITS',1,1,INSTAT,72,'EVENT',2,2,noch,
     '            noco,5,CO,OPTION2,STRING,0.3,0.,ERROR,*9999)

              ELSE IF(ABBREV(OPTION(Input_Choice),'RETURN',2)) THEN
                CALL INPUT_MODE(72,1,'CHOICE','REQUEST',ERROR,*9999)
                CALL INPUT_MODE(71,1,'CHOICE','REQUEST',ERROR,*9999)
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('UNITS')
      RETURN
 9999 CALL ERRORS('UNITS',ERROR)
      CALL EXITS('UNITS')
      RETURN 1
      END

      SUBROUTINE ZERO_SPARSE_MATRIX(ISC_MATRIX,ISR_MATRIX,
     '  MATRIX_M,MATRIX_N,ROW,SPARSITY_TYPE,MATRIX,ERROR,*)
	
C#### Subroutine: ZERO_SPARSE_MATRIX
C###  Description:
C###    ZERO_SPARSE_MATRIX initialises a matrix to all zeros, taking 
C###    into account how the matrix is stored (ie sparsely).
C###    If ROW=0 then entire matrix is zeroed, else
C###    just the particular row.
C     MATRIX_M = Number of rows in the matrix
C              = NYT(1,nc,nx) - GK,GQ
C              = NOT(1,1,nr,nx) - GKK
C     MATRIX_N = Number of columns in the matrix
C              = NYT(2,nc,nx) - GK,GQ
C              = NOT(2,1,nr,nx) - GKK

      IMPLICIT NONE
!     Parameter List
      INTEGER ISC_MATRIX(*),ISR_MATRIX(*),
     '  MATRIX_M,MATRIX_N,ROW,SPARSITY_TYPE
      REAL*8 MATRIX(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER START,FINISH,nz

      CALL ENTERS('ZERO_SPARSE_MATRIX',*9999)

      IF(SPARSITY_TYPE.EQ.1) THEN
        IF(ROW.EQ.0) THEN
          START = ISC_MATRIX(1)-(MATRIX_M+1)
          FINISH = ISC_MATRIX(MATRIX_M+1)-1-(MATRIX_M+1)
        ELSE
          START = ISC_MATRIX(ROW)-(MATRIX_M+1)
          FINISH = ISC_MATRIX(ROW+1)-1-(MATRIX_M+1)
        ENDIF
C GMH 14/9/95 Do we use a BLAS call for this?  
        DO nz=START,FINISH
          MATRIX(nz)=0.0D0
        ENDDO
      ELSE IF(SPARSITY_TYPE.EQ.2) THEN
        ERROR='>>Not implemented yet'
        GOTO 9999
      ELSE IF(SPARSITY_TYPE.EQ.0) THEN
        IF(ROW.EQ.0) THEN
          START = 1
          FINISH = MATRIX_M*MATRIX_N
C GMH 14/9/95 Do we use a BLAS call for this?  
          DO nz=START,FINISH
            MATRIX(nz)=0.0D0
          ENDDO
        ELSE
          START = ROW
          FINISH = ROW+(MATRIX_N-1)*MATRIX_M
C         Consecutive numbers are in a column, we want to do
C           one value per column
          DO nz=START,FINISH,MATRIX_M
            MATRIX(nz)=0.0D0
          ENDDO
        ENDIF
      ELSE
        ERROR='>>Invalid sparsity type'
        GOTO 9999
      ENDIF

      CALL EXITS('ZERO_SPARSE_MATRIX')
      RETURN
 9999 CALL ERRORS('ZERO_SPARSE_MATRIX',ERROR)
      CALL EXITS('ZERO_SPARSE_MATRIX')
      RETURN 1
      END
  
Module FE02
===========

C 21/2/97 LC archived section from routine :
C
C#### Subroutine: DXIDNU
C###  Description:
C###    DXIDNU evaluates derivatives (DXIXN) of Xi- wrt
C###    undeformed Nu(fibre)-coords, and their inverse DXNXI.
C###    This routine assumes that GL and GU are the undeformed metric 
C###    tensors of the Xi coordinate system and that XG contains 
C###    the Gauss pt position, interpolated material axis
C###    orientations, and derivatives with respect to Xi coordinates.

C old way of handling fibre/sheets
C**** 13Sep88: Nu(1) lies in the Xi(1)-Xi(2) plane at angle eta (fibre
C****          angle) to Xi(1). If JTYP9=0 fibre angle is 0 by default.
C****          Nu(2) is orthog to the fibre coord & lies in the Xi(1)-
C****          Xi(2) plane & the remaining Nu coord is orthog to
C****          this plane. The fibre angle eta is eta(1), the angle
C****          between the Nu(1) (fibre) coord & the Xi(1) axis.
C****          The base vectors of the Nu coords have unit length.
C****          Whether the derivs are of deformed or undeformed nu-
C****          coords is determined by the metric tensors with respect
C****          to Xi (GL and GU).
C      REAL*8 C1,C2,DETERM,DXDNU(3,3),DXDT(3,3),ETA1,ETA2,
C     '  RD,RDSQ,RGU33,RK1,RK2
C      IF(NITB.EQ.1) THEN
C        DXINU(1,1)=DSQRT(GL(1,1)) !arclength along 1D element
C
C      ELSE IF(NITB.EQ.2) THEN
C        IF(JTYP9.EQ.0) THEN
C          ETA1=0.0D0
C          ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))
C          C1=1.0D0
C          C2=DCOS(ETA2)
C        ELSE IF(JTYP9.GE.1) THEN !fibre angle defined
C          IF(JTYP12.EQ.1) THEN !fibre angle defined wrt Xi(1) coord
C            ETA1=XG(NJ_LOC(NJL_FIBR,1),1)
C            ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))-ETA1
C          ELSE IF(JTYP12.EQ.2) THEN !fibre angle defined wrt Xi(2) coord
C            ETA2=-XG(NJ_LOC(NJL_FIBR,1),1)
C            ETA1=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))+ETA2
C          ENDIF
C          C1=DCOS(ETA1)
C          C2=DCOS(ETA2)
C        ENDIF
C        RK1=DSQRT(GL(1,1))*C1
C        RK2=DSQRT(GL(2,2))*C2
C        RDSQ=(GL(1,1)*GL(2,2)-GL(1,2)*GL(1,2))
C        RD=DSQRT(RDSQ)
C        RGU33=DSQRT(GU(3,3))
C        DXINU(1,1)=(GL(2,2)*RK1-GL(1,2)*RK2)/RDSQ
C        DXINU(2,1)=(GL(1,1)*RK2-GL(1,2)*RK1)/RDSQ
C        DXINU(3,1)= 0.0D0
C        DXINU(1,2)=-RK2/RD
C        DXINU(2,2)= RK1/RD
C        DXINU(3,2)= 0.0D0
C        DXINU(1,3)= 0.0D0
C        DXINU(2,3)= 0.0D0
C        DXINU(3,3)= RGU33
C
C      ELSE IF(NITB.EQ.3) THEN
C        IF(JTYP9.EQ.0) THEN
C          ETA1=0.0D0
C          ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))
C          C1=1.0D0
C          C2=DCOS(ETA2)
C        ELSE IF(JTYP9.GE.1) THEN !fibre angle defined
C          IF(JTYP12.EQ.1) THEN !fibre angle defined wrt Xi(1) coord
C            ETA1=XG(NJ_LOC(NJL_FIBR,1),1)
C            ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))-ETA1
C          ELSE IF(JTYP12.EQ.2) THEN !fibre angle defined wrt Xi(2) coord
C            ETA2=-XG(NJ_LOC(NJL_FIBR,1),1)
C            ETA1=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))+ETA2
C          ENDIF
C          C1=DCOS(ETA1)
C          C2=DCOS(ETA2)
C        ENDIF
C        RK1=DSQRT(GL(1,1))*C1
C        RK2=DSQRT(GL(2,2))*C2
C        RDSQ=(GL(1,1)*GL(2,2)-GL(1,2)*GL(1,2))
C        RD=DSQRT(RDSQ)
C        RGU33=DSQRT(GU(3,3))
C        DXINU(1,1)=(GL(2,2)*RK1-GL(1,2)*RK2)/RDSQ
C        DXINU(2,1)=(GL(1,1)*RK2-GL(1,2)*RK1)/RDSQ
C        DXINU(3,1)= 0.0D0
C        DXINU(1,2)=-RK2/RD
C        DXINU(2,2)= RK1/RD
C        DXINU(3,2)= 0.0D0
C        DXINU(1,3)=RGU33*(GL(1,2)*GL(2,3)-GL(2,2)*GL(1,3))/RDSQ
C        DXINU(2,3)=RGU33*(GL(1,2)*GL(1,3)-GL(1,1)*GL(2,3))/RDSQ
C        DXINU(3,3)=RGU33
C      ENDIF
C      CALL INVERT(NITB,DXINU,DNUXI,DETERM)

C 21/2/97 removed section from : 
      SUBROUTINE GAUSS2(IBT,INP,nb,NGAP,PG,WG,XIG,ERROR,*)

C#### Subroutine: GAUSS2
C###  Description:
C###    GAUSS2 defines the Gaussian quadrature coords XIGG and weights 
C###    WG and evaluates the basis function Gauss point array PG for
C###    Simplex/Serendipity/Sector  type  basis function nb.

C gmh 3/8/94  Never reached - GAUSS1 called instead.
C       Use standard Gauss quadrature to integrate over the element
C        LMAX=JIDINT(DSQRT(DBLE(NGAP(1))))
C        ng=0
C        DO L=1,LMAX
C          XI(2)=0.5D0+D(L,LMAX)
C          WEIGHT=W(L,LMAX)
C          DO K=1,LMAX
C            XI(1)=0.5D0+D(K,LMAX)
C            ng=ng+1
C            XIG(1,ng)=XI(1)
C            XIG(2,ng)=XI(2)
C            WG(ng)=WEIGHT*W(K,LMAX)
C            DO nn=1,NNT(nb)
C              DO nu=1,NUT(nb)
C                PG(nn,nu,ng)=PSI5(nu,nn,XI)
C              ENDDO !nu
C            ENDDO !nn
C          ENDDO !K
C        ENDDO !L



C LC 21/2/97 removed section from routine 
C
C#### Subroutine: GKSINI
C###  Description:
C###    GKSINI clears workstations.

C cpb 7/2/96 Initialise segment arrays
C      CALL RESET(TRANS)
C      DO iw=1,2*NJT-3
C        IF(IWKS(iw).GT.0) THEN
C          CALL ACWK(iw,1,ERROR,*9999)
C          IF(ISAXES(iw).NE.0) THEN
C            CALL DELETE_SEGMENT(ISAXES(iw),ISEG,iw,ERROR,*9999)
C          ENDIF
C          IF(ISGAUS(iw).NE.0) THEN
C            CALL DELETE_SEGMENT(ISGAUS(iw),ISEG,iw,ERROR,*9999)
C          ENDIF
C          IF(ISRULE(iw).NE.0) THEN
C            CALL DELETE_SEGMENT(ISRULE(iw),ISEG,iw,ERROR,*9999)
C          ENDIF
C          DO np=1,NPM
C            IF(ISNONO(iw,np).NE.0) THEN
C              CALL DELETE_SEGMENT(ISNONO(iw,np),ISEG,iw,ERROR,*9999)
C            ENDIF
C          ENDDO
C          IF(ISGRID(iw).NE.0) THEN
C            CALL DELETE_SEGMENT(ISGRID(iw),ISEG,iw,ERROR,*9999)
C          ENDIF
C          DO ne=1,NEM
C            IF(ISELNO(iw,ne).NE.0) THEN
C              CALL DELETE_SEGMENT(ISELNO(iw,ne),ISEG,iw,ERROR,*9999)
C            ENDIF
C            IF(ISMATE(iw,ne).NE.0) THEN
C              CALL DELETE_SEGMENT(ISMATE(iw,ne),ISEG,iw,ERROR,*9999)
C            ENDIF
C          ENDDO
C          DO nr=1,NRM
C            IF(ISLINE(iw,nr).NE.0) THEN
C              CALL DELETE_SEGMENT(ISLINE(iw,nr),ISEG,iw,ERROR,*9999)
C            ENDIF
C            IF(ISSCAL(iw,nr).NE.0) THEN
C              CALL DELETE_SEGMENT(ISSCAL(iw,nr),ISEG,iw,ERROR,*9999)
C            ENDIF
C            IF(ISSURF(iw,nr).NE.0) THEN
C              CALL DELETE_SEGMENT(ISSURF(iw,nr),ISEG,iw,ERROR,*9999)
C            ENDIF
C          ENDDO
C          IF(ISLINO(iw).NE.0) THEN
C            CALL DELETE_SEGMENT(ISLINO(iw),ISEG,iw,ERROR,*9999)
C          ENDIF
C          DO ne=1,NEM
C            DO nofibr=1,NRM
C              IF(ISFIBR(iw,ne,nofibr).NE.0) THEN
C                CALL DELETE_SEGMENT(ISFIBR(iw,ne,nofibr),ISEG,iw,ERROR,
C     '            *9999)
C              ENDIF
C            ENDDO
C          ENDDO
C          CALL DAWK(iw,1,ERROR,*9999)
C        ENDIF
C      ENDDO


C LC 21/2/97 archived section from routine :
C
C#### Subroutine: LINSEG
C###  Description:
C###    LINSEG defines the various parameters associated with 
C###    global line segments nl=1,NLT, and element line segments 
C###    NAE=1,NLE(nb).
C
C***      Lagrange/Hermite
C
c cpb 7/8/95 Rewritting this code
C          NAE=0
C         Loop thru the edges for each Xi direction. 
C         ie when NI1=1 the Xi(1)-direction edges are generated
C                       and MI(1)=1, MI(2)=2, MI(3)=3
C            when NI1=2 the Xi(2)-direction edges are generated
C                       and MI(1)=2, MI(2)=3, MI(3)=1
C            when NI1=3 the Xi(3)-direction edges are generated
C                       and MI(1)=3, MI(2)=1, MI(3)=2
C         Note: M(ni) records the position index for an edge (or node 
C           on the edge) corresponding to Xi(ni).
C         ie when NI1=1 M(1) gives the node posit along Xi(1)
C                       M(2)   "    "  edge posit rel to the Xi(2) dir.
C                       M(3)   "    "  edge posit rel to the Xi(3) dir.
C          IF(NITB.EQ.1) THEN
C            NAE=NAE+1
C            DO n1=1,IBT(2,1,nb)+1
C              DO nn=1,NNT(nb)
C                IF(INP(nn,1,nb).EQ.n1) NNL(n1,NAE,nb)=nn
C              ENDDO !nn
C            ENDDO !n1
C          ELSE
C            DO ni1=1,NITB
C             Define dirs MI(1),MI(2) & MI(3) cyclically shiftd from NI1
C             note that MI(3) must be treated different for the 2D case 
C              MI(1)=ni1
C              MI(2)=MJ(ni1+1,NITB-1)
C              IF(NITB.EQ.2) THEN
C                MI(3)=3
C                M3=0
C                MM3=1
C              ELSE IF(NITB.EQ.3) THEN
C                MI(3)=MJ(ni1+2,NITB-1)
C                M3=IBT(2,MI(3),nb)
C                MM3=M3
C              ENDIF
C              IF(DOP) THEN
C                WRITE(OP_STRING,'('' MI(1)='',I1,'' MI(2)='',I1,'
C     '            //''' MI(3)='',I1)') MI(1),MI(2),MI(3)
C                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C              ENDIF
C             Pick up the edges in the MI(3) (=cyclic NI1+2) direction
C              DO n3=1,M3+1,MM3
C                M(MI(3))=n3
C               Pick up the edges in the MI(2) (=cyclic NI1+1) direction
C                DO n2=1,IBT(2,MI(2),nb)+1,IBT(2,MI(2),nb)
C                  M(MI(2))=n2
C                  NAE=NAE+1
C                 Pick up the nodes along the MI(1) (=NI1) edge
C                  DO n1=1,IBT(2,MI(1),nb)+1
C                    IF(DOP) THEN
C                      WRITE(OP_STRING,'(/'' N1='',I1,'' N2='',I1,'
C     '                  //''' N3='',I1)') n1,n2,n3
C                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                    ENDIF
C                    M(MI(1))=n1
C                    nn=1
C                    FOUND=.FALSE.
C                    DO WHILE(nn.LE.NNT(nb).AND..NOT.FOUND)
C                      CONT=.TRUE.
C                      DO ni=1,NITB
C                        IF(INP(nn,ni,nb).NE.M(ni)) CONT=.FALSE.
C                      ENDDO !ni
C                      IF(CONT) THEN
C                        NNL(n1,NAE,nb)=nn
C                        FOUND=.TRUE.
C                        IF(DOP) THEN
C                          WRITE(OP_STRING,'('' For basis type nb='','
C     '                      //'I2,'' element edge NAE='',I2,'' has '
C     '                      //'element node NNL('',I1,'',NAE,nb)='','
C     '                      //'I2)') nb,NAE,n1,nn
C                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                        ENDIF
C                      ELSE
C                        nn=nn+1
C                      ENDIF
C                    ENDDO !nn
C                  ENDDO !n1
C                ENDDO !n2
C              ENDDO !n3
C            ENDDO !ni1
C          ENDIF
C          NLE(nb)=NAE

C LC 21/2/97 archived section from :
C
C#### Subroutine: MAT_VEC_NG
C###  Description:
C###    MAT_VEC_NG calculates direction cosines of undeformed 
C###    material vectors at Gauss point ng.
C###    This routine assumes XG contains Gauss pt coordinates, 
C###    microstructural orientations, and derivatives wrt Xi.
C
C old MPN 25-Apr-96
CC**** Delta is angle between h (wall normal) and v (normal to wall in
CC**** yz-plane).
CC**** Phi is angle between v and r.
CC**** Note: If Beta>0 A_VECTOR is calculated with Beta terms included.
C      REAL*8 ALFA,ARG1,ARG2,BETA,GAMA,cosh_LAMDA,
C     '  cos_ALFA,cos_BETA,cos_GAMA,cos_DELTA,cos_PHI,cos_MU,
C     '  cos_THETA,coth_LAMDA,cot_MU,
C     '  D1,D2,D3,DELTA,dLAMDA_dXi1,dLAMDA_dXi2,dMU_dXi1,dMU_dXi2,
C     '  dTHETA_dXi1,dX_dXI1(3),PHI,RLAMDA,RMU,sinh_LAMDA,
C     '  sin_ALFA,sin_BETA,sin_GAMA,sin_DELTA,sin_PHI,sin_MU,sin_THETA,
C     '  SUMSQ,THETA,tanh_LAMDA,tan_MU
C      IF(NITB.EQ.1) THEN    !1d
C      	IF(ITYP10(1).EQ.1) THEN  !r.c. coordinates
C          DO nj=1,3
C      	    dX_dXI1(nj)=XG(nj,2)
C          ENDDO
C        ELSE IF(ITYP10(1).EQ.2) THEN !cylindrical polar coordinates
C          THETA=XG(2,1)
C          sin_THETA=DSIN(THETA)
C          cos_THETA=DCOS(THETA)
C          dX_dXI1(1)=XG(1,2)*cos_THETA - XG(1,1)*sin_THETA*XG(2,2)
C          dX_dXI1(2)=XG(1,2)*sin_THETA + XG(1,1)*cos_THETA*XG(2,2)
C          dX_dXI1(3)=XG(3,2)
C	ENDIF
C      	SUMSQ=0.0D0
C      	DO nj=1,3
C      	  SUMSQ=SUMSQ+dX_dXI1(nj)*dX_dXI1(nj)
C      	ENDDO
C      	SUMSQ=DSQRT(SUMSQ)
C      	DO nj=1,3
C      	  IF(SUMSQ.GT.1.D-12) A_VECTOR(nj)=dX_dXI1(nj)/SUMSQ
C      	  B_VECTOR(nj)=0.0d0
C      	  C_VECTOR(nj)=0.0d0
C      	ENDDO
C      ELSE IF(NITB.EQ.2) THEN   !2D
C      	IF(ITYP10(1).EQ.1) THEN  !r.c. coordinates only so far
C      	  nj=NJ_LOC(NJL_FIBR,1)  !Fibre angle position
C         ALFA=XG(nj,1)          !Fibre angle
C         dX_dXI1(1)=XG(1,2)
C         dX_dXI1(2)=XG(2,2)
C         PHI=DATAN2(dX_dXI1(2),dX_dXI1(1)) !Angle betw x and xi
C	  A_VECTOR(1)=DCOS(ALFA+PHI)
C	  A_VECTOR(2)=DSIN(ALFA+PHI)
C	  A_VECTOR(3)=0.0D0
C	  B_VECTOR(1)=-DSIN(ALFA+PHI)
C	  B_VECTOR(2)=DCOS(ALFA+PHI)
C	  B_VECTOR(3)=0.0D0
C	  C_VECTOR(1)=0.0D0
C	  C_VECTOR(2)=0.0D0
C	  C_VECTOR(3)=0.0D0
C        ENDIF
C      ELSE IF(NITB.EQ.3) THEN   !3D
C      	IF(ITYP10(1).EQ.1) THEN  !r.c. coordinates
C
C!  Alfa (fibre angle)  : rotation in x,y plane about cross axis
C      	  nj=NJ_LOC(NJL_FIBR,1)
C          IF(nj.EQ.0) THEN
C            WRITE(OP_STRING,
C     '        '('' >>WARNING !!! Assuming fibre angle of zero'')')
C            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C            ALFA       = 0.0D0
C          ELSE
C            ALFA       = XG(nj,1)
C          ENDIF 
C      	  IF(DOP) THEN
C      	    WRITE(OP_STRING,'(''   ALFA ='',E12.3)') ALFA
C      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      	  ENDIF
C      	  sin_ALFA   = DSIN(ALFA)
C      	  cos_ALFA   = DCOS(ALFA)
C
C!  Beta (imbrication angle)  : rotation of fibre axis about sheet axis
C      	  nj=NJ_LOC(NJL_FIBR,2) 
C          IF(nj.EQ.0) THEN
C            WRITE(OP_STRING,
C     '        '('' >>WARNING !!! Assuming imbrication angle of zero'')')
C            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C            BETA       = 0.0D0
C          ELSE
C            BETA       = XG(nj,1)
C          ENDIF 
C      	  IF(DOP) THEN
C      	    WRITE(OP_STRING,'(''   BETA ='',E12.3)') BETA
C      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      	  ENDIF
C      	  sin_BETA   = DSIN(BETA)
C      	  cos_BETA   = DCOS(BETA)
C
C!  Gamma (sheet angle)  : rotation of sheet axis about fibre axis
C      	  nj=NJ_LOC(NJL_FIBR,3) 
C          IF(nj.EQ.0) THEN
C            WRITE(OP_STRING,
C     '        '('' >>WARNING !!! Assuming sheet angle of zero'')')
C            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
C            GAMA       = 0.0D0
C          ELSE
C            GAMA       = XG(nj,1)
C          ENDIF 
C      	  IF(DOP) THEN
C      	    WRITE(OP_STRING,'(''  GAMMA ='',E12.3)') GAMA
C      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      	  ENDIF
C      	  sin_GAMA   = DSIN(GAMA)
C      	  cos_GAMA   = DCOS(GAMA)
CC GBS Removed 22-MAR-1996
Cc          WRITE(OP_STRING,
Cc     '      '('' >>WARNING !!! Glen has hacked this to be orthog'')')
Cc          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
Cc      	  sin_GAMA   = 1.0D0
Cc      	  cos_GAMA   = 0.0D0
C      	  A_VECTOR(1) =  cos_ALFA*cos_BETA + sin_ALFA*sin_BETA*sin_GAMA
C      	  A_VECTOR(2) =  sin_ALFA*cos_BETA - cos_ALFA*sin_GAMA*sin_BETA
C      	  A_VECTOR(3) = -cos_GAMA*sin_BETA
C
C      	  B_VECTOR(1) =  cos_ALFA*sin_BETA - sin_ALFA*sin_GAMA*cos_BETA
C      	  B_VECTOR(2) =  sin_ALFA*sin_BETA + cos_ALFA*sin_GAMA*cos_BETA
C      	  B_VECTOR(3) =  cos_GAMA*cos_BETA
C
C 	  C_VECTOR(1) =  sin_ALFA*cos_GAMA
C  	  C_VECTOR(2) = -cos_ALFA*cos_GAMA
C  	  C_VECTOR(3) =  sin_GAMA
CC GMH this is hacked in so that we have orthogonal axes
CC GBS Removed 22-MAR-1996
Cc          CALL CROSS(A_VECTOR,B_VECTOR,C_VECTOR)
C      
C      	ELSE IF(ITYP10(1).EQ.4) THEN !prolate coordinates
C
C!  Lamda
C      	  RLAMDA     = XG(1,1)
C      	  IF(DOP) THEN
C      	    WRITE(OP_STRING,'('' RLAMDA ='',E12.3)') RLAMDA
C      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      	  ENDIF
C      	  dLAMDA_dXi1= XG(1,2)
C      	  dLAMDA_dXi2= XG(1,4)
C      	  sinh_LAMDA = DSINH(RLAMDA)
C      	  cosh_LAMDA = DCOSH(RLAMDA)
C      	  IF(DABS(cosh_LAMDA).GT.1.0D-6) THEN
C      	    tanh_LAMDA = sinh_LAMDA/cosh_LAMDA
C      	  ELSE
C      	    tanh_LAMDA = 1.0D8
C      	  ENDIF
C      	  IF(DABS(sinh_LAMDA).GT.1.0D-6) THEN
C      	    coth_LAMDA = cosh_LAMDA/sinh_LAMDA
C      	  ELSE
C      	    coth_LAMDA = 1.0D8
C      	  ENDIF
C
C!  Mu	  
C      	  RMU        = XG(2,1)
C      	  IF(DOP) THEN
C      	    WRITE(OP_STRING,'(''    RMU ='',E12.3)') RMU
C      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      	  ENDIF
C      	  dMU_dXi1   = XG(2,2)
C      	  dMU_dXi2   = XG(2,4)
C      	  sin_MU     = DSIN(RMU)
C      	  cos_MU     = DCOS(RMU)
C      	  IF(DABS(cos_MU).GT.1.0D-6) THEN
C      	    tan_MU = sin_MU/cos_MU
C      	  ELSE
C      	    tan_MU = 1.0D8
C      	  ENDIF
C      	  IF(DABS(sin_MU).GT.1.0D-6) THEN
C      	    cot_MU     = cos_MU/sin_MU
C      	  ELSE
C      	    cot_MU = 1.0D8
C      	  ENDIF
C
C!  Theta
C      	  THETA      = XG(3,1)
C      	  IF(DOP) THEN
C      	    WRITE(OP_STRING,'(''  THETA ='',E12.3)') THETA
C      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      	  ENDIF
C      	  dTHETA_dXi1= XG(3,2)
C      	  sin_THETA  = DSIN(THETA)
C      	  cos_THETA  = DCOS(THETA)
C
C!  Alfa (fibre angle)
C      	  nj=NJ_LOC(NJL_FIBR,1) 
C      	  ALFA       = XG(nj,1)
C      	  IF(DOP) THEN
C      	    WRITE(OP_STRING,'(''   ALFA ='',E12.3)') ALFA
C      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      	  ENDIF
C      	  sin_ALFA   = DSIN(ALFA)
C      	  cos_ALFA   = DCOS(ALFA)
C
C!  Beta (imbrication angle)
C      	  nj=NJ_LOC(NJL_FIBR,2) 
C      	  BETA       = XG(nj,1)
C      	  IF(DOP) THEN
C      	    WRITE(OP_STRING,'(''   BETA ='',E12.3)') BETA
C      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      	  ENDIF
C      	  sin_BETA   = DSIN(BETA)
C      	  cos_BETA   = DCOS(BETA)
C
C!  Gamma (sheet angle)
C      	  nj=NJ_LOC(NJL_FIBR,3) 
C      	  GAMA       = XG(nj,1)
C      	  IF(DOP) THEN
C      	    WRITE(OP_STRING,'(''   GAMA ='',E12.3)') GAMA
C      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      	  ENDIF
C      	  sin_GAMA   = DSIN(GAMA)
C      	  cos_GAMA   = DCOS(GAMA)
C
C!  Phi (eqtn 8 in "Laminar structure of the heart II")
C      	  ARG1=coth_LAMDA*dLAMDA_dXi1+cot_MU*dMU_dXi1
C      	  ARG2=dTHETA_dXi1
C      	  PHI =DATAN2(ARG1,ARG2)
C      	  IF(PHI.GT.PI/2.0D0) THEN
C      	    PHI=PHI-PI
C      	  ELSE IF(PHI.LT.-PI/2.0D0) THEN
C      	    PHI=PHI+PI
C      	  ENDIF
C      	  IF(DOP) THEN
C      	    WRITE(OP_STRING,'(''    PHI ='',E12.3)') PHI
C      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      	  ENDIF
C      	  sin_PHI   = DSIN(PHI)
C      	  cos_PHI   = DCOS(PHI)
C
C!  Delta (eqtn 12 in "Laminar structure of the heart II")
C      	  ARG1=(tan_MU*dLAMDA_dXi2+tanh_LAMDA*dMU_dXi2)*cos_PHI
C      	  ARG2=tanh_LAMDA*dLAMDA_dXi2-tan_MU*dMU_dXi2
C      	  DELTA=DATAN2(ARG1,ARG2)
C      	  IF(DELTA.GT.PI/2.0D0) THEN
C      	    DELTA=DELTA-PI
C      	  ELSE IF(DELTA.LT.-PI/2.0D0) THEN
C      	    DELTA=DELTA+PI
C      	  ENDIF
C      	  IF(DOP) THEN
C      	    WRITE(OP_STRING,'(''  DELTA ='',E12.3)') DELTA
C      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      	  ENDIF
C      	  sin_DELTA  = DSIN(DELTA)
C      	  cos_DELTA  = DCOS(DELTA)
C
C!  A vector (eqtn 3 in "Laminar structure of the heart II")
C      	  IF(DABS(BETA).LT.1.0D-5) THEN !no imbrication angle
C      	    A_VECTOR(1)=-sin_ALFA*cos_DELTA
C      	    A_VECTOR(2)=-cos_ALFA*sin_PHI-sin_ALFA*sin_DELTA*cos_PHI
C      	    A_VECTOR(3)=-cos_ALFA*cos_PHI+sin_ALFA*sin_DELTA*sin_PHI
C      	  ELSE !include imbrication angle
C      	    D1=-(sin_ALFA*cos_BETA-cos_ALFA*sin_BETA*sin_GAMA)*cos_DELTA
C     '	       + sin_BETA*cos_GAMA*sin_DELTA
C      	    D2=  cos_ALFA*cos_BETA+sin_ALFA*sin_BETA*sin_GAMA
C      	    D3=-(sin_ALFA*cos_BETA-cos_ALFA*sin_BETA*sin_GAMA)*sin_DELTA
C     '	       - sin_BETA*cos_GAMA*cos_DELTA
C      	    A_VECTOR(1)= D1
C      	    A_VECTOR(2)=-D2*sin_PHI+D3*cos_PHI
C      	    A_VECTOR(3)=-D2*cos_PHI-D3*sin_PHI
C      	  ENDIF
C      	  CALL ROTATION(1,THETA,A_VECTOR,ERROR,*9999)
C
C!  B vector (eqtn 4 in "Laminar structure of the heart II")
C      	  ARG1=cos_ALFA*sin_GAMA*sin_DELTA-cos_GAMA*cos_DELTA
C      	  B_VECTOR(1)=-cos_ALFA*sin_GAMA*cos_DELTA - cos_GAMA*sin_DELTA
C      	  B_VECTOR(2)= sin_ALFA*sin_GAMA*sin_PHI - ARG1*cos_PHI
C      	  B_VECTOR(3)= sin_ALFA*sin_GAMA*cos_PHI + ARG1*sin_PHI
C      	  CALL ROTATION(1,THETA,B_VECTOR,ERROR,*9999)
C
C!  C vector (eqtn 5 in "Laminar structure of the heart II")
C      	  CALL CROSS(A_VECTOR,B_VECTOR,C_VECTOR)
C      	ENDIF
C      ENDIF



C LC 21/2/97 removed from routine:
C
C#### Subroutine: MAT_VEC_XI
C###  Description:
C###    MAT_VEC_XI calculates direction cosines of material vectors at 
C###    XI in element ne.
C###    This routine assumes that XE contains undeformed 
C###    element vertex coordinates, and microstructural material angles.
C
CC**** Note: XE must be passed in (calculated in previous call to XPXE).
CC**** A_VECTOR is fibre angle vector (in sheet)
CC**** B_VECTOR is sheet angle vector (in sheet orthog to fibres)
CC**** C_VECTOR   is normal to sheet
CC**** Alfa is fibre angle
CC**** Beta is imbrication angle
CC**** Gama is sheet angle
CC**** Delta is angle between h (wall normal) and v (normal to wall in
CC**** yz-plane)
CC**** Phi is angle between v and r
CC**** GAMA,PHI & DELTA are also returned.
CC**** Note: If Beta>0 A_VECTOR is calculated with Beta terms included.
C      REAL*8 ARG1,ARG2,cosh_LAMDA,
C     '  cos_ALFA,cos_BETA,cos_GAMA,cos_DELTA,cos_PHI,cos_MU,
C     '  coth_LAMDA,cot_MU,
C     '  D1,D2,D3,dLAMDA_dXi1,dLAMDA_dXi2,dMU_dXi1,dMU_dXi2,
C     '  dTHETA_dXi1,
C     '  PXI,RLAMDA,RMU,sinh_LAMDA,sin_ALFA,sin_BETA,sin_GAMA,sin_DELTA,
C     '  sin_PHI,sin_MU,THETA,tanh_LAMDA,tan_MU
C      IF(ITYP10(1).EQ.4) THEN !prolate coordinates
C        IF(DOP) THEN
C          WRITE(OP_STRING,'('' Xi coords: '',3E12.3)') (XI(I),I=1,3)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C!  Lamda
C        nb=NBJ(1)
C        RLAMDA     = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C     '    XE(1,1))
C        IF(DOP) THEN
C          WRITE(OP_STRING,'('' RLAMDA ='',E12.3)') RLAMDA
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C        dLAMDA_dXi1= PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,XI,
C     '    XE(1,1))
C        dLAMDA_dXi2= PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,4,XI,
C     '    XE(1,1))
C        sinh_LAMDA = DSINH(RLAMDA)
C        cosh_LAMDA = DCOSH(RLAMDA)
C        IF(DABS(cosh_LAMDA).GT.1.0D-6) THEN
C          tanh_LAMDA = sinh_LAMDA/cosh_LAMDA
C        ELSE
C          tanh_LAMDA = 1.0D8
C        ENDIF
C        IF(DABS(sinh_LAMDA).GT.1.0D-6) THEN
C          coth_LAMDA = cosh_LAMDA/sinh_LAMDA
C        ELSE
C          coth_LAMDA = 1.0D8
C        ENDIF
C
C!  Mu
C        nb=NBJ(2)
C        RMU        = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C     '    XE(1,2))
C        IF(DOP) THEN
C          WRITE(OP_STRING,'(''    RMU ='',E12.3)') RMU
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C        dMU_dXi1   = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,XI,
C     '    XE(1,2))
C        dMU_dXi2   = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,4,XI,
C     '    XE(1,2))
C        sin_MU     = DSIN(RMU)
C        cos_MU     = DCOS(RMU)
C        IF(DABS(cos_MU).GT.1.0D-6) THEN
C          tan_MU = sin_MU/cos_MU
C        ELSE
C          tan_MU = 1.0D8
C        ENDIF
C        IF(DABS(sin_MU).GT.1.0D-6) THEN
C          cot_MU     = cos_MU/sin_MU
C        ELSE
C          cot_MU = 1.0D8
C        ENDIF
C
C!  Theta
C        nb=NBJ(3)
C        THETA = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C     '    XE(1,3))
C        IF(DOP) THEN
C          WRITE(OP_STRING,'(''  THETA ='',E12.3)') THETA
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C        dTHETA_dXi1= PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,XI,
C     '    XE(1,3))
C
C!  Alfa (fibre angle)
C      	nj=NJ_LOC(NJL_FIBR,1) 
C        nb=NBJ(nj)
C        ALFA = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C     '    XE(1,nj))
C        IF(DOP) THEN
C          WRITE(OP_STRING,'(''   ALFA ='',E12.3)') ALFA
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C        sin_ALFA   = DSIN(ALFA)
C        cos_ALFA   = DCOS(ALFA)
C
C!  Beta (imbrication angle)
C      	nj=NJ_LOC(NJL_FIBR,2) 
C        nb=NBJ(nj)
C        IF(nb.GT.0) THEN
C          BETA=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C     '      XE(1,nj))
C        ELSE
C          BETA=0.0D0
C        ENDIF
C        IF(DOP) THEN
C          WRITE(OP_STRING,'(''   BETA ='',E12.3)') BETA
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C        sin_BETA   = DSIN(BETA)
C        cos_BETA   = DCOS(BETA)
C
C!  Gama (sheet angle)
C      	nj=NJ_LOC(NJL_FIBR,3) 
C        nb=NBJ(nj)
C        GAMA = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C     '    XE(1,nj))
C        IF(DOP) THEN
C          WRITE(OP_STRING,'(''   GAMA ='',E12.3)') GAMA
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C        sin_GAMA   = DSIN(GAMA)
C        cos_GAMA   = DCOS(GAMA)
C
C!  Phi  (eqtn 8 in "Laminar structure of the heart II")
C        ARG1=coth_LAMDA*dLAMDA_dXi1+cot_MU*dMU_dXi1
C        ARG2=dTHETA_dXi1
C        PHI =DATAN2(ARG1,ARG2)
C        IF(PHI.GT.PI/2.0D0) THEN
C          PHI=PHI-PI
C        ELSE IF(PHI.LT.-PI/2.0D0) THEN
C          PHI=PHI+PI
C        ENDIF
C        IF(DOP) THEN
C          WRITE(OP_STRING,'(''    PHI ='',E12.3)') PHI
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C        sin_PHI   = DSIN(PHI)
C        cos_PHI   = DCOS(PHI)
C
C!  Delta (eqtn 12 in "Laminar structure of the heart II")
C        ARG1=(tan_MU*dLAMDA_dXi2+tanh_LAMDA*dMU_dXi2)*cos_PHI
C        ARG2=tanh_LAMDA*dLAMDA_dXi2-tan_MU*dMU_dXi2
C        DELTA=DATAN2(ARG1,ARG2)
C        IF(DELTA.GT.PI/2.0D0) THEN
C          DELTA=DELTA-PI
C        ELSE IF(DELTA.LT.-PI/2.0D0) THEN
C          DELTA=DELTA+PI
C        ENDIF
C        IF(DOP) THEN
C          WRITE(OP_STRING,'(''  DELTA ='',E12.3)') DELTA
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C        sin_DELTA  = DSIN(DELTA)
C        cos_DELTA  = DCOS(DELTA)
C
C!  A vector (eqtn 3 in "Laminar structure of the heart II")
C        IF(DABS(BETA).LT.1.0D-5) THEN !no imbrication angle
C          A_VECTOR(1)=-sin_ALFA*cos_DELTA
C          A_VECTOR(2)=-cos_ALFA*sin_PHI-sin_ALFA*sin_DELTA*cos_PHI
C          A_VECTOR(3)=-cos_ALFA*cos_PHI+sin_ALFA*sin_DELTA*sin_PHI
C        ELSE !include imbrication angle
C          D1=-(sin_ALFA*cos_BETA-cos_ALFA*sin_BETA*sin_GAMA)*cos_DELTA
C     '       + sin_BETA*cos_GAMA*sin_DELTA
C          D2=  cos_ALFA*cos_BETA+sin_ALFA*sin_BETA*sin_GAMA
C          D3=-(sin_ALFA*cos_BETA-cos_ALFA*sin_BETA*sin_GAMA)*sin_DELTA
C     '       - sin_BETA*cos_GAMA*cos_DELTA
C          A_VECTOR(1)= D1
C          A_VECTOR(2)=-D2*sin_PHI+D3*cos_PHI
C          A_VECTOR(3)=-D2*cos_PHI-D3*sin_PHI
C        ENDIF
C        CALL ROTATION(1,THETA,A_VECTOR,ERROR,*9999)
C        IF(DOP) THEN
C          WRITE(OP_STRING,'('' a_vector:'',3E12.3)') 
C     '      (A_VECTOR(I),I=1,3)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C!  B vector (eqtn 4 in "Laminar structure of the heart II")
C        ARG1=cos_ALFA*sin_GAMA*sin_DELTA-cos_GAMA*cos_DELTA
C        B_VECTOR(1)=-cos_ALFA*sin_GAMA*cos_DELTA - cos_GAMA*sin_DELTA
C        B_VECTOR(2)= sin_ALFA*sin_GAMA*sin_PHI - ARG1*cos_PHI
C        B_VECTOR(3)= sin_ALFA*sin_GAMA*cos_PHI + ARG1*sin_PHI
C        CALL ROTATION(1,THETA,B_VECTOR,ERROR,*9999)
C        IF(DOP) THEN
C          WRITE(OP_STRING,'('' b_vector:'',3E12.3)') 
C     '      (B_VECTOR(I),I=1,3)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C!  C vector (eqtn 5 in "Laminar structure of the heart II")
C        CALL CROSS(A_VECTOR,B_VECTOR,C_VECTOR)
C        IF(DOP) THEN 
C          WRITE(OP_STRING,'('' c_vector:'',3E12.3)') 
C     '      (C_VECTOR(I),I=1,3)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDIF
C      ENDIF


C LC 21/2/97 removed section from routine :
C
C#### Subroutine: TOFFEL
C###  Description:
C###    TOFFEL evaluates Christoffel symbol of 2nd kind CHTOFF(ia,ib,ic)
C###    and, if CURVE is .true., the curvature tensor DBM(ia,ib,ic)
C###    in Xi-coordinate system at current Gauss point.
C**** X3G(nj,nd) are 3rd derivs of Xj (see subroutine X3XG)
C**** G(nj,1) is tangent to surface in Xi(1) direction
C**** G(nj,2) is    "     "    "     " Xi(2)     "
C**** G(nj,3) is    "     "    "      "Xi(3)     "
C**** DG(nj,ib,ic) are cpts of base vector G(nj,ib) derivs wrt Xi(ic)
C
C old MPN 24-Apr-96: DZX and DZXX functions in FE01 now used
C        IF(ITYP10(nr).EQ.2) THEN !cylindrical polar coords
C          RAD=XG(1,1)
C          THETA=XG(2,1)
C          SINT=DSIN(THETA)
C          COST=DCOS(THETA)
C          DZX(1,1)=COST
C          DZX(1,2)=-RAD*SINT
C          DZX(2,1)=SINT
C          DZX(2,2)=RAD*COST
C          DZX(3,3)=1.0D0
C          DZXX(1,2,1)=-SINT
C          DZXX(2,2,1)=COST
C          DZXX(1,1,2)=-SINT
C          DZXX(1,2,2)=-RAD*COST
C          DZXX(2,1,2)=COST
C          DZXX(2,2,2)=-RAD*SINT
C
C        ELSE IF(ITYP10(nr).EQ.3) THEN !spherical polar coords
C          RAD  =XG(1,1)
C          THETA=XG(2,1)
C          PHI  =XG(3,1)
C          SINT=DSIN(THETA)
C          COST=DCOS(THETA)
C          SINP=DSIN(PHI)
C          COSP=DCOS(PHI)
C          DZX(1,1)=COST*COSP
C          DZX(1,2)=-RAD*SINT*COSP
C          DZX(1,3)=-RAD*COST*SINP
C          DZX(2,1)=SINT*COSP
C          DZX(2,2)=RAD*COST*COSP
C          DZX(2,3)=-RAD*SINT*SINP
C          DZX(3,1)=SINP
C          DZX(3,3)=RAD*COSP
C          DZXX(1,2,1)=-SINT*COSP
C          DZXX(1,3,1)=-COST*SINP
C          DZXX(2,2,1)=COST*COSP
C          DZXX(2,3,1)=-SINT*SINP
C          DZXX(3,3,1)=COSP
C          DZXX(1,1,2)=-SINT*COSP
C          DZXX(1,2,2)=-RAD*COST*COSP
C          DZXX(1,3,2)=RAD*SINT*SINP
C          DZXX(2,1,2)=COST*COSP
C          DZXX(2,2,2)=-RAD*SINT*COSP
C          DZXX(2,3,2)=-RAD*COST*SINP
C          DZXX(1,1,3)=-COST*SINP
C          DZXX(1,2,3)=RAD*SINT*SINP
C          DZXX(1,3,3)=-RAD*COST*COSP
C          DZXX(2,1,3)=-SINT*SINP
C          DZXX(2,2,3)=-RAD*COST*SINP
C          DZXX(2,3,3)=-RAD*SINT*COSP
C          DZXX(3,1,3)=COSP
C          DZXX(3,3,3)=-RAD*SINP
C
C        ELSE IF(ITYP10(nr).EQ.4) THEN !prolate spheroidal coords
C          LAMDA=XG(1,1)
C          mu=XG(2,1)
C          THETA=XG(3,1)
C          SINHL=DSINH(LAMDA)
C          COSHL=DCOSH(LAMDA)
C          SINM=DSIN(MU)
C          COSM=DCOS(MU)
C          SINT=DSIN(THETA)
C          COST=DCOS(THETA)
C          DZX(1,1)=FOCUS*SINHL*COSM
C          DZX(1,2)=-FOCUS*COSHL*SINM
C          DZX(1,3)=0.0D0
C          DZX(2,1)=FOCUS*COSHL*SINM*COST
C          DZX(2,2)=FOCUS*SINHL*COSM*COST
C          DZX(2,3)=-FOCUS*SINHL*SINM*SINT
C          DZX(3,1)=FOCUS*COSHL*SINM*SINT
C          DZX(3,2)=FOCUS*SINHL*COSM*SINT
C          DZX(3,3)=FOCUS*SINHL*SINM*COST
C          DZXX(1,1,1)=FOCUS*COSHL*COSM
C          DZXX(1,2,1)=-FOCUS*SINHL*SINM
C          DZXX(2,1,1)=FOCUS*SINHL*SINM*COST
C          DZXX(2,2,1)=FOCUS*COSHL*COSM*COST
C          DZXX(2,3,1)=-FOCUS*COSHL*SINM*SINT
C          DZXX(3,1,1)=FOCUS*SINHL*SINM*SINT
C          DZXX(3,2,1)=FOCUS*COSHL*COSM*SINT
C          DZXX(3,3,1)=FOCUS*COSHL*SINM*COST
C          DZXX(1,1,2)=-FOCUS*SINHL*SINM
C          DZXX(1,2,2)=-FOCUS*COSHL*COSM
C          DZXX(2,1,2)=FOCUS*COSHL*COSM*COST
C          DZXX(2,2,2)=-FOCUS*SINHL*SINM*COST
C          DZXX(2,3,2)=-FOCUS*SINHL*COSM*SINT
C          DZXX(3,1,2)=FOCUS*COSHL*COSM*SINT
C          DZXX(3,2,2)=-FOCUS*SINHL*SINM*SINT
C          DZXX(3,3,2)=FOCUS*SINHL*COSM*COST
C          DZXX(2,1,3)=-FOCUS*COSHL*SINM*SINT
C          DZXX(2,2,3)=-FOCUS*SINHL*COSM*SINT
C          DZXX(2,3,3)=-FOCUS*SINHL*SINM*COST
C          DZXX(3,1,3)=FOCUS*COSHL*SINM*COST
C          DZXX(3,2,3)=FOCUS*SINHL*COSM*COST
C          DZXX(3,3,3)=-FOCUS*SINHL*SINM*SINT
C
C        ENDIF !ityp10


      SUBROUTINE CALC_IMAGE_DATA(NBJ,NELIST,NKJE,NPF,NPNE,NRE,
     '  NVJE,SE,XA,XE,XP,ZD,ERROR,*)

C#### Subroutine: CALC_IMAGE_DATA
C###  Description:
C###    CALC_IMAGE_DATA defines a reduced data set ZD for an image 
C###    based on which pixels lie inside of, or near, the defined 
C###    quadratic elements.  A relatively "efficient" algorithm is 
C###    used whereby pixels are compared to a bilinear element roughly 
C###    surrounding the actual element.

C**** In its present form, this routine will only work with quadratic
C**** elements.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:pics01.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NELIST(0:NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,NCTR,ne,nolist
      REAL*8 QSCAL,XI1,XI2,XPIX,XSCAL,XX1,XX2,XX3,XX4,YPIX,YSCAL,YY1,
     '  YY2,YY3,YY4
      LOGICAL LINCLUDE(12,12),LMEMBER
c      LOGICAL LINCLUDE(512,512),LMEMBER

      CALL ENTERS('CALC_IMAGE_DATA',*9999)

C ... Mask image to determine which pixels lie inside bilinear
C     elements which completely surround respective curved elements.
C     accordingly, there are NEELEM(0,nr) bilinear elems corresponding
C     to NEELEM(0,nr) curved (probably biquadratic elements).
C     if XPIX,YPIX lies in any of these bilinear elements,
C        set LINCLUDE(I,J)=.TRUE.


      XSCAL=DBLE(XMAX-XMIN)/DBLE(IMGX)   !pixel scaling
      YSCAL=DBLE(YMAX-YMIN)/DBLE(IMGY)

      DO j=1,IMGY
        DO i=1,IMGX
          LINCLUDE(j,i)=.FALSE.
        ENDDO   !end I
      ENDDO     !end J

      QSCAL=1.2d0
      DO nolist=1,NELIST(0)
        ne=NELIST(nolist)
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '    SE(1,1,ne),XA,XE,XP,ERROR,*9999)

C ...   XX1...YY4 are nodes of a bilinear element which surround the
C       the element defined by XE

        XX1=XE(9,1)-QSCAL*(XE(9,1)-XE(1,1))
        XX2=XE(7,1)+QSCAL*(XE(3,1)-XE(7,1))
        XX3=XE(3,1)-QSCAL*(XE(3,1)-XE(7,1))
        XX4=XE(1,1)+QSCAL*(XE(9,1)-XE(1,1))
        YY1=XE(9,2)-QSCAL*(XE(9,2)-XE(1,2))
        YY2=XE(7,2)+QSCAL*(XE(3,2)-XE(7,2))
        YY3=XE(3,2)-QSCAL*(XE(3,2)-XE(7,2))
        YY4=XE(1,2)+QSCAL*(XE(9,2)-XE(1,2))

        DO j=1,IMGY
          DO i=1,IMGX

            XPIX=DBLE(i)*XSCAL
            YPIX=DBLE(IMGY+1-j)*YSCAL
            CALL BIL_XI_FNDR(XX1,XX2,XX3,XX4,YY1,YY2,YY3,YY4,
     '        XPIX,YPIX,XI1,XI2,LMEMBER,ERROR,*9999)

            IF(LMEMBER) LINCLUDE(j,i)=.TRUE.

          ENDDO   !end I
        ENDDO     !end J
      ENDDO       !end NOLIST

      NCTR=0
      DO j=1,IMGY
        DO i=1,IMGX
          IF(LINCLUDE(j,i)) THEN
            NCTR=NCTR+1
            ZD(NCTR,1)=DBLE(i)*XSCAL
            ZD(NCTR,2)=DBLE(IMGY+1-j)*YSCAL
          ENDIF
        ENDDO
      ENDDO
      NDT=NCTR    !number of points held in ZD

      CALL EXITS('CALC_IMAGE_DATA')
      RETURN
 9999 CALL ERRORS('CALC_IMAGE_DATA',ERROR)
      CALL EXITS('CALC_IMAGE_DATA')
      RETURN 1
      END


      SUBROUTINE CALC_STRIPE_PTS(NT_STRIPE_DATA,NDP,WD,ZD,GROUP_NAME,
     '  ERROR,*)

C**** This routine finds the intersections of a group of bezier
C**** polylines and subdivides the resulting "chunks" into points based
C**** on arc length.  Input NT_STRIPE_DATA is the total number of data
C**** points to be generated from the group of polylines.  In the case
C**** where a value of 0 for parameter NT_STRIPE_DATA is supplied,
C**** discretization parameters from the previous call are used.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:fgbez00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grou00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:plin00.cmn'
      INCLUDE 'cmiss$reference:str00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER NDP(*),NT_STRIPE_DATA
      REAL*8 WD(NJM,*),ZD(NJM,*)
      CHARACTER ERROR*(*),GROUP_NAME*(*)
!     Local Variables
      INTEGER I,IBOUND,ICTR(80),IFAIL,IFRINGE(200),IIMAX,INTYPE,IP,
     '  IPRINT,IRANK(200),ISTATE(2),IV,IW(5),J,K1,K2,KK,
     '  LDM,LH,LIW,LW,M1,M2,MAXCAL,N1,N2,NCTLPTS(80),nd,NDUM,
     '  NDUM2,nj,NOGRPL,NOPTS,np,NPOLYLINES,NSEGS(80)
      REAL*8 ACN,ARCLENT,ARCLENT1,ARCLENT2,BCN,BL(2),BU(2),
     '  CHUNK1(200,4),CHUNK2(200,16),DELS,DELTA(2),DIST,ETA,F,FEST,G(2),
     '  HESD(2),HESL(5),PTS(1040,2),SLENC,SPACID,STEPMX,TAB1(200,4),
     '  TAB2(200,4),W(3800),   !Maximum #pts is 3800
     '  WW(18),X(2),XFA,XFB,XKSIB,XLIST(80,15,2),
     '  XTOL
      CHARACTER ORDER*1
      LOGICAL LOCSCH
      EXTERNAL E04JBQ,FUNCT1,MONIT

      CALL ENTERS('CALC_STRIPE_PTS',*9999)

      NOGRPL=0
c     loop through all groups to findout which nogrpl corresponds to
c     name group_name
      DO I=1,NTGRPL
        IF(LAGRPL(I).EQ.GROUP_NAME) THEN
          NOGRPL=I
        ELSE
        ENDIF
      ENDDO


      NPOLYLINES=LIGRPL(0,NOGRPL)       !Number of polylines

      II=1
      DO I=1,NPOLYLINES
        NSEGS(I)=NT_PLIN_SECTIONS(LIGRPL(I,NOGRPL)) !# segments in Ith polyline
        NCTLPTS(I)=NSEGS(I)*3+1  !# control pts in Ith polyline
        DO J=1,NCTLPTS(I)
          PTS(II,1)=PLIN_DATA(1,J,LIGRPL(I,NOGRPL))
          PTS(II,2)=PLIN_DATA(2,J,LIGRPL(I,NOGRPL))
          II=II+1
        ENDDO      !end J
      ENDDO        !end I
      IIMAX=II-1    !the combined total of nodes and c.pts over all polylines
C ... Initialize XLIST(I,J,K)
      DO I=1,80
        DO J=1,15
          XLIST(I,J,1)=9.9999D0
        ENDDO    !end J
      ENDDO      !end I

      II=0
      DO I=1,80
         ICTR(I)=0
      ENDDO

      DO I=1,NPOLYLINES-1
        DO J=I+1,NPOLYLINES
          DO K1=1,NSEGS(I)
            DO K2=1,NSEGS(J)
              NDUM=0
              DO KK=1,I
                IF(KK.GT.1) NDUM=NDUM+NCTLPTS(KK-1)
              ENDDO      !end KK
              IF(K1.GT.1) NDUM=NDUM+(K1-1)*3
              XA(1)=PTS(NDUM+1,1)
              XA(2)=PTS(NDUM+2,1)
              XA(3)=PTS(NDUM+3,1)
              XA(4)=PTS(NDUM+4,1)
              YA(1)=PTS(NDUM+1,2)
              YA(2)=PTS(NDUM+2,2)
              YA(3)=PTS(NDUM+3,2)
              YA(4)=PTS(NDUM+4,2)

              NDUM=0
              DO KK=1,J
                IF(KK.GT.1) NDUM=NDUM+NCTLPTS(KK-1)
              ENDDO     !end KK
              IF(K2.GT.1) NDUM=NDUM+(K2-1)*3
              XB(1)=PTS(NDUM+1,1)
              XB(2)=PTS(NDUM+2,1)
              XB(3)=PTS(NDUM+3,1)
              XB(4)=PTS(NDUM+4,1)
              YB(1)=PTS(NDUM+1,2)
              YB(2)=PTS(NDUM+2,2)
              YB(3)=PTS(NDUM+3,2)
              YB(4)=PTS(NDUM+4,2)

C ...         parameters needed for E04JBF
              IWHICHF=1    !in common block to E04JBF
              KTYP26=3
              KTYP27=1
              np=2
              IPRINT=-1       !MONIT is called at the last iteration
              LOCSCH=.TRUE.
              INTYPE=0
              MAXCAL=560
              ETA=0.5D0
              XTOL=0.000001D0
              STEPMX=20.0D0
              FEST=0.0D0
              DELTA(1)=1.0D-4
              DELTA(2)=1.0D-4
              IBOUND=0
              BL(1)=0.00001D0
              BL(2)=0.00001D0
              BU(1)=0.99999D0
              BU(2)=0.99999D0
              X(1)=0.5D0
              X(2)=0.5D0
              LH=5   !can be 1
              LIW=5
              LW=18
              IFAIL=1

              CALL E04JBF(np,FUNCT1,MONIT,IPRINT,LOCSCH,INTYPE,E04JBQ,
     '                    MAXCAL,ETA,XTOL,STEPMX,FEST,DELTA,IBOUND,
     '                    BL,BU,X,HESL,LH,HESD,ISTATE,F,G,IW,LIW,
     '                    WW,LW,IFAIL)

c             ictr(i) gives the current number of intersections located
c             in bezier polyline i

              IF(F.LT.1.0D-04) THEN      !intersection found
                ICTR(I)=ICTR(I)+1
                XLIST(I,ICTR(I),1)=X(1)
                XLIST(I,ICTR(I),2)=K1
                ICTR(J)=ICTR(J)+1
                XLIST(J,ICTR(J),1)=X(2)
                XLIST(J,ICTR(J),2)=K2
              ELSE                       !no intersection found
              ENDIF

            ENDDO    !end K2
          ENDDO      !end K1
        ENDDO        !end J
      ENDDO          !end I


C ... compute the total combined arc length of polylines
      ARCLENT=0.0D0
      DO I=1,NPOLYLINES
        DO J=1,NSEGS(I)
          NDUM=0
          DO KK=1,I
             IF(KK.GT.1) NDUM=NDUM+NCTLPTS(KK-1)
          ENDDO          !end KK
          NDUM=NDUM+3*(J-1)

          XA(1)=PTS(NDUM+1,1)
          XA(2)=PTS(NDUM+2,1)
          XA(3)=PTS(NDUM+3,1)
          XA(4)=PTS(NDUM+4,1)
          YA(1)=PTS(NDUM+1,2)
          YA(2)=PTS(NDUM+2,2)
          YA(3)=PTS(NDUM+3,2)
          YA(4)=PTS(NDUM+4,2)

          NSIMP=100
          XFA=0.0D0
          XFB=1.0D0
          CALL SIMPS(XA,YA,NSIMP,XFA,XFB,ACDIFF)
          ARCLENT=ARCLENT+ACDIFF
         ENDDO      !end J
       ENDDO        !end I


c ... Compute the ideal spacing based on total number
c     of data points prescribed.  recall NT_STRIPE_DATA was passed
c     into subroutine as input.

      IF(NT_STRIPE_DATA.NE.0) SPACID=ARCLENT/NT_STRIPE_DATA


C CPB 21/3/94 I Think TAB would be better suited as an integer array


c ... Define d.p. array TAB1(100,4) to hold all
c     unranked points in spatial order
c     TAB1 is a table with columns representing:
c     1) index of polyline
c     2) bezier segment number
c     3) XKSI
c     4) type (intersection:5, head endpoint of curve:6,
c              tail endpoint of curve:7)
c              note: endpoints may be excluded from analysis

c ... First, compute array tab1 containing intersections
      II=0
      DO I=1,NPOLYLINES
        DO J=1,15
          IF(XLIST(I,J,1).NE.9.9999D0) THEN
            II=II+1
            TAB1(II,1)=DBLE(I)     !curve number
            NDUM=0
            DO KK=1,I
              IF(KK.GT.1) NDUM=NDUM+NSEGS(KK-1)
            ENDDO   !end KK
            TAB1(II,2)=NDUM+XLIST(I,J,2)   !bezier segment no.
            TAB1(II,3)=XLIST(I,J,1)        !XKSI
            TAB1(II,4)=5.0D0               !type
3001        FORMAT(I4,4(F12.6,1X))
          ELSE
          ENDIF
        ENDDO   !end J
      ENDDO     !end I


c ... Second, put in the endpoints of the bezier polylines
      DO I=1,NPOLYLINES

        II=II+1
        TAB1(II,1)=DBLE(I)
        TAB1(II,2)=0.0D0
        DO KK=1,I
          IF(KK.GT.1) TAB1(II,2)=TAB1(II,2)+DBLE(NSEGS(KK-1))
        ENDDO   !end KK
        TAB1(II,2)=DBLE(TAB1(II,2))+1.0D0
        TAB1(II,3)=0.0D0
        TAB1(II,4)=6.0D0
        II=II+1
        TAB1(II,1)=DBLE(I)
        TAB1(II,2)=0.0D0
        DO KK=1,I
          TAB1(II,2)=TAB1(II,2)+DBLE(NSEGS(KK))
        ENDDO  !end KK
        TAB1(II,3)=1.0D0
        TAB1(II,4)=7.0D0
      ENDDO    !end I


c ... Now rank the elements of  using M01DEF
      LDM=200
      M1=1
      M2=II
      N1=1
      N2=4
      ORDER='A'
      IFAIL=0
      CALL M01DEF(TAB1,LDM,M1,M2,N1,N2,ORDER,IRANK,IFAIL)

c ... sort matrix tab1, putting sorted rows in matrix tab2

      DO I=1,II
        TAB2(IRANK(I),1)=TAB1(I,1)
        TAB2(IRANK(I),2)=TAB1(I,2)
        TAB2(IRANK(I),3)=TAB1(I,3)
        TAB2(IRANK(I),4)=TAB1(I,4)
      ENDDO !end I


c ... Now deterDMIN1(e array CHUNK1(i,4) which identifies
c     the cartesian coordinates of the bezier segment or segments
c     in which the Ith chunk lies.
c     CHUNK1(i,1)= global index of bezier segment in which head endpoint
c                  of chunk i lies
c     CHUNK1(i,2)= XKSI of head endpoint
c     CHUNK1(i,3)= global index of bezier segment in which tail endpoint
c                  of chunk i lies
c     CHUNK1(i,4)= XKSI of tail endpoint

      IP=0        !We proceed with assumption that
      IV=0        !II=total number of rows in TAB2
300   IP=IP+1
      IV=IV+1
      IF(IV.GE.II) GO TO 311
      IF(TAB2(IV,4).EQ.7) IV=IV+1
      CHUNK1(IP,1)=TAB2(IV,2)
      CHUNK1(IP,2)=TAB2(IV,3)
      CHUNK1(IP,3)=TAB2(IV+1,2)
      CHUNK1(IP,4)=TAB2(IV+1,3)
      GO TO 300
311   IP=IP-1     !total number of chunks


C CBP 21/3/94 There is alot of assigning integers to doubles and
C checking of doubles against integers and vise versa. This code 
C should be looked at !!!!!!!!!!!!!


c ... Determine array CHUNK2(200,16) which identifies
c     nodes and control points corresponding to the
c     bezier segment in which the chunk's 2 endpoints lie.
c     format x1 x2 x3 x4 y1 y2 y3 y4 x1 x2 ... etc.
c            .... 8 pairs of x and y ...
c     1st  8 columns: nodes of bez. seg. containing "head" point
c     2nd  8 columns: nodes of bez. seg. containing "tail" point


      DO I=1,IP    !1 through total number of chunks

c ..    Identify the index of the node which is 1 less than
c       the head node of a given bezier segment
        NDUM=0       !no. accum. c-pts of all previous b-segments
        NDUM2=0      !no. accum. b-segments in all p-lines to present
        DO KK=1,NPOLYLINES
          IF(KK.GT.1) NDUM=NDUM+NCTLPTS(KK-1)
          NDUM2=NDUM2+NSEGS(KK)
          IF(NDUM2.GE.CHUNK1(I,1)) GO TO 310
        ENDDO
310     NDUM=NDUM+3*(CHUNK1(I,1)-NDUM2+NSEGS(KK)-1)

        CHUNK2(I,1)=PTS(NDUM+1,1)
        CHUNK2(I,2)=PTS(NDUM+2,1)
        CHUNK2(I,3)=PTS(NDUM+3,1)
        CHUNK2(I,4)=PTS(NDUM+4,1)
        CHUNK2(I,5)=PTS(NDUM+1,2)
        CHUNK2(I,6)=PTS(NDUM+2,2)
        CHUNK2(I,7)=PTS(NDUM+3,2)
        CHUNK2(I,8)=PTS(NDUM+4,2)

c ..    Identify the index of the node which is 1 less than
c       the tail node of a given bezier segment
        NDUM=0
        NDUM2=0
        DO KK=1,NPOLYLINES

          IF(KK.GT.1) NDUM=NDUM+NCTLPTS(KK-1)
          NDUM2=NDUM2+NSEGS(KK)
          IF(NDUM2.GE.CHUNK1(I,3)) GO TO 312
        ENDDO
312     NDUM=NDUM+3*(CHUNK1(I,3)-NDUM2+NSEGS(KK)-1)

        CHUNK2(I,9) =PTS(NDUM+1,1)
        CHUNK2(I,10)=PTS(NDUM+2,1)
        CHUNK2(I,11)=PTS(NDUM+3,1)
        CHUNK2(I,12)=PTS(NDUM+4,1)
        CHUNK2(I,13)=PTS(NDUM+1,2)
        CHUNK2(I,14)=PTS(NDUM+2,2)
        CHUNK2(I,15)=PTS(NDUM+3,2)
        CHUNK2(I,16)=PTS(NDUM+4,2)
      ENDDO    !end I
3002  FORMAT(I4,16(F12.6,1X))

c ... Determine array IFRINGE(I) which contains 1 or 0
c     depending on whether chunk is a fringe or an
c     internal segment
c                        IFRINGE(I)=1 ......fringe
c                        IFRINGE(I)=0 ......internal segment
      IFRINGE(1)=1
      IFRINGE(IP)=1
      DO I=2,IP-1
        DIST=((CHUNK2(I,1)-CHUNK2(I-1,9))**2+
     '        (CHUNK2(I,5)-CHUNK2(I-1,13))**2)**.5
        IF(DIST.LT.0.01) THEN
          IFRINGE(I)=0 !This condition will regard last chunk in polyline as
!                       internal; then correction is made when next chunk is
!                       considered; see below
        ELSE
          IFRINGE(I)=1
          IFRINGE(I-1)=1    !Here is where correction is made
        ENDIF
      ENDDO     !end i

c ... For each chunk, determine if it straddles a node.
c     if so, integrate the lengths ST1 and ST2, the two
c     connecting pieces. ST=ST1+ST2.  if not, compute
c     ST directly.

      II=0               !data point counter
      IWHICHF=2          !in common block to FUNCT1
      KTYP26=4
      KTYP27=1
      DO I=1,IP
        IF(IFRINGE(I).EQ.1) GO TO 400
        IF(CHUNK1(I,1).EQ.CHUNK1(I,3)) THEN  !no straddle

          XKSIA=CHUNK1(I,2)
          XKSIB=CHUNK1(I,4)
          NSIMP=50

          XA(1)=CHUNK2(I,1)
          XA(2)=CHUNK2(I,2)
          XA(3)=CHUNK2(I,3)
          XA(4)=CHUNK2(I,4)
          YA(1)=CHUNK2(I,5)
          YA(2)=CHUNK2(I,6)
          YA(3)=CHUNK2(I,7)
          YA(4)=CHUNK2(I,8)

          CALL SIMPS(XA,YA,NSIMP,XKSIA,XKSIB,ACDIFF)

c ...     Determine datapoint interval
          IF(NT_STRIPE_DATA.EQ.0) THEN  !2nd routine call
            NOPTS=NOPTS_PRIOR(I)
          ELSE                           !1st routine call
            NOPTS=ACDIFF/SPACID+0.5D0
            NOPTS_PRIOR(I)=NOPTS         !prep for 2nd routine call
          ENDIF
          IF(NOPTS.GT.0) THEN
            DELS=ACDIFF/NOPTS            !spacing for this chunk
          ELSE
            DELS=0.0D0                   !loop will be skipped anyway
          ENDIF

c ...     Do a loop from 0,NOPTS:
c         for each value of s, determine the value
c         of XKSI which produces it (using E04JBF)
c         compute x,y value of that point: define ZD([1,2],I)

          XKSIA=CHUNK1(I,2)    !constrained to be at end of chunk
          BL(2)=0.00001D0      !contrain other 1 variable
          BU(2)=0.99999D0

c ...     place point at "head" end of bezier
          DIST=0.0D0
          IF(I.EQ.1) THEN
            DIST=10.0D0
          ELSE
            DIST=((CHUNK2(I,1)-CHUNK2(I-1,9))**2+
     '           (CHUNK2(I,5)-CHUNK2(I-1,13))**2)**0.5D0
          ENDIF
          IF(DIST.GT.0.001D0) THEN   !chunk is first in polyline
!                                    (criterion had been .01)
            II=II+1
            ZD(1,II)=CHUNK2(I,1)
            ZD(2,II)=CHUNK2(I,5)
          ELSE
          ENDIF

          DO J=0,NOPTS
            NSIMP=40
            SLEN=J*DELS

            XA(1)=CHUNK2(I,1)
            XA(2)=CHUNK2(I,2)
            XA(3)=CHUNK2(I,3)
            XA(4)=CHUNK2(I,4)
            YA(1)=CHUNK2(I,5)
            YA(2)=CHUNK2(I,6)
            YA(3)=CHUNK2(I,7)
            YA(4)=CHUNK2(I,8)

C ...       INPUTS TO E04JBF
            np=2
            IPRINT=-1       !MONIT is called at the last iteration
            LOCSCH=.TRUE.
            INTYPE=0
            MAXCAL=560
            ETA=0.5D0
            XTOL=0.0D0
            STEPMX=2.0D0
            FEST=0.0D0
            DELTA(1)=1.0D-4
            DELTA(2)=1.0D-4
            IBOUND=0
            BL(1)=CHUNK1(I,2)
            BU(1)=CHUNK1(I,4)
            X(1)=(CHUNK1(I,2)+CHUNK1(I,4))/2.0D0      !guess for XKSIB
            LH=5   !can be 1
            LIW=5
            LW=18
            IFAIL=1

            CALL E04JBF(np,FUNCT1,MONIT,IPRINT,LOCSCH,INTYPE,E04JBQ,
     '                  MAXCAL,ETA,XTOL,STEPMX,FEST,DELTA,IBOUND,
     '                  BL,BU,X,HESL,LH,HESD,ISTATE,F,G,IW,LIW,
     '                  WW,LW,IFAIL)

            XKSIB=X(1)      !the converged value of XKSIB
            IF(DABS(F).GT.1.0001D0) THEN
              WRITE(OP_STRING,'(''F doesnt vanish'',7(F10.5,1X))')
     '          F,XKSIA,X(1),BL(1),
     '          BU(1),SLEN,ACDIFF
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
C1005        FORMAT('F doesnt vanish',7(F10.5,1X))

c ...       for computed value of xksi,compute x,y
            II=II+1
            IF(XKSIB.GT.1.0D0.OR.XKSIB.LT.0.0D0) PRINT*,
     '        'XKSIB is out of the proper range of [0,1]'

            ACN=XA(1)*(1.0D0-XKSIB)**3
            ACN=ACN+XA(2)*3.0D0*XKSIB*(1.0D0-XKSIB)**2
            ACN=ACN+XA(3)*3.0D0*XKSIB*XKSIB*(1.0D0-XKSIB)
            ACN=ACN+XA(4)*XKSIB**3
            ZD(1,II)=ACN           !X global pt.

            BCN=YA(1)*(1.0D0-XKSIB)**3
            BCN=BCN+YA(2)*3.0D0*XKSIB*(1.0D0-XKSIB)**2
            BCN=BCN+YA(3)*3.0D0*XKSIB*XKSIB*(1.0D0-XKSIB)
            BCN=BCN+YA(4)*XKSIB**3
            ZD(2,II)=BCN           !Y global pt.
            W(II)=DELS             !weighting associated with point
          ENDDO !end J

        ELSE                       !chunk straddles an intersection

c ...     compute the chunk's arc length.  we will compute
c         the arclength in 2 parts

          NSIMP=50
          XKSIA=CHUNK1(I,2)        !endpoints
          XKSIB=1.0D0                !of first part

          XA(1)=CHUNK2(I,1)        !control pts of bezier segment
          XA(2)=CHUNK2(I,2)        !in which first part of
          XA(3)=CHUNK2(I,3)        !chunk lies
          XA(4)=CHUNK2(I,4)
          YA(1)=CHUNK2(I,5)
          YA(2)=CHUNK2(I,6)
          YA(3)=CHUNK2(I,7)
          YA(4)=CHUNK2(I,8)

          CALL SIMPS(XA,YA,NSIMP,XKSIA,XKSIB,ARCLENT1)
c ...     This completes integration of first part of chunk

          NSIMP=50
          XKSIA=0.0D0               !endpoints
          XKSIB=CHUNK1(I,4)         !of second part
          XA(1)=CHUNK2(I,9)         !control pts of bezier segment
          XA(2)=CHUNK2(I,10)        !in which second part
          XA(3)=CHUNK2(I,11)        !of chunk lies
          XA(4)=CHUNK2(I,12)
          YA(1)=CHUNK2(I,13)
          YA(2)=CHUNK2(I,14)
          YA(3)=CHUNK2(I,15)
          YA(4)=CHUNK2(I,16)

          CALL SIMPS(XA,YA,NSIMP,XKSIA,XKSIB,ARCLENT2)
c ...     This completes integration of the second chunk
c         arclent=arclent1+arclent2


c ...     Determine datapoint interval
          IF(NT_STRIPE_DATA.EQ.0) THEN            !2nd routine call
            NOPTS=NOPTS_PRIOR(I)
          ELSE                                    !1st routine call
            NOPTS=(ARCLENT1+ARCLENT2)/SPACID+0.50D0
            NOPTS_PRIOR(I)=NOPTS                  !Prep for 2nd routine call
          ENDIF
          IF(NOPTS.GT.0) THEN
            DELS=(ARCLENT1+ARCLENT2)/NOPTS        !Spacing for this chunk
          ELSE
            DELS=0.0D0                            !Loop will be skipped anyway
          ENDIF


c ...     Do a loop from 0,nopts
c         for each value of S, determine the value of XKSI
c         which produces it (using E04JBF)
c         compute x,y value of that point: define ZD([1,2],II)

          XKSIA=CHUNK1(I,2)
          XKSIB=XKSIA

c ...     put point at "head" end of bezier segment
          DIST=0.
          IF(I.EQ.1) THEN
            DIST=10.0D0
          ELSE
            DIST=((CHUNK2(I,1)-CHUNK2(I-1,9))**2+
     '            (CHUNK2(I,5)-CHUNK2(I-1,13))**2)**0.5D0
          ENDIF
          IF(DIST.GT.0.001D0) THEN       !chunk is first in polyline
!                                        (criterion had been .01)
            II=II+1
            ZD(1,II)=CHUNK2(I,1)
            ZD(2,II)=CHUNK2(I,5)
          ELSE
          ENDIF

          DO J=0,NOPTS
            SLENC=J*DELS                       !Cumulative length
            IF(SLENC.LE.ARCLENT1) THEN        !First chunk
              NSIMP=40
              SLEN=SLENC
              XA(1)=CHUNK2(I,1)      !Control pts of bezier segment
              XA(2)=CHUNK2(I,2)      !in which first part of
              XA(3)=CHUNK2(I,3)      !chunk lies
              XA(4)=CHUNK2(I,4)
              YA(1)=CHUNK2(I,5)
              YA(2)=CHUNK2(I,6)
              YA(3)=CHUNK2(I,7)
              YA(4)=CHUNK2(I,8)

C ...         Inputs to E04JBF
              np=2
              IPRINT=-1       !MONIT is called at the last iteration
              LOCSCH=.TRUE.
              INTYPE=0
              MAXCAL=560
              ETA=0.5D0
              XTOL=0.0D0
              STEPMX=2.0D0
              FEST=0.0D0
              DELTA(1)=1.0D-4
              DELTA(2)=1.0D-4
              IBOUND=0
              BL(1)=0.00001D0
              BU(1)=0.99999D0
              X(1)=(XKSIA+1.0D0)/2.0D0               !Guess for XKSIB
              LH=5   !Can be 1
              LIW=5
              LW=18
              IFAIL=1

              CALL E04JBF(np,FUNCT1,MONIT,IPRINT,LOCSCH,INTYPE,E04JBQ,
     '                    MAXCAL,ETA,XTOL,STEPMX,FEST,DELTA,IBOUND,
     '                    BL,BU,X,HESL,LH,HESD,ISTATE,F,G,IW,LIW,
     '                    WW,LW,IFAIL)


              XKSIB=X(1)      !The converged value of XKSIB
c ...         for computed value of XKSI, compute x,y
              II=II+1

              ACN=XA(1)*(1.0D0-XKSIB)**3
              ACN=ACN+XA(2)*3.0D0*XKSIB*(1.0D0-XKSIB)**2
              ACN=ACN+XA(3)*3.0D0*XKSIB*XKSIB*(1.0D0-XKSIB)
              ACN=ACN+XA(4)*XKSIB**3
              ZD(1,II)=ACN                  !X global pt.

              BCN=YA(1)*(1.0D0-XKSIB)**3
              BCN=BCN+YA(2)*3.0D0*XKSIB*(1.0D0-XKSIB)**2
              BCN=BCN+YA(3)*3.0D0*XKSIB*XKSIB*(1.0D0-XKSIB)
              BCN=BCN+YA(4)*XKSIB**3
              ZD(2,II)=BCN                  !Y global pt
              W(II)=DELS

            ELSE                   !Second chunk
              SLEN=SLENC-ARCLENT1
              XKSIA=0.0D0
              NSIMP=40
              XA(1)=CHUNK2(I,9)    !Control points of bezier segment
              XA(2)=CHUNK2(I,10)   !in which first part of
              XA(3)=CHUNK2(I,11)   !chunk lies
              XA(4)=CHUNK2(I,12)
              YA(1)=CHUNK2(I,13)
              YA(2)=CHUNK2(I,14)
              YA(3)=CHUNK2(I,15)
              YA(4)=CHUNK2(I,16)

C ...         INPUTS TO E04JBF
              np=2    !NP = 4? See above!!
              IPRINT=-1       !MONIT is called at the last iteration
              LOCSCH=.TRUE.
              INTYPE=0
              MAXCAL=560
              ETA=0.5D0
              XTOL=0.0D0
              STEPMX=2.0D0
              FEST=0.0D0
              DELTA(1)=1.0D-4
              DELTA(2)=1.0D-4
              IBOUND=0
              BL(1)=0.00001D0
              BU(1)=0.99999D0
              X(1)=XKSIB/2.0D0        !Guess for XKSIB
              LH=5   !Can be 1
              LIW=5
              LW=18
              IFAIL=1
              CALL E04JBF(np,FUNCT1,MONIT,IPRINT,LOCSCH,INTYPE,E04JBQ,
     '                    MAXCAL,ETA,XTOL,STEPMX,FEST,DELTA,IBOUND,
     '                    BL,BU,X,HESL,LH,HESD,ISTATE,F,G,IW,LIW,
     '                    WW,LW,IFAIL)


              XKSIB=X(1)      !The converged value of xksib

c ...         For computed value of XKSI, compute x,y
              II=II+1
              IF(XKSIB.GT.1.0D0.OR.XKSIB.LT.0.0D0) PRINT*,
     '            'XKSIB is out of the proper range of [0,1]'
              ACN=XA(1)*(1.0D0-XKSIB)**3
              ACN=ACN+XA(2)*3.0D0*XKSIB*(1.0D0-XKSIB)**2
              ACN=ACN+XA(3)*3.0D0*XKSIB*XKSIB*(1.0D0-XKSIB)
              ACN=ACN+XA(4)*XKSIB**3
              ZD(1,II)=ACN                  !X global pt.

              BCN=YA(1)*(1.0D0-XKSIB)**3
              BCN=BCN+YA(2)*3.0D0*XKSIB*(1.0D0-XKSIB)**2
              BCN=BCN+YA(3)*3.0D0*XKSIB*XKSIB*(1.0D0-XKSIB)
              BCN=BCN+YA(4)*XKSIB**3
              ZD(2,II)=BCN                  !Y global pt.
              W(II)=DELS

            ENDIF
          ENDDO    !END J
        ENDIF
 400  ENDDO        !END I; All chunks have been analyzed


C ... Set value of NDT, the dimension of the ZD array
      IF(NT_STRIPE_DATA.NE.0) THEN      !1st call
        NDT=II                  !parameter to be passed back to FE24
        NDT_PRIOR=NDT           !prep for 2nd call
      ELSE
        NDT=NDT_PRIOR           !2nd call
      ENDIF

!     Use 1 for weights
      DO nd=1,NDT
        NDP(nd)=ND
        DO nj=1,NJT
          WD(nj,nd)=1.0D0
        ENDDO
      ENDDO

 9998 CALL EXITS('CALC_STRIPE_PTS')
      RETURN
 9999 CALL ERRORS('CALC_STRIPE_PTS',ERROR)
      CALL EXITS('CALC_STRIPE_PTS')
      RETURN 1
      END


      SUBROUTINE CALC_NVNEP(IBT,INP,NBJ,NEELEM,NJE,NJP,NPE,NPNODE,
     '  nr,NVNE,NVNP,NXI,XP,ERROR,*)

C**** Calculates the solution mapping arrays NVNE and NVNP
!new MPN 13/12/93: using nv in ZP to store multiple theta's for nodes 
!                  that lie on axis in polar or prolate coords
!                  Assumes that not more than 4 of the
!                  nodes of any element lie on the axis (ie not single
!                  element, FULL prolate or sphere etc)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
!     Parameter List
      INTEGER IBT(2,NIM,*),INP(NNM,NIM,*),NBJ(NCM,NJM,*),
     '  NEELEM(0:NEM,0:*),
     '  NJE(*),NJP(*),NPE(NNM,NBM,*),NPNODE(0:NPM,0:*),
     '  nr,NVNE(NNM,NBM,*),NVNP(0:NVM,NJM,NPM,*),NXI(-NIM:NIM,0:*)
      REAL*8 XP(NKM,NVM,NJM,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nc,NC1,NC2,NC3,NC4,ne,NE1,NE3,NESTART,
     '  nj,njj,NJJ1,NJJ2,nn,NONODE,NOELEM,np,NPE1,NPE2,
     '  NPE3,NPE4,NPE5,NPE6,NPE7,NPE8,nv,nv1,nv2,nv3,nv4
      REAL*8 AXIS_COORD1,AXIS_COORD2,TOL
      LOGICAL CHECKED(100),FOUND,LASTXI1,LASTXI3 
      DATA TOL /1.0d-6/

      CALL ENTERS('CALC_NVNEP',*9999)

      nc=1 ! temporary cpb 23/11/94
      nv=1 ! temporary cpb 23/11/94

!     Initialise NVNP and NVNE nv mapping arrays
      DO NONODE=1,NPNODE(0,nr)
        np=NPNODE(NONODE,nr)
        DO NJJ1=1,3  ! geometry/fibres/field
          DO NJJ2=1,NJ_LOC(NJJ1,0)
            nj=NJ_LOC(NJJ1,NJJ2)
            NVNP(0,nj,np,nc)=1
          ENDDO
        ENDDO
      ENDDO
      DO NOELEM=1,NEELEM(0,nr)
        ne=NEELEM(NOELEM,nr)
        DO nb=1,NBT
          DO nn=1,NNT(nb)
            NVNE(nn,nb,ne)=1
          ENDDO
        ENDDO
      ENDDO

      IF(ITYP11(1).GT.1.AND.NJE(1).GT.1) THEN
        IF(ITYP11(1).EQ.2) THEN       !cyl. polar
          nj=2   !theta
          njj=1  !radius       
          AXIS_COORD1=0.0d0    !radius=0 is on axis
          AXIS_COORD2=0.0d0    !radius=0 is on axis
        ELSE IF(ITYP11(1).EQ.3) THEN  !sph. polar
          nj=2   !theta
          njj=3  !phi
          AXIS_COORD1= PI/2.0d0 !phi= pi/2 is on axis
          AXIS_COORD2=-PI/2.0d0 !phi=-pi/2 is on axis
        ELSE IF(ITYP11(1).EQ.4) THEN  !prolate sph.
          nj=3   !theta
          njj=2  !mu
          AXIS_COORD1=0.0d0    !mu= 0 is on axis
          AXIS_COORD2=0.0d0    !mu=pi is on axis
        ELSE IF(ITYP11(1).EQ.5) THEN  !oblate sph.
        ENDIF
        nb=NBJ(1,nj,NEELEM(1,nr))  !geometric basis for theta in elem 1
        IF(NIT(nb).EQ.3) THEN
          CALL ASSERT(NEELEM(0,nr).LT.100,
     '      ' >>Increase size of CHECKED array to NEM',ERROR,*9999)
          DO NOELEM=1,NEELEM(0,nr) !initialise 
            ne=NEELEM(NOELEM,nr)
            CHECKED(ne)=.FALSE.
          ENDDO

C ***     Determine element neighbours
          CALL NENXI(IBT,INP,NBJ,NEELEM,NPE,NXI,ERROR,*9999)

C ***     Calc nv maps for global nodes (NVNP) and elem nodes (NVNE)
          NESTART=1
          DO WHILE(NESTART.NE.0)
            NESTART=0
!           Loop through all unchecked elements to find any elements
!           with 2 nodes repeated for Xi1=0,1
            NOELEM=0
            FOUND=.FALSE.
            DO WHILE(NOELEM.LT.NEELEM(0,nr).AND..NOT.FOUND)
              NOELEM=NOELEM+1
              IF(.NOT.CHECKED(NOELEM)) THEN
                ne=NEELEM(NOELEM,nr)
                NPE1=NPE(1,nb,ne)
                NPE2=NPE(2,nb,ne)
                NPE3=NPE(3,nb,ne)
                NPE4=NPE(4,nb,ne)
                NPE5=NPE(5,nb,ne)
                NPE6=NPE(6,nb,ne)
                NPE7=NPE(7,nb,ne)
                NPE8=NPE(8,nb,ne)
!               Check for zero radius (cyl or sph) or mu (prolate) 
!               at repeated element nodes
                IF(   NPE1.EQ.NPE2          !1st & 2nd nodes in common
     '            .AND.(DABS(XP(1,nv,njj,NPE1)-AXIS_COORD1).LT.TOL
     '              .OR.DABS(XP(1,nv,njj,NPE1)-AXIS_COORD2).LT.TOL)
     '            .AND.NPE3.NE.NPE4
     '            .OR.NPE3.EQ.NPE4          !3rd & 4th nodes in common 
     '            .AND.(DABS(XP(1,nv,njj,NPE3)-AXIS_COORD1).LT.TOL
     '              .OR.DABS(XP(1,nv,njj,NPE3)-AXIS_COORD2).LT.TOL)
     '            .AND.NPE1.NE.NPE2
     '            .OR.NPE5.EQ.NPE6          !5th & 6th nodes in common
     '            .AND.(DABS(XP(1,nv,njj,NPE5)-AXIS_COORD1).LT.TOL
     '              .OR.DABS(XP(1,nv,njj,NPE5)-AXIS_COORD2).LT.TOL)
     '            .AND.NPE7.NE.NPE8
     '            .OR.NPE7.EQ.NPE8          !7th & 8th nodes in common
     '            .AND.(DABS(XP(1,nv,njj,NPE7)-AXIS_COORD1).LT.TOL
     '              .OR.DABS(XP(1,nv,njj,NPE7)-AXIS_COORD2).LT.TOL)
     '            .AND.NPE5.NE.NPE6) THEN
                  DO WHILE(NXI(-3,ne).NE.0)
                    ne=NXI(-3,ne) !move to innermost element (-Xi3 dirn)
                  ENDDO
                  NE3=NE
                  DO WHILE(NXI(-1,ne).NE.0.AND.NXI(-1,ne).NE.NE3)
                    ne=NXI(-1,ne) !move to first element in -Xi1 dirn
                  ENDDO           !and stop if have moved completely
                                  !around circumference
                  IF(NXI(-1,ne).EQ.0) THEN !have boundary
                    NESTART=NE
                  ELSE                     !complete ring
                    NESTART=NE3
                  ENDIF
                  FOUND=.TRUE.
                ENDIF
              ENDIF
            ENDDO !noelem while loop
20          IF(NESTART.NE.0) THEN  !unchecked elem with repeated nodes
              NE3=NESTART 
              LASTXI3=.FALSE.
              DO WHILE(.NOT.LASTXI3) !moves through elems in +Xi3 dirn
                NE1=NE3 
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' NE1='',I4,'' NE3='',I4)') NE1,NE3
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
                nv1=1
                nv2=1
                nv3=1
                nv4=1
                LASTXI1=.FALSE.
                DO WHILE(.NOT.LASTXI1) !moves around elems in +Xi1 dirn
                  NPE1=NPE(1,nb,NE1)   
                  NPE2=NPE(2,nb,NE1)
                  NPE3=NPE(3,nb,NE1)
                  NPE4=NPE(4,nb,NE1)   !global nodes numbers of
                  NPE5=NPE(5,nb,NE1)   !element nodes
                  NPE6=NPE(6,nb,NE1)
                  NPE7=NPE(7,nb,NE1)
                  NPE8=NPE(8,nb,NE1)
                  IF(NPE1.EQ.NPE2          !1st & 2nd nodes in common
     '              .AND.(DABS(XP(1,nv,njj,NPE1)-AXIS_COORD1).LT.TOL
     '                .OR.DABS(XP(1,nv,njj,NPE1)-AXIS_COORD2).LT.TOL)
     '              .AND.NPE3.NE.NPE4) THEN  
                    NVNE(1,nb,ne1)=nv1
                    nv1=nv1+1
                    NVNE(2,nb,ne1)=nv1
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' nv1='',I4)') nv1
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDIF
                  IF(NPE3.EQ.NPE4          !3rd & 4th nodes in common
     '              .AND.(DABS(XP(1,nv,njj,NPE3)-AXIS_COORD1).LT.TOL
     '                .OR.DABS(XP(1,nv,njj,NPE3)-AXIS_COORD2).LT.TOL)
     '              .AND.NPE1.NE.NPE2) THEN  
                    NVNE(3,nb,ne1)=nv2
                    nv2=nv2+1
                    NVNE(4,nb,ne1)=nv2
                  ENDIF
                  IF(NPE5.EQ.NPE6          !5th & 6th nodes in common
     '              .AND.(DABS(XP(1,nv,njj,NPE5)-AXIS_COORD1).LT.TOL
     '                .OR.DABS(XP(1,nv,njj,NPE5)-AXIS_COORD2).LT.TOL)
     '              .AND.NPE7.NE.NPE8) THEN  
                    NVNE(5,nb,ne1)=nv3
                    nv3=nv3+1
                    NVNE(6,nb,ne1)=nv3
                  ENDIF
                  IF(NPE7.EQ.NPE8          !7th & 8th nodes in common
     '              .AND.(DABS(XP(1,nv,njj,NPE7)-AXIS_COORD1).LT.TOL
     '                .OR.DABS(XP(1,nv,njj,NPE7)-AXIS_COORD2).LT.TOL)
     '              .AND.NPE5.NE.NPE6) THEN  
                    NVNE(7,nb,ne1)=nv4
                    nv4=nv4+1
                    NVNE(8,nb,ne1)=nv4
                  ENDIF
                  DO NOELEM=1,NEELEM(0,nr)  !Check off current element
                    IF(NE1.EQ.NEELEM(NOELEM,nr)) CHECKED(NOELEM)=.TRUE.
                  ENDDO
                  IF(NXI(1,NE1).EQ.0.OR.NXI(1,NE1).EQ.NE3) THEN
                    LASTXI1=.TRUE.  !no +Xi1 adj elems or back to 1st
                  ELSE
                    NE1=NXI(1,NE1)  !move to adjacent elem in +Xi1 dirn
                  ENDIF
                ENDDO  !.NOT.LASTXI1
                IF(nv1.GT.1) THEN          !1st & 2nd nodes in common
                  CALL ASSERT(nv1.LE.NVM,'>>Increase NVMX',ERROR,*9999)
                  IF(DOP) THEN
                    WRITE(OP_STRING,'('' NE1='',I4,'' NE3='',I4)') 
     '                NE1,NE3
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                  IF(NXI(1,NE1).EQ.NE3) THEN
                    NVNE(2,nb,ne1)=1
                    NVNP(0,nj,npe1,nc)=nv1-1
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' NVNP changed to'',I4)') nv1-1
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE
                    NVNP(0,nj,npe1,nc)=nv1
                  ENDIF
                ENDIF
                IF(nv2.GT.1) THEN          !3rd & 4th nodes in common
                  CALL ASSERT(nv2.LE.NVM,'>>Increase NVMX',ERROR,*9999)
                  IF(NXI(1,NE1).EQ.NE3) THEN
                    NVNE(4,nb,ne1)=1
                    NVNP(0,nj,npe3,nc)=nv2-1
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' NVNP changed to'',I4)') nv2-1
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE
                    NVNP(0,nj,npe3,nc)=nv2
                  ENDIF
                ENDIF
                IF(nv3.GT.1) THEN          !5th & 6th nodes in common
                  CALL ASSERT(nv3.LE.NVM,'>>Increase NVMX',ERROR,*9999)
                  IF(NXI(1,NE1).EQ.NE3) THEN
                    NVNE(6,nb,ne1)=1
                    NVNP(0,nj,npe5,nc)=nv3-1
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' NVNP changed to'',I4)') nv3-1
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE
                    NVNP(0,nj,npe5,nc)=nv3
                  ENDIF
                ENDIF
                IF(nv4.GT.1) THEN          !7th & 8th nodes in common
                  CALL ASSERT(nv4.LE.NVM,'>>Increase NVMX',ERROR,*9999)
                  IF(NXI(1,NE1).EQ.NE3) THEN
                    NVNE(8,nb,ne1)=1
                    NVNP(0,nj,npe7,nc)=nv4-1
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' NVNP changed to'',I4)') nv4-1
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ENDIF
                  ELSE
                    NVNP(0,nj,NPE7,nc)=nv4
                  ENDIF
                ENDIF
                NE3=NXI(3,NE3)  !move to next element in +Xi3
                IF(NE3.EQ.0.OR.NE3.EQ.NESTART) LASTXI3=.TRUE.
              ENDDO  !.NOT.LASTXI3
            ENDIF
          ENDDO !NESTART.NE.0
        ENDIF
      ENDIF

! Fill in other NVNP values (not used yet)
      DO NONODE=1,NPNODE(0,nr)
        np=NPNODE(NONODE,nr)
        DO NJJ1=1,3  ! geometry/fibres/field
          DO NJJ2=1,NJ_LOC(NJJ1,0)
            nj=NJ_LOC(NJJ1,NJJ2)
            DO nv=1,NVNP(0,nj,np,nc)
              NVNP(nv,nj,np,nc)=nv
            ENDDO
          ENDDO
        ENDDO
      ENDDO


      CALL EXITS('CALC_NVNEP')
      RETURN
 9999 CALL ERRORS('CALC_NVNEP',ERROR)
      CALL EXITS('CALC_NVNEP')
      RETURN 1
      END

      SUBROUTINE CPCP2(ISPL,NQE,NUM_CP,CP_LOCAL,NUM_CPSUB,CPSUB,ERROR,*)

C#### Subroutine: CPCP2
C###  Description:
C###    CPCP2 returns in CPSUB the control points obtained by
C###    refining the CP_LOCAL control point mesh.  Cubics are divided into 
C###    halves, resulting in 7 independent points.  Linear Lagrange are 
C###    divided into 6 equal parts, returning 7 points.
C     ISPL=1 : Lagrange/Hermite tensor product
C     ISPL=2 : B-spline tensor product

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:bspln00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER ISPL,NQE(NSM,NBFM),NUM_CP(2,*),NUM_CPSUB(*)
      REAL*8 CP_LOCAL(3,4,NXM),CPSUB(3,8,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER n1,N1OFF,n2,N2OFF,ni,ni1,ni2,nj,NKN,NOFF(2),ns1,ns2,
     '  NSOFF(2)
      REAL*8 A(6,7)

C     matrix coeffs:     linear                     cubic
      DATA A      /1.0000D0,0.0000D0,  1.000D0,0.000D0,0.000D0,0.000D0,
     '             0.8333D0,0.1667D0,  0.500D0,0.500D0,0.000D0,0.000D0,
     '             0.6667D0,0.3333D0,  0.250D0,0.500D0,0.250D0,0.000D0,
     '             0.5000D0,0.5000D0,  0.125D0,0.375D0,0.375D0,0.125D0,
     '             0.3333D0,0.6667D0,  0.000D0,0.250D0,0.500D0,0.250D0,
     '             0.1667D0,0.8333D0,  0.000D0,0.000D0,0.500D0,0.500D0,
     '             0.0000D0,1.0000D0,  0.000D0,0.000D0,0.000D0,1.000D0/

      IF(ISPL.EQ.1) THEN
        DO nj=1,NJT
          N1OFF=NUM_CP(1,nj)-2
          N2OFF=NUM_CP(2,nj)-2
          NUM_CPSUB(1)=7
          NUM_CPSUB(2)=7
          DO ns2=1,NUM_CPSUB(2)
            DO ns1=1,NUM_CPSUB(1)
              CPSUB(nj,ns1,ns2)=0.0D0
              DO n2=1,NUM_CP(2,nj)
                DO n1=1,NUM_CP(1,nj)
                  CPSUB(nj,ns1,ns2)=CPSUB(nj,ns1,ns2)+
     '              CP_LOCAL(nj,n1,n2)*A(n1+N1OFF,ns1)*A(n2+N2OFF,ns2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE IF(ISPL.EQ.2) THEN
        DO ni=1,2
          NOFF(ni)=NQE(1,ni)-1-MBSPL(ni)
C         find the knots in the refined knot vector where ne starts
          NSOFF(ni)=1
          DO NKN=1,NTKN(ni)-1
            IF(BKNOT(NKN,ni).LT.BKNOT(NKN+1,ni)) THEN
              IF(NKN.EQ.NQE(1,ni))GOTO 20
              NSOFF(ni)=NSOFF(ni)+6
            ELSE
              NSOFF(ni)=NSOFF(ni)+1
            ENDIF
          ENDDO
 20       NSOFF(ni)=NSOFF(ni)-(MBSPL(ni)+1)/2-1
        ENDDO
        NUM_CPSUB(1)=7
        NUM_CPSUB(2)=7
        DO ns2=1,NUM_CPSUB(2)
          DO ns1=1,NUM_CPSUB(1)
            DO nj=1,NJT
              CPSUB(nj,ns1,ns2)=0.0D0
              DO n2=1,NUM_CP(2,nj)
                DO n1=1,NUM_CP(1,nj)
                  CPSUB(nj,ns1,ns2)=CPSUB(nj,ns1,ns2)+CP_LOCAL(nj,n1,n2)
     '            *ALPHA(n1+NOFF(1),ns1+NSOFF(1),1)
     '            *ALPHA(n2+NOFF(2),ns2+NSOFF(2),2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF(DOP) THEN
        DO nj=1,NJT
          DO ni2=1,NUM_CPSUB(2)
            WRITE(OP_STRING,
     '        '('' CPSUB('',I1,'',..,'',I1,''): '',7E11.3)')
     '        nj,ni2,(CPSUB(nj,ni1,ni2),ni1=1,NUM_CPSUB(1))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('CPCP2')
      RETURN
 9999 CALL ERRORS('CPCP2',ERROR)
      CALL EXITS('CPCP2')
      RETURN 1
      END


      SUBROUTINE FEMREINI(IBT,IDO,INP,NDET,NDIPOLES,NEELEM,NEL,NFF,NHP,
     '  NKE,NKEF,NKH,NLL,NLLINE,NNF,NNL,NONY,NP_INTERFACE,NPL,
     '  NPNE,NPNODE,NPNY,NVHE,NVHP,NVJE,NW,NYNE,NYNP,DET,DL,FEXT,
     '  PE,XA,XP,YP,ZA,ZC,ZP,FIX,ERROR,*)

C#### Subroutine: FEMREINI
C###  Description:
C###    FEMREINI initialises all arrays and variables.

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:anal00.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:fsklib.inc'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NDET(NBFM,0:NNM),NDIPOLES(NRM,NXM),NEELEM(0:NE_R_M,0:NRM),
     '  NEL(0:NELM,NLM),
     '  NFF(6,NEM),NHP(NPM,0:NRM,NXM),NKE(NKM,NNM,NBFM,NEM),
     '  NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM),
     '  NLL(12,NEM),
     '  NLLINE(0:NL_R_M,0:NRM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),NP_INTERFACE(0:NPM,0:3),
     '  NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM),
     '  NVJE(NNM,NBFM,NJM,NEM),
     '  NW(NEM,3),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 DET(NBFM,0:NNM,NGM,6),DL(3,NLM),FEXT(NIFEXTM,NGM,NEM),
     '  PE(2,NEM),
     '  XA(NAM,NJM,NEM),XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),ZC(NJM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER I,iw,J,k,MAXNITW,na,nb,nc,ne,ng,nh,ni,niy,
     '  nj,nk,nl,nn,no,np,nr,nrc,nv,nx,ny

      CALL ENTERS('FEMREINI',*9999)

      DO nr=0,NRM
        DO ne=0,NE_R_M
          NEELEM(ne,nr)=0
        ENDDO
        DO np=0,NP_R_M
          NPNODE(np,nr)=0
        ENDDO
        NET(nr)=0
        NPT(nr)=0
      ENDDO
      DO nx=1,NXM
        DO nr=1,NRM
          DO nrc=1,NRCM
            DO ny=1,NYM
              DO no=0,NOYM
                NONY(no,ny,nrc,nr,nx)=0
              ENDDO
            ENDDO
          ENDDO
        ENDDO !nr
      ENDDO !nx

      DO nb=1,NBFM
        NFE(nb)=0
        DO ni=1,NIM
          IBT(1,ni,nb)=0
          IBT(2,ni,nb)=0
          IBT(3,ni,nb)=0
          DO nn=1,NNM
            INP(nn,ni,nb)=1
            DO nk=1,NKM
              IDO(nk,nn,ni,nb)=1
            ENDDO
          ENDDO
        ENDDO
        DO nn=1,NNM
          DO nk=1,NKM
            IDO(nk,nn,0,nb)=0
          ENDDO
        ENDDO
        DO J=1,12
          DO I=1,4
            NNL(I,J,nb)=0
          ENDDO
        ENDDO
        DO J=1,6
          DO I=0,17
            NNF(I,J,nb)=0
          ENDDO !i
          DO i=1,16
            DO k=0,4
              NKEF(k,i,j,nb)=0
            ENDDO !k
          ENDDO !i
        ENDDO !j
      ENDDO !nb

      DO ne=1,NEM
        DO I=1,6
          NFF(I,ne)=0
        ENDDO
        DO I=1,12
          NLL(I,ne)=0
        ENDDO
      ENDDO

      DO nl=1,NLM
        DO I=1,5
          DO nj=0,3
            NPL(I,nj,nl)=0
          ENDDO
        ENDDO
        DO I=0,NELM
          NEL(I,nl)=0
        ENDDO
        DO I=1,3
          DL(I,nl)=0.0d0
        ENDDO
      ENDDO
      DO nl=0,NL_R_M
        DO nr=0,NRM
          NLLINE(nl,nr)=0
        ENDDO
      ENDDO

      DO ne=1,NEM
        DO nb=1,NBFM
          DO nn=1,NNM
            DO nk=1,NKM
              NKE(nk,nn,nb,ne)=0
            ENDDO
            NPNE(nn,nb,ne)=0
          ENDDO
        ENDDO
      ENDDO

      DO nb=1,NBFM
        DO nn=1,NNM
          NKT(nn,nb)=0
        ENDDO
        NNT(nb)=0
        NKT(0,nb)=0
        NIT(nb)=0
        NST(nb)=0
      ENDDO

      DO np=1,NPM
        DO nj=1,NJM
          DO nv=1,NVM
            DO nk=1,NKM
              XP(nk,nv,nj,np)=0.0d0
            ENDDO
          ENDDO
        ENDDO
        DO nc=1,NCM
          DO nh=1,NHM
            DO nv=1,NVM
              DO nk=1,NKM
                ZP(nk,nv,nh,np,nc)=0.0d0
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO ne=1,NEM
        DO nj=1,NJM
          DO na=1,NAM
            XA(na,nj,ne)=0.0d0
          ENDDO !na
        ENDDO !nj
        DO nj=1,NJM
          ZC(nj,ne)=0.0d0
        ENDDO !nj
        DO nc=1,NCM
          DO nh=1,NHM
            DO na=1,NAM
              ZA(na,nh,nc,ne)=0.0d0
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO nr=1,NRM
        DO nc=1,NCM
          DO nrc=0,NRCM
            DO np=1,NPM
              DO nh=1,NHM
          DO nv=1,NVM
        	  DO nk=1,NKM
      		    NYNP(nk,nv,nh,np,nrc,nc,nr)=0
      		  ENDDO !nk
      		ENDDO !nv
      	      ENDDO !nh
      	    ENDDO !np
      	  ENDDO !nrc
        ENDDO !nc
      ENDDO !nr
      DO ne=1,NEM
        DO nc=1,NCM
          DO nrc=1,NRCM
            DO nh=1,NHM
              DO na=1,NAM
                NYNE(na,nh,nrc,nc,ne)=0
              ENDDO !na
            ENDDO !nh
          ENDDO !nrc
        ENDDO !nc
      ENDDO !ne

      DO nrc=0,NRCM
        DO ny=1,NYM
          DO i=0,6
            NPNY(i,ny,nrc)=0
          ENDDO !i
        ENDDO !ny
      ENDDO !nrc

      IF(USE_BEM.NE.0) THEN
        DO nb=1,NBFM
          DO nn=0,NNM
            DO ng=1,NGM
              DO I=1,6
                DET(nb,nn,ng,I)=1.0d0
              ENDDO
            ENDDO
            NDET(nb,nn)=1
          ENDDO
        ENDDO
      ENDIF

      DO ne=1,NEM
        DO ng=1,NGM
          DO ni=1,NIFEXTM
            FEXT(ni,ng,ne)=0.0d0
          ENDDO !ni
        ENDDO !ng
      ENDDO !ne

      DO ne=1,NEM
        NW(ne,1)=1
        PE(1,ne)=0.0d0
        PE(2,ne)=0.0d0
        DO nb=1,NBFM
          DO nj=1,NJM
            DO nn=1,NNM
              NVJE(nn,nb,nj,ne)=0
            ENDDO
          ENDDO
          DO nh=1,NHM
            DO nn=1,NNM
              NVHE(nn,nb,nh,ne)=0
            ENDDO !nn
          ENDDO !nh
        ENDDO
      ENDDO

      DO nx=1,NXM
        DO niy=1,NIYM
          DO ny=1,NYM
            YP(ny,niy,nx)=0.0d0
          ENDDO !ny
        ENDDO !niy
      ENDDO! nx

      DO nx=1,NXM
        DO i=1,NIYFIXM
          DO ny=1,NYM
            FIX(ny,i,nx)=.FALSE.
          ENDDO !ny
        ENDDO !i
      ENDDO !nx

      DO nr=1,NRM
        DO np=1,NPM
          DO nc=1,NCM
            DO nh=1,NHM
              NKH(nh,np,nc)=0
              NVHP(nh,np,nc)=0
            ENDDO
          ENDDO
          DO nx=1,NXM
            NHP(np,nr,nx)=0
          ENDDO
        ENDDO
      ENDDO

      DO np=0,NP_R_M
        DO I=0,3
          NP_INTERFACE(np,I)=0
        ENDDO
      ENDDO

      DO nr=1,NRM
        DO nx=1,NXM
          NDIPOLES(nr,nx)=0
        ENDDO
      ENDDO

      DO nr=1,NRT
        ANAL_CHOICE(nr)=0
      ENDDO

      MAXNITW=10 ! max # windows that can be parsed
      DO iw=1,MAXNITW
        IF(IWKG(iw).EQ.1) THEN
          IWKG(iw)=0
          CALL CLOSE_WS(iw,ERROR,*9999)
        ENDIF
      ENDDO !iw

      CALL EXITS('FEMREINI')
      RETURN
 9999 CALL ERRORS('FEMREINI',ERROR)
      CALL EXITS('FEMREINI')
      RETURN 1
      END


      SUBROUTINE GAUSS3(IBT,nb,NGAP,PG,WG,XIG,ERROR,*)

C#### Subroutine: GAUSS3
C###  Description:
C###    GAUSS3 defines the Gaussian quadrature coords XIGG and weights 
C###    WG and evaluates the basis function Gauss point array PG for
C###    B-spline tensor product type basis function nb.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),nb,NGAP(NIM)
      REAL*8 PG(NSM,NUM,NGM),WG(NGM),XIG(NIM,NGM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,i1,i2,i3,j,k,ng,ng1,ng2,ng3,ni,nk,nu
      REAL*8 D(5,5),PSI3,W(5,5),XI(3),XIGG(5,5,5,3)

      DATA D/5*0.0D0,-0.288675134594813D0,0.288675134594813D0,3*0.0D0,
     '      -0.387298334620741D0,0.0D0,0.387298334620741D0,2*0.0D0,
     '      -0.430568155797026D0, -0.169990521792428D0,
     '       0.169990521792428D0,  0.430568155797026D0,0.0D0,
     '      -0.453089922969332D0, -0.269234655052841D0,0.0D0,
     '       0.269234655052841D0,0.453089922969332D0/
      DATA W/1.0D0,4*0.0D0,0.5D0,0.5D0,3*0.0D0,0.277777777777778D0,
     '       0.444444444444444D0,0.277777777777778D0,2*0.0D0,
     '       0.173927422568727D0,0.326072577431273D0,
     '       0.326072577431273D0,0.173927422568727D0,0.0D0,
     '       0.118463442528094D0,0.239314335249683D0,
     '       0.284444444444444D0,0.239314335249683D0,
     '       0.118463442528094D0/

      CALL ENTERS('GAUSS3',*9999)
      ng1=NGAP(1)
      ng2=1
      ng3=1
      IF(NIT(nb).GT.1) ng2=NGAP(2)
      IF(NIT(nb).GT.2) ng3=NGAP(3)
      DO 20 k=1,ng3
      DO 20 j=1,ng2
      DO 20 i=1,ng1
        XIGG(i,j,k,1)=0.5D0+D(i,ng1)
        XIGG(i,j,k,2)=0.5D0+D(j,ng2)
        XIGG(i,j,k,3)=0.5D0+D(k,ng3)
        ng=i+(j-1+(k-1)*ng2)*ng1
        WG(ng)=W(i,ng1)*W(j,ng2)*W(k,ng3)
        DO ni=1,NIT(nb)
          XI(ni)=XIGG(i,j,k,ni)
          XIG(ni,ng)=XI(ni)
        ENDDO
        nk=0
        DO 10 i3=1,IBT(2,3)+1
        DO 10 i2=1,IBT(2,2)+1
        DO 10 i1=1,IBT(2,1)+1
          nk=nk+1
        DO 10 nu=1,NUT(nb)
          PG(nk,nu,ng)=PSI3(i1,i2,i3,NIT(nb),nu,XI)
 10     CONTINUE
 20   CONTINUE

      NBASEF(nb,0)=1 !number of children in family
      NBASEF(nb,1)=nb !parent number
      NFBASE(1,nb)=nb !family number of global basis number nbbem
      NFBASE(2,nb)=1  !local child number in family

      CALL EXITS('GAUSS3')
      RETURN
 9999 CALL ERRORS('GAUSS3',ERROR)
      CALL EXITS('GAUSS3')
      RETURN 1
      END
                  

      SUBROUTINE GAUSSPWB(ALIM,BLIM,N,NNSING,WEIGHT,ABSCIS,ICALL,
     '  ERROR,*)

C#### Subroutine: GAUSSPWB
C###  Description:
C###    GAUSSPWB is called from GAUSS10.
C###    This subroutine returns the weights and abscissae for a
C###    Gauss-Legendre integration scheme of order N. If the scheme of
C###    order N is not available then the next highest order scheme is
C###    returned and N is changed to the appropriate value. Upon 
C###    entrance A,B contain the lower and upper limits of the 
C###    integration interval, unchanged on exit. On exit the arrays 
C###    ABSCIS and WEIGHT contain the abscissae and weights
C###    respectively. NNSING refers to the node number (for singular
C###    basis functions only)

C**** ICALL=10 if called from Gauss10
C**** This subroutine could be easily modified
C**** to generate quadrature schemes by altering the first param in the
C**** call to the NAG routine. (see the NAG manuals)
C**** NNSING indicates location of physical singularity.  If NNSING=1
C**** then singularity is at first local node and if NNSING=NNT then
C**** singularity is at last local node.

      IMPLICIT NONE
      INTEGER NGMAX
      PARAMETER (NGMAX = 64)
!     Parameter List
      INTEGER N,NNSING,ICALL
      REAL*8 ALIM,BLIM,ABSCIS(NGMAX),WEIGHT(NGMAX)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,IFAIL,IND(NGMAX)
      REAL*8 ADP,BDP,ABSCISDP(NGMAX),WEIGHTDP(NGMAX)
C ***  IND(n) An array telling which gauss scheme will be generated.
C ***  i.e. if a scheme is requested that is not present then the next
C ***  highest available scheme is used.
      EXTERNAL D01BAZ

      DATA IND/1,2,3,4,5,6,2*8,2*10,2*12,2*14,2*16,4*20,4*24,8*32,16*48,
     '      16*64/

      CALL ENTERS('GAUSSPWB',*9999)
      ADP=ALIM
      BDP=BLIM
      IFAIL=1
      IF(ICALL.EQ.10)THEN
C ***   Standard Gauss quadrature (low, mid and high order schemes)

        CALL D01BBF(D01BAZ,ADP,BDP,0,IND(N),WEIGHTDP,ABSCISDP,IFAIL)

      ELSEIF(ICALL.EQ.11)THEN
C ***   Calculates the quadrature scheme for an integrand containing an
C ***   r**(-1/2) type singularity
C ***     General integrand for this NAG routine is of the form
C ***        (b-xi)**c * (xi-a)**d * f(xi)
C ***     Want to have an integrand of the form
C ***        (1-xi)**(-0.5) * f(xi)    or     xi**(-0.5) * f(xi)
C ***       for the squareroot singularity at the corner of the cavity.
C ***     The type of singularity depends on the value of NNSING

        IF(NNSING.EQ.1)THEN !Second type of singularity
          CALL D01BCF(1,ADP,BDP,0.0d0,-0.5d0,IND(N),WEIGHTDP,
     '      ABSCISDP,IFAIL)
        ELSE !First type of singularity
          CALL D01BCF(1,ADP,BDP,-0.5d0,0.0d0,IND(N),WEIGHTDP,
     '      ABSCISDP,IFAIL)
        ENDIF
      ENDIF

C *** This next line changes the order of the gauss point scheme wanted
C *** to the order of the gauss scheme that is given, if the required
C *** scheme is not avaliable.
      N=IND(N)
      DO I=1,N
        ABSCIS(I)=ABSCISDP(N+1-I)
        WEIGHT(I)=WEIGHTDP(N+1-I)
      ENDDO

      CALL EXITS('GAUSSPWB')
      RETURN
9999  CALL ERRORS('GAUSSPWB',ERROR)
      CALL EXITS('GAUSSPWB')
      RETURN 1
      END


      SUBROUTINE GLOBALH(IBT,INP,NBH,NBJ,NEELEM,NHE,NHP,NJP,
     '  NKH,NKJ,NPE,NP_INTERFACE,NPNODE,NPNY,NQE,nr,NVNE,NVNP,NW,nx,NXI,NYNE,
     '  NYNP,XP,ERROR,*)

C**** Calc.s nodal arrays from element arrays NHE(ne),NBH(nc,nh,ne)
C**** and NPE(nn,nb,ne):
C**** Arrays extended to have nc dependence by AJP 17-12-91.
C#### NHP(np) is number of dependent variables at node np.
C#### NKH(nc,nh,np) is number of derivatives for dependent variable nh
C###  and equation number nc at node np.
C###    For FEM problems nc=1 (only one equation generated per dep var).
C###    For BEM problems equations are generated for the dependent 
C###    variable and its normal derivative(s).
C###    For BEM problems :
C###               nc=1 corresponds to the current dep. variable
C###               nc=2       "      "  its normal derivative
C###               (on the lowest corner elem if np is at a corner)
C###               (nc=3 corresponds to the normal derivative on the
C###                next highest corner elem if appropriate and
C###                nc=4 corresponds to the normal derivative on the
C###                highest corner element if appropriate)
C#### NYNP(nk,nv,nh,np,nc,nx) for problem nx is mapping between 
C###   nk,nv,nh,np,nc and ny
C###  NYNP(nk,nv,nh,np,nc) is the global mapping between 
C###   nk,nc,nh,np and ny for use in coupled hermite problems where 
C###   derivative equations are only generated for one region at nodes
C###   on an interface.
C#### NYNE(na,nc,nh,ne) is mapping between na,nh,nc,ne and ny
C#### NVNP(0,nh,np,nc) is the mapping between nh,np and nr 
C###   and the maximum nc value.  This array is set up in
C###   IPEQUA or DEEQUA and possibly modified in DECORN.
C###   NVNP(nv=1..nvmax,nh,np,nc) identifies the element associated with
C###   the current nv value at node np. Thus it identifies which element
C###   a normal derivative is associated with. 
C#### NYT(nc,nr,nx) is the total number of ny values for nr and problem
C###   type nx.
C**** New AJP 13-1-94
C**** KTYP93(nc) indicates whether cross derivatives are to be included in the
C**** dependent variable interpolation if bicubic hermite elements are used
C**** (see IPINI9).  If they are not (KTYP93(nc)=1) then NKH is altered from
C**** 4 to 3 for the appropriate nodes.
C**** AJP 7-6-94.  Adjusted NKH back to 4 for both ktyp93 cases. 
C#### NPNY is the reverse mapping to NYNP
C###  NPNY(1,ny)=nk
C###  NPNY(2,ny)=nv
C###  NPNY(3,ny)=nh
C###  NPNY(4,ny)=np
C###  NPNY(5,ny)=nc
                 
      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:bem000.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp90.cmn'
!     Parameter List
      INTEGER IBT(2,NIM,*),INP(NNM,NIM,*),NBH(NCM,NHM,*),NBJ(NCM,NJM,*),
     '  NCNE(NNM,NBM,NEM,*),NCNP(0:NCM,NJM,NPM,*),
     '  NEELEM(0:NEM,0:*),NHE(*),NHP(*),NJP(*),
     '  NKH(NCM,NHM,*),NKJ(NCM,NJM,*),NPE(NNM,NBM,*),NPNODE(0:NPM,0:*),
     '  NP_INTERFACE(0:NPM,0:*),NPNY(5,*),NQE(NSM,NBM,*),nr,NW(NEM,*),
     '  nx,NXI(-NIM:NIM,0:*),
     '  NYNE(NAM,NVM,NHM,*),NYNP(NKM,NVM,NHM,NPM,*)
      REAL*8 XP(NKM,NCM,NJM,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,na,nb,nc,ncr,ne,nh,nj,nk,NKK,nn,NOELEM,
     '  NONODE,np,NR1,NR2,ns,nv,NY(0:16),NY_MAX,NY_ZERO
      LOGICAL ALL_DERIVATIVES,NO_EQUATION
                 
      CALL ENTERS('GLOBALH',*9999)

      nv=1 ! temporary cpb 23/11/94

C *** Set up NHP and NKH arrays for non-FE50 problems
C *** (FE50 problems set up NHP & NKH in IPINI5)
      IF(ITYP1(nr).NE.5) THEN

        DO NONODE=1,NPNODE(0,nr)
          np=NPNODE(NONODE,nr)
          DO nc=1,NCNP(0,1,np,nr) !Added AJP 21-1-92
            DO NOELEM=1,NEELEM(0,nr)
              ne=NEELEM(NOELEM,nr)
              nb=NBH(1,1,ne)
              ns=0
!news AJP 20-1-94
              IF(NB.GT.0.AND.NB.LE.NBM)THEN !AJP 17-6-94
!               Added so when this is called from Refine the 
!               arrays can stillbe set up even if NBH hasn't 
!               been defined yet
!newe
                DO nn=1,NNT(nb)
                  IF(NPE(nn,nb,ne).EQ.NP) THEN
                    IF(NHE(ne).GT.NHP(np)) NHP(np)=NHE(ne)
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' np='',I4,'' ne='',I4,'
     '                  //''' nn='',I2,'' NHP(np)='',I1)') np,ne,nn,
     '                  NHP(np)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                    DO nh=1,NHE(ne)
                      nb=NBH(nc,nh,ne)
                      IF(ITYP2(nr).EQ.1.AND.NW(ne,1).EQ.3) THEN 
                        !linear elastic beam requires derivs for all nh
                        NKH(nc,nh,np)=2  
                      ELSE IF(ITYP2(nr).EQ.1.AND.NW(ne,1).EQ.6) THEN 
                        !lin. elast plate may require derivs for all nh
                        IF(NJT.EQ.3.AND.NBH(nc,3,ne).GT.0) THEN 
                          !3D & bending incl
                          IF(NKT(nn,NBH(nc,3,ne)).EQ.4) THEN
                            NKH(nc,nh,np)=4  !transv. defl.n is bicubic
                          ELSE
                            IF(NKT(nn,nb).GT.NKH(nc,nh,np)) 
     '                        NKH(nc,nh,np)=NKT(nn,nb)
                          ENDIF
                        ELSE
                          IF(NKT(nn,nb).GT.NKH(nc,nh,np)) 
     '                      NKH(nc,nh,np)=NKT(nn,nb)
                        ENDIF
                      ELSE
                        IF(NKT(nn,nb).GT.NKH(nc,nh,np)) 
     '                    NKH(nc,nh,np)=NKT(nn,nb)
                        IF(JTYP2.EQ.1) THEN
                          DO nk=1,NKT(nn,nb)
                            ns=NS+1
                            NKK=NQE(ns,nb,ne)
                            IF(NKK.GT.NKH(nc,nh,np)) NKH(nc,nh,np)=NKK
                          ENDDO !End of NK loop
                        ENDIF
                      ENDIF
                      IF(DOP) THEN
                        WRITE(OP_STRING,'('' nh='',I1,'' nb='',I2,'
     '                     //''' nc='',I1,'' nkh(nc,nh,np)='',I2)') 
     '                     nh,nb,nc,NKH(nc,nh,np)
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF
                    ENDDO ! NH
                    GO TO 10 !don't need to complete loops once NP found
                  ENDIF ! NPE decision
                ENDDO ! NN
              ENDIF ! NB decision
 10           CONTINUE
            ENDDO ! NE
          ENDDO ! NC
        ENDDO ! NP
      ENDIF 

!new MPN: Separate out BEM/coupled problems from everything else
      IF((ITYP1(nr).NE.9).AND. !not BEM
     '   (KTYP90.EQ.0)) THEN   !not coupled FE/BE problems

!old MPN 15-Nov-94: moved to GLOBALJ
!        !Compute NCNE and NCNP mapping arrays
!        CALL CALC_NCNEP(IBT,INP,NBH,NBJ,NCNE,NCNP,NEELEM,NHE,NHP,
!     '    NPE,NPNODE,nr,NXI,XP,ERROR,*9999)      

        ncr=1 ! Temporary cpb 13/11/94

C ***   Set NKH(2..,nh,np)=NKH(1,nh,np)  (Needs to be generalised)
        DO NONODE=1,NPNODE(0,nr)
          np=NPNODE(NONODE,nr)
          DO nh=1,NHP(np)
            DO nc=2,NCNP(0,nh,np,nr)
              NKH(nc,nh,np)=NKH(1,nh,np)
            ENDDO ! NC
          ENDDO ! NH
        ENDDO ! NP

C ***   Compute NYNP and NPNY mapping arrays
        NY_ZERO=0
        DO NONODE=1,NPNODE(0,NR)
          NP=NPNODE(NONODE,NR)
          DO nh=1,NHP(NP) !Set up NYNP and NPNY
            DO nv=1,NCNP(0,nh,np,NR)
C ***         Set NKH(2..,nh,np)=NKH(1,nh,np)  (Could generalise later)
              IF(NC.GT.1) NKH(nc,NH,NP)=NKH(1,NH,NP)
              DO nk=1,NKH(nc,NH,NP) !global derivatives
                NY_ZERO=NY_ZERO+1
                NYNP(nk,nv,nh,np,nc)=NY_ZERO
                NPNY(1,ny_zero)=nk
                NPNY(2,ny_zero)=nv
                NPNY(3,ny_zero)=nh
                NPNY(4,ny_zero)=np
                NPNY(5,ny_zero)=nc
              ENDDO !End of NK loop
            ENDDO !End of nv loop
          ENDDO !End of NH loop
        ENDDO !End of nonode loop (i.e. NP loop)

C!!!     Above code could be replaced with following call but NRLIST
C!!!     must be passed through
c        CALL CALC_NYNP(NCNP,NHP,NKH,NPNODE,NPNY,NRLIST,nx,NYNP,
c     '    ERROR,*9999)

c cpb 23/11/94 this should go to calc_nynp

C ***   Compute NYNE mapping array
        DO NOELEM=1,NEELEM(0,nr)
          ne=NEELEM(NOELEM,nr)
          DO nh=1,NHE(ne)
C ***       Set NBH(2..,nh,ne)=NBH(1,nh,ne)  (Needs to be generalised)
            nb=NBH(1,nh,ne)
            IF(NB.GT.0.AND.NB.LE.NBM)THEN !AJP 17-6-94
              DO nn=1,NNT(nb)
                DO nc=2,NCNE(nn,nb,ne,nr)
                  NBH(nc,nh,ne)=NBH(1,nh,ne)
                ENDDO
              ENDDO
! NOTE:     This will need to be generalised if nv dependence is 
!           required for element dof - MPN
              DO na=1,NAT(nb)
                NY_ZERO=NY_ZERO+1
                NYNE(na,nv,nh,ne)=NY_ZERO
              ENDDO
            ENDIF
          ENDDO ! NH
        ENDDO ! NE

        NYT(1,nr,nx)=NY_ZERO
        CALL ASSERT(NYT(1,nr,nx).LE.NYM,'>Increase NYMX',ERROR,*9999)
!newe
      ELSE !BEM or coupled FE/BE problems

        !Set up NY(nc) to be the maximum ny in the
        !regions with smaller nr values
        NY_MAX=0
        DO nc=0,NCM              
          NY(nc)=0
        ENDDO
        DO NR2=1,NR-1
          DO NONODE=1,NPNODE(0,NR2) 
            np=NPNODE(NONODE,NR2)
            DO ncr=0,NCRM
              DO nh=1,NHM
                DO nv=1,NVM
                  DO nk=1,NKM 
                    NY_MAX=NYNP(nk,nv,nh,np,nc)
                    IF(NY(nc).LT.NY_MAX) NY(nc)=NY_MAX
                  ENDDO
                ENDDO
              ENDDO
            ENDDO ! ncr
          ENDDO
        ENDDO
        DO nc=1,NCM
          NYT(nc,nr,nx)=NY(nc) !Initial value
        ENDDO
        NY_ZERO=NY(0)

        DO NONODE=1,NPNODE(0,nr) !to compute new arrays
          np=NPNODE(NONODE,nr)
          NO_EQUATION=.FALSE.

          IF(KTYP91.EQ.2) THEN 
            !Derivative equations generated in one region only
            IF(NP_INTERFACE(np,0).GT.1)THEN
              !Node is an interface node
              NR1=NP_INTERFACE(np,1) !Assumes interface is between only
              NR2=NP_INTERFACE(np,2) !two regions.
              IF((ITYP4(NR1).EQ.2).AND.(ITYP4(NR2).EQ.2).AND.
     '          (NR.NE.NR1))THEN
                !Interface between 2 be regions and not in the first
                !region
                NO_EQUATION=.TRUE.
              ELSE IF((ITYP4(NR1).EQ.1).AND.(ITYP4(NR2).EQ.2).AND.
     '          (NR.NE.NR1))THEN  
                !Interface between an fe region (nr1) and a be region
                !(nr2) and in the be region
                NO_EQUATION=.TRUE.
              ELSE IF((ITYP4(NR1).EQ.2).AND.(ITYP4(NR2).EQ.1).AND.
     '          (NR.NE.NR2))THEN
                !Interface between a be region (nr1) and an fe region
                !(nr2) and in the be region
                NO_EQUATION=.TRUE.
              ENDIF
            ENDIF
          ENDIF
!         !Check all geometric derivatives at node NP.  If they
!         !are all zero then node is part of a distorted element
!         !so ny value should not be assigned for these derivatives.
!         !Currently only distinguishing between the cases when 
!         !all derrivatives are used and no derivatives are used.
!         !Only do this if this is a BE region.
          IF(ITYP4(nr).EQ.2)THEN !BE region
            ALL_DERIVATIVES=.FALSE.
            nj=0
            DO WHILE((NJ.LT.NJP(np)).AND.(.NOT.ALL_DERIVATIVES))
              nj=NJ+1
              nk=1
              DO WHILE((NK.LT.NKJ(1,nj,np)).AND.(.NOT.ALL_DERIVATIVES))
                nk=NK+1
                IF(DABS(XP(nk,1,nj,np)).GT.RDELTA)ALL_DERIVATIVES=.TRUE.
              ENDDO
            ENDDO
          ELSE
            ALL_DERIVATIVES=.TRUE.
          ENDIF

          DO nh=1,NHP(np) !Set up NYNP and NPNY
            DO nc=1,NCNP(0,nh,np,nr)
              IF(.NOT.ALL_DERIVATIVES)THEN
                NKH(nc,nh,np)=1
              ENDIF
! AJP 24-6-94 KTYP93 is not region dependent.  What happens when we set
! up nynp for a coupled 3d FE-BE problem - we will miss out the cross
! derivatives in the FE region as well.  Is this what we want ??
              DO nk=1,MAX(NKH(nc,nh,np)-KTYP93(nc),1) !global derivatives
                IF(NC.EQ.1)THEN
                  IF(NO_EQUATION)THEN
                    IF(NK.EQ.1)THEN
                      NY_ZERO=NY_ZERO+1
                    ENDIF
                  ELSE
                    NY_ZERO=NY_ZERO+1
                  ENDIF

C cpb 22/11/94                  NYNP(nk,nc,nh,0,np,nr)=NY_ZERO

c PJH 10/3/94     Added to be more generally available
                  NPNY(1,ny_zero)=nk
                  NPNY(2,ny_zero)=nv
                  NPNY(3,ny_zero)=nh
                  NPNY(4,ny_zero)=np
                  NPNY(5,ny_zero)=nc
                  !Gives the row number of the equations when not all
                  !the equations at the interface will be constructed.
                ENDIF
                IF(NC.LE.2)THEN
                  NY(nc)=NY(nc)+1
                  NYNP(nk,nv,nh,np,nc)=NY(nc)
                  NPNY(1,ny(nc))=nk
                  NPNY(2,ny(nc))=nv
                  NPNY(3,ny(nc))=nh
                  NPNY(4,ny(nc))=np
                  NPNY(5,ny(nc))=nc
                ELSE
                  !NC greater than 2 corresponds to normal derivative
	          !values at other corner elements.
                  !Equations for these are
                  !assembled in the normal derivative matrix
    	          !(ie nc=2 matrix)
                  !in column(s) next to the equation for the first
                  !normal derivative.
                  NY(2)=NY(2)+1

                  ncr=1 ! Temporary cpb 13/11/94

                  NYNP(nk,nv,nh,np,nc)=NY(2)
                  NPNY(1,ny(2))=nk
                  NPNY(2,ny(2))=nv
                  NPNY(3,ny(2))=nh
                  NPNY(4,ny(2))=np
                  NPNY(5,ny(2))=nc
                ENDIF
              ENDDO !End of NK loop
            ENDDO !End of NC loop
          ENDDO !End of NH loop
        ENDDO !End of nonode loop (i.e. NP loop)

        DO NOELEM=1,NEELEM(0,nr) !Set up NYNE
          ne=NEELEM(NOELEM,nr)
          DO nh=1,NHE(ne)
C            DO nc=1,NCNP(0,nh,np,nr) !HUH??? *****what is NP??*****
            DO nc=1,NCM
              nb=NBH(nc,nh,ne)
! New AJP 4-2-94
              IF(NB.GT.0.AND.NB.LE.NBT)THEN
! Added so when this is called from Refine the arrays can still
! be set up even if NBH hasn't been defined yet
! End new
                DO na=1,NAT(NBH(nc,nh,ne))
                  IF(NC.LE.2)THEN
                    NY(nc)=NY(nc)+1
                    NYNE(na,nv,nh,ne)=NY(nc)
                  ELSE IF(NC.GT.2)THEN
                    NY(2)=NY(2)+1
                    NYNE(na,nv,nh,ne)=NY(2)
                  ENDIF
                ENDDO
              ENDIF
            ENDDO !NC loop
          ENDDO !NH loop
        ENDDO !NOELEM loop (ie NE loop)
        DO nc=1,NCT(nr)
          NYT(nc,nr,nx)=NY(nc)-NYT(nc,nr,nx) !Total # of ny vals for nr
c         NY(nc)=0
          CALL ASSERT(NYT(nc,nr,nx).LE.NYM,'>Increase NYMX',ERROR,*9999)
        ENDDO
      ENDIF !BEM or coupled FE/BE decision

      CALL EXITS('GLOBALH')
      RETURN
 9999 CALL ERRORS('GLOBALH',ERROR)
      CALL EXITS('GLOBALH')
      RETURN 1
      END


C KAT 2002-01-30: Old version
      SUBROUTINE NENXI(IBT,INP,NBJ,NEELEM,NENP,NPNE,NXI,ERROR,*)

C#### Subroutine: NENXI
C###  Description:
C###    NENXI finds elements surrounding element ne.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),INP(NNM,NIM,NBFM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),
     '  NPNE(NNM,NBFM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER MAX_LOC_ELEM
      PARAMETER (MAX_LOC_ELEM=300)
      INTEGER i,i1,i2,i3,j,LOC_ELEM_LIST(0:MAX_LOC_ELEM),
     '  MP111,MP112,MP121,MP122,MP211,MP212,MP221,
     '  MP222,nb,nb_out,ne,ne_loc,nee,ni,NITB,NITBHOLD,nn,noelem_loc,
     '  noelem,noelem1,noloc,nn_loc,np_loc,NP111,NP112,NP121,NP122,
     '  NP211,NP212,NP221,NP222,nr,nrr,NUM_COLLAPSED,XI_COLLAPSED
      LOGICAL SECTOR,ELEM_FOUND

      CALL ENTERS('NENXI',*9999)

      DO i=-NIM,NIM
        DO j=0,NEIM
          NXI(i,j,0)=0
        ENDDO
      ENDDO

      DO nr=1,NRT
C LKC 22-JUL-1999 Initialise in case NEELEM(0,nr)=0
        ne=0
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
!          IF(DOP) THEN
!!           write(*,*) 'searching for neighbour of ',ne
!            WRITE(OP_STRING,'('' Searching for neighbour of ='',I5)')
!     '        ne
!            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
!          ENDIF
                  
          DO i=-NIM,NIM
            DO j=0,NEIM
              NXI(i,j,ne)=0
            ENDDO
          ENDDO
          i1=1
          i2=1
          i3=1
          nb=NBJ(1,ne)
          nb_out=NBJ(1,ne)
          NITB=NIT(nb)
          NITBHOLD=NITB
          NP111=0
          NP211=0
          NP121=0
          NP221=0
          NP112=0
          NP212=0
          NP122=0
          NP222=0
          DO nn=1,NNT(nb)
            i1=INP(nn,1,nb)
            IF(NITB.GE.2) I2=INP(nn,2,nb)
            IF(NITB.EQ.3) I3=INP(nn,3,nb)
            IF(i1.EQ.1.AND.i2.EQ.1.AND.i3.EQ.1) NP111=NPNE(nn,nb,ne)
            IF(i1.EQ.2.AND.i2.EQ.1.AND.i3.EQ.1) NP211=NPNE(nn,nb,ne)
            IF(i1.EQ.1.AND.I2.EQ.2.AND.i3.EQ.1) NP121=NPNE(nn,nb,ne)
            IF(i1.EQ.2.AND.I2.EQ.2.AND.i3.EQ.1) NP221=NPNE(nn,nb,ne)
            IF(i1.EQ.1.AND.I2.EQ.1.AND.i3.EQ.2) NP112=NPNE(nn,nb,ne)
            IF(i1.EQ.2.AND.I2.EQ.1.AND.i3.EQ.2) NP212=NPNE(nn,nb,ne)
            IF(i1.EQ.1.AND.I2.EQ.2.AND.i3.EQ.2) NP122=NPNE(nn,nb,ne)
            IF(i1.EQ.2.AND.I2.EQ.2.AND.i3.EQ.2) NP222=NPNE(nn,nb,ne)
          ENDDO
          IF(NIM.GT.1) THEN
            IF(IBT(1,1,nb).EQ.3.AND.IBT(1,2,nb).EQ.3) THEN 
C Hermite simplex
              IF(NKT(1,nb).EQ.1) THEN !Apex at node 1
                NP211=NP111
              ELSE                    !Apex at node 3
                NP221=NP121
              ENDIF
C CPB 29/5/97 Adding sector elements
            ELSE
              SECTOR=.FALSE.
              NUM_COLLAPSED=0
              DO ni=1,NIT(nb)
                IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) THEN
                  SECTOR=.TRUE.                  
                  NUM_COLLAPSED=NUM_COLLAPSED+1
                  XI_COLLAPSED=ni
                ENDIF
              ENDDO
              IF(SECTOR) THEN
                IF(NUM_COLLAPSED.EQ.1) THEN
                  IF(NITB.EQ.2) THEN
                    IF(IBT(1,XI_COLLAPSED,nb).EQ.5) THEN
                      IF(XI_COLLAPSED.EQ.1) THEN
                        NP211=NP111
                      ELSE
                        NP121=NP111
                      ENDIF
                    ELSE
                      IF(XI_COLLAPSED.EQ.1) THEN
                        NP221=NP121
                      ELSE
                        NP221=NP211
                      ENDIF
                    ENDIF
                  ELSE                   
                    IF(IBT(1,XI_COLLAPSED,nb).EQ.5) THEN
                      IF(XI_COLLAPSED.EQ.1) THEN
                        IF(IBT(3,1,nb).EQ.2) THEN
                          NP211=NP111
                          NP212=NP112
                        ELSE IF(IBT(3,1,nb).EQ.3) THEN
                          NP211=NP111
                          NP221=NP121
                        ENDIF 
                      ELSE IF(XI_COLLAPSED.EQ.2) THEN
                        IF(IBT(3,1,nb).EQ.1) THEN
                          NP211=NP111
                          NP122=NP112
                        ELSE IF(IBT(3,1,nb).EQ.3) THEN
                          NP121=NP111
                          NP221=NP211
                        ENDIF 
                      ELSE
                        IF(IBT(3,1,nb).EQ.1) THEN
                          NP112=NP111
                          NP122=NP121
                        ELSE IF(IBT(3,1,nb).EQ.2) THEN
                          NP112=NP111
                          NP212=NP211
                        ENDIF 
                      ENDIF
                    ELSE
                      IF(XI_COLLAPSED.EQ.1) THEN
                        IF(IBT(3,1,nb).EQ.2) THEN
                          NP221=NP121
                          NP222=NP122
                        ELSE IF(IBT(3,1,nb).EQ.3) THEN
                          NP212=NP112
                          NP222=NP122
                        ENDIF 
                      ELSE IF(XI_COLLAPSED.EQ.2) THEN
                        IF(IBT(3,1,nb).EQ.1) THEN
                          NP221=NP211
                          NP222=NP212
                        ELSE IF(IBT(3,1,nb).EQ.3) THEN
                          NP122=NP112
                          NP222=NP212
                        ENDIF 
                      ELSE
                        IF(IBT(3,1,nb).EQ.1) THEN
                          NP212=NP211
                          NP222=NP221
                        ELSE IF(IBT(3,1,nb).EQ.2) THEN
                          NP122=NP121
                          NP222=NP221
                        ENDIF 
                      ENDIF
                    ENDIF
                  ENDIF
                ELSE IF(NUM_COLLAPSED.EQ.2) THEN
                  IF(IBT(1,XI_COLLAPSED,nb).EQ.5) THEN
                    IF(IBT(3,XI_COLLAPSED,nb).EQ.1) THEN
                      NP211=NP111
                      NP112=NP111
                      NP122=NP111
                    ELSE IF(IBT(3,XI_COLLAPSED,nb).EQ.2) THEN
                      NP211=NP111
                      NP112=NP111
                      NP212=NP111
                    ELSE
                      NP211=NP111
                      NP121=NP111
                      NP221=NP111
                    ENDIF
                  ELSE
                    IF(IBT(3,XI_COLLAPSED,nb).EQ.1) THEN
                      NP221=NP211
                      NP212=NP211
                      NP222=NP211
                    ELSE IF(IBT(3,XI_COLLAPSED,nb).EQ.2) THEN
                      NP221=NP121
                      NP122=NP121
                      NP222=NP121
                    ELSE
                      NP212=NP112
                      NP122=NP112
                      NP222=NP112
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
!         IF(DOP) THEN
!	    WRITE(OP_STRING,'('' ne='',I5,'' NPijk:'',8I5)')
!    '        ne,NP111,NP211,NP121,NP221,NP112,NP212,NP122,NP222
!           CALL WRITES(IODI,OP_STRING,ERROR,*9999)
!         ENDIF


CC LKC 8-NOV-1999
CC   Solve the problem of neighbouring elements by searching first 
CC   in the current region then search the other regions starting from 1. 
CC   The order of the region
CC   numbers may/will effect the neighbours which are found.
CC
CC
CCC LKC 17-JUL-1999 We really only want to loop over adjacent elements in 
CCC  the same region otherwise interfaces elements may be considered to 
CCC  be adjacent.
CCC
CCC news mpn 14-Sep-95
CCC         loop over all regions for adjacent elements
CCC          DO nrr=1,NRT
CCC            DO noelem1=1,NEELEM(0,nrr)
CCC              nee=NEELEM(noelem1,nrr)
CCCC old          DO nee=1,NET(nr)
CCC          DO noelem1=1,NEELEM(0,nr)
CCC            nee=NEELEM(noelem1,nr)
CC
CC CS 18-OCT-1999 LKC comments not true. We do want to know about 
CC  adjacent elements in different regions eg pressure coupling. 
CC  Reverting back for now till a solution is agreed upon.
CC rgb lets solve the problem of element looping by instead
CC     going over all elements surrounding an element
          nrr=nr
          LOC_ELEM_LIST(0)=0
          DO noloc=1,MAX_LOC_ELEM
            LOC_ELEM_LIST(noloc)=0
          ENDDO
          DO nn_loc=1,NNT(nb_out)
            np_loc=NPNE(nn_loc,nb_out,ne)
            DO noelem_loc=1,NENP(np_loc,0,nrr)
              ne_loc=NENP(np_loc,noelem_loc,nrr)
              noloc=0
              ELEM_FOUND=.FALSE.
              DO WHILE((.NOT.ELEM_FOUND).AND.
     '          (noloc.LE.LOC_ELEM_LIST(0)))
                noloc=noloc+1
                IF(LOC_ELEM_LIST(noloc).EQ.ne_loc) THEN
                  ELEM_FOUND=.TRUE.
                ENDIF
              ENDDO
              IF((.NOT.ELEM_FOUND).AND.(ne_loc.NE.ne)) THEN
                IF((LOC_ELEM_LIST(0)+1).GT.MAX_LOC_ELEM) THEN
                  ERROR='Too many maximum local elements, increase'
     '              //' MAX_LOC_ELEM'
                  GOTO 9999
                ENDIF
                LOC_ELEM_LIST(0)=LOC_ELEM_LIST(0)+1
                LOC_ELEM_LIST(LOC_ELEM_LIST(0))=ne_loc
              ENDIF
            ENDDO
          ENDDO
          
          DO noelem1=1,LOC_ELEM_LIST(0)
            nee=LOC_ELEM_LIST(noelem1)    

!            IF(DOP) THEN
!!             write(*,*) 'current nr checking elem:',nee            
!              WRITE(OP_STRING,'(''  Checking element ='',I5)')
!     '          nee
!              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
!            ENDIF
      
            i1=1
            i2=1
            i3=1
            nb=NBJ(1,nee)
            NITB=NIT(nb)
            IF(NITB.EQ.NITBHOLD) THEN
              MP111=0
              MP211=0
              MP121=0
              MP221=0
              MP112=0
              MP212=0
              MP122=0
              MP222=0
              DO nn=1,NNT(nb)
                i1=INP(nn,1,nb)
                IF(NITB.GE.2) i2=INP(nn,2,nb)
                IF(NITB.EQ.3) i3=INP(nn,3,nb)
                IF(i1.EQ.1.AND.i2.EQ.1.AND.i3.EQ.1)
     '            MP111=NPNE(nn,nb,nee)
                IF(i1.EQ.2.AND.i2.EQ.1.AND.i3.EQ.1)
     '            MP211=NPNE(nn,nb,nee)
                IF(i1.EQ.1.AND.i2.EQ.2.AND.i3.EQ.1)
     '            MP121=NPNE(nn,nb,nee)
                IF(i1.EQ.2.AND.i2.EQ.2.AND.i3.EQ.1)
     '            MP221=NPNE(nn,nb,nee)
                IF(i1.EQ.1.AND.i2.EQ.1.AND.i3.EQ.2)
     '            MP112=NPNE(nn,nb,nee)
                IF(i1.EQ.2.AND.i2.EQ.1.AND.i3.EQ.2)
     '            MP212=NPNE(nn,nb,nee)
                IF(i1.EQ.1.AND.i2.EQ.2.AND.i3.EQ.2)
     '            MP122=NPNE(nn,nb,nee)
                IF(i1.EQ.2.AND.i2.EQ.2.AND.i3.EQ.2)
     '            MP222=NPNE(nn,nb,nee)
              ENDDO
              IF(NIM.GT.1) THEN
                IF(IBT(1,1,nb).EQ.3.AND.IBT(1,2,nb).EQ.3) THEN
C                 Hermite simplex
                  IF(NKT(1,nb).EQ.1) THEN !Apex at node 1
                    MP211=MP111
                  ELSE !Apex at node 3
                    MP221=MP121
                  ENDIF
C CPB 29/5/97 Adding sector elements
                ELSE
                  SECTOR=.FALSE.
                  NUM_COLLAPSED=0
                  DO ni=1,NIT(nb)
                    IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) THEN
                      SECTOR=.TRUE.                  
                      NUM_COLLAPSED=NUM_COLLAPSED+1
                      XI_COLLAPSED=ni
                    ENDIF
                  ENDDO
                  IF(SECTOR) THEN
                    IF(NUM_COLLAPSED.EQ.1) THEN
                      IF(NITB.EQ.2) THEN
                        IF(IBT(1,XI_COLLAPSED,nb).EQ.5) THEN
                          IF(XI_COLLAPSED.EQ.1) THEN
                            MP211=MP111
                          ELSE
                            MP121=MP111
                          ENDIF
                        ELSE
                          IF(XI_COLLAPSED.EQ.1) THEN
                            MP221=MP121
                          ELSE
                            MP221=MP211
                          ENDIF
                        ENDIF
                      ELSE                   
                        IF(IBT(1,XI_COLLAPSED,nb).EQ.5) THEN
                          IF(XI_COLLAPSED.EQ.1) THEN
                            IF(IBT(3,1,nb).EQ.2) THEN
                              MP211=MP111
                              MP212=MP112
                            ELSE IF(IBT(3,1,nb).EQ.3) THEN
                              MP211=MP111
                              MP221=MP121
                            ENDIF 
                          ELSE IF(XI_COLLAPSED.EQ.2) THEN
                            IF(IBT(3,1,nb).EQ.1) THEN
                              MP211=MP111
                              MP122=MP112
                            ELSE IF(IBT(3,1,nb).EQ.3) THEN
                              MP121=MP111
                              MP221=MP211
                            ENDIF 
                          ELSE
                            IF(IBT(3,1,nb).EQ.1) THEN
                              MP112=MP111
                              MP122=MP121
                            ELSE IF(IBT(3,1,nb).EQ.2) THEN
                              MP112=MP111
                              MP212=MP211
                            ENDIF 
                          ENDIF
                        ELSE
                          IF(XI_COLLAPSED.EQ.1) THEN
                            IF(IBT(3,1,nb).EQ.2) THEN
                              MP221=MP121
                              MP222=MP122
                            ELSE IF(IBT(3,1,nb).EQ.3) THEN
                              MP212=MP112
                              MP222=MP122
                            ENDIF 
                          ELSE IF(XI_COLLAPSED.EQ.2) THEN
                            IF(IBT(3,1,nb).EQ.1) THEN
                              MP221=MP211
                              MP222=MP212
                            ELSE IF(IBT(3,1,nb).EQ.3) THEN
                              MP122=MP112
                              MP222=MP212
                            ENDIF 
                          ELSE
                            IF(IBT(3,1,nb).EQ.1) THEN
                              MP212=MP211
                              MP222=MP221
                            ELSE IF(IBT(3,1,nb).EQ.2) THEN
                              MP122=MP121
                              MP222=MP221
                            ENDIF 
                          ENDIF
                        ENDIF
                      ENDIF
                    ELSE IF(NUM_COLLAPSED.EQ.2) THEN
                      IF(IBT(1,XI_COLLAPSED,nb).EQ.5) THEN
                        IF(IBT(3,XI_COLLAPSED,nb).EQ.1) THEN
                          MP211=MP111
                          MP112=MP111
                          MP122=MP111
                        ELSE IF(IBT(3,XI_COLLAPSED,nb).EQ.2) THEN
                          MP211=MP111
                          MP112=MP111
                          MP212=MP111
                        ELSE
                          MP211=MP111
                          MP121=MP111
                          MP221=MP111
                        ENDIF
                      ELSE
                        IF(IBT(3,XI_COLLAPSED,nb).EQ.1) THEN
                          MP221=MP211
                          MP212=MP211
                          MP222=MP211
                        ELSE IF(IBT(3,XI_COLLAPSED,nb).EQ.2) THEN
                          MP221=MP121
                          MP122=MP121
                          MP222=MP121
                        ELSE
                          MP212=MP112
                          MP122=MP112
                          MP222=MP112
                        ENDIF
                      ENDIF
                    ENDIF !NUM_COLLAPSED
                  ENDIF !SECTOR
                ENDIF
              ENDIF !NIM
              IF(NITB.EQ.1) THEN
                IF((NP111.EQ.MP211).AND.(NP111.NE.0))THEN
                  NXI(-1,0,NE)=NXI(-1,0,NE)+1
                  NXI(-1,NXI(-1,0,ne),ne)=nee                 
!                  WRITE(*,*) '& ', NXI(-1,0,NE)

                ENDIF
                IF((NP211.EQ.MP111).AND.(NP211.NE.0))THEN 
                  NXI(1,0,NE)=NXI(1,0,NE)+1
                  NXI( 1,NXI(1,0,ne),ne)=nee
!                  WRITE(*,*) '&& ',NXI(1,0,NE)

                ENDIF
                IF(((NP111.EQ.MP111).AND.(nee.NE.ne))
     '            .AND.(NP111.NE.0))THEN 
                  NXI(-1,0,NE)=NXI(-1,0,NE)+1
                  NXI(-1,NXI(-1,0,ne),ne)=nee
!                  WRITE(*,*) '&&& ',NXI(-1,0,NE)

                ENDIF
                IF(((NP211.EQ.MP211).AND.(nee.NE.ne))
     '            .AND.(NP111.NE.0))THEN 
                  NXI(1,0,NE)=NXI(1,0,NE)+1
                  NXI(1,NXI(1,0,ne),ne)=nee
!                  WRITE(*,*) '&&&& ',NXI(1,0,NE)

                ENDIF
              ELSE IF(NITB.EQ.2) THEN
                IF(((NP111.EQ.MP211).AND.(NP121.EQ.MP221))
     '            .AND.(NP111.NE.0))THEN
                  NXI(-1,0,NE)=NXI(-1,0,NE)+1
                  NXI(-1,NXI(-1,0,ne),ne)=nee
                ENDIF
                IF(((NP211.EQ.MP111).AND.(NP221.EQ.MP121))
     '            .AND.(NP211.NE.0))THEN
                  NXI(1,0,NE)=NXI(1,0,NE)+1
                  NXI(1,NXI(1,0,ne),ne)=nee
                ENDIF
                IF(((NP111.EQ.MP121).AND.(NP211.EQ.MP221))
     '            .AND.(NP111.NE.0))THEN
                  NXI(-2,0,NE)=NXI(-2,0,NE)+1
                  NXI(-2,NXI(-2,0,ne),ne)=nee
                ENDIF
                IF(((NP121.EQ.MP111).AND.(NP221.EQ.MP211))
     '            .AND.(NP121.NE.0))THEN
                  NXI(2,0,NE)=NXI(2,0,NE)+1
                  NXI( 2,NXI(2,0,ne),ne)=nee
                ENDIF
              ELSE IF(NITB.EQ.3) THEN
!               IF(DOP) THEN
!	          WRITE(OP_STRING,'('' nee='',I5,'' MPijk:'',8I5)')
!    '              nee,MP111,MP211,MP121,MP221,MP112,MP212,MP122,MP222
!                 CALL WRITES(IODI,OP_STRING,ERROR,*9999)
!               ENDIF
                IF(((NP111.EQ.MP211).AND.(NP121.EQ.MP221).AND.
     '            (NP112.EQ.MP212).AND.(NP122.EQ.MP222))
     '            .AND.(NP111.NE.0))THEN
                  NXI(-1,0,ne)=1
                  NXI(-1,1,ne)=nee
                ENDIF
                IF(((NP211.EQ.MP111).AND.(NP221.EQ.MP121).AND.
     '            (NP212.EQ.MP112).AND.(NP222.EQ.MP122)) 
     '            .AND.(NP211.NE.0))THEN
                  NXI( 1,0,ne)=1
                  NXI( 1,1,ne)=nee
                ENDIF
                IF(((NP111.EQ.MP121).AND.(NP211.EQ.MP221).AND.
     '            (NP112.EQ.MP122).AND.(NP212.EQ.MP222)) 
     '            .AND.(NP111.NE.0)) THEN
                  NXI(-2,0,ne)=1
                  NXI(-2,1,ne)=nee
                ENDIF
                IF(((NP121.EQ.MP111).AND.(NP221.EQ.MP211).AND.
     '            (NP122.EQ.MP112).AND.(NP222.EQ.MP212)) 
     '            .AND.(NP121.NE.0))THEN
                  NXI(2,0,ne)=1
                  NXI(2,1,ne)=nee
                ENDIF
                IF(((NP111.EQ.MP112).AND.(NP211.EQ.MP212).AND.
     '            (NP121.EQ.MP122).AND.(NP221.EQ.MP222)) 
     '            .AND.(NP111.NE.0))THEN
                  NXI(-3,0,ne)=1
                  NXI(-3,1,ne)=nee
                ENDIF
                IF(((NP112.EQ.MP111).AND.(NP212.EQ.MP211).AND.
     '            (NP122.EQ.MP121).AND.(NP222.EQ.MP221))
     '            .AND.(NP112.NE.0))THEN
                  NXI(3,0,ne)=1
                  NXI(3,1,ne)=nee
                ENDIF
              ENDIF !NITB
            ENDIF            
      	  ENDDO ! noelem (ne)

C*** Loop over all regions (except the current)
          DO nrr=1,NRT
            IF(nrr.NE.nr) THEN !do not check the current region
               
              LOC_ELEM_LIST(0)=0
              DO noloc=1,MAX_LOC_ELEM
                LOC_ELEM_LIST(noloc)=0
              ENDDO
              DO nn_loc=1,NNT(nb_out)
                np_loc=NPNE(nn_loc,nb_out,ne)
                DO noelem_loc=1,NENP(np_loc,0,nrr)
                  ne_loc=NENP(np_loc,noelem_loc,nrr)
                  noloc=0
                  ELEM_FOUND=.FALSE.
                  DO WHILE((.NOT.ELEM_FOUND).AND.
     '              (noloc.LE.LOC_ELEM_LIST(0)))
                    noloc=noloc+1
                    IF(LOC_ELEM_LIST(noloc).EQ.ne_loc) THEN
                      ELEM_FOUND=.TRUE.
                    ENDIF
                  ENDDO
                  IF((.NOT.ELEM_FOUND).AND.(ne_loc.NE.ne)) THEN
                    IF((LOC_ELEM_LIST(0)+1).GT.MAX_LOC_ELEM) THEN
                      ERROR='Too many maximum local elements, increase'
     '                  //' MAX_LOC_ELEM'
                      GOTO 9999
                    ENDIF
                    LOC_ELEM_LIST(0)=LOC_ELEM_LIST(0)+1
                    LOC_ELEM_LIST(LOC_ELEM_LIST(0))=ne_loc
                  ENDIF
                ENDDO
              ENDDO

              DO noelem1=1,LOC_ELEM_LIST(0)
                nee=LOC_ELEM_LIST(noelem1)              
              
!                IF(DOP) THEN
!!                write(*,*) 'checking elem:',nee
!                  WRITE(OP_STRING,'(''  Checking element ='',I5)')
!     '              nee
!                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
!                ENDIF

                i1=1
                i2=1
                i3=1
                nb=NBJ(1,nee)
                NITB=NIT(nb)
                IF(NITB.EQ.NITBHOLD) THEN
                  MP111=0
                  MP211=0
                  MP121=0
                  MP221=0
                  MP112=0
                  MP212=0
                  MP122=0
                  MP222=0
                  DO nn=1,NNT(nb)
                    i1=INP(nn,1,nb)
                    IF(NITB.GE.2) i2=INP(nn,2,nb)
                    IF(NITB.EQ.3) i3=INP(nn,3,nb)
                    IF(i1.EQ.1.AND.i2.EQ.1.AND.i3.EQ.1)
     '                MP111=NPNE(nn,nb,nee)
                    IF(i1.EQ.2.AND.i2.EQ.1.AND.i3.EQ.1)
     '                MP211=NPNE(nn,nb,nee)
                    IF(i1.EQ.1.AND.i2.EQ.2.AND.i3.EQ.1)
     '                MP121=NPNE(nn,nb,nee)
                    IF(i1.EQ.2.AND.i2.EQ.2.AND.i3.EQ.1)
     '                MP221=NPNE(nn,nb,nee)
                    IF(i1.EQ.1.AND.i2.EQ.1.AND.i3.EQ.2)
     '                MP112=NPNE(nn,nb,nee)
                    IF(i1.EQ.2.AND.i2.EQ.1.AND.i3.EQ.2)
     '                MP212=NPNE(nn,nb,nee)
                    IF(i1.EQ.1.AND.i2.EQ.2.AND.i3.EQ.2)
     '                MP122=NPNE(nn,nb,nee)
                    IF(i1.EQ.2.AND.i2.EQ.2.AND.i3.EQ.2)
     '                MP222=NPNE(nn,nb,nee)
                  ENDDO
                  IF(NIM.GT.1) THEN
                    IF(IBT(1,1,nb).EQ.3.AND.IBT(1,2,nb).EQ.3) THEN
C                     Hermite simplex
                      IF(NKT(1,nb).EQ.1) THEN !Apex at node 1
                        MP211=MP111
                      ELSE !Apex at node 3
                        MP221=MP121
                      ENDIF
C CPB 29/5/97 Adding sector elements
                    ELSE
                      SECTOR=.FALSE.
                      NUM_COLLAPSED=0
                      DO ni=1,NIT(nb)
                        IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) THEN
                          SECTOR=.TRUE.                  
                          NUM_COLLAPSED=NUM_COLLAPSED+1
                          XI_COLLAPSED=ni
                        ENDIF
                      ENDDO
                      IF(SECTOR) THEN
                        IF(NUM_COLLAPSED.EQ.1) THEN
                          IF(NITB.EQ.2) THEN
                            IF(IBT(1,XI_COLLAPSED,nb).EQ.5) THEN
                              IF(XI_COLLAPSED.EQ.1) THEN
                                MP211=MP111
                              ELSE
                                MP121=MP111
                              ENDIF
                            ELSE
                              IF(XI_COLLAPSED.EQ.1) THEN
                                MP221=MP121
                              ELSE
                                MP221=MP211
                              ENDIF
                            ENDIF
                          ELSE                   
                            IF(IBT(1,XI_COLLAPSED,nb).EQ.5) THEN
                              IF(XI_COLLAPSED.EQ.1) THEN
                                IF(IBT(3,1,nb).EQ.2) THEN
                                  MP211=MP111
                                  MP212=MP112
                                ELSE IF(IBT(3,1,nb).EQ.3) THEN
                                  MP211=MP111
                                  MP221=MP121
                                ENDIF 
                              ELSE IF(XI_COLLAPSED.EQ.2) THEN
                                IF(IBT(3,1,nb).EQ.1) THEN
                                  MP211=MP111
                                  MP122=MP112
                                ELSE IF(IBT(3,1,nb).EQ.3) THEN
                                  MP121=MP111
                                  MP221=MP211
                                ENDIF 
                              ELSE
                                IF(IBT(3,1,nb).EQ.1) THEN
                                  MP112=MP111
                                  MP122=MP121
                                ELSE IF(IBT(3,1,nb).EQ.2) THEN
                                  MP112=MP111
                                  MP212=MP211
                                ENDIF 
                              ENDIF
                            ELSE
                              IF(XI_COLLAPSED.EQ.1) THEN
                                IF(IBT(3,1,nb).EQ.2) THEN
                                  MP221=MP121
                                  MP222=MP122
                                ELSE IF(IBT(3,1,nb).EQ.3) THEN
                                  MP212=MP112
                                  MP222=MP122
                                ENDIF 
                              ELSE IF(XI_COLLAPSED.EQ.2) THEN
                                IF(IBT(3,1,nb).EQ.1) THEN
                                  MP221=MP211
                                  MP222=MP212
                                ELSE IF(IBT(3,1,nb).EQ.3) THEN
                                  MP122=MP112
                                  MP222=MP212
                                ENDIF 
                              ELSE
                                IF(IBT(3,1,nb).EQ.1) THEN
                                  MP212=MP211
                                  MP222=MP221
                                ELSE IF(IBT(3,1,nb).EQ.2) THEN
                                  MP122=MP121
                                  MP222=MP221
                                ENDIF 
                              ENDIF
                            ENDIF
                          ENDIF
                        ELSE IF(NUM_COLLAPSED.EQ.2) THEN
                          IF(IBT(1,XI_COLLAPSED,nb).EQ.5) THEN
                            IF(IBT(3,XI_COLLAPSED,nb).EQ.1) THEN
                              MP211=MP111
                              MP112=MP111
                              MP122=MP111
                            ELSE IF(IBT(3,XI_COLLAPSED,nb).EQ.2) THEN
                              MP211=MP111
                              MP112=MP111
                              MP212=MP111
                            ELSE
                              MP211=MP111
                              MP121=MP111
                              MP221=MP111
                            ENDIF
                          ELSE
                            IF(IBT(3,XI_COLLAPSED,nb).EQ.1) THEN
                              MP221=MP211
                              MP212=MP211
                              MP222=MP211
                            ELSE IF(IBT(3,XI_COLLAPSED,nb).EQ.2) THEN
                              MP221=MP121
                              MP122=MP121
                              MP222=MP121
                            ELSE
                              MP212=MP112
                              MP122=MP112
                              MP222=MP112
                            ENDIF
                          ENDIF
                        ENDIF !NUM_COLLAPSED
                      ENDIF !SECTOR
                    ENDIF
                  ENDIF !NIM.GT.1
                  IF(NITB.EQ.1) THEN
                    IF((NP111.EQ.MP211).AND.(NP111.NE.0))THEN
                      NXI(-1,0,NE)=NXI(-1,0,NE)+1
                      NXI(-1,NXI(-1,0,ne),ne)=nee
!                      write(*,*) '*', NXI(-1,0,NE)
                    ENDIF
                    IF((NP211.EQ.MP111).AND.(NP211.NE.0))THEN 
                      NXI(1,0,NE)=NXI(1,0,NE)+1
                      NXI( 1,NXI(1,0,ne),ne)=nee
!                      write(*,*) '**', NXI(1,0,NE)
                    ENDIF
                    IF(((NP111.EQ.MP111).AND.(nee.NE.ne))
     '                .AND.(NP111.NE.0))THEN 
                      NXI(-1,0,NE)=NXI(-1,0,NE)+1
                      NXI(-1,NXI(-1,0,ne),ne)=nee
!                      write(*,*) '***', NXI(-1,0,NE)
                    ENDIF
                    IF(((NP211.EQ.MP211).AND.(nee.NE.ne))
     '                .AND.(NP111.NE.0))THEN 
                      NXI(1,0,NE)=NXI(1,0,NE)+1
                      NXI(1,NXI(1,0,ne),ne)=nee
!                      write(*,*) '****', NXI(1,0,NE)                     
                    ENDIF
                    
                  ELSE IF(NITB.EQ.2) THEN
                    IF(((NP111.EQ.MP211).AND.(NP121.EQ.MP221))
     '                .AND.(NP111.NE.0))THEN
                      NXI(-1,0,NE)=NXI(-1,0,NE)+1
                      NXI(-1,NXI(-1,0,ne),ne)=nee
                    ENDIF
                    IF(((NP211.EQ.MP111).AND.(NP221.EQ.MP121))
     '                .AND.(NP211.NE.0))THEN
                      NXI(1,0,NE)=NXI(1,0,NE)+1
                      NXI(1,NXI(1,0,ne),ne)=nee
                    ENDIF
                    IF(((NP111.EQ.MP121).AND.(NP211.EQ.MP221))
     '                .AND.(NP111.NE.0))THEN
                      NXI(-2,0,NE)=NXI(-2,0,NE)+1
                      NXI(-2,NXI(-2,0,ne),ne)=nee
                    ENDIF
                    IF(((NP121.EQ.MP111).AND.(NP221.EQ.MP211))
     '                .AND.(NP121.NE.0))THEN
                      NXI(2,0,NE)=NXI(2,0,NE)+1
                      NXI( 2,NXI(2,0,ne),ne)=nee
                    ENDIF
                  ELSE IF(NITB.EQ.3) THEN
!               IF(DOP) THEN
!	          WRITE(OP_STRING,'('' nee='',I5,'' MPijk:'',8I5)')
!    '              nee,MP111,MP211,MP121,MP221,MP112,MP212,MP122,MP222
!                 CALL WRITES(IODI,OP_STRING,ERROR,*9999)
!               ENDIF
                    IF(((NP111.EQ.MP211).AND.(NP121.EQ.MP221).AND.
     '                (NP112.EQ.MP212).AND.(NP122.EQ.MP222))
     '                .AND.(NP111.NE.0))THEN
                      NXI(-1,0,ne)=1
                      NXI(-1,1,ne)=nee
                    ENDIF
                    IF(((NP211.EQ.MP111).AND.(NP221.EQ.MP121).AND.
     '                (NP212.EQ.MP112).AND.(NP222.EQ.MP122)) 
     '                .AND.(NP211.NE.0))THEN
                      NXI( 1,0,ne)=1
                      NXI( 1,1,ne)=nee
                    ENDIF
                    IF(((NP111.EQ.MP121).AND.(NP211.EQ.MP221).AND.
     '                (NP112.EQ.MP122).AND.(NP212.EQ.MP222)) 
     '                .AND.(NP111.NE.0)) THEN
                      NXI(-2,0,ne)=1
                      NXI(-2,1,ne)=nee
                    ENDIF
                    IF(((NP121.EQ.MP111).AND.(NP221.EQ.MP211).AND.
     '                (NP122.EQ.MP112).AND.(NP222.EQ.MP212)) 
     '                .AND.(NP121.NE.0))THEN
                      NXI(2,0,ne)=1
                      NXI(2,1,ne)=nee
                    ENDIF
                    IF(((NP111.EQ.MP112).AND.(NP211.EQ.MP212).AND.
     '                (NP121.EQ.MP122).AND.(NP221.EQ.MP222)) 
     '                .AND.(NP111.NE.0))THEN
                      NXI(-3,0,ne)=1
                      NXI(-3,1,ne)=nee
                    ENDIF
                    IF(((NP112.EQ.MP111).AND.(NP212.EQ.MP211).AND.
     '                (NP122.EQ.MP121).AND.(NP222.EQ.MP221))
     '                .AND.(NP112.NE.0))THEN
                      NXI(3,0,ne)=1
                      NXI(3,1,ne)=nee
                    ENDIF
                  ENDIF !NITB
                ENDIF !NITB
              ENDDO ! noelem (nee)
              
            ENDIF !NR.EQ.nee
          ENDDO ! nrr

        ENDDO !nrr

C LKC 22-JUL-1999 if ne=0 then no elements in a region
C        NXI(0,0,ne)=1
C        NXI(0,1,ne)=ne
        IF(ne.NE.0) THEN
          NXI(0,0,ne)=1
          NXI(0,1,ne)=ne
        ENDIF
        
!        IF(DOP) THEN
!C$          call mp_setlock()
!          WRITE(OP_STRING,'('' NXI(...,'',I5,''):'',7(X,I5))')
!     '      ne,(NXI(I,1,ne),I=-NIM,NIM)
!          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
!C$          call mp_unsetlock()
!        ENDIF
      ENDDO ! nr

      IF(DOP) THEN
        DO nr=1,NRT
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            WRITE(OP_STRING,'('' NXI(...,'',I5,''):'',7(X,I5))')
     '        ne,(NXI(I,1,ne),I=-NIM,NIM)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
      ENDIF
      
      CALL EXITS('NENXI')
      RETURN
 9999 CALL ERRORS('NENXI',ERROR)
      CALL EXITS('NENXI')
      RETURN 1
      END

      
Module FE03
=========== 

C LC 24/2/97 Archived section from routine 
C
C
C#### Subroutine: PROJ_ORTHOG
C###  Description:
C###    <HTML>
C###    Finds the closest point on a given element to a given
C###  world coordinate (ZD).  Starts from the given XI position.  If
C###  a suitable projection cannot be found, FOUND is set to false.
C###  The scalar distance between the projection and the point is
C###  returned in DISTANCE
C
C GMH 29/10/95 Not being done correctly
C          XI(1)=XID(1,ND)
C          XI(2)=XID(2,ND)
C          XI(3)=XID(3,ND) !new AAY 30 May 1991
C          IN_ELEM = (LD(ND).EQ.NE) !new AAY 14 Aug 1991
C          ITMAX=15
C          IF(ITYP10(nr).EQ.1) THEN !rectanglar cartesian
C            IT=0
C            IF((NEW.AND.LD(ND).EQ.0).OR.
C     '        (.NOT.NEW.AND.IN_ELEM)) THEN
C !not found yet
C              SQND = SQMAX
C              CALL CLOS31(IBT,IDO,INP,IT,ITMAX,NBJ(1,NE),
C     '          SQND,XE,XI,ZD(1,ND),ERROR,*9999)
C              IF(SQND.GE.SQMAX) LD(nd)=0 !AAY 1 May 94
C            ELSE
C              XI(1)= -2.d0 !so it won't be put in xid
C            ENDIF
C          ENDIF !ityp10(nr)=1



      INTEGER FUNCTION NYPJK(NJO,nk,NKJ,np,NPO)

C**** Returns the global degree of freedom number ny associated with
C**** derivative nk, optimisation dimension NJO and global node np
C**** for fitting problems.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:fit000.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER NJO,nk,NKJ(NJM,NPM),np,NPO(0:*)
!     Local Variables
      INTEGER n1,nj1,njj,nkk,np1,ny

      ny=0
      DO n1=1,NPO(0)
        np1=NPO(n1)
        DO njj=1,NL_FIT(0,1)
          nj1=NL_FIT(njj,1)
          DO nkk=1,NKJ(nj1,np1)
            ny=ny+1
            IF((np.EQ.np1).AND.(NJO.EQ.nj1).AND.(nk.EQ.nkk)) THEN
              NYPJK=ny
              RETURN
            ENDIF
          ENDDO
        ENDDO 
      ENDDO
      NYPJK=0

      RETURN
      END


      INTEGER FUNCTION NYF(NB,NJ,NK,NP,NPNODE)

C**** Returns the global degree of freedom number NY associated with
C**** derivative NK, geometric coordinate NJ and global node NP.
C**** This is the NY number under which Fourier coefficients are stored.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER NB,NJ,NK,NP,NPNODE(0:NPM,0:*)
!     Local Variables
      INTEGER N1,NJ1,NK1,NP1,NR

      NYF=0
      DO NJ1=1,NJT
        DO NR=1,NRT
          DO N1=1,NPNODE(0,NR)
            NP1=NPNODE(N1,NR)
            DO NK1=1,NKT(0,NB)
              NYF=NYF+1
              IF((NP.EQ.NP1).AND.(NJ.EQ.NJ1).AND.(NK.EQ.NK1)) THEN
                RETURN
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END


      INTEGER FUNCTION NYPHK(nh,nk,NKH,np,NPO)

C**** Returns the global degree of freedom number ny associated with
C**** derivative nk, optimisation variable nh and global node np for 
C**** fitting problems.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:fit000.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER nh,nk,NKH(NHM,NPM,*),np,NPO(0:*)
!     Local Variables
      INTEGER n1,n1h,n1hh,n1k,N1P,ny

C CPB 19/4/94 Changed parameter list but I can't find where it is
C called from.

      ny=0
      NYPHK=0
      DO n1=1,NPO(0)
        N1P=NPO(n1)
C CPB 19/4/94 Adding NJ_FIT
C        DO n1h=NJO0,NJO1
        DO n1hh=1,NJ_FIT(0)
          n1h=NJ_FIT(n1hh)
          DO n1k=1,NKH(n1h,N1P,1)
            ny=ny+1
            IF((np.EQ.N1P).AND.(nh.EQ.n1h).AND.(nk.EQ.n1k)) THEN
              NYPHK=ny
              RETURN
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NYPHK=0

      RETURN
      END

     
      SUBROUTINE CLOS12(IT,ITMAX,NPL,nr,DL,SQ,TOL,VMAX,XP,XID,ZD)

C#### Subroutine: CLOS12
C###  Description:
C###    CLOS12 finds the XI-coordinates at the closest approach of a 
C###    line to a data point with coordinates XD using a modified 
C###    Newton algorithm.

C**** ITMAX is the maximum number of iterations.
C**** TOL is the required tolerance of the solution.
C**** VMAX denotes the maximum step length per iteration.
C**** NOTE: Works only for 1-D elements in cylindrical polar
C**** coordinates

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IT,ITMAX,NPL(20),nr
      REAL*8 DL(3),SQ,TOL,VMAX,XID,XP(NKM,NVM,NJM,NPM),ZD(NJM)
!     Local Variables
      CHARACTER ERROR*10
      INTEGER it1,it2,nj
      REAL*8 D2SQXI,D2XXI(3),DSQXI,DSQXIV,DXI,DXXI(3),DZ(3),DZXI(3),
     '  D2ZXI(3),OCO,OSO,OOX,PLXI,SQLIN,SQOLD,TOL2,V,W,X(3),XCO,
     '  XIDLIN,XOO,XSO,Z(3)

C     CALL ENTERS('CLOS12',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'(A)') ' >CLOS12 1-D cylindrical polar'
      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      TOL2=TOL**2
      DO nj=1,2
        X(nj)=PLXI(nj,NPL,nr,1,DL,XID,XP)
      ENDDO
      XOO=X(1)
      OCO=DCOS(X(2)) 

      OSO=DSIN(X(2))
      OOX=X(3)
      XCO=XOO*OCO
      XSO=XOO*OSO
      Z(1)=XCO
      Z(2)=XSO
      Z(3)=OOX
      DZ(1)=Z(1)-ZD(1)
      DZ(2)=Z(2)-ZD(2)
      SQ=DZ(1)**2+DZ(2)**2
      IF(DOP) THEN
        WRITE(OP_STRING,'(/'' XD(nj)='',2(E12.6,4X),/'' XID(ni)='','
     '    //'(E12.6,4X),/''  X(nj)='',2(E12.6,4X),'
     '    //'/''      SQ='',E12.6)')
     '    (ZD(nj),nj=1,2),XID,(X(nj),nj=1,2),SQ    
      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      DO 6 it1=1,ITMAX
        DO 2 nj=1,2
          DXXI(nj)=PLXI(nj,NPL,nr,2,DL,XID,XP)
          D2XXI(nj)=PLXI(nj,NPL,nr,3,DL,XID,XP)
    2   CONTINUE
        DZXI(1)=DXXI(1)*OCO-DXXI(2)*XSO
        DZXI(2)=DXXI(1)*OSO+DXXI(2)*XCO
        D2ZXI(1)=D2XXI(1)*OCO-D2XXI(2)*XSO
     '                               -2.d0*DXXI(1)*DXXI(2)*OSO
     '           -DXXI(2)*DXXI(2)*XCO
        D2ZXI(2)=D2XXI(1)*OSO-D2XXI(2)*XCO
     '                               +2.d0*DXXI(1)*DXXI(2)*OCO
     '           -DXXI(2)*DXXI(2)*XSO                   
                   
        DSQXI=DZXI(1)*DZ(1)+DZXI(2)*DZ(2)+DZXI(3)*DZ(3)
        D2SQXI=DZXI(1)*DZXI(1)+D2ZXI(1)*DZ(1)
     '        +DZXI(2)*DZXI(2)+D2ZXI(2)*DZ(2)
        W=1.d0
        IF(DABS(DSQXI).GT.D2SQXI*VMAX) THEN
          V=-DSIGN(VMAX,DSQXI)
        ELSE
          V=-DSQXI/D2SQXI
        ENDIF
        DSQXIV=DSQXI*V
        IF(DOP) THEN
          WRITE(OP_STRING,'(''    IT1='',I4,10X,''DSQXI(ni)='','
     '      //'(E12.6,4X),/18X,''D2SQXI(MI,ni)='',(E12.6,4X),/26X,'
     '      //'''V(ni)='',(E12.6,4X))') it1,DSQXI,D2SQXI,V
          CALL WRITES(IODI,OP_STRING,ERROR,*9999) 
        ENDIF
C ***   Performs line search
        DO 4 it2=1,ITMAX
          XIDLIN=XID+V*W
          DO 3 nj=1,2
            X(nj)=PLXI(nj,NPL,nr,1,DL,XIDLIN,XP)
    3     CONTINUE
          XOO=X(1)
          OCO=DCOS(X(2))
          OSO=DSIN(X(2))
          OOX=X(3)
          XCO=XOO*OCO
          XSO=XOO*OSO
          Z(1)=XCO
          Z(2)=XSO
          DZ(1)=Z(1)-ZD(1)         
          DZ(2)=Z(2)-ZD(2)
          SQLIN=DZ(1)**2+DZ(2)**2
          IF(DOP) THEN
            WRITE(OP_STRING,'(8X,''IT2='',I4,14X,''W='',E12.6,/21X,'
     '        //'''XIDLIN(ni)='',(E12.6,4X),/26X,''X(nj)='','
     '        //'2(E12.6,4X),/26X,''SQLIN='',E12.6)') 
     '        it2,W,XIDLIN,(X(nj),nj=1,2),SQLIN
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(SQLIN.LE.SQ) GOTO 5
          W=(DSQXIV*W*W)/(DSQXIV*W+SQ-SQLIN)
    4   CONTINUE
    5   SQOLD=SQ
        SQ=SQLIN
        DXI=V*W
        XID=XID+DXI
        IF(DOP) THEN
          WRITE(OP_STRING,'(24X,''XID(ni)='',(E12.6,4X))') XID    
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IT=it1
        IF(((DXI**2)/(1.d0+XID**2).LE.TOL2)
     '      .AND.((SQOLD-SQ)/(1.d0+SQ).LE.TOL)) GOTO 7
    6 CONTINUE

C   7 CALL EXITS('CLOS12')
    7 CONTINUE
 9999 RETURN
      END

                    
      SUBROUTINE CLOS13(IT,ITMAX,NPL,nr,DL,SQ,TOL,VMAX,XP,XID,ZD)

C#### Subroutine: CLOS13
C###  Description:
C###    CLOS13 finds the XI-coordinates at the closest approach of a 
C###    line to a data point with coordinates XD using a modified 
C###    Newton algorithm.

C**** ITMAX is the maximum number of iterations.
C**** TOL is the required tolerance of the solution.
C**** VMAX denotes the maximum step length per iteration.
C**** NOTE: Works only for 1-D elements in spherical polar
C**** coordinates

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IT,ITMAX,NPL(20),nr
      REAL*8 DL(3),SQ,TOL,VMAX,XP(NKM,NVM,NJM,NPM),XID,ZD(NJM)
!     Local Variables
      CHARACTER ERROR*10
      INTEGER it1,it2,nj
      REAL*8 D2SQXI,D2XXI(3),D2ZXI(3),DSQXI,DSQXIV,DXI,DXXI(3),DZ(3),
     '  DZXI(3),OCC,OCO,OCS,OOC,OOS,OSC,OSO,OSS,PLXI,SQLIN,SQOLD,TOL2,
     '  V,W,X(3),XCC,XCS,XIDLIN,XOC,XOO,XOS,XSC,XSS,Z(3)

C     CALL ENTERS('CLOS13',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'(A)') ' >CLOS13 1-D spherical polar' 
       CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      TOL2=TOL**2
      DO nj=1,NJT
        X(nj)=PLXI(nj,NPL,nr,1,DL,XID,XP)
      ENDDO
      XOO=X(1)
      OCO=DCOS(X(2))
      OSO=DSIN(X(2))
      OOC=DCOS(X(3))
      OOS=DSIN(X(3))
      OCC=OCO*OOC
      OCS=OCO*OOS
      OSC=OSO*OOC
      OSS=OSO*OOS
      XOC=XOO*OOC
      XOS=XOO*OOS
      XCC=XOO*OCC
      XCS=XOO*OCS
      XSC=XOO*OSC
      XSS=XOO*OSS
      Z(1)=XCC
      Z(2)=XSC
      Z(3)=XOS
      DZ(1)=Z(1)-ZD(1)
      DZ(2)=Z(2)-ZD(2)
      DZ(3)=Z(3)-ZD(3)
      SQ=DZ(1)**2+DZ(2)**2+DZ(3)**2
      IF(DOP) THEN
        WRITE(OP_STRING,'(/'' XD(nj)='',3(E12.6,4X),/''XID(ni)='','
     '    //'(E12.6,4X),/''  X(nj)='',3(E12.6,4X),'   
   
     '    //'/''     SQ='',E12.6)')
     '    (ZD(nj),nj=1,NJT),XID,(X(nj),nj=1,NJT),SQ
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF    
      DO 6 it1=1,ITMAX
        DO nj=1,NJT
          DXXI(nj)=PLXI(nj,NPL,nr,2,DL,XID,XP)
          D2XXI(nj)=PLXI(nj,NPL,nr,3,DL,XID,XP)
        ENDDO
        DZXI(1)=DXXI(1)*OCC-DXXI(2)*XSC-DXXI(3)*XCS
        DZXI(2)=DXXI(1)*OSC+DXXI(2)*XCC-DXXI(3)*XSS
        DZXI(3)=DXXI(1)*OOS            +DXXI(3)*XOC
        D2ZXI(1)=D2XXI(1)*OCC-D2XXI(2)*XSC-D2XXI(3)*XCS
     '                               -2.d0*DXXI(1)*DXXI(2)*OSC
     '           -DXXI(2)*DXXI(2)*XCC+2.d0*DXXI(2)*DXXI(3)*XSS
     '           -DXXI(3)*DXXI(3)*XCC-2.d0*DXXI(2)*DXXI(3)*OCS
        D2ZXI(2)=D2XXI(1)*OSC+D2XXI(2)*XCC-D2XXI(3)*XSS
     '                               +2.d0*DXXI(1)*DXXI(2)*OCC
     '           -DXXI(2)*DXXI(2)*XSC-2.d0*DXXI(2)*DXXI(3)*XCS
     '           -DXXI(3)*DXXI(3)*XSC-2.d0*DXXI(3)*DXXI(1)*OSS
        D2ZXI(3)=D2XXI(1)*OOS             +D2XXI(3)*XOC
     '           +DXXI(3)*DXXI(3)*XOS+2.d0*DXXI(3)*DXXI(1)*OOC
        DSQXI=DZXI(1)*DZ(1)+DZXI(2)*DZ(2)+DZXI(3)*DZ(3)
        D2SQXI=DZXI(1)*DZXI(1)+D2ZXI(1)*DZ(1)
     '        +DZXI(2)*DZXI(2)+D2ZXI(2)*DZ(2)
     '        +DZXI(3)*DZXI(3)+D2ZXI(3)*DZ(3)
        W=1.d0
        IF(DABS(DSQXI).GT.D2SQXI*VMAX) THEN
          V=-DSIGN(VMAX,DSQXI)
        ELSE
          V=-DSQXI/D2SQXI    
            
        ENDIF
        DSQXIV=DSQXI*V
        IF(DOP) THEN
          WRITE(OP_STRING,'(''    IT1='',I4,10X,''DSQXI(ni)='','
     '      //'(E12.6,4X),/18X,''D2SQXI(MI,ni)='',(E12.6,4X),'
     '      //'/26X,''V(ni)='',(E12.6,4X))') it1,DSQXI,D2SQXI,V
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF  
        DO 4 it2=1,ITMAX
          XIDLIN=XID+V*W
          DO nj=1,NJT
            X(nj)=PLXI(nj,NPL,nr,1,DL,XIDLIN,XP)
          ENDDO
          XOO=X(1)
          OCO=DCOS(X(2))
          OSO=DSIN(X(2))
          OOC=DCOS(X(3))
          OOS=DSIN(X(3))
          OCC=OCO*OOC
          OCS=OCO*OOS
          OSC=OSO*OOC
          OSS=OSO*OOS
          XOC=XOO*OOC
          XOS=XOO*OOS
          XCC=XOO*OCC
          XCS=XOO*OCS
          XSC=XOO*OSC
          XSS=XOO*OSS
          Z(1)=XCC
          Z(2)=XSC         
        
          Z(3)=XOS
          DZ(1)=Z(1)-ZD(1)
          DZ(2)=Z(2)-ZD(2)
          DZ(3)=Z(3)-ZD(3)
          SQLIN=DZ(1)**2+DZ(2)**2+DZ(3)**2
          IF(DOP) THEN
            WRITE(OP_STRING,'(8X,''IT2='',I4,14X,''W='',E12.6,/21X,'
     '        //'''XIDLIN(ni)='',(E12.6,4X),/26X,''X(nj)='','
     '        //'3(E12.6,4X),/26X,''SQLIN='',E12.6)') 
     '        it2,W,XIDLIN,(X(nj),nj=1,NJT),SQLIN
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(SQLIN.LE.SQ) GOTO 5
          W=(DSQXIV*W*W)/(DSQXIV*W+SQ-SQLIN)
    4   CONTINU
      DZ(1)=Z(1)-ZD(1)
      DZ(2)=Z(2)-ZD(2)
      SQ=DZ(1)**2+DZ(2)**2
      IF(DOP) THEN
        WRITE(OP_STRING,'(/'' XD(nj)='',2(E12.6,4X),/''XID(1)='','
     '    //'(E12.6,4X),/''  X(nj)='',2(E12.6,4X),'
     '    //'/''     SQ='',E12.6)')
     '    (ZD(nj),nj=1,2),XID,(X(nj),nj=1,2),SQ  
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      DO 6 it1=1,ITMAX
        DO nj=1,2        
       
           DXXI(nj)=PLXI(nj,NPL,nr,2,DL,XID,XP)
          D2XXI(nj)=PLXI(nj,NPL,nr,3,DL,XID,XP)
        ENDDO
        DZXI(1)=DXXI(1)*SCO-DXXI(2)*CSO
        DZXI(2)=DXXI(1)*CSC+DXXI(2)*SCC-DXXI(3)*SSS
        D2ZXI(1)=D2XXI(1)*SCO-D2XXI(2)*CSO
     '           +DXXI(1)*DXXI(1)*CCO-2.d0*DXXI(1)*DXXI(2)*SSO
     '           -DXXI(2)*DXXI(2)*CCO
        D2ZXI(2)=D2XXI(1)*CSC+D2XXI(2)*SCC-D2XXI(3)*SSS
     '           +DXXI(1)*DXXI(1)*SSC+2.d0*DXXI(1)*DXXI(2)*CCC
     '           -DXXI(2)*DXXI(2)*SSC-2.d0*DXXI(2)*DXXI(3)*SCS
     '           -DXXI(3)*DXXI(3)*SSC-2.d0*DXXI(3)*DXXI(1)*CSS
        DSQXI=DZXI(1)*DZ(1)+DZXI(2)*DZ(2)+DZXI(3)*DZ(3)
        D2SQXI=DZXI(1)*DZXI(1)+D2ZXI(1)*DZ(1)
     '        +DZXI(2)*DZXI(2)+D2ZXI(2)*DZ(2)
        W=1.d0
        IF(DABS(DSQXI).GT.D2SQXI*VMAX) THEN
          V=-DSIGN(VMAX,DSQXI)
        ELSE
          V=-DSQXI/D2SQXI
        ENDIF
        DSQXIV=DSQXI*V
        IF(DOP) THEN
          WRITE(OP_STRING,'(''    IT1='',I4,10X,''DSQXI(ni)='','
     '      //'(E12.6,4X),/18X,''D2SQXI(MI,ni)='',(E12.6,4X),'
     '      //'/26X,''V(ni)='',(E12.6,4X))')   it1,DSQXI,D2SQXI,V
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF   

C ***   Performs line search    
        DO 4 it2=1,ITMAX
          XIDLIN=XID+V*W
          DO nj=1,2
            X(nj)=PLXI(nj,NPL,nr,1,DL,XIDLIN,XP)
          ENDDO
          COO=DCOSH(X(1))*FOCUS
          SOO=DSINH(X(1))*FOCUS
          OCO=DCOS(X(2))
          OSO=DSIN(X(2))
          OOC=DCOS(X(3))
          OOS=DSIN(X(3))
          CCO=COO*OCO
          CSO=COO*OSO
          SCO=SOO*OCO
          SSO=SOO*OSO
          CCC=COO*OCO*OOC
          CSC=COO*OSO*OOC
          CSS=COO*OSO*OOS
          SCC=SOO*OCO*OOC
          SCS=SOO*OCO*OOS
          SSC=SOO*OSO*OOC
          SSS=SOO*OSO*OOS
          Z(1)=CCO
          Z(2)=SSC
          DZ(1)=Z(1)-ZD(1)
          DZ(2)=Z(2)-ZD(2)
          SQLIN=DZ(1)**2+DZ(2)**2
          IF(DOP) THEN
            WRITE(OP_STRING,'(8X,''IT2='',I4,14X,''W='',E12.6,/21X,'
     '        //'''XIDLIN(ni)='',(E12.6,4X),/26X,''X(nj)='','
     '        //'2(E12.6,4X),/26X,''SQLIN='',E12.6)')        
                            
     '        it2,W,XIDLIN,(X(nj),nj=1,2),SQLIN
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF   
          IF(SQLIN.LE.SQ) GO TO 5
          W=(DSQXIV*W*W)/(DSQXIV*W+SQ-SQLIN)
    4   CONTINUE
    5   SQOLD=SQ
        SQ=SQLIN
        DXI=V*W
        XID=XID+DXI
        IF(DOP) THEN
          WRITE(OP_STRING,'(24X,''XID(ni)='',(E12.6,4X))') XID
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF     
        IT=it1
        IF(((DXI**2)/(1.d0+XID**2).LE.TOL2)
     '      .AND.((SQOLD-SQ)/(1.d0+SQ).LE.TOL)) GO TO 7
    6 CONTINUE

C   7 CALL EXITS('CLOS14')
    7 CONTINUE
 9999 RETURN
      END      


      SUBROUTINE CLOS22(IBT,IDO,INP,IT,ITMAX,NPF,SQ,TOL,VMAX,XE,XID,ZD)

C#### Subroutine: CLOS22
C###  Description:
C###    CLOS22 finds the XI-coordinates at the closest approach of a 
C###    face to a data point with coordinates XD using a modified 
C###    Newton algorithm.

C**** ITMAX is the maximum number of iterations.
C**** TOL is the required tolerance of the solution.
C**** VMAX denotes the maximum step length per iteration.
C**** NOTE: Works only for 2-D elements in cylindrical polar
C**** coordinates.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IT,ITMAX,NPF(15)
      REAL*8 SQ,TOL,VMAX,XE(NSM,NJM),XID(NIM),ZD(NJM)
!     Local Variables
      CHARACTER ERROR*10
      INTEGER it1,it2,MI,nb,ni,nj
      REAL*8 D2SQXI(2,2),D2XXI(3,2,2),D2ZXI(3,2,2),DELTA,DET,DSQXI(2),
     '  DSQXIV,DXI(2),DXXI(3,2),DZ(3),DZXI(3,2),OCO,OOX,OSO,PXI,
     '  SQLIN,SQOLD,TOL2,V(2),V2,VMAX2,W,X(3),XCO,XIDLIN(2),XOO,XSO,Z(3)

C     CALL ENTERS('CLOS22',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'(A)') ' >CLOS22 2-D cylindrical polar'
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      TOL2=TOL**2
      VMAX2=VMAX**2
      nb=NPF(10)
      DO nj=1,NJT
        X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),NPF(9+nj),1,XID,
     '    XE(1,nj))
      ENDDO
      XOO=X(1)
      OCO=DCOS(X(2))
      OSO=DSIN(X(2))
      OOX=X(3)
      XCO=XOO*OCO
      XSO=XOO*OSO
      Z(1)=XCO
      Z(2)=XSO
      Z(3)=OOX
      DZ(1)=Z(1)-ZD(1)
      DZ(2)=Z(2)-ZD(2)
      DZ(3)=Z(3)-ZD(3)
      SQ=DZ(1)**2+DZ(2)**2+DZ(3)**2
      IF(DOP) THEN
        WRITE(OP_STRING,'(/'' ZD(nj)='',3(E12.6,4X),/''XID(ni)='','
     '    //'2(E12.6,4X),/''  Z(nj)='',3(E12.6,4X),'
     '    //'/''     SQ='',E12.6)')
     '    (ZD(nj),nj=1,NJT),(XID(ni),ni=1,2),(Z(nj),nj=1,NJT),SQ
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF    
      DO 6 it1=1,ITMAX
        DO nj=1,NJT
          DXXI(nj,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      NPF(9+nj),2,XID,XE(1,nj))
          DXXI(nj,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      NPF(9+nj),4,XID,XE(1,nj))
          D2XXI(nj,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      NPF(9+nj),3,XID,XE(1,nj))
          D2XXI(nj,1,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      NPF(9+nj),6,XID,XE(1,nj))
          D2XXI(nj,2,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      NPF(9+nj),5,XID,XE(1,nj))
        ENDDO
        DZXI(1,1)=DXXI(1,1)*OCO-DXXI(2,1)*XSO
        DZXI(1,2)=DXXI(1,2)*OCO-DXXI(2,2)*XSO
        DZXI(2,1)=DXXI(1,1)*OSO+DXXI(2,1)*XCO
        DZXI(2,2)=DXXI(1,2)*OSO+DXXI(2,2)*XCO
        DZXI(3,1)=DXXI(3,1)
        DZXI(3,2)=DXXI(3,2)
        D2ZXI(1,1,1)=D2XXI(1,1,1)*OCO-D2XXI(2,1,1)*XSO
     '           -DXXI(2,1)*DXXI(2,1)*XCO-2.d0*DXXI(1,1)*DXXI(2,1)*OSO
        D2ZXI(1,1,2)=D2XXI(1,1,2)*OCO-D2XXI(2,1,2)*XSO
     '           -DXXI(2,1)*DXXI(2,2)*XCO-2.d0*DXXI(1,1)*DXXI(2,2)*OSO
        D2ZXI(1,2,2)=D2XXI(1,2,2)*OCO-D2XXI(2,2,2)*XSO
     '           -DXXI(2,2)*DXXI(2,2)*XCO-2.d0*DXXI(1,2)*DXXI(2,2)*OSO
        D2ZXI(2,1,1)=D2XXI(1,1,1)*OSO-D2XXI(2,1,1)*XCO
     '           -DXXI(2,1)*DXXI(2,1)*XSO+2.d0*DXXI(1,1)*DXXI(2,1)*OCO
        D2ZXI(2,1,2)=D2XXI(1,1,2)*OSO-D2XXI(2,1,2)*XCO
     '           -DXXI(2,1)*DXXI(2,2)*XSO+2.d0*DXXI(1,1)*DXXI(2,2)*OCO
        D2ZXI(2,2,2)=D2XXI(1,2,2)*OSO-D2XXI(2,2,2)*XCO
     '           -DXXI(2,2)*DXXI(2,2)*XSO+2.d0*DXXI(1,2)*DXXI(2,2)*OCO
        D2ZXI(3,1,1)=D2XXI(3,1,1)
        D2ZXI(3,1,2)=D2XXI(3,1,2)
        D2ZXI(3,2,2)=D2XXI(3,2,2)
        DSQXI(1)=DZXI(1,1)*DZ(1)+DZXI(2,1)*DZ(2)+DZXI(3,1)*DZ(3)
        DSQXI(2)=DZXI(1,2)*DZ(1)+DZXI(2,2)*DZ(2)+DZXI(3,2)*DZ(3)
        D2SQXI(1,1)=DZXI(1,1)*DZXI(1,1)+D2ZXI(1,1,1)*DZ(1)
     '             +DZXI(2,1)*DZXI(2,1)+D2ZXI(2,1,1)*DZ(2)
     '             +DZXI(3,1)*DZXI(3,1)+D2ZXI(3,1,1)*DZ(3)
        D2SQXI(1,2)=DZXI(1,1)*DZXI(1,2)+D2ZXI(1,1,2)*DZ(1)
     '             +DZXI(2,1)*DZXI(2,2)+D2ZXI(2,1,2)*DZ(2)
     '             +DZXI(3,1)*DZXI(3,2)+D2ZXI(3,1,2)*DZ(3)
        D2SQXI(2,2)=DZXI(1,2)*DZXI(1,2)+D2ZXI(1,2,2)*DZ(1)
     '             +DZXI(2,2)*DZXI(2,2)+D2ZXI(2,2,2)*DZ(2)
     '             +DZXI(3,2)*DZXI(3,2)+D2ZXI(3,2,2)*DZ(3)
        DET=D2SQXI(1,1)*D2SQXI(2,2)-D2SQXI(1,2)**2
        IF(DABS(DET).LE.TOL) THEN
          V(1)=-DSQXI(1)
          V(2)=-DSQXI(2)
        ELSE
          V(1)=-(D2SQXI(2,2)*DSQXI(1)-D2SQXI(1,2)*DSQXI(2))/DET
          V(2)=-(DSQXI(2)+D2SQXI(1,2)*V(1))/D2SQXI(2,2)
          DSQXIV=DSQXI(1)*V(1)+DSQXI(2)*V(2)
          DELTA=DSQRT((DSQXI(1)**2+DSQXI(2)**2)*(V(1)**2+V(2)**2))*TOL
          IF(DABS(DSQXIV).LE.DELTA) THEN
            V(1)=-DSQXI(1)
            V(2)=-DSQXI(2)
          ELSE IF(DSQXIV.GT.DELTA) THEN
            V(1)=-V(1)
            V(2)=-V(2)
          ENDIF
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,'(''    IT1='',I4,10X,''DSQXI(ni)='','
     '      //'2(E12.6,4X),/18X,''D2SQXI(MI,ni)='',2(E12.6,4X),/32X,'
     '      //'2(E12.6,4X),/26X,''V(ni)='',2(E12.6,4X))')  
     '      it1,(DSQXI(ni),ni=1,2),((D2SQXI(MI,ni),MI=1,2),ni=1,2),
     '      (V(ni),ni=1,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF 
        W=1.d0
        V2=V(1)**2+V(2)**2
        IF(V2.GE.VMAX2) W=VMAX/DSQRT(V2)
        DO 4 it2=1,ITMAX
          XIDLIN(1)=XID(1)+V(1)*W
          XIDLIN(2)=XID(2)+V(2)*W
          DO 3 nj=1,NJT
            X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),NPF(9+nj),1,
     '        XIDLIN,XE(1,nj))
 3        CONTINUE
          XOO=X(1)
          OCO=DCOS(X(2))
          OSO=DSIN(X(2))
          OOX=X(3)
          XCO=XOO*OCO
          XSO=XOO*OSO
          Z(1)=XCO
          Z(2)=XSO
          Z(3)=OOX
          DZ(1)=Z(1)-ZD(1)
          DZ(2)=Z(2)-ZD(2)
          DZ(3)=Z(3)-ZD(3)
          SQLIN=DZ(1)**2+DZ(2)**2+DZ(3)**2
          IF(DOP) THEN
            WRITE(OP_STRING,'(8X,''IT2='',I4,14X,''W='',E12.6,/21X,'
     '        //'''XIDLIN(ni)='',2(E12.6,4X),/26X,''Z(nj)='','
     '        //'3(E12.6,4X),/26X,''SQLIN='',E12.6)') 
     '        it2,W,(XIDLIN(ni),ni=1,2),(Z(nj),nj=1,NJT),SQLIN
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF    
          IF(SQLIN.LE.SQ) GOTO 5
          W=(DSQXIV*W*W)/(DSQXIV*W+SQ-SQLIN)
    4   CONTINUE
    5   SQOLD=SQ
        SQ=SQLIN
        DXI(1)=V(1)*W
        DXI(2)=V(2)*W
        XID(1)=XID(1)+DXI(1)
        XID(2)=XID(2)+DXI(2)
        IF(DOP) THEN
          WRITE(OP_STRING,'(24X,''XID(ni)='',2(E12.6,4X))')
     '      (XID(ni),ni=1,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF  
        IT=it1
        IF(((DXI(1)**2+DXI(2)**2)/(1.d0+XID(1)**2+XID(2)**2).LE.TOL2)
     '    .AND.((SQOLD-SQ)/(1.d0+SQ).LE.TOL)) GOTO 7
    6 CONTINUE
C   7 CALL EXITS('CLOS22')
    7 CONTINUE
 9999 RETURN
      END


      SUBROUTINE CLOS23(IBT,IDO,INP,IT,ITMAX,NPF,SQ,TOL,VMAX,XE,XID,ZD)

C#### Subroutine: CLOS23
C###  Description:
C###    CLOS23 finds the XI-coordinates at the closest approach of a 
C###    face to a data point with coordinates XD using a modified 
C###    Newton algorithm.

C**** ITMAX is the maximum number of iterations.
C**** TOL is the required tolerance of the solution.
C**** VMAX denotes the maximum step length per iteration.
C**** NOTE: Works only for 2-D elements in spherical polar
C**** coordinates.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IT,ITMAX,NPF(15)
      REAL*8 SQ,TOL,VMAX,XE(NSM,NJM),XID(NIM),ZD(NJM)
!     Local Variables
      CHARACTER ERROR*10
      INTEGER it1,it2,MI,nb,ni,nj
      REAL*8 D2SQXI(2,2),D2XXI(3,2,2),D2ZXI(3,2,2),DELTA,DET,DSQXI(2),
     '  DSQXIV,DXI(2),DXXI(3,2),DZ(3),DZXI(3,2),OCC,OCO,OCS,OOC,
     '  OOS,OSC,OSO,OSS,PXI,SQLIN,SQOLD,TOL2,V(2),V2,VMAX2,W,X(3),XCC,
     '  XCS,XIDLIN(2),XOC,XOO,XOS,XSC,XSS,Z(3)

C     CALL ENTERS('CLOS23',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'(A)') ' >CLOS23 2-D spherical polar'
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      TOL2=TOL**2
      VMAX2=VMAX**2
      nb=NPF(10)
      DO 1 nj=1,NJT
        X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),NPF(9+nj),1,XID,
     '    XE(1,nj))
 1    CONTINUE
      XOO=X(1)
      OCO=DCOS(X(2))
      OSO=DSIN(X(2))
      OOC=DCOS(X(3))
      OOS=DSIN(X(3))
      OCC=OCO*OOC
      OCS=OCO*OOS
      OSC=OSO*OOC
      OSS=OSO*OOS
      XOC=XOO*OOC
      XOS=XOO*OOS
      XCC=XOO*OCC
      XCS=XOO*OCS
      XSC=XOO*OSC
      XSS=XOO*OSS
      Z(1)=XCC
      Z(2)=XSC
      Z(3)=XOS
      DZ(1)=Z(1)-ZD(1)
      DZ(2)=Z(2)-ZD(2)
      DZ(3)=Z(3)-ZD(3)
      SQ=DZ(1)**2+DZ(2)**2+DZ(3)**2
      IF(DOP) THEN
        WRITE(OP_STRING,'(/'' ZD(nj)='',3(E12.6,4X),/''XID(ni)='','
     '    //'2(E12.6,4X),/''  Z(nj)='',3(E12.6,4X),'
     '    //'/''     SQ='',E12.6)')
     '    (ZD(nj),nj=1,NJT),(XID(ni),ni=1,2),(Z(nj),nj=1,NJT),SQ
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF  
      DO 6 it1=1,ITMAX
        DO 2 nj=1,NJT
          DXXI(nj,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      NPF(9+nj),2,XID,XE(1,nj))
          DXXI(nj,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      NPF(9+nj),4,XID,XE(1,nj))
          D2XXI(nj,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      NPF(9+nj),3,XID,XE(1,nj))
          D2XXI(nj,1,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      NPF(9+nj),6,XID,XE(1,nj))
          D2XXI(nj,2,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      NPF(9+nj),5,XID,XE(1,nj))
    2   CONTINUE
        DZXI(1,1)=DXXI(1,1)*OCC-DXXI(2,1)*XSC-DXXI(3,1)*XCS
        DZXI(1,2)=DXXI(1,2)*OCC-DXXI(2,2)*XSC-DXXI(3,2)*XCS
        DZXI(2,1)=DXXI(1,1)*OSC+DXXI(2,1)*XCC-DXXI(3,1)*XSS
        DZXI(2,2)=DXXI(1,2)*OSC+DXXI(2,2)*XCC-DXXI(3,2)*XSS
        DZXI(3,1)=DXXI(1,1)*OOS              +DXXI(3,1)*XOC
        DZXI(3,2)=DXXI(1,2)*OOS              +DXXI(3,2)*XOC
        D2ZXI(1,1,1)=D2XXI(1,1,1)*OCC-D2XXI(2,1,1)*XSC-D2XXI(3,1,1)*XCS
     '                                   -2.d0*DXXI(1,1)*DXXI(2,1)*OSC
     '           -DXXI(2,1)*DXXI(2,1)*XCC+2.d0*DXXI(2,1)*DXXI(3,1)*XSS
     '           -DXXI(3,1)*DXXI(3,1)*XCC-2.d0*DXXI(3,1)*DXXI(1,1)*OCS
        D2ZXI(1,1,2)=D2XXI(1,1,2)*OCC-D2XXI(2,1,2)*XSC-D2XXI(3,1,2)*XCS
     '                                   -2.d0*DXXI(1,1)*DXXI(2,2)*OSC
     '           -DXXI(2,1)*DXXI(2,2)*XCC+2.d0*DXXI(2,1)*DXXI(3,2)*XSS
     '           -DXXI(3,1)*DXXI(3,2)*XCC-2.d0*DXXI(3,1)*DXXI(1,2)*OCS
        D2ZXI(1,2,2)=D2XXI(1,2,2)*OCC-D2XXI(2,2,2)*XSC-D2XXI(3,2,2)*XCS
     '                                   -2.d0*DXXI(1,2)*DXXI(2,2)*OSC
     '           -DXXI(2,2)*DXXI(2,2)*XCC+2.d0*DXXI(2,2)*DXXI(3,2)*XSS
     '           -DXXI(3,2)*DXXI(3,2)*XCC-2.d0*DXXI(3,2)*DXXI(1,2)*OCS
        D2ZXI(2,1,1)=D2XXI(1,1,1)*OSC+D2XXI(2,1,1)*XCC-D2XXI(3,1,1)*XSS
     '                                   +2.d0*DXXI(1,1)*DXXI(2,1)*OCC
     '           -DXXI(2,1)*DXXI(2,1)*XSC-2.d0*DXXI(2,1)*DXXI(3,1)*XCS
     '           -DXXI(3,1)*DXXI(3,1)*XSC-2.d0*DXXI(3,1)*DXXI(1,1)*OSS
        D2ZXI(2,1,2)=D2XXI(1,1,2)*OSC+D2XXI(2,1,2)*XCC-D2XXI(3,1,2)*XSS
     '                                   +2.d0*DXXI(1,1)*DXXI(2,2)*OCC
     '           -DXXI(2,1)*DXXI(2,2)*XSC-2.d0*DXXI(2,1)*DXXI(3,2)*XCS
     '           -DXXI(3,1)*DXXI(3,2)*XSC-2.d0*DXXI(3,1)*DXXI(1,2)*OSS
        D2ZXI(2,2,2)=D2XXI(1,2,2)*OSC+D2XXI(2,2,2)*XCC-D2XXI(3,2,2)*XSS
     '                                   +2.d0*DXXI(1,2)*DXXI(2,2)*OCC
     '           -DXXI(2,2)*DXXI(2,2)*XSC-2.d0*DXXI(2,2)*DXXI(3,2)*XCS
     '           -DXXI(3,2)*DXXI(3,2)*XSC-2.d0*DXXI(3,2)*DXXI(1,2)*OSS
        D2ZXI(3,1,1)=D2XXI(1,1,1)*OOS                 +D2XXI(3,1,1)*XOC
     '           -DXXI(3,1)*DXXI(3,1)*XOS+2.d0*DXXI(3,1)*DXXI(1,1)*OOC
        D2ZXI(3,1,2)=D2XXI(1,1,2)*OOS                 +D2XXI(3,1,2)*XOC
     '           -DXXI(3,1)*DXXI(3,2)*XOS+2.d0*DXXI(3,1)*DXXI(1,2)*OOC
        D2ZXI(3,2,2)=D2XXI(1,2,2)*OOS                 +D2XXI(3,2,2)*XOC
     '           -DXXI(3,2)*DXXI(3,2)*XOS+2.d0*DXXI(3,2)*DXXI(1,2)*OOC
        DSQXI(1)=DZXI(1,1)*DZ(1)+DZXI(2,1)*DZ(2)+DZXI(3,1)*DZ(3)
        DSQXI(2)=DZXI(1,2)*DZ(1)+DZXI(2,2)*DZ(2)+DZXI(3,2)*DZ(3)
        D2SQXI(1,1)=DZXI(1,1)*DZXI(1,1)+D2ZXI(1,1,1)*DZ(1)
     '             +DZXI(2,1)*DZXI(2,1)+D2ZXI(2,1,1)*DZ(2)
     '             +DZXI(3,1)*DZXI(3,1)+D2ZXI(3,1,1)*DZ(3)
        D2SQXI(1,2)=DZXI(1,1)*DZXI(1,2)+D2ZXI(1,1,2)*DZ(1)
     '             +DZXI(2,1)*DZXI(2,2)+D2ZXI(2,1,2)*DZ(2)
     '             +DZXI(3,1)*DZXI(3,2)+D2ZXI(3,1,2)*DZ(3)
        D2SQXI(2,2)=DZXI(1,2)*DZXI(1,2)+D2ZXI(1,2,2)*DZ(1)
     '             +DZXI(2,2)*DZXI(2,2)+D2ZXI(2,2,2)*DZ(2)
     '             +DZXI(3,2)*DZXI(3,2)+D2ZXI(3,2,2)*DZ(3)
        DET=D2SQXI(1,1)*D2SQXI(2,2)-D2SQXI(1,2)**2
        IF(DABS(DET).LE.TOL) THEN
          V(1)=-DSQXI(1)
          V(2)=-DSQXI(2)
        ELSE
          V(1)=-(D2SQXI(2,2)*DSQXI(1)-D2SQXI(1,2)*DSQXI(2))/DET
          V(2)=-(DSQXI(2)+D2SQXI(1,2)*V(1))/D2SQXI(2,2)
          DSQXIV=DSQXI(1)*V(1)+DSQXI(2)*V(2)
          DELTA=DSQRT((DSQXI(1)**2+DSQXI(2)**2)*(V(1)**2+V(2)**2))*TOL
          IF(DABS(DSQXIV).LE.DELTA) THEN
            V(1)=-DSQXI(1)
            V(2)=-DSQXI(2)
          ELSE IF(DSQXIV.GT.DELTA) THEN
            V(1)=-V(1)
            V(2)=-V(2)
          ENDIF
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,'(''    IT1='',I4,10X,''DSQXI(ni)='','
     '      //'2(E12.6,4X),/18X,''D2SQXI(MI,ni)='',2(E12.6,4X),/32X,'
     '      //'2(E12.6,4X),/26X,''V(ni)='',2(E12.6,4X))') 
     '      it1,(DSQXI(ni),ni=1,2),((D2SQXI(MI,ni),MI=1,2),ni=1,2),
     '      (V(ni),ni=1,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF   
        W=1.d0
        V2=V(1)**2+V(2)**2
        IF(V2.GE.VMAX2) W=VMAX/DSQRT(V2)
        DO 4 it2=1,ITMAX
          XIDLIN(1)=XID(1)+V(1)*W
          XIDLIN(2)=XID(2)+V(2)*W
          DO 3 nj=1,NJT
            X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),NPF(9+nj),1,
     '                XIDLIN,XE(1,nj))
 3        CONTINUE
          XOO=X(1)
          OCO=DCOS(X(2))
          OSO=DSIN(X(2))
          OOC=DCOS(X(3))
          OOS=DSIN(X(3))
          OCC=OCO*OOC
          OCS=OCO*OOS
          OSC=OSO*OOC
          OSS=OSO*OOS
          XOC=XOO*OOC
          XOS=XOO*OOS
          XCC=XOO*OCC
          XCS=XOO*OCS
          XSC=XOO*OSC
          XSS=XOO*OSS
          Z(1)=XCC
          Z(2)=XSC
          Z(3)=XOS
          DZ(1)=Z(1)-ZD(1)
          DZ(2)=Z(2)-ZD(2)
          DZ(3)=Z(3)-ZD(3)
          SQLIN=DZ(1)**2+DZ(2)**2+DZ(3)**2
          IF(DOP) THEN
            WRITE(OP_STRING,'(8X,''IT2='',I4,14X,''W='',E12.6,/21X,'
     '        //'''XIDLIN(ni)='',2(E12.6,4X),/26X,''Z(nj)='','
     '        //'3(E12.6,4X),/26X,''SQLIN='',E12.6)') 
     '        it2,W,(XIDLIN(ni),ni=1,2),(Z(nj),nj=1,NJT),SQLIN
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF  
          IF(SQLIN.LE.SQ) GOTO 5
          W=(DSQXIV*W*W)/(DSQXIV*W+SQ-SQLIN)
    4   CONTINUE
    5   SQOLD=SQ
        SQ=SQLIN
        DXI(1)=V(1)*W
        DXI(2)=V(2)*W
        XID(1)=XID(1)+DXI(1)
        XID(2)=XID(2)+DXI(2)
        IF(DOP) THEN
          WRITE(OP_STRING,'(24X,''XID(ni)='',2(E12.6,4X))')
     '      (XID(ni),ni=1,2) 
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IT=it1
        IF(((DXI(1)**2+DXI(2)**2)/(1.d0+XID(1)**2+XID(2)**2).LE.TOL2)
     '    .AND.((SQOLD-SQ)/(1.d0+SQ).LE.TOL)) GO TO 7
    6 CONTINUE

C   7 CALL EXITS('CLOS23')
    7 CONTINUE
 9999 RETURN
      END

      SUBROUTINE CONT_PLANE_X(NBJ,NJE,NKE,NLS_ICON_ELEM,NLS_ICON_X,
     '  NPNE,NPF,NQE,NLS_CON_PSI,NLS_CON_XI,SE,XA,XE,XP,NCONT,
     '  NORMAL,P,NELEM,XI,SQMAX,ERROR,*) 
                                                                        
C**** AAY 10 Dec 94

C#### Subroutine: CONT_PLANE_X
C###  Description:
C###    CONT_PLANE_X finds the XI-coordinates at the intersection of a 
C###    contour to a plane passing through point P with normal NORMAL. 
C###    Xi coords are put in XI, element number in NELEM if the distance
C###    is > SQMAX the projection is rejected.
                   
      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NCONT,NELEM,NJE(NEM),NKE(NKM,NNM,NBFM,NEFM),
     '  NLS_ICON_ELEM(1000,100),NLS_ICON_X(100),NPNE(NNM,NBFM,NEFM),
     '  NPF(15,NFM),NQE(NSM,NBFM,NEFM)
      REAL*8 NORMAL(*),NLS_CON_PSI(16,2,1000,3,100),
     '  NLS_CON_XI(3,2,1000,100),P(*),
     '  SE(NSM,NBFM,NEFM),SQMAX,XA(NAM,NJM,NQM),XE(NSM,NJM),
     '  XI(3),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER n,nb,ne,ni,nj,ns,nseg
      REAL*8 DIST,DIST_MIN,S,X(3),Z(3,2),ZDIFF,Z_VAL(2)

      CALL ENTERS('CONT_PLANE_X',*9999)

      ne=0 !element number
      DIST_MIN=1.0d5 !for chosing the closest intersection
C     do for all line segments
      DO NSEG=1,NLS_ICON_X(NCONT)
        !find the distances from the line to the image
        IF(NLS_ICON_ELEM(NSEG,NCONT).NE.NE) THEN ! element has changed
          ne=NLS_ICON_ELEM(NSEG,NCONT)
          CALL XPXE(NBJ(1,NE),NJE(NE),NKE(1,1,1,NE),NPNE(1,1,NE),
     '      NPF(1,1),NQE(1,1,NE),SE(1,1,NE),XA,XE,XP,ERROR,*9999)
        ENDIF
        DO n=1,2 !each end of line segment
          DO nj=1,NJT
            nb=NBJ(nj,ne)
            X(nj)=0.0d0
            DO ns=1,NST(nb)
              X(nj)=X(nj)+NLS_CON_PSI(ns,n,NSEG,nj,NCONT)*XE(ns,nj)
            ENDDO !nn
          ENDDO !nj
          CALL XZ(JTYP3,X,Z(1,n))
          Z_VAL(n)=0.0d0
          DO nj=1,NJT
            Z_VAL(n)=Z_VAL(n)+(NORMAL(nj)*(Z(nj,n)-P(nj)))
          ENDDO !nj
        ENDDO !n
        IF(Z_VAL(1)*Z_VAL(2).LE.0.0d0)THEN !intersect
          ZDIFF=Z_VAL(2)-Z_VAL(1)
          IF(ABS(ZDIFF).LT.1.0d-5) THEN
            S=0.5d0
          ELSE
            S=(-Z_VAL(1))/ZDIFF
          ENDIF
          DO nj=1,NJT
            Z(nj,1)=Z(nj,1)+S*(Z(nj,2)-Z(nj,1))
          ENDDO !nj
          DIST=0.0d0
          DO nj=1,NJT
            DIST=DIST+(P(nj)-Z(nj,1))**2
          ENDDO
          IF(DIST.LT.DIST_MIN.AND.DIST.LT.SQMAX)THEN
            NELEM=ne
            DO ni=1,3
              XI(ni)=NLS_CON_XI(NI,1,NSEG,NCONT)+S*
     '          (NLS_CON_XI(ni,2,NSEG,NCONT)-
     '          NLS_CON_XI(ni,1,NSEG,NCONT))
            ENDDO !ni
            DIST_MIN=DIST
          ENDIF
        ENDIF !intersect
      ENDDO !nseg

 9998 CALL EXITS('CONT_PLANE_X')
      RETURN                                                            
 9999 CALL ERRORS('CONT_PLANE_X',ERROR)
      CALL EXITS('CONT_PLANE_X')
      RETURN 1
      END                                                               

      
      SUBROUTINE FITGAU(IBT,IDISP,IDO,INP,ISC_GKK,ISC2_GKK,ISR_GKK,
     '  ISR2_GKK,IWK1,IWK2,LGE,LN,NBH,NBJ,NKE,NKH,NKJ,NONY,NPF,NPNE,
     '  NPNODE,NPNY, NQE,nr,NRE,NVHE,NVJE,NVHP,NVJP,nx,NYNE,NYNO,NYNP,
     '  NYNR,CONY,CYNO,ER,ES,GK,GKK,GR,GRR,PG,SE,WG,WK1,WU,XA,XE,XIG,XO,
     '  XP,YG,FIX,ERROR,*)

C#### Subroutine: FITGAU
C###  Description:
C###    Fits field variables to Gauss point data

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:fit000.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDISP(10),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISC_GKK(NZ_ISC_M),ISC2_GKK(NZ_ISC2_M),
     '  ISR_GKK(NZ_ISR_M),ISR2_GKK(NZ_ISR2_M),IWK1(5*NOM),
     '  IWK2(8*NOM),LGE(NHM*NSM,NRCM),LN(0:*),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NKE(NKM,NNM,NBFM,NEFM),NKH(NHM,NPM,NCM),
     '  NKJ(NJM,NPM),NONY(0:NOYM,NYM,NRCM,0:NRM),NPF(15,NFM),
     '  NPNE(NNM,NBFM,NEFM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM),NQE(NSM,NBFM,NEFM),nr,NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEFM),NVJE(NNM,NBFM,NJM,NEFM),
     '  NVHP(NHM,NPM,NCM),NVJP(NJM,NPM),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM)
      REAL*8 CONY(0:NOYM,NYM,NRCM,0:NRM),CYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),GK(NYROWM,NZ_GK_M),
     '  GKK(NZ_GKK_M),GR(NYROWM),GRR(NOM),PG(NSM,NUM,NGM,NBM),
     '  SE(NSM,NBFM,NEFM),WG(NGM,NBM),WK1(4*NOM),WU(0:10,*),
     '  XA(NAM,NJM,NQM),XE(NSM,NJM),XIG(NIM,NGM,NBM),XO(NOM),
     '  XP(NKM,NVM,NJM,NPM),YG(NGM,NJM,NEM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,5)
!     Local Variables
      INTEGER i,IFAIL,l,last_nj,nb,nbb,ne,ng,NGTOT,nhs1,
     '  nhs2,NHST(2),nj,nj0,njj,nk,no1,no2,no_nynr1,no_nynr2,noy1,
     '  noy2,np,ns1,ns2,nv,ny1,ny2,nyo1,nzz
      REAL*8 co1,co2,EDG,PXI,RESID1,RPMIN,SAED,SMED,SQED,VALUE
      LOGICAL ABORT(4),FIRSTA,GROW,LBLOCK,UPDATE
      DATA ABORT/2*.TRUE.,.FALSE.,.TRUE./,GROW/.TRUE./,LBLOCK/.TRUE./

      CALL ENTERS('FITGAU',*9999)

      FIRSTA=.TRUE. !temporary cpb 11/5/95
 
      nj0=NLH_FIT(1,1,1) !is the first fit variable nj
      last_nj=nj0     !satisfy FTNCHEK
      DO njj=1,NUM_FIT(0) !Loop over the number of fit variables
        nj=NLH_FIT(1,1,njj) !is the nj currently being fitted
        IF(nj.EQ.nj0) THEN !first fit variable is being fitted
          UPDATE=.TRUE.
        ELSE
          UPDATE=.FALSE.
C***      Check if the GK matrix needs to be updated
          DO l=1,LN(0)
            ne=LN(l)
            nb=NBJ(nj,ne)
            nbb=NBJ(last_nj,ne)
            IF(nb.NE.nbb) UPDATE=.TRUE.
          ENDDO !l (ne)
        ENDIF
        last_nj=nj

C*** Calculate solution mapping arrays for the current fit variable

        CALL GLOBALF(nj,NKJ,NONY,NPNODE,nr,NVJP,
     '    nx,NYNO,NYNP,NYNR,CONY,CYNO,FIX,ERROR,*9999)

        IF(NOT(2,nr,nx).EQ.0) THEN
          ERROR=' >>The number of unknowns is zero'
          GOTO 9999
        ENDIF

C*** Initialise variables

        DO no_nynr1=1,NYNR(0,1,1,nr)
          ny1=NYNR(no_nynr1,1,1,nr)
          GR(ny1)=0.d0
          IF(UPDATE) THEN
            DO no_nynr2=1,NYNR(0,2,1,nr)
              ny2=NYNR(no_nynr2,2,1,nr)
              GK(ny1,ny2)=0.d0
            ENDDO !no_nynr2
          ENDIF !UPDATE
        ENDDO !no_nynr1

        DO no1=1,NOT(1,nr,nx)
          GRR(no1)=0.d0
        ENDDO !no1

        IF(UPDATE) THEN
          DO nzz=1,NZ_GKK_M
            GKK(nzz)=0.d0
          ENDDO !nz
        ENDIF

        DO l=1,LN(0) !Loop over elements in the fit
          ne=LN(l)
          nb=NBJ(nj,ne)
          CALL MELGEF(LGE,NBH(1,1,ne),ne,NHST,njj,NKH,NPNE(1,1,ne),
     '      nr,NVHE(1,1,nj,ne),nx,NYNE,NYNP,ERROR,*9999)
          IF(DOP) THEN
            FORMAT='(/'' Element'',I5,'', Number of variables: '','
     '        //'''NHST(1)='',I3,'', NHST(2)='',I3)'
            WRITE(OP_STRING,FORMAT) ne,NHST(1),NHST(2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            FORMAT='('' LGE(1..,1): '',14I5,:(/13X,14I5))'
            WRITE(OP_STRING,FORMAT) (LGE(nhs1,1),nhs1=1,NHST(1))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            FORMAT='('' LGE(1..,2): '',14I5,:(/13X,14I5))'
            WRITE(OP_STRING,FORMAT) (LGE(nhs1,2),nhs1=1,NHST(2))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF            
          DO nhs1=1,NHST(1)
            ER(nhs1)=0.d0
            IF(UPDATE) THEN
              DO nhs2=1,NHST(2)
                ES(nhs1,nhs2)=0.d0
              ENDDO !nhs2
            ENDIF !UPDATE
          ENDDO !nhs1
          CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),NQE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA,XE,XP,ERROR,*9999)
          CALL YGER(nb,ER,PG,SE(1,1,ne),YG(1,NG_FIT(1,njj),ne),
     '      ERROR,*9999)
          IF(UPDATE) THEN
            CALL YGES(nb,ES,PG,SE(1,1,ne),WG,WU(0,ne),ERROR,*9999)
          ENDIF !UPDATE
          IF(DOP) THEN
            WRITE(OP_STRING,'(/'' Element '',I4,'' rhs vector '
     '        //'and stiffness matrix:'')') ne
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO ns1=1,NHST(1)
              WRITE(OP_STRING,'('' ER('',I2,'')    = '',D12.4)') 
     '          ns1,ER(ns1)
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' ES('',I2,'',1..)= '',6D12.4,'
     '          //':(/13X,6D12.4))') ns1,(ES(ns1,ns2),ns2=1,NHST(2))
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF

C*** Assemble element stiffness matrix into global system.

          DO nhs1=1,NHST(1)
            ny1=IABS(LGE(nhs1,1))
            GR(ny1)=GR(ny1)+ER(nhs1)
            IF(UPDATE) THEN
              DO nhs2=1,NHST(2)
                ny2=IABS(LGE(nhs2,2))
                GK(ny1,ny2)=GK(ny1,ny2)+ES(nhs1,nhs2)
              ENDDO !nhs2
            ENDIF !UPDATE
          ENDDO !nhs1

        ENDDO !l (ne)

        IF(UPDATE) THEN
          WRITE(OP_STRING,'('' Assembly of global matrix and RHS '','
     '      //'''complete'')')
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,'('' Assembly of RHS vector complete'')')
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF !UPDATE

        IF(DOP) THEN
          WRITE(OP_STRING,
     '	    '(/'' Global load vector GR & stiffness matrix GK:'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NYNR(0,1,1,nr)='',I5,'
     '      //''', NYNR(0,2,1,nr)='',I5)') NYNR(0,1,1,nr),NYNR(0,2,1,nr)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO no_nynr1=1,NYNR(0,1,1,nr)
            ny1=NYNR(no_nynr1,1,1,nr)
            WRITE(OP_STRING,'('' GR('',I5,'')    = '',D12.4)') 
     '        ny1,GR(ny1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' GK('',I5,'',1..)= '',6D12.4,'
     '        //':(/16X,6D12.4))') ny1,(GK(ny1,NYNR(no_nynr2,2,1,nr)),
     '        no_nynr2=1,NYNR(0,2,1,nr))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO !no_nynr1
        ENDIF
 
C*** Generate the reduced system of equations

        IF(UPDATE) nzz=0
        DO no_nynr1=1,NYNR(0,1,1,nr) !loop of equations (rows) of GK
          ny1=NYNR(no_nynr1,1,1,nr) !equation #
          DO noy1=1,NONY(0,ny1,1,nr) 
            no1=NONY(noy1,ny1,1,nr)
            co1=CONY(noy1,ny1,1,nr)
            GRR(no1)=GR(ny1)*co1
            IF(UPDATE) THEN
              DO no_nynr2=1,NYNR(0,2,1,nr)
                ny2=NYNR(no_nynr2,2,1,nr)
                DO noy2=1,NONY(0,ny2,2,nr)
                  no2=NONY(noy2,ny2,2,nr)
                  co2=CONY(noy2,ny2,2,nr)
                  IF(DABS(GK(ny1,ny2)).GT.1.0d-6) THEN
                    nzz=nzz+1
                    IF(nzz.LE.NZ_GKK_M) GKK(nzz)=GK(ny1,ny2)*co1*co2
                    IF(FIRSTA) THEN
                      IF(nzz.LE.NZ_ISR_M) ISR_GKK(nzz)=no1
                      IF(nzz.LE.NZ_ISC_M) ISC_GKK(nzz)=no2
                    ELSE
                      IF(nzz.LE.NZ_ISR2_M) ISR2_GKK(nzz)=no1
                      IF(nzz.LE.NZ_ISC2_M) ISC2_GKK(nzz)=no2
                    ENDIF
                  ENDIF
                ENDDO !noy2
              ENDDO !no_nynr2
            ENDIF !UPDATE
          ENDDO !noy1
        ENDDO !no_nynr1

        IF(UPDATE) THEN
          NZZT(nr)=nzz
          CALL ASSERT(NZZT(nr).LE.NZ_GKK_M,'>> Increase NZ_GKK_MX',
     '      ERROR,*9999)
          IF(FIRSTA) THEN
            CALL ASSERT(NZZT(nr).LE.NZ_ISC_M,'>> Increase NZ_ISC_MX',
     '        ERROR,*9999)
            CALL ASSERT(NZZT(nr).LE.NZ_ISR_M,'>> Increase NZ_ISR_MX',
     '        ERROR,*9999)
          ELSE
            CALL ASSERT(NZZT(nr).LE.NZ_ISC2_M,'>> Increase NZ_ISC2_MX',
     '        ERROR,*9999)
            CALL ASSERT(NZZT(nr).LE.NZ_ISR2_M,'>> Increase NZ_ISR2_MX',
     '        ERROR,*9999)
          ENDIF
        ENDIF !UPDATE

        IF(DOP) THEN
          WRITE(OP_STRING,'('' NOT(1,nr,nx)='',I5,'
     '      //''', NOT(2,nr,nx)='',I5)') NOT(1,nr,nx),NOT(2,nr,nx)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(/'' Global rhs vector GRR:'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(/'' GRR: '',6D12.4,/:(6X,6D12.4))') 
     '      (GRR(no1),no1=1,NOT(1,nr,nx))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '(/'' Global stiffness matrix GKK (unfactorised):'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NZZT(nr)='',I6)') NZZT(nr)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' ISR: '',6(I5,X),/:(6X,6(I5,X)))') 
     '      (ISR_GKK(nzz),nzz=1,NZZT(nr))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' ISC: '',6(I5,X),/:(6X,6(I5,X)))')
     '      (ISC_GKK(nzz),nzz=1,NZZT(nr))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' GKK: '',6D12.4,/:(6X,6D12.4))') 
     '      (GKK(nzz),nzz=1,NZZT(nr))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        
        IF(UPDATE) THEN
          IF(FIRSTA) THEN
            IFAIL=111
            CALL F01BRF(NOT(1,nr,nx),NZZT(nr),GKK,NZ_ISC_M,ISR_GKK,
     '        NZ_ISR_M,ISC_GKK,0.1d0,IWK1,IWK2,WK1,LBLOCK,GROW,ABORT,
     '        IDISP,IFAIL)
            IF(IFAIL.EQ.0) THEN
              WRITE(OP_STRING,'('' F01BRF has completed LU '','
     '          //'''decomposition:'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              DO I=1,10
                WRITE(OP_STRING,'('' IDISP('',I2,'')='',I10)') I,
     '            IDISP(I)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDDO
              WRITE(OP_STRING,'('' F01BRF: WK1(1)='',D12.4)') WK1(1)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL ASSERT(NZM.GT.IDISP(2),'>>NZM too small',ERROR,*9999)
            ELSE IF(IFAIL.NE.0) THEN
              WRITE(ERROR,'('' >> IFAIL='',I3,'' in Nag routine '','
     '          //'''F01BRF'')') IFAIL
              GOTO 9999
            ENDIF
            FIRSTA=.FALSE.
            
          ELSE IF(.NOT.FIRSTA) THEN
            IFAIL=111
            CALL F01BSF(NOT(1,nr,nx),NZZT(nr),GKK,NZ_ISC_M,ISR2_GKK,
     '        ISC2_GKK,ISC_GKK,IWK1,IWK2,WK1,GROW,1.0d-4,RPMIN,
     '        ABORT(1),IDISP,IFAIL)
            IF(IFAIL.EQ.0) THEN
              WRITE(OP_STRING,'('' F01BSF has completed LU '','
     '          //'''decomposition:''/,'' WK1(1)='',D12.4,'
     '		//''' RPMIN='',D12.4)') WK1(1),RPMIN
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ELSE IF(ifail.NE.0) THEN
              WRITE(ERROR,'('' >> IFAIL='',I3,'' in Nag routine '','
     '          //'''F01BSF'')') IFAIL
              GOTO 9999
            ENDIF
          ENDIF
        ENDIF !UPDATE

        CALL F04AXF(NOT(1,nr,nx),GKK,NZ_ISC_M,ISC_GKK,IWK1,GRR,WK1,1,
     '    IDISP,RESID1)
        IFAIL=0
        IF(IFAIL.EQ.0) THEN
          WRITE(OP_STRING,'('' F04AXF has completed forward and back'','
     '      //''' substitution'')')
      	   CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE IF(IFAIL.NE.0) THEN
          WRITE(OP_STRING,'('' IFAIL='',I3,'' in Nag routine '','
     '      //'''F04AXF'')') IFAIL
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          GO TO 9999
        ENDIF
      
        DO no1=1,NOT(1,nr,nx)
          XO(no1)=GRR(no1)
        ENDDO

        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' Fitted values XO(no):'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          FORMAT='('' XO('',I5,'')= '',D12.4)'
          DO no1=1,NOT(1,nr,nx)
            WRITE(OP_STRING,FORMAT) no1,XO(no1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF

        WRITE(OP_STRING,'(/'' Fitted nodal values XP(nk,nv,nj,np):'')')
      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

        DO no1=1,NOT(2,nr,nx)
          DO nyo1=1,NYNO(0,no1,2,nr)
            ny1=NYNO(nyo1,no1,2,nr)
            co1=CYNO(nyo1,no1,2,nr)
            IF(NPNY(0,ny1,0).EQ.1) THEN
              nk=NPNY(1,ny1,0)
              nv=NPNY(2,ny1,0)
              np=NPNY(4,ny1,0)
            ELSE
              ERROR='>> Element dofs not implemented'
              GOTO 9999
            ENDIF
C MPN 17-Aug-95: don't want to add increment
C old       XP(nk,nv,nj,np)=XP(nk,nv,nj,np)+XO(no1)*co1 !adds increment
            XP(nk,nv,nj,np)=XO(no1)*co1
            WRITE(OP_STRING,'('' XP('',I1,'','',I2,'','',I1,'','',I4,'
     '        //''') = '',D12.6)') nk,nv,nj,np,XP(nk,nv,nj,np)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDDO !nyo1
        ENDDO !no1                                        

C ***   Calculate RMS errors for fit
        SMED=0.d0
        SAED=0.d0
        SQED=0.d0
        NGTOT=0
        DO l=1,LN(0)
          ne=LN(l)
          nb=NBJ(nj,ne)
          CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF,NPNE(1,1,ne),
     '      NQE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),SE(1,1,ne),XA,XE,
     '      XP,ERROR,*9999)
          DO ng=1,NGT(nb)
            VALUE=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '        XIG(1,ng,nb),XE(1,nj))
            EDG=YG(ng,NG_FIT(1,njj),ne)-VALUE
            SMED=SMED+EDG
            SAED=SAED+DABS(EDG)
            SQED=SQED+EDG**2
            NGTOT=NGTOT+1
          ENDDO !ng
        ENDDO !l (ne)
        IF(NGTOT.GT.1) THEN
          WRITE(OP_STRING,
     '     '(/'' Average error           : '',D12.6,'' +/- '',D12.6,'
     '    //'/'' Average absolute error  : '',D12.6,'' +/- '',D12.6,'
     '    //'/'' Root mean squared error : '',D12.6)') 
     '      SMED/NGTOT,
     '      DSQRT((SQED-SMED**2/NGTOT)/(NGTOT-1)),
     '      SAED/NGTOT,DSQRT((SQED-SAED**2/NGTOT)/(NGTOT-1)),
     '      DSQRT(SQED/NGTOT)
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ENDDO !njj
  
      CALL EXITS('FITGAU')
      RETURN
 9999 CALL ERRORS('FITGAU',ERROR)
      CALL EXITS('FITGAU')
      RETURN 1
      END


      SUBROUTINE FITGEO(IBT,IDO,INP,IWORK,LD,LDR,LGE,
     '  LN,NXI,NBJ,NDDL,NDLT,NFF,NJE,NJP,NKE,NKJ,NLL,NNL,NONY,
     '  NPE,NPF,NPL,NPO,NQE,nr,nx,
     '  CONY,DL,EDD,ER,ES,SCALE,SE,SQ,VE,
     '  WD,WDL,XA,XE,XID,XIDL,XP,YP,ZD,ZDL,
     '  GKK,GRR,WK1,WK2,FIX,ERROR,*)

C**** Fits geometric variables to data points (nonlinear).
C *** Note that RWORK has been deleted. Check that its replacement WK1
C *** is dimensioned large enough.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b21.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:fit000.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
!     Parameter List
      INTEGER IBT(2,NIM,*),IDO(NKM,0:NIM,*),INP(NNM,NIM,*),IWORK(*),
     '  LD(*),LDR(*),LGE(*),LN(0:*),
     '  NBJ(NJM,*),NDDL(NEM,*),NDLT(*),NFF(6,*),NJE(*),NJP(*),
     '  NKE(NKM,NNM,NBM,*),NKJ(NJM,*),NLL(12,*),NNL(4,12,*),
     '  NONY(0:NOYM,*),NPE(NNM,NBM,*),NPF(12,*),NPL(20,*),NPO(0:*),
     '  NQE(NSM,NBM,*),nr,nx,NXI(-NIM:NIM,0:*)
      REAL*8 CONY(0:NOYM,*),DL(3,*),EDD(*),ER(*),
     '  ES(NVM,*),SCALE(*),SE(NSM,NBM,*),SQ(*),VE(NSM,NKM,*),
     '  WD(NJM,*),WDL(NJM,*),XA(NAM,NJM,*),XE(NSM,*),XID(NIM,*),
     '  XIDL(NIM,*),XP(NKM,NJM,*),YP(NYM,*),ZD(NJM,*),
     '  ZDL(NJM,*)
      REAL*8 GKK(*),GRR(*),WK1(*),WK2(*)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,*)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,INFO
      REAL*8 ACC,X(4),XI(3),Z(4)
      CHARACTER CHAR1*100,CHAR2*100,CHAR3*100,CFROMI*100,CFROML*100,
     '  CFROMR*100,CIOT*1,CJET*1,CJOT*1,CNJ*4,CNK*4,CNO*4,CNP*4,CNY*4,
     '  CSC*12,CXO*12,FNAME*8,LFTYPE*8
      LOGICAL EXIST,FIRST,LFNO,LFXI,NODE,SETXID

C     COMMON /VA13BD/ SETXID,IPRINT,LP,MAXFUN,MODE,NFUN   !put in an include
C      EXTERNAL FUNC                                      !file if needed again

      CALL ENTERS('FITGEO',*9999)
        FORMAT='($,'' Specify the maximum number of'//
     '    ' function calls  [unlimited]: '',I4)'
        CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,IMAX,
     '    INFO,ERROR,*9999)
C        MAXFUN=IDATA(1) !appears useless - uncomment if needed MPN 20/1/92
        RDEFLT(1)=1.0E-4
        CHAR1=CFROMR(RDEFLT(1),'(E12.6)')
        CALL TRIM(CHAR1,IBEG1,IEND1)
        FORMAT='($,'' Specify the relative accuracy for termination '//
     '    '['//CHAR1(IBEG1:IEND1)//' ]: '',E12.6)'
        CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
        ACC=RDATA(1)

        FORMAT='(/'' Optimisation begins ...'')'
        WRITE(OP_STRING,FORMAT)
      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C       CALL VA13AD(FUNC,IBT,IDO,INP,IO5,LD,LDR,LN,NXI,NBJ,NDDL,NDLT,
C    '    NFF,NJE,NJP,NKE,NKJ,NONY,NPE,NPF,NPL,NPO,NQE,
C    '    CONY,DL,WK1,SE,SQ,SS,VE,WD,WDL,XA,XE,XID,
C    '    XIDL,XO,XP,ZD,ZDL,SCALE,ACC,FIX,ERROR,*9999)
        FORMAT='(/$,'' Do you want to save the current data values'//
     '    ' on a disk file [N]? '',A)'
        CALL AINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,ADATA,ANO,
     '    INFO,ERROR,*9999)
        IF(ADATA(1).EQ.'Y') THEN
          CALL TRIM(FOPTM,IBEG,IEND)
          FORMAT='(/$,'' Specify the name of the data file ['//
     '      FOPTM(IBEG:IEND)//'.dat]: '',A)'
          CALL CINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,CDATA,FOPTM,8,
     '      INFO,ERROR,*9999)
          FNAME=CDATA(1)
C          CALL OPENF(9,'DISK',FNAME//' IODATA '//IOMODE,'NEW',
C     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
C          CALL IODATA('WRITE','GEOMETRY',9,NJT+JTYP9+JTYP11,WD,ZD,
C     '      ERROR,*9999)
        ENDIF

      CALL EXITS('FITGEO')
      RETURN
 9999 CALL TRIM(ERROR,IBEG,IEND)
      ERROR=ERROR(IBEG:IEND)//' >FITGEO'
      RETURN 1
      END


      SUBROUTINE FITFOU(IBT,ICN,ICN2,IDO,INP,IRN,IRN2,IWK1,IWK2,LN,
     '  NBJ,NDDL,NDLT,NJE,NKE,NKJ,NONY,NPF,NPNE,NPO,NQE,nr,NRE,
     '  NVNE,NVNP,nx,NYNO,CONY,ER,ES,GK,PG,SE,WD,WDL,WG,WU,XA,XE,
     '  XID,XIDL,
     '  XP,YP,ZD,ZDL,GKK,GRR,WK1,XO,ERROR,*)

C**** Fits Fourier motion coefficients to data points (KTYP8=5).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b21.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:fit000.cmn'
      INCLUDE 'cmiss$reference:four00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:moti00.cmn'
!     Parameter List
      INTEGER IBT(2,NIM,*),ICN(*),ICN2(*),IDO(NKM,0:NIM,*),
     '  INP(NNM,NIM,*),IRN(*),IRN2(*),IWK1(*),IWK2(*),LN(0:*),
     '  NBJ(NJM,*),NDDL(NEM,*),NDLT(*),NJE(*),NKE(NKM,NNM,NBFM,*),
     '  NKJ(NJM,*),NONY(0:NOYM,NYM,*),NPF(15,*),NPNE(NNM,NBFM,*),
     '  NPO(0:*),
     '  NQE(NSM,NBFM,*),nr,NRE(*),NVNE(NNM,NBFM,NJM,*),NVNP(NJM,NPM,*),
     '  nx,NYNO(0:NYOM,NOOPM,*)
C!!!    note gk dimensioned to NOM here due to space problems AAY 3/12/90
      REAL*8 CONY(0:NOYM,NYM,*),ER(*),ES(NHM*NSM,*),GK(NOM,*),
     '  PG(NSM,NUM,NGM,*),SE(NSM,NBFM,*),WD(NJM,*),WDL(NJM,*),
     '  WG(NGM,*),WU(0:10,*),XA(NAM,NJM,*),XE(NSM,*),XID(NIM,*),
     '  XIDL(NIM,*),XP(NKM,NVM,NJM,*),YP(NYM,*),ZD(NJM,*),ZDL(NJM,*)
      REAL*8 GKK(*),GRR(*),WK1(*),XO(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,IBEG,IDISP(10),IEND,IFAIL,l,n,n1k,n1n,n1o,n1oy,N1P,
     '  n1s,N1Y,n2k,n2n,n2o,n2oy,N2P,n2s,N2Y,nb,ne,nj,NJO,nk,no,noy,np,
     '  nrc,ny,NYPJK,nz
      REAL*8 C1O,C2O,CO,RESID1,RPMIN
      LOGICAL ABORT(4),FIRSTS,GROW,LBLOCK
      DATA ABORT/2*.TRUE.,.FALSE.,.TRUE./,GROW/.TRUE./,LBLOCK/.TRUE./

      CALL ENTERS('FITFOU',*9999)
      CALL ASSERT(OMEGA.GT.0.0D0,'>>no Fourier basis function defined',
     '  ERROR,*9999)

C      NJO=NJO0
      NJO=NJ_FIT(1,1)
      !calculate optimisation dofs no=NONY, total NOT(nr,nx) & coeffs CONY
      CALL NY_NO(NJO,NKJ,NONY,NPO,nr,nx,NYNO,CONY,ERROR,*9999)

      DO n1o=1,NOT(1,nr,nx)
        GRR(n1o)=0.D0
        IF(CALCGSM)THEN
          DO n2o=1,NOT(1,nr,nx)
            GK(n1o,n2o)=0.0D0
          ENDDO
        ENDIF
      ENDDO
      IF(CALCGSM)THEN
        DO nz=1,NZM
          GKK(nz)=0.0D0
        ENDDO
      ENDIF

      DO l=1,LN(0)
        ne=LN(l)
        nb=NBJ(NJO,ne) !Fourier basis
        CALL XPXE(NBJ(1,ne),NJE(ne),NKE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),NQE(1,1,ne),NRE(ne),NVNE(1,1,1,ne),NVNP,
     '    SE(1,1,ne),XA,XE,XP,ERROR,*9999)
        CALL ZDERF(IBT,IDO,INP,NBJ,NDDL,NDLT(ne),ne,NJO,ER,SE(1,1,ne),
     '    WD,WDL,XE,XID,XIDL,ZD,ZDL,ERROR,*9999)
        IF(CALCGSM)THEN
          CALL ZDESF(IBT(1,1,nb),IDO(1,0,nb),INP(1,1,nb),nb,
     '      NDLT(ne),ne,NJO,ES,PG,SE(1,1,ne),WDL,WG,WU(1,ne),
     '      XIDL,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,'(/'' Element '',I2,'' rhs vector '
     '        //'and stiffness matrix:'')') ne
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
	    DO n1s=1,NST(nb)
              WRITE(OP_STRING,'('' ER= '',E12.6,'' ES= '','
     '          //'(T23,8(E12.6,1X)))')
     '          ER(n1s),(ES(n1s,n2s),n2s=1,NST(nb))
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF
        ENDIF
        n1s=0
        DO n1n=1,NNT(NBO)
          N1P=NPNE(n1n,NBO,ne)
          DO n1k=1,NKT(n1n,nb)
            N1Y=NYPJK(NJO,n1k,NKJ,N1P,NPO)
            n1s=n1s+1
            IF(DOP) THEN
              WRITE(OP_STRING,*) ' N1N=',n1n,' N1P=',N1P,' N1K=',n1k,
     '          ' N1Y=',N1Y
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            DO n1oy=1,NONY(0,N1Y,1)
              n1o=NONY(n1oy,N1Y,1)
              C1O=CONY(n1oy,N1Y,1)
              IF(n1o.GT.0) THEN
                GRR(n1o)=GRR(n1o)+ER(n1s)*C1O
                IF(CALCGSM)THEN
                  n2s=0
                  DO n2n=1,NNT(NBO)
                    N2P=NPNE(n2n,NBO,ne)
                    DO n2k=1,NKT(n2n,nb)
                      N2Y=NYPJK(NJO,n2k,NKJ,N2P,NPO)
                      n2s=n2s+1
                      DO n2oy=1,NONY(0,N2Y,2)
                        n2o=NONY(n2oy,N2Y,2)
                        C2O=CONY(n2oy,N2Y,2)
                        IF(n2o.GT.0) THEN
                          GK(n2o,n1o)=GK(n2o,n1o)+ES(n2s,n1s)*C2O*C1O
                        ENDIF
                      ENDDO
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      IF(CALCGSM)THEN
        nz=0
        DO n1o=1,NOT(1,nr,nx)
          DO n2o=1,NOT(2,nr,nx)
            IF(DABS(GK(n1o,n2o)).GT.1.0D-6) THEN
              nz=nz+1
              GKK(nz)=GK(n1o,n2o)
              IRN(nz)=n1o
              ICN(nz)=n2o
            ENDIF
          ENDDO
        ENDDO
        NZZT(nr)=nz
      ENDIF

      WRITE(OP_STRING,'('' Assembly of global matrix complete'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(/'' NOT(1:2,nr,nx)='',2I4,'' NZZT(nr)='',I6)')
     '  (NOT(nrc,nr,nx),nrc=1,2),NZZT(nr)
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C     IF(DOP) THEN
C       WRITE(OP_STRING,'(/A)')
C    '    ' Global rhs vector GRR and stiffness matrix GK:'
C      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C       DO n1o=1,NOT(1,nr,nx)
C         WRITE(OP_STRING,'('' GRR= '',D12.4,'' GK= '',
C    '      (T23,8(E12.4,1X)))') GRR(n1o),(GK(n1o,n2o),n2o=1,NOT(1,nr,nx))
C      	  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C       ENDDO
C     ENDIF

      IF(DOP) THEN
C       FORMAT='('' Global rhs vector GRR and stiffness matrix GKK:'')'
C       WRITE(OP_STRING,FORMAT)
C      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C       FORMAT='('' GRR:'',8D12.4/,(4X,8D12.4))'
C       WRITE(OP_STRING,FORMAT) (GRR(no),no=1,NOT(1,nr,nx))
C      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C       FORMAT='('' IRN: '',20I5/,(6X,20I5))'
C       WRITE(OP_STRING,FORMAT) (IRN(nz),nz=1,NZZT(nr))
C      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C       FORMAT='('' ICN: '',20I5/,(6X,20I5))'
C       WRITE(OP_STRING,FORMAT) (ICN(nz),nz=1,NZZT(nr))
C      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C       FORMAT='('' GKK:'',8D12.4/,(4X,8D12.4))'
C       WRITE(OP_STRING,FORMAT) (GKK(nz),nz=1,NZZT(nr))
C      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(CALCGSM)THEN
        FIRSTS=.TRUE.
        IF(FIRSTS) THEN
          IFAIL=111
          CALL F01BRF(NOT(1,nr,nx),NZZT(nr),GKK,NZM,IRN,NZM,ICN,0.5d0,
     '      IWK1,IWK2,WK1,LBLOCK,GROW,ABORT,IDISP,IFAIL)
          IF(IFAIL.EQ.0) THEN
            WRITE(OP_STRING,'('' F01BRF has completed LU '','
     '        //'''decomposition:''/,'' IDISP(1)='',I5,'
     '        //''' IDISP(2)='',I5, '' IDISP(3)='',I5,'
     '        //''' IDISP(4)='',I5,'' IDISP(5)='',I5,/'
     '        //''' IDISP(6)='',I5,'' IDISP(7)='',I5, '
     '        //''' IDISP(8)='',I5, '' IDISP(9)='',I5,'
     '        //''' IDISP(10)='',I5)') (IDISP(I),I=1,10)
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' F01BRF: WK1(1)='',D12.4)') WK1(1)
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL ASSERT(NZM.GT.IDISP(2),'>>NZM too small',ERROR,*9999)
          ELSE IF(ifail.NE.0) THEN
            WRITE(OP_STRING,'('' IFAIL='',I3,'' in Nag routine '','
     '        //'''F01BRF'')') IFAIL
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            GO TO 9999
          ENDIF
        ELSE IF(.NOT.FIRSTS) THEN
          IFAIL=111
          CALL F01BSF(NOT(1,nr,nx),NZZT(nr),GKK,NZM,IRN2,ICN2,ICN,
     '      IWK1,IWK2,WK1,GROW,1.0d-4,RPMIN,ABORT(1),IDISP,IFAIL)
          IF(IFAIL.EQ.0) THEN
            WRITE(OP_STRING,'('' F01BSF has completed LU '','
     '        //'''decomposition:''/,'' WK1(1)='',D12.4,'' RPMIN='','
     '        //'D12.4)') WK1(1),RPMIN
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(ifail.NE.0) THEN
            WRITE(OP_STRING,'('' IFAIL='',I3,'' in Nag routine '','
     '        //'''F01BSF'')') IFAIL
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            GO TO 9999
          ENDIF
        ENDIF
        FIRSTS=.FALSE.
      ENDIF
      CALL F04AXF(NOT(1,nr,nx),GKK,NZM,ICN,IWK1,GRR,WK1,1,IDISP,RESID1)
      IF(IFAIL.EQ.0) THEN
        WRITE(OP_STRING,'('' F04AXF has completed forward and back'','
     '    //''' substitution'')')
      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ELSE IF(ifail.NE.0) THEN
        WRITE(OP_STRING,'('' IFAIL='',I3,'' in Nag routine FAXF'')')
     '    IFAIL
       	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        GO TO 9999
      ENDIF

      DO no=1,NOT(1,nr,nx)
        XO(no)=GRR(no)
      ENDDO
      WRITE(OP_STRING,'(/'' Fitted values XO(no):'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      FORMAT='('' XO('',I4,'')= '',D12.4)'
      DO no=1,NOT(1,nr,nx)
        WRITE(OP_STRING,FORMAT) no,XO(no)
      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDDO

      FORMAT='(/'' Fitted Fourier coefficients YP(ny,1):'')'
      WRITE(OP_STRING,FORMAT)
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      FORMAT='('' YP('',I4,'',1) = '',(E12.4)))'
      ny=0
      DO n=1,NPO(0)
        np=NPO(n)
        DO nj=1,NJT
          ny=(np-1)*NJT*NKT(0,nb)+(nj-1)*NKT(0,nb)
          DO nk=1,NKT(0,nb)
            ny=ny+1 !these NYs for YP (nodal dof)
            IF(nj.EQ.NJG)THEN !update nodal values
              N2Y=NYPJK(NJO,nk,NKJ,np,NPO) !optimisation dof
              DO noy=1,NONY(0,N2Y,2)
                no=NONY(noy,N2Y,2)
                CO=CONY(noy,N2Y,2)
                IF(no.GT.0) THEN
                  YP(ny,MOTION_IY)=XO(no)*CO
                ENDIF
              ENDDO
              WRITE(OP_STRING,FORMAT) ny,YP(ny,MOTION_IY)
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      CALL EXITS('FITFOU')
      RETURN
 9999 CALL TRIM(ERROR,IBEG,IEND)
      ERROR=ERROR(IBEG:IEND)//' >FITFOU'
      RETURN 1
      END


      SUBROUTINE IOHESS(COMAND,IUNIT,NOT,XO,WK1,ERROR,*)

C**** Handles I/O of Hessian matrix for nonlinear optimisation.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
!     Parameter List
      INTEGER IUNIT,NOT
      REAL*8 WK1(*),XO(*)
      CHARACTER COMAND*(*),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,NOT1
      CHARACTER ACCESS*10,CFROMI*10,FORM*11,IFMT*16,RFMT*16
      LOGICAL OPENED

      CALL ENTERS('IOHESS',*9999)
C     IF(DOP) THEN
C       WRITE(OP_STRING,FMT='(A)') ' >IOHESS'
C       CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C     ENDIF
C     INQUIRE(UNIT=IUNIT,OPENED=OPENED)
C     IF(OPENED) THEN
C       INQUIRE(UNIT=IUNIT,FORM=FORM)
C       IF(FORM.EQ.'FORMATTED') THEN
C         INQUIRE(UNIT=IUNIT,ACCESS=ACCESS)
C         IF(ACCESS.EQ.'SEQUENTIAL') THEN
C           REWIND(UNIT=IUNIT)
C           NPEREC=5
C           NWORKT=NOT*(NOT+1)/2
C           IFMT='(I10)'
C           RFMT='('//
C    '        CTRIM(CFROMI(NPEREC,'(I10)'),IBEG,IEND)(IBEG:IEND)//
C    '        'E24.16)'
C           IF(COMAND.EQ.'READ') THEN
C             NOT1=NOT
C             READ(UNIT=IUNIT,FMT=IFMT,IOSTAT=IOSTAT) NOT
C             IF((IOSTAT.EQ.0).AND.(NOT.EQ.NOT1)) THEN
C               READ(UNIT=IUNIT,FMT=RFMT,IOSTAT=IOSTAT)
C    '            (XO(NO),NO=1,NOT)
C               READ(UNIT=IUNIT,FMT=RFMT,IOSTAT=IOSTAT)
C    '            (WK1(NWORK),NWORK=1,NWORKT)
C               IF(IOSTAT.NE.0) THEN
C                 ERROR=' Read error'
C                 GOTO 9999
C               ENDIF
C             ENDIF
C           ELSE IF(COMAND.EQ.'WRITE') THEN
C             WRITE(UNIT=IUNIT,FMT=IFMT) NOT
C             WRITE(UNIT=IUNIT,FMT=RFMT) (XO(NO),NO=1,NOT)
C             WRITE(UNIT=IUNIT,FMT=RFMT) (WK1(NWORK),NWORK=1,NWORKT)
C           ELSE
C             ERROR=' Command error: COMAND='//COMAND
C             GOTO 9999
C           ENDIF
C         ELSE
C             ERROR=' File '//
C    '          CTRIM(CFROMI(IUNIT,'(I10)'),IBEG,IEND)(IBEG:IEND)//
C    '          ' is not connected for sequential access'
C             GOTO 9999
C         ENDIF
C       ELSE
C         ERROR=' File '//
C    '      CTRIM(CFROMI(IUNIT,'(I10)'),IBEG,IEND)(IBEG:IEND)//
C    '      ' is not connected for formatted i/o'
C         GOTO 9999
C       ENDIF
C     ELSE
C       ERROR=' File '//
C    '    CTRIM(CFROMI(IUNIT,'(I10)'),IBEG,IEND)(IBEG:IEND)//
C    '    ' is not opened'
C         GOTO 9999
C     ENDIF

      CALL EXITS('IOHESS')
      RETURN
 9999 CALL TRIM(ERROR,IBEG,IEND)
      ERROR=ERROR(1:IEND)//' >IOHESS'
      RETURN 1
      END


      SUBROUTINE IPDATA(IBT,IDO,INP,IWORK,LD,LDR,LGE,
     '  LN,NXI,NBJ,NCNE,NCNP,NDDL,NDLT,NEELEM,
     '  NFF,NJE,NJP,NKE,NKJ,NLL,NNL,NONY,
     '  NPE,NPF,NPL,NPNODE,NPO,NQE,
     '  CONY,DL,EDD,ER,ES,GKK,GRR,SCALE,SE,SQ,VE,
     '  WD,WDL,XA,XE,XID,XIDL,XO,XP,YP,ZD,ZDL,FIX,ERROR,*)

C**** Reads data points from keyboard.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b21.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:fit000.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
!     Parameter List
      INTEGER IBT(2,NIM,*),IDO(NKM,0:NIM,*),INP(NNM,NIM,*),IWORK(*),
     '  LD(*),LDR(*),LGE(*),LN(0:*),
     '  NBJ(NJM,*),
     '  NCNE(NNM,NBM,NEM,*),NCNP(0:NCM,NJM,NPM,*),
     '  NDDL(NEM,*),NDLT(*),NEELEM(0:NEM,0:*),
     '  NFF(6,*),NJE(*),NJP(*),
     '  NKE(NKM,NNM,NBM,*),NKJ(NJM,*),NLL(12,*),NNL(4,12,*),
     '  NONY(0:NOYM,*),NPE(NNM,NBM,*),NPF(12,*),NPL(20,*),
     '  NPNODE(0:NPM,0:*),
     '  NPO(0:*),NQE(NSM,NBM,*),NXI(-NIM:NIM,0:*)
      REAL*8 CONY(0:NOYM,*),DL(3,*),EDD(*),ER(*),
     '  ES(NVM,*),SCALE(*),SE(NSM,NBM,*),SQ(*),VE(NSM,NKM,*),WD(NJM,*),
     '  WDL(NJM,*),XA(NAM,NJM,*),XE(NSM,*),XID(NIM,*),XIDL(NIM,*),
     '  XP(NKM,NJM,*),YP(NYM,*),ZD(NJM,*),ZDL(NJM,*)
      REAL*8 GKK(*),GRR(*),XO(*)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,*)
!     Local Variables
      INTEGER IBEG,IBEG0,IBEG1,ICOORD,IEND,IEND0,IEND1,INFO,LASTL,LFMAX,
     '  NB,ND,NE,NI,NJ,NJJ
      REAL*8 X(4),XI(3),Z(4)
      CHARACTER CFROMI*100,CFROML*100,CFROMR*100,CHAR1*100,CHAR2*100,
     '  CHAR3*100,CIOT*1,CJET*1,CJOT*1,CNJ*4,CNK*4,CNO*4,CNP*4,CNY*4,
     '  CSC*12,CXO*12,FNAME*8,LFTYPE*8
      LOGICAL EXIST,FIRST,LFNO,LFXI,NODE,SETXID

C     COMMON /VA13BD/ SETXID,IPRINT,LP,MAXFUN,MODE,NFUN  !put in an include
C      EXTERNAL FUNC                                     !file if needed again

      CALL ENTERS('IPDATA',*9999)

        CALL TRIM(LFTYPE,IBEG0,IEND0)
        FORMAT='(/'' Do you wish to specify the '//LFTYPE(IBEG0:IEND0)
     '    //' associated with each data point [N]? '',A)'
        CALL AINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,ADATA,ANO,
     '    INFO,ERROR,*9999)
        IF(ADATA(1).EQ.'Y') THEN
          LFNO=.TRUE.
          FORMAT='(/'' Do you wish to specify the '//
     '      LFTYPE(IBEG0:IEND0)//
     '      ' coordinate for each data point [N]? '',A)'
          CALL AINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,ADATA,ANO,
     '      INFO,ERROR,*9999)
          IF(ADATA(1).EQ.'Y') THEN
            LFXI=.TRUE.
            DO NI=1,NIOT
              RDEFLT(NI)=0.50D0
            ENDDO
          ELSE
            LFXI=.FALSE.
          ENDIF
        ELSE
          LFNO=.FALSE.
          LFXI=.FALSE.
        ENDIF
        IF((KTYP1.EQ.1).OR.(.NOT.LFXI)) THEN
          FORMAT='(/'' Specify the coordinate system in which the'//
     '      ' data are defined (Default is 1):'''//
     '      '/''   (1) Rectangular Cartesian'''//
     '      '/''   (2) Cylindrical polar'''//
     '      '/''   (3) Spherical polar'''//
     '      '/''   (4) Prolate spheroidal'''//
     '      '/''   (5) Oblate spheroidal'''//
     '      '/''    '',I1)'
          CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IONE,1,5,
     '      INFO,ERROR,*9999)
          ICOORD=IDATA(1)
        ENDIF
        FORMAT='(/'' How many data points do you want to specify'//
     '    ' [0]? '',I4)'
        CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IZERO,0,NDM,
     '    INFO,ERROR,*9999)
        NDT=IDATA(1)
        LASTL=LD(1)
        DO 3070 ND=1,NDT
          CHAR1=CFROMI(ND,'(I10)')
          CALL TRIM(CHAR1,IBEG1,IEND1)
          FORMAT='(/'' For data point '//CHAR1(IBEG1:IEND1)//':'')'
          CALL INOUT(IOTYPE,IVDU,IFILE,FORMAT,ERROR,*9999)
          IF(KTYP1.EQ.2) THEN
            FORMAT='('' Enter the '//CJOT//' field coordinates'//
     '        ' [0..]: ''/1X,'//CJOT//'(E12.6,4X))'
            CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,NJOT,RDATA,RZERO,
     '        -RMAX,RMAX,INFO,ERROR,*9999)
C cpb 7/4/94 adding nj_loc
C           DO NJ=NJT+1,NJT+NJTOT
            DO NJJ=1,NJ_LOC(NJL_FIEL,0)
              NJ=NJ_LOC(NJL_FIEL,NJJ)
C              ZD(NJ,ND)=RDATA(NJ-NJT)
              ZD(NJ,ND)=RDATA(NJ_LOC(NJL_GEOM,NJJ))
              WD(NJ,ND)=1.0D0
           ENDDO
          ENDIF
          IF((KTYP1.EQ.1).OR.(.NOT.LFXI)) THEN
            FORMAT='('' Enter the '//CJET//' geometric coordinates'//
     '        ' [0..]: ''/1X,'//CJET//'(E12.6,4X))'
            CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,NJT,RDATA,RZERO,
     '        -RMAX,RMAX,INFO,ERROR,*9999)
            CALL XZ(ICOORD,RDATA,ZD(1,ND))
            DO NJ=1,NJT
              WD(NJ,ND)=1.0D0
            ENDDO
          ENDIF
          IF(LFNO) THEN
            CHAR1=CFROMI(LASTL,'(I10)')
            CALL TRIM(CHAR1,IBEG1,IEND1)
            FORMAT='('' Enter the initial '//LFTYPE(IBEG0:IEND0)//
     '        ' number ['//CHAR1(IBEG1:IEND1)//']: '',I4)'
            CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,LASTL,0,NLM,
     '        INFO,ERROR,*9999)
            LD(ND)=IDATA(1)
            LASTL=IDATA(1)
          ELSE
            LD(ND)=1
          ENDIF
          IF(LFXI) THEN
            FORMAT='('' Enter the '//CIOT//' '//
     '        LFTYPE(IBEG0:IEND0)//' coordinates [midpoint]: ''/1X,'//
     '        CIOT//'(E12.6,4X))'
            CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,NIOT,RDATA,RDEFLT,
     '        -RMAX,RMAX,INFO,ERROR,*9999)
            DO 3050 NI=1,NIOT
              XID(NI,ND)=RDATA(NI)
 3050       CONTINUE
          ELSE
            CALL XCOORD(IBT,IDO,INP,NBJ,NCNE,NCNP,NE,NEELEM,
     '        NJE,NJP,NKE,NPE,NPF,NPNODE,NQE,
     '        SE,XA,XE,XI,XP,ZD(1,ND),ERROR,*9999)
            LD(ND)=NE
            NDLT(NE)=NDLT(NE)+1
            NDDL(NE,NDLT(NE))=ND
            NB=NBJ(1,NE)
            DO 3060 NI=1,NIT(NB)
              XID(NI,ND)=XI(NI)
 3060       CONTINUE
          ENDIF
 3070   CONTINUE
        LFMAX=NET(1)
        FORMAT='(/'' Data was read from the terminal'')'
        WRITE(OP_STRING,FORMAT)
      	CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      CALL EXITS('IPDATA')
      RETURN
 9999 CALL TRIM(ERROR,IBEG,IEND)
      ERROR=ERROR(IBEG:IEND)//' >IPDATA'
      RETURN 1
      END


      SUBROUTINE LNLXI(INP,LN,NPE,NXI)

C**** Note: This is now redundant, having been superceded by NENXI
C****  - PJH 4-Jan-1991
C**** The elements surrounding each element LN(L),L=1,LN(0)
C**** are returned in LXI (NXI(-1,L) & NXI(1,L) for 1D;  NXI(-2..2,L)
C**** for 2D).If an adjacent line or face does not exist NXI(i,l) is
C**** set to 0.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER INP(NNM,NIM,*),LN(0:*),NPE(NNM,NBM,*),NXI(-NIM:NIM,0:*)
!     Local Variables
      CHARACTER ERROR*10 
      INTEGER L,LS,N11,N12,N21,N22,NB,NBS,NE,NES,NI,NN,NP11,NP11S,
     '  NP12,NP12S,NP21,NP21S,NP22,NP22S

C     CALL ENTERS('LNLXI',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'(A)') ' >LNLXI'
      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      DO 50 L=1,LN(0)
        DO 5 NI=-NIM,NIM
          NXI(NI,L)=0
 5      CONTINUE
C****   Note temporary nb
        NB=1
        NE=LN(L)
        IF(NIT(NB).EQ.1) THEN
C         NP1=NPE(1,NB,NE)
C         IF(NPE(2,NB,NE).EQ.1) NP2=NPE(7,NB,NE)
C         IF(NPE(2,NB,NE).EQ.2) NP2=NPE(8,NB,NE)
C         IF(NPE(2,NB,NE).EQ.3) NP2=NPE(9,NB,NE)
C         IF(NPE(2,NB,NE).EQ.4) NP2=NPE(7,NB,NE)
C         DO 30 LS=1,LN(0)
C           IF(LS.NE.L) THEN
C             NLS=LN(LS)
C             NP1S=NPL(6,NLS)
C             IF(NPL(2,NLS).EQ.1) NP2S=NPL(7,NLS)
C             IF(NPL(2,NLS).EQ.2) NP2S=NPL(8,NLS)
C             IF(NPL(2,NLS).EQ.3) NP2S=NPL(9,NLS)
C             IF(NPL(2,NLS).EQ.4) NP2S=NPL(7,NLS)
C             IF(NP1.EQ.NP2S) NXI(-1,L)=-LS
C             IF(NP2.EQ.NP1S) NXI( 1,L)=-LS
C           ENDIF
C30       CONTINUE
        ELSE IF(NIT(NB).EQ.2) THEN
          DO 10 NN=1,NNT(NB)
            IF((INP(NN,1,NB).EQ.1).AND.(INP(NN,2,NB).EQ.1)) N11=NN
            IF((INP(NN,1,NB).EQ.1).AND.(INP(NN,2,NB).EQ.2)) N12=NN
            IF((INP(NN,1,NB).EQ.2).AND.(INP(NN,2,NB).EQ.1)) N21=NN
            IF((INP(NN,1,NB).EQ.2).AND.(INP(NN,2,NB).EQ.2)) N22=NN
 10       CONTINUE
          NP11=NPE(N11,NB,NE)
          NP12=NPE(N12,NB,NE)
          NP21=NPE(N21,NB,NE)
          NP22=NPE(N22,NB,NE)
          DO 20 LS=1,LN(0)
            IF(LS.NE.L) THEN
              NES=LN(LS)
C****         Note temporary nb
              NBS=1
              DO 15 NN=1,NNT(NBS)
                IF((INP(NN,1,NBS).EQ.1).AND.(INP(NN,2,NBS).EQ.1)) N11=NN
                IF((INP(NN,1,NBS).EQ.1).AND.(INP(NN,2,NBS).EQ.2)) N12=NN
                IF((INP(NN,1,NBS).EQ.2).AND.(INP(NN,2,NBS).EQ.1)) N21=NN
                IF((INP(NN,1,NBS).EQ.2).AND.(INP(NN,2,NBS).EQ.2)) N22=NN
 15           CONTINUE
              NP11S=NPE(N11,NBS,NES)
              NP12S=NPE(N12,NBS,NES)
              NP21S=NPE(N21,NBS,NES)
              NP22S=NPE(N22,NBS,NES)
              IF(NP11.EQ.NP12S.AND.NP21.EQ.NP22S) NXI(-2,L)=LS
              IF(NP11.EQ.NP21S.AND.NP12.EQ.NP22S) NXI(-1,L)=LS
              IF(NP21.EQ.NP11S.AND.NP22.EQ.NP12S) NXI( 1,L)=LS
              IF(NP12.EQ.NP11S.AND.NP22.EQ.NP21S) NXI( 2,L)=LS
            ENDIF
 20       CONTINUE
        ENDIF
C       DO 40 NI=1,NIM
C         IF(NXI(-NI,L).EQ.0) NXI(-NI,L)=-L
C         IF(NXI( NI,L).EQ.0) NXI( NI,L)=-L
C40     CONTINUE
 50   CONTINUE

C     CALL EXITS('LNLXI')
 9999 RETURN
      END


      SUBROUTINE NEWXID(ITMAX,LD,LN,NBJ_face,
     '  NDDL,NDLT,NJP,NKE,NPF,NPNE,NPL,NQE,NRE,NVJE,NXI,
     '  SE,SQ,XA,XE,XID,XP,ZD,ERROR,*)

C#### Subroutine: NEWXID
C###  Description:
C###    NEWXID calculates Xi coordinates corresponding to closest 
C###    approach of a data point to a region defined by XP.

C**** If any Xi lie outside [0,1] the region is redefined by XP for a
C**** segment, if it exists, corresponding to the Xi coordinate most
C**** violated in the previous segment.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
!     Parameter List
      INTEGER ITMAX,LD(NDM),LN(0:NEM),NBJ_face(NJM),NDDL(NEM,NDEM),
     '  NDLT(NEM),NJP(NPM),NKE(NKM,NNM,NBFM,NEFM),NPF(15,NFM),
     '  NPNE(NNM,NBFM,NEFM),NPL(20,NLM),NQE(NSM,NBFM,NEFM),NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEFM),NXI(-NIM:NIM,0:NEM)
      REAL*8 SE(NSM,NBFM,NEFM),SQ(NDM),XA(NAM,NJM,NQM),XE(NSM,NJM),
     '  XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IT,l,l1,L2,l3,LMAX,nd,ndf1,NDL,ndl1,ne,
     '  NF1,NF2,ni,NIOUT(4),
     '  nj,njj,njj1,njj2,NJOT,nl,NL1,nou,NOUT,nr
      REAL*8 ABSOUT(8)

      CALL ENTERS('NEWXID',*9999)
C GMH 30/8/95 I dont think that IT is used at the moment,
C             but for now I will just initialise it to
C             stop FTNCHEK from spewing.
      IT=0
      IF(DOP) THEN
        WRITE(OP_STRING,'(A)') ' >NEWXID'
      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
CGMH29/8/95 Unused      VMAX=0.25D0
      DO 9 l1=1,LN(0)
        IF(DOP) THEN
          WRITE(OP_STRING,'('' L1='',I4,'' LN(L1)='',I4)') l1,LN(l1)
      	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(LN(l1).GT.0) THEN
          IF(DOP) THEN
            WRITE(OP_STRING,'('' LXI='',5I4)')
     '        ((NXI(ni,LMAX),ni=-2,2),LMAX=1,LN(0))
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          NF1=LN(l1)
          NJOT=3
          DO 5 ndf1=1,NDLT(NF1)
            nd=NDDL(NF1,ndf1)
            L2=l1
            DO 4 l3=1,LN(0)
              NF2=LN(L2)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' NDF1='',I4,'' nd='',I4,'
     '            //''' L2='',I4,'' LN(L2)='',I4)') ndf1,nd,L2,NF2
      		CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              ne=NPF(6,nf2)
              nr=NRE(ne)
              njj=0
              DO njj1=1,2  !geometry/fibres
                DO njj2=1,NJ_LOC(njj1,0)
                  nj=NJ_LOC(njj1,njj2)
                  njj=njj+1
                  NBJ_face(nj)=NPF(9+njj,nf2)
                ENDDO ! njj2
              ENDDO ! njj1
              CALL XPXE(NBJ_face,NKE(1,1,1,NF2),NPF(1,NF2),
     '          NPNE(1,1,NF2),NQE(1,1,NF2),NRE(nf2),NVJE(1,1,1,ne),
     '          SE(1,1,NF2),XA,XE,XP,ERROR,*9999)
c Note: PJH 21Nov94 The following commented calls need updating
              IF(ITYP10(nr).EQ.1) THEN
c               CALL CLOS21(IBT,IDO,INP,IT,ITMAX,NJOT,NKE(1,1,1,NF2),
c    '            NPF(1,NF2),NPNE(1,1,NF2),NQE(1,1,NF2),
c    '            SE(1,1,NF2),SQ(nd),TOL,VMAX,XE,XID(1,nd),ZD(1,nd))
              ELSE IF(ITYP10(nr).EQ.2) THEN
c               CALL CLOS22(IBT,IDO,INP,IT,ITMAX,NPF(1,NF2),SQ(nd),
c    '            TOL,VMAX,XE,XID(1,nd),ZD(1,nd))
              ELSE IF(ITYP10(nr).EQ.3) THEN
c               CALL CLOS23(IBT,IDO,INP,IT,ITMAX,NJOT,NPF(1,NF2),
c    '            SQ(nd),TOL,VMAX,XE,XID(1,nd),ZD(1,nd))
              ELSE IF(ITYP10(nr).EQ.4) THEN
c               CALL CLOS24(IBT,IDO,INP,IT,ITMAX,NJOT,NKE(1,1,1,NF2),
c    '            NPF(1,NF2),NPNE(1,1,NF2),NQE(1,1,NF2),
c    '            SE(1,1,NF2),SQ(nd),TOL,VMAX,XE,XID(1,nd),ZD(1,nd))
              ENDIF
              IF(IT.GE.ITMAX) THEN
                WRITE(OP_STRING,
     '            '('' WARNING: Convergence not reached in CLOS2X''/'
     '            //''' nd='',I5,'' ZD='',(E12.6,1X)/'
     '            //''' LD='',I4,'' XID'',(E12.6,1X),''SQ='',E12.6/)')
     '            nd,(ZD(nj,nd),nj=1,NJOT),
     '            LD(nd),(XID(ni,nd),ni=1,2),SQ(nd)
      		CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
              NOUT=0
              DO 1 ni=1,2
                IF(XID(ni,nd).LT.0.d0) THEN
                  NOUT=NOUT+1
                  NIOUT(NOUT)=-ni
                  ABSOUT(NOUT)=-XID(ni,nd)
                ELSE IF(XID(ni,nd).GT.1.d0) THEN
                  NOUT=NOUT+1
                  NIOUT(NOUT)=ni
                  ABSOUT(NOUT)=XID(ni,nd)-1.d0
                ENDIF
    1         CONTINUE
              IF(NOUT.EQ.0) THEN
                LD(nd)=LN(L2)
                GO TO 5
              ENDIF
              CALL RSORT(NOUT,ABSOUT,NIOUT)
              DO 2 nou=NOUT,1,-1
                IF(ABS(NXI(NIOUT(nou),L2)).NE.L2) GO TO 3
                IF(nou.EQ.1) THEN
                  LD(nd)=LN(L2)
                  GO TO 5
                ENDIF
    2         CONTINUE
    3         IF(DOP) THEN
                WRITE(OP_STRING,'('' LXI='',7I4)') 
     '            (NXI(ni,L2),ni=-2,2)
      		CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              L2=ABS(NXI(NIOUT(nou),L2))
              ni=ABS(NIOUT(nou))
              IF(XID(ni,nd).LT.0.d0) THEN
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' XID(ni,nd)='',E12.4)') XID(ni,nd)
      		  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                XID(ni,nd)=DMAX1(XID(ni,nd)+1.d0,0.75D0)
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' Changed to '',E12.4)') XID(ni,nd)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ELSE IF(XID(ni,nd).GT.1.d0) THEN
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' XID(ni,nd)='',E12.4)') XID(ni,nd)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                XID(ni,nd)=DMIN1(XID(ni,nd)-1.d0,0.25D0)
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' Changed to '',E12.4)') XID(ni,nd)
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF
              IF(DOP) THEN
                WRITE(OP_STRING,'('' L2='',I2,'' NIOUT='',I2,'
     '            //''' XID(ni,nd)='',E12.6)') L2,NIOUT(nou),XID(ni,nd)
      		CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
    4       CONTINUE
    5     CONTINUE
        ELSE IF(LN(l1).LT.0) THEN
          NL1=ABS(LN(l1))
          NJOT=NJP(NPL(6,NL1))
          DO 8 ndl1=1,NDLT(NL1)
            nd=NDDL(NL1,ndl1)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' L1='',I4,'' NL1='',I4,'' NJOT='','
     '          //'I4,'' NDL1='',I4,'' nd='',I4,'' LD(nd)='',I4)')
     '          l1,NL1,NJOT,ndl1,nd,LD(nd)    
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            L2=l1
            DO 7 l3=1,LN(0)
!GMH29/8/95 Unused              NL2=ABS(LN(L2))
c Note: PJH 21Nov94 The following commented calls need updating
              IF(ITYP10(nr).EQ.1) THEN
c               CALL CLOS11(IT,ITMAX,NPL(1,NL2),DL(1,NL2),
c    '                      SQ(nd),TOL,VMAX,XP,XID(1,nd),ZD(1,nd))
              ELSE IF(ITYP10(nr).EQ.2) THEN
c               CALL CLOS12(IT,ITMAX,NPL(1,NL2),nr,DL(1,NL2),
c    '                      SQ(nd),TOL,VMAX,XP,XID(1,nd),ZD(1,nd))
              ELSE IF(ITYP10(nr).EQ.3) THEN
c               CALL CLOS13(IT,ITMAX,NPL(1,NL2),nr,DL(1,NL2),
c    '                      SQ(nd),TOL,VMAX,XP,XID(1,nd),ZD(1,nd))
              ELSE IF(ITYP10(nr).EQ.4) THEN
c               CALL CLOS14(IT,ITMAX,NPL(1,NL2),nr,DL(1,NL2),
c    '                      SQ(nd),TOL,VMAX,XP,XID(1,nd),ZD(1,nd))
              ENDIF
              IF(IT.GE.ITMAX) THEN
                WRITE(OP_STRING,
     '            '('' WARNING: Convergence not reached in CLOS1X''/'
     '            //''' nd='',I5,'' ZD='',(E12.6,1X)/'
     '            //''' LD='',I4,'' XID'',(E12.6,1X),''SQ='',E12.6/)')
     '          nd,(ZD(nj,nd),nj=1,NJOT),
     '          LD(nd),(XID(ni,nd),ni=1,1),SQ(nd)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
              NOUT=0
              IF(XID(1,nd).LT.0.d0) THEN
                NOUT=NOUT+1
                NIOUT(NOUT)=-1
                ABSOUT(NOUT)=-XID(1,nd)
              ELSE IF(XID(1,nd).GT.1.d0) THEN
                NOUT=NOUT+1
                NIOUT(NOUT)=1
                ABSOUT(NOUT)=XID(1,nd)-1.d0
              ELSE
                LD(nd)=LN(L2)
                GO TO 8
              ENDIF
              IF(ABS(NXI(NIOUT(NOUT),L2)).NE.L2) GO TO 6
              LD(nd)=LN(L2)
              GO TO 7
    6         L2=ABS(NXI(NIOUT(NOUT),L2))
              ni=ABS(NIOUT(NOUT))
              IF(XID(ni,nd).LT.0.d0) XID(ni,nd)=DMAX1(XID(ni,nd)+1.d0,
     '          0.75D0)
              IF(XID(ni,nd).GT.1.d0) XID(ni,nd)=DMIN1(XID(ni,nd)-1.d0,
     '          0.25D0)
    7       CONTINUE
    8     CONTINUE
        ENDIF
    9 CONTINUE
      DO nl=1,NLM
        NDLT(nl)=0
      ENDDO
      DO nd=1,NDT
        nl=ABS(LD(nd))
        NDLT(nl)=NDLT(nl)+1
        NDL=NDLT(nl)
        NDDL(nl,NDL)=nd
      ENDDO
      IF(DOP) THEN
        DO nd=1,NDT
          IF(LD(nd).GT.0) THEN
            WRITE(OP_STRING,'(''     nd='',I4,4X,''LD(nd)='',I4,9X,'
     '        //'''XID(ni,nd)='',2(E12.6,4X),4X,''SQ(nd)='',E12.6)')
     '        nd,LD(nd),(XID(ni,nd),ni=1,2),SQ(nd)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ELSE IF(LD(nd).LT.0) THEN
            WRITE(OP_STRING,'(''     nd='',I4,4X,''LD(nd)='',I4,9X,'
     '        //'''XID(1,nd)='',E12.6,4X,''SQ(nd)='',E12.6)')
     '        nd,LD(nd),XID(1,nd),SQ(nd)  
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO
        DO l=1,LN(0)
          nl=ABS(LN(l))
          WRITE(OP_STRING,'('' L='',I4,'' nl='',I4,'
     '      //''' NDLT(nl)='',I4)') l,nl,NDLT(nl)  
      	  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          IF(NDLT(nl).GT.0) THEN
            WRITE(OP_STRING,'(12X,''NDDL(nl,ndl)='',10I4,'
     '        //'16(/24X,10I4))') (NDDL(nl,NDL),NDL=1,NDLT(nl))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('NEWXID')
      RETURN
 9999 CALL ERRORS('NEWXID',ERROR)
      CALL EXITS('NEWXID')
      RETURN 1
      END

      SUBROUTINE NY_NO(NJO,NKJ,NONY,NPO,nr,nx,NYNO,CONY,ERROR,*)

C#### Subroutine: NY_NO
C###  Description:
C###    NY_NO calculates optimisation dofs no from mesh dofs ny for 
C###    region nr.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:fit000.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER NJO,NKJ(NJM,NPM),NONY(0:NOYM,NYM,NRCM),NPO(0:*),nr,nx,
     '  NYNO(0:NYOM,NOOPM,NRCM)
      REAL*8 CONY(0:NOYM,NYM,NRCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER n,nj,njj,nk,no,noy,np,nrc,ny

      CALL ENTERS('NY_NO',*9999)
      nrc=1 ! temporary
      ny=0
      no=0
      DO n=1,NPO(0)              !loop over all nodes assoc with data
        np=NPO(n)
        DO njj=1,NUM_FIT(0)       !loop over all fit variables nj
          nj=NLH_FIT(1,1,njj)
          DO nk=1,NKJ(nj,np)
            ny=ny+1
            DO noy=1,NONY(0,ny,nrc)    !If NONY(0,ny,nrc)>0 (ny in fit)
              IF(nj.EQ.NJO) THEN !.. & nj=current variable njo
                no=no+1          !.. then increment no
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' np='',I4,'' nj='',I1,'' nk='','
     '              //'I1,'' ny='',I4,'' noy='',I4,'' no='',I4)')
     '              np,nj,nk,ny,noy,no
      		  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                NONY(noy,ny,nrc)=no    !.. & record
                CONY(noy,ny,nrc)=1.d0 !.. & set coupling coeff to 1
                NYNO(0,no,nrc)=1
                NYNO(1,no,nrc)=ny
              ELSE
                NONY(noy,ny,nrc)=0
                CONY(noy,ny,nrc)=0.d0
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      NOT(1,1,nr,nx)=no                     !.. & record total
      IF(DOP) THEN
        WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5)') NOT(1,1,nr,nx)
      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('NY_NO')
      RETURN
 9999 CALL ERRORS('NY_NO',ERROR)
      CALL EXITS('NY_NO')
      RETURN 1
      END
      

      SUBROUTINE PROJ21(IT,ITMAX,NBJ,NE,SQ,XE,XI,ZD,ZDD,ERROR,*)

C**** Finds the XI-coordinates at the closest approach of a 2D element to a
C**** line connecting data point ZD with ZDD. Uses e04lbf.
C**** NOTE: Works only for 2-D elements in rectangular cartesian coords.
C**** Note: E04LBF calls PROJ_FUNC etc which need IBT, which is passed 
C**** via a common block edata00 as IBT2. This fudge could be eliminated
C**** by using a different Nag routine (E04UPF ?) which allows user parameters
C**** to be passed through.

C
C cpb 7/10/94 edata00.cmn has also be passed to archive. If this routine
C is ever reused it will need to be done through e04upf to avoid using
C edata00
C 

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:edata00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IT,ITMAX,NBJ(*),NE
      REAL*8 SQ,XE(NSM,*),XI(*),ZD(*),ZDD(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBOUND,IFAIL,ISTATE(3),IPRINT,IWORK(2),LH,LIW,LW,NJ,NS,
     '  NUM_INDEP
      REAL*8 BOUND_LOWER(3),BOUND_UPPER(3),ETA,HESL(3),HESD(3),
     '  FDERIV(3),STEPMX,WORK(24),XTOL
      EXTERNAL PROJ_FUNCT,PROJ_HESS,PROJ_MONIT

      DATA NUM_INDEP/3/,IPRINT/20/,ETA/0.5/,XTOL/1.0D-05/,STEPMX/5.0D0/,
     '  IBOUND/0/,BOUND_LOWER/0.0D0,0.0D0,0.0D0/,
     '  BOUND_UPPER/1.0D0,1.0D0,1.0D+6/,
     '  LH/3/,LIW/2/,LW/24/

      CALL ENTERS('PROJ21',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'(A)') ' >PROJ21 2-D rectangular cartesian'
      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
C     copy data into common block for proj_funct and proj_hess
      DO NJ=1,NJT
        U1(NJ)=ZDD(NJ)
        U2(NJ)=ZD(NJ)
        DO NS=1,NST(NBJ(NJ))
          XENE(NS,NJ)=XE(NS,NJ)
        ENDDO
      ENDDO
      IFAIL=1
      CALL E04LBF(NUM_INDEP,PROJ_FUNCT,PROJ_HESS,PROJ_MONIT,IPRINT,
     '  ITMAX,ETA,XTOL,
     '  STEPMX,IBOUND,BOUND_LOWER,BOUND_UPPER,XI,HESL,LH,HESD,ISTATE,
     '  SQ,FDERIV,IWORK,LIW,WORK,LW,IFAIL)
      IF(IFAIL.NE.0.AND.DOP)THEN
        WRITE(OP_STRING,*)' E04LBF returns IFAIL=',IFAIL
      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IT=NITER
      IF(XI(3).LT.1.0D-02)IT=ITMAX !wrong element flag

 9998 CALL EXITS('PROJ21')
      RETURN
 9999 CALL ERRORS('PROJ21',ERROR)
      CALL EXITS('PROJ21')
      RETURN 1
      END


      SUBROUTINE PROJ_FUNCT(IFLAG,NUM_INDEP,XI,F,FDERIV,IW,LIW,W,LW)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:edata00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IFLAG,IW(*),LIW,LW,NUM_INDEP
      REAL*8 F,FDERIV(*),W(*),XI(*)
!     Local Variables
      INTEGER NB,NJ
      REAL*8 DUDT(3),DVDS1(3),DVDS2(3),PFXI,U(3),V(3)

C     calculate position on line at t=xi(3)
      DO NJ=1,NJT
        U(NJ)=U1(NJ)+XI(3)*(U2(NJ)-U1(NJ))
        DUDT(NJ)=(U2(NJ)-U1(NJ))
      ENDDO
C     calculate position on element at xi=xi(2),xi(3)
      DO NJ=1,NJT
        NB=NBJ2(NJ,N_E)
        V(NJ)=PFXI(IBT2(1,1,NB),IDO2(1,0,NB),INP2(1,1,NB),NAN2,
     '    NB,1,XENE(1,NJ),XI)
        DVDS1(NJ)=PFXI(IBT2(1,1,NB),IDO2(1,0,NB),INP2(1,1,NB),NAN2,
     '    NB,2,XENE(1,NJ),XI)
        DVDS2(NJ)=PFXI(IBT2(1,1,NB),IDO2(1,0,NB),INP2(1,1,NB),NAN2,
     '    NB,4,XENE(1,NJ),XI)
      ENDDO
      F=0.0D0
      FDERIV(1)=0.0D0
      FDERIV(2)=0.0D0
      FDERIV(3)=0.0D0
      DO NJ=1,NJT
        F=F+(U(NJ)-V(NJ))*(U(NJ)-V(NJ))
        FDERIV(1)=FDERIV(1)-2*(U(NJ)-V(NJ))*DVDS1(NJ)
        FDERIV(2)=FDERIV(2)-2*(U(NJ)-V(NJ))*DVDS2(NJ)
        FDERIV(3)=FDERIV(3)+2*(U(NJ)-V(NJ))*DUDT(NJ)
      ENDDO

      RETURN
      END


      SUBROUTINE PROJ_HESS(IFLAG,NUM_INDEP,XI,FHESL,LH,FHESD,IW,LIW,W,
     '  LW)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:edata00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IFLAG,IW(*),NUM_INDEP,LH,LIW,LW
      REAL*8 FHESL(*),FHESD(*),W(*),XI(*)
!     Local Variables
      INTEGER NB,NJ
      REAL*8 DUDT(3),DVDS1(3),DVDS2(3),D2VDS1(3),D2VDS2(3),D2VDS12(3),
     '  PFXI,U(3),V(3)

C     calculate position on line at t=xi(3)
      DO NJ=1,NJT
        U(NJ)=U1(NJ)+XI(3)*(U2(NJ)-U1(NJ))
        DUDT(NJ)=(U2(NJ)-U1(NJ))
      ENDDO
C     calculate position on element at xi=xi(2),xi(3)
      DO NJ=1,NJT
        NB=NBJ2(NJ,N_E)
        V(NJ)=PFXI(IBT2(1,1,NB),IDO2(1,0,NB),INP2(1,1,NB),NAN2,
     '    NB,1,XENE(1,NJ),XI)
        DVDS1(NJ)=PFXI(IBT2(1,1,NB),IDO2(1,0,NB),INP2(1,1,NB),NAN2,
     '    NB,2,XENE(1,NJ),XI)
        DVDS2(NJ)=PFXI(IBT2(1,1,NB),IDO2(1,0,NB),INP2(1,1,NB),NAN2,
     '    NB,4,XENE(1,NJ),XI)
        D2VDS1(NJ)=PFXI(IBT2(1,1,NB),IDO2(1,0,NB),INP2(1,1,NB),NAN2,
     '    NB,3,XENE(1,NJ),XI)
        D2VDS2(NJ)=PFXI(IBT2(1,1,NB),IDO2(1,0,NB),INP2(1,1,NB),NAN2,
     '    NB,5,XENE(1,NJ),XI)
        D2VDS12(NJ)=PFXI(IBT2(1,1,NB),IDO2(1,0,NB),INP2(1,1,NB),NAN2,
     '    NB,6,XENE(1,NJ),XI)
      ENDDO
      FHESL(1)=0.0D0
      FHESL(2)=0.0D0
      FHESL(3)=0.0D0
      DO NJ=1,NJT
        FHESL(1)=FHESL(1)-2*((U(NJ)-V(NJ))*D2VDS12(NJ)-DVDS2(NJ)
     '    *DVDS1(NJ))
        FHESL(2)=FHESL(2)-2*DVDS1(NJ)*DUDT(NJ)
        FHESL(3)=FHESL(3)-2*DVDS2(NJ)*DUDT(NJ)
      ENDDO
      FHESD(1)=0.0D0
      FHESD(2)=0.0D0
      FHESD(3)=0.0D0
      DO NJ=1,NJT
        FHESD(1)=FHESD(1)-2*((U(NJ)-V(NJ))*D2VDS1(NJ)-DVDS1(NJ)
     '    *DVDS1(NJ))
        FHESD(2)=FHESD(2)-2*((U(NJ)-V(NJ))*D2VDS2(NJ)-DVDS2(NJ)
     '    *DVDS2(NJ))
        FHESD(3)=FHESD(3)+2*DUDT(NJ)*DUDT(NJ)
      ENDDO

      RETURN
      END


      SUBROUTINE PROJ_MONIT(NUM_INDEP,XI,F,FDERIV,ISTATE,GPJNRM,COND,
     '  POSDEF,NUMITER,NF,IW,LIW,W,LW)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:edata00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER ISTATE(*),IW(*),LIW,LW,NF,NUM_INDEP,NUMITER
      REAL*8 XI(*),F,FDERIV(*),GPJNRM,COND,W(*)
      LOGICAL POSDEF
!     Local Variables
      CHARACTER ERROR*10
      
      IF(DOP)THEN
        WRITE(OP_STRING,*)' NITER=',NUMITER
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*)' XI=',XI(1),XI(2),XI(3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*)' F=',F
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*)' COND=',COND
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*)' FDERIV=',FDERIV(1),FDERIV(2),FDERIV(3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      NITER=NUMITER

 9999 RETURN
      END


      SUBROUTINE XPES(IBT,IDO,INP,NBJ,NDDL,NDLT,ne,NJE,NJV,
     '  ER,ES,SE,WD,WDL,XID,XIDL,ZD,ZDL,ERROR,*)

C#### Subroutine: XPES
C###  Description:
C###    Evaluates element stiffness matrix ES(ms,ns) and r.h.s. ER(ns)
C###    in calculation of least squares fit of linear field variables,
C###    defined by nodal values XP(nk,nv,nj,np), to the set of data
C###    values XD(nj,nd) with weights WD(nj,nd) at local coordinate
C###    values XID(ni,nd), where nj=NJE(1)+njv.
C**** NDLT         is  no. data points within current element ne
C**** NDDL(ne,nde)  "  global data pt no. of element data pt nde
C**** ZDL(nh,nde)   "  r.c. coords & value of   "      "      "
C**** WDL(nh,nde)   "  weighting factor for     "      "      "
C**** XIDL(ni,nde)  "  Xi-coordinate of         "      "      "

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),NBJ(NJM),
     '  NDDL(NEM,NDEM),NDLT,ne,NJE,NJV
      REAL*8 ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),SE(NSM,NBFM),WD(NJM,NDM),
     '  WDL(NHM,NDEM),XID(NIM,NDM),XIDL(NIM,NDEM),ZD(NJM,NDM),
     '  ZDL(NHM,NDEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nde,NITB,nj,nk1,nk2,nn1,nn2,ns1,ns2
      REAL*8 PSI1,SYM1,SYM2

      CALL ENTERS('XPES',*9999)

C CPB 19/4/94 Can't find where this subroutine is called from ???

      IF(DOP) THEN
        WRITE(OP_STRING,'(A)') ' >XPES'
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
C CPB 19/4/94 This nj location needs to be generalised
      nj=NJE+NJV
      nb=NBJ(nj)
      NITB=NIT(nb)
      CALL ZDZDL(1,NDDL,NDLT,NITB,ne,nx,WD,WDL,XID,XIDL,ZD,ZDL,
     '  ERROR,*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'(I3,4E11.3)') NDDL(ne,1),ZDL(3,1),WDL(3,1),
     '    XIDL(1,1),XIDL(2,1)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      ns1=0
      DO nn1=1,NNT(nb)
        DO nk1=1,NKT(nn1,nb)
          ns1=ns1+1
          SYM1=0.d0
          DO nde=1,NDLT
            SYM1=SYM1+PSI1(IBT,IDO,INP,nb,1,nk1,nn1,XIDL(1,nde))
     '        *ZDL(nj,nde)*WDL(nj,nde)
          ENDDO
          ER(ns1)=ER(ns1)+SYM1*SE(ns1,nb)
          ns2=0
          DO nn2=1,NNT(nb)
            DO nk2=1,NKT(nn2,nb)
              ns2=ns2+1
              SYM2=0.d0
              DO nde=1,NDLT
                SYM2=SYM2+PSI1(IBT,IDO,INP,nb,1,nk1,nn1,XIDL(1,nde))
     '            *PSI1(IBT,IDO,INP,nb,1,nk2,nn2,XIDL(1,nde))
     '            *WDL(nj,nde)
              ENDDO
              ES(ns1,ns2)=ES(ns1,ns2)+SYM2*SE(ns1,nb)*SE(ns2,nb)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      CALL EXITS('XPES')
      RETURN
 9999 CALL ERRORS('XPES',ERROR)
      CALL EXITS('XPES')
      RETURN 1
      END


      SUBROUTINE ZDERF(IBT,IDO,INP,NBJ,NDDL,NDLT,ne,NJO,
     '  ER,SE,WD,WDL,XE,XID,XIDL,ZD,ZDL,ERROR,*)

C**** Evaluates element r.h.s. ER(ns) for Fourier fits
C**** in calculation of least squares fit of linear field variables,
C**** defined by nodal values XP(nk,nv,nj,np), to the set of data
C**** values XD(nj,nd) with weights WD(nj,nd) at local coordinate
C**** values XID(ni,nd), where nj=NJE(1)+NJO.
C**** NDLT         is  no. data points within current element ne
C**** NDDL(ne,nde)  "  global data pt no. of element data pt nde
C**** ZDL(nj,nde)   "  r.c. coords & value of   "      "      "
C**** WDL(nj,nde)   "  weighting factor for     "      "      "
C**** XIDL(ni,nde)  "  Xi-coordinate of         "      "      "

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:fit000.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
!     Parameter List
      INTEGER  IBT(2,NIM,*),IDO(NKM,0:NIM,*),INP(NNM,NIM,*),
     '  NBJ(NJM,*),NDDL(NEM,*),NDLT,ne,NJO
      REAL*8 ER(*),SE(NSM,*),WD(NJM,*),WDL(NJM,*),XE(NSM,*),
     '  XID(NIM,*),XIDL(NIM,*),ZD(NJM,*),ZDL(NJM,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER n1s,nb,NBG,nde,NITB,nj,nk1,nn1,ns1,NSA
      REAL*8 PSI1,PXI,SUM1,X(3),Z(3)

      CALL ENTERS('ZDERF',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'(A)') ' >ZDERF'
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      nb=NBJ(NJO,ne) !Fourier basis
      NITB=NIT(nb)
      DO n1s=1,NST(nb)
        ER(n1s)=0.0D0
      ENDDO
      CALL ZDZDL(NDDL,NDLT,NITB,ne,WD,WDL,XID,XIDL,ZD,ZDL,ERROR,*9999)
      IF(KTYP58.EQ.2)THEN  !displacement
        DO nde=1,NDLT
          DO nj=1,NJT
            NBG=NBJ(nj,ne)
            X(nj)=PXI(IBT(1,1,NBG),IDO(1,0,NBG),INP(1,1,NBG),NBG,1,
     '        XIDL(1,nde),XE(1,nj))
          ENDDO
          CALL XZ(ITYP10(1),X,Z)
          ZDL(NJO,nde)=ZDL(NJO,nde)-Z(NJ_FIT(1,1))
          IF(DOP) THEN
            WRITE(OP_STRING,'('' nde,zdl(njo,nde):'',I5,E13.4)')
     '        nde,ZDL(NJO,nde)
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF
      ns1=0
      DO nn1=1,NNT(nb)
        DO nk1=1,NKT(nn1,nb)
          ns1=ns1+1
          SUM1=0.0D0
          DO nde=1,NDLT
            SUM1=SUM1+PSI1(IBT(1,1,nb),IDO(1,0,nb),INP(1,1,nb),nb,1,
     '        nk1,nn1,XIDL(1,nde))*ZDL(NJO,nde)*WDL(NJO,nde)
          ENDDO
          NSA=(nn1-1)*NKT(0,NBO)+MOD(ns1-1,NKT(0,NBO))+1  !spatial dof (for SE)
          ER(ns1)=ER(ns1)+SUM1*SE(NSA,NBO)
        ENDDO
      ENDDO

      CALL EXITS('ZDERF')
      RETURN
 9999 CALL ERRORS('ZDERF',ERROR)
      CALL EXITS('ZDERF')
      RETURN 1
      END


      SUBROUTINE ZDESF(IBT,IDO,INP,nb,NDLT,ne,NJO,ES,PG,SE,WDL,WG,WU,
     '  XIDL,ERROR,*)

C**** Evaluates element stiffness matrix ES(ms,ns)
C**** in calculation of least squares fit of FOURIER field variables,

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:fit000.cmn'
      INCLUDE 'cmiss$reference:four00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
!     Parameter List
      INTEGER IBT(2,*),IDO(NKM,0:*),INP(NNM,*),nb,NDLT,ne,NJO
      REAL*8 ES(NHM*NSM,*),PG(NSM,NUM,NGM,*),SE(NSM,*),WDL(NJM,*),
     '  WG(NGM,*),WU(0:*),XIDL(NIM,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER na1,na2,nde,ng,NITB,nk1,nk2,NKF1,NKF2,nn1,nn2,ns1,ns2,
     '  NSA,NSB,nu
      REAL*8 F,PRECALC(44,4,4000),PSI1,SMOOTH,SUM

      CALL ENTERS('ZDESF',*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'(A)') ' >ZDESF'
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*) ' nb,NUT(NBO)=',nb,NUT(NBO)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*) ' NST(nb)=',NST(nb)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      NITB=NIT(nb)
      DO ns1=1,NST(nb)
        DO ns2=1,ns1
          ES(ns1,ns2)=0.0D0
        ENDDO
      ENDDO
      DO nn1=1,NNT(nb)
        DO nk1=1,NKT(nn1,nb)
          DO nde=1,NDLT
            PRECALC(nk1,nn1,nde)=PSI1(IBT,IDO,INP,nb,1,nk1,nn1,
     '        XIDL(1,nde))
          ENDDO
        ENDDO
      ENDDO
      ns1=0
      DO nn1=1,NNT(NBO)              !loop over nodes in element
        NKF1=0
        DO na1=1,IBT(2,NIT(nb))      !loop over terms in Fourier series
          DO nk1=1,NKT(nn1,NBO)          !loop over spatial derivs per node
            NKF1=NKF1+1              !nodal dof for nb (time basis)
            ns1=ns1+1                !element dof for nb
            NSA=(nn1-1)*NKT(0,NBO)+nk1 !element dof for NBO (geom. basis)
            ns2=0
            DO nn2=1,NNT(NBO)
              NKF2=0
              DO na2=1,IBT(2,NIT(nb))
                DO nk2=1,NKT(nn2,NBO)
                  NKF2=NKF2+1
                  ns2=ns2+1
                  NSB=(nn2-1)*NKT(0,NBO)+nk2
                  SMOOTH=0.0D0
                  IF(ktyp12.NE.0.AND.na1.EQ.na2)THEN !calculate smoothing
                    DO ng=1,NGT(NBO)
                      DO nu=2,NUT(NBO)
                        SMOOTH=SMOOTH+PG(NSA,nu,ng,NBO)
     '                          *PG(NSB,nu,ng,NBO)
     '                          *WU(nu-1)*WG(ng,NBO)
                      ENDDO
                    ENDDO
                    IF(na1.EQ.1)THEN
                      F=2.0D0*PI
                    ELSE
                      F=PI
                    ENDIF
                    SMOOTH=SMOOTH*F
                  ENDIF
                  SUM=0.0D0
                  DO nde=1,NDLT
                    SUM=SUM+PRECALC(NKF1,nn1,nde)*PRECALC(NKF2,nn2,nde)
     '                *WDL(NJO,nde)
                  ENDDO
                  ES(ns1,ns2)=(SMOOTH+SUM)*SE(NSA,NBO)*SE(NSB,NBO)
                  ES(ns2,ns1)=ES(ns1,ns2)
                ENDDO
              ENDDO
            ENDDO
            IF(DOP) THEN
              WRITE(OP_STRING,*)' finished ns=',ns1
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      CALL EXITS('ZDESF')
      RETURN
 9999 CALL ERRORS('ZDESF',ERROR)
      CALL EXITS('ZDESF')
      RETURN 1
      END

Module FE04
===========

      SUBROUTINE DIVH1(IBT,IDO,IDRN,INP,NBH,NBJ,NCO,ne,NEELEM,NHE,NJE,
     '  NJP,NKE,NKJ,NPF,NPNE,NPNODE,NPSTART,NQE,nr,NRE,NVJE,NVJP,NW,nx,
     '  CE,SE,XA,XE,XII,XP,NOCROSS,ERROR,*)

C#### Subroutine: DIVH1
C###  Description:
C###    DIVH1 subdivides Hermite tensor product element ne in the
C###    Xi(IDRN) direction at Xi(IDRN)=XII.

C**** IDR2 & IDR3 are the two orthogonal directions.
C**** NPSTART is first node number for new nodes
C**** NNIP(i1,i2,i3) records the element vertex number nn for each
C****   set of node position indices I1,I2,I3 where these are in the
C****   IDRN,IDR2 & IDR3 directions, respectively.
C**** NPNE1(nn) & NPNE2(nn) are temporary storage of global node numbers
C****   for two elements created by the subdivision.
C**** IV(nn,idrn) are vertices of new element corresponding to newly
C****   created nodes
C**** DSDXI(ni) are derivatives of arc-length wrt Xi

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),IDRN,
     '  INP(NNM,NIM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NCO(NEM),ne,
     '  NEELEM(0:NE_R_M,0:NRM),NHE(NEM),NJE(NEM),NJP(NPM),
     '  NKE(NKM,NNM,NBFM,NEFM),NKJ(NJM,NPM),NPF(15,NFM),
     '  NPNE(NNM,NBFM,NEFM),NPNODE(0:NP_R_M,0:NRM),NPSTART,
     '  NQE(NSM,NBFM,NEFM),nr,NRE(NEM),NW(NEM,2),nx,
     '  NVJE(NNM,NBFM,NJM,NEFM),NVJP(NJM,NPM)
      REAL*8 CE(NMM,NEM),SE(NSM,NBFM,NEFM),XA(NAM,NJM,NQM),XE(NSM,NJM),
     '  XII,XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL NOCROSS
!     Local Variables
      INTEGER i1,i2,i3,IDR(5),IDR2,IDR3,IK,k,m,n,nb,nb1,nc,
     '  ne_NEW,nh,nhx,ni,nj,njj,nk,nm,nn,nn_list(0:16),NNIP(4,4,4),
     '  NNTB,NODE,noelem,nonode,np,np1,np2,NPNE1(64),NPNE2(64),
     '  NPTNEW,NPTOLD,NP_TEST,nr2,ns,NTOT,nu,NUMI1,NUMI2,NUMI3,
     '  nv,nv1,nv2
      REAL*8 AA,DIFF,DSDXI(3),G1,G3,PXI,R,RC,RR,RRC,SLX,SMX,
     '  REFINE_GETXIPOS,SUM,X(11,6),XI(3),XS(4)
      LOGICAL EXIST,FOUND 

      DATA IDR/1,2,3,1,2/

      CALL ENTERS('DIVH1',*9999)

      nc=1 ! temporary cpb 22/11/94
      NPTNEW=NPSTART-1
      ne_NEW=NET(0)+1
      CALL ASSERT(ne_NEW.LE.NEM,'>>Too many elements -Increase NEMX',
     '  ERROR,*9999)

      nb1=NBJ(1,ne)
      NNTB=NNT(nb1)
      IF(NIT(nb1).EQ.1) THEN
        IDR2=0
        IDR3=0
      ELSEIF(NIT(nb1).EQ.2) THEN
        IDR3=0
        IF(IDRN.EQ.1) THEN
          IDR2=2
        ELSE
          IDR2=1
        ENDIF      
      ELSE 
        IDR2=IDR(IDRN+1)
        IDR3=IDR(IDRN+2)
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,*) ' ne= ',ne,' nb1= ',nb1,' NNT(nb1)= ',
     '	  NNT(nb1)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

C cpb 2/12/96 Fixing REFINE
      IF(NIT(nb1).EQ.1) THEN
        NUMI2=1
        NUMI3=1
      ELSE IF(NIT(nb1).EQ.2) THEN
        NUMI3=1
      ENDIF
      IF(NUMI3.NE.1) THEN
        IF(IBT(1,IDR3,nb1).EQ.1) THEN !Lagrange
          NUMI3=IBT(2,IDR3,nb1)+1
        ELSE IF(IBT(1,IDR3,nb1).EQ.2) THEN !Hermite
          NUMI3=2
        ELSE
          ERROR='>>Invalid basis type'
          GOTO 9999
        ENDIF
      ENDIF
      IF(NUMI2.NE.1) THEN
        IF(IBT(1,IDR2,nb1).EQ.1) THEN !Lagrange
          NUMI2=IBT(2,IDR2,nb1)+1
        ELSE IF(IBT(1,IDR2,nb1).EQ.2) THEN !Hermite
          NUMI2=2
        ELSE
          ERROR='>>Invalid basis type'
          GOTO 9999
        ENDIF
      ENDIF
      IF(IBT(1,IDRN,nb1).EQ.1) THEN !Lagrange
        NUMI1=IBT(2,IDRN,nb1)+1
      ELSE IF(IBT(1,IDRN,nb1).EQ.2) THEN !Hermite
        NUMI1=2
      ELSE
        ERROR='>>Invalid basis type'
        GOTO 9999
      ENDIF
C*** Find NNIP(i1,i2,i3), the element vertex numbers nn
C     IBT2=IBT(2,IDRN,nb1)
      DO i3=1,NUMI3
        DO i2=1,NUMI2
          DO i1=1,NUMI1
            NNIP(i1,i2,i3)=0
            IF(NUMI2.EQ.1.AND.NUMI3.EQ.1) THEN !1 Xi direction only
              DO nn=1,NNTB
                IF(INP(nn,IDRN,nb1).EQ.i1) NNIP(i1,i2,i3)=nn
              ENDDO !nn
            ELSE IF(NUMI3.EQ.1) THEN !2 Xi directions only              
              DO nn=1,NNTB
                IF(INP(nn,IDRN,nb1).EQ.i1.AND.
     '            INP(nn,IDR2,nb1).EQ.i2) NNIP(i1,i2,i3)=nn
              ENDDO !nn
            ELSE !3 Xi directions
              DO nn=1,NNTB
                IF(INP(nn,IDRN,nb1).EQ.i1.AND.
     '            INP(nn,IDR2,nb1).EQ.i2.AND.
     '            INP(nn,IDR3,nb1).EQ.i3) NNIP(i1,i2,i3)=nn
              ENDDO !nn
            ENDIF
            CALL ASSERT(NNIP(i1,i2,i3).NE.0,
     '        '>>Could not find local node',ERROR,*9999)
            IF(DOP) THEN
      	      WRITE(OP_STRING,'('' NNIP('',I1,'','',I1,'','
     '          //''',I1,'')= '',i2)') i1,i2,i3,NNIP(i1,i2,i3)
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO !i1
        ENDDO !i2
      ENDDO !i3

      CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '  NQE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),SE(1,1,ne),XA,XE,XP,
     '  ERROR,*9999)

C Create temporary version info in new element
      DO nb=1,NBFT
        DO nn=1,NNT(nb)
          DO nj=1,NJ_LOC(0,0,nr)
            NVJE(nn,nb,nj,ne_NEW)=NVJE(nn,nb,nj,ne)
          ENDDO !nj
        ENDDO !nn
      ENDDO !nb

C *** Loop over Xi directions
      m=0
      DO i3=1,NUMI3
        IF(NUMI3.EQ.1) THEN
          XI(3)=0.0d0
        ELSE
          XI(IDR3)=DBLE(i3-1)/DBLE(NUMI3-1)
        ENDIF
        DO i2=1,NUMI2
          IF(NUMI2.EQ.1) THEN
            XI(2)=0.0d0
          ELSE
            XI(IDR2)=DBLE(i2-1)/DBLE(NUMI2-1)
          ENDIF
          DO i1=1,NUMI1-1
            XI(IDRN)=REFINE_GETXIPOS(i1,IBT,IDRN,nb1,XII)
            DO nj=1,NJT !to find coords of proposed new node position
              nb=NBJ(nj,ne)
              XS(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '          nb,1,XI,XE(1,nj))
            ENDDO
            IF(DOP) THEN
              WRITE(OP_STRING,FMT='('' XI: '',3E12.3)') 
     '          (XI(ni),ni=1,3)
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,FMT='('' XS: '',3E12.3)')
     '		(XS(nj),nj=1,NJT)
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            m=m+1

            DO np=1,NPTNEW

C ***         Search for existing nodes
              IK=0
              DO nj=1,NJT
                DIFF= DABS(XS(nj)-XP(1,1,nj,np))
                IF((DIFF.LT.(1.d-3)* DABS(XS(nj)).OR.DIFF.LT.1.d-6).
     '          OR.(ITYP10(1).eq.2.and.nj.EQ.2.AND.(DABS(XP(1,1,1,np)).
     '            LT.1.d-6.OR. DABS(DIFF-2.0D0*PI).LT.1.d-6)).
     '          OR.(ITYP10(1).eq.3.and.nj.EQ.2.AND.(DABS(XP(1,1,1,np)).
     '            LT.1.d-6.OR. DABS(DIFF-2.0D0*PI).LT.1.d-6)).
     '          OR.(ITYP10(1).eq.3.and.nj.EQ.3.AND.(DABS(XP(1,1,1,np)).
     '            LT.1.d-6.OR. DABS(DIFF-2.0D0*PI).LT.1.d-6)).
     '          OR.(ITYP10(1).eq.4.and.nj.EQ.3.AND.(DABS(XP(1,1,2,np)).
     '            LT.1.d-6.OR. DABS(XP(1,1,2,np)-PI).LT.1.d-6.
     '          OR. DABS(DIFF-2.0D0*PI).LT.1.d-5))) THEN
                    IK=IK+1
                ENDIF
C               IF(DOP) THEN
C                 WRITE(IO4,'('' nj='',i1,'' ik='',i1,
C                   '' XP(1,1,2,np)'',E11.3,'' diff='',E11.3)')
C    '              nj,IK,XP(1,1,2,np),DIFF
C               ENDIF
              ENDDO
              IF(IK.EQ.NJT) THEN
                EXIST=.TRUE.
                NODE=np
                GO TO 41
              ELSE
                EXIST=.FALSE.
              ENDIF
            ENDDO

 41         IF(.NOT.EXIST) THEN !no existing node found
              NPTNEW=NPTNEW+1
              CALL ASSERT(NPTNEW.LE.NPM,
     '          '>>Too many nodes - increase NPM',ERROR,*9999)
              NPNODE(0,nr)=NPNODE(0,nr)+1
              NPNODE(0,0)=NPNODE(0,0)+1
              NPNODE(NPNODE(0,nr),nr)=NPTNEW
              NPNODE(NPNODE(0,0),0)=NPTNEW
              NODE=NPTNEW
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(NODE,.FALSE.,ERROR,*9999)
              NPT(nr)=NPTNEW
              NPT( 0)=NPTNEW

C ***         Calculate nodal coords & derivs wrt Xi-directions
              DO nj=1,NJ_LOC(0,0,nr)
                nb=NBJ(nj,ne)
!AJP 12/7/96 Check on nb added becuase if a fibre field has been
!defined in another region nj_loc(0,0) may exceed what is needed for 
!the current region and nb may not be defined.
                IF(nb.GT.0) THEN 
                  DO nu=1,NUT(nb)
                    X(nu,nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                INP(1,1,nb),nb,nu,XI,XE(1,nj))
                  ENDDO
                ELSE
                  X(1,nj)=0.0d0 !since this is used later
                ENDIF !nb
              ENDDO

C ***         If an angle is 2*pi set it to zero
              IF(ITYP10(1).GE.2) THEN
                IF(DABS(X(1,2)-2.0D0*PI).LT.1.d-5) X(1,2)=0.0D0
                IF(ITYP10(1).GE.3) THEN
                  IF(DABS(X(1,3)-2.0D0*PI).LT.1.d-5) X(1,3)=0.0D0
                ENDIF
              ENDIF

              IF(DOP) THEN
                DO nj=1,NJT
                  WRITE(OP_STRING,'('' EXIST='',L1,'' NODE='',I4)')
     '              EXIST,NODE
      		  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c                  WRITE(OP_STRING,'('' X(nu,'',I1,''): '',10E12.3)')
c     '              nj,(X(nu,nj),nu=1,NUT(NBJ(nj,ne)))
c      		  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(*,'('' X(1,'',I1,''): '',10E12.3)')
     '              nj,X(1,nj)
                ENDDO
              ENDIF

              NPTOLD=NPNODE(NPNODE(0,nr)-1,nr)
              NJP(NODE)=NJP(NPTOLD)
              DO nj=1,NJ_LOC(0,0,nr)
                NKJ(nj,NODE)=NKJ(nj,NPTOLD)
              ENDDO

              IF(NBI(nb1).EQ.3) THEN !global line derivatives
C ***           Calculate derivatives of s* wrt Xi
C ***           Note: Only correct for bicubic at present
                IF(IDRN.EQ.1) THEN
                  IF(i3.EQ.1) THEN
                    IF(i2.EQ.1) THEN
                      DSDXI(1)=0.25D0*(SE( 2,nb1,ne)+SE( 6,nb1,ne))
                      DSDXI(2)=0.50D0*(SE( 3,nb1,ne)+SE( 7,nb1,ne))
                    ELSE IF(i2.EQ.2) THEN
                      DSDXI(1)=0.25D0*(SE(10,nb1,ne)+SE(14,nb1,ne))
                      DSDXI(2)=0.50D0*(SE(11,nb1,ne)+SE(15,nb1,ne))
                    ENDIF
                  ELSE IF(i3.EQ.2) THEN
                    IF(i2.EQ.1) THEN
                      DSDXI(1)=0.25D0*(SE(18,nb1,ne)+SE(22,nb1,ne))
                      DSDXI(2)=0.50D0*(SE(19,nb1,ne)+SE(23,nb1,ne))
                    ELSE IF(i2.EQ.2) THEN
                      DSDXI(1)=0.25D0*(SE(26,nb1,ne)+SE(30,nb1,ne))
                      DSDXI(2)=0.50D0*(SE(27,nb1,ne)+SE(31,nb1,ne))
                    ENDIF
                  ENDIF
                ELSE IF(IDRN.EQ.2) THEN
                  IF(i2.EQ.1) THEN
                    IF(i3.EQ.1) THEN
                      DSDXI(1)=0.50D0*(SE( 2,nb1,ne)+SE(10,nb1,ne))
                      DSDXI(2)=0.25D0*(SE( 3,nb1,ne)+SE(11,nb1,ne))
                    ELSE IF(i3.EQ.2) THEN
                      DSDXI(1)=0.50D0*(SE( 6,nb1,ne)+SE(14,nb1,ne))
                      DSDXI(2)=0.25D0*(SE( 7,nb1,ne)+SE(15,nb1,ne))
                    ENDIF
                  ELSE IF(i2.EQ.2) THEN
                    IF(i3.EQ.1) THEN
                      DSDXI(1)=0.50D0*(SE(18,nb1,ne)+SE(26,nb1,ne))
                      DSDXI(2)=0.25D0*(SE(19,nb1,ne)+SE(27,nb1,ne))
                    ELSE IF(i3.EQ.2) THEN
                      DSDXI(1)=0.50D0*(SE(22,nb1,ne)+SE(30,nb1,ne))
                      DSDXI(2)=0.25D0*(SE(23,nb1,ne)+SE(31,nb1,ne))
                    ENDIF
                  ENDIF
                ENDIF !idrn
                IF(DOP) THEN
      		  WRITE(OP_STRING,'('' DSDXI(ni): '',3E12.3)')
     '		    (DSDXI(ni),ni=1,NIT(nb1))
      		  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF

C cpb 5/12/96 re-adding arc-length scaling
c cpb 8/10/94 changing around nbi(nb) = 4 and 5
C             ELSE IF(NBI(nb1).EQ.4) THEN !arc-length derivatives
              ELSE IF(NBI(nb1).EQ.5.OR.NBI(nb1).EQ.6) 
     '            THEN !arc-length or ave. arc-length derivatives
C ***           Calculate derivatives of arc-length wrt Xi
                DO ni=1,NIT(nb1)
                  nu=1+ni*(1+ni)/2
                  SUM=0.0D0
                  IF(ITYP10(1).EQ.1) THEN
                    DO nj=1,NJT
                      SUM=SUM+X(nu,nj)**2
                    ENDDO
                  ELSE IF(ITYP10(1).EQ.2) THEN
                    R=X(1,1)
                    RR=R*R
                    SUM=SUM+X(nu,1)**2+RR*X(nu,2)**2
                    IF(NJT.GT.2) SUM=SUM+X(nu,3)**2
                  ELSE IF(ITYP10(1).EQ.3) THEN
                    R=X(1,1)
                    RR=R*R
                    RC=R*DCOS(X(1,3))
                    RRC=RC*RC
                    SUM=SUM+X(nu,1)**2+RRC*X(nu,2)**2+RR*X(nu,3)**2
                  ELSE IF(ITYP10(1).EQ.4) THEN
                    AA=FOCUS*FOCUS
                    SLX=DSINH(X(1,1))
                    SMX=DSIN(X(1,2))
                    G1=AA*(SLX*SLX+SMX*SMX)
                    G3=AA* SLX*SLX*SMX*SMX
                    SUM=SUM+G1*(X(nu,1)**2+X(nu,2)**2)
                    IF(NJT.GT.2) SUM=SUM+G3*X(nu,3)**2
                  ENDIF
                  DSDXI(ni)=DSQRT(SUM)
                ENDDO !ni
                IF(DOP) THEN
      		  WRITE(OP_STRING,'('' DSDXI(ni): '',3E12.3)')
     '		    (DSDXI(ni),ni=1,NIT(nb1))
      		  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF !dop
              ENDIF !nbi

C cpb 5/12/96 re-adding arc-length scaling
c cpb 8/10/94 changing around nbi(nb) = 4 and 5
C             IF(NBI(nb1).EQ.3.OR.NBI(nb1).EQ.4) THEN
              IF(NBI(nb1).EQ.3.OR.NBI(nb1).EQ.5.OR.NBI(nb1).EQ.6) THEN
C ***           Store nodal derivatives wrt arc-length
                DO nj=1,NJ_LOC(0,0,nr)
                  XP(1,1,nj,NODE)=X(1,nj)
                  nb=NBJ(nj,ne)
                  DO nk=2,NKT(0,nb)
                    IF(nk.EQ.2) THEN
                      IF(DABS(DSDXI(1)).GT.1.d-6) THEN
                        XP(2,1,nj,NODE)=X(2,nj)/DSDXI(1)
                      ELSE
                        WRITE(OP_STRING,
     '                    '('' Warning: Arc length deriv is zero'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ENDIF
                    ELSE IF(nk.EQ.3) THEN
                      IF(DABS(DSDXI(2)).GT.1.d-6) THEN
                        XP(3,1,nj,NODE)=X(4,nj)/DSDXI(2)
                      ELSE
                        WRITE(OP_STRING,
     '                    '('' Warning: Arc length deriv is zero'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ENDIF
                    ELSE IF(nk.EQ.4) THEN
                      IF(NOCROSS) THEN  !New AJP 14-6-93
                        XP(4,1,nj,NODE)=0.0D0
                      ELSE
                        IF(DABS(DSDXI(1)).GT.1.d-6.AND.DABS(DSDXI(2))
     '                    .GT.1.d-6) THEN
                          XP(4,1,nj,NODE)=X(6,nj)/(DSDXI(1)*DSDXI(2))
                        ELSE
                          WRITE(OP_STRING,
     '                      '('' Warning: Arc length deriv is zero'')')
      		          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ENDIF
                      ENDIF
                    ELSE IF(nk.EQ.5) THEN
                      IF(DABS(DSDXI(3)).GT.1.d-6) THEN
                        XP(5,1,nj,NODE)=X(7,nj)/DSDXI(3)
                      ELSE
                        WRITE(OP_STRING,
     '                    '('' Warning: Arc length deriv is zero'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ENDIF
                    ELSE IF(nk.EQ.6) THEN
                      IF(DABS(DSDXI(2)).GT.1.d-6.AND.DABS(DSDXI(3))
     '                  .GT.1.d-6) THEN
                        XP(6,1,nj,NODE)=X(9,nj)/(DSDXI(2)*DSDXI(3))
                      ELSE
                        WRITE(OP_STRING,
     '                    '('' Warning: Arc length deriv is zero'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ENDIF
                    ELSE IF(nk.EQ.7) THEN
                      IF(DABS(DSDXI(3)).GT.1.d-6.AND.DABS(DSDXI(1))
     '                  .GT.1.d-6) THEN
                        XP(7,1,nj,NODE)=X(10,nj)/(DSDXI(3)*DSDXI(1))
                      ELSE
                        WRITE(OP_STRING,
     '                    '('' Warning: Arc length deriv is zero'')')
      			CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ENDIF
                    ENDIF
                  ENDDO !nk
                ENDDO !nj

              ELSE
                DO nj=1,NJ_LOC(0,0,nr)
                  XP(1,1,nj,NODE)=X(1,nj)
                ENDDO
              ENDIF !nbi

C Check for multiple versions of fibre/sheet angles in nodes at either 
C end of the line on the original element which contains the new NODE.

              IF(DOP) THEN
                WRITE(OP_STRING,'(/'' IDRN='',I1,'
     '            //''' i1='',I1,'' i2='',I1,'' i3='',I1)') 
     '            IDRN,i1,i2,i3
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              nb=NBJ(1,ne)
              np1=NPNE(NNIP(1,i2,i3),nb,ne)
              np2=NPNE(NNIP(2,i2,i3),nb,ne)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' NODE='',I6,'' has adj element'
     '            //' vertices '',2I2,'' --> global nodes '',2I6)')
     '            NODE,NNIP(1,i2,i3),NNIP(2,i2,i3),np1,np2
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' #versions for np1 (nj=1..): '','
     '            //'10I3)') (NVJP(nj,np1),nj=1,NJ_LOC(0,0,nr))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' #versions for np2 (nj=1..): '','
     '            //'10I3)') (NVJP(nj,np2),nj=1,NJ_LOC(0,0,nr))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF !dop

              DO njj=1,NJ_LOC(NJL_GEOM,0,nr) !loop over geometric coords
                nj=NJ_LOC(NJL_GEOM,njj,nr)
                NVJP(nj,NODE)=NVJP(nj,np1)
              ENDDO !nj
              
C Deal with special case of node on axis for prolate coords
              IF(ITYP10(nr).EQ.4) THEN !prolate coords
                IF(DABS(XP(1,1,2,NODE)).LT.1.0d-6) THEN !on axis
                  nj=3 !for theta coordinate
                  IF(NVJP(nj,NODE).GT.1) THEN !multiple versions
                    write(*,'('' New node on axis with mult thetas'')')
                    DO nv=1,NVJP(nj,NODE)
                      DO nk=1,NKJ(nj,NODE)                  
                        XP(nk,nv,nj,NODE)=XII *XP(nk,nv,nj,np1)
     '                    +(1.d0-XII)*XP(nk,nv,nj,np2)
                      ENDDO !nk
                    ENDDO !nv
                  ENDIF !nvjp>1
                ENDIF !mu=0 (on axis)
              ENDIF !ityp10=4 (prolate)
        
              DO njj=1,NJ_LOC(NJL_FIEL,0,nr) !loop over field coords
                nj=NJ_LOC(NJL_FIEL,njj,nr)
                NVJP(nj,NODE)=NVJP(nj,np1)
              ENDDO !nj

              DO njj=1,NJ_LOC(NJL_FIBR,0,nr) !loop over fibre coords
                nj=NJ_LOC(NJL_FIBR,njj,nr)
                NVJP(nj,NODE)=MAX(NVJP(nj,np1),NVJP(nj,np2))
                IF(NVJP(nj,NODE).GT.1) THEN !multiple versions
                  DO nv=1,NVJP(nj,NODE)
                    DO nk=1,NKJ(nj,NODE)
                      IF(DOP) THEN
                        WRITE(OP_STRING,'('' nj='',I1,'' nv='',I1,'
     '                    //''' nk='',I1)') nj,nv,nk
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF
                      IF(NVJP(nj,np1).EQ.NVJP(nj,np2)) THEN
                        XP(nk,nv,nj,NODE)=XII *XP(nk,nv,nj,np1)
     '                             +(1.d0-XII)*XP(nk,nv,nj,np2)
                      ELSE IF(NVJP(nj,np1).EQ.1) THEN
                        XP(nk,nv,nj,NODE)=XII *XP(nk, 1,nj,np1)
     '                             +(1.d0-XII)*XP(nk,nv,nj,np2)
                      ELSE IF(NVJP(nj,np2).EQ.1) THEN
                        XP(nk,nv,nj,NODE)=XII *XP(nk,nv,nj,np1)
     '                             +(1.d0-XII)*XP(nk, 1,nj,np2)
                      ELSE
                        ERROR='Incompatible versions'
                        GO TO 9999  
                      ENDIF !nvjp
                      IF(DOP) THEN
                        WRITE(OP_STRING,'('' XP(nk,nv,nj,NODE)='','
     '                    //'E12.3)') XP(nk,nv,nj,NODE)
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF
                    ENDDO !nk
                  ENDDO !nv

                ENDIF !nvjp>1
              ENDDO !nj

            ELSE !Node exits, but possibly only in another region
              FOUND=.FALSE.
              nonode=1
              DO WHILE((nonode.LE.NPNODE(0,nr)).AND.(.NOT.FOUND))
                NP_TEST=NPNODE(nonode,nr)
                IF(np_test.eq.np)THEN !Node np is in current region
                  FOUND=.TRUE.
                ELSE
                  nonode=nonode+1
                ENDIF
              ENDDO
              IF(.NOT.FOUND)THEN !Update arrays to include np
                NPNODE(0,nr)=NPNODE(0,nr)+1
                NPNODE(NPNODE(0,nr),nr)=np
                IF(np.GT.NPT(nr))THEN
                  NPT(nr)=np
                ENDIF
! New AJP 21-1-94
                !Node exists but in another region.
                !Update current region arrays.
                !Firstly find other region number of node np
                FOUND=.FALSE.
                nr2=1
                DO WHILE ((nr2.LE.NRT).AND.(.NOT.FOUND))
                  IF(nr2.ne.nr)THEN
                    nonode=1
                    DO WHILE ((nonode.LE.NPNODE(0,nr2)).AND.
     '                        (.NOT.FOUND))
                      NP_TEST=NPNODE(nonode,nr2)
                      IF(np.EQ.NP_TEST)THEN
                        FOUND=.TRUE.
                      ELSE
                        nonode=nonode+1
                      ENDIF
                    ENDDO
                    IF(.NOT.FOUND)nr2=nr2+1
                  ELSE
                    nr2=nr2+1
                  ENDIF
                ENDDO
                IF(nr2.GT.NRT)THEN
                  WRITE(OP_STRING,
     '              '('' Cant find node in other regions ?'')')
        	  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GOTO 9999
                ENDIF
              ENDIF !.not.found
            ENDIF !.not.exist
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(NODE,.FALSE.,ERROR,*9999)
C*** Calculate the element arrays for the node
            IF(IBT(1,IDRN,nb).EQ.1.AND.IBT(2,IDRN,nb).EQ.2) 
     '        THEN !Quadratic Lagrange
              IF(XII.EQ.0.5d0) THEN
                IF(i1.EQ.1) THEN
                  NPNE1(NNIP(2,i2,i3))=NODE
                ELSE
                  NPNE2(NNIP(2,i2,i3))=NODE
                ENDIF
              ELSE IF(XII.GT.0.5d0) THEN
                NPNE2(NNIP(i1,i2,i3))=NODE
                IF(i1.EQ.1) NPNE1(NNIP(3,i2,i3))=NODE
              ELSE
                NPNE1(NNIP(i1+1,i2,i3))=NODE
                IF(i1.EQ.2) NPNE2(NNIP(1,i2,i3))=NODE
              ENDIF
            ELSE IF(IBT(1,IDRN,nb).EQ.1.AND.IBT(2,IDRN,nb).EQ.3) 
     '          THEN !Cubic Lagrange
              IF(XII.EQ.(1.0d0/3.0d0)) THEN
                IF(i1.EQ.3) THEN
                  NPNE2(NNIP(2,i2,i3))=NODE
                ELSE
                  NPNE1(NNIP(i1+1,i2,i3))=NODE
                ENDIF
              ELSE IF(XII.EQ.(2.0d0/3.0d0)) THEN
                IF(i1.EQ.1) THEN
                  NPNE1(NNIP(3,i2,i3))=NODE
                ELSE
                  NPNE2(NNIP(i1,i2,i3))=NODE
                ENDIF
              ELSE IF(XII.GT.(2.0d0/3.0d0)) THEN
                NPNE2(NNIP(i1,i2,i3))=NODE
                IF(i1.EQ.1) NPNE1(NNIP(4,i2,i3))=NODE
              ELSE IF(XII.GT.(1.0d0/3.0d0)) THEN
                IF(i1.LE.2) NPNE1(NNIP(i1+2,i2,i3))=NODE 
                IF(i1.GE.2) NPNE2(NNIP(i1-1,i2,i3))=NODE 
              ELSE
                NPNE1(NNIP(i1+1,i2,i3))=NODE
                IF(i1.EQ.3) NPNE2(NNIP(1,i2,i3))=NODE
              ENDIF
            ELSE
              NPNE1(NNIP(2,i2,i3))=NODE
              NPNE2(NNIP(1,i2,i3))=NODE
            ENDIF

            IF(DOP) THEN
              WRITE(OP_STRING,'('' NODE='',I6,'' exists'')') NODE
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF !dop

C Update version #s for elements ne & ne_NEW (PJH 4Jul96)
            DO nj=1,NJ_LOC(0,0,nr)
              nb=NBJ(nj,ne)
              nn_list(0)=0
!AJP 12/7/96 Check on nb added becuase if a fibre field has been
!defined in another region nj_loc(0,0) may exceed what is needed for 
!the current region and nb may not be defined.
              IF(nb.GT.0) THEN 
                DO nn=1,NNT(nb) !store vertices coincident with NODE
                  IF(NPNE(nn,nb,ne).EQ.NODE) THEN !vertex nn is at NODE
                    nn_list(0)=nn_list(0)+1
                    nn_list(nn_list(0))=nn
                  ENDIF 
                ENDDO !nn
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' nj='',I1,'' Element vertices'
     '              //' coinciding with node'',I6,'' are '',5I3)') 
     '              nj,NODE,(nn_list(n),n=1,nn_list(0))
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF !dop
              ENDIF
              IF(nn_list(0).GT.0.AND.NVJP(nj,NODE).GT.1) THEN !multiple
C                                                             !versions at NODE 
                nv1=NVJE(nn_list(1),nb,nj,ne) !1st vertex
                nv2=NVJE(nn_list(2),nb,nj,ne) !2nd vertex
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' current multiple versions:'
     '              //' nv1='',I2,'' nv2='',I2)') nv1,nv2
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF !dop
                NVJP(nj,NODE)=NVJP(nj,NODE)+1 !new version# total
                nv=NVJP(nj,NODE)              !new version#
                CALL ASSERT(nv.LE.NVM,
     '            '>>Too many versions -Increase NVMX',ERROR,*9999)
                NVJE(nn_list(2),nb,nj,ne)=nv  !new version for vertex 2
                NVJE(nn_list(1),nb,nj,ne_NEW)=nv !///y for v1 at ne_NEW
                DO nk=1,NKJ(nj,NODE)
C PJH 4Jul96 Note that this doesn't handle nkt>1 case yet
                  XP(nk,nv,nj,NODE)= XII *XE(nn_list(1),nj)
     '                        +(1.d0-XII)*XE(nn_list(2),nj)
                ENDDO !nk
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' Update vertices '',I2,'
     '              //''' & '',I2,'' of elements '',I5,'' &'''
     '              //',I5,'',respec, to version'',I2)')
     '              nn_list(2),nn_list(1),ne,ne_NEW,nv
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' Value(nk=1) is '',E12.3)') 
     '              XP(1,nv,nj,NODE)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF !dop
              ENDIF !nvjp
            ENDDO !nj

          ENDDO !i1

C*** Fill in the old nodes into the old and new elements
          IF(IBT(1,IDRN,nb).EQ.1.AND.IBT(2,IDRN,nb).EQ.2) 
     '      THEN !Quadratic Lagrange
            IF(XII.EQ.0.5d0) THEN
              NPNE1(NNIP(1,i2,i3))=NPNE(NNIP(1,i2,i3),nb,ne)
              NPNE1(NNIP(3,i2,i3))=NPNE(NNIP(2,i2,i3),nb,ne)
              NPNE2(NNIP(1,i2,i3))=NPNE(NNIP(2,i2,i3),nb,ne)
              NPNE2(NNIP(3,i2,i3))=NPNE(NNIP(3,i2,i3),nb,ne)
            ELSE IF(XII.GT.0.5d0) THEN
              NPNE1(NNIP(1,i2,i3))=NPNE(NNIP(1,i2,i3),nb,ne)
              NPNE1(NNIP(2,i2,i3))=NPNE(NNIP(2,i2,i3),nb,ne)
              NPNE2(NNIP(3,i2,i3))=NPNE(NNIP(3,i2,i3),nb,ne)
            ELSE
              NPNE1(NNIP(1,i2,i3))=NPNE(NNIP(1,i2,i3),nb,ne)
              NPNE2(NNIP(2,i2,i3))=NPNE(NNIP(2,i2,i3),nb,ne)
              NPNE2(NNIP(3,i2,i3))=NPNE(NNIP(3,i2,i3),nb,ne)
            ENDIF                  
          ELSE IF((IBT(1,IDRN,nb).EQ.1.OR.
     '        IBT(1,IDRN,nb).EQ.5.OR.IBT(1,IDRN,nb).EQ.6).AND.
     '        IBT(2,IDRN,nb).EQ.3) THEN !Cubic Lagrange
            IF(XII.EQ.(1.0d0/3.0d0)) THEN
              NPNE1(NNIP(1,i2,i3))=NPNE(NNIP(1,i2,i3),nb,ne)
              NPNE1(NNIP(4,i2,i3))=NPNE(NNIP(2,i2,i3),nb,ne)
              NPNE2(NNIP(1,i2,i3))=NPNE(NNIP(2,i2,i3),nb,ne)
              NPNE2(NNIP(3,i2,i3))=NPNE(NNIP(3,i2,i3),nb,ne)
              NPNE2(NNIP(4,i2,i3))=NPNE(NNIP(4,i2,i3),nb,ne)
            ELSE IF(XII.EQ.(2.0d0/3.0d0)) THEN
              NPNE1(NNIP(1,i2,i3))=NPNE(NNIP(1,i2,i3),nb,ne)
              NPNE1(NNIP(2,i2,i3))=NPNE(NNIP(2,i2,i3),nb,ne)
              NPNE1(NNIP(4,i2,i3))=NPNE(NNIP(3,i2,i3),nb,ne)
              NPNE2(NNIP(1,i2,i3))=NPNE(NNIP(3,i2,i3),nb,ne)
              NPNE2(NNIP(4,i2,i3))=NPNE(NNIP(4,i2,i3),nb,ne)
            ELSE IF(XII.GT.(2.0d0/3.0d0)) THEN
              NPNE1(NNIP(1,i2,i3))=NPNE(NNIP(1,i2,i3),nb,ne)
              NPNE1(NNIP(2,i2,i3))=NPNE(NNIP(2,i2,i3),nb,ne)
              NPNE1(NNIP(3,i2,i3))=NPNE(NNIP(3,i2,i3),nb,ne)
              NPNE2(NNIP(4,i2,i3))=NPNE(NNIP(4,i2,i3),nb,ne)
            ELSE IF(XII.GT.(1.0d0/3.0d0)) THEN
              NPNE1(NNIP(1,i2,i3))=NPNE(NNIP(1,i2,i3),nb,ne)
              NPNE1(NNIP(2,i2,i3))=NPNE(NNIP(2,i2,i3),nb,ne)
              NPNE2(NNIP(3,i2,i3))=NPNE(NNIP(3,i2,i3),nb,ne)
              NPNE2(NNIP(4,i2,i3))=NPNE(NNIP(4,i2,i3),nb,ne)
            ELSE
              NPNE1(NNIP(1,i2,i3))=NPNE(NNIP(1,i2,i3),nb,ne)
              NPNE2(NNIP(2,i2,i3))=NPNE(NNIP(2,i2,i3),nb,ne)
              NPNE2(NNIP(3,i2,i3))=NPNE(NNIP(3,i2,i3),nb,ne)
              NPNE2(NNIP(4,i2,i3))=NPNE(NNIP(4,i2,i3),nb,ne)
            ENDIF
          ELSE
            NPNE1(NNIP(1,i2,i3))=NPNE(NNIP(1,i2,i3),nb,ne)
            NPNE2(NNIP(2,i2,i3))=NPNE(NNIP(2,i2,i3),nb,ne)
          ENDIF

        ENDDO !i2
      ENDDO !i3

      IF(DOP) THEN
        WRITE(OP_STRING,'('' NPNE1: '',27I3)') (NPNE1(nn),nn=1,NNTB)
      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NPNE2: '',27I3)') (NPNE2(nn),nn=1,NNTB)
      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      DO nb=1,NBFT
        DO nn=1,NNT(nb)
          NPNE(nn,nb,ne)=NPNE1(nn)
        ENDDO
      ENDDO !nb

      IF(IDRN.LT.3.AND.NKT(0,nb1).GT.1) THEN
C ***   Update scaling factors in current element & new element

        DO nb=1,NBFT
          IF(NKT(0,nb).GT.1) THEN
            IF(IBT(1,IDRN,nb).EQ.2) THEN
              IF(NKT(0,nb).LE.2) THEN
                DO nn=1,NNT(nb)
                  ns=(nn-1)*NKT(0,nb)+2
                  SE(ns,nb,ne)=SE(ns,nb,ne)/2.d0
                  IF(DOP) THEN
      		    WRITE(IO4,'('' ne='',I4,'' nb='',i1,'' ns='''
     '		      //',i2,'//' '' SE='',E12.3)') 
     '                ne,nb,ns,SE(ns,nb,ne)
      		    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDDO
              ELSE IF(NKT(0,nb).EQ.4) THEN
                IF(NIT(nb).EQ.2) THEN
                  NTOT=1
                ELSE IF(NIT(nb).EQ.3) THEN
                  NTOT=2
                ENDIF
                DO n=1,NTOT
                  k=16*(n-1)
                  IF(IDRN.EQ.1) THEN
C ***               Xi(1) derivs
                    SE( 2+k,nb,ne_NEW)=0.25D0*(SE( 2+k,nb,ne)
     '                +SE( 6+K,nb,ne))
                    SE( 6+k,nb,ne_NEW)=0.50D0* SE( 6+k,nb,ne)
                    SE(10+k,nb,ne_NEW)=0.25D0*(SE(10+k,nb,ne)
     '                +SE(14+k,nb,ne))
                    SE(14+k,nb,ne_NEW)=0.50D0* SE(14+k,nb,ne)
                    SE( 6+k,nb,ne)    =0.25D0*(SE( 2+k,nb,ne)
     '                +SE( 6+k,nb,ne))
                    SE( 2+k,nb,ne)    =0.50D0* SE( 2+k,nb,ne)
                    SE(14+k,nb,ne)    =0.25D0*(SE(10+k,nb,ne)
     '                +SE(14+k,nb,ne))
                    SE(10+k,nb,ne)    =0.50D0* SE(10+k,nb,ne)
C ***               Xi(2) derivs
                    SE( 3+k,nb,ne_NEW)=0.50D0*(SE( 3+k,nb,ne)
     '                +SE( 7+k,nb,ne))
                    SE( 7+k,nb,ne_NEW)=      SE( 7+k,nb,ne)
                    SE(11+k,nb,ne_NEW)=0.50D0*(SE(11+k,nb,ne)
     '                +SE(15+k,nb,ne))
                    SE(15+k,nb,ne_NEW)=      SE(15+k,nb,ne)
                    SE( 7+k,nb,ne)    =0.50D0*(SE( 3+k,nb,ne)
     '                +SE( 7+k,nb,ne))
                    SE( 3+k,nb,ne)    =      SE( 3+k,nb,ne)
                    SE(15+k,nb,ne)    =0.50D0*(SE(11+k,nb,ne)
     '                +SE(15+k,nb,ne))
                    SE(11+k,nb,ne)    =      SE(11+k,nb,ne)
                  ELSE IF(IDRN.EQ.2) THEN
C ***               Xi(1) derivs
                    SE( 2+k,nb,ne_NEW)=0.50D0*(SE( 2+k,nb,ne)
     '                +SE(10+k,nb,ne))
                    SE( 6+k,nb,ne_NEW)=0.50D0*(SE( 6+k,nb,ne)
     '                +SE(14+k,nb,ne))
                    SE(10+k,nb,ne_NEW)=      SE(10+k,nb,ne)
                    SE(14+k,nb,ne_NEW)=      SE(14+k,nb,ne)
                    SE(10+k,nb,ne)    =0.50D0*(SE( 2+k,nb,ne)
     '                +SE(10+k,nb,ne))
                    SE(14+k,nb,ne)    =0.50D0*(SE( 6+k,nb,ne)
     '                +SE(14+k,nb,ne))
                    SE( 2+k,nb,ne)    =      SE( 2+k,nb,ne)
                    SE( 6+k,nb,ne)    =      SE( 6+k,nb,ne)
C ***               Xi(2) derivs
                    SE( 3+k,nb,ne_NEW)=0.25D0*(SE( 3+k,nb,ne)
     '                +SE(11+k,nb,ne))
                    SE( 7+k,nb,ne_NEW)=0.25D0*(SE( 7+k,nb,ne)
     '                +SE(15+k,nb,ne))
                    SE(11+k,nb,ne_NEW)=0.50D0* SE(11+k,nb,ne)
                    SE(15+k,nb,ne_NEW)=0.50D0* SE(15+k,nb,ne)
                    SE(11+k,nb,ne)    =0.25D0*(SE( 3+k,nb,ne)
     '                +SE(11+k,nb,ne))
                    SE(15+k,nb,ne)    =0.25D0*(SE( 7+k,nb,ne)
     '                +SE(15+k,nb,ne))
                    SE( 3+k,nb,ne)    =0.50D0* SE( 3+k,nb,ne)
                    SE( 7+k,nb,ne)    =0.50D0* SE( 7+k,nb,ne)
                  ELSE IF(IDRN.EQ.3) THEN
                    write(OP_STRING,*) '!!!!! additions needed'
      		    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
C ***             Xi(1),Xi(2) derivs
                  SE( 4+k,nb,ne_NEW)=SE( 2+k,nb,ne_NEW)
     '                              *SE( 3+k,nb,ne_NEW)
                  SE( 8+k,nb,ne_NEW)=SE( 6+k,nb,ne_NEW)
     '                              *SE( 7+k,nb,ne_NEW)
                  SE(12+k,nb,ne_NEW)=SE(10+k,nb,ne_NEW)
     '                              *SE(11+k,nb,ne_NEW)
                  SE(16+k,nb,ne_NEW)=SE(14+k,nb,ne_NEW)
     '                              *SE(15+k,nb,ne_NEW)
                  SE( 4+k,nb,ne)  =SE( 2+k,nb,ne)  *SE( 3+k,nb,ne)
                  SE( 8+k,nb,ne)  =SE( 6+k,nb,ne)  *SE( 7+k,nb,ne)
                  SE(12+k,nb,ne)  =SE(10+k,nb,ne)  *SE(11+k,nb,ne)
                  SE(16+k,nb,ne)  =SE(14+k,nb,ne)  *SE(15+k,nb,ne)
                ENDDO
              ELSE
                WRITE(OP_STRING,'('' SEs not calculated'')')
      		CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ns=0
                DO nn=1,NNT(nb)
                  DO nk=1,NKT(nn,nb)
C OLD AJP 25-5-93    ns=nk+(nn-1)*NKT(0,nb)
                    ns=ns+1
                    SE(ns,nb,ne_NEW)=SE(ns,nb,ne)
                  ENDDO
                ENDDO
              ENDIF !nkt
              IF(DOP) THEN
                DO ns=1,NST(nb)
                  WRITE(OP_STRING,'('' ne='',I4,'' nb='',i1,'' ns='','
     '		    //'I2,'//' '' SE='',E12.3)') ne,nb,ns,SE(ns,nb,ne)
      		  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDDO
              ENDIF !dop
            ENDIF !ibt
          ENDIF !nkt
        ENDDO !nb
      ENDIF !idrn

C *** Form other arrays for new element ne_NEW
      NJE(ne_NEW)=NJE(ne)
      NHE(ne_NEW)=NHE(ne)
      NW(ne_NEW,1)=NW(ne,1)
      NCO(ne_NEW)=NCO(ne)
      NRE(ne_NEW)=NRE(ne)
      DO nm=1,NMM
        CE(nm,ne_NEW)=CE(nm,ne)
      ENDDO !nm
      DO nj=1,NJ_LOC(0,0,nr)
        NBJ(nj,ne_NEW)=NBJ(nj,ne)
      ENDDO !nj

      IF(CALL_EQUA.OR.CALL_FIT.OR.CALL_OPTI) THEN
        DO nhx=1,NHE(ne)
          nh=NH_LOC(nhx,nx)
          DO nc=1,NCM
            NBH(nh,nc,ne_NEW)=NBH(nh,nc,ne)
          ENDDO !nc
        ENDDO !nhx
      ENDIF !call_equa etc

      DO nb=1,NBFT
        DO nn=1,NNT(nb)
          NPNE(nn,nb,ne_NEW)=NPNE2(nn)
        ENDDO !nn
      ENDDO !nb

      DO nb=1,NBFT
        DO nn=1,NNT(nb)
          DO nk=1,NKT(nn,nb)
            NKE(nk,nn,nb,ne_NEW)=NKE(nk,nn,nb,ne)
          ENDDO !nk
        ENDDO !nn
      ENDDO !nb

      NPSTART=NPT(0)+1
      NET(0)=ne_NEW   !is new highest element # in mesh
      NET(nr)=ne_NEW  !is new highest element # in region nr
      CALL ASSERT(NET(nr).LE.NEM,
     '  '>>Too many elements - increase NEM',ERROR,*9999)
      NEELEM(0,nr)=NEELEM(0,nr)+1 !is new #elements in region nr
      NEELEM(NEELEM(0,nr),nr)=ne_NEW !is new element number

!     !Update element list in each group
      CALL UPGREL(ne,ne_NEW,ERROR,*9999)

! Check that the global nodes of new element ne_NEW are in the 
! node list for the current region
      nb=NBJ(1,ne_NEW)
      DO nn=1,NNT(nb)
        np=NPNE(nn,nb,ne_NEW)
        EXIST=.FALSE.
        DO nonode=1,NPNODE(0,nr)
          IF(NPNODE(nonode,nr).eq.np) EXIST=.TRUE.
        ENDDO
        IF(.NOT.EXIST) THEN
          NPNODE(0,nr)=NPNODE(0,nr)+1
          NPNODE(NPNODE(0,nr),nr)=np
        ENDIF
      ENDDO

C *** Diagnostic output
      IF(DOP) THEN
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nb=1,NBFT
            WRITE(OP_STRING,'('' NPNE(nn,'',i1,'','',I4,'
     '        //'''): '',20I4)') nb,ne,(NPNE(nn,nb,ne),nn=1,NNT(nb))
      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
        DO np=1,NPT(nr)
          DO nj=1,NJP(np)
            DO nv=1,NVJP(nj,np)
              WRITE(OP_STRING,'('' XP(nk=..,nv='',I2,'',nj='',I1,'
     '          //''',np='',I4,''): '',10E11.4)')
     '          nv,nj,np,(XP(nk,nv,nj,np),nk=1,NKJ(nj,np))
      	      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !nv
          ENDDO !nj
        ENDDO !np
      ENDIF !dop

      CALL EXITS('DIVH1')
      RETURN
 9999 CALL ERRORS('DIVH1',ERROR)
      CALL EXITS('DIVH1')
      RETURN 1
      END

C 24/2/97 LC removed from routine :

      SUBROUTINE BASIS1(IBT,IDO,INP,nb,NGAP,PG,WG,XIG,ERROR,*)

C#### Subroutine: BASIS1
C###  Description: 
C###    BASIS1 inputs parameters for Lagrange/Hermite tensor 
C###    product basis functions and Gauss-Legendre or Gauss-Lobatto 
C###    quadrature.


C cpb 13/7/95 Adding quadratic Hermite interpolation and this now needs
C to be calculated after INP
C      IF(IOTYPE.NE.3) THEN
C        nk=0
C        IF(NIT(nb).EQ.1) THEN
C          DO ni1=1,IBT(1,1)
C            nk=nk+1
C            IDO(nk,1)=ni1
C          ENDDO
C        ELSE IF(NIT(nb).EQ.2) THEN
C          DO ni2=1,IBT(1,2)
C            DO ni1=1,IBT(1,1)
C              nk=nk+1
C              IDO(nk,1)=ni1
C              IDO(nk,2)=ni2
C            ENDDO
C          ENDDO
C        ELSE IF(NIT(nb).EQ.3) THEN
C          DO ni3=1,IBT(1,3)
C            DO ni2=1,IBT(1,2)
C              DO ni1=1,IBT(1,1)
C                IF(nk.LT.7) THEN
C                  nk=nk+1
C                  IDO(nk,1)=ni1
C                  IDO(nk,2)=ni2
C                  IDO(nk,3)=ni3
C                ENDIF
C              ENDDO
C            ENDDO
C          ENDDO
C        ENDIF
C        NKT(0,nb)=nk
C        NST(nb)=0
C        DO nn=1,NNT(nb) !AJP 25-5-93  All nkt(nn,nb) the same for now
C          NKT(nn,nb)=NKT(0,nb)
C          NST(nb)=NST(nb)+NKT(nn,nb)
C        ENDDO !nn
C      ENDIF

c cpb 13/7/95 Need to calc NKT (and hence NST) after INP
C      IF(NST(nb).GT.0) THEN
C      IF(IOTYPE.NE.3) THEN
c cpb 19/7/95 Generalising INP calculation
C        nn=0
C        IF(NIT(nb).EQ.1) THEN
C          DO n1=1,IBT(2,1)+1
C            nn=nn+1
C            INP(nn,1)=n1
C          ENDDO
C        ELSE IF(NIT(nb).EQ.2) THEN
C          DO n2=1,IBT(2,2)+1
C            DO n1=1,IBT(2,1)+1
C              nn=nn+1
C              INP(nn,1)=n1
C              INP(nn,2)=n2
C            ENDDO
C          ENDDO
C        ELSE IF(NIT(nb).EQ.3) THEN
C          DO n3=1,IBT(2,3)+1
C            DO n2=1,IBT(2,2)+1
C              DO n1=1,IBT(2,1)+1
C                nn=nn+1
C                INP(nn,1)=n1
C                INP(nn,2)=n2
C                INP(nn,3)=n3
C              ENDDO
C            ENDDO
C          ENDDO
C        ENDIF
C        CALL ASSERT(nn.LE.NNM,' >>Need to increase NNM',ERROR,*9999)
C        ENDIF

c cpb 14/7/95 Old NKT,IDO calculation
C        nk=0
C        IF(NIT(nb).EQ.1) THEN
C          DO ni1=1,IBT(1,1)
C            nk=nk+1
C            IDO(nk,1)=ni1
C          ENDDO
C        ELSE IF(NIT(nb).EQ.2) THEN
C          DO ni2=1,IBT(1,2)
C            DO ni1=1,IBT(1,1)
C              nk=nk+1
C              IDO(nk,1)=ni1
C              IDO(nk,2)=ni2
C            ENDDO
C          ENDDO
C        ELSE IF(NIT(nb).EQ.3) THEN
C          DO ni3=1,IBT(1,3)
C            DO ni2=1,IBT(1,2)
C              DO ni1=1,IBT(1,1)
C                IF(nk.LT.7) THEN
C                  nk=nk+1
C                  IDO(nk,1)=ni1
C                  IDO(nk,2)=ni2
C                  IDO(nk,3)=ni3
C                ENDIF
C              ENDDO
C            ENDDO
C          ENDDO
C        ENDIF !NIT

C LC 24/2/97 LC removed section from :
      SUBROUTINE BASIS2(IBT,IDO,INP,nb,NGAP,PG,WG,XIG,ERROR,*)


C#### Subroutine: BASIS2
C###  Description:
C###    BASIS2 inpust parameters for 2D or 3D Simplex/Serendipity/
C###    Sector basis functions and Gauss-Legendre or triangular 
C###    geometry quadrature.

c cpb 15/7/95 Generalising sectors
C        IBT(1,1)=5 !Sector Xi1 direction
C        IBT(1,2)=5 !Sector Xi2 direction
C        IBT(2,1)=SECT(1) !as input in Xi1 direction
C        IBT(2,2)=SECT(2) !as input in Xi2 direction
C
Cc        INP(1,1)=1 !Xi1 direction
Cc        INP(1,2)=1 !Xi2 direction
Cc        INP(2,1)=1 !Xi1 direction Xi1=0
Cc        INP(2,2)=2 !Xi2 direction
Cc        INP(3,1)=2 !Xi1 direction Xi1=1/3 
Cc        INP(3,2)=2 !Xi2 direction
Cc        INP(4,1)=3 !Xi1 direction Xi1=2/3 
Cc        INP(4,2)=2 !Xi2 direction
Cc        INP(5,1)=4 !Xi1 direction Xi1=1
Cc        INP(5,2)=2 !Xi2 direction
Cc        NNT(nb)=5  !5 nodes in cubic sector element
C
C        INP(1,1)=1 !set node 1 specially
C        INP(1,2)=1
C        n1=SECT(1)+1
C        IF(n1.EQ.5) n1=2
C        n2=SECT(2)+1
C        IF(n2.EQ.5) n2=2
C        nn=2
C        DO k2=2,n2
C          DO k1=1,n1
C            INP(nn,1)=k1
C            INP(nn,2)=k2
C            nn=nn+1
C          ENDDO
C        ENDDO
C        NNT(nb)=n1*(n2-1)+1

C PJH   9Aug91 - use default INP
C       IF(NST(nb).GT.0) THEN
C         CHAR1=' '
C         k=0
C         DO nn=1,NNT(nb)
C           DO ni=1,NIT(nb)
C             k=k+1
C             IDEFLT(k)=INP(nn,ni)
C             WRITE(CHAR1(k:k),'(I1)') IDEFLT(k)
C           ENDDO
C         ENDDO
C         CALL TRIM(CHAR1,IBEG,IEND)
C         FORMAT='($,'' Enter the node position indices ['//
C    '      CHAR1(IBEG:IEND)//']: '',40I2)'
C         IF(IOTYPE.EQ.3) THEN
C           k=0
C           DO nn=1,NNT(nb)
C             DO ni=1,NIT(nb)
C               k=k+1
C               IDATA(k)=INP(nn,ni)
C             ENDDO
C           ENDDO
C         ENDIF
C         CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,k,
C    '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,21,
C    '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C         IF(iotype.ne.3) THEN
C           k=0
C           DO nn=1,NNT(nb)
C             DO ni=1,NIT(nb)
C               k=k+1
C               INP(nn,ni)=IDATA(k)
C             ENDDO
C           ENDDO
C         ENDIF
C       ENDIF

C LC 24/2/97 removed all except call line :

c cpb 16/7/95 Generalising Sector
C        IF(NBSC(1,nb).ne.3) THEN
          CALL GAUSS2(IBT,INP,nb,NGAP(1,nb),PG,WG,XIG,ERROR,*9999)
C        ELSE !Sector
CC gmh 12/8/94 Set up IDO if it is cubic hermite in xi3 direction
CC             this bit taken from BASIS1
C
Cc cpb 14/7/95 Temporary nn loop
C          DO nn=1,NNT(nb)
C            IF(NIT(nb).EQ.2) THEN
C              nk=0
C              n1=1
C              IF(IBT(2,1).EQ.4) n1=2
C              n2=1
C              IF(IBT(2,2).EQ.4) n2=2
C              ni3=1
C              DO ni2=1,n2
C                DO ni1=1,n1
C                  IF(nk.LT.7) THEN
C                    nk=nk+1
C                    IDO(nk,nn,1)=ni1
C                    IDO(nk,nn,2)=ni2
C                    IDO(nk,nn,3)=ni3
C                  ENDIF
C                ENDDO
C              ENDDO
C            ELSE IF(NIT(nb).EQ.3) THEN 
C              nk=0
C              n1=1
C              IF(IBT(2,1).EQ.4) n1=2
C              n2=1
C              IF(IBT(2,2).EQ.4) n2=2
C              DO ni3=1,IBT(1,3)
C                DO ni2=1,n2
C                  DO ni1=1,n1
C                    IF(nk.LT.7) THEN
C                      nk=nk+1
C                      IDO(nk,nn,1)=ni1
C                      IDO(nk,nn,2)=ni2
C                      IDO(nk,nn,3)=ni3
C                    ENDIF
C                  ENDDO
C                ENDDO
C              ENDDO
C            ENDIF
C          ENDDO !nn
C          NKT(0,nb)=nk
C          NST(nb)=0
C          DO nn=1,NNT(nb)
C            NKT(nn,nb)=NKT(0,nb)
C            NST(nb)=NST(nb)+NKT(nn,nb)
C          ENDDO 
C          DO nn=1,NNT(nb)
C            DO nk=1,NKT(0,nb)
C              IDO1=IDO(nk,nn,1)
C              IDO2=1
C              IDO3=1
C              IF(NIT(nb).GE.2) IDO2=IDO(nk,nn,2)
C              IF(NIT(nb).EQ.3) IDO3=IDO(nk,nn,3)
C              IF(IDO1.EQ.1.AND.IDO2.EQ.1.AND.IDO3.EQ.1) IDO(nk,nn,0)=1
C              IF(IDO1.EQ.2.AND.IDO2.EQ.1.AND.IDO3.EQ.1) IDO(nk,nn,0)=2
C              IF(IDO1.EQ.1.AND.IDO2.EQ.2.AND.IDO3.EQ.1) IDO(nk,nn,0)=4
C              IF(IDO1.EQ.1.AND.IDO2.EQ.1.AND.IDO3.EQ.2) IDO(nk,nn,0)=7
C              IF(IDO1.EQ.2.AND.IDO2.EQ.2.AND.IDO3.EQ.1) IDO(nk,nn,0)=6
C              IF(IDO1.EQ.1.AND.IDO2.EQ.2.AND.IDO3.EQ.2) IDO(nk,nn,0)=9
C              IF(IDO1.EQ.2.AND.IDO2.EQ.1.AND.IDO3.EQ.2) IDO(nk,nn,0)=10
C            ENDDO !nk
C          ENDDO !nn
CC gmh 3/8/94  Sector is just cubic * linear lagrange, so we can
CC             just use GAUSS1 to set the points.
CC	      Set the NGAPs to 'standard' form ie 
CC             NGAP1=NGAP2=sqrt(num points input)
C          NGAP(3,nb)=NGAP(2,nb)
C          NGAP(2,nb)=JIDINT(DSQRT(DBLE(NGAP(1,nb))))
C          NGAP(1,nb)=NGAP(2,nb)
C          CALL GAUSS1(IBT,IDO,INP,nb,NGAP(1,nb),PG,WG,XIG,ERROR,*9999)
C        ENDIF


C LC 24/2/97 archived entire routine : 

C Old BASIS6 - replaced AJP 31-5-93
C      SUBROUTINE BASIS6(IBT,IDO,INP,nb,NGAP,PG,WG,XIG,ERROR,*)
C
CC#### Subroutine: BASIS6
CC###  Description:C
CC**** Singular Boundary Element basis function routine.
CC**** Input of parameter for 1-d basis functions designed to integrate
CC**** particular singular integrals.
CC**** Currently sets up basis function to handle r**(-1/2) type integrands
CC**** See description in GAUSSPWB (FE90) for further information.
C
C      IMPLICIT NONE
C      INCLUDE 'cmiss$reference:b00.cmn'
C      INCLUDE 'cmiss$reference:b01.cmn'
C      INCLUDE 'cmiss$reference:bem000.cmn'
C      INCLUDE 'cmiss$reference:cbdi02.cmn'
C      INCLUDE 'cmiss$reference:geom00.cmn'
C      INCLUDE 'cmiss$reference:inout00.cmn'
C      INCLUDE 'cmiss$reference:ityp00.cmn'
C      INCLUDE 'cmiss$reference:jtyp00.cmn'
C      INCLUDE 'cmiss$reference:sing00.cmn'
C!     Parameter List
C      INTEGER IBT(3,NIM),IDO(NKM,0:*),INP(NNM,NIM),nb,NGAP(NIM,NBM)
C      REAL*8 PG(NSM,NUM,NGM),WG(NGM),XIG(NIM,NGM)
C      CHARACTER ERROR*(*)
C!     Local Variables
C      INTEGER IB,IBEG,ICHAR,IDO1,IDO2,IDO3,IEND,INFO,K,n1,n2,ni,
C     '  ni1,ni2,ni3,nk,nn,NOQUES
C      CHARACTER ANS,BUFFER*80,CHAR1*100,CHAR2*3,CHAR3*3,
C     '  CONTYP*9,EMAP*9,FNAME*8,FPLOT*8,LFTYPE*8,TITLE*79,TPLOT*79
C      LOGICAL FILEIP,MORE
C
C      CALL ENTERS('BASIS6',*9999)
C      FILEIP=.FALSE.
C      NOQUES=0
C      INFO=2
C
C      DO K=1,100
C        CHAR1(K:K)=' '
C      ENDDO
C      FORMAT='($,'' Enter the number of Xi-coordinates [1]: '',I1)'
C      IF(IOTYPE.EQ.3) IDATA(1)=NIT(nb)
C      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NJT,
C     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C      IF(iotype.ne.3) NIT(nb)=IDATA(1)
C
C      IF(iotype.ne.3) THEN
C        NUT(nb)=NIT(nb)*NIT(nb)+2
C        NGT(nb)=1
C        NNT(nb)=1
C        DO ni=1,NIT(nb)
C          NGAP(ni,nb)=1
C          IBT(1,ni)=1
C          IBT(2,ni)=0
C          IBT(3,ni)=0
C        ENDDO !ni
C      ENDIF
C
C      DO ni=1,NIT(nb)
C        WRITE(CHAR1,'(I1)') ni
C        FORMAT='(/'' The interpolant in the Xi('//
C     '    CHAR1(1:1)//') direction is [1]: '''//
C     '    '/''   (1) Linear Lagrange'''//
C     '    '/''   (2) Quadr. Lagrange'''//
C     '    '/''   (3) Cubic  Lagrange'''//
C     '    '/''   (4) Cubic  Hermite'''//
C     '    '/$,''    '',I1)'
C        IF(IOTYPE.EQ.3) THEN
C          IF(IBT(1,ni).EQ.1) THEN
C            IB=IBT(2,ni)
C          ELSE IF(IBT(1,ni).EQ.2 ) THEN
C            IB=4
C          ENDIF
C          IDATA(1)=IB
C        ENDIF
C        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(iotype.ne.3) IB=IDATA(1)
C
C        IF(IB.EQ.1) THEN
C          IDEFLT(1)=2
C          IDEFLT(2)=10
C        ELSE IF(IB.GT.1) THEN
C          IDEFLT(1)=3
C          IDEFLT(2)=20
C        ENDIF
C        WRITE(CHAR2,'(I3)') IDEFLT(1)
C        WRITE(CHAR3,'(I3)') IDEFLT(2)
C        FORMAT='($,'' Enter the number of Gauss points in the Xi('//
C     '    CHAR1(1:1)//') direction for the''/'' 
C     '  low and high order schemes'
C     '      //'['//CHAR2(1:3)//','//CHAR3(1:3)//']: '',I3,I3)'
C        IF(IOTYPE.EQ.3) THEN
C          IDATA(1)=NGAP(ni,NBQ(NBBT(nb),nb))
C          IDATA(2)=NGAP(ni,NBQ(1,nb))
C        ENDIF
C        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,2,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,64,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(iotype.ne.3) THEN
C          NBQ(1,nb)=IDATA(1) !Low order
C          NBQ(2,nb)=IDATA(2) !High order
C          IF(NBQ(1,1).GT.NBQ(2,1)) THEN !Entered in the wrong order
C            NBQ(1,nb)=IDATA(2)
C            NBQ(2,nb)=IDATA(1)
C          ENDIF
C          NGAP(ni,nb)=NBQ(2,nb) !High order scheme
C
C          NGT(nb)=NGT(nb)*NGAP(ni,nb)
C          CALL ASSERT(NGT(nb).LE.NGM,' >>Need to increase NGM',
C     '      ERROR,*9999)
C          IF(IB.LE.3) THEN
C            IBT(1,ni)=1
C            IBT(2,ni)=IB
C            NNT(nb)=NNT(nb)*(IB+1)
C          ELSE
C            IBT(1,ni)=2
C            IBT(2,ni)=1
C            NNT(nb)=NNT(nb)*2
C          ENDIF
C        ENDIF
C      ENDDO
C
C      IF(iotype.ne.3) THEN
C        nk=0
C        IF(NIT(nb).EQ.1) THEN
C          DO ni1=1,IBT(1,1)
C            nk=nk+1
C            IDO(nk,1)=ni1
C          ENDDO
C        ELSE IF(NIT(nb).EQ.2) THEN
C          DO ni2=1,IBT(1,2)
C            DO ni1=1,IBT(1,1)
C              nk=nk+1
C              IDO(nk,1)=ni1
C              IDO(nk,2)=ni2
C            ENDDO
C          ENDDO
C        ENDIF
C        NKT(0,nb)=nk
C        NST(nb)=0
C        DO nn=1,NNT(nb) !AJP 25-5-93  All nkt(nn,nb) the same for now
C          NKT(nn,nb)=NKT(0,nb)
C          NST(nb)=NST(nb)+NKT(nn,nb)
C        ENDDO
C      ENDIF
C
C      IF(NST(nb).GT.0) THEN
C        IF(iotype.ne.3) THEN
C          nn=0
C          IF(NIT(nb).EQ.1) THEN
C            DO n1=1,IBT(2,1)+1
C              nn=nn+1
C              INP(nn,1)=n1
C            ENDDO
C          ELSE IF(NIT(nb).EQ.2) THEN
C            DO n2=1,IBT(2,2)+1
C              DO n1=1,IBT(2,1)+1
C                nn=nn+1
C                INP(nn,1)=n1
C                INP(nn,2)=n2
C              ENDDO
C            ENDDO
C          ENDIF
C          CALL ASSERT(nn.LE.NNM,' >>Need to increase NNM',ERROR,*9999)
C        ENDIF
C
C        K=0
C        DO nn=1,NNT(nb)
C          DO ni=1,NIT(nb)
C            K=K+1
C            IDEFLT(K)=INP(nn,ni)
C            WRITE(CHAR1(k:k),'(I1)') IDEFLT(k)
C          ENDDO
C        ENDDO
C        CALL TRIM(CHAR1,IBEG,IEND)
C        FORMAT='($,'' Enter the node position indices ['//
C     '    CHAR1(IBEG:IEND)//']: '',40I2)'
C        IF(IOTYPE.EQ.3) THEN
C          K=0
C          DO nn=1,NNT(nb)
C            DO ni=1,NIT(nb)
C              K=K+1
C              IDATA(K)=INP(nn,ni)
C            ENDDO
C          ENDDO
C        ENDIF
C        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,K,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,40,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(iotype.ne.3) THEN
C          K=0
C          DO nn=1,NNT(nb)
C            DO ni=1,NIT(nb)
C              K=K+1
C              INP(nn,ni)=IDATA(K)
C            ENDDO
C          ENDDO
C        ENDIF
C
C        IF(NKT(0,nb).GT.1) THEN
C          K=0
C          DO nk=1,NKT(0,nb)
C            DO ni=1,NIT(nb)
C              K=K+1
C              IDEFLT(K)=IDO(nk,ni)
C              WRITE(CHAR1(K:K),'(I1)') IDEFLT(K)
C            ENDDO
C          ENDDO
C          CALL TRIM(CHAR1,IBEG,IEND)
C          FORMAT='($,'' Enter the derivative order indices ['//
C     '      CHAR1(IBEG:IEND)//']: '',40I2)'
C          IF(IOTYPE.EQ.3) THEN
C            K=0
C            DO nk=1,NKT(0,nb)
C              DO ni=1,NIT(nb)
C                K=K+1
C                IDATA(K)=IDO(nk,ni)
C              ENDDO
C            ENDDO
C          ENDIF
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,K,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,40,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(iotype.ne.3) THEN
C            K=0
C            DO nk=1,NKT(0,nb)
C              DO ni=1,NIT(nb)
C                K=K+1
C                IDO(nk,ni)=IDATA(K)
C              ENDDO
C            ENDDO
C          ENDIF
C        ENDIF
C        DO nk=1,NKT(0,nb)
C          IDO1=IDO(nk,1)
C          IDO2=1
C          IDO3=1
C          IF(NIT(nb).GE.2) IDO2=IDO(nk,2)
C          IF(NIT(nb).EQ.3) IDO3=IDO(nk,3)
C          IF(IDO1.EQ.1.AND.IDO2.EQ.1.AND.IDO3.EQ.1) IDO(nk,0)=1
C          IF(IDO1.EQ.2.AND.IDO2.EQ.1.AND.IDO3.EQ.1) IDO(nk,0)=2
C          IF(IDO1.EQ.1.AND.IDO2.EQ.2.AND.IDO3.EQ.1) IDO(nk,0)=4
C          IF(IDO1.EQ.1.AND.IDO2.EQ.1.AND.IDO3.EQ.2) IDO(nk,0)=7
C          IF(IDO1.EQ.2.AND.IDO2.EQ.2.AND.IDO3.EQ.1) IDO(nk,0)=6
C          IF(IDO1.EQ.1.AND.IDO2.EQ.2.AND.IDO3.EQ.2) IDO(nk,0)=9
C          IF(IDO1.EQ.2.AND.IDO2.EQ.1.AND.IDO3.EQ.2) IDO(nk,0)=10
C        ENDDO
C      ENDIF
C
C      CALL EXITS('BASIS6')
C      RETURN
C 9999 CALL ERRORS('BASIS6',ERROR)
C      CALL EXITS('BASIS6')
C      RETURN 1
C      END

C removed entire routine : OLD BASIS7 Replace GBS 25-APR-1996

C      SUBROUTINE BASIS7(IBT,IDO,INP,nb,NGAP,PG,WG,XIG,ERROR,*)

C#### Subroutine: BASIS7
C###  Description:
C###    BASIS7 inputs parameters for Lagrange/Hermite tensor product 
C###    basis functions and Gauss-Legendre or Gauss-Lobatto quadrature.
C**** INP(nn,ni,nb) gives the index for element node nn in each Xi-dir.
C**** Thus: INP(nn,1..,nb) = 1,2,2 indicates that node nn is the first
C**** node in the Xi(1) dir and second in the Xi(2) & Xi(3) directions.

C      IMPLICIT NONE
C      INCLUDE 'cmiss$reference:b00.cmn'
C      INCLUDE 'cmiss$reference:b01.cmn'
C      INCLUDE 'cmiss$reference:cbdi02.cmn'
C      INCLUDE 'cmiss$reference:geom00.cmn'
C      INCLUDE 'cmiss$reference:inout00.cmn'
C!     Parameter List
C      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,
C     '  NGAP(NIM,NBM)
C      REAL*8 PG(NSM,NUM,NGM),WG(NGM),XIG(NIM,NGM)
C      CHARACTER ERROR*(*)
C!     Local Variables
C      INTEGER IB,IBEG,ICHAR,IDO1,IDO2,IDO3,IEND,INFO,k,n1,n2,n3,
C     '  ng,ni,ni1,ni2,ni3,nk,nn,NOQUES
C      CHARACTER CHAR1*100,CHAR2*1
C      LOGICAL FILEIP

C      CALL ENTERS('BASIS7',*9999)
C      FILEIP=.FALSE.
C      NOQUES=0
C      INFO=2
C
C      DO k=1,100
C        CHAR1(k:k)=' '
C      ENDDO
C      FORMAT='($,'' Enter the number of Xi-coordinates [1]: '',I1)'
C      IF(IOTYPE.EQ.3) IDATA(1)=NIT(nb)
C      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NJT,
C     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C      IF(iotype.ne.3) NIT(nb)=IDATA(1)
C      IF(IOTYPE.EQ.1.AND.NIT(nb).EQ.3) THEN
C        WRITE(OP_STRING,'(''>>Remember to define 2D face basis '','
C     '    //'''functions'')')
C      	CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      ENDIF
C
C      IF(iotype.ne.3) THEN
C        NUT(nb)=NIT(nb)*NIT(nb)+2
C        NGT(nb)=1
C        NNT(nb)=1
C        DO ni=1,NIT(nb)
C          NGAP(ni,nb)=1
C          IBT(1,ni)=1
C          IBT(2,ni)=0
C          IBT(3,ni)=0
C        ENDDO
C      ENDIF
C
C      NOQUES=0
C      DO ni=1,NIT(nb)
C        WRITE(CHAR1,'(I1)') ni
C        FORMAT='(/'' The interpolant in the Xi('//
C     '    CHAR1(1:1)//') direction is [1]: '''//
C     '    '/''   (1) Linear Lagrange'''//
C     '    '/''   (2) Quadr. Lagrange'''//
C     '    '/''   (3) Cubic  Lagrange'''//
C     '    '/''   (4) Cubic  Hermite'''//
C     '    '/$,''    '',I1)'
C        IF(IOTYPE.EQ.3) THEN
C          IF(IBT(1,ni).EQ.1) THEN
C            IB=IBT(2,ni)
C          ELSE IF(IBT(1,ni).EQ.2) THEN
C            IB=4
C          ENDIF
C          IDATA(1)=IB
C        ENDIF
C        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(iotype.ne.3) IB=IDATA(1)
C
C        IF(nb.GT.1.AND.NIT(nb).EQ.NIT(nb-1)) THEN
C          IDEFLT(1)=NGAP(ni,nb-1)
C        ELSE
C          IF(IB.EQ.1) THEN
C            IDEFLT(1)=2
C          ELSE IF(IB.GT.1) THEN
C            IDEFLT(1)=3
C          ENDIF
C        ENDIF
C        WRITE(CHAR2,'(I1)') IDEFLT(1)
C        FORMAT='($,'' Enter # Gauss points in Xi('
C     '    //CHAR1(1:1)//') direction ['//CHAR2(1:1)
C     '    //'] (incl. bdry): '',I1)'
C        IF(IOTYPE.EQ.3) IDATA(1)=NGAP(ni,nb)
C        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,11,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(iotype.ne.3) THEN
C          NGAP(ni,nb)=IDATA(1)
C          NGT(nb)=NGT(nb)*NGAP(ni,nb)
C          CALL ASSERT(NGT(nb).LE.NGM,' >>Need to increase NGM',
C     '      ERROR,*9999)
C          IF(IB.LE.3) THEN
C            IBT(1,ni)=1
C            IBT(2,ni)=IB
C            NNT(nb)=NNT(nb)*(IB+1)
C          ELSE
C            IBT(1,ni)=2
C            IBT(2,ni)=1
C            NNT(nb)=NNT(nb)*2
C          ENDIF
C          CALL ASSERT(NNT(nb).LE.NNM,' >>Need to increase NNM',
C     '      ERROR,*9999)
C        ENDIF
C      ENDDO
C
C      IF(iotype.ne.3) THEN
Cc cpb 14/7/95 Temporary nn loop
C        DO nn=1,NNT(nb)
C          nk=0
C          IF(NIT(nb).EQ.1) THEN
C            DO ni1=1,IBT(1,1)
C              nk=nk+1
C              IDO(nk,nn,1)=ni1
C            ENDDO
C          ELSE IF(NIT(nb).EQ.2) THEN
C            DO ni2=1,IBT(1,2)
C              DO ni1=1,IBT(1,1)
C                nk=nk+1
C                IDO(nk,nn,1)=ni1
C                IDO(nk,nn,2)=ni2
C              ENDDO
C            ENDDO
C          ELSE IF(NIT(nb).EQ.3) THEN
C            DO ni3=1,IBT(1,3)
C              DO ni2=1,IBT(1,2)
C                DO ni1=1,IBT(1,1)
C                  IF(nk.LT.7) THEN
C                    nk=nk+1
C                    IDO(nk,nn,1)=ni1
C                    IDO(nk,nn,2)=ni2
C                    IDO(nk,nn,3)=ni3
C                  ENDIF
C                ENDDO
C              ENDDO
C            ENDDO
C          ENDIF
C        ENDDO !nn
C        NKT(0,nb)=nk
C        NST(nb)=0
C        DO nn=1,NNT(nb) !AJP 25-5-93  All nkt(nn,nb) the same for now
C          NKT(nn,nb)=NKT(0,nb)
C          NST(nb)=NST(nb)+NKT(nn,nb)
C        ENDDO
C      ENDIF
C
C      IF(NST(nb).GT.0) THEN
C        IF(iotype.ne.3) THEN
C          nn=0
C          IF(NIT(nb).EQ.1) THEN
C            DO n1=1,IBT(2,1)+1
C              nn=nn+1
C              INP(nn,1)=n1
C            ENDDO
C          ELSE IF(NIT(nb).EQ.2) THEN
C            DO n2=1,IBT(2,2)+1
C              DO n1=1,IBT(2,1)+1
C                nn=nn+1
C                INP(nn,1)=n1
C                INP(nn,2)=n2
C              ENDDO
C            ENDDO
C          ELSE IF(NIT(nb).EQ.3) THEN
C            DO n3=1,IBT(2,3)+1
C              DO n2=1,IBT(2,2)+1
C                DO n1=1,IBT(2,1)+1
C                  nn=nn+1
C                  INP(nn,1)=n1
C                  INP(nn,2)=n2
C                  INP(nn,3)=n3
C                ENDDO
C              ENDDO
C            ENDDO
C          ENDIF
C          CALL ASSERT(nn.LE.NNM,' >>Need to increase NNM',ERROR,*9999)
C        ENDIF
C
C        k=0
C        DO nn=1,NNT(nb)
C          DO ni=1,NIT(nb)
C            k=k+1
C            IDEFLT(k)=INP(nn,ni)
C            WRITE(CHAR1(k:k),'(I1)') IDEFLT(k)
C          ENDDO
C        ENDDO
C        CALL TRIM(CHAR1,IBEG,IEND)
C        FORMAT='($,'' Enter the node position indices ['//
C     '    CHAR1(IBEG:IEND)//']: '',40I2)'
C        IF(IOTYPE.EQ.3) THEN
C          k=0
C          DO nn=1,NNT(nb)
C            DO ni=1,NIT(nb)
C              k=k+1
C              IDATA(k)=INP(nn,ni)
C            ENDDO
C          ENDDO
C        ENDIF
C        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,k,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,40,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(iotype.ne.3) THEN
C          k=0
C          DO nn=1,NNT(nb)
C            DO ni=1,NIT(nb)
C              k=k+1
C              INP(nn,ni)=IDATA(k)
C            ENDDO
C          ENDDO
C        ENDIF
C
C        IF(NKT(0,nb).GT.1) THEN
C          k=0
C          DO nn=1,NNT(nb)
C            DO nk=1,NKT(nn,nb)
C              DO ni=1,NIT(nb)
C                k=k+1
C                IDEFLT(k)=IDO(nk,nn,ni)
C                WRITE(CHAR1(k:k),'(I1)') IDEFLT(k)
C              ENDDO
C            ENDDO
C          ENDDO !nn
C          CALL TRIM(CHAR1,IBEG,IEND)
C          FORMAT='($,'' Enter the derivative order indices ['//
C     '      CHAR1(IBEG:IEND)//']: '',40I2)'
C          IF(IOTYPE.EQ.3) THEN
C            k=0
C            DO nn=1,NNT(nb)
C              DO nk=1,NKT(nn,nb)
C                DO ni=1,NIT(nb)
C                  k=k+1
C                  IDATA(k)=IDO(nk,nn,ni)
C                ENDDO
C              ENDDO
C            ENDDO !nn
C          ENDIF
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,k,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,40,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(iotype.ne.3) THEN
C            k=0
C            DO nn=1,NNT(nb)
C              DO nk=1,NKT(nn,nb)
C                DO ni=1,NIT(nb)
C                  k=k+1
C                  IDO(nk,nn,ni)=IDATA(k)
C                ENDDO
C              ENDDO
C            ENDDO !nn
C          ENDIF
C        ENDIF
C        DO nn=1,NNT(nb)
C          DO nk=1,NKT(nn,nb)
C            IDO1=IDO(nk,nn,1)
C            IDO2=1
C            IDO3=1
C            IF(NIT(nb).GE.2) IDO2=IDO(nk,nn,2)
C            IF(NIT(nb).EQ.3) IDO3=IDO(nk,nn,3)
C            IF(IDO1.EQ.1.AND.IDO2.EQ.1.AND.IDO3.EQ.1) IDO(nk,nn,0)=1
C            IF(IDO1.EQ.2.AND.IDO2.EQ.1.AND.IDO3.EQ.1) IDO(nk,nn,0)=2
C            IF(IDO1.EQ.1.AND.IDO2.EQ.2.AND.IDO3.EQ.1) IDO(nk,nn,0)=4
C            IF(IDO1.EQ.1.AND.IDO2.EQ.1.AND.IDO3.EQ.2) IDO(nk,nn,0)=7
C            IF(IDO1.EQ.2.AND.IDO2.EQ.2.AND.IDO3.EQ.1) IDO(nk,nn,0)=6
C            IF(IDO1.EQ.1.AND.IDO2.EQ.2.AND.IDO3.EQ.2) IDO(nk,nn,0)=9
C            IF(IDO1.EQ.2.AND.IDO2.EQ.1.AND.IDO3.EQ.2) IDO(nk,nn,0)=10
C          ENDDO !nk
C        ENDDO !nn
C      ENDIF
C
C      CALL GAUSS7(IBT,IDO,INP,nb,NGAP(1,nb),PG,WG,XIG,ERROR,*9999)
C
C      IF(DOP) THEN
C        DO ng=1,NGT(nb)
C          WRITE(OP_STRING,'('' XIG(ni,'',I2,''): '',3E11.3)')
C     '       ng,(XIG(ni,ng),ni=1,NIT(nb))
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDDO
C      ENDIF
C
C      CALL EXITS('BASIS7')
C      RETURN
C 9999 CALL ERRORS('BASIS7',ERROR)
C      CALL EXITS('BASIS7')
C      RETURN 1
C      END

Module FE05
===========


C removed section from routine :

      SUBROUTINE FIBRE2(INDEX,IBT,IDO,INP,IW,NAN,NBJ,nr,XE,XG,
     '  ERROR,*)

C#### Subroutine: FIBRE2
C###  Description:
C###    FIBRE2 draws fitted fibre field vectors (DXIF) on curvilinear 
C###    plane defined by constant Xi(3) in element ne at increments of 
C###    DXI1 and DXI2 in Xi(1) and Xi(2) directions, respectively.

C old MPN 30-Apr-96: old way of calculating material axes
C          nb=NBJ(NJ_LOC(NJL_FIBR,1))
C          IF(NNT(nb).GT.0) THEN
C            IF(JTYP9.GT.0) THEN !fibre angle field defined
C              ETA=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C     '          XE(1,NJ_LOC(NJL_FIBR,1)))
C            ENDIF
C          ELSE IF(NAT(nb).GT.0) THEN
C            ETA=XA(1,NJ_LOC(NJL_FIBR,1),ne)
C          ENDIF
C          CALL XMG(IBT,IDO,INP,NBJ,NJE,GXL,XE,XI,ERROR,*9999)
C          RWX=GXL(1,1)*GXL(2,2)-GXL(1,2)**2
C          IF(DOP) THEN
C            WRITE(OP_STRING,'('' GXL:'',3E11.3)') 
C     '        GXL(1,1),GXL(2,2),GXL(1,2)
C      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C            WRITE(OP_STRING,'('' RWX= '',E11.3)') RWX
C      	    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
C          IF(JTYP12.EQ.1) THEN
C            ETA1=ETA
C            ETA2=DACOS(GXL(1,2)/DSQRT(GXL(1,1)*GXL(2,2)))-ETA1
C          ELSE IF(JTYP12.EQ.2) THEN
C            ETA2=ETA
C            ETA1=DACOS(GXL(1,2)/DSQRT(GXL(1,1)*GXL(2,2)))-ETA2
C          ENDIF
C          C1=DCOS(ETA1)
C          C2=DCOS(ETA2)
C          DXINU(1,1)=( GXL(2,2)*DSQRT(GXL(1,1))*C1
C     '                -GXL(1,2)*DSQRT(GXL(2,2))*C2)/RWX
C          DXINU(2,1)=(-GXL(1,2)*DSQRT(GXL(1,1))*C1
C     '                +GXL(1,1)*DSQRT(GXL(2,2))*C2)/RWX
C          IF(DABS(C2).GT.1.d-6) THEN
C            RK=DSQRT(GXL(1,1)/GXL(2,2))*C1/C2
C            DXINU(1,2)=-1.d0/DSQRT(GXL(1,1)+RK*(RK*GXL(2,2)-
C     '        2.d0*GXL(1,2)))
C            DXINU(2,2)=-RK*DXINU(1,2)
C          ELSE
C            DENOM=DSQRT(RWX*GXL(1,1))
C            DXINU(1,2)=-GXL(1,2)/DENOM
C            DXINU(2,2)= GXL(1,1)/DENOM
C          ENDIF


C removed routine from :

      SUBROUTINE SHEET1(INDEX,IBT,IDO,INP,IW,NAN,NBJ,ne,NEELEM,
     '  NKE,NPF,NPNE,NQE,NRE,NVJE,NXI,
     '  DXI2,DXI3,SE,THETA,XA,XE,XG,XP,ERROR,*)

C#### Subroutine: SHEET1
C###  Description:
C###    SHEET1 draws fitted sheet angles on constant Xi_1,X plane.

C CPB 30/3/93 Not sure about the nj locations for this commented section
C         nb=NBJ(NJT+1)
C         IF(NNT(nb).GT.0) THEN
C           ETA=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C    '        XE(1,NJT+2))
C         ELSE IF(NAT(nb).GT.0) THEN
C           ETA=XA(1,NJT+2,ne)
C         ENDIF
C         CALL XMG(IBT,IDO,INP,NBJ,NJE,GL,XE,XI,ERROR,*9999)
C         G=GL(1,1)*GL(2,2)-GL(1,2)**2
C         IF(DOP) THEN
C           WRITE(IO4,'('' GL:'',3E11.3)') GL(1,1),GL(2,2),GL(1,2)
C           WRITE(IO4,'('' G= '',E11.3)') G
C         ENDIF
C         IF(JTYP12.EQ.1) THEN
C           ETA1=ETA
C           ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))-ETA1
C         ELSE IF(JTYP12.EQ.2) THEN
C           ETA2=ETA
C           ETA1=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))-ETA2
C         ENDIF
C         C1=DCOS(ETA1)
C         C2=DCOS(ETA2)
C         DXINU(1,1)=( GL(2,2)*DSQRT(GL(1,1))*C1
C    '                -GL(1,2)*DSQRT(GL(2,2))*C2)/G
C         DXINU(2,1)=(-GL(1,2)*DSQRT(GL(1,1))*C1
C    '                +GL(1,1)*DSQRT(GL(2,2))*C2)/G
C         IF(DABS(C2).GT.1.d-6) THEN
C           RK=DSQRT(GL(1,1)/GL(2,2))*C1/C2
C           DXINU(1,2)=-1.d0/DSQRT(GL(1,1)+RK*(RK*GL(2,2)-2.d0*GL(1,2)))
C           DXINU(2,2)=-RK*DXINU(1,2)
C         ELSE
C           DENOM=DSQRT(G*GL(1,1))
C           DXINU(1,2)=-GL(1,2)/DENOM
C           DXINU(2,2)= GL(1,1)/DENOM
C         ENDIF
C         SCALE=DSQRT(GL(1,1))
C         DX1=DXIF*SCALE*DXINU(1,1)
C         DX2=DXIF*SCALE*DXINU(2,1)
C         IF(DOP) THEN
C           WRITE(IO4,'('' DXIF='',E11.3,'' DX1='',E11.3,
C    '        '' DX2='',E11.3)') DXIF,DX1,DX2
C         ENDIF
C         XI(1)=XIP1-DX1/2.d0
C         XI(2)=XIP2-DX2/2.d0
C         IF(DOP) THEN
C           WRITE(IO4,'('' XI(1)='',E11.3,'' XI(2)='',E11.3)') XI(1),XI(2)
C         ENDIF
C         IF(PROJEC(1:2).EQ.'XI') THEN
C           Z(1)=XI(1)
C           Z(2)=XI(2)
C         ELSE
C           DO nj=1,NJT
C             nb=NBJ(nj)
C             X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C    '          XE(1,nj))
C           ENDDO
C           IF(DOP) THEN
C             WRITE(IO4,'('' X(1)='',E11.3,'' X(2)='',E11.3)') X(1),X(2)
C           ENDIF
C           IF(IW.LE.3) THEN
C             CALL XZ(ITYP10(1),X,Z)
C           ELSE IF(IW.EQ.4) THEN
C             Z(2)=X(2)
C             Z(3)=X(3)
C           ENDIF
C         ENDIF
C         XX(1)=Z(1)
C         YY(1)=Z(2)
C         ZZ(1)=Z(3)
C         XI(1)=XIP1+DX1/2.d0
C         XI(2)=XIP2+DX2/2.d0
C
C         IF(PROJEC(1:2).EQ.'XI') THEN
C           Z(1)=XI(1)
C           Z(2)=XI(2)
C         ELSE
C           DO nj=1,NJT
C             nb=NBJ(nj)
C             X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
C    '          XE(1,nj))
C           ENDDO
C           IF(IW.LE.3) THEN
C             CALL XZ(ITYP10(1),X,Z)
C           ELSE IF(IW.EQ.4) THEN
C             Z(2)=X(2)
C             Z(3)=X(3)
C           ENDIF
C         ENDIF
C         XX(2)=Z(1)
C         YY(2)=Z(2)
C         ZZ(2)=Z(3)
C         DO n=1,2
C           Z3D(1,n)=XX(n)
C           Z3D(2,n)=YY(n)
C           Z3D(3,n)=ZZ(n)
C         ENDDO
C         CALL POLYLINE(INDEX,IW,2,Z3D,ERROR,*9999)
C       ENDDO
C     ENDDO


      REAL*8 FUNCTION PS1()

C#### Function: PS1
C###  Type: REAL*8
C###  Description:
C###    PS1 evaluates 1D linear B-spline polynomial coefficients.

      IMPLICIT NONE
!     Parameter List

      PS1=1.0d0

      RETURN
      END


      REAL*8 FUNCTION PS2(K,N,SP,SP1)

C#### Function: PS2
C###  Type: REAL*8
C###  Description:
C###    PS2 evaluates 1D quadratic B-spline polynomial coefficients.
C**** K is polynomial term index
C**** N is basis function  index

      IMPLICIT NONE
!     Parameter List
      INTEGER K,N
      REAL*8 SP,SP1

      GO TO (10,20,30),N
 10     GO TO (11,12,13),K
 11       PS2=SP/(1.0D0+SP)
          RETURN
 12       PS2=-2.0D0*SP/(1.0D0+SP)
          RETURN
 13       PS2=SP/(1.0D0+SP)
          RETURN
 20     GO TO (21,22,23),K
 21       PS2=1.0D0/(1.0D0+SP)
          RETURN
 22       PS2=1.0D0/(1.0D0+SP)
          RETURN
 23       PS2=-(SP/(1.0D0+SP)+1.0D0/(1.0D0+SP1))
          RETURN
 30     GO TO (31,32,33),K
 31       PS2=0.0D0
          RETURN
 32       PS2=0.0D0
          RETURN
 33       PS2=1.0D0/(1.0D0+SP1)
          RETURN
      END


      REAL*8 FUNCTION PS3()

C#### Function: PS3
C###  Type: REAL*8
C###  Description:
C###    PS3 evaluates 1D cubic B-spline polynomial coefficients.
C**** K is polynomial term index
C**** N is basis function  index

      IMPLICIT NONE
!     Parameter List

      PS3=1.0d0

      RETURN
      END

      

      REAL*8 FUNCTION PSE(K1,K2,K3,N1,N2,N3,IBT,NITB)

C#### Function: PSE
C###  Type: REAL*8
C###  Description:
C###    PSE evaluates tensor product B-spline basis functions.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),K1,K2,K3,N1,N2,N3,NITB
!     Local Variables
      INTEGER K(3),N(3),ni
      REAL*8 PS1,PS2,PS3

      K(1)=K1
      K(2)=K2
      K(3)=K3
      N(1)=N1
      N(2)=N2
      N(3)=N3
      PSE=1.0D0
      DO ni=1,NITB
        IF(IBT(2,ni).EQ.1) PSE=PSE*PS1()
        IF(IBT(2,ni).EQ.2) PSE=PSE*PS2(K(ni),N(ni),1.0d0,1.0d0)
        IF(IBT(2,ni).EQ.3) PSE=PSE*PS3()
      ENDDO
      RETURN
      END


      REAL*8 FUNCTION PSI3(I1,I2,I3,NITB,nu,XI)

C#### Function: PSI3
C###  Type: REAL*8
C###  Description:
C###    PSI3 evaluates tensor product B-spline basis functions at XI.

      IMPLICIT NONE
!     Parameter List
      INTEGER I1,I2,I3,NITB,nu
      REAL*8 XI(3)
!     Local Variables
      INTEGER I(3),IPU(11,3),ni
      REAL*8 PSL

      DATA IPU/1,2,3,1,1,2,1,1,2,1,2,
     '         1,1,1,2,3,2,1,1,1,2,2,
     '         1,1,1,1,1,1,2,3,2,2,2/

      I(1)=I1
      I(2)=I2
      I(3)=I3
      PSI3=1.0D0
      DO ni=1,NITB
        PSI3=PSI3*PSL(I(ni),IPU(nu,ni),ni,XI(ni))
      ENDDO

      RETURN
      END


      REAL*8 FUNCTION PSL(I,K,ni,XI)

C#### Function: PSL
C###  Type: REAL*8
C###  Description:
C###    PSL evaluates linear,quadratic or cubic spline basis functions.
C**** I is index for the B-spline
C**** K is partial derivative index
C     Note: this is not a good way of doing this- all the B-splines are
C     evaluated and just one is returned. But what the hell...

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:bspln00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER I,K,ni
      REAL*8 XI         
!     Local Variables
      INTEGER i1,J,JMAX,left,MO,N
      PARAMETER (JMAX=20)
      REAL*8 B(JMAX+3),DL_LOCAL(4),DR(4),SAVED,TERM
      CHARACTER ERROR*10

C     find the knot to the left of XI
      DO left=1,NTKN(ni)
        IF(BKNOT(left,ni).GT.XI) GOTO 10
      ENDDO
 10   left=MAX(1,left-1)
      MO=MBSPL(ni)+1
      B(1)=1.0D0
      DO i1=2,6
        B(i1)=0.0D0
      ENDDO
      IF(DOP) THEN
        WRITE(OP_STRING,*)' bknot(..,1):',(BKNOT(i1,1),i1=1,NTKN(1))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*)' bknot(..,2):',(BKNOT(i1,2),i1=1,NTKN(2))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*)' psl: left,i,k,ni,xi,mo=',
     '    left,i,k,ni,xi,mo
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      DO J=1,MO-K
        DR(J)=BKNOT(left+J,ni)-XI
        DL_LOCAL(J)=XI-BKNOT(left+1-J,ni)
        SAVED=0.0D0
        DO N=1,J
          TERM=B(N)/(DR(N)+DL_LOCAL(J+1-N))
          B(N)=SAVED+DR(N)*TERM
          SAVED=DL_LOCAL(J+1-N)*TERM
        ENDDO
        B(J+1)=SAVED
      ENDDO
      IF(DOP) THEN
        WRITE(OP_STRING,*)' b=',(B(i1),i1=1,MO-K+1)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
C
      IF(K.EQ.1) THEN
        PSL=B(I)
      ELSE IF(K.EQ.2) THEN
        PSL=(MO-1)*(-B(I+1)/(BKNOT(left+MO,ni)-BKNOT(left+1,ni))
     '              +B(I+1)/(BKNOT(left+MO-1,ni)-BKNOT(left,ni)))
      ELSE IF(K.EQ.3) THEN
        PSL=(MO-1)*((MO-2)*(-B(I+2)/(BKNOT(left+MO,ni)-BKNOT(left+2,ni))
     '                      +B(I+1)/(BKNOT(left+MO-1,ni)
     '                      -BKNOT(left+1,ni)))/(BKNOT(left+MO,ni)
     '                      -BKNOT(left+1,ni))+
     '              (MO-2)*(-B(I+1)/(BKNOT(left+MO-1,ni)
     '                      -BKNOT(left+1,ni))+B(I)
     '                      /(BKNOT(left+MO-2,ni)-BKNOT(left,ni))
     '                     )/(BKNOT(left+MO-1,ni)-BKNOT(left,ni)))
      ENDIF

 9999 RETURN
      END


Module FE06
===========

      SUBROUTINE SHEET3(INDEX,IBT,IDO,INP,NAN,NBJ,ne,nr,IW,
     '  XE,XG,XID,ZD,ERROR,*)

C#### Subroutine: SHEET3
C###  Description:
C###    SHEET3 draws measured sheet angles on Xi_3,X plane 
C###    (CALC_SHEET=.TRUE.).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b14.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),IW,NAN(NIM,NAM,NBFM),NBJ(NJM,NEM),ne,nr
      REAL*8 XE(NSM,NJM),XG(NJM,NUM),XID(NIM),ZD(NJM)
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 GAMA,GAMA_proj,DELTA,DX,
     '  FITTED_ALFA,
     '  FITTED_A_VECTOR(3),FITTED_B_VECTOR(3),FITTED_C_VECTOR(3),
     '  F_VECTOR(3),NORM_LINE(3),
     '  POINTS(3,2),RAD,C_xfv(3),X

      CALL ENTERS('SHEET3',*9999)

      CALL ASSERT(.FALSE.,'>> Old code - needs updating',
     '  ERROR,*9999)

      CALL ASSERT(IW.EQ.13,
     '  ' Incorrect workstation ID: sb 13 for sheets',ERROR,*9999)

C new MAT_VEC_XI does not pass back
C FITTED_ALFA,FITTED_BETA,FITTED_GAMA,DELTA,PHI any more.
      CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ(1,ne),nr,
     '  FITTED_A_VECTOR,FITTED_B_VECTOR,FITTED_C_VECTOR,
     '  XE,XG,XID,ERROR,*9999)
     

C!     Calculate PHI and DELTA (ignore all fitted quantities)
C      CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ(1,ne),nr,
C     '  FITTED_A_VECTOR,FITTED_B_VECTOR,FITTED_C_VECTOR,
C     '  FITTED_ALFA,FITTED_BETA,FITTED_GAMA,DELTA,PHI,
C     '  XE,XID,ERROR,*9999)

      RAD=DSQRT(ZD(2)**2+ZD(3)**2) !this needs fixing!!!!!!!!!!!!!
      X=ZD(1)
      GAMA=ZD(6)
      !h_vector in xfv coords
C      H_xfv(1)=DSIN(DELTA)
C      H_xfv(2)=0.d0
C      H_xfv(3)=DCOS(DELTA)
      !GAMA vector in xfv coords

!      Ian Le Grice 28-3-92
!      Define line of intersection of myocardial sheet and
!       (v,x)-plane in x,g,v coordinate system

      F_VECTOR(1)= 0.d0
      F_VECTOR(2)= 1.d0
      F_VECTOR(3)= 0.d0

      C_xfv(1)= DCOS(FITTED_ALFA)*DCOS(GAMA)*DCOS(DELTA)
     '          -DSIN(GAMA)*DSIN(DELTA)
      C_xfv(2)= DSIN(FITTED_ALFA)*DCOS(GAMA)
      C_xfv(3)=-DCOS(FITTED_ALFA)*DCOS(GAMA)*DSIN(DELTA)
     '          +DSIN(GAMA)*DCOS(DELTA)

      NORM_LINE(1)= F_VECTOR(2)*C_xfv(3)-F_VECTOR(3)*C_xfv(2)
      NORM_LINE(2)= F_VECTOR(3)*C_xfv(1)-F_VECTOR(1)*C_xfv(3)
      NORM_LINE(3)= F_VECTOR(1)*C_xfv(2)-F_VECTOR(2)*C_xfv(1)

      GAMA_proj=DATAN2(-NORM_LINE(1),NORM_LINE(3))

!      GAMA_xfv(1)=-DCOS(FITTED_ALFA)*DSIN(GAMA)*DCOS(DELTA)+DCOS(GAMA)*DSIN(DELTA)
!      GAMA_xfv(3)= DCOS(FITTED_ALFA)*DSIN(GAMA)*DSIN(DELTA)+DCOS(GAMA)*DCOS(DELTA)
!      ABS_GAMA_xv=DSQRT(GAMA_xfv(1)**2+GAMA_xfv(3)**2)
!      !unit projection of GAMA vector in xv-plane
!      GAMA_xv(1)=GAMA_xfv(1)/ABS_GAMA_xv
!      GAMA_xv(2)=GAMA_xfv(3)/ABS_GAMA_xv
!      GAMA_proj=DATAN2(-GAMA_xv(1),GAMA_xv(2))

!      WRITE(*,'('' rad='',E12.3,'' x='',E12.3,'' gama='',E12.3,
!     '  '' gama_proj='',E12.3)') RAD,X,GAMA,GAMA_PROJ

      DX=0.02D0*FOCUS
      POINTS(1,1)=X  -DX*DSIN(GAMA_proj)
      POINTS(2,1)=RAD+DX*DCOS(GAMA_proj)
      POINTS(1,2)=X  +DX*DSIN(GAMA_proj)
      POINTS(2,2)=RAD-DX*DCOS(GAMA_proj)
      CALL POLYLINE(INDEX,IW,2,POINTS,ERROR,*9999)

      CALL EXITS('SHEET3')
      RETURN
 9999 CALL ERRORS('SHEET3',ERROR)
      CALL EXITS('SHEET3')
      RETURN 1
      END

Module FE07
===========

C 24/2/97 LC  removed 2 section from :

C#### Subroutine: GENSOL
C###  Description: 
C###    GENSOL calls solution routines.

c cpb 1/5/95 Replace Fourier analysis with Quasi-static analysis
C          ELSE IF(ITYP5(nr,nx).EQ.4) THEN !problem is Fourier analysis
C            CALL TRIM(FILE02,IBEG,IEND)
C            IUNIT=IOFILE2
C            CALL OPENF(IUNIT,'DISK',FILE02(IBEG:IEND)//'.fourier',
C     '        'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
C            NSTEP=DLOG10(F1/F0)/dFreq
C            WRITE(IUNIT,'('' Fourier analysis: Freq range '',D12.3,'
C     '	      //''' to '',D12.3,'' Log step = '',D12.3,'
C     '        //''' Number steps = '',I4)') F0,F1,dFreq,NSTEP+1
C            WRITE(IUNIT,'('' Number of nodes = '',I5)') NPT(nr)
C            DO istep=0,NSTEP
C              F=F0*10.0D0**(istep*dFreq)
C              OMEGA=2.0D0*PI*F
C              WRITE(OP_STRING,'('' Fourier analysis: Freq='',D12.3,'
C     '	        //''' omega='',D12.3)')F,OMEGA
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C              WRITE(IUNIT,'(/'' Freq='',D12.3,'' omega='',D12.3)') 
C     '          F,OMEGA
C              IF(istep.EQ.0) THEN
C                IP=1
C              ELSE
C                IP=2
C              ENDIF
C
CC cpb 12/3/95 Solve6 to fe_archive
C
C              CALL SOLVE6(IP,IBT,IDO,INP,LGE,
C     '          NBH,NBJ,NEELEM,NHE,NJE,NKE,NKH(1,1,1,nr),
C     '          NONY(0,1,1,nr),NPF,NPNE,NQE,
C     '          nr,NRE,NVJE,NVJP,NW,nx,NYNE,NYNP,NYNR(0,0,1,nr),
C     '          CE,CG,CONY(0,1,1,nr),CP,
C     '          ED,EM,ER,ES,GD,GK,GM,GR,OMEGA,
C     '          PG,RG,SE,VE,WG,XA,XE,XG,XP,YG,
C     '          YP,ZE,
C     '          ZG,WK1,FIX,
C     '          SYMM,GKC,GKKC,GRRC,XOC,ERROR,*9999)
C              DO no_nynr=1,NYNR(0,1,1,nr)
C                ny=NYNR(no_nynr,1,1,nr)
C                DO noy=1,NONY(0,ny,nrc,nr)
C                  no=NONY(noy,ny,nrc,nr)
C                  XOC_ABS=CDABS(XOC(no))
C                  XOC_ANG=DATAN2(DIMAG(XOC(no)),DREAL(XOC(no)))
C                  WRITE(IUNIT,'('' ny ='',I4,'' XOC('',I4,'')= '
C     '              //''',2D12.4,'' abs(XOC)='',D12.4,'' ang(XOC)='
C     '              //''',D12.4)')ny,no,XOC(no),XOC_ABS,XOC_ANG
C                ENDDO
C              ENDDO
C              !put data into spreadsheet
c             no=2
c             NOROW=istep+1
c             DATA(NOROW,1)=OMEGA
c             DATA(NOROW,2)=CDABS(XOC(no))
c             DATA(NOROW,3)=DATAN2(DIMAG(XOC(no)),DREAL(XOC(no)))
C            ENDDO
C            CALL CLOSEF(IUNIT,ERROR,*9999)
c           NT_DAT=NSTEP+1
c           NT_COL=3
c           NAME_COL(1)='frequency (rads/s)'
c           NAME_COL(2)='magnitude'
c           NAME_COL(3)='phase'


C 24/2/97 LC removed section from : 

C#### Subroutine: SOLVE9
C###  Description:
C###    SOLVE9 finds solution of an unsymmetric system of linear 
C###    equations no=1,NOT(1,1,0,nx) for coupled FEM/BEM problems.

C*** Solve reordered system
c cpb 18/5/95 Old GMRES solver
CC       level of fill-in
C        DO no1=1,NOT(1,1,0,nx)
C          FILLIN_GMRES(no1)=0 !diagonal only
CC          FILLIN_GMRES(no1)=1 !no fill-in
C        ENDDO
CC       NOPCNV needs to be called only if the sparsity pattern has been
CC         changed.  It calculates the symbolic factorization
C        IF(UPDATE_MATRIX) THEN
C          CALL NOPCNV(
C     '      IWK1_GMRES, !working integer array dimensioned to at least 
CC                        the number of non-zeros
C     '      FILLIN_GMRES, !level of fill-in for incomplete LU 
CC                          factorization (in) dimensioned to the 
CC                          number of rows
C     '      NOT(1,1,0,nx), !number of rows
C     '      ISC_GKK, !the sparsity struct for the matrix (in/unchanged)
C     '      IS_PRECON_GMRES, !the sparsity structure for the incomplete 
CC                             LU factorization (out)
C     '      IPIVOT_GMRES, !pivot locations, dimensioned to the number
CC                          of rows
C     '      IWK2_GMRES, !working array, dimensioned to the no of rows
C     '      NZ_GMRES_M !dimensioned size of the working array for 
CC                       storing the incomplete LU factorization (WK1).  
C     '      )
C        ENDIF
CC       need to be called for every solution
C        NUM_ORTHOG=NUM_GMRES_ORTHOG(nx) !number of orthogonalisations
C        ABSOLUTE_TOLERANCE=GMRES_ABSTOL !needs to be set each time
C        MAX_ITER=MAX_GMRES_ITERS !maximum number of iterations, 
CC                                  needs to be set each time
C        CALL NRMLZ1(
C     '    NOT(1,1,0,nx), !number of rows
C     '    GKK, !the matrix
C     '    GRR, !right hand side
C     '    ISC_GKK, !the sparsity structure for the matrix (in/unchanged)
C     '    ABSOLUTE_TOLERANCE, !absolute tolerance
C     '    XO !solution vector
C     '    )
C        CALL SOLVE_GMRES(
C     '    NOT(1,1,0,nx), !number of rows
C     '    GKK, !the matrix
C     '    PRECON_GMRES, !working array for storing the pre-conditioner
C     '    GRR, !right hand side
C     '    ISC_GKK, !the sparsity structure for the matrix (in/unchanged)
C     '    IS_PRECON_GMRES, !the sparsity structure for the incomplete 
CC                            LU factorization (out)
C     '    NUM_ORTHOG, !number of previous search dirs that are saved
C     '    ORTHOG_GMRES, !working array for storing the previous search
CC                        directions
C     '    ABSOLUTE_TOLERANCE, !absolute tolerance
C     '    IPIVOT_GMRES, !pivot locations, dimensioned to the number 
CC                         of rows
C     '    WK1_GMRES, !working array dimensioned to the number of rows
C     '    XO, !solution vector
C     '    WK2_GMRES, !working array dimensioned to the number of rows
C     '    WK3_GMRES, !working array dimensioned to the number of rows
C     '    QR_GMRES, !working array for storing the QR factorization,
CC                     dimensioned at least (NUM_ORTHOG+1)*(NUM_ORTHG+1)
C     '    C_GMRES, !working array for storing the cosines, dimensioned
CC                    at least NUM_ORTHOG+1
C     '    S_GMRES, !working array for storing the sines, dimensioned
CC                    at least NUM_ORTHOG+1
C     '    G_GMRES, !working array for storing the error norms, 
CC                    dimensioned at least NUM_ORTHOG+1
C     '    Y_GMRES, !working array for storing orthogonalization 
CC                    coefficients, dimensioned at least NUM_ORTHOG+1
C     '    MAX_ITER !maximum number of iterations
C     '  )



C CPB 27/3/96 This is the old assemble3 before I rewrote it.
      SUBROUTINE ASSEMBLE3(IBT,IDO,INP,LGE,NBH,NBJ,NEELEM,NHE,NHP,NJE,
     '  NKE,NKH,NPF,NPNE,NPNODE,NQE,nr,NVHE,NVHP,
     '  NVJE,NVJP,NW,nx,NYNE,NYNP,NYNR,NZNY,
     '  CE,CG,CP,ED,EM,ER,ES,GD,GK,GM,GR,PG,RG,SE,VE,WG,XA,XE,XG,
     '  XP,YG,YP,ZA,ZE,ZG,ZP,DYNAM1,DYNAM2,FIRST,UPDATE_MATRIX,
     '  UPDATE_VECTOR,ERROR,*)

C#### Subroutine: ASSEMBLE3
C###  Description:
C###    ASSEMBLE3 assembles the global unreduced matrices GK, 
C###    GD, GM, etc. for time dependent FEM problems.

C**** If DYNAM1 is true the problem will use the GD matrix
C**** If DYNAM2 is true the problem will use the GM matrix
C**** FIRST indicates whether or not the routine is being called for
C****   the first time. If FIRST is true the (temporary) matrix sparsity
C****   patterns are calculated. If the routine is not being called
C****   for the first time FIRST should be false.
C**** If UPDATE_MATRIX is true the GK (and GD, GM matrices) will be
C****   calculated otherwise just the R.H.S. vector (GR) will be
C****   calculated.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LGE(NHM*NSM,NRCM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NHE(NEM),NHP(NPM),NJE(NEM),
     '  NKE(NKM,NNM,NBFM,NEFM),NKH(NHM,NPM,NCM),NPF(15,NFM),
     '  NPNE(NNM,NBFM,NEFM),NPNODE(0:NP_R_M,0:NRM),
     '  NQE(NSM,NBFM,NEFM),nr,NVHE(NNM,NBFM,NHM,NEFM),NVHP(NHM,NPM,NCM),
     '  NVJE(NNM,NBFM,NJM,NEFM),NVJP(NJM,NPM),NW(NEM,2),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM),NZNY(NYM,*)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CP(NMM,NPM),ED(NHM*NSM,NHM*NSM),
     '  EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),
     '  GD(NZ_GD_M),GK(NZ_GK_M),GM(NZ_GM_M),GR(NYROWM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEFM),
     '  VE(NSM,NKM,NEFM),
     '  WG(NGM,NBM),XA(NAM,NJM,NQM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),YG(NGM,NJM,NEM),YP(NYM,NIYM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      LOGICAL DYNAM1,DYNAM2,FIRST,UPDATE_MATRIX,UPDATE_VECTOR
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nc,ne,nhs1,nhs2,NHST(2),noelem,no_nynr1,no_nynr2,
     '  ny1,ny2,nz,nzy

      CALL ENTERS('ASSEMBLE3',*9999)

      nc=1 !just want the GK matrix variables

      IF(UPDATE_MATRIX) THEN

C***  Initialise matrices for this region

        IF(FIRST) THEN
          nz=0
          NZT(1,nx)=0
          DO no_nynr1=1,NYNR(0,1,nc) !loop over local rows
            ny1=NYNR(no_nynr1,1,nc) !is local row number
            DO no_nynr2=1,NYNR(0,2,nc) !loop over local columns
              ny2=NYNR(no_nynr2,2,nc) !is local column number
              NZNY(ny1,ny2)=0 !temporary sparse storage scheme
            ENDDO !no_nynr2
          ENDDO !no_nynr1
          DO nz=1,NZ_GK_M
            GK(nz)=0.0d0
          ENDDO
          IF(DYNAM1) THEN
            NZT(3,nx)=0
            DO nz=1,NZ_GD_M
              GD(nz)=0.0d0
            ENDDO
          ENDIF
          IF(DYNAM2) THEN
            NZT(4,nx)=0
            DO nz=1,NZ_GM_M
              GM(nz)=0.0d0
            ENDDO
          ENDIF
        ENDIF !first

        DO no_nynr1=1,NYNR(0,1,nc) !loop over local rows
          ny1=NYNR(no_nynr1,1,nc) !is local row number
          GR(ny1)=0.0d0
          DO no_nynr2=1,NYNR(0,2,nc) !loop over local columns
            ny2=NYNR(no_nynr2,2,nc) !is local column number
            nz=NZNY(ny1,ny2)
            GK(nz)=0.0d0
            IF(DYNAM1) THEN
              GD(nz)=0.0d0
            ENDIF
            IF(DYNAM2) THEN
              GM(nz)=0.0d0
            ENDIF
          ENDDO !no_nynr2
        ENDDO !no_nynr1

C***  Find element stiffness matrices

        DO noelem=1,NEELEM(0,nr)    !main element loop
          ne=NEELEM(noelem,nr)
          IF(NW(ne,1).GT.0) THEN
            CALL MELGE(LGE,NBH(1,1,ne),1,ne,NHE(ne),NHST,NKH,
     '        NPNE(1,1,ne),nr,NVHE(1,1,1,ne),nx,NYNE,NYNP,ERROR,*9999)
            IF(IWRIT4(nr,nx).GE.5) THEN
              FORMAT='(/'' Element '',I4,'', Number of variables: '','
     '          //'''NHST(1)='',I3,'', NHST(2)='',I3)'
              WRITE(OP_STRING,FORMAT) ne,NHST(1),NHST(2)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              FORMAT='('' LGE(1..,1): '',14I5,/:(13X,14I5))'
              WRITE(OP_STRING,FORMAT) (LGE(nhs1,1),nhs1=1,NHST(1))
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              FORMAT='('' LGE(1..,2): '',14I5,/:(13X,14I5))'
              WRITE(OP_STRING,FORMAT) (LGE(nhs2,2),nhs2=1,NHST(2))
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF            

C*** Initialize element arrays

            DO nhs1=1,NHST(1)
              ER(nhs1)=0.0d0
              DO nhs2=1,NHST(2)
                ES(nhs1,nhs2)=0.0d0
                IF(DYNAM1) THEN
                  ED(nhs1,nhs2)=0.0d0
                ENDIF
                IF(DYNAM2) THEN
                  EM(nhs1,nhs2)=0.0d0
                ENDIF
              ENDDO
            ENDDO

            IF(ITYP1(nr,nx).EQ.3) THEN !partial differential equation
              CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),
     '          NPNE(1,1,ne),NQE(1,1,ne),nr,NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA,XE,XP,ERROR,*9999)
              CALL XPES30(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),ne,
     '          NHE(ne),NJE(ne),NPNE(1,1,ne),nr,nx,
     '          CE(1,ne),CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '          VE(1,1,ne),WG,XE,XG,YG(1,1,ne),ZE,ZG,
     '          UPDATE_MATRIX,UPDATE_VECTOR,ERROR,*9999)
            ELSE IF(ITYP1(nr,nx).EQ.4) THEN !linear elasticity
              CALL XPES40(NBH(1,1,ne),NBJ(1,ne),
     '          NHE(ne),NJE(ne),NKE(1,1,1,ne),NPF,NPNE(1,1,ne),
     '          NQE(1,1,ne),nr,NVJE(1,1,1,ne),NW(ne,1),nx,
     '          CE(1,ne),CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '          VE(1,1,ne),WG,XA,XE,XG,XP,YG(1,1,ne),
     '          UPDATE_MATRIX,ERROR,*9999)
            ENDIF

            IF(IWRIT4(nr,nx).GE.5) THEN
              WRITE(OP_STRING,
     '          '(/'' Element load vector ER & stiffness matrix ES:'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              DO nhs1=1,NHST(1)
                WRITE(OP_STRING,'('' ER('',I2,'')='',D12.4,'' ES: '','
     '            //'6D12.4,/(25X,6D12.4))') 
     '            nhs1,ER(nhs1),(ES(nhs1,nhs2),nhs2=1,NHST(2))
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDDO !nhs1
              IF(DYNAM1) THEN
                WRITE(OP_STRING,
     '            '(/'' Element damping matrix ED:'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                DO nhs1=1,NHST(1)
                  WRITE(OP_STRING,'('' nhs1=  '',12X,'' ED: '','
     '              //'6D12.4,/(25X,6D12.4))') nhs1, 
     '              (ED(nhs1,nhs2),nhs2=1,NHST(2))
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDDO !nhs1
              ENDIF !dynam1
              IF(DYNAM2) THEN
                WRITE(OP_STRING,
     '            '(/'' Element mass matrix EM:'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                DO nhs1=1,NHST(1)
                  WRITE(OP_STRING,'('' nhs1=  '',12X,'' EM: '','
     '              //'6D12.4,/(25X,6D12.4))') nhs1, 
     '              (EM(nhs1,nhs2),nhs2=1,NHST(2))
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDDO !nhs1
              ENDIF !dynam2
            ENDIF !iwrit

C*** Assemble element matrices into global matrices

            DO nhs1=1,NHST(1)
              ny1=IABS(LGE(nhs1,1))
              GR(ny1)=GR(ny1)+ER(nhs1)
              DO nhs2=1,NHST(2)
                ny2=IABS(LGE(nhs2,2))
                IF(NZNY(ny1,ny2).EQ.0) THEN !increm & record position
                  nz=nz+1                   !of new dof
                  NZT(1,nx)=NZT(1,nx)+1
                  NZNY(ny1,ny2)=nz
                  nzy=nz
                ELSE                        !get existing 1D position
                  nzy=NZNY(ny1,ny2)
                ENDIF
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' nhs1='',I3,'' ny1='',I5,'
     '              //''' nhs2='',I3,'' ny2='',I5,'' nzy='',I7)') 
     '              nhs1,ny1,nhs2,ny2,nzy
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                IF(nzy.LE.NZ_GK_M) THEN
                  GK(nzy)=GK(nzy)+ES(nhs1,nhs2)   !global stiff.n matrix
                ENDIF
                IF(DYNAM1) THEN                                
                  IF(nzy.LE.NZ_GD_M) THEN
                    GD(nzy)=GD(nzy)+ED(nhs1,nhs2) !global damping matrix
                  ENDIF
                ENDIF
                IF(DYNAM2) THEN                                
                  IF(nzy.LE.NZ_GM_M) THEN
                    GM(nzy)=GM(nzy)+EM(nhs1,nhs2) !global mass matrix
                  ENDIF
                ENDIF
              ENDDO !nhs2
            ENDDO !nhs1
          ENDIF !nw
        ENDDO !noelem (ne)
        IF(nz.GT.NZ_GK_M) THEN
          ERROR='>>Increase NZ_GK_MX'
          GOTO 9999
        ENDIF
        IF(DYNAM1) THEN
      	  NZT(3,nx)=NZT(1,nx)
          IF(nz.GT.NZ_GD_M) THEN
            ERROR='>>Increase NZ_GD_MX'
            GOTO 9999
          ENDIF
        ENDIF
        IF(DYNAM2) THEN
      	  NZT(4,nx)=NZT(1,nx)
          IF(nz.GT.NZ_GM_M) THEN
            ERROR='>>Increase NZ_GM_MX'
            GOTO 9999
          ENDIF
        ENDIF
	
      ELSE  !just calculate RHS vector
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(NW(ne,1).GT.0) THEN

            CALL MELGE(LGE,NBH(1,1,ne),nc,ne,NHE(ne),NHST,NKH,
     '        NPNE(1,1,ne),nr,NVHE(1,1,1,ne),nx,NYNE,NYNP,ERROR,*9999)

            DO nhs1=1,NHST(1)        !initialize element arrays
              ER(nhs1)=0.0d0
            ENDDO

            CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),NQE(1,1,ne),nr,NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA,XE,XP,ERROR,*9999)
            CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '        nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
            CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '        SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)

            IF(ITYP1(nr,nx).EQ.3) THEN !partial differential equation
              CALL XPES30(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),ne,
     '          NHE(ne),NJE(ne),NPNE(1,1,ne),nr,nx,
     '          CE(1,ne),CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '          VE(1,1,ne),WG,XE,XG,YG(1,1,ne),ZE,ZG,
     '          UPDATE_MATRIX,UPDATE_VECTOR,ERROR,*9999)

            ELSE IF(ITYP1(nr,nx).EQ.4) THEN !linear elasticity
              CALL XPES40(NBH(1,1,ne),NBJ(1,ne),
     '          NHE(ne),NJE(ne),NKE(1,1,1,ne),NPF,NPNE(1,1,ne),
     '          NQE(1,1,ne),nr,NVJE(1,1,1,ne),NW(ne,1),nx,
     '          CE(1,ne),CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '          VE(1,1,ne),WG,XA,XE,XG,XP,YG(1,1,ne),
     '          UPDATE_MATRIX,ERROR,*9999)
            ENDIF

            IF(IWRIT4(nr,nx).GE.5) THEN
              WRITE(OP_STRING,
     '          '(/'' Element load vector ER:'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              DO nhs1=1,NHST(1)
                WRITE(OP_STRING,'('' ER('',I2,'')='',D12.4)') nhs1,
     '            ER(nhs1)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDDO !nhs1
            ENDIF

            DO nhs1=1,NHST(1)
              ny1=IABS(LGE(nhs1,1))
              GR(ny1)=GR(ny1)+ER(nhs1)
            ENDDO

          ENDIF !nw(ne,1)>0
        ENDDO !noelem
      ENDIF !update matrix

      IF(IWRIT4(nr,nx).GE.4) THEN
        IF(UPDATE_MATRIX) THEN
          WRITE(OP_STRING,'(/'' Global load vector GR:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NYNR(0,1,1)='',I5)') NYNR(0,1,1)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' GR:'',8(X,D12.4),/(4X,8(X,D12.4)))')
     '      (GR(NYNR(no_nynr1,1,1)),no_nynr1=1,NYNR(0,1,1))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(/'' Global stiffness matrix GK:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NYNR(0,2,1)='',I5)') NYNR(0,2,1)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DO no_nynr1=1,NYNR(0,1,1)
            ny1=NYNR(no_nynr1,1,1)
            WRITE(OP_STRING,'('' ny1='',I5,'' NZNY:'',16I7,'
     '        //'/:(16X,16I7))') ny1,(NZNY(ny1,NYNR(no_nynr2,2,1)),
     '        no_nynr2=1,NYNR(0,2,1)) 
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDDO
          WRITE(OP_STRING,'('' GK:'',8(X,D12.4),/(4X,8(X,D12.4)))')
     '      (GK(nz),nz=1,NZT(1,nx))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          IF(DYNAM1) THEN
            WRITE(OP_STRING,
     '	      '(/'' Global damping matrix GD:'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' GD:'',8(X,D12.4),/(4X,8(X,D12.4)))')
     '        (GD(nz),nz=1,NZT(3,nx))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF !dynam1
          IF(DYNAM2) THEN
            WRITE(OP_STRING,
     '	      '(/'' Global mass matrix GM:'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' GM:'',8(X,D12.4),/(4X,8(X,D12.4)))')
     '        (GM(nz),nz=1,NZT(4,nx))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE !just want R.H.S. vector
          WRITE(OP_STRING,'(/'' Global load vector GR:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NYNR(0,1,1)='',I5)') NYNR(0,1,1)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' GR:'',8(X,D12.4),/(4X,8(X,D12.4)))')
     '      (GR(NYNR(no_nynr1,1,1)),no_nynr1=1,NYNR(0,1,1))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF    

      CALL EXITS('ASSEMBLE3')
      RETURN
 9999 CALL ERRORS('ASSEMBLE3',ERROR)
      CALL EXITS('ASSEMBLE3')
      RETURN 1
      END


      SUBROUTINE ASSEMBLE8(ERROR,*)

C#### Subroutine: ASSEMBLE8
C###  Description:
C###    ASSEMBLE8 assembles the global unreduced matrices GK, 
C###    GM, etc. for Laplace transform problems.

      IMPLICIT NONE
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('ASSEMBLE8',*9999)

      CALL EXITS('ASSEMBLE8')
      RETURN
 9999 CALL ERRORS('ASSEMBLE8',ERROR)
      CALL EXITS('ASSEMBLE8')
      RETURN 1
      END


 
      SUBROUTINE FUNCT1(IFLAG,N,XC,FC,GC,IWK,LIW,W,LW)

C**** Evaluates objective function FC for NAG optimization routine E04JBF.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:data00.cmn'
      INCLUDE 'cmiss$reference:fgbez00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:opti00.cmn'
      INCLUDE 'cmiss$reference:xif00.cmn'
!     Parameter List
      INTEGER IFLAG,IWK(*),LIW,LW,N
      REAL*8 FC,GC(*),W(*),XC(*)
!     Local Variables
      INTEGER I,IBEG,IEND,KK,KK1,NOCO,NOOPTI,NP,NTCOD,
     '  NTCOQUD(16),INTWORK(1)
      REAL*8 ACN,B(4),BCN,CL20,CL21,CL22,F01,F02,F11,F12,S01,S02,S1,
     '  S11,S12,S2,VL20,VL21,VL22,W1,XDIFFERENCE,XI1,XI2,XKA,XKB,XKSIB,
     '  XNPB,YNPB,REALWORK(1)
      CHARACTER COD(16)*20,COQUD(16,1)*1,ERROR*(MXCH),STRING*(MXCH)
      LOGICAL END

      CALL ENTERS('FUNCT1',*9999)
      
      WRITE(OP_STRING,'('' FUNCT1 NOT IN OPERATION'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      
C CPB 20/7/93 SUBROUTINE NOT NEEDED ANYMORE
C 
C      NOCO=1
C
C      IF(KTYP26.EQ.1.AND.KTYP27.EQ.2) THEN
CC **    Optimize material params by minimizing sum of squared reaction diffs
C        DO NOCO=1,6
C          NTCOQUD(NOCO)=0
C        ENDDO
C        DO NOOPTI=1,NTOPTI
C          PAOPTI(NOOPTI)=XC(NOOPTI)
C        ENDDO
C        WRITE(OP_STRING,'('' FUNCT1 XC: '',10E12.3)')
C     '	  (PAOPTI(NOOPTI),NOOPTI=1,
C     '    NTOPTI)
C        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C
C        IF(KTYP28.EQ.0) THEN
C          COD(1)='FEM'
C          COD(2)='SOLVE'
C          COD(3)='STEP'
C          COD(4)='1'
C          COD(5)='UPDATE'
C          COD(6)='1'
C          NTCOD=6
C          NOCO=1
C          CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,
C     '      END,STRING,INTWORK,REALWORK,ERROR,*9999)
C        ENDIF
C
C        COD(1)='FEM'
C        COD(2)='EVALUATE'
C        COD(3)='OBJECTIVE'
C        NTCOD=3
C        NOCO=1
C        CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,END,
C     '    STRING,INTWORK,REALWORK,ERROR,*9999)
C        FC=FUNC
C
C      ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.1) THEN
CC **    Objective function is minimum area of trapezoids
C        DO NOCO=1,6
C          NTCOQUD(NOCO)=0
C        ENDDO
C        DO NOOPTI=1,NTOPTI
C          PAOPTI(NOOPTI)=XC(NOOPTI)
C        ENDDO
C        WRITE(OP_STRING,'('' FUNCT1 XC: '',10E12.3)')
C     '	  (PAOPTI(NOOPTI),NOOPTI=1,
C     '    NTOPTI)
C       CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C
C        COD(1)='FEM'
C        COD(2)='DEFINE'
C        COD(3)='POLYLINE'
C        COQUD(3,1)='c'
C        COD(4)='WITH'
C        COD(5)='OPT'
C       COD(6)='ENDPOINT'
C       COD(7)='30'
C       NTCOD=7
C       NTCOQUD(3)=1
C       NOCO=1
C       CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,END,
C    '    STRING,INTWORK,REALWORK,ERROR,*9999)
C
C       COD(1)='FEM'
C       COD(2)='DEFINE'
C       COD(3)='FIT'
C       COQUD(3,1)='r'
C       COD(4)='GEOM'
C       NTCOD=4
C       NTCOQUD(3)=1
C       NOCO=1
C       CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,END,
C    '    STRING,INTWORK,REALWORK,ERROR,*9999)
C
C       COD(1)='FEM'
C       COD(2)='FIT'
C       COD(3)='GEOM'
C       NTCOD=3
C       NOCO=1
C       CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,END,
C    '    STRING,INTWORK,REALWORK,ERROR,*9999)
C
C       COD(1)='FEM'
C       COD(2)='UPDATE'
C       COD(3)='NODES'
C       NTCOD=3
C       NOCO=1
C       CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,END,
C    '    STRING,INTWORK,REALWORK,ERROR,*9999)
C
C       COD(1)='FEM'
C       COD(2)='LIST'
C       COD(3)='DATA'
C       COD(4)='ERROR'
C       NTCOD=4
C       NOCO=1
C       CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,END,
C    '    STRING,INTWORK,REALWORK,ERROR,*9999)
C
C       FUNC=SQED
C       FC=FUNC
C
C       IF(UPVUOP) THEN
C         COD(1)='FEM'
C         COD(2)='DEFINE'
C         COD(3)='LINE'
C         COQUD(3,1)='s'
C         NTCOD=3
C         NTCOQUD(3)=1
C         NOCO=1
C         CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,
C    '      END,STRING,INTWORK,REALWORK,ERROR,*9999)
C
C         COD(1)='FEM'
C         COD(2)='DEFINE'
C         COD(3)='DATA'
C         COQUD(3,1)='s'
C         NTCOD=3
C         NTCOQUD(3)=1
C         NOCO=1
C         CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,
C    '      END,STRING,INTWORK,REALWORK,ERROR,*9999)
C
C         COD(1)='FEM'
C         COD(2)='DEFINE'
C         COD(3)='DATA'
C         COQUD(3,1)='s'
C         COD(4)='PROJECTIONS'
C         NTCOD=4
C         NTCOQUD(3)=1
C         NOCO=1
C          CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,
C    '      END,STRING,INTWORK,REALWORK,ERROR,*9999)
C       ENDIF
C
C     ELSE IF(KTYP26.LE.2.AND.KTYP27.EQ.4) THEN
C **    Objective function is hydrostatic pressure condition
C       DO NOCO=1,6
C         NTCOQUD(NOCO)=0
C       ENDDO
C       DO NOOPTI=1,NTOPTI
C         PAOPTI(NOOPTI)=XC(NOOPTI)
C       ENDDO
C       WRITE(OP_STRING,'('' FUNCT1 XC: '',10E12.3)')
C    '	  (PAOPTI(NOOPTI),NOOPTI=1,
C    '    NTOPTI)
C       CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C
C       COD(1)='FEM'
C       COD(2)='SOLVE'
C       NTCOD=2
C       NOCO=1
C       CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,END,
C    '    STRING,INTWORK,REALWORK,ERROR,*9999)
C
C       COD(1)='FEM'
C       COD(2)='EVALUATE'
C       COD(3)='OBJECTIVE'
C
C       NOCO=1
C       CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,END,
C    '    STRING,INTWORK,REALWORK,ERROR,*9999)
C       FC=FUNC
C
C       IF(UPVUOP) THEN
C         COD(1)='FEM'
C         COD(2)='UPDATE'
C         COD(3)='MESH'
C         NTCOD=3
C         NOCO=1
C         CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,
C    '      END,STRING,INTWORK,REALWORK,ERROR,*9999)
C       ENDIF
C
C     ELSE IF(KTYP26.EQ.3.AND.KTYP27.EQ.1) THEN !Dave's stripe calcs.
C
C       np=2
C        xksia=xc(1)
C        xksib=xc(2)
C        if (xksia.lt.1.0d-6) xksia=1.0d-6
C        if (xksia.gt.0.999999d0) xksia=0.999999d0     !Ask Peter how
C        b(1)=(1.0d0-xksia)**3                            !to deal with this
C        b(2)=3.0d0*xksia*(1.-xksia)**2                   !underflow problem
C        b(3)=3.0d0*xksia*xksia*(1.-xksia)
C        b(4)=xksia*xksia*xksia
C        xxa=0.0d0
C        yya=0.0d0
C        do i=1,4
C          xxa=xxa+b(i)*xa(i)
C          yya=yya+b(i)*ya(i)
C        end do   !end i
C
C        if (xksib.lt.1.0d-6) xksib=1.0d-6
C        if (xksib.gt.0.999999d0) xksib=0.999999d0
C        b(1)=(1.0d0-xksib)**3
C        b(2)=3.0d0*xksib*(1.0d0-xksib)**2
C       b(3)=3.0d0*xksib*xksib*(1.0d0-xksib)
C       b(4)=xksib*xksib*xksib
C       xxb=0.0d0
C       yyb=0.0d0
C       do i=1,4
C         xxb=xxb+b(i)*xb(i)
C         yyb=yyb+b(i)*yb(i)
C       end do    !end i
C
C       fc=((yya-yyb)**2+(xxa-xxb)**2)**0.5d0
C
C     ELSE IF(KTYP26.EQ.4.AND.KTYP27.EQ.1) THEN !Dave's stripe calcs.
C       np=2
C ...   xksia passed in with common block
C       xksib=xc(1)
C       acdiff=0.0d0
C       acn=0.0d0
C       bcn=0.0d0
C       w1=0.0d0
C       xka=0.0d0
C       do kk=1,nsimp+1      !Simpson's Rule
C         xka=xksia+(xksib-xksia)*(kk-1.0d0)/dble(nsimp)
C
C         acn=xa(1)*(-3.0d0*(1.0d0-xka)**2)
C         acn=acn+xa(2)*(3.0d0*(1.0d0-xka)*(1.0d0-3.0d0*xka))
C         acn=acn+xa(3)*(3.0d0*xka*(2.0d0-3.0d0*xka))
C         acn=acn+xa(4)*3.0d0*xka*xka
C
C         bcn=ya(1)*(-3.0d0*(1.-xka)**2)
C         bcn=bcn+ya(2)*(3.0d0*(1.-xka)*(1.0d0-3.0d0*xka))
C         bcn=bcn+ya(3)*(3.0d0*xka*(2.0d0-3.0d0*xka))
C         bcn=bcn+ya(4)*3.0d0*xka*xka
C
C         w1=2.0d0                 !presume kk is odd
C         kk1=kk/2
C         kk1=kk1*2
C         if (kk1.eq.kk) w1=4.   !kk is even
C         if (kk.eq.1.or.kk.eq.nsimp+1) w1=1.    !kk is endpt.
C         acdiff=acdiff+((acn*acn+bcn*bcn)**0.5d0)*w1*(xksib-xksia)
C    '           /dble(nsimp)/3.0d0
C       end do !end kk
C
C       fc=1.0d0+(slen-acdiff)**2     !exact length - estimated length
C
C 
C     ELSE IF(KTYP26.EQ.4.AND.KTYP27.EQ.2) THEN !xi calculation in
C                                               !cubic-linear element
C ...   Find x_endo         !x coordinate of endocardial point
C                           !corresponding to trial value of xi1
C       XI1=XC(1)
C       S01=1.0D0-3.0D0*XI1*XI1+2.0D0*XI1*XI1*XI1  !evaluate hermite basis function
C       S02=3.0D0*XI1*XI1-2.0D0*XI1*XI1*XI1
C       S11=XI1-2.0D0*XI1*XI1+XI1*XI1*XI1
C       S12=-(XI1*XI1-XI1*XI1*XI1)
C       F01=XNODVAL(1,1)                  !Nodal values
C       F02=XNODVAL(3,1)
C       F11=XNODVAL(2,1)
C       F12=XNODVAL(4,1)
C       X_ENDO=S01*F01+S02*F02+S11*F11+S12*F12
C
C ...   Find y_endo         !y coordinate of endocardial point
C                           !corresponding to trial value of xi1
C       F01=XNODVAL(1,2)
C       F02=XNODVAL(3,2)
C       F11=XNODVAL(2,2)
C       F12=XNODVAL(4,2)
C       Y_ENDO=S01*F01+S02*F02+S11*F11+S12*F12
C
C ...   Find x_epi          !x coordinate of epicardial point
C                           !corresponding to trial value of xi1
C       F01=XNODVAL(5,1)
C       F02=XNODVAL(7,1)
C       F11=XNODVAL(6,1)
C       F12=XNODVAL(8,1)
C       X_EPI=S01*F01+S02*F02+S11*F11+S12*F12
C
C ...   Find y_epi          !y coordinate of epicardial point
C                           !corresponding to trial value of xi1
C       F01=XNODVAL(5,2)
C       F02=XNODVAL(7,2)
C       F11=XNODVAL(6,2)
C       F12=XNODVAL(8,2)
C       Y_EPI=S01*F01+S02*F02+S11*F11+S12*F12
C
C ...   Calculate (XNPB,YNPB) the nearest point to (XPXIF,YPXIF) on the
C       line drawn between (X_ENDO,Y_ENDO) and (X_EPI,Y_EPI)
C       XDIFFERENCE=DABS(X_EPI-X_ENDO)
C       IF(XDIFFERENCE.LT.1.0D-8) THEN      !Avoid singular case
C         XNPB=X_EPI
C         YNPB=YPXIF
C       ELSE                                !Nonsingular case
C         S1=X_ENDO-X_EPI
C         S2=Y_ENDO-Y_EPI
C         XNPB=(XPXIF*S1*S1+X_ENDO*S2*S2-S1*S2*(Y_ENDO-YPXIF))
C    '          /(S1*S1+S2*S2)
C         YNPB=Y_ENDO+S2*(XNPB-X_ENDO)/S1
C       ENDIF
C
C       FC=(XNPB-XPXIF)**2+(YNPB-YPXIF)**2      !value of obj. function
C
C
C     ELSE IF(KTYP26.EQ.4.AND.KTYP27.EQ.3) THEN !xi calculation in
C                                               !Lagrange quadratic element
C       XI1=XC(1)
C       XI2=XC(2)
C
C       VL20=2.0D0*(XI1-0.5D0)*(XI1-1.0D0)     !shape factors
C       VL21=-4.0D0*(XI1-1.0D0)*XI1
C       VL22=2.0D0*XI1*(XI1-0.5D0)
C       CL20=2.0D0*(XI2-0.5D0)*(XI2-1.0D0)
C       CL21=-4.0D0*(XI2-1.0D0)*XI2
C       CL22=2.0D0*XI2*(XI2-0.5D0)
C
C ...   Pass array XNODVAL containing nodal values via common
C
C       XNPB=VL20*CL20*XNODVAL(1,1)+VL21*CL20*XNODVAL(2,1)+
C    '       VL22*CL20*XNODVAL(3,1)+VL20*CL21*XNODVAL(4,1)+
C    '       VL21*CL21*XNODVAL(5,1)+VL22*CL21*XNODVAL(6,1)+
C    '       VL20*CL22*XNODVAL(7,1)+VL21*CL22*XNODVAL(8,1)+
C    '       VL22*CL22*XNODVAL(9,1)
C
C       YNPB=VL20*CL20*XNODVAL(1,2)+VL21*CL20*XNODVAL(2,2)+
C    '       VL22*CL20*XNODVAL(3,2)+VL20*CL21*XNODVAL(4,2)+
C    '       VL21*CL21*XNODVAL(5,2)+VL22*CL21*XNODVAL(6,2)+
C    '       VL20*CL22*XNODVAL(7,2)+VL21*CL22*XNODVAL(8,2)+
C    '       VL22*CL22*XNODVAL(9,2)
C
C       FC=(XNPB-XPXIF)**2+(YNPB-YPXIF)**2      !value of obj. function
C
C
C     ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.6) THEN
C **    Optimize geometric params by minimizing fluid interface residual
C       DO NOCO=1,6
C         NTCOQUD(NOCO)=0
C       ENDDO
C       DO NOOPTI=1,NTOPTI
C         PAOPTI(NOOPTI)=XC(NOOPTI)
C       ENDDO
C       IF(DOP) THEN
C         WRITE(OP_STRING,
C    '    '('' FUNCT1 XC: '',10E12.3)') (PAOPTI(NOOPTI),NOOPTI=1,NTOPTI)
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C       ENDIF
C
C       COD(1)='FEM'
C       COD(2)='UPDATE'
C       COD(3)='NODE'
C       COD(4)='INTERFACE'
C       COD(5)='FLUX'
C       NTCOD=5
C       NOCO=1
C       IF(DOP)THEN
C         WRITE(OP_STRING,'('' fem update node interface flux..'')')
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C       ENDIF
C       CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,END,
C    '    STRING,INTWORK,REALWORK,ERROR,*9999)
C
C       COD(1)='FEM'
C       COD(2)='SOLVE'
C       COD(3)='FOR'
C       COD(4)='1'
C       NTCOD=4
C       NOCO=1
C       IF(DOP) THEN
C         WRITE(OP_STRING,'('' fem solve for 1..'')')
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C       ENDIF  
C       CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,END,
C    '    STRING,INTWORK,REALWORK,ERROR,*9999)
C
C       COD(1)='FEM'
C       COD(2)='SOLVE'
C       COD(3)='FOR'
C       COD(4)='2'
C       NTCOD=4
C       NOCO=1
C       IF(DOP) THEN
C         WRITE(OP_STRING,'('' fem solve for 2..'')')
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C       ENDIF
C       CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,END,
C    '    STRING,INTWORK,REALWORK,ERROR,*9999)
C
C       COD(1)='FEM'
C       COD(2)='EVALUATE'
C       COD(3)='OBJECTIVE'
C       NTCOD=3
C       NOCO=1
C       IF(DOP) THEN
C         WRITE(OP_STRING,'('' fem evaluate objective..'')')
C         CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C       ENDIF
C       CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,END,
C    '    STRING,INTWORK,REALWORK,ERROR,*9999)
C       FC=FUNC
C
C       IF(UPVUOP) THEN
C         COD(1)='FEM'
C         COD(2)='UPDATE'
C         COD(3)='NODE'
C         COD(4)='INTERFACE'
C         COD(5)='INCREMENT'
C         NTCOD=5
C         NOCO=1
C         WRITE(OP_STRING,'('' fem update node interface position..'')')
C         CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C         CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,
C    '      END,STRING,INTWORK,REALWORK,ERROR,*9999)
C
C         COD(1)='FEM'
C         COD(2)='DEFINE'
C         COD(3)='LINE'
C         COQUD(3,1)='s'
C         NTCOD=3
C         NOCO=1
C         NTCOQUD(3)=1
C         CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,
C    '      END,STRING,INTWORK,REALWORK,ERROR,*9999)
C       ENDIF
C
C     ENDIF

      CALL EXITS('FUNCT1')
      RETURN
 9999 CALL ERRORS('FUNCT1',ERROR)
      CALL TRIM(ERROR,IBEG,IEND)
      WRITE(OP_STRING,'(A)') ' '//ERROR(IBEG:IEND)
      CALL WRITES(IOER,OP_STRING,ERROR,*9998)
      ERROR(1:)=' '
      CALL EXITS('FUNCT1')
 9998 RETURN
      END


      SUBROUTINE FUNCT3(N,X,FVEC,IFLAG)

C**** Evaluates residuals FVEC for NAG routine C05NCF from variables X.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp100.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:opti00.cmn'
!     Parameter List
      INTEGER IFLAG,INTWORK(1),N
      REAL*8 X(*),FVEC(*),REALWORK(1)
!     Local Variables
      INTEGER IBEG,IEND,NOCO,NOOPTI,NTCOD,NTCOQUD(16)
      CHARACTER COD(16)*20,COQUD(16,1)*1,ERROR*(MXCH),STRING*(MXCH)
      LOGICAL END

      CALL ENTERS('FUNCT3',*9999)

      WRITE(OP_STRING,'('' FUNCT1 NOT IN OPERATION'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      
C CPB 20/7/93 SUBROUTINE NOT NEEDED ANYMORE
C     NOCO=1
C     IF(KTYP27.EQ.3) THEN
C       DO NOCO=1,6
C         NTCOQUD(NOCO)=0
C       ENDDO
C       DO NOOPTI=1,NTOPTI
C         PAOPTI(NOOPTI)=X(NOOPTI)
C       ENDDO
C       WRITE(OP_STRING,*)(' X(',NOOPTI,')=',X(NOOPTI),NOOPTI=1,NTOPTI)
C       CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C
C       COD(1)='FEM'
C       COD(2)='SOLVE'
C       COD(3)='FOR'
C       COD(4)='1'
C       NTCOD=4
C       NOCO=1
C       CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,END,
C    '    STRING,INTWORK,REALWORK,ERROR,*9999)
C       COD(1)='FEM'
C       COD(2)='SOLVE'
C       COD(3)='FOR'
C       COD(4)='2'
C       NTCOD=4
C       NOCO=1
C       CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,END,
C    '    STRING,INTWORK,REALWORK,ERROR,*9999)
C       COD(1)='FEM'
C       COD(2)='EVALUATE'
C       COD(3)='OBJECTIVE'
C       COD(4)='FOR'
C       COD(5)='2' !Always put BE region number here.
C       NTCOD=5
C       NOCO=1
C       CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,END,
C    '    STRING,INTWORK,REALWORK,ERROR,*9999)
C       DO NOOPTI=1,NTOPTI
C CPB 20/7/93 NOT NEEDED NOW THAT THE ARRAYS ARE NOT IN THE COMMON BLOCK
C          FVEC(NOOPTI)=RESID(NOOPTI)
C         WRITE(OP_STRING,*)' RESID(',NOOPTI,')=',RESID(NOOPTI)
C         CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C       ENDDO
C       IF(KTYP100.EQ.2) THEN
C         WRITE(OP_STRING,*)' THETA_SAT=',THETA_SAT
C         CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C       ENDIF
C       IF(UPVUOP) THEN
C         COD(1)='FEM'
C         COD(2)='UPDATE'
C         COD(3)='MESH'
C         NTCOD=3
C         NOCO=1
C         CALL FEM(ISEG_TEMP,NOCO,NTCOD,NTCOQUD,COD,COQUD,CSEG_TEMP,
C    '      END,STRING,INTWORK,REALWORK,ERROR,*9999)
C       ENDIF
C     ENDIF

      CALL EXITS('FUNCT3')
      RETURN
 9999 CALL ERRORS('FUNCT3',ERROR)
      CALL TRIM(ERROR,IBEG,IEND)
      WRITE(OP_STRING,'(A)') ' '//ERROR(IBEG:IEND)
      CALL WRITES(IOER,OP_STRING,ERROR,*9998)
      ERROR(1:)=' '
      CALL EXITS('FUNCT3')
 9998 RETURN
      END


      SUBROUTINE GEOMIN(INTWORK,REALWORK,ERROR,*)

C#### Subroutine: GEOMIN
C###  Description:
C###    GEOMIN performs nonlinear geometric fitting.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:data00.cmn'
      INCLUDE 'cmiss$reference:fit001.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER INTWORK(*)
      REAL*8 REALWORK(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER MXSG
      PARAMETER (MXSG=2000)
      INTEGER ISEG(MXSG),KOUNT
      REAL*8 SQEDOLD
      CHARACTER CSEG(MXSG)*20,STRING*(MXCH)
      LOGICAL CONTINUE,END

      CALL ENTERS('GEOMIN',*9999)
      FIRSTA=.TRUE.
      IF(KTYP27.EQ.5) THEN
        KOUNT=0
        SQED=0.0D0
        CONTINUE=.TRUE.
        DO WHILE(CONTINUE)
          KOUNT=KOUNT+1
          SQEDOLD=SQED
          CO(1)='FEM'
          CO(2)='READ'
          CO(3)='GEOMIN'
          COQU(3,1)='COM'
          COQU(3,2)='DOC'
          NTCO=3
          NTCOQU(1)=0
          NTCOQU(2)=0
          NTCOQU(3)=2
          noco=1
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)
          IF(DABS(SQED-SQEDOLD).GT.0.02D0*SQEDOLD.AND.KOUNT.LE.10) THEN
            CONTINUE=.TRUE.
          ELSE
            CONTINUE=.FALSE.
          ENDIF
        ENDDO
      ENDIF

 9998 CALL EXITS('GEOMIN')
      RETURN
 9999 CALL ERRORS('GEOMIN',ERROR)
      CALL EXITS('GEOMIN')
      RETURN 1
      END




C CPB 29/3/96 This is the old assemble3 before I rewrote it.
      SUBROUTINE MARCH1(IBT,IDISP,IDO,INP,ISC_GKK,ISC2_GKK,ISR_GKK,
     '  ISR2_GKK,IWK1,IWK2,LGE,NBH,NBJ,NEELEM,NHE,NHP,NJE,NKE,NKH,NKJ,
     '  NONY,NPF,NP_INTERFACE,NPNE,NPNODE,NPNY,NQE,nr,
     '  NVHE,NVHP,NVJE,NVJP,
     '  NW,nx,NYNE,NYNO,NYNP,NYNR,NZNY,CE,CG,CONY,CP,
     '  CYNO,ED,EM,ER,ES,GD,GK,GM,GR,PG,RG,SE,VE,WG,
     '  XA,XE,XG,XO,XP,YG,YP,ZA,ZE,ZG,ZP,GKK,GRR,WK1,
     '  FIX,ERROR,*)

C#### Subroutine: MARCH1
C###  Description:
C###    MARCH1 performs time integration of linear or nonlinear 
C###    (in terms other than transient) probs (+ nonlinear b.c.s) by 
C###    linear,quadratic or cubic algorithm.

C**** (KTYP22=1,2 or 3) with fixed or automatic time steps (KTYP23=1,2).
C**** Error is calculated & time step adjusted if KTYP23=2.
C**** For 2nd order problems initial accel.ns are calculated from d.e.
C**** with  known initial velocities and displacements (cubic only).
C**** YP(ny,1) is solution vector at time T+DT (returned from SOLVE3).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b08.cmn'
      INCLUDE 'cmiss$reference:b10.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:head00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:marc00.cmn'
      INCLUDE 'cmiss$reference:moti00.cmn'
      INCLUDE 'cmiss$reference:outp00.cmn'
      INCLUDE 'cmiss$reference:time01.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDISP(10),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISC_GKK(NISC_GKKM),ISC2_GKK(NISC2_GKKM),
     '  ISR_GKK(NISR_GKKM),ISR2_GKK(NISR2_GKKM),IWK1(5*NOM),IWK2(8*NOM),
     '  LGE(NHM*NSM,NRCM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NHE(NEM),NHP(NPM),NJE(NEM),
     '  NKE(NKM,NNM,NBFM,NEFM),NKH(NHM,NPM,NCM),NKJ(NJM,NPM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM),NPF(15,NFM),
     '  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEFM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),NQE(NSM,NBFM,NEFM),
     '  nr,NVHE(NNM,NBFM,NHM,NEFM),NVHP(NHM,NPM,NCM),
     '  NVJE(NNM,NBFM,NJM,NEFM),NVJP(NJM,NPM),NW(NEM,2),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),NZNY(NZ_NZNY_M)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CONY(0:NOYM,NYM,NRCM,0:NRM),
     '  CP(NMM,NPM),CYNO(0:NYOM,NOOPM,NRCM,0:NRM),ED(NHM*NSM,NHM*NSM),
     '  EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),
     '  GD(NZ_GD_M),GK(NZ_GK_M),GM(NZ_GM_M),GR(NYROWM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEFM),
     '  VE(NSM,NKM,NEFM),WG(NGM,NBM),XA(NAM,NJM,NQM),XE(NSM,NJM),
     '  XG(NJM,NUM),XO(NOM),XP(NKM,NVM,NJM,NPM),YG(NGM,NJM,NEM),
     '  YP(NYM,NIYM),ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      REAL*8 GKK(NZ_GKK_M),GRR(NOM),WK1(4*NOM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER i,IBEG,IBEG1,IBEG5,IEND,IEND1,IEND5,INCR,IP,
     '  k,KFAC,KOUNT,n,nb,nc,ne,nh,nhx,
     '  no_coeffs,noelem,nonode,no_unit,np,NP1,NP2,
     '  NSTEP,NT_COEFFS,NT_UNIT,NWRIT,ny,NY_UNIT(20),
     '  no_nynr,no_nynr2,VERSION
      REAL*8 BDRY_CONC,C1,C2,Conc_N2,d_Conc_N2_dt,d_Conc_N2_dVol,
     '  ERR,FLOW,Flow_N2_integral,Peak_flow,PF1,
     '  Prev_Conc_N2,SNORM,SUM,T,T_cycle,T_end_breath_hold,
     '  T_end_expiration,T_end_inspiration,T_expire,Tidal_volume,
     '  TOL,VOL_CH
      CHARACTER CHAR*5,FILE*100
      LOGICAL CONVERGED,DYNAM1,DYNAM2,FILEIP,FIRST,
     '  FLOW_CALC,LINEAR,LUNG_PRINT,OUTPUT,UPDATE_MATRIX,UPDATE_VECTOR
      SAVE NSTEP

      CALL ENTERS('MARCH1',*9999)

c     timea=secnds(0.0D0)

      IP=1
      NWRIT=1
      INCR=0
      DT=TINCR
      IF(ITYP6(nr,nx).EQ.1) THEN      !linear equations
        LINEAR=.TRUE.
      ELSE IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear equations or bcs
        LINEAR=.FALSE.
      ENDIF
      FILEIP=.FALSE.
      FIRST=.TRUE.
      FLOW_CALC=.FALSE.
      UPDATE_MATRIX=.TRUE.
      UPDATE_VECTOR=.TRUE.

C*** Determine whether or not the problem needs the GD and GM matrices.
C*** If the problem needs GD DYNAM1 should be set the true, false
C*** otherwise. If the problem needs GM DYNAM2 should be set to true
C*** false otherwise.

C!!!cpb 9/12/94 Only Navier-Stokes covered at the moment.

      IF(ITYP2(nr,nx).EQ.3.OR.ITYP2(nr,nx).EQ.5) THEN !Advec-Diffusion
        DYNAM1=.TRUE.  !Use GD                        !or Navier-Stokes
        DYNAM2=.FALSE. !Don't use GM
      ELSE
        DYNAM1=.TRUE.  !Use GD
        DYNAM2=.TRUE.  !Use GM
      ENDIF

      IF(RESTART) THEN
        T=T_RESTART
        FIRST=.FALSE.
        UPDATE_MATRIX=.FALSE.
      ELSE IF(.NOT.RESTART) THEN !perform initial tasks
        T=TSTART
        NSTEP=0
        
        CALL IOHIST(IUNIT,NIYLIST,NPNY,NRLIST,NUMTIMEDATA,NYNR,
     '    TOTALNQ,TOTALNY,TIME,YP,YPMAX,YPMIN,COMMAND,FILEFORMAT,
     '    FILENAME,'OPEN',.TRUE.,.FALSE.,ERROR,*9999)


C*** Put initial conditions into current solUTION YP(ny,1) & previous
C*** solution YP(ny,8)
        DO no_nynr=1,NYNR(0,0,1,nr)
          ny=NYNR(no_nynr,0,1,nr)
          IF(.NOT.FIX(ny,1)) YP(ny,1)=YP(ny,3)
          YP(ny,8)=YP(ny,1)
        ENDDO

        IF(ITYP2(nr,nx).EQ.5.AND.ITYP3(nr,nx).EQ.2) THEN !Lung model
          !Write headings for acinus file
          CALL TRIM(FILE02,IBEG,IEND)
          CALL OPENF(IOFILE4,'DISK',FILE02(IBEG:IEND)//'.acinus','NEW',
     '      'SEQUEN','FORMATTED',160,ERROR,*9999)
          IOFI = IOFILE4
          IF(CALL_MOTI) THEN
            CALL OPMOTI(IBT,NBH,NEELEM,NHP,NKH,
     '        NPNODE,nr,FIX,YP,ERROR,*9999)
          ENDIF
          Bdry_Conc=YP(1,1) !to use later when YP(1,1) has been reset
        ENDIF

      ENDIF !.not.restart

C***  Main time loop

 10   CONTINUE !loop start

        NSTEP=NSTEP+1

        IF(ITYP2(nr,nx).EQ.5.AND.ITYP3(nr,nx).EQ.2) THEN !Lung model
          IF(KTYP58(nr).EQ.3) THEN !update flow using Fourier basis
            FLOW_CALC=.TRUE.
            NT_COEFFS=IBT(2,NIT(NB_MOTION),NB_MOTION)
C           calculate flow
            FLOW=0.0d0
            DO no_coeffs=1,NT_COEFFS
              FLOW=FLOW+PF1(no_coeffs,1,T+0.5d0*DT)*
     '          FLOW_COEFFS(no_coeffs)
            ENDDO
            WRITE(OP_STRING,'('' Flow at time '',D12.3,'
     '        //''' is '',D12.3)') T+0.5d0*DT,FLOW
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C           calculate volume change from original undeformed mesh
            VOL_CH=0.0d0
            DO no_coeffs=1,NT_COEFFS
              VOL_CH=VOL_CH+PF1(no_coeffs,-1,T+0.5d0*DT)
     '          *FLOW_COEFFS(no_coeffs)
            ENDDO
            WRITE(OP_STRING,'('' Volume change at time '',D12.3,'
     '        //''' is '',D12.3)')T+0.5d0*DT,VOL_CH
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL MESH_FLOW(NEELEM,NPNE,nr,CE,
     '        FLOW,VOL_CH,XP,ERROR,*9999)
            IF(FLOW.GE.0.0d0) THEN !Bdry condition at mouth of model
              FIX(1,1)=.TRUE.
              YP(1,1)=Bdry_Conc
              YP(1,4)=Bdry_Conc
            ELSE !change to zero flux b.c. when flow is out of lung
              FIX(1,1)=.FALSE.
              FIX(1,2)=.TRUE.
              YP(1,2) = 0.0d0
            ENDIF

          ELSE IF(KTYP58(nr).EQ.4) THEN !update flow using Lung 
            FLOW_CALC=.TRUE.        !..gas flow coeffs
C           calculate flow & volume change from original undef. mesh
            T_end_inspiration=FLOW_COEFFS(2)
            T_end_breath_hold=FLOW_COEFFS(2)+FLOW_COEFFS(3)
            T_end_expiration =FLOW_COEFFS(2)+FLOW_COEFFS(3)
     '                                      +FLOW_COEFFS(4)
            Tidal_volume=FLOW_COEFFS(1)/23.0D0 !mm^3/acinus
            Peak_flow=Tidal_volume*PI/(2.0D0*FLOW_COEFFS(2))
c           T_cycle=DMOD(T,FLOW_COEFFS(2)+FLOW_COEFFS(3)+FLOW_COEFFS(4))
            T_cycle=T+0.5d0*DT
            IF(T_cycle.LE.T_end_inspiration) THEN      !inspiration
              FLOW  = Peak_flow*DSIN(T_cycle*PI/FLOW_COEFFS(2))
              VOL_CH= VOL_CH+FLOW*DT
              WRITE(OP_STRING,'('' Time = '',D12.3,'
     '          //''' Inspiration: '','' Flow = '',D12.3,'
     '          //''' Volume change = '',D12.3)')
     '          T+0.5D0*DT,FLOW,VOL_CH
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ELSE IF(T_cycle.LE.T_end_breath_hold) THEN !breath-hold
              FLOW  = 0.0D0
              VOL_CH= Peak_flow*FLOW_COEFFS(2)/PI*2.0D0
              WRITE(OP_STRING,'('' Time = '',D12.3,'
     '          //''' Breath-hold: '','' Flow = '',D12.3,'
     '          //''' Volume change = '',D12.3)')
     '          T+0.5D0*DT,FLOW,VOL_CH
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ELSE IF(T_cycle.LE.T_end_expiration ) THEN !expiration
              T_cycle=T+0.5d0*DT
              T_expire=T_cycle-T_end_breath_hold  
              FLOW  =-Peak_flow*DSIN(T_expire*PI/FLOW_COEFFS(4))
              VOL_CH=VOL_CH+FLOW*DT
              WRITE(IOOP,'('' Time = '',D12.3,'' Expiration : '','
     '          //''' Flow = '',D12.3,'' Volume change = '',D12.3)')
     '          T+0.5D0*DT,FLOW,VOL_CH
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL MESH_FLOW(NEELEM,NPNE,nr,CE,FLOW,VOL_CH,XP,ERROR,*9999)
            IF(FLOW.GE.0.0D0) THEN !Bdry condition at mouth of model
              FIX(1,1)=.TRUE.
              YP(1,1)=Bdry_Conc
            ELSE !change to zero flux b.c. when flow is out of lung
              FIX(1,1)=.FALSE.
              ne=NEELEM(0,nr)+2
              FIX(ne,1)=.TRUE.
              YP(ne,1) = 0.0d0
            ENDIF

            IF(T_cycle.GE.T_end_inspiration.AND.
     '         T_cycle-T_end_inspiration.LE.DT) THEN
              WRITE(IOFILE4,'(/'' End inspiration:''/,1X,15(''=''))')
              LUNG_PRINT=.TRUE.
            ELSE IF(T_cycle.GE.T_end_breath_hold.AND.
     '              T_cycle-T_end_breath_hold.LE.DT) THEN
              WRITE(IOFILE4,'(/'' End breath_hold:''/,1X,14(''=''))')
              LUNG_PRINT=.TRUE.
            ELSE IF(DABS(T_end_expiration-T_cycle).LE.1.5d0*DT) THEN
              WRITE(IOFILE4,'(/'' End expiration:''/,1X,15(''=''))')
              LUNG_PRINT=.TRUE.
            ELSE
              LUNG_PRINT=.FALSE.
            ENDIF
            IF(LUNG_PRINT) THEN
              WRITE(IOFILE4,'(/'' Time = '',D12.4)') T
              WRITE(IOFILE4,'(1X,3A5,4A12)') 'Elem','nd#1','nd#2',
     '		'Gen #','Volume','Conc 1','Conc 2'
              nc=1 !TEMPORARY AJP 17-12-91
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                nb=NBH(NH_LOC(1,nx),nc,ne)
                NP1=NPNE(1,nb,ne)
                NP2=NPNE(2,nb,ne)
                C1=YP(NP1,1)
                C2=YP(NP2,1)
                WRITE(IOFILE4,'(1X,3I5,4D12.4)') 
     '            ne,NP1,NP2,CE(1,ne),CE(4,ne),C1,C2
              ENDDO
            ENDIF !lung_print
          ENDIF !ktyp58(nr)
        ENDIF !ityp2=5 & ityp3=2

C*** Check if the stiffness matrices needs to be (re)calculated and if
C*** so (re)calculate them.

        IF(UPDATE_MATRIX.AND.NSTEP.EQ.1) THEN !First time so assemble
          CALL ASSEMBLE3(IBT,IDO,INP,LGE,NBH,NBJ,NEELEM,NHE,NHP,NJE,
     '      NKE,NKH,NPF,NPNE,NPNODE,NQE,nr,NVHE,
     '      NVHP,NVJE,NVJP,NW,nx,NYNE,NYNP,NYNR(0,0,1,nr),NZNY,
     '      CE,CG,CP,ED,EM,ER,ES,GD,GK,GM,GR,PG,RG,SE,VE,WG,
     '      XA,XE,XG,XP,YG,YP,ZA,ZE,ZG,ZP,DYNAM1,DYNAM2,FIRST,
     '      UPDATE_MATRIX,UPDATE_VECTOR,ERROR,*9999) 
        ELSE
          IF(IP.EQ.1) THEN
            IF(FLOW_CALC) THEN
              IP=1 !recalculate stiffness matrix
              UPDATE_MATRIX=.TRUE.
              CALL ASSEMBLE3(IBT,IDO,INP,LGE,NBH,NBJ,NEELEM,NHE,NHP,NJE,
     '          NKE,NKH,NPF,NPNE,NPNODE,NQE,nr,NVHE,
     '          NVHP,NVJE,NVJP,NW,nx,NYNE,NYNP,NYNR(0,0,1,nr),NZNY,
     '          CE,CG,CP,ED,EM,ER,ES,GD,GK,GM,GR,PG,RG,
     '          SE,VE,WG,XA,XE,XG,XP,YG,YP,ZA,ZE,ZG,ZP,DYNAM1,DYNAM2,
     '          FIRST,UPDATE_MATRIX,UPDATE_VECTOR,ERROR,*9999) 
            ELSE
              IP=3 !use previous stiffness matrix
              UPDATE_MATRIX=.FALSE.
              UPDATE_VECTOR=.FALSE.
            ENDIF !flow_calc
          ENDIF !ip=1
        ENDIF !update_matrix & nstep=1

        IF(NWRIT.EQ.IWRIT1(nr,nx)) THEN
          OUTPUT=.TRUE.
          NWRIT =1
        ELSE
          OUTPUT=.FALSE.
          NWRIT=NWRIT+1
        ENDIF

C*** ?????

        IF(T_cycle.GE.T_end_inspiration.AND.
     '     T_cycle-T_end_inspiration.LE.DT) THEN
C          CALL GLOBALH(NONY,NP_INTERFACE,NPNY,nr,nx,
C     '      NYNE,NYNO,NYNP,NYNR,CONY,CYNO,FIX,ERROR,*9999)
        ENDIF

        no_unit=0

C*** Adjust any incremental boundary conditions

        DO no_nynr=1,NYNR(0,0,1,nr) !loop over global variables of GK 
          ny=NYNR(no_nynr,0,1,nr) !global variable #
          IF(FIX(ny,2)) THEN !incremental bc
            YP(ny,1)=YP(ny,1)+YP(ny,2)
          ENDIF
        ENDDO !no_nynr

C*** Calculate any initial accelerations etc.

        IF(IP.EQ.1) THEN
          INIT=.TRUE.
          IF(KTYP22.EQ.3) THEN !cubic time-stepping algorithm
C!!! See backup version of march1
          ENDIF
        ENDIF

        A1=DT*THETA(1)
        IF(KTYP22.GE.2) A2=DT**2/2.0d0*THETA(2)        
        IF(KTYP22.EQ.3) A3=DT**3/6.0d0*THETA(3)

        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' A1='',D10.3,'' A2='',D10.3,'' A3='','
     '      //'D10.3)') A1,A2,A3
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO i=1,7
            WRITE(OP_STRING,'(/'' YP(ny,'',I2,''):'',10D11.3,'
     '	      //':(/11X,10D11.3))') i,(YP(NYNR(no_nynr2,0,1,nr),i),
     '        no_nynr2=1,NYNR(0,0,1,nr))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)    
          ENDDO
        ENDIF

C*** Solve the problem for the current time step, iterating on the
C*** solution process if necessary.

        CONVERGED=.FALSE.
        KOUNT=0
        DO WHILE(.NOT.CONVERGED)
          KOUNT=KOUNT+1
          CALL SOLVE3(IDISP,IP,ISC_GKK,ISC2_GKK,ISR_GKK,ISR2_GKK,
     '      IWK1,IWK2,NONY(0,1,1,nr),
     '      NPNY,nr,nx,NYNE,NYNO(0,1,1,nr),NYNP,
     '      NYNR(0,0,1,nr),NZNY,
     '      CONY(0,1,1,nr),CYNO(0,1,1,nr),GD,GK,GKK,GM,GR,GRR,WK1,
     '      XO,YP,DYNAM1,DYNAM2,FIX,ERROR,*9999)
          IF(LINEAR) THEN
            CONVERGED=.TRUE.
          ELSE IF(.NOT.LINEAR) THEN
C!!!cpb 9/12/94 This needs to be looked at.
            SUM=0.0D0
            IF(KOUNT.GT.1) THEN !check for convergence
              DO no_nynr=1,NYNR(0,0,1,nr) !Loop over global vars of GK
                ny=NYNR(no_nynr,0,1,nr) !global variable #
                SUM=SUM+(YP(ny,1)-YP(ny,15))**2
              ENDDO
              IF(SUM.LT.1.0D-6) CONVERGED=.TRUE.
            ENDIF
            WRITE(OP_STRING,'('' Nonlinear iteration: Kount='',i3,'
     '	      //''' sum='',e12.3)')K OUNT,SUM
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            !Note: yp(1) on entry to solve3 is an estimate of new time 
            !step solution and on exit is time derivative
            IF(.NOT.CONVERGED) THEN
              DO no_nynr=1,NYNR(0,0,1,nr) !Loop over global vars of GK
                ny=NYNR(no_nynr,0,1,nr) !global variable #
                YP(ny,15)=YP(ny,1) !temporary storage of slope
                YP(ny,1)=YP(ny,4)+DT*YP(ny,1) !is new time sol estimate
              ENDDO
            ENDIF !not.converged
          ENDIF !linear
        ENDDO !not.converged

C*** Adjust any time varying parameters.

        IF(ITYP2(nr,nx).EQ.5.AND.ITYP3(nr,nx).EQ.2) THEN !Lung model
          IF(KTYP58(nr).EQ.4) THEN !update flow using Lung flow coeffs
            IF(T_cycle.LE.T_end_inspiration) THEN      !inspiration
              Conc_N2=0.0D0 !Nitrogen conc during inspiration at node 1
              Flow_N2_integral=0.0D0
              d_Conc_N2_dVol_max=0.0D0
            ELSE IF(T_cycle.LE.T_end_expiration ) THEN !expiration
              Prev_Conc_N2=Conc_N2      !stores previous time step conc
              Conc_N2=Bdry_Conc-YP(1,1) !is latest [N2] at node 1
              d_Conc_N2_dt=(Conc_N2-Prev_Conc_N2)/DT !rate change [N2]
              d_Conc_N2_dVol=-d_Conc_N2_dt/FLOW
              Flow_N2_integral=Flow_N2_integral+FLOW*Conc_N2
              WRITE(OP_STRING,'(''  [N2]      = '',D12.3)') Conc_N2
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' d[N2]/dt   = '',D12.3)')
     '          d_Conc_N2_dt
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' d[N2]/dvol = '',D12.3)')
     '		d_Conc_N2_dVol
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Flow*[N2]  = '',D12.3)')
     '          FLOW*Conc_N2
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Integral   = '',D12.3)')
     '		FLOW_N2_integral
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              IF(d_Conc_N2_dVol.GT.d_Conc_N2_dVol_max) THEN
                d_Conc_N2_dVol_max=d_Conc_N2_dVol
                Time_Max_d_Conc_N2_dVol=T_cycle
              ENDIF
            ENDIF !t_cycle
          ENDIF !ktyp58(nr)=4
        ENDIF !ityp2 & ityp3

C***    Write output

        IF(OUTPUT) THEN
          WRITE(OP_STRING,'(/'' Solution at time T+DT='',D11.4,'
     '      //''' with DT='',D11.4)') T+DT,DT
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(59,'(/'' Solution at time T+DT='',D11.4,'
     '      //''' with DT='',D11.4)')  T+DT,DT
        ENDIF
        IF(OUTPUT) THEN
          CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '      nr,NVHP,nx,NYNE,NYNP,
     '      YP,ZA,ZP,ERROR,*9999)
          CALL ZPOP(4,NBH,1,NEELEM,NHE,NHP,NKH,NPNODE,nr,
     '      NVHP,nx,NYNE,NYNP,YP(1,4),ZA,ZP,FIX(1,4),ERROR,*9999)
        ENDIF

C CPB 26/5/95 NEW HISTORY FILE FORMAT
C        WRITE(IOFILE2,'('' YP(ny,1,nx) at t='',D11.4,'' :'')') T
C        WRITE(IOFILE2,'(:(1X,10(D12.5,X)))') (YP(NYNR(no_nynr2,0,1,nr),
C     '    1),no_nynr2=1,NYNR(0,0,1,nr))

        WRITE(IOFILE2,'(/'' YP(ny,1) at t='',D11.4,'' :'')') T
        WRITE(IOFILE2,'(X,5D13.5,/:(X,5D13.5)))') 
     '    ((YP(NYNR(no_nynr,0,nc,nr),1),no_nynr=1,NYNR(0,0,nc,nr)),
     '      nc=1,NCT(nr))

C       WRITE(IOFILE2,'('' YP(ny,5,nx) at t='',D11.4,'' :'')') T
C       WRITE(IOFILE2,'(:(1X,10(D12.5,X)))') (YP(NYNR(no_nynr2,0,1,nr),
C     '    5),no_nynr2=1,NYNR(0,0,1,nr))

        WRITE(IOFILE3,'(I5,9D12.4)') NSTEP,T+DT,
     '    (YP(NODE_HISTORY(N),1),N=1,NODE_HISTORY(0))

C*** Adjust time stepping parameters if required

        IF(KTYP23.EQ.2) THEN
          KFAC=1
          DO K=1,KTYP22+1
            KFAC=KFAC*K
          ENDDO
          ERR=0.0d0
          DO no_nynr=1,NYNR(0,0,1,nr) !Loop over global variables of GK
            ny=NYNR(no_nynr,0,1,nr) !global variable #
            YP(ny,14)=(YP(ny,13)-YP(ny,14))*DT**KTYP22/KFAC
            ERR=ERR+YP(ny,14)**2
          ENDDO
          ERR=DSQRT(ERR)
          SNORM=0.0D0
          DO no_nynr=1,NYNR(0,0,1,nr) !Loop over global variables of GK
            ny=NYNR(no_nynr,0,1,nr) !global variable #
            SNORM=SNORM+YP(ny,1)**2
          ENDDO
          SNORM=DSQRT(SNORM)
          TOL=SNORM
          IF(ERR.GT.TOL) THEN
            DT=DT/2.0d0
            INCR=0
          ELSE IF(ERR.LT.TOL/2.0d0) THEN
            INCR=INCR+1
            IF(INCR.GT.1) THEN
              DT=DT*1.250d0
              INCR=0
            ENDIF
          ENDIF
        ENDIF

C*** Transfer current solution to previous solution and increment time

        T=T+DT
        DO no_nynr2=1,NYNR(0,0,1,nr) !loop over global variables of GK
          ny=NYNR(no_nynr2,0,1,nr) !global variable #
          YP(ny,5)=YP(ny,1)  !is solution vector at time T
        ENDDO

	IF(DOP) THEN
	  WRITE(OP_STRING,'('' TIMES : CURRENT ='',D10.4,'
     '	    //''' DELTA='',D10.4,'' INIT='',D10.4,'' FINAL='',D10.4)')
     '	    T,DT,TSTART,TFINISH
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)    
	ENDIF

C*** Check for numeric keypad entry interrupt

        CALL GETSTR2(ERROR,*9998)

C*** Check if all the current time is less than the finish time and
C*** repeat solution if so.

      IF(T.LT.TFINISH) GOTO 10
                       
      CALL CLOSEF(IOFILE2,ERROR,*9999)
      CALL CLOSEF(IOFILE3,ERROR,*9999)
      CALL CLOSEF(IOFILE4,ERROR,*9999)

      IF(FILEIP) THEN
        DO no_unit=1,NT_UNIT
          CALL CLOSEF(NY_UNIT(no_unit),ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('MARCH1')
      RETURN

 9998 T_RESTART=T
      CALL EXITS('MARCH1')
      RETURN

 9999 IF(FILEIP) THEN
        DO no_unit=1,NT_UNIT
          CLOSE(UNIT=no_unit)
        ENDDO
      ENDIF
      CALL ERRORS('MARCH1',ERROR)
      CALL EXITS('MARCH1')
      RETURN 1
      END


      SUBROUTINE MARCH4(IBT,IDO,INP,ISCONO,ISCONT,ISEG,ISELNO,ISFIBR,
     '  ISFIEL,ISGAUS,ISHIST,ISINCR,ISISOC,ISLINE,ISLINO,ISNONO,ISREAC,
     '  ISSECT,ISSTRE,ISSTRM,ISSURF,ISVELO,ITHRES,IWK,
     '  MXI,NAN,NBH,NBJ,NCNE,NCNP,NEELEM,NELIST,NGAP,NHE,NHP,NJE,
     '  NJP,NKE,NKH,NKJ,NLL,NLLIST,NPE,NPF,NPL,NPLIST,NPNODE,NQE,
     '  nr,NRE,NTCOVA,NTIW,NW,nx,NXI,NYNE,NYNP,A,CE,CG,COVA,
     '  CP,CSEG,DL,PG,SE,THRES,VE,XA,XE,XF,XG,XGRC,XIG,
     '  XP,YG,YP,ZA,ZE,ZF,ZG,ZP,FIX,ERROR,*)

C**** Performs time integration of cardiac activation equations
C**** with variable time step algorithm.
C**** YP(ny,1) is solution vector at time T+DT
C**** YP(ny,3) is incremental boundary conditions
C****       4  "  solution vector at time T
C****       5  "  reaction vector at time T+DT
C**** FIX(ny,5) is .true. for input from a file (FILE07)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b08.cmn'
      INCLUDE 'cmiss$reference:b10.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:gks000.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp30.cmn'
      INCLUDE 'cmiss$reference:suben00.cmn'
      INCLUDE 'cmiss$reference:time02.cmn'      
!     Parameter List
      INTEGER IBT(2,NIM,*),IDO(NKM,0:NIM,*),INP(NNM,NIM,*),
     '  ISCONO(NHM,*),ISCONT(NHM,NEM,*),ISEG(*),ISELNO(NWM,*),
     '  ISFIBR(NWM,NEM,*),ISFIEL(NWM,*),ISGAUS(NWM,NGM,*),ISHIST(0:*),
     '  ISINCR(*),ISISOC(*),ISLINE(NWM,*),ISLINO(*),ISNONO(NWM,*),
     '  ISREAC(*),ISSECT(*),ISSTRE(NEM,*),ISSTRM(NEM,*),ISSURF(NWM,*),
     '  ISVELO(NEM,*),ITHRES(3,NGM,*),IWK(*),MXI(2,*),
     '  NAN(NIM,NAM,*),NBH(NHM,NCM,*),NBJ(NJM,*),
     '  NCNE(NNM,NBM,NEM,*),NCNP(NHM,0:NCM,NPM,*),
     '  NEELEM(0:NEM,0:*),NELIST(0:*),
     '  NGAP(NIM,*),NHE(*),NHP(*),NJE(*),NJP(*),
     '  NKE(NKM,NNM,NBM,*),NKH(NHM,NCM,*),NKJ(NJM,*),
     '  NLL(12,*),NLLIST(0:*),
     '  NPE(NNM,NBM,*),NPF(12,*),NPL(20,*),NPLIST(0:*),
     '  NPNODE(0:NPM,0:*),NQE(NSM,NBM,*),nr,NRE(*),
     '  NTCOVA(*),NTIW,NW(NEM,*),nx,NXI(-NIM:NIM,0:*),
     '  NYNE(NAM,NHM,NCM,*),NYNP(NKM,NHM,0:NCM,NPM,0:*)
      REAL*8 A(NSM,*),CE(NMM,*),CG(NMM,*),
     '  COVA(NEM,*),CP(NMM,*),
     '  DL(3,*),PG(NSM,NUM,NGM,*),SE(NSM,NBM,*),THRES(3,NGM,*),
     '  VE(NSM,NKM,*),
     '  XA(NAM,NJM,*),XE(NSM,*),XF(NSM,*),XG(NJM,*),XGRC(NJM,*),
     '  XIG(NIM,NGM,*),
     '  XP(NKM,NJM,*),YG(NGM,NJM,*),YP(NYM,*),
     '  ZA(NAM,NHM,NCM,*),ZE(NSM,*),ZF(NSM,*),ZG(NHM,*),
     '  ZP(NKM,NHM,NCM,*)
      CHARACTER CSEG(*)*(*),ERROR*(*)
      LOGICAL FIX(NYM,*)
!     Local Variables
      INTEGER NEIGM,NEMAX,NGMAX,NHISMX,NSAMX
      PARAMETER (NEIGM=27,NEMAX=24,NGMAX=125,NHISMX=300,NSAMX=50)
      INTEGER I,IBEG,IBEG1,IBEG5,IEND,IEND1,IEND5,IHANDLE,
     '  IME(2000),IMG(2000),INCR,INE(2000),ING(2000),IP,ISTATUS,
     '  LE(NEMAX),LG(NEMAX,NGMAX),MAXIT,MSA,MSA1,MSA2,
     '  N,NA,NB,NC,ND,NE,NEDUM,NEHIS(NHISMX),NEND,NG,NG1,NGDUM,
     '  NGHIS(NHISMX),
     '  NH,NHIS,NITB,NJ,NK,NLE1,NLET,NLGT(NEMAX),NMAX,NNEIG,
     '  NOELEM,NONODE,NO_UNIT,NP,NPHIS,NPTTT,
     '  NSE1,NSE2,NSE3,NSTEP,NT_UNIT,NUAT,NUATT,
     '  NWRIT,NY,NY_UNIT(20)
      REAL*8 ATIME(600),CPULAS,CPUTOT,DUM1,ELAPSED,SSTA(NSAMX),
     '  STA(NGMAX,NEMAX),
     '  T,TA(NGMAX,NEMAX,NEIGM),THIS(NHISMX),THRES3(125,24),TIMER
      CHARACTER CFROMI*5,CHAR*5,FILE*100
      LOGICAL EFLAG,FILEIP,LINEAR,OUTPUT

      CALL ENTERS('MARCH4',*9999)
      NC=1 !Temporary AJP 17-12-91
      IP=1
      NWRIT=1
      NSTEP=0
      INCR=0
      T=TSTART
      DT=TINCR
      IF(ITYP6(nr).EQ.1) THEN !linear equations
        LINEAR=.TRUE.
      ELSE                    !nonlinear equations
        LINEAR=.FALSE.
      ENDIF
      FILEIP=.FALSE.
      NO_UNIT=0
      DO NONODE=1,NPNODE(0,NR)
        NP=NPNODE(NONODE,NR)
        CHAR=CFROMI(NP,'(I5)')
        CALL TRIM(CHAR,IBEG1,IEND1)
        DO NH=1,NHP(NP)
          DO NK=1,NKH(NH,NC,NP)
            ny=NYNP(nk,nh,1,np,nr)
            IF(NK.EQ.1.AND.FIX(NY,5)) THEN
              NO_UNIT=NO_UNIT+1
              NY_UNIT(NO_UNIT)=20+NO_UNIT
              CALL TRIM(FILE07,IBEG5,IEND5)
              FILE=FILE07(IBEG5:IEND5)//'.iphist_'
     '          //CHAR(IBEG1:IEND1)
              FILEIP=.TRUE.
              CALL TRIM(FILE,IBEG,IEND)
              CALL OPENF(NY_UNIT(NO_UNIT),'DISK',FILE(IBEG:IEND),'OLD',
     '          'SEQUEN','FORMATTED',132,ERROR,*9999)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' FILE='',A)') FILE(IBEG:IEND)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NT_UNIT=NO_UNIT
      CALL TRIM(FILE02,IBEG,IEND)
      CALL OPENF(IOFILE2,'DISK',FILE02(IBEG:IEND)//'.history','NEW',
     '  'SEQUEN','FORMATTED',132,ERROR,*9999)

      IHANDLE=0
      ISTATUS=TIMER(IHANDLE,T_CPU)
      nuat=0
      nuatt=0
      NB=NBJ(1,1)
      NITB=NIT(NB)
      IF(NGT(NB).EQ.3**NITB) THEN
        NMAX=3
      ELSE IF(NGT(NB).EQ.5**NITB) THEN
        NMAX=5
      ELSE
        ERROR='>>>Incorrect # of Gauss points'
        GO TO 9999
      ENDIF
      CALL CPCG(1,NBJ(1,1),NPE(1,1,1),CE(1,1),CG,CP,PG,ERROR,*9999)
      IF(KTYP23.EQ.2) THEN !variable time step model
        MSA1=NINT(CG(3,1))  !starting value of # of points to activate
        MSA2=NINT(CG(4,1))  !upper limit of # of points activated at once
      ENDIF
      NSE1=NINT(CG(5,1))  !# of Gauss points along a search in Xi1,2 dir.
      NSE3=NINT(CG(6,1))  !# of Gauss points along a search in Xi3 dir.
!      namx=1 !?????
!      NSUBLV=NINT(CG(7,1))*NAMX**2 !# of faster LV subendo Gauss points
!      NSUBRV=NINT(CG(8,1))*NAMX**2 !# of faster RV subendo Gauss points
      VSUBEN=CG(9,1)      !speedup factor in subendocardial layers
!     NSUBEN=NSUBLV !!!!temporary
c     read(57,*)irvsuben     !flag for including single layer rv Purkinje
      NGTOP_LAYER=NMAX*NMAX*(NMAX-1)+1  !Gauss point no of top layer
      NGBOT_LAYER=NMAX*NMAX             !Gauss point no of bottom layer
      NSE2=(2*NSE1+1)*(2*NSE1+1)
      IF(NITB.GT.2) THEN
        IF(NSE3.EQ.1) THEN
          NSE2=NSE2+18
        ELSE
          NSE2=NSE2*3
        ENDIF
      ENDIF
      NSE2=NSE2-1
      IF(NSE2.GT.NEIGM) THEN
      	WRITE(OP_STRING,'('' NSE2='',I4,'' > MAX='',I4)') 
     '    NSE2,NEIGM
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	CALL EXIT
      ENDIF
      WRITE(OP_STRING,'('' NSE1='',I3,'' NSE2='',I3,'' NSE3='',I3)')
     '  NSE1,NSE2,NSE3
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      NHIS=0
      MAXIT=0
      WRITE(OP_STRING,'(//,
     '  '' ***INITIAL CONDITIONS*** - ISOCHRONE NUMBER 1'',/)')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(IOFILE2,'('' ***INITIAL CONDITIONS*** - '','
     '	//'''ISOCHRONE NUMBER 1'')')
      DO NOELEM=1,NEELEM(0,NR)
        NE=NEELEM(NOELEM,NR)
      	NB=NBJ(1,NE)
      	MAXIT=MAXIT+NGT(NB)
        IF(NITB.EQ.2) THEN
          WRITE(OP_STRING,'('' ne='',I4,'' ITHRES: '',25I1)')
     '      NE,(ITHRES(1,NG,NE),NG=1,NGT(NB))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(IOFILE2,'('' ne='',I4,'' ITHRES: '',25I1)')
     '      NE,(ITHRES(1,NG,NE),NG=1,NGT(NB))
          WRITE(IOFILE2,'((9E12.3))') (THRES(1,NG,NE),NG=1,NGT(NB))
        ELSE IF(NITB.EQ.3.AND.NET(nr).EQ.8) THEN
          WRITE(OP_STRING,'('' ne='',I4)')NE
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' ITHRES: '',25I1)')
     '      (ITHRES(1,NG,NE),NG=1,NGT(NB))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(IOFILE2,'('' ne='',I4)')NE
          WRITE(IOFILE2,'('' ITHRES: '',25I1)')(ITHRES(1,NG,NE),NG=1,
     '	    NGT(NB))
          WRITE(IOFILE2,'((9E12.3))') (THRES(1,NG,NE),NG=1,NGT(NB))
        ENDIF
        DO NG=1,NGT(NB)
          IF(THRES(1,NG,NE).GT.1.0D-5) THEN
            IF(KTYP23.EQ.2) THEN
              NHIS=NHIS+1
              IF(KTYP31.EQ.1) THEN
                THIS(NHIS)=THRES(1,NG,NE)
              ELSE
                THIS(NHIS)=-THRES(1,NG,NE)
              ENDIF
              NGHIS(NHIS)=NG
      	      NEHIS(NHIS)=NE
	            DO NH=1,NHIS-1
                IF(THIS(NH).GT.THIS(NHIS)) THEN
                  DUM1=THIS(NHIS)
                  NGDUM=NGHIS(NHIS)
                  NEDUM=NEHIS(NHIS)
                  THIS(NHIS)=THIS(NH)
                  NGHIS(NHIS)=NGHIS(NH)
                  NEHIS(NHIS)=NEHIS(NH)
                  THIS(NH)=DUM1
                  NGHIS(NH)=NGDUM
                  NEHIS(NH)=NEDUM
                ENDIF
              ENDDO
            ELSE
              THRES(2,NG,NE)=THRES(1,NG,NE)
            ENDIF
            THRES(1,NG,NE) = 0.0D0
          ELSE
C           IF(ITHRES(1,NG,NE).EQ.1.AND.KTYP31.EQ.2) THEN
            IF(ITHRES(1,NG,NE).EQ.1.AND.KTYP23.EQ.2.AND.KTYP31.EQ.2) 
     '        THEN
              WRITE(OP_STRING,'('' ITHRES=1,THRES=0,NE='',I3,'
     '          //''' NG='',I4)') NE,NG
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              NHIS=NHIS+1
              IF(KTYP31.EQ.1) THEN
                THIS(NHIS)=THRES(1,NG,NE)
              ELSE
                THIS(NHIS)=-THRES(1,NG,NE)
              ENDIF
              NGHIS(NHIS)=NG
              NEHIS(NHIS)=NE
              DO NH=1,NHIS-1
                IF(THIS(NH).GT.THIS(NHIS)) THEN
                  DUM1=THIS(NHIS)
                  NGDUM=NGHIS(NHIS)
                  NEDUM=NEHIS(NHIS)
                  THIS(NHIS)=THIS(NH)
                  NGHIS(NHIS)=NGHIS(NH)
                  NEHIS(NHIS)=NEHIS(NH)
                  THIS(NH)=DUM1
                  NGHIS(NH)=NGDUM
                  NEHIS(NH)=NEDUM
                ENDIF
              ENDDO
              ITHRES(1,NG,NE)=0
            ENDIF
          ENDIF
          IF(KTYP23.EQ.2) THEN
            STA(NG,NE)=10000.0D0
            IF(ITHRES(1,NG,NE).EQ.1) THEN
              NUAT=NUAT+1
              NUATT=NUATT+1
              ING(NUAT)=NG
              INE(NUAT)=NE
              IMG(NUAT)=0
              IME(NUAT)=0
              ATIME(NUAT)=THRES(2,NG,NE)
            ENDIF
            DO NNEIG=1,NEIGM
              TA(NG,NE,NNEIG)=10000.0D0
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      NPTTT=MAXIT
      NPHIS=0.8*NPTTT
C     WRITE(OP_STRING,'('' NPHIS='',I6)')NPHIS
      MAXIT=MAXIT+50
      NEND=0
      DO I=5,MSA2
        NEND=NEND+I
      ENDDO
      NPTTT=NPTTT-NEND 	
      IF(KTYP23.EQ.2) THEN
        NLET=0
        DO NA=1,NUAT
          EFLAG=.FALSE.
          DO NLE1=1,NLET
            IF(INE(NA).EQ.LE(NLE1)) THEN
              EFLAG=.TRUE.
              LG(NLE1,NLGT(NLE1)+1)=ING(NA)
              NLGT(NLE1)=NLGT(NLE1)+1
            ENDIF
          ENDDO
          IF(.NOT.EFLAG) THEN
            NLET=NLET+1
            LE(NLET)=INE(NA)   !List of initially active elements
            NLGT(NLET)=1
            LG(NLET,1)=ING(NA) !List of initially active gauss points
          ENDIF
        ENDDO
      ELSE
        NUAT=1
        NUATT=1
      ENDIF

      IF(KTYP31.EQ.2) THEN
      	MSA=MSA1-1
      ELSE
      	MSA=MSA1
      ENDIF
 10   CONTINUE                   !loop start
      IF(NUATT.GT.NPTTT.AND.MSA.GT.5) THEN
        MSA=MSA-1              !decrements # Gauss pts which can be activated
      ELSE IF(MSA.LT.MSA2.AND.DT.NE.TINCR) THEN
        MSA=MSA+1              !increments # Gauss pts ...
      ENDIF
      NSTEP=NSTEP+1
      WRITE(OP_STRING,'(//,'' **ISOCHRONE NUMBER '',I4,''**''/)')
     '	NSTEP+1
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
 11   CONTINUE
      IF(NWRIT.EQ.IWRIT1) THEN
        OUTPUT=.TRUE.
        NWRIT =1
      ELSE
        OUTPUT=.FALSE.
        NWRIT=NWRIT+1
      ENDIF
      
      IF(KTYP23.EQ.1) THEN      !constant time step model
        CALL TFRONT(IBT,INP,ITHRES,NBJ,NCNE,NCNP,NEELEM,NJE,NKE,NPE,NPF,
     '    NQE,NRE,NSE1,NSE2,NSE3,NXI,
     '    CE,CG,CP,PG,SE,T,THRES,VE,XA,XE,XG,XIG,
     '    XP,YG,ERROR,*9999)
      ELSE IF(KTYP23.EQ.2) THEN !variable time step model
      	  CALL FFRONT(IBT,IME,IMG,INE,ING,INP,ITHRES,LE,LG,MSA,
     '    NBJ,NEELEM,NEHIS,NEIGM,NEMAX,NGHIS,NGMAX,NHIS,NJE,NKE,NLET,
     '    NLGT,NMAX,NPE,NPF,NPNODE,NQE,NRE,NSE1,NSE2,NSE3,NUAT,
     '    NUATT,NW,nxI,CE,CG,CP,PG,SE,SSTA,STA,T,THIS,THRES,
     '    TA,VE,XA,XE,XG,XIG,XP,ZE,ZG,ERROR,*9999)
      ENDIF
      
      IF(OUTPUT) THEN
        WRITE(OP_STRING,'(/'' Solution at time T+DT='',E11.4,'
     '	  //''' with DT='',E11.4,/)') T+DT,DT
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(IOFILE2,'('' ISOCHRONE NUMBER '',I4)') NSTEP+1
        IF(NITB.EQ.3.AND.NET(1).EQ.24) THEN
          IF(NMAX.EQ.3) THEN     !3*3*3 Gauss point elements
           DO N=1,3
	            DO NG1= 3,1,-1
     	          WRITE(IOFILE2,'(6(1X,12(I1)))')
     '            ((ITHRES(1,NG,NE),NG=1+(NG1-1)*3,3+(NG1-1)*3),
     '            NE=12+4*N,(N-1)*4+13,-1)
     '            ,((ITHRES(1,NG,NE),NG=10+(NG1-1)*3,12+(NG1-1)*3),
     '            NE=12+4*N,(N-1)*4+13,-1)
     '            ,((ITHRES(1,NG,NE),NG=19+(NG1-1)*3,21+(NG1-1)*3),
     '            NE=12+4*N,(N-1)*4+13,-1)
     '            ,((ITHRES(1,NG,NE),NG=1+(NG1-1)*3,3+(NG1-1)*3),
     '             NE=4*N,(N-1)*4+1,-1)
     '            ,((ITHRES(1,NG,NE),NG=10+(NG1-1)*3,12+(NG1-1)*3),
     '            NE=4*N,(N-1)*4+1,-1)
     '            ,((ITHRES(1,NG,NE),NG=19+(NG1-1)*3,21+(NG1-1)*3),
     '            NE=4*N,(N-1)*4+1,-1)
	              WRITE(OP_STRING,'(6(1X,12(I1)))')
     '            ((ITHRES(1,NG,NE),NG=1+(NG1-1)*3,3+(NG1-1)*3),
     '            NE=12+4*N,(N-1)*4+13,-1)
     '            ,((ITHRES(1,NG,NE),NG=10+(NG1-1)*3,12+(NG1-1)*3),
     '            NE=12+4*N,(N-1)*4+13,-1)
     '            ,((ITHRES(1,NG,NE),NG=19+(NG1-1)*3,21+(NG1-1)*3),
     '            NE=12+4*N,(N-1)*4+13,-1)
     '           ,((ITHRES(1,NG,NE),NG=1+(NG1-1)*3,3+(NG1-1)*3),
     '            NE=4*N,(N-1)*4+1,-1)
     '           ,((ITHRES(1,NG,NE),NG=10+(NG1-1)*3,12+(NG1-1)*3),
     '            NE=4*N,(N-1)*4+1,-1)
     '            ,((ITHRES(1,NG,NE),NG=19+(NG1-1)*3,21+(NG1-1)*3),
     '            NE=4*N,(N-1)*4+1,-1)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDDO
            ENDDO
          ELSE IF(NMAX.EQ.5) THEN !5*5*5 Gauss point elements
       	    DO N=1,3
	            DO NG1=NMAX,1,-1
	              WRITE(IOFILE2,'(5(1X,20(I1)))')
     ' 	          ((ITHRES(1,NG,NE),NG=1+(NG1-1)*5,5+(NG1-1)*5),
     '            NE=12+4*N,(N-1)*4+13,-1)
     '            ,((ITHRES(1,NG,NE),NG=26+(NG1-1)*5,30+(NG1-1)*5),
     '            NE=12+4*N,(N-1)*4+13,-1)
     ' 	          ,((ITHRES(1,NG,NE),NG=51+(NG1-1)*5,55+(NG1-1)*5),
     '            NE=12+4*N,(N-1)*4+13,-1)
     ' 	          ,((ITHRES(1,NG,NE),NG=76+(NG1-1)*5,80+(NG1-1)*5),
     '            NE=12+4*N,(N-1)*4+13,-1)
     '	          ,((ITHRES(1,NG,NE),NG=101+(NG1-1)*5,105+(NG1-1)*5),
     '            NE=12+4*N,(N-1)*4+13,-1)
	              IF(DOP) THEN
                  WRITE(OP_STRING,'(5(1X,20(I1)))')
     '	            ((ITHRES(1,NG,NE),NG=1+(NG1-1)*5,5+(NG1-1)*5),
     '              NE=12+4*N,(N-1)*4+13,-1)
     '	            ,((ITHRES(1,NG,NE),NG=26+(NG1-1)*5,30+(NG1-1)*5),
     '              NE=12+4*N,(N-1)*4+13,-1)
     '	            ,((ITHRES(1,NG,NE),NG=51+(NG1-1)*5,55+(NG1-1)*5),
     '              NE=12+4*N,(N-1)*4+13,-1)
     '	            ,((ITHRES(1,NG,NE),NG=76+(NG1-1)*5,80+(NG1-1)*5),
     '              NE=12+4*N,(N-1)*4+13,-1)
     '	            ,((ITHRES(1,NG,NE),NG=101+(NG1-1)*5,105+(NG1-1)*5),
     '              NE=12+4*N,(N-1)*4+13,-1)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO
c           WRITE(IOFILE2,'(/)')
            IF(DOP) THEN
              WRITE(OP_STRING,'(/)')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
      	    DO N=1,3
	            DO NG1= NMAX,1,-1
	              WRITE(IOFILE2,'(5(1X,20(I1)))')
     '	          ((ITHRES(1,NG,NE),NG=1+(NG1-1)*5,5+(NG1-1)*5),
     '            NE=4*N,(N-1)*4+1,-1)
     '	          ,((ITHRES(1,NG,NE),NG=26+(NG1-1)*5,30+(NG1-1)*5),
     '            NE=4*N,(N-1)*4+1,-1)
     '	          ,((ITHRES(1,NG,NE),NG=51+(NG1-1)*5,55+(NG1-1)*5),
     '            NE=4*N,(N-1)*4+1,-1)
     '	          ,((ITHRES(1,NG,NE),NG=76+(NG1-1)*5,80+(NG1-1)*5),
     '            NE=4*N,(N-1)*4+1,-1)
     '	          ,((ITHRES(1,NG,NE),NG=101+(NG1-1)*5,105+(NG1-1)*5),
     '            NE=4*N,(N-1)*4+1,-1)
	              IF(DOP) THEN
                  WRITE(OP_STRING,'(5(1X,20(I1)))')
     '	            ((ITHRES(1,NG,NE),NG=1+(NG1-1)*5,5+(NG1-1)*5),
     '              NE=4*N,(N-1)*4+1,-1)
     '	            ,((ITHRES(1,NG,NE),NG=26+(NG1-1)*5,30+(NG1-1)*5),
     '              NE=4*N,(N-1)*4+1,-1)
     '	            ,((ITHRES(1,NG,NE),NG=51+(NG1-1)*5,55+(NG1-1)*5),
     '              NE=4*N,(N-1)*4+1,-1)
     '	            ,((ITHRES(1,NG,NE),NG=76+(NG1-1)*5,80+(NG1-1)*5),
     '              NE=4*N,(N-1)*4+1,-1)
     '	            ,((ITHRES(1,NG,NE),NG=101+(NG1-1)*5,105+(NG1-1)*5),
     '              NE=4*N,(N-1)*4+1,-1)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO
	        ENDIF
        ENDIF

        DO NOELEM=1,NEELEM(0,NR)
          NE=NEELEM(NOELEM,NR)
          NB=NBJ(1,NE)
          IF(NITB.EQ.2) THEN
            IF(KTYP23.EQ.2) THEN      !variable time step model
              WRITE(OP_STRING,'('' ne='',I4,'' ITHRES: '',27I1)')
     '          NE,(ITHRES(1,NG,NE),NG=1,NGT(NB))
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              WRITE(IOFILE2,'('' ne='',I4,'' ITHRES: '',27I1)')
     '          NE,(ITHRES(1,NG,NE),NG=1,NGT(NB))
              WRITE(IOFILE2,'((9E12.3))') (THRES(1,NG,NE),NG=1,NGT(NB))
            ELSE                      !fixed timestep
              WRITE(OP_STRING,'('' ne='',I4,'' ITHRES: '',27I1)')
     '          NE,(ITHRES(1,NG,NE),NG=1,NGT(NB))
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              WRITE(IOFILE2,'('' ne='',I4,'' ITHRES: '',27I1)')
     '          NE,(ITHRES(1,NG,NE),NG=1,NGT(NB))
              DO NG=1,NGT(NB)
                IF(ITHRES(1,NG,NE).EQ.0) THEN
                  THRES3(NG,NE)=0.0D0
                ELSE
                  THRES3(NG,NE)=T+DT-(THRES(1,NG,NE))
                ENDIF
              ENDDO
              WRITE(IOFILE2,'((9E12.3))') (THRES3(NG,NE),NG=1,NGT(NB))
            ENDIF
          ELSE IF(NITB.EQ.3.AND.NET(nr).EQ.8) THEN
            WRITE(OP_STRING,'('' ne='',I4)')NE
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' ITHRES: '',25I1)')(ITHRES(1,NG,NE),
     '	      NG=1,NGT(NB))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(IOFILE2,'('' ne='',I4)')NE
            WRITE(IOFILE2,'('' ITHRES: '',25I1)')(ITHRES(1,NG,NE),NG=1,
     '        NGT(NB))
            WRITE(IOFILE2,'((9E12.3))') (THRES(1,NG,NE),NG=1,NGT(NB))
          ENDIF
c           IF(NUATT.GT.NPHIS) THEN
c             WRITE(IOFILE2,'('' ne='',I4,'' ITHRES: '',27I1)')
c    '          NE,(ITHRES(1,NG,NE),NG=1,NGT(NB))
c             WRITE(IOFILE2,'((9E12.4))') (THRES(1,NG,NE),NG=1,NGT(NB))
c           ENDIF
        ENDDO

      ENDIF

      IF(UPVU) THEN
        CALL UPVIEW(IBT,IDO,INP,ISCONO,ISCONT,ISEG,ISELNO,ISFIBR,
     '    ISFIEL,ISGAUS,ISHIST,ISINCR,ISISOC,ISLINE,ISLINO,ISNONO,
     '    ISREAC,ISSECT,ISSTRE,ISSTRM,ISSURF,ISVELO,
     '    ITHRES,IWK,MXI,NAN,NBH,NBJ,NCNE,NCNP,
     '    NEELEM,NELIST,NGAP,NHE,NHP,NJE,NJP,
     '    NKE,NKH,NKJ,NLL,NLLIST,NPE,NPF,NPL,NPLIST,
     '    NPNODE,NQE,NRE,NTCOVA,
     '    NTIW,NW,NYNE,NYNP,
     '    CE,CG,COVA,CP,CSEG,DL,PG,SE,VE,XA,XE,XF,XG,XIG,XP,YG,YP,
     '    ZA,ZE,ZF,ZG,ZP,FIX,ERROR,*9999)
      ENDIF

      T=T+DT
      IF(DOP) THEN
        WRITE(OP_STRING,'('' TIMES : CURRENT ='',E10.4,'' DELTA='',
     '    E10.4,'' INIT='',E10.4,'' FINAL='',E10.4)')
     '	  T,DT,TSTART,TFINISH
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

C ***   Check whether all Gauss pts are fully active
      IF(KTYP31.EQ.1) THEN
        NUAT=0
        DO NOELEM=1,NEELEM(0,NR)
          NE=NEELEM(NOELEM,NR)
          NB=NBJ(1,NE)
          DO NG=1,NGT(NB)
            IF(ITHRES(1,NG,NE).EQ.1) NUAT=NUAT+1
          ENDDO
        ENDDO
      ENDIF
     
      IF(NSTEP.LT.MAXIT) THEN
	IF(KTYP31.EQ.2) THEN
	  IF(NSTEP.GT.20) THEN
            IF(NUAT.GT.0.AND.T.GT.TFINISH) GO TO 10
	  ELSE
	    IF(T.GT.TFINISH) GO TO 10
	  ENDIF
	ELSE
	  IF(NUAT.GT.0.AND.T.LT.TFINISH) GO TO 10
	ENDIF
      ELSE
	WRITE(OP_STRING,'('' Exceeded max iterations ='',I8)') MAXIT
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      ELAPSED=TIMER(IHANDLE,T_CPU)
      CPULAS=ELAPSED-CPUTOT
      CPUTOT=ELAPSED
      WRITE(OP_STRING,
     '  '(//,'' Solution took '',I8,'' iterations and '','
     '	//'E11.4,'' s of cpu time'')') NSTEP,CPUTOT
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(59,'(//,'' Solution took '',I8,'' iterations and '','
     '	//'E11.4,'' s of cpu time'')') NSTEP,CPUTOT
      CALL CLOSEF(IOFILE2,ERROR,*9999)

      IF(NITB.EQ.3.AND.NET(nr).EQ.24) THEN
        IF(KTYP31.EQ.1) THEN
          CALL TRIM(FILE02,IBEG,IEND)
          CALL OPENF(IOFILE2,'DISK',FILE02(IBEG:IEND)//'.ipdat','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
          WRITE(IOFILE2,
     '      '('' Model epicardial Gauss point activation times'')')
          ND=0
          DO NE=13,24
            CALL XPXE(NBJ(1,NE),NCNE(1,1,ne,NRE(ne)),NCNP,
     '        NJE(NE),NKE(1,1,1,NE),NPE(1,1,NE),
     '        NPF(1,1),NQE(1,1,NE),NRE(ne),SE(1,1,NE),XA,XE,XP,
     '        ERROR,*9999)
            DO NG=19,27
              CALL XEXG(NBJ(1,NE),NG,NJE(NE),PG,VE(1,1,NE),XE,XG,ERROR,
     '          *9999)
              XG(4,1)=THRES(1,NG,NE)
              ND=ND+1
              WRITE(IOFILE2,'(I3,4E12.4,'' 1 1 1 1'')')
     '		ND,(XG(NJ,1),NJ=1,4)
            ENDDO
          ENDDO
          CALL CLOSEF(IOFILE2,ERROR,*9999)
        ENDIF
      ENDIF
      IF(FILEIP) THEN
        DO NO_UNIT=1,NT_UNIT
          CALL CLOSEF(NY_UNIT(NO_UNIT),ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('MARCH4')
      RETURN
 9999 IF(FILEIP) THEN
        DO NO_UNIT=1,NT_UNIT
          CLOSE(UNIT=NO_UNIT)
        ENDDO
      ENDIF
      CALL ERRORS('MARCH4',ERROR)
      CALL EXITS('MARCH4')
      RETURN 1
      END


C***  The following archived version of MARCH4 is complete and 
C***  identical to the pre-stripped version on 9-8-01.  Exporting and 
C***  output options have been moved to other routines (by PM).  The 
C***  non-archived version will now have large 'commented out' regions
C***  removed, particularly those for outputing data.  MHT 9-8-01.

      SUBROUTINE MARCH4(BC_POINTS,BRANCH,CALCULATED,CONECT,CQ,ER,ES,
     '  GKK,GRR,ISC_GKK,ISR_GKK,LGE,NBJ,NEELEM,NENQ,NHQ,NPNE,NQNE,NQS,
     '  nr,nx,NXQ,NYNQ,NYQNR,XP,XIP,XQ,YQ,TIME_VALUES,
     '  TIME_VARAIBLE_NAMES,NTIME_POINTS,NTIME_INTERP,ERROR,*)

C#### Subroutine: MARCH4
C###  Description:
C###    MARCH4 performs time integration of linear or nonlinear 
C###    equations using explicit or implicit finite differences.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbfe01.cmn'
      INCLUDE 'cmiss$reference:coro00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grid00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'                             
      INCLUDE 'cmiss$reference:ktyp00.cmn'                             
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
      INCLUDE 'cmiss$reference:mach00.inc'
      INCLUDE 'cmiss$reference:time02.cmn'
      INCLUDE 'cmiss$reference:tol00.cmn'

!     Parameter List
      INTEGER BC_POINTS(3,3,0:NQM),ISC_GKK(NISC_GKKM),
     '  ISR_GKK(NISR_GKKM),LGE(NHM*NSM,NRCM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),NHQ(NRM),
     '  NPNE(NNM,NBFM,NEM),NQNE(NEM,NQEM),NQS(NEM),nr,nx,
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM),NYNQ(NHM,NQM,0:NRCM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM),
     '  NTIME_POINTS(NTIMEVARSM),NTIME_INTERP(NTIMEVARSM)              
      REAL*8 CQ(NMM,NQM),ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),GKK(NZ_GKK_M),
     '  GRR(NOM),TIME,XQ(NJM,NQM),XP(NKM,NVM,NJM,NPM),XIP(NIM,NPM),
     '  YQ(NYQM,NIQM),TIME_VALUES(2,0:NTIMEPOINTSM+1,NTIMEVARSM)
      CHARACTER ERROR*(*),TIME_VARAIBLE_NAMES(NTIMEVARSM)*(*)         
      LOGICAL CALCULATED(NQM)
!     Local Variables
      INTEGER ADJACENT,BIFUR_COUNT,CONECT(-1:1,0:2,NQM),COUNT,CURRENT,
     '  I,ii,j,jj,LABLED,DATA_COUNT,FILE_COUNT,OLD_CURENT,nb1,nb2,ne1,
     '  ne2,ne3,ne_current,nj,nj1,nj2,njj,np1,np2,nq,nqq,nq1,nq1_pre,
     '  nq2,nq2_pre,nq3,nq_old,no_bc_points,no_ne,no_nq,
     '  no_nynr,NUM_NODES,ny,ny_p,ny_r,ny_v,POINTS(3),TIME_STEPS,
     '  PATH_ARRAY(60),IBEG,IEND 
      REAL*8 FA,FB,GRAD_DIF,LAMBDA,LAMBDA1,LEN_COUNT,STEP,VEL,TRACE1,
     '  TRACE2,TRACE3,XI_DIST1,XI_DIST2, XI_SUM1,XI_SUM2
      REAL TIME_START2(1),TIME_STOP(1),TOT_BITIME,TOT_BOTIME,TOT_GRTIME
      CHARACTER FMT*100
      INTEGER*4 WORK_PTR
      LOGICAL BRANCH(NQM),DISCONT,DONE,FOUND,UPDATE_MATRIX,
     '  UPDATE_VECTOR,HALF_TIME_STEP,PRINT_FILE

      CALL ENTERS('MARCH4',*9999)

      FMT='(4E18.10)'
      IF(ITYP16(nr,nx).EQ.4) THEN !Lax-Wendroff  
        IF(ITYP3(nr,nx).EQ.1) THEN !flow in elastic tube  
          DO nq=1,NQT
            CALCULATED(nq)=.FALSE.
          ENDDO
          CALCULATED(NQ_START(nr))=.TRUE.!ipinit NPS 4/2/97
          CURRENT=NXQ(1,1,NQ_START(nr),1)
          IF(CURRENT.EQ.0) THEN
            CURRENT=NXQ(-1,1,NQ_START(nr),1)
          ENDIF
          CONECT(1,0,NQ_START(nr))=1
          CONECT(-1,1,CURRENT)=NQ_START(nr)
          CONECT(-1,0,CURRENT)=1
          CONECT(1,1,NQ_START(nr))=CURRENT
          OLD_CURENT=NQ_START(nr)
          NUM_NODES=1
          DO WHILE (CURRENT.NE.NQ_START(nr))
            ADJACENT=0
            DO i=-1,1,2
              DO j=1,NXQ(i,0,CURRENT,1)
                ADJACENT=ADJACENT+1
                POINTS(ADJACENT)=NXQ(i,j,CURRENT,1)
                IF (CALCULATED(NXQ(i,j,CURRENT,1))) THEN
!                  CONECT(-1,1,CURRENT)=POINTS(ADJACENT)
                  LABLED=ADJACENT
                ENDIF
              ENDDO
            ENDDO                  
            CONECT(-1,1,CURRENT)=OLD_CURENT !NEW 28/5/99
            CONECT(0,1,CURRENT)=CURRENT
            CONECT(-1,0,CURRENT)=1
            COUNT=0
            DO i=1,ADJACENT
              IF (i.NE.LABLED) THEN
                IF(.NOT.CALCULATED(POINTS(I))) THEN
                  COUNT=COUNT+1
                  CONECT(1,COUNT,CURRENT)=POINTS(I)
                ENDIF
              ENDIF
            ENDDO !i
            IF(ADJACENT.EQ.2) THEN
              IF(.NOT.(CALCULATED(POINTS(1))
     '        .AND.CALCULATED(POINTS(2)))) THEN
C a normal position
                CONECT(1,0,CURRENT)=1                
              ENDIF
            ENDIF !ADJACENT=2
            IF(ADJACENT.EQ.3) THEN  !recalculates NXQ to remove
              CONECT(1,0,CURRENT)=2  !conectivity at diastole
              DO i=1,ADJACENT        !grid points at a bifurcation
                IF(i.NE.LABLED) THEN 
                  nq1=POINTS(i)
                  DO ii=1,ADJACENT
                    IF((ii.NE.i).AND.(ii.NE.LABLED)) THEN
                      nq2=POINTS(ii)
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
              i=-3
              FOUND=.FALSE.
              DO WHILE((.NOT.FOUND).AND.(I.LE.1))
                I=I+2
                COUNT=0
                DO WHILE((COUNT.LT.NXQ(i,0,nq1,1)).AND.(.NOT.FOUND))
                  COUNT=COUNT+1
                  IF(NXQ(i,count,nq1,1).EQ.nq2) THEN
                    FOUND=.TRUE.
                  ENDIF
                ENDDO
              ENDDO
              IF(FOUND) THEN
                IF(COUNT.GE.2) THEN
                  NXQ(i,0,nq1,1)=NXQ(i,0,nq1,1)-1
                ELSE IF(COUNT.EQ.1) THEN
                  NXQ(i,1,nq1,1)=NXQ(i,NXQ(i,0,nq1,1),nq1,1)
                  NXQ(i,0,nq1,1)=NXQ(i,0,nq1,1)-1
                ENDIF
              ENDIF
              i=-3
              FOUND=.FALSE.
              DO WHILE((.NOT.FOUND) .AND.(I.LT.1))
                I=I+2
                COUNT=0
                DO WHILE((COUNT.LT.NXQ(i,0,nq2,1)).AND.(.NOT.FOUND)) 
                  COUNT=COUNT+1
                  IF(NXQ(i,count,nq2,1).EQ.nq1) THEN
                    FOUND=.TRUE.
                  ENDIF
                ENDDO
              ENDDO
              IF(FOUND) THEN
                IF(COUNT.GE.2) THEN
                  NXQ(i,0,nq2,1)=NXQ(i,0,nq2,1)-1
                ELSE IF(COUNT.EQ.1) THEN
                  NXQ(i,1,nq2,1)=NXQ(i,NXQ(i,0,nq2,1),nq2,1)
                  NXQ(i,0,nq2,1)=NXQ(i,0,nq2,1)-1
                ENDIF
              ENDIF !FOUND
            ENDIF !ADJACENT=3
            CALCULATED(CURRENT)=.TRUE.
            NUM_NODES=NUM_NODES+1
            IF((ADJACENT.EQ.2).AND.CALCULATED(POINTS(1))
     '        .AND.CALCULATED(POINTS(2))) THEN
              FOUND=.FALSE.
              DO WHILE ((.NOT.FOUND).AND.(CURRENT.NE.NQ_START(nr)))
                CURRENT=CONECT(-1,1,CURRENT)
                IF(CONECT(1,0,CURRENT).GT.1) THEN
                  IF(.NOT.CALCULATED(CONECT(1,2,CURRENT))) THEN
                    FOUND=.TRUE.
                    OLD_CURENT=CURRENT !NEW 28/5/99
                    CURRENT=CONECT(1,2,CURRENT)
                  ENDIF 
                ENDIF
              ENDDO
            ELSE IF(ADJACENT.EQ.1) THEN         
              FOUND=.FALSE.
              DO WHILE ((.NOT.FOUND).AND.(CURRENT.NE.NQ_START(nr)))
                CURRENT=CONECT(-1,1,CURRENT)
                IF(CONECT(1,0,CURRENT).GT.1) THEN
                  IF (.NOT.CALCULATED(CONECT(1,2,CURRENT))) THEN
                    FOUND=.TRUE.
                    OLD_CURENT=CURRENT !NEW 28/5/99
                    CURRENT=CONECT(1,2,CURRENT)
                  ENDIF 
                ENDIF
              ENDDO
            ELSE
              OLD_CURENT=CURRENT           !NEW 28/5/99
              CURRENT=CONECT(1,1,CURRENT)
            ENDIF !ADJACENT=1
          ENDDO !ADJACENT=2
C Calculating the grid points at bifurcations and the 
C start and end points where boundary condition need to be applied
          no_bc_points=0
          DO nq=NQR(1,nr),NQR(2,nr) !New BC_POINTS array NPS 29/5/99
            ADJACENT=0
            DO ii=-1,1,2
              DO nqq=1,NXQ(ii,0,nq,1)
                ADJACENT=ADJACENT+1
              ENDDO
            ENDDO
            CALL ASSERT(ADJACENT.LE.3,'trifurcation',ERROR,*9999)
            IF((ADJACENT.NE.2)) THEN !normal pt or terminal pt
              no_bc_points=no_bc_points+1
              BC_POINTS(1,1,no_bc_points)=nq
              IF(ADJACENT.NE.3) THEN !not bifurcation
                BRANCH(no_bc_points)=.FALSE.
              ELSE !is bifurcation
                BRANCH(no_bc_points)=.TRUE.
                BC_POINTS(2,1,no_bc_points)=CONECT(1,1,nq) 
                BC_POINTS(3,1,no_bc_points)=CONECT(1,2,nq)
C Determine the a1,b1,c1 for the bifurcation
                DO jj=1,3 !segments
                  nq1=NXQ(1,1,BC_POINTS(jj,1,no_bc_points),1)
                  nq2=NXQ(-1,1,BC_POINTS(jj,1,no_bc_points),1)
                  IF (nq1.NE.0) THEN
                    ne1=NENQ(1,nq1)
                  ELSE
                    ne1=0
                  ENDIF
                  IF (nq2.NE.0) THEN
                    ne2=NENQ(1,nq2)
                  ELSE
                    ne2=0
                  ENDIF
                  ne_current=NENQ(1,BC_POINTS(jj,1,no_bc_points)) 
                  IF(ne1.EQ.ne_current) THEN
                    BC_POINTS(jj,2,no_bc_points)=
     '                NXQ(1,1,BC_POINTS(jj,1,no_bc_points),1) 
                  ELSE                  
                    BC_POINTS(jj,2,no_bc_points)=
     '                NXQ(-1,1,BC_POINTS(jj,1,no_bc_points),1) 
                  ENDIF
                ENDDO !jj
C this determines the a2,b2,c2 for the bifurcation
                DO jj=1,3 !segments
                  nq1=CONECT(-1,1,BC_POINTS(jj,2,no_bc_points))
                  nq2=BC_POINTS(jj,1,no_bc_points)
                  IF (nq1.eq.nq2) THEN
                    BC_POINTS(jj,3,no_bc_points)=nq2 
                  ELSE
                    BC_POINTS(jj,3,no_bc_points)=
     '                BC_POINTS(jj,2,no_bc_points)     
                  ENDIF
                ENDDO !jj
C this determines the half space step grid points a1-a2,b1-b2,c1-c2 
C for the bifurcation

c This is the place to put the connectivity for the branch points
c need to record the colection of three branch points around
c a bifurcation NPS 17/11/96
              ENDIF !ADJACENT.NE.3
            ENDIF !ADJACENT.NE.2
          ENDDO !nq
          BC_POINTS(1,1,0)=no_bc_points
        ENDIF !ITYP3(nr,nx)=1
      ENDIF !ITYP16(nr,nx)=4
C This section of code calculates the conectivity of a 1D branching 
C network in a 3d host mesh, NXQ can not be used directly because
C adjacent elemnts may not have consistant Xi directions.
c      CALL OPENF(IFILE,'DISK','cq.data','UNKNOWN',
c     '  'SEQUEN','FORMATTED',160,ERROR,*9999)         
c      DO nq=NQR(1,nr),NQR(2,nr) !ipinit NPS 4/2/97
c        IF(CONECT(1,0,nq).eq.1) THEN
c          NQ_NEXT=CONECT(1,1,nq)
c          R1=CQ(4,nq)
c          R2=CQ(4,nq_next)
c          IF(((R2/R1).LT.0.99D0).or.((R2/R1).GT.1.001D0)) THEN
c            WRITE(IFILE,*) 'NQ,nq_next,nenq(1,NQ),nenq(1,NQ_next)'
c     '        ,NQ,nq_next,nenq(1,NQ),nenq(1,NQ_next),r2/r1
c            write(IFILE,*) '-----------'
c          ENDIF
c        ENDIF
c      ENDDO 
      DO no_nq=1,BC_POINTS(1,1,0) !now loops have been added this should
        nq=BC_POINTS(1,1,no_nq) ! be rechecked NPS 29/5/99
        IF(BRANCH(no_nq)) THEN
          DO i=1,2
            DISCONT=.FALSE.
            nq1=conect(1,I,nq)
            trace1=YQ(NYNQ(1,nq1,0),2)
            ne1=INT(CQ(8,nq1))
            nq2=conect(1,1,nq1)
            IF(nq2.NE.0) THEN
              trace2=YQ(NYNQ(1,nq2,0),2)
              ne2=INT(cq(8,nq2))
              nq3=CONECT(1,1,nq2)
              trace3=YQ(NYNQ(1,nq3,0),2)
              ne3=INT(cq(8,nq3))
              DO WHILE((CONECT(1,0,nq3).EQ.1).AND.(.NOT.DISCONT))
                IF((ne2.NE.ne1).OR.(ne2.NE.ne3).OR.DISCONT) THEN
                  GRAD_DIF=TRACE1-2.0d0*TRACE2+TRACE3
                  IF((DABS(GRAD_DIF).GT.0.5d0).OR.DISCONT) THEN
                    DISCONT=.TRUE.
                  ENDIF
                ENDIF !ne2
C!!! Something uninitialized here
                IF (((DABS(TRACE1-TRACE2)).GT.2.50D0).OR.
     '            ((DABS(TRACE2-TRACE3)).GT.2.50D0)) THEN
                   DISCONT=.TRUE.
                ENDIF
C the above code removes discontinuities in the trace do to linear 
C interpolation of the host mesh. Once a full C1 conintous host 
C mechanics mesh is implimented it can be removed. NPS 12/11/98
                nq1=nq2
                trace1=trace2
                ne1=ne2
                nq2=nq3
                trace2=trace3
                ne2=ne3
                nq3=conect(1,1,nq2)
                trace3=YQ(NYNQ(1,nq3,0),2)
                ne3=INT(cq(8,nq3))
              ENDDO
              IF(DISCONT) THEN
                nq1=conect(1,I,nq)
                TRACE1=YQ(NYNQ(1,nq1,0),2)
                DO WHILE(CONECT(1,0,nq1).EQ.1)
                  nq1=conect(1,1,nq1)
                  YQ(NYNQ(1,nq1,0),2)=TRACE1                 
                  YQ(NYNQ(4,nq1,0),2)=TRACE1                 
                ENDDO
              ENDIF ! DISCONT
            ENDIF ! (nq2.NE.0) 
          ENDDO !i
        ENDIF !BRANCH(no_nq)
      ENDDO !no_nq
      DO no_nq=1,BC_POINTS(1,1,0) ! This code removes dicontinuities
        nq=BC_POINTS(1,1,no_nq)   ! in the lambda vessel stretch values
        IF(BRANCH(no_nq)) THEN ! between vessel segments NPS 26/11/99
          DO I=1,2             ! checking both branches
            nq1=conect(1,I,nq)
            LAMBDA1=YQ(NYNQ(3,nq1,0),2)
            DO WHILE(CONECT(1,0,nq1).EQ.1) !while not end of element
              nq1=conect(1,1,nq1)
              YQ(NYNQ(3,nq1,0),2)=LAMBDA1              
              YQ(NYNQ(6,nq1,0),2)=LAMBDA1              
            ENDDO
          ENDDO !I=1,2 
        ENDIF !BRANCH(no_nq)
      ENDDO !no_nq
      CALL CLOSEF(IFILE,ERROR,*9999)
C PM 26-JUL-01
C      TIME=T0                                                        
	TIME=T_RESTART
      WORK_PTR=0
      PRINT_FILE=.TRUE.
      TOT_BITIME=0.0
      TOT_BOTIME=0.0
      TOT_GRTIME=0.0
C      FILE_COUNT=1
C PM 26-JUL-01
C	TIME_STEPS=INT((T1-T0)/TINCR)
	TIME_STEPS=INT((TFINISH-TSTART+1.0d-6)/DT)
C PM 26-JUL-01
C      CORO_FILE(1)='1'
C      CORO_FILE(2)='2'
C      CORO_FILE(3)='3'
C      CORO_FILE(4)='4'
C      CORO_FILE(5)='5'
C      CORO_FILE(6)='6'
C      CORO_FILE(7)='7'
C      CORO_FILE(8)='8'
C      CORO_FILE(9)='9'
C      CORO_FILE(10)='b'
C      CORO_FILE(11)='c'
C      CORO_FILE(12)='d'
C      CORO_FILE(13)='e'
C      CORO_FILE(14)='f'
C      CORO_FILE(15)='g'
C      CORO_FILE(16)='h'
C      CORO_FILE(17)='i'
C      CORO_FILE(18)='j'
C      CORO_FILE(19)='k'
C      CORO_FILE(20)='l'
C      CORO_FILE(21)='m'
      DO I=1,50
        PATH_ARRAY(I)=1     
      ENDDO
C      CALL OPENF(IFILE,'DISK','path.data','UNKNOWN',
C     '  'SEQUEN','FORMATTED',160,ERROR,*9999)
C      DO I=1,50
C        READ(IFILE,*) PATH_ARRAY(I)     
C      ENDDO
C      CALL CLOSEF(IFILE,ERROR,*9999)
C PM 26-JUL-01
C      TIME_STEPS=INT((T1-T0)/TINCR)
      UPDATE_VECTOR=.TRUE.
      IF(ITYP16(nr,nx).NE.4) THEN ! not an explicite scheme
        CALL ALLOCATE_MEMORY(NOQT(1,1,nr,nx)*NOQT(2,1,nr,nx),0,CHARTYPE,
     '    WORK_PTR,MEM_INIT,ERROR,*9999)
      ELSE
        CALL ALLOCATE_MEMORY(1,0,CHARTYPE,
     '    WORK_PTR,MEM_INIT,ERROR,*9999)
      ENDIF
      UPDATE_MATRIX=.TRUE.
      CALL CPU_TIMER(CPU_USER,TIME_START2)
C PM 26-JUL-01 : to define boundary conditions from a file
	IF(KTYP3_INIT(nx).EQ.3) THEN
	  CALL STRING_TRIM(FILE03,IBEG,IEND)	
        CALL OPENF(IOFILE1,'TERM',FILE03(IBEG:IEND),'UNKNOWN',
     '    'SEQUEN','UNFORMATTED',132,ERROR,*9999)
      ENDIF
      DO I=1,TIME_STEPS
        IF(ITYP16(nr,nx).eq.4) THEN
          HALF_TIME_STEP=.TRUE.
          UPDATE_VECTOR=.TRUE.
C PM 26-JUL-01
C          TIME=TIME+(TINCR/2.0d0) !half time step
           TIME=TIME+(DT/2.0d0)    !half time step
           CALL ASSEMBLE9(BC_POINTS,BRANCH,CQ,CONECT,ER,ES,GKK,GRR,
     '       HALF_TIME_STEP,ISC_GKK,ISR_GKK,LGE,1,1,nr,nx,NHQ,
     '       NQ_START(nr),NXQ,NYNQ,NYQNR,UPDATE_MATRIX,UPDATE_VECTOR,
     '       %VAL(WORK_PTR),TIME,TOT_BITIME,TOT_BOTIME,TOT_GRTIME,XQ,YQ,
     '       TIME_VALUES,TIME_VARAIBLE_NAMES,NTIME_POINTS,NTIME_INTERP,
     '       ERROR,*9999)
C no call to solve is needed for an explicite scheme, in this
C case the right hand side is written directly into YQ
          HALF_TIME_STEP=.FALSE.
C PM 26-JUL-01
C          TIME=TIME+(TINCR/2.0d0) !full time step
          TIME=TIME+(DT/2.0d0)                                          
        ENDIF
        DO no_nynr=1,NYQNR(0,0,1,nr) 
          ny=NYQNR(no_nynr,0,1,nr)
          YQ(ny,8)=YQ(ny,1) !temporary storage for previous time step
        ENDDO
        UPDATE_VECTOR=.TRUE.
        CALL ASSEMBLE9(BC_POINTS,BRANCH,CQ,CONECT,ER,ES,GKK,GRR,
     '    HALF_TIME_STEP,ISC_GKK,ISR_GKK,LGE,1,1,nr,nx,NHQ,
     '    NQ_START(nr),
     '    NXQ,NYNQ,NYQNR,UPDATE_MATRIX,UPDATE_VECTOR,%VAL(WORK_PTR),
     '    TIME,TOT_BITIME,TOT_BOTIME,TOT_GRTIME,XQ,YQ,
     '    TIME_VALUES,TIME_VARAIBLE_NAMES,NTIME_POINTS,
     '    NTIME_INTERP,             
     '    ERROR,*9999)                                                
C no call to solve is needed for an explicite scheme, in this
C case the right hand side is written directly into YQ
        HALF_TIME_STEP=.FALSE.
        IF(ITYP16(nr,nx).eq.4) THEN !Lax-Wendroff  
          IF(ITYP3(nr,nx).eq.1) THEN !flow in elastic tube 
            DATA_COUNT=DATA_COUNT+1
            IF(DATA_COUNT.GE.((T1-T0)/(20*TINCR))) THEN
C              WRITE(OP_STRING,'(''TIME'',F12.8)') TIME 
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DATA_COUNT=0
              IF (DATA_FILE_COUNT.EQ.0) THEN
                FILE_COUNT=FILE_COUNT+1
C PM 26-JUL-01 : This is now done using exgrid
C                CALL OPENF(IOFILE1,'DISK','vel.data'//
C     '            CORO_FILE(FILE_COUNT),'NEW',
C     '            'SEQUEN','FORMATTED',132,ERROR,*9999)
C                CALL OPENF(IOFILE2,'DISK','pres.data'//
C     '            CORO_FILE(FILE_COUNT),'NEW',
C     '            'SEQUEN','FORMATTED',132,ERROR,*9999)
C                CALL OPENF(IOFILE3,'DISK','flow.data'//
C     '            CORO_FILE(FILE_COUNT),'NEW',
C     '            'SEQUEN','FORMATTED',132,ERROR,*9999)
C                CALL OPENF(IOFILE4,'DISK','rad.data'//
C     '            CORO_FILE(FILE_COUNT),'NEW',
C     '            'SEQUEN','FORMATTED',132,ERROR,*9999)
C                CALL OPENF(IOFILE5,'DISK',CORO_FILE(FILE_COUNT)//
C     '            '.exnode',
C     '            'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
C                CALL OPENF(IOFILE6,'DISK',CORO_FILE(FILE_COUNT)//
C     '            'v.exnode',
C     '            'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
                PRINT_FILE=.TRUE.
              ELSE IF(PRINT_FILE) THEN !only one file per time step
C MHT 9-8-1 unused
c                IND1=(DATA_FILE_COUNT/10)+1                
c                IND2=(MOD(DATA_FILE_COUNT,10))+1
C                CALL OPENF(IOFILE1,'DISK','veld.data'//
C     '            FILE_STRING(IND1:IND1)//
C     '            FILE_STRING(IND2:IND2),'NEW',
C     '            'SEQUEN','FORMATTED',132,ERROR,*9999)
C                CALL OPENF(IOFILE2,'DISK','presd.data'//
C     '            FILE_STRING(IND1:IND1)//
C     '            FILE_STRING(IND2:IND2),'NEW',
c     '            'SEQUEN','FORMATTED',132,ERROR,*9999)
C                CALL OPENF(IOFILE3,'DISK','flowd.data'//
C     '            FILE_STRING(IND1:IND1)//
C     '            FILE_STRING(IND2:IND2),'NEW',
C     '            'SEQUEN','FORMATTED',132,ERROR,*9999)
C                CALL OPENF(IOFILE4,'DISK','radd.data'//
C     '            FILE_STRING(IND1:IND1)//
C     '            FILE_STRING(IND2:IND2),'NEW',
C     '            'SEQUEN','FORMATTED',132,ERROR,*9999)
C                CALL OPENF(IOFILE5,'DISK', FILE_STRING(IND1:IND1)//
C     '            FILE_STRING(IND2:IND2)//
C     '            'd.exnode',
C     '            'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
C                CALL OPENF(IOFILE6,'DISK', FILE_STRING(IND1:IND1)//
C     '            FILE_STRING(IND2:IND2)//
C     '            'vd.exnode',
C     '            'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)
              ENDIF
              IF(PRINT_FILE) THEN
                PRINT_FILE=.FALSE.
                nq=NQ_START(nr)
                LAMBDA=(YQ(NYNQ(3,nq,0),9)*(T1-TIME)/(T1-T0))
     '                +(((TIME-T0)/(T1-T0))*
     '                YQ(NYNQ(3,nq,0),2))
                COUNT=1
                LEN_COUNT=0.0d0
                nj=NJ_LOC(NJL_FIEL,1,nr)
                DONE=.FALSE.
                BIFUR_COUNT=1
                DO while(.NOT.DONE)
                  ny=NYNQ(3,nq,0)
                  WRITE(IOFILE1,FMT) LEN_COUNT,YQ(ny,1),
     '              (-1.0d0*YQ((ny+3),1)),LAMBDA
                  IF(conect(1,1,nq).EQ.0) THEN
                    DONE=.TRUE.
                  ELSE
                    nq_old=nq
                    COUNT=COUNT+1
                    IF (CONECT(1,0,nq).GT.1)THEN
                      COUNT=COUNT-1          
                      ne1=NENQ(1,CONECT(1,1,nq))
                      ne2=NENQ(1,CONECT(1,2,nq))
                      nb1=NBJ(1,ne1)
                      nb2=NBJ(1,ne2)
                      IF(XP(1,1,nj,NPNE(1,nb1,ne1)).GT.
     '                  XP(1,1,nj,NPNE(NNT(nb1),nb1,ne1))) THEN
                        np1=NPNE(NNT(nb1),nb1,ne1)
                      ELSE
                        np1=NPNE(1,nb1,ne1)
                      ENDIF
                      IF(XP(1,1,nj,NPNE(1,nb2,ne2)).GT.
     '                  XP(1,1,nj,NPNE(NNT(nb2),nb2,ne2))) THEN
                        np2=NPNE(NNT(nb2),nb2,ne2)
                      ELSE
                        np2=NPNE(1,nb2,ne2)
                      ENDIF
                      IF(PATH_ARRAY(BIFUR_COUNT).EQ.1) THEN
                        IF(XP(1,1,nj,np1).GT.XP(1,1,nj,np2)) THEN
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE1) THEN
                            nq=CONECT(1,1,nq)
                          ELSE
                            nq=CONECT(1,2,nq)
                          ENDIF
                        ELSE
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE2) THEN
                            nq=CONECT(1,1,nq)
                          ELSE
                            nq=CONECT(1,2,nq)
                          ENDIF
                        ENDIF
                      ELSE
                        IF(XP(1,1,nj,np1).GT.XP(1,1,nj,np2)) THEN
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE1) THEN
                            nq=CONECT(1,2,nq)
                          ELSE
                            nq=CONECT(1,1,nq)
                          ENDIF
                        ELSE
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE2) THEN
                            nq=CONECT(1,2,nq)
                          ELSE
                            nq=CONECT(1,1,nq)
                          ENDIF
                        ENDIF 
                      ENDIF
                      BIFUR_COUNT=BIFUR_COUNT+1
                    ELSE
                      nq=CONECT(1,1,nq)
                    ENDIF
                    STEP=0.0d0
C PM 26-JUL-01
C                    LAMBDA=(YQ(NYNQ(3,nq,0),9)*(T1-TIME)/(T1-T0))
C     '                +(((TIME-T0)/(T1-T0))*
C     '                YQ(NYNQ(3,nq,0),2))
                    LAMBDA=(YQ(NYNQ(3,nq,0),9)*
     '                (TFINISH-TIME)/(TFINISH-TSTART))
     '                +(((TIME-TSTART)/(TFINISH-TSTART))*
     '                YQ(NYNQ(3,nq,0),2))
                    DO njj=1,3
                      STEP=STEP+((XQ(njj,nq)-XQ(njj,nq_old))**2.0d0)    
                    ENDDO
                    STEP=STEP*LAMBDA
                    LEN_COUNT=LEN_COUNT+(STEP**0.5d0)
                  ENDIF
                ENDDO
                nq=NQ_START(nr)
                COUNT=1
                LEN_COUNT=0.0d0
                nj=NJ_LOC(NJL_FIEL,1,nr)
                DONE=.FALSE.
                BIFUR_COUNT=1
                DO while(.NOT.DONE)                 
                  ny=NYNQ(1,nq,0)
                  WRITE(IOFILE2,FMT) LEN_COUNT,YQ(ny,1),YQ(ny,2),
     '              YQ(ny+3,1)
                  FA=PI*(YQ(NYNQ(2,nq,0),1)**2.0d0)*
     '              YQ(NYNQ(3,nq,0),1)  
                  FB=-PI*(YQ(NYNQ(5,nq,0),1)**2.0d0)*
     '              YQ(NYNQ(6,nq,0),1)
                  WRITE(IOFILE3,FMT) LEN_COUNT,FA,FB,0.0D0
                  IF(conect(1,1,nq).EQ.0) THEN
                    FA=PI*(YQ(NYNQ(2,nq,0),1)**2.0d0)*
     '                YQ(NYNQ(3,nq,0),1)  
C MHT 9-8-1 P1 set but not used
c                    P1=YQ(ny,1)-YQ(NYNQ(2,nq,0),4)*FA
                    FB=-PI*(YQ(NYNQ(5,nq,0),1)**2.0d0)*
     '                YQ(NYNQ(6,nq,0),1)  
C MHT 9-8-1 P2 set but not used
c                    P2=YQ((ny+3),1)+YQ(NYNQ(5,nq,0),4)*FB
C PM 26-JUL-01                    
C                    WRITE(IOFILE2,FMT) LEN_COUNT,P1,YQ(ny,2),
C     '                P2
C                    WRITE(IOFILE3,FMT) LEN_COUNT,
C     '                YQ(NYNQ(3,nq,0),4),YQ(NYNQ(3,nq,0),4)
                    DONE=.TRUE.
                  ELSE
                    nq_old=nq
                    COUNT=COUNT+1
                    IF (CONECT(1,0,nq).GT.1)THEN
                      COUNT=COUNT-1          
                      ne1=NENQ(1,CONECT(1,1,nq))
                      ne2=NENQ(1,CONECT(1,2,nq))
                      nb1=NBJ(1,ne1)
                      nb2=NBJ(1,ne2)
                      IF(XP(1,1,nj,NPNE(1,nb1,ne1)).GT.
     '                  XP(1,1,nj,NPNE(NNT(nb1),nb1,ne1))) THEN
                        np1=NPNE(NNT(nb1),nb1,ne1)
                      ELSE
                        np1=NPNE(1,nb1,ne1)
                      ENDIF
                      IF(XP(1,1,nj,NPNE(1,nb2,ne2)).GT.
     '                  XP(1,1,nj,NPNE(NNT(nb2),nb2,ne2))) THEN
                        np2=NPNE(NNT(nb2),nb2,ne2)
                      ELSE
                        np2=NPNE(1,nb2,ne2)
                      ENDIF
                      IF(PATH_ARRAY(BIFUR_COUNT).EQ.1) THEN
                        IF(XP(1,1,nj,np1).GT.XP(1,1,nj,np2)) THEN
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE1) THEN
                            nq=CONECT(1,1,nq)
                          ELSE
                            nq=CONECT(1,2,nq)
                          ENDIF
                        ELSE
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE2) THEN
                            nq=CONECT(1,1,nq)
                          ELSE
                            nq=CONECT(1,2,nq)
                          ENDIF
                        ENDIF
                      ELSE
                        IF(XP(1,1,nj,np1).GT.XP(1,1,nj,np2)) THEN
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE1) THEN
                            nq=CONECT(1,2,nq)
                          ELSE
                            nq=CONECT(1,1,nq)
                          ENDIF
                        ELSE
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE2) THEN
                            nq=CONECT(1,2,nq)
                          ELSE
                            nq=CONECT(1,1,nq)
                          ENDIF
                        ENDIF 
                      ENDIF
                      BIFUR_COUNT=BIFUR_COUNT+1
                    ELSE
                      nq=CONECT(1,1,nq)
                    ENDIF
                    STEP=0.0d0
C PM 26-JUL-01
C                   LAMBDA=(YQ(NYNQ(3,nq,0),9)*(T1-TIME)/(T1-T0))
C     '                +(((TIME-T0)/(T1-T0))*
C     '                YQ(NYNQ(3,nq,0),2))
                   LAMBDA=(YQ(NYNQ(3,nq,0),9)*
     '               (TFINISH-TIME)/(TFINISH-TSTART))
     '                +(((TIME-TSTART)/(TFINISH-TSTART))*
     '                YQ(NYNQ(3,nq,0),2))
                    DO njj=1,3
                      STEP=STEP+((XQ(njj,nq)-XQ(njj,nq_old))**2.0d0)
                    ENDDO
                    STEP=STEP*LAMBDA
                    LEN_COUNT=LEN_COUNT+(STEP**0.5d0)
                  ENDIF
                ENDDO
                nq=NQ_START(nr)
                COUNT=1
                LEN_COUNT=0.0d0
                nj=NJ_LOC(NJL_FIEL,1,nr)
                DONE=.FALSE.
                BIFUR_COUNT=1
                DO while(.NOT.DONE)
                  ny=NYNQ(2,nq,0)
                  YQ(ny,6)=1.0d0 !labels grid points along a path
C                  WRITE(IOFILE4,*) LEN_COUNT,YQ(ny,1),YQ((ny+3),1)
                  IF(conect(1,1,nq).EQ.0) THEN
                    DONE=.TRUE.
                  ELSE
                    nq_old=nq
                    COUNT=COUNT+1
                    IF (CONECT(1,0,nq).GT.1)THEN
                      COUNT=COUNT-1          
                      ne1=NENQ(1,CONECT(1,1,nq))
                      ne2=NENQ(1,CONECT(1,2,nq))
                      nb1=NBJ(1,ne1)
                      nb2=NBJ(1,ne2)
                      IF(XP(1,1,nj,NPNE(1,nb1,ne1)).GT.
     '                  XP(1,1,nj,NPNE(NNT(nb1),nb1,ne1))) THEN
                        np1=NPNE(NNT(nb1),nb1,ne1)
                      ELSE
                        np1=NPNE(1,nb1,ne1)
                      ENDIF
                      IF(XP(1,1,nj,NPNE(1,nb2,ne2)).GT.
     '                  XP(1,1,nj,NPNE(NNT(nb2),nb2,ne2))) THEN
                        np2=NPNE(NNT(nb2),nb2,ne2)
                      ELSE
                        np2=NPNE(1,nb2,ne2)
                      ENDIF
                      IF(PATH_ARRAY(BIFUR_COUNT).EQ.1) THEN
                        IF(XP(1,1,nj,np1).GT.XP(1,1,nj,np2)) THEN
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE1) THEN
                            nq=CONECT(1,1,nq)
                          ELSE
                            nq=CONECT(1,2,nq)
                          ENDIF
                        ELSE
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE2) THEN
                            nq=CONECT(1,1,nq)
                          ELSE
                            nq=CONECT(1,2,nq)
                          ENDIF
                        ENDIF
                      ELSE
                        IF(XP(1,1,nj,np1).GT.XP(1,1,nj,np2)) THEN
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE1) THEN
                            nq=CONECT(1,2,nq)
                          ELSE
                            nq=CONECT(1,1,nq)
                          ENDIF
                        ELSE
                          IF(NENQ(1,CONECT(1,1,nq)).EQ.NE2) THEN
                            nq=CONECT(1,2,nq)
                          ELSE
                            nq=CONECT(1,1,nq)
                          ENDIF
                        ENDIF 
                      ENDIF
                      BIFUR_COUNT=BIFUR_COUNT+1
                    ELSE
                      nq=CONECT(1,1,nq)
                    ENDIF
                    STEP=0.0d0
C PM 26-JUL-01
C                    LAMBDA=(YQ(NYNQ(3,nq,0),9)*(T1-TIME)/(T1-T0))
C     '                +(((TIME-T0)/(T1-T0))*
C     '                YQ(NYNQ(3,nq,0),2))
                    LAMBDA=(YQ(NYNQ(3,nq,0),9)*(TFINISH-TIME)/
     '                (TFINISH-TSTART))
     '                +(((TIME-TSTART)/(TFINISH-TSTART))*
     '                YQ(NYNQ(3,nq,0),2))
                    DO njj=1,3
                      STEP=STEP+((XQ(njj,nq)-XQ(njj,nq_old))**2.0d0)
                    ENDDO
                    STEP=STEP*LAMBDA
                    LEN_COUNT=LEN_COUNT+(STEP**0.5d0)
                  ENDIF
                ENDDO   
C MHT 9-8-1 unused
C                OFFSET=10000
C PM 26-JUL-01              
C                WRITE(IOFILE5,'(1X,''Group name: grid1'')') 
C                WRITE(IOFILE5,'(1X,''#Fields=5'')') 
C                WRITE(IOFILE5,'(1X,''1) coordinates, coordinate, '
C     '            //'rectangular cartesian, #Components=3'')')
C                WRITE(IOFILE5,'(1X,''  x.  Value index= 1, '
C     '            //'#Derivatives=0'')')
C                WRITE(IOFILE5,'(1X,''  y.  Value index= 2, '
C     '            //'#Derivatives=0'')')
C                WRITE(IOFILE5,'(1X,''  z.  Value index= 3, '
C     '            //'#Derivatives=0'')')               
C                WRITE(IOFILE5,'(1X,''2) pressure, field, '
C     '            //'rectangular cartesian, #Components=1'')')
C                WRITE(IOFILE5,'(1X,''  pressure.  '
C     '            //'Value index= 4, #Derivatives=0'')')
C                WRITE(IOFILE5,'(1X,''3) radius, field, '
C     '            //'rectangular cartesian, #Components=1'')')
C                WRITE(IOFILE5,'(1X,''  radius.  '
C     '            //'Value index= 5, #Derivatives=0'')')
C                WRITE(IOFILE5,'(1X,''4) velocity, field, '
C     '            //'rectangular cartesian, #Components=1'')')
C                WRITE(IOFILE5,'(1X,''  velocity.  '
C     '            //'Value index= 6, #Derivatives=0'')')
C                WRITE(IOFILE5,'(1X,''5) path, field, '
C     '            //'rectangular cartesian, #Components=1'')')
C                WRITE(IOFILE5,'(1X,''  path.  '
C     '            //'Value index= 7, #Derivatives=0'')')
                DO nq=NQR(1,nr),NQR(2,nr) 
C                  WRITE(IOFILE5,'(1X,''Node: '',I7)') nq+OFFSET
                  ny_p=NYNQ(1,nq,0)
                  ny_r=NYNQ(2,nq,0)
                  ny_v=NYNQ(3,nq,0)
                  VEL=YQ(ny_v,1)
                  IF(DABS(VEL).LT.LOOSE_TOL) THEN
                    VEL=LOOSE_TOL
                  ENDIF
c MHT 16-06-99 PATH set but not used, 9-8-1 trace set not used
c                  PATH=MAX(YQ(ny_r,6),LOOSE_TOL)
c                  TRACE=YQ(NYNQ(1,nq,0),2)
C PM 26-JUL-01
C                  WRITE(IOFILE5,'(1X,7E13.5)') 
C     '              XQ(1,nq),XQ(2,nq),XQ(3,nq),
C     '              YQ(ny_p,1),YQ(ny_r,1),VEL,TRACE
                ENDDO
C                WRITE(IOFILE6,'(1X,''Group name: grid1'')') 
C                WRITE(IOFILE6,'(1X,''#Fields=5'')') 
C                WRITE(IOFILE6,'(1X,''1) coordinates, coordinate, '
C     '            //'rectangular cartesian, #Components=3'')')
C                WRITE(IOFILE6,'(1X,''  x.  Value index= 1, '
C     '            //'#Derivatives=0'')')
C                WRITE(IOFILE6,'(1X,''  y.  Value index= 2, '
C     '            //'#Derivatives=0'')')
C                WRITE(IOFILE6,'(1X,''  z.  Value index= 3, '
C     '            //'#Derivatives=0'')')               
C                WRITE(IOFILE6,'(1X,''2) pressure, field, '
C     '            //'rectangular cartesian, #Components=1'')')
C                WRITE(IOFILE6,'(1X,''  pressure.  '
C     '            //'Value index= 4, #Derivatives=0'')')
C                WRITE(IOFILE6,'(1X,''3) radius, field, '
C     '            //'rectangular cartesian, #Components=1'')')
C                WRITE(IOFILE6,'(1X,''  radius.  '
C     '            //'Value index= 5, #Derivatives=0'')')
C                WRITE(IOFILE6,'(1X,''4) velocity, field, '
C     '            //'rectangular cartesian, #Components=1'')')
C                WRITE(IOFILE6,'(1X,''  velocity.  '
C     '            //'Value index= 6, #Derivatives=0'')')
C                WRITE(IOFILE6,'(1X,''5) path, field, '
C     '            //'rectangular cartesian, #Components=1'')')
C                WRITE(IOFILE6,'(1X,''  path.  '
C     '            //'Value index= 7, #Derivatives=0'')')
                DO nq=NQR(1,nr),NQR(2,nr) 
C                  WRITE(IOFILE6,'(1X,''Node: '',I7)') nq+OFFSET
                  ny_p=NYNQ(4,nq,0)
                  ny_r=NYNQ(5,nq,0)
                  ny_v=NYNQ(6,nq,0)                 
                  VEL=YQ(ny_v,1)
                  IF(DABS(VEL).LT.LOOSE_TOL) THEN
                    VEL=LOOSE_TOL
                  ENDIF
C MHT 9-8-1 TRACE set not used
c                  TRACE=YQ(NYNQ(1,nq,0),2)
                  LAMBDA=YQ(NYNQ(3,nq,0),2)

C                  WRITE(IOFILE6,'(1X,7E13.5)') 
C     '              XQ(1,nq),XQ(2,nq),XQ(3,nq),
C     '              YQ(ny_p,1),YQ(ny_r,1),VEL,LAMBDA
                ENDDO
C                CALL CLOSEF(IOFILE1,ERROR,*9999) 
C                CALL CLOSEF(IOFILE2,ERROR,*9999) 
C                CALL CLOSEF(IOFILE3,ERROR,*9999) 
C                CALL CLOSEF(IOFILE4,ERROR,*9999) 
C                CALL CLOSEF(IOFILE5,ERROR,*9999) 
C               CALL CLOSEF(IOFILE6,ERROR,*9999) 
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO !time
C PM 26-JUL-01
	T_RESTART=TIME
      WRITE(OP_STRING,'(''TIME'',F12.8,'' seconds'')') TIME             
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)                           
      CALL CPU_TIMER(CPU_USER,TIME_STOP)
C MHT 9-8-1 unused
c      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
C PM 26-JUL-01
C      WRITE(OP_STRING,'(/''total CPU time for grid '
C     '  //'solution: '',D15.7,'' s'')') ELAPSED_TIME
C      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      WRITE(OP_STRING,'(/''total CPU time for grid '
C     '  //'assemble: '',D15.7,'' s'')') TOT_GRTIME
C      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      WRITE(OP_STRING,'(/''total CPU time for grid '
C     '  //'bifurcation solution: '',D15.7,'' s'')') TOT_BITIME
C      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      WRITE(OP_STRING,'(/''total CPU time for grid '
C     '  //'black box solution: '',D15.7,'' s'')') TOT_BOTIME
C      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      IF(ITYP16(nr,nx).eq.4) THEN !Lax-Wendroff  
        IF(ITYP3(nr,nx).eq.1) THEN !flow in elastic tube 
          nj=NJ_LOC(NJL_FIEL,2,nr) !allocates grid point information
          nj1=NJ_LOC(NJL_FIEL,3,nr) ! to nodal fields
          nj2=NJ_LOC(NJL_FIEL,1,nr)
          DO no_ne=1,NEELEM(0,nr)
            ne1=NEELEM(no_ne,nr)
            nb1=NBJ(1,ne1)
            np1=NPNE(1,nb1,ne1)
            np2=NPNE(NNT(nb1),nb1,ne1)
            nq1=NQNE(ne1,1)
            nq1_pre=CONECT(-1,1,NQ1)
            nq2=NQNE(ne1,NQET(NQS(ne1)))
            nq2_pre=CONECT(-1,1,NQ2)
            XI_DIST1=(DABS(XIP(1,NP1)-YQ(NYNQ(1,nq1,0),10)))
     '        +(DABS(XIP(2,NP1)-YQ(NYNQ(2,nq1,0),10)))
     '        +(DABS(XIP(3,NP1)-YQ(NYNQ(3,nq1,0),10)))
            XI_DIST2=(DABS(XIP(1,NP2)-YQ(NYNQ(1,nq2,0),10)))
     '        +(DABS(XIP(2,NP2)-YQ(NYNQ(2,nq2,0),10)))
     '        +(DABS(XIP(3,NP2)-YQ(NYNQ(3,nq2,0),10)))
            XI_SUM1=DABS(YQ(NYNQ(1,nq1,0),10))
     '        +DABS(YQ(NYNQ(2,nq1,0),10))
     '        +DABS(YQ(NYNQ(3,nq1,0),10))
            XI_SUM2=DABS(YQ(NYNQ(1,nq2,0),10))
     '        +DABS(YQ(NYNQ(2,nq2,0),10))
     '        +DABS(YQ(NYNQ(3,nq2,0),10))
C could use radius for this loop as well
            IF((XI_SUM2+XI_SUM1).GT.LOOSE_TOL) THEN
              IF((XI_DIST1+XI_DIST2).GT.LOOSE_TOL) THEN
                np2=NPNE(1,nb1,ne1) 
                np1=NPNE(NNT(nb1),nb1,ne1) 
              ENDIF
            ENDIF
            IF(NQ1_PRE.EQ.0) THEN
              ny_p=NYNQ(3,nq1,0)
              ny_r=NYNQ(2,nq1,0)
              ny_v=NYNQ(6,nq1,0) !now doing venous pressure 
              IF(DABS(YQ(ny_p,1)).GT.LOOSE_TOL) THEN
                XP(1,1,nj,np1)=YQ(ny_p,1) 
              ELSE
                XP(1,1,nj,np1)=LOOSE_TOL
              ENDIF
              XP(1,1,nj2,np1)=YQ(ny_r,3) !radius left unchanged
              IF(DABS(YQ(ny_v,1)).GT.LOOSE_TOL) THEN
                XP(1,1,nj1,np1)=YQ(ny_v,1)
              ELSE
                XP(1,1,nj1,np1)=LOOSE_TOL
              ENDIF
            ELSE  
              IF(CONECT(1,0,NQ1_pre).LT.2) THEN 
C not a distal bifuraction
                ny_p=NYNQ(3,nq1,0)
                ny_r=NYNQ(2,nq1,0)
                ny_v=NYNQ(6,nq1,0) !now doing venous pressue
                IF(DABS(YQ(ny_p,1)).GT.LOOSE_TOL) THEN
                  XP(1,1,nj,np1)=YQ(ny_p,1) 
                ELSE
                  XP(1,1,nj,np1)=LOOSE_TOL
                ENDIF
C                if(YQ(ny_r,1).lt.XP(1,1,nj2,np1)) then
C                  write(*,*) ne1,np1,nq1,nq1_pre,CONECT(1,0,NQ1_pre)
C                  ny_p=ny_p
C                endif
                XP(1,1,nj2,np1)=YQ(ny_r,3) !radius left unchanged
                IF(DABS(YQ(ny_v,1)).GT.LOOSE_TOL) THEN
                  XP(1,1,nj1,np1)=YQ(ny_v,1)
                ELSE
                  XP(1,1,nj1,np1)=LOOSE_TOL
                ENDIF
              ENDIF
            ENDIF
            IF(NQ2_PRE.EQ.0) THEN
              ny_p=NYNQ(3,nq2,0)
              ny_r=NYNQ(2,nq2,0)
              ny_v=NYNQ(6,nq2,0) !now doing venous pressure
              IF(DABS(YQ(ny_p,1)).GT.LOOSE_TOL) THEN
                XP(1,1,nj,np2)=YQ(ny_p,1) 
              ELSE
                XP(1,1,nj,np2)=LOOSE_TOL
              ENDIF
              XP(1,1,nj2,np2)=YQ(ny_r,3) !radius left unchanged
              IF(DABS(YQ(ny_v,1)).GT.LOOSE_TOL) THEN
                XP(1,1,nj1,np2)=YQ(ny_v,1)
              ELSE
                XP(1,1,nj1,np2)=LOOSE_TOL
              ENDIF
            ELSE
              IF(CONECT(1,0,NQ2_pre).LT.2) THEN 
C not a distal bifuraction
                ny_p=NYNQ(3,nq2,0)
                ny_r=NYNQ(2,nq2,0)
                ny_v=NYNQ(6,nq2,0)  !now doing  VENOUS PRESSURE 
                IF(DABS(YQ(ny_p,1)).GT.LOOSE_TOL) THEN
                  XP(1,1,nj,np2)=YQ(ny_p,1) 
                ELSE
                  XP(1,1,nj,np2)=LOOSE_TOL
                ENDIF
                XP(1,1,nj2,np2)=YQ(ny_r,3)
                IF(DABS(YQ(ny_v,1)).GT.LOOSE_TOL) THEN
                  XP(1,1,nj1,np2)=YQ(ny_v,1)
                ELSE
                  XP(1,1,nj1,np2)=LOOSE_TOL
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDIF
      CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)

      IF(DATA_FILE_COUNT.NE.0) THEN
        DATA_FILE_COUNT=DATA_FILE_COUNT+1
      ENDIF
 
      CALL EXITS('MARCH4')
      RETURN
 9999 CALL ERRORS('MARCH4',ERROR)
      CALL EXITS('MARCH4')
      RETURN 1
      END



      SUBROUTINE MARCH19(ADAMS_IWORK,ICQS,ICQS_SPATIAL,IICQS_SPATIAL,
     '  IRCQS_SPATIAL,LD,NBJ,NDDATA,nq,NRLIST,NXLIST,ADAMS_WORK,RCQS,
     '  RCQS_SPATIAL,WD,XID,YQS,ZD,RHS_ROUTINE,ERROR,*)

C**** Subroutine: MARCH19
C***  Description:
C***    MARCH19 is used to solve time dependent cellular models from
C***    FE19. ZD(3,nd) is used to store each solution variable to be
C***    passed into IOSIGN to be written to the signal file specified
C***    in the .ipsolv file, which also specifies whether the derived
C***    values are to be written to the file. NDDATA(0,0) stores the 
C***    total number of output variables (treated as data points).

C***  Created by David Nickerson, May 1999

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:adam00.cmn'
c      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cell00.cmn'
      INCLUDE 'cmiss$reference:cell02.cmn'
      INCLUDE 'cmiss$reference:cell_reserved.inc'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
c      INCLUDE 'cmiss$reference:grid00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
c      INCLUDE 'cmiss$reference:iwrit00.cmn'
c      INCLUDE 'cmiss$reference:ktyp30.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
      INCLUDE 'cmiss$reference:maqloc00.inc'
c      INCLUDE 'cmiss$reference:marc00.cmn'
      INCLUDE 'cmiss$reference:nqloc00.inc'
      INCLUDE 'cmiss$reference:sign00.cmn'
c      INCLUDE 'cmiss$reference:solv00.cmn'
      INCLUDE 'cmiss$reference:time02.cmn'
c      INCLUDE 'cmiss$reference:tol00.cmn'

!     Parameter list
      INTEGER ADAMS_IWORK(ADAMS_LIWORK),ICQS(NQIM),
     '  ICQS_SPATIAL(NQISVM,NQM),IICQS_SPATIAL(0:NQISVM),
     '  IRCQS_SPATIAL(0:NQRSVM),LD(NDM),NBJ(NJM,NEM),
     '  NDDATA(0:NDM,0:NRM),nq,NRLIST(0:NRM),NXLIST(0:NXM)
      REAL*8 ADAMS_WORK(ADAMS_LWORK),RCQS(NQRM),
     '  RCQS_SPATIAL(NQRSVM,NQM),WD(NJM,NDM),XID(NIM,NDM),
     '  YQS(NIQSM,NQM),ZD(NJM,NDM)
      EXTERNAL RHS_ROUTINE
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER IFAIL,nd,NUMTIMEDATA,NUMTIMEDATA1,nxc,nx,SIZES(10)
      REAL*8 ENDT,SIGNALMAX(9),SIGNALMIN(9),STARTT,T,TIME
      LOGICAL ENDFILE,STIFF_EQNS
      CHARACTER FILEFORMAT*6

      LOGICAL EXTEND_INTERVAL
      PARAMETER(EXTEND_INTERVAL=.TRUE.)

      INTEGER L_DY
      PARAMETER(L_DY=100)
      REAL*8 DY(L_DY)

      CALL ENTERS('MARCH19',*9999)
 
      !set the tabulation time step
      DT = TINCR
      !Get class information
      nxc=NXLIST(1)
      CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
      CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '  ERROR,*9999)
      CALL ASSERT(CALL_CELL,'>>Need to define cell parameters',
     '  ERROR,*9999)
      CALL ASSERT(CALL_CELL_MATE,
     '  '>>Need to define cell material parameters',ERROR,*9999)
C *** Initialise the model
      IF (CELL_SPATIALLY_VARYING) CALL CELL_ASSIGN_SPATIAL(ICQS,
     '  ICQS_SPATIAL,IICQS_SPATIAL,IRCQS_SPATIAL,nq,RCQS,RCQS_SPATIAL,
     '  ERROR,*9999)
      SIZES(NUM_EQN)=CELL_NUM_STATE
      SIZES(NUM_PARAM)=CELL_NUM_PARAMETERS
      SIZES(NUM_AII)=CELL_NUM_AII
      SIZES(NUM_AIO)=CELL_NUM_AIO
      SIZES(NUM_ARI)=CELL_NUM_ARI
      SIZES(NUM_ARO)=CELL_NUM_ARO
      SIZES(NUM_CONTROL)=CELL_NUM_CONTROL
      SIZES(NUM_DERIVED)=CELL_NUM_DERIVED
      SIZES(NUM_MODEL)=CELL_NUM_MODEL
      SIZES(NUM_PROTOCOL)=CELL_NUM_PROTOCOL
      !get the solution signal file type
      IF(BINTIMEFILE.EQ.1) THEN
        FILEFORMAT='BINARY'
      ELSE
        FILEFORMAT='ASCII'
      ENDIF
C *** Set the number of data points required
      IF(OUTPUT_DERIVED) THEN
        NDDATA(0,0)=CELL_NUM_STATE+CELL_NUM_DERIVED
        NDDATA(0,1)=NDDATA(0,0)
      ELSE
        NDDATA(0,0)=CELL_NUM_STATE
        NDDATA(0,1)=CELL_NUM_STATE
      ENDIF
C *** Check that the maximum number of data points is greater than the
C     number of "signals" required
      CALL ASSERT(NDDATA(0,0).LE.NDM,'>>Increase NDM',ERROR,*9999)
      DO nd=1,NDDATA(0,0)
        WD(1,nd)=1.0d0
        WD(2,nd)=1.0d0
        WD(3,nd)=1.0d0
        ZD(1,nd)=0.0d0
        ZD(2,nd)=0.0d0
        XID(1,nd)=0.0d0
        XID(2,nd)=0.0d0
        LD(nd)=1
        NDDATA(nd,1)=nd
      ENDDO
      SIGNAL_REGNAME(1,IOFILE1)=''
      SIGNAL_NUMREGIONS(IOFILE1)=1
      SIGNAL_NUMELEC(1,IOFILE1)=NDDATA(0,0)
C *** Model specific initialisation - is this required or should
C     it all be taken care of in ipini3 or RHS routine or somewhere ????
      IF(ITYP3(NRLIST(1),nx).EQ.6) THEN !LR2
        CALL LR_INIT_CELL(ICQS,RCQS,YQS,ERROR,*9999)
      ELSEIF(ITYP3(NRLIST(1),nx).EQ.9) THEN !HH
        !CALL HH_INIT_CELL() - needs fixing
      ELSE
        CALL ASSERT(.FALSE.,'Model not implemented',ERROR,*9999)
      ENDIF
C *** open the signal file
      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,T0,WD,XID,ZD,'WRITE',FILEFORMAT,
     '  FILE02,'OPEN',ENDFILE,.TRUE.,ERROR,*9999)
C *** Write the electrode data
      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA1,0,SIGNALMAX,
     '  SIGNALMIN,T0,WD,XID,ZD,'WRITE',FILEFORMAT,
     '  FILE02,'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)
cC *** Write initial values to file
c      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
c     '  SIGNALMIN,T0,WD,XID,ZD,'WRITE',FILEFORMAT,
c     '  FILE02,'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999) 
C *** Write initial values to file
      DO nd=1,CELL_NUM_STATE
        ZD(3,nd)=YQS(CELL_STATE_OFFSET-1+nd,nq)
      ENDDO
C *** Add derived values if required
      IF(OUTPUT_DERIVED) THEN
        DO nd=CELL_NUM_STATE+1,NDDATA(0,0)
          ZD(3,nd)=YQS(CELL_DERIVED_OFFSET-1+nd-CELL_NUM_STATE,nq)
        ENDDO
      ENDIF
      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,T0,WD,XID,ZD,'WRITE',FILEFORMAT,
     '  FILE02,'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999) 
      IFAIL=1
      STIFF_EQNS=.FALSE.
C *** Solve the model
      DO T=T0,T1,DT
        IF(IFAIL.EQ.1) THEN
          STARTT=T
          ENDT=T+DT
        ELSE
          STARTT=ADAMS_WORK(4)
          ENDT=T+DT
        ENDIF
        CALL ADAMS(ICQS(CELL_AII_OFFSET),ICQS(CELL_AIO_OFFSET),
     '    ICQS(CELL_CONTROL_OFFSET),ADAMS_ERROR_TYPE,IFAIL,
     '    ADAMS_IWORK,ADAMS_LIWORK,ADAMS_LWORK,ADAMS_MAX_ITERS,
     '    ADAMS_MAX_ORDER,ICQS(CELL_MODEL_OFFSET),CELL_NUM_STATE,
     '    SIZES,ICQS(CELL_VARIANT_OFFSET),ADAMS_ABS_ERR,
     '    RCQS(CELL_ARI_OFFSET),RCQS(CELL_ARO_OFFSET),
     '    YQS(CELL_DERIVED_OFFSET,1),DY,ADAMS_MAX_STEP,
     '    RCQS(CELL_PARAMETERS_OFFSET),RCQS(CELL_PROTOCOL_OFFSET),
     '    ADAMS_REL_ERR,STARTT,ENDT,ADAMS_WORK,YQS(CELL_STATE_OFFSET,1),
     '    EXTEND_INTERVAL,STIFF_EQNS,ADAMS_USE_ROUND_CTRL,RHS_ROUTINE,
     '    ERROR)
        IF(IFAIL.NE.2) GOTO 9999
        IF(STIFF_EQNS) THEN
          WRITE(OP_STRING,'('' >>WARNING: Equations appear to be '
     '      //'stiff'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
C ***   Write values to signal file
        DO nd=1,CELL_NUM_STATE
          ZD(3,nd)=YQS(CELL_STATE_OFFSET-1+nd,1)
        ENDDO
C ***   Add derived values if required
        IF(OUTPUT_DERIVED) THEN
          DO nd=CELL_NUM_STATE+1,NDDATA(0,0)
            ZD(3,nd)=YQS(CELL_DERIVED_OFFSET-1+nd-CELL_NUM_STATE,1)
          ENDDO
        ENDIF
        TIME = T+DT !required due to FTNCHEK error
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,FILE02,
     '    'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999) 
      ENDDO
C *** Close the signal file
      CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,T1,WD,XID,ZD,'CLOSE',FILEFORMAT,
     '  FILE02,' ',ENDFILE,.TRUE.,ERROR,*9999)

      T_RESTART=T
      
      CALL EXITS('MARCH19')
      RETURN
 9999 CALL ERRORS('MARCH19',ERROR)
      CALL EXITS('MARCH19')
      RETURN 1
      END


      SUBROUTINE MODAL(IP,IBT,IDO,INP,IWK1,LGE,NBH,NBJ,
     '  NEELEM,NHE,NHP,NJE,NKE,NKH,NONY,NPF,NPNE,NPNODE,NQE,nr,
     '  NVHE,NVHP,NVJE,NVJP,NW,nx,NYNE,NYNP,NYNR,
     '  CE,CG,CONY,CP,ED,EM,ER,ES,GD,GK,GM,GR,PG,
     '  RG,SE,VE,WG,XA,XE,XG,XP,YG,YP,ZA,ZE,ZG,ZP,
     '  EIGV,GKK,GMM,GRR,WK1,WK2,WK3,WK4,XO,DYNAM,FIX,FNY,SYMM,ERROR,*)

C#### Subroutine: MODAL
C###  Description:
C###    MODAL finds eigenvalues and vectors of standard eigenproblem 
C###    if DYNAM is .false. or generalized eigenproblem if DYNAM is 
C###    .true. for square matrices NOT(nr,nx)*NOT(nr,nx).

C**** IP controls LINEAR/FIRST/UPDATE parameters as follows:
C**** IP=1 : linear=T; first=T; update=T  (linear eqtns with update)
C****  " 2     "    T    "   T         F  (  "      "   w.out  "   )
C****  " 3     "    F    "   T         T  (nonlin eqtns)
C**** SYMM indicates matrix symmetry. 
C**** IWRIT4(nr,nx) controls diagnostic output.
C**** Element stiffness matrix ES & load vector ER for linear/nonlinear
C**** problems are generated by XPES/ZEES, respec.,  and assembled into
C**** the global stiffness matrix GKK & load vector GRR.
C**** Boundary element stiffness matrices are assembled directly into a
C**** matrix GK.
C**** The element mass, damping and stiffness matrices EM,ED & ES are
C**** assembled into global matrices GM,GD & GK.
C**** The mapping coeffs NONY(noy,ny,nrc),CONY(noy,ny,nrc),
C**** noy=1,NONY(0,ny,nrc)
C**** are calculated to reduce the system of equations ny=1,NYNR to the
C**** system no=1,NOT(nr,nx) by removing constraints, if UPDATE=.true.
C**** NONY(0,ny,nrc) is 0,1 when FIX(ny,1) is .true.,.false., resp.
C**** YP(ny,1) contains essential b.c.s, defined by FIX(ny,1), on entry
C****                   solution (increm. if non.l. or dynam), on exit.
C****                or rhs vector if IP is 5 or 6.
C**** YP(ny,2)    "     nodal flux or stress bc.s defined by FIX(ny,2).
C**** YP(ny,3)    "     real part of eigenvalue
C**** YP(ny,4)    "     imag part of eigenvalue
C**** YP(ny,5)    "     reaction forces for linear solution, on exit.
C**** LGE(nhs,nrc) is location in global system of elem. var. nhs
C**** (=1,NHST(nrc)) for rows (nrc=1) and columns (nrc=2)

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:moda00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IP,IWK1(5*NOM),LGE(NHM*NSM,NRCM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NHE(NEM),
     '  NHP(NPM),NJE(NEM),NKE(NKM,NNM,NBFM,NEFM),NKH(NHM,NPM,NCM),
     '  NONY(0:NOYM,NYM,NRCM),NPF(15,NFM),
     '  NPNE(NNM,NBFM,NEFM),NPNODE(0:NP_R_M,0:NRM),
     '  NQE(NSM,NBFM,NEFM),nr,  
     '  NVHE(NNM,NBFM,NHM,NEFM),NVHP(NHM,NPM,NCM),
     '  NVJE(NNM,NBFM,NJM,NEFM),NVJP(NJM,NPM),
     '  NW(NEM,2),nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNR(0:NY_R_M,0:NRCM,NCM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CONY(0:NOYM,NYM,NRCM),CP(NMM,NPM),
     '  ED(NHM*NSM,NHM*NSM),EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),
     '  ES(NHM*NSM,NHM*NSM),GD(NYM,*),GK(NYM,*),GM(NYM,*),GR(NYROWM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEFM),
     '  VE(NSM,NKM,NEFM),WG(NGM,NBM),XA(NAM,NJM,NQM),XE(NSM,NJM),
     '  XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),YG(NGM,NJM,NEM),YP(NYM,NIYM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),
     '  ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      REAL*8 EIGV(NOM,NTM),GKK(NOM,*),GMM(NOM,*),GRR(NOM),
     '  WK1(4*NOM),WK2(3*NOM),WK3(NOM),WK4(NOM,7),XO(NOM)
      CHARACTER ERROR*(*)
      LOGICAL DYNAM,FIX(NYM,NIYFIXM),FNY(NZ_FNY_M),SYMM
!     Local Variables
      INTEGER IFAIL,MO,nc,ne,no,no1,no2,noelem,
     '  no_nynr,no_nynr1,no_nynr2,noy,noy1,noy2,
     '  nhs,nhs1,nhs2,NHST(2),nrc,ny,ny1,ny2
      REAL*8 AA,ALB,BB,CO1,CO2,EPS,EPS1,TOL,UB,X02AJF
      CHARACTER CFROMI*4,FORMAT*132
      LOGICAL FIRST,LINEAR,LOGIC(6,3),UPDATE_MATRIX,UPDATE_VECTOR
      DATA LOGIC/2*.TRUE.,4*.FALSE.,3*.TRUE.,3*.FALSE.,.TRUE.,.FALSE.,
     '           2*.TRUE.,2*.FALSE./

      CALL ENTERS('MODAL',*9999)

      ERROR='>>Recode MODAL'
      IF(ERROR.NE.' ') GOTO 9999

      nc=1 !temporary cpb 1/12/94

      LINEAR=LOGIC(IP,1)
      FIRST =LOGIC(IP,2)
      UPDATE_MATRIX=LOGIC(IP,3)
      IF(IWRIT4(nr,nx).GT.0) THEN
        FORMAT='(/'' >MODAL   linear='',L1,'' first='',L1,'' symm='''//
     '         ',L1,'' dynam='',L1,'' update='',L1)'
        WRITE(OP_STRING,FORMAT) LINEAR,FIRST,SYMM,DYNAM,UPDATE_MATRIX
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      no=0
      DO nrc=1,2
        DO no_nynr=1,NYNR(0,nrc,1)
          ny=NYNR(no_nynr,nrc,1)
          IF(FIX(ny,1)) THEN
            NONY(0,ny,nrc)=0
          ELSE
            NONY(0,ny,nrc)=1
          ENDIF
          DO noy=1,NONY(0,ny,nrc)
            no=no+1
            NONY(noy,ny,nrc)=no
            CONY(noy,ny,nrc)=1.0d0
          ENDDO
        ENDDO
      ENDDO
      NOT(1,1,nr,nx)=no
      DO no1=1,NOT(1,1,nr,nx)
        GRR(no1)=0.0d0
        DO no2=1,NOT(2,1,nr,nx)
          GKK(no1,no2)=0.0d0
          GMM(no1,no2)=0.0d0
        ENDDO
      ENDDO

      DO no_nynr=1,NYNR(0,1,1)
        ny=NYNR(no_nynr,1,1)
        IF(FIX(ny,2)) THEN
          GR(ny)=YP(ny,2)
        ELSE
          GR(ny)=0.0d0
        ENDIF
      ENDDO

      IF(UPDATE_MATRIX) THEN
        DO no_nynr1=1,NYNR(0,1,1)
          ny1=NYNR(no_nynr1,1,1)
          DO no_nynr2=1,NYNR(0,2,1)
            ny2=NYNR(no_nynr2,2,1)
            GK(ny1,ny2)=0.0d0
            IF(DYNAM) THEN
              GM(ny1,ny2)=0.0d0
              GD(ny1,ny2)=0.0d0
            ENDIF
          ENDDO
        ENDDO

        DO noelem=1,NEELEM(0,nr)    !main element loop
          ne=NEELEM(noelem,nr)
          IF(NW(ne,1).GE.0) THEN
            CALL MELGE(LGE,NBH(1,1,ne),nc,ne,NHE(ne),NHST,NKH,
     '        NPNE(1,1,ne),nr,NVHE(1,1,1,ne),nx,NYNE,NYNP,ERROR,*9999)

            DO nhs1=1,NHST(1)    !Initialize element arrays
              ER(nhs1)=0.0d0
              DO nhs2=1,NHST(2)
                ES(nhs1,nhs2)=0.0d0
              ENDDO
            ENDDO

            IF(LINEAR) THEN
              IF(NW(ne,1).LE.20) THEN
                IF(ITYP1(nr,nx).EQ.3) THEN
                  CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),
     '              NPNE(1,1,ne),NQE(1,1,ne),nr,NVJE(1,1,1,ne),
     '              SE(1,1,ne),XA,XE,XP,ERROR,*9999)
                  CALL XPES30(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),ne,
     '              NHE(ne),NJE(ne),NPNE(1,1,ne),nr,nx,
     '              CE(1,ne),CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '              VE(1,1,ne),WG,XE,XG,YG(1,1,ne),ZE,ZG,
     '              UPDATE_MATRIX,UPDATE_VECTOR,ERROR,*9999)
                ELSE IF(ITYP1(nr,nx).EQ.4) THEN
                  CALL XPES40(NBH(1,1,ne),NBJ(1,ne),
     '              NHE(ne),NJE(ne),NKE(1,1,1,ne),
     '              NPF,NPNE(1,1,ne),NQE(1,1,ne),nr,
     '              NVJE(1,1,1,ne),NW(ne,1),nx,
     '              CE(1,ne),CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '              VE(1,1,ne),WG,XA,XE,XG,XP,YG(1,1,ne),
     '              UPDATE_MATRIX,ERROR,*9999)
                ENDIF
              ELSE
                IF(ITYP1(nr,nx).EQ.3) THEN
                ELSE IF(ITYP1(nr,nx).EQ.4) THEN
                ENDIF
              ENDIF
            ELSE IF(.NOT.LINEAR) THEN
            ENDIF

            IF(iwrit4(nr,nx).GE.5.and.NW(ne,1).LE.20) THEN
              WRITE(OP_STRING,
     '          '(/'' Element load vector ER & stiffness matrix ES:'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              DO nhs1=1,NHST(1)
                WRITE(OP_STRING,'('' ER('',I2,'')='',D12.4,'' ES: '','
     '            //'8D12.4/,(25X,8D12.4))') 
     '            nhs1,ER(nhs1),(ES(nhs1,nhs2),nhs2=1,NHST(2))
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDDO
              IF(DYNAM) THEN
                WRITE(OP_STRING,'(/'' Element matrices EM & ED:'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                DO nhs1=1,NHST(1)
                  WRITE(OP_STRING,'('' nhs1='',I2,'' EM: '',8D12.4/,'
     '              //'(12X,8D12.4))')nhs1,(EM(nhs1,nhs2),nhs2=1,
     '              NHST(2))
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDDO
                DO nhs1=1,NHST(1)
                  WRITE(OP_STRING,'('' nhs1='',I2,'' ED: '',8D12.4/,'
     '              //'(12X,8D12.4))')nhs1,(ED(nhs1,nhs2),nhs2=1,
     '              NHST(2))
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDDO
              ENDIF
            ENDIF
                  
            DO nhs1=1,NHST(1)
              ny1=LGE(nhs1,1)
              GR(ny1)=GR(ny1)+ER(nhs1)
              DO nhs2=1,NHST(2)
                ny2=LGE(nhs2,2)
                GK(ny1,ny2)=GK(ny1,ny2)+ES(nhs1,nhs2)
                IF(DYNAM) THEN
                  GM(ny1,ny2)=GM(ny1,ny2)+EM(nhs1,nhs2)
                  GD(ny1,ny2)=GD(ny1,ny2)+ED(nhs1,nhs2)
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO

        IF(IWRIT4(nr,nx).GE.4) THEN
          WRITE(OP_STRING,
     '      '(/'' Global load vector GR & stiffness matrix GK:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DO no_nynr1=1,NYNR(0,1,1)
            ny1=NYNR(no_nynr1,1,1)
            WRITE(OP_STRING,'('' GR('',I4,'')='',D12.4,'' GK: '','
     '	      //'D12.4'//'/,(25X,8D12.4))') ny1,GR(ny1),
     '	      (GK(ny1,NYNR(no_nynr2,2,1)),no_nynr2=1,NYNR(0,2,1))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDDO
          IF(DYNAM) THEN
            WRITE(OP_STRING,'(/'' Global matrices GM & GD:'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            DO no_nynr1=1,NYNR(0,1,1)
              ny1=NYNR(no_nynr1,1,1)
              WRITE(OP_STRING,'('' NY1='',I2,'' GM: '',8D12.4'//
     '          '/,(12X,8D12.4))') ny1,
     '	      (GM(ny1,NYNR(no_nynr2,2,1)),no_nynr2=1,NYNR(0,2,1))
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDDO
            DO no_nynr1=1,NYNR(0,1,1)
              ny1=NYNR(no_nynr1,1,1)
              WRITE(OP_STRING,'('' NY1='',I2,'' GD: '',8D12.4'//
     '          '/,(12X,8D12.4))') ny1,
     '	      (GD(ny1,NYNR(no_nynr2,2,1)),no_nynr2=1,NYNR(0,2,1))
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF
        ENDIF

      ELSE IF(.NOT.UPDATE_MATRIX) THEN
        DO noelem=1,NEELEM(0,nr)    !main element loop
          ne=NEELEM(noelem,nr)
          IF(NW(ne,1).GE.0) THEN
            CALL MELGE(LGE,NBH(1,1,ne),1,ne,NHE(ne),NHST,NKH,
     '        NPNE(1,1,ne),nr,NVHE(1,1,1,ne),nx,NYNE,NYNP,ERROR,*9999)

            CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),NQE(1,1,ne),nr,NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA,XE,XP,ERROR,*9999)
            CALL YPZP(4,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '        nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
            CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '        SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
            CALL XPES30(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),ne,
     '        NHE(ne),NJE(ne),NPNE(1,1,ne),nr,nx,
     '        CE(1,ne),CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '        VE(1,1,ne),WG,XE,XG,YG(1,1,ne),ZE,ZG,
     '        UPDATE_MATRIX,UPDATE_VECTOR,ERROR,*9999)

            DO nhs=1,NHST(1)
              ny=IABS(LGE(nhs,1))
              GR(ny)=GR(ny)+ER(nhs)
            ENDDO
          ENDIF
        ENDDO
      ENDIF

      nrc=2 !temporary ajp 30/11/94
      DO no_nynr1=1,NYNR(0,1,1)
        ny1=NYNR(no_nynr1,1,1)
        DO noy1=1,NONY(0,ny1,1)
          no1=NONY(noy1,ny1,1)
          CO1=CONY(noy1,ny1,1)
          DO no_nynr2=1,NYNR(0,2,1)
            ny2=NYNR(no_nynr2,2,1)
            AA=GK(ny1,ny2)
            IF(DYNAM) THEN
              BB=GM(ny1,ny2)
            ELSE
              IF(ny1.EQ.ny2) THEN
                BB=1.0D0
              ELSE
                BB=0.0D0
              ENDIF
            ENDIF
            IF(NONY(0,ny2,2).ne.0) THEN
              DO noy2=1,NONY(0,ny2,2)
                no2=NONY(noy2,ny2,2)
                CO2=CONY(noy2,ny2,2)
                GKK(no1,no2)=GKK(no1,no2)+AA*CO1*CO2
                GMM(no1,no2)=GMM(no1,no2)+BB*CO1*CO2
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      IF(IWRIT4(nr,nx).GE.3) THEN
        WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I4)') NOT(1,1,nr,nx)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/'' Global stiffness matrix GKK:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        DO no1=1,NOT(1,1,nr,nx)
          WRITE(OP_STRING,'(/'' GKK('',I4,'',.)='',8D12.4,'
     '	    //'(/25X,8D12.4))') no1,(GKK(no1,no2),no2=1,NOT(2,1,nr,nx))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'(/'' Global mass matrix GMM:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        DO no1=1,NOT(1,1,nr,nx)
          WRITE(OP_STRING,'(/'' GMM('',I4,'',.)='',8D12.4,'
     '	    //'(/25X,8D12.4))') no1,(GMM(no1,no2),no2=1,NOT(2,1,nr,nx))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

C Old code
C      IFAIL=1
C      CALL F02BJF(NOT,GKK,NOM,GMM,NOM,-1.D0,WK1,WK2,XO,.TRUE.,
C     '  EIGV,NOM,IWK1,IFAIL)
C      IF(ifail.ne.0) THEN
C        ERROR=' IFAIL='//CFROMI(IFAIL,'(I4)')//' in F02BJF'
C        GOTO 9999
C      ENDIF
C
C      WRITE(OP_STRING,'(/'' Eigenvalues and eigenvectors:'')')
C      no=0
C      DO no_nynr=1,NYNR(0,1,1)
C        ny=NYNR(no_nynr,1,1)
C        IF(.NOT.FIX(ny,1)) THEN
C          no=no+1
C          IF(DABS(XO(no)).GT.1.0D-6) THEN
C            ALFA=XO(no)
C            YP(ny,3)=WK1(no)/ALFA !is real component of eigenvalue
C            YP(ny,4)=WK2(no)/ALFA !is imag component of eigenvalue
C            EVALR=WK1(no)/ALFA
C            EVALI=WK2(no)/ALFA
C            WRITE(OP_STRING,'(1X,I4,'') Real: '',D12.4,'' Imaginary: '',
C     '        D12.4)') no,EVALR,EVALI
C            WRITE(OP_STRING,'(1X,10D12.4)') (EIGV(MO,no),MO=1,NOT(nr,nx))
C          ELSE
C            WRITE(OP_STRING,'(1X,I4,'') Eigenvalue is singular'')') no
C          ENDIF
C        ENDIF
C      ENDDO
C
C New code to calculate only a small number of the eigenvalues
C AJP 5-7-91
      IFAIL=1
      !Reduce from general Kx=lamda.Mx to standard eigenproblem Pz=lamda.z
      CALL F01AEF(NOT(1,1,nr,nx),GKK,NOM,GMM,NOM,WK1,IFAIL)
      IF(ifail.ne.0) THEN
        ERROR=' IFAIL='//CFROMI(IFAIL,'(I4)')//' in F01AEF'
        GOTO 9999
      ENDIF
      !Householder reduction of matrix to tridiagonal form
      CALL F01AGF(NOT(1,1,nr,nx),TOL,GKK,NOM,GRR,WK2,WK3)
      IF(KTYP16.EQ.1) THEN !Lowest evalues wanted
        ALB=0.0D0
        UB=10.0D0
      ELSE IF(KTYP16.EQ.2) THEN !Highest evalues wanted
        ALB=0.0D0   !May need changing
        UB=100.0D0  !May need changing
      ENDIF
      EPS=X02AJF() !Machine precision
      EPS1=0.0D0
      IF(KTYP17.GT.NTM) THEN
        ERROR=' no. of evalues req. exceeds NTM ='//CFROMI(NTM,'(I4)')
        GOTO 9999
      ENDIF
1000  IFAIL=1
      !Solve tridiagonal problem - gives final eigenvalues in XO
      CALL F02BEF(NOT(1,1,nr,nx),GRR,ALB,UB,EPS,EPS1,WK2,WK3,NTM,NTMAX,
     '  XO,EIGV,NOM,IWK1,WK4,FNY,IFAIL)
      IF(IFAIL.EQ.1) THEN !More than NTM evalues in given interval.
        UB=UB/2 !Since only a few evalues reqd, reduce interval size and repeat.
        WRITE(OP_STRING,*)
     '    ' WARNING - decreased interval size for evalue search'
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GOTO 1000
      ENDIF
      IF(ifail.ne.0) THEN
        ERROR=' IFAIL='//CFROMI(IFAIL,'(I4)')//' in F02BEF'
        GOTO 9999
      ENDIF
      IF(NTMAX.LT.KTYP17) THEN !Not enough evalues in interval
        UB=UB+UB/2.0d0 !Inc interval size to ensure sufficient evalues
        WRITE(OP_STRING,*)
     '    ' WARNING - increased interval size for evalue search'
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GOTO 1000
      ENDIF
      !Return eigenvectors of matrix P
      CALL F01AHF(NOT(1,1,nr,nx),1,NTMAX,GKK,NOM,WK2,EIGV,NOM)
      !Return eigenvectors of Kx=lamda.Mx
      CALL F01AFF(NOT(1,1,nr,nx),1,NTMAX,GMM,NOM,WK1,EIGV,NOM)

      WRITE(OP_STRING,'(/'' Eigenvalues and eigenvectors:'')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      no=0
      DO no_nynr=1,NYNR(0,1,1)
        ny=NYNR(no_nynr,1,1)
        IF(.NOT.FIX(ny,1)) THEN
          no=no+1
          IF(no.LE.NTMAX) THEN
            YP(ny,3)=XO(no)
            YP(ny,4)=0.0D0    !Evalues are real
            WRITE(OP_STRING,'(1X,I4,1X,D12.4)') no,XO(no)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(1X,10D12.4)')(EIGV(MO,no),MO=1,
     '        NOT(1,1,nr,nx))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDDO

      CALL EXITS('MODAL')
      RETURN
 9999 CALL ERRORS('MODAL',ERROR)
      CALL EXITS('MODAL')
      RETURN 1
      END

      SUBROUTINE MONIT(N,XC,FC,GC,ISTATE,GPJNRM,COND,POSDEF,NITER,nf,
     '  IW,LIW,W,LW)

C#### Subroutine: MONIT
C###  Description:
C###    MONIT monitors minimisation process for NAG routine E04JBF.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER ISTATE(*),IW(*),LIW,LW,N,nf,NITER
      REAL*8 COND,FC,GC(*),GPJNRM,W(*),XC(*)
      LOGICAL POSDEF
!     Local Variables
      CHARACTER ERROR*10

      WRITE(OP_STRING,'('' Number of iterations performed  = '',I4)')
     '	NITER
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Number of function evaluations  = '',I4)')
     '	nf
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,
     '  '('' Norm of projected gradient vector = '',D12.3)') GPJNRM
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,
     '  '('' Condition # of projected Hessian  = '',D12.3)') COND
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

 9999 RETURN
      END

      SUBROUTINE SOLVE1(IP,IBT,IDO,INP,LGE,NBH,NBJ,
     '  NEELEM,NHE,NHP,NJE,NKE,NKH,NONY,NPE,NPF,NPNODE,NQE,nr,
     '  NRE,NVNE,NVNP,NW,nx,NYNE,NYNP,A,B,CE,CG,CONY,CP,ED,EM,ER,ES,GD,
     '  GK,GM,GR,PAOPTI,PBOPTI,PG,RG,SE,VE,WG,XA,XE,XG,
     '  XP,YG,YP,ZA,ZE,ZG,ZP,GKK,GRR,WK1,WK2,XO,DYNAM,FIX,
     '  SYMM,ERROR,*)

C**** Finds solution of a system of linear equations no=1,NOT(1,nr,nx).
C**** IP controls LINEAR/FIRST/UPDATE parameters as follows:
C**** IP=1 : linear=T; first=T; update=T  (linear eqtns with update)
C****  " 2     "    T    "   T         F  (  "      "   w.out  "   )
C****  " 3     "    F    "   T         T  (nonlin eqtns)
C**** SYMM indicates matrix symmetry. IWRIT4 controls diagnostic output.
C**** DYNAM is .true. when SOLVE1 called in time integration procedures.
C**** Element stiffness matrix ES & load vector ER for linear/nonlinear
C**** problems are generated by XPES/ZEES, respec.,  and assembled into
C**** the global stiffness matrix GKK & load vector GRR.
C**** In dynamic problems (DYNAM =.true.) the element mass, damping and
C**** stiffness matrices EM,ED & ES are assembled into global matrices
C**** GM,GD & GK which are then mapped into the solution matrix GKK.
C**** Note: These contain 2nd,1st & 0th order time-deriv coeffs,respec.
C**** The mapping coeffs NONY(noy,ny,nrc),CONY(noy,ny,nrc),
C**** noy=1,NONY(0,ny,nrc)
C**** are calc.d to reduce the system of equations ny=1,NYT(nrc,1,nx) to
C**** system no=1,NOT(1,nr,nx) by removing constraints, if UPDATE=.true.
C**** NONY(0,ny,nrc) is 0,1 when FIX(ny,1) is .true.,.false., resp.
C**** Note: In dynamic problems GKK is reformed at each time step since
C**** time step may change, so linear eqtn solution is done in one step
C**** either by Crout for non-symmetric problems or Cholesky for symm.
C**** KTYP22=0 : Static solution of initial accelerations.
C****   "    1 : first  order (Crank-Nicholson) time integration.
C****   "    2 : second order (Newmark; Wood-Zienkiewicz etc)
C****   "    3 : third  order (Houbolt; Hilbert-Hughes etc)
C**** YP(ny,1) contains essential b.c.s, defined by FIX(ny,1), on entry
C****                   solution (increm. if non.l. or dynam), on exit.
C****                or rhs vector if IP is 5 or 6.
C**** YP(ny,2)    "     nodal flux or stress bc.s defined by FIX(ny,2).
C**** YP(ny,3)    "     incremental b.c.s             "      FIX(ny,3).
C**** YP(ny,5)    "     reaction forces for linear solution, on exit.
C**** LGE(nhs,nrc) is location in global system of elem. var. nhs
C****   (=1,NHST(nrc))
C**** 23-9-92 AJP
C**** If FIX(ny,3) is true for a linear problem then we have a 
C**** mixed bc of the form aT+b(dT/dn)=c,
C**** where a,b and c are stored in YP(ny,1), YP(ny,2) and YP(ny,3) 
C**** respectively. dt/dn is treated as an integrated flux value.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b08.cmn'
      INCLUDE 'cmiss$reference:b10.cmn'
      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:coup00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'    
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp90.cmn'
      INCLUDE 'cmiss$reference:ktyp100.cmn'
      INCLUDE 'cmiss$reference:opti00.cmn'
!     Parameter List
      INTEGER IBT(2,NIM,*),IDO(NKM,0:NIM,*),INP(NNM,NIM,*),IP,
     '  LGE(NHM*NSM,*),NBH(NHM,NCM,*),NBJ(NJM,*),
     '  NEELEM(0:NEM,0:*),
     '  NHE(*),NHP(*),NJE(*),NKE(NKM,NNM,NBFM,*),NKH(NHM,NCM,*),
     '  NONY(0:NOYM,NYM,*),NPE(NNM,NBFM,*),
     '  NPF(12,*),
     '  NPNODE(0:NPM,0:*),NQE(NSM,NBFM,*),
     '  nr,NRE(*),NVNE(NNM,NBFM,NJM,*),NVNP(NJM,NPM,*),
     '  NW(NEM,*),nx,NYNE(NAM,NHM,NRCM,NCM,*),
     '  NYNP(NKM,NVM,NHM,NPM,NRCM,*)
      REAL*8 A(NSM,*),B(NSM,*),CE(NMM,*),CG(NMM,*),CONY(0:NOYM,NYM,*),
     '  CP(NMM,*),ED(NVM,*),EM(NVM,*),ER(*),ES(NVM,*),
     '  GD(NYM,*),GK(NYM,*),GM(NYM,*),GR(*),
     '  PAOPTI(*),PBOPTI(*),
     '  PG(NSM,NUM,NGM,*),RG(*),
     '  SE(NSM,NBFM,*),VE(NSM,NKM,*),WG(NGM,*),
     '  XA(NAM,NJM,*),XE(NSM,*),XG(NJM,*),XP(NKM,NVM,NJM,*),
     '  YG(NGM,NJM,*),YP(NYM,*),
     '  ZA(NAM,NHM,NCM,*),ZE(NSM,*),ZG(NHM,*),
     '  ZP(NKM,NVM,NHM,NPM,*)
      REAL*8 GKK(NOM,*),GRR(*),WK1(*),WK2(*),XO(*)
      CHARACTER ERROR*(*)
      LOGICAL DYNAM,FIX(NYM,*),SYMM
!     Local Variables
      INTEGER IFAIL,nc,ne,nj,no,no1,no2,noelem,nonode2,noopti,
     '  noy,noy1,noy2,NP1,NP2,NP3,nrc,nhs,nhs1,nhs2,NHST(2),nv,ny,
     '  ny1,ny2,NY3,NY4  
      REAL*8 AA,AB,BB,CO,CO1,CO2,SF,SUM,XMOVE,XFIX
      CHARACTER FORMAT*500
      LOGICAL FIRST,LINEAR,LOGIC(6,3),MIXED,UPDATE
      DATA LOGIC/2*.TRUE.,4*.FALSE.,3*.TRUE.,3*.FALSE.,.TRUE.,.FALSE.,
     '           2*.TRUE.,2*.FALSE./

      CALL ENTERS('SOLVE1',*9999)
                 
      nc=1 ! temporary cpb 22/11/94
      nrc=2 ! temporary ajp 30/11/94
      nv=1 ! temporary cpb 22/11/94

      LINEAR=LOGIC(IP,1)
      FIRST =LOGIC(IP,2)
      UPDATE=LOGIC(IP,3)
      IF(IWRIT4.GT.0) THEN
        FORMAT='(/'' >SOLVE1   linear='',L1,'' first='',L1,'' symm='''//
     '    ',L1,'' update='',L1)'
        WRITE(OP_STRING,FORMAT) LINEAR,FIRST,SYMM,UPDATE
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL ASSERT(NYT(2,1,nx).GT.0,
     '  '>>no initial conditions for this region',ERROR,*9999)
      MIXED=.FALSE.

      IF(NOT(1,nr,nx).EQ.0) THEN
        ERROR=' >>The number of unknowns is zero'
        GO TO 9999
      ENDIF
      DO no1=1,NOT(1,nr,nx)
        GRR(no1)=0.0d0
        DO no2=1,NOT(2,nr,nx)
          GKK(no1,no2)=0.0d0
        ENDDO
      ENDDO

! Note: this should be moved to upopti - but check with AJP first.

      IF(KTYP27.EQ.3) THEN !Coupled saturated-unsaturated flow.
C ***   Update nodal coordinates from optimisation params
C ***   As well as updating interface nodal coordinates we need
C ***   to ensure that elements in the FE region are not becoming too 
C ***   distorted by updating the coordinates of all other nodes in 
C ***   each FE element (except for the nodes around the fixed bdy).
C ***   NPRAD(I,noopti,nr), I=1,..,NPRAD(0,noopti,nr) are the nodes
C ***   on the same radial line as the current (noopti) free surface 
C ***   node which are in the FE region and not on any fixed surface.
C ***   The node on the fixed surface is identified 
C ***   by NPRAD(-1,noopti,nr).  
C ***   This node number is needed in the mesh updating.
C ***   The array NPRAD is set up in IPOPTI
        DO noopti=1,NTOPTI
          NP1=NPJOIN(noopti,1) !Appropriate free surface node number
          NP3=NPRAD(-1,noopti,nr) !Appropriate fixed surface node number
          SF=PAOPTI(noopti)/PBOPTI(noopti)
          DO nj=1,NJT
            XFIX=XP(1,nv,nj,NP3)
            XMOVE=XP(1,nv,nj,NP1)
            IF(DABS(XFIX-XMOVE).GT.0.001D0) THEN
              !If XFIX-XMOVE=0 then no update necessary.
              DO nonode2=1,NPRAD(0,noopti,nr)
                NP2=NPRAD(nonode2,noopti,nr) !Node numbers on same radial line
                !Updating nodal coordinates on same radial line as current
                !free surface node.  Need to ensure that this node remains
                !in the same position relative to the free surface node
                !and the nodes on the fixed surface.  Therefore need to
                !use the coordinates of both the free surface node and the
                !fixed surface node.
                XP(1,nv,nj,NP2)=((SF-1.0d0)*XFIX*XMOVE+
     '            (XFIX-SF*XMOVE)*XP(1,nv,nj,NP2))/(XFIX-XMOVE)
              ENDDO !End of loop over radial line nodes.
            ENDIF
            !Updating free surface nodal coordinates
            XP(1,nv,nj,NP1)=SF*XMOVE
          ENDDO
          PBOPTI(noopti)=PAOPTI(noopti)
          IF(KTYP100.EQ.2) THEN !Flow around cavities
            IF(ITYP10(1).EQ.1) THEN !Cartesians
              ny=NYNP(1,1,1,np3,nrc,nc)
              YP(ny,1)=XP(1,nv,NJT,NP1)-XP(1,nv,NJT,NP3)
            ENDIF
          ENDIF
        ENDDO

      ENDIF

! Update RHS vector GR from flux boundary conditions
      DO ny=1,NYT(nrc,1,nx)
        IF(FIX(ny,2)) THEN
          GR(ny)=YP(ny,2)
        ELSEIF(FIX(ny,3).AND.(ITYP6(nr).EQ.1))THEN !Mixed bc and linear
          GR(ny)=YP(ny,3)/YP(ny,2) !c/b
          MIXED=.TRUE.
        ELSE
          GR(ny)=0.0D0
        ENDIF
      ENDDO

      IF(UPDATE) THEN
        DO ny1=1,NYT(1,1,nx)
          DO ny2=1,NYT(2,1,nx)
            GK(ny1,ny2)=0.0D0
            IF(DYNAM) THEN
              GM(ny1,ny2)=0.0D0
              GD(ny1,ny2)=0.0D0
            ENDIF
          ENDDO
        ENDDO

        DO noelem=1,NEELEM(0,nr) !is main element loop
          ne=NEELEM(noelem,nr)
          IF(NW(ne,1).GT.0) THEN

            CALL MELGE(LGE,NBH(1,1,ne),nc,NHE(ne),NHST,NKH,NPE(1,1,ne),
     '        NVNE(1,1,1,ne),NYNE(1,1,1,nc,ne),NYNP(1,1,1,1,1,nc),
     '        ERROR,*9999)

            IF(IWRIT4.GT.1) THEN
              FORMAT='(/'' Element'',I3,'' Number of variables '
     '          //'NHST='',2(I3)/,'' LGE(1..,1): '',18I4/,(6X,18I4),'
     '          //''' LGE(1..,2):  '',18I4/,(6X,18,I4))'
              WRITE(OP_STRING,FORMAT) ne,NHST(1),NHST(2),(LGE(nhs,1),
     '          nhs=1,NHST(1)),(LGE(nhs,2),nhs=1,NHST(2))
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF            

            DO nhs1=1,NHST(1)
              ER(nhs1)=0.0d0
              DO nhs2=1,NHST(2)
                ES(nhs1,nhs2)=0.0d0
              ENDDO
            ENDDO

            IF(NW(ne,1).LE.20) THEN !finite element
              IF(ITYP1(nr).EQ.3) THEN !partial differential equation
                CALL XPXE(NBJ(1,ne),NJE(ne),NKE(1,1,1,ne),NPE(1,1,ne),
     '            NPF(1,1),NQE(1,1,ne),nr,NVNE(1,1,1,ne),NVNP,
     '            SE(1,1,ne),XA,XE,XP,ERROR,*9999)
                CALL XPES30(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),ne,
     '            NHE(ne),NJE(ne),NPE(1,1,ne),nr,
     '            CE(1,ne),CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '            VE(1,1,ne),WG,XE,XG,ZE,ZG,UPDATE,ERROR,*9999)
              ELSE IF(ITYP1(nr).EQ.4) THEN !linear elasticity
                CALL XPES40(NBH(1,1,ne),NBJ(1,ne),
     '            NHE(ne),NJE(ne),NKE(1,1,1,ne),NPE(1,1,ne),NPF,
     '            NQE(1,1,ne),nr,NVNE(1,1,1,ne),NVNP,NW(ne,1),CE(1,ne),
     '            CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '            VE(1,1,ne),WG,XA,XE,XG,XP,YG(1,1,ne),UPDATE,
     '            ERROR,*9999)
              ENDIF
            ELSE !boundary element
              ERROR=' >>Bdry elements should use SOLVE4'
              GO TO 9999
            ENDIF

            IF(IWRIT4.GE.3) THEN
              WRITE(OP_STRING,
     '          '(/'' Element load vector ER & stiffness matrix ES:'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              DO nhs1=1,NHST(1)
                WRITE(OP_STRING,'('' ER('',I2,'')='',E12.4,'' ES: '','
     '            //'8E12.4/,(25X,8E12.4))') 
     '            nhs1,ER(nhs1),(ES(nhs1,nhs2),nhs2=1,NHST(2))
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDDO
              IF(DYNAM) THEN
                WRITE(OP_STRING,'(/'' Element matrices EM & ED:'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                DO nhs1=1,NHST(1)
                  WRITE(OP_STRING,
     '              '('' nhs1='',I2,'' EM: '',8E12.4/,(12X,8E12.4))')
     '              nhs1,(EM(nhs1,nhs2),nhs2=1,NHST(2))
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDDO
                DO nhs1=1,NHST(1)
                  WRITE(OP_STRING,'('' nhs1='',I2,'' ED: '',8E12.4/,'
     '              //'(12X,8E12.4))')nhs1,(ED(nhs1,nhs2),nhs2=1,
     '              NHST(2))
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDDO
              ENDIF
            ENDIF

            DO nhs1=1,NHST(1)
              ny1=IABS(LGE(nhs1,1))
              GR(ny1)=GR(ny1)+ER(nhs1)
              DO nhs2=1,NHST(2)
                ny2=IABS(LGE(nhs2,2))
                GK(ny1,ny2)=GK(ny1,ny2)+ES(nhs1,nhs2)
                IF(DYNAM) THEN
                  GM(ny1,ny2)=GM(ny1,ny2)+EM(nhs1,nhs2)
                  GD(ny1,ny2)=GD(ny1,ny2)+ED(nhs1,nhs2)
                ENDIF
              ENDDO
            ENDDO

          ENDIF
        ENDDO

        IF(IWRIT4.GE.4) THEN
          WRITE(OP_STRING,'(/'' Global load vector GR & stiffness'
     '      //' matrix GK:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DO ny1=1,NYT(1,1,nx)
            WRITE(OP_STRING,'('' GR('',I4,'')='',E12.4,'' GK: '','
     '	      //'8E12.4'//'/,(25X,8E12.4))') 
     '        ny1,GR(ny1),(GK(ny1,ny2),ny2=1,NYT(2,1,nx))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDDO
          IF(DYNAM) THEN
            WRITE(OP_STRING,'(/'' Global matrices GM & GD:'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            DO ny1=1,NYT(1,1,nx)
              WRITE(OP_STRING,'('' NY1='',I2,'' GM: '',8E12.4'
     '          //'/,(12X,8E12.4))') 
     '          ny1,(GM(ny1,ny2),ny2=1,NYT(2,1,nx))
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDDO
            DO ny1=1,NYT(1,1,nx)
              WRITE(OP_STRING,'('' NY1='',I2,'' GD: '',8E12.4'
     '          //'/,(12X,8E12.4))') 
     '          ny1,(GD(ny1,ny2),ny2=1,NYT(2,1,nx))
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF
        ENDIF

      ELSE IF(.NOT.UPDATE) THEN
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(NW(ne,1).GT.0) THEN

            CALL MELGE(LGE,NBH(1,1,ne),nc,NHE(ne),NHST,NKH,NPE(1,1,ne),
     '        NVNE(1,1,1,ne),NYNE(1,1,1,nc,ne),NYNP(1,1,1,1,1,nc),
     '        ERROR,*9999)

            CALL XPXE(NBJ(1,ne),NJE(ne),NKE(1,1,1,ne),NPE(1,1,ne),
     '        NPF(1,1),NQE(1,1,ne),nr,NVNE(1,1,1,ne),NVNP,SE(1,1,ne),
     '        XA,XE,XP,ERROR,*9999)
            CALL YPZP(4,NBH,NEELEM,
     '        NHE,NHP,NKH,NPNODE,nr,NVNP,NYNE,NYNP,
     '        YP,ZA,ZP,ERROR,*9999)
            CALL ZPZE(1,NBH(1,1,ne),NHE(ne),NKE(1,1,1,ne),NPE(1,1,ne),
     '        NPF(1,1),nr,NVNE(1,1,1,ne),NVNP,NW(ne,1),SE(1,1,ne),
     '        ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
            CALL XPES30(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),ne,
     '        NHE(ne),NJE(ne),NPE(1,1,ne),nr,
     '        CE(1,ne),CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '        VE(1,1,ne),WG,XE,XG,ZE,ZG,UPDATE,ERROR,*9999)
            DO nhs=1,NHST(1)
              ny=IABS(LGE(nhs,1))
              GR(ny)=GR(ny)+ER(nhs)
            ENDDO

          ENDIF
        ENDDO
      ENDIF

      IF(KTYP90.EQ.6) THEN ! Linear boundary layer CPB/AJP 26/7/94
        DO ny1=1,NYT(nrc,1,nx)
          DO noy1=1,NONY(0,ny1,1)
            no1=NONY(noy1,ny1,1)
            CO1=CONY(noy1,ny1,1)
            IF(noy1.EQ.1) THEN
              NY3=ny1
            ELSE
              NY3=1
              DO WHILE(NY3.LE.NYT(nrc,1,nx).AND..NOT.(NONY(1,NY3,1).EQ.
     '          no1.and.ny3.ne.ny1))
                NY3=NY3+1
              ENDDO
              IF(NONY(1,NY3,1).ne.no1) THEN
                ERROR='>>Coupled no not found'
                GOTO 9999
              ENDIF
            ENDIF
            BB=GR(NY3)
            GRR(no1)=GRR(no1)+BB*CO1
            DO ny2=1,NYT(nrc,1,nx)
              IF(NONY(0,ny2,2).EQ.0) THEN
                AA=GK(NY3,ny2)
                GRR(no1)=GRR(no1)-AA*CO1*YP(ny2,1)
              ELSE
                DO noy2=1,NONY(0,ny2,2)
                  no2=NONY(noy2,ny2,2)
                  CO2=CONY(noy2,ny2,2)
                  IF(noy2.EQ.1) THEN
                    NY4=ny2
                  ELSE
                    NY4=1
                    DO WHILE(NY4.LE.NYT(nrc,1,nx).AND..NOT.
     '                (NONY(1,NY4,2).eq.no2.and.ny4.ne.ny2))
                      NY4=NY4+1
                    ENDDO
                    IF(NONY(1,NY4,2).ne.no2) THEN
                      ERROR='>>Coupled no not found'
                      GOTO 9999
                    ENDIF
                  ENDIF
                  AA=GK(NY3,NY4)
C cpb 27/7/94 Do we want both coupling coefficients for the noy2=1 case?
                  GKK(no1,no2)=GKK(no1,no2)+AA*CO1*CO2
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO ny1=1,NYT(nrc,1,nx)
          DO noy1=1,NONY(0,ny1,1)
            no1=NONY(noy1,ny1,1)
            CO1=CONY(noy1,ny1,1)
            BB=GR(ny1)
            DO ny=1,NYT(nrc,1,nx)
              IF(KTYP22.GE.1) BB=BB-GK(ny1,ny)*YP(ny,10)
              IF(KTYP22.GE.2) BB=BB-GD(ny1,ny)*YP(ny,11)
              IF(KTYP22.EQ.3) BB=BB-GM(ny1,ny)*YP(ny,12)
            ENDDO
            GRR(no1)=GRR(no1)+BB*CO1
            DO ny2=1,NYT(nrc,1,nx)
              IF(.NOT.DYNAM) THEN
                AA=GK(ny1,ny2)
              ELSE IF(DYNAM) THEN
                IF(KTYP22.EQ.0) THEN
                  AA=GM(ny1,ny2)
                ELSE IF(KTYP22.EQ.1) THEN
                  AA=GD(ny1,ny2)+A1*GK(ny1,ny2)
                ELSE IF(KTYP22.EQ.2) THEN
                  AA=GM(ny1,ny2)+A1*GD(ny1,ny2)+A2*GK(ny1,ny2)
                ELSE IF(KTYP22.EQ.3) THEN
                  AA=A1*GM(ny1,ny2)+A2*GD(ny1,ny2)+A3*GK(ny1,ny2)
                ENDIF
              ENDIF
              IF(FIX(ny2,3).AND.(ITYP6(nr).EQ.1))THEN !Mixed bc and linear
                AB=YP(ny2,1)/YP(ny2,2) !a/b
              ELSE
                AB=0
              ENDIF
              IF(NONY(0,ny2,2).EQ.0) THEN
                GRR(no1)=GRR(no1)-AA*CO1*YP(ny2,1)
              ELSE
                DO noy2=1,NONY(0,ny2,2)
                  no2=NONY(noy2,ny2,2)
                  CO2=CONY(noy2,ny2,2)
                  GKK(no1,no2)=GKK(no1,no2)+AA*CO1*CO2
                  IF(ny1.EQ.ny2)GKK(no1,no2)=GKK(no1,no2)+AB
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
     
      IF(IWRIT4.GE.4) THEN
        WRITE(OP_STRING,
     '	  '(/'' Global load vector GRR & stiffness matrix GKK:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NOT(1,nr,nx)='',I4)') NOT(1,nr,nx)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        DO no1=1,NOT(1,nr,nx)
          WRITE(OP_STRING,'(/'' GRR('',I4,'')='',D12.4,'' GKK: '','
     '	    //'8D12.4,(/25X,8D12.4))') no1,GRR(no1),(GKK(no1,no2),
     '	    no2=1,NOT(2,nr,nx))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      IFAIL=1
      IF(ITYP2(nr).EQ.3.OR.ITYP6(nr).EQ.2) THEN 
C ***   advection term or nonlin eqn
C ***   Approximate solution of non-symmetric GKK by Crout
        CALL F04ARF(GKK,NOM,GRR,NOT(1,nr,nx),XO,WK1,IFAIL)
        IF(ifail.ne.0) THEN
          WRITE(OP_STRING,*)'IFAIL= ',IFAIL
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ERROR=' IFAIL<>0 in F04ARF'
          GOTO 9999
        ENDIF
C        CALL NL_LINEAR_FULL(1,2,GKK,NOM,GRR,NOT(1,nr,nx),XO,GD,NYM,
C     '     IWK1,WK1,WK2,IFAIL,ERROR,*9999)

      ELSE IF(ITYP2(nr).ne.3.AND.ITYP6(nr).EQ.1) THEN 
C ***   linear symmetric eqns
C ***   Accurate solution of symmetric GKK by Cholesky
        CALL F04ASF(GKK,NOM,GRR,NOT(1,nr,nx),XO,WK1,WK2,IFAIL)
        IF(ifail.ne.0) THEN
          WRITE(OP_STRING,*)'IFAIL= ',IFAIL
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ERROR=' IFAIL<>0 in F04ASF'
          GOTO 9999
        ENDIF
C         CALL NL_LINEAR_FULL(1,1,GKK,NOM,GRR,NOT(1,nr,nx),XO,GD,NYM,
C     '     IWK1,WK1,WK2,IFAIL,ERROR,*9999)
      ENDIF

      IF(IWRIT4.GE.4) THEN
        FORMAT='(/'' Solution values:'')'
        WRITE(OP_STRING,FORMAT)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        FORMAT='('' XO('',I4,'')= '',D12.6)'
        DO no=1,NOT(1,nr,nx)
          WRITE(OP_STRING,FORMAT) no,XO(no)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      IF(KTYP90.EQ.6) THEN ! Linear boundary layer CPB/AJP 26/7/94
        DO ny=1,NYT(nrc,1,nx)
          DO noy=1,NONY(0,ny,1)
            no=NONY(noy,ny,1)
            CO=CONY(noy,ny,1)
            IF(noy.EQ.1) THEN
              IF(.NOT.FIX(ny,1).OR.((MIXED).AND..NOT.
     '          (FIX(ny,3).AND.(ITYP6(nr).EQ.1)))) YP(ny,1)=XO(no)
            ELSE
              ny2=1
              DO WHILE(ny2.LE.NYT(nrc,1,nx).AND..NOT.(NONY(1,ny2,2).EQ.
     '          no.and.ny2.ne.ny))
                ny2=ny2+1
              ENDDO
              IF(NONY(1,ny2,2).eq.no) THEN
                YP(ny2,1)=YP(ny,1)/CO+YP(ny2,1) ! u2 = u1/(1/a) + b
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ELSE
        DO ny=1,NYT(nrc,1,nx)
          SUM=0.0D0
          DO noy=1,NONY(0,ny,2)
            no=NONY(noy,ny,2)
            CO=CONY(noy,ny,2)
            SUM=SUM+XO(no)*CO
          ENDDO
          IF(.NOT.FIX(ny,1).OR.((MIXED).AND..NOT.
     '      (FIX(ny,3).AND.(ITYP6(nr).EQ.1)))) YP(ny,1)=SUM
        ENDDO
      ENDIF

      DO ny1=1,NYT(1,1,nx)
        SUM=0.0D0
        DO ny2=1,NYT(2,1,nx)
          SUM=SUM+GK(ny1,ny2)*YP(ny2,1)
        ENDDO
        YP(ny1,5)=SUM-GR(ny1)
        IF(FIX(ny1,3)) YP(ny1,5)=YP(ny1,5)+YP(ny1,3)/YP(ny1,2)
      ENDDO

      CALL EXITS('SOLVE1')
      RETURN
 9999 CALL ERRORS('SOLVE1',ERROR)
      CALL EXITS('SOLVE1')
      RETURN 1
      END


      SUBROUTINE SOLVE1C(IP,IBT,IDO,INP,LGE,
     '  NBH,NBJ,NEELEM,NHE,NHP,NJE,NKE,NKH,
     '  NONY,NP1OPT,NPE,NPF,
     '  NPNODE,NQE,nr,NRE,NVNE,NVNP,NW,nx,NYNE,NYNP,
     '  A,B,CE,CG,CONY,CP,ED,EM,ER,ES,
     '  GD,GK,GM,GR,PAOPTI,PBOPTI,PG,RG,SE,VE,
     '  WG,XA,XE,XG,XP,YG,YP,ZA,ZE,ZG,ZP,
     '  GKK,GRR,WK1,WK2,XO,DYNAM,FIX,SYMM,ERROR,*)

C**** Finds solution of a complex system of linear equations ny=1,NOT(nr).
C**** IP controls LINEAR/FIRST/UPDATE parameters as follows:
C**** IP=1 : linear=T; first=T; update=T  (linear eqtns with update)
C****  " 2     "    T    "   T         F  (  "      "   w.out  "   )
C****  " 3     "    F    "   T         T  (nonlin eqtns)
C**** SYMM indicates matrix symmetry. IWRIT4 controls diagnostic output.
C**** DYNAM is .true. when SOLVE1C called in time integration procedures.
C**** Element stiffness matrix ES & load vector ER for linear/nonlinear
C**** problems are generated by XPES/ZEES, respec.,  and assembled into
C**** the global stiffness matrix GKK & load vector GRR.
C**** Boundary element stiffness matrices are assembled directly into a
C**** matrix GK.
C**** In dynamic problems (DYNAM =.true.) the element mass, damping and
C**** stiffness matrices EM,ED & ES are assembled into global matrices
C**** GM,GD & GK which are then mapped into the solution matrix GKK.
C**** Note: These contain 2nd,1st & 0th order time-deriv coeffs,respec.
C**** The mapping coeffs NONY(noy,ny,nrc),CONY(noy,ny,nrc),
C**** noy=1,NONY(0,ny,nrc)
C**** are calc.d to reduce system of equations ny=1,NYT(nrc,1,nx) to the
C**** system no=1,NOT(nrc,nr,nx) by removing constraints, if UPDATE=.true.
C**** NONY(0,ny,nrc) is 0,1 when FIX(ny,1) is .true.,.false., resp.
C**** Note: In dynamic problems GKK is reformed at each time step since
C**** time step may change, so linear eqtn solution is done in one step
C**** either by Crout for non-symmetric problems or Cholesky for symm.
C**** KTYP22=0 : Static solution of initial accelerations.
C****   "    1 : first  order (Crank-Nicholson) time integration.
C****   "    2 : second order (Newmark; Wood-Zienkiewicz etc)
C****   "    3 : third  order (Houbolt; Hilbert-Hughes etc)
C**** YP(ny,1) contains essential b.c.s, defined by FIX(ny,1), on entry
C****                   solution (increm. if non.l. or dynam), on exit.
C****                or rhs vector if IP is 5 or 6.
C**** YP(ny,2)    "     nodal flux or stress bc.s defined by FIX(ny,2).
C**** YP(ny,3)    "     incremental b.c.s             "      FIX(ny,3).
C**** YP(ny,5)    "     reaction forces for linear solution, on exit.
C**** LGE(nhs,nrc) is location in global system of elem. var. nhs
C****   (=1,NHST(nrc))

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b08.cmn'
      INCLUDE 'cmiss$reference:b10.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:opti00.cmn'
      INCLUDE 'cmiss$reference:solv00.cmn'
!     Parameter List
      INTEGER IBT(2,NIM,*),IDO(NKM,0:NIM,*),INP(NNM,NIM,*),IP,
     '  LGE(NHM*NSM,*),NBH(NHM,NCM,*),NBJ(NJM,*),
     '  NEELEM(0:NEM,0:*),
     '  NHE(*),NHP(*),NJE(*),NKE(NKM,NNM,NBFM,*),NKH(NHM,NCM,*),
     '  NONY(0:NOYM,NYM,*),NP1OPT(*),NPE(NNM,NBFM,*),NPF(12,*),
     '  NPNODE(0:NPM,0:*),NQE(NSM,NBFM,*),nr,NRE(*),
     '  NVNE(NNM,NBFM,NJM,*),NVNP(NJM,NPM,*),
     '  NW(NEM,*),nx,NYNE(NAM,NHM,NRCM,NCM,*),
     '  NYNP(NKM,NVM,NHM,NPM,NRCM,*)
      REAL*8 A(NSM,*),B(NSM,*),CE(NMM,*),CG(NMM,*),CONY(0:NOYM,NYM,*),
     '  CP(NMM,*),ED(NVM,*),EM(NVM,*),ER(*),ES(NVM,*),
     '  GD(NYM,*),GK(NYM,*),GM(NYM,*),GR(*),
     '  PAOPTI(*),PBOPTI(*),
     '  PG(NSM,NUM,NGM,*),RG(*),
     '  SE(NSM,NBFM,*),VE(NSM,NKM,*),WG(NGM,*),
     '  XA(NAM,NJM,*),XE(NSM,*),XG(NJM,*),XP(NKM,NVM,NJM,*),
     '  YG(NGM,NJM,*),YP(NYM,*),ZA(NAM,NHM,NCM,*),ZE(NSM,*),
     '  ZG(NHM,*),ZP(NKM,NVM,NHM,NPM,*)
      REAL*8 GKK(NOM,*),GRR(*),WK1(*),WK2(*),XO(*)
      CHARACTER ERROR*(*)
      LOGICAL DYNAM,FIX(NYM,*),SYMM
!     Local Variables
      INTEGER IFAIL,IREPLY,nc,ne,no,no1,no2,noelem,noopti,
     '  noy,noy1,noy2,NP1,nhs,nhs1,nhs2,NHST(2),nrc,nv,ny,ny1,ny2
      REAL*8 AA,BB,CO,CO1,CO2,DPHOUT,MULT(3),PHIIN,PHIOUT,SUM,TOL,
     '  TRANLS
      CHARACTER FORMAT*500
      LOGICAL FIRST,LINEAR,LOGIC(6,3),UPDATE
      DATA LOGIC/2*.TRUE.,4*.FALSE.,3*.TRUE.,3*.FALSE.,.TRUE.,.FALSE.,
     '           2*.TRUE.,2*.FALSE./
      DATA MULT /0.2618D0,1.3090D0,1.5708D0/, TOL /1.0D-5/

      CALL ENTERS('SOLVE1C',*9999)

      nv=1 ! Temporary cpb 24/11/94
      nc=1 ! Temporary cpb 1/12/94

      CALL ASSERT(NOM*NOT(1,nr,nx).LE.NZM,'>>NZM too small',ERROR,*9999)
      LINEAR=LOGIC(IP,1)
      FIRST =LOGIC(IP,2)
      UPDATE=LOGIC(IP,3)
      IF(IWRIT4.GT.0) THEN
        FORMAT='(/'' >SOLVE1C   linear='',L1,'' first='',L1,'
     '    //''' symm='''//',L1,'' update='',L1)'
        WRITE(OP_STRING,FORMAT) LINEAR,FIRST,SYMM,UPDATE
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      no=0
      DO nrc=1,2
        DO ny=1,NYT(nrc,1,nx)
          IF(FIX(ny,1)) THEN
            NONY(0,ny,nrc)=0
          ELSE
            NONY(0,ny,nrc)=1
          ENDIF
          DO noy=1,NONY(0,ny,nrc)
            no=no+1
            NONY(noy,ny,nrc)=no
            CONY(noy,ny,nrc)=1.0D0
          ENDDO
        ENDDO
      ENDDO
      NOT(2,nr,nx)=no
      IF(NOT(2,nr,nx).EQ.0) THEN
        ERROR=' >>The number of unknowns is zero'
        GO TO 9999
      ENDIF
      DO no1=1,NOT(1,nr,nx)
        GRR(no1)=0.0d0
        DO no2=1,NOT(2,nr,nx)
          GKK(no1,no2)=0.0d0
        ENDDO
      ENDDO

      nrc=2 ! temporary ajp 30/11/94
      DO ny=1,NYT(nrc,1,nx)
        IF(FIX(ny,2)) THEN
          GR(ny)=YP(ny,2)
        ELSE
          GR(ny)=0.0D0
        ENDIF
      ENDDO

      IF(KTYP27.EQ.4) THEN
C ***   Update nodal coordinates from optimisation params
        DO noopti=1,NTOPTI
          NP1=NP1OPT(noopti)
          XP(1,nv,2,NP1)=PAOPTI(noopti)/PBOPTI(noopti)*XP(1,nv,2,NP1)
          PBOPTI(noopti)=PAOPTI(noopti)
        ENDDO
      ENDIF

      IF(UPDATE) THEN
        DO ny1=1,NYT(1,1,nx)
          DO ny2=1,NYT(2,1,nx)
            GK(ny1,ny2)=0.0D0
            IF(DYNAM) THEN
              GM(ny1,ny2)=0.0D0
              GD(ny1,ny2)=0.0D0
            ENDIF
          ENDDO
        ENDDO

        DO noelem=1,NEELEM(0,nr) !is main element loop
          ne=NEELEM(noelem,nr)
          IF(NW(ne,1).GT.0) THEN

            CALL MELGE(LGE,NBH(1,1,ne),nc,NHE(ne),NHST,NKH,NPE(1,1,ne),
     '        NVNE(1,1,1,ne),NYNE(1,1,1,nc,ne),NYNP(1,1,1,1,1,nc),
     '        ERROR,*9999)

            DO nhs1=1,NHST(1)    !Initialize element arrays
              ER(nhs1)=0.0d0
              DO nhs2=1,NHST(2)
                ES(nhs1,nhs2)=0.0d0
              ENDDO
            ENDDO

            IF(NW(ne,1).LE.20) THEN !finite element
              IF(ITYP1(nr).EQ.3) THEN !partial differential equation
                CALL XPXE(NBJ(1,ne),NJE(ne),NKE(1,1,1,ne),NPE(1,1,ne),
     '            NPF(1,1),NQE(1,1,ne),nr,NVNE(1,1,1,ne),NVNP,
     '            SE(1,1,ne),XA,XE,XP,ERROR,*9999)
                CALL XPES30(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),ne,
     '            NHE(ne),NJE(ne),NPE(1,1,ne),nr,
     '            CE(1,ne),CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '            VE(1,1,ne),WG,XE,XG,ZE,ZG,UPDATE,ERROR,*9999)
              ELSE IF(ITYP1(nr).EQ.4) THEN !linear elasticity
                CALL XPES40(NBH(1,1,ne),NBJ(1,ne),
     '            NHE(ne),NJE(ne),NKE(1,1,1,ne),
     '            NPE(1,1,ne),NPF,NQE(1,1,ne),nr,NVNE(1,1,1,ne),
     '            NVNP,NW(ne,1),
     '            CE(1,ne),CG,CP,ED,EM,ER,ES,
     '            PG,RG,SE(1,1,ne),VE(1,1,ne),WG,
     '            XA,XE,XG,XP,YG(1,1,ne),
     '            UPDATE,ERROR,*9999)
              ENDIF
            ELSE !boundary element
              ERROR=' >>Bdry elements should use SOLVE4'
              GO TO 9999
            ENDIF

            IF(IWRIT4.GE.3) THEN
              WRITE(OP_STRING,
     '          '(/'' Element load vector ER & stiffness matrix ES:'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              DO nhs1=1,NHST(1)
                WRITE(OP_STRING,'('' ER('',I2,'')='',E12.4,'' ES: '','
     '            //'8E12.4/,(25X,8E12.4))') 
     '            nhs1,ER(nhs1),(ES(nhs1,nhs2),nhs2=1,NHST(2))
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDDO
              IF(DYNAM) THEN
                WRITE(OP_STRING,'(/'' Element matrices EM & ED:'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                DO nhs1=1,NHST(1)
                  WRITE(OP_STRING,'('' nhs1='',I2,'' EM: '',8E12.4/,'
     '              //'(12X,8E12.4))') nhs1,(EM(nhs1,nhs2),nhs2=1,
     '              NHST(2))
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDDO
                DO nhs1=1,NHST(1)
                  WRITE(OP_STRING,'('' nhs1='',I2,'' ED: '',8E12.4/,'
     '              //'(12X,8E12.4))') nhs1,(ED(nhs1,nhs2),nhs2=1,
     '              NHST(2))
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDDO
              ENDIF
            ENDIF

            DO nhs1=1,NHST(1)
              ny1=IABS(LGE(nhs1,1))
              GR(ny1)=GR(ny1)+ER(nhs1)
              DO nhs2=1,NHST(2)
                ny2=IABS(LGE(nhs2,2))
                GK(ny1,ny2)=GK(ny1,ny2)+ES(nhs1,nhs2)
                IF(DYNAM) THEN
                  GM(ny1,ny2)=GM(ny1,ny2)+EM(nhs1,nhs2)
                  GD(ny1,ny2)=GD(ny1,ny2)+ED(nhs1,nhs2)
                ENDIF
              ENDDO
            ENDDO

          ENDIF
        ENDDO

        IF(IWRIT4.GE.4) THEN
          WRITE(OP_STRING,'(/'' Global load vector GR & stiffness'
     '      //' matrix GK:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DO ny1=1,NYT(1,1,nx)
            WRITE(OP_STRING,'('' GR('',I4,'')='',E12.4,'' GK: '','
     '	      //'8E12.4'//'/,(25X,8E12.4))') 
     '        ny1,GR(ny1),(GK(ny1,ny2),ny2=1,NYT(2,1,nx))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDDO
          IF(DYNAM) THEN
            WRITE(OP_STRING,'(/'' Global matrices GM & GD:'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            DO ny1=1,NYT(1,1,nx)
              WRITE(OP_STRING,'('' NY1='',I2,'' GM: '',8E12.4'
     '          //'/,(12X,8E12.4))') 
     '          ny1,(GM(ny1,ny2),ny2=1,NYT(2,1,nx))
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDDO
            DO ny1=1,NYT(1,1,nx)
              WRITE(OP_STRING,'('' NY1='',I2,'' GD: '',8E12.4'
     '          //'/,(12X,8E12.4))') 
     '          ny1,(GD(ny1,ny2),ny2=1,NYT(2,1,nx))
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF
        ENDIF

      ELSE IF(.NOT.UPDATE) THEN
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(NW(ne,1).GT.0) THEN

            CALL MELGE(LGE,NBH(1,1,ne),nc,NHE(ne),NHST,NKH,NPE(1,1,ne),
     '        NVNE(1,1,1,ne),NYNE(1,1,1,nc,ne),NYNP(1,1,1,1,1,nc),
     '        ERROR,*9999)

            CALL XPXE(NBJ(1,ne),NJE(ne),NKE(1,1,1,ne),NPE(1,1,ne),
     '        NPF(1,1),NQE(1,1,ne),nr,NVNE(1,1,1,ne),NVNP,SE(1,1,ne),
     '        XA,XE,XP,ERROR,*9999)
            CALL YPZP(4,NBH,NEELEM,
     '        NHE,NHP,NKH,NPNODE,nr,NVNP,NYNE,NYNP,
     '        YP,ZA,ZP,ERROR,*9999)
            CALL ZPZE(1,NBH(1,1,ne),NHE(ne),NKE(1,1,1,ne),NPE(1,1,ne),
     '        NPF(1,1),nr,NVNE(1,1,1,ne),NVNP,NW(ne,1),SE(1,1,ne),
     '        ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
            CALL XPES30(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),ne,
     '        NHE(ne),NJE(ne),NPE(1,1,ne),nr,
     '        CE(1,ne),CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '        VE(1,1,ne),WG,XE,XG,ZE,ZG,UPDATE,ERROR,*9999)

            DO nhs=1,NHST(1)
              ny=IABS(LGE(nhs,1))
              GR(ny)=GR(ny)+ER(nhs)
            ENDDO

          ENDIF
        ENDDO
      ENDIF

      DO ny1=1,NYT(1,1,nx)
        DO noy1=1,NONY(0,ny1,1)
          no1=NONY(noy1,ny1,1)
          CO1=CONY(noy1,ny1,1)
          BB=GR(ny1)
          DO ny=1,NYT(2,1,nx)
            IF(KTYP22.GE.1) BB=BB-GK(ny1,ny)*YP(ny,10)
            IF(KTYP22.GE.2) BB=BB-GD(ny1,ny)*YP(ny,11)
            IF(KTYP22.EQ.3) BB=BB-GM(ny1,ny)*YP(ny,12)
          ENDDO
          GRR(no1)=GRR(no1)+BB*CO1
          DO ny2=1,NYT(2,1,nx)
            IF(.NOT.DYNAM) THEN
              AA=GK(ny1,ny2)
            ELSE IF(DYNAM) THEN
              IF(KTYP22.EQ.0) THEN
                AA=GM(ny1,ny2)
              ELSE IF(KTYP22.EQ.1) THEN
                AA=GD(ny1,ny2)+A1*GK(ny1,ny2)
              ELSE IF(KTYP22.EQ.2) THEN
                AA=GM(ny1,ny2)+A1*GD(ny1,ny2)+A2*GK(ny1,ny2)
              ELSE IF(KTYP22.EQ.3) THEN
                AA=A1*GM(ny1,ny2)+A2*GD(ny1,ny2)+A3*GK(ny1,ny2)
              ENDIF
            ENDIF
            IF(NONY(0,ny2,2).EQ.0) THEN
              GRR(no1)=GRR(no1)-AA*CO1*YP(ny2,1)
            ELSE
              DO noy2=1,NONY(0,ny2,2)
                no2=NONY(noy2,ny2,2)
                CO2=CONY(noy2,ny2,2)
                GKK(no1,no2)=GKK(no1,no2)+AA*CO1*CO2
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      IF(IWRIT4.GE.4) THEN
        WRITE(OP_STRING,
     '    '(/'' Global load vector GRR & stiffness matrix GKK:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NOT(1,nr,nx)='',I4)') NOT(1,nr,nx)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        DO no1=1,NOT(1,nr,nx)
          WRITE(OP_STRING,'(/'' GRR('',I4,'')='',D12.4,'' GKK: '','
     '	    //'8D12.4,(/25X,8D12.4))') no1,GRR(no1),(GKK(no1,no2),
     '	    no2=1,NOT(2,nr,nx))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      IFAIL=1
      CALL F04ATF(GKK,NOM,GRR,NOT(1,nr,nx),XO,GD,NOM,WK1,WK2,IFAIL)
      IF(ifail.ne.0) THEN
        WRITE(OP_STRING,*)'IFAIL= ',IFAIL
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ERROR=' IFAIL<>0 in F04ATF'
        GOTO 9999
      ENDIF

      IF(IWRIT4.GE.4) THEN
        FORMAT='(/'' Solution values:'')'
        WRITE(OP_STRING,FORMAT)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        FORMAT='('' XO('',I4,'')= '',D12.6)'
        DO no=1,NOT(1,nr,nx)
          WRITE(OP_STRING,FORMAT) no,XO(no)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

C *** Write solutions to the YP array
      DO ny=1,NYT(nrc,1,nx)
        SUM=0.0D0
        DO noy=1,NONY(0,ny,2)
          no=NONY(noy,ny,2)
          CO=CONY(noy,ny,2)
          SUM=SUM+XO(no)*CO
        ENDDO
        IF(.NOT.FIX(ny,1)) YP(ny,1)=SUM
      ENDDO

      DO ny1=1,NYT(1,1,nx)
        SUM=0.0D0
        DO ny2=1,NYT(2,1,nx)
          SUM=SUM+GK(ny1,ny2)*YP(ny2,1)
        ENDDO
        YP(ny1,5)=SUM-GR(ny1)
      ENDDO

C***  This is specific to the silencer with 3 outlet nodes as silencer1 (M.P.N.)
C***  This code is very crude and needs to be generalised to remain permanent.
      WRITE(OP_STRING,*)
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,*)' Calculate: 1) Pole Parameters B*S/(i*z0),D'
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,*)
     '  '            2) Pole Parameters A,C*z0/(i*S) & Trans. Loss'
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,*)'            3) Nothing ?'
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      READ(*,*)IREPLY

      IF(IREPLY.EQ.1) THEN
C***    Find inlet and outlet nodes
        INLET(0)=0
        OUTLET(0)=0
        DO ny=1,NYT(nrc,1,nx)
          IF(FIX(ny,2).AND.(DABS(YP(ny,2)).GT.1.E-5)) THEN
            INLET(0)=INLET(0)+1
            INLET(INLET(0))=ny
          ELSE IF(FIX(ny,1).AND.(YP(ny,1).EQ.0.0D0)) THEN
            OUTLET(0)=OUTLET(0)+1
            OUTLET(OUTLET(0))=ny
          ENDIF
        ENDDO

C***    Find average value for the inlet potential
        PHIIN=0.0D0
        DO no=1,INLET(0)
          PHIIN=PHIIN+YP(INLET(no),1)
        ENDDO
        IF(INLET(0).ne.0) PHIIN=PHIIN/DBLE(INLET(0))
        WRITE(OP_STRING,*)' PHIIN=',PHIIN
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C***    Find average value for the outlet velocity
C***    WARNING: MULT must contain the values of the multipliers of the natural
C***    boundary conditions (obtained from the integration of the linear basis
C***    functions.) Here MULT outlet contains 3 nodes hence MULT has 3 values.
        DPHOUT=0.0D0
        DO no=1,OUTLET(0)
          DPHOUT=DPHOUT+YP(OUTLET(no),5)/MULT(no)
        ENDDO
        IF(OUTLET(0).ne.0) DPHOUT=DPHOUT/DBLE(OUTLET(0))
        WRITE(OP_STRING,*)' DPHOUT=',DPHOUT
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C***   Calculate two of the four pole parameters (ie B and D)
        IF(DABS(DPHOUT).GT.TOL) THEN
          BPOLE=-CE(1,1)*PHIIN/DPHOUT
          DPOLE=-1.0D0/DPHOUT
          WRITE(OP_STRING,*)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,*)' B*S/(i*z0)=',BPOLE
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,*)' D=',DPOLE
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,*)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,*)' Zero Outlet Velocity!'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF((IREPLY.EQ.2).AND.((DABS(BPOLE).GT.TOL).OR.
     '  (DABS(DPOLE).GT.TOL))) THEN
C***    Find average value for the inlet potential
        PHIIN=0.0D0
        DO no=1,INLET(0)
          PHIIN=PHIIN+YP(INLET(no),1)
        ENDDO
        IF(INLET(0).ne.0) PHIIN=PHIIN/DBLE(INLET(0))
        WRITE(OP_STRING,*)' PHIIN=',PHIIN
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C***    Find average value for the outlet potential
        PHIOUT=0.0D0
        DO no=1,OUTLET(0)
          PHIOUT=PHIOUT+YP(OUTLET(no),1)
        ENDDO
        IF(OUTLET(0).ne.0) PHIOUT=PHIOUT/DBLE(OUTLET(0))
        WRITE(OP_STRING,*)' PHIOUT=',PHIOUT
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C***   Calculate two of the four pole parameters (ie A and C)
        IF(DABS(PHIOUT).GT.TOL) THEN
          APOLE=PHIIN/PHIOUT
          CPOLE=-1.0D0/(PHIOUT*CE(1,1))
          WRITE(OP_STRING,*)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,*)' A=',APOLE
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,*)' B*S/(i*z0)=',BPOLE
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,*)' C*z0/(i*S)=',CPOLE
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,*)' D=',DPOLE
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C***      Calculate the transmission loss
          TRANLS=20*DLOG10(0.5*DSQRT((APOLE+DPOLE)**2+(BPOLE+CPOLE)**2))
          WRITE(OP_STRING,*)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,*)' Transmission Loss = ',TRANLS
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,*)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,*)' Zero Outlet Potential!'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF((IREPLY.EQ.2).AND.((DABS(BPOLE).LE.TOL).AND.
     '  (DABS(DPOLE).LE.TOL))) THEN
        WRITE(OP_STRING,*)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*)' Need to calculate B*S/(i*z0)'
     '    //' and D first.'
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      ENDIF

      CALL EXITS('SOLVE1C')
      RETURN
 9999 CALL ERRORS('SOLVE1C',ERROR)
      CALL EXITS('SOLVE1C')
      RETURN 1
      END


      SUBROUTINE SOLVE2(IP,IBT,ICN,ICN2,IDO,INP,IRN,IRN2,IWK1,IWK2,
     '  LGE,NAN,NBH,NBJ,NEELEM,NFF,NGAP,
     '  NHE,NHP,NJE,NKE,NKH,
     '  NLL,NNF,NNL,NONY,NPE,NPF,NPL,
     '  NPNODE,NQE,nr,NRE,NVNE,NVNP,NW,nx,NYNE,
     '  NYNP,CE,CG,CONY,CP,D_RE,D_RI3,D_TG,D_ZG,
     '  ER,ES,FEXT,GK,GR,PF,PG,RE1,RE2,RG,
     '  SE,VE,WG,XA,XE,XF,XG,XP,YP,ZA,ZE,ZE1,ZF,ZG,ZG1,ZP,
     '  GKK,GRR,WK1,FIX,FNY,ERROR,*)

C**** Finds solution of system of unsymm sparse linear equations
C****   no=1,NOT(nrc,nr,nx).
C**** IP controls LINEAR/FIRST/UPDATE parameters as follows:
C**** IP=1 : linear=T; first=T; update=T  (linear eqtns with update)
C****  " 2     "    T    "   T         F  (  "      "   w.out  "   )
C****  " 3     "    F    "   T         T  (nonlin eqtns,load increm)
C****  " 4     "    F    "   F         T  (  "      "   full Newton)
C****  " 5     "    F    "   F         F  (  "      "   modified nr)
C****  " 6     "    F    "   F         F  (  "      "   quasi-Newtn)
C**** IWRIT4 controls diagnostic output.
C**** Element stiffness matrix ES & load vector ER for nonlinear
C**** problems are generated by ZEES, and assembled into
C**** the global stiffness matrix GK & load vector GR, and then into the
C**** constraint reduced matrices (non-zeros only) GKK and GRR.
C**** The mapping coeffs NONY(noy,ny,nrc),CONY(noy,ny,nrc),
C**** noy=1,NONY(0,ny,nrc)
C**** are calc.d to reduce system of equations ny=1,NYT(nrc,1,nx) to the
C**** system no=1,NOT(nrc,nr,nx) by removing constraints, if UPDATE=.true.
C**** NONY(0,ny,nrc) is 0,1 when FIX(ny,1) is .true.,.false., resp.
C**** YP(ny,1) contains essential b.c.s, defined by FIX(ny,1), on entry
C****                   incremental solution, on exit.
C****                or rhs vector if IP is 5 or 6.
C**** YP(ny,2)    "     nodal flux or stress bc.s defined by FIX(ny,2).
C**** YP(ny,3)    "     incremental b.c.s             "      FIX(ny,3).
C**** LGE(nhs,nrc) is location in global system of elem. var. nhs (=1,
C****   NHST(nrc))
C**** IRN(nz) is row index of the non-zero stored in GKK(nz),nz=1,NZT.
C**** ICN(nz) is col   "      "      "        "       "        "
C**** FNY(ny1,ny2) is .true. if corresponding global value is non-zero.
      
      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b08.cmn'
      INCLUDE 'cmiss$reference:b10.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:gen000.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
!     Parameter List
      INTEGER IBT(2,NIM,*),ICN(*),ICN2(*),IDO(NKM,0:NIM,*),
     '  INP(NNM,NIM,*),IP,IRN(*),IRN2(*),IWK1(*),IWK2(*),LGE(NHM*NSM,*),
     '  NAN(NIM,NAM,*),NBH(NHM,NCM,*),NBJ(NJM,*),
     '  NEELEM(0:NEM,0:*),NFF(6,*),NGAP(NIM,*),NHE(*),
     '  NHP(*),NJE(*),NKE(NKM,NNM,NBFM,*),NKH(NHM,NCM,*),
     '  NLL(12,*),NNF(0:14,6,*),NNL(4,12,*),
     '  NONY(0:NOYM,NYM,*),NPE(NNM,NBFM,*),NPF(12,*),NPL(20,*),
     '  NPNODE(0:NPM,0:*),NQE(NSM,NBFM,*),
     '  nr,NRE(*),
     '  NVNE(NNM,NBFM,NJM,*),NVNP(NJM,NPM,*),
     '  NW(NEM,*),nx,NYNE(NAM,NHM,NRCM,NCM,*),
     '  NYNP(NKM,NVM,NHM,NPM,NRCM,*)
      REAL*8 CE(NMM,*),CG(NMM,*),CONY(0:NOYM,NYM,*),CP(NMM,*),
     '  D_RE(NSM,NHM,*),D_RI3(*),D_TG(3,3,*),D_ZG(NHM,NUM,*),
     '  ER(*),ES(NVM,*),FEXT(8,NGM,*),
     '  GK(NYM,*),GR(*),PF(2,*),PG(NSM,NUM,NGM,*),RE1(NSM,*),
     '  RE2(NSM,*),RG(*),SE(NSM,NBFM,*),VE(NSM,NKM,*),WG(NGM,*),
     '  XA(NAM,NJM,*),XE(NSM,*),XF(NSM,*),XG(NJM,*),XP(NKM,NVM,NJM,*),
     '  YP(NYM,*),ZA(NAM,NHM,NCM,*),ZE(NSM,*),ZE1(NSM,*),ZF(NSM,*),
     '  ZG(NHM,*),ZG1(NHM,*),ZP(NKM,NVM,NHM,NPM,*)
      REAL*8 GKK(*),GRR(*),WK1(*)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,*),FNY(NYM,*)
!     Local Variables
      INTEGER IDISP(10),IFAIL,nc,ne,no,no1,no2,noelem,noy,noy1,noy2,
     '  nhs1,nhs2,NHST(2),nrc,ny,ny1,ny2,nz
      REAL*8 AA,CO,CO1,CO2,RESID1,RPMIN,SUM
      CHARACTER FORMAT*500
      LOGICAL ABORT(4),ELEM,GROW,LBLOCK,LINEAR,LOGIC(6,3),UPDATE
      DATA LOGIC/2*.TRUE.,4*.FALSE.,3*.TRUE.,3*.FALSE.,.TRUE.,.FALSE.,
     '           2*.TRUE.,2*.FALSE./
      DATA ABORT/2*.TRUE.,.FALSE.,.TRUE./,GROW/.TRUE./,LBLOCK/.TRUE./

      CALL ENTERS('SOLVE2',*9999)

      nc=1 ! temporary cpb 1/12/94

      LINEAR=LOGIC(IP,1)
      UPDATE=LOGIC(IP,3)
      IF(IWRIT4.GT.0) THEN
        FORMAT='(/'' >SOLVE2  linear='',L1,'' first='',L1,'
     '    //''' update='',L1)'
        WRITE(OP_STRING,FORMAT) LINEAR,FIRSTS,UPDATE
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      IF(UPDATE) THEN
        IF(FIRSTS) THEN
          DO ny1=1,NYT(1,1,nx)
            DO ny2=1,NYT(2,1,nx)
              FNY(ny1,ny2)=.FALSE.
            ENDDO
          ENDDO
          no=0
          DO nrc=1,2
            DO ny=1,NYT(nrc,1,nx)
              IF(FIX(ny,1)) THEN
                NONY(0,ny,nrc)=0
              ELSE
                NONY(0,ny,nrc)=1
              ENDIF
              DO noy=1,NONY(0,ny,nrc)
                no=no+1
                NONY(noy,ny,nrc)=no
                CONY(noy,ny,nrc)=1.0D0
              ENDDO
            ENDDO
          ENDDO
          NOT(2,nr,nx)=no
        ENDIF
        CALL ASSERT(NOT(2,nr,nx).GT.0,
     '    '>>Number of unknowns is zero',ERROR,*9999)
        CALL ASSERT(NOT(2,nr,nx).LE.NOM,
     '    '>>NOM too small (sb >NOT(nr,nx))',ERROR,*9999)

        DO ny1=1,NYT(1,1,nx)
          GR(ny1)=0.0D0 !ny1 not correct for GR
          DO ny2=1,NYT(2,1,nx)
            GK(ny1,ny2)=0.0D0
          ENDDO
        ENDDO

        DO noelem=1,NEELEM(0,nr) !is main element loop
          ne=NEELEM(noelem,nr)
          IF(NW(ne,1).GT.0) THEN

            CALL MELGE(LGE,NBH(1,1,ne),nc,NHE(ne),NHST,NKH,NPE(1,1,ne),
     '        NVNE(1,1,1,ne),NYNE(1,1,1,nc,ne),NYNP(1,1,1,1,1,nc),
     '        ERROR,*9999)

            ELEM=.FALSE.
            DO nhs1=1,NHST(1)   !Initialize element arrays
              ER(nhs1)=0.0d0
              ny=IABS(LGE(nhs1,1))
              IF(.NOT.FIX(ny,1)) ELEM=.TRUE.
              DO nhs2=1,NHST(2)
                ES(nhs1,nhs2)=0.0D0
              ENDDO
            ENDDO

            IF(ELEM) THEN
              IF(LINEAR) THEN
                ERROR=' >>Linear equations should use SOLVE1'
                GO TO 9999

              ELSE IF(.NOT.LINEAR) THEN
                CALL XPXE(NBJ,NJE,NKE(1,1,1,ne),NPE(1,1,ne),NPF(1,1),
     '            NQE(1,1,ne),nr,NVNE(1,1,1,ne),NVNP,SE(1,1,ne),XA,XE,
     '            XP,ERROR,*9999)
                CALL ZPZE(1,NBH,NHE,NKE(1,1,1,ne),
     '            NPE(1,1,ne),NPF(1,1),nr,NVNE(1,1,1,ne),NVNP,NW,
     '            SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
                CALL ZEES(IBT,IDO,INP,LGE,NAN,NBH(1,1,ne),NBJ(1,ne),
     '            ne,NFF(1,ne),NHE(ne),NJE(ne),
     '            NNF,NPE,NPF,NPL,nr,NW(ne,1),
     '            CE(1,ne),CG,CP,D_RE,D_RI3,D_TG,D_ZG,ES,FEXT(1,1,ne),
     '            PF(1,ne),PG,RE1,RE2,RG,SE,VE(1,1,ne),WG,
     '            XE,XG,ZE,ZE1,ZG,ZG1,FIX,ERROR,*9999)

                IF(iwrit4.ge.3.and.NW(ne,1).LE.20) THEN
                  WRITE(OP_STRING,
     '              '(/'' Element load vector ER & stiffness matrix'
     '              //' ES:'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  DO nhs1=1,NHST(1)
                    WRITE(OP_STRING,'('' ER('',I2,'')='',E12.4,'
     '                //''' ES: '',8E12.4/,(25X,8E12.4))') 
     '                nhs1,ER(nhs1),(ES(nhs1,nhs2),nhs2=1,NHST(2))
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDDO
                ENDIF

                DO nhs1=1,NHST(1)
                  ny1=IABS(LGE(nhs1,1))
                  GR(ny1)=GR(ny1)+ER(nhs1)
                  DO nhs2=1,NHST(2)
                    ny2=IABS(LGE(nhs2,2))
                    FNY(ny1,ny2)=.TRUE.
                    GK(ny1,ny2)=GK(ny1,ny2)+ES(nhs1,nhs2)
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
          ENDIF
        ENDDO

        IF(IWRIT4.GE.4) THEN
          WRITE(OP_STRING,'(/'' Global rhs vector GR & stiffness'
     '      //' matrix GK:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DO ny1=1,NYT(1,1,nx)
            WRITE(OP_STRING,'('' GR('',I2,'')='',E12.4,'' GK: '','
     '	      //'8E12.4'//'/,(25X,8E12.4))') 
     '        ny1,GR(ny1),(GK(ny1,ny2),ny2=1,NYT(2,1,nx))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF
      ENDIF

      IF(.NOT.LINEAR) THEN
        CALL ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBJ,
     '    NEELEM,NFF,NHE,NHP,NJE,
     '    NKE,NKH,NNF,NPE,NPF,
     '    NPNODE,NQE,nr,NRE,NVNE,NVNP,NW,NYNE,
     '    NYNP,CE,CG,CP,FEXT,PF,PG,RE1,RG,GR,SE,VE,WG,
     '    XA,XE,XG,XP,YP,ZA,ZE,ZG,ZP,FIX,ERROR,*9999)
        DO ny=1,NYT(1,1,nx)
          GR(ny)=-GR(ny)
        ENDDO
      ENDIF

C **  Calculate constrained-reduced 1D global matrix GKK
      DO no=1,NOT(1,nr,nx)
        GRR(no)=0.0D0
      ENDDO

      nz=0
      DO ny1=1,NYT(1,1,nx)
        DO noy1=1,NONY(0,ny1,1)
          no1=NONY(noy1,ny1,1)
          CO1=CONY(noy1,ny1,1)
          GRR(no1)=GRR(no1)+GR(ny1)*CO1
          DO ny2=1,NYT(2,1,nx)
            AA=GK(ny1,ny2)
            IF(NONY(0,ny2,2).EQ.0) THEN
              GRR(no1)=GRR(no1)-AA*CO1*YP(ny2,1)
            ELSE IF(UPDATE) THEN
              DO noy2=1,NONY(0,ny2,2)
                no2=NONY(noy2,ny2,2)
                CO2=CONY(noy2,ny2,2)
                IF(FNY(ny1,ny2)) THEN
                  nz=nz+1
                  GKK(nz)=AA*CO1*CO2
                  IF(FIRSTS) THEN
                    IRN(nz) =no1
                    IRN2(nz)=no1
                    ICN(nz) =no2
                    ICN2(nz)=no2
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NZT=nz
      IF(NZT.GT.NZ_GKK_M) THEN
        FORMAT='('' Increase array dimensions: NZT='',I6,'
     '    //''' NZ_GKK_M='',I6)'
        WRITE(OP_STRING,FORMAT) NZT,NZ_GKK_M
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      IF(IWRIT4.GE.4) THEN
        WRITE(OP_STRING,'(/'' NOT(1,nr,nx)='',I4,'' NZT='',I6)') 
     '    NOT(1,nr,nx),NZT
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '(/'' Global rhs vector GRR & stiffness matrix GKK:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '(/'' GRR: '',8D12.4/,(4X,8D12.4))') (GRR(no),no=1,
     '    NOT(1,nr,nx))
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '(/'' IRN: '',  20I5/,(6X,20I5))')   (IRN(nz),nz=1,NZT)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '(/'' ICN: '',  20I5/,(6X,20I5))')   (ICN(nz),nz=1,NZT)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '(/'' GKK: '',8D12.4/,(4X,8D12.4))') (GKK(nz),nz=1,NZT)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      IF(UPDATE) THEN
        IF(FIRSTS) THEN
          IFAIL=111
          CALL F01BRF(NOT(1,nr,nx),NZT,GKK,NZ_GKK_M,IRN,NZ_GKK_M,ICN,
     '      0.3D0,IWK1,IWK2,WK1,LBLOCK,GROW,ABORT,IDISP,IFAIL)
          IF(IWRIT4.GT.0) THEN
            WRITE(OP_STRING,'('' F01BRF: WK1(1)='',D12.4)') WK1(1)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          IF(ifail.ne.0) THEN
            WRITE(OP_STRING,'('' IFAIL='',I3,'' in F01BRF'''
     '	      //')') IFAIL
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            GO TO 9999
          ENDIF
          FIRSTS=.FALSE.
        ELSE
          IFAIL=111
          CALL F01BSF(NOT(1,nr,nx),NZT,GKK,NZ_GKK_M,IRN2,ICN2,ICN,
     '      IWK1,IWK2,WK1,GROW,1.D-4,RPMIN,ABORT(1),IDISP,IFAIL)
          IF(IWRIT4.GT.0) THEN
            WRITE(OP_STRING,'('' F01BSF: WK1(1)='',D12.4,'
     '	      //''' RPMIN='',D12.4)') WK1(1),RPMIN
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          IF(ifail.ne.0) THEN
            WRITE(OP_STRING,'('' IFAIL='',I3,'' in F01BSF'''
     '	      //')') IFAIL
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            GO TO 9999
          ENDIF
          IF(WK1(1).GT.1.0D4) THEN
            FIRSTS=.TRUE.
          ELSE
            FIRSTS=.FALSE.
          ENDIF
        ENDIF
      ENDIF
      CALL F04AXF(NOT(1,nr,nx),GKK,NZ_GKK_M,ICN,IWK1,GRR,WK1,1,IDISP,
     '  RESID1)
      IF(ifail.ne.0) THEN
        WRITE(OP_STRING,'('' IFAIL='',I3,'' in F04AXF'')') IFAIL
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GO TO 9999
      ENDIF

      IF(IWRIT4.GE.4) THEN
        FORMAT='(/'' Solution values:'')'
        WRITE(OP_STRING,FORMAT)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        FORMAT='('' GRR('',I4,'')= '',D12.6)'
        DO no=1,NOT(1,nr,nx)
          WRITE(OP_STRING,FORMAT) no,GRR(no)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      DO ny=1,NYT(1,1,nx)
        SUM=0.0D0
        DO noy=1,NONY(0,ny,1)
          no=NONY(noy,ny,1)
          CO=CONY(noy,ny,1)
          SUM=SUM+GRR(no)*CO
        ENDDO
        IF(.NOT.FIX(ny,1)) YP(ny,1)=SUM
      ENDDO

 9998 CALL EXITS('SOLVE2')
      RETURN
 9999 CALL ERRORS('SOLVE2',ERROR)
      CALL EXITS('SOLVE2')
      RETURN 1
      END


C CPB 28/3/96 This is the old assemble3 before I rewrote it.
      SUBROUTINE SOLVE3(IDISP,IP,ISC_GKK,ISC2_GKK,ISR_GKK,ISR2_GKK,
     '  IWK1,IWK2,NONY,NPNY,nr,nx,NYNE,NYNO,NYNP,NYNR,NZNY,
     '  CONY,CYNO,GD,GK,GKK,GM,GR,GRR,WK1,XO,YP,DYNAM1,
     '  DYNAM2,FIX,ERROR,*)

C#### Subroutine: SOLVE3
C###  Description:
C###    SOLVE3 solves advection diffusion equations with constant 
C###    time step.  

C**** Returns solution and not increment as in previous SOLVE3. 
C**** Still uses outdated spares structures.  
C**** Matrices are assembled in ASSEMBLE3 and just reorganised here.
C**** Other comments from old SOLVE3 (some of which may not be true):
C****
C**** Finds solution of system of unsymm sparse linear eqns no=1,
C**** NOT(nrc,nr,nx) in time dependent problems.
C**** IP controls LINEAR/FIRST/UPDATE parameters as follows:
C**** IP=1: linear=T; first=T; update=T  (first call       with  update)
C****  " 2    "    T    "   T         F  ( not used here               )
C****  " 3    "    T    "   F         F  (subsequent calls w.out update)
C**** DYNAM is .true. when SOLVE3 called in time integration procedures.
C**** Element stiffness matrix ES & load vector ER for linear/nonlinear
C**** problems are generated by XPES/ZEES, respec.,  and assembled into
C**** the global stiffness matrix GK & load vector GR, and then into the
C**** constraint reduced matrices (non-zeros only) GKK and GRR.
C**** In dynamic problems (DYNAM =.true.) the element mass, damping and
C**** stiff. matrices EM,ED & ES are assembled into 1D global matrices
C**** GM,GD & GK which are then mapped into the solution matrix GKK.
C**** Note: These contain 2nd,1st & 0th order time-deriv coeffs,respec.
C**** The mapping coeffs NONY(noy,ny,nrc),CONY(noy,ny,nrc),
C**** noy=1,NONY(0,ny,nrc) are calc.d to reduce system of equations
C**** ny=1,NYNR(no_nynr,nrc,nc) to the system no=1,NOT(nrc,nr,nx) 
C**** by removing constraints, if UPDATE=.true.
C**** NONY(0,ny,nrc) is 0,1 when FIX(ny,1) is .true.,.false., resp.
C**** Note: In dynamic problems with constant time step GK is constant
C**** so F01BRF is used to factorize matrix on 1st time step only.
C**** KTYP22=0 : Static solution of initial accelerations.
C****   "    1 : first  order (Crank-Nicholson) time integration.
C****   "    2 : second order (Newmark; Wood-Zienkiewicz etc)
C****   "    3 : third  order (Houbolt; Hilbert-Hughes etc)
C**** YP(ny,1) on entry: contains b.c.s, defined by FIX(ny,1),
C****          (or current estimate of soln for nonlinear case)
C****          on exit : solution.
C****                    or rhs vector if IP is 5 or 6.
C**** YP(ny,5) is previous solution
C**** LGE(nhs,nrc) is location in global system of elem. var. nhs (=1,
C**** NHST(nrc))

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b08.cmn'
      INCLUDE 'cmiss$reference:b12.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
!     Parameter List
      INTEGER IDISP(10),IP,ISC_GKK(NISC_GKKM),ISC2_GKK(NISC2_GKKM),
     '  ISR_GKK(NISR_GKKM),ISR2_GKK(NISR2_GKKM),IWK1(5*NOM),
     '  IWK2(8*NOM),NONY(0:NOYM,NYM,NRCM),NPNY(0:6,NYM,0:NRCM),
     '  nr,nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNR(0:NY_R_M,0:NRCM,NCM),
     '  NZNY(NYM,*)
      REAL*8 CONY(0:NOYM,NYM,NRCM),CYNO(0:NYOM,NOOPM,NRCM),GD(NZ_GD_M),
     '  GK(NZ_GK_M),GKK(NZ_GKK_M),GM(NZ_GM_M),GR(NYROWM),
     '  GRR(NOM),WK1(4*NOM),XO(NOM),YP(NYM,NIYM)
      LOGICAL DYNAM1,DYNAM2,FIX(NYM,NIYFIXM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER GETNYR,IFAIL,
     '  no1,no2,no_nynr1,no_nynr2,
     '  noy1,noy2,ny1,ny2,ny3,nyo2,nz,nzz
      REAL*8 AA,BB,co1,co2,RESID1,RPMIN,SUM
      LOGICAL ABORT(4),FIRST,GROW,LBLOCK,LINEAR,LOGIC(3,3),
     '  UPDATE_MATRIX
      DATA LOGIC/3*.TRUE.,2*.TRUE.,.FALSE.,.TRUE.,2*.FALSE./

      CALL ENTERS('SOLVE3',*9999)

      IF(NOT(2,1,nr,nx).EQ.0) THEN
        ERROR=' >>The number of unknowns is zero'
        GOTO 9999
      ENDIF

      LINEAR=LOGIC(IP,1)
      FIRST =LOGIC(IP,2)
      UPDATE_MATRIX=LOGIC(IP,3)
      IF(IWRIT4(nr,nx).GT.0) THEN
        FORMAT='(/'' >SOLVE3   linear='',L1,'' first='',L1,'
     '    //''' update_matrix='',L1)'
        WRITE(OP_STRING,FORMAT) LINEAR,FIRST,UPDATE_MATRIX
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

C*** Initialise matrices

      DO no1=1,NOT(1,1,nr,nx)
        GRR(no1)=0.0d0
      ENDDO !no1

C*** Update RHS vector GR from flux boundary conditions

      DO no_nynr1=1,NYNR(0,0,2) !loop over global variables for nc=2
        ny1=NYNR(no_nynr1,0,2)  !is global (flux) variable number
        IF(FIX(ny1,1)) THEN  !flux is set as a b.c.
          ny2=GETNYR(2,NPNY,nr,1,0,ny1,NYNE,NYNP) !is global (flux) row number
          DO noy1=1,NONY(0,ny2,1)
            no1=NONY(noy1,ny2,1)
            co1=CONY(noy1,ny2,1)
            GRR(no1)=GRR(no1)+YP(ny1,1)*co1 !GRR value to applied 
          ENDDO !noy1                       !flux bc
        ENDIF !fix(ny1,1)
      ENDDO !no_nynr1

C***  Calculate constrained-reduced 1D global matrix GKK

      nzz=0
      NZZT(1,nr,nx)=0
      DO no_nynr1=1,NYNR(0,1,1) !loop over equations [rows] of GK
        ny1=NYNR(no_nynr1,1,1) !equation # 
        DO noy1=1,NONY(0,ny1,1)
          no1=NONY(noy1,ny1,1) !solution var # attached to ny1
          co1=CONY(noy1,ny1,1) !row coupling coefficient
          BB=GR(ny1)
          DO no_nynr2=1,NYNR(0,2,1) !loop over local variables of GK
            ny2=NYNR(no_nynr2,2,1) !local variable #
            ny3=GETNYR(1,NPNY,nr,0,2,ny2,NYNE,NYNP)
            !global variable # for local var ny2
            nz=NZNY(ny1,ny2)
            IF(nz.GT.0) THEN
C cpb 9/12/94 Old iy storage scheme
C             IF(KTYP22.GE.1) BB=BB-GK(nz)*YP(ny2,10)
C             IF(KTYP22.GE.2) BB=BB-GD(nz)*YP(ny2,11)
C             IF(KTYP22.EQ.3) BB=BB-GM(nz)*YP(ny2,12)
              IF(KTYP22.GE.1) BB=BB-GK(nz)*YP(ny3,5) !iy=5 previous sol
            ENDIF
          ENDDO !no_nynr2
          GRR(no1)=GRR(no1)+BB*co1
          DO no_nynr2=1,NYNR(0,2,1) !loop over local variables of GK
            ny2=NYNR(no_nynr2,2,1) !local variable #
            ny3=GETNYR(1,NPNY,nr,0,2,ny2,NYNE,NYNP)
            !global variable # for local var ny2
            nz=NZNY(ny1,ny2)
            IF(nz.GT.0) THEN
              IF(.NOT.DYNAM1.AND..NOT.DYNAM2) THEN
                AA=GK(nz)
              ELSE
                IF(KTYP22.EQ.0) THEN
                  AA=GM(nz)
                ELSE IF(KTYP22.EQ.1) THEN
                  AA=GD(nz)+A1*GK(nz)
                ELSE IF(KTYP22.EQ.2) THEN
                  AA=GM(nz)+A1*GD(nz)+A2*GK(nz)
                ELSE IF(KTYP22.EQ.3) THEN
                  AA=A1*GM(nz)+A2*GD(nz)+A3*GK(nz)
                ENDIF
              ENDIF
              IF(FIX(ny3,1)) THEN !variable is set as a bc
                GRR(no1)=GRR(no1)-AA*co1*YP(ny3,1)
              ELSE
                DO noy2=1,NONY(0,ny3,2) !solution vars attached to ny2
                  no2=NONY(noy2,ny3,2) !solution var #
                  co2=CONY(noy2,ny3,2) !variable coupling coeff
                  nzz=nzz+1
                  NZZT(1,nr,nx)=NZZT(1,nr,nx)+1
                  IF(nzz.LE.NZ_GKK_M) THEN
                    GKK(nzz)=AA*co1*co2
                  ENDIF
                  IF(FIRST) THEN
                    IF(nzz.LE.NISR_GKKM) THEN
                      ISR_GKK(nzz) =no1
                    ENDIF
                    IF(nzz.LE.NISR2_GKKM) THEN
                      ISR2_GKK(nzz)=no1
                    ENDIF
                    IF(nzz.LE.NISC_GKKM) THEN
                      ISC_GKK(nzz) =no2
                    ENDIF
                    IF(nzz.LE.NISC2_GKKM) THEN
                      ISC2_GKK(nzz)=no2
                    ENDIF
                  ENDIF
                ENDDO
C!!!cpb 9/12/94 additive constant may not be correct. Needs to be
C!!!checked.
C                co3=CYNO(0,no2,2) !is the additive variable constant
C                GRR(no1)=GRR(no1)-AA*co3 !adjust R.H.S. vector
C                                         !for the add const.
              ENDIF !fix(ny3,1)
            ENDIF !nz>0
          ENDDO !no_nynr2
        ENDDO !noy1
      ENDDO !no_nynr1
      IF(NZZT(1,nr,nx).GT.NZ_GKK_M) THEN
        ERROR='>>Increase NZ_GKK_MX'
        GOTO 9999
      ENDIF
      IF(NZZT(1,nr,nx).GT.NISR_GKKM) THEN
        ERROR='>>Increase NISR_GKKMX'
        GOTO 9999
      ENDIF
      IF(NZZT(1,nr,nx).GT.NISR2_GKKM) THEN
        ERROR='>>Increase NISR2_GKKMX'
        GOTO 9999
      ENDIF
      IF(NZZT(1,nr,nx).GT.NISC_GKKM) THEN
        ERROR='>>Increase NISC_GKKMX'
        GOTO 9999
      ENDIF
      IF(NZZT(1,nr,nx).GT.NISC2_GKKM) THEN
        ERROR='>>Increase NISC2_GKKMX'
        GOTO 9999
      ENDIF

C      IF(KTYP4.NE.0) THEN !Output global matrices
C        CALL WRITE_SOL_MATRIX(nr,nx,GKK,GRR,ERROR,*9999)
C      ENDIF

      IF(IWRIT4(nr,nx).GE.3) THEN
        WRITE(OP_STRING,
     '    '(/'' Global rhs vector GRR & stiffness matrix GKK:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NOT(1,1,nr,nx)='',I5,'
     '    //''', NOT(2,1,nr,nx)='','
     '    //'I5,'', NZZT(1,nr,nx)='',I6)') NOT(1,1,nr,nx),
     '    NOT(2,1,nr,nx),NZZT(1,nr,nx)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '(/'' GRR: '',8D12.4,/:(6X,8D12.4))') (GRR(no1),no1=1,
     '    NOT(1,1,nr,nx))
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' ISR: '',10(I5,X),/:(6X,10(I5,X)))') (ISR_GKK(nz),nz=1,
     '    NZZT(1,nr,nx))
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' ISC: '',10(I5,X),/:(6X,10(I5,X)))') (ISC_GKK(nz),nz=1,
     '    NZZT(1,nr,nx))
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' GKK: '',8D11.3,/:(6X,8D11.3))') (GKK(nz),nz=1,
     '    NZZT(1,nr,nx))
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      IF(FIRST) THEN
        IFAIL=111
        CALL F01BRF(NOT(1,1,nr,nx),NZZT(1,nr,nx),GKK,NZ_GKK_M,ISR_GKK,
     '    NZ_GKK_M,ISC_GKK,0.3d0,IWK1,IWK2,WK1,LBLOCK,GROW,ABORT,
     '    IDISP,IFAIL)
        IF(IWRIT4(nr,nx).GE.2) THEN
          WRITE(OP_STRING,'(/'' F01BRF: WK1(1)='',D12.4)') WK1(1)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(IFAIL.NE.0) THEN
          WRITE(OP_STRING,'('' IFAIL='',I3,'' in F01BRF'')') IFAIL
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          GOTO 9999
        ENDIF
        FIRST=.FALSE.
      ELSE IF(.NOT.FIRST) THEN
        IFAIL=111
        CALL F01BSF(NOT(1,1,nr,nx),NZZT(1,nr,nx),GKK,NZ_GKK_M,ISR2_GKK,
     '    ISC2_GKK,ISC_GKK,IWK1,IWK2,WK1,GROW,1.0d-4,RPMIN,
     '    ABORT(1),IDISP,IFAIL)
        IF(IWRIT4(nr,nx).GE.2) THEN
          WRITE(OP_STRING,'(/'' F01BSF: WK1(1)='',D12.4,'' RPMIN='','
     '	    //'D12.4)')WK1(1),RPMIN
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(IFAIL.NE.0) THEN
          WRITE(OP_STRING,'('' IFAIL='',I3,'' in F01BSF'')') IFAIL
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          GOTO 9999
        ENDIF
      ENDIF
      CALL F04AXF(NOT(1,1,nr,nx),GKK,NZ_GKK_M,ISC_GKK,IWK1,GRR,WK1,1,
     '  IDISP,RESID1)
      IF(IFAIL.NE.0) THEN
        WRITE(OP_STRING,'('' IFAIL='',I3,'' in F04AXF'')') IFAIL
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
        GOTO 9999
      ENDIF

C*** Put solution values into XO (for completness as XO should contain
C*** the solution vector)

      DO no1=1,NOT(1,1,nr,nx)
        XO(no1)=GRR(no1)
      ENDDO

      IF(IWRIT4(nr,nx).GE.2) THEN
        FORMAT='(/'' Solution values:'')'
        WRITE(OP_STRING,FORMAT)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        FORMAT='('' XO('',I5,'')= '',D13.6)'
        DO no1=1,NOT(1,1,nr,nx)
          WRITE(OP_STRING,FORMAT) no1,XO(no1)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

C*** Put solution values into YP

      DO no2=1,NOT(2,1,nr,nx)
        co1=CYNO(0,no2,2)
        DO nyo2=1,NYNO(0,no2,2)
          ny2=NYNO(nyo2,no2,2)
          co2=CYNO(nyo2,no2,2)
C!!!cpb 9/12/94 Other KTYP22's will need to be added. See backup
C!!! version of march1
          IF(KTYP22.EQ.1) THEN
            YP(ny2,1)=TINCR*(co2*XO(no2)+co1)+YP(ny2,5)
          ENDIF
        ENDDO !nyo2
      ENDDO !no2

C!!! What is this piece of code calculating ???
C!!! It was calculating the "residuals" i.e. GK*incremental soln - GR ??
      DO no_nynr1=1,NYNR(0,1,1) !Loop over the equations/rows of GK 
        ny1=NYNR(no_nynr1,1,1) !equation/row #
        SUM=0.0d0
        DO no_nynr2=1,NYNR(0,2,1) !Loop over the local columns of GK
          ny2=NYNR(no_nynr2,2,1) !local variable/column #
          ny3=NYNR(no_nynr2,0,1) !global variable/column # for local
                                 !variable ny2
          nz=NZNY(ny1,ny2)
          SUM=SUM+GK(nz)*YP(ny3,1)
        ENDDO
        YP(ny1,4)=SUM-GR(ny1) !iy=4 stores the residuals
      ENDDO

      CALL EXITS('SOLVE3')
      RETURN
 9999 CALL ERRORS('SOLVE3',ERROR)
      CALL EXITS('SOLVE3')
      RETURN 1
      END


      SUBROUTINE SOLVE5(IP,IBT,ICN,ICN2,IDO,INP,IRN,IRN2,IWK1,IWK2,LGE,
     '  NBH,NBJ,NEELEM,
     '  NHE,NHP,NJE,NKE,NKH,NONY,NP1OPT,NPE,NPF,NPNODE,
     '  NQE,nr,NRE,NVNE,NVNP,NW,nx,NYNE,NYNP,NZNY,
     '  CE,CG,CONY,CP,ED,EM,ER,ES,GK,GR,
     '  PAOPTI,PBOPTI,PG,RG,SE,VE,WG,
     '  XA,XE,XG,XP,YG,YP,ZA,ZE,ZG,ZP,
     '  GKK,GRR,WK1,WK2,FIX,SYMM,ERROR,*)

C**** Added 21-6-91.
C**** Finds solution of a system of SPARSE symmetric linear equations 
C**** no=1,NOT(nrc,nr,nx) for static problems only.
C**** Only the nonzeroes of the global matrix are stored, using the array
C**** NZNY(ny,ny) to indicate which entries are nonzero (see solve3).
C**** It is possible to just have 2 1d arrays rather than NZNY, but this
C**** is not yet implemented.
C**** IRN(nz) is row index of the non-zero stored in GKK(nz),nz=1,NZT.
C**** ICN(nz) is col   "      "      "        "       "        "
C**** YP(ny,1) contains essential b.c.s, defined by FIX(ny,1), on entry
C****                   incremental solution, on exit.
C****                or rhs vector if IP is 5 or 6.
C**** YP(ny,2)    "     nodal flux or stress bc.s defined by FIX(ny,2).
C**** YP(ny,3)    "     incremental b.c.s             "      FIX(ny,3).
C**** LGE(nhs,nrc) is location in global system of elem. var nhs (=1,
C****   NHST(nrc))
C
C Old comments form solve2 (from which solve5 was created).
C**** The mapping coeffs NONY(noy,ny,nrc),CONY(noy,ny,nrc),
C**** noy=1,NONY(0,ny,nrc)
C**** are calc.d to reduce system of equations ny=1,NYT(nrc,1,nx) to the
C**** system no=1,NOT(nrc,nr,nx) by removing constraints, if UPDATE=.true.
C**** NONY(0,ny,nrc) is 0,1 when FIX(ny,1) is .true.,.false., resp.
C**** IRN(nz) is row index of the non-zero stored in GKK(nz),nz=1,NZT.
C**** ICN(nz) is col   "      "      "        "       "        "
C**** FNY(ny1,ny2) is .true. if corresponding global value is non-zero.
C**** IP controls LINEAR/FIRST/UPDATE parameters as follows:
C**** IP=1 : linear=T; first=T; update=T  (linear eqtns with update)
C****  " 2     "    T    "   T         F  (  "      "   w.out  "   )
C****  " 3     "    F    "   T         T  (nonlin eqtns)
C**** SYMM indicates matrix symmetry. IWRIT4 controls diagnostic output.
C**** DYNAM is .true. when SOLVE5 called in time integration procedures.
C**** Element stiffness matrix ES & load vector ER for linear/nonlinear
C**** problems are generated by XPES/ZEES, respec.,  and assembled into
C**** the global stiffness matrix GKK & load vector GRR.
C**** KTYP22=0 : Static solution of initial accelerations.
C****   "    1 : first  order (Crank-Nicholson) time integration.
C****   "    2 : second order (Newmark; Wood-Zienkiewicz etc)
C****   "    3 : third  order (Houbolt; Hilbert-Hughes etc)
C**** YP(ny,1) contains essential b.c.s, defined by FIX(ny,1), on entry
C****                   solution (increm. if non.l. or dynam), on exit.
C****                or rhs vector if IP is 5 or 6.
C**** YP(ny,2)    "     nodal flux or stress bc.s defined by FIX(ny,2).
C**** YP(ny,3)    "     incremental b.c.s             "      FIX(ny,3).
C**** YP(ny,5)    "     reaction forces for linear solution, on exit.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b08.cmn'
      INCLUDE 'cmiss$reference:b10.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp100.cmn'
      INCLUDE 'cmiss$reference:opti00.cmn'
!     Parameter List
      INTEGER IBT(2,NIM,*),ICN(*),ICN2(*),IDO(NKM,0:NIM,*),
     '  INP(NNM,NIM,*),IP,
     '  IRN(*),IRN2(*),IWK1(*),IWK2(*),LGE(NHM*NSM,*),
     '  NBH(NHM,NCM,*),NBJ(NJM,*),
     '  NEELEM(0:NEM,0:*),
     '  NHE(*),NHP(*),NJE(*),
     '  NKE(NKM,NNM,NBFM,*),NKH(NHM,NCM,*),NONY(0:NOYM,NYM,*),
     '  NP1OPT(*),NPE(NNM,NBFM,*),NPF(12,*),
     '  NPNODE(0:NPM,0:*),NQE(NSM,NBFM,*),nr,
     '  NRE(*),NVNE(NNM,NBFM,NJM,*),NVNP(NJM,NPM,*),NW(NEM,*),
     '  nx,NYNE(NAM,NHM,NRCM,NCM,*),
     '  NYNP(NKM,NVM,NHM,NPM,NRCM,*),NZNY(NYM,*)
      REAL*8 CE(NMM,*),CG(NMM,*),CONY(0:NOYM,NYM,*),
     '  CP(NMM,*),ED(NVM,*),EM(NVM,*),ER(*),ES(NVM,*),
     '  GK(*),GR(*),PAOPTI(*),PBOPTI(*),
     '  PG(NSM,NUM,NGM,*),RG(*),
     '  SE(NSM,NBFM,*),VE(NSM,NKM,*),WG(NGM,*),
     '  XA(NAM,NJM,*),XE(NSM,*),XG(NJM,*),XP(NKM,NVM,NJM,*),
     '  YG(NGM,NJM,*),YP(NYM,*),ZA(NAM,NHM,NCM,*),ZE(NSM,*),
     '  ZG(NHM,*),ZP(NKM,NVM,NHM,NPM,*)
      REAL*8 GKK(*),GRR(*),WK1(*),WK2(*)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,*),SYMM
!     Local Variables
      INTEGER IFAIL,INFORM(4),nc,ne,nj,no,no1,no2,noelem,NOITS(2),
     '  nonode2,
     '  noopti,noy,noy1,noy2,NP1,NP2,NP3,nhs,nhs1,nhs2,NHST(2),
     '  nrc,nv,ny,ny1,ny2,nz,NZY,NZZ,NZZT
      REAL*8 AA,ACC(2),BB,CO,CO1,CO2,DROPTL,SF,SUM,XFIX,XMOVE
      CHARACTER FORMAT*500
      LOGICAL ABORT(3),ELEM,FIRST,LINEAR,LOGIC(6,3),UPDATE
      DATA LOGIC/2*.TRUE.,4*.FALSE.,3*.TRUE.,3*.FALSE.,.TRUE.,.FALSE.,
     '           2*.TRUE.,2*.FALSE./
      DATA ABORT/3*.FALSE./

      CALL ENTERS('SOLVE5',*9999)

      nc=1 ! temporary cpb 22/11/94
      nrc=2 ! temporary ajp 30/11/94
      nv=1 ! temporary cpb 22/11/94

      LINEAR=LOGIC(IP,1)
      FIRST =LOGIC(IP,2)
      UPDATE=LOGIC(IP,3)
      IF(IWRIT4.GT.0) THEN
        FORMAT='(/'' >SOLVE5   linear='',L1,'' first='',L1,'' symm='''
     '    //',L1,'' update='',L1)'
        WRITE(OP_STRING,FORMAT) LINEAR,FIRST,SYMM,UPDATE
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      IF(KTYP27.EQ.3) THEN !Coupled saturated-unsaturated flow.
C ***   Update nodal coordinates from optimisation params
C ***   As well as updating interface nodal coordinates we need
C ***   to ensure that elements in the FE region are not becoming too distorted
C ***   by updating the coordinates of all other nodes in each FE element
C ***   (except for the nodes around the fixed boundary).
C ***   NPRAD(I,noopti,nr), I=1,..,NPRAD(0,noopti,nr) are the nodes
C ***   on the same radial line as the current (noopti) free surface node
C ***   which are in the FE region and not on any fixed surface.
C ***   The node on the fixed surface is identified by
C ***   NPRAD(-1,noopti,nr).  This node number is needed in the mesh updating.
C ***   The array NPRAD is set up in IPOPTI

        IF(KTYP100.EQ.1) THEN !Flow from cavities
          DO noopti=1,NTOPTI
            NP1=NPJOIN(noopti,1) !Appropriate free surface node number
            NP3=NPRAD(-1,noopti,nr) !Appropriate fixed surface node number
            SF=PAOPTI(noopti)/PBOPTI(noopti)
            DO nj=1,NJT
              XFIX=XP(1,nv,nj,NP3)
              XMOVE=XP(1,nv,nj,NP1)
              IF(DABS(XFIX-XMOVE).GT.0.001D0) THEN
                !If XFIX-XMOVE=0 then no update necessary.
                DO nonode2=1,NPRAD(0,noopti,nr)
                  NP2=NPRAD(nonode2,noopti,nr) !Node numbers on same radial line
                  !Updating nodal coordinates on same radial line as current
                  !free surface node.  Need to ensure that this node remains
                  !in the same position relative to the free surface node
                  !and the nodes on the fixed surface.  Therefore need to
                  !use the coordinates of both the free surface node and the
                  !fixed surface node.
                  XP(1,nv,nj,NP2)=((SF-1.0d0)*XFIX*XMOVE+
     '                         (XFIX-SF*XMOVE)*XP(1,nv,nj,NP2))
     '                        /(XFIX-XMOVE)
                ENDDO !End of loop over radial line nodes.
              ENDIF
              !Updating free surface nodal coordinates
              XP(1,nv,nj,NP1)=SF*XMOVE
            ENDDO
            PBOPTI(noopti)=PAOPTI(noopti)
          ENDDO

        ELSE IF(KTYP100.EQ.2) THEN !Flow around cavities
          !AJP 26-7-91.
          !Need to increment Theta_sat.  Use the last noopti variable to
          !do this and impose symmetry on the XP coordinates for the last
          !noopti value.
          DO noopti=1,NTOPTI-1
            NP1=NPJOIN(noopti,1) !Appropriate free surface node number
            NP3=NPRAD(-1,noopti,nr) !Appropriate fixed surface node number
            SF=PAOPTI(noopti)/PBOPTI(noopti)
            DO nj=1,NJT
              XFIX=XP(1,nv,nj,NP3)
              XMOVE=XP(1,nv,nj,NP1)
              IF(DABS(XFIX-XMOVE).GT.0.001D0) THEN
                !If XFIX-XMOVE=0 then no update necessary.
                DO nonode2=1,NPRAD(0,noopti,nr)
                  NP2=NPRAD(nonode2,noopti,nr) !Node numbers on same radial line
                  !Updating nodal coordinates on same radial line as current
                  !free surface node.  Need to ensure that this node remains
                  !in the same position relative to the free surface node
                  !and the nodes on the fixed surface.  Therefore need to
                  !use the coordinates of both the free surface node and the
                  !fixed surface node.
                  XP(1,nv,nj,NP2)=((SF-1.0d0)*XFIX*XMOVE+
     '                         (XFIX-SF*XMOVE)*XP(1,nv,nj,NP2))
     '                        /(XFIX-XMOVE)
                ENDDO !End of loop over radial line nodes.
              ENDIF
              !Updating free surface nodal coordinates
              XP(1,nv,nj,NP1)=SF*XMOVE
            ENDDO
            PBOPTI(noopti)=PAOPTI(noopti)
            IF(ITYP10(1).EQ.1) THEN !Cartesians
              ny=NYNP(1,1,1,np3,nrc,nc)
              YP(ny,1)=XP(1,nv,NJT,NP3)-XP(1,nv,NJT,NP1)
            ENDIF
          ENDDO
          !Increment theta_sat
          THETA_SAT=PAOPTI(NTOPTI)/PBOPTI(NTOPTI)*THETA_SAT
          PBOPTI(NTOPTI)=PAOPTI(NTOPTI)
          !Impose symmetry on the appropriate XP variables
          NP1=NPJOIN(NTOPTI,1) !Appropriate free surface node number
          NP3=NPRAD(-1,NTOPTI,nr) !Appropriate fixed surface node number
          XP(1,nv,1,NP1)=-XP(1,nv,1,NPJOIN(1,1))
          XP(1,nv,1,NP3)=-XP(1,nv,1,NPRAD(-1,1,nr))
          DO nj=2,NJT
            XP(1,nv,nj,NP1)=XP(1,nv,nj,NPJOIN(1,1))
            XP(1,nv,nj,NP3)=XP(1,nv,nj,NPRAD(-1,1,nr))
          ENDDO
          DO nonode2=1,NPRAD(0,NTOPTI,nr)
            NP2=NPRAD(nonode2,NTOPTI,nr)
            XP(1,nv,1,NP2)=-XP(1,nv,1,NPRAD(nonode2,1,nr))
            DO nj=2,NJT
              XP(1,nv,nj,NP2)=XP(1,nv,nj,NPRAD(nonode2,1,nr))
            ENDDO
          ENDDO
          IF(ITYP10(1).EQ.1) THEN !Cartesians
            ny=NYNP(1,1,1,np3,nrc,nc)
            YP(ny,1)=XP(1,nv,NJT,NP3)-XP(1,nv,NJT,NP1)
          ENDIF
        ENDIF !Ktyp100 loop.
      ELSE IF(KTYP27.EQ.4) THEN
C ***   Update nodal coordinates from optimisation params
        DO noopti=1,NTOPTI
          NP1=NP1OPT(noopti)
          XP(1,nv,2,NP1)=PAOPTI(noopti)/PBOPTI(noopti)*XP(1,nv,2,NP1)
          PBOPTI(noopti)=PAOPTI(noopti)
        ENDDO
      ENDIF

      DO ny=1,NYT(1,1,nx)
        IF(FIX(ny,2)) THEN
          GR(ny)=YP(ny,2)
        ELSE
          GR(ny)=0.0D0
        ENDIF
      ENDDO

      IF(UPDATE) THEN           !Update global GK
        DO nz=1,NZ_GK_M           !initialize global arrays
          GK(nz)=0.0D0
        ENDDO
        IF(FIRST) THEN
          DO ny1=1,NYT(1,1,nx)
            DO ny2=1,NYT(2,1,nx)
              NZNY(ny1,ny2)=0
            ENDDO
          ENDDO
          no=0
          DO nrc=1,2
            DO ny=1,NYT(1,1,nx)
              IF(FIX(ny,1)) THEN
                NONY(0,ny,nrc)=0
              ELSE
                NONY(0,ny,nrc)=1
              ENDIF
              DO noy=1,NONY(0,ny,nrc)
                no=no+1
                NONY(noy,ny,nrc)=no
                CONY(noy,ny,nrc)=1.0D0
              ENDDO
            ENDDO
          ENDDO
          NOT(2,nr,nx)=no
          CALL ASSERT(NOT(2,nr,nx).GT.0,'>>Number of unknowns is zero',
     '      ERROR,*9999)
          CALL ASSERT(NOT(2,nr,nx).LE.NOM,
     '      '>>NOM too small (sb >NOT(nr,nx))',ERROR,*9999)
          CALL ASSERT(NYT(2,1,nx).GT.0,
     '      '>>no initial conditions for this region',ERROR,*9999)
        ENDIF

        nz=0
        DO noelem=1,NEELEM(0,nr) !is main element loop
          ne=NEELEM(noelem,nr)
          IF(NW(ne,1).GT.0) THEN

            CALL MELGE(LGE,NBH(1,1,ne),nc,NHE(ne),NHST,NKH,NPE(1,1,ne),
     '        NVNE(1,1,1,ne),NYNE(1,1,1,nc,ne),NYNP(1,1,1,1,1,nc),
     '        ERROR,*9999)

            IF(IWRIT4.GT.1) THEN                
              WRITE(OP_STRING,FORMAT) ne,NHST(1),(LGE(nhs,1),
     '          nhs=1,NHST(1))
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF

            ELEM=.FALSE.
            DO nhs1=1,NHST(1)   !Initialize element arrays
              ER(nhs1)=0.0D0
              DO nhs2=1,NHST(2)
                ES(nhs1,nhs2)=0.0D0
              ENDDO
            ENDDO

            IF(NW(ne,1).LE.20) THEN !finite element
              IF(ITYP1(nr).EQ.3) THEN !partial differential equation
                CALL XPXE(NBJ(1,ne),NJE(ne),NKE(1,1,1,ne),NPE(1,1,ne),
     '            NPF(1,1),NQE(1,1,ne),nr,NVNE(1,1,1,ne),NVNP,
     '            SE(1,1,ne),XA,XE,XP,ERROR,*9999)
                CALL XPES30(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),ne,
     '            NHE(ne),NJE(ne),NPE(1,1,ne),nr,
     '            CE(1,ne),CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '            VE(1,1,ne),WG,XE,XG,ZE,ZG,UPDATE,ERROR,*9999)
              ELSE IF(ITYP1(nr).EQ.4) THEN !linear elasticity
                CALL XPES40(NBH(1,1,ne),NBJ(1,ne),
     '            NHE(ne),NJE(ne),NKE(1,1,1,ne),NPE(1,1,ne),NPF,
     '            NQE(1,1,ne),nr,NVNE(1,1,1,ne),NVNP,NW(ne,1),
     '            CE(1,ne),CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '            VE(1,1,ne),WG,XA,XE,XG,XP,YG(1,1,ne),
     '            UPDATE,ERROR,*9999)
              ENDIF
            ELSE !boundary element
              ERROR=' >>Bdry elements should use SOLVE4'
              GO TO 9999
            ENDIF

            IF(IWRIT4.GE.3) THEN
              WRITE(OP_STRING,
     '          '(/'' Element load vector ER & stiffness matrix ES:'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              DO nhs1=1,NHST(1)
                WRITE(OP_STRING,'('' ER('',I2,'')='',E12.4,'
     '		  //''' ES: '',8E12.4/,(25X,8E12.4))') 
     '		  nhs1,ER(nhs1),(ES(nhs1,nhs2),nhs2=1,NHST(1))
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF
            DO nhs1=1,NHST(1)
              ny1=IABS(LGE(nhs1,1))
              GR(ny1)=GR(ny1)+ER(nhs1)
              DO nhs2=1,NHST(2)
                ny2=IABS(LGE(nhs2,2))
                IF(NZNY(ny1,ny2).EQ.0) THEN !increm & record pos.n of new dof
                  nz=nz+1
                  NZNY(ny1,ny2)=nz
                  NZY=nz
                ELSE                        !get existing 1D position
                  NZY=NZNY(ny1,ny2)
                ENDIF
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' nhs1='',I3,'' NY1='',I5,'
     '		    //''' nhs2='',I3,'' NY2='',I5,'' NZY='',I7)') 
     '		    nhs1,ny1,nhs2,ny2,NZY
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                GK(NZY)=GK(NZY)+ES(nhs1,nhs2)   !global stiffness matrix
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        NZT=nz
        IF(NZT.GT.NZ_GKK_M) THEN
          FORMAT='('' Increase array dimensions: NZT='',I6,'
     '      //''' NZ_GKK_M='',I6)'
          WRITE(OP_STRING,FORMAT) NZT,NZ_GKK_M
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

        IF(IWRIT4.GE.4) THEN
          WRITE(OP_STRING,'(/'' Global 1D arrays:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(/'' NOT(1,nr,nx)='',I4,'' NZT='',I6)') 
     '      NOT(1,nr,nx),NZT
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DO ny1=1,NYT(1,1,nx)
            WRITE(OP_STRING,'(/'' NZNY('',I2,'',ny2):''/,(4X,20I5))')
     '        ny1,(NZNY(ny1,ny2),ny2=1,NYT(2,1,nx))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDDO
          WRITE(OP_STRING,'(/'' GR: '',8E12.4/,(4X,8E12.4))') 
     '      (GR(ny),ny=1,NYT(1,1,nx))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(/'' GK: '',8E12.4/,(4X,8E12.4))') 
     '      (GK(nz),nz=1,NZT)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(.NOT.UPDATE) THEN
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(NW(ne,1).GT.0) THEN

            CALL MELGE(LGE,NBH(1,1,ne),nc,NHE(ne),NHST,NKH,NPE(1,1,ne),
     '        NVNE(1,1,1,ne),NYNE(1,1,1,nc,ne),NYNP(1,1,1,1,1,nc),
     '        ERROR,*9999)

            DO nhs1=1,NHST(1)        !initialize element arrays
              ER(nhs1)=0.0D0
            ENDDO

            IF(NW(ne,1).LE.20) THEN !finite element
              IF(ITYP1(nr).EQ.3) THEN !partial differential equation
                CALL XPXE(NBJ(1,ne),NJE(ne),NKE(1,1,1,ne),NPE(1,1,ne),
     '            NPF(1,1),NQE(1,1,ne),nr,NVNE(1,1,1,ne),NVNP,
     '            SE(1,1,ne),XA,XE,XP,ERROR,*9999)
                CALL YPZP(4,NBH,NEELEM,
     '            NHE,NHP,NKH,NPNODE,nr,NVNP,NYNE,NYNP,
     '            YP,ZA,ZP,ERROR,*9999)
                CALL ZPZE(1,NBH(1,1,ne),NHE(ne),NKE(1,1,1,ne),
     '            NPE(1,1,ne),NPF(1,1),nr,NVNE(1,1,1,ne),NVNP,NW(ne,1),
     '            SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '            ERROR,*9999)
                CALL XPES30(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),ne,
     '            NHE(ne),NJE(ne),NPE(1,1,ne),nr,
     '            CE(1,ne),CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '            VE(1,1,ne),WG,XE,XG,ZE,ZG,UPDATE,ERROR,*9999)
              ELSE IF(ITYP1(nr).EQ.4) THEN !linear elasticity
                CALL XPES40(NBH(1,1,ne),NBJ(1,ne),
     '            NHE(ne),NJE(ne),NKE(1,1,1,ne),NPE(1,1,ne),NPF,
     '            NQE(1,1,ne),nr,NVNE(1,1,1,ne),NVNP,NW(ne,1),
     '            CE(1,ne),CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '            VE(1,1,ne),WG,XA,XE,XG,XP,YG(1,1,ne),
     '            UPDATE,ERROR,*9999)
              ENDIF
            ELSE !boundary element
              ERROR=' >>Bdry elements should use SOLVE4'
              GO TO 9999
            ENDIF

            DO nhs=1,NHST(1)
              ny=IABS(LGE(nhs,1))
              GR(ny)=GR(ny)+ER(nhs)
            ENDDO

          ENDIF
        ENDDO
      ENDIF

C **  Calculate 1D global constraint-reduced matrix GKK
C **  (upper triangle only since GKK is symmetric)
      NZZ=0
      DO no=1,NOT(1,nr,nx)
        GRR(no)=0.0D0
      ENDDO
      DO ny1=1,NYT(1,1,nx)
        DO noy1=1,NONY(0,ny1,1)
          no1=NONY(noy1,ny1,1)
          CO1=CONY(noy1,ny1,1)
          BB=GR(ny1)
          DO ny2=1,NYT(2,1,nx)
            nz=NZNY(ny1,ny2)
            IF(nz.GT.0) THEN
              IF(KTYP22.GE.1) BB=BB-GK(nz)*YP(ny2,10)
            ENDIF
          ENDDO
          GRR(no1)=GRR(no1)+BB*CO1
          DO ny2=1,NYT(2,1,nx)
            nz=NZNY(ny1,ny2)
            IF(nz.GT.0) THEN
              AA=GK(nz)
              IF(NONY(0,ny2,2).EQ.0) THEN
                GRR(no1)=GRR(no1)-AA*CO1*YP(ny2,1)
                IF(IWRIT4.GT.1)
     '            WRITE(OP_STRING,'('' ny1='',I3,'' ny2='',I3,'
     '              //''' GRR('',I2,'')='',E12.3,'' GK('',I5,'')='','
     '              //'E12.3)') ny1,ny2,no1,GRR(no1),nz,GK(nz)
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ELSE
                IF(ny2.GE.ny1) THEN !only need to loop over upper triangle
                  DO noy2=1,NONY(0,ny2,2)
                    no2=NONY(noy2,ny2,2)
                    CO2=CONY(noy2,ny2,2)
                    NZZ=NZZ+1
                    GKK(NZZ)=AA*CO1*CO2
                    IF(FIRST) THEN
                      IRN(NZZ) =no1
                      IRN2(NZZ)=no1
                      ICN(NZZ) =no2
                      ICN2(NZZ)=no2
                    ENDIF
                  ENDDO
                ENDIF !End of loop over upper triangle only
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NZZT=NZZ
      IF(2*NZZT.GT.NZ_GKK_M) THEN
        FORMAT='('' Increase array dimensions: NZZT='',I6,'
     '    //''' NZ_GKK_M='',I6)'
        WRITE(OP_STRING,FORMAT) NZZT,NZ_GKK_M
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        FORMAT=
     '  '('' (Note:For sparse symmetric solver, NZ_GKK_M sb 2*NZZT) '')'
        WRITE(OP_STRING,FORMAT)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      IF(IWRIT4.GE.4) THEN
        WRITE(OP_STRING,'(/'' NOT(1,nr,nx)='',I4,'' NZZT='',I6)') 
     '    NOT(1,nr,nx),NZZT
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '(/'' Global rhs vector GRR & stiffness matrix GKK:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '(/'' GRR: '',8D12.4/,(4X,8D12.4))') (GRR(no),no=1,
     '    NOT(1,nr,nx))
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '(/'' IRN: '',  20I5/,(6X,20I5))')   (IRN(nz),nz=1,NZZT)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '(/'' ICN: '',  20I5/,(6X,20I5))')   (ICN(nz),nz=1,NZZT)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '(/'' GKK: '',8D12.4/,(4X,8D12.4))') (GKK(nz),nz=1,NZZT)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      IFAIL=111
      DROPTL=0.1d0
      CALL F01MAF(NOT(1,nr,nx),NZZT,GKK,NZ_GKK_M,IRN,NZ_GKK_M,ICN,
     '  DROPTL,0.8d0,WK1,IWK1,IWK2,ABORT,INFORM,IFAIL)
      IF(ifail.ne.0) THEN
        WRITE(OP_STRING,'('' IFAIL='',I3,'' in F01MAF'')') IFAIL
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GO TO 9999
      ENDIF
      ACC(1)=1.0D-6
      NOITS(1)=100
      IFAIL=111
      CALL F04MAF(NOT(1,nr,nx),NZZT,GKK,NZ_GKK_M,IRN,NZ_GKK_M,ICN,GRR,
     '  ACC,NOITS,WK1,WK2,IWK1,INFORM,IFAIL)
      IF(ifail.ne.0) THEN
        WRITE(OP_STRING,'('' IFAIL='',I3,'' in F04MAF'')') IFAIL
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        GO TO 9999
      ENDIF

      IF(IWRIT4.GE.4) THEN
        FORMAT='(/'' Solution values:'')'
        WRITE(OP_STRING,FORMAT)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        FORMAT='('' GRR('',I4,'')= '',D12.6)'
        DO no=1,NOT(1,nr,nx)
          WRITE(OP_STRING,FORMAT) no,GRR(no)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      DO ny=1,NYT(1,1,nx)
        SUM=0.0D0
        DO noy=1,NONY(0,ny,1)
          no=NONY(noy,ny,1)
          CO=CONY(noy,ny,1)
          SUM=SUM+GRR(no)*CO
        ENDDO
        IF(.NOT.FIX(ny,1)) YP(ny,1)=SUM
      ENDDO

      DO ny1=1,NYT(1,1,nx)
        SUM=0.0D0
        DO ny2=1,NYT(2,1,nx)
          nz=NZNY(ny1,ny2)
          SUM=SUM+GK(nz)*YP(ny2,1)
        ENDDO
        YP(ny1,5)=SUM-GR(ny1)
      ENDDO
C
      CALL EXITS('SOLVE5')
      RETURN
 9999 CALL ERRORS('SOLVE5',ERROR)
      CALL EXITS('SOLVE5')
      RETURN 1
      END


c cpb 25/10/94 This routine is not called anywhere and should be
C archived ?

      SUBROUTINE SOLVE6(IP,IBT,IDO,INP,LGE,
     '  NBH,NBJ,NEELEM,NHE,NJE,NKE,NKH,NONY,NPF,NP_INTERFACE,NPNE,
     '  NPNY,NQE,nr,NRE,NVHE,NVJE,NVJP,NW,nx,NXI,NYNE,NYNP,NYNR,
     '  CE,CG,CONY,CP,ED,EM,ER,ES,GD,GK,GM,GR,OMEGA,PG,
     '  RG,SE,VE,WG,XA,XE,XG,XP,YG,YP,ZE,ZG,
     '  WK1,FIX,SYMM,GKC,GKKC,GRRC,XOC,ERROR,*)

C**** Solves complex eqtns for freq domain analysis
C**** IP controls LINEAR/FIRST/UPDATE parameters as follows:
C**** IP=1 : linear=T; first=T; update=T  (linear eqtns with update)
C****  " 2     "    T    "   T         F  (  "      "   w.out  "   )
C****  " 3     "    F    "   T         T  (nonlin eqtns)
C**** SYMM indicates matrix symmetry. IWRIT4 controls diagnostic output.
C**** Element stiffness matrix ES & load vector ER for linear/nonlinear
C**** problems are generated by XPES/ZEES, respec.,  and assembled into
C**** the global stiffness matrix GKK & load vector GRR.
C**** The element mass, damping and
C**** stiffness matrices EM,ED & ES are assembled into global matrices
C**** GM,GD & GK which are then mapped into the complex matrix GKKC.
C**** The mapping coeffs NONY(noy,ny,nrc),CONY(noy,ny,nrc),
C**** noy=1,NONY(0,ny,nrc)
C**** are calc.d to reduce system of equations ny=1,NYNR(no_nynr,nrc,nc)
C**** to the system no=1,NOT(nrc,nr,nx) by removing constraints, 
C**** if UPDATE=.true.
C**** NONY(0,ny,nrc) is 0,1 when FIX(ny,1) is .true.,.false., resp.
C**** YP(ny,1) contains essential b.c.s, defined by FIX(ny,1), on entry
C****                   solution (increm. if non.l. or dynam), on exit.
C****                or rhs vector if IP is 5 or 6.
C**** YP(ny,2)    "     nodal flux or stress bc.s defined by FIX(ny,2).
C**** YP(ny,3)    "     incremental b.c.s             "      FIX(ny,3).
C**** YP(ny,5)    "     reaction forces for linear solution, on exit.
C**** LGE(nhs,nrc) is location in global system of elem. var. nhs (=1,
C****   NHST(nrc))

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:acti00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b08.cmn'
      INCLUDE 'cmiss$reference:b10.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:iwrit00.cmn'
      INCLUDE 'cmiss$reference:ktyp50.cmn'
!     Parameter List
      INTEGER IBT(2,NIM,*),IDO(NKM,0:NIM,*),INP(NNM,NIM,*),IP,
     '  LGE(NHM*NSM,*),NBH(NHM,NCM,*),NBJ(NJM,*),
     '  NEELEM(0:NEM,0:*),
     '  NHE(*),NJE(*),NKE(NKM,NNM,NBFM,*),NKH(NHM,NPM,*),
     '  NONY(0:NOYM,NYM,*),NPF(15,*),NP_INTERFACE(0:NPM,0:*),
     '  NPNE(NNM,NBFM,*),NPNY(0:5,NYM,0:*),
     '  NQE(NSM,NBFM,*),
     '  nr,NRE(*),
     '  NVHE(NNM,NBFM,NHM,*),NVJE(NNM,NBFM,NJM,*),NVJP(NJM,*),
     '  NW(NEM,*),nx,NXI(-NIM:NIM,0:*),NYNE(NAM,NHM,0:NRCM,NCM,*),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,*),NYNR(0:NYM,0:NRCM,*)
      REAL*8 CE(NMM,*),CG(NMM,*),
     '  CONY(0:NOYM,NYM,*),CP(NMM,*),
     '  ED(NHM*NSM,*),EM(NHM*NSM,*),ER(*),ES(NHM*NSM,*),
     '  GD(NYM,*),GK(NYM,*),GM(NYM,*),OMEGA,
     '  PG(NSM,NUM,NGM,*),RG(*),
     '  SE(NSM,NBFM,*),VE(NSM,NKM,*),WG(NGM,*),
     '  XA(NAM,NJM,*),XE(NSM,*),XG(NJM,*),XP(NKM,NVM,NJM,*),
     '  YG(NGM,NJM,*),YP(NYM,*),
     '  ZE(NSM,*),ZG(NHM,*)
      REAL*8 GR(*),WK1(*)
      COMPLEX*16 GKC(NCYM,*)
      COMPLEX*16 GKKC(NCYM,*),GRRC(*),XOC(*)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,*),SYMM
!     Local Variables
      INTEGER ID,IFAIL,nc,ne,NE1,no,no1,no2,noelem,
     '  no_nynr,no_nynr1,no_nynr2,
     '  noy,noy1,noy2,nhs1,nhs2,NHST(2),nrc,ny,ny1,ny2
      REAL*8 AA,BB,CO,CO1,CO2,DD,DETI,DETR,SUM,XOC_ABS,XOC_ANG
      COMPLEX*16 EKC(8,8)
      CHARACTER CFROMI*1,CHAR1*1,FORMAT*100
      LOGICAL FIRST,INLIST,LINEAR,LOGIC(6,3),UPDATE
      DATA LOGIC/2*.TRUE.,4*.FALSE.,3*.TRUE.,3*.FALSE.,.TRUE.,.FALSE.,
     '           2*.TRUE.,2*.FALSE./

      CALL ENTERS('SOLVE6',*9999)

      ERROR='>> RECODE SOLVE6'
      IF(ERROR.NE.' ') GOTO 9999

      nc=1 !temporary cpb 1/12/94

      CALL ASSERT(NOM*NOT(1,nr,nx).LE.NZM,'>>NZM too small',ERROR,*9999)
      LINEAR=LOGIC(IP,1)
      FIRST =LOGIC(IP,2)
      UPDATE=LOGIC(IP,3)
      IF(IWRIT4.GT.0) THEN
        FORMAT='(/'' >SOLVE6   linear='',L1,'' first='',L1,'' symm='''
     '    //',L1,'' update='',L1)'
        WRITE(OP_STRING,FORMAT) LINEAR,FIRST,SYMM,UPDATE
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      no=0
      DO nrc=1,2
        DO no_nynr=1,NYNR(0,1,1)
          ny=NYNR(no_nynr,1,1)
          IF(FIX(ny,1)) THEN
            NONY(0,ny,nrc)=0
          ELSE
            NONY(0,ny,nrc)=1
          ENDIF
          DO noy=1,NONY(0,ny,nrc)
            no=no+1
            NONY(noy,ny,nrc)=no
            CONY(noy,ny,nrc)=1.0D0
          ENDDO
        ENDDO
      ENDDO
      NOT(2,nr,nx)=no
      CALL ASSERT(NOT(2,nr,nx).GT.0,'>>The number of unknowns is zero',
     '  ERROR,*9999)
      CALL ASSERT(NOT(2,nr,nx).LE.NCYM,
     '  '>>NCYM too small (sb >NOT(2,nr,nx))',ERROR,*9999)
      DO no1=1,NOT(1,nr,nx)
        GRRC(no1)=DCMPLX(0.0d0,0.0d0)
        DO no2=1,NOT(2,nr,nx)
          GKKC(no1,no2)=DCMPLX(0.0D0,0.0D0)
        ENDDO
      ENDDO

      IF(UPDATE) THEN
        DO no_nynr=1,NYNR(0,1,1)
          ny=NYNR(no_nynr,1,1)
          IF(FIX(ny,2)) THEN
            GR(ny)=YP(ny,2)
          ELSE
            GR(ny)=0.0D0
          ENDIF
        ENDDO

        DO no_nynr1=1,NYNR(0,1,1)
          ny1=NYNR(no_nynr1,1,1)
          DO no_nynr2=1,NYNR(0,2,1)
            ny2=NYNR(no_nynr2,2,1)
            GK(ny1,ny2)=0.0D0
            GM(ny1,ny2)=0.0D0
            GD(ny1,ny2)=0.0D0
          ENDDO
        ENDDO

        DO noelem=1,NEELEM(0,nr) !is main element loop
          ne=NEELEM(noelem,nr)
          IF(NW(ne,1).GT.0) THEN

            CALL MELGE(LGE,NBH(1,1,ne),nc,ne,NHE(ne),NHST,NKH,
     '        NP_INTERFACE,NPNE(1,1,ne),NPNY,nr,NVHE(1,1,1,ne),
     '        NW(1,1),NXI,NYNE,NYNP,NYNR,ERROR,*9999)

            DO nhs1=1,NHST(1)    !Initialize element arrays
              ER(nhs1)=0.0d0
              DO nhs2=1,NHST(2)
                ES(nhs1,nhs2)=0.0d0
                EM(nhs1,nhs2)=0.0d0
                ED(nhs1,nhs2)=0.0d0
              ENDDO
            ENDDO

            IF(NW(ne,1).LE.20) THEN !finite element
              IF(ITYP1(nr).EQ.3) THEN !partial differential equation
                CALL XPXE(NBJ(1,ne),NJE(ne),NKE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),NQE(1,1,ne),nr,NVJE(1,1,1,ne),NVJP,
     '            SE(1,1,ne),XA,XE,XP,ERROR,*9999)
                CALL XPES30(IBT,IDO,INP,NBH(1,1,ne),NBJ(1,ne),ne,
     '            NHE(ne),NJE(ne),NPNE(1,1,ne),nr,
     '            CE(1,ne),CG,CP,ED,EM,ER,ES,PG,RG,SE(1,1,ne),
     '            VE(1,1,ne),WG,XE,XG,ZE,ZG,UPDATE,ERROR,*9999)
              ELSE IF(ITYP1(nr).EQ.4) THEN !linear elasticity
                CALL XPES40(NBH(1,1,ne),NBJ(1,ne),
     '            NHE(ne),NJE(ne),NKE(1,1,1,ne),
     '            NPF,NPNE(1,1,ne),NQE(1,1,ne),nr,NVJE(1,1,1,ne),
     '            NVJP,NW(ne,1),
     '            CE(1,ne),CG,CP,ED,EM,ER,ES,
     '            PG,RG,SE(1,1,ne),VE(1,1,ne),WG,
     '            XA,XE,XG,XP,YG(1,1,ne),
     '            UPDATE,ERROR,*9999)
              ENDIF
            ELSE !boundary element
              ERROR=' >>Bdry elements should use SOLVE4'
              GO TO 9999
            ENDIF

            IF(IWRIT4.GE.3) THEN
              WRITE(OP_STRING,
     '          '(/'' Element load vector ER & stiffness matrix ES:'
     '          //' (Element '',I5,'')'')') ne
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              DO nhs1=1,NHST(1)
                WRITE(OP_STRING,'('' ER('',I2,'')='',E12.4,'
     '		  //''' ES: '',8E12.4/,(25X,8E12.4))') 
     '		  nhs1,ER(nhs1),(ES(nhs1,nhs2),nhs2=1,NHST(2))
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDDO
              WRITE(OP_STRING,'(/'' Element matrices EM & ED:'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              DO nhs1=1,NHST(1)
                WRITE(OP_STRING,'('' nhs1='',I2,'' EM: '',8E12.4/,'
     '            //'(12X,8E12.4))') nhs1,(EM(nhs1,nhs2),nhs2=1,NHST(2))
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDDO
              DO nhs1=1,NHST(1)
                WRITE(OP_STRING,'('' nhs1='',I2,'' ED: '','
     '            //'8E12.4/,(12X,8E12.4))') 
     '            nhs1,(ED(nhs1,nhs2),nhs2=1,NHST(2))
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF

            DO nhs1=1,NHST(1)
              ny1=IABS(LGE(nhs1,1))
              GR(ny1)=GR(ny1)+ER(nhs1)
              DO nhs2=1,NHST(2)
                ny2=IABS(LGE(nhs2,2))
                GK(ny1,ny2)=GK(ny1,ny2)+ES(nhs1,nhs2)
                GM(ny1,ny2)=GM(ny1,ny2)+EM(nhs1,nhs2)
                GD(ny1,ny2)=GD(ny1,ny2)+ED(nhs1,nhs2)
              ENDDO
            ENDDO

          ENDIF
        ENDDO

        IF(IWRIT4.GE.4) THEN
          WRITE(OP_STRING,
     '      '(/'' Global load vector GR & stiffness matrix GK:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DO no_nynr1=1,NYNR(0,1,1)
            ny1=NYNR(no_nynr1,1,1)
            WRITE(OP_STRING,'('' GR('',I4,'')='',E12.4,'' GK: '','
     '	      //'8E12.4'//'/,(25X,8E12.4))') ny1,GR(ny1),
     '	      (GK(ny1,NYNR(no_nynr2,2,1)),no_nynr2=1,NYNR(0,2,1))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDDO
          WRITE(OP_STRING,'(/'' Global matrices GM & GD:'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DO no_nynr1=1,NYNR(0,1,1)
            ny1=NYNR(no_nynr1,1,1)
            WRITE(OP_STRING,'('' NY1='',I2,'' GM: '',8E12.4'
     '        //'/,(12X,8E12.4))') ny1,
     '	      (GM(ny1,NYNR(no_nynr2,2,1)),no_nynr2=1,NYNR(0,2,1))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDDO
          DO no_nynr1=1,NYNR(0,1,1)
            ny1=NYNR(no_nynr1,1,1)
            WRITE(OP_STRING,'('' NY1='',I2,'' GD: '',8E12.4'
     '        //'/,(12X,8E12.4))') ny1,
     '	      (GD(ny1,NYNR(no_nynr2,2,1)),no_nynr2=1,NYNR(0,2,1))
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF

      ENDIF !update

      IF(KTYP59.EQ.4) THEN !compute active cell terms
          DO no_nynr1=1,NYNR(0,1,1)
            ny1=NYNR(no_nynr1,1,1)
            DO no_nynr2=1,NYNR(0,2,1)
              ny2=NYNR(no_nynr2,2,1)
            GKC(ny1,ny2)=DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
        DO noelem=1,NEELEM(0,nr) !is main element loop
          ne=NEELEM(noelem,nr)
          IF(INLIST(ne,NE_ACTI(1),NE_ACTI(0),NE1)) THEN !elem is active
            CALL MELGE(LGE,NBH(1,1,ne),nc,ne,NHE(ne),NHST,NKH,
     '        NP_INTERFACE,NPNE(1,1,ne),NPNY,nr,NVHE(1,1,1,ne),
     '        NW(1,1),NXI,NYNE,NYNP,NYNR,ERROR,*9999)
            CALL ACTIVE(NBH(1,1,ne),EKC,OMEGA,ERROR,*9999)
            DO nhs1=1,NHST(1)
              ny1=IABS(LGE(nhs1,1))
              DO nhs2=1,NHST(2)
                ny2=IABS(LGE(nhs2,2))
                GKC(ny1,ny2)=GKC(ny1,ny2)-EKC(nhs1,nhs2)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

      DO no_nynr1=1,NYNR(0,1,1)
        ny1=NYNR(no_nynr1,1,1)
        DO noy1=1,NONY(0,ny1,1)
          no1=NONY(noy1,ny1,1)
          CO1=CONY(noy1,ny1,1)
          BB=GR(ny1)
          GRRC(no1)=GRRC(no1)+DCMPLX(BB*CO1)
          DO no_nynr2=1,NYNR(0,2,1)
            ny2=NYNR(no_nynr2,2,1)
            AA=-OMEGA**2*GM(ny1,ny2)+GK(ny1,ny2)
            DD= OMEGA*GD(ny1,ny2)
            IF(NONY(0,ny2,2).EQ.0) THEN
              GRRC(no1)=GRRC(no1)-DCMPLX(AA*CO1*YP(ny2,1))
            ELSE
              DO noy2=1,NONY(0,ny2,2)
                no2=NONY(noy2,ny2,2)
                CO2=CONY(noy2,ny2,2)
                GKKC(no1,no2)=GKKC(no1,no2)+DCMPLX(AA*CO1*CO2,
     '            DD*CO1*CO2)
                IF(KTYP59.EQ.4) THEN !add active cell terms
                  IF(DABS(CO1*CO2).GT.1.0D-6) THEN
                    GKKC(no1,no2)=GKKC(no1,no2)+DCMPLX(GKC(ny1,ny2))
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      IF(IWRIT4.GE.4) THEN
        WRITE(OP_STRING,
     '    '(/'' Global load vector GRRC & stiffness matrix GKKC:'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NOT(1,nr,nx)='',I4)') NOT(1,nr,nx)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        DO no1=1,NOT(1,nr,nx)
          WRITE(OP_STRING,'(/'' GRRC('',I4,'')='',2D12.4,'
     '	    //''' GKKC: '',8D12.4,(/25X,8D12.4))') 
     '	    no1,GRRC(no1),(GKKC(no1,no2),no2=1,NOT(2,nr,nx))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      IFAIL=1
      CALL F03AHF(NOT(1,nr,nx),GKKC,NCYM,DETR,DETI,ID,WK1,IFAIL)
      IF(ifail.ne.0) THEN
        CHAR1=CFROMI(IFAIL,'(I1)')
        ERROR=' IFAIL='//CHAR1//' in F01BNF'
        GOTO 9999
      ENDIF
      CALL F04AKF(NOT(1,nr,nx),1,GKKC,NCYM,WK1,GRRC,NCYM)
      DO no=1,NOT(1,nr,nx)
        XOC(no)=GRRC(no)
      ENDDO

      IF(IWRIT4.GE.1) THEN
        FORMAT='(/'' Solution values:'')'
        WRITE(OP_STRING,FORMAT)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        FORMAT='('' XOC('',I4,'')= '',2D12.4,'
     '    //''' abs(XOC)='',D12.4,'' ang(XOC)='',D12.4)'
        DO no=1,NOT(1,nr,nx)
          WRITE(OP_STRING,FORMAT) no,XOC(no),XOC_ABS,XOC_ANG
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

C *** Write solutions to the YP array
      DO no_nynr=1,NYNR(0,1,1)
        ny=NYNR(no_nynr,1,1)
        SUM=0.0D0
        DO noy=1,NONY(0,ny,2)
          no=NONY(noy,ny,2)
          CO=CONY(noy,ny,2)
          SUM=SUM+XOC(no)*CO
        ENDDO
        IF(.NOT.FIX(ny,1)) YP(ny,1)=SUM
      ENDDO

C **  Needs attention!!!!!!
      DO no_nynr1=1,NYNR(0,1,1)
        ny1=NYNR(no_nynr1,1,1)
        SUM=0.0D0
        DO no_nynr2=1,NYNR(0,2,1)
          ny2=NYNR(no_nynr2,2,1)
          SUM=SUM+GK(ny1,ny2)*YP(ny2,1)
        ENDDO
        YP(ny1,5)=SUM-GR(ny1)
      ENDDO

      CALL EXITS('SOLVE6')
      RETURN
 9999 CALL ERRORS('SOLVE6',ERROR)
      CALL EXITS('SOLVE6')
      RETURN 1
      END
 

Module FE09
===========

      REAL*8 FUNCTION PYTHAG(A,B)

C#### Function: PYTHAG
C###  Description:
C###    PYTHAG finds DSQRT(A**2+B**2) without overflow or destructive
C###    underflow

      IMPLICIT NONE
!     Parameter List
      REAL*8 A,B
!     Local Variables
      REAL*8 P,R,S,T,U
C
      P = DMAX1(DABS(A),DABS(B))
      IF(P.NE.0.0d0) THEN
        R = (DMIN1(DABS(A),DABS(B))/P)**2        
 10     T = 4.0d0 + R
        IF(T.EQ.4.0d0) GOTO 20
        S = R/T
        U = 1.0d0 + 2.0d0*S
        P = U*P
        R = (S/U)**2*R
        GOTO 10
      ENDIF
 20   PYTHAG=P
C
      RETURN
      END


      REAL*8 FUNCTION TIMER(IHANDLE,T_FLAG)

C#### Function: TIMER
C###  Type: REAL*8
C###  Description: 
C###    TIMER returns elapsed time (if T_FLAG=T_ELAPSED) or 
C###    CPU time (if T_FLAG=T_CPU) in seconds since initialisation.
C###    To initialises timer call TIMER with T_FLAG=T_INITIALISE.
C###    NOTE: ETIME uses single precision real numbers.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:time02.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER IHANDLE,T_FLAG
!     Local Variables
      INTEGER CERROR(50),CERRLEN,ERR
      REAL*8 RETURN_TIME
      CHARACTER ERROR*255

      ERR=0
      CALL CTIMER(RETURN_TIME,IHANDLE,T_FLAG,ERR,CERROR)
      TIMER = RETURN_TIME
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,OP_STRING(1))
        OP_STRING(1)(CERRLEN:)=' '
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF

 9999 RETURN
      END




