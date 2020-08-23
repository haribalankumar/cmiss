      SUBROUTINE DEWIND(NPNODE,XP,ZD,STRING,ERROR,*)

C#### Subroutine: DEWIND
C###  Description:
C###    DEWIND sets world coordinate dimensions.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'disp00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'post00.cmn'
!     Parameter List
      INTEGER NPNODE(0:NP_R_M,0:NRM)
      REAL*8 XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE,iw,IWK(6),
     '  N3CO,nd,nj,noiw,nonode,np,nr,NTIW,nv
c     INTEGER INSTAT,noch
      REAL*8 SMAX,X(3),Z(3),ZZMAX(3),ZZMIN(3)
      CHARACTER CHAR1*1,FILE*(MXCH),STATUS*3
      LOGICAL ABBREV,ALL_REGIONS,CALCU,CBBREV,FILIO,FIRST_TIME,
     '  GENER,MOUSE

      CALL ENTERS('DEWIND',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

        WRITE(CHAR1,'(I1)') 2*NJT-3

C#### Command: FEM define window<;(c/l/p/r/w)[c]><;FILENAME[$current]><;(PATH/example)[$current]>
C###  Parameter:      <on (WS#/xy/yz/xz)[xz]>
C###  Description:
C###    World coordinate window dimensions are read from or written to
C###    a file FILENAME.ipwind in the directory
C###    specified by PATH with $current specifing the current default file..
C###  Parameter:      <for (postscript/encapsulated)[encapsulated]>
C###    Specify whether the format of the window output is
C###    postscript or encapsulated postscript.
CC###  Parameter:      <with (colour/black_white)[black_white]>
CC###    Specify whether the window is colour or balck and white.
CC###  Parameter:      <in (portrait/landscape)[portrait]>
CC###    Specify whether orientation of the window is portrait or landscape.

        OP_STRING(1)=STRING(1:IEND)//'<;(c/l/p/r/w)[c]>'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<on (WS#/xy/yz/xz)[xy]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

        OP_STRING(1)=STRING(1:IEND)//' <for (postscript/encapsulated)>'
C        OP_STRING(2)=BLANK(1:15)
C     '    //'<with (colour/black_white)[black_white]>'
C        OP_STRING(3)=BLANK(1:15)
C     '    //'<in (portrait/landscape)[portrait]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEWIND',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set Use_GraphicsM=1',ERROR,*9999)
        nv=1 ! Temporary MPN 12-Nov-94
        IPFILE=1 !is input file version number on 24-Jan-1990
        CALL PARSE_QUALIFIERS(' CDLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(NTCOQU(noco).EQ.0) THEN
          CALCU=.TRUE.
        ELSE IF(CALCU) THEN !to define FILE for postscript file
          CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        ENDIF
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        EPS=.FALSE.
        POSTSCRIPT=.FALSE.
        IF(CBBREV(CO,'ON',1,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'XY',1,noco+1,NTCO,N3CO).OR
     '      .CBBREV(CO,'YX',1,noco+1,NTCO,N3CO)) THEN
            NTIW=1
            IWK(1)=1
          ELSE IF(CBBREV(CO,'YZ',1,noco+1,NTCO,N3CO).
     '         OR.CBBREV(CO,'ZY',1,noco+1,NTCO,N3CO)) THEN
            NTIW=1
            IWK(1)=2
          ELSE IF(CBBREV(CO,'XZ',1,noco+1,NTCO,N3CO).
     '         OR.CBBREV(CO,'ZX',1,noco+1,NTCO,N3CO)) THEN
            NTIW=1
            IWK(1)=3
          ELSE
            CALL PARSIL(CO(N3CO+1),6,NTIW,IWK,ERROR,*9999)
          ENDIF
        ELSE
          NTIW=1
          IWK(1)=1
        ENDIF

        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'ENCAPSULATED',1)) THEN
            EPS=.TRUE.
            POSTSCRIPT=.FALSE.
C            DEWIND_POSTSCRIPT=.TRUE.
          ELSE IF(ABBREV(CO(N3CO+1),'POSTSCRIPT',1)) THEN
            EPS=.FALSE.
            POSTSCRIPT=.TRUE.
C            DEWIND_POSTSCRIPT=.TRUE.
          ENDIF
        ENDIF

        IF(POSTSCRIPT.OR.EPS) THEN
C          IF(CBBREV(CO,'WITH',1,noco+1,NTCO,N3CO)) THEN
C            IF(ABBREV(CO(N3CO+1),'COLOUR',1)) THEN
C              COLOURPS=.TRUE.
C            ELSE
C              COLOURPS=.FALSE.
C            ENDIF
C          ELSE
C            COLOURPS=.FALSE.
C          ENDIF
C          IF(CBBREV(CO,'IN',1,noco+1,NTCO,N3CO)) THEN
C            IF(ABBREV(CO(N3CO+1),'LANDSCAPE',1)) THEN
C              PORTRAIT=.FALSE.
C            ELSE
C              PORTRAIT=.TRUE.
C            ENDIF
C          ELSE
C            PORTRAIT=.TRUE.
C          ENDIF
          IF(NJT.LE.2) THEN ! Use GKS postscript
            NTIW=1
            IWK(1)=15
          ELSE ! 3D so use PHIGS postscript
            NTIW=1
            IWK(1)=16
          ENDIF
C CPB 5/8/93 Adding in all the postscript workstation types. Note: PHIGS
C and GKS use the same workstation numbers for the different postscript
C formats.
C
C 21/11/96
C          IF(POSTSCRIPT) THEN
C            IF(COLOURPS) THEN
C              IF(PORTRAIT) THEN
C                PSTYPE=273678398 ! ie x1050003e - colour a4 portrait ps
C              ELSE
C                PSTYPE=5242942   ! ie x0050003e - colour a4 landscape ps
C              ENDIF
C            ELSE
C              IF(PORTRAIT) THEN
C                PSTYPE=273678397 ! ie x1050003d - b/w a4 portrait ps
C              ELSE
C                PSTYPE=5242941   ! ie x0050003d - b/w a4 landscape ps
C              ENDIF
C            ENDIF
C          ELSE
C            IF(COLOURPS) THEN
C              IF(PORTRAIT) THEN
C                PSTYPE=273678402 ! ie x10500042 - colour a4 portrait eps
C              ELSE
C                PSTYPE=5242946   ! ie x00500042 - colour a4 landscape eps
C              ENDIF
C            ELSE
C              IF(PORTRAIT) THEN
C                PSTYPE=273678401 ! ie x10500041 - b/w a4 portrait eps
C              ELSE
C                PSTYPE=5242945   ! ie x00500041 - b/w a4 landscape eps
C              ENDIF
C            ENDIF
C          ENDIF
        ENDIF

        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          ALL_REGIONS=.FALSE. !when parse_regions not called
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'wind',STATUS,
     '        ERR,ERROR,*9999)
            CALL IPWIND(ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO

        ELSE IF(CALCU) THEN
C cpb 7/3/95 Adding data points into screen calculation
          DO nj=1,NJT
            ZZMAX(nj)=-1.0d6
            ZZMIN(nj)= 1.0d6
          ENDDO
          IF(NDT.GT.0) THEN
            DO nd=1,NDT
              DO nj=1,NJT
                X(nj)=ZD(nj,nd)
              ENDDO !nj
              CALL XZ(ITYP10(1),X,Z)
              DO nj=1,NJT
                IF(Z(nj).LT.ZZMIN(nj)) ZZMIN(nj)=Z(nj)
                IF(Z(nj).GT.ZZMAX(nj)) ZZMAX(nj)=Z(nj)
              ENDDO !nj
            ENDDO !nd
          ENDIF
          IF(NPT(1).GT.0) THEN
            DO nr=1,NRT
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  X(nj)=XP(1,nv,nj,np)
                ENDDO
                CALL XZ(ITYP10(1),X,Z)
                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  IF(Z(nj).LT.ZZMIN(nj)) ZZMIN(nj)=Z(nj)
                  IF(Z(nj).GT.ZZMAX(nj)) ZZMAX(nj)=Z(nj)
                ENDDO !nj
              ENDDO !nonode
            ENDDO
          ELSE IF(NDT.EQ.0) THEN
            ZZMIN(1)=-1.0d0
            ZZMIN(2)=-1.0d0
            ZZMIN(3)=-1.0d0
            ZZMAX(1)= 1.0d0
            ZZMAX(2)= 1.0d0
            ZZMAX(3)= 1.0d0
          ENDIF

          IF(NJT.EQ.1) THEN !only 1 global coordinate defined
            ZZMIN(2)=-1.0d0
            ZZMAX(2)=1.0d0
          ENDIF

          SMAX=-1.0d6
          DO nj=1,NJT
            IF((ZZMAX(nj)-ZZMIN(nj)).GT.SMAX)
     '        SMAX=ZZMAX(nj)-ZZMIN(nj)
          ENDDO !nj

          XMIN=REAL((ZZMIN(1)+ZZMAX(1))/2.0d0-0.75d0*SMAX)
          YMIN=REAL((ZZMIN(2)+ZZMAX(2))/2.0d0-0.75d0*SMAX)
          XMAX=REAL((ZZMIN(1)+ZZMAX(1))/2.0d0+0.75d0*SMAX)
          YMAX=REAL((ZZMIN(2)+ZZMAX(2))/2.0d0+0.75d0*SMAX)
          IF(NJT.EQ.3) THEN
            ZMIN=REAL((ZZMIN(3)+ZZMAX(3))/2.0d0-0.75d0*SMAX)
            ZMAX=REAL((ZZMIN(3)+ZZMAX(3))/2.0d0+0.75d0*SMAX)
            DIAG=SQRT((XMAX-XMIN)**2+(YMAX-YMIN)**2+(ZMAX-ZMIN)**2)
          ELSE
            ZMIN=0.0
            ZMAX=0.0
            DIAG=SQRT((XMAX-XMIN)**2+(YMAX-YMIN)**2)
          ENDIF
          IF(DIAG.LT.1.0e-6) THEN !revert to standard window
            XMIN=-1.0
            XMAX= 1.0
            YMIN=-1.0
            YMAX= 1.0
            ZMIN=-1.0
            ZMAX= 1.0
          ENDIF

          SCREEN_WIDTH=49.0

        ENDIF

        IF(DOP) THEN
          WRITE(OP_STRING,*) ' Window: ',
     '      XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,*) ' Diag: ',DIAG
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,*) ' Screen_width= ',SCREEN_WIDTH
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        DO noiw=1,NTIW !to define workstation IWK(noiw)
          iw=IWK(noiw)
          IF(IWKS(iw).EQ.1) THEN !workstation is already defined
            CALL CLOSE_WS(iw,ERROR,*9999) !..so close first
            IWKS(iw)=0
          ENDIF
          CALL SETUP(iw,ERROR,*9999)

          IF(DOP) THEN                               !temporary
            WRITE(OP_STRING,*) ' in dewind:'
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*) ' iw,njt=',iw,NJT
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

          CALL ASSERT(WINDOW_TYPE(1:3).NE.'VWS','CODE MOVED TO ARCHIVE',
     '      ERROR,*9999)

        ENDDO
      ENDIF

      CALL EXITS('DEWIND')
      RETURN
 9999 CALL ERRORS('DEWIND',ERROR)
      CALL EXITS('DEWIND')
      RETURN 1
      END



