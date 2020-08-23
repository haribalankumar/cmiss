      SUBROUTINE DRPLIN(ISEG,ISPLIN,PAOPTI,WD,ZD,CSEG,STRING,ERROR,*)

C#### Subroutine: DRPLIN
C###  Description:
C###    DRPLIN draws polyline segment.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'draw00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'plin00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISPLIN(NWM,NGRSEGM)
      REAL*8 PAOPTI(*),WD(NJM,NDM),ZD(NJM,NDM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER NTLM
      PARAMETER(NTLM=20)
      INTEGER i,IBEG,IBEG1,IEND,IEND1,IFROMC,INDEX,INDEX_PLIN,INDEX_OLD,
     '  INDEX_POLYLINE,INSTAT,iw,IWK(6),j,n,N3CO,nd,nj,nl,
     '  NOPLIN,no_point,nopts,nolines,NO_SECTIONS,NTIW,NT_POINT,NTPTS
C DPN 17-12-97 - add new arrays to pass through to BEZIER
      INTEGER ISL2BE_LOC(NTLM),ISL3BE_LOC(NTLM),ISN2BE_LOC(NTLM),
     '  ISN3BE_LOC(NTLM)
      REAL*8 ALFA,ETA,PL(3,21),PTS(3,20),RFROMC,WEIGHT,XBEZ(4),XI,
     '  XPOINTS(20),XPTS(20),YBEZ(4),YPOINTS(20),YPTS(20)
      CHARACTER FILE*100,STATUS*3,TYPE*10,CHAR*2
      LOGICAL ABBREV,CALCU,CBBREV,FILIO,MOUSE,OPTIM,SEGME

      CALL ENTERS('DRPLIN',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM draw polyline;l/p/r/w<;FILENAME>[default]<;example>
C###  Description:
C###    Draws a polyline.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME>['//FILE00(IBEG1:IEND1)//']<;example>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw polyline;s
C###  Parameter:    <index ID#[1]>
C###  Parameter:    <on WS_ID#[1]>
C###    Specify the workstation (GX window) to draw the
C###    polyline on.
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Description:
C###    ??

        OP_STRING(1)=STRING(1:IEND)//';s <index ID#[1]>'
        OP_STRING(2)=BLANK(1:15)//'<on WS_ID#[1]>'
        OP_STRING(3)=BLANK(1:15)//'<rgb=RGB[black]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw polyline;m
C###  Parameter:    <(line/bezier)[line]>
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to draw the
C###    polyline on.
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Description:
C###    Draws a polyline specified by mouse.

        OP_STRING(1)=STRING(1:IEND)//';m <(line/bezier)[line]>'
        OP_STRING(2)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<rgb=RGB[black]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw polyline;c
C###  Parameter:    <index ID#[1]>
C###  Parameter:    <with optimisation_parameters>
C###  Parameter:    <endpoint WEIGHT#>
C###  Description:
C###    Draws a calculated polyline.

        OP_STRING(1)=STRING(1:IEND)//';c <index ID#[1]>'
        OP_STRING(2)=BLANK(1:15)//'<with optimisation_parameters>'
        OP_STRING(3)=BLANK(1:15)//'<endpoint WEIGHT#>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRPLIN',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL CHECKQ('CDLMPRSW',noco,1,CO,COQU,STRING,*1)
        FILIO=.FALSE.
        MOUSE=.FALSE.
        SEGME=.FALSE.
        CALCU=.FALSE.
        IF(ABBREV(COQU(noco,1),'P',1)) THEN
          FILIO=.TRUE.
          IOTYPE=1
          STATUS='NEW'
          BACKUP=.FALSE.
        ELSE IF(ABBREV(COQU(noco,1),'R',1)) THEN
          FILIO=.TRUE.
          IOTYPE=2
          STATUS='OLD'
        ELSE IF(ABBREV(COQU(noco,1),'W',1)) THEN
          FILIO=.TRUE.
          IOTYPE=3
          STATUS='NEW'
        ELSE IF(ABBREV(COQU(noco,1),'L',1)) THEN
          FILIO=.TRUE.
          IOTYPE=4
          STATUS='OLD'
        ELSE IF(ABBREV(COQU(noco,1),'M',1)) THEN
          MOUSE=.TRUE.
          IOTYPE=3
          STATUS='NEW'
        ELSE IF(ABBREV(COQU(noco,1),'S',1)) THEN
          SEGME=.TRUE.
        ELSE IF(ABBREV(COQU(noco,1),'C',1)) THEN
          CALCU=.TRUE.
        ENDIF
        IF(MOUSE.OR.SEGME) THEN
          CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
        ENDIF

        IF(CBBREV(CO,'INDEX',1,noco+1,NTCO,N3CO)) THEN
          INDEX_PLIN=IFROMC(CO(N3CO+1))
        ELSE
          INDEX_PLIN=1
        ENDIF

        IF(SEGME) THEN
        ELSE IF(MOUSE) THEN
          IF(CBBREV(CO,'BEZIER',2,noco+1,NTCO,N3CO)) THEN
            TYPE='BEZIER'
          ELSE
            TYPE='LINE'
          ENDIF
        ELSE IF(CALCU) THEN
          IF(CBBREV(CO,'WITH',2,noco+1,NTCO,N3CO)) THEN
            OPTIM=.TRUE.
          ELSE
            OPTIM=.FALSE.
          ENDIF
          IF(CBBREV(CO,'ENDPOINT',2,noco+1,NTCO,N3CO)) THEN
            WEIGHT=RFROMC(CO(N3CO+1))
          ELSE
            WEIGHT=1.0D0
          ENDIF
        ENDIF

        IF(FILIO) THEN
          CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.ipplin',STATUS,
     '      'DIRECT','FORMATTED',132,ERROR,*9999)
          CALL IPPLIN(ERROR,*9999)
C !!!! Note: need to define noplin_index
          CALL CLOSEF(IFILE,ERROR,*9999)
          IF(BACKUP) THEN
            CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.ipplin','OLD',
     '        'DIRECT','FORMATTED',132,ERROR,*9999)
            IOTYPE=2
            CALL IPPLIN(ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
          ENDIF

        ELSE IF(SEGME) THEN
          IW=IWK(1)
          IF(ADD) THEN
            NTPLIN=NTPLIN+1
          ELSE IF(NTPLIN.EQ.0) THEN
            NTPLIN=1
          ENDIF
          NOPLIN=NOPLIN_INDEX(INDEX_PLIN) !is segment no. of INDEX_PLIN
          IF(INDEX_PLIN_TYPE(INDEX_PLIN).EQ.1) THEN      !Piecewise linear
            NT_POINT=NT_PLIN_SECTIONS(INDEX_PLIN)+1
            CALL ACWK(iw,1,ERROR,*9999)
            CALL OPEN_SEGMENT(ISPLIN(iw,NOPLIN),ISEG,iw,'polyline_line',
     '        INDEX,INDEX_OLD,INDEX_PLIN,1,CSEG,ERROR,*9999)
            CALL POLYLINE(INDEX,iw,NT_POINT,PLIN_DATA(1,1,INDEX_PLIN),
     '        ERROR,*9999)
            CALL CLOSE_SEGMENT(ISPLIN(iw,NOPLIN),iw,ERROR,*9999)
            CALL DAWK(iw,1,ERROR,*9999)

          ELSE IF(INDEX_PLIN_TYPE(INDEX_PLIN).EQ.2) THEN !Bezier
            CALL ACWK(iw,1,ERROR,*9999)
            CALL OPEN_SEGMENT(ISPLIN(iw,NOPLIN),ISEG,iw,
     '        'polyline_bezier',
     '        INDEX,INDEX_OLD,INDEX_PLIN,1,CSEG,ERROR,*9999)
            DO nl=1,NT_PLIN_SECTIONS(INDEX_PLIN)
              DO i=1,4
                XBEZ(i)=PLIN_DATA(1,3*(nl-1)+i,INDEX_PLIN)
                YBEZ(i)=PLIN_DATA(2,3*(nl-1)+i,INDEX_PLIN)
              ENDDO
              CALL BEZIER_POINTS(PL,XBEZ(1),YBEZ(1),ERROR,*9999)
              CALL POLYLINE(1,iw,21,PL,ERROR,*9999)
C              CHAR=CFROMI(INDEX_PLIN,'(I2)')
              WRITE(CHAR,'(I2)') INDEX_PLIN
              CALL STRING_TRIM(CHAR,IBEG,IEND)
              CALL TEXT(1,iw,CHAR(IBEG:IEND),PL(1,1),ERROR,*9999)
            ENDDO
            CALL CLOSE_SEGMENT(ISPLIN(iw,NOPLIN),iw,ERROR,*9999)
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF

        ELSE IF(MOUSE) THEN
          iw=IWK(1)
          IF(TYPE(1:4).EQ.'LINE') THEN
            NT_PLIN=NT_PLIN+1  !increment polyline total
            NTPLIN=NTPLIN+1    !increment segment count
            INDEX_PLIN=NT_PLIN !is latest polyline
            INDEX_PLIN_TYPE(INDEX_PLIN)=1
            NOPLIN_INDEX(INDEX_PLIN)=NTPLIN !is segment no. of INDEX_PLIN
            CALL ACWK(iw,0,ERROR,*9999)
            CALL OPEN_SEGMENT(ISPLIN(iw,NTPLIN),ISEG,iw,'polyline_line',
     '        INDEX,INDEX_OLD,INDEX_PLIN,1,CSEG,ERROR,*9999)
            WRITE(OP_STRING,'('' >>Locate start of line'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL LOCATOR(INSTAT,0.0D0,XPTS(1),0.0D0,
     '        YPTS(1),ERROR,*9999)
            XSEGMENT_DATA(1,NTSG)=REAL(XPTS(1))                 !Node
            YSEGMENT_DATA(1,NTSG)=REAL(YPTS(1))                 !Node
            PLIN_DATA(1,1,INDEX_PLIN)=XPTS(1)
            PLIN_DATA(2,1,INDEX_PLIN)=YPTS(1)
            WRITE(OP_STRING,'('' >>Locate polyline points  '//
     '        '(use 2nd button to end)'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            NTPTS=1
            DO WHILE(INSTAT.EQ.1.AND.NTPTS.LT.18)
              CALL LOCATOR(INSTAT,
     '          XPTS(1),XPTS(2),YPTS(1),YPTS(2),ERROR,*9999)
              IF(INSTAT.EQ.1) THEN
                NTPTS=NTPTS+1
                PTS(1,1)=XPTS(1)
                PTS(2,1)=YPTS(1)
                PTS(1,2)=XPTS(2)
                PTS(2,2)=YPTS(2)
                CALL POLYLINE(1,iw,2,PTS,ERROR,*9999)
                XSEGMENT_DATA(NTPTS,NTSG)=REAL(XPTS(2))         !Node
                YSEGMENT_DATA(NTPTS,NTSG)=REAL(YPTS(2))         !Node
                PLIN_DATA(1,NTPTS,INDEX_PLIN)=XPTS(2)
                PLIN_DATA(2,NTPTS,INDEX_PLIN)=YPTS(2)
                XPTS(1)=XPTS(2)
                YPTS(1)=YPTS(2)
              ENDIF
            ENDDO
            CALL CLOSE_SEGMENT(ISPLIN(iw,NTPLIN),iw,ERROR,*9999)
            CALL DAWK(iw,0,ERROR,*9999)
            ISEGMENT_DATA(1,NTSG)=NTPTS-1
            NT_PLIN_SECTIONS(INDEX_PLIN)=NTPTS-1

          ELSE IF(TYPE(1:6).EQ.'BEZIER') THEN
            NT_PLIN=NT_PLIN+1  !increment polyline total
            NTPLIN=NTPLIN+1    !increment segment count
            INDEX_PLIN=NT_PLIN !is latest polyline
            INDEX_PLIN_TYPE(INDEX_PLIN)=2
            NOPLIN_INDEX(INDEX_PLIN)=NTPLIN !is segment no. of INDEX_PLIN
            CALL ACWK(iw,0,ERROR,*9999)
            CALL OPEN_SEGMENT(ISPLIN(iw,NTPLIN),ISEG,iw,
     '        'polyline_bezier',INDEX,INDEX_OLD,INDEX_PLIN,1,CSEG,
     '        ERROR,*9999)
            WRITE(OP_STRING,'('' >>Locate start of Bezier curve'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL LOCATOR(INSTAT,0.0D0,XPTS(1),
     '        0.0D0,YPTS(1),
     '        ERROR,*9999)
            XSEGMENT_DATA(1,NTSG)=REAL(XPTS(1)) !Records x-position of first node
            YSEGMENT_DATA(1,NTSG)=REAL(YPTS(1)) !Records y-position of first node
            XPOINTS(1)=XPTS(1) !First node
            YPOINTS(1)=YPTS(1) !First node
            WRITE(OP_STRING,'('' >>Locate points on Bezier curve  '//
     '        '(use 2nd button to end)'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            NTPTS=0
            NO_SECTIONS=0
            DO WHILE(INSTAT.EQ.1.AND.NTPTS.LT.18)
              CALL LOCATOR(INSTAT,
     '          XPTS(1),XPTS(2),YPTS(1),YPTS(2),ERROR,*9999)
              IF(INSTAT.EQ.1) THEN
                NO_SECTIONS=NO_SECTIONS+1
                NTPTS=NTPTS+3
                PTS(1,1)=XPTS(1)
                PTS(2,1)=YPTS(1)
                PTS(1,2)=XPTS(2)
                PTS(2,2)=YPTS(2)
                CALL POLYLINE(1,iw,2,PTS,ERROR,*9999)
                !control points initially at nodes
                XPOINTS(NTPTS-1)=XPOINTS(NTPTS-2)
                YPOINTS(NTPTS-1)=YPOINTS(NTPTS-2)
                XPOINTS(NTPTS  )=XPTS(2)
                YPOINTS(NTPTS  )=YPTS(2)
                !second node
                XPOINTS(NTPTS+1)=XPTS(2)
                YPOINTS(NTPTS+1)=YPTS(2)
                XPTS(1)=XPTS(2)
                YPTS(1)=YPTS(2)
              ENDIF
            ENDDO
            CALL CLOSE_SEGMENT(ISPLIN(iw,NTPLIN),iw,ERROR,*9999)
            CALL DAWK(iw,0,ERROR,*9999)
C DPN 17-12-97 - initialise arrays to pass through to BEZIER(), fixing
C                ftnchk error
            DO nolines = 1,NTLM
              ISL2BE_LOC(nolines)=0
              ISL3BE_LOC(nolines)=0
              ISN2BE_LOC(nolines)=0
              ISN3BE_LOC(nolines)=0
            ENDDO
C DPN 17-12-97 - update call to BEZIER()
            CALL BEZIER(INDEX,INDEX_PLIN,ISPLIN(iw,NTPLIN),ISEG,
     '        ISL2BE_LOC,ISL3BE_LOC,ISN2BE_LOC,ISN3BE_LOC,iw,
     '        NO_SECTIONS,XPOINTS,YPOINTS,CSEG,ERROR,*9999)
            NTPTS=NTPTS+1
            NT_PLIN_SECTIONS(INDEX_PLIN)=NO_SECTIONS
            IF(DOP) THEN
              WRITE(OP_STRING,
     '          '('' INDEX_PLIN='',I2,'' NT_PLIN_SECTIONS='',I2)')
     '          INDEX_PLIN,NT_PLIN_SECTIONS(INDEX_PLIN)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            DO nopts=1,NTPTS
              PLIN_DATA(1,nopts,INDEX_PLIN)=XPOINTS(nopts)
              PLIN_DATA(2,nopts,INDEX_PLIN)=YPOINTS(nopts)
              IF(DOP) THEN
                WRITE(OP_STRING,
     '            '('' PLIN_DATA(nj,'',I2,'','',I2,''):'',3E12.3)')
     '            nopts,INDEX_PLIN,(PLIN_DATA(nj,nopts,INDEX_PLIN),
     '            nj=1,NJT)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
          ENDIF

        ELSE IF(CALCU) THEN
          IF(OPTIM) THEN
            ALFA=PAOPTI(1) !defined by optimiser (initially by DEOPTI)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' ALFA='',E12.3)') ALFA
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE
            ALFA=1.0D0
          ENDIF
          IF(TYPE(1:4).EQ.'LINE') THEN
            NT_POINT=NT_PLIN_SECTIONS(INDEX_PLIN)+1
            nd=0
            DO no_point=1,NT_POINT-1
              DO n=0,10 !to create 11 points on line
                ETA=0.1D0*DBLE(n)
                XI=(ALFA+(1.0D0-ALFA)*ETA)*ETA
                nd=nd+1
                DO nj=1,NJT
                  ZD(nj,nd)=
     '              (1.d0-XI)*PLIN_DATA(nj,no_point,INDEX_PLIN)
     '                    +XI*PLIN_DATA(nj,no_point+1,INDEX_PLIN)
                  IF(n.EQ.0.OR.n.EQ.10) THEN
                    WD(nj,nd)=WEIGHT !for endpoints only
                  ELSE
                    WD(nj,nd)=1.0D0
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
            NDT=nd

          ELSE IF(TYPE(1:6).EQ.'BEZIER') THEN
            nd=0
            DO nl=1,NT_PLIN_SECTIONS(INDEX_PLIN)
              DO i=1,4
                XBEZ(i)=PLIN_DATA(1,3*(nl-1)+i,INDEX_PLIN)
                YBEZ(i)=PLIN_DATA(2,3*(nl-1)+i,INDEX_PLIN)
              ENDDO
              CALL BEZIER_POINTS(PL,XBEZ(1),YBEZ(1),ERROR,*9999)
              DO j=1,21
                nd=nd+1
                DO nj=1,NJT
                  ZD(nj,nd)=PL(nj,j)
                  IF(j.EQ.1.OR.j.EQ.21) THEN
                    WD(nj,nd)=WEIGHT !for endpoints only
                  ELSE
                    WD(nj,nd)=1.0D0
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
            NDT=nd
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('DRPLIN')
      RETURN
 9999 CALL ERRORS('DRPLIN',ERROR)
      CALL EXITS('DRPLIN')
      RETURN 1
      END


