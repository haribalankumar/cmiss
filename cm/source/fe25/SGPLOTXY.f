      SUBROUTINE SGPLOTXY(POLYTYPE,INDEX,ISEG,ISEGNUM,ISPLOTXY,
     '  iw,NPTS,D_PTS,X,Y,CSEG,ERROR,*)

C#### Subroutine: SGPLOTXY
C###  Description:
C###    SGPLOTXY plots the vector X versus the vector Y. This routine is
C###    a generic plotting routine and can handle both linear and log
C###    scales. SGPLOTXY creates segment ISPLOTXY. The common block file
C###    plot02.cmn must be included when using this routine.
C****
C****   Arguments
C****   =========
C****
C****   POLYTYPE (input) EXTERNAL
C****   Specifies whether a POLYLINE_DYNAM or POLYMARKER_DYNAM is used.
C****
C****   INDEX    (input) INTEGER
C****   The POLYTYPE index (colour) for each series.
C****
C****   ISEG     (input) INTEGER array
C****   Specifies whether the segment is not yet created/created but not
C****   visible/ created and visible.
C****
C****   ISEGNUM  (input) INTEGER
C****   Specifies the graphical segment number. This routine has been
C****   designed to use only two segments to simplify the procedure.
C****   If ISEGNUM = 1, the segment containing the x-y data is created,
C****   if ISEGNUM = 2, the segment containing the axes and labels is
C****   created. NOTE : The nature of this routine means that multiple
C****   plots can be created on segment 1.
C****
C****   ISPLOTXY (input) INTEGER
C****   Specifies the graphical segment number.
C****
C****   iw       (input) INTEGER
C****   Specifies the workstaion number.
C****
C****   NPTS     (input) INTEGER
C****   The number of data points for the arrays X and Y.
C****
C****   D_PTS    (workspace) REAL*8 array, dimension (3,NPTS)
C****   Workspace array for storing GX data points.
C****
C****   X        (input) DOUBLE PRECISION array, dimension (NPTS)
C****
C****   Y        (input) DOUBLE PRECISION array, dimension (NPTS)
C****
C****   Common block
C****   ============
C****
C****   XMAX     (input) REAL*8
C****   Specifies the maximum x-axis value.
C****
C****   XMIN     (input) REAL*8
C****   Specifies the minimum x-axis value.
C****
C****   YMAX     (input) REAL*8
C****   Specifies the maximum y-axis value.
C****
C****   YMIN     (input) REAL*8
C****   Specifies the minimum y-axis value.
C****
C****   NEXTPLOT (input) CHARACTER*7
C****   Specifies whether ISEGNUM = 1 should be created.
C****   = 'REPLACE':  a new plot will be created on ISEGNUM = 1.
C****   = 'ADD'    :  the current plot is overlayed on the visible
C****                 ISEGNUM = 1.
C****   NOTE : If NEXTPLOT = 'ADD' the axis bounds should not be changed.
C****
C****   TITLE    (input) CHARACTER*80
C****   Specifies the title of the plot.
C****
C****   XLABEL   (input) CHARACTER*80
C****   Specifies the x-label of the plot.
C****
C****   YLABEL   (input) CHARACTER*80
C****   Specifies the y-label of the plot.
C****
C****   XSCALE   (input) CHARACTER *6
C****   Specifies the scale of the x-axis.
C****   = 'LINEAR':  a linear scale.
C****   = 'LOG'   :  a logrithmic scale.
C****
C****   YSCALE   (input) CHARACTER *6
C****   Specifies the scale of the y-axis.
C****   = 'LINEAR':  a linear scale.
C****   = 'LOG'   :  a logrithmic scale.
C****
C****  ISHOLD    (input) LOGICAL
C****  Specifies whether the current plot properties are held.
C****  = .TRUE. :   the current plot segment (ISEGNUM = 1) is held open
C****               for subsequent plots.
C****  = .FALSE.:   the current plot segment (ISEGNUM = 1) is closed,
C****               such that subsequent calls to SGPLOTXY will erase
C****               the current segment.
C****  For plotting only one data set on the current axis NEXTPLOT =
C****  'REPLACE' and ISHOLD = .FALSE. such that the segment is both
C****  opened and closed. For plotting multiple data sets on the
C****  current axis, for the first call NEXTPLOT = 'REPLACE' and ISHOLD =
C****  .TRUE.. For the remaining data sets NEXTPLOT = 'ADD' to ensure
C****  the current segment is not reopened. ISHOLD = .FALSE. on the
C****  final data set to ensure that the current segment is closed.
C****
C****  XGRID     (input) LOGICAL
C****  Specifies the x-grid properties
C****  = .TRUE. :  adds the x-grid lines to the current axis.
C****  = .FALSE.:  no x-grid lines are created.
C****
C****  YGRID     (input) LOGICAL
C****  Specifies the y-grid properties
C****  = .TRUE. :  adds the y-grid lines to the current axis.
C****  = .FALSE.:  no y-grid lines are created.

      IMPLICIT NONE
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'plot02.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISEGNUM,ISPLOTXY,iw,NPTS
      REAL*8 D_PTS(3,*),X(*),Y(*)
      CHARACTER CSEG(*)*(*),ERROR*(*)
      EXTERNAL POLYTYPE
!     Local Variables
      INTEGER i,IBEG,IDATA,IEND,INDEX_AXES,INDEX_GRID,INDEX_OLD,SIGFIG
      INTEGER*4 R_PTS_PTR
      REAL*8 CPT(2),DELTA,MAXVAL(2),MINVAL(2),ONE,RANGE,RDATA,ZERO
      CHARACTER CHAR*10,FORMAT_STR
      LOGICAL LOOP,MINOR
      PARAMETER (ONE=1.0d0,ZERO=0.0d0)
!     Functions
      INTEGER CM_FLOOR,ILOG10,INDEX_POLYLINE

      CALL ENTERS('SGPLOTXY',*9999)

      CALL ASSERT((XMAX.NE.XMIN).AND.(YMAX.NE.YMIN),'>>Non valid data',
     '  ERROR,*9999)

      MINVAL(1)=XMIN
      MINVAL(2)=YMIN
      MAXVAL(1)=XMAX
      MAXVAL(2)=YMAX

      IF(ISEGNUM.EQ.2) THEN
        CALL OPEN_SEGMENT(ISPLOTXY,ISEG,iw,'PLOTXY',INDEX,INDEX_OLD,1,1,
     '    CSEG,ERROR,*9999)
        ! Title
        CPT(1)=0.50d0
        CPT(2)=1.05d0
        CALL STRING_TRIM(TITLE,IBEG,IEND)
        CALL TEXT(0,iw,TITLE(IBEG:IEND),CPT,ERROR,*9999)
        ! X-label
        CPT(2)=-0.10d0
        CALL STRING_TRIM(XLABEL,IBEG,IEND)
        CALL TEXT(0,iw,XLABEL(IBEG:IEND),CPT,ERROR,*9999)
        ! Y-label
        CPT(1)=-0.10d0
        CPT(2)=0.50d0
        CALL STRING_TRIM(YLABEL,IBEG,IEND)
        CALL TEXT(0,iw,YLABEL(IBEG:IEND),CPT,ERROR,*9999)

        ! Axes
        INDEX_AXES=INDEX_POLYLINE(0,'SOLID','WIDTH2','BLACK')
        INDEX_GRID=INDEX_POLYLINE(0,'DOTTED','WIDTH1','GREY')

        DO i=1,2 ! Loop over x and y axis
          IF((XSCALE(1:6).EQ.'LINEAR'.AND.i.EQ.1).OR.
     '      (YSCALE(1:6).EQ.'LINEAR'.AND.i.EQ.2)) THEN
            RANGE=(MAXVAL(i)-MINVAL(i))/10.0d0
            RDATA=10**DBLE(ILOG10(RANGE))
            IF(RANGE.GE.3.50d0*RDATA) THEN
              IF(RANGE.GE.7.50d0*RDATA) THEN
                DELTA=10.0d0*RDATA
              ELSE
                DELTA=5.0d0*RDATA
              ENDIF
            ELSE
              IF(RANGE.GE.1.50d0*RDATA) THEN
                DELTA=2.0d0*RDATA
              ELSE
                DELTA=RDATA
              ENDIF
            ENDIF
            RDATA=DBLE(CM_FLOOR(MINVAL(i)/DELTA)+1)*DELTA
            ! Determine format
            SIGFIG=ILOG10(DELTA)
            IF(SIGFIG.GE.3.OR.SIGFIG.LE.-3) THEN
              FORMAT_STR='E'
            ELSEIF(SIGFIG.GE.0) THEN
              FORMAT_STR='I'
            ELSEIF(SIGFIG.GE.-2) THEN
              FORMAT_STR='F'
            ENDIF
            DO WHILE(RDATA.LE.MAXVAL(i))
              D_PTS(i,1)=(RDATA-MINVAL(i))/(MAXVAL(i)-MINVAL(i))
              D_PTS(i,2)=D_PTS(i,1)
              D_PTS(3-i,1)=ZERO
              IF((XGRID.AND.i.EQ.2).OR.
     '          (YGRID.AND.i.EQ.1)) THEN
                ! Grid
                D_PTS(3-i,2)=ONE
                CALL POLYLINE(INDEX_GRID,iw,2,D_PTS,ERROR,*9999)
              ENDIF
              ! Tics
              D_PTS(3-i,2)=0.01d0
              CALL POLYLINE(INDEX_AXES,iw,2,D_PTS,ERROR,*9999)
              D_PTS(3-i,1)=0.99d0
              D_PTS(3-i,2)=ONE
              CALL POLYLINE(INDEX_AXES,iw,2,D_PTS,ERROR,*9999)
              CPT(i)=D_PTS(i,1)
              CPT(3-i)=-0.05d0
              CALL SGPLOTXYLABEL(IBEG,IEND,2,RDATA,FORMAT_STR,CHAR,
     '          ERROR,*9999)
              CALL TEXT(0,iw,CHAR(IBEG:IEND),CPT,ERROR,*9999)
              RDATA=RDATA+DELTA
            ENDDO
          ELSEIF((XSCALE(1:3).EQ.'LOG'.AND.i.EQ.1).OR.
     '      (YSCALE(1:3).EQ.'LOG'.AND.i.EQ.2))THEN
            IDATA=ILOG10(MAXVAL(i))
            IF(DLOG10(MAXVAL(i)/MINVAL(i)).GT.5.0d0) THEN
              ! Plot major decades
              MINOR=.FALSE.
              RDATA=10.0d0**DBLE(IDATA)
            ELSE
              ! Plot minor decades
              MINOR=.TRUE.
              RDATA=CM_FLOOR(MAXVAL(i)/10.0d0**DBLE(IDATA))*10.0d0**
     '          DBLE(IDATA)
            ENDIF
            LOOP=.TRUE.
            DO WHILE(LOOP)
              D_PTS(i,1)=DLOG10(RDATA/MINVAL(i))/DLOG10(MAXVAL(i)
     '            /MINVAL(i))
              IF(D_PTS(i,1).GE.ZERO) THEN
                D_PTS(i,2)=D_PTS(i,1)
                D_PTS(3-i,1)=ZERO
                IF((XGRID.AND.i.EQ.2).OR.
     '             (YGRID.AND.i.EQ.1)) THEN
                  ! Grid
                  D_PTS(3-i,2)=ONE
                  CALL POLYLINE(INDEX_GRID,iw,2,D_PTS,ERROR,*9999)
                ENDIF

                ! Tics
                IF((.NOT.MINOR.AND.(MOD(IDATA,5).EQ.0)).OR.
     '            (MINOR.AND.(DABS(RDATA-10.0d0
     '            **DBLE(IDATA))
     '            /(DABS(RDATA)+ONE).LE.CONVERG_TOL))) THEN
                  DELTA=0.01d0
                  CPT(i)=D_PTS(i,1)
                  CPT(3-i)=-0.05d0
                  CALL SGPLOTXYLABEL(IBEG,IEND,1,DBLE(IDATA),'E',CHAR,
     '              ERROR,*9999)
                  CALL TEXT(0,iw,CHAR(IBEG:IEND),CPT,ERROR,*9999)
                  IF(MINOR) IDATA=IDATA-1
                ELSE
                  DELTA=0.005d0
                ENDIF
                D_PTS(3-i,2)=DELTA
                CALL POLYLINE(INDEX_AXES,iw,2,D_PTS,ERROR,*9999)
                D_PTS(3-i,1)=ONE-DELTA
                D_PTS(3-i,2)=ONE
                CALL POLYLINE(INDEX_AXES,iw,2,D_PTS,ERROR,*9999)
                IF(MINOR) THEN
                  RDATA=RDATA-10.0d0**DBLE(IDATA)
                ELSE
                  IDATA=IDATA-1
                  RDATA=10.0d0**DBLE(IDATA)
                ENDIF
              ELSE
                LOOP=.FALSE.
              ENDIF
            ENDDO
          ENDIF
        ENDDO

        ! Draw bounding box
        D_PTS(1,1)=ZERO
        D_PTS(2,1)=ZERO
        D_PTS(1,2)=ONE
        D_PTS(2,2)=ZERO
        CALL POLYLINE(INDEX_AXES,iw,2,D_PTS,ERROR,*9999)
        D_PTS(1,1)=ONE
        D_PTS(2,1)=ONE
        CALL POLYLINE(INDEX_AXES,iw,2,D_PTS,ERROR,*9999)
        D_PTS(1,2)=ZERO
        D_PTS(2,2)=ONE
        CALL POLYLINE(INDEX_AXES,iw,2,D_PTS,ERROR,*9999)
        D_PTS(1,1)=ZERO
        D_PTS(2,1)=ZERO
        CALL POLYLINE(INDEX_AXES,iw,2,D_PTS,ERROR,*9999)
        CALL CLOSE_SEGMENT(ISPLOTXY,iw,ERROR,*9999)
      ELSE
        IF(NEXTPLOT(1:7).EQ.'REPLACE') THEN
          CALL OPEN_SEGMENT(ISPLOTXY,ISEG,iw,'PLOTXY',INDEX,INDEX_OLD,1,
     '      1,CSEG,ERROR,*9999)
        ENDIF
        ! Draw Plot
        IF(XSCALE(1:6).EQ.'LINEAR') THEN
          DO i=1,NPTS
            D_PTS(1,i)=(X(i)-MINVAL(1))/(MAXVAL(1)-MINVAL(1))
          ENDDO
        ELSEIF(XSCALE(1:3).EQ.'LOG') THEN
          DO i=1,NPTS
           D_PTS(1,i)=DLOG10(X(i)/MINVAL(1))/DLOG10(MAXVAL(1)/MINVAL(1))
          ENDDO
        ENDIF
        IF(YSCALE(1:6).EQ.'LINEAR') THEN
          DO i=1,NPTS
            D_PTS(2,i)=(Y(i)-MINVAL(2))/(MAXVAL(2)-MINVAL(2))
          ENDDO
        ELSEIF(YSCALE(1:3).EQ.'LOG') THEN
          DO i=1,NPTS
            D_PTS(2,i)=DLOG10(Y(i)/MINVAL(2))/DLOG10(MAXVAL(2)/
     '        MINVAL(2))
          ENDDO
        ENDIF
        R_PTS_PTR=0
        CALL ALLOCATE_MEMORY(NPTS*2,0,SPTYPE,R_PTS_PTR,MEM_INIT,
     '    ERROR,*9999)
        CALL POLYTYPE(INDEX,iw,NPTS,NPTS,%VAL(R_PTS_PTR),D_PTS,
     '    ERROR,*9999)
        CALL FREE_MEMORY(R_PTS_PTR,ERROR,*9999)

        IF(.NOT.ISHOLD) THEN
          CALL CLOSE_SEGMENT(ISPLOTXY,iw,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('SGPLOTXY')
      RETURN
 9999 CALL ERRORS('SGPLOTXY',ERROR)
      CALL EXITS('SGPLOTXY')
      RETURN 1
      END

