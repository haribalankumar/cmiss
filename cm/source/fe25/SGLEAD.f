      SUBROUTINE SGLEAD(INDEX,ISEG,ISLEAD,nlead,nolead,NTSIG,FACTOR,
     '  XSIG,YSIG,XMAX,XMIN,YMAX,YMIN,CSEG,
     '  PTITLE,PLOTAXESSCALE,PLOTLABELLEAD,PLOTTITLE,ERROR,*)

C#### Subroutine: SGLEAD
C###  Description:
C###    SGLEAD creates segment ISLead. iw is 9.
C###    If nlead=0 ISLEAD is segment containing axes and labels.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'lead00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISLEAD,nlead,nolead,NTSIG
      REAL*8 FACTOR,XSIG(*),YSIG(*),XMAX,XMIN,YMAX,YMIN
      CHARACTER CSEG(*)*(*),PTITLE*(*),ERROR*(*)
      LOGICAL PLOTAXESSCALE,PLOTLABELLEAD,PLOTTITLE
!     Local Variables
      INTEGER i,IBEG,IEND,INDEX_OLD,iw,NCENTR,nosig
      INTEGER*4 R_PTS_PTR
      REAL*8 YRANGE
      REAL*8 PTS(3,8192) !1 kHz sampling frequency for 8 seconds
      CHARACTER CHAR*10,TITLE*15,ERROR1*255

      CALL ENTERS('SGLEAD',*9999)

      iw=9
      CALL OPEN_SEGMENT(ISLEAD,ISEG,iw,'LEAD',INDEX,INDEX_OLD,
     '  nolead,1,CSEG,ERROR,*9999)

      YRANGE=SQRT((YMAX-YMIN)**2)

      IF(nlead.EQ.0) THEN !Draw axes
C       X axis
        PTS(1,1)=0.0d0
        PTS(2,1)=0.5d0
        PTS(1,2)=1.0d0
        PTS(2,2)=0.5d0
        CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
C       Y axis
        PTS(1,1)=0.0d0
        PTS(2,1)=0.0d0
        PTS(1,2)=0.0d0
        PTS(2,2)=1.0d0
        CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
        PTS(2,1)=0.50d0
        PTS(2,2)=0.53d0
C       X axis tics + scale
        DO i=1,10
          PTS(1,1)=DBLE(i)/10.0d0
          PTS(1,2)=PTS(1,1)
          CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
C MPN 21Mar2002: added title, labellead, axesscale, vmin, vmax options.
C old         IF(MOD(I,2).EQ.0) THEN
          IF(PLOTAXESSCALE.AND.MOD(I,2).EQ.0) THEN
            WRITE(CHAR,'(E9.2)') XMIN+PTS(1,1)*(XMAX-XMIN)
            CALL STRING_TRIM(CHAR,IBEG,IEND)
            CALL TEXT(1,iw,CHAR(IBEG:IEND),PTS(1,1),ERROR,*9999)
          ENDIF
        ENDDO
        IF(DOP) THEN
          WRITE(OP_STRING,*) ' ymin,ymax=',ymin,ymax
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
C       Y axis tics + scale
        DO i=1,5
          PTS(1,1)=0.0d0
          PTS(2,1)=DBLE(i-1)/4.0d0
          PTS(1,2)=0.02d0
          PTS(2,2)=PTS(2,1)
          CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
C MPN 21Mar2002: added title, labellead, axesscale, vmin, vmax options.
          IF(PLOTAXESSCALE) THEN
            WRITE(CHAR,'(E9.2)') YMIN+PTS(2,1)*YRANGE
            CALL STRING_TRIM(CHAR,IBEG,IEND)
            PTS(1,1)=-0.01d0
            CALL TEXT(1,iw,CHAR(IBEG:IEND),PTS(1,1),ERROR,*9999)
          ENDIF
        ENDDO
C MPN 21Mar2002: added title, labellead, axesscale, vmin, vmax options.
        IF(PLOTTITLE) THEN
          PTS(1,1)=0.5d0
          PTS(2,1)=-0.10d0
          CALL TEXT(1,iw,PTITLE,PTS(1,1),ERROR,*9999)
        ENDIF

      ELSE IF(nlead.GT.0) THEN !Draw time plots
        DO nosig=1,NTSIG
          PTS(1,nosig)=XSIG(nosig)
          PTS(2,nosig)=(YSIG(nosig)-YMIN)/YRANGE*FACTOR
        ENDDO
        R_PTS_PTR=0
        CALL ALLOCATE_MEMORY(NTSIG*2,1,SPTYPE,R_PTS_PTR,MEM_INIT,ERROR,
     '    *9999)
        CALL POLYLINE_DYNAM(INDEX,iw,NTSIG,NTSIG,%VAL(R_PTS_PTR),
     '    PTS,ERROR,*9999)
        CALL FREE_MEMORY(R_PTS_PTR,ERROR,*9999)
        NCENTR=NTSIG/2
        TITLE=LEADTITLE(nlead)
        CALL STRING_TRIM(TITLE,IBEG,IEND)
C MPN 21Mar2002: added title, labellead, axesscale, vmin, vmax options.
        IF(PLOTLABELLEAD) THEN
          PTS(1,1)=XSIG(NCENTR)
          PTS(2,1)=(YSIG(NCENTR)-YMIN)/YRANGE*FACTOR
          CALL TEXT(1,iw,TITLE(IBEG:IEND),PTS(1,1),ERROR,*9999)
        ENDIF
      ENDIF

      CALL CLOSE_SEGMENT(ISLEAD,iw,ERROR,*9999)

      CALL EXITS('SGLEAD')
      RETURN
 9999 CALL ERRORS('SGLEAD',ERROR)
      IF(R_PTS_PTR.NE.0) THEN
        CALL FREE_MEMORY(R_PTS_PTR,ERROR1,*1111)
      ENDIF
 1111 CALL EXITS('SGLEAD')
      RETURN 1
      END


