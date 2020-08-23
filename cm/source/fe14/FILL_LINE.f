      SUBROUTINE FILL_LINE(IBUNDLE,iw,NTPTS,D_PTS,ERROR,*)

C#### Subroutine: FILL_LINE
C###  Description:
C###    FILL_LINE draws a polyline on iw with index IBUNDLE.
C###    It calculates the colours better when drawing 1-d fields.
C###    It uses the same colour scale as FILL_AREA

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'disp00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBUNDLE,iw,NTPTS
      REAL*8 D_PTS(3,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj,nopts,NJ1,NJ2
      REAL R_PTS(NTPTS,2)

      CALL ENTERS('FILL_LINE',*9999)

      IF(NTPTS.EQ.0) GOTO 9998

      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' iw='',I2,'' NTPTS='',I6)') iw,NTPTS
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nopts=1,NTPTS
          WRITE(OP_STRING,
     '      '('' D_PTS(nj,'',I6,''): '',3E12.3)')
     '      nopts,(D_PTS(nj,nopts),nj=1,NJT)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      IF(IW.EQ.1) THEN
        NJ1=1
        NJ2=2
      ELSE IF(IW.EQ.2) THEN
        NJ1=2
        NJ2=3
      ELSE IF(IW.EQ.3) THEN
        NJ1=1
        NJ2=3
      ELSE IF(IW.EQ.7) THEN
        NJ1=2
        NJ2=3
      ELSE IF(IW.EQ.9) THEN
        NJ1=1
        NJ2=2
      ELSE IF(IW.EQ.10) THEN
        NJ1=1
        NJ2=2
      ELSE IF(IW.EQ.11) THEN
        NJ1=1
        NJ2=2
        DO nopts=1,NTPTS
          CALL MAP10(1,D_PTS(2,nopts))
        ENDDO
      ELSE IF(IW.EQ.12) THEN
        NJ1=1
        NJ2=2
      ELSE IF(IW.EQ.13) THEN
        NJ1=1
        NJ2=2
      ENDIF

      IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
        DO nopts=1,NTPTS
          R_PTS(nopts,1)=REAL(D_PTS(NJ1,nopts))   !is real x
          R_PTS(nopts,2)=REAL(D_PTS(NJ2,nopts))   !is real y
        ENDDO

        CALL FILL_LINE_GX(IBUNDLE,NTPTS,R_PTS,ERROR,*9999)
      ENDIF

 9998 CALL EXITS('FILL_LINE')
      RETURN
 9999 CALL ERRORS('FILL_LINE',ERROR)
      CALL EXITS('FILL_LINE')
      RETURN 1
      END


