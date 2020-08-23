      SUBROUTINE POLYMARKER_DYNAM(IBUNDLE,iw,LDRPTS,NTPTS,R_PTS,D_PTS,
     '  ERROR,*)

C#### Subroutine: POLYMARKER_DYNAM
C###  Description:
C###    POLYMARKER_DYNAM draws a polymarker on iw with index IBUNDLE. It
C###    differs from POLYMARKER in that the REAL data array for GX is
C###    dynamically allocated outside the subroutine.
C###  See-Also: PoLYMARKER

C**** D_PTS(1..3,nopts) contains the REAL*8 3D coords of each point.
C**** IF iw=4 coords are curvilinear coords, else rect. cart.
C**** If IBUNDLE is 0 the primitive will use the previously
C**** defined polymarker index.
C**** NOTE: If iw is 15 or 16 (postscript) IBUNDLE is reset to be black.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'disp00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBUNDLE,iw,LDRPTS,NTPTS
      REAL R_PTS(LDRPTS,2)
      REAL*8 D_PTS(3,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj,nopts,NJ1,NJ2,INDEX

      CALL ENTERS('POLYMARKER_DYNAM',*9999)

      CALL ASSERT(NTPTS.LE.LDRPTS,'>>NTPTS > LDRPTS',ERROR,*9999)

      IF(NTPTS.EQ.0) GOTO 9998

      INDEX=IBUNDLE

      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' iw='',I2,'' NTPTS='',I6)') iw,NTPTS
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nopts=1,NTPTS
          WRITE(OP_STRING,'('' D_PTS(nj,'',I6,''): '',3E12.3)')
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

        CALL POLYMARKER_GX(INDEX,NTPTS,R_PTS(1,1),R_PTS(1,2),
     '    ERROR,*9999)
      ENDIF

 9998 CALL EXITS('POLYMARKER_DYNAM')
      RETURN
 9999 CALL ERRORS('POLYMARKER_DYNAM',ERROR)
      CALL EXITS('POLYMARKER_DYNAM')
      RETURN 1
      END

