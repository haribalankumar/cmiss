      SUBROUTINE SGSCAL(INDEX,ISEG,ISSCAL,iw,NOSCAL,COLOUR,CSEG,ERROR,*)

C#### Subroutine: SGSCAL
C###  Description:
C###    SGSCAL creates segment ISSCAL for scale.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'colo00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'map000.cmn'
      INCLUDE 'scal00.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISSCAL,iw,NOSCAL
      CHARACTER CSEG(*)*(*),ERROR*(*)
      LOGICAL COLOUR
!     Local Variables
      INTEGER IBEG,IEND,INDEX_FIELD,INDEX_OLD,no_scale,NT_SCALE
      REAL*8 PT(3),PTS(3,4),TIC(3,2),ZVAL
      CHARACTER CHAR12*12

      CALL ENTERS('SGSCAL',*9999)
      CALL OPEN_SEGMENT(ISSCAL,ISEG,iw,'SCAL',INDEX,INDEX_OLD,
     '  NOSCAL,1,CSEG,ERROR,*9999)

      IF(COLOUR_WS) THEN
        NT_SCALE=25
      ELSE
        NT_SCALE=16
      ENDIF

      DO no_scale=1,NT_SCALE
        ZVAL=ZMINI+DBLE(no_scale-1)/DBLE(NT_SCALE-1)*(ZMAXI-ZMINI)
        IF(iw.EQ.4) THEN
          IF(DABS(ZVAL).EQ.0.0D0) THEN
            CHAR12='  0.0 '
          ELSE IF(DABS(ZVAL).LT.0.001D0) THEN
            CHAR12='<0.001'
          ELSE IF(DABS(ZVAL).GE.0.001D0.AND.DABS(ZVAL).LT.10.0D0) THEN
            WRITE(CHAR12,'(F6.3)') ZVAL
          ELSE IF(DABS(ZVAL).GE.10.0D0.AND.DABS(ZVAL).LT.100.0D0) THEN
            WRITE(CHAR12,'(F6.2)') ZVAL
          ELSE IF(DABS(ZVAL).GE.100.0D0.AND.DABS(ZVAL).LT.1.0D3) THEN
            WRITE(CHAR12,'(F6.1)') ZVAL
          ELSE IF(DABS(ZVAL).GE.1.0D3.AND.DABS(ZVAL).LT.1.0D4) THEN
            WRITE(CHAR12,'(F6.0)') ZVAL
          ELSE
            CHAR12='>10000'
          ENDIF
        ELSE
          WRITE(CHAR12,'(E12.3)') ZVAL
        ENDIF
        IF(NJT.EQ.2) THEN
          PTS(1,1)=DBLE(XMIN)+DBLE(no_scale-1)/DBLE(NT_SCALE-1)
     '                       *DBLE(XMAX-XMIN)                !x1
          PTS(2,1)=DBLE(YMAX)-0.05d0*DBLE(YMAX-YMIN)         !y1
          PTS(1,2)=PTS(1,1)+DBLE(XMAX-XMIN)/DBLE(NT_SCALE-1) !x2
          PTS(2,2)=PTS(2,1)                                  !y2=y1
          PTS(1,3)=PTS(1,2)                                  !x3=x2
          PTS(2,3)=DBLE(YMAX)                                !y3
          PTS(1,4)=PTS(1,1)                                  !x4=x1
          PTS(2,4)=PTS(2,3)                                  !y4=y3
C PJH 2Feb95 PTS(1,1)=DBLE(XMAX)-0.22D0*DBLE(XMAX-XMIN)
C         PTS(2,1)=DBLE(YMIN)+(0.1D0+0.8D0*DBLE(no_scale-1)/
C    '      DBLE(NT_SCALE-1))*DBLE(YMAX-YMIN)
C         PTS(1,2)=PTS(1,1)+0.05D0*DBLE(XMAX-XMIN)
C         PTS(2,2)=PTS(2,1)
C         PTS(1,3)=PTS(1,2)
C         PTS(2,3)=PTS(2,2)+0.8D0/DBLE(NT_SCALE-1)*DBLE(YMAX-YMIN)
C         PTS(1,4)=PTS(1,1)
C         PTS(2,4)=PTS(2,3)
C         PT(1)=PTS(1,2)+0.07D0*DBLE(XMAX-XMIN)
C         PT(2)=0.5D0*(PTS(2,1)+PTS(2,4))
        ELSE IF(NJT.EQ.3) THEN
          IF(iw.EQ.1) THEN
            PTS(1,1)=DBLE(XMAX)-0.22D0*DBLE(XMAX-XMIN)
            PTS(3,1)=DBLE(YMIN)+(0.1D0+0.8D0*DBLE(no_scale-1)/
     '        DBLE(NT_SCALE-1))*DBLE(YMAX-YMIN)
            PTS(1,2)=PTS(1,1)+0.05D0*DBLE(XMAX-XMIN)
            PTS(3,2)=PTS(3,1)
            PTS(1,3)=PTS(1,2)
            PTS(3,3)=PTS(3,2)+0.8D0/DBLE(NT_SCALE-1)*DBLE(YMAX-YMIN)
            PTS(1,4)=PTS(1,1)
            PTS(3,4)=PTS(3,3)
            PT(1)=PTS(1,2)+0.07D0*DBLE(XMAX-XMIN)
            PT(3)=0.5D0*(PTS(3,1)+PTS(3,4))
          ELSE IF(iw.EQ.2) THEN
          ELSE IF(iw.EQ.3) THEN
          ELSE IF(iw.EQ.4) THEN
            PTS(1,1)=-0.9D0+1.7D0*DBLE(no_scale-1)/DBLE(NT_SCALE-1)
            PTS(1,2)=PTS(1,1)+1.6D0/DBLE(NT_SCALE-1)
            PTS(1,3)=PTS(1,2)
            PTS(1,4)=PTS(1,1)
            PTS(2,1)=0.6D0
            PTS(2,2)=PTS(2,1)
            PTS(2,3)=0.7D0
            PTS(2,4)=PTS(2,3)
            PT(1)=0.5D0*(PTS(1,1)+PTS(1,2))
            PT(2)=0.75D0
            TIC(1,1)=PT(1)
            TIC(1,2)=PT(1)
            TIC(2,1)=0.705D0
            TIC(2,2)=0.72D0
          ENDIF
        ENDIF
        IF(COLOUR) THEN
! PJH 9Jan96 INDEX_FIELD=249-NINT((ZVAL-ZMINI)/ZDIFF*232.0d0)
C       Keep index between 45 and 200
          INDEX_FIELD=200-NINT((ZVAL-ZMINI)/ZDIFF*155.0d0)
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' ZVAL='',E12.3,'' INDEX_FIELD='',I3)')
     '        ZVAL,INDEX_FIELD
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE
          INDEX_FIELD=16-NINT((ZVAL-ZMINI)/ZDIFF*15.0D0)
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' ZVAL='',E12.3,'' INDEX_FIELD='',I2)')
     '        ZVAL,INDEX_FIELD
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

        PROJEC='RECTANGULAR'
        CALL FILL_AREA(INDEX_FIELD,iw,4,PTS,ERROR,*9999)

        IF(iw.NE.4) THEN
          IF(no_scale.EQ.1) THEN
            PT(1)=DBLE(XMIN)+0.055d0*DBLE(XMAX-XMIN)
            PT(2)=DBLE(YMAX)-0.075d0*DBLE(YMAX-YMIN)
            CALL TEXT(INDEX,iw,CHAR12(1:12),PT,ERROR,*9999)
          ELSE IF(no_scale.EQ.NT_SCALE) THEN
            PT(1)=DBLE(XMAX)-0.090d0*DBLE(XMAX-XMIN)
            PT(2)=DBLE(YMAX)-0.075d0*DBLE(YMAX-YMIN)
            CALL TEXT(INDEX,iw,CHAR12(1:12),PT,ERROR,*9999)
          ENDIF
        ELSE IF(iw.EQ.4) THEN
          IF(MOD(no_scale,2).EQ.0) THEN !label every 2nd fill-area
            CALL STRING_TRIM(CHAR12,IBEG,IEND)
            CALL TEXT(INDEX,iw,CHAR12(IBEG:IEND),PT,ERROR,*9999)
            CALL POLYLINE(1,iw,2,TIC,ERROR,*9999)
          ENDIF
        ENDIF

      ENDDO !no_scale

      CALL CLOSE_SEGMENT(ISSCAL,iw,ERROR,*9999)

      CALL EXITS('SGSCAL')
      RETURN
 9999 CALL ERRORS('SGSCAL',ERROR)
      CALL EXITS('SGSCAL')
      RETURN 1
      END


