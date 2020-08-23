      SUBROUTINE SGHIST(INDEX,ISEG,ISHIST,np,NTHIST,CSEG,TMAX,TMIN,
     '  XHIST,YHIST,YMAX,YMIN,ERROR,*)

C#### Subroutine: SGHIST
C###  Description:
C###    SGHIST creates segment ISHIST. iw is 10.
C###    If np=0 ISHIST is segment containing axes and labels.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISHIST,np,NTHIST
      REAL*8 TMAX,TMIN,XHIST(*),YHIST(*),YMAX,YMIN,YTOT,TTOT
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER i,IBEG,IEND,INDEX_OLD,iw,NCENTR,nohist
      REAL*8 PTS(3,5000)
      CHARACTER CHAR*10

      CALL ENTERS('SGHIST',*9999)

      iw=10
      CALL OPEN_SEGMENT(ISHIST,ISEG,iw,'HIST',INDEX,INDEX_OLD,
     '  np,1,CSEG,ERROR,*9999)

      IF(DOP) THEN
        WRITE(OP_STRING,'('' Index = '',I5,'', Index_old = '',I5)')
     '    INDEX,INDEX_OLD
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NTHIST = '',I5)') NTHIST
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' TMIN = '',D11.4,'', TMAX = '',D11.4)')
     '    TMIN,TMAX
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' YMIN = '',D11.4,'', YMAX = '',D11.4)')
     '    YMIN,YMAX
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      YTOT=DABS(YMIN)+DABS(YMAX)
      TTOT=DABS(TMIN)+DABS(TMAX)
      IF(np.EQ.0) THEN !Draw axes
C       X axis
        IF(YMIN.LE.0.0d0) THEN
          PTS(1,1)=0.0d0
          PTS(2,1)=-YMIN/YTOT
          PTS(1,2)=1.0d0
          PTS(2,2)=PTS(2,1)
        ELSE
          PTS(1,1)=0.0d0
          PTS(2,1)=0.0d0
          PTS(1,2)=1.0d0
          PTS(2,2)=0.0d0
        ENDIF
        CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
C       Y axis
        PTS(1,1)=0.0d0
        PTS(2,1)=0.0d0
        PTS(1,2)=0.0d0
        PTS(2,2)=1.0d0
        CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
C       X tic marks
        IF(YMIN.LE.0.0d0) THEN
          PTS(2,1)=-0.03d0-YMIN/YTOT
          PTS(2,2)=0.01d0-YMIN/YTOT
        ELSE
          PTS(2,1)=-0.03d0
          PTS(2,2)=0.01d0
        ENDIF
        DO i=1,10
          PTS(1,1)=DBLE(i)/10.0d0
          PTS(1,2)=PTS(1,1)
          CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
          IF(MOD(I,2).EQ.0) THEN
            WRITE(CHAR,'(F7.2)') PTS(1,1)*(TMIN+TMAX)-TMIN
            CALL STRING_TRIM(CHAR,IBEG,IEND)
            IF(YMIN.LE.0.0d0) THEN
              PTS(2,1)=-0.08d0-YMIN/YTOT
            ELSE
              PTS(2,1)=-0.08d0
            ENDIF
            CALL TEXT(1,iw,CHAR(IBEG:IEND),PTS(1,1),ERROR,*9999)
            IF(YMIN.LE.0.0d0) THEN
              PTS(2,1)=-0.03d0-YMIN/YTOT
              PTS(2,2)=0.01d0-YMIN/YTOT
            ELSE
              PTS(2,1)=-0.03d0
              PTS(2,2)=0.01d0
            ENDIF
          ENDIF
        ENDDO
C       Y tic marks
        PTS(1,1)=-0.02d0
        PTS(1,2)=0.01d0
        DO i=0,5
          PTS(2,1)=DBLE(i)/5.0d0
          PTS(2,2)=PTS(2,1)
          CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)
          WRITE(CHAR,'(F7.2)') PTS(2,1)*YTOT+YMIN
          CALL STRING_TRIM(CHAR,IBEG,IEND)
          PTS(1,1)=-0.07d0
          CALL TEXT(1,iw,CHAR(IBEG:IEND),PTS(1,1),ERROR,*9999)
          PTS(1,1)=-0.02d0
        ENDDO
C       Plot lables
        PTS(1,1)=0.5d0
        PTS(2,1)=-0.2d0
        CALL TEXT(1,iw,'Nodal time history',PTS(1,1),ERROR,*9999)
        PTS(1,1)=1.05d0
        IF(YMIN.LE.0.0d0) THEN
          PTS(2,1)=-YMIN/YTOT
        ELSE
          PTS(2,1)=0.0d0
        ENDIF
        CALL TEXT(1,iw,'t',PTS(1,1),ERROR,*9999)

      ELSE IF(np.GT.0) THEN !Draw time plots
        DO nohist=1,NTHIST
          PTS(1,nohist)=(XHIST(nohist)-TMIN)/TTOT
          PTS(2,nohist)=(YHIST(nohist)-YMIN)/YTOT
        ENDDO
        CALL POLYLINE(INDEX,iw,NTHIST,PTS,ERROR,*9999)
        NCENTR=NTHIST/2
        WRITE(CHAR,'(I5)') np
        CALL STRING_TRIM(CHAR,IBEG,IEND)
        PTS(1,1)=(XHIST(NCENTR)-TMIN)/TTOT
        PTS(2,1)=(YHIST(NCENTR)-YMIN)/YTOT
        CALL TEXT(1,iw,CHAR(IBEG:IEND),PTS(1,1),ERROR,*9999)
      ENDIF

      CALL CLOSE_SEGMENT(ISHIST,iw,ERROR,*9999)

      CALL EXITS('SGHIST')
      RETURN
 9999 CALL ERRORS('SGHIST',ERROR)
      CALL EXITS('SGHIST')
      RETURN 1
      END


