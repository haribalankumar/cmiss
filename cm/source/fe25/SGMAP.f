      SUBROUTINE SGMAP(INDEX,ISEG,ISMAP,iw,CSEG,ERROR,*)

C#### Subroutine: SGMAP
C###  Description:
C###    SGMAP creates segment ISMAP.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'map000.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISMAP,iw
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,INDEX_OLD,ne1,ne2,nodx,NTDX
      REAL*8 X(3),XL(3,20)
      CHARACTER CHAR*5

      CALL ENTERS('SGMAP',*9999)
      CALL OPEN_SEGMENT(ISMAP,ISEG,iw,'MAP',INDEX,INDEX_OLD,
     '  1,1,CSEG,ERROR,*9999)

      PROJEC=MAP_PROJEC

      IF(PROJEC(1:6).EQ.'HAMMER'.OR.PROJEC(1:5).EQ.'POLAR') THEN
        IF(iw.EQ.4) THEN
          WRITE(CHAR,'(F5.0)') RMAP
          CALL STRING_TRIM(MAP_PROJEC,IBEG,IEND)
          X(1)=-0.95D0
          X(2)=0.95D0
          PROJEC='RECTANGULAR'
          CALL TEXT(1,iw,MAP_PROJEC(IBEG:IEND)//' projection'//
     '      '   Rotation='//CHAR(1:4)//' degrees',X,ERROR,*9999)
          PROJEC=MAP_PROJEC
C old 21/3/97
C       ELSE IF(iw.EQ.5) THEN
C         COLOUR=254
        ENDIF
        IF(PROJEC(1:6).EQ.'HAMMER'.AND.RMU.GT.0.0D0) THEN
C ***     Draw line at theta=0
          NTDX=20
          DO nodx=1,NTDX
            XL(2,nodx)=DBLE(nodx-1)/DBLE(NTDX-1)*RMU/180.0D0*PI
            XL(3,nodx)=0.0D0
          ENDDO
          CALL POLYLINE(INDEX,iw,NTDX,XL,ERROR,*9999)
C ***     Draw line at theta=2*pi
          DO nodx=1,NTDX
            XL(2,nodx)=DBLE(nodx-1)/DBLE(NTDX-1)*RMU/180.0D0*PI
            XL(3,nodx)=1.9999D0*PI
          ENDDO
          CALL POLYLINE(INDEX,iw,NTDX,XL,ERROR,*9999)
        ENDIF
        IF(LABEL) THEN
          X(2)=0.68D0*PI
          X(3)=1.5D0*PI
          CALL TEXT(1,iw,'LAD',X,ERROR,*9999)
          X(2)=0.68D0*PI
          X(3)=PI
          CALL TEXT(1,iw,'LV',X,ERROR,*9999)
          X(2)=0.68D0*PI
          X(3)=0.5D0*PI
          CALL TEXT(1,iw,'PDA',X,ERROR,*9999)
        ENDIF

      ELSE IF(PROJEC(1:11).EQ.'RECTANGULAR') THEN
C ***   Draw box
        XL(1,1)=-1.0D0*XI1MAP
        XL(2,1)=-1.0D0*XI2MAP
        XL(1,2)=       XI1MAP
        XL(2,2)=-1.0D0*XI2MAP
        XL(1,3)=       XI1MAP
        XL(2,3)=       XI2MAP
        XL(1,4)=-1.0D0*XI1MAP
        XL(2,4)=       XI2MAP
        XL(1,5)=-1.0D0*XI1MAP
        XL(2,5)=-1.0D0*XI2MAP
        CALL POLYLINE(INDEX,iw,5,XL,ERROR,*9999)

      ELSE IF(PROJEC(1:2).EQ.'XI') THEN
C ***   Draw element boundaries
        DO ne1=1,NET_XI1
          MXI1=ne1
          DO ne2=1,NET_XI2
            MXI2=ne2
            XL(1,1)=0.0D0
            XL(2,1)=0.0D0
            XL(1,2)=1.0D0
            XL(2,2)=0.0D0
            XL(1,3)=1.0D0
            XL(2,3)=1.0D0
            XL(1,4)=0.0D0
            XL(2,4)=1.0D0
            XL(1,5)=0.0D0
            XL(2,5)=0.0D0
            CALL POLYLINE(INDEX,iw,5,XL,ERROR,*9999)
          ENDDO
        ENDDO
      ENDIF

      CALL CLOSE_SEGMENT(ISMAP,iw,ERROR,*9999)

      CALL EXITS('SGMAP')
      RETURN
 9999 CALL ERRORS('SGMAP',ERROR)
      CALL EXITS('SGMAP')
      RETURN 1
      END


