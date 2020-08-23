      SUBROUTINE SGLINE(INDEX,ISEG,ISLINE,ISLINO,iw,NLLIST,NOLINE,NPL,
     '  nr,nx,CSEG,DEFORM,DL,SOLID,STATIC,XP,ZP,ERROR,*)

C#### Subroutine: SGLINE
C###  Description:
C###    SGLINE creates line segment ISLINE containing lines
C###    and line segment ISLINO containing line numbers.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'map000.cmn'
      INCLUDE 'time00.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISLINE,ISLINO,iw,NLLIST(0:NLM),NOLINE,
     '  NPL(5,0:3,NLM),nr,nx
      REAL*8 DL(3,NLM),XP(NKM,NVM,NJM,NPM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),DEFORM*(*),ERROR*(*),SOLID*(*)
      LOGICAL STATIC
!     Local Variables
      INTEGER i,IBEG,IEND,INDEX_OLD,nl,NMID,nolist,NTDX
      REAL*8 XL(3,20)
      CHARACTER CHAR4*4,CLABEL*52

      CALL ENTERS('SGLINE',*9999)

C *** Define line segments
      CLABEL='LINE '//DEFORM(1:3)//' '//SOLID(1:3)
      IF(STATIC) THEN
        CLABEL(43:52)='          '
      ELSE
        WRITE(CLABEL(43:52),'(E10.3)') TIME
      ENDIF
      CALL OPEN_SEGMENT(ISLINE,ISEG,iw,CLABEL,INDEX,INDEX_OLD,
     '  NOLINE,1,CSEG,ERROR,*9999)

      IF(iw.EQ.4) PROJEC=MAP_PROJEC
      IF(DEFORM(1:8).EQ.'DEFORMED') THEN
        NTDX=20
      ELSE
        NTDX=0
      ENDIF

      DO nolist=1,NLLIST(0)
        nl=NLLIST(nolist)
        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' Region nr='',I2,'' Line nl= '',I6)')
     '      nr,nl
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(iw.NE.4) THEN
          IF(DEFORM(1:8).EQ.'DEFORMED') THEN
            CALL ZPXL(1,NPL(1,0,nl),nr,NTDX,nx,DL(1,nl),XL,ZP,
     '        ERROR,*9999)
          ELSE
            CALL XPXL(DEFORM,1,NPL(1,0,nl),nr,NTDX,DL(1,nl),XL,XP,
     '        ERROR,*9999)
          ENDIF
        ELSE IF(iw.EQ.4.AND.NPL(1,0,nl).LE.2) THEN
          IF(DEFORM(1:8).EQ.'DEFORMED') THEN
            CALL ZPXL(0,NPL(1,0,nl),nr,NTDX,nx,DL(1,nl),XL,ZP,
     '      ERROR,*9999)
          ELSE
            CALL XPXL(DEFORM,0,NPL(1,0,nl),nr,NTDX,DL(1,nl),XL,XP,
     '        ERROR,*9999)
          ENDIF
        ELSE IF(iw.EQ.5.OR.iw.EQ.6) THEN
          IF(PROJEC(1:6).EQ.'HAMMER') THEN
            IF(DEFORM(1:8).EQ.'DEFORMED') THEN
              CALL ZPXL(0,NPL(1,0,nl),nr,NTDX,nx,DL(1,nl),XL,ZP,
     '          ERROR,*9999)
            ELSE
              CALL XPXL(DEFORM,0,NPL(1,0,nl),nr,NTDX,DL(1,nl),XL,XP,
     '          ERROR,*9999)
            ENDIF
          ENDIF
        ENDIF
        IF(NJT.EQ.1) THEN
          DO i=1,20
            XL(2,i)=0.0d0
          ENDDO
        ENDIF
        CALL POLYLINE(INDEX,iw,NTDX,XL,ERROR,*9999)
      ENDDO

      CALL CLOSE_SEGMENT(ISLINE,iw,ERROR,*9999)

C *** Define line numbers
      IF(iw.LE.2.OR.iw.EQ.4.AND.NOLINE.EQ.1) THEN
        CALL OPEN_SEGMENT(ISLINO,ISEG,iw,'LINO',INDEX,INDEX_OLD,
     '    1,2,CSEG,ERROR,*9999)

        DO nolist=1,NLLIST(0)
          nl=NLLIST(nolist)
          IF(iw.LE.2) THEN
            IF(DEFORM(1:8).EQ.'DEFORMED') THEN
              CALL ZPXL(1,NPL(1,0,nl),nr,NTDX,nx,DL(1,nl),XL,ZP,
     '          ERROR,*9999)
            ELSE
              CALL XPXL(DEFORM,1,NPL(1,0,nl),nr,NTDX,DL(1,nl),XL,XP,
     '          ERROR,*9999)
            ENDIF
          ELSE IF(iw.EQ.4) THEN
            IF(DEFORM(1:8).EQ.'DEFORMED') THEN
              CALL ZPXL(0,NPL(1,0,nl),nr,NTDX,nx,DL(1,nl),XL,ZP,
     '          ERROR,*9999)
            ELSE
              CALL XPXL(DEFORM,0,NPL(1,0,nl),nr,NTDX,DL(1,nl),XL,XP,
     '          ERROR,*9999)
            ENDIF
          ENDIF
          IF(NTDX.EQ.2) THEN
            XL(1,1)=(XL(1,1)+XL(1,2))/2.0D0
            XL(2,1)=(XL(2,1)+XL(2,2))/2.0D0
            IF(NJT.EQ.3) THEN
              XL(3,1)=(XL(3,1)+XL(3,2))/2.0D0
            ENDIF
          ELSE
            NMID=NTDX/2
            XL(1,1)=XL(1,NMID)
            XL(2,1)=XL(2,NMID)
            IF(NJT.EQ.3) THEN
              XL(3,1)=XL(3,NMID)
            ENDIF
          ENDIF
          WRITE(CHAR4,'(I4)') nl
          CALL STRING_TRIM(CHAR4,IBEG,IEND)
C CPB 14/3/94 Don't draw line numbers as they are invisible anyway
          CALL TEXT(1,iw,CHAR4(IBEG:IEND),XL(1,1),ERROR,*9999)
        ENDDO

        CALL CLOSE_SEGMENT(ISLINO,iw,ERROR,*9999)
      ENDIF

      CALL EXITS('SGLINE')
      RETURN
 9999 CALL ERRORS('SGLINE',ERROR)
      CALL EXITS('SGLINE')
      RETURN 1
      END


