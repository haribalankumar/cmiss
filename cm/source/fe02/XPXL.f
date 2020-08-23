      SUBROUTINE XPXL(TYPE,IP,NPL,nr,NTDX,DL,XL,XP,ERROR,*)

C#### Subroutine: XPXL
C###  Description:
C###    XPXL creates world coordinate polyline array
C###    XL(nj,nodx),nodx=1,ntdx from XP(nk,nv,nj,np).
C**** If NTDX=0 on entry,
C**** NTDX is set to 2 when straight line plotting is
C**** appropriate, or is set to 20 when the line may be curved.
C**** If IP=0 XL is returned in same coordinate system as XP,
C**** otherwise XL is returned in rectangular cartesian coordinates.
C**** (took out XI3,YP 30/10/90 AAY)
C**** TYPE is either UNDEFORMED geometry lines or FIELD variable lines.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'bspln00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IP,NPL(5,0:3),nr,NTDX
      REAL*8 DL(3),XL(3,20),XP(NKM,NVM,NJM,NPM)
      CHARACTER TYPE*(*),ERROR*(*)
!     Local Variables
      INTEGER I,j,nj,nk,nodx
      REAL*8 PLXI,XI,XX(3),ZZ(3)

      CALL ENTERS('XPXL',*9999)

      IF(NTDX.LE.2) THEN
        IF(ITYP10(nr).EQ.1.AND.NPL(1,1).LE.1.AND.NPL(1,2).LE.1
     '    .AND.NPL(1,3).LE.1) THEN
          NTDX=2
        ELSE IF(ITYP10(nr).EQ.1.AND.NPL(1,1).EQ.5.
     '    AND.MBSPL(1).EQ.1) THEN
          NTDX=2
        ELSE
          NTDX=20
        ENDIF
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,*) 'NPL nj=0: ',(NPL(I,1),I=1,5)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO I=1,2
          DO j=1,NJT
            IF(TYPE(1:10).EQ.'UNDEFORMED') THEN
              nj=NJ_LOC(NJL_GEOM,j,nr)
            ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
              nj=NJ_LOC(NJL_FIEL,j,nr)
            ENDIF
            WRITE(OP_STRING,'('' XP(nk,1,'',I1,'','',I4,''): '','
     '        //'4E11.3)') nj,NPL(1+I,1),(XP(nk,1,nj,NPL(1+I,1)),nk=1,4)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
        WRITE(OP_STRING,*) 'DL: ',(DL(I),I=1,2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      DO nodx=1,NTDX
        XI=DBLE(nodx-1)/DBLE(NTDX-1)
        DO j=1,NJT
          IF(TYPE(1:10).EQ.'UNDEFORMED') THEN
            nj=NJ_LOC(NJL_GEOM,j,nr)
          ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
            nj=NJ_LOC(NJL_FIEL,j,nr)
          ENDIF
          XX(j)=PLXI(nj,NPL,nr,1,DL,XI,XP)
        ENDDO
        IF(IP.EQ.0) THEN
          CALL XZ(1,XX,ZZ)
        ELSE
          CALL XZ(ITYP10(nr),XX,ZZ)
        ENDIF

        DO j=1,NJT
          XL(j,nodx)=ZZ(j)
        ENDDO
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' XL(j,'',I2,''): '',3E11.4)') nodx,
     '      (XL(j,nodx),j=1,NJT)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
      ENDDO

      CALL EXITS('XPXL')
      RETURN
 9999 CALL ERRORS('XPXL',ERROR)
      CALL EXITS('XPXL')
      RETURN 1
      END


