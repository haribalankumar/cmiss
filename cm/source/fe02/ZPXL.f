      SUBROUTINE ZPXL(IP,NPL,nr,NTDX,nx,DL,XL,ZP,ERROR,*)

C#### Subroutine: ZPXL
C###  Description:
C###    ZPXL creates world coordinate polyline array
C###    XL(nj,nodx),nodx=1,ntdx from ZP(nk,1,nh,np,nc).

C**** If NTDX=0 on entry, NTDX is set =2 when straight line plotting is
C**** appropriate, or is set to 20 when the line may be curved.
C**** If IP=0 XL is returned in same coordinate system as XP,
C**** otherwise XL is returned in rectangular cartesian coordinates.
C**** (took out XI3,YP 30/10/90 AAY)
C**** NOTE: This routine was changed on 23/7/96 to use
C**** ZP(nk,1,nh,np,nc) and not ZP(nk,1,nj,np,nc)

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'bspln00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IP,NPL(5,0:3),nr,NTDX,nx
      REAL*8 DL(3),XL(3,20),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,nc,nh,njj,nk,nodx
      REAL*8 PLXIZ,XI,XX(3),ZZ(3)

      CALL ENTERS('ZPXL',*9999)
C      nc=1 !temporary cpb 23/11/94
      IF(NTDX.EQ.0) THEN
        IF(ITYP10(nr).EQ.1.AND.NPL(1,1).LE.1.AND.NPL(1,2).LE.1
     '    .AND.NPL(1,3).LE.1) THEN
          NTDX=2
        ELSE IF(ITYP10(nr).EQ.1.AND.NPL(1,1).EQ.5
     '    .AND.MBSPL(1).EQ.1) THEN
          NTDX=2
        ELSE
          NTDX=20
        ENDIF
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        nc=1
        WRITE(OP_STRING,*) 'NPL nj=0: ',(NPL(I,0),I=1,5)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO I=1,2
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nh=NH_LOC(njj,nx)
            WRITE(OP_STRING,'('' ZP(nk,1,'',I1,'','',I4,'
     '        //''',nc): '',4E11.3)')
     '        nh,NPL(1+I,1),(ZP(nk,1,nh,NPL(1+I,1),nc),nk=1,4)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO
        WRITE(OP_STRING,*) 'DL: ',(DL(I),I=1,2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      DO nodx=1,NTDX
        XI=DBLE(nodx-1)/DBLE(NTDX-1)
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nh=NH_LOC(njj,nx)
          XX(njj)=PLXIZ(nh,njj,NPL,nr,1,DL,XI,ZP)
        ENDDO
C newe
        IF(IP.EQ.0) THEN
          CALL XZ(1,XX,ZZ)
        ELSE
          CALL XZ(ITYP10(nr),XX,ZZ)
        ENDIF

        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          XL(njj,nodx)=ZZ(njj)
        ENDDO
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' XL(nj,'',I2,''): '',3E11.4)') nodx,
     '      (XL(NH_LOC(njj,nx),nodx),njj=1,NJ_LOC(NJL_GEOM,0,nr))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
      ENDDO

      CALL EXITS('ZPXL')
      RETURN
 9999 CALL ERRORS('ZPXL',ERROR)
      CALL EXITS('ZPXL')
      RETURN 1
      END


