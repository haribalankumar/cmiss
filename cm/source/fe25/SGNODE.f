      SUBROUTINE SGNODE(INDEX,ISEG,ISNONO,iw,NHP,NKH,np,nr,nx,NYNP,
     '  FIX,XWC,YWC,ZWC,CSEG,ERROR,*)

C#### Subroutine: SGNODE
C###  Description:
C###    SGNODE creates new node segment ISNONO(iw,np).

      IMPLICIT NONE
      INCLUDE 'call00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'map000.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISNONO,iw,NHP(NPM),NKH(NHM,NPM,NCM),np,nr,
     '  nx,NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 XWC,YWC,ZWC
      CHARACTER CSEG(*)*(*),ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER IBEG,IEND,INDEX_OLD,nc,nh,nhx,nk,nv,ny
      REAL*8 DELTA,PTS(3,2)
      CHARACTER CHAR*4
      LOGICAL FOUND

      CALL ENTERS('SGNODE',*9999)

      CALL OPEN_SEGMENT(ISNONO,ISEG,iw,'NONO',INDEX,INDEX_OLD,
     '  np,1,CSEG,ERROR,*9999)

      WRITE(CHAR,'(I4)') np
      CALL STRING_TRIM(CHAR,IBEG,IEND)
      PTS(1,1)=XWC
      PTS(2,1)=YWC
      PTS(3,1)=ZWC
      IF(iw.EQ.4) THEN
        PROJEC=MAP_PROJEC
        IF(PROJEC(1:2).EQ.'XI') THEN
c         MXI1=MXI(1,?) !bottom left coords
c         MXI2=MXI(2,?) !..for Xi map projection
        ELSE
          CALL ZX(ITYP10(nr),PTS,PTS)
        ENDIF
      ENDIF
      IF(iw.NE.5.and.iw.NE.6) THEN
        CALL TEXT(INDEX,iw,CHAR(IBEG:IEND),PTS,ERROR,*9999)
      ELSE
        CALL POLYMARKER(INDEX,iw,1,PTS,ERROR,*9999)
      ENDIF

C *** Draw line segments to indicate boundary conditions
      IF(CALL_INIT) THEN
        IF(NHP(np).GT.0.AND.ITYP10(nr).EQ.1) THEN
          DELTA=0.006D0*DBLE(DIAG)
          IF(iw.EQ.1) THEN
            nc=1 !Temporary AJP 18-12-91
            nv=1
C            NYTOT=0
C            DO N=1,np-1
C              DO nhx=1,NHP(N)
C                nh=NH_LOC(nhx,nx)
C                NYTOT=NYTOT+NKH(nh,N,nc)
C              ENDDO
C            ENDDO

            FOUND=.FALSE.
            DO nhx=1,NHP(np)
              nh=NH_LOC(nhx,nx)
              DO nk=1,NKH(nh,np,nc)
                ny=NYNP(nk,nv,nh,np,0,nc,nr)
                IF(FIX(ny,1)) FOUND=.TRUE.
                IF(FIX(ny,2)) FOUND=.TRUE.
              ENDDO
            ENDDO

            IF(NHP(np).EQ.1) THEN !only one dependent variable
              IF(FOUND) THEN
                PTS(1,1)=XWC-DELTA
                PTS(1,2)=XWC+DELTA
                PTS(2,1)=YWC-1.7D0*DELTA
                PTS(2,2)=PTS(2,1)
                PTS(3,1)=ZWC-1.7D0*DELTA
                PTS(3,2)=PTS(3,1)
                CALL POLYLINE(5,iw,2,PTS,ERROR,*9999)
              ENDIF
            ELSE IF(NHP(np).EQ.2) THEN !two dependent variables
              IF((ITYP1(nr,nx).EQ.5.AND.FOUND).OR !finite elasticity
     '          .(ITYP1(nr,nx).NE.5.AND.FOUND)) THEN !other
                PTS(1,1)=XWC-1.4D0*DELTA
                PTS(1,2)=PTS(1,1)
                PTS(2,1)=YWC-DELTA
                PTS(2,2)=YWC+DELTA
                PTS(3,1)=ZWC-DELTA
                PTS(3,2)=ZWC+DELTA
                CALL POLYLINE(5,iw,2,PTS,ERROR,*9999)
              ENDIF
              IF((ITYP1(nr,nx).EQ.5.AND.FOUND).OR !finite elasticity
     '          .(ITYP1(nr,nx).NE.5.AND.FOUND)) THEN !other
                PTS(1,1)=XWC-DELTA
                PTS(1,2)=XWC+DELTA
                PTS(2,1)=YWC-1.7D0*DELTA
                PTS(2,2)=PTS(2,1)
                PTS(3,1)=ZWC-1.7D0*DELTA
                PTS(3,2)=PTS(3,1)
                CALL POLYLINE(5,iw,2,PTS,ERROR,*9999)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF !call_init

      CALL CLOSE_SEGMENT(ISNONO,iw,ERROR,*9999)

      CALL EXITS('SGNODE')
      RETURN
 9999 CALL ERRORS('SGNODE',ERROR)
      CALL EXITS('SGNODE')
      RETURN 1
      END


