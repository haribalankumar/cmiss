      SUBROUTINE SGAXES(INDEX,ISAXES,ISEG,iw,CSEG,ERROR,*)

C#### Subroutine: SGAXES
C###  Description:
C###    SGAXES creates axes segment ISAXES.

      IMPLICIT NONE
!     Parameter List
      INCLUDE 'axes00.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INTEGER INDEX,ISAXES,ISEG(*),iw
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER INDEX_OLD,naxis,nj,nox,noy,noz,np
      REAL*8 AXES(3,2,3),  !AXES(nj,np,naxis) defines pts for each axis
     '     TEXT_POS(3,3),  !TEXT_POS(nj,naxis) defines text positions
     '     TICS(3,2),      !TICS(nj,nt) defines tic polylines
     '     XXDIF,XXMAX,XXMIN,YYDIF,YYMAX,YYMIN,ZZDIF,ZZMAX,ZZMIN

      CALL ENTERS('SGAXES',*9999)
      CALL OPEN_SEGMENT(ISAXES,ISEG,iw,'AXES',INDEX,INDEX_OLD,
     '  1,1,CSEG,ERROR,*9999)

      XXMIN=DBLE(XMIN)+0.01d0*DBLE(XMAX-XMIN)
      XXMAX=DBLE(XMAX)-0.01d0*DBLE(XMAX-XMIN)
      YYMIN=DBLE(YMIN)+0.01d0*DBLE(YMAX-YMIN)
      YYMAX=DBLE(YMAX)-0.01d0*DBLE(YMAX-YMIN)
      ZZMIN=DBLE(ZMIN)+0.01d0*DBLE(ZMAX-ZMIN)
      ZZMAX=DBLE(ZMAX)-0.01d0*DBLE(ZMAX-ZMIN)
      XXDIF=XXMAX-XXMIN
      YYDIF=YYMAX-YYMIN
      ZZDIF=ZZMAX-ZZMIN
      DO nj=1,3
        DO np=1,2
          DO naxis=1,3
            AXES(nj,np,naxis)=0.0d0
          ENDDO
        ENDDO
      ENDDO

      IF(IWKT(iw).EQ.1) THEN !GKS
        AXES(1,1,1)=XXMIN
        AXES(1,2,1)=XXMAX
        AXES(2,1,2)=YYMIN
        AXES(2,2,2)=YYMAX
        AXES(3,1,3)=ZZMIN
        AXES(3,2,3)=ZZMAX
        TEXT_POS(1,1)=XXMAX-0.01D0*XXDIF
        TEXT_POS(2,1)=0.02D0*YYDIF
        TEXT_POS(3,1)=0.02D0*ZZDIF
        TEXT_POS(1,2)=0.02D0*XXDIF
        TEXT_POS(2,2)=YYMAX-0.01D0*YYDIF
        TEXT_POS(3,2)=0.02D0*ZZDIF
        TEXT_POS(1,3)=0.02D0*XXDIF
        TEXT_POS(2,3)=0.02D0*YYDIF
        TEXT_POS(3,3)=ZZMAX-0.01D0*ZZDIF

      ELSE IF(IWKT(iw).EQ.2) THEN !PHIGS
        AXES(1,2,1)=0.3D0*DBLE(XMAX)
        AXES(2,2,2)=0.3D0*DBLE(YMAX)
        AXES(3,2,3)=0.3D0*DBLE(ZMAX)
        DO nj=1,3
          DO naxis=1,3
            TEXT_POS(nj,naxis)=AXES(nj,2,naxis)
          ENDDO
        ENDDO
      ENDIF

      CALL POLYLINE(INDEX,iw,2,AXES(1,1,1),ERROR,*9999)
      CALL POLYLINE(INDEX,iw,2,AXES(1,1,2),ERROR,*9999)
      CALL TEXT(1,iw,'X',TEXT_POS(1,1),ERROR,*9999)
      CALL TEXT(1,iw,'Y',TEXT_POS(1,2),ERROR,*9999)
      IF(NJT.EQ.3) THEN
        CALL POLYLINE(INDEX,iw,2,AXES(1,1,3),ERROR,*9999)
        CALL TEXT(1,iw,'Z',TEXT_POS(1,3),ERROR,*9999)
      ENDIF

C     tics on the x axis
      IF(NTX.GT.0) THEN
        TICS(2,1)=0.0D0
        TICS(2,2)=YYDIF/60.0D0
        TICS(3,1)=0.0D0
        TICS(3,2)=ZZDIF/60.0D0
        DO nox=1,NTX
          TICS(1,1)=XLIST(nox)
          TICS(1,2)=XLIST(nox)
          CALL POLYLINE(INDEX,iw,2,TICS,ERROR,*9999)
        ENDDO
      ENDIF
C     tics on the Y axis
      IF(NTY.GT.0) THEN
        TICS(1,1)=0.0D0
        TICS(1,2)=XXDIF/60.0D0
        TICS(3,1)=0.0D0
        TICS(3,2)=ZZDIF/60.0D0
        DO noy=1,NTY
          TICS(2,1)=YLIST(noy)
          TICS(2,2)=YLIST(noy)
          CALL POLYLINE(INDEX,iw,2,TICS,ERROR,*9999)
        ENDDO
      ENDIF
C     tics on the Z axis
      IF(NTZ.GT.0) THEN
        TICS(2,1)=0.0D0
        TICS(2,2)=YYDIF/60.0D0
        TICS(1,1)=0.0D0
        TICS(1,2)=XXDIF/60.0D0
        DO noz=1,NTZ
          TICS(3,1)=XLIST(noz)
          TICS(3,2)=XLIST(noz)
          CALL POLYLINE(INDEX,iw,2,TICS,ERROR,*9999)
        ENDDO
      ENDIF

      CALL CLOSE_SEGMENT(ISAXES,iw,ERROR,*9999)

      CALL EXITS('SGAXES')
      RETURN
 9999 CALL ERRORS('SGAXES',ERROR)
      CALL EXITS('SGAXES')
      RETURN 1
      END


