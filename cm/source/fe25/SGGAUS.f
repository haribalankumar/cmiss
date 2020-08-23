      SUBROUTINE SGGAUS(INDEX,ISEG,ISGAUS,iw,MXI,ne,ng,NITB,nr,
     '  CSEG,TYPE,X,XG,XIG,DEFORM,ERROR,*)

C#### Subroutine: SGGAUS
C###  Description:
C###    SGGAUS creates Gauss point segment ISGAUS.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'map000.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISGAUS,iw,MXI(2),ne,ng,NITB,nr
      REAL*8 X(*),XG(NJM,NUM),XIG(NIM)
      CHARACTER CSEG(*)*(*),ERROR*(*),TYPE*(*)
      LOGICAL DEFORM
!     Local Variables
      INTEGER INDEX_OLD,nj
      REAL*8 AXES(3,2),A_VECTOR(3),B_VECTOR(3),C_VECTOR(3),DAXES,Z(3)
      CHARACTER CHAR3*3

      CALL ENTERS('SGGAUS',*9999)
      WRITE(CHAR3,'(I3)') ng
      CALL OPEN_SEGMENT(ISGAUS,ISEG,iw,'GA'//CHAR3(1:3),INDEX,INDEX_OLD,
     '  ne,1,CSEG,ERROR,*9999)

      IF(iw.EQ.4) THEN
        CALL ASSERT(.NOT.DEFORM,'>>>ERROR: Cannot draw deformed Gauss'
     '    //' pts on the projection window',ERROR,*9999)
        PROJEC=MAP_PROJEC
        IF(PROJEC(1:2).EQ.'XI') THEN
          MXI1=MXI(1) !bottom left coords
          MXI2=MXI(2) !..for Xi map projection
          Z(1)=XIG(1) !Xi coordinates
          Z(2)=XIG(2) !..of Gauss point
        ELSE
          DO nj=1,3
            Z(nj)=X(nj)
          ENDDO
        ENDIF
      ELSE
        IF(DEFORM) THEN
          CALL XZ(ITYP11(nr),X,Z)
        ELSE
          CALL XZ(ITYP10(nr),X,Z)
        ENDIF
      ENDIF

      IF(TYPE(1:5).EQ.'POINT') THEN
        CALL POLYMARKER(INDEX,iw,1,Z,ERROR,*9999)

      ELSE IF(TYPE(1:4).EQ.'AXES') THEN
        CALL ASSERT(.NOT.DEFORM,'>>>ERROR: Not implemented',ERROR,*9999)

        DO nj=1,NJT
          AXES(nj,1)=Z(nj)
        ENDDO
        DAXES=DBLE(DIAG)/40.0D0

        IF(NJT.EQ.2) THEN
          AXES(1,2)=AXES(1,1)+DAXES
          AXES(2,2)=AXES(2,1)
          CALL POLYLINE(INDEX,iw,2,AXES,ERROR,*9999)
          AXES(1,2)=AXES(1,1)
          AXES(2,2)=AXES(2,1)+DAXES
          CALL POLYLINE(INDEX,iw,2,AXES,ERROR,*9999)

        ELSE IF(NJT.EQ.3) THEN
          CALL MAT_VEC_NG(NITB,nr,A_VECTOR,B_VECTOR,C_VECTOR,
     '      XG,ERROR,*9999)
          DO nj=1,NJT
            AXES(nj,2)=AXES(nj,1)+DAXES*A_VECTOR(nj)
          ENDDO
          CALL POLYLINE(INDEX,iw,2,AXES,ERROR,*9999)
          DO nj=1,NJT
            AXES(nj,2)=AXES(nj,1)+DAXES*B_VECTOR(nj)
          ENDDO
          CALL POLYLINE(INDEX,iw,2,AXES,ERROR,*9999)
          DO nj=1,NJT
            AXES(nj,2)=AXES(nj,1)+DAXES*C_VECTOR(nj)
          ENDDO
          CALL POLYLINE(INDEX,iw,2,AXES,ERROR,*9999)
        ENDIF
      ENDIF

      CALL CLOSE_SEGMENT(ISGAUS,iw,ERROR,*9999)

      CALL EXITS('SGGAUS')
      RETURN
 9999 CALL ERRORS('SGGAUS',ERROR)
      CALL EXITS('SGGAUS')
      RETURN 1
      END


