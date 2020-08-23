      SUBROUTINE SGCLOC(INDEX,ISCLOC,ISEG,iw,CSEG,THETA,ERROR,*)

C#### Subroutine: SGCLOC
C###  Description:
C###    SGCLOC creates clock segment ISCLOC.
C###    THETA (0..2*pi) gives hand position, with zero being vertical.

      IMPLICIT NONE
      INCLUDE 'cloc00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER INDEX,ISCLOC,ISEG(*),iw
      REAL*8 THETA
      CHARACTER CSEG(*)*(*),ERROR*(*)
!    Local Variables
      INTEGER INDEX_OLD,nj
      REAL*8 PTS(3,2),RADIUS

      CALL ENTERS('SGCLOC',*9999)
      CALL OPEN_SEGMENT(ISCLOC,ISEG,iw,'CLOC',INDEX,INDEX_OLD,
     '  1,1,CSEG,ERROR,*9999)

      CALL CIRCLE('POLYLINE',iw,R_CLOCK,X_CLOCK,20,ERROR,*9999)
      DO nj=1,NJT
        PTS(nj,1)=X_CLOCK(nj) !is clock centre
      ENDDO
      RADIUS=0.9D0*R_CLOCK
      IF(NJT.EQ.2) THEN
        PTS(1,2)=PTS(1,1)+RADIUS*DSIN(THETA)
        PTS(2,2)=PTS(2,1)+RADIUS*DCOS(THETA)
        PTS(3,2)=0.0D0
      ELSE IF(NJT.EQ.3) THEN
        IF(iw.EQ.1) THEN
          PTS(1,2)=PTS(1,1)+RADIUS*DSIN(THETA)
          PTS(2,2)=PTS(2,1)
          PTS(3,2)=PTS(3,1)+RADIUS*DCOS(THETA)
        ELSE IF(iw.EQ.2) THEN
          PTS(1,2)=PTS(1,1)
          PTS(2,2)=PTS(2,1)+RADIUS*DSIN(THETA)
          PTS(3,2)=PTS(3,1)+RADIUS*DCOS(THETA)
        ELSE IF(iw.EQ.3) THEN
c         CALL ZZ(Z,Z,TRANS)
          PTS(1,2)=PTS(1,1)+RADIUS*DSIN(THETA)
          PTS(2,2)=PTS(2,1)+RADIUS*DCOS(THETA)
          PTS(3,2)=PTS(3,1)
        ENDIF
      ENDIF
      CALL POLYLINE(INDEX,iw,2,PTS,ERROR,*9999)

      CALL CLOSE_SEGMENT(ISCLOC,iw,ERROR,*9999)
      CALL EXITS('SGCLOC')
      RETURN
 9999 CALL ERRORS('SGCLOC',ERROR)
      CALL EXITS('SGCLOC')
      RETURN 1
      END


