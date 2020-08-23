      SUBROUTINE CIRCLE(TYPE,iw,RADIUS,Z,NLINES,ERROR,*)

C#### Subroutine: CIRCLE
C###  Description:
C###    CIRCLE draws polyline (TYPE='POLYLINE') or fill-area
C###    TYPE='FILL AREA') circle on workstation iw with radius RADIUS
C###    about the point Z(nj) with NLINES lines (max 100).

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'trans00.cmn'
!     Parameter List
      INTEGER iw,NLINES
      REAL*8 RADIUS,Z(3)
      CHARACTER ERROR*(*),TYPE*(*)
!     Local Variables
      INTEGER i,nj
      REAL*8 PTS(3,100),THETA

      CALL ENTERS('CIRCLE',*9999)

      IF(DOP) THEN
        WRITE(OP_STRING,'('' Radius='',E12.4,'
     '    //''' Position='',3(E12.4,X))') RADIUS,(Z(nj),nj=1,NJT)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      DO i=1,NLINES+1
        THETA=2.0D0*PI*DBLE(i-1)/DBLE(NLINES)
        IF(NJT.EQ.2) THEN
          PTS(1,i)=Z(1)+RADIUS*DCOS(THETA)
          PTS(2,i)=Z(2)+RADIUS*DSIN(THETA)
          PTS(3,i)=0.0D0
        ELSE IF(NJT.EQ.3) THEN
          IF(iw.EQ.1) THEN
            PTS(1,i)=Z(1)+RADIUS*DCOS(THETA)
            PTS(2,i)=0.0D0
            PTS(3,i)=Z(3)+RADIUS*DSIN(THETA)
          ELSE IF(iw.EQ.2) THEN
            PTS(1,i)=0.0D0
            PTS(2,i)=Z(2)+RADIUS*DCOS(THETA)
            PTS(3,i)=Z(3)+RADIUS*DSIN(THETA)
          ELSE IF(iw.EQ.3) THEN
            CALL ZZ(Z,Z,TRANS)
            PTS(1,i)=Z(1)+RADIUS*DCOS(THETA)
            PTS(2,i)=Z(2)+RADIUS*DSIN(THETA)
            PTS(3,i)=0.0D0
          ENDIF
        ENDIF
      ENDDO

      IF(TYPE(1:8).EQ.'POLYLINE') THEN
        CALL POLYLINE(1,iw,NLINES+1,PTS,ERROR,*9999)
      ELSE IF(TYPE(1:9).EQ.'FILL AREA') THEN
        CALL FILL_AREA(1,iw,NLINES+1,PTS,ERROR,*9999)
      ENDIF

      CALL EXITS('CIRCLE')
      RETURN
 9999 CALL ERRORS('CIRCLE',ERROR)
      CALL EXITS('CIRCLE')
      RETURN 1
      END


