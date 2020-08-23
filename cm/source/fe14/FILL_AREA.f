      SUBROUTINE FILL_AREA(IBUNDLE,iw,NTPTS,PTS,ERROR,*)

C#### Subroutine: FILL_AREA
C###  Description:
C###    FILL_AREA draws a fill-area on iw with index IBUNDLE.

C**** PTS(1..3,nopts) is a real*8 array containing the 3D coords of
C**** each point.
C**** If the IBUNDLE is 0 the primitive will use the previously
C**** defined fill-area index.
C**** NOTE: If iw is 15 or 16 (postscript) IBUNDLE is reset to be black.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'disp00.cmn'
!     Parameter List
      INTEGER I,IBUNDLE,iw,NTPTS
      REAL*8 PTS(3,40)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nopts,NJ1,NJ2
      REAL R_PTS(NTPTS,2)

      CALL ENTERS('FILL_AREA',*9999)

      IF(NTPTS.EQ.0) GOTO 9998

      IF(DOP) THEN
        WRITE(OP_STRING,'(/'' PTS(1,1..):'',4E12.3)')
     '    (PTS(1,nopts),nopts=1,4)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' PTS(2,1..):'',4E12.3)')
     '    (PTS(2,nopts),nopts=1,4)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' PTS(3,1..):'',4E12.3)')
     '    (PTS(3,nopts),nopts=1,4)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(IW.EQ.1) THEN
        NJ1=1
        NJ2=2
      ELSE IF(IW.EQ.2) THEN
        NJ1=2
        NJ2=3
      ELSE IF(IW.EQ.3) THEN
        NJ1=1
        NJ2=3
      ENDIF

      DO I=1,NTPTS
        R_PTS(I,1)=REAL(PTS(NJ1,I))
        R_PTS(I,2)=REAL(PTS(NJ2,I))
      ENDDO

      IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN

        CALL FILL_AREA_GX(IBUNDLE,NTPTS,R_PTS,ERROR,*9999)

      ENDIF

 9998 CALL EXITS('FILL_AREA')
      RETURN
 9999 CALL ERRORS('FILL_AREA',ERROR)
      CALL EXITS('FILL_AREA')
      RETURN 1
      END


