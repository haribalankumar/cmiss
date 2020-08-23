      SUBROUTINE VECTOR(INDEX,iw,CENTRE,DIRECTION,SCALE,ERROR,*)

C#### Subroutine: VECTOR
C###  Description:
C###    VECTOR draws a vector given by the centre position and the
C###    direction at a give scale.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'disp00.cmn'
!     Parameter List
      INTEGER INDEX,iw
      REAL*8 CENTRE(3),DIRECTION(3),SCALE
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER NJ1,NJ2
      REAL R_CENTRE(2),R_DIRECTION(2)

      CALL ENTERS('VECTOR',*9999)

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

      IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
        R_CENTRE(1)=REAL(CENTRE(NJ1))   !is real x
        R_CENTRE(2)=REAL(CENTRE(NJ2))   !is real y
        R_DIRECTION(1)=REAL(DIRECTION(NJ1)*SCALE)   !is real x
        R_DIRECTION(2)=REAL(DIRECTION(NJ2)*SCALE)   !is real y
        
        CALL VECTOR_GX(INDEX,R_CENTRE,R_DIRECTION,ERROR,*9999)
      ENDIF

      CALL EXITS('VECTOR')
      RETURN
 9999 CALL ERRORS('VECTOR',ERROR)
      CALL EXITS('VECTOR')
      RETURN 1
      END


