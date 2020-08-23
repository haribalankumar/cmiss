      SUBROUTINE SETTRANSMAT(PHI,THETA,TRANSLATION,TRANSMAT,ERROR,*)

C#### Subroutine: SETTRANSMAT
C###  Description:
C###    SETTRANSMAT sets the transformation matrix for a general
C###    rotation by THETA in 2D or THETA, PHI in 3D + a translation (i.e
C###    a rotation by theta about the z axis and the a rotation by phi
C###    about the x axis).

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      REAL*8 PHI,THETA,TRANSLATION(3),TRANSMAT(3,4)
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('SETTRANSMAT',*9999)

      IF(NJT.EQ.2) THEN
        TRANSMAT(1,1)=DCOS(THETA)
        TRANSMAT(2,1)=-DSIN(THETA)
        TRANSMAT(1,2)=DSIN(THETA)
        TRANSMAT(2,2)=DCOS(THETA)
        TRANSMAT(1,3)=TRANSLATION(1)
        TRANSMAT(2,3)=TRANSLATION(2)
      ELSE IF(NJT.EQ.3) THEN
        TRANSMAT(1,1)=DCOS(THETA)
        TRANSMAT(2,1)=-DSIN(THETA)
        TRANSMAT(3,1)=0.0d0
        TRANSMAT(1,2)=DSIN(THETA)*DCOS(PHI)
        TRANSMAT(2,2)=DCOS(THETA)*DCOS(PHI)
        TRANSMAT(3,2)=-DSIN(PHI)
        TRANSMAT(1,3)=DSIN(THETA)*DSIN(PHI)
        TRANSMAT(2,3)=DCOS(THETA)*DSIN(PHI)
        TRANSMAT(3,3)=DCOS(PHI)
        TRANSMAT(1,4)=TRANSLATION(1)
        TRANSMAT(2,4)=TRANSLATION(2)
        TRANSMAT(3,4)=TRANSLATION(3)
      ENDIF

      CALL EXITS('SETTRANSMAT')
      RETURN
 9999 CALL ERRORS('SETTRANSMAT',ERROR)
      CALL EXITS('SETTRANSMAT')
      RETURN 1
      END


