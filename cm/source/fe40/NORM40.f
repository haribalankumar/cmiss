      SUBROUTINE NORM40(XG,XN,ERROR,*)

C#### Subroutine: NORM40
C###  Description:
C###    NORM40 finds the normal vector XN(nj) in the cartesian
C###    reference frame. Note: XG(i,j) must be turned around.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      REAL*8 XG(NJM,NUM),XN(NJM)
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 XGT1(3),XGT2(3)

      CALL ENTERS('NORM40',*9999)
      IF(ITYP10(1).EQ.1) THEN
        XGT1(1)=XG(1,2)
        XGT1(2)=XG(2,2)
        XGT1(3)=XG(3,2)
        XGT2(1)=XG(1,4)
        XGT2(2)=XG(2,4)
        XGT2(3)=XG(3,4)
      ELSE IF(ITYP10(1).EQ.2) THEN
        XGT1(1)=XG(1,2)*DCOS(XG(2,1))-XG(1,1)*DSIN(XG(2,1))*XG(2,2)
        XGT1(2)=XG(1,2)*DSIN(XG(2,1))+XG(1,1)*DCOS(XG(2,1))*XG(2,2)
        XGT1(3)=XG(3,2)
        XGT2(1)=XG(1,4)*DCOS(XG(2,1))-XG(1,1)*DSIN(XG(2,1))*XG(2,4)
        XGT2(2)=XG(1,4)*DSIN(XG(2,1))+XG(1,1)*DCOS(XG(2,1))*XG(2,4)
        XGT2(3)=XG(3,4)
      ELSE IF(ITYP10(1).EQ.3) THEN
        XGT1(1)=XG(1,2)*DCOS(XG(3,1))*DCOS(XG(2,1))
     '         -XG(1,1)*DSIN(XG(3,1))*DCOS(XG(2,1))*XG(3,2)
     '         -XG(1,1)*DCOS(XG(3,1))*DSIN(XG(2,1))*XG(2,2)
        XGT1(2)=XG(1,2)*DCOS(XG(3,1))*DSIN(XG(2,1))
     '         -XG(1,1)*DSIN(XG(3,1))*DSIN(XG(2,1))*XG(3,2)
     '         +XG(1,1)*DCOS(XG(3,1))*DCOS(XG(2,1))*XG(2,2)
        XGT1(3)=XG(1,2)*DSIN(XG(3,1))+XG(1,1)*DCOS(XG(3,1))*XG(3,2)
        XGT2(1)=XG(1,4)*DCOS(XG(3,1))*DCOS(XG(2,1))
     '         -XG(1,1)*DSIN(XG(3,1))*DCOS(XG(2,1))*XG(3,4)
     '         -XG(1,1)*DCOS(XG(3,1))*DSIN(XG(2,1))*XG(2,4)
        XGT2(2)=XG(1,4)*DCOS(XG(3,1))*DSIN(XG(2,1))
     '         -XG(1,1)*DSIN(XG(3,1))*DSIN(XG(2,1))*XG(3,4)
     '         +XG(1,1)*DCOS(XG(3,1))*DCOS(XG(2,1))*XG(2,4)
        XGT2(3)=XG(1,4)*DSIN(XG(3,1))+XG(1,1)*DCOS(XG(3,1))*XG(3,4)
      ENDIF
C
C  Finds normal components XN(nj)
C

C     nr = 1 in NJ_LOC is temporary RGB - replacing NJE
      IF(NJ_LOC(NJL_GEOM,0,1).EQ.2) THEN
        XN(1)=XGT1(2)
        XN(2)=-XGT1(1)

      ELSE IF(NJ_LOC(NJL_GEOM,0,1).EQ.3) THEN
        XN(1)=XGT1(2)*XGT2(3)-XGT1(3)*XGT2(2)
        XN(2)=XGT1(3)*XGT2(1)-XGT1(1)*XGT2(3)
        XN(3)=XGT1(1)*XGT2(2)-XGT1(2)*XGT2(1)
      ENDIF

      CALL EXITS('NORM40')
      RETURN
 9999 CALL ERRORS('NORM40',ERROR)
      CALL EXITS('NORM40')
      RETURN 1
      END


