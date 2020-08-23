      SUBROUTINE ZETG44(CE,FX,FY,FYX,FZ,ZE,ERROR,*)

C#### Subroutine: ZETG44
C###  Description:
C###    ZETG44 is for link elements only. Basis function is assumed
C###    to be 1D linear in x and z, and 1D cubic Hermite in y.

C**** CE(1) is stiffness in kN/m in X-direction
C**** CE(2) is stiffness in kN/m in Y-direction
C**** CE(3) is stiffness in kN/m in Z-direction
C**** CE(4) is bending stiffness in x-y plane
C**** FX,FY,FZ are forces generated in X,Y,Z directions in kN.
C**** FYX is moment generated in x-y plane in kNm.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      REAL*8 CE(NMM),FX,FY,FYX,FZ,ZE(NSM,NHM)
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 SX,SY,SYX,SZ

      CALL ENTERS('ZETG44',*9999)
      SX=ZE(2,1)-ZE(1,1)
      SY=ZE(3,2)-ZE(1,2)
      SZ=ZE(2,3)-ZE(1,3)
      SYX=ZE(4,2)-ZE(2,2)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,*) ' SX=',SX,' SY=',SY,' SZ=',SZ,' SYX=',SYX
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      FX=CE(1)*SX
      FY=CE(2)*SY
      FZ=CE(3)*SZ
      FYX=CE(4)*SYX
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,*) ' FX=',FX,' FY=',FY,' FZ=',FZ,' FYX=',FYX
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZETG44')
      RETURN
 9999 CALL ERRORS('ZETG44',ERROR)
      CALL EXITS('ZETG44')
      RETURN 1
      END


