      SUBROUTINE DXRCDX(ICOORD,NJT,DXRCX,X,ERROR,*)

C#### Subroutine: DXRCDX
C###  Description:
C###    <HTML>
C###    DXRCDX calculates dXRC(i)/dX(j) at X, where XRC(i) are
C###    rectangular Cartesian and X(j) are curvilinear coordinates
C###    defined by ICOORD.
C###    <PRE>
C###    ICOORD=1 for rectangular cartesian coordinates;
C###           2 for cylindrical polar coordinates;
C###           3 for spherical polar coordinates;
C###           4 for prolate spheriodal coordinates; and
C###           5 for oblate spheroidal coordinates.
C###    </PRE></HTML>

      IMPLICIT NONE
      INCLUDE 'b14.cmn'
!     Parameter List
      INTEGER ICOORD,NJT
      REAL*8  DXRCX(3,3),X(3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j
      REAL*8 COSHX1,COSX2,COSX3,SINHX1,SINX2,SINX3

      CALL ENTERS('DXRCDX',*9999)

      GO TO (100,200,300,400,500),ICOORD
        ERROR='>> Unknown coordinate system'
        GO TO 9999
 100  CONTINUE !rectangular cartesian already
        DO j=1,NJT
          DO i=1,NJT
            IF(i.EQ.j) THEN
              DXRCX(i,j)=1.0D0
            ELSE
              DXRCX(i,j)=0.0D0
            ENDIF
          ENDDO
        ENDDO
        GO TO 1000
 200  CONTINUE !cylindrical polar
        COSX2=DCOS(X(2))
        SINX2=DSIN(X(2))
        DXRCX(1,1)=COSX2
        DXRCX(2,1)=SINX2
        DXRCX(1,2)=-X(1)*SINX2
        DXRCX(2,2)=X(1)*COSX2
        IF(NJT.EQ.3) THEN
          DXRCX(3,1)=0.0D0
          DXRCX(3,2)=0.0D0
          DXRCX(1,3)=0.0D0
          DXRCX(2,3)=0.0D0
          DXRCX(3,3)=1.0D0
        ENDIF
        GO TO 1000
 300  CONTINUE !spherical polar
        COSX2=DCOS(X(2))
        SINX2=DSIN(X(2))
        COSX3=DCOS(X(3))
        SINX3=DSIN(X(3))
        DXRCX(1,1)=COSX2*COSX3
        DXRCX(2,1)=SINX2*COSX3
        DXRCX(3,1)=SINX3
        DXRCX(1,2)=-X(1)*SINX2*COSX3
        DXRCX(2,2)=X(1)*COSX2*COSX3
        DXRCX(3,2)=0.0D0
        DXRCX(1,3)=-X(1)*COSX2*SINX3
        DXRCX(2,3)=-X(1)*SINX2*SINX3
        DXRCX(3,3)=X(1)*COSX3
        GO TO 1000
 400  CONTINUE !prolate spheriodal
        COSHX1=DCOSH(X(1))
        SINHX1=DSINH(X(1))
        COSX2=DCOS(X(2))
        SINX2=DSIN(X(2))
        COSX3=DCOS(X(3))
        SINX3=DSIN(X(3))
        DXRCX(1,1)=FOCUS*SINHX1*COSX2
        DXRCX(2,1)=FOCUS*COSHX1*SINX2*COSX3
        DXRCX(3,1)=FOCUS*COSHX1*SINX2*SINX3
        DXRCX(1,2)=-FOCUS*COSHX1*SINX2
        DXRCX(2,2)=FOCUS*SINHX1*COSX2*COSX3
        DXRCX(3,2)=FOCUS*SINHX1*COSX2*SINX3
        DXRCX(1,3)=0.0D0
        DXRCX(2,3)=-FOCUS*SINHX1*SINX2*SINX3
        DXRCX(3,3)=FOCUS*SINHX1*SINX2*COSX3
        GO TO 1000
 500  CONTINUE !oblate spheroidal
        COSHX1=DCOSH(X(1))
        SINHX1=DSINH(X(1))
        COSX2=DCOS(X(2))
        SINX2=DSIN(X(2))
        COSX3=DCOS(X(3))
        SINX3=DSIN(X(3))
        DXRCX(1,1)=FOCUS*SINHX1*COSX2*COSX3
        DXRCX(2,1)=FOCUS*COSHX1*SINX2
        DXRCX(3,1)=FOCUS*SINHX1*COSX2*SINX3
        DXRCX(1,2)=-FOCUS*COSHX1*SINX2*COSX3
        DXRCX(2,2)=FOCUS*SINHX1*COSX2
        DXRCX(3,2)=-FOCUS*COSHX1*SINX2*SINX3
        DXRCX(1,3)=-FOCUS*COSHX1*COSX2**SINX3
        DXRCX(2,3)=0.0D0
        DXRCX(3,3)=FOCUS*COSHX1*COSX2*COSX3

 1000 CALL EXITS('DXRCDX')
      RETURN
 9999 CALL ERRORS('DXRCDX',ERROR)
      CALL EXITS('DXRCDX')
      RETURN 1
      END


