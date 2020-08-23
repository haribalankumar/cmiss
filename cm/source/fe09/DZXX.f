      REAL*8 FUNCTION DZXX(ICOORD,i,j,k,X)

C#### Function: DZXX
C###  Type: REAL*8
C###  Description:
C###    <HTML>
C###    DZXX calculates D2Z(i)/[DX(j).DX(k)] at X, where Z(i) are
C###    rectangular Cartesian and X(k) are curvilinear coordinates
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
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER ICOORD,i,j,k
      REAL*8  X(3)
      CHARACTER ERROR*10

      DZXX=0.0d0
      IF(ICOORD.EQ.1) THEN !rectangular cartesian coords
C       second derivs all zero

      ELSE IF(ICOORD.EQ.2) THEN !cylindrical polar coords
        IF(k.EQ.1) THEN
          IF(j.EQ.2) THEN
            IF(i.EQ.1) THEN
              DZXX=-DSIN(X(2))
            ELSE IF(i.EQ.2) THEN
              DZXX=DCOS(X(2))
            ENDIF
          ENDIF
        ELSE IF(k.EQ.2) THEN
          IF(j.EQ.1) THEN
            IF(i.EQ.1) THEN
              DZXX=-DSIN(X(2))
            ELSE IF(i.EQ.2) THEN
              DZXX=DCOS(X(2))
            ENDIF
          ELSE IF(j.EQ.2) THEN
            IF(i.EQ.1) THEN
              DZXX=-X(1)*DCOS(X(2))
            ELSE IF(i.EQ.2) THEN
              DZXX=-X(1)*DSIN(X(2))
            ENDIF
          ENDIF
        ENDIF

      ELSE IF(ICOORD.EQ.3) THEN !spherical polar coords
        IF(k.EQ.1) THEN
          IF(j.EQ.2) THEN
            IF(i.EQ.1) THEN
              DZXX=-DSIN(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.2) THEN
              DZXX=DCOS(X(2))*DCOS(X(3))
            ENDIF
          ELSE IF(j.EQ.3) THEN
            IF(i.EQ.1) THEN
              DZXX=-DCOS(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.2) THEN
              DZXX=-DSIN(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=DCOS(X(3))
            ENDIF
          ENDIF
        ELSE IF(k.EQ.2) THEN
          IF(j.EQ.1) THEN
            IF(i.EQ.1) THEN
              DZXX=-DSIN(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.2) THEN
              DZXX=DCOS(X(2))*DCOS(X(3))
            ENDIF
          ELSE IF(j.EQ.2) THEN
            IF(i.EQ.1) THEN
              DZXX=-X(1)*DCOS(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.2) THEN
              DZXX=-X(1)*DSIN(X(2))*DCOS(X(3))
            ENDIF
          ELSE IF(j.EQ.3) THEN
            IF(i.EQ.1) THEN
              DZXX=X(1)*DSIN(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.2) THEN
              DZXX=-X(1)*DCOS(X(2))*DSIN(X(3))
            ENDIF
          ENDIF
        ELSE IF(k.EQ.3) THEN
          IF(j.EQ.1) THEN
            IF(i.EQ.1) THEN
              DZXX=-DCOS(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.2) THEN
              DZXX=-DSIN(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=DCOS(X(3))
            ENDIF
          ELSE IF(j.EQ.2) THEN
            IF(i.EQ.1) THEN
              DZXX=X(1)*DSIN(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.2) THEN
              DZXX=-X(1)*DCOS(X(2))*DSIN(X(3))
            ENDIF
          ELSE IF(j.EQ.3) THEN
            IF(i.EQ.1) THEN
              DZXX=-X(1)*DCOS(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.2) THEN
              DZXX=-X(1)*DSIN(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=-X(1)*DSIN(X(3))
            ENDIF
          ENDIF
        ENDIF

      ELSE IF(ICOORD.EQ.4) THEN !prolate spheroidal coords
        IF(k.EQ.1) THEN
          IF(j.EQ.1) THEN
            IF(i.EQ.1) THEN
              DZXX=FOCUS*DCOSH(X(1))*DCOS(X(2))
            ELSE IF(i.EQ.2) THEN
              DZXX=FOCUS*DSINH(X(1))*DSIN(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=FOCUS*DSINH(X(1))*DSIN(X(2))*DSIN(X(3))
            ENDIF
          ELSE IF(j.EQ.2) THEN
            IF(i.EQ.1) THEN
              DZXX=-FOCUS*DSINH(X(1))*DSIN(X(2))
            ELSE IF(i.EQ.2) THEN
              DZXX=FOCUS*DCOSH(X(1))*DCOS(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=FOCUS*DCOSH(X(1))*DCOS(X(2))*DSIN(X(3))
            ENDIF
          ELSE IF(j.EQ.3) THEN
            IF(i.EQ.2) THEN
              DZXX=-FOCUS*DCOSH(X(1))*DSIN(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=FOCUS*DCOSH(X(1))*DSIN(X(2))*DCOS(X(3))
            ENDIF
          ENDIF
        ELSE IF(k.EQ.2) THEN
          IF(j.EQ.1) THEN
            IF(i.EQ.1) THEN
              DZXX=-FOCUS*DSINH(X(1))*DSIN(X(2))
            ELSE IF(i.EQ.2) THEN
              DZXX=FOCUS*DCOSH(X(1))*DCOS(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=FOCUS*DCOSH(X(1))*DCOS(X(2))*DSIN(X(3))
            ENDIF
          ELSE IF(j.EQ.2) THEN
            IF(i.EQ.1) THEN
              DZXX=-FOCUS*DCOSH(X(1))*DCOS(X(2))
            ELSE IF(i.EQ.2) THEN
              DZXX=-FOCUS*DSINH(X(1))*DSIN(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=-FOCUS*DSINH(X(1))*DSIN(X(2))*DSIN(X(3))
            ENDIF
          ELSE IF(j.EQ.3) THEN
            IF(i.EQ.2) THEN
              DZXX=-FOCUS*DSINH(X(1))*DCOS(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=FOCUS*DSINH(X(1))*DCOS(X(2))*DCOS(X(3))
            ENDIF
          ENDIF
        ELSE IF(k.EQ.3) THEN
          IF(j.EQ.1) THEN
            IF(i.EQ.2) THEN
              DZXX=-FOCUS*DCOSH(X(1))*DSIN(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=FOCUS*DCOSH(X(1))*DSIN(X(2))*DCOS(X(3))
            ENDIF
          ELSE IF(j.EQ.2) THEN
            IF(i.EQ.2) THEN
              DZXX=-FOCUS*DSINH(X(1))*DCOS(X(2))*DSIN(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=FOCUS*DSINH(X(1))*DCOS(X(2))*DCOS(X(3))
            ENDIF
          ELSE IF(j.EQ.3) THEN
            IF(i.EQ.2) THEN
              DZXX=-FOCUS*DSINH(X(1))*DSIN(X(2))*DCOS(X(3))
            ELSE IF(i.EQ.3) THEN
              DZXX=-FOCUS*DSINH(X(1))*DSIN(X(2))*DSIN(X(3))
            ENDIF
          ENDIF
        ENDIF

      ELSE IF(ICOORD.EQ.5) THEN !oblate spheroidal coords
        WRITE(OP_STRING,'('' 2nd derivs not implemented for oblate '
     '    //'spheroidal coords'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

 9999 RETURN
      END


