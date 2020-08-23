      SUBROUTINE DIPOLE_SOLVE(CHOICE_ANAL,DIPOLE_POSITION,
     '  DIPOLE_ORIENTATION,ERROR,*)

C#### Subroutine: DIPOLE_SOLVE
C###  Description:
C###    DIPOLE_SAVE calculates the coefficients of the Legendre
C###    polynomials in the series expansion of the solutions for a
C###    dipole inside nspheres concentric spheres.

      IMPLICIT NONE
      INCLUDE 'anal00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'mesh03.cmn'
!     Parameter list
      INTEGER CHOICE_ANAL
      REAL*8 DIPOLE_POSITION(3),DIPOLE_ORIENTATION(3)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER i,ICOL,IFAIL,IPIV(20),
     '  j,NSERIES,NSERIES_MAX,NUM_VARS,nospheres,
     '  MSERIES,MSERIES_MAX,row
      REAL*8 A(20,20),B(20)

      CALL ENTERS('DIPOLE_SOLVE',*9999)

C***  Set up matrix equation for coefficients
      CALL ASSERT(NSPHERES.LE.10,' >>A matrix not large enough',
     '  ERROR,*9999)
      IF(NJT.EQ.3.AND.CHOICE_ANAL.EQ.15) THEN !eccentric,multiple
        NSERIES_MAX=ANAL_N
        MSERIES_MAX=1
      ELSE
        NSERIES_MAX=1
        MSERIES_MAX=0
      ENDIF
      DO NSERIES=1,NSERIES_MAX
        DO MSERIES=0,MSERIES_MAX
          DO i=1,20
            B(i)=0.0d0
            DO j=1,20
              A(i,j)=0.0d0
            ENDDO !j
          ENDDO !i

C cpb 26/1/96 Adding 2D dipole solution
          IF(NJT.EQ.2) THEN
            A(1,1)=MESH3_RAD(1)
            A(1,2)=-MESH3_RAD(1)
            A(1,3)=-1.0d0/MESH3_RAD(1)
            A(2,1)=SIGMA(1)
            A(2,2)=-SIGMA(2)
            A(2,3)=SIGMA(2)/(MESH3_RAD(1)**2)
            B(2)=(SIGMA(1)-SIGMA(2))/(SIGMA(1)*2.0d0*PI*MESH3_RAD(1)**2)
            ICOL=2
            DO nospheres=2,NSPHERES-1
              A(2*nospheres-1,ICOL)=MESH3_RAD(nospheres)
              A(2*nospheres-1,ICOL+1)=1.0d0/MESH3_RAD(nospheres)
              A(2*nospheres-1,ICOL+2)=-MESH3_RAD(nospheres)
              A(2*nospheres-1,ICOL+3)=-1.0d0/MESH3_RAD(nospheres)
              A(2*nospheres,ICOL)=SIGMA(nospheres)
              A(2*nospheres,ICOL+1)=-SIGMA(nospheres)/
     '          (MESH3_RAD(nospheres)**2)
              A(2*nospheres,ICOL+2)=-SIGMA(nospheres+1)
              A(2*nospheres,ICOL+3)=SIGMA(nospheres+1)/
     '          (MESH3_RAD(nospheres)**2)
              B(2*nospheres)=(SIGMA(nospheres)-SIGMA(nospheres+1))/
     '          (SIGMA(1)*2.0d0*PI*MESH3_RAD(nospheres)**2)
              ICOL=ICOL+2
            ENDDO !nospheres
            A(2*NSPHERES-1,2*NSPHERES-2)=MESH3_RAD(NSPHERES)
            A(2*NSPHERES-1,2*NSPHERES-1)=1.0d0/MESH3_RAD(NSPHERES)
            A(2*NSPHERES-1,2*NSPHERES)=-1.0d0/MESH3_RAD(NSPHERES)
            A(2*NSPHERES,2*NSPHERES-2)=SIGMA(NSPHERES)
            A(2*NSPHERES,2*NSPHERES-1)=-SIGMA(NSPHERES)/
     '        (MESH3_RAD(NSPHERES)**2)
            B(2*NSPHERES)=SIGMA(NSPHERES)/
     '        (SIGMA(1)*2.0d0*PI*MESH3_RAD(NSPHERES)**2)
            NUM_VARS=2*NSPHERES
          ELSE IF(NJT.EQ.3) THEN
            IF(CHOICE_ANAL.NE.15) THEN !not eccentric
              A(1,1)=MESH3_RAD(1)
              A(1,2)=-MESH3_RAD(1)
              A(1,3)=-1.0d0/(MESH3_RAD(1)*MESH3_RAD(1))
              A(2,1)=SIGMA(1)
              A(2,2)=-SIGMA(2)
              A(2,3)=2.0d0*SIGMA(2)/(MESH3_RAD(1)**3)
              B(2)=(SIGMA(2)-SIGMA(1))/(SIGMA(1)*2.0d0*PI*
     '          MESH3_RAD(1)**3)
              ICOL=2
              DO nospheres=2,NSPHERES-1
                A(2*nospheres-1,ICOL)=MESH3_RAD(nospheres)
                A(2*nospheres-1,ICOL+1)=1.0d0/(MESH3_RAD(nospheres)**2)
                A(2*nospheres-1,ICOL+2)=-MESH3_RAD(nospheres)
                A(2*nospheres-1,ICOL+3)=-1.0d0/(MESH3_RAD(nospheres)**2)
                A(2*nospheres,ICOL)=SIGMA(nospheres)
                A(2*nospheres,ICOL+1)=-2.0d0*SIGMA(nospheres)/
     '            (MESH3_RAD(nospheres)**3)
                A(2*nospheres,ICOL+2)=-SIGMA(nospheres+1)
                A(2*nospheres,ICOL+3)=2.0d0*SIGMA(nospheres+1)/
     '            (MESH3_RAD(nospheres)**3)
                B(2*nospheres)=(SIGMA(nospheres+1)-SIGMA(nospheres))/
     '            (SIGMA(1)*2.0d0*PI*MESH3_RAD(nospheres)**3)
                ICOL=ICOL+2
              ENDDO !nospheres
              A(2*NSPHERES-1,2*NSPHERES-2)=MESH3_RAD(NSPHERES)
              A(2*NSPHERES-1,2*NSPHERES-1)=1.0d0/
     '          (MESH3_RAD(NSPHERES)**2)
              A(2*NSPHERES-1,2*NSPHERES)=-1.0d0/(MESH3_RAD(NSPHERES)**2)
              A(2*NSPHERES,2*NSPHERES-2)=SIGMA(NSPHERES)
              A(2*NSPHERES,2*NSPHERES-1)=-2.0d0*SIGMA(NSPHERES)/
     '          (MESH3_RAD(NSPHERES)**3)
              B(2*NSPHERES)=-SIGMA(NSPHERES)/
     '          (SIGMA(1)*2.0d0*PI*MESH3_RAD(NSPHERES)**3)
              NUM_VARS=2*NSPHERES
            ELSE !eccentric
C              DIPOLE_F=DIPOLE_POSITION(3) !z distance
C              CALL ASSERT(DIPOLE_F.NE.0.0D0,
C     '          '>>Cannot be a centric dipole',ERROR,*9999)
C              CALL ASSERT(DIPOLE_POSITION(1).EQ.0.0D0,
C     '          '>>Dipole must be on the z axis',ERROR,*9999)
C              CALL ASSERT(DIPOLE_POSITION(2).EQ.0.0D0,
C     '          '>>Dipole must be on the z axis',ERROR,*9999)
C              CALL ASSERT(DIPOLE_ORIENTATION(2).EQ.0.0D0,
C     '          '>>Must be an x dipole',ERROR,*9999)
C              CALL ASSERT(DIPOLE_ORIENTATION(3).EQ.0.0D0,
C     '          '>>Must be on the z axis',ERROR,*9999)
CC           Inner sphere - dipole and r^n term
C              A(1,1)=MESH3_RAD(1)**NSERIES
C              A(1,2)=-1.0d0/(MESH3_RAD(1)**(NSERIES+1))
C              A(1,3)=-MESH3_RAD(1)**NSERIES
C              B(1)  =-(DIPOLE_F**(NSERIES-1))/(MESH3_RAD(1)**(NSERIES+1))
C              A(2,1)=SIGMA(1)*DBLE(NSERIES)*(MESH3_RAD(1)**(NSERIES-1))
C              A(2,2)=SIGMA(2)*DBLE(NSERIES+1)/(MESH3_RAD(1)**(NSERIES+2))
C              A(2,3)=-SIGMA(2)*DBLE(NSERIES)*(MESH3_RAD(1)**(NSERIES-1))
C              B(2)  =SIGMA(1)*DBLE(NSERIES+1)*(DIPOLE_F**(NSERIES-1))/
C     '          (MESH3_RAD(1)**(NSERIES+2))
C              ICOL=2
C              DO nospheres=2,NSPHERES-1
C                A(2*nospheres-1,ICOL)  =1.0d0/(MESH3_RAD(nospheres)**
C     '            (NSERIES+1))
C                A(2*nospheres-1,ICOL+1)=MESH3_RAD(nospheres)**NSERIES
C                A(2*nospheres-1,ICOL+2)=-1.0d0/(MESH3_RAD(nospheres)**
C     '          (NSERIES+1))
C                A(2*nospheres-1,ICOL+3)=-MESH3_RAD(nospheres)**NSERIES
C                A(2*nospheres,ICOL)    =-SIGMA(nospheres)*DBLE(NSERIES+1)/
C     '            (MESH3_RAD(nospheres)**(NSERIES+2))
C                A(2*nospheres,ICOL+1)  =SIGMA(nospheres)*DBLE(NSERIES)*
C     '            (MESH3_RAD(nospheres)**(NSERIES-1))
C                A(2*nospheres,ICOL+2) =SIGMA(nospheres+1)*DBLE(NSERIES+1)/
C     '            (MESH3_RAD(nospheres)**(NSERIES+2))
C                A(2*nospheres,ICOL+3)  =-SIGMA(nospheres+1)*DBLE(NSERIES)*
C     '            (MESH3_RAD(nospheres)**(NSERIES-1))
C                ICOL=ICOL+2
C              ENDDO !nospheres
C              A(2*NSPHERES-1,ICOL)=-SIGMA(NSPHERES)*DBLE(NSERIES+1)/
C     '          (MESH3_RAD(NSPHERES)**(NSERIES+2))
C              A(2*NSPHERES-1,ICOL+1)=SIGMA(NSPHERES)*DBLE(NSERIES)*
C     '          (MESH3_RAD(NSPHERES)**(NSERIES-1))
C              NUM_VARS=2*NSPHERES-1
C             columns are B1,A2,B2,A3,B3 etc
              IF(DIPOLE_POSITION(1).NE.0.0d0) THEN
                CALL ASSERT(DIPOLE_POSITION(2).EQ.0.0D0,
     '            '>>Dipole must be on the x axis',ERROR,*9999)
                CALL ASSERT(DIPOLE_POSITION(3).EQ.0.0D0,
     '            '>>Dipole must be on the x axis',ERROR,*9999)
                DIPOLE_F=DIPOLE_POSITION(1) !x distance
                DIPOLE_AXIS=1
                DIPOLE_STRENGTH(1)=-DIPOLE_ORIENTATION(3)
                DIPOLE_STRENGTH(2)=DIPOLE_ORIENTATION(2)
                DIPOLE_STRENGTH(3)=DIPOLE_ORIENTATION(1)
              ELSE IF (DIPOLE_POSITION(2).NE.0.0d0) THEN
                CALL ASSERT(DIPOLE_POSITION(3).EQ.0.0D0,
     '            '>>Dipole must be on the y axis',ERROR,*9999)
                CALL ASSERT(DIPOLE_POSITION(1).EQ.0.0D0,
     '            '>>Dipole must be on the y axis',ERROR,*9999)
                DIPOLE_F=DIPOLE_POSITION(2) !y distance
                DIPOLE_AXIS=2
                DIPOLE_STRENGTH(1)=DIPOLE_ORIENTATION(1)
                DIPOLE_STRENGTH(2)=-DIPOLE_ORIENTATION(3)
                DIPOLE_STRENGTH(3)=DIPOLE_ORIENTATION(2)
              ELSE IF (DIPOLE_POSITION(3).NE.0.0d0) THEN
                CALL ASSERT(DIPOLE_POSITION(1).EQ.0.0D0,
     '            '>>Dipole must be on the z axis',ERROR,*9999)
                CALL ASSERT(DIPOLE_POSITION(2).EQ.0.0D0,
     '            '>>Dipole must be on the z axis',ERROR,*9999)
                DIPOLE_F=DIPOLE_POSITION(3) !z distance
                DIPOLE_AXIS=3
                DIPOLE_STRENGTH(1)=DIPOLE_ORIENTATION(1)
                DIPOLE_STRENGTH(2)=DIPOLE_ORIENTATION(2)
                DIPOLE_STRENGTH(3)=DIPOLE_ORIENTATION(3)
              ELSE
                DIPOLE_F=0.0d0
                DIPOLE_AXIS=3 !Z axis is the easiest
                DIPOLE_STRENGTH(1)=DIPOLE_ORIENTATION(1)
                DIPOLE_STRENGTH(2)=DIPOLE_ORIENTATION(2)
                DIPOLE_STRENGTH(3)=DIPOLE_ORIENTATION(3)
              ENDIF
c              CALL ASSERT(DIPOLE_F.NE.0.0D0,
c     '          '>>Cannot be a centric dipole',ERROR,*9999)
C              CALL ASSERT(DIPOLE_ORIENTATION(2).EQ.0.0D0,
C     '          '>>Must be an x dipole',ERROR,*9999)
C              CALL ASSERT(DIPOLE_ORIENTATION(3).EQ.0.0D0,
C     '          '>>Must be on the z axis',ERROR,*9999)
C              We need a series from n=0 to infinity
              row=1
C             Do the potential terms
              DO nospheres=1,NSPHERES-1
                IF(nospheres.NE.1) THEN
                  A(row,(nospheres*2)-2)=1.0d0/
     '              (MESH3_RAD(nospheres)**(NSERIES+1))
                ENDIF
                A(row,(nospheres*2)-1)=MESH3_RAD(nospheres)**NSERIES
                A(row,(nospheres*2)-0)=-1.0d0/
     '            (MESH3_RAD(nospheres)**(NSERIES+1))
                A(row,(nospheres*2)+1)=-MESH3_RAD(nospheres)**NSERIES
C               Evaluate the RHS
                IF(nospheres.EQ.1) THEN
                  B(row)=-(DIPOLE_F**(NSERIES-1))/
     '              (MESH3_RAD(nospheres)**(NSERIES+1))
                ENDIF
                row=row+1
              ENDDO
C             Do the derivative terms
              DO nospheres=1,NSPHERES
                IF(nospheres.NE.1) THEN
                  A(row,(nospheres*2)-2)=-SIGMA(nospheres)*
     '              DBLE(NSERIES+1)/(MESH3_RAD(nospheres)**(NSERIES+2))
                ENDIF
                A(row,(nospheres*2)-1)=SIGMA(nospheres)*DBLE(NSERIES)*
     '            (MESH3_RAD(nospheres)**(NSERIES-1))
                IF(nospheres.NE.NSPHERES) THEN
                  A(row,(nospheres*2)+0)=SIGMA(nospheres+1)*
     '              DBLE(NSERIES+1)/
     '              (MESH3_RAD(nospheres)**(NSERIES+2))
                  A(row,(nospheres*2)+1)=-SIGMA(nospheres+1)*
     '              DBLE(NSERIES)*
     '              (MESH3_RAD(nospheres)**(NSERIES-1))
                ENDIF
C               Evaluate the RHS
                IF(nospheres.EQ.1) THEN
                  B(row)=(SIGMA(1)*DBLE(NSERIES+1)*
     '              DIPOLE_F**(NSERIES-1))/
     '              (MESH3_RAD(nospheres)**(NSERIES+2))
                ENDIF
                row=row+1
              ENDDO
              NUM_VARS=2*NSPHERES-1
            ENDIF !ANAL_CHOICE.NE.15
          ENDIF !NJT=2/3

C***  Solve Ax=B for the series coefficients DIPOLE_COEFF
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(DIPOLE_SOLVE_1)
            WRITE(OP_STRING,'('' Matrix of '
     '        //'coefficents for a dipole:'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            FORMAT='('' A('',I2,'',... )= '',/(25X,8D12.4))'
            DO nospheres=1,NUM_VARS
              WRITE(OP_STRING,FORMAT) nospheres,
     '          (A(nospheres,J),J=1,NUM_VARS)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDDO
            FORMAT='('' B('',I4,'')= '',D12.6)'
            DO nospheres=1,NUM_VARS
              WRITE(OP_STRING,FORMAT) nospheres,B(nospheres)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDDO !nospheres
CC$OMP END CRITICAL(DIPOLE_SOLVE_1)
          ENDIF
          IFAIL=1
          CALL DGETRF(NUM_VARS,NUM_VARS,A,20,IPIV,IFAIL)
          IF(IFAIL.EQ.0) THEN
            CALL DGETRS('N',NUM_VARS,1,A,20,IPIV,B,20,IFAIL)
            IF(IFAIL.NE.0) THEN
              WRITE(ERROR,'('' >>INFO='',I5,'' in DGETRS'')') IFAIL
              GOTO 9999
            ENDIF
          ELSE
            WRITE(ERROR,'('' >>INFO='',I5,'' in DGETRS'')') IFAIL
            GOTO 9999
          ENDIF
          DO i=1,NUM_VARS
            DIPOLE_COEFF(i,NSERIES*(MSERIES_MAX+1)+MSERIES)=B(i)
          ENDDO
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(DIPOLE_SOLVE_2)
            WRITE(OP_STRING,'('' Solution coefficients'
     '        //'for a dipole (series '',I4,''):'')') NSERIES
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            FORMAT='('' DIPOLE_COEFF ('',I4,'')= '',D12.6)'
            DO nospheres=1,NUM_VARS
              WRITE(OP_STRING,FORMAT) nospheres,
     '          DIPOLE_COEFF(nospheres,NSERIES*2+MSERIES)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDDO !nospheres
CC$OMP END CRITICAL(DIPOLE_SOLVE_2)
          ENDIF
        ENDDO !MSERIES
      ENDDO !NSERIES

      CALL EXITS('DIPOLE_SOLVE')
      RETURN
 9999 CALL ERRORS('DIPOLE_SOLVE',ERROR)
      CALL EXITS('DIPOLE_SOLVE')
      RETURN 1
      END


