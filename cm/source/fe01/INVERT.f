      SUBROUTINE INVERT(n,A,B,AA)

C#### Subroutine: INVERT
C###  Description:
C###    INVERT returns the inverse of matrix A as B and det(A) as AA.
C###    Matrix A may be no larger than 3*3 (N=3). Note that in
C###    both cases A and B are dimensioned to A(3,3) and B(3,3).

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER n
      REAL*8 A(3,3),AA,B(3,3)
!     Local Variables
      CHARACTER ERROR*10
C      INTEGER i,j,M(5)
      REAL*8 MDET(3)

C MPN 29Jun2000: replaced by XERO_TOL in tol00.cmn
C      PARAMETER(ZERO_TOLERANCE=1.0d-15)

C      DATA M/1,2,3,1,2/

C     CALL ENTERS('INVERT',*9999)

      IF(N.EQ.1) THEN !1*1 matrix (for completeness - do not delete!)
        AA=A(1,1)
        IF(DABS(AA).GT.ZERO_TOL) THEN
          B(1,1)=1.0d0/AA
        ELSE
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' >>Warning: A(1,1) is zero in INVERT'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
C news MPN 29Jun2000: initialise B, coz computations still carry on!
          B(1,1)=0.0d0
C newe
        ENDIF

      ELSE IF(N.EQ.2) THEN !2*2 matrix
        AA=A(1,1)*A(2,2)-A(1,2)*A(2,1)
        IF(DABS(AA).GT.ZERO_TOL) THEN
          B(1,1)= A(2,2)/AA
          B(1,2)=-A(1,2)/AA
          B(2,1)=-A(2,1)/AA
          B(2,2)= A(1,1)/AA
        ELSE
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,
     '      '('' >>Warning: Zero determinant in 2*2 INVERT'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
C news MPN 29Jun2000: initialise B, coz computations still carry on!
          B(1,1)=0.0d0
          B(1,2)=0.0d0
          B(2,1)=0.0d0
          B(2,2)=0.0d0
C newe
        ENDIF

      ELSE IF(N.EQ.3) THEN !3*3 matrix
C        AA=DET(A)
        MDET(1)=A(2,2)*A(3,3)-A(2,3)*A(3,2)
        MDET(2)=A(3,2)*A(1,3)-A(3,3)*A(1,2)
        MDET(3)=A(1,2)*A(2,3)-A(1,3)*A(2,2)
        AA=MDET(1)*A(1,1)+MDET(2)*A(2,1)+MDET(3)*A(3,1)
        IF(DABS(AA).GT.ZERO_TOL) THEN
          B(1,1)=MDET(1)/AA
          B(1,2)=MDET(2)/AA
          B(1,3)=MDET(3)/AA
          B(2,1)=(A(2,3)*A(3,1)-A(3,3)*A(2,1))/AA
          B(2,2)=(A(3,3)*A(1,1)-A(1,3)*A(3,1))/AA
          B(2,3)=(A(1,3)*A(2,1)-A(2,3)*A(1,1))/AA
          B(3,1)=(A(2,1)*A(3,2)-A(3,1)*A(2,2))/AA
          B(3,2)=(A(3,1)*A(1,2)-A(1,1)*A(3,2))/AA
          B(3,3)=(A(1,1)*A(2,2)-A(2,1)*A(1,2))/AA
C          DO i=1,3
C            DO j=1,3
C              B(i,j)=(A(M(j+1),M(i+1))*A(M(j+2),M(i+2))
C     '               -A(M(j+2),M(i+1))*A(M(j+1),M(i+2)))/AA
C            ENDDO !j
C          ENDDO !i
        ELSE
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,
     '      '('' >>Warning: Zero determinant in 3*3 INVERT'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
C news MPN 29Jun2000: initialise B, coz computations still carry on!
          B(1,1)=0.0d0
          B(1,2)=0.0d0
          B(1,3)=0.0d0
          B(2,1)=0.0d0
          B(2,2)=0.0d0
          B(2,3)=0.0d0
          B(3,1)=0.0d0
          B(3,2)=0.0d0
          B(3,3)=0.0d0
C newe
        ENDIF
      ELSE
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,
     '    '('' >>Warning: Matrix larger than 3x3 - cannot INVERT!'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

C     CALL EXITS('INVERT')

 9999 RETURN
      END


