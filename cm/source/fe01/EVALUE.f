      SUBROUTINE EVALUE(NDIM,A,EVAL,ERROR,*)

C#### Subroutine: EVALUE
C###  Description:
C###    EVALUE solves quadratic(NDIM=2) or cubic(NDIM=3) characteristic
C###    polynomial for symmetric matrix A. The eigenvalues are returned
C###    in EVAL(j),j=1,NDIM.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER NDIM
      REAL*8 A(3,3),EVAL(NDIM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i
      REAL*8 ANG,B2,B3,C1,C2,D,DET,Q,Q3,R,RI1,RI2,RI3,RI4,RQ,TEMP,
     '  THETA

      REAL*8 ZERO_TOLERANCE
      PARAMETER(ZERO_TOLERANCE=1.0d-10)

      CALL ENTERS('EVALUE',*9999)

c cpb 1/6/96 Ordering the eigenvalues for largest to smallest absolute
c value

      IF(NDIM.EQ.1) THEN
        EVAL(1)=A(1,1)
      ELSE IF(NDIM.EQ.2) THEN
        IF(DABS(A(1,2)).GT.ZERO_TOLERANCE) THEN
          RI1=A(1,1)+A(2,2)
          RI2=A(1,1)*A(2,2)-A(1,2)**2
          B2=RI1/2.0d0
          C1=RI1*RI1
          C2=4.0D0*RI2
          IF(C2.GT.C1) THEN
            ERROR='>>Complex roots found in quadratic equation'
            GOTO 9999
          ENDIF
          B3=DSQRT(C1-C2)/2.0d0
          EVAL(1)=B2+B3
          EVAL(2)=B2-B3
        ELSE
          EVAL(1)=A(1,1)
          EVAL(2)=A(2,2)
        ENDIF
        IF(DABS(EVAL(2)).GT.DABS(EVAL(1))) THEN
          TEMP=EVAL(1)
          EVAL(1)=EVAL(2)
          EVAL(2)=TEMP
        ENDIF
      ELSE IF(NDIM.EQ.3) THEN
        RI1=A(1,1)+A(2,2)+A(3,3)
        RI2=A(1,1)*A(2,2)+A(2,2)*A(3,3)+A(3,3)*A(1,1)
     '    -(A(1,2)**2+A(2,3)**2+A(3,1)**2)
        RI3=DET(A)
        RI4=RI1/3.0d0
        Q=RI4*RI4-RI2/3.0d0
        R=RI4*(RI4*RI4-RI2/2.0D0)+RI3/2.0d0
        Q3=Q*Q*Q
        D=R*R-Q3
        IF(DABS(D).GT.ZERO_TOLERANCE) THEN
          ERROR='>>Complex roots found in solution of cubic equation'
          GOTO 9999
        ENDIF
        RQ=DSQRT(DABS(Q))
C GMH.  1/8/95 We know that Q3 is greater than or equal to zero.
C       This means we can cut out the test for equal eivalues,
C       and just use the general case.  If r=0 then the roots will be
C       equal.
C       If Q3 is equal to zero then r must equal zero, so
C       according to the formula, theta = 0
        IF(DABS(Q).LT.ZERO_TOLERANCE) THEN
          THETA=0.0d0
        ELSE
          THETA=DACOS(R/DSQRT(DABS(Q3)))/3.0d0
        ENDIF
        ANG=2.0d0*PI/3.0d0
        EVAL(1)=2.0d0*RQ*DCOS(THETA)+RI4
        EVAL(2)=2.0d0*RQ*DCOS(THETA+ANG)+RI4
        EVAL(3)=2.0d0*RQ*DCOS(THETA+2.0d0*ANG)+RI4
C        IF(D.LT.(-1.0d0*ZER_TOLERANCE)) THEN
C          THETA=DACOS(R/DSQRT(DABS(Q3)))/3.0d0
C          ANG=2.0d0*PI/3.0d0
C          DO i=1,3
C            EVAL(i)=2.0d0*RQ*DCOS(THETA+(i-1)*ANG)+RI4
C          ENDDO
C        ELSE
C          DIFF=(RI4+2.0d0*RQ)*(RI4-RQ)**2-RI3
C          IF(DABS(DIFF).GT.ZERO_TOLERANCE) RQ=-1.0d0*RQ
C          EVAL(1)=RI4+2.0d0*RQ
C          EVAL(2)=RI4-RQ
C          EVAL(3)=EVAL(2)
C        ENDIF
        DO i=1,2
          IF(DABS(EVAL(3)).GT.DABS(EVAL(i))) THEN
            TEMP=EVAL(i)
            EVAL(i)=EVAL(3)
            EVAL(3)=TEMP
          ENDIF
        ENDDO !i
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' EVAL: '',3D12.4)') (EVAL(i),i=1,NDIM)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('EVALUE')
      RETURN
 9999 CALL ERRORS('EVALUE',ERROR)
      CALL EXITS('EVALUE')
      RETURN 1
      END


