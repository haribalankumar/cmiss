      REAL*8 FUNCTION M_PLGNDR_SIN(L,M,X)

C#### Function: M_PLGNDR_SIN
C###  Type: REAL*8
C###  Description:
C###    PLGNDR computes the associated Legendre
C###    polynomial m*Pl,m(x) (l is the subscript
C###    and m is the superscript) divided by sin theta,
C###    where x=cos theta.  Here l and m are
C###    integers satisfying 1 <= M <= L while X lies in the range
C###    -1 <= X <= 1 (see numerical recipes, 2nd edition p247).
C###  See-Also: CHEBYSHEV2,PLGNDR

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter list
      INTEGER L,M
      REAL*8 X
!     Local variables
      INTEGER i,ll
      REAL*8 FACT,PLL,PMM,PMMP1,SOMX2
      CHARACTER ERROR*10

      IF(M.LT.0.OR.ABS(X).GT.1) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' >>ERROR: Invalid arguments in '
     '    //'M_PLGNDR_SIN'')')
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
        !How does one write out error messages correctly in functions ??
      ELSE
        IF(M.GT.L) THEN
          M_PLGNDR_SIN=0.0d0
        ELSE
          IF(M.EQ.0) THEN
            M_PLGNDR_SIN=0.0d0
          ELSE
            PMM=1.0d0 !Compute Pm,m
            SOMX2=DSQRT((1.0d0-X)*(1.0d0+X))
            FACT=1.0d0
C       We must remove the '*SOMX2' from one of the multiplications
            PMM=-PMM*FACT
            IF(M.GT.1) THEN
              FACT=FACT+2.0d0
              DO i=1,M
                PMM=-PMM*FACT*SOMX2
                FACT=FACT+2.0d0
              ENDDO
            ENDIF
            IF(L.EQ.M) THEN
              M_PLGNDR_SIN=DBLE(M)*PMM
            ELSE
              PMMP1=X*(2*M+1)*PMM !Compute Pm+1,m
              IF(L.EQ.M+1) THEN
                M_PLGNDR_SIN=DBLE(M)*PMMP1
              ELSE !Compute Pl,m l>m+1
                DO ll=M+2,L
                  PLL=(X*(2*ll-1)*PMMP1-(ll+M-1)*PMM)/(ll-M)
                  PMM=PMMP1
                  PMMP1=PLL
                ENDDO
                M_PLGNDR_SIN=DBLE(M)*PLL
              ENDIF
            ENDIF
          ENDIF !m=0
        ENDIF !M.GT.L
      ENDIF
 9999 RETURN
      END


