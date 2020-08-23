      REAL*8 FUNCTION PLGNDR(L,M,X)

C#### Function: PLGNDR
C###  Type: REAL*8
C###  Description:
C###    PLGNDR computes the associated Legendre polynomial Pl,m(x) (l is
C###    the subscript and m is the superscript).  Here l and m are
C###    integers satisfying 0 <= M <= L while X lies in the range
C###    -1 <= X <= 1 (see numerical recipes, 2nd edition p247).
C###  See-Also: CHEBYSHEV2,M_PLGNDR_SIN

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
        WRITE(OP_STRING,'('' >>ERROR: Invalid arguments in PLGNDR'')')
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
        !How does one write out error messages correctly in functions ??
      ELSE
        IF(M.GT.L) THEN
          PLGNDR=0.0d0
        ELSE
          PMM=1.0d0 !Compute Pm,m
          IF(M.GT.0) THEN
            SOMX2=DSQRT((1.0d0-X)*(1.0d0+X))
            FACT=1.0d0
            DO i=1,M
              PMM=-PMM*FACT*SOMX2
              FACT=FACT+2.0d0
            ENDDO
          ENDIF
          IF(L.EQ.M) THEN
            PLGNDR=PMM
          ELSE
            PMMP1=X*(2*M+1)*PMM !Compute Pm+1,m
            IF(L.EQ.M+1) THEN
              PLGNDR=PMMP1
            ELSE !Compute Pl,m l>m+1
              DO ll=M+2,L
                PLL=(X*(2*ll-1)*PMMP1-(ll+M-1)*PMM)/(ll-M)
                PMM=PMMP1
                PMMP1=PLL
              ENDDO
              PLGNDR=PLL
            ENDIF
          ENDIF
        ENDIF !M.GT.L
      ENDIF
 9999 RETURN
      END

