      SUBROUTINE ISORTP(n,IDATA,PIVOT)

C#### Subroutine: ISORTP
C###  Description:
C###    ISORTP sorts N integer IDATA values into a non-decreasing
C###    sequence using a bubble sort algorithm. It also returns
C###    a list of the final pivot vector used to reorder the array.

      IMPLICIT NONE
!     Parameter List
      INTEGER IDATA(*),n,PIVOT(*)
!     Local Variables
      INTEGER FLAG,i,ITEMP,j,k,PTEMP

      DO i=1,n
        PIVOT(i)=i
      ENDDO

      IF(n.LE.1) RETURN
      FLAG=n
      DO 2 i=1,n
        k=FLAG-1
        FLAG=0
        DO 1 j=1,k
          IF(IDATA(j).GT.IDATA(j+1)) THEN
            ITEMP=IDATA(j)
            IDATA(j)=IDATA(j+1)
            IDATA(j+1)=ITEMP

            PTEMP=PIVOT(j)
            PIVOT(j)=PIVOT(j+1)
            PIVOT(j+1)=PTEMP

            FLAG=j
          ENDIF
 1      CONTINUE
        IF(FLAG.EQ.0) THEN
          GO TO 3
        ENDIF
 2    CONTINUE

 3    RETURN
      END


