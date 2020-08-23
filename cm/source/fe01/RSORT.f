      SUBROUTINE RSORT(n,RDATA,INDEX)

C#### Subroutine: RSORT
C###  Description:
C###    RSORT sorts N DATA values into a non-decreasing sequence
C###    using a bubble sort algorithm.
C**** INDEX(N) is rearranged in the same order as DATA(N).

      IMPLICIT NONE
!     Parameter List
      INTEGER INDEX(*),n
      REAL*8 RDATA(*)
!     Local Variables
      INTEGER FLAG,i,ITEMP,j,k
      REAL*8 TEMP

C     CALL ENTERS('RSORT',*9999)
      IF(N.LE.1) THEN
        GO TO 9999
      ENDIF
      FLAG=n
      DO 2 i=1,n
        k=FLAG-1
        FLAG=0
        DO 1 j=1,k
          IF(RDATA(j).GT.RDATA(j+1)) THEN
            TEMP=RDATA(j)
            RDATA(j)=RDATA(j+1)
            RDATA(j+1)=TEMP
            ITEMP=INDEX(j)
            INDEX(j)=INDEX(j+1)
            INDEX(j+1)=ITEMP
            FLAG=j
          ENDIF
 1      CONTINUE
        IF(FLAG.EQ.0) THEN
          GO TO 9999
        ENDIF
 2    CONTINUE

C9999 CALL EXITS('RSORT')
 9999 RETURN
      END


