      SUBROUTINE IHEAPSORT(N,RA)

C#### Subroutine: IHEAPSORT
C###  Description:
C###    Numerical recipes integer heap sort algorithm. Use when N > 50

      IMPLICIT NONE
!     Parameter List
      INTEGER N,RA(N)
!     Local Variables
      INTEGER I,IR,J,L
      INTEGER RRA

      IF(N.LT.2) RETURN

      L=N/2+1
      IR=N
 10   CONTINUE
      IF(L.GT.1) THEN
        L=L-1
        RRA=RA(L)
      ELSE
        RRA=RA(IR)
        RA(IR)=RA(1)
        IR=IR-1
        IF(IR.EQ.1) THEN
          RA(1)=RRA
          RETURN
        ENDIF
      ENDIF
      I=L
      J=L+L
 20   IF(J.LE.IR) THEN
        IF(J.LT.IR) THEN
          IF(RA(J).LT.RA(J+1)) J=J+1
        ENDIF
        IF(RRA.LT.RA(J)) THEN
          RA(I)=RA(J)
          I=J
          J=J+J
        ELSE
          J=IR+1
        ENDIF
        GOTO 20
      ENDIF
      RA(I)=RRA
      GOTO 10
      END


