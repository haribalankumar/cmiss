      LOGICAL FUNCTION FOUND_LIST(FOUND,REF)

C#### Function: FOUND_LIST
C###  Description:
C###    Constructs a search list of the tetrahedral ordering
C###    to find the referneced tetrahedral. Genmesh function.
CC JMB 16-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
!     Parameter List
      INTEGER FOUND(0:N_GM),REF
!     Local Variables
      INTEGER I,NFOUND

      NFOUND=FOUND(0)

      DO I=1,NFOUND
        IF(FOUND(I).EQ.REF)THEN
          FOUND_LIST=.FALSE.
          GOTO 10
        ENDIF !found
      ENDDO !i
      NFOUND=NFOUND+1
      FOUND(NFOUND)=REF
      ! Update
      FOUND(0)=NFOUND
      FOUND_LIST=.TRUE.

   10 RETURN
      END


C      REAL*8 FUNCTION IND(A,B,C,COORD)
C
CC#### Function: IND
CC###  Description:
CC###    Genmesh function.
CCC JMB 27-NOV-2001
C
C      IMPLICIT NONE
C!     Parameter List
C      REAL*8 A(3),B(3),C(3),COORD(3)
C
C      IND=
C     '  ((B(2)-A(2))*(C(3)-B(3))-(C(2)-B(2))*(B(3)-A(3)))*
C     '  (B(1)-COORD(1))+
C     '  ((B(3)-A(3))*(C(1)-B(1))-(C(3)-B(3))*(B(1)-A(1)))*
C     '  (B(2)-COORD(2))+
C     '  ((B(1)-A(1))*(C(2)-B(2))-(C(1)-B(1))*(B(2)-A(2)))*
C     '  (B(3)-COORD(3))
C
C      RETURN
C      END

