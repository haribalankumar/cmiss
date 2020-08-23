      REAL*8 FUNCTION IND(A,B,C,COORD)

C#### Function: IND
C###  Description:
CC JMB 27-NOV-2001

      IMPLICIT NONE
!     Parameter List
      REAL*8 A(3),B(3),C(3),COORD(3)

      IND=(B(2)- A(2))*(C(3)-C(2))-(C(2)-B(2))
     '  *(B(3)- A(3))*(B(1)-COORD(1))+(B(3)-A(3))
     '  *(C(1)- B(1))-(C(3)-B(3))*(B(1)-A(1))
     '  *(B(2)- COORD(2))+(B(1)-A(1))*(C(2)-B(2))
     '  -(C(1)- B(1))*(B(2)-A(2))*(B(3)-COORD(3))

      RETURN
      END


