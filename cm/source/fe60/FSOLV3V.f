      SUBROUTINE FSOLV3V(A,B,ERROR,*)

C#### Subroutine FSOLV3V
C###  Description:
C###    FSOLV3V solves a 3x3 system of equations quickly

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      REAL*8 A(3,3),B(3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i
      REAL*8 PIVOT1,PIVOT2,PIVOT3,PIVMAX,TEMP1,TEMP2,TEMP3

      CALL ENTERS('FSOLV3V',*9999)

      PIVMAX=MAX(ABS(A(1,1)),ABS(A(2,1)),ABS(A(3,1)))
      CALL ASSERT(PIVMAX.NE.(0.d0),'Singular matrix',ERROR,*9999)
      IF(PIVMAX.EQ.ABS(A(2,1)))THEN !swap rows 1 and 2
        DO i=1,3
          TEMP1=A(1,i)
          A(1,i)=A(2,i)
          A(2,i)=TEMP1
        ENDDO !i
        TEMP1=B(1)
        B(1)=B(2)
        B(2)=TEMP1
      ELSE IF(PIVMAX.EQ.ABS(A(3,1)))THEN !swap rows 1 and 3
        DO i=1,3
          TEMP1=A(1,i)
          A(1,i)=A(3,i)
          A(3,i)=TEMP1
        ENDDO !i
        TEMP1=B(1)
        B(1)=B(3)
        B(3)=TEMP1
      ENDIF
      PIVOT1=A(1,1)
      PIVOT2=A(2,1)
      PIVOT3=A(3,1)
      DO i=1,3
        A(1,i)=A(1,i)/PIVOT1
        A(2,i)=A(2,i)-PIVOT2*A(1,i)
        A(3,i)=A(3,i)-PIVOT3*A(1,i)
      ENDDO !i
      B(1)=B(1)/PIVOT1
      B(2)=B(2)-PIVOT2*B(1)
      B(3)=B(3)-PIVOT3*B(1)
      PIVMAX=MAX(ABS(A(2,2)),ABS(A(3,2)))
      CALL ASSERT(PIVMAX.NE.(0.d0),'Singular matrix',ERROR,*9999)
      IF(PIVMAX.EQ.ABS(A(3,2)))THEN !swap rows 2 and 3
        DO i=2,3
          TEMP1=A(2,i)
          A(2,i)=A(3,i)
          A(3,i)=TEMP1
        ENDDO !i
        TEMP1=B(2)
        B(2)=B(3)
        B(3)=TEMP1
      ENDIF
      PIVOT2=A(2,2)
      PIVOT3=A(3,2)
      DO i=2,3
        A(2,i)=A(2,i)/PIVOT2
        A(3,i)=A(3,i)-PIVOT3*A(2,i)
      ENDDO !i
      B(2)=B(2)/PIVOT2
      B(3)=B(3)-PIVOT3*B(2)
      TEMP3=B(3)/A(3,3)
      TEMP2=(B(2)-TEMP3*A(2,3))/A(2,2)
      TEMP1=(B(1)-A(1,3)*TEMP3-A(1,2)*TEMP2)/A(1,1)

      B(1)=TEMP1
      B(2)=TEMP2
      B(3)=TEMP3

      CALL EXITS('FSOLV3V')
      RETURN
 9999 CALL ERRORS('FSOLV3V',ERROR)
      CALL EXITS('FSOLV3V')
      RETURN 1
      END


