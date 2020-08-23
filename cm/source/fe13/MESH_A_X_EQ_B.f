      SUBROUTINE MESH_A_X_EQ_B(MATRIX,VECTOR,SOLUTION,ERROR,*)

C#### Subroutine: MESH_A_X_EQ_B
C###  Description:

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
      REAL*8 MATRIX(3,3),VECTOR(3),SOLUTION(3)
      CHARACTER ERROR*(*)
      !Local variables
      INTEGER i,j,k,pivot_row
      REAL*8 A(3,4),max,pivot_value,TEMP(4)

      CALL ENTERS('MESH_A_X_EQ_B',*9999)

      DO i=1,3
        DO j=1,3
          A(i,j)=MATRIX(i,j)
        ENDDO !j
        A(i,4)=VECTOR(i)
      ENDDO !i
      
      DO k=1,2
        max=0.d0
        DO i=k,3
          IF(DABS(A(i,k)).GT.max)THEN
            max=DABS(A(i,k))
            pivot_row=i
          ENDIF
        ENDDO !i
        IF(pivot_row.NE.k)THEN
          DO j=1,4
            TEMP(j)=A(k,j)
            A(k,j)=A(pivot_row,j)
            A(pivot_row,j)=TEMP(j)
          ENDDO !j
        ENDIF
        pivot_value=A(k,k)
        DO j=1,4
          A(k,j)=A(k,j)/pivot_value
        ENDDO !j
        DO i=k+1,3
          DO j=k+1,4
            A(i,j)=A(i,j)-A(i,k)*A(k,j)
          ENDDO
          A(i,k)=0.d0
        ENDDO
      ENDDO !N
      A(3,4)=A(3,4)/A(3,3)
      A(2,4)=A(2,4)-A(3,4)*A(2,3)
      A(1,4)=A(1,4)-A(3,4)*A(1,3)-A(2,4)*A(1,2)
      DO i=1,3
        SOLUTION(i)=A(i,4)
      ENDDO

      CALL EXITS('MESH_A_X_EQ_B')
      RETURN
 9999 CALL ERRORS('MESH_A_X_EQ_B',ERROR)
      CALL EXITS('MESH_A_X_EQ_B')
      RETURN 1
      END



