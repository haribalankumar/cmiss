      SUBROUTINE MINI_LINEAR_SYSTEM(MatrixSize,NonZeros,SparseCol,
     & SparseRow,Solution,SparseVal,RHS,ERROR,*)

C#### Subroutine MINI_LINEAR_SYSTEM
C###  Description:
C###  This function is uses UMFPACK to solve a linear system of 
C###  equations. It allows us to solve small linear problems within
C###  large problems without using unnecessarily large amounts of memory
C###  in the CMISS structure. 
C###  Created by ARC: 2010-01-10 by ARC

C###  Input:
C###  MatrixSize -> the size of the matrix.
C###  NonZeros -> The number of non-zero entries.
C###  SparseCol -> Define the column indices and Ax values for UMFPACK. 
C###  SparseRow ->Define the vector Ap for UMFPACK.
C###  SparseVal -> The values of non-zero entries.
C###  RHS -> The vector on the RHS of the linear system

C###  Output:
C###  Solution -> Solution vector
      IMPLICIT NONE
      include 'ptr00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'solver.inc'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'mach00.inc'


!     Parameter List
      INTEGER MatrixSize,NonZeros,SparseCol(NonZeros),
     &  SparseRow(MatrixSize+1)
      REAL*8 Solution(MatrixSize),SparseVal(Nonzeros),RHS(MatrixSize)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER OUTPUTCODE 
      REAL*8 UMF_PARAM(40)
      
       CALL ENTERS('MINI_LINEAR_SYSTEM',*9999)

       CALL UMFPACK4_FREE(ISOLV1_PTR(2),.TRUE.,RSOLV1_PTR(2),
     &        .TRUE.,ERROR,*9999)
        UMF_PARAM(1)=0.1d0
       CALL UMFPACK4_FACTOR(SparseVal,Matrixsize,NonZeros,SparseCol,
     &  SparseRow,ISOLV1_PTR(2),RSOLV1_PTR(2),UMF_PARAM,    
     &        SOLV_UMF_CONTROL(1,2),SOLV_UMF_INFO(1,2),
     &        SOLV_ANORM(2),
     &     ERROR,*9999) 
      CALL UMFPACK4_SOLVE(SparseVal,Matrixsize,NonZeros, SparseCol,
     &    SparseRow,RHS,SOLUTION,RSOLV1_PTR(2), 
     &    SOLV_UMF_CONTROL(1,2), SOLV_UMF_INFO(1,2), SOLV_ANORM(2),
     &   0,ERROR,*9999)  
       CALL UMFPACK4_FREE(ISOLV1_PTR(2),.TRUE.,RSOLV1_PTR(2),
     &        .TRUE.,ERROR,*9999) 


      CALL EXITS('MINI_LINEAR_SYSTEM')
      RETURN
 9999 CALL ERRORS('MINI_LINEAR_SYSTEM',ERROR)
      CALL EXITS('MINI_LINEAR_SYSTEM')
      RETURN 1
      END
