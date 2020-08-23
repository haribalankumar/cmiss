      SUBROUTINE LADDERSOL_MATRIX(NonZeros,MatrixSize,submatrixsize,
     & SparseCol,SparseRow,SparseVal,RHS,Pin,Pout,ERROR,*)
C#### Subroutine LADDERSOL_MATRIX
C###  CREATED: JAN 2010 by ARC

C###  Description:
C###  Sets up ladder matrix entries that are independent of iteration
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'pulm00.cmn'
      !INPUT AND OUTPUT PARAMETER LIST
      INTEGER NonZeros,MatrixSize,submatrixsize,SparseRow(MatrixSize+1),
     &  SparseCol(Nonzeros)
      REAL*8 SparseVal(Nonzeros),RHS(MatrixSize),Pin,Pout
      CHARACTER ERROR*(*)
      !Local variables
      INTEGER count1,count,i,j

      CALL ENTERS('LADDERSOL_MATRIX',*9999)
C###  Define the vector Ap (for UMFPACK) this vector starts at zero and then each entry is the cumulative total of non-zero entries as you step through columns from L to right
      SparseRow(1)=1
      SparseRow(2)=3
      SparseRow(3)=7
      SparseRow(4)=9
      SparseRow(5)=13
      DO i=2,num_symm_gen-1
         SparseRow(6+4*(i-2))=SparseRow(6+4*(i-2)-1)+2+i
         SparseRow(7+4*(i-2))=SparseRow(6+4*(i-2))+3+i
         SparseRow(8+4*(i-2))=SparseRow(7+4*(i-2))+2+i
         SparseRow(9+4*(i-2))=SparseRow(8+4*(i-2))+3+i
      ENDDO
      DO i=1,num_symm_gen
        SparseRow(4*num_symm_gen-3+i)=SparseRow(submatrixsize+i)+3
      ENDDO
      SparseRow(MatrixSize+1)=SparseRow(MatrixSize)+num_symm_gen+1
      
C###  DEFINE THE RHS VECTOR (size=Matrix size)- These are the BCS      
      DO count1=1,MatrixSize
         RHS(count1)=0.d0
      ENDDO
      RHS(1)=-Pin
      RHS(3)=Pout

C###  Define the column indices and Ax values for UMFPACK. These step through the rowss and give the column index and the value. (start at column 1)

C...  ENTRIES THAT ARE INDEPENDENT OF ITERATION      
      DO i=1,Nonzeros
         SparseCol(i)=0
         SparseVal(i)=0.d0
      ENDDO
C...  CONSERVATION OF FLOW
      SparseCol(NonZeros)=MatrixSize
      SparseVal(NonZeros)=1.d0
      DO i=1,num_symm_gen
         SparseCol(NonZeros-i)=MatrixSize-i
         SparseVal(NonZeros-i)=-2.d0**(num_symm_gen+1-i)
         SparseCol(NonZeros-num_symm_gen-1-3*(i-1))=MatrixSize-i
      ENDDO
C...  CAPILLARY CONNECTIONS
C...  Prior to final generation..      
      DO i=1,num_symm_gen-1
        SparseCol(NonZeros-4*num_symm_gen+3*(i-1))=1+4*(i-1)
        SparseVal(NonZeros-4*num_symm_gen+3*(i-1))=1.d0
        SparseCol(NonZeros-4*num_symm_gen+3*(i-1)+1)=3+4*(i-1)
        SparseVal(NonZeros-4*num_symm_gen+3*(i-1)+1)=-1.d0
      ENDDO
C...  FINAL GENERATION
      SparseCol(NonZeros-num_symm_gen-2)=MatrixSize-num_symm_gen-1
      SparseVal(NonZeros-num_symm_gen-2)=-1.d0
      SparseCol(NonZeros-num_symm_gen-3)=MatrixSize-num_symm_gen-3
      SparseVal(NonZeros-num_symm_gen-3)=1.d0
C...  First generation
      SparseCol(1)=1
      Sparseval(1)=-1.d0
      SparseCol(2)=MatrixSize
      SparseCol(3)=1
      SparseVal(3)=1.d0
      SparseCol(4)=2
      SparseVal(4)=-1.d0
      SparseCol(5)=submatrixsize+1
      SparseCol(6)=MatrixSize
      SparseCol(7)=3
      SparseVal(7)=1.d0
      SparseCol(8)=MatrixSize
      SparseCol(9)=3
      SparseVal(9)=-1.d0
      SparseCol(10)=4
      SparseVal(10)=1.d0
      SparseCol(11)=submatrixsize+1
      SparseCol(12)=MatrixSize
      count=12
      DO i=2,num_symm_gen-1
      SparseCol(count+1)=4*i-6
      SparseVal(count+1)=1.d0
      SparseCol(count+2)=4*i-3
      SparseVal(count+2)=-1.d0
        DO j=2,i
          SparseCol(count+3+(j-2))=submatrixsize+j-1
        ENDDO
      count=count+i+2
      SparseCol(count)=MatrixSize
      SparseCol(count+1)=5+4*(i-2)
      SparseVal(count+1)=1.d0
      SparseCol(count+2)=6+4*(i-2)
      SparseVal(count+2)=-1.d0
        DO j=2,i+1
          SparseCol(count+3+(j-2))=submatrixsize+j-1
        ENDDO
      count=count+i+3
      SparseCol(count)=MatrixSize
      SparseCol(count+1)=4*i-4
      SparseVal(count+1)=-1.d0
      SparseCol(count+2)=4*i-1
      SparseVal(count+2)=1.d0
         DO j=2,i
          SparseCol(count+3+(j-2))=submatrixsize+j-1
        ENDDO      
      count=count+i+2
      SparseCol(count)=MatrixSize
      SparseCol(count+1)=4*i-1
      SparseVal(count+1)=-1.d0
      SparseCol(count+2)=4*i
      SparseVal(count+2)=1.d0
        DO j=2,i+1
          SparseCol(count+3+(j-2))=submatrixsize+j-1
        ENDDO
      count=count+i+3
      SparseCol(count)=MatrixSize
      ENDDO

C###  EXIT SUBROUTINE 'LADDERSOL_MATRIX'
      CALL EXITS('LADDERSOL_MATRIX')
      RETURN
 9999 CALL ERRORS('LADDERSOL_MATRIX',ERROR)
      CALL EXITS('LADDERSOL_MATRIX')
      RETURN 1
      END

