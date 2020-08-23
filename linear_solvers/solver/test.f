
      SUBROUTINE RAND_WEIGHT(A,LDA,N,B,SPARSE_A,ISC_A,ISR_A,ERROR,*)

C#### Subroutine: RAND_WEIGHT
C###  Description:
C###    Randomly weight each equation.
C###  Written by Stuart Norris 22/05/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,SPARSE_A,ISC_A(*),ISR_A(*)
      REAL*8 A(*),B(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,IJ,ISEED
      REAL*8 WEIGHT(N)
!     Subroutines
      INTEGER IRAND
      REAL*8  RAND
      EXTERNAL IRAND,RAND,SRAND,DSCAL


      ISEED=2
      CALL SRAND(ISEED)
      DO I=1,N
        WEIGHT(I)=SIGN(1000.0D0*RAND()+1.0D0,RAND()-0.5D0)
      ENDDO

C     Weight dense arrays
      IF(SPARSE_A.EQ.0) THEN
        DO I=1,N
          CALL DSCAL(N,WEIGHT(I),A(I),LDA)
          B(I)=WEIGHT(I)*B(I)          
        ENDDO

C     Weight sparse arrays
      ELSE IF(SPARSE_A.EQ.1) THEN
        DO I=1,N
          DO IJ=ISR_A(I),ISR_A(I+1)-1
            A(IJ)=WEIGHT(I)*A(IJ)
          ENDDO
          B(I)=B(I)*WEIGHT(I)
        ENDDO

      ELSE
        ERROR='>>Sparsity type not implemented'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('RAND_WEIGHT',ERROR)
      RETURN 1
      END


      SUBROUTINE MIX_ARRAYS(A,LDA,N,NZA,B,SPARSE_A,ISC_A,ISR_A,ERROR,*)

C#### Subroutine: MIX_ARRAYS
C###  Description:
C###    Mix up the array ordering, so as to test pivoting in the
C###    solvers.
C###  Written by Stuart Norris 17/05/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,NZA,SPARSE_A,ISC_A(*),ISR_A(*)
      REAL*8 A(*),B(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,IROW,ITMP,ISEED,IJ_TMP,IJ_A,NJ
      INTEGER ISR_TMP(N+1),ISC_TMP(NZA),IDX(N)
      REAL*8 A_TMP(NZA),B_TMP(N)
!     Subroutines
      INTEGER IRAND
      EXTERNAL IRAND,SRAND,DCOPY


      ISEED=2
      CALL SRAND(ISEED)
C     Scramble up the index to get our final ordering
      DO I=1,N
        IDX(I)=I
      ENDDO
      DO I=1,N
        IROW=MOD(IRAND(),N)+1
        IF(IROW.NE.I) THEN
          ITMP=IDX(I)
          IDX(I)=IDX(IROW)
          IDX(IROW)=ITMP
        ENDIF
      ENDDO

C     Load the data
      IF(SPARSE_A.EQ.0) THEN
        DO I=1,N
          CALL DCOPY(N,A(IDX(I)),LDA,A_TMP(I),LDA)
          B_TMP(I)=B(IDX(I))
        ENDDO
        CALL DCOPY(NZA,A_TMP(1),1,A(1),1)
        CALL DCOPY(N,B_TMP(1),1,B(1),1)

      ELSE IF(SPARSE_A.EQ.1) THEN

C       Create the mixed array
        ISR_TMP(1)=1
        DO I=1,N
          NJ=ISR_A(IDX(I)+1)-ISR_A(IDX(I))

          ISR_TMP(I+1)=ISR_TMP(I)+NJ
          DO J=1,NJ
            IJ_TMP=ISR_TMP(I)+J-1
            IJ_A=ISR_A(IDX(I))+J-1

            ISC_TMP(IJ_TMP)=ISC_A(IJ_A)
            A_TMP(IJ_TMP)=A(IJ_A)
          ENDDO
          B_TMP(I)=B(IDX(I))
        ENDDO

C       Copy back to the original array
        DO I=1,N+1
          ISR_A(I)=ISR_TMP(I)
        ENDDO
        DO I=1,NZA
          ISC_A(I)=ISC_TMP(I)
        ENDDO
        CALL DCOPY(NZA,A_TMP(1),1,A(1),1)
        CALL DCOPY(N,B_TMP(1),1,B(1),1)

      ELSE
        ERROR='>>Sparsity type not implemented'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('MIX_ARRAYS',ERROR)
      RETURN 1
      END


      SUBROUTINE WRITE_MATRIX_DIMS(SPARSE_A,N,NZA,ANORM,BNORM,ERROR,*)

C#### Subroutine: WRITE_MATRIX_DIMS
C###  Description:
C###    Write out the matrix dimensions.
C###  Written by Stuart Norris 07/08/02

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER SPARSE_A,N,NZA
      REAL*8 ANORM,BNORM
      CHARACTER ERROR*(*)
!     Local Variables
!     External routines


      IF(SPARSE_A.EQ.0) THEN
        WRITE(OP_STRING,'('' Dense Linear System'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ELSE
        WRITE(OP_STRING,'('' Sparse Linear System'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      WRITE(OP_STRING,7000) ' No. Equations      : ',N
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,7000) ' No. Elements       : ',NZA
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,7200) ' Norm of A          : ',ANORM
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,7200) ' Norm of RHS        : ',BNORM
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'()')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      RETURN

 7000 FORMAT(A,I18)
 7100 FORMAT(A,7X,F10.5,A)
 7200 FORMAT(A,7X,1P,E11.4)

 9999 CALL ERRORS('WRITE_MATRIX_DIMS',ERROR)
      RETURN 1
      END


      SUBROUTINE CHECK_POS_DEF(A,LDA,N,SPARSE_A,ISC_A,ISR_A,OUTPUTCODE,
     '  ERROR,*)

C#### Subroutine: 
C###  Description:
C###    Determines if a matrix is positive definite, a condition that
C###    required by some direct and iterative linear solvers.
C###
C###    This is a bit difficult, so at the moment we just check for 
C###     diagonal dominance.
C###  Written by Stuart Norris 24/05/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,SPARSE_A,ISC_A(*),ISR_A(*),OUTPUTCODE
      REAL*8 A(*)
      CHARACTER ERROR*(*)
!     Local Variables


      CALL CHECK_DIAG_DOM(A,LDA,N,SPARSE_A,ISC_A,ISR_A,OUTPUTCODE,
     '  ERROR,*9999)

      RETURN

 9999 CALL ERRORS('CHECK_POS_DEF',ERROR)
      RETURN 1
      END


      SUBROUTINE CHECK_DIAG_DOM(A,LDA,N,SPARSE_A,ISC_A,ISR_A,OUTPUTCODE,
     '  ERROR,*)

C#### Subroutine: CHECK_DIAG_DOM
C###  Description:
C###    Determines if a matrix is diagonally dominant, a condition that
C###    required by some direct and iterative linear solvers.
C###  Written by Stuart Norris 26/02/02

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER LDA,N,SPARSE_A,ISC_A(*),ISR_A(*),OUTPUTCODE
      REAL*8 A(*)
      CHARACTER ERROR*(*)
!     Local Parameters
      REAL*8 EPS
      PARAMETER(EPS=1.0D-8)
!     Local Variables
      INTEGER I,J,NELEM,NDIAG,NZERO
      REAL*8 SUM,DIAG,RATIO_MIN,RATIO_MAX,RATIO_MEAN,RATIO
      LOGICAL HAVE_DIAG
!     Functions
      REAL*8 DLAMCH


      WRITE(OP_STRING,'('' Check Matrix Diagonal Dominance: '')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      CALL WRITES(IOOP,' ',ERROR,*9999)
      NDIAG=0
      NELEM=0
      NZERO=0
      RATIO_MEAN=0.0D0
      RATIO_MIN=0.0D0
      RATIO_MAX=0.0D0
      HAVE_DIAG=.TRUE.

      DO I=1,N

C       Non-sparse matricies
        IF(SPARSE_A.EQ.0) THEN
          SUM=0.0D0
          DO J=1,I-1
            SUM=SUM+ABS(A(I+LDA*(J-1)))
          ENDDO
          DIAG=ABS(A(I+LDA*(I-1)))
          DO J=I+1,N
            SUM=SUM+ABS(A(I+LDA*(J-1)))
          ENDDO

C       Sparse matricies: CMISS compressed row
        ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN
          SUM=0.0D0
          DIAG=0.0D0
          HAVE_DIAG=.FALSE.
          DO J=ISR_A(I),ISR_A(I+1)-1
            IF(ISC_A(J).EQ.I) THEN
              DIAG=ABS(A(J))
              HAVE_DIAG=.TRUE.
            ELSE
              SUM=SUM+ABS(A(J))
            ENDIF
          ENDDO

C       Other non supported format
        ELSE
          ERROR='>>Sparsity type not implemented'
          GOTO 9999
        ENDIF

C       Increase diagonal by tolerence
        IF(DIAG.GT.EPS) THEN
          DIAG=DIAG*(1.0D0+EPS)
        ELSE IF(DIAG.NE.0.0D0) THEN
          DIAG=DIAG+EPS
        ENDIF

C       Calculate ratio, and update
        IF(SUM.NE.0.0D0) THEN
          RATIO=DIAG/SUM
        ELSE
          RATIO=DLAMCH('O')
        ENDIF
        IF(I.EQ.1) THEN
          RATIO_MIN=RATIO
          RATIO_MAX=RATIO
        ENDIF
        RATIO_MIN=MIN(RATIO_MIN,RATIO)
        RATIO_MAX=MAX(RATIO_MAX,RATIO)
        RATIO_MEAN=RATIO_MEAN*(DBLE(I-1)/DBLE(I)) + RATIO/DBLE(I)

C       Process output
        IF(OUTPUTCODE.GE.3) THEN
          IF(.NOT.HAVE_DIAG) THEN
            WRITE(OP_STRING,'('' No Diagonal Element:     '',I9))') I
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(DIAG.EQ.0.0D0) THEN
            WRITE(OP_STRING,'('' Zero on Diagonal:        '',I9))') I
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(SUM.GT.DIAG) THEN
            WRITE(OP_STRING,'('' Not Diagonally Dominant: '','
     '        //'I9,2(1X,D11.4))') I,SUM,DIAG
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

        IF(.NOT.HAVE_DIAG) THEN
          NELEM=NELEM+1
        ELSE IF(DIAG.EQ.0.0D0) THEN
          NZERO=NZERO+1
        ELSE IF(SUM.GT.DIAG) THEN
          NDIAG=NDIAG+1
        ENDIF

      ENDDO

C     Summary
      IF(NELEM.GT.0 .OR. NDIAG.GT.0) THEN
        WRITE(OP_STRING,'('' Matrix is Not Diagonally Dominant'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(5X,I7,'' non dominant rows'')') NDIAG
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        IF(SPARSE_A.NE.0) THEN
          WRITE(OP_STRING,'(5X,I7,'' missing diagonals'')') NELEM
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        WRITE(OP_STRING,'(5X,I7,'' zeros on diagonals'')') NZERO
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(5X,I7,'' total rows'')') N
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ELSE
        WRITE(OP_STRING,'('' Matrix is Diagonally Dominant'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      WRITE(OP_STRING,'('' Maximum dominance'',2X,1P,D11.4)') RATIO_MAX
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Minimum dominance'',2X,1P,D11.4)') RATIO_MIN
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Mean dominance'',5X,1P,D11.4)') RATIO_MEAN
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      CALL WRITES(IOOP,' ',ERROR,*9999)

      RETURN

 9999 CALL ERRORS('CHECK_DIAG_DOM',ERROR)
      RETURN 1
      END


      SUBROUTINE CHECK_SYMM(A,LDA,N,NZA,SPARSE_A,ISC_A,ISR_A,OUTPUTCODE,
     '  ERROR,*)

C#### Subroutine: CHECK_SYMM
C###  Description:
C###    Checks if a matrix is symmetric.
C###  Written by Stuart Norris 28/02/02

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER LDA,N,NZA,SPARSE_A,ISC_A(*),ISR_A(*),OUTPUTCODE
      REAL*8 A(*)
      CHARACTER ERROR*(*)
!     Local Parameters
      REAL*8 EPS
      PARAMETER(EPS=1.0D-8)
!     Local Variables
      INTEGER I,J,IJ,NELEM,NROWS,NDIFF,NFILL,NZERO
      INTEGER IWORK(N),ITC_A(N+1),ITR_A(NZA),IDX_A(NZA)
      REAL*8 DIFF,AII,AIJ,AJI
      LOGICAL LROW


      WRITE(OP_STRING,'('' Check Matrix Symmetry: '')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      NELEM=0
      NDIFF=0
      NROWS=0
      NFILL=0
      NZERO=0
      AIJ=0.0D0

C     Non-sparse matricies
      IF(SPARSE_A.EQ.0) THEN

        DO I=1,N-1
          LROW=.FALSE.

          AII=A(I+LDA*(I-1))
          IF(AII.NE.0.0D0) NELEM=NELEM+1

          DO J=I+1,N
            AIJ=A(I+LDA*(J-1))
            AJI=A(J+LDA*(I-1))

            IF(AIJ.NE.0.0D0) NELEM=NELEM+1
            IF(AJI.NE.0.0D0) NELEM=NELEM+1

            IF(ABS(AIJ).GT.EPS) THEN
              DIFF=ABS(AIJ-AJI)/ABS(AIJ)
            ELSE
              DIFF=ABS(AIJ-AJI)
            ENDIF

            IF(DIFF.GT.EPS) THEN
              NDIFF=NDIFF+1
              LROW=.TRUE.
            ENDIF
          ENDDO

          IF(LROW) NROWS=NROWS+1
        ENDDO

C     Sparse matricies: CMISS compressed row
      ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

C       Generate a compressed column index for use in the factorisation.

        CALL TRANS_INDEX(N,NZA,ISC_A,ISR_A,ITC_A,ITR_A,IDX_A,IWORK,
     '    ERROR,*9999)

        DO I=1,N-1
          LROW=.FALSE.

C         Set the column index
          IF(I.EQ.1) THEN
            DO J=1,N
              IWORK(J)=0
            ENDDO
          ELSE
            DO J=ITC_A(I-1),ITC_A(I)-1
              IWORK(ITR_A(J))=0
            ENDDO
          ENDIF
          DO J=ITC_A(I),ITC_A(I+1)-1
            IWORK(ITR_A(J))=IDX_A(J)
          ENDDO

          DO IJ=ISR_A(I),ISR_A(I+1)-1
            J=ISC_A(IJ)

C           Above the diagonal
            IF(J.GT.I) THEN

C             Symmetric fill
              IF(IWORK(J).NE.0) THEN
                AIJ=A(IJ)
                AJI=A(IWORK(J))

                IF(AIJ.EQ.0.0D0) NZERO=NZERO+1
                IF(AJI.EQ.0.0D0) NZERO=NZERO+1

                IF(ABS(AIJ).GT.EPS) THEN
                  DIFF=ABS(AIJ-AJI)/ABS(AIJ)
                ELSE
                  DIFF=ABS(AIJ-AJI)
                ENDIF

                IF(DIFF.GT.EPS) THEN
                  NDIFF=NDIFF+2
                  LROW=.TRUE.
                ENDIF

C             Asymmetric fill
              ELSE
                IF(AIJ.EQ.0.0D0) THEN
                  NZERO=NZERO+1
                ELSE
                  NDIFF=NDIFF+1
                  LROW=.TRUE.
                ENDIF
                NFILL=NFILL+1
              ENDIF

            ELSE IF(J.EQ.I) THEN
              AII=A(IJ)
              IF(AII.EQ.0.0D0) NZERO=NZERO+1
            ENDIF
          ENDDO

          IF(LROW) NROWS=NROWS+1
        ENDDO

        NELEM=NZA

      ELSE
         ERROR='>>Sparsity type not implemented'
         GOTO 9999
      ENDIF

C     Summary
      IF(NDIFF.GT.0) THEN
        WRITE(OP_STRING,'('' Matrix is Not Symmetric'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(I12,'' asymmetric elements'')') NDIFF
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        IF(SPARSE_A.NE.0) THEN
          WRITE(OP_STRING,'(I12,'' asymmetric fill elements'')') NFILL
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(I12,'' zeroed elements'')') NZERO
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        WRITE(OP_STRING,'(I12,'' total elements'')') NELEM
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(I12,'' asymmetric rows'')') NROWS
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(I12,'' total rows'')') N
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ELSE
        WRITE(OP_STRING,'('' Matrix is Symmetric'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      CALL WRITES(IOOP,' ',ERROR,*9999)

      RETURN

 9999 CALL ERRORS('CHECK_SYMM',ERROR)
      RETURN 1
      END


      SUBROUTINE CHECK_FILL(A,LDA,N,SPARSE_A,ISC_A,ISR_A,OUTPUTCODE,
     '  ERROR,*)

C#### Subroutine: CHECK_FILL
C###  Description:
C###    Check the fill pattern of a matrix
C###  Written by Stuart Norris 14/05/02

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER LDA,N,SPARSE_A,ISC_A(*),ISR_A(*),OUTPUTCODE
      REAL*8 A(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,IJ,NROWS,NCOLS,NZERO,NELEM,NPHYS,NZROW,NERR
      INTEGER*8 NTOTAL
      REAL*8 AIJ,SUM,PELEM,PZERO,PFILL,COL(N)
!     External routines
      REAL*8 DASUM
      EXTERNAL DASUM


      WRITE(OP_STRING,'('' Check Matrix Fill Pattern: '')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      NROWS=0
      NCOLS=0
      NZERO=0
      NELEM=0
      NPHYS=0
      NZROW=0
      NERR=0

C     Non-sparse matricies
      IF(SPARSE_A.EQ.0) THEN

        DO I=1,N
          SUM=DASUM(N,A(I),LDA)
          IF(SUM.EQ.0.0D0) THEN
            NROWS=NROWS+1
          ENDIF

          SUM=DASUM(N,A((I-1)*LDA+1),1)
          IF(SUM.EQ.0.0D0) THEN
            NCOLS=NCOLS+1
          ENDIF

          DO J=1,N
            AIJ=A(J+(I-1)*LDA)
            IF(AIJ.EQ.0.0D0) NZERO=NZERO+1
          ENDDO
        ENDDO
        NPHYS=N**2
        NELEM=NPHYS-NZERO

C     Sparse matricies: CMISS compressed row
      ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

        DO J=1,N
          COL(J)=0.0D0
        ENDDO

        DO I=1,N
          SUM=0.0D0
          DO IJ=ISR_A(I),ISR_A(I+1)-1
            AIJ=ABS(A(IJ))

            IF(AIJ.EQ.0.0D0) THEN
              NZERO=NZERO+1
            ELSE
              NELEM=NELEM+1
            ENDIF

            SUM=SUM+ABS(A(IJ))
            COL(ISC_A(IJ))=COL(ISC_A(IJ))+ABS(A(IJ))
          ENDDO

          IF(ISR_A(I+1)-ISR_A(I).LT.0) THEN
            NERR=NERR+1
          ELSE IF(ISR_A(I+1)-ISR_A(I).EQ.0) THEN
            NZROW=NZROW+1
          ELSE IF(SUM.EQ.0.0D0) THEN
            NROWS=NROWS+1
          ENDIF
        ENDDO

        DO J=1,N
          IF(COL(J).EQ.0.0D0) NCOLS=NCOLS+1
        ENDDO
        NPHYS=ISR_A(N+1)-1

      ELSE
         ERROR='>>Sparsity type not implemented'
         GOTO 9999
      ENDIF

      NTOTAL=N
      NTOTAL=NTOTAL**2
      PELEM=(DBLE(NELEM)/DBLE(NPHYS))*100.0D0
      PZERO=(DBLE(NZERO)/DBLE(NPHYS))*100.0D0
      PFILL=(DBLE(NPHYS)/DBLE(NTOTAL))*100.0D0

C     Summary
      WRITE(OP_STRING,'(I18,'' zero sum rows'')') NROWS
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(I18,'' zero sum columns'')') NCOLS
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(I18,'' total rows'')') N
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      IF(SPARSE_A.NE.0) THEN
        WRITE(OP_STRING,'(I18,'' erroniously indexed rows'')') NERR
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(I18,'' zero element rows'')') NZROW
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      CALL WRITES(IOOP,' ',ERROR,*9999)

      WRITE(OP_STRING,'(I18,'' nonzero elements'')') NELEM
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(I18,'' zero elements'')') NZERO
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(I18,'' physical memory'')') NPHYS
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(I18,'' matrix size'')') NTOTAL
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      CALL WRITES(IOOP,' ',ERROR,*9999)

      WRITE(OP_STRING,'(8X,F9.5,''% nonzero elements'')') PELEM
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(8X,F9.5,''% zero elements'')') PZERO
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(8X,F9.5,''% fill'')') PFILL
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      CALL WRITES(IOOP,' ',ERROR,*9999)

      RETURN

 9999 CALL ERRORS('CHECK_FILL',ERROR)
      RETURN 1
      END


      SUBROUTINE CHECK_SYMM_FILL(A,LDA,N,NZA,SPARSE_A,ISC_A,ISR_A,
     '  OUTPUTCODE,ERROR,*)

C#### Subroutine: CHECK_SYMM_FILL
C###  Description:
C###    Checks if a matrix has matrix structure
C###  Written by Stuart Norris 24/05/02

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER LDA,N,NZA,SPARSE_A,ISC_A(*),ISR_A(*),OUTPUTCODE
      REAL*8 A(*)
      CHARACTER ERROR*(*)
!     Local Parameters
      REAL*8 EPS
      PARAMETER(EPS=1.0D-8)
!     Local Variables
      INTEGER I,J,IJ
      INTEGER NELEM,NZERO,NNSYM,NVSYM,NROWS,NDZERO,NDIAG,NPHYS
      INTEGER IWORK(N),ITC_A(N+1),ITR_A(NZA),IDX_A(NZA)
      REAL*8 AIJ,AJI
      LOGICAL LROWS,LDIAG


      WRITE(OP_STRING,'('' Check Matrix Structural Symmetry: '')')
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      NELEM=0
      NZERO=0
      NNSYM=0
      NVSYM=0
      NROWS=0
      NDZERO=0
      NDIAG=0

C     Non-sparse matricies
      IF(SPARSE_A.EQ.0) THEN

        DO I=1,N
          LROWS=.FALSE.

          DO J=1,N
            AIJ=A(I+(J-1)*LDA)
            AJI=A(J+(I-1)*LDA)

            IF(AIJ.EQ.0.0D0) THEN
              NZERO=NZERO+1
            ELSE IF(AJI.EQ.0.0D0) THEN
              NNSYM=NNSYM+1
              LROWS=.TRUE.
            ENDIF

            IF(I.EQ.J) THEN
              IF(AIJ.EQ.0.0) NDZERO=NDZERO+1
            ENDIF
          ENDDO

          IF(LROWS) NROWS=NROWS+1
        ENDDO
        NPHYS=N**2
        NELEM=NPHYS-NZERO

C     Sparse matricies: CMISS compressed row
      ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN

C       Generate a compressed column index for use in the factorisation.

        CALL TRANS_INDEX(N,NZA,ISC_A,ISR_A,ITC_A,ITR_A,IDX_A,IWORK,
     '    ERROR,*9999)

        DO I=1,N-1
          LROWS=.FALSE.
          LDIAG=.FALSE.

C         Set the column index
          IF(I.EQ.1) THEN
            DO J=1,N
              IWORK(J)=0
            ENDDO
          ELSE
            DO J=ITC_A(I-1),ITC_A(I)-1
              IWORK(ITR_A(J))=0
            ENDDO
          ENDIF
          DO J=ITC_A(I),ITC_A(I+1)-1
            IWORK(ITR_A(J))=IDX_A(J)
          ENDDO

          DO IJ=ISR_A(I),ISR_A(I+1)-1
            J=ISC_A(IJ)
            IF(A(IJ).EQ.0.0D0) NZERO=NZERO+1

C           Off diagonal
            IF(J.NE.I) THEN

              IF(IWORK(J).EQ.0) THEN
                NNSYM=NNSYM+1
                IF(A(IJ).NE.0.0D0) NVSYM=NVSYM+1
                LROWS=.TRUE.
              ENDIF

C           On Diagonal
            ELSE
              IF(A(IJ).EQ.0.0D0) NDZERO=NDZERO+1
              LDIAG=.TRUE.
            ENDIF
          ENDDO

          IF(.NOT.LDIAG) NDIAG=NDIAG+1
          IF(LROWS) NROWS=NROWS+1
        ENDDO
        NELEM=NZA-NZERO
        NPHYS=NZA

      ELSE
         ERROR='>>Sparsity type not implemented'
         GOTO 9999
      ENDIF

C     Summary
      IF(NVSYM.GT.0 .OR. NNSYM.GT.0) THEN
        WRITE(OP_STRING,'('' Matrix isnt Structurally Symmetric'')')
      ELSE
        WRITE(OP_STRING,'('' Matrix is Structurally Symmetric'')')
      ENDIF
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

      WRITE(OP_STRING,'(I12,'' asymmetric nonzeros'')') NVSYM
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      IF(SPARSE_A.NE.0) THEN
        WRITE(OP_STRING,'(I12,'' asymmetric elements'')') NNSYM
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      CALL WRITES(IOOP,' ',ERROR,*9999)

      WRITE(OP_STRING,'(I12,'' zero elements'')') NZERO
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(I12,'' zero diagonals'')') NDZERO
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      IF(SPARSE_A.NE.0) THEN
        WRITE(OP_STRING,'(I12,'' missing diagonal elements'')') NDIAG
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      WRITE(OP_STRING,'(I12,'' total nonzeros'')') NELEM
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(I12,'' physical memory'')') NPHYS
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(I12,'' matrix size'')') N**2
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      CALL WRITES(IOOP,' ',ERROR,*9999)      

      WRITE(OP_STRING,'(I12,'' asymmetric rows'')') NROWS
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(I12,'' total rows'')') N
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      CALL WRITES(IOOP,' ',ERROR,*9999)

      RETURN

 9999 CALL ERRORS('CHECK_SYMM_FILL',ERROR)
      RETURN 1
      END
