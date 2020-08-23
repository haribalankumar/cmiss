      SUBROUTINE INIT_SPARSE_MATRIX(COLLIST1,COLLIST2,ISC,ISR,
     '  LDA_COLLIST2,LDA_ROWLIST2,LIST_TYPE,M,NZMAX,NZTOT,
     '  ROWLIST1,ROWLIST2,SPARSENESS,MATRIX,ERROR,*)

C#### Subroutine: INIT_SPARSE_MATRIX
C###  Description:
C###    INIT_SPARSE_MATRIX initialises (sets to zero) a sparse matrix.

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA_COLLIST2,LDA_ROWLIST2,COLLIST1(0:*),
     '  COLLIST2(0:LDA_COLLIST2,*),ISC(*),ISR(*),LIST_TYPE,
     '  M,NZMAX,NZTOT,ROWLIST1(0:*),ROWLIST2(0:LDA_ROWLIST2,*),
     '  SPARSENESS
      REAL*8 MATRIX(NZMAX)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER col,i,j,noclist1,noclist2,norlist1,norlist2,nz,row

      CALL ENTERS('INIT_SPARSE_MATRIX',*9999)

      IF(SPARSENESS.EQ.0) THEN !No sparsity
        IF(LIST_TYPE.EQ.0) THEN
C$OMP PARALLEL DO
C$OMP&  PRIVATE(nz),
C$OMP&  SHARED(NZTOT, MATRIX)
          DO nz=1,NZTOT
            MATRIX(nz)=0.0d0
          ENDDO !nz
C$OMP END PARALLEL DO
        ELSE IF(LIST_TYPE.EQ.1) THEN
C$OMP PARALLEL DO
C$OMP&  PRIVATE(norlist1,i,noclist1,j,nz),
C$OMP&  SHARED(ROWLIST1, COLLIST1, MATRIX)
          DO norlist1=1,ROWLIST1(0)
            i=ROWLIST1(norlist1)
            DO noclist1=1,COLLIST1(0)
              j=COLLIST1(noclist1)
              nz=i+(j-1)*M
              MATRIX(nz)=0.0d0
            ENDDO !noclist (j)
          ENDDO !norlist (i)
C$OMP END PARALLEL DO
        ELSE IF(LIST_TYPE.EQ.2) THEN
C$OMP PARALLEL DO
C$OMP&  PRIVATE(norlist1,row,norlist2,i,noclist1,col,noclist2,j,nz),
C$OMP&  SHARED(ROWLIST1,ROWLIST2,COLLIST1,COLLIST2,MATRIX)
          DO norlist1=1,ROWLIST1(0)
            row=ROWLIST1(norlist1)
            DO norlist2=1,ROWLIST2(0,row)
              i=ROWLIST2(norlist2,row)
              DO noclist1=1,COLLIST1(0)
                col=COLLIST1(noclist1)
                DO noclist2=1,COLLIST2(0,col)
                  j=COLLIST2(noclist2,col)
                  nz=i+(j-1)*M
                  MATRIX(nz)=0.0d0
                ENDDO !noclist2
              ENDDO !noclist1
            ENDDO !norlist2
          ENDDO !norlist1
C$OMP END PARALLEL DO
        ELSE
          ERROR='>>Invalid list type'
          GOTO 9999
        ENDIF
      ELSE IF(SPARSENESS.EQ.1) THEN !compressed-row sparsity
        IF(LIST_TYPE.EQ.0) THEN
C$OMP PARALLEL DO
C$OMP&  PRIVATE(nz),
C$OMP&  SHARED(NZTOT, MATRIX)
          DO nz=1,NZTOT
            MATRIX(nz)=0.0d0
          ENDDO !nz
C$OMP END PARALLEL DO
        ELSE IF(LIST_TYPE.EQ.1) THEN
C$OMP PARALLEL DO
C$OMP&  PRIVATE(norlist1,i,nz),
C$OMP&  SHARED(ROWLIST1, ISR, MATRIX)
          DO norlist1=1,ROWLIST1(0)
            i=ROWLIST1(norlist1)
            DO nz=ISR(i),ISR(i+1)-1
              MATRIX(nz)=0.0d0
            ENDDO !nz
          ENDDO !norlist (i)
C$OMP END PARALLEL DO
        ELSE IF(LIST_TYPE.EQ.2) THEN
C$OMP PARALLEL DO
C$OMP&  PRIVATE(norlist1,row,norlist2,i,nz),
C$OMP&  SHARED(ROWLIST1, ROWLIST2, ISR, MATRIX)
          DO norlist1=1,ROWLIST1(0)
            row=ROWLIST1(norlist1)
            DO norlist2=1,ROWLIST2(0,row)
              i=ROWLIST2(norlist2,row)
              DO nz=ISR(i),ISR(i+1)-1
                MATRIX(nz)=0.0d0
              ENDDO !nz
            ENDDO !norlist2
          ENDDO !norlist1
C$OMP END PARALLEL DO
        ELSE
          ERROR='>>Invalid list type'
          GOTO 9999
        ENDIF
      ELSE IF(SPARSENESS.EQ.2.OR.SPARSENESS.EQ.4.OR.SPARSENESS.EQ.5)
     '  THEN !row-column sparsity
C$OMP PARALLEL DO
C$OMP&  PRIVATE(nz),
C$OMP&  SHARED(NZTOT, MATRIX)
        DO nz=1,NZTOT
          MATRIX(nz)=0.0d0
        ENDDO !nz
C$OMP END PARALLEL DO
      ELSE IF(SPARSENESS.EQ.3) THEN !compressed-column sparsity
        IF(LIST_TYPE.EQ.0) THEN
C$OMP PARALLEL DO
C$OMP&  PRIVATE(nz),
C$OMP&  SHARED(NZTOT, MATRIX)
          DO nz=1,NZTOT
            MATRIX(nz)=0.0d0
          ENDDO !nz
C$OMP END PARALLEL DO
        ELSE IF(LIST_TYPE.EQ.1) THEN
C$OMP PARALLEL DO
C$OMP&  PRIVATE(noclist1,j,nz),
C$OMP&  SHARED(COLLIST1, ISC, MATRIX)
          DO noclist1=1,COLLIST1(0)
            j=COLLIST1(noclist1)
            DO nz=ISC(j),ISC(j+1)-1
              MATRIX(nz)=0.0d0
            ENDDO !nz
          ENDDO !norlist (j)
C$OMP END PARALLEL DO
        ELSE IF(LIST_TYPE.EQ.2) THEN
C$OMP PARALLEL DO
C$OMP&  PRIVATE(noclist1,col,noclist2,j,nz),
C$OMP&  SHARED(COLLIST1, COLLIST2, ISC, MATRIX)
          DO noclist1=1,COLLIST1(0)
            col=COLLIST1(noclist1)
            DO noclist2=1,COLLIST2(0,col)
              j=COLLIST2(noclist2,col)
              DO nz=ISC(j),ISC(j+1)-1
                MATRIX(nz)=0.0d0
              ENDDO !nz
            ENDDO !noclist2
          ENDDO !noclist1
C$OMP END PARALLEL DO
        ENDIF
      ELSE
        ERROR='>>Unknown sparsity'
        GOTO 9999
      ENDIF

      CALL EXITS('INIT_SPARSE_MATRIX')
      RETURN
 9999 CALL ERRORS('INIT_SPARSE_MATRIX',ERROR)
      CALL EXITS('INIT_SPARSE_MATRIX')
      RETURN 1
      END


