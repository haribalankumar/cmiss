
      SUBROUTINE SOLVER_DUMMY(A,LDA,N,SPARSE_A,PIVOT_A,IPIVOT_A,
     '  ISC_A,ISR_A,ERROR,*)

C#### Subroutine: SOLVER_DUMMY
C###  Description:
C###    Dummy routine
C###  Written by Stuart Norris 25/03/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,SPARSE_A,IPIVOT_A(*),ISC_A(*),ISR_A(*)
      REAL*8 A(*)
      CHARACTER ERROR*(*)
      LOGICAL PIVOT_A
!     Local Variables


      RETURN
 9999 CALL ERRORS('SOLVER_DUMMY',ERROR)
      RETURN 1
      END


      SUBROUTINE SET_ZERO(A,LDA,N,B,SPARSE_A,ISC_A,ISR_A,SETZERO,ERROR,
     '  *)

C#### Subroutine: SET_ZERO
C###  Description:
C###    Set one variable to zero in a system of equations
C###  Written by Stuart Norris 15/07/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,SPARSE_A,ISC_A(*),ISR_A(*),SETZERO
      REAL*8 A(*),B(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IJ


      IF(SETZERO.GT.N) THEN
        WRITE(ERROR,7000) '>>SETZERO > N (',SETZERO,' > ',N,')'
        GOTO 9999

      ELSE IF(SETZERO.GT.0) THEN
        IF(SPARSE_A.EQ.0) THEN
          DO IJ=SETZERO,LDA*N,LDA
            A(IJ)=0.0D0
          ENDDO
          A(SETZERO+(SETZERO-1)*LDA)=1.0D0
        ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN
          DO IJ=ISR_A(SETZERO),ISR_A(SETZERO+1)-1
            IF(ISC_A(IJ).NE.SETZERO) THEN
              A(IJ)=0.0D0
            ELSE
              A(IJ)=1.0D0
            ENDIF
          ENDDO
        ELSE
          ERROR='>>Sparsity type not implemented'
          GOTO 9999
        ENDIF

        B(SETZERO)=0.0D0
      ENDIF

      RETURN

 7000 FORMAT(A,2(I20,A))

 9999 CALL ERRORS('SET_ZERO',ERROR)
      RETURN 1
      END


      SUBROUTINE KLUDGE_PIV(A,LDA,N,NZA,B,SPARSE_A,ISC_A,ISR_A,LEVEL,
     '  ERROR,*)

C#### Subroutine: 
C###  Description:
C###    A kludgy reordering of the matrix, for testing various full
C###    pivoting strategies.
C###  Written by Stuart Norris 23/05/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,NZA,SPARSE_A,ISC_A(*),ISR_A(*),LEVEL
      REAL*8 A(*),B(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ILOW,IHIGH,I,J,II,IJ,IP,NIJ,JJ
      INTEGER IDX(N),ISR_T(N+1),ISC_T(NZA),JDX(N)
      REAL*8 D(NZA),E(N),BTMP


C     Dense arrays
      IF(SPARSE_A.EQ.0) THEN

        ILOW=1
        IHIGH=N
        DO I=1,N
          II=I+(I-1)*LDA
          IF(A(II).EQ.0.0D0) THEN
            IDX(IHIGH)=I
            IHIGH=IHIGH-1
          ELSE
            IDX(ILOW)=I
            ILOW=ILOW+1
          ENDIF
        ENDDO
        IF(IHIGH-ILOW.EQ.1) THEN
          WRITE(*,*) 'IHIGH-ILOW.NE.1!'
        ENDIF

        DO I=1,N
          E(I)=B(IDX(I))
        ENDDO
        DO J=1,N
          DO I=1,N
            IJ=I+(J-1)*LDA
            IP=IDX(I)+(IDX(J)-1)*LDA
            D(IJ)=A(IP)
          ENDDO
        ENDDO

        CALL DCOPY(N,E,1,B,1)
        CALL DCOPY(NZA,D,1,A,1)

C       Assumes symmetric fill patterns
        IF(LEVEL.GT.1) THEN
          DO I=1,N
            IDX(I)=0
          ENDDO
          DO I=ILOW+1,N
            DO J=1,ILOW
              IJ=I+(J-1)*LDA
              IF(A(IJ).NE.0.0D0 .AND. IDX(I).EQ.0) THEN
                IDX(I)=J
                IDX(J)=I
              ENDIF
            ENDDO
          ENDDO
          DO I=ILOW+1,N
            CALL DSWAP(N,A(I),LDA,A(IDX(I)),LDA)
            BTMP=B(I)
            B(I)=B(IDX(I))
            B(IDX(I))=BTMP
          ENDDO
        ENDIF

C     Sparse arrays
      ELSE IF(SPARSE_A.EQ.1) THEN

        ILOW=1
        IHIGH=N
        DO I=1,N
          DO II=ISR_A(I),ISR_A(I+1)-1
            J=ISC_A(II)
            IF(I.EQ.J) THEN
              IF(A(II).EQ.0.0D0) THEN
                IDX(IHIGH)=I
                JDX(I)=IHIGH
                IHIGH=IHIGH-1
              ELSE
                IDX(ILOW)=I
                JDX(I)=ILOW
                ILOW=ILOW+1
              ENDIF
              GOTO 100
            ENDIF
          ENDDO
          IDX(IHIGH)=I
          JDX(I)=IHIGH
          IHIGH=IHIGH-1
 100      CONTINUE
        ENDDO
        IF(ILOW-IHIGH.NE.1) THEN
          WRITE(*,*) 'ILOW-IHIGH.NE.1!'
        ENDIF

        DO I=1,N
          E(I)=B(IDX(I))
        ENDDO
        CALL DCOPY(N,E,1,B,1)

        ISR_T(1)=1
        DO I=1,N
          NIJ=ISR_A(IDX(I)+1)-ISR_A(IDX(I))
          ISR_T(I+1)=ISR_T(I)+NIJ
          DO J=1,NIJ
            IJ=ISR_T(I)+J-1
            IP=ISR_A(IDX(I))+J-1

            D(IJ)=A(IP)
            ISC_T(IJ)=JDX(ISC_A(IP))
          ENDDO
        ENDDO


        IF(LEVEL.LT.2) THEN
          CALL DCOPY(NZA,D,1,A,1)
          CALL ICOPY(N+1,ISR_T,1,ISR_A,1)
          CALL ICOPY(NZA,ISC_T,1,ISC_A,1)

C       Assumes symmetric fill patterns
        ELSE !IF(LEVEL.GT.1) THEN
          print *, 'kludge x 2'
          DO I=1,N
            IDX(I)=I
          ENDDO
          DO I=ILOW,N

            II=I
            JJ=I

C           DO J=1,ILOW-1
            DO IJ=ISR_T(I),ISR_T(I+1)-1
              J=ISC_T(IJ)
              IF(J.LT.ILOW) THEN
                IF(D(IJ).NE.0.0D0 .AND. IDX(J).EQ.J) THEN
c                  IDX(I)=J
c                  IDX(J)=I
c                  goto 1010
                  jj=j
                  ii=i
                ENDIF
              ENDIF
            ENDDO
            IF(JJ.EQ.I) WRITE(*,*) 'NODIAG: ',I,N
c 1010       continue
            IDX(ii)=jj
            IDX(jj)=ii
c            print *, i,idx(i),idx(idx(i))
          ENDDO

          DO I=ILOW,N
            BTMP=B(I)
            B(I)=B(IDX(I))
            B(IDX(I))=BTMP
          ENDDO

          ISR_A(1)=1
          DO I=1,N
            NIJ=ISR_T(IDX(I)+1)-ISR_T(IDX(I))
            ISR_A(I+1)=ISR_A(I)+NIJ
            DO J=1,NIJ
              IJ=ISR_A(I)+J-1
              IP=ISR_T(IDX(I))+J-1

              A(IJ)=D(IP)
              ISC_A(IJ)=ISC_T(IP)
            ENDDO
          ENDDO

        ENDIF

      ELSE
        ERROR='>>Sparsity type not implemented'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('KLUDGE_PIV',ERROR)
      RETURN 1
      END


      SUBROUTINE Examine(A,LDA,N,B,SPARSE_A,ISC_A,ISR_A,ERROR,*)

C#### Subroutine: Examine
C###  Description:
C###    Randomly weight each equation.
C###  Written by Stuart Norris 22/05/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,SPARSE_A,ISC_A(*),ISR_A(*)
      REAL*8 A(*),B(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,IJ,J,IMAX(N),NMAX(N),IPIV(N),JPIV(N),IPIVMAX,JPIVMAX
      REAL*8 DMAX
!     Subroutines


      DO I=1,N
        IMAX(I)=0
        NMAX(I)=0
        IPIV(I)=0
        JPIV(I)=0
      ENDDO

C     Examine dense arrays
      IF(SPARSE_A.EQ.0) THEN
        DO I=1,N
          DMAX=0.0D0
          DO J=1,N
            IJ=I+(J-1)*LDA
            IF(ABS(A(IJ)).GT.DMAX) THEN
              DMAX=ABS(A(IJ))
              IMAX(I)=J
              NMAX(I)=1
            ELSE IF(ABS(A(IJ)).EQ.DMAX) THEN
              NMAX(I)=NMAX(I)+1
            ENDIF
          ENDDO

          IF(DMAX.EQ.0.0D0) THEN
            IMAX(I)=0
            NMAX(I)=0
          ENDIF

          IF(NMAX(I).EQ.1) THEN
            J=IMAX(I)
            JPIV(J)=JPIV(J)+1
          ELSE IF(NMAX(I).GT.1) THEN
            DO J=1,N
              IJ=I+(J-1)*LDA
              IF(ABS(A(IJ)).EQ.DMAX) THEN
                JPIV(J)=JPIV(J)+1
              ENDIF
            ENDDO
          ENDIF
        ENDDO

C     Examine sparse arrays
      ELSE IF(SPARSE_A.EQ.1) THEN
        DO I=1,N
          DMAX=0.0D0
          DO IJ=ISR_A(I),ISR_A(I+1)-1
            IF(ABS(A(IJ)).GT.DMAX) THEN
              J=ISC_A(IJ)
              DMAX=ABS(A(IJ))
              IMAX(I)=J
              NMAX(I)=1
            ELSE IF(ABS(A(IJ)).EQ.DMAX) THEN
              NMAX(I)=NMAX(I)+1
            ENDIF
          ENDDO

          IF(DMAX.EQ.0.0D0) THEN
            IMAX(I)=0
            NMAX(I)=0
          ENDIF

          IF(NMAX(I).EQ.1) THEN
            J=IMAX(I)
            JPIV(J)=JPIV(J)+1
          ELSE IF(NMAX(I).GT.1) THEN
            DO IJ=ISR_A(I),ISR_A(I+1)-1
              IF(ABS(A(IJ)).EQ.DMAX) THEN
                J=ISC_A(IJ)
                JPIV(J)=JPIV(J)+1
              ENDIF
            ENDDO
          ENDIF
        ENDDO

      ELSE
        ERROR='>>Sparsity type not implemented'
        GOTO 9999
      ENDIF

      DO I=1,N
        IF(IMAX(I).EQ.0) THEN
          WRITE(*,'(A,2(1X,I9))') 'Examine: singular row ',I,N
        ELSE
          IPIV(IMAX(I))=I
        ENDIF
      ENDDO

      IPIVMAX=0
      JPIVMAX=0

      DO I=1,N
C       IF(IPIV(I).EQ.0) THEN
        IF(JPIV(I).EQ.0) THEN
          WRITE(*,'(A,2(1X,I9))') 'Examine: no pivot ',I,N
        ENDIF
        IPIVMAX=MAX(IPIV(I),IPIVMAX)
        JPIVMAX=MAX(JPIV(I),JPIVMAX)
      ENDDO

      WRITE(*,*) 'MAXES: ',IPIVMAX,JPIVMAX

      RETURN

 9999 CALL ERRORS('Examine',ERROR)
      RETURN 1
      END


      SUBROUTINE TRANS_INDEX(N,NZA,ISC,ISR,ITC,ITR,IDX,IWORK,ERROR,*)

C#### Subroutine: TRANS_INDEX
C###  Description:
C###    Generate a compressed column-row index of a compressed row-column 
C###    array. Note that no error checking is done, and so we may end up
C###    stepping off the end of an array.
C###  Written by Stuart Norris 07/05/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NZA
      INTEGER ISC(NZA),ISR(N+1),IDX(NZA),ITC(N+1),ITR(NZA),IWORK(N)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,K,IJ


      DO I=1,N+1
        ITC(I)=0
      ENDDO
      DO I=1,N
        IWORK(I)=0
      ENDDO

      DO I=1,NZA
        J=ISC(I)+1
        IF(J.GT.N+1 .OR. J.LT.2) THEN
          ERROR='>>ERROR in ISC'
          GOTO 9999
        ENDIF
        ITC(J)=ITC(J)+1
      ENDDO

      ITC(1)=1
      DO I=2,N+1
        ITC(I)=ITC(I)+ITC(I-1)
      ENDDO

      DO I=1,N
        DO IJ=ISR(I),ISR(I+1)-1
          J=ISC(IJ)
          IF(J.GT.N .OR. J.LT.1) THEN
            ERROR='>>ERROR in ISC'
            GOTO 9999
          ENDIF
          K=ITC(J)+IWORK(J)
          IF(K.GT.NZA .OR. K.LT.1) THEN
            ERROR='>>ERROR in K'
            GOTO 9999
          ENDIF
          IDX(K)=IJ
          IWORK(J)=IWORK(J)+1

          ITR(K)=I
        ENDDO
      ENDDO

      RETURN

 9999 CALL ERRORS('TRANS_INDEX',ERROR)
      RETURN 1
      END


      SUBROUTINE TRANS_INDEX_SYMM(N,NZA,ISC,ISR,IDX,IWORK,ERROR,*)

C#### Subroutine: TRANS_INDEX_SYMM
C###  Description:
C###    Generate a compressed column-row index of a compressed row-column 
C###    array for a symmetric matrix. Note that no error checking is done,
C###    and so we may end up stepping off the end of an array, not to
C###    mention even worse crap that may occur when the matrix isn't actually
C###    symmetric.
C###  Written by Stuart Norris 10/05/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NZA
      INTEGER ISC(NZA),ISR(N+1),IDX(NZA),IWORK(N)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,K,IJ


      DO I=1,N
        IWORK(I)=0
      ENDDO

      DO I=1,N
        DO IJ=ISR(I),ISR(I+1)-1
          J=ISC(IJ)
          IF(J.GT.N .OR. J.LT.1) THEN
            ERROR='>>ERROR in ISC'
            GOTO 9999
          ENDIF
          K=ISR(J)+IWORK(J)
          IF(K.GT.NZA .OR. K.LT.1) THEN
            ERROR='>>ERROR in K'
            GOTO 9999
          ENDIF
          IDX(K)=IJ
          IWORK(J)=IWORK(J)+1
        ENDDO
      ENDDO

      RETURN

 9999 CALL ERRORS('TRANS_INDEX_SYMM',ERROR)
      RETURN 1
      END


      SUBROUTINE WRITE_BOEING(A,N,NZA,SPARSE_A,ISC_A,ISR_A,FILENM,
     '  ERROR,*)

C#### Subroutine: WRITE_SPARSE
C###  Description:
C###    Write a sparse array out in Boeing format
C###  Written by Stuart Norris 21/08/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NZA,SPARSE_A,ISC_A(*),ISR_A(*)
      REAL*8 A(*)
      CHARACTER*(*) FILENM,ERROR
!     Local constants
      INTEGER UNIT
      PARAMETER(UNIT=10)
!     Local Variables
      INTEGER I,J,IERR,NZ,NP
      INTEGER TOTCRD,PTRCRD,INDCRD,VALCRD,RHSCRD
      INTEGER NROW,NCOL,NNZERO,NELTVL
      CHARACTER TITLE*72,KEY*8,MXTYPE*3
      CHARACTER PTRFMT*20,INDFMT*20,VALFMT*20,RHSFMT*20
!     Functions
      INTEGER LEN_TRIM
      EXTERNAL LEN_TRIM


      CALL ASSERT(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4,
     '  '>>Only sparse matricies dealt with',ERROR,*9999)

      OPEN(UNIT=UNIT,FILE=FILENM,FORM='FORMATTED',STATUS='UNKNOWN',
     '  IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to open '''//FILENM(:LEN_TRIM(FILENM))//''''
        GOTO 9999
      ENDIF

      IF(SPARSE_A.EQ.1) THEN
        WRITE(*,*) 'Dumping sparse array'
        NP=N+1
        NZ=ISR_A(N+1)-1

        TITLE='Matrix from soltest'
        KEY=' '
        WRITE(UNIT,'(A72,A8)') TITLE,KEY

        PTRCRD=NP/10
        IF(MOD(NP,10).NE.0) PTRCRD=PTRCRD+1
        INDCRD=NZ/10
        IF(MOD(NZ,10).NE.0) INDCRD=INDCRD+1
        VALCRD=NZ/5
        IF(MOD(NZ,5).NE.0) VALCRD=VALCRD+1
        RHSCRD=0
        TOTCRD=PTRCRD+INDCRD+VALCRD+RHSCRD+3
        WRITE(UNIT,'(5(I14))') TOTCRD,PTRCRD,INDCRD,VALCRD,RHSCRD

        MXTYPE='RUA'
        NROW=N
        NCOL=N
        NNZERO=NZ
        NELTVL=0
        WRITE(UNIT,'(A3,11X,4(I14))') MXTYPE,NROW,NCOL,NNZERO,NELTVL

        PTRFMT='(10I8)'
        INDFMT='(10I8)'
        VALFMT='(5E16.8)'
        RHSFMT=' '
        WRITE(UNIT,'(2(A16),2(A20))') PTRFMT,INDFMT,VALFMT,RHSFMT

        DO I=1,PTRCRD
          WRITE(UNIT,PTRFMT) (ISR_A(J),J=((I-1)*10+1),MIN(I*10,NP))
        ENDDO

        DO I=1,INDCRD
          WRITE(UNIT,INDFMT) (ISC_A(J),J=((I-1)*10+1),MIN(I*10,NZ))
        ENDDO

        DO I=1,VALCRD
          WRITE(UNIT,VALFMT) (A(J),J=((I-1)*5+1),MIN(I*5,NZ))
        ENDDO

      ENDIF

      CLOSE(UNIT,IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to close '''//FILENM(:LEN_TRIM(FILENM))//''''
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('WRITE_BOEING',ERROR)
      RETURN 1
      END


      SUBROUTINE WRITE_TRIPLET(A,LDA,N,NZA,B,SPARSE_A,ISC_A,ISR_A,
     '  FILENM,ERROR,*)

C#### Subroutine: WRITE_TRIPLET
C###  Description:
C###    Write a sparse array out in triplet form.
C###  Written by Stuart Norris 16/07/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,NZA,SPARSE_A,ISC_A(*),ISR_A(*)
      REAL*8 A(*),B(*)
      CHARACTER*(*) FILENM,ERROR
!     Local constants
      INTEGER UNIT
      PARAMETER(UNIT=10)
!     Local Variables
      INTEGER I,J,IJ,IERR
!     Functions
      INTEGER LEN_TRIM
      EXTERNAL LEN_TRIM


      OPEN(UNIT=UNIT,FILE=FILENM,FORM='FORMATTED',STATUS='UNKNOWN',
     '  IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to open '''//FILENM(:LEN_TRIM(FILENM))//''''
        GOTO 9999
      ENDIF

C     Dense matrix
      IF(SPARSE_A.EQ.0) THEN
        WRITE(*,*) 'Writing dense array'

        WRITE(UNIT,7000,IOSTAT=IERR) N,NZA
        IF(IERR.NE.0) THEN
          ERROR='>>Unable to write N and NZA'
          GOTO 9999
        ENDIF

        DO I=1,N
          DO J=1,N
            WRITE(UNIT,7100,IOSTAT=IERR) I-1,J-1,A(I+(J-1)*LDA)
            IF(IERR.NE.0) THEN
              ERROR='>>Unable to write I,J,A'
              GOTO 9999
            ENDIF
          ENDDO
        ENDDO

C     Compressed row
      ELSE IF(SPARSE_A.EQ.1) THEN
        WRITE(*,*) 'Writing sparse array'

        WRITE(UNIT,7000,IOSTAT=IERR) N,NZA
        IF(IERR.NE.0) THEN
          ERROR='>>Unable to write N and NZA'
          GOTO 9999
        ENDIF

        DO I=1,N
          DO IJ=ISR_A(I),ISR_A(I+1)-1
            WRITE(UNIT,7100,IOSTAT=IERR) I-1,ISC_A(IJ)-1,A(IJ)
            IF(IERR.NE.0) THEN
              ERROR='>>Unable to write I,J,A'
              GOTO 9999
            ENDIF
          ENDDO
        ENDDO

C     Unknown format

      ELSE
        ERROR='>>Unknown sparsity pattern'
        GOTO 9999
      ENDIF

      CLOSE(UNIT,IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to close '''//FILENM(:LEN_TRIM(FILENM))//''''
        GOTO 9999
      ENDIF

      RETURN

 7000 FORMAT(I12,1X,I12)
 7100 FORMAT(2(I12,1X),F10.6)

 9999 CALL ERRORS('WRITE_TRIPLET',ERROR)
      RETURN 1
      END


      SUBROUTINE WRITE_SOLN2D(X,N,NX,NY,FILENM,ERROR,*)

C#### Subroutine: WRITE_SOLN2D
C###  Description:
C###    Write out the solution from a calculation of a 2D FV/FD system.
C###  Written by Stuart Norris 16/07/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NX,NY
      REAL*8 X(2:NX-1,2:NY-1)
      CHARACTER*(*) FILENM,ERROR
!     Local constants
      INTEGER UNIT
      PARAMETER(UNIT=11)
!     Local Variables
      INTEGER I,J,IERR
!     Functions
      INTEGER LEN_TRIM
      EXTERNAL LEN_TRIM


      OPEN(UNIT=UNIT,FILE=FILENM,FORM='FORMATTED',STATUS='UNKNOWN',
     '  IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to open '''//FILENM(:LEN_TRIM(FILENM))//''''
        GOTO 9999
      ENDIF

      WRITE(*,*) 'Writing sparse array'

      ERROR='>>Unable to write I,J,X'
      DO I=1,NX
        WRITE(UNIT,7100,ERR=9999) I,1,-1.0D0
      ENDDO
      WRITE(UNIT,'()',ERR=9999)

      DO J=2,NY-1
        WRITE(UNIT,7100,ERR=9999) 1,J,-1.0D0
        DO I=2,NX-1
          WRITE(UNIT,7100,ERR=9999) I,J,REAL(X(I,J))
        ENDDO
        WRITE(UNIT,7100,ERR=9999) NX,J,1.0D0
        WRITE(UNIT,'()',ERR=9999)
      ENDDO

      DO I=1,NX
        WRITE(UNIT,7100,ERR=9999) I,NY,1.0D0
      ENDDO
      WRITE(UNIT,'()',ERR=9999)
      ERROR=' '

      CLOSE(UNIT,IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to close '''//FILENM(:LEN_TRIM(FILENM))//''''
        GOTO 9999
      ENDIF

      RETURN

 7000 FORMAT(I12,1X,I12)
 7100 FORMAT(2(I12,1X),E12.5E2)

 9999 CALL ERRORS('WRITE_SOLN2D',ERROR)
      RETURN 1
      END


      SUBROUTINE READ_TRIPLET_DIM(N,NZA,FILENM,ERROR,*)

C#### Subroutine: READ_TRIPLET_DIM
C###  Description:
C###    Read the dimensions of a sparse array from a triplet file.
C###  Written by Stuart Norris 16/07/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NZA
      CHARACTER*(*) FILENM,ERROR
!     Local constants
      INTEGER UNIT
      PARAMETER(UNIT=10)
!     Local Variables
      INTEGER IERR
!     Functions
      INTEGER LEN_TRIM
      EXTERNAL LEN_TRIM


      OPEN(UNIT=UNIT,FILE=FILENM,FORM='FORMATTED',STATUS='OLD',
     '  IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to open '''//FILENM(:LEN_TRIM(FILENM))//''''
        GOTO 9999
      ENDIF

      READ(UNIT,*,IOSTAT=IERR) N,NZA
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to read N and NZA'
        GOTO 9999
      ENDIF

      CLOSE(UNIT,IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to close '''//FILENM(:LEN_TRIM(FILENM))//''''
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('READ_TRIPLET_DIM',ERROR)
      RETURN 1
      END


      SUBROUTINE READ_TRIPLET(A,N,NMAX,MMAX,ISC_A,ISR_A,FILENM,ERROR,*)

C#### Subroutine: READ_TRIPLET
C###  Description:
C###    Read a sparse array from a triplet file.
C###  Written by Stuart Norris 16/07/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NMAX,MMAX,ISC_A(*),ISR_A(*)
      REAL*8 A(*)
      CHARACTER*(*) FILENM,ERROR
!     Local constants
      INTEGER UNIT
      PARAMETER(UNIT=10)
!     Local Variables
      INTEGER I,J,IJ,IERR,NZA,IROW
      REAL*8 R
!     Functions
      INTEGER LEN_TRIM
      EXTERNAL LEN_TRIM


      OPEN(UNIT=UNIT,FILE=FILENM,FORM='FORMATTED',STATUS='OLD',
     '  IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to open '''//FILENM(:LEN_TRIM(FILENM))//''''
        GOTO 9999
      ENDIF

      WRITE(*,*) 'Reading sparse array'

      READ(UNIT,*,IOSTAT=IERR) N,NZA
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to read N and NZA'
        GOTO 9999
      ENDIF

      IF(N.GT.NMAX) THEN
        WRITE(ERROR,'(A,2(1X,I5))') '>>N > NMAX ::',N,NMAX
        GOTO 9999
      ENDIF
      IF(NZA.GT.MMAX) THEN
        WRITE(ERROR,'(A,2(1X,I5))') '>>NZA > MMAX ::',NZA,MMAX
        GOTO 9999
      ENDIF

      ISR_A(1)=1
      IROW=1
      DO IJ=1,NZA
        READ(UNIT,*,IOSTAT=IERR) I,J,R
        IF(IERR.NE.0) THEN
          ERROR='>>Unable to write I,J,A'
          GOTO 9999
        ENDIF
        A(IJ)=R
        ISC_A(IJ)=J+1
        IF(I+1.GT.IROW) THEN
          ISR_A(I+1)=IJ
          IROW=I+1
        ENDIF
      ENDDO
      ISR_A(N+1)=NZA+1

      CLOSE(UNIT,IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to close '''//FILENM(:LEN_TRIM(FILENM))//''''
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('READ_TRIPLET',ERROR)
      RETURN 1
      END


      SUBROUTINE WRITE_SPARSE(A,LDA,N,NZA,B,SPARSE_A,ISC_A,ISR_A,FILENM,
     '  ERROR,*)

C#### Subroutine: WRITE_SPARSE
C###  Description:
C###    Write a sparse array out to be viewed by VSM.
C###  Written by Stuart Norris 09/04/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,NZA,SPARSE_A,ISC_A(*),ISR_A(*)
      REAL*8 A(*),B(*)
      CHARACTER*(*) FILENM,ERROR
!     Local constants
      INTEGER UNIT
      PARAMETER(UNIT=10)
!     Local Variables
      INTEGER I,J,IERR,NZ
!     Functions
      INTEGER LEN_TRIM
      EXTERNAL LEN_TRIM


      OPEN(UNIT=UNIT,FILE=FILENM,FORM='UNFORMATTED',STATUS='UNKNOWN',
     '  IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to open '''//FILENM(:LEN_TRIM(FILENM))//''''
        GOTO 9999
      ENDIF

C     Dense matrix
      IF(SPARSE_A.EQ.0) THEN
        WRITE(*,*) 'Dumping dense array'

        WRITE(UNIT,IOSTAT=IERR) N
        IF(IERR.NE.0) THEN
          ERROR='>>Unable to write N'
          GOTO 9999
        ENDIF

        WRITE(UNIT,IOSTAT=IERR) (I*N,I=0,N)
        IF(IERR.NE.0) THEN
          ERROR='>>Unable to write ISR_A'
          GOTO 9999
        ENDIF

        WRITE(UNIT,IOSTAT=IERR) ((I-1,I=1,N),J=1,N)
        IF(IERR.NE.0) THEN
          ERROR='>>Unable to write ISC_A'
          GOTO 9999
        ENDIF
        WRITE(UNIT,IOSTAT=IERR) ((A(I+(J-1)*LDA),J=1,N),I=1,N)
        IF(IERR.NE.0) THEN
          ERROR='>>Unable to write A'
          GOTO 9999
        ENDIF

C     Compressed row
      ELSE IF(SPARSE_A.EQ.1) THEN
        WRITE(*,*) 'Dumping sparse array'
        NZ=ISR_A(N+1)-1

        WRITE(UNIT,IOSTAT=IERR) N
        IF(IERR.NE.0) THEN
          ERROR='>>Unable to write N'
          GOTO 9999
        ENDIF

        WRITE(UNIT,IOSTAT=IERR) (ISR_A(I)-1,I=1,N+1)
        IF(IERR.NE.0) THEN
          ERROR='>>Unable to write ISR_A'
          GOTO 9999
        ENDIF

        WRITE(UNIT,IOSTAT=IERR) (ISC_A(I)-1,I=1,NZ)
        IF(IERR.NE.0) THEN
          ERROR='>>Unable to write ISC_A'
          GOTO 9999
        ENDIF

        WRITE(UNIT,IOSTAT=IERR) (A(I),I=1,NZ)
        IF(IERR.NE.0) THEN
          ERROR='>>Unable to write A'
          GOTO 9999
        ENDIF

C     Row-column sparse arrays
      ELSE IF(SPARSE_A.EQ.2) THEN
        CALL WRITE_SPARSE_RC(A,N,NZA,SPARSE_A,ISC_A,ISR_A,UNIT,
     '    ERROR,*9999)

C     Compressed column sparse arrays
      ELSE IF(SPARSE_A.EQ.3) THEN
        ERROR='>>Unable to cope with Compressed-column arrays'
        GOTO 9999

C     Fancy Row-column sparse arrays
      ELSE IF(SPARSE_A.GE.4 .AND. SPARSE_A.LE.6) THEN
        CALL WRITE_SPARSE_RC(A,N,NZA,SPARSE_A,ISC_A,ISR_A,UNIT,
     '    ERROR,*9999)

C     Unknown format

      ELSE
        ERROR='>>Unknown sparsity pattern'
        GOTO 9999
      ENDIF

      WRITE(UNIT,IOSTAT=IERR) (B(I),I=1,N)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to write B'
        GOTO 9999
      ENDIF

      CLOSE(UNIT,IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to close '''//FILENM(:LEN_TRIM(FILENM))//''''
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('WRITE_SPARSE',ERROR)
      RETURN 1
      END


      SUBROUTINE WRITE_MA28AD_FACTORED(A,LDA,N,B,ISR_A,ISC_LU,ISR_LU,
     '  IKEEP,FILENAME,ERROR,*)

C#### Subroutine: WRITE_MA28AD_FACTORED
C###  Description:
C###    Write a sparse array factorised by Harwell's ma28ad()
C###  Written by Stuart Norris 02/05/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,LDA,ISR_A(*),ISC_LU(*),ISR_LU(*),IKEEP(N,*)
      REAL*8 A(*),B(*)
      CHARACTER FILENAME*(*),ERROR*(*)
!     Local Variables
      INTEGER I,NZA,SPARSE_A
      INTEGER LEN_TRIM
      EXTERNAL LEN_TRIM


      SPARSE_A=1

      ISR_LU(1)=1
      DO I=1,N
        ISR_LU(I+1)=ISR_LU(I)+IKEEP(I,1)
      ENDDO
      NZA=ISR_LU(I)-1

      WRITE(*,*) 'Dumping factored system to ''',
     '  FILENAME(:LEN_TRIM(FILENAME)),''''
      CALL WRITE_SPARSE(A,LDA,N,NZA,B,SPARSE_A,ISC_LU,ISR_LU,FILENAME,
     '  ERROR,*9999)

      RETURN

 9999 CALL ERRORS('WRITE_MA28AD_FACTORED',ERROR)
      RETURN 1
      END


      SUBROUTINE WRITE_SPARSE_RC(A,N,NZ,SPARSE_A,ISC_A,ISR_A,UNIT,
     '  ERROR,*)

C#### Subroutine: 
C###  Description:
C###    Sort a sparse array, so that it can be written out in compressed
C###    row format.
C###  Written by Stuart Norris 18/04/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NZ,SPARSE_A,ISC_A(*),ISR_A(*),UNIT
      REAL*8 A(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,K,IERR,OFFSET
      INTEGER LSR_A(N+1),IDX(NZ)


C     Row-column sparsity
      IF(SPARSE_A.EQ.2 .OR. SPARSE_A.EQ.4 .OR. SPARSE_A.EQ.6) THEN
        K=1
        DO I=1,N
          LSR_A(I)=K
          DO J=1,NZ
            IF(ISR_A(J).EQ.I) THEN
              IDX(K)=J
              K=K+1
            ENDIF
          ENDDO
        ENDDO
        LSR_A(N+1)=K
        OFFSET=0

C     Row-column sparsity, with all ROW+COL data in ISC_A
      ELSE IF(SPARSE_A.EQ.5) THEN
        K=1
        DO I=1,N
          LSR_A(I)=K
          DO J=1,NZ
            IF(ISC_A(J).EQ.I) THEN
              IDX(K)=J
              K=K+1
            ENDIF
          ENDDO
        ENDDO
        LSR_A(N+1)=K
        OFFSET=NZ

C     Compressed column data
C     ELSE IF(SPARSE_A.EQ.3) THEN

      ELSE
        ERROR='>>Unknown sparsity pattern'
        GOTO 9999
      ENDIF

C     Write out the data

      WRITE(UNIT,IOSTAT=IERR) N
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to write N'
        GOTO 9999
      ENDIF

      WRITE(UNIT,IOSTAT=IERR) (LSR_A(I)-1,I=1,N+1)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to write ISR_A'
        GOTO 9999
      ENDIF

      WRITE(UNIT,IOSTAT=IERR) (ISC_A(IDX(I)+OFFSET)-1,I=1,NZ)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to write ISC_A'
        GOTO 9999
      ENDIF

      WRITE(UNIT,IOSTAT=IERR) (A(IDX(I)),I=1,NZ)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to write A'
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('WRITE_SPARSE_RC',ERROR)
      RETURN 1
      END


      SUBROUTINE READ_SPARSE_DIM(N,NZA,FILENM,ERROR,*)

C#### Subroutine: READ_SPARSE_DIM
C###  Description:
C###    Read in the dimensions of a sparse array written by WRITE_SPARSE
C###  Written by Stuart Norris 11/04/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NZA
      CHARACTER ERROR*(*),FILENM*(*)
!     Local constants
      INTEGER UNIT
      PARAMETER(UNIT=10)
!     Local Variables
      INTEGER I,IERR,ISR_A
!     Functions
      INTEGER LEN_TRIM
      EXTERNAL LEN_TRIM


      OPEN(UNIT=UNIT,FILE=FILENM,FORM='UNFORMATTED',STATUS='OLD',
     '  IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to open '''//FILENM(:LEN_TRIM(FILENM))//''''
        GOTO 9999
      ENDIF

      READ(UNIT,IOSTAT=IERR) N
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to read N'
        GOTO 9999
      ENDIF

      READ(UNIT,IOSTAT=IERR) (ISR_A,I=1,N+1)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to read ISR_A'
        GOTO 9999
      ENDIF
      NZA=ISR_A

      CLOSE(UNIT,IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to close '''//FILENM(:LEN_TRIM(FILENM))//''''
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('READ_SPARSE_DIM',ERROR)
      RETURN 1
      END


      SUBROUTINE READ_SPARSE(A,N,NMAX,MMAX,B,ISC_A,ISR_A,FILENM,ERROR,*)

C#### Subroutine: READ_SPARSE
C###  Description:
C###    Read in a sparse array written by WRITE_SPARSE
C###  Written by Stuart Norris 10/04/02

      IMPLICIT NONE
!     Parameter List
      INTEGER N,NMAX,MMAX,ISC_A(*),ISR_A(*)
      REAL*8 A(*),B(*)
      CHARACTER ERROR*(*),FILENM*(*)
!     Local constants
      INTEGER UNIT
      PARAMETER(UNIT=10)
!     Local Variables
      INTEGER I,IJ,IERR,NZ,ISCTMP(MMAX)
      REAL*8 ATMP(MMAX),BTMP(NMAX)
!     Functions
      INTEGER LEN_TRIM
      EXTERNAL LEN_TRIM


      OPEN(UNIT=UNIT,FILE=FILENM,FORM='UNFORMATTED',STATUS='OLD',
     '  IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to open '''//FILENM(:LEN_TRIM(FILENM))//''''
        GOTO 9999
      ENDIF

      READ(UNIT,IOSTAT=IERR) N
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to read N'
        GOTO 9999
      ENDIF
      IF(N.GT.NMAX) THEN
        WRITE(ERROR,'(A,2(1X,I5))') '>>N > NMAX ::',N,NMAX
        GOTO 9999
      ENDIF

c      READ(UNIT,IOSTAT=IERR) (ISR_A(I),I=1,N+1)
      READ(UNIT,IOSTAT=IERR) (ISCTMP(I),I=1,N+1)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to read ISR_A'
        GOTO 9999
      ENDIF
c      DO I=1,N+1
c        ISR_A(I)=ISR_A(I)+1
c      ENDDO
C$OMP PARALLEL DO
      DO I=1,N+1
        ISR_A(I)=ISCTMP(I)+1
      ENDDO
C$OMP END PARALLEL DO

      NZ=ISR_A(N+1)-1
      IF(NZ.GT.MMAX) THEN
        WRITE(ERROR,'(A,2(1X,I5))') '>>NZ > MMAX ::',NZ,MMAX
        GOTO 9999        
      ENDIF

      READ(UNIT,IOSTAT=IERR) (ISCTMP(I),I=1,NZ)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to read ISC_A'
        GOTO 9999
      ENDIF
c      DO I=1,NZ
c        ISC_A(I)=ISC_A(I)+1
c      ENDDO
C$OMP PARALLEL DO
      DO I=1,N
        DO IJ=ISR_A(I),ISR_A(I+1)-1
          ISC_A(IJ)=ISCTMP(IJ)+1
        ENDDO
      ENDDO
C$OMP END PARALLEL DO

      READ(UNIT,IOSTAT=IERR) (ATMP(I),I=1,NZ)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to read A'
        GOTO 9999
      ENDIF
C$OMP PARALLEL DO
      DO I=1,N
        DO IJ=ISR_A(I),ISR_A(I+1)-1
          A(IJ)=ATMP(IJ)
        ENDDO
      ENDDO
C$OMP END PARALLEL DO

      READ(UNIT,IOSTAT=IERR) (BTMP(I),I=1,N)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to read B'
        GOTO 9999
      ENDIF
C$OMP PARALLEL DO
      DO I=1,N
        B(I)=BTMP(I)
      ENDDO
C$OMP END PARALLEL DO

      CLOSE(UNIT,IOSTAT=IERR)
      IF(IERR.NE.0) THEN
        ERROR='>>Unable to close '''//FILENM(:LEN_TRIM(FILENM))//''''
        GOTO 9999
      ENDIF

      RETURN

 9999 CALL ERRORS('READ_SPARSE',ERROR)
      RETURN 1
      END


      SUBROUTINE WRITE_SYSTEM(TRANS,UNIT,A,B,LDA,N,L,M,SPARSE_A,PIVOT_A,
     '  ISC_A,ISR_A,ERROR,*)

C#### Subroutine: WRITE_SYSTEM
C###  Description:
C###    Writes a system of equations to the specified unit.
C###  Written by Stuart Norris 21/03/02

      IMPLICIT NONE
!     Parameter List
      INTEGER UNIT,LDA,N,M,L,SPARSE_A,ISC_A(*),ISR_A(*)
      REAL*8 A(*),B(*)
      CHARACTER ERROR*(*),TRANS*(*)
      LOGICAL PIVOT_A
!     Local Variables
      INTEGER I,J,K
      REAL*8 X(N)
      CHARACTER*32 FORM*(32),TYPE*(7),ORIENT*(17)


      CALL ASSERT(SPARSE_A.EQ.0 .OR. SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4,
     '  '>>Sparsity type not implemented',ERROR,*9999)

      IF(SPARSE_A.EQ.0) THEN
        TYPE='  dense'
      ELSE
        TYPE=' sparse'
      ENDIF
      IF(TRANS(1:1).EQ.'T') THEN
        ORIENT=' transpose system'
      ELSE
        ORIENT=' system'
      ENDIF

      WRITE(UNIT,'()')
      WRITE(UNIT,'(4(A,I3),2A)') 'System: ',L,' x ',M,' of a ',N,' x ',
     '  N,TYPE,ORIENT
      WRITE(UNIT,'()')
C     WRITE(FORM,'(A,I3,A)') '(I3,1P,',M,'(1X,E12.5),A,E12.5)'
      WRITE(FORM,'(A,I3,A)') '(I3,1P,',M,'(1X,E9.2),A,E9.2)'

      DO I=1,L

C       Non-sparse matricies, with and without pivoting
        IF(SPARSE_A.EQ.0) THEN
          IF(TRANS(1:1).EQ.'T') THEN
            CALL DCOPY(N,A((I-1)*LDA+1),1,X,1)
          ELSE
            CALL DCOPY(N,A(I),LDA,X,1)
          ENDIF

C       Sparse matricies: CMISS compressed row
        ELSE IF(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) THEN
          DO J=1,N
            X(J)=0.0D0
          ENDDO
          IF(TRANS(1:1).EQ.'T') THEN
            DO K=1,N
C             CALL GET_IJ_MAT(J,K,I,LDA,N,SPARSE_A,PIVOT_A,IPIVOT_A,
C    '          ISC_A,ISR_A,ERROR,*9999)
C             IF(J.NE.0) X(K)=A(J)
              DO J=ISR_A(K),ISR_A(K+1)-1
                IF(ISC_A(J).EQ.I) THEN
                  X(K)=A(J)
                  GOTO 100
                ENDIF
              ENDDO
 100          CONTINUE
            ENDDO
          ELSE
            DO J=ISR_A(I),ISR_A(I+1)-1
              X(ISC_A(J))=A(J)
            ENDDO
          ENDIF

        ELSE
          ERROR='>>Sparsity type not implemented'
          GOTO 9999
        ENDIF

        WRITE(UNIT,FORM) I,(X(J),J=1,M),' : ',B(I)
      ENDDO
      WRITE(UNIT,'()')

      RETURN

 9999 CALL ERRORS('WRITE_SYSTEM',ERROR)
      RETURN 1
      END


      SUBROUTINE SPARSE_TO_DENSE(A,D,LDA,N,SPARSE_A,ISC_A,ISR_A,ERROR,*)

C#### Subroutine: SPARSE_TO_DENSE
C###  Description:
C###    Convert the sparse matrix in A() into a dense matrix in D.
C###  Written by Stuart Norris 11/04/02

      IMPLICIT NONE
!     Parameter List
      INTEGER LDA,N,SPARSE_A,ISC_A(*),ISR_A(*)
      REAL*8 A(*),D(LDA,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,IJ


      IF(SPARSE_A.EQ.0) THEN
        ERROR='>>Trying to convert a dense matrix to a dense matrix'
        GOTO 9999
      ENDIF

C     test case A
C$OMP PARALLEL DO PRIVATE(I,J)
      DO J=1,N
        DO I=1,N
          D(I,J)=0.0D0
        ENDDO
      ENDDO
C$OMP END PARALLEL DO

cC     test case B
c      DO J=1,N
cC$OMP PARALLEL DO PRIVATE(I)
c        DO I=1,N
c          D(I,J)=0.0D0
c        ENDDO
cC$OMP END PARALLEL DO
c      ENDDO

C$OMP PARALLEL DO PRIVATE(I,J,IJ)
      DO I=1,N
        DO IJ=ISR_A(I),ISR_A(I+1)-1
          J=ISC_A(IJ)
          D(I,J)=A(IJ)
        ENDDO
      ENDDO
C$OMP END PARALLEL DO

      RETURN

 9999 CALL ERRORS('SPARSE_TO_DENSE',ERROR)
      RETURN 1
      END
