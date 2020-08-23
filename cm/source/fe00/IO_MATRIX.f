      SUBROUTINE IO_MATRIX(COMAND,DEVICE,ISC_MATRIX,ISR_MATRIX,
     '  IUNIT,MATRIX_M,MATRIX_N,MATRIX_MAXM,NZTOT,MMAX,NMAX,NZMAX,
     '  ROWLIST,SPARSITY_TYPE,
     '  MATRIX,FILEFORMAT,INOUTTYPE,ALLROWS,ERROR,*)

C#### Subroutine: IO_MATRIX
C###  Description:
C###    IO_MATRIX reads and writes a single matrix, allowing for the
C###    way the matrix is stored (ie sparsely, and what scheme).

C***  Note: The matrices will be read/writen according to the
C***        FILEFORMAT and IOOUTTYPE. In addition the matrix sparsity
C***        (if required)
C***  MATRIX_M    = Number of rows in the matrix
C***  MATRIX_N    = Number of columns in the matrix
C***  MATRIX_MAXM = The leading dimension of the Matrix

      IMPLICIT NONE
      INCLUDE 'binf00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'mach00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER ISC_MATRIX(*),ISR_MATRIX(*),IUNIT,MATRIX_M,MATRIX_N,
     '  MATRIX_MAXM,MMAX,NMAX,NZMAX,NZTOT,ROWLIST(0:*),SPARSITY_TYPE
      REAL*8 MATRIX(*)
      CHARACTER COMAND*(*),DEVICE*(*),ERROR*(*),
     '  FILEFORMAT*(*),INOUTTYPE*(*)
      LOGICAL ALLROWS
!     Local Variables
      INTEGER CERROR(50),col,col_loop,column,ERR,NUMINDICES(1),
     '  NUMVALUES(1),nz,NZEROTOT(1),row,row_loop,SPARSITYCODE(1),
     '  SPARSEARRAYSIZE,TIMECODE(1)
      REAL*8 DUMMY_VALUE,KBYTESAVING,PERCENTSAVING,TIME(1)
      CHARACTER FMT*500,
     '  LINE*132,MODE*4,SPARSITY*6,SPARSITY_STRING*18
      LOGICAL ISBINFILEOPEN


C!!! LKC 17-JAN-2000
C!!!   There is a problem with the same variable (NZ_DUMMY)
C!!!   being passed in as NZMAX and NZTOT

      CALL ENTERS('IO_MATRIX',*9999)

C***  Reading

      IF(COMAND(1:4).EQ.'READ') THEN
C       Must read from disk
C!!!! AJP 26-3-97.  It appears that no check has been made on the
C!!!! size of the array in which the matrix is to be stored (only the
C!!!! leading dimension is pasted in so a full check cannot be made
C!!!! here).
        CALL ASSERT(DEVICE(1:4).EQ.'DISK','Invalid device for read',
     '    ERROR,*9999)
        IF(FILEFORMAT(1:5).EQ.'ASCII') THEN

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C          READ(IUNIT,'('' Number of rows='',I5,'
C     '      //''', Number of columns='',I5)') MATRIX_M,MATRIX_N
          FMT='('' Number of rows='',I5,'
     '      //''', Number of columns='',I5)'
          READ(IUNIT,FMT) MATRIX_M,MATRIX_N

          CALL ASSERT(MATRIX_M.LE.MMAX,'ERROR, increase M dimension',
     '      ERROR,*9999)
          CALL ASSERT(MATRIX_N.LE.NMAX,'ERROR, increase N dimension',
     '      ERROR,*9999)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C          READ(IUNIT,'('' Number of values = '',I10)') NZTOT
C          READ(IUNIT,'('' Number of values stored = '',I10)') NZTOT
          FMT='('' Number of values = '',I10)'
          READ(IUNIT,FMT) NZTOT
          FMT='('' Number of values stored = '',I10)'
          READ(IUNIT,FMT) NZTOT

          CALL ASSERT(NZTOT.LE.NZMAX,'ERROR, increase NZMAX',
     '      ERROR,*9999)

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C          READ(IUNIT,'('' Matrix is stored as a '',A6,'' matrix'')')
C     '      SPARSITY
          FMT='('' Matrix is stored as a '',A6,'' matrix'')'
          READ(IUNIT,FMT) SPARSITY

C         Read in sparsity patterns
          IF(SPARSITY(1:6).EQ.'sparse') THEN
            READ(IUNIT,'(1X,A)') LINE
            IF(LINE(1:4).EQ.'Both') THEN
              MODE='BOTH'
            ELSE IF(LINE(1:4).EQ.'Only') THEN
              MODE='ONLY'
            ELSE
              ERROR='>>Invalid File Format'
              GOTO 9999
            ENDIF

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C            READ(IUNIT,'('' Sparsity pattern used is '',A17)')
C     '        SPARSITY_STRING
            FMT='('' Sparsity pattern used is '',A17)'
            READ(IUNIT,FMT) SPARSITY_STRING
            IF(SPARSITY_STRING(1:14).EQ.'Compressed Row') THEN
              SPARSITY_TYPE=1
              READ(IUNIT,'(A)') LINE
              READ(IUNIT,'(A)') LINE

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C              READ(IUNIT,'(1X,7(1X,I10),:/(1X,7(1X,I10)))')
C     '          (ISR_MATRIX(nz),nz=1,MATRIX_M+1)
              FMT='(1X,7(1X,I10),:/(1X,7(1X,I10)))'
              READ(IUNIT,FMT) (ISR_MATRIX(nz),nz=1,MATRIX_M+1)
              READ(IUNIT,'(A)') LINE
              CALL ASSERT(NZTOT.EQ.ISR_MATRIX(MATRIX_M+1)-1,
     '          '>>Invalid row indicies',ERROR,*9999)

              DO row_loop=1,MATRIX_M

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(IUNIT,'('' Row '',I7,'':'',8(1X,I7),'
C     '            //':/(13X,8(1X,I7)))') row,(ISC_MATRIX(nz),
C     '            nz=ISR_MATRIX(row_loop),ISR_MATRIX(row_loop+1)-1)
                FMT='('' Row '',I7,'':'',8(1X,I7),:/(13X,8(1X,I7)))'
                READ(IUNIT,FMT) row,(ISC_MATRIX(nz),
     '            nz=ISR_MATRIX(row_loop),ISR_MATRIX(row_loop+1)-1)
              ENDDO !row_loop
            ELSE IF(SPARSITY_STRING(1:10).EQ.'Row Column') THEN
              SPARSITY_TYPE=2
              READ(IUNIT,'(A)') LINE
              READ(IUNIT,'(A)') LINE

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C              READ(IUNIT,'(1X,9(1X,I7),:/(1X,9(1X,I7)))')
C     '          (ISR_MATRIX(nz),nz=1,NZTOT)
C              READ(IUNIT,'(A)') LINE
C              READ(IUNIT,'(1X,9(1X,I7),:/(1X,9(1X,I7)))')
C     '          (ISC_MATRIX(nz),nz=1,NZTOT)

              FMT='(1X,9(1X,I7),:/(1X,9(1X,I7)))'
              READ(IUNIT,FMT) (ISR_MATRIX(nz),nz=1,NZTOT)
              READ(IUNIT,'(A)') LINE
              READ(IUNIT,'(1X,9(1X,I7),:/(1X,9(1X,I7)))')
     '          (ISC_MATRIX(nz),nz=1,NZTOT)


            ELSE IF(SPARSITY_STRING(1:17).EQ.'Compressed Column') THEN
              SPARSITY_TYPE=3
              READ(IUNIT,'(A)') LINE
              READ(IUNIT,'(A)') LINE
              READ(IUNIT,'(1X,7(1X,I10),:/(1X,7(1X,I10)))')
     '          (ISC_MATRIX(nz),nz=1,MATRIX_N+1)
              READ(IUNIT,'(A)') LINE
              CALL ASSERT(NZTOT.EQ.ISC_MATRIX(MATRIX_N+1)-1,
     '          '>>Invalid column indicies',ERROR,*9999)
              DO col_loop=1,MATRIX_N
C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(IUNIT,'('' Col '',I7,'':'',8(1X,I7),'
C     '            //':/(13X,8(1X,I7)))') col,(ISR_MATRIX(nz),
C     '            nz=ISC_MATRIX(col_loop),ISC_MATRIX(col_loop+1)-1)
                FMT='('' Col '',I7,'':'',8(1X,I7),:/(13X,8(1X,I7)))'
                READ(IUNIT,FMT) col,(ISR_MATRIX(nz),
     '            nz=ISC_MATRIX(col_loop),ISC_MATRIX(col_loop+1)-1)
              ENDDO !col_loop
            ELSE IF(SPARSITY_STRING(1:17).EQ.'Sorted Row Column') THEN
              SPARSITY_TYPE=4
              READ(IUNIT,'(A)') LINE
              READ(IUNIT,'(A)') LINE
              READ(IUNIT,'(1X,9(1X,I7),:/(1X,9(1X,I7)))')
     '          (ISR_MATRIX(nz),nz=1,NZTOT)
              READ(IUNIT,'(A)') LINE
              READ(IUNIT,'(1X,9(1X,I7),:/(1X,9(1X,I7)))')
     '          (ISC_MATRIX(nz),nz=1,NZTOT)
            ELSE IF(SPARSITY_STRING(1:18).EQ.'Umfpack Row Column') THEN
              SPARSITY_TYPE=5
              READ(IUNIT,'(A)') LINE
              READ(IUNIT,'(A)') LINE
              READ(IUNIT,'(1X,9(1X,I7),:/(1X,9(1X,I7)))')
     '          (ISC_MATRIX(nz),nz=1,2*NZTOT)
            ELSE
              ERROR='>>Invalid sparsity type'
              GOTO 9999
            ENDIF !sparsity_type
          ELSE
C***        If the matrix is stored in the file as a full matrix then
C***        how to store it in memory when it is read back in is not
C***        determined (could be stored in column major (fortran) format
C***        or as a continous block of memory (as in sparseness=0). To
C***        resolve this the way it is stored in memory is passed into
C***        IO_MATRIX in the SPARSITY_TYPE variable. If SPARSITY_TYPE
C***        =-1 then the matrix is stored in the normal fortran column
C***        major mode (ie. MATRIX_MAXM is used to index). If the type
C***        =0 the block storage is used (ie. MATRIX_M is used to index)
C***        If neither is true block storage defaults.
CC AJPs 191297
C old            IF(SPARSITY_TYPE.NE.-1.OR.SPARSITY_TYPE.NE.0) THEN
            IF(SPARSITY_TYPE.NE.-1.AND.SPARSITY_TYPE.NE.0) THEN
CC AJPe
              SPARSITY_TYPE=0
            ENDIF
            MODE='BOTH'
          ENDIF


C*** Read in array values
          IF(MODE(1:4).EQ.'ONLY'.AND.INOUTTYPE.NE.'SPARSITY') THEN
            WRITE(OP_STRING,'('' >>Warning: Can not read matrix '
     '        //'values'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(MODE(1:4).EQ.'BOTH'.AND.INOUTTYPE.EQ.'SPARSITY') THEN
C***        Skip over values
            READ(IUNIT,'(A)') LINE
            IF(SPARSITY_TYPE.EQ.-1.OR.SPARSITY_TYPE.EQ.0) THEN
              DO row_loop=1,MATRIX_M

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(IUNIT,'('' Row '',I7,'':'',4D16.8,'
C     '            //':/(13X,4D16.8))') row,(DUMMY_VALUE,col=1,MATRIX_N)
                FMT='('' Row '',I7,'':'',4D16.8,:/(13X,4D16.8))'
                READ(IUNIT,FMT) row,(DUMMY_VALUE,col=1,MATRIX_N)
              ENDDO !row_loop
            ELSE IF(SPARSITY_TYPE.EQ.1) THEN
              DO row_loop=1,MATRIX_M

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(IUNIT,'('' Row '',I7,'':'',4D16.8,'
C     '            //':/(13X,4D16.8))') row,(DUMMY_VALUE,
C     '            nz=ISR_MATRIX(row_loop),ISR_MATRIX(row_loop+1)-1)
                FMT='('' Row '',I7,'':'',4D16.8,:/(13X,4D16.8))'
                READ(IUNIT,FMT) row,(DUMMY_VALUE,
     '            nz=ISR_MATRIX(row_loop),ISR_MATRIX(row_loop+1)-1)
              ENDDO !row_loop
            ELSE IF(SPARSITY_TYPE.EQ.2.OR.SPARSITY_TYPE.EQ.4) THEN
              READ(IUNIT,'(1X,5D16.8,:/(1X,5D16.8))') (DUMMY_VALUE,
     '          nz=1,NZTOT)
            ELSE IF(SPARSITY_TYPE.EQ.3) THEN
              DO col_loop=1,MATRIX_N

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(IUNIT,'('' Col '',I7,'':'',4D16.8,'
C     '            //':/(13X,4D16.8))') col,(DUMMY_VALUE,
C     '            nz=ISC_MATRIX(col_loop),ISC_MATRIX(col_loop+1)-1)
                FMT='('' Col '',I7,'':'',4D16.8,:/(13X,4D16.8))'
                READ(IUNIT,FMT) col,(DUMMY_VALUE,
     '            nz=ISC_MATRIX(col_loop),ISC_MATRIX(col_loop+1)-1)
              ENDDO !col_loop
            ELSE IF(SPARSITY_TYPE.EQ.5) THEN
              READ(IUNIT,'(1X,5D16.8,:/(1X,5D16.8))') (DUMMY_VALUE,
     '          nz=1,NZTOT)
            ENDIF
          ELSE IF(MODE(1:4).EQ.'BOTH'.AND.INOUTTYPE.NE.'SPARSITY') THEN
            READ(IUNIT,'(A)') LINE
            IF(SPARSITY_TYPE.EQ.-1) THEN ! No sparsity (fortran)
              DO row_loop=1,MATRIX_M

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(IUNIT,'('' Row '',I7,'':'',4D16.8,'
C     '            //':/(13X,4D16.8))') row,(MATRIX(row+(col-1)*
C     '            MATRIX_MAXM),col=1,MATRIX_N)
                FMT='('' Row '',I7,'':'',4D16.8,'
     '            //':/(13X,4D16.8))'
                READ(IUNIT,FMT) row,(MATRIX(row+(col-1)*
     '            MATRIX_MAXM),col=1,MATRIX_N)
              ENDDO !row_loop

            ELSE IF(SPARSITY_TYPE.EQ.0) THEN !No sparsity (compact)
              DO row_loop=1,MATRIX_M

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(IUNIT,'('' Row '',I7,'':'',4D16.8,'
C     '            //':/(13X,4D16.8))') row,(MATRIX(row+(col-1)*
C     '            MATRIX_M),col=1,MATRIX_N)
                FMT='('' Row '',I7,'':'',4D16.8,:/(13X,4D16.8))'

C LKC 18-APR-2002 is this really wrong?? It appears so anyway,
C   and the MFI array reads in wrong but all other ones seem ok...
C                READ(IUNIT,FMT) row,(MATRIX(row+(col-1)*
C     '            MATRIX_M),col=1,MATRIX_N)
                READ(IUNIT,FMT) row,(MATRIX(row+(col-1)*
     '            MATRIX_MAXM),col=1,MATRIX_N)

              ENDDO !row_loop
            ELSE IF(SPARSITY_TYPE.EQ.1) THEN !Compressed Row
              DO row_loop=1,MATRIX_M

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(IUNIT,'('' Row '',I7,'':'',4D16.8,'
C     '            //':/(13X,4D16.8))') row,(MATRIX(nz),
C     '            nz=ISR_MATRIX(row_loop),ISR_MATRIX(row_loop+1)-1)
                FMT='('' Row '',I7,'':'',4D16.8,:/(13X,4D16.8))'
                READ(IUNIT,FMT) row,(MATRIX(nz),
     '            nz=ISR_MATRIX(row_loop),ISR_MATRIX(row_loop+1)-1)

              ENDDO !row_loop
            ELSE IF(SPARSITY_TYPE.EQ.2.OR.SPARSITY_TYPE.EQ.4) THEN !Row-Column
              READ(IUNIT,'(1X,5D16.8,:/(1X,5D16.8))') (MATRIX(nz),
     '          nz=1,NZTOT)
            ELSE IF(SPARSITY_TYPE.EQ.3) THEN !Compressed Column
              DO col_loop=1,MATRIX_N

C LKC 26-JUN-2001 This is invalid FORTRAN - but now works
C                READ(IUNIT,'('' Col '',I7,'':'',4D16.8,'
C     '            //':/(13X,4D16.8))') col,(MATRIX(nz),
C     '            nz=ISC_MATRIX(col_loop),ISC_MATRIX(col_loop+1)-1)
                FMT='('' Col '',I7,'':'',4D16.8,:/(13X,4D16.8))'
                READ(IUNIT,FMT) col,(MATRIX(nz),
     '            nz=ISC_MATRIX(col_loop),ISC_MATRIX(col_loop+1)-1)
              ENDDO !col_loop
            ELSE IF(SPARSITY_TYPE.EQ.5) THEN
              !Row-Column
              READ(IUNIT,'(1X,5D16.8,:/(1X,5D16.8))') (MATRIX(nz),
     '          nz=1,NZTOT)
            ELSE
              ERROR='>>Invalid sparsity type'
              GOTO 9999
            ENDIF
          ENDIF !array/all_inout

        ELSE IF(FILEFORMAT(1:6).EQ.'BINARY') THEN

C         Read number of indicies
          CALL BINREADFILE(IUNIT,INTTYPE,1,NUMINDICES,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          CALL ASSERT(NUMINDICES(1).EQ.2,
     '      '>>Number of indices <> 2, not implemented',ERROR,*9999)
C         Read the index values
          CALL BINREADFILE(IUNIT,INTTYPE,NUMINDICES(1),INTDATA,
     '      REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          MATRIX_M=INTDATA(1)
          MATRIX_N=INTDATA(2)
C         Read the size of the matrix

C LKC 17-JAN=2000 Why are we setting this value before and after?
          NZEROTOT(1)=NZTOT
          CALL BINREADFILE(IUNIT,INTTYPE,1,NZEROTOT,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          NZTOT=NZEROTOT(1)
          CALL ASSERT(MATRIX_M.LE.MMAX,'ERROR, increase M dimension',
     '      ERROR,*9999)
          CALL ASSERT(MATRIX_N.LE.NMAX,'ERROR, increase N dimension',
     '      ERROR,*9999)
C LKC 17-JAN-2000 Wrong error message?
C          CALL ASSERT(NZTOT.LE.NZMAX,'ERROR, increase N dimension',
C     '      ERROR,*9999)
          CALL ASSERT(NZTOT.LE.NZMAX,'ERROR, increase NZ Max',
     '      ERROR,*9999)


C         Read the number of values stored
          CALL BINREADFILE(IUNIT,INTTYPE,1,NUMVALUES,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          IF((INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '      INOUTTYPE(1:9).EQ.'ALL_INOUT').AND.NUMVALUES(1).EQ.0) THEN
            ERROR='>>File does not contain array values'
            GOTO 9999
          ENDIF
C         Read the timecode
          CALL BINREADFILE(IUNIT,INTTYPE,1,TIMECODE,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          IF(TIMECODE(1).EQ.1) THEN
C           Read the time
            CALL BINREADFILE(IUNIT,DPTYPE,1,INTDATA,REAL4DATA,
     '        TIME,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          ENDIF
C         Read the sparsity type
          CALL BINREADFILE(IUNIT,INTTYPE,1,SPARSITYCODE,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          IF(SPARSITYCODE(1).GT.0) THEN
            SPARSITY_TYPE=SPARSITYCODE(1)
          ELSE
C***        If the matrix is stored in the file as a full matrix then
C***        how to store it in memory when it is read back in is not
C***        determined (could be stored in column major (fortran) format
C***        or as a continous block of memory (as in sparseness=0). To
C***        resolve this the way it is stored in memory is passed into
C***        IO_MATRIX in the SPARSITY_TYPE variable. If SPARSITY_TYPE
C***        =-1 then the matrix is stored in the normal fortran column
C***        major mode (ie. MATRIX_MAXM is used to index). If the type
C***        =0 the block storage is used (ie. MATRIX_M is used to index)
C***        If neither is true block storage defaults.
CC AJPs            IF(SPARSITY_TYPE.NE.-1.OR.SPARSITY_TYPE.NE.0) THEN
            IF(SPARSITY_TYPE.NE.-1.AND.SPARSITY_TYPE.NE.0) THEN
CC AJPe

              SPARSITY_TYPE=0
            ENDIF
          ENDIF
          IF(SPARSITY_TYPE.GT.0) THEN
C           Read the sparsity information
            IF(INOUTTYPE(1:6).EQ.'ARRAYS') THEN
              ERROR='>>Cannot just read the values for a sparse matrix'
              GOTO 9999
            ENDIF
            IF(SPARSITY_TYPE.EQ.1) THEN !compressed row
C             Read the row offsets
              CALL BINREADFILE(IUNIT,INTTYPE,MATRIX_M+1,ISR_MATRIX,
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '          *9999)
C             Read the column numbers
              CALL BINREADFILE(IUNIT,INTTYPE,ISR_MATRIX(MATRIX_M+1)-1,
     '          ISC_MATRIX,REAL4DATA,REAL8DATA,CHARDATA,
     '          LOGDATA,SINTDATA,ERROR,*9999)
              NZTOT=ISR_MATRIX(MATRIX_M+1)-1
            ELSE IF(SPARSITY_TYPE.EQ.2.OR.SPARSITY_TYPE.EQ.4) THEN !row-column
C             Read the number of non-zeros
              NZEROTOT(1)=NZTOT
              CALL BINREADFILE(IUNIT,INTTYPE,1,NZEROTOT,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              NZTOT=NZEROTOT(1)
C             Read the row numbers
              CALL BINREADFILE(IUNIT,INTTYPE,NZTOT,ISR_MATRIX,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C             Read the column numbers
              CALL BINREADFILE(IUNIT,INTTYPE,NZTOT,ISC_MATRIX,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
            ELSE IF(SPARSITY_TYPE.EQ.3) THEN !compressed column
C             Read the column offsets
              CALL BINREADFILE(IUNIT,INTTYPE,MATRIX_N+1,ISC_MATRIX,
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '          *9999)
C             Read the row numbers
              CALL BINREADFILE(IUNIT,INTTYPE,ISC_MATRIX(MATRIX_N+1)-1,
     '          ISR_MATRIX,REAL4DATA,REAL8DATA,CHARDATA,
     '          LOGDATA,SINTDATA,ERROR,*9999)
              NZTOT=ISC_MATRIX(MATRIX_N+1)-1
            ELSE IF(SPARSITY_TYPE.EQ.5) THEN
              !umfpack row-column
C             Read the number of non-zeros
              NZEROTOT(1)=NZTOT
              CALL BINREADFILE(IUNIT,INTTYPE,1,NZEROTOT,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              NZTOT=NZEROTOT(1)
C             Read the row and column numbers
              CALL BINREADFILE(IUNIT,INTTYPE,2*NZTOT,ISC_MATRIX,
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
            ELSE
              ERROR='>>Unknown sparsity format'
              GOTO 9999
            ENDIF
          ENDIF
          IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '      INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
C           Read the values
            IF(SPARSITY_TYPE.EQ.-1) THEN
              CALL ASSERT(NUMVALUES(1).EQ.MATRIX_M*MATRIX_N,
     '          '>>Invalid number of matrix values',ERROR,*9999)
              DO column=1,MATRIX_N
                CALL BINREADFILE(IUNIT,DPTYPE,MATRIX_M,INTDATA,
     '            REAL4DATA,MATRIX((column-1)*MATRIX_MAXM+1),
     '            CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              ENDDO !column
            ELSE IF(SPARSITY_TYPE.EQ.0) THEN
              CALL BINREADFILE(IUNIT,DPTYPE,NUMVALUES(1),INTDATA,
     '          REAL4DATA,MATRIX,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
            ELSE IF(SPARSITY_TYPE.EQ.1) THEN
              CALL BINREADFILE(IUNIT,DPTYPE,NUMVALUES(1),INTDATA,
     '          REAL4DATA,MATRIX,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
            ELSE IF(SPARSITY_TYPE.EQ.2.OR.SPARSITY_TYPE.EQ.4) THEN
              CALL BINREADFILE(IUNIT,DPTYPE,NUMVALUES(1),INTDATA,
     '          REAL4DATA,MATRIX,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
            ELSE IF(SPARSITY_TYPE.EQ.3) THEN
              CALL BINREADFILE(IUNIT,DPTYPE,NUMVALUES(1),INTDATA,
     '          REAL4DATA,MATRIX,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
            ELSE IF(SPARSITY_TYPE.EQ.5) THEN
              CALL BINREADFILE(IUNIT,DPTYPE,NUMVALUES(1),INTDATA,
     '          REAL4DATA,MATRIX,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
            ENDIF
          ELSE
            IF(NUMVALUES(1).GT.0) THEN
C             Skip the values
              CALL BINSKIPFILE(IUNIT,NUMVALUES(1)*DPSIZE,ERROR,*9999)
            ENDIF
          ENDIF

        ENDIF

C***  Writing

      ELSE IF(COMAND.EQ.'WRITE') THEN

        IF(FILEFORMAT(1:5).EQ.'ASCII') THEN

          IF(DEVICE(1:4).EQ.'DISK') THEN

            WRITE(IUNIT,'('' Number of rows='',I5,'
     '        //''', Number of columns='',I5)') MATRIX_M,
     '        MATRIX_N
CC AJPs 191297
C old            WRITE(IUNIT,'('' Number of values = '',I10)') NZTOT
            IF(SPARSITY_TYPE.EQ.-1) THEN
              WRITE(IUNIT,'('' Number of values = '',I10)')
     '          MATRIX_M*MATRIX_N
              WRITE(IUNIT,'('' Number of values stored = '',I10)')
     '          MATRIX_M*MATRIX_N
            ELSE
              WRITE(IUNIT,'('' Number of values = '',I10)') NZTOT
              WRITE(IUNIT,'('' Number of values stored = '',I10)') NZTOT
            ENDIF
CC AJPe
            IF(SPARSITY_TYPE.GT.0) THEN
              WRITE(IUNIT,'('' Matrix is stored as a sparse '
     '          //'matrix'')')
              IF(INOUTTYPE(1:6).EQ.'ARRAYS') THEN
C***            It is an error to just write out the values of a sparse
C***            matrix without the sparsity pattern as the matrix will
C***            not be able to be read back in again.
                ERROR='>>Cannot just write the values of a sparse '
     '            //'matrix'
                GOTO 9999
              ENDIF
              IF(INOUTTYPE(1:8).EQ.'SPARSITY') THEN
                WRITE(IUNIT,'('' Only sparsity is contained'')')
              ELSE
                WRITE(IUNIT,'('' Both sparsity and values are '
     '            //'contained'')')
              ENDIF
              IF(SPARSITY_TYPE.EQ.1) THEN
                WRITE(IUNIT,'('' Sparsity pattern used is '
     '            //'Compressed Row'')')
C               Write out sparsity patterns
                WRITE(IUNIT,'('' Sparsity indices:'')')
                WRITE(IUNIT,'('' Row offsets:'')')
                CALL WRITE_LONG(INTTYPE,1,1,IUNIT,MATRIX_M+1,
     '            7,7,ISR_MATRIX,%VAL(0),'(1X,7(1X,I10))',
     '            '(1X,7(1X,I10))',ERROR,*9999)
                WRITE(IUNIT,'('' Column numbers:'')')
                DO row=1,MATRIX_M
                  WRITE(LINE,'(I7)') row
                  CALL WRITE_LONG(INTTYPE,ISR_MATRIX(row),1,IUNIT,
     '              ISR_MATRIX(row+1)-1,8,8,ISC_MATRIX,%VAL(0),
     '              '('' Row '//LINE(1:7)//':'',8(1X,I7))',
     '              '(13X,8(1X,I7))',
     '              ERROR,*9999)
                ENDDO !row
              ELSE IF(SPARSITY_TYPE.EQ.2) THEN
                WRITE(IUNIT,'('' Sparsity pattern used is '
     '            //'Row Column    '')')
C               Write out sparsity patterns
                WRITE(IUNIT,'('' Sparsity indices:'')')
                WRITE(IUNIT,'('' Row numbers:'')')
                CALL WRITE_LONG(INTTYPE,1,1,IUNIT,NZTOT,9,9,
     '            ISR_MATRIX,%VAL(0),'(1X,9(1X,I7))','(1X,9(1X,I7))',
     '            ERROR,*9999)
                WRITE(IUNIT,'('' Column numbers:'')')
                CALL WRITE_LONG(INTTYPE,1,1,IUNIT,NZTOT,
     '            9,9,ISC_MATRIX,%VAL(0),
     '            '(1X,9(1X,I7))','(1X,9(1X,I7))',
     '            ERROR,*9999)
              ELSE IF(SPARSITY_TYPE.EQ.3) THEN
                WRITE(IUNIT,'('' Sparsity pattern used is '
     '            //'Compressed Column'')')
C               Write out sparsity patterns
                WRITE(IUNIT,'('' Sparsity indices:'')')
                WRITE(IUNIT,'('' Column offsets:'')')
                CALL WRITE_LONG(INTTYPE,1,1,IUNIT,MATRIX_N+1,
     '            7,7,ISC_MATRIX,%VAL(0),'(1X,7(1X,I10))',
     '           '(1X,7(1X,I10))',ERROR,*9999)
                WRITE(IUNIT,'('' Row numbers:'')')
                DO column=1,MATRIX_N
                  WRITE(LINE,'(I7)') column
                  CALL WRITE_LONG(INTTYPE,ISC_MATRIX(column),1,IUNIT,
     '              ISC_MATRIX(column+1)-1,8,8,ISR_MATRIX,%VAL(0),
     '              '('' Col '//LINE(1:7)//':'',8(1X,I7))',
     '              '(13X,8(1X,I7))',
     '              ERROR,*9999)
                ENDDO !row
              ELSE IF(SPARSITY_TYPE.EQ.4) THEN
                WRITE(IUNIT,'('' Sparsity pattern used is '
     '            //'Sorted Row Column'')')
C               Write out sparsity patterns
                WRITE(IUNIT,'('' Sparsity indices:'')')
                WRITE(IUNIT,'('' Row numbers:'')')
                CALL WRITE_LONG(INTTYPE,1,1,IUNIT,NZTOT,9,9,
     '            ISR_MATRIX,%VAL(0),'(1X,9(1X,I7))','(1X,9(1X,I7))',
     '            ERROR,*9999)
                WRITE(IUNIT,'('' Column numbers:'')')
                CALL WRITE_LONG(INTTYPE,1,1,IUNIT,NZTOT,
     '            9,9,ISC_MATRIX,%VAL(0),
     '            '(1X,9(1X,I7))','(1X,9(1X,I7))',
     '            ERROR,*9999)
              ELSE IF(SPARSITY_TYPE.EQ.5) THEN
                WRITE(IUNIT,'('' Sparsity pattern used is '
     '            //'Umfpack Row Column    '')')
C               Write out sparsity patterns
                WRITE(IUNIT,'('' Sparsity indices:'')')
                WRITE(IUNIT,'('' Row numbers:'')')
                CALL WRITE_LONG(INTTYPE,1,1,IUNIT,2*NZTOT,9,9,
     '            ISC_MATRIX,%VAL(0),'(1X,9(1X,I7))','(1X,9(1X,I7))',
     '            ERROR,*9999)
              ENDIF
            ELSE
              WRITE(IUNIT,'('' Matrix is stored as a full   '
     '          //'matrix'')')
            ENDIF

C           Write out array values
            IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '        INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
              WRITE(IUNIT,'('' Matrix values:'')')
              IF(SPARSITY_TYPE.EQ.-1) THEN
                DO row=1,MATRIX_M
                  WRITE(LINE,'(I7)') row
                  CALL WRITE_LONG(DPTYPE,row,MATRIX_MAXM,IUNIT,
     '              row+(MATRIX_N-1)*MATRIX_MAXM,
     '              4,4,%VAL(0),MATRIX,
     '              '('' Row '//LINE(1:7)//':'',4D16.8)',
     '              '(13X,4D16.8))',
     '              ERROR,*9999)
                ENDDO !row
              ELSE IF(SPARSITY_TYPE.EQ.0) THEN
                DO row=1,MATRIX_M
                  WRITE(LINE,'(I7)') row
                  CALL WRITE_LONG(DPTYPE,row,MATRIX_M,IUNIT,
     '              row+(MATRIX_N-1)*MATRIX_M,
     '              4,4,%VAL(0),MATRIX,
     '              '('' Row '//LINE(1:7)//':'',4D16.8)',
     '              '(13X,4D16.8))',
     '              ERROR,*9999)
                ENDDO !row
              ELSE IF(SPARSITY_TYPE.EQ.1) THEN !compressed row
                DO row=1,MATRIX_M
                  WRITE(LINE,'(I7)') row
                  CALL WRITE_LONG(DPTYPE,ISR_MATRIX(row),1,
     '              IUNIT,ISR_MATRIX(row+1)-1,4,4,%VAL(0),MATRIX,
     '              '('' Row '//LINE(1:7)//':'',4D16.8)',
     '              '(13X,4D16.8))',
     '              ERROR,*9999)
                ENDDO !row
              ELSE IF(SPARSITY_TYPE.EQ.2.OR.SPARSITY_TYPE.EQ.4) THEN !rowcolumn
                CALL WRITE_LONG(DPTYPE,1,1,IUNIT,
     '            NZTOT,
     '            5,5,%VAL(0),MATRIX,
     '            '(1X,5D16.8)',
     '            '(1X,5D16.8)',
     '            ERROR,*9999)
              ELSE IF(SPARSITY_TYPE.EQ.3) THEN !compressed column
                DO column=1,MATRIX_N
                  WRITE(LINE,'(I7)') column
                  CALL WRITE_LONG(DPTYPE,ISC_MATRIX(column),1,
     '              IUNIT,ISC_MATRIX(column+1)-1,4,4,%VAL(0),MATRIX,
     '              '('' Col '//LINE(1:7)//':'',4D16.8)',
     '              '(13X,4D16.8))',
     '              ERROR,*9999)
                ENDDO !row
              ELSE IF(SPARSITY_TYPE.EQ.5) THEN !rowcolumn
                CALL WRITE_LONG(DPTYPE,1,1,IUNIT,NZTOT,5,5,%VAL(0),
     '            MATRIX,'(1X,5D16.8)','(1X,5D16.8)',ERROR,*9999)
              ELSE
                ERROR='>>Invalid sparsity type'
                GOTO 9999
              ENDIF             !sparsity
            ENDIF               !array/all_inout

          ELSE IF(DEVICE(1:4).EQ.'TERM') THEN

C           Write out sparsity patterns
            WRITE(OP_STRING,'('' Number of rows='',I5,'
     '        //''', Number of columns='',I5)') MATRIX_M,MATRIX_N
            CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
            IF(SPARSITY_TYPE.GT.0) THEN
              IF(SPARSITY_TYPE.EQ.1) THEN !compressed row
                WRITE(OP_STRING,'('' Matrix is stored as a compressed '
     '            //'row sparse matrix'')')
                CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                SPARSEARRAYSIZE=ISR_MATRIX(MATRIX_M+1)-1
              ELSE IF(SPARSITY_TYPE.EQ.2) THEN !row-column
                WRITE(OP_STRING,'('' Matrix is stored as a row-column '
     '            //'sparse matrix'')')
                CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                SPARSEARRAYSIZE=2*NZTOT
              ELSE IF(SPARSITY_TYPE.EQ.3) THEN !compressed column
                WRITE(OP_STRING,'('' Matrix is stored as a compressed '
     '            //'column sparse matrix'')')
                CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                SPARSEARRAYSIZE=ISC_MATRIX(MATRIX_N+1)-1
              ELSE IF(SPARSITY_TYPE.EQ.4) THEN !row-column
                WRITE(OP_STRING,'('' Matrix is stored as a sorted '
     '            //'row-column sparse matrix'')')
                CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                SPARSEARRAYSIZE=2*NZTOT
              ELSE IF(SPARSITY_TYPE.EQ.5) THEN !umfpack row-column
                WRITE(OP_STRING,'('' Matrix is stored as a Umfpack '
     '            //'row-column sparse matrix'')')
                CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                SPARSEARRAYSIZE=2*NZTOT
              ENDIF
              WRITE(OP_STRING,'('' Sparsity array size = '',I10)')
     '          SPARSEARRAYSIZE
              CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Number of non-zeros = '',I10)')
     '          NZTOT
              CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Full matrix size    = '',I10)')
     '          MATRIX_M*MATRIX_N
              CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Percentage sparsity =    '',F6.2,'
     '          //'''%'')') DBLE(NZTOT)/DBLE(MATRIX_M*MATRIX_N)*100.0d0
              CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
              KBYTESAVING=DBLE(8*MATRIX_M*MATRIX_N-4*SPARSEARRAYSIZE-
     '          8*NZTOT)/DBLE(1024)
              PERCENTSAVING=DBLE(1024.0d0*KBYTESAVING)/
     '          DBLE(8*MATRIX_M*MATRIX_N)*100.0d0
              WRITE(OP_STRING,'('' kilobytes saved     = '',F10.3)')
     '          KBYTESAVING
              CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Percentage saving   =    '',F6.1,'
     '          //'''%'')') PERCENTSAVING
              CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
            ELSE
              WRITE(OP_STRING,'('' Matrix is stored as a full '
     '          //'matrix'')')
              CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
            ENDIF

            IF(INOUTTYPE(1:8).EQ.'SPARSITY'.OR.
     '        INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
              IF(SPARSITY_TYPE.GT.0) THEN
                WRITE(OP_STRING,'(/'' Sparsity indices:'')')
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(SPARSITY_TYPE.EQ.1) THEN
                WRITE(OP_STRING,'('' Row offsets:'')')
                CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                IF(ALLROWS) THEN
                  CALL WRITE_LONG(INTTYPE,1,1,IUNIT,MATRIX_M+1,
     '              7,7,ISR_MATRIX,%VAL(0),
     '              '(1X,7(1X,I10))','(1X,7(1X,I10))',ERROR,*9999)
                ELSE
C GMH 31/5/96  Cannot split up - too complicated
                  WRITE(OP_STRING,'(1X,7(1X,I10),:/(1X,7(1X,I10)))')
     '              (ISR_MATRIX(ROWLIST(nz)),nz=1,ROWLIST(0))
                  CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                ENDIF
                WRITE(OP_STRING,'('' Column numbers:'')')
                CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                IF(ALLROWS) THEN
                  DO row=1,MATRIX_M
                    WRITE(LINE,'(I7)') row
                    CALL WRITE_LONG(INTTYPE,ISR_MATRIX(row),1,IUNIT,
     '                ISR_MATRIX(row+1)-1,8,8,ISC_MATRIX,%VAL(0),
     '                '('' Row '//LINE(1:7)//':'',8(1X,I7))',
     '                '(13X,8(1X,I7))',ERROR,*9999)
                  ENDDO !row
                ELSE
                  DO row_loop=1,ROWLIST(0)
                    row=ROWLIST(row_loop)
                    WRITE(OP_STRING,'('' Row '',I7,'':'',8(1X,I7),'
     '                //':/(13X,8(1X,I7)))') row,(ISC_MATRIX(nz),
     '                nz=ISR_MATRIX(row),ISR_MATRIX(row+1)-1)
                    CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                  ENDDO !row
                ENDIF
              ELSE IF(SPARSITY_TYPE.EQ.2.OR.SPARSITY_TYPE.EQ.4) THEN
                WRITE(OP_STRING,'('' Row numbers:'')')
                CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                CALL WRITE_LONG(INTTYPE,1,1,IUNIT,NZTOT,8,8,
     '            ISR_MATRIX,%VAL(0),
     '            '(1X,8(1X,I7))','(1X,8(1X,I7))',ERROR,*9999)
                WRITE(OP_STRING,'('' Column numbers:'')')
                CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                CALL WRITE_LONG(INTTYPE,1,1,IUNIT,
     '            NZTOT,
     '            8,8,ISC_MATRIX,%VAL(0),
     '            '(1X,8(1X,I7))',
     '            '(1X,8(1X,I7))',
     '            ERROR,*9999)
              ELSE IF(SPARSITY_TYPE.EQ.3) THEN
                WRITE(OP_STRING,'('' Column offsets:'')')
                CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                IF(ALLROWS) THEN
                  CALL WRITE_LONG(INTTYPE,1,1,IUNIT,MATRIX_N+1,
     '              7,7,ISC_MATRIX,%VAL(0),
     '              '(1X,7(1X,I10))','(1X,7(1X,I10))',
     '              ERROR,*9999)
                ELSE
                  WRITE(OP_STRING,'(1X,7(1X,I10),:/(1X,7(1X,I10)))')
     '              (ISC_MATRIX(ROWLIST(nz)),nz=1,ROWLIST(0))
                  CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                ENDIF
                WRITE(OP_STRING,'('' Row numbers:'')')
                CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                IF(ALLROWS) THEN
                  DO column=1,MATRIX_N
                    WRITE(LINE,'(I7)') column
                    CALL WRITE_LONG(INTTYPE,ISC_MATRIX(column),1,
     '                IUNIT,ISC_MATRIX(column+1)-1,8,8,ISR_MATRIX,
     '                %VAL(0),'('' Col '//LINE(1:7)//':'',8(1X,I7))',
     '                '(13X,8(1X,I7))',ERROR,*9999)
                  ENDDO !column
                ELSE
                  DO col_loop=1,ROWLIST(0)
                    column=ROWLIST(col_loop)
                    WRITE(OP_STRING,'('' Col '',I7,'':'',8(1X,I7),'
     '                //':/(13X,8(1X,I7)))') column,(ISR_MATRIX(nz),
     '                nz=ISC_MATRIX(column),ISC_MATRIX(column+1)-1)
                    CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                  ENDDO !row
                ENDIF
              ELSE IF(SPARSITY_TYPE.EQ.5) THEN
                WRITE(OP_STRING,'('' Row and column numbers:'')')
                CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                CALL WRITE_LONG(INTTYPE,1,1,IUNIT,2*NZTOT,8,8,
     '            ISR_MATRIX,%VAL(0),'(1X,8(1X,I7))',
     '            '(1X,8(1X,I7))',ERROR,*9999)
                WRITE(OP_STRING,'('' Column numbers:'')')
              ENDIF
            ENDIF

C           Write out array values
            IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '        INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
              WRITE(OP_STRING,'(/'' Matrix values:'')')
              CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
              IF(SPARSITY_TYPE.EQ.-1) THEN
                IF(ALLROWS) THEN
                  DO row=1,MATRIX_M
                    WRITE(LINE,'(I7)') row
                    CALL WRITE_LONG(DPTYPE,row,MATRIX_MAXM,IUNIT,
     '                row+(MATRIX_N-1)*MATRIX_MAXM,
     '                5,5,%VAL(0),MATRIX,
     '                '('' Row '//LINE(1:7)//':'',5D12.4)',
     '                '(13X,5D12.4)',
     '                ERROR,*9999)
                  ENDDO !row
                ELSE
                  DO row_loop=1,ROWLIST(0)
                    row=ROWLIST(row_loop)
                    WRITE(OP_STRING,'('' Row '',I7,'':'',5D12.4,'
     '                //':/(13X,5D12.4))') row,
     '                (MATRIX(row+(col-1)*MATRIX_MAXM),col=1,MATRIX_N)
                    CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                  ENDDO !row
                ENDIF
              ELSE IF(SPARSITY_TYPE.EQ.0) THEN !not sparse
                IF(ALLROWS) THEN
                  DO row=1,MATRIX_M
                    WRITE(LINE,'(I7)') row
                    CALL WRITE_LONG(DPTYPE,row,MATRIX_M,IUNIT,
     '                row+(MATRIX_N-1)*MATRIX_M,
     '                5,5,%VAL(0),MATRIX,
     '                '('' Row '//LINE(1:7)//':'',5D12.4)',
     '                '(13X,5D12.4))',
     '                ERROR,*9999)
                  ENDDO !row
                ELSE
                  DO row_loop=1,ROWLIST(0)
                    row=ROWLIST(row_loop)
                    WRITE(OP_STRING,'('' Row '',I7,'':'',5D12.4,'
     '                //':/(13X,5D12.4))') row,
     '                (MATRIX(row+(col-1)*MATRIX_M),col=1,MATRIX_N)
                    CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                  ENDDO !row_loop
                ENDIF
              ELSE IF(SPARSITY_TYPE.EQ.1) THEN !compressed row
                IF(ALLROWS) THEN
                  DO row=1,MATRIX_M
                    WRITE(LINE,'(I7)') row
                    CALL WRITE_LONG(DPTYPE,ISR_MATRIX(row),1,
     '                IUNIT,ISR_MATRIX(row+1)-1,5,5,%VAL(0),MATRIX,
     '                '('' Row '//LINE(1:7)//':'',5D12.4)',
     '                '(13X,5D12.4))',ERROR,*9999)
                  ENDDO !row
                ELSE
                  DO row_loop=1,ROWLIST(0)
                    row=ROWLIST(row_loop)
                    WRITE(OP_STRING,'('' Row '',I7,'':'',5D12.4,'
     '                //':/(13X,5D12.4))') row,(MATRIX(nz),
     '                nz=ISR_MATRIX(row),ISR_MATRIX(row+1)-1)
                    CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                  ENDDO         !row_loop
                ENDIF
              ELSE IF(SPARSITY_TYPE.EQ.2.
     '            OR.SPARSITY_TYPE.EQ.4) THEN !row-column
                CALL WRITE_LONG(DPTYPE,1,1,IUNIT,
     '            NZTOT,
     '            5,5,%VAL(0),MATRIX,
     '            '(1X,5(1X,D11.4))',
     '            '(1X,5(1X,D11.4))',
     '            ERROR,*9999)
              ELSE IF(SPARSITY_TYPE.EQ.3) THEN !compressed column
                IF(ALLROWS) THEN
                  DO column=1,MATRIX_N
                    WRITE(LINE,'(I7)') column
                    CALL WRITE_LONG(DPTYPE,ISC_MATRIX(column),1,
     '                IUNIT,ISC_MATRIX(column+1)-1,5,5,%VAL(0),MATRIX,
     '                '('' Col '//LINE(1:7)//':'',5D12.4)',
     '                '(13X,5D12.4))',ERROR,*9999)
                  ENDDO !column
                ELSE
                  DO col_loop=1,ROWLIST(0)
                    column=ROWLIST(col_loop)
                    WRITE(OP_STRING,'('' Col '',I7,'':'',5D12.4,'
     '                //':/(13X,5D12.4))') column,(MATRIX(nz),
     '                nz=ISC_MATRIX(column),ISC_MATRIX(column+1)-1)
                    CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
                  ENDDO !col_loop
                ENDIF
              ELSE IF(SPARSITY_TYPE.EQ.5) THEN !row-column
                CALL WRITE_LONG(DPTYPE,1,1,IUNIT,NZTOT,5,5,
     '            %VAL(0),MATRIX,
     '            '(1X,5(1X,D11.4))','(1X,5(1X,D11.4))',ERROR,*9999)
              ELSE
                ERROR='>>Invalid sparsity type'
                GOTO 9999
              ENDIF !sparsity
            ENDIF !array/all_inout
          ENDIF

        ELSE IF(FILEFORMAT(1:6).EQ.'BINARY') THEN

          CALL ASSERT(DEVICE(1:4).EQ.'DISK',
     '      '>>Invalid device for binary write',ERROR,*9999)

C         Write the number of indicies
          INTDATA(1)=2
          CALL BINWRITEFILE(IUNIT,INTTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C         Write the index values
          INTDATA(1)=MATRIX_M
          INTDATA(2)=MATRIX_N
          CALL BINWRITEFILE(IUNIT,INTTYPE,2,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C         Write out the size of the matrix
          IF(SPARSITY_TYPE.EQ.-1) THEN
            NZTOT=MATRIX_M*MATRIX_N
          ENDIF
          NZEROTOT(1)=NZTOT
          CALL BINWRITEFILE(IUNIT,INTTYPE,1,NZEROTOT,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          NZTOT=NZEROTOT(1)
C         Write the number of values stored
          IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '      INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
            NUMVALUES(1)=NZTOT
          ELSE
            NUMVALUES(1)=0
          ENDIF
          CALL BINWRITEFILE(IUNIT,INTTYPE,1,NUMVALUES,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C         Write the timecode (not implemented at the moment)
          INTDATA(1)=0
          CALL BINWRITEFILE(IUNIT,INTTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C         Write the time (for future development)
C         Write the sparsity format
          IF(SPARSITY_TYPE.EQ.-1) THEN
            INTDATA(1)=0
          ELSE
            INTDATA(1)=SPARSITY_TYPE
          ENDIF
          CALL BINWRITEFILE(IUNIT,INTTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
C         Write the sparsity information
          IF(SPARSITY_TYPE.GT.0) THEN
            IF(SPARSITY_TYPE.EQ.1) THEN
C             Write the row offsets
              CALL BINWRITEFILE(IUNIT,INTTYPE,MATRIX_M+1,ISR_MATRIX,
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '          *9999)
C             Write the column numbers
              CALL BINWRITEFILE(IUNIT,INTTYPE,NZTOT,
     '          ISC_MATRIX,REAL4DATA,REAL8DATA,
     '          CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
            ELSE IF(SPARSITY_TYPE.EQ.2.OR.SPARSITY_TYPE.EQ.4) THEN
C             Write the number of non zeros
              NZEROTOT(1)=NZTOT
              CALL BINWRITEFILE(IUNIT,INTTYPE,1,NZEROTOT,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              NZTOT=NZEROTOT(1)
C             Write the row numbers
              CALL BINWRITEFILE(IUNIT,INTTYPE,NZTOT,ISR_MATRIX,
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
C             Write the column numbers
              CALL BINWRITEFILE(IUNIT,INTTYPE,NZTOT,ISC_MATRIX,
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
            ELSE IF(SPARSITY_TYPE.EQ.3) THEN
C             Write the column offsets
              CALL BINWRITEFILE(IUNIT,INTTYPE,MATRIX_N+1,ISC_MATRIX,
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,
     '          *9999)
C             Write the row numbers
              CALL BINWRITEFILE(IUNIT,INTTYPE,NZTOT,
     '          ISR_MATRIX,REAL4DATA,REAL8DATA,
     '          CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
            ELSE IF(SPARSITY_TYPE.EQ.5) THEN
C             Write the number of non zeros
              NZEROTOT(1)=NZTOT
              CALL BINWRITEFILE(IUNIT,INTTYPE,1,NZEROTOT,REAL4DATA,
     '          REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              NZTOT=NZEROTOT(1)
C             Write the row and column numbers
              CALL BINWRITEFILE(IUNIT,INTTYPE,2*NZTOT,ISC_MATRIX,
     '          REAL4DATA,REAL8DATA,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
            ELSE
              ERROR='>>Invalid sparsity type'
              GOTO 9999
            ENDIF
          ENDIF
C         Write the values
          IF(INOUTTYPE(1:6).EQ.'ARRAYS'.OR.
     '      INOUTTYPE(1:9).EQ.'ALL_INOUT') THEN
            IF(SPARSITY_TYPE.EQ.-1) THEN
              DO column=1,MATRIX_N
                CALL BINWRITEFILE(IUNIT,DPTYPE,MATRIX_M,INTDATA,
     '            REAL4DATA,MATRIX((column-1)*MATRIX_MAXM+1),
     '            CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
              ENDDO !column
            ELSE IF(SPARSITY_TYPE.EQ.0) THEN
              CALL BINWRITEFILE(IUNIT,DPTYPE,NZTOT,INTDATA,
     '          REAL4DATA,MATRIX,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
            ELSE IF(SPARSITY_TYPE.EQ.1) THEN
              CALL BINWRITEFILE(IUNIT,DPTYPE,NZTOT,INTDATA,
     '          REAL4DATA,MATRIX,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
            ELSE IF(SPARSITY_TYPE.EQ.2.OR.SPARSITY_TYPE.EQ.4) THEN
              CALL BINWRITEFILE(IUNIT,DPTYPE,NZTOT,INTDATA,
     '          REAL4DATA,MATRIX,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
            ELSE IF(SPARSITY_TYPE.EQ.3) THEN
              CALL BINWRITEFILE(IUNIT,DPTYPE,NZTOT,INTDATA,
     '          REAL4DATA,MATRIX,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
            ELSE IF(SPARSITY_TYPE.EQ.5) THEN
              CALL BINWRITEFILE(IUNIT,DPTYPE,NZTOT,INTDATA,
     '          REAL4DATA,MATRIX,CHARDATA,LOGDATA,SINTDATA,
     '          ERROR,*9999)
            ELSE
              ERROR='>>Invalid sparsity type'
              GOTO 9999
            ENDIF
          ENDIF

        ENDIF
      ENDIF !read/write

      CALL EXITS('IO_MATRIX')
      RETURN
 9999 CALL ERRORS('IO_MATRIX',ERROR)
      CALL EXITS('IO_MATRIX')
      IF(ISBINFILEOPEN(IUNIT)) CALL BINARYCLOSEFILE(IUNIT,ERR,CERROR)
      RETURN 1
      END


