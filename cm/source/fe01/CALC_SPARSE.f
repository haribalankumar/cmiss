      SUBROUTINE CALC_SPARSE(ISCMAX,ISRMAX,ISC,ISR,M,N,NZTOT,
     '  SPARSENESS,WORK_ARRAY,ERROR,*)

C#### Subroutine: CALC_SPARSE
C###  Description:
C###    CALC_SPARSE calculates sparsity patterns from a image of
C###    the sparsity pattern from the work array

C KAT 2Feb99: Transposed WORK_ARRAY to reduce page swapping.

      IMPLICIT NONE
!     Parameter List
      INTEGER ISCMAX,ISRMAX,ISC(ISCMAX),ISR(ISRMAX),M,N,NZTOT,
     '  SPARSENESS
      LOGICAL*1 WORK_ARRAY(N,M)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,nzz

      CALL ENTERS('CALC_SPARSE',*9999)

      IF(SPARSENESS.EQ.1) THEN !Compressed row format

C LKC 16-JAN-1999 Change to error statements for additional info output
C        CALL ASSERT((N+1).LE.ISRMAX,'>>Increase ISRMAX',ERROR,*9999)
        IF(N+1.GT.ISRMAX) THEN
          WRITE(ERROR,'(''Increase ISRMAX ( '',I8,'')'')') ISRMAX
          GOTO 9999
        ENDIF

CC       Initalise sparsity arrays
CCC$DOACROSS local(i) share(M, ISR)
C        DO i=1,M+1
C          ISR(i)=1
C        ENDDO !i
CCC$DOACROSS local(nz) share(ISCMAX, ISC)
C        DO nz=1,ISCMAX
C          ISC(nz)=0
C        ENDDO !nz
C       Calculate sparsity pattern
        nzz=1
        DO i=1,M
          ISR(i)=nzz
          DO j=1,N
            IF(WORK_ARRAY(j,i)) THEN
C***          The element is set so add it to the sparsity pattern
              IF(nzz.LE.ISCMAX) THEN
C***            Set the column numbers
                ISC(nzz)=j
              ENDIF
              nzz=nzz+1
            ENDIF
          ENDDO !j
        ENDDO !i
        ISR(M+1)=nzz
        NZTOT=nzz-1

C LKC 16-JAN-1999 Change to error statements for additional info output
C        CALL ASSERT(NZTOT.LE.ISCMAX,'>>Increase ISCMAX',ERROR,*9999)
        IF(NZTOT.GT.ISCMAX) THEN
          WRITE(ERROR,'(''Increase ISCMAX ( '',I8,'')'')') ISCMAX
          GOTO 9999
        ENDIF

      ELSE IF(SPARSENESS.EQ.2.OR.SPARSENESS.EQ.4) THEN !row-column format
CC       Initalise sparsity arrays
CCC$DOACROSS local(nz) share(ISRMAX, ISR)
C        DO nz=1,ISRMAX
C          ISR(nz)=0
C        ENDDO !nz
CCC$DOACROSS local(nz) share(ISCMAX, ISC)
C        DO nz=1,ISCMAX
C          ISC(nz)=0
C        ENDDO !nz
C       Calculate sparsity pattern
        nzz=0
        DO i=1,M
          DO j=1,N
            IF(WORK_ARRAY(j,i)) THEN
C***          The element is set so add it to the sparsity pattern
              nzz=nzz+1
              IF(nzz.LE.ISRMAX) THEN
                ISR(nzz)=i
              ENDIF
              IF(nzz.LE.ISCMAX) THEN
                ISC(nzz)=j
              ENDIF
            ENDIF
          ENDDO !j
        ENDDO !i
        NZTOT=nzz

C LKC 16-JAN-1999 Change to error statements for additional info output
C        CALL ASSERT(nzz.LE.ISRMAX,'>>Increase ISRMAX',ERROR,*9999)
C        CALL ASSERT(nzz.LE.ISCMAX,'>>Increase ISCMAX',ERROR,*9999)
        IF(nzz.GT.ISRMAX) THEN
          WRITE(ERROR,'(''Increase ISRMAX ( '',I8,'')'')') ISRMAX
          GOTO 9999
        ENDIF
        IF(nzz.GT.ISCMAX) THEN
          WRITE(ERROR,'(''Increase ISCMAX ( '',I8,'')'')') ISCMAX
          GOTO 9999
        ENDIF


      ELSE IF(SPARSENESS.EQ.3) THEN !Compressed column format

C LKC 16-JAN-1999 Change to error statements for additional info output
C        CALL ASSERT((N+1).LE.ISCMAX,'>>Increase ISCMAX',ERROR,*9999)
        IF(N+1.GT.ISCMAX) THEN
          WRITE(ERROR,'(''Increase ISCMAX ( '',I8,'')'')') ISCMAX
          GOTO 9999
        ENDIF

CC       Initalise sparsity arrays
CCC$DOACROSS local(i) share(N, ISC)
C        DO i=1,N+1
C          ISC(i)=1
C        ENDDO !i
CCC$DOACROSS local(nz) share(ISRMAX, ISR)
C        DO nz=1,ISRMAX
C          ISR(nz)=0
C        ENDDO !nz
C       Calculate sparsity pattern
        nzz=1
        DO j=1,N
          ISC(j)=nzz
          DO i=1,M
            IF(WORK_ARRAY(j,i)) THEN
C***          The element is set so add it to the sparsity pattern
              IF(nzz.LE.ISRMAX) THEN
C***            Set the column numbers
                ISR(nzz)=i
              ENDIF
              nzz=nzz+1
            ENDIF
          ENDDO !i
        ENDDO !j
        ISC(N+1)=nzz
        NZTOT=nzz-1
C LKC 16-JAN-1999 Change to error statements for additional info output
C        CALL ASSERT(NZTOT.LE.ISRMAX,'>>Increase ISRMAX',ERROR,*9999)
        IF(NZTOT.GT.ISRMAX) THEN
          WRITE(ERROR,'(''Increase ISRMAX ( '',I8,'')'')') ISRMAX
          GOTO 9999
        ENDIF

      ELSE IF(SPARSENESS.EQ.5) THEN
        !Umfpack row-column format
        nzz=0
        DO i=1,M
          DO j=1,N
            IF(WORK_ARRAY(j,i)) THEN
C***          The element is set so add it to the sparsity pattern
              nzz=nzz+1
              IF(nzz.LE.ISCMAX) THEN
                ISC(nzz)=i
              ENDIF
            ENDIF
          ENDDO !j
        ENDDO !i
        NZTOT=nzz
        DO i=1,M
          DO j=1,N
            IF(WORK_ARRAY(j,i)) THEN
              nzz=nzz+1
              IF(nzz.LE.ISCMAX) THEN
                ISC(nzz)=j
              ENDIF
            ENDIF
          ENDDO !j
        ENDDO !i
C LKC 16-JAN-1999 Change to error statements for additional info output
C        CALL ASSERT(nzz.LE.ISCMAX,'>>Increase ISCMAX',ERROR,*9999)
        IF(nzz.GT.ISCMAX) THEN
          WRITE(ERROR,'(''Increase ISCMAX ( '',I8,'')'')') ISCMAX
          GOTO 9999
        ENDIF

      ELSE
        ERROR='>>Unknown sparsity format'
        GOTO 9999
      ENDIF

      CALL EXITS('CALC_SPARSE')
      RETURN
 9999 CALL ERRORS('CALC_SPARSE',ERROR)
      CALL EXITS('CALC_SPARSE')
      RETURN 1
      END


