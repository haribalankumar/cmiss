      SUBROUTINE OPSTFMAT(COLLIST,ISC,ISR,IUNIT,M,N,NZMATRIX,ROWLIST,
     '  SPARSENESS,MATRIX,VECTOR,MATNAME,VECTNAME,NOLIST,OPMATRIX,
     '  OPVECTOR,ERROR,*)

C#### Subroutine: OPSTFMAT
C###  Description:
C###    OPSTFMAT outputs global stiffness matrices and/or global
C###    stiffness vectors to a screen unit in a standard way.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER COLLIST(0:*),ISC(*),ISR(*),IUNIT,M,N,NZMATRIX,
     '  ROWLIST(0:*),SPARSENESS
      REAL*8 MATRIX(*),VECTOR(NYROWM)
      CHARACTER ERROR*(*),MATNAME*3,VECTNAME*3
      LOGICAL NOLIST,OPMATRIX,OPVECTOR
!     Local Variables
      INTEGER i,j,noclist,norlist
      CHARACTER COL*5,NAME*3,ROW*5

      CALL ENTERS('OPSTFMAT',*9999)

      IF(SPARSENESS.EQ.0) THEN !No sparsity
        IF(NOLIST) THEN
          DO i=1,M
            IF(OPVECTOR) THEN
              WRITE(OP_STRING,'(1X,A3,''('',I5,'')    ='',1X,D11.4)')
     '          VECTNAME,i,VECTOR(i)
              CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
            ENDIF
            IF(OPMATRIX) THEN
C cpb 27/6/96 Adding long output for matrix
C              WRITE(OP_STRING,'(1X,A3,''('',I5,'',1..)='',5(1X,D11.4),'
C     '          //'/:(16X,5(1X,D11.4)))') MATNAME,i,(MATRIX(i+(j-1)*N),
C     '          j=1,N)
C              CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
              WRITE(NAME,'(A3)') MATNAME
              WRITE(ROW,'(I5)') i
              CALL WRITE_LONG(DPTYPE,i,M,IUNIT,i+(N-1)*M,5,5,%VAL(0),
     '          MATRIX,
     '          '('' '//NAME//'('//ROW//',1..)='',5(1X,D11.4))',
     '          '(16X,5(1X,D11.4))',ERROR,*9999)
            ENDIF
          ENDDO !i
        ELSE
          DO norlist=1,ROWLIST(0)
            i=ROWLIST(norlist)
            IF(OPVECTOR) THEN
              WRITE(OP_STRING,'(1X,A3,''('',I5,'')    ='',1X,D11.4)')
     '          VECTNAME,i,VECTOR(i)
              CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
            ENDIF
            IF(OPMATRIX) THEN
              WRITE(NAME,'(A3)') MATNAME
              WRITE(ROW,'(I5)') i
              CALL WRITE_LONG_IDX(DPTYPE,N,COLLIST,IUNIT,5,5,
     '          %VAL(0),MATRIX(i),
     '          '('' '//NAME//'('//ROW//',1..)='',5(1X,D11.4))',
     '          '(16X,5(1X,D11.4))',ERROR,*9999)
            ENDIF
          ENDDO !i
        ENDIF
      ELSE IF(SPARSENESS.EQ.1) THEN !Compressed row sparsity
        IF(NOLIST) THEN
          DO i=1,M
            IF(OPVECTOR) THEN
              WRITE(OP_STRING,'(1X,A3,''('',I5,'')    ='',1X,D11.4)')
     '          VECTNAME,i,VECTOR(i)
              CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
            ENDIF
            IF(OPMATRIX) THEN
              WRITE(NAME,'(A3)') MATNAME
              WRITE(ROW,'(I5)') i
              CALL WRITE_LONG(DPTYPE,ISR(i),1,IUNIT,ISR(i+1)-1,
     '          5,5,%VAL(0),MATRIX,'('' '//NAME//'('//ROW//',1..)='','
     '          //'5(1X,D11.4))','(16X,5(1X,D11.4))',ERROR,*9999)
            ENDIF
          ENDDO !i
        ELSE
          DO norlist=1,ROWLIST(0)
            i=ROWLIST(norlist)
            IF(OPVECTOR) THEN
              WRITE(OP_STRING,'(1X,A3,''('',I5,'')    ='',1X,D11.4)')
     '          VECTNAME,i,VECTOR(i)
              CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
            ENDIF
            IF(OPMATRIX) THEN
              WRITE(NAME,'(A3)') MATNAME
              WRITE(ROW,'(I5)') i
              CALL WRITE_LONG(DPTYPE,ISR(i),1,IUNIT,ISR(i+1)-1,
     '          5,5,%VAL(0),MATRIX,'('' '//NAME//'('//ROW//',1..)='','
     '          //'5(1X,D11.4))','(16X,5(1X,D11.4))',ERROR,*9999)
            ENDIF
          ENDDO !i
        ENDIF
      ELSE IF(SPARSENESS.EQ.2.OR.SPARSENESS.EQ.4) THEN !Row-column sparsity
        IF(OPVECTOR) THEN
          IF(NOLIST) THEN
            WRITE(NAME,'(A3)') VECTNAME
            CALL WRITE_LONG(DPTYPE,1,1,IUNIT,M,6,6,%VAL(0),VECTOR,
     '        '('' '//NAME//':'',6(1X,D11.4))',
     '        '(5X,6(1X,D11.4))',ERROR,*9999)
          ELSE
            WRITE(OP_STRING,
     '            '(1X,A3,'':'',6(1X,D11.4),:/(5X,6(1X,D11.4)))')
     '        VECTNAME,(VECTOR(ROWLIST(norlist)),norlist=1,ROWLIST(0))
            CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        IF(OPMATRIX) THEN
          WRITE(NAME,'(A3)') MATNAME
          CALL WRITE_LONG(DPTYPE,1,1,IUNIT,NZMATRIX,6,6,
     '        %VAL(0),MATRIX,
     '        '('' '//NAME//':'',6(1X,D11.4))',
     '        '(5X,6(1X,D11.4))',ERROR,*9999)
        ENDIF
      ELSE IF(SPARSENESS.EQ.3) THEN !Compressed column sparsity
        IF(NOLIST) THEN
          DO i=1,M
            IF(OPVECTOR) THEN
              WRITE(OP_STRING,'(1X,A3,''('',I5,'')    ='',1X,D11.4)')
     '          VECTNAME,i,VECTOR(i)
              CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
          IF(OPMATRIX) THEN
            DO j=1,N
              WRITE(NAME,'(A3)') MATNAME
              WRITE(COL,'(I5)') j
              CALL WRITE_LONG(DPTYPE,ISC(j),1,IUNIT,ISC(j+1)-1,
     '          5,5,%VAL(0),MATRIX,'('' '//NAME//'(1..,'//COL//')='','
     '          //'5(1X,D11.4))','(16X,5(1X,D11.4))',ERROR,*9999)
            ENDDO !i
          ENDIF
        ELSE
          DO norlist=1,ROWLIST(0)
            i=ROWLIST(norlist)
            IF(OPVECTOR) THEN
              WRITE(OP_STRING,'(1X,A3,''('',I5,'')    ='',1X,D11.4)')
     '          VECTNAME,i,VECTOR(i)
              CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO !norlist
          IF(OPMATRIX) THEN
            DO noclist=1,ROWLIST(0)
              j=ROWLIST(noclist)
              WRITE(NAME,'(A3)') MATNAME
              WRITE(COL,'(I5)') j
              CALL WRITE_LONG(DPTYPE,ISC(j),1,IUNIT,ISC(j+1)-1,
     '          5,5,%VAL(0),MATRIX,'('' '//NAME//'(1..,'//COL//')='','
     '          //'5(1X,D11.4))','(16X,5(1X,D11.4))',ERROR,*9999)
            ENDDO !j
          ENDIF
        ENDIF
      ELSE IF(SPARSENESS.EQ.5) THEN !Umfpack row-column sparsity
        IF(OPVECTOR) THEN
          IF(NOLIST) THEN
            WRITE(NAME,'(A3)') VECTNAME
            CALL WRITE_LONG(DPTYPE,1,1,IUNIT,M,6,6,%VAL(0),VECTOR,
     '        '('' '//NAME//':'',6(1X,D11.4))',
     '        '(5X,6(1X,D11.4))',ERROR,*9999)
          ELSE
            WRITE(OP_STRING,
     '        '(1X,A3,'':'',6(1X,D11.4),:/(5X,6(1X,D11.4)))')
     '        VECTNAME,(VECTOR(ROWLIST(norlist)),norlist=1,ROWLIST(0))
            CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        IF(OPMATRIX) THEN
          WRITE(NAME,'(A3)') MATNAME
          CALL WRITE_LONG(DPTYPE,1,1,IUNIT,NZMATRIX,6,6,
     '        %VAL(0),MATRIX,
     '        '('' '//NAME//':'',6(1X,D11.4))',
     '        '(5X,6(1X,D11.4))',ERROR,*9999)
        ENDIF
      ELSE
        ERROR='>>Unkown sparsity'
        GOTO 9999
      ENDIF

      CALL EXITS('OPSTFMAT')
      RETURN
 9999 CALL ERRORS('OPSTFMAT',ERROR)
      CALL EXITS('OPSTFMAT')
      RETURN 1
      END


