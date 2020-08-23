      SUBROUTINE EIGENPROBLEM(AMMAX,ANMAX,BMMAX,BNMAX,EVECMMAX,
     '  M,MAXNUMEIGEN,NUMEIGEN,NUMEIGENREQUIRED,OUTPUTCODE,
     '  PROBLEM_TYPE,SOLVER_TYPE,A,B,EIGENVALUES,EIGENVECTORS,
     '  ERROR,*)

C#### Subroutine: EIGENPROBLEM
C###  Description:
C###    EIGENPROBLEM solves either an eigenproblem Ax=lambda.x or a
C###    generalised eigenproblem Ax=lambda.Bx. The problem type is
C###    determined by the PROBLEM_TYPE with PROBLEM_TYPE=1 meaning
C###    a standard eigenproblem and PROBLEM_TYPE=2 meaning a generalised
C###    eigenproblem. How the problem is solved is determined by
C###    SOLVER_TYPE.

C**** (re)created by Martin Buist 29 April 1997

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER AMMAX,ANMAX,BMMAX,BNMAX,EVECMMAX,M,
     '  MAXNUMEIGEN,NUMEIGEN,NUMEIGENREQUIRED,OUTPUTCODE,
     '  PROBLEM_TYPE,SOLVER_TYPE
      REAL*8 A(AMMAX,ANMAX),B(BMMAX,BNMAX),EIGENVALUES(MAXNUMEIGEN),
     '  EIGENVECTORS(EVECMMAX,MAXNUMEIGEN)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER NB
      PARAMETER (NB=32) ! Optimal blocksize for DSYTRD from ILAENV.
      INTEGER i,INFO,no,nt,ntt,SMALLEST_LOC,j
      REAL ELAPSED_TIME,TIME_START1(1),TIME_START2(1),TIME_STOP(1)
      REAL*8 DNRM2,SMALLEST,TEMP,VEC_NORM,
     '  SUM,W(M),WORKD((NB+2)*M)

      CALL ENTERS('EIGENPROBLEM',*9999)

      CALL CPU_TIMER(CPU_USER,TIME_START1)

      IF(PROBLEM_TYPE.EQ.1) THEN !standard eigenproblem
        ERROR='>>Not implemented yet'
        GOTO 9999
      ELSE IF(PROBLEM_TYPE.EQ.2) THEN !generalised eigenproblem
        IF(SOLVER_TYPE.EQ.1) THEN
          INFO=0
          CALL DSYGV(1,'V','L',M,A,AMMAX,B,BMMAX,W,WORKD,(NB+2)*M,INFO)
          IF(INFO.NE.0) THEN
            WRITE(ERROR,'('' >>Non-zero return from DSYGV'')')
            GOTO 9999
          ENDIF
          NUMEIGEN=NUMEIGENREQUIRED
          DO j=1,NUMEIGEN
            EIGENVALUES(j)=W(j)
          ENDDO
          DO i=1,M
            DO j=1,NUMEIGEN
              EIGENVECTORS(i,j)=A(i,j)
            ENDDO
          ENDDO
        ELSE
          ERROR='>>Invalid SOLVER_TYPE'
          GOTO 9999
        ENDIF
      ELSE
        ERROR='>>Invalid PROBLEM_TYPE'
        GOTO 9999
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_START2)
      ! This seems to sort so that the eigenvalues are in ascending order.
      ! Would the order of ascending magnitude be better?
      ! (For positive-definite matrices they are all positive anyway.)
      DO nt=1,NUMEIGEN-1
        SMALLEST=EIGENVALUES(nt)
        SMALLEST_LOC=nt
        DO ntt=nt,NUMEIGEN
          IF(EIGENVALUES(ntt).LT.SMALLEST) THEN
            SMALLEST=EIGENVALUES(ntt)
            SMALLEST_LOC=ntt
          ENDIF
        ENDDO !ntt
        IF(SMALLEST_LOC.NE.nt) THEN
          TEMP=EIGENVALUES(nt)
          EIGENVALUES(nt)=EIGENVALUES(SMALLEST_LOC)
          EIGENVALUES(SMALLEST_LOC)=TEMP
          DO no=1,M
            TEMP=EIGENVECTORS(no,nt)
            EIGENVECTORS(no,nt)=EIGENVECTORS(no,SMALLEST_LOC)
            EIGENVECTORS(no,SMALLEST_LOC)=TEMP
          ENDDO !no
        ENDIF
      ENDDO !nt
      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
      IF(OUTPUTCODE.GE.1) THEN
        WRITE(OP_STRING,'(/'' CPU time to sort eigenvalues : '','
     '    //'D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_START2)
      DO nt=1,NUMEIGEN
        VEC_NORM=DNRM2(M,EIGENVECTORS(1,nt),1)
        IF(VEC_NORM.GE.1.0d-10) THEN
          ! Eigenvectors are normalized, but, even if all eigenvalues have
          ! multiplicity 1, eigenvectors are not unique in their sign.
          ! Try to make them as positive as possible so as to reduce the
          ! variation in results.
          SUM=0.d0
          DO i=1,M
            SUM=SUM+EIGENVECTORS(i,nt)
          ENDDO ! i
          IF(SUM.LT.0.d0) THEN
            VEC_NORM=-VEC_NORM
          ENDIF
          CALL DSCAL(M,1.0d0/VEC_NORM,EIGENVECTORS(1,nt),1)
        ELSE
          WRITE(ERROR,'('' >>Eigenvector '',I5,'' is zero'')') nt
          GOTO 9999
        ENDIF
      ENDDO !nt
      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START2(1)
      IF(OUTPUTCODE.GE.1) THEN
        WRITE(OP_STRING,'(/'' CPU time to normalise eigenvectors '
     '    //': '',D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(OUTPUTCODE.GE.2) THEN
        WRITE(OP_STRING,'(/'' Eigenvalues :'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(5(1X,D13.5),:/(5(1X,D13.5)))')
     '    (EIGENVALUES(nt),nt=1,NUMEIGEN)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/'' Eigenvectors :'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        DO nt=1,NUMEIGEN
          WRITE(OP_STRING,
     '      '(1X,I5,'')'',5(1X,D13.5),:/(7X,5(1X,D13.5)))')
     '      nt,(EIGENVECTORS(no,nt),no=1,M)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDDO !nt
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)
      IF(OUTPUTCODE.GT.1) THEN
        WRITE(OP_STRING,'(/'' Total CPU time for eigenproblem '
     '    //'solution : '',D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF


      CALL EXITS('EIGENPROBLEM')
      RETURN
 9999 CALL ERRORS('EIGENPROBLEM',ERROR)
      CALL EXITS('EIGENPROBLEM')
      RETURN 1
      END

