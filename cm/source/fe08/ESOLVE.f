      SUBROUTINE ESOLVE(VECT,N,A,LDA,EVALUES,EVECTORS,ERROR,*)

C#### Subroutine: ESOLVE
C###  Description:
C###    ESOLVE is designed to buffer the calls to the lapack solvers
C###    to make it easier change solvers. This routine finds the
C###    eigenvalues and/or eigenvectors of a symmetric matrix, eg. a
C###    symmetric stress tensor. The current solver is DSYEV.
C###    IF VECT is true, eigenvectors are calculated.
C**** Created by Martin Buist, 30 April 1997

      IMPLICIT NONE

      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'mach00.inc'

!     Parameter List
      INTEGER LDA,N
      REAL*8 A(LDA,LDA),EVALUES(LDA),EVECTORS(LDA,LDA)
      CHARACTER ERROR*(*)
      LOGICAL VECT
!     Local Variables
      INTEGER i,j,IFAIL
C KAT 14Jan00:
      INTEGER*4 WORKSPACE_PTR
C      REAL*8 WORKSPACE(LDA*3)

      CALL ENTERS('ESOLVE',*9999)

      WORKSPACE_PTR=0

      DO i=1,N
        DO j=1,N
          EVECTORS(i,j)=A(i,j)
        ENDDO
      ENDDO

      CALL ALLOCATE_MEMORY(LDA*3,1,DPTYPE,WORKSPACE_PTR,MEM_INIT,
     '  ERROR,*9999)

      IFAIL=0
      IF(VECT) THEN
        CALL DSYEV('V','L',N,EVECTORS,LDA,EVALUES,%VAL(WORKSPACE_PTR),
     '    LDA*3,IFAIL)
      ELSE
        CALL DSYEV('N','L',N,EVECTORS,LDA,EVALUES,%VAL(WORKSPACE_PTR),
     '    LDA*3,IFAIL)
      ENDIF

      CALL FREE_MEMORY(WORKSPACE_PTR,ERROR,*9999)

      IF(IFAIL.NE.0) THEN
        WRITE(OP_STRING,'(''IFAIL non-zero in return from DSYEV'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        GOTO 9999
      ENDIF

      CALL EXITS('ESOLVE')
      RETURN
 9999 CALL ERRORS('ESOLVE',ERROR)
      CALL EXITS('ESOLVE')
      RETURN 1
      END


