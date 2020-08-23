      SUBROUTINE EIGEN1(n,IA,A,EVAL,EVEC,ERROR,*)

C#### Subroutine: EIGEN1
C###  Description:
C###    EIGEN1 finds N eigenvalues and eigenvectors for REAL*8 symmetric
C###    matrix A(IA,IA) (IA<5) using LAPACK routine DSYEV.
C###    The eigenvalues are returned in EVAL(j),j=1,N.

CC###    matrix A(IA,IA) (IA<5) using Nag routine F02ABF.  The eigenvalues
CC###    are returned in EVAL(j),j=1,N.
      IMPLICIT NONE
!     Parameter List
      INTEGER IA,n
      REAL*8 A(IA,IA),EVAL(IA),EVEC(IA,IA)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IFAIL,ni,mi
      REAL*8 WORK(15)

      CALL ENTERS('EIGEN1',*9999)
      CALL ASSERT(IA.LE.5,'EIGEN1 is only dimensioned up to 5',
     '  ERROR,*9999)
      IFAIL=1
C      CALL F02ABF(A,IA,n,EVAL,EVEC,IA,WORK,IFAIL)
      DO ni=1,n
        DO mi=1,n
          EVEC(ni,mi)=A(ni,mi)
        ENDDO
      ENDDO
C MLB 19/3/97
C This may not give evectors as accurately as NAG
        CALL DSYEV('V','L',n,EVEC,IA,EVAL,WORK,15,IFAIL)
      IF(ifail.NE.0) THEN
        ERROR='IFAIL=1 in DSYEV'
        GO TO 9999
      ENDIF

      CALL EXITS('EIGEN1')
      RETURN
 9999 CALL ERRORS('EIGEN1',ERROR)
      CALL EXITS('EIGEN1')
      RETURN 1
      END


