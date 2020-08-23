      SUBROUTINE CALC_LATTICE_WEIGHTS(nq,SUPPORT,M,XQ,ERROR,*)

C#### Subroutine: CALC_LATTICE_WEIGHTS
C###  Description:
C###    CALC_LATTICE_WEIGHTS calls the Moore-Penrose 
C###    pseudo inverse subroutine for a matrix A which is made
C###    up of the x offsets of supporting grid points
C###    around the grid point of interest.
C###    The inverse of this matrix A then gives the weights
C###    to be used in finite difference methods.      
C###    These nine weights are in order du/dx1, du/dx2, du/dx3,
C###    d^2u/dx1^2, d^2u/dx2^2, d^2u/dx3^2, d^2u/dx1dx2,
C###    d^2u/dx1dx3, d^2u/dx2dx3.  

      IMPLICIT NONE

      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER nq,SUPPORT(0:NQGM)
      REAL*8 M(9,NQGM),XQ(NJM,NQM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,nq1,nqi
      REAL*8 A(NQGM,9)

      CALL ENTERS('CALC_LATTICE_WEIGHTS',*9999)

      DO i=1,NQGM
        DO j=1,9
          A(i,j)=0.0d0
          M(j,i)=0.0d0
        ENDDO
      ENDDO

C     Construct the A matrix for nq
      DO i=1,SUPPORT(0)
        nq1=SUPPORT(i)
        IF(nq1.EQ.nq) nqi=i

        A(i,1)=XQ(1,nq1)-XQ(1,nq)
        IF(NJT.GE.2) A(i,2)=XQ(2,nq1)-XQ(2,nq)
        IF(NJT.EQ.3) A(i,3)=XQ(3,nq1)-XQ(3,nq)
        A(i,4)=A(i,1)*A(i,1)*0.5d0
        A(i,5)=A(i,2)*A(i,2)*0.5d0
        A(i,6)=A(i,3)*A(i,3)*0.5d0
        A(i,7)=A(i,1)*A(i,2)
        A(i,8)=A(i,1)*A(i,3)
        A(i,9)=A(i,2)*A(i,3)
      ENDDO

      CALL MP_PSEUDO_INVERSE(A,M,NQGM,9,SUPPORT(0),9,ERROR,*9999)

C     Calculate the weights for grid point nq.
      DO j=1,9
        M(j,nqi)=0.0d0
        DO i=1,SUPPORT(0)
          IF(i.NE.nqi) M(j,nqi)=M(j,nqi)-M(j,i)
        ENDDO
      ENDDO

      CALL EXITS('CALC_LATTICE_WEIGHTS')
      RETURN
 9999 CALL ERRORS('CALC_LATTICE_WEIGHTS',ERROR)
      CALL EXITS('CALC_LATTICE_WEIGHTS')
      RETURN 1
      END
