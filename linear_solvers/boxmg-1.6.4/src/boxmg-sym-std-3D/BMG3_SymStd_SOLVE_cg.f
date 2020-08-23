      SUBROUTINE BMG3_SymStd_SOLVE_cg( 
     &         q, qf, ii, jj, kk, abd, bbd, nabd1, nabd2,
     &         BMG_IOFLAG, BMG_iPARMS 
     &         )

C***BEGIN PROLOGUE  BMG3_SymStd_SOLVE_cg
C***DESCRIPTION
C             BMG3_SymStd_SOLVE_cg does a direct solve on 
C             the coarsest grid. it uses the linpack routine spbsl.
C***PARAMETERS
C***INPUT
C   QF        Refer to BMG3D.
 
C   II        Number of grid points in x direction, including
C             two fictitious points.
C   JJ        Number of grid points in y direction, including
C             two fictitious points.
C   KK        Number of grid points in z direction, including
C             two fictitious points.
C   ABD       Refer to BMG3D.
C   BBD       Refer to BMG3D.
C   NABD1     Refer to BMG3D.
C   ABD       Refer to BMG3D.
C   NABD1     Refer to BMG3D.
C***OUTPUT
C   Q         Refer to BMG3D
C
C***ROUTINES CALLED  DPBSL
C***CHANGES
C
C  1999/12/06  - indented the code, made do-loops end on enddo (M.Berndt)
C  2000/02/27  - removed dimension "*", cleaned up declarations (JDM)
C
C***END PROLOGUE  BMG3_SymStd_SOLVE_cg
C
C ==========================================================================

      IMPLICIT NONE
      
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER ii, jj, kk, nabd1, nabd2
      REAL*8  abd(nabd1,nabd2), bbd(nabd2), q(ii,jj,kk), qf(ii,jj,kk)

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)
      INTEGER  BMG_iPARMS(NBMG_iPARMS)

C ----------------------------
C     Local Declarations
C
      integer i, i1, i2, ibw, j, j1, k, kt, k1, n, info

C =========================================================================

c     
c     direct solve on coarsest grid
c     

      i1=ii-1
      j1=jj-1
      k1=kk-1

      i2=i1-1
      n=i2*(j1-1)*(kk-2)
      ibw=i2*j1+1

      kt=0
      do k=2,k1
         do j=2,j1
            do i=2,i1
               kt=kt+1
               bbd(kt)=qf(i,j,k)
            enddo
         enddo
      enddo

C -------------------------------------------------------
C     Solve using the LAPACK routine DPBTRS
C -------------------------------------------------------

      CALL DPBTRS ('U', KT, IBW, 1, ABD, NABD1, BBD, NABD2, INFO) 

      IF (INFO .NE. 0) THEN

         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'Coarse grid solve failed!'
            WRITE(*,510) 'INFO = ', INFO
         END IF

         BMG_iPARMS(id_BMG3_Ext_Err_Code) = INFO
         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,20)
         RETURN

      ENDIF


      kt=0
      do k=2,k1
         do j=2,j1
            do i=2,i1
               kt=kt+1
               q(i,j,k)=bbd(kt)
            enddo
         enddo
      enddo



C ========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_SOLVE_cg.f',/,5X,A)
 510  FORMAT (5X,A,1X,I3)

C ===========================================

      RETURN
      END
