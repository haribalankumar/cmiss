      SUBROUTINE BMG3_SymStd_SETUP_cg_LU( 
     &                SO, ii, jj, kk, NStncl, abd, nabd1, nabd2,
     &                BMG_IOFLAG, BMG_iPARMS
     &                )

C***BEGIN PROLOGUE  BMG3_SymStd_SETUP_cg_LU
C***DESCRIPTION
C             BMG3_SymStd_SETUP_cg_LU sets up the matrix on the 
C             coarsest grid, and using a linpack routine it forms 
C             the l-u decomposition of the matrix.
C***PARAMETERS
C***INPUT
C   SO        Refer to BMG3D.
C   II        Number of grid points in x direction, including
C             two fictitious points.
C   JJ        Number of grid points in y direction, including
C             two fictitious points.
C   KK        Number of grid points in z direction, including
C             two fictitious points.
C   NABD1     Leading dimension of ABD.
C***INPUT/OUTPUT
C   ABD       Refer to BMG3D.
C
C***ROUTINES CALLED  DPBFA
C***CHANGES
C
C  1999/12/07  - code indented, do loops end on enddo (M.Berndt)
C  2000/02/27  - removed dimension "*", cleaned up declarations (JDM)
C
C***END PROLOGUE  BMG3_SymStd_SETUP_cg_LU
C
C ==========================================================================

      IMPLICIT NONE

C ----------------------------
C     Includes
C 
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER ii, info, jj, kk, nabd1, nabd2, NStncl
      REAL*8  abd(nabd1,nabd2), so(ii,jj,kk,NStncl)
      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)
      INTEGER BMG_iPARMS(NBMG_iPARMS)

C ----------------------------
C     Local Declarations
C
      INTEGER i, i1, i2, ibw, j, j1, k, k1, kl, BC

C =========================================================================

      BC = BMG_iPARMS(id_BMG3_BC)

C -------------------------------------------------------
C     Copy the operator on the coarsest grid into ABD 
C -------------------------------------------------------

      info = 0

      i1=ii-1
      j1=jj-1
      k1=kk-1

      i2=i1-1
      ibw=i2*j1+1


      IF ( NStncl.EQ.14 ) THEN 

C$OMP PARALLEL DO
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k,kl)
C$OMP& SHARED(abd,i1,j1,k1,so)
         DO kl=1, (i1-1)*(j1-1)*(k1-1)
            i = mod((kl-1),(i1-1))+2
            j = mod((kl-1)/(i1-1),(j1-1))+2
            k = (kl-1)/((i1-1)*(j1-1))+2
            !
            abd(ibw+1,kl)    = so(i,j,k,kp)
            abd(ibw,kl)      = -so(i,j,k,kpw)
            abd(ibw-i1+3,kl) = -so(i+1,j,k,kpnw)
            abd(ibw-i1+2,kl) = -so(i,j,k,kps)
            abd(ibw-i1+1,kl) = -so(i,j,k,kpsw)
            !
            abd(ibw-(j1-2)*i2+2,kl) = -so(i+1,j+1,k,kbne)
            abd(ibw-(j1-2)*i2+1,kl) = -so(i,j+1,k,kbn)
            abd(ibw-(j1-2)*i2,kl)   = -so(i,j+1,k,kbnw)
            abd(ibw-(j1-1)*i2+2,kl) = -so(i+1,j,k,kbe)
            abd(ibw-(j1-1)*i2+1,kl) = -so(i,j,k,kb)
            abd(ibw-(j1-1)*i2,kl)   = -so(i,j,k,kbw)
            !
            abd(3,kl) = -so(i+1,j,k,kbse)
            abd(2,kl) = -so(i,j,k,kbs)
            abd(1,kl) = -so(i,j,k,kbsw)
            !
         ENDDO
C$OMP END PARALLEL DO

      ELSE IF ( NStncl.EQ.4 ) THEN 
         
C$OMP PARALLEL DO
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k,kl)
C$OMP& SHARED(abd,i1,j1,k1,so)
         DO kl=1, (i1-1)*(j1-1)*(k1-1)
            i = mod((kl-1),(i1-1))+2
            j = mod((kl-1)/(i1-1),(j1-1))+2
            k = (kl-1)/((i1-1)*(j1-1))+2
            !
            abd(ibw+1,kl)    = so(i,j,k,kp)
            abd(ibw,kl)      = -so(i,j,k,kpw)
            abd(ibw-i1+3,kl) = rZERO
            abd(ibw-i1+2,kl) = -so(i,j,k,kps)
            abd(ibw-i1+1,kl) = rZERO
            !
            abd(ibw-(j1-2)*i2+2,kl) = rZERO
            abd(ibw-(j1-2)*i2+1,kl) = rZERO
            abd(ibw-(j1-2)*i2,kl)   = rZERO
            abd(ibw-(j1-1)*i2+2,kl) = rZERO
            abd(ibw-(j1-1)*i2+1,kl) = -so(i,j,k,kb)
            abd(ibw-(j1-1)*i2,kl)   = rZERO
            !
            abd(3,kl) = rZERO
            abd(2,kl) = rZERO
            abd(1,kl) = rZERO
            !
         ENDDO
C$OMP END PARALLEL DO

      ELSE
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,500)
            WRITE(*,510) 'NEED: NStncl = 4 or 14 '
            WRITE(*,510) 'HAVE: NStncl = ', NStncl
         END IF
         
         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,4)
         RETURN
         
      ENDIF

      IF( BC.EQ.BMG_BCs_indef_nonper ) THEN
         ABD(ibw+1,1) = 2.D0*ABD(ibw+1,1)
      ENDIF         

C -------------------------------------------------------
C     Factor using the LAPACK routine DPBTRF
C -------------------------------------------------------

      CALL DPBTRF('U', NABD2, IBW, ABD, NABD1, INFO) 

      IF (INFO .NE. 0) THEN

         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'Coarse grid Cholesky decomposition failed!'
            WRITE(*,510) 'INFO = ', INFO
         END IF

         BMG_iPARMS(id_BMG3_Ext_Err_Code) = INFO
         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,5)
         RETURN

      ENDIF


C ========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_SETUP_cg_LU.f',/,5X,A)
 510  FORMAT (5X,A,1X,I3)

C ===========================================

      RETURN
      END
