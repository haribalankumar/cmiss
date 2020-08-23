      SUBROUTINE BMG3_SymStd_residual(
     &                kg, NOG, ifd, q, qf, so, RES, ii, jj, kk, NStncl,
     &                BMG_IJK_MAP
     &                )

C ==========================================================================
C
C   BMG3_SymStd_updown.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_updown performs all the necessary tasks associated
C   with moving "down" to a coarser grid or moving "up" to 
C   a finer grid.  It is necessary because f77 does not support
C   recursion.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   2000/02/17  - written (JDM)
C   2000/02/23  - pulled out relaxation into a separate 
C                 driver routine  (M.Berndt)
C   2000/03/05  - eliminated "magic number" indexing of IGRD
C               - eliminated ijst, ijstc pointer confusion
C
C ==================================================================
C   INPUT:
C ========================
C
C
C
C ==================================================================
C   OUTPUT:
C ===========================
C
C
C
C ==================================================================
C   LOCAL:
C ========================
C
C
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER ifd, ii, jj, kg, kk, NOG, NStncl, 
     &        BMG_IJK_MAP(3*(ii-2)*(jj-2)*(kk-2))
      REAL*8  q(ii,jj,kk), qf(ii,jj,kk), RES(ii,jj,kk),
     &        SO(ii,jj,kk,NStncl)

C ----------------------------
C     Local Declarations
C
      integer I, I1, I2, J, J1, J2, K, K1, K2, Kl
      INTEGER maxIJK, maxIJ, ibmg_start
      REAL*8  starttime, endtime, OMP_GET_WTIME
      EXTERNAL OMP_GET_WTIME

C ==========================================================================

      I2=II-2
      J2=JJ-2
      K2=KK-2

      I1 = II-1
      J1 = JJ-1
      K1 = KK-1

      maxIJ  = I2*J2
      maxIJK = I2*J2*K2
 
      starttime = OMP_GET_WTIME()

      IF( kg.lt.NOG .or. IFD.ne.1 ) THEN

C$OMP PARALLEL DO 
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,J,K)
C$OMP& SHARED(i1,j1,k1,RES,so,qf,q)
         DO k=2,k1
            DO j=2,j1
               DO i=2,i1

                  RES(i,j,k) = qf(i,j,k)
     &                 + so( i ,j,k,kpw)*q(i-1,j,k) ! *
     &                 + so(i+1,j,k,kpw)*q(i+1,j,k)           
     &                 + so(i, j ,k,kps)*q(i,j-1,k) ! *
     &                 + so(i,j+1,k,kps)*q(i,j+1,k)
     &                 + so( i , j ,k,kpsw)*q(i-1,j-1,k) ! *
     &                 + so(i+1,j+1,k,kpsw)*q(i+1,j+1,k)
     &                 + so(i,j, k ,kb)*q(i,j,k-1) ! *
     &                 + so(i,j,k+1,kb)*q(i,j,k+1)
     &                 + so(i, j , k ,kbs)*q(i,j-1,k-1) ! *
     &                 + so(i,j+1,k+1,kbs)*q(i,j+1,k+1)
     &                 + so( i , j , k ,kbsw)*q(i-1,j-1,k-1) ! *
     &                 + so(i+1,j+1,k+1,kbsw)*q(i+1,j+1,k+1)
     &                 + so( i ,j+1,k,kpnw)*q(i-1,j+1,k)
     &                 + so(i+1, j ,k,kpnw)*q(i+1,j-1,k)
     &                 + so( i ,j+1, k ,kbnw)*q(i-1,j+1,k-1)
     &                 + so(i+1, j ,k+1,kbnw)*q(i+1,j-1,k+1)
     &                 + so(i,j+1, k ,kbn)*q(i,j+1,k-1)
     &                 + so(i, j ,k+1,kbn)*q(i,j-1,k+1)
     &                 + so(i+1,j+1, k ,kbne)*q(i+1,j+1,k-1)
     &                 + so( i , j ,k+1,kbne)*q(i-1,j-1,k+1)
     &                 + so(i+1,j, k ,kbe)*q(i+1,j,k-1)
     &                 + so( i ,j,k+1,kbe)*q(i-1,j,k+1)
     &                 + so(i+1, j , k ,kbse)*q(i+1,j-1,k-1)
     &                 + so( i ,j+1,k+1,kbse)*q(i-1,j+1,k+1)
     &                 + so(i+1,j,k+1,kbw )*q(i+1,j,k+1)
     &                 + so( i ,j, k ,kbw )*q(i-1,j,k-1)    
     &                 - so(i,j,k,kp)*q(i,j,k)

               ENDDO
            ENDDO
         ENDDO
C$OMP END PARALLEL DO

      ELSE

C$OMP PARALLEL DO 
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,J,K)
C$OMP& SHARED(i1,j1,k1,RES,so,qf,q)
         DO k=2,k1
            DO j=2,j1
               DO i=2,i1

                  RES(i,j,k) = qf(i,j,k)
     &                 + so(i,j,k,kpw)*q(i-1,j,k)
     &                 + so(i,j+1,k,kps)*q(i,j+1,k)
     &                 + so(i+1,j,k,kpw)*q(i+1,j,k)
     &                 + so(i,j,k,kps)*q(i,j-1,k)
     &                 + so(i,j,k,kb)*q(i,j,k-1)
     &                 + so(i,j,k+1,kb)*q(i,j,k+1)
     &                 - so(i,j,k,kp)*q(i,j,k)
                  
               ENDDO
            ENDDO
         ENDDO
C$OMP END PARALLEL DO

      ENDIF

      endtime = OMP_GET_WTIME()

      PRINT *,kg,NOG,' Computing Residual Time = ', endtime-starttime

C ==========================================================================

      RETURN
      END
