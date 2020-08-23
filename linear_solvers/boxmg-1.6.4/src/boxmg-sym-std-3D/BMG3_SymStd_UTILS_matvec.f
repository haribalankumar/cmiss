      SUBROUTINE BMG3_SymStd_UTILS_matvec(
     &                NOG, Q, QF, SO, Nx, Ny, Nz, NStncl, BMG_IJK_MAP
     &                )

C ==========================================================================
C
C   BMG3_SymStd_UTILS_matvec.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_UTILS_matvec performs a matrix multiplication that
C   is used within the pcg routine called BMG3_SymStd_SOLVE_pcg.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   2003/07/29  - written from residual calculation. (TMA)
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
      INTEGER IFD, kg, Nx, Ny, Nz, NOG, NStncl, BMG_IJK_MAP(*)
      REAL*8  Q(Nx,Ny,Nz), QF(Nx,Ny,Nz), SO(Nx,Ny,Nz,NStncl)

C ----------------------------
C     Local Declarations
C
      INTEGER i, j, k, kk, Nxm2, Nym2, Nzm2, ibmg_start, maxXYZ

C ==========================================================================

      Nxm2 = Nx-2
      Nym2 = Ny-2
      Nzm2 = Nz-2

      maxXYZ = Nxm2*Nym2*Nzm2
 
      ibmg_start = BMG_IJK_MAP(NOG)

      IF( NStncl.EQ.14 ) THEN

C$OMP PARALLEL DO 
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k,kk)
C$OMP& SHARED(BMG_IJK_MAP,ibmg_start,maxXYZ,Nxm2,Nym2,Nzm2,Q,QF,SO)
         DO kk = 0, maxXYZ-1

            i = BMG_IJK_MAP(ibmg_start+3*kk  )
            j = BMG_IJK_MAP(ibmg_start+3*kk+1)
            k = BMG_IJK_MAP(ibmg_start+3*kk+2)

            QF(i,j,k) =  SO(i,j,k,kp)*Q(i,j,k)
     &           - SO(i,j,k,kpw)*Q(i-1,j,k)
     &           - SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &           - SO(i,j+1,k,kps)*Q(i,j+1,k)
     &           - SO(i+1,j+1,k,kpsw) *Q(i+1,j+1,k)
     &           - SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &           - SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &           - SO(i,j,k,kps)*Q(i,j-1,k)
     &           - SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &           - SO(i,j,k,kb)*Q(i,j,k-1)
     &           - SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &           - SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &           - SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &           - SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &           - SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &           - SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &           - SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &           - SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &           - SO(i,j,k+1,kb)*Q(i,j,k+1)
     &           - SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &           - SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &           - SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &           - SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &           - SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &           - SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &           - SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &           - SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
            
         ENDDO
C$OMP END PARALLEL DO

      ELSE

C$OMP PARALLEL DO 
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k,kk)
C$OMP& SHARED(BMG_IJK_MAP,ibmg_start,maxXYZ,Nxm2,Nym2,Nzm2,Q,QF,SO)
         DO kk = 0, maxXYZ-1

            i = BMG_IJK_MAP(ibmg_start+3*kk  )
            j = BMG_IJK_MAP(ibmg_start+3*kk+1)
            k = BMG_IJK_MAP(ibmg_start+3*kk+2)

            QF(i,j,k) =  SO(i,j,k,kp)*Q(i,j,k)
     &           - SO(i,j,k,kpw)*Q(i-1,j,k)
     &           - SO(i,j+1,k,kps)*Q(i,j+1,k)
     &           - SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &           - SO(i,j,k,kps)*Q(i,j-1,k)
     &           - SO(i,j,k,kb)*Q(i,j,k-1)
     &           - SO(i,j,k+1,kb)*Q(i,j,k+1)

         ENDDO
C$OMP END PARALLEL DO

      ENDIF

C ==========================================================================

      RETURN
      END
