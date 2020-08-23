      SUBROUTINE BMG3_SymStd_COPY_rV_23_xz( 
     &                       Q, Q2, iPL, Nx, Ny, Nz 
     &                       )

C =======================================================================
C 
C   BMG3_SymStd_COPY_rV_23_xz
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_COPY_rV_23_xz copies the a REAL 2D vector to the
C   iPL{th}-(x,z) plane of a REAL 3D vector.
C
C =======================================================================

      IMPLICIT   NONE

C ---------------------------
C     Argument Declarations:
C
      INTEGER  iPL, Nx, Ny, Nz
      REAL*8   Q(Nx,Ny,Nz), Q2(Nx,Nz)

C ---------------------------
C     Local Declarations:
C
      INTEGER  i, k, kk, Nxm2, Nzm2, maxXZ

C ========================================================================

      Nxm2 = Nx-2
      Nzm2 = Nz-2

      maxXZ = Nxm2*Nzm2

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,k,kk)
C$OMP& SHARED(iPL,maxXZ,Nxm2,Q,Q2)
      DO kk=0,maxXZ-1
         i=mod(kk,Nxm2)+2
         k=(kk/Nxm2)+2
         Q(i,iPL,k) = Q2(i,k)
      ENDDO
C$OMP END PARALLEL DO

C =====================

      RETURN
      END

