      SUBROUTINE BMG3_SymStd_COPY_rV_23_yz( 
     &                       Q, Q2, iPL, Nx, Ny, Nz 
     &                       )

C =======================================================================
C 
C   BMG3_SymStd_COPY_rV_23_yz
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_COPY_rV_23_yz copies a REAL 2D vector to the
C   iPL{th}-(x,y) plane of a REAL 3D vector.
C
C =======================================================================

      IMPLICIT   NONE

C ---------------------------
C     Argument Declarations:
C
      INTEGER  iPL, Nx, Ny, Nz
      REAL*8   Q(Nx,Ny,Nz), Q2(Ny,Nz)

C ---------------------------
C     Local Declarations:
C
      INTEGER  j, k, kk, Nym2, Nzm2, maxYZ

C ========================================================================

      Nym2 = Ny-2
      Nzm2 = Nz-2

      maxYZ = Nym2*Nzm2

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(j,k,kk)
C$OMP& SHARED(iPL,maxYZ,Nym2,Q,Q2)
      DO kk=0,maxYZ-1
         j=mod(kk,Nym2)+2
         k=(kk/Nym2)+2
         Q(iPL,j,k) = Q2(j,k)
      ENDDO
C$OMP END PARALLEL DO

C =====================

      RETURN
      END

