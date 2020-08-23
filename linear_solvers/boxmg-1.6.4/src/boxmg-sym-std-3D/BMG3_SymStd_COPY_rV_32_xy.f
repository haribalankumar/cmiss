      SUBROUTINE BMG3_SymStd_COPY_rV_32_xy( 
     &                       Q, Q2, iPL, Nx, Ny, Nz 
     &                       )

C =======================================================================
C 
C   BMG3_SymStd_COPY_rV_32_xy
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_COPY_rV_32_xy copies the a REAL 2D vector to the
C   iPL{th}-(x,y) plane of a REAL 3D vector.
C
C =======================================================================

      IMPLICIT   NONE

C ---------------------------
C     Argument Declarations:
C
      INTEGER  iPL, Nx, Ny, Nz
      REAL*8   Q(Nx,Ny,Nz), Q2(Nx,Ny)

C ---------------------------
C     Local Declarations:
C
      INTEGER  i, j, kk, Nxm2, Nym2, maxXY

C ========================================================================

      Nxm2 = Nx-2
      Nym2 = Ny-2
      
      maxXY = Nxm2*Nym2

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,kk)
C$OMP& SHARED(iPL,maxXY,Nxm2,Q,Q2)
      DO kk=0,maxXY-1
         i=mod(kk,Nxm2)+2
         j=(kk/Nxm2)+2
         Q2(i,j) = Q(i,j,iPL)
      ENDDO
C$OMP END PARALLEL DO

C =====================

      RETURN
      END

