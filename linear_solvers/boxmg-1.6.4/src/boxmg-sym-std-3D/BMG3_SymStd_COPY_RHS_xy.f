      SUBROUTINE BMG3_SymStd_COPY_RHS_xy( 
     &                       SO, Q, QF, QF2, iPL, Nx, Ny, Nz, NStncl
     &                       )

C =======================================================================
C 
C   BMG3_SymStd_COPY_RHS_xy
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_COPY_RHS_xy creates the right hand side vector, QF,
C   from the iPL{th}-(x,y) plane of the current 3D vector.
C
C =======================================================================

      IMPLICIT   NONE

C ---------------------------
C     Includes
C
      INCLUDE 'BMG_stencils.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER  iPL, Nx, Ny, Nz, NStncl
      REAL*8   Q(Nx,Ny,Nz), QF(Nx,Ny,Nz), QF2(Nx,Ny),
     &         SO(Nx,Ny,Nz,NStncl)

C ---------------------------
C     Local Declarations:
C
      INTEGER  i, j, kk, maxXY, Nxm2, Nym2

C ========================================================================

      Nxm2 = Nx-2
      Nym2 = Ny-2

      maxXY  = Nxm2*Nym2

      IF ( NStncl.EQ.14 ) THEN
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,kk)
C$OMP& SHARED(iPL,maxXY,Nxm2,Q,QF,QF2,SO)
         DO kk=0,maxXY-1
            i = mod(kk,Nxm2)+2
            j = (kk/Nxm2)+2
            QF2(i,j) = QF(i,j,iPL)
     &           + SO(i,j,iPL,kb)*Q(i,j,iPL-1)
     &           + SO(i,j,iPL,kbw)*Q(i-1,j,iPL-1)
     &           + SO(i,j+1,iPL,kbnw)*Q(i-1,j+1,iPL-1)
     &           + SO(i,j+1,iPL,kbn)*Q(i,j+1,iPL-1)
     &           + SO(i+1,j+1,iPL,kbne)*Q(i+1,j+1,iPL-1)
     &           + SO(i+1,j,iPL,kbe)*Q(i+1,j,iPL-1)
     &           + SO(i+1,j,iPL,kbse)*Q(i+1,j-1,iPL-1)
     &           + SO(i,j,iPL,kbs)*Q(i,j-1,iPL-1)
     &           + SO(i,j,iPL,kbsw)*Q(i-1,j-1,iPL-1)
     &           + SO(i,j,iPL+1,kbe)*Q(i-1,j,iPL+1)
     &           + SO(i,j+1,iPL+1,kbse)*Q(i-1,j+1,iPL+1)
     &           + SO(i,j+1,iPL+1,kbs)*Q(i,j+1,iPL+1)
     &           + SO(i+1,j+1,iPL+1,kbsw)*Q(i+1,j+1,iPL+1)
     &           + SO(i+1,j,iPL+1,kbw)*Q(i+1,j,iPL+1)
     &           + SO(i,j,iPL+1,kb)*Q(i,j,iPL+1)
     &           + SO(i+1,j,iPL+1,kbnw)*Q(i+1,j-1,iPL+1)
     &           + SO(i,j,iPL+1,kbn)*Q(i,j-1,iPL+1)
     &           + SO(i,j,iPL+1,kbne)*Q(i-1,j-1,iPL+1)
         ENDDO
C$OMP END PARALLEL DO
      ELSE
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,kk)
C$OMP& SHARED(iPL,maxXY,Nxm2,Q,QF,QF2,SO)
         DO kk=0,maxXY-1
            i = mod(kk,Nxm2)+2
            j = (kk/Nxm2)+2
            QF2(i,j) = QF(i,j,iPL)
     &           + SO(i,j,iPL,kb)*Q(i,j,iPL-1)
     &           + SO(i,j,iPL+1,kb)*Q(i,j,iPL+1)
         ENDDO
C$OMP END PARALLEL DO
         !
      ENDIF

C =====================

      RETURN
      END

