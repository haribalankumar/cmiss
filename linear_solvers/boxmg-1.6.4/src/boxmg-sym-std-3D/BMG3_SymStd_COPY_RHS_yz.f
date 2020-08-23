      SUBROUTINE BMG3_SymStd_COPY_RHS_yz( 
     &                       SO, Q, QF, QF2, iPL, Nx, Ny, Nz, NStncl
     &                       )

C =======================================================================
C 
C   BMG3_SymStd_COPY_RHS_yz
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_COPY_RHS_yz creates the right hand side vector, QF,
C   from the iPL{th}-(y,z) plane of the current 3D vector.
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
      REAL*8   Q(Nx,Ny,Nz), QF(Nx,Ny,Nz), QF2(Ny,Nz),
     &         SO(Nx,Ny,Nz,NStncl)

C ---------------------------
C     Local Declarations:
C
      INTEGER  j, k, kk, maxYZ, Nxm2, Nym2, Nzm2
C ========================================================================

      Nym2 = Ny-2
      Nzm2 = Nz-2

      maxYZ = Nym2*Nzm2

      IF ( NStncl.EQ.14 ) THEN
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(j,k,kk)
C$OMP& SHARED(iPL,maxYZ,Nym2,Q,QF,QF2,SO)
         DO kk=0,maxYZ-1
            j = mod(kk,Nym2)+2
            k = (kk/Nym2)+2
            !
            QF2(j,k) = QF(iPL,j,k)
     &           + SO(iPL,j+1,k,kpnw)*Q(iPL-1,j+1,k)
     &           + SO(iPL,j,k,kpw)*Q(iPL-1,j,k)
     &           + SO(iPL,j,k,kpsw)*Q(iPL-1,j-1,k)
     &           + SO(iPL,j+1,k,kbnw)*Q(iPL-1,j+1,k-1)
     &           + SO(iPL,j,k,kbw)*Q(iPL-1,j,k-1)
     &           + SO(iPL,j,k,kbsw)*Q(iPL-1,j-1,k-1)
     &           + SO(iPL,j+1,k+1,kbse)*Q(iPL-1,j+1,k+1)
     &           + SO(iPL,j,k+1,kbe)*Q(iPL-1,j,k+1)
     &           + SO(iPL,j,k+1,kbne)*Q(iPL-1,j-1,k+1)
     &           + SO(iPL+1,j+1,k,kpsw)*Q(iPL+1,j+1,k)
     &           + SO(iPL+1,j,k,kpw)*Q(iPL+1,j,k)
     &           + SO(iPL+1,j,k,kpnw)*Q(iPL+1,j-1,k)
     &           + SO(iPL+1,j+1,k,kbne)*Q(iPL+1,j+1,k-1)
     &           + SO(iPL+1,j,k,kbe)*Q(iPL+1,j,k-1)
     &           + SO(iPL+1,j,k,kbse)*Q(iPL+1,j-1,k-1)
     &           + SO(iPL+1,j+1,k+1,kbsw)*Q(iPL+1,j+1,k+1)
     &           + SO(iPL+1,j,k+1,kbw)*Q(iPL+1,j,k+1)
     &           + SO(iPL+1,j,k+1,kbnw)*Q(iPL+1,j-1,k+1)
            !
         ENDDO
C$OMP END PARALLEL DO
         !
      ELSE
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(j,k,kk)
C$OMP& SHARED(iPL,maxYZ,Nym2,Q,QF,QF2,SO)
         DO kk=0,maxYZ-1
            j = mod(kk,Nym2)+2
            k = (kk/Nym2)+2
            !
            QF2(j,k) = QF(iPL,j,k)
     &           + SO(iPL,j,k,kpw)*Q(iPL-1,j,k)
     &           + SO(iPL+1,j,k,kpw)*Q(iPL+1,j,k)
         ENDDO
C$OMP END PARALLEL DO
         !
      ENDIF

C =====================

      RETURN
      END

