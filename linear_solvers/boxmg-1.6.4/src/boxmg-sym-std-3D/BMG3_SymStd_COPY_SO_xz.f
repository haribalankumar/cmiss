      SUBROUTINE BMG3_SymStd_COPY_SO_xz(
     &                SO, SO_xz, iPL, Nx, Ny, Nz,
     &                NStncl_3D, NStncl_2D,
     &                BMG_IOFLAG, BMG_iPARMS
     &                )


C =======================================================================
C 
C   BMG3_SymStd_COPY_SO_xz
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_COPY_SO_xz copies the iPL{th}-(x,y) plane of the 
C   3D stencil into SO_xz.
C
C =======================================================================

      IMPLICIT   NONE

      INCLUDE    'BMG_constants.h'
      INCLUDE    'BMG_stencils.h'
      INCLUDE    'BMG_parameters.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER  Nx, iPL, Ny, Nz, NStncl_2D, NStncl_3D
      REAL*8   SO(Nx,Ny,Nz,NStncl_3D), SO_xz(Nx,Nz,NStncl_2D)

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)
      INTEGER  BMG_iPARMS(NBMG_iPARMS)

C ---------------------------
C     Local Declarations:
C
      INTEGER  i, k, kk, maxXZ

C ========================================================================

      maxXZ  = Nx*Nz

      IF ( NStncl_3D.EQ.14 .AND. NStncl_2D.EQ.5 ) THEN
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,k,kk)
C$OMP& SHARED(iPL,maxXZ,Nx,SO,SO_xz)
         DO kk=0,maxXZ-1
            i = mod(kk,Nx)+1
            k = (kk/Nx)+1
            SO_xz(i,k,ko)  = SO(i,iPL,k,kp)
            SO_xz(i,k,kw)  = SO(i,iPL,k,kpw)
            SO_xz(i,k,ks)  = SO(i,iPL,k,kb)
            SO_xz(i,k,ksw) = SO(i,iPL,k,kbw)
            SO_xz(i,k,knw) = SO(i,iPL,k,kbe)
         ENDDO
C$OMP END PARALLEL DO
         !
      ELSEIF ( NStncl_3D.EQ.4 .AND. NStncl_2D.EQ.3 ) THEN
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,k,kk)
C$OMP& SHARED(iPL,maxXZ,Nx,SO,SO_xz)
         DO kk=0,maxXZ-1
            i = mod(kk,Nx)+1
            k = (kk/Nx)+1
            SO_xz(i,k,ko) = so(i,iPL,k,kp)
            SO_xz(i,k,kw) = so(i,iPL,k,kpw)
            SO_xz(i,k,ks) = so(i,iPL,k,kb)
         ENDDO
C$OMP END PARALLEL DO
         !
      ELSE
         !
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'Inconsitent number of stencil entries: '
            WRITE(*,510) 'HAVE: NStncl_3D = ', NStncl_3D
            WRITE(*,520) 'HAVE: NStncl_2D = ', NStncl_2D
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,23)
         RETURN
         !
      ENDIF

C ========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_COPY_SO_xz.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END
