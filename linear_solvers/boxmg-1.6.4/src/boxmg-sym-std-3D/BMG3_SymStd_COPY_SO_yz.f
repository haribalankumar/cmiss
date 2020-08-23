      SUBROUTINE BMG3_SymStd_COPY_SO_yz(
     &                SO, SO_yz, iPL, Nx, Ny, Nz,
     &                NStncl_3D, NStncl_2D,
     &                BMG_IOFLAG, BMG_iPARMS
     &                )


C =======================================================================
C 
C   BMG3_SymStd_COPY_SO_yz
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_COPY_SO_yz copies the iPL{th}-(y,z) plane of the 
C   3D stencil into SO_yz.
C
C =======================================================================

      IMPLICIT   NONE

      INCLUDE    'BMG_constants.h'
      INCLUDE    'BMG_stencils.h'
      INCLUDE    'BMG_parameters.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER  iPL, NStncl_2D, NStncl_3D, Nx, Ny, Nz
      REAL*8   SO(Nx,Ny,Nz,NStncl_3D), SO_yz(Ny,Nz,NStncl_2D)

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)
      INTEGER  BMG_iPARMS(NBMG_iPARMS)

C ---------------------------
C     Local Declarations:
C
      INTEGER  j, k, kk, maxYZ

C ========================================================================

      maxYZ  = Ny*Nz

      IF ( NStncl_3D.EQ.14 .AND. NStncl_2D.EQ.5 ) THEN
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(j,k,kk)
C$OMP& SHARED(iPL,maxYZ,Ny,SO,SO_yz)
         DO kk=0,maxYZ-1
            j = mod(kk,Ny)+1
            k = (kk/Ny)+1
            SO_yz(j,k,ko)  = SO(iPL,j,k,kp)
            SO_yz(j,k,kw)  = SO(iPL,j,k,kps)
            SO_yz(j,k,ks)  = SO(iPL,j,k,kb)
            SO_yz(j,k,ksw) = SO(iPL,j,k,kbs)
            SO_yz(j,k,knw) = SO(iPL,j,k,kbn)
         ENDDO
C$OMP END PARALLEL DO
         !
      ELSEIF ( NStncl_3D.EQ.4 .AND. NStncl_2D.EQ.3 ) THEN
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(j,k,kk)
C$OMP& SHARED(iPL,maxYZ,Ny,SO,SO_yz)
         DO kk=1,maxYZ-1
            j = mod(kk,Ny)+1
            k = (kk/Ny)+1
            SO_yz(j,k,ko) = SO(iPL,j,k,kp)
            SO_yz(j,k,kw) = SO(iPL,j,k,kps)
            SO_yz(j,k,ks) = SO(iPL,j,k,kb)
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

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,22)
         RETURN
         !
      ENDIF

C ========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_COPY_SO_yz.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END
