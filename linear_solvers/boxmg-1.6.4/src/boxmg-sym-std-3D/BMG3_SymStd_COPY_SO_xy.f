      SUBROUTINE BMG3_SymStd_COPY_SO_xy(
     &                SO, SO_xy, iPL, Nx, Ny, Nz,
     &                NStncl_3D, NStncl_2D,
     &                BMG_IOFLAG, BMG_iPARMS  
     &                )


C =======================================================================
C 
C   BMG3_SymStd_COPY_SO_xy
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_COPY_SO_xy copies the iPL{th}-(x,y) plane of the 
C   3D stencil into SO_xy.
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
      REAL*8   SO(Nx,Ny,Nz,NStncl_3D), SO_xy(Nx,Ny,NStncl_2D)

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)
      INTEGER  BMG_iPARMS(NBMG_iPARMS)

C ---------------------------
C     Local Declarations:
C
      INTEGER  i, j, kk, maxXY

C ========================================================================

      maxXY  = Nx*Ny

      IF ( NStncl_3D.EQ.14 .AND. NStncl_2D.EQ.5 ) THEN
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,kk)
C$OMP& SHARED(iPL,maxXY,Nx,SO,SO_xy)
         DO kk=0,maxXY-1
            i = mod(kk,Nx)+1
            j = (kk/Nx)+1
            SO_xy(i,j,ko)  = SO(i,j,iPL,kp)
            SO_xy(i,j,kw)  = SO(i,j,iPL,kpw)
            SO_xy(i,j,ks)  = SO(i,j,iPL,kps)
            SO_xy(i,j,ksw) = SO(i,j,iPL,kpsw)
            SO_xy(i,j,knw) = SO(i,j,iPL,kpnw)
         ENDDO
C$OMP END PARALLEL DO
         !
      ELSEIF ( NStncl_3D.EQ.4 .AND. NStncl_2D.EQ.3 ) THEN
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,kk)
C$OMP& SHARED(iPL,maxXY,Nx,SO,SO_xy)
         DO kk=0,maxXY-1
            i = mod(kk,Nx)+1
            j = (kk/Nx)+1
            SO_xy(i,j,ko) = SO(i,j,iPL,kp)
            SO_xy(i,j,kw) = SO(i,j,iPL,kpw)
            SO_xy(i,j,ks) = SO(i,j,iPL,kps)
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
         
         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,21)
         RETURN

         !
      ENDIF

C ========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_COPY_SO_xy.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END
