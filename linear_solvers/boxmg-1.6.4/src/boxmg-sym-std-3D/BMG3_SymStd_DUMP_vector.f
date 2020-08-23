      SUBROUTINE BMG3_SymStd_DUMP_vector( 
     &                BMG_IOFLAG, Q, Nx, Ny, Nz, kg, NOG,
     &                FILE_dir, FILE_base, GHOST
     &                )

C ==========================================================================
C
C   BMG3_SymStd_DUMP_vector.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_DUMP_vector.f outputs a component of a vector
C   interactively to the screen or dumps the entire vector into
C   an ASCII text file. Thus it is intended for debugging on
C   reasonably size grids.
C
C   - assumes that Q(p_U) has been used in the call, where
C     p_U is returned from BMG3_SymStd_GET_pointers.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2003/08/27 (JDM)
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
      INCLUDE 'BMG_parameters.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER   kg, NOG, Nx, Ny, Nz
      REAL*8    Q(Nx,Ny,Nz)
      LOGICAL   BMG_IOFLAG(NBMG_IOFLAG), GHOST
      CHARACTER FILE_base*(*), FILE_dir*(*)

C --------------------------
C     Local Declarations:
C
      INTEGER   i, j, k
      INTEGER   IEND, ISTART, JEND, JSTART, KEND, KSTART
      CHARACTER ANS*1, FILEo*(30)

C ==========================================================================

      IF (.FALSE.) THEN 

         WRITE(*,25) '**** INTERACTIVE STENCIL EVALUATION *****'
         WRITE(*,26) 'Finest grid index = ', NOG
         WRITE(*,28) 'Current grid index = ', kg
         WRITE(*,28) 'Index ranges include ghost points ... '
         WRITE(*,29) '[1,Nx] = [1,',Nx,']'
         WRITE(*,29) '[1,Ny] = [1,',Ny,']'
         WRITE(*,29) '[1,Nz] = [1,',Nz,']'
         WRITE(*,*) 

 10      CONTINUE

            WRITE(*,40) 'Enter index value (i,j,k) .... '
            READ(*,*) i, j, k

            IF ( i.LT.1 .OR. i.GT.Nx 
     &           .OR. j.LT.1 .OR. j.GT.Ny 
     &           .OR. k.LT.1 .OR. k.GT.Nz ) THEN
               WRITE(*,*) 'error index out of bounds'
               GOTO 10
            ENDIF

            WRITE(*,50) Q(i,j,k)

 20      CONTINUE

         WRITE(*,45) 'Another stencil (Y/N) .... '
         READ(*,'(A)') ANS
      
         IF (INDEX('Yy',ANS).NE.0) THEN
            GOTO 10
         ELSE IF (INDEX('Nn',ANS).EQ.0) THEN
            WRITE(*,*) 'Error .... '
            GOTO 20
         ENDIF

      ELSE

         IF ( GHOST ) THEN 
            ISTART = 1
            JSTART = 1
            KSTART = 1
            IEND   = Nx
            JEND   = Ny
            KEND   = Nz
         ELSE
            ISTART = 2
            JSTART = 2
            KSTART = 2
            IEND   = Nx-1
            JEND   = Ny-1
            KEND   = Nz-1
         ENDIF

         IF ( FILE_dir.EQ.'' ) THEN
            FILEo = FILE_base//'.'//CHAR(ICHAR('0')+kg) 
         ELSE
            FILEo = FILE_dir//FILE_base//'.'//CHAR(ICHAR('0')+kg) 
         ENDIF

         OPEN ( 10, FILE=FILEo, STATUS='UNKNOWN' ) 

         DO k=KSTART, KEND
            DO j=JSTART, JEND
               DO i=ISTART, IEND
                  WRITE(10,700) i, j, k, Q(i,j,k)
               ENDDO
            ENDDO
         ENDDO
         
         CLOSE( 10 )

      ENDIF


 25   FORMAT (/,/,6X,A,/)
 26   FORMAT (/,8X,A,I2)
 28   FORMAT (8X,A,I2,/)
 29   FORMAT (10X,A,I5,A)

 30   FORMAT (/,4X,A)
 40   FORMAT (4X,A,$)
 45   FORMAT (/,4X,A,$)

 50   FORMAT (/,24X,1P,E12.5)
 52   FORMAT (/,8X,1P,3(E12.5,4X))
 55   FORMAT (8X,1P,3(E12.5,4X))
 60   FORMAT (24X,1P,E12.5,/)
 62   FORMAT (8X,1P,3(E12.5,4X),/)


 70   FORMAT (/,24X,1P,E12.5)
 75   FORMAT (16X,1P,2(E12.5,4X))
 80   FORMAT (24X,1P,E12.5)
 85   FORMAT (24X,1P,E12.5,/)

 700  FORMAT (I2,2(',',1X,I2),1P,',',1X,E16.8)

      RETURN
      END
         
