      SUBROUTINE BMG2_SymStd_DUMP_stencil( 
     &                BMG_IOFLAG, SO, Nx, Ny, IFD, NStncl, kg, NOG
     &                )

C ==========================================================================
C
C   BMG2_SymStd_DUMP_stencil.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_DUMP_stencil.f outputs the stencil interactively
C   to the screen (or SOMEDAY could dump the entire stencil into 
C   an ASCII text file). Thus it is intended for debugging the 
C   stencil generation on reasonable size grids. 
C
C   - assumes that SO(p_SO) has been used in the call, where
C     p_SO is returned by BMG2_SymStd_GET_pointers.
C   - assumes that you have provided the value of NStncl
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2003/07/08 (JDM)
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
      INTEGER IFD, kg, NOG, NStncl, Nx, Ny, Nz
      REAL*8  SO(Nx,Ny,NStncl)
      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)

C --------------------------
C     Local Declarations:
C
      INTEGER   i, j
      CHARACTER ANS*1, FILEo*30

C ==========================================================================

      IF ( BMG_IOFLAG(iBMG2_OUT_STENCIL_TTY) ) THEN 

         WRITE(*,25) '**** INTERACTIVE STENCIL EVALUATION *****'
         WRITE(*,26) 'Finest grid index = ', NOG
         WRITE(*,28) 'Current grid index = ', kg
         WRITE(*,28) 'Index ranges include ghost points ... '
         WRITE(*,29) '[1,Nx] = [1,',Nx,']'
         WRITE(*,29) '[1,Ny] = [1,',Ny,']'
         WRITE(*,*) 

 10      CONTINUE

            WRITE(*,40) 'Enter index value (i,j) .... '
            READ(*,*) i, j

            IF ( i.LT.1 .OR. i.GT.Nx 
     &           .OR. j.LT.1 .OR. j.GT.Ny ) THEN 
               WRITE(*,*) 'error index out of bounds'
               GOTO 10
            ENDIF

            IF (NStncl.EQ.3) THEN
         
               WRITE(*,50) -SO(i,j+1,ks)
               WRITE(*,55) -SO(i,j,kw),SO(i,j,ko),-SO(i+1,j,kw)
               WRITE(*,60) -SO(i,j,ks)
               
            ELSEIF (NStncl.EQ.5) THEN
                  
               WRITE(*,52) -SO(i,j+1,knw),
     &                     -SO(i,j+1,ks),
     &                     -SO(i+1,j+1,ksw)
               WRITE(*,55) -SO(i,j,kw),
     &                      SO(i,j,ko),
     &                     -SO(i+1,j,kw)
               WRITE(*,62) -SO(i,j,ksw),
     &                     -SO(i,j,ks),
     &                     -SO(i+1,j,knw)
               
            ENDIF

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

         FILEo = 'output/'
     &           //'Stencil.'//CHAR(ICHAR('0')+kg) 

         OPEN ( 10, FILE=FILEo, STATUS='UNKNOWN' ) 

         IF ( NStncl.EQ.3 ) THEN

            DO i=2, Nx-1
               DO j=2, Ny-1
                  
                  WRITE(10,700) i, j,
     &                          -SO(i,j+1,ks),
     &                          -SO(i,j,kw),
     &                           SO(i,j,ko),
     &                          -SO(i+1,j,kw),
     &                          -SO(i,j,ks)

               ENDDO
            ENDDO

         ELSEIF ( NStncl.EQ.5 ) THEN

            DO j=2, Ny-1
               DO i=2, Nx-1

                  WRITE(10,800) i, j, 
     &                          -SO(i,j+1,knw),
     &                          -SO(i,j+1,ks),
     &                          -SO(i+1,j+1,ksw),
     &                          -SO(i,j,kw),
     &                           SO(i,j,ko),
     &                          -SO(i+1,j,kw),
     &                          -SO(i,j,ksw),
     &                          -SO(i,j,ks),
     &                          -SO(i+1,j,knw)

               ENDDO
            ENDDO
            
         ENDIF

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


 700  FORMAT (I2,',',1X,I2,1P,5(',',1X,E14.7))
 800  FORMAT (I2,',',1X,I2,1P,9(',',1X,E14.7))
 
      RETURN
      END
         
