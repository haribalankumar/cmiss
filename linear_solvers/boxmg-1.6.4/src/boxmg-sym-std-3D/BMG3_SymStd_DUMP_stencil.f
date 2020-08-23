      SUBROUTINE BMG3_SymStd_DUMP_stencil( 
     &                BMG_IOFLAG, SO, Nx, Ny, Nz, IFD, NStncl, kg, NOG
     &                )

C ==========================================================================
C
C   BMG3_SymStd_DUMP_stencil.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_DUMP_stencil.f outputs the stencil interactively
C   to the screen or dumps the entire stencil into an ASCII text
C   file. Thus it is intended for debugging the stencil generation
C   on reasonable size grids.
C
C   - assumes that SO(p_SO) has been used in the call, where
C     p_SO is returned from BMG3key.
C   - assumes that you have provided the value of NStncl
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2003/01/30 (JDM)
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
      REAL*8  SO(Nx,Ny,Nz,NStncl)
      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)

C --------------------------
C     Local Declarations:
C
      INTEGER   i, j, k
      CHARACTER ANS*1, FILEo*30

C ==========================================================================

      IF ( BMG_IOFLAG(iBMG3_OUT_STENCIL_TTY) ) THEN 

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

            IF (NStncl.EQ.4) THEN

               WRITE(*,*)
               WRITE(*,50) -SO(i,j+1,k,kps)
               WRITE(*,*)
         
               WRITE(*,50) -SO(i,j,k+1,kb)
               WRITE(*,55) -SO(i,j,k,kpw),SO(i,j,k,kp),-SO(i+1,j,k,kpw)
               WRITE(*,60) -SO(i,j,k,kb)
               
               WRITE(*,*)
               WRITE(*,50) -SO(i,j,k,kps)
               WRITE(*,*)
               
            ELSEIF (NStncl.EQ.14) THEN
                  
               WRITE(*,52) -SO(i,j+1,k+1,kbse),
     &                     -SO(i,j+1,k+1,kbs),
     &                     -SO(i+1,j+1,k+1,kbsw)
               WRITE(*,55) -SO(i,j+1,k,kpnw),
     &                     -SO(i,j+1,k,kps),
     &                     -SO(i+1,j+1,k,kpsw)
               WRITE(*,62) -SO(i,j+1,k,kbnw),
     &                     -SO(i,j+1,k,kbn),
     &                     -SO(i+1,j+1,k,kbne)
         
               WRITE(*,52) -SO(i,j,k+1,kbe),
     &                     -SO(i,j,k+1,kb),
     &                     -SO(i+1,j,k+1,kbw)
               WRITE(*,55) -SO(i,j,k,kpw),
     &                      SO(i,j,k,kp),
     &                     -SO(i+1,j,k,kpw)
               WRITE(*,62) -SO(i,j,k,kbw),
     &                     -SO(i,j,k,kb),
     &                     -SO(i+1,j,k,kbe)
               
               WRITE(*,52) -SO(i,j,k+1,kbne),
     &                     -SO(i,j,k+1,kbn),
     &                     -SO(i+1,j,k+1,kbnw)
               WRITE(*,55) -SO(i,j,k,kpsw),
     &                     -SO(i,j,k,kps),
     &                     -SO(i+1,j,k,kpnw)
               WRITE(*,62) -SO(i,j,k,kbsw),
     &                     -SO(i,j,k,kbs),
     &                     -SO(i+1,j,k,kbse)

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

         IF (NStncl.EQ.4) THEN
            
            DO i=2, Nx-1
               DO j=2, Ny-1
                  DO k=2, Nz-1
                     WRITE(10,700) i, j, k,
     &                            -SO(i,j+1,k,kps),
     &                            -SO(i,j,k+1,kb),
     &                            -SO(i,j,k,kpw),
     &                             SO(i,j,k,kp),
     &                            -SO(i+1,j,k,kpw),
     &                            -SO(i,j,k,kb),
     &                            -SO(i,j,k,kps)
                  ENDDO
               ENDDO
            ENDDO
        
         ELSE
             
            DO k=2, Nz-1
               DO j=2, Ny-1
                  DO i=2, Nx-1
                     WRITE(10,800) i, j, k, 'N',
     &                               -SO(i,j+1,k+1,kbse),
     &                               -SO(i,j+1,k+1,kbs),
     &                               -SO(i+1,j+1,k+1,kbsw),
     &                               -SO(i,j+1,k,kpnw),
     &                               -SO(i,j+1,k,kps),
     &                               -SO(i+1,j+1,k,kpsw),
     &                               -SO(i,j+1,k,kbnw),
     &                               -SO(i,j+1,k,kbn),
     &                               -SO(i+1,j+1,k,kbne)
         
                     WRITE(10,800) i, j, k, 'O',
     &                               -SO(i,j,k+1,kbe),
     &                               -SO(i,j,k+1,kb),
     &                               -SO(i+1,j,k+1,kbw),
     &                               -SO(i,j,k,kpw),
     &                               SO(i,j,k,kp),
     &                               -SO(i+1,j,k,kpw),
     &                               -SO(i,j,k,kbw),
     &                               -SO(i,j,k,kb),
     &                               -SO(i+1,j,k,kbe)
               
                     WRITE(10,800) i, j, k, 'S',
     &                               -SO(i,j,k+1,kbne),
     &                               -SO(i,j,k+1,kbn),
     &                               -SO(i+1,j,k+1,kbnw),
     &                               -SO(i,j,k,kpsw),
     &                               -SO(i,j,k,kps),
     &                               -SO(i+1,j,k,kpnw),
     &                               -SO(i,j,k,kbsw),
     &                               -SO(i,j,k,kbs),
     &                               -SO(i+1,j,k,kbse)

                  ENDDO
               ENDDO
            ENDDO
               
   
         ENDIF

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

 700  FORMAT (I2,2(',',1X,I2),1P,7(',',1X,E14.7))
 800  FORMAT (I2,2(',',1X,I2),1X,A,1P,9(',',1X,E12.4))

      RETURN
      END
         
