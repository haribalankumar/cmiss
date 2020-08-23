      SUBROUTINE BMG2_SymStd_DUMP_restrict( 
     &                BMG_IOFLAG, CI, KFG, KCG, NOG,
     &                Nx, Ny, Nx_c, Ny_c
     &                )

C ==========================================================================
C
C   BMG2_SymStd_DUMP_restrict.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_DUMP_restrict.f outputs the interpolation weights
C   interactively to the screen or dumps the entire array
C   of interpolation weights into an ASCII text file. Thus it is 
C   intended for debugging the stencil generation on reasonable size 
C   grids.
C
C   - assumes that CI(p_CIC) has been used in the call, where
C     p_CIC is returned from BMG3_SymStd_GET_pointers.
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
      INTEGER  KCG, KFG, NOG, Nx, Nx_c, Ny, Ny_c

      REAL*8   CI(Nx_c,Ny_c,8)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C --------------------------
C     Local Declarations:
C
      INTEGER   ic, jc, IO_n
      CHARACTER ANS*1, FILEo*30

C ==========================================================================


      IF ( BMG_IOFLAG(iBMG2_OUT_RESTRICT_TTY) ) THEN

         WRITE(*,125) '**** INTERACTIVE RESTRICTION EVALUATION *****'
         WRITE(*,126) 'Finest grid index = ', NOG
         WRITE(*,126) 'Fine grid index   = ', KFG
         WRITE(*,128) 'Coarse grid index = ', KCG
         WRITE(*,127) 'Coarse-grid index ranges ' 
     &                //'(including ghost points) ... '
         WRITE(*,129) '[1,Nx_c] = [1,',Nx_c,']'
         WRITE(*,129) '[1,Ny_c] = [1,',Ny_c,']'
         WRITE(*,*) 
            
 10      CONTINUE   ! <<<<<<< BEGIN: Data entry

            WRITE(*,140) 'Enter the coarse-grid coordinates ... '
            READ(*,*) ic, jc

            IF ( ic.LT.2 .OR. ic.GT.Nx_c-1
     &               .OR. jc.LT.2 .OR. jc.GT.Ny_c-1 ) THEN

               WRITE(*,*) 'error index out of bounds'
               GOTO 10

            ENDIF

 30      CONTINUE      ! <<<<<<< END: Data entry

         WRITE(*,152) CI(ic,jc+1,LSE),
     &                CI(ic,jc+1,LB),
     &                CI(ic+1,jc+1,LSW)

         WRITE(*,155) CI(ic,jc,LR),
     &                rONE,
     &                CI(ic+1,jc,LL)

         WRITE(*,162) CI(ic,jc,LNE),
     &                CI(ic,jc,LA),
     &                CI(ic+1,jc,LNW)
            
 40      CONTINUE

            WRITE(*,145) 'Another restriction (Y/N) .... '
            READ(*,'(A)') ANS

            IF (INDEX('Yy',ANS).NE.0) THEN
               GOTO 10
            ELSE IF (INDEX('Nn',ANS).EQ.0) THEN
               WRITE(*,*) 'Error .... '
               GOTO 40
            ENDIF

 60      CONTINUE

      ELSE

         IO_n = 10

         FILEo = 'output/'
     &           //'Restrict.'//CHAR(ICHAR('0')+KCG) 

         OPEN ( IO_n, FILE=FILEo, STATUS='UNKNOWN' ) 
         
         DO jc=2, Ny_c-1
            DO ic=2, Nx_c-1

               WRITE(IO_n,800) ic, jc,
     &                         CI(ic,jc+1,LSE),
     &                         CI(ic,jc+1,LB),
     &                         CI(ic+1,jc+1,LSW),
     &                         CI(ic,jc,LR),
     &                         rONE,
     &                         CI(ic+1,jc,LL),
     &                         CI(ic,jc,LNE),
     &                         CI(ic,jc,LA),
     &                         CI(ic+1,jc,LNW)
 
            ENDDO
         ENDDO

         CLOSE( IO_n )

      ENDIF


 125  FORMAT (/,/,6X,A,/)
 126  FORMAT (/,8X,A,I2)
 127  FORMAT (8X,A,/)
 128  FORMAT (8X,A,I2,/)
 129  FORMAT (10X,A,I5,A)

 130  FORMAT (/,4X,A)
 140  FORMAT (4X,A,$)
 145  FORMAT (/,4X,A,$)

 150  FORMAT (/,24X,1P,E12.5)
 152  FORMAT (/,8X,1P,3(E12.5,4X))
 155  FORMAT (8X,1P,3(E12.5,4X))
 160  FORMAT (24X,1P,E12.5,/)
 162  FORMAT (8X,1P,3(E12.5,4X),/)


 170  FORMAT (/,24X,1P,E12.5)
 175  FORMAT (16X,1P,2(E12.5,4X))
 180  FORMAT (24X,1P,E12.5)
 185  FORMAT (24X,1P,E12.5,/)

 800  FORMAT (I2,',',1X,I2,1P,9(',',1X,E14.7))

      RETURN
      END
         
