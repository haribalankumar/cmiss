      SUBROUTINE OPMODA(LIST,nr,nx,EIGVAL,EIGVEC,VALUES,VECTORS,ERROR,*)

C#### Subroutine: OPMODA
C###  Description:
C###    OPMODA outputs modal values i.e. eigenvalues and eigenvectors.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
C      INCLUDE 'cmiss$reference:moda00.cmn'
!     Parameter List
      INTEGER LIST(0:NLISTM),nr,nx
      REAL*8 EIGVAL(NTM,2),EIGVEC(NOM,NTM,2)
      CHARACTER ERROR*(*)
      LOGICAL VALUES,VECTORS
!     Local Variables
      INTEGER mode,no,no_mode
      REAL*8 ANG_FREQUENCY,EVALR,FREQUENCY

      CALL ENTERS('OPMODA',*9999)

      WRITE(OP_STRING,'(/'' Region '',I2,'' (nx='',I1,'')'')') nr,nx
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      DO no_mode=1,LIST(0)
        mode=LIST(no_mode)
        WRITE(OP_STRING,'(/'' Mode # '',I5)') mode
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(VALUES) THEN
          EVALR=EIGVAL(mode,1)
          ANG_FREQUENCY=DSQRT(EVALR)
          FREQUENCY=ANG_FREQUENCY/(2.0d0*PI)
          WRITE(OP_STRING,'(/''   Mode Frequency    = '',D12.5,'
     '      //''' kHz'')') FREQUENCY
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( ''   Angular Frequency = '',D12.5,'
     '      //''' rad/s'')') ANG_FREQUENCY*1000.0d0
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( ''   Eigenvalue        = '',D12.5)') EVALR
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(VECTORS) THEN
          WRITE(OP_STRING,'(/''   Mode shape (eigenvector):'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C rgb 9/9/1999 If NOT is too large then the write to
C              OP_STRING will break.
          IF(NOT(2,1,nr,nx).GT.200) THEN
            DO no=1,NOT(2,1,nr,nx)
              WRITE(OP_STRING,'(D12.5)') EIGVEC(no,mode,1)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO
          ELSE
            WRITE(OP_STRING,'(4X,5(1X,D12.5),:/(4X,5(1X,D12.5)))')
     '        (EIGVEC(no,mode,1),no=1,NOT(2,1,nr,nx))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDDO !no_mode (mode)

      CALL EXITS('OPMODA')
      RETURN
 9999 CALL ERRORS('OPMODA',ERROR)
      CALL EXITS('OPMODA')
      RETURN 1
      END


