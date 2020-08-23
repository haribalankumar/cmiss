      SUBROUTINE OPEIGE(EIGV,FIX,YP,ERROR,*)

C#### Subroutine: OPEIGE
C###  Description:
C###    OPEIGE outputs eigenvalues.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'moda00.cmn'
!     Parameter List
      REAL*8 EIGV(NOM,NTM,2),YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER MO,no,nr,nrc,nx,ny
      REAL*8 EVALI,EVALR

      CALL ENTERS('OPEIGE',*9999)

      nrc=2 !Temporary AJP 30-11-94
      nx=1 !Temporary

      WRITE(OP_STRING,'(/'' Eigenvalues and eigenvectors:'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      no=0
      nr=1 !needs fixing
      NTMAX=NTM !common variable was not set
      DO ny=1,NYT(nrc,1,nx)
        IF(.NOT.FIX(ny,1)) THEN
          no=no+1
          IF(no.LE.NTMAX) THEN
            EVALR=YP(ny,3)
            EVALI=YP(ny,4)
            WRITE(OP_STRING,'(/1X,I4,'') Real: '',D12.4,'
     '        //''' Imaginary: '',D12.4)') no,EVALR,EVALI
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(1X,10D12.4)')
     '        (EIGV(MO,no,2),MO=1,NOT(nrc,1,nr,nx))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDDO

      CALL EXITS('OPEIGE')
      RETURN
 9999 CALL ERRORS('OPEIGE',ERROR)
      CALL EXITS('OPEIGE')
      RETURN 1
      END


