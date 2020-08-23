      SUBROUTINE OPESTFMAT(NHST,IUNIT,MATRIX,VECTOR,MATNAME,VECTNAME,
     '  OPMATRIX,OPVECTOR,ERROR,*)

C#### Subroutine: OPESTFMAT
C###  Description:
C###    OPESTFMAT outputs element stiffness matrices and/or element
C###    stiffness vectors to a screen unit in a standard way.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NHST(2),IUNIT
      REAL*8 MATRIX(NHM*NSM,NHM*NSM),VECTOR(NHM*NSM)
      CHARACTER ERROR*(*),MATNAME*2,VECTNAME*2
      LOGICAL OPMATRIX,OPVECTOR
!     Local Variables
      INTEGER nhs1,nhs2

      CALL ENTERS('OPESTFMAT',*9999)

      DO nhs1=1,NHST(1)
        IF(OPVECTOR) THEN
          WRITE(OP_STRING,'(1X,A2,''('',I2,'')    ='',1X,D11.4)')
     '      VECTNAME,nhs1,VECTOR(nhs1)
          CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
        ENDIF
        IF(OPMATRIX) THEN
          WRITE(OP_STRING,'(1X,A2,''('',I2,'',1..)='',6(1X,D11.4),'
     '      //':/(12X,6(1X,D11.4)))') MATNAME,nhs1,(MATRIX(nhs1,nhs2),
     '      nhs2=1,NHST(2))
          CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
        ENDIF
      ENDDO !nhs1

      CALL EXITS('OPESTFMAT')
      RETURN
 9999 CALL ERRORS('OPESTFMAT',ERROR)
      CALL EXITS('OPESTFMAT')
      RETURN 1
      END


