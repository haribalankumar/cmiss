      REAL*8 FUNCTION HYPGREEN(IGREN,IMAT,CE,RAD,RJAC,XN,XNO,XR)

C#### Function: HYPGREEN
C###  Type: REAL*8
C###  Description:
C###    HYPGREEN evaluates the second derviative of the Green's function
C###    for the hypersingular integral equations. RJAC is the Jacobian
C###    obtained from any element subdivision or polar coordinate
C###    transformation.  RJAC=1 if no transformation has been used. If
C###    the generalised Laplace equation is being solved then IMAT
C###    determines whether the conductivity needs to be included in the
C###    HYPGREEN expression (depends on the coefficient of HYPGREEN in
C###    the calling routine).  IMAT=0 means no material properties to be
C###    included and IMAT=1 means to include them.  This allows
C###    different regions to be coupled correctly.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IGREN,IMAT
      REAL*8 CE(NMM),RJAC,RAD,XN(NJM),XNO(*),XR(NJM)
!     Local Variables
      INTEGER nj
      REAL*8 NDOTNO,RDOTN,RDOTNO
      CHARACTER ERROR*10

      NDOTNO=0.0d0
      RDOTN=0.0d0
      RDOTNO=0.0d0
      DO nj=1,NJT
        NDOTNO=NDOTNO+XN(nj)*XNO(nj)
        RDOTN=RDOTN+XR(nj)*XN(nj)
        RDOTNO=RDOTNO+XR(nj)*XNO(nj)
      ENDDO
      IF(IGREN.EQ.1) THEN !Laplace eqtn 2D
        HYPGREEN=1.0d0/(2.0d0*PI*RAD*RAD)*(NDOTNO-2.0d0*RDOTN*RDOTNO/
     '    (RAD*RAD))*RJAC
        !no minus sign.  Found after great pain AJP, CPB 20-12-93
      ELSE IF(IGREN.EQ.2) THEN !Laplace eqtn 3D
        HYPGREEN=1.0d0/(4.0d0*PI*RAD*RAD*RAD)*
     '    (NDOTNO-3.0d0*RDOTN*RDOTNO/(RAD*RAD))*RJAC
      ELSE IF(IGREN.EQ.7) THEN !Generalised Laplace eqtn 2D
        !Depending on the coefficient of HYPGREEN the
        !conductivity is included or not (if the coefficient is just a
        !potential then it is not needed, but if it is the normal
        !derivative then it does need to be included).
        IF(IMAT.EQ.0) THEN
          HYPGREEN=1.0d0/(2.0d0*PI*RAD*RAD)*
     '      (NDOTNO-2.0d0*RDOTN*RDOTNO/(RAD*RAD))*RJAC
        ELSE
          HYPGREEN=1.0d0/(2.0d0*PI*RAD*RAD*CE(1))*
     '      (NDOTNO-2.0d0*RDOTN*RDOTNO/(RAD*RAD))*RJAC
        ENDIF
      ELSE IF(IGREN.EQ.8) THEN !Generalised Laplace eqtn 3D
        !Depending on the coefficient of HYPGREEN the
        !conductivity is included or not (if the coefficient is just a
        !potential then it is not needed, but if it is the normal
        !derivative then it does need to be included).
        IF(IMAT.EQ.0) THEN
          HYPGREEN=1.0d0/(4.0d0*PI*RAD*RAD*RAD)*
     '      (NDOTNO-3.0d0*RDOTN*RDOTNO/(RAD*RAD))*RJAC
        ELSE
          HYPGREEN=1.0d0/(4.0d0*PI*RAD*RAD*RAD*CE(1))*
     '      (NDOTNO-3.0d0*RDOTN*RDOTNO/(RAD*RAD))*RJAC
        ENDIF
      ELSEIF(IGREN.EQ.9) THEN
C       !Green's fn for Poisson eqtn, 2d cartesian, special source term
        !Source term is -div(k.grad(g)), so use Laplace fundamental soln
        HYPGREEN=1.0d0/(2.0d0*PI*RAD*RAD)*(NDOTNO-2.0d0*RDOTN*RDOTNO/
     '    (RAD*RAD))*RJAC
      ELSE IF(IGREN.EQ.10) THEN !Laplace eqtn 3D
C       !Green's fn for Poisson eqtn, 3d cartesian, special source term
        !Source term is -div(k.grad(g)), so use Laplace fundamental soln
        HYPGREEN=1.0d0/(4.0d0*PI*RAD*RAD*RAD)*
     '    (NDOTNO-3.0d0*RDOTN*RDOTNO/(RAD*RAD))*RJAC
      ELSE
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' >>ERROR:IGREN incorrect'')')
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
9999  RETURN
      END


