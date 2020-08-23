      SUBROUTINE PETROV_PREP(NITB,CG,CONT0,DNORM,DZXI,GU,GUxDZ,GU2DZ,
     '  MDDEN1,MDDEN2,MGAL,MPET,MATERIALPET,
     '  NEEDDERIV,ERROR,*)

C#### Subroutine: PETROV_PREP
C###  Description:
C###    PETROV_PREP calculates multipliers necessary in the calculation
C###    of derivative weights when using Petrov-Galerkin FEM.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'load00.cmn'
      INCLUDE 'petg00.inc'
!     Parameter List
      INTEGER NITB
      REAL*8 CG(NMM),CONT0(2),DNORM,DZXI(3),GU(3,3),GUxDZ(3),
     '  GU2DZ(3),MDDEN1,MDDEN2,MGAL,MPET
      CHARACTER ERROR*(*)
      LOGICAL MATERIALPET,NEEDDERIV
!     Local Variables
      INTEGER mi,ni
      REAL*8 CONTDEN,CONTEXP1,CONTEXP2,CONTNUM,DAVG,DENOM1,DENOM2,
     '  DENOM3,DET,DINVAVG,D2NORM,GL(3,3),MDCONT,PECLET,PECLETMAX,PPET,
     '  SUM,TEMP,TEMP1,XINORM

      DATA PECLETMAX/0d0/

      CALL ENTERS('PETROV_PREP',*9999)

C     Estimate of coupling in direction of propagation.
C      IF(NITB.EQ.3) THEN
C        DAVG=DMAX1(GU(1,1),GU(2,2),GU(3,3))
C      ELSE IF(NITB.EQ.2) THEN
C        DAVG=DMAX1(GU(1,1),GU(2,2))
C      ELSE IF(NITB.EQ.1) THEN
C        DAVG=GU(1,1)
C      ENDIF
      SUM=0.0d0
      DO ni=1,NITB
        SUM=SUM+GU(ni,ni)
      ENDDO
      DAVG=SUM/NITB
C     Estimate of coupling reciprocal in direction of propagation.
      CALL INVERT(NITB,GU,GL,DET)
      SUM=0.0d0
      DO ni=1,NITB
        SUM=SUM+GL(ni,ni)
      ENDDO
      DINVAVG=SUM/NITB

      IF(MATERIALPET) THEN !deriv weights in material coord dirn.
C       Norm of weight derivative direction.
        D2NORM=0.0d0
        DO ni=1,NITB
          D2NORM=D2NORM+GUxDZ(ni)*GUxDZ(ni)
        ENDDO
CC       Estimate of coupling in direction of propagation.
C        DAVG=1.0d0/DINVAVG
      ELSE !deriv weights in element xi coord based dirn
C       Norm of first xi derivatives.
        XINORM=0.0d0
        DO ni=1,NITB
          XINORM=XINORM+DZXI(ni)*DZXI(ni)
        ENDDO
      ENDIF !deriv weight direction

C     Calculate multipliers for Galerkin and derivative weights.
      TEMP1=(1.0d0-PETSMOOTH)*CG(2)*CG(2)
      TEMP=PETSMOOTH*CG(1)*CG(1)
      PPET=PETRATIO*CG(2)*FACTOR
      IF(MATERIALPET) THEN
C***    For material coord based limited derivative weights:
        DENOM1=DSQRT(TEMP1*DNORM+TEMP)
        DENOM2=DSQRT(TEMP1*D2NORM+TEMP*DAVG)
        DENOM3=PETLIMIT*PPET*DENOM1+DENOM2
C        DENOM3=PETLIMIT*PPET+DAVG
        PECLET=CG(2)*FACTOR*CONT0(2)/DSQRT(DAVG)
C        PECLET=CG(2)*FACTOR*DENOM1*CONT0(2)/DENOM2
        CONTEXP1=DEXP(-PECLET*CONT0(1))
        CONTEXP2=DEXP(-PECLET)
        CONTNUM=1.0d0-CONTEXP1
        CONTDEN=1.0d0-CONTEXP2
        MDCONT=1.0d0!+PECLET*(CONT0(1)*CONTEXP1/CONTNUM-CONTEXP2/CONTDEN)
        IF(CONTDEN.GT.0.0d0) THEN
          CONT0(1)=CONTNUM/CONTDEN
        ELSE
          CONT0(1)=0.0D0
        ENDIF
        MGAL=1.0d0
        IF(DENOM2.NE.0.0d0) THEN
          MPET=CONT0(1)*CG(2)*PPET*DENOM1/(DENOM3*DENOM2)
        ELSE !can't choose a direction for derivative
          MPET=0.0d0
        ENDIF
        IF(NEEDDERIV) THEN
          IF(DENOM2.NE.0.0d0) THEN
            MDDEN1=(MDCONT/DENOM1-PETLIMIT*PPET/DENOM3)*TEMP1/DENOM1
            MDDEN2=(-1.0d0/DENOM3-MDCONT/DENOM2)*TEMP1/DENOM2
C            MDDEN1=0.0d0
C            MDDEN0=-TEMP1/(DENOM2*DENOM2)
          ELSE !can't choose a direction for derivative
            MDDEN1=0.0d0
            MDDEN2=0.0d0
          ENDIF
C         Find temporary vector, GU2DZ = GU * GU * DZXI.
          DO ni=1,NITB
            SUM=0.0d0
            DO mi=1,NITB
              SUM=SUM+GUxDZ(mi)*GU(mi,ni)
            ENDDO !mi
            GU2DZ(ni)=SUM
          ENDDO !ni
        ENDIF !NEEDDERIV
C      ELSE IF(MATERIALPET) THEN
CC***    For material coord based full derivative weights:
C        PPET=PPET*CG(2)
C        DENOM1=DSQRT(TEMP1*DNORM+TEMP)
C        DENOM2=DSQRT((TEMP1*D2NORM+TEMP*DAVG)*DAVG)
C        DENOM3=PPET*FACTOR*DENOM1+DENOM2
C        PECLET=CG(2)*FACTOR*CONT0(2)/DSQRT(DAVG)
C        IF(PECLET.NE.0.0d0) THEN
C          CONT0(1)=(1.0d0-DEXP(-PECLET*CONT0(1)))/(1.0d0-DEXP(-PECLET))
C        ELSE
C          CONT0(1)=0.0D0
C        ENDIF
C        MGAL=DENOM2/DENOM3
C        IF(DENOM2.NE.0.0d0) THEN
C          MPET=CONT0(1)*PPET/DENOM3
C        ELSE !can't choose a direction for derivative
C          MPET=0.0d0
C        ENDIF
C        IF(NEEDDERIV) THEN
C          IF(DENOM2.NE.0.0d0) THEN
C            MDGAL1=-TEMP1*PPET*FACTOR/(DENOM1*DENOM3)
C            MDGAL2=TEMP1*DAVG*(1.0d0/DENOM2-1.0d0/DENOM3)/DENOM2
C            MDDEN1=-TEMP1*PPET*FACTOR/(DENOM1*DENOM3)
C            MDDEN2=-TEMP1*DAVG/(DENOM2*DENOM3)
C          ELSE !can't choose a direction for derivative
C            MDGAL1=0.0d0
C            MDGAL2=0.0d0
C            MDDEN1=0.0d0
C            MDDEN2=0.0d0
C          ENDIF
CC         Find temporary vector, GU2DZ = GU * GU * DZXI.
C          DO ni=1,NITB
C            SUM=0.0d0
C            DO mi=1,NITB
C              SUM=SUM+GUxDZ(mi)*GU(mi,ni)
C            ENDDO !mi
C            GU2DZ(ni)=SUM
C          ENDDO !ni
C        ENDIF !NEEDDERIV
      ELSE
C***    For element local coord based derivative weights:
        DENOM1=DSQRT(TEMP1*XINORM+TEMP*DINVAVG)
        DENOM2=DSQRT(TEMP1*DNORM+TEMP)
C        DENOM3=(DINVAVG*PETLIMIT*PPET+1.0d0)*DENOM1
        DENOM3=PETLIMIT*PPET*DENOM1+DENOM2
        PECLET=CG(2)*FACTOR*CONT0(2)/DSQRT(DAVG)
        IF(PECLET.NE.0.0d0) THEN
          CONT0(1)=(1.0d0-DEXP(-PECLET*CONT0(1)))/(1.0d0-DEXP(-PECLET))
        ELSE
          CONT0(1)=0.0D0
        ENDIF
        MGAL=1.0d0
        IF(DENOM3.NE.0.0d0) THEN
          MPET=CONT0(1)*CG(2)*PPET/DENOM3
        ELSE !can't choose a direction for derivative
          MPET=0.0d0
        ENDIF
        IF(NEEDDERIV) THEN
          IF(DENOM1.NE.0.0d0.AND.DENOM2.NE.0.0d0) THEN
C            MDDEN0=-TEMP1*(DINVAVG*PETLIMIT*PPET+1.0d0)/(DENOM3*DENOM1)
C            MDDEN1=0.0d0
            MDDEN1=-TEMP1*PETLIMIT*PPET/(DENOM3*DENOM1)
            MDDEN2=-TEMP1/(DENOM3*DENOM2)
          ELSE !can't choose a direction for derivative
            MDDEN1=0.0d0
            MDDEN2=0.0d0
          ENDIF
        ENDIF !NEEDDERIV
      ENDIF !deriv weight direction

      IF(DOP) THEN
        DENOM1=DSQRT(DNORM)
        DENOM2=DSQRT(D2NORM)
        IF(DENOM2.GT.0d0) THEN
          PECLET=CG(2)*FACTOR*DENOM1/DENOM2
          IF(PECLET.GT.PECLETMAX) THEN
            PECLETMAX=PECLET
            WRITE(OP_STRING,'(''Largest Peclet Number='',E12.3)') PECLET
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('PETROV_PREP')
      RETURN
 9999 CALL ERRORS('PETROV_PREP',ERROR)
      CALL EXITS('PETROV_PREP')
      RETURN 1
      END



