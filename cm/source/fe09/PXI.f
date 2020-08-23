
      REAL*8 FUNCTION PXI(IBT,IDO,INP,nb,nu,XI,XE)

C#### Function: PXI
C###  Type: REAL*8
C###  Description:
C###    PXI interpolates nodal array XE at XI with basis functions as
C###    originally defined (jtyp5=1) or transformed to monomial form
C###    (jtyp5=2).

C#### Variable: nu
C###  Type: INTEGER
C###  Set_up: PXI
C###  Description:
C###    <HTML>
C###    nu is the calculated derivative loop variable.
C###    <PRE>
C###    nu=1 function
C###    nu=2 derivative wrt Xi1
C###    nu=3 derivative wrt Xi1^2
C###    nu=4 derivative wrt Xi2
C###    nu=5 derivative wrt Xi2^2
C###    nu=6 derivative wrt Xi1 Xi2
C###    nu=7 derivative wrt Xi3
C###    nu=8 derivative wrt Xi3^2
C###    nu=9 derivative wrt Xi1 Xi3
C###    nu=10 derivative wrt Xi2 Xi3
C###    nu=11 derivative wrt Xi1 XI2 Xi3
C###    </PRE>
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'b14.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,nu
      REAL*8 XE(NSM),XI(3)
!     Local Variables
      INTEGER i1,i2,i3,IP1,IP2,IP3,ni,nk,nn,ns
      REAL*8 COSHX,CSS,D,DES,PSI1,PSI2,PSI2_HERMITE,PSI2_XI,
     '  PSI5,RAD,SINHX,SS,SUM1,SUM2,THETA,XL(4)
      LOGICAL SECTOR
      SAVE RAD,COSHX,SINHX

      IF(JTYP5.EQ.1) THEN      !conventional form of interpolant
        PXI=0.0d0
        SECTOR=.FALSE.
        DO ni=1,NIT(nb)
          IF(IBT(1,ni).EQ.5.OR.IBT(1,ni).EQ.6) SECTOR=.TRUE.
        ENDDO !ni
        IF(SECTOR) THEN
          ns=0
          DO nn=1,NNT(nb)
            DO nk=1,NKT(nn,nb)
              ns=ns+1
              PXI=PXI+PSI5(IBT,IDO,INP,nb,nu,nk,nn,XI)*XE(ns)
            ENDDO !nk
          ENDDO !nn
        ELSE
          IF(IBT(1,1).LE.2.OR.IBT(1,1).EQ.9) THEN
            ns=0
            DO nn=1,NNT(nb)
              DO nk=1,NKT(nn,nb)
                ns=ns+1
                PXI=PXI+PSI1(IBT,IDO,INP,nb,nu,nk,nn,XI)*XE(ns)
              ENDDO
            ENDDO
          ELSE IF(IBT(1,1).EQ.3.AND.IBT(2,1).EQ.4) THEN
!         !Special hermite simplex.  Added AJP 5-6-93
            ns=0
            DO nn=1,NNT(nb)
              DO nk=1,NKT(nn,nb)
                ns=ns+1
                PXI=PXI+PSI2_HERMITE(IDO,INP,nb,nu,nk,nn,XI)*XE(ns)
              ENDDO
            ENDDO
          ELSE IF(IBT(1,1).EQ.3.AND.IBT(3,1).EQ.2) THEN
!         !Simplex using Xi parametric coords  Added CS 16-9-97
            DO nn=1,NNT(nb)
              PXI=PXI+PSI2_XI(nb,nu,nn,XI)*XE(nn)
            ENDDO
          ELSE IF(IBT(1,1).EQ.3) THEN
            XL(2)=XI(1)
            XL(3)=XI(2)
            XL(4)=XI(3)
            XL(1)=1.0D0-XL(2)-XL(3)-XL(4)
            DO nn=1,NNT(nb)
              PXI=PXI+PSI2(IBT,INP,nb,nu,nn,XL)*XE(nn)
            ENDDO
C          ELSE IF(IBT(1,1).EQ.5) THEN
C            nk=0
C            DO i3=1,IBT(2,3)+1
C              DO i2=1,IBT(2,2)+1
C                DO i1=1,IBT(2,1)+1
C                  nk=nk+1
C                  PXI=PXI+PSI3(i1,i2,i3,NIT(nb),nu,XI)*XE(nk)
C                ENDDO
C              ENDDO
C            ENDDO
          ENDIF
        ENDIF
      ELSE IF(JTYP5.EQ.2) THEN !monomial form of interpolant
        IP1=IBT(2,1)
        IP2=IBT(2,2)
        IP3=IBT(2,3)
        ns=1
        PXI=0.0D0
        DO i3=0,IP3
          SUM2=0.0D0
          DO i2=0,IP2
            SUM1=XE(ns)
            DO i1=1,IP1
              ns=ns+1
              SUM1=SUM1+XE(ns)*XI(1)**IP1
            ENDDO
            SUM2=SUM2+SUM1*XI(2)**IP2
          ENDDO
          PXI=PXI+SUM2*XI(3)**IP3
        ENDDO
      ENDIF

      IF(JTYP10.GE.2) THEN
        IF(ITYP10(1).EQ.2) THEN
          RAD=DSQRT(PXI)
          IF(nu.EQ.1) THEN
            PXI=RAD
          ELSE
            PXI=PXI/(2.0D0*RAD)
          ENDIF
        ELSE IF(ITYP10(1).EQ.3) THEN
          RAD=PXI**(1.0D0/3.0D0)
          IF(nu.EQ.1) THEN
            PXI=RAD
          ELSE
            PXI=PXI/(3.0D0*RAD*RAD)
          ENDIF
        ELSE IF(ITYP10(1).EQ.4) THEN
          IF(JTYP10.EQ.2) THEN
            SS=PXI/(FOCUS*FOCUS)
            SINHX=DSQRT(SS)
            COSHX=DSQRT(1.0D0+SS)
            IF(nu.EQ.1) THEN
              PXI=DLOG(COSHX+SINHX)
            ELSE
              PXI=PXI/(2.0D0*FOCUS*FOCUS*SINHX*COSHX)
            ENDIF
          ELSE IF(JTYP10.EQ.3) THEN
            CSS=PXI/FOCUS**3
            DES=CSS*CSS-4.0D0/27.0D0
            IF(DES.GT.0.0D0) THEN
              D=((CSS+DSQRT(DES))/2.0D0)**(1.0D0/3.0D0)
              COSHX=D+1.0D0/(3.0D0*D)
            ELSE
              THETA=DACOS(CSS*DSQRT(27.0D0)/2.0D0)
              COSHX=2.0D0/DSQRT(3.0D0)*DCOS(THETA/3.0D0)
            ENDIF
C            SINHX=DSQRT(COSHX*COSHX-1.0D0)
C !!! Added abs to avoid problems with small negative numbers
C !!! AJP 11/7/95
            SINHX=DSQRT(DABS(COSHX*COSHX-1.0D0))
            IF(nu.EQ.1) THEN
              PXI=DLOG(COSHX+SINHX)
            ELSE
              PXI=PXI/((3.0D0*COSHX*COSHX-1.0D0)*SINHX)/FOCUS**3
            ENDIF
          ENDIF
        ELSE IF(ITYP10(1).EQ.5) THEN
        ENDIF
      ENDIF

      RETURN
      END


