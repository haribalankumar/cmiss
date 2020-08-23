      REAL*8 FUNCTION PSI8(NAN,nb,nu,XI)

C#### Function: PSI8
C###  Type: REAL*8
C###  Description:
C###    <HTML>
C###    PSI8 evaluates auxiliary basis functions at Xi.
C###    <PRE>
C###    NABTYP(nb)=1 Legendre polynomials
C###               2 Fourier  basis functions
C###               3 Pressure basis functions
C###    </PRE> </HTML>

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NAN(NIM),nb,nu
      REAL*8 XI(3)
!     Local Variables
      INTEGER IPU(11,3),ni
      REAL*8 PAF,PAL0,PAL1,PAL2,PAP4,PH3

      DATA IPU/1,2,3,1,1,2,1,1,2,1,2,
     '         1,1,1,2,3,2,1,1,1,2,2,
     '         1,1,1,1,1,1,2,3,2,2,2/

      PSI8=1.0D0
      DO ni=1,NIT(nb)
        IF(NABTYP(nb).EQ.1) THEN      !Legendre basis
          IF(NAN(ni).EQ.0) PSI8=PSI8*PAL0(IPU(nu,ni))
          IF(NAN(ni).EQ.1) PSI8=PSI8*PAL1(IPU(nu,ni),XI(ni))
          IF(NAN(ni).GE.2) PSI8=PSI8*PAL2(IPU(nu,ni),NAN(ni),XI(ni))
        ELSE IF(NABTYP(nb).EQ.2) THEN !Fourier  basis
          PSI8=PSI8*PAF(IPU(nu,ni),NAN(ni),XI(ni))
        ELSE IF(NABTYP(nb).EQ.3) THEN !Pressure basis
          IF(NAN(ni).LT.0) PSI8=0.0d0
          IF(NAN(ni).EQ.0) PSI8=PSI8*PAL0(IPU(nu,ni))
          IF(NAN(ni).EQ.1) PSI8=PSI8*PAL1(IPU(nu,ni),XI(ni))
          IF(NAN(ni).EQ.2) PSI8=PSI8*PAL2(IPU(nu,ni),NAN(ni),XI(ni))
! new MPN 18Mar98: changed cubic Hermite fns
          IF(NAN(ni).EQ.3) PSI8=PSI8*PH3(1,1,IPU(nu,ni),XI(ni))
          IF(NAN(ni).EQ.4) PSI8=PSI8*PAP4(IPU(nu,ni),XI(ni))
          IF(NAN(ni).EQ.5) PSI8=PSI8*PH3(1,2,IPU(nu,ni),XI(ni))
          IF(NAN(ni).EQ.6) PSI8=PSI8*PH3(2,2,IPU(nu,ni),XI(ni))
          IF(NAN(ni).EQ.7) PSI8=PSI8*PH3(2,1,IPU(nu,ni),XI(ni))
! old
c          IF(NAN(ni).EQ.3) PSI8=PSI8*PH3(2,1,IPU(nu,ni),XI(ni))
c          IF(NAN(ni).EQ.4) PSI8=PSI8*PAP4(IPU(nu,ni),XI(ni))
c          IF(NAN(ni).EQ.5) PSI8=PSI8*PH3(1,2,IPU(nu,ni),XI(ni))
c          IF(NAN(ni).EQ.6) PSI8=PSI8*PH3(2,2,IPU(nu,ni),XI(ni))
        ENDIF
      ENDDO !ni

      RETURN
      END


