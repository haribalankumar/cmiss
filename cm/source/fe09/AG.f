      REAL*8 FUNCTION AG(nb,nhx,nr,ns,PPG,TG,ZG)

C#### Function: AG
C###  Type: REAL*8
C###  Description:
C###    AG evaluates integrand of the domain integral in the virtual
C###    work equation for all 4 coordinate systems and all 8 equation
C###    types.  Note: PPG(ns,1+nj) has derivative of basis function
C###    wrt x(nj), where x coord type is as defined by DXIX.

C cpb 18/4/97 JTYP9 replaced with NJ_LOC
C**** 19Sep88: Changed IF(JTYP9.EQ.0) to IF(KTYP53(nr).EQ.1) and
C****          changed IF(JTYP9.EQ.1) to IF(KTYP53(nr).GE.2) since
C****          KTYP53(nr) now defines whether stresses in the
C****          constitutive law are referred to theta (KTYP53(nr)=1)
C****          or Nu (KTYP53(nr)>1).

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER nb,nhx,nr,ns
      REAL*8 PPG(64,4),TG(3,3),ZG(NHM,NUM)
!     Local Variables
      INTEGER mi,ni,NITB,NU1(0:3)
      REAL*8 CLZ,CMZ,CSLZ,CSMZ,CZ,
     '  D1,D2,D3,DLA,DMA,DRA,DTA,E1,E2,G1,
     '  PA,PGA,PGG,RA,RZ,SLZ,SMZ,SZ,TA
      DATA NU1/1,2,4,7/

      AG=0.0d0
      NITB=NIT(nb)
      IF(ITYP10(nr).EQ.1) THEN       !rectangular cartesian coords
        DO mi=1,NITB
          DO ni=1,NITB
            AG=AG+TG(mi,ni)*ZG(nhx,NU1(ni))*PPG(ns,1+mi)
          ENDDO !ni
        ENDDO !mi

      ELSE IF(ITYP10(nr).EQ.2) THEN  !cylindrical polar coords
        RZ=ZG(1,1)
        PGG=PPG(ns,1)
        DO ni=1,NITB
          D1=ZG(1,NU1(ni))
          D2=ZG(2,NU1(ni))
          D3=ZG(3,NU1(ni))
          DO mi=1,NITB
            DRA=ZG(1,NU1(mi)) !check this shells etc (KTYP51(nr).GT.3)
            DTA=ZG(2,NU1(mi)) !  "     "   "     "    "     "     "
            PGA=PPG(ns,1+mi)
            IF(nhx.EQ.1) THEN
              AG=AG+TG(mi,ni)*(D1*PGA+D2*RZ*DTA*PGG)
            ELSE IF(nhx.EQ.2) THEN
              AG=AG+TG(mi,ni)*(-D1*DTA/RZ*PGG+D2*(PGA-DRA/RZ*PGG))
            ELSE IF(nhx.EQ.3) THEN
              AG=AG+TG(mi,ni)*D3*PGA
            ENDIF
          ENDDO !mi
        ENDDO !ni

      ELSE IF(ITYP10(nr).EQ.3) THEN  !spherical polar coords
        RZ=ZG(1,1)
        CZ=DCOS(ZG(3,1))
        SZ=DSIN(ZG(3,1))
        PGG=PPG(ns,1)
        DO ni=1,NITB
          D1=ZG(1,NU1(ni))
          D2=ZG(2,NU1(ni))
          D3=ZG(3,NU1(ni))
          DO mi=1,NITB
            RA=ZG(1,NU1(mi))
            TA=ZG(2,NU1(mi))
            PA=ZG(3,NU1(mi))
            PGA=PPG(ns,1+mi)
            IF(nhx.EQ.1) THEN
              AG=AG+TG(mi,ni)*(D1*PGA+RZ*(D2*CZ*CZ*TA+D3*PA)*PGG)
            ELSE IF(nhx.EQ.2) THEN
              AG=AG+TG(mi,ni)*(TA*(D3*SZ/CZ-D1/RZ)*PGG
     '                        +D2*(PGA-(RA/RZ-SZ/CZ*PA)*PGG))
            ELSE IF(nhx.EQ.3) THEN
              AG=AG+TG(mi,ni)*(D3*(PGA-RA/RZ*PGG)
     '                       -(D1*PA/RZ+D2*CZ*SZ*TA)*PGG)
            ENDIF
          ENDDO !mi
        ENDDO !ni

      ELSE IF(ITYP10(nr).EQ.4) THEN  !prolate spheroidal coords
        SLZ=DSINH(ZG(1,1))
        SMZ=DSIN (ZG(2,1))
        CLZ=DSQRT(1.0d0+SLZ*SLZ)
        CMZ=DSQRT(1.0d0-SMZ*SMZ)
        CSLZ=CLZ/SLZ
        CSMZ=CMZ/SMZ
        G1=SLZ*SLZ+SMZ*SMZ
        E1=CLZ*SLZ/G1
        E2=CMZ*SMZ/G1
        PGG=PPG(ns,1)
        DO ni=1,NITB
          D1=ZG(1,NU1(ni))
          D2=ZG(2,NU1(ni))
          D3=ZG(3,NU1(ni))
          DO mi=1,NITB
            DLA=ZG(1,NU1(mi))
            DMA=ZG(2,NU1(mi))
            DTA=ZG(3,NU1(mi))
            PGA=PPG(ns,1+mi)
            IF(nhx.EQ.1) THEN
              AG=AG+TG(mi,ni)*(D1*(PGA-(E1*DLA+E2*DMA)*PGG)
     '             +D2*(E1*DMA-E2*DLA)*PGG+D3*E1*SMZ*SMZ*DTA*PGG)
            ELSE IF(nhx.EQ.2) THEN
              AG=AG+TG(mi,ni)*(D1*(E2*DLA-E1*DMA)*PGG
     '             +D2*(PGA-(E1*DLA+E2*DMA)*PGG)+D3*E2*SLZ*SLZ*DTA*PGG)
            ELSE IF(nhx.EQ.3) THEN
              AG=AG+TG(mi,ni)*(-(D1*CSLZ+D2*CSMZ)*DTA*PGG
     '             +D3*(PGA-(CSLZ*DLA+CSMZ*DMA)*PGG))
            ENDIF
          ENDDO !mi
        ENDDO !ni
      ELSE IF(ITYP10(nr).EQ.5) THEN  !oblate spheroidal coords
      ENDIF

      RETURN
      END


