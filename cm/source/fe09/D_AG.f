      REAL*8 FUNCTION D_AG(PARAMTYPE,nb,nhx,nr,ns,
     '  D_TG,D_ZG,PPG,TG,ZG)

C#### Function: D_AG
C###  Type: REAL*8
C###  Description:
C###    D_AG evaluates the derivative wrt PARAMTYPE parameters, of the
C###    integrand of the domain integral in the virtual work
C###    equation for all 4 coordinate systems and all 8 equation types.
C###    Note: PPG(ns,1+nj) has derivative of basis function wrt x(nj),
C###    where x coord type is as defined by DXIX.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER nb,nhx,nr,ns
      REAL*8 D_TG(3,3),D_ZG(NHM,NUM),PPG(64,4),TG(3,3),ZG(NHM,NUM)
      CHARACTER PARAMTYPE*(*)
!     Local Variables
      INTEGER mi,ni,NITB,NU1(0:3)
      REAL*8 AG,CLZ,CMZ,CSLZ,CSMZ,CZ,D1,D2,D3,D_CLZ,D_CMZ,D_CSLZ,
     '  D_CSMZ,D_CZ,D_D1,D_D2,D_D3,D_DLA,D_DMA,D_DRA,D_DTA,D_E1,
     '  D_E2,D_G1,D_PA,D_RA,D_RZ,D_SLZ,D_SMZ,D_SZ,D_TA,DLA,DMA,DRA,
     '  DTA,E1,E2,G1,PA,PGA,PGG,RA,RZ,SLZ,SMZ,SZ,TA
      DATA NU1/1,2,4,7/

      D_AG=0.0d0
      IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
        D_AG=AG(nb,nhx,nr,ns,PPG,D_TG,ZG)

      ELSEIF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
        NITB=NIT(nb)

        IF(ITYP10(nr).EQ.1) THEN       !rectangular cartesian coords
          DO mi=1,NITB
            DO ni=1,NITB
              D_AG=D_AG+(D_TG(mi,ni)*ZG(nhx,NU1(ni))
     '                  +TG(mi,ni)*D_ZG(nhx,NU1(ni)))*PPG(ns,1+mi)
            ENDDO
          ENDDO

        ELSE IF(ITYP10(nr).EQ.2) THEN  !cylindrical polar coords
          RZ  =  ZG(1,1)
          D_RZ=D_ZG(1,1)
          PGG=PPG(ns,1)
          DO ni=1,NITB
            D1  =  ZG(1,NU1(ni))
            D_D1=D_ZG(1,NU1(ni))
            D2  =  ZG(2,NU1(ni))
            D_D2=D_ZG(2,NU1(ni))
            D3  =  ZG(3,NU1(ni))
            D_D3=D_ZG(3,NU1(ni))
            DO mi=1,NITB
              DRA  =  ZG(1,NU1(mi)) !check shells etc(KTYP51(nr).GT.3)
              D_DRA=D_ZG(1,NU1(mi))
              DTA  =  ZG(2,NU1(mi)) !  "    "     "    "     "     "
              D_DTA=D_ZG(2,NU1(mi))
              PGA=PPG(ns,1+mi)
              IF(nhx.EQ.1) THEN
                D_AG=D_AG+D_TG(mi,ni)*(D1*PGA+D2*RZ*DTA*PGG)
     '                   +TG(mi,ni)*(D_D1*PGA
     '                     +(D_D2*RZ*DTA+D2*D_RZ*DTA+D2*RZ*D_DTA)*PGG)
              ELSE IF(nhx.EQ.2) THEN
                D_AG=D_AG+D_TG(mi,ni)*(-D1*DTA/RZ*PGG
     '                                 +D2*(PGA-DRA/RZ*PGG))
     '                   +TG(mi,ni)*(D_D2*PGA
     '                            +(-D_D1*DTA/RZ
     '                              -D1*D_DTA/RZ
     '                              +D1*DTA*D_RZ/(RZ*RZ)
     '                              -D_D2*DRA/RZ
     '                              -D2*D_DRA/RZ
     '                              +D2*DRA*D_RZ/(RZ*RZ))*PGG)
              ELSE IF(nhx.EQ.3) THEN
                D_AG=D_AG+(D_TG(mi,ni)*D3+TG(mi,ni)*D_D3)*PGA
              ENDIF
            ENDDO
          ENDDO

        ELSE IF(ITYP10(nr).EQ.3) THEN  !spherical polar coords
          RZ  =  ZG(1,1)
          D_RZ=D_ZG(1,1)
          CZ  = DCOS(ZG(3,1))
          D_CZ=-DSIN(ZG(3,1))*D_ZG(3,1)
          SZ  = DSIN(ZG(3,1))
          D_SZ= DCOS(ZG(3,1))*D_ZG(3,1)
          PGG=PPG(ns,1)
          DO ni=1,NITB
            D1  =  ZG(1,NU1(ni))
            D_D1=D_ZG(1,NU1(ni))
            D2  =  ZG(2,NU1(ni))
            D_D2=D_ZG(2,NU1(ni))
            D3  =  ZG(3,NU1(ni))
            D_D3=D_ZG(3,NU1(ni))
            DO mi=1,NITB
              RA  =  ZG(1,NU1(mi))
              D_RA=D_ZG(1,NU1(mi))
              TA  =  ZG(2,NU1(mi))
              D_TA=D_ZG(2,NU1(mi))
              PA  =  ZG(3,NU1(mi))
              D_PA=D_ZG(3,NU1(mi))
              PGA=PPG(ns,1+mi)
              IF(nhx.EQ.1) THEN
                D_AG=D_AG+D_TG(mi,ni)*(D1*PGA+
     '                                 RZ*(D2*CZ*CZ*TA+D3*PA)*PGG)
     '                   +TG(mi,ni)*(D_D1*PGA
     '                             +(D_RZ*D2*CZ*CZ*TA
     '                              +RZ*D_D2*CZ*CZ*TA
     '                              +RZ*D2*2.0d0*CZ*D_CZ*TA
     '                              +RZ*D2*CZ*CZ*D_TA
     '                              +D_RZ*D3*PA
     '                              +RZ*D_D3*PA
     '                              +RZ*D3*D_PA)*PGG)
              ELSE IF(nhx.EQ.2) THEN
                D_AG=D_AG+D_TG(mi,ni)*(TA*(D3*SZ/CZ-D1/RZ)*PGG
     '                         +D2*(PGA-(RA/RZ-SZ/CZ*PA)*PGG))
     '                   +TG(mi,ni)*(D_D2*PGA
     '                             +(D_TA*D3*SZ/CZ
     '                              +TA*D_D3*SZ/CZ
     '                              +TA*D3*D_SZ/CZ
     '                              -TA*D3*SZ*D_CZ/(CZ*CZ)
     '                              -D_TA*D1/RZ
     '                              -TA*D_D1/RZ
     '                              +TA*D1*D_RZ/(RZ*RZ)
     '                              -D_D2*RA/RZ
     '                              -D2*D_RA/RZ
     '                              +D2*RA*D_RZ/(RZ*RZ)
     '                              +D_D2*SZ/CZ*PA
     '                              +D2*D_SZ/CZ*PA
     '                              -D2*SZ*D_CZ/(CZ*CZ)*PA
     '                              +D2*SZ/CZ*D_PA)*PGG)
              ELSE IF(nhx.EQ.3) THEN
                D_AG=D_AG+D_TG(mi,ni)*(D3*(PGA-RA/RZ*PGG)
     '                         -(D1*PA/RZ+D2*CZ*SZ*TA)*PGG)
     '                   +TG(mi,ni)*(D_D3*PGA
     '                            +(-D_D3*RA/RZ
     '                              -D3*D_RA/RZ
     '                              +D3*RA*D_RZ/(RZ*RZ)
     '                              -D_D1*PA/RZ
     '                              -D1*D_PA/RZ
     '                              +D1*PA*D_RZ/(RZ*RZ)
     '                              -D_D2*CZ*SZ*TA
     '                              -D2*D_CZ*SZ*TA
     '                              -D2*CZ*D_SZ*TA
     '                              -D2*CZ*SZ*D_TA)*PGG)
              ENDIF
            ENDDO
          ENDDO

        ELSE IF(ITYP10(nr).EQ.4) THEN  !prolate spheroidal coords
          SLZ  =DSINH(ZG(1,1))
          D_SLZ=DCOSH(ZG(1,1))*D_ZG(1,1)
          SMZ  =DSIN (ZG(2,1))
          D_SMZ=DCOS (ZG(2,1))*D_ZG(2,1)
          CLZ  = DSQRT(1.0d0+SLZ*SLZ)
          D_CLZ= SLZ*D_SLZ/DSQRT(1.0d0+SLZ*SLZ)
          CMZ  = DSQRT(1.0d0-SMZ*SMZ)
          D_CMZ=-SMZ*D_SMZ/DSQRT(1.0d0-SMZ*SMZ)
          CSLZ  =CLZ/SLZ
          D_CSLZ=D_CLZ/SLZ-CLZ*D_SLZ/(SLZ*SLZ)
          CSMZ  =CMZ/SMZ
          D_CSMZ=D_CMZ/SMZ-CMZ*D_SMZ/(SMZ*SMZ)
          G1  =SLZ*SLZ+SMZ*SMZ
          D_G1=2.0d0*(SLZ*D_SLZ+SMZ*D_SMZ)
          E1  =CLZ*SLZ/G1
          D_E1=D_CLZ*SLZ/G1+CLZ*D_SLZ/G1-CLZ*SLZ*D_G1/(G1*G1)
          E2  =CMZ*SMZ/G1
          D_E2=D_CMZ*SMZ/G1+CMZ*D_SMZ/G1-CMZ*SMZ*D_G1/(G1*G1)
          PGG=PPG(ns,1)
          DO ni=1,NITB
            D1  =  ZG(1,NU1(ni))
            D_D1=D_ZG(1,NU1(ni))
            D2  =  ZG(2,NU1(ni))
            D_D2=D_ZG(2,NU1(ni))
            D3  =  ZG(3,NU1(ni))
            D_D3=D_ZG(3,NU1(ni))
            DO mi=1,NITB
              DLA  =  ZG(1,NU1(mi))
              D_DLA=D_ZG(1,NU1(mi))
              DMA  =  ZG(2,NU1(mi))
              D_DMA=D_ZG(2,NU1(mi))
              DTA  =  ZG(3,NU1(mi))
              D_DTA=D_ZG(3,NU1(mi))
              PGA=PPG(ns,1+mi)
              IF(nhx.EQ.1) THEN
                D_AG=D_AG+D_TG(mi,ni)*(D1*(PGA-(E1*DLA+E2*DMA)*PGG)
     '             +D2*(E1*DMA-E2*DLA)*PGG+D3*E1*SMZ*SMZ*DTA*PGG)
     '                   +TG(mi,ni)*(D_D1*PGA
     '                            +(-D_D1*E1*DLA-D1*D_E1*DLA-D1*E1*D_DLA
     '                              -D_D1*E2*DMA-D1*D_E2*DMA-D1*E2*D_DMA
     '                              +D_D2*E1*DMA+D2*D_E1*DMA+D2*E1*D_DMA
     '                              -D_D2*E2*DLA-D2*D_E2*DLA-D2*E2*D_DLA
     '                              +D_D3*E1*SMZ*SMZ*DTA
     '                              +D3*D_E1*SMZ*SMZ*DTA
     '                              +D3*E1*2.0d0*SMZ*D_SMZ*DTA
     '                              +D3*E1*SMZ*SMZ*D_DTA)*PGG)
              ELSE IF(nhx.EQ.2) THEN
                D_AG=D_AG+D_TG(mi,ni)*(D1*(E2*DLA-E1*DMA)*PGG
     '             +D2*(PGA-(E1*DLA+E2*DMA)*PGG)+D3*E2*SLZ*SLZ*DTA*PGG)
     '                   +TG(mi,ni)*(D_D2*PGA
     '                             +(D_D1*E2*DLA+D1*D_E2*DLA+D1*E2*D_DLA
     '                              -D_D1*E1*DMA-D1*D_E1*DMA-D1*E1*D_DMA
     '                              -D_D2*E1*DLA-D2*D_E1*DLA-D2*E1*D_DLA
     '                              -D_D2*E2*DMA-D2*D_E2*DMA-D2*E2*D_DMA
     '                              +D_D3*E2*SLZ*SLZ*DTA
     '                              +D3*D_E2*SLZ*SLZ*DTA
     '                              +D3*E2*2.0d0*SLZ*D_SLZ*DTA
     '                              +D3*E2*SLZ*SLZ*D_DTA)*PGG)
              ELSE IF(nhx.EQ.3) THEN
                D_AG=D_AG+D_TG(mi,ni)*(-(D1*CSLZ+D2*CSMZ)*DTA*PGG
     '             +D3*(PGA-(CSLZ*DLA+CSMZ*DMA)*PGG))
     '                   +TG(mi,ni)*(D_D3*PGA
     '                +(-D_D1*CSLZ*DTA-D1*D_CSLZ*DTA-D1*CSLZ*D_DTA
     '                  -D_D2*CSMZ*DTA-D2*D_CSMZ*DTA-D2*CSMZ*D_DTA
     '                  -D_D3*CSLZ*DLA-D3*D_CSLZ*DLA-D3*CSLZ*D_DLA
     '                  -D_D3*CSMZ*DMA-D3*D_CSMZ*DMA-D3*CSMZ*D_DMA)*PGG)
              ENDIF
            ENDDO
          ENDDO
        ELSE IF(ITYP10(nr).EQ.5) THEN  !oblate spheroidal coords
        ENDIF
      ENDIF

      RETURN
      END

C FE60 Functions
C ==============

