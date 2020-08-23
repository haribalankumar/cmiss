      SUBROUTINE ZEZG(JP,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*)

C#### Subroutine: ZEZG
C###  Description:
C###    ZEZG evaluates the Gauss point array ZG(nhx,nu) from the
C###    element array ZE(ns,nhx) at the current Gauss point ng.
C###    JP=0: return 1st derivs wrt Xi
C###    JP=1: return 1st derivs multiplied by dXIX
C###    JP=2: return 1st & 2nd derivs wrt Xi
C###    JP=3: return 1st & 2nd derivs multiplied by dXIX

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER JP,NBH(NHM),ng,NHE,nx
      REAL*8 DXIX(3,3),PG(NSM,NUM,NGM,NBM),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,nb,nh,nhx,ni,NSTNAT,nu,NU1(0:3),NU2(3,3)
C      INTEGER NU2(3,3)
      REAL*8 COSHZ,CSS,D,DES,DZDXI(3),RAD,SINHZ,SS,THETA,ZG_tmp(11)
!     Functions
      REAL*8 DDOT

      DATA NU1/1,2,4,7/
      DATA NU2/3,6,9,6,5,10,9,10,8/

      CALL ENTERS('ZEZG',*9999)

      DO nhx=1,NHE
        nh=NH_LOC(nhx,nx)
        nb=NBH(nh)
        NSTNAT=NST(nb)+NAT(nb)

C       The value of the dependent variable
        ZG(nhx,1)=DDOT(NSTNAT,PG(1,1,ng,nb),1,ZE(1,nhx),1)

        IF(JP.EQ.0) THEN      !return 1st derivs wrt Xi
          DO ni=1,NIT(nb)
            ZG(nhx,NU1(ni))=
     '        DDOT(NSTNAT,PG(1,NU1(ni),ng,nb),1,ZE(1,nhx),1)
          ENDDO

        ELSE IF(JP.EQ.1) THEN !return 1st derivs multiplied by dXIX
C MPN 18-Jul-95 Initialising ZG
C !!! KAT 21Nov97:  This is done solely to make example_721 run successfully
C               The fact that it uses these results indicates that the
C               example must be giving the wrong answer.
          DO nu=2,NUT(nb)
            ZG(nhx,nu)=0.0d0
          ENDDO
C         1st derivatives wrt Xi
          DO ni=1,NIT(nb)
            DZDXI(ni)=
     '        DDOT(NSTNAT,PG(1,NU1(ni),ng,nb),1,ZE(1,nhx),1)
          ENDDO
C         1st derivatives wrt X
          IF(NIT(nb).EQ.1) THEN
            ZG(nhx,2)=DZDXI(1)*DXIX(1,1)
          ELSE IF(NIT(nb).EQ.2) THEN
            ZG(nhx,2)=DZDXI(1)*DXIX(1,1) + DZDXI(2)*DXIX(2,1)
            ZG(nhx,4)=DZDXI(1)*DXIX(1,2) + DZDXI(2)*DXIX(2,2)
          ELSE IF(NIT(nb).EQ.3) THEN
            ZG(nhx,2)=DZDXI(1)*DXIX(1,1) +
     '        DZDXI(2)*DXIX(2,1) + DZDXI(3)*DXIX(3,1)
            ZG(nhx,4)=DZDXI(1)*DXIX(1,2) +
     '        DZDXI(2)*DXIX(2,2) + DZDXI(3)*DXIX(3,2)
            ZG(nhx,7)=DZDXI(1)*DXIX(1,3) +
     '        DZDXI(2)*DXIX(2,3) + DZDXI(3)*DXIX(3,3)
          ENDIF !nit

        ELSE IF(JP.EQ.2) THEN !return 1st & 2nd derivs wrt Xi
          DO nu=2,NUT(nb)
            ZG(nhx,nu)=DDOT(NSTNAT,PG(1,nu,ng,nb),1,ZE(1,nhx),1)
          ENDDO

        ELSE IF(JP.EQ.3) THEN !return 1st & 2nd derivs mult by dXIX (PJH 21/7/98)
C         Initialise array
          DO nu=2,NUT(nb)
            ZG(nhx,nu)=0.0d0
          ENDDO
C         1st derivatives wrt Xi
          DO ni=1,NIT(nb)
            DZDXI(ni)=
     '        DDOT(NSTNAT,PG(1,NU1(ni),ng,nb),1,ZE(1,nhx),1)
          ENDDO
C         1st derivatives wrt X
          IF(NIT(nb).EQ.1) THEN
            ZG(nhx,2)=DZDXI(1)*DXIX(1,1)
          ELSE IF(NIT(nb).EQ.2) THEN
            ZG(nhx,2)=DZDXI(1)*DXIX(1,1) + DZDXI(2)*DXIX(2,1)
            ZG(nhx,4)=DZDXI(1)*DXIX(1,2) + DZDXI(2)*DXIX(2,2)
          ELSE IF(NIT(nb).EQ.3) THEN
            ZG(nhx,2)=DZDXI(1)*DXIX(1,1) +
     '        DZDXI(2)*DXIX(2,1) + DZDXI(3)*DXIX(3,1)
            ZG(nhx,4)=DZDXI(1)*DXIX(1,2) +
     '        DZDXI(2)*DXIX(2,2) + DZDXI(3)*DXIX(3,2)
            ZG(nhx,7)=DZDXI(1)*DXIX(1,3) +
     '        DZDXI(2)*DXIX(2,3) + DZDXI(3)*DXIX(3,3)
C           2nd derivatives wrt Xi
            ZG(nhx, 3)=DDOT(NSTNAT,PG(1, 3,ng,nb),1,ZE(1,nhx),1)
            ZG(nhx, 5)=DDOT(NSTNAT,PG(1, 5,ng,nb),1,ZE(1,nhx),1)
            ZG(nhx, 6)=DDOT(NSTNAT,PG(1, 6,ng,nb),1,ZE(1,nhx),1)
            ZG(nhx, 8)=DDOT(NSTNAT,PG(1, 8,ng,nb),1,ZE(1,nhx),1)
            ZG(nhx, 9)=DDOT(NSTNAT,PG(1, 9,ng,nb),1,ZE(1,nhx),1)
            ZG(nhx,10)=DDOT(NSTNAT,PG(1,10,ng,nb),1,ZE(1,nhx),1)
            ZG(nhx,11)=DDOT(NSTNAT,PG(1,11,ng,nb),1,ZE(1,nhx),1)

C!!! Why recalculate all the dot products?
C new RJD 6/11/98
            ZG_tmp(3)=DDOT(NSTNAT,PG(1, 3,ng,nb),1,ZE(1,nhx),1)
            ZG_tmp(5)=DDOT(NSTNAT,PG(1, 5,ng,nb),1,ZE(1,nhx),1)
            ZG_tmp(6)=DDOT(NSTNAT,PG(1, 6,ng,nb),1,ZE(1,nhx),1)
            ZG_tmp(8)=DDOT(NSTNAT,PG(1, 8,ng,nb),1,ZE(1,nhx),1)
            ZG_tmp(9)=DDOT(NSTNAT,PG(1, 9,ng,nb),1,ZE(1,nhx),1)
            ZG_tmp(10)=DDOT(NSTNAT,PG(1,10,ng,nb),1,ZE(1,nhx),1)
            ZG_tmp(11)=DDOT(NSTNAT,PG(1,11,ng,nb),1,ZE(1,nhx),1)

C           2nd derivatives wrt Nu coords (RJD 8/10/98)
            ZG(nhx, 3)=ZG_tmp(3)*DXIX(1,1)*DXIX(1,1)   !11 deriv
     '        +ZG_tmp(6)*DXIX(1,1)*DXIX(2,1)
     '        +ZG_tmp(9)*DXIX(1,1)*DXIX(3,1)
     '        +ZG_tmp(6)*DXIX(2,1)*DXIX(1,1)
     '        +ZG_tmp(5)*DXIX(2,1)*DXIX(2,1)
     '        +ZG_tmp(10)*DXIX(2,1)*DXIX(3,1)
     '        +ZG_tmp(9)*DXIX(3,1)*DXIX(1,1)
     '        +ZG_tmp(10)*DXIX(3,1)*DXIX(2,1)
     '        +ZG_tmp(8)*DXIX(3,1)*DXIX(3,1)
            ZG(nhx, 5)=ZG_tmp(3)*DXIX(1,2)*DXIX(1,2)   !22 deriv
     '        +ZG_tmp(6)*DXIX(1,2)*DXIX(2,2)
     '        +ZG_tmp(9)*DXIX(1,2)*DXIX(3,2)
     '        +ZG_tmp(6)*DXIX(2,2)*DXIX(1,2)
     '        +ZG_tmp(5)*DXIX(2,2)*DXIX(2,2)
     '        +ZG_tmp(10)*DXIX(2,2)*DXIX(3,2)
     '        +ZG_tmp(9)*DXIX(3,2)*DXIX(1,2)
     '        +ZG_tmp(10)*DXIX(3,2)*DXIX(2,2)
     '        +ZG_tmp(8)*DXIX(3,2)*DXIX(3,2)
            ZG(nhx, 6)=ZG_tmp(3)*DXIX(1,1)*DXIX(1,2)   !12 deriv
     '        +ZG_tmp(6)*DXIX(1,1)*DXIX(2,2)
     '        +ZG_tmp(9)*DXIX(1,1)*DXIX(3,2)
     '        +ZG_tmp(6)*DXIX(2,1)*DXIX(1,2)
     '        +ZG_tmp(5)*DXIX(2,1)*DXIX(2,2)
     '        +ZG_tmp(10)*DXIX(2,1)*DXIX(3,2)
     '        +ZG_tmp(9)*DXIX(3,1)*DXIX(1,2)
     '        +ZG_tmp(10)*DXIX(3,1)*DXIX(2,2)
     '        +ZG_tmp(8)*DXIX(3,1)*DXIX(3,2)
            ZG(nhx, 8)=ZG_tmp(3)*DXIX(1,3)*DXIX(1,3)   !33 deriv
     '        +ZG_tmp(6)*DXIX(1,3)*DXIX(2,3)
     '        +ZG_tmp(9)*DXIX(1,3)*DXIX(3,3)
     '        +ZG_tmp(6)*DXIX(2,3)*DXIX(1,3)
     '        +ZG_tmp(5)*DXIX(2,3)*DXIX(2,3)
     '        +ZG_tmp(10)*DXIX(2,3)*DXIX(3,3)
     '        +ZG_tmp(9)*DXIX(3,3)*DXIX(1,3)
     '        +ZG_tmp(10)*DXIX(3,3)*DXIX(2,3)
     '        +ZG_tmp(8)*DXIX(3,3)*DXIX(3,3)
            ZG(nhx, 9)=ZG_tmp(3)*DXIX(1,1)*DXIX(1,3)   !13 deriv
     '        +ZG_tmp(6)*DXIX(1,1)*DXIX(2,3)
     '        +ZG_tmp(9)*DXIX(1,1)*DXIX(3,3)
     '        +ZG_tmp(6)*DXIX(2,1)*DXIX(1,3)
     '        +ZG_tmp(5)*DXIX(2,1)*DXIX(2,3)
     '        +ZG_tmp(10)*DXIX(2,1)*DXIX(3,3)
     '        +ZG_tmp(9)*DXIX(3,1)*DXIX(1,3)
     '        +ZG_tmp(10)*DXIX(3,1)*DXIX(2,3)
     '        +ZG_tmp(8)*DXIX(3,1)*DXIX(3,3)
            ZG(nhx,10)=ZG_tmp(3)*DXIX(1,2)*DXIX(1,3)   !23 deriv
     '        +ZG_tmp(6)*DXIX(1,2)*DXIX(2,3)
     '        +ZG_tmp(9)*DXIX(1,2)*DXIX(3,3)
     '        +ZG_tmp(6)*DXIX(2,2)*DXIX(1,3)
     '        +ZG_tmp(5)*DXIX(2,2)*DXIX(2,3)
     '        +ZG_tmp(10)*DXIX(2,2)*DXIX(3,3)
     '        +ZG_tmp(9)*DXIX(3,2)*DXIX(1,3)
     '        +ZG_tmp(10)*DXIX(3,2)*DXIX(2,3)
     '        +ZG_tmp(8)*DXIX(3,2)*DXIX(3,3)
          ENDIF !nit
        ENDIF !jp
      ENDDO

C KAT 14Nov97:  Old and slow:
C      DO nhx=1,NHE
C        nh=NH_LOC(nhx,nx)
C        nb=NBH(nh)
CC MPN 18-Jul-95 Initialising ZG
C        DO nu=1,NUT(nb)
C          ZG(nhx,nu)=0.0d0
C        ENDDO
C        DO ni=0,NIT(nb)
C          SUM=0.0D0
C          DO ns=1,NST(nb)+NAT(nb)
C            IF(jp.eq.0.or.ni.EQ.0) THEN
C              PGG=PG(ns,NU1(ni),ng,nb)
C            ELSE IF(JP.GT.0) THEN
C              IF(NIT(nb).EQ.1) THEN
C                PGG=PG(ns,2,ng,nb)*DXIX(1,1)
C              ELSE IF(NIT(nb).EQ.2) THEN
C                PGG=PG(ns,2,ng,nb)*DXIX(1,ni) +
C     '              PG(ns,4,ng,nb)*DXIX(2,ni)
C              ELSE IF(NIT(nb).EQ.3) THEN
C                PGG=PG(ns,2,ng,nb)*DXIX(1,ni) +
C     '              PG(ns,4,ng,nb)*DXIX(2,ni) +
C     '              PG(ns,7,ng,nb)*DXIX(3,ni)
CC              ELSE
CC                PGG=PGX(nb,ni,ns,DXIX,PG(1,1,ng,nb))
C              ENDIF
C            ENDIF
C            SUM=SUM+PGG*ZE(ns,nhx)
C          ENDDO
C          ZG(nhx,NU1(ni))=SUM
C        ENDDO
C      ENDDO

      IF(JTYP10.GE.2) THEN !isochoric elements
        nb=NBH(NH_LOC(1,nx))
        IF(ITYP11(1).EQ.2) THEN      !cylindrical polar coords
          RAD=DSQRT(ZG(1,1))
          ZG(1,1)=RAD
        ELSE IF(ITYP11(1).EQ.3) THEN !spherical polar coords
          RAD=ZG(1,1)**(1.0D0/3.0D0)
          ZG(1,1)=RAD
        ELSE IF(ITYP11(1).EQ.4) THEN !prolate-spheroidal coords
          IF(JTYP10.EQ.2) THEN
            SS=ZG(1,1)/(FOCUS*FOCUS)
            SINHZ=DSQRT(SS)
            COSHZ=DSQRT(1.0D0+SS)
          ELSE IF(JTYP10.EQ.3) THEN
            CSS=ZG(1,1)/FOCUS**3
            DES=CSS*CSS-4.0D0/27.0D0
            IF(DES.GT.0.0D0) THEN
              D=((CSS+DSQRT(DES))/2.0D0)**(1.0D0/3.0D0)
              COSHZ=D+1.0d0/(3.0d0*D)
            ELSE
              THETA=DACOS(CSS*DSQRT(27.0D0)/2.0D0)
              COSHZ=2.0D0/DSQRT(3.0D0)*DCOS(THETA/3.0D0)
            ENDIF
            SINHZ=DSQRT(COSHZ*COSHZ-1.0D0)
          ENDIF
          ZG(1,1)=DLOG(COSHZ+SINHZ)
        ENDIF !coord system

        DO ni=1,NIT(nb)
          IF(ITYP11(1).EQ.2) THEN
            ZG(1,NU1(ni))=ZG(1,NU1(ni))/(2.0D0*RAD)
          ELSE IF(ITYP11(1).EQ.3) THEN
            ZG(1,NU1(ni))=ZG(1,NU1(ni))/(3.0D0*RAD*RAD)
          ELSE IF(ITYP11(1).EQ.4) THEN
            IF(JTYP10.EQ.2) THEN
              ZG(1,NU1(ni))=ZG(1,NU1(ni))/
     '          (2.0D0*FOCUS*FOCUS*SINHZ*COSHZ)
            ELSE IF(JTYP10.EQ.3) THEN
              ZG(1,NU1(ni))=ZG(1,NU1(ni))/((3.0D0*COSHZ*COSHZ-1.0D0)
     '          *SINHZ)/FOCUS**3
            ENDIF
          ENDIF
        ENDDO !ni

        IF(JP.GT.1) THEN ! need second derivs
          DO ni=1,NIT(nb)
            DO mi=1,ni
              nu=NU2(ni,mi)
              IF(ITYP11(1).EQ.2) THEN
                ZG(1,nu)=-ZG(1,nu)/(4.0D0*RAD**3)
              ELSE IF(ITYP11(1).EQ.3) THEN
                ZG(1,nu)=-ZG(1,nu)*2.0D0/(9.0D0*RAD**5)
              ELSE IF(ITYP11(1).EQ.4) THEN
                IF(JTYP10.EQ.2) THEN
                  ZG(1,nu)=-ZG(1,nu)*(COSHZ*COSHZ+SINHZ*SINHZ)
     '              /(4.0D0*FOCUS**4*SINHZ**3*COSHZ**3)
                ELSE IF(JTYP10.EQ.3) THEN
                  ZG(1,nu)=-ZG(1,nu)*(2.0D0+9.0D0*SINHZ*SINHZ)*COSHZ/
     '              (FOCUS**6*(3.0D0*COSHZ*COSHZ-1.0D0)**3*SINHZ**3)
                ENDIF
              ENDIF
            ENDDO !mi
          ENDDO !ni
        ENDIF !JP
      ENDIF !isochoric elements

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO nhx=1,NHE
          nh=NH_LOC(nhx,nx)
          WRITE(OP_STRING,'('' ZG(nhx='',I1,'',nu=1..): '',6D12.4,'
     '      //'/(19X,6D12.4))') nhx,(ZG(nhx,nu),nu=1,NUT(NBH(nh)))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZEZG')
      RETURN
 9999 CALL ERRORS('ZEZG',ERROR)
      CALL EXITS('ZEZG')
      RETURN 1
      END


