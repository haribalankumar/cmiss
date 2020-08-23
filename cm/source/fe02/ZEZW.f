      SUBROUTINE ZEZW(D2X,IP,IBT,IDO,INP,NAN,NBH,NHTOT,nr,nx,
     '  DXIX,ZE,ZW,XI,ERROR,*)

C#### Subroutine: ZEZW
C###  Description:
C###    ZEZW interpolates ZE(nhx,1) and first derivs ZE(nhx,NU1(ni))
C###    at XI for Lagrange/Hermite tensor product basis functions with
C###    auxiliary element parameters and returns ZW(nhx,NU1(ni)). ZW
C###    is corrected for JTYP10=2. Derivatives are wrt Xi if IP=0 or
C###    wrt nu if IP=1.

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
      INTEGER D2X,IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IP,
     '  NAN(NIM,NAM,NBFM),NBH(NHM),NHTOT,nr,nx
      REAL*8 DXIX(3,3),XI(3),ZE(NSM,NHM),ZW(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER j,k,mi,nb,nh,nhx,ni,NITB,nu,NU1(0:3),NU2(3,3)
      REAL*8 COSHZ,CSS,D,DES,DZHX(3),PFXI,RAD,SINHZ,SS,THETA

      DATA NU1/1,2,4,7/
      DATA NU2/3,6,9,6,5,10,9,10,8/

      CALL ENTERS('ZEZW',*9999)

      NITB=NIT(NBH(NH_LOC(1,nx)))
      DO ni=0,NITB
        DO nhx=1,NHTOT
          nh=NH_LOC(nhx,nx)
          nb=NBH(nh)
!new      MPN 19/7/93 array overrun for NNT in PFXI if nb=0
          IF(nb.GT.0) THEN
            ZW(nhx,NU1(ni))=PFXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        NAN(1,1,nb),nb,NU1(ni),ZE(1,nhx),XI)
          ELSE
            ZW(nhx,NU1(ni))=0.0d0 !should never get into here!!!
            WRITE(OP_STRING,'('' >>Warning: nb is zero in ZEZW'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO !nhx
      ENDDO !ni

      IF(D2X.EQ.1) THEN !calculate second derivatives
        DO j=1,NITB
          DO k=1,NITB
            DO nhx=1,NHTOT
              nh=NH_LOC(nhx,nx)
              nb=NBH(nh)
              IF(nb.GT.0) THEN
                ZW(nhx,NU2(j,k))=PFXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '            INP(1,1,nb),NAN(1,1,nb),nb,NU2(j,k),ZE(1,nhx),XI)
              ELSE
                ZW(nhx,NU2(j,k))=0.0d0 !should never get into here!!!
                WRITE(OP_STRING,'('' >>Warning: nb is zero in ZEZW'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF(IP.EQ.1) THEN
        IF(D2X.EQ.1) THEN
          WRITE(OP_STRING,'('' >>Warning: second derivatives not'
     '      //' supported'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        DO nhx=1,NHTOT
          nh=NH_LOC(nhx,nx)
          DO mi=1,NITB
            DZHX(mi)=0.0D0
            DO ni=1,NITB
              DZHX(mi)=DZHX(mi)+ZW(nhx,NU1(ni))*DXIX(ni,mi)
            ENDDO !ni
          ENDDO !mi
          DO mi=1,NITB
            ZW(nhx,NU1(mi))=DZHX(mi)
          ENDDO !mi
        ENDDO !nhx
      ENDIF

      IF(JTYP10.GE.2) THEN
        IF(D2X.EQ.1) THEN
          WRITE(OP_STRING,'('' >>Warning: second derivatives not'
     '      //' supported'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        nb=NBH(NH_LOC(1,nx))
        IF(ITYP10(nr).EQ.2) THEN
          RAD=DSQRT(ZW(1,1))
          ZW(1,1)=RAD
        ELSE IF(ITYP10(nr).EQ.3) THEN
          RAD=ZW(1,1)**(1.0D0/3.0D0)
          ZW(1,1)=RAD
        ELSE IF(ITYP10(nr).EQ.4) THEN
          IF(JTYP10.EQ.2) THEN
            SS=ZW(1,1)/(FOCUS*FOCUS)
            SINHZ=DSQRT(SS)
            COSHZ=DSQRT(1.0D0+SS)
          ELSE IF(JTYP10.EQ.3) THEN
            CSS=ZW(1,1)/FOCUS**3
            DES=CSS*CSS-4.0D0/27.0D0
            IF(DES.GT.0.0D0) THEN
              D=((CSS+DSQRT(DES))/2.0D0)**(1.0D0/3.0D0)
              COSHZ=D+1.0D0/(3.0D0*D)
            ELSE
              THETA=DACOS(CSS*DSQRT(27.0D0)/2.0D0)
              COSHZ=2.0D0/DSQRT(3.0D0)*DCOS(THETA/3.0D0)
            ENDIF
            SINHZ=DSQRT(COSHZ*COSHZ-1.0D0)
          ENDIF
          ZW(1,1)=DLOG(COSHZ+SINHZ)
        ENDIF
        DO ni=1,NIT(nb)
          IF(ITYP10(nr).EQ.2) THEN
            ZW(1,NU1(ni))=ZW(1,NU1(ni))/(2.0D0*RAD)
          ELSE IF(ITYP10(nr).EQ.3) THEN
            ZW(1,NU1(ni))=ZW(1,NU1(ni))/(3.0D0*RAD*RAD)
          ELSE IF(ITYP10(nr).EQ.4) THEN
            IF(JTYP10.EQ.2) THEN
              ZW(1,NU1(ni))=ZW(1,NU1(ni))/
     '          (2.0D0*FOCUS*FOCUS*SINHZ*COSHZ)
            ELSE IF(JTYP10.EQ.3) THEN
              ZW(1,NU1(ni))=ZW(1,NU1(ni))/((3.0D0*COSHZ*COSHZ-1.0D0)*
     '          SINHZ)/FOCUS**3
            ENDIF
          ENDIF
        ENDDO !ni
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO nhx=1,NHTOT
          nh=NH_LOC(nhx,nx)
          nb=NBH(nh)
          WRITE(OP_STRING,'('' ZW(nhx='',I1,'',nu=1..):'',6D12.4,'
     '      //'/(18X,6D12.4))') nhx,(ZW(nhx,nu),nu=1,NUT(nb))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nhx
CC$      call mp_unsetlock()
      ENDIF !dop

      CALL EXITS('ZEZW')
      RETURN
 9999 CALL ERRORS('ZEZW',ERROR)
      CALL EXITS('ZEZW')
      RETURN 1
      END


