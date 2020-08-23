      SUBROUTINE XEXW(D2X,IBT,IDO,INP,NAN,NBJ,nr,XE,XW,XI,ERROR,*)

C#### Subroutine: XEXW
C###  Description:
C###    XEXW interpolates XE(nj,1) and first derivatives XE(nj,NU1(ni))
C###    at XI for Lagrange/Hermite tensor product basis functions
C###    and returns XW(nj,NU1(ni)). XW is corrected for JTYP10=2
C###    (isochoric interpolation).

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
     '  INP(NNM,NIM,NBFM),NAN(NIM,NAM,NBFM),NBJ(NJM),nr
C JPP 3Apr2003 Parameterize XI with NIM
C     REAL*8 XI(3),XE(NSM,NJM),XW(NJM,NUM)
      REAL*8 XI(NIM),XE(NSM,NJM),XW(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER j,k,nb,nj,njj1,njj2,ni,NITB,nu,NU1(0:3),NU2(3,3)
      REAL*8 COSHZ,CSS,D,DES,PFXI,RAD,SINHZ,SS,THETA

      DATA NU1/1,2,4,7/
      DATA NU2/3,6,9,6,5,10,9,10,8/

      CALL ENTERS('XEXW',*9999)

      NITB=NIT(NBJ(1))
      DO ni=0,NITB
        DO njj1=1,3 !geometry/fibres/field
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            nb=NBJ(nj)
!new        MPN 19/7/93 array overrun for NNT in PFXI if nb=0
            IF(nb.GT.0) THEN
              XW(nj,NU1(ni))=PFXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '          NAN(1,1,nb),nb,NU1(ni),XE(1,nj),XI)
            ELSE
              XW(nj,NU1(ni))=0.0d0 !should never get into here!!!
            ENDIF
          ENDDO !njj2
        ENDDO !njj1
      ENDDO !ni

      IF(D2X.EQ.1) THEN !calculate second derivatives
        DO j=1,NITB
          DO k=1,NITB
            DO njj1=1,3 !geometry/fibres/field
              DO njj2=1,NJ_LOC(njj1,0,nr)
                nj=NJ_LOC(njj1,njj2,nr)
                nb=NBJ(nj)
                IF(nb.GT.0) THEN
                  XW(nj,NU2(j,k))=PFXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '              INP(1,1,nb),NAN(1,1,nb),nb,NU2(j,k),XE(1,nj),XI)
                ELSE
                  XW(nj,NU2(j,k))=0.0d0 !should never get into here
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF(JTYP10.GE.2) THEN
        IF(D2X.EQ.1) THEN
          WRITE(OP_STRING,'('' >>Warning: second derivatives not'
     '      //' updated for isochoric interpolation'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        nb=NBJ(1)
        IF(ITYP10(nr).EQ.2) THEN
          RAD=DSQRT(XW(1,1))
          XW(1,1)=RAD
        ELSE IF(ITYP10(nr).EQ.3) THEN
          RAD=XW(1,1)**(1.0D0/3.0D0)
          XW(1,1)=RAD
        ELSE IF(ITYP10(nr).EQ.4) THEN
          IF(JTYP10.EQ.2) THEN
            SS=XW(1,1)/(FOCUS*FOCUS)
            SINHZ=DSQRT(SS)
            COSHZ=DSQRT(1.0D0+SS)
          ELSE IF(JTYP10.EQ.3) THEN
            CSS=XW(1,1)/FOCUS**3
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
          XW(1,1)=DLOG(COSHZ+SINHZ)
        ENDIF
        DO ni=1,NIT(nb)
          IF(ITYP10(nr).EQ.2) THEN
            XW(1,NU1(ni))=XW(1,NU1(ni))/(2.0D0*RAD)
          ELSE IF(ITYP10(nr).EQ.3) THEN
            XW(1,NU1(ni))=XW(1,NU1(ni))/(3.0D0*RAD*RAD)
          ELSE IF(ITYP10(nr).EQ.4) THEN
            IF(JTYP10.EQ.2) THEN
              XW(1,NU1(ni))=XW(1,NU1(ni))/
     '          (2.0D0*FOCUS*FOCUS*SINHZ*COSHZ)
            ELSE IF(JTYP10.EQ.3) THEN
              XW(1,NU1(ni))=XW(1,NU1(ni))/((3.0D0*COSHZ*COSHZ-1.0D0)*
     '          SINHZ)/FOCUS**3
            ENDIF
          ENDIF
        ENDDO !ni
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO njj1=1,3
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            WRITE(OP_STRING,'('' XW(nj='',I1,'',nu=1..): '',6D12.4,'
     '        //'/(18X,6D12.4))') nj,(XW(nj,nu),nu=1,NUT(NBJ(nj)))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO !njj2
        ENDDO !njj1
CC$      call mp_unsetlock()
      ENDIF !dop

      CALL EXITS('XEXW')
      RETURN
 9999 CALL ERRORS('XEXW',ERROR)
      CALL EXITS('XEXW')
      RETURN 1
      END


