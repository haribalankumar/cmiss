      SUBROUTINE XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*)

C#### Subroutine: XEXG
C###  Description:
C###    XEXG evaluates Gauss point array XG from element node array
C###    XE at current Gauss point ng. XG is corrected for JTYP10=2
C###    (isochoric interpolation).

C#### Variable: XG(nj,nu)
C###  Type: REAL*8
C###  Set_up: XEXG
C###  Description:
C###    <HTML> <PRE>
C###    XG(nj, 1) is Xj
C###    XG(nj, 2) is d(Xj)/d(Xi1)         !1st deriv wrt Xi1
C###    XG(nj, 3) is d2(Xj)/d(Xi1)2
C###    XG(nj, 4) is d(Xj)/d(Xi2)         !1st deriv wrt Xi2
C###    XG(nj, 5) is d2(Xj)/d(Xi2)2
C###    XG(nj, 6) is d2(Xj)/d(Xi1)d(Xi2)
C###    XG(nj, 7) is d(Xj)/d(Xi3)         !1st deriv wrt Xi3
C###    XG(nj, 8) is d2(Xj)/d(Xi3)2
C###    XG(nj, 9) is d2(Xj)/d(Xi1)d(Xi3)
C###    XG(nj,10) is d2(Xj)/d(Xi2)d(Xi3)
C###    </PRE> </HTML>

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
C$    INCLUDE 'mp00.cmn'
!     Parameter List
      INTEGER NBJ(NJM),ng,nr
      REAL*8 PG(NSM,NUM,NGM,NBM),XE(NSM,NJM),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,na,nb,ni,nj,njj1,njj2,nk,NSTNAT,nu,NU2(3,3)
      REAL*8 COSHX,CSS,D,DDOT,DES,RAD,SINHX,SS,SUM,SUMM,THETA
      DATA NU2/3,6,9,6,5,10,9,10,8/

      CALL ENTERS('XEXG',*9999)

      DO njj1=1,3,2   !Loop over geometry and field
        DO njj2=1,NJ_LOC(njj1,0,nr)
          nj=NJ_LOC(njj1,njj2,nr)
          nb=NBJ(nj)
          IF(nb.GT.0) THEN
            NSTNAT=NST(nb)+NAT(nb)
            DO nu=1,NUT(nb)
              IF(NNT(nb).GT.0) THEN
                XG(nj,nu)=DDOT(NSTNAT,PG(1,nu,ng,nb),1,XE(1,nj),1)
              ELSE
                SUM=0.0D0
                IF(NAT(nb).GT.0) THEN
                  DO nk=1,NKT(0,nb)
                    SUMM=0.0D0
                    DO na=1,NAT(nb)
                      SUMM=SUMM+XE(na,nj)
                    ENDDO
                    SUM=SUM+PG(nk,nu,ng,nb)*SUMM
                  ENDDO !nk
                ENDIF
                XG(nj,nu)=SUM
              ENDIF

CC$    IF(mp_my_threadnum().EQ.1) THEN
CC$    call mp_setlock()
C      WRITE(OP_STRING,'('' PG(nk,nu,ng,nb)'',6D12.4)') PG(nk,nu,ng,nb)
C      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$    call mp_unsetlock()
CC$    ENDIF


            ENDDO !nu
          ENDIF !nb>0
        ENDDO !njj2
      ENDDO !njj1

C     Fibre angle interpolation
      IF(NJ_LOC(NJL_FIBR,0,nr).GT.0) THEN !fibre angle defined
        DO njj1=1,NJ_LOC(NJL_FIBR,0,nr)
          nj=NJ_LOC(NJL_FIBR,njj1,nr)
          nb=NBJ(nj)
          NSTNAT=NST(nb)+NAT(nb)
          IF(nb.GT.0) THEN
            DO nu=1,NUT(nb)
              IF(nb.GT.0.AND.NNT(nb).GT.0) THEN
                XG(nj,nu)=DDOT(NSTNAT,PG(1,nu,ng,nb),1,XE(1,nj),1)
              ELSE
                XG(nj,nu)=0.0d0
              ENDIF
            ENDDO !nu
          ENDIF !nb
        ENDDO !njj1
      ENDIF !fibres

C     Isochoric interpolation
      IF(JTYP10.GE.2) THEN
        nb=NBJ(1)

        IF(ITYP10(nr).EQ.2) THEN      !cyl polar
          RAD=DSQRT(XG(1,1))
          XG(1,1)=RAD
        ELSE IF(ITYP10(nr).EQ.3) THEN !sph polar
          RAD=XG(1,1)**(1.0D0/3.0D0)
          XG(1,1)=RAD
        ELSE IF(ITYP10(nr).EQ.4) THEN !prolate sph
          IF(JTYP10.EQ.2) THEN
            SS=XG(1,1)/(FOCUS*FOCUS)
            SINHX=DSQRT(SS)
            COSHX=DSQRT(1.0D0+SS)
          ELSE IF(JTYP10.EQ.3) THEN
            CSS=XG(1,1)/FOCUS**3
            DES=CSS*CSS-4.0D0/27.0D0
            IF(DES.GT.0.0D0) THEN
              D=((CSS+DSQRT(DES))/2.0D0)**(1.0D0/3.0D0)
              COSHX=D+1.0D0/(3.0D0*D)
            ELSE
              THETA=DACOS(CSS*DSQRT(27.0D0)/2.0D0)
              COSHX=2.0D0/DSQRT(3.0D0)*DCOS(THETA/3.0D0)
            ENDIF
            SINHX=DSQRT(COSHX*COSHX-1.D0)
          ENDIF
          XG(1,1)=DLOG(COSHX+SINHX)
        ELSE IF(ITYP10(nr).EQ.5) THEN !oblate sph
        ENDIF

        DO ni=1,NIT(nb)
          nu=1+ni*(ni+1)/2
          IF(ITYP10(nr).EQ.2) THEN      !cyl polar
            XG(1,nu)=XG(1,nu)/(2.0D0*RAD)
          ELSE IF(ITYP10(nr).EQ.3) THEN !sph polar
            XG(1,nu)=XG(1,nu)/(3.0D0*RAD*RAD)
          ELSE IF(ITYP10(nr).EQ.4) THEN !prolate sph
            IF(JTYP10.EQ.2) THEN
              XG(1,nu)=XG(1,nu)/(2.0D0*FOCUS*FOCUS*SINHX*COSHX)
            ELSE IF(JTYP10.EQ.3) THEN
              XG(1,nu)=XG(1,nu)/((3.0D0*COSHX*COSHX-1.0D0)*SINHX)
     '                               /FOCUS**3
            ENDIF
          ELSE IF(ITYP10(nr).EQ.5) THEN !oblate sph
          ENDIF
        ENDDO !ni

        DO ni=1,NIT(nb)
          DO mi=1,ni
            nu=NU2(ni,mi)
            IF(ITYP10(nr).EQ.2) THEN      !cyl polar
              XG(1,nu)=-XG(1,nu)/(4.0D0*RAD**3)
            ELSE IF(ITYP10(nr).EQ.3) THEN !sph polar
              XG(1,nu)=-XG(1,nu)*2.0D0/(9.0D0*RAD**5)
            ELSE IF(ITYP10(nr).EQ.4) THEN !prolate sph
              IF(JTYP10.EQ.2) THEN
                XG(1,nu)=-XG(1,nu)*(COSHX*COSHX+SINHX*SINHX)
     '            /(4.0d0*FOCUS**4*SINHX**3*COSHX**3)
              ELSE IF(JTYP10.EQ.3) THEN
                XG(1,nu)=-XG(1,nu)*(2.0D0+9.0D0*SINHX*SINHX)*COSHX/
     '            (FOCUS**6*(3.0D0*COSHX*COSHX-1)**3*SINHX**3)
              ENDIF
            ELSE IF(ITYP10(nr).EQ.5) THEN !oblate sph
            ENDIF
          ENDDO !mi
        ENDDO !ni
      ENDIF !jtyp10

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO njj1=1,3
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            WRITE(OP_STRING,
     '        '('' XG(nj='',I1,'',nu=1..): '',4(1X,D12.4),'
     '        //'/(18X,4(1X,D12.4)))') nj,(XG(nj,nu),nu=1,NUT(NBJ(nj)))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO !njj2
        ENDDO !njj1
CC$      call mp_unsetlock()
      ENDIF !dop

      CALL EXITS('XEXG')
      RETURN
 9999 CALL ERRORS('XEXG',ERROR)
      CALL EXITS('XEXG')
      RETURN 1
      END


