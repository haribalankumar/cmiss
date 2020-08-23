      SUBROUTINE XFXG(NBJF,NEAXT,NIEF,ng,nr,PG,XDF,XGF,ERROR,*)

C#### Subroutine: XFXG
C###  Description:
C###    XFXG evaluates Gauss point array XGF from face node array
C###    XDF at current Gauss point ng. XGF is corrected for JTYP10=2
C###    (isochoric interpolation).

C#### Variable: XGF(nj,nu,neax)
C###  Type: REAL*8
C###  Set_up: XFXG
C###  Description:
C###    <HTML>
C###    <P> Gauss point derivative nu of geometric variable nj at a face
C###    calculated from adjacent element neax. </P>

C###    <PRE>
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
!     Parameter List
      INTEGER NBJF(NJM),NEAXT,NIEF(0:2),ng,nr
      REAL*8 PG(NSM,NUM,NGM,NBM),XDF(NSFM,2,2,NJM),XGF(NJM,NUM,2)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,neax,ni_f,nj,njj,njj1,njj2,NSTNAT,nu,NU1(3)
      REAL*8 COSHX,CSS,D,DDOT,DES,RAD,SINHX,SS,THETA

      DATA NU1/2,4,7/
C      DATA NU2/3,6,9,6,5,10,9,10,8/

      CALL ENTERS('XFXG',*9999)

C     Interpolate geometry and its first derivs
      DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
        nj=NJ_LOC(NJL_GEOM,njj,nr)
        nb=NBJF(nj)
        NSTNAT=NST(nb)+NAT(nb)
C       Gauss point positions (needed if not rectangular Cartesian)
        IF(ITYP10(nr).NE.1) THEN
          XGF(nj,1,1)=DDOT(NSTNAT,PG(1,1,ng,nb),1,XDF(1,1,1,nj),1)
        ENDIF
C       first derivatives wrt out of face xi
        DO neax=1,NEAXT
          XGF(nj,NU1(NIEF(0)),neax)=
     '      DDOT(NSTNAT,PG(1,1,ng,nb),1,XDF(1,2,neax,nj),1)
        ENDDO !neax
C       first derivatives wrt in face xi
        DO ni_f=1,NIT(nb)
          XGF(nj,NU1(NIEF(ni_f)),1)=
     '      DDOT(NSTNAT,PG(1,NU1(ni_f),ng,nb),1,XDF(1,1,1,nj),1)
        ENDDO !ni_f
        IF(NEAXT.EQ.2) THEN
          IF(ITYP10(nr).NE.1) THEN
            XGF(nj,1,2)=XGF(nj,1,1)
          ENDIF
          DO ni_f=1,NIT(nb)
            XGF(nj,NU1(NIEF(ni_f)),2)=XGF(nj,NU1(NIEF(ni_f)),1)
          ENDDO !ni_f
        ENDIF
      ENDDO !njj

C     Interpolate fibre angles
      DO njj=1,NJ_LOC(NJL_FIBR,0,nr)
        nj=NJ_LOC(NJL_FIBR,njj,nr)
        nb=NBJF(nj)
        DO neax=1,NEAXT
          XGF(nj,1,neax)=
     '      DDOT(NST(nb)+NAT(nb),PG(1,1,ng,nb),1,XDF(1,1,neax,nj),1)
        ENDDO !neax
      ENDDO !njj


C     Isochoric interpolation
      IF(JTYP10.GE.2) THEN
        nb=NBJF(1)

        DO neax=1,NEAXT
          IF(ITYP10(nr).EQ.2) THEN !cyl polar
            RAD=DSQRT(XGF(1,1,neax))
            XGF(1,1,neax)=RAD
          ELSE IF(ITYP10(nr).EQ.3) THEN !sph polar
            RAD=XGF(1,1,neax)**(1.0D0/3.0D0)
            XGF(1,1,neax)=RAD
          ELSE IF(ITYP10(nr).EQ.4) THEN !prolate sph
            IF(JTYP10.EQ.2) THEN
              SS=XGF(1,1,neax)/(FOCUS*FOCUS)
              SINHX=DSQRT(SS)
              COSHX=DSQRT(1.0D0+SS)
            ELSE IF(JTYP10.EQ.3) THEN
              CSS=XGF(1,1,neax)/FOCUS**3
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
            XGF(1,1,neax)=DLOG(COSHX+SINHX)
          ELSE IF(ITYP10(nr).EQ.5) THEN !oblate sph
          ENDIF

          DO ni_f=1,NIT(nb)
            nu=NU1(ni_f)
            IF(ITYP10(nr).EQ.2) THEN !cyl polar
              XGF(1,nu,neax)=XGF(1,nu,neax)/(2.0D0*RAD)
            ELSE IF(ITYP10(nr).EQ.3) THEN !sph polar
              XGF(1,nu,neax)=XGF(1,nu,neax)/(3.0D0*RAD*RAD)
            ELSE IF(ITYP10(nr).EQ.4) THEN !prolate sph
              IF(JTYP10.EQ.2) THEN
                XGF(1,nu,neax)=XGF(1,nu,neax)/
     '            (2.0D0*FOCUS*FOCUS*SINHX*COSHX)
              ELSE IF(JTYP10.EQ.3) THEN
                XGF(1,nu,neax)=XGF(1,nu,neax)/
     '            ((3.0D0*COSHX*COSHX-1.0D0)*SINHX)/FOCUS**3
              ENDIF
            ELSE IF(ITYP10(nr).EQ.5) THEN !oblate sph
            ENDIF
          ENDDO !ni_f

C       Second derivatives not used.
C        DO ni=1,NIT(nb)
C          DO mi=1,ni
C            nu=NU2(ni,mi)
C            IF(ITYP10(nr).EQ.2) THEN      !cyl polar
C              XG(1,nu)=-XG(1,nu)/(4.0D0*RAD**3)
C            ELSE IF(ITYP10(nr).EQ.3) THEN !sph polar
C              XG(1,nu)=-XG(1,nu)*2.0D0/(9.0D0*RAD**5)
C            ELSE IF(ITYP10(nr).EQ.4) THEN !prolate sph
C              IF(JTYP10.EQ.2) THEN
C                XG(1,nu)=-XG(1,nu)*(COSHX*COSHX+SINHX*SINHX)
C     '            /(4.0d0*FOCUS**4*SINHX**3*COSHX**3)
C              ELSE IF(JTYP10.EQ.3) THEN
C                XG(1,nu)=-XG(1,nu)*(2.0D0+9.0D0*SINHX*SINHX)*COSHX/
C     '            (FOCUS**6*(3.0D0*COSHX*COSHX-1)**3*SINHX**3)
C              ENDIF
C            ELSE IF(ITYP10(nr).EQ.5) THEN !oblate sph
C            ENDIF
C          ENDDO !mi
C        ENDDO !ni
        ENDDO !neax
      ENDIF !jtyp10

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO njj1=1,3
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            WRITE(OP_STRING,'('' XGF(nj='',I1,'',nu=1..,neax='','//
     '        'I1,''): '',4(1X,D12.4),/(18X,4(1X,D12.4)))')
     '        nj,neax,(XGF(nj,nu,neax),nu=1,NUT(NBJF(nj)))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO !njj2
        ENDDO !njj1
CC$      call mp_unsetlock()
      ENDIF !dop

      CALL EXITS('XFXG')
      RETURN
 9999 CALL ERRORS('XFXG',ERROR)
      CALL EXITS('XFXG')
      RETURN 1
      END


