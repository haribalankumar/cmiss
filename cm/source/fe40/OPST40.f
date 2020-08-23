      SUBROUTINE OPST40(nb,NBH,NBJ,ne,NHE,nw,nx,
     '  CG,ENERGY_ELEMENT,PG,TYPE,WG,XE,XG,YG,ZE,ZG,
     '  PRINCIPAL,COORDS,ERROR,*)

C#### Subroutine: OPST40
C###  Description:
C###    If TYPE='WRITES':o/p stress,strain and SE at Gauss points.
C###    IF TYPE='STRESS':puts 1st principal stresses in YG(1,ng,ne);
C###      2nd principal stresses in YG(2,ng,ne);
C###      principal angles in YG(3,ng,ne); and
C###      strain energy in YG(4,ng,ne)
C###    If growth law is defined (KTYP60=1 or 2), material density in
C###      YG(5,ng,ne)

C**** Note: ENERGY_GAUSS is strain energy at Gauss point
C****       ENERGY_ELEMENT is total for current element
C**** XG are undeformed Xj coords and derivs wrt Xi
C**** ZG are displacements and derivs wrt Nu coords

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbst02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp60.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'press00.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER nb,NBH(NHM),NBJ(NJM),ne,NHE,nw,nx
      REAL*8 CG(NMM,NGM),ENERGY_ELEMENT,PG(NSM,NUM,NGM,NBM),
     '  XE(NSM,NJM),WG(NGM,NBM),XG(NJM,NUM),
     '  YG(NIYGM,NGM),ZE(NSM,NHM),ZG(NHM,NUM)
      LOGICAL PRINCIPAL
      CHARACTER COORDS*(*),ERROR*(*),TYPE*(*)
!     Local Variables
      INTEGER i,IFAIL,j,k,l,mi,ng,nh,nhx,ni,nj,njj1,njj2,
     '  nr,ns,NU1(0:3)
      REAL*8 AXL(3,3),AXIAL_STRAIN,AZ,AZL(3,3),
     '  DS0,DS1,DSDX1,DX1,DX2,
     '  DXINU(3,3),DXIX(3,3),EAU,EAUV,EG(3,3),EGNU(3,3),
     '  EIV,EIV1,EIV2,
     '  ENERGY_GAUSS,ETA1,ETA2,EVAL(3),EVEC(3),EVECT(3,3),
     '  GL(3,3),GU(3,3),
     '  PGX,PHI,PPG(9,2,16),PV,RG,RGX,RGZ,RL,RT(3,3),SZ,
     '  TENSION,TG(3,3),TGNU(3,3),WK1_LOCAL(10)
      CHARACTER FORMAT*132
      DATA NU1/1,2,4,7/

      CALL ENTERS('OPST40',*9999)
C GBS 10-NOV-1994
      nr=1   !Temporary
      ENERGY_ELEMENT=0.0d0

!Cable or truss elements
      IF(nw.EQ.1) THEN
        IF(TYPE(1:6).EQ.'WRITES') THEN
          WRITE(OP_STRING,
     '      '(/'' Element '',I3,'' (Cable or truss):'')') ne
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

        !Calculate strain & tension (in kN)
        CALL ZETG41(CG(1,1),AXIAL_STRAIN,SZ,TENSION,
     '    XE,ZE,ERROR,*9999)
        ENERGY_GAUSS=0.5d0*TENSION*AXIAL_STRAIN
        ENERGY_ELEMENT=ENERGY_ELEMENT+ENERGY_GAUSS
        IF(DABS(TENSION).LT.PRSTMIN) PRSTMIN=DABS(TENSION)
        IF(DABS(TENSION).GT.PRSTMAX) PRSTMAX=DABS(TENSION)
        IF(TYPE(1:6).EQ.'WRITES') THEN
          WRITE(OP_STRING,
     '      '('' Strain= '',E11.3,'' Tension = '',E11.3,'' kN'')')
     '      AXIAL_STRAIN,TENSION
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(TYPE(1:6).EQ.'STRESS') THEN
          DO ng=1,NGT(NBH(NH_LOC(1,nx)))
c           YG(1,ng)=TENSION !stores 1st principal stress
            YG(1,ng)=0.0d0 !stores 1st principal stress
            YG(2,ng)=0.0d0 !stores 2nd principal stress
            YG(3,ng)=0.0d0 !stores principal angle
            YG(4,ng)=ENERGY_GAUSS !stores strain energy
          ENDDO
        ENDIF

!Batten elements
      ELSE IF(nw.EQ.2) THEN
        IF(TYPE(1:6).EQ.'WRITES') THEN
          WRITE(OP_STRING,'(/'' Element '',I3,'' (Batten):''/)') ne
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        DO ng=1,NGT(NBH(NH_LOC(1,nx)))
          CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
          CALL XGMG(0,1,NBJ(1),nr,DXIX,GL,GU,RG,XG,ERROR,*9999)
          CALL ZETG42(NBH,ng,CG(1,ng),DS0,DS1,DSDX1,EIV1,EIV2,PG,
     '      PV,RL,XE,ZE,ERROR,*9999)
          IF(TYPE(1:6).EQ.'WRITES') THEN
            WRITE(OP_STRING,'('' Gauss point '',I2,'' EIV1='',E11.4,'
     '        //''' kNm  EIV2='',E11.4,'' kN  PV='',E11.4,'' kN'')')
     '        ng,EIV1,EIV2,PV
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO

!Beam elements
      ELSE IF(nw.EQ.3) THEN
        IF(TYPE(1:6).EQ.'WRITES') THEN
          WRITE(OP_STRING,'(/'' Element '',I3,'' (Beam):''/)') ne
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        DO ng=1,NGT(NBH(NH_LOC(1,nx)))
          CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
          CALL XGMG(0,1,NBJ(1),nr,DXIX,GL,GU,RG,XG,ERROR,*9999)
          DX1=DXIX(1,1)
          DX2=DX1*DX1

          !Calculate axial force (in kN) & bending moment (in kNm)
          CALL ZETG43(NBH,ng,AXIAL_STRAIN,CG(1,ng),DX2,EAU,EAUV,EIV,
     '      PG,XE,ZE,ERROR,*9999)
          IF(TYPE(1:6).EQ.'WRITES') THEN
            WRITE(OP_STRING,
     '        '('' Gauss point ''   ,I2,''  Axial strain  =''  ,E11.4,'
     '          //'/16X,'' Axial force   =''   ,E11.4,'' kN'','
     '          //'/16X,'' Bending moment='',E11.4,'' kNm'')')
     '          NG,AXIAL_STRAIN,EAU,EIV
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO

!Membrane elements
      ELSE IF(nw.EQ.5) THEN
        IF(TYPE(1:6).EQ.'WRITES') THEN
          WRITE(OP_STRING,'(/'' Element '',I3,'' (Membrane):'')') ne
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        DO ng=1,NGT(NBH(Nh_LOC(1,nx)))

          !Calculate DXINU: derivs of Xi wrt Nu (fibre) coords
          !(unless no fibres are defined; then stresses wrt Xi coords)
          CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
          CALL XGMG(1,0,NBJ(1),nr,DXINU,GL,GU,RGX,XG,ERROR,*9999)
          IF(JTYP12.LE.1) THEN
            ETA1=XG(NJ_LOC(NJL_FIBR,1,nr),1)
            IF(TYPE(1:6).EQ.'WRITES') THEN
              WRITE(OP_STRING,
     '          '(/'' Gauss point '',I2,'' Fibre angle='',F6.1,'
     '          //''' degrees wrt Xi(1) coord'')') ng,ETA1*180.0d0/PI
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE IF(JTYP12.EQ.2) THEN
            ETA2=XG(NJ_LOC(NJL_FIBR,1,nr),1)
            IF(TYPE(1:6).EQ.'WRITES') THEN
              WRITE(OP_STRING,'(/'' Gauss point '',I2,'
     '          //''' Fibre angle='',F6.1,'
     '          //''' degrees wrt Xi(2) coord'')') ng,ETA2*180.0d0/PI
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
          IF(TYPE(1:6).EQ.'WRITES') THEN
            WRITE(OP_STRING,
     '        '('' Pressure = '',E11.3,'' kPa'')') PRESS(ng,ne)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          DO nb=1,NBFT
            IF(nb.EQ.NBH(NH_LOC(1,nx)).OR.nb.EQ.NBH(NH_LOC(2,nx))
     '        .OR.nb.EQ.NBH(NH_LOC(3,nx))) THEN
              DO ns=1,NST(nb)
                DO ni=1,2
                  PPG(nb,ni,ns)=PGX(nb,ni,ns,DXINU,PG(1,1,ng,nb))
                ENDDO
              ENDDO
            ENDIF
          ENDDO

          !Calculate XG,ZG & metric tensors - all wrt Nu coords
          CALL ZEAZ45(NBJ,ng,nr,AXL,AZ,AZL,PG,PPG,
     '      XE,XG,ZE,ZG,ERROR,*9999)
          IF(NG.EQ.1.AND.DOP) THEN
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              nb=NBJ(nj)
              IF(TYPE(1:6).EQ.'WRITES') THEN
                WRITE(OP_STRING,'('' XE(ns,'',I1,''): '',9E12.4)')
     '            nj,(XE(ns,nj),ns=1,NST(nb))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
            IF(TYPE(1:6).EQ.'WRITES') THEN
              WRITE(OP_STRING,'('' DXINU:'',4E11.3)')
     '          DXINU(1,1),DXINU(1,2),DXINU(2,1),DXINU(2,2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            DO nhx=1,NHE
              nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
              nb=NBH(nh)
              IF(TYPE(1:6).EQ.'WRITES') THEN
                WRITE(OP_STRING,'('' ZE(ns,'',I1,''): '',9E12.4)')
     '            nhx,(ZE(ns,nhx),ns=1,NST(nb))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
            DO nhx=1,NHE
              nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
              nb=NBH(nh)
              IF(TYPE(1:6).EQ.'WRITES') THEN
                WRITE(OP_STRING,
     '            '('' ZG('',I1,'',NU1(ni)): '',8E15.8)')
     '            nh,(ZG(nh,NU1(ni)),ni=0,NIT(nb))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
            IF(TYPE(1:6).EQ.'WRITES') THEN
              WRITE(OP_STRING,
     '          '('' RGX='',E11.3,'' RGZ='',E11.3)') RGX,RGZ
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF

          !Calculate stress & strain tensors TG,EG & e.values & vectors
          CALL AZTG45(nw,AXL,AZL,CG(1,ng),EG,TG,ERROR,*9999)
          IF(TYPE(1:6).EQ.'WRITES') THEN
            WRITE(OP_STRING,
     '        '('' Strain tensor referred to fibre coords:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,
     '        '('' EG(1,1)='',E12.4,'' EG(2,2)='',E12.4,'
     '        //''' EG(1,2)='',E12.4)') EG(1,1),EG(2,2),EG(1,2)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,
     '        '('' Stress tensor referred to fibre coords '
     '        //'(/unit undeformed length):'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,
     '        '('' TG(1,1)='',E12.4,'' TG(2,2)='',E12.4,'
     '        //''' TG(1,2)='',E12.4,'' kN/m'')')
     '        TG(1,1),TG(2,2),TG(1,2)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL EVALUE(2,TG,EVAL,ERROR,*9999)
          CALL EVECTR(2,TG,EVAL(1),EVEC,ERROR,*9999)
          PHI=DATAN2(EVEC(2),EVEC(1))
          IF(TYPE(1:6).EQ.'WRITES') THEN
            WRITE(OP_STRING,
     '        '('' Princ. stresses:'',2E12.4,'' kN/m  PHI = '','
     '        //'F6.1,'' degrees'')') EVAL(1),EVAL(2),PHI*180.0d0/PI
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(EVAL(1).GT.PRSTMAX) PRSTMAX=EVAL(1)
          IF(DABS(EVAL(1)).LT.PRSTMIN) PRSTMIN=DABS(EVAL(1))
          IF(TYPE(1:6).EQ.'STRESS') THEN
            YG(1,ng)=EVAL(1) !stores 1st principal stress
            YG(2,ng)=EVAL(2) !stores 2nd principal stress
            YG(3,ng)=PHI     !stores principal angle
          ENDIF
        ENDDO

!3D stress
      ELSE IF(nw.EQ.9) THEN
        IF(TYPE(1:6).EQ.'WRITES') THEN
          WRITE(OP_STRING,'(/'' Element '',I3,'' (3D):'')') ne
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

        DO ng=1,NGT(NBH(NH_LOC(1,nx)))
          CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
!PJH 19Jul95
!         CALL XGMG(1,NIT(NBJ(1)),NBJ(1),nr,DXINU,GL,GU,RG,XG,
!    '      ERROR,*9999)
! Call XGMG with ip=0 to give DXINU as derivs wrt r.c. ref coords
C GMH 4/1/96 Call with ip=1 so we get stress wrt fibre coords
          CALL XGMG(1,NIT(NBJ(1)),NBJ(1),nr,DXINU,GL,GU,RG,XG,
     '      ERROR,*9999)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            DO ni=1,3
              WRITE(OP_STRING,'('' DXINU(ni='',I1,'',nj):'',3E11.3)')
     '          ni,(DXINU(ni,nj),nj=1,3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
CC$          call mp_unsetlock()
          ENDIF

!PJH 19Jul95
!         DO ni=1,3
!           DO nj=1,3
!             SUM1=0.d0
!             SUM2=0.d0
!             DO k=1,3
!               SUM1=XG(nj,NU1(k))*DXINU(k,ni)
!               SUM2=XG(ni,NU1(k))*DXINU(k,nj)
!             ENDDO !k
!             dX_dNu(ni,nj)=0.5d0*(SUM1+SUM2)
!           ENDDO !nj
!         ENDDO !ni

          !Calc stress & strain tensors TG,EG & e.values & vectors
          CALL ZEZG(1,NBH,ng,NHE,nx,DXINU,PG,ZE,ZG,ERROR,*9999)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              WRITE(OP_STRING,'('' XG('',I1,'',NU1(ni)): '',8E15.8)')
     '          nj,(XG(nj,NU1(ni)),ni=0,NIT(NBJ(nj)))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
            DO nhx=1,NHE
              nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
              WRITE(OP_STRING,'('' ZG('',I1,'',NU1(ni)): '',8E15.8)')
     '          nh,(ZG(nh,NU1(ni)),ni=0,NIT(NBH(nh)))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
CC$          call mp_unsetlock()
          ENDIF !dop

! Christoffels only needed if strains are referred to non-orthogonal
! cartesian coords. For the isotropic case EG and TG are referred to
! the reference r.c. coords and CHTOFF & dX_dNu not used. (PJH 19Jul95)
!          CALL TOFFEL(ITYP10(nr),NBJ(1),nr,CHTOFF,DBM,GU,XG,X3G,
!     '      .FALSE.,ERROR,*9999)
C GMh 5/1/96 Calc stress and strain wrt material coords
          CALL MAT_VEC_NG(3,nr,RT(1,1),RT(1,2),RT(1,3),XG,
     '      ERROR,*9999)
          CALL CGS9(nr,nw,nx,CG(1,ng),EGNU,RT,TGNU,ZG,ERROR,*9999)
          ENERGY_GAUSS=0.d0
          DO i=1,3
            DO j=1,3
              ENERGY_GAUSS=ENERGY_GAUSS+0.5d0*(TGNU(i,j)*EGNU(i,j))
            ENDDO !j
          ENDDO !i
          ! GR shouldn't this be integrated over the element volume?
          ENERGY_ELEMENT=ENERGY_ELEMENT+ENERGY_GAUSS
C GMH 5/1/96 We must now rotate to get the stress tensor in terms of
C            reference coordinates
          DO i=1,3
            DO j=1,3
              IF(NJ_LOC(njl_fibr,0,nr).NE.0) THEN
                TG(i,j)=0.0D0
                EG(i,j)=0.0D0
                DO k=1,3
                  DO l=1,3
                    TG(i,j)=TG(i,j)+RT(k,i)*TGNU(k,l)*RT(l,j)
                    EG(i,j)=EG(i,j)+RT(k,i)*EGNU(k,l)*RT(l,j)
                  ENDDO !l
                ENDDO !k
              ELSE
                TG(i,j)=TGNU(i,j)
                EG(i,j)=EGNU(i,j)
              ENDIF
            ENDDO !j
          ENDDO !i

          IF(TYPE(1:6).EQ.'WRITES') THEN
            FORMAT='(/'' ng='',I2,'' X(j,0):     '',6E11.4)'
            WRITE(OP_STRING,FMT=FORMAT) ng,
     '        ((XG(NJ_LOC(njj1,njj2,nr),1),
     '        njj2=1,NJ_LOC(njj1,0,nr)),njj1=1,2)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            FORMAT='( '' ====='', '' Z(j,0):     '',6E11.4)'
            WRITE(OP_STRING,FMT=FORMAT)
     '        (ZG(NH_LOC(nhx,nx),1),nhx=1,NHE)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            FORMAT='(7X,''Strains and stresses wrt Nu coords:'')'
            WRITE(OP_STRING,FMT=FORMAT)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO mi=1,3
              IF(mi.EQ.1) THEN
                FORMAT='(7X,''EG(1,i):'',E11.4,20X,''TG(1,i):'',E11.4)'
              ELSE IF(mi.EQ.2) THEN
                FORMAT='(10X,  ''2    '',2E11.4,12X,  ''2    '',2E11.4)'
              ELSE IF(mi.EQ.3) THEN
                FORMAT='(10X,  ''3    '',3E11.4,1X,   ''3    '',3E11.4)'
              ENDIF
              WRITE(OP_STRING,FMT=FORMAT)
     '          (EGNU(mi,ni),ni=1,mi),(TGNU(mi,ni),ni=1,mi)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO !mi
          ENDIF !type

          IF(PRINCIPAL) THEN
            IFAIL=0
C            CALL F02ABF(TG,3,3,EVAL,EVECT,3,WK1_LOCAL,IFAIL)
            DO ni=1,3
              DO mi=1,3
                EVECT(ni,mi)=TG(ni,mi)
              ENDDO
            ENDDO
C MLB 19/3/97
C This may not give evectors as accurately as NAG
            CALL DSYEV('V','L',3,EVECT,3,EVAL,WK1_LOCAL,10,IFAIL)
            CALL ASSERT(IFAIL.EQ.0,
     '        'Could not find eivalues and eivectors',ERROR,*9999)
            IF(TYPE(1:6).EQ.'WRITES') THEN
              FORMAT='(7X,''Principal stresses'//
     '          ' and eigenvectors wrt reference coords:'')'
              WRITE(OP_STRING,FMT=FORMAT)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              DO mi=1,3
                IF(mi.EQ.1) THEN
                  FORMAT='(7X,''PST(1):'',D11.4,'' VCT(1,i):'',3D11.4)'
                ELSE IF(mi.EQ.2) THEN
                  FORMAT='(7X,''    2  '',D11.4,''     2    '',3D11.4)'
                ELSE IF(mi.EQ.3) THEN
                  FORMAT='(7X,''    3  '',D11.4,''     3    '',3D11.4)'
                ENDIF
                WRITE(OP_STRING,FMT=FORMAT)
     '            EVAL(mi),(EVECT(mi,ni),ni=1,3)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO
              WRITE(OP_STRING,
     '          '('' Strain energy = '',E12.4,'' kJ/m^2'')')
     '          ENERGY_GAUSS
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
C            CALL EVALUE(3,TG,EVAL,ERROR,*9999)
C            CALL EVECTR(3,TG,EVAL(1),EVEC,ERROR,*9999)
C GMH 10/7/95 Must check that both arguments are not zero
C            IF((EVEC(2).NE.0.0D0).OR.(EVEC(1).NE.0.0D0)) THEN
C              PHI=DATAN2(EVEC(2),EVEC(1))
C            ELSE
C              PHI=0.0D0
C            ENDIF
C            IF(TYPE(1:6).EQ.'WRITES') THEN
C              WRITE(OP_STRING,
C     '          '('' Princ. stresses:'',3E12.4,'' kN/m  PHI = '','
C     '          //'F6.1,'' degrees'')')
C     '          EVAL(1),EVAL(2),EVAL(3),PHI*180.d0/PI
C              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C              WRITE(OP_STRING,
C     '          '('' Strain energy = '',E12.4,'' kJ/m^2'')')
C     '          ENERGY_GAUSS
C              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C            ENDIF !type
            IF(EVAL(3).GT.PRSTMAX) PRSTMAX=EVAL(3)
            IF(EVAL(1).LT.PRSTMIN) PRSTMIN=EVAL(1)
            IF(TYPE(1:6).EQ.'STRESS') THEN
              CALL ASSERT(NJM.GE.4,
     '          'NIYGM must be >=4 to store principal stresses in YG',
     '          ERROR,*9999)
              YG(1,ng)=EVAL(1) !stores 1st principal stress
              YG(2,ng)=EVAL(2) !stores 2nd principal stress
              YG(3,ng)=EVAL(3) !stores 3rd principal stress
              YG(4,ng)=ENERGY_GAUSS !stores strain energy
            ENDIF
          ELSE
            IF(TYPE(1:6).EQ.'STRESS') THEN
              CALL ASSERT(NJM.GE.10,
     '          'NIYGM must be >=10 to store stress tensor in YG',
     '          ERROR,*9999)
              YG(1,ng)=TG(1,1) !stores stress tensor
              YG(2,ng)=TG(2,2)
              YG(3,ng)=TG(1,2)
              YG(4,ng)=TG(1,3)
              YG(5,ng)=TG(2,1)
              YG(6,ng)=TG(2,3)
              YG(7,ng)=TG(3,1)
              YG(8,ng)=TG(3,2)
              YG(9,ng)=TG(3,3)
              YG(10,ng)=ENERGY_GAUSS !stores strain energy
            ENDIF
          ENDIF
        ENDDO
        IF(TYPE(1:6).EQ.'WRITES') THEN
          WRITE(OP_STRING,
     '      '('' Total S.E. for element = '',E12.4,'' kJ/m^2'')')
     '      ENERGY_ELEMENT
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

!Plane stress and plane strain
      ELSE IF(nw.EQ.11.OR.nw.EQ.12) THEN
        IF(TYPE(1:6).EQ.'WRITES') THEN
          IF(nw.EQ.11) THEN
            WRITE(OP_STRING,
     '        '(/'' Element '',I3,'' (Plane stress):'')') ne
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(nw.EQ.12) THEN
            WRITE(OP_STRING,
     '        '(/'' Element '',I3,'' (Plane strain):'')') ne
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

        DO ng=1,NGT(NBH(NH_LOC(1,nx)))
          CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
          CALL XGMG(1,NIT(NBJ(1)),NBJ(1),nr,DXINU,GL,GU,RG,XG,
     '      ERROR,*9999)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            DO ni=1,2
              WRITE(OP_STRING,'('' DXINU(ni='',I1,'',nj):'',3E11.3)')
     '          ni,(DXINU(ni,nj),nj=1,2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
CC$          call mp_unsetlock()
          ENDIF

C cpb 14/5/96 Old way
C          dX_dNu(1,1)=XG(1,NU1(1))*DXINU(1,1)+XG(1,NU1(2))*DXINU(2,1)
C          dX_dNu(2,1)=XG(2,NU1(1))*DXINU(1,1)+XG(2,NU1(2))*DXINU(2,1)
C          dX_dNu(1,2)=XG(1,NU1(1))*DXINU(1,2)+XG(1,NU1(2))*DXINU(2,2)
C          dX_dNu(2,2)=XG(2,NU1(1))*DXINU(1,2)+XG(2,NU1(2))*DXINU(2,2)
C          IF(JTYP12.LE.1) THEN
C            IF(NJ_LOC(NJL_FIBR,1,nr).GT.0) THEN
C              ETA1=XG(NJ_LOC(NJL_FIBR,1,nr),1)
C            ELSE
C              ETA1=0.d0
C            ENDIF
C            IF(TYPE(1:6).EQ.'WRITES') THEN
C              WRITE(OP_STRING,
C     '          '(/'' Gauss point '',I2,'' Fibre angle='',F6.1,'
C     '          //''' degrees wrt Xi(1) coord'')') ng,ETA1*180.0d0/PI
C              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C            ENDIF
C          ELSE IF(JTYP12.EQ.2) THEN
C            IF(NJ_LOC(NJL_FIBR,1,nr).GT.0) THEN
C              ETA2=XG(NJ_LOC(NJL_FIBR,1,nr),1)
C            ELSE
C              ETA2=0.d0
C            ENDIF
C            IF(TYPE(1:6).EQ.'WRITES') THEN
C              WRITE(OP_STRING,
C     '          '(/'' Gauss point '',I2,'' Fibre angle='',F6.1,'
C     '          //''' degrees wrt Xi(2) coord'')') ng,ETA2*180.0d0/PI
C              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C            ENDIF
C          ENDIF

          IF(TYPE(1:6).EQ.'WRITES') THEN
            WRITE(OP_STRING,'(/'' Gauss point '',I2)') ng
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

          IF(KTYP60.GT.0) THEN !Growth law: Output density and Young's modulus
            WRITE(OP_STRING,
     '        '('' Density = '',E12.4,'' Youngs modulus = '','
     '        //'E12.4)') CG(5,ng),CG(1,ng)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

          !Calculate stress & strain tensors TG,EG & e.values & vectors
          CALL ZEZG(1,NBH,ng,NHE,nx,DXINU,PG,ZE,ZG,ERROR,*9999)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              WRITE(OP_STRING,'('' XG('',I1,'',NU1(ni)): '',8E15.8)')
     '          nj,(XG(nj,NU1(ni)),ni=0,NIT(NBJ(nj)))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
            DO nhx=1,NHE
              nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
              WRITE(OP_STRING,'('' ZG('',I1,'',NU1(ni)): '',8E15.8)')
     '          nh,(ZG(nh,NU1(ni)),ni=0,NIT(NBH(nh)))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
CC$          call mp_unsetlock()
          ENDIF
          CALL MAT_VEC_NG(2,nr,RT(1,1),RT(1,2),RT(1,3),XG,
     '      ERROR,*9999)
          CALL CGS11(nr,nw,nx,CG(1,ng),EGNU,RT,TGNU,ZG,ERROR,*9999)
C CPB 13/5/95 Rotate the stress tensor to get back into reference
C coordinates
          DO i=1,2
            DO j=1,2
              IF((NJ_LOC(njl_fibr,0,nr).NE.0).OR.
     '          (COORDS.EQ.'Reference')) THEN
                TG(i,j)=0.0D0
                EG(i,j)=0.0D0
                DO k=1,2
                  DO l=1,2
                    TG(i,j)=TG(i,j)+RT(i,k)*TGNU(k,l)*RT(j,l)
                    EG(i,j)=EG(i,j)+RT(i,k)*EGNU(k,l)*RT(j,l)
                  ENDDO !l
                ENDDO !k
              ELSE
                TG(i,j)=TGNU(i,j)
                EG(i,j)=EGNU(i,j)
              ENDIF
            ENDDO !j
          ENDDO !i
          ENERGY_GAUSS=0.5d0*(TG(1,1)*EG(1,1)+TG(2,2)*EG(2,2)+
     '      TG(1,2)*EG(1,2))
          ENERGY_ELEMENT=ENERGY_ELEMENT+ENERGY_GAUSS*WG(ng,nb)
          IF(TYPE(1:6).EQ.'WRITES') THEN
            WRITE(OP_STRING,
     '        '('' EG(1,1)='',D12.4,'' EG(2,2)='',D12.4,'
     '        //''' EG(1,2)='',D12.4)') EG(1,1),EG(2,2),EG(1,2)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,
     '        '('' TG(1,1)='',D12.4,'' TG(2,2)='',D12.4,'
     '        //''' TG(1,2)='',D12.4,'' kN/m'')')
     '        TG(1,1),TG(2,2),TG(1,2)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL EVALUE(2,TG,EVAL,ERROR,*9999)
          CALL EVECTR(2,TG,EVAL(1),EVEC,ERROR,*9999)
          PHI=DATAN2(EVEC(2),EVEC(1))
          IF(TYPE(1:6).EQ.'WRITES') THEN
            WRITE(OP_STRING,
     '        '('' Princ. stresses:'',2D12.4,'' kN/m  PHI = '','
     '        //'F6.1,'' degrees'')') EVAL(1),EVAL(2),PHI*180.d0/PI
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,
     '        '('' Strain energy = '',D12.4,'' kJ/m^2'')')
     '        ENERGY_GAUSS
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(EVAL(1).GT.PRSTMAX) PRSTMAX=EVAL(1)
          IF(DABS(EVAL(2)).LT.PRSTMIN) PRSTMIN=DABS(EVAL(2))
          IF(TYPE(1:6).EQ.'STRESS') THEN
            IF(PRINCIPAL) THEN
              YG(1,ng)=EVAL(1)      !stores 1st principal stress
              YG(2,ng)=EVAL(2)      !stores 2nd principal stress
              YG(3,ng)=PHI          !stores principal angle
              YG(4,ng)=ENERGY_GAUSS !stores strain energy
            ELSE
              YG(1,ng)=TG(1,1) !stores stress tensor
              YG(2,ng)=TG(2,2)
              YG(3,ng)=TG(1,2)
              YG(4,ng)=ENERGY_GAUSS !stores strain energy
            ENDIF
          ENDIF
        ENDDO
        IF(TYPE(1:6).EQ.'WRITES') THEN
          WRITE(OP_STRING,
     '      '('' Total S.E. for element = '',D12.4,'' kJ/m^2'')')
     '      ENERGY_ELEMENT
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF !nw

      CALL EXITS('OPST40')
      RETURN
 9999 CALL ERRORS('OPST40',ERROR)
      CALL EXITS('OPST40')
      RETURN 1
      END


