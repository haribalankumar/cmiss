      SUBROUTINE ZERE55(INP,NBH,ne,NHE,nr,nx,
     '  CG,PG,RE,WG,ZE,ZEREF,ZG,ERROR,*)

C#### Subroutine: ZERE55
C###  Description:
C###    ZERE55 calculates constant-volume elements for cavity loading.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER INP(NNM,NIM,NBFM),NBH(NHM),ne,NHE,nr,nx
      REAL*8 CG(NMM,NGM),PG(NSM,NUM,NGM,NBM),RE(NSM,NHM),
     '  WG(NGM,NBM),ZE(NSM,NHM),ZEREF(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,IXI1,IXI2,IXI3,nb,ng,nh,nhx,nk,nn,ns
      REAL*8 deltaMU,DXIX(3,3),GZ,GZL(3,3),GZU(3,3),
     '  kstif,press_current,press_init,
     '  press_sign,RGZ,RGZREF,RWGZ,RWGZREF,VOLDEF,VOLUND
      LOGICAL FOUNDnn

      CALL ENTERS('ZERE55',*9999)

      VOLDEF=0.0d0  !  Deformed element volume
      VOLUND=0.0d0  !Undeformed element volume

      DO i=1,3
        DO j=1,3
          DXIX(i,j)=0.0d0
        ENDDO
      ENDDO

      DO nhx=1,NH_LOC(0,nx)
        nh=NH_LOC(nhx,nx)
        DO ns=1,NST(NBH(nh))+NAT(NBH(nh))
          RE(ns,nh)=0.0d0
        ENDDO !ns
      ENDDO !nh

      CALL ASSERT(KTYP52(nr).EQ.2,'Must be incompressible',ERROR,*9999)

C *** Main Gauss point loop
      DO ng=1,NGT(NBH(NH_LOC(1,nx)))
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(/'' Gauss pt '',I3)') NG
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' ------------''/)')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF

C!!! requires a different reference state than the ventricular wall
CC ***   Interpolate Gauss pt geometric var.s XG and derivs wrt Xi
C        CALL XEXG(NBJ,ng,nr,PG,VE,XE,XG,ERROR,*9999)
CC ***   Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
CC ***   derivs (DXIX) of Xi wrt Xj (JP=0) coords.
C        CALL XGMG(0,NIT(NBJ(1)),NBJ(1),nr,DXIX,GXL,GXU,RGX(ng),XG,
C     '    ERROR,*9999)
C        RWG=RGX(ng)*WG(ng,NBH(NH_LOC(1,nx)))
C        IF(JTYP4.EQ.2) RWG=RWG*2.0d0*PI*XG(1,1)    !cyl symm about x
C        IF(JTYP4.EQ.3) RWG=RWG*2.0d0*PI*XG(2,1)    !cyl symm about y
C        IF(JTYP4.EQ.4) RWG=RWG*4.0d0*PI*XG(1,1)**2 !spherical symmetry

C ***   Interpolate dependent var.s ZG and derivs wrt Xi from ZEREF
        CALL ZEZG(0,NBH,ng,NHE,nx,DXIX,PG,ZEREF,ZG,ERROR,*9999)
C ***   Calculate deformed metric tensors wrt Xi (GZL,GZU)
        CALL ZGMG(NBH(NH_LOC(1,nx)),nr,GZ,GZL,GZU,ZG,ERROR,*9999)
        RGZREF=DSQRT(GZ)
        RWGZREF=RGZREF*WG(ng,NBH(NH_LOC(1,nx)))
        IF(JTYP4.EQ.2) RWGZREF=RWGZREF*2.0d0*PI*ZG(1,1)    !cyl symm (x)
        IF(JTYP4.EQ.3) RWGZREF=RWGZREF*2.0d0*PI*ZG(2,1)    !cyl symm (y)
        IF(JTYP4.EQ.4) RWGZREF=RWGZREF*4.0d0*PI*ZG(1,1)**2 !sph symmetry
        press_init=ZG(NH_LOC(0,nx),1) !initial cav. press.

C ***   Interpolate dependent var.s ZG and derivs wrt Xi from ZE
        CALL ZEZG(0,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
C ***   Calculate deformed metric tensors wrt Xi (GZL,GZU)
        CALL ZGMG(NBH(NH_LOC(1,nx)),nr,GZ,GZL,GZU,ZG,ERROR,*9999)
        RGZ=DSQRT(GZ)
        RWGZ=RGZ*WG(ng,NBH(NH_LOC(1,nx)))
        IF(JTYP4.EQ.2) RWGZ=RWGZ*2.0d0*PI*ZG(1,1)    !cyl symm about x
        IF(JTYP4.EQ.3) RWGZ=RWGZ*2.0d0*PI*ZG(2,1)    !cyl symm about y
        IF(JTYP4.EQ.4) RWGZ=RWGZ*4.0d0*PI*ZG(1,1)**2 !spherical symmetry
        press_current=ZG(NH_LOC(0,nx),1) !current cav. press.

C ***   pressure-displ constraints on lambda,mu,theta
        kstif=CG(1,ng) !GENERALISE!!!
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' >>mu stiffness='',D12.4,'
     '      //''' initial cavity pressure='',D12.4,'
     '      //''' current cavity pressure='',D12.4,'
     '      //''' RWGZREF='',D12.4,'' RWGZ='',D12.4)')
     '      kstif,press_init,press_current,RWGZREF,RWGZ
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Residual calcs:'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF !DOP
        DO nhx=1,NJ_LOC(NJL_GEOM,0,nr) !loop over geom vars
          nh=NH_LOC(nhx,nx)
          nb=NBH(nh) !basis fn for nh
          DO IXI3=1,2
            DO IXI2=1,2
              IF(IXI2.EQ.1) THEN
                press_sign= 1.0d0
              ELSE IF(IXI2.EQ.2) THEN
                press_sign=-1.0d0
              ENDIF
              DO IXI1=1,2
                ns=0
                nn=1
                FOUNDnn=.FALSE.
                DO WHILE(.NOT.FOUNDnn.AND.nn.LE.NNT(nb))
                  IF(INP(nn,1,nb).EQ.IXI1.AND.INP(nn,2,nb).EQ.IXI2
     '              .AND.INP(nn,3,nb).EQ.IXI3) THEN
                    FOUNDnn=.TRUE.
                    DO nk=1,NKT(nn,nb)
                      ns=ns+1
                      deltaMU=ZE(ns,nhx)-ZEREF(ns,nhx) !delta MU;deriv
                      RE(ns,nh)=RE(ns,nh)+(kstif*deltaMU+press_sign*
     '                  (press_current-press_init))*
     '                  PG(ns,1,ng,nb)*RWGZREF
                      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                      call mp_setlock()
                        WRITE(OP_STRING,'('' nh='',I2,'' ns='',I2,'
     '                    //''' deltaMU='',D12.4,'
     '                    //''' PG(ns,1,ng,nb):'',D12.4,'
     '                    //''' Resid='',D12.4)')
     '                    nh,ns,deltaMU,PG(ns,1,ng,nb),(kstif*deltaMU+
     '                    press_sign*(press_current-press_init))
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                      call mp_unsetlock()
                      ENDIF !DOP
                    ENDDO !nk
                  ELSE
                    DO nk=1,NKT(nn,nb)
                      ns=ns+1
                    ENDDO !nk
                    nn=nn+1
                  ENDIF
                ENDDO !while...
              ENDDO !IXI1
            ENDDO !IXI2
          ENDDO !IXI3
        ENDDO !nhx

        VOLUND=VOLUND+RWGZREF
        VOLDEF=VOLDEF+RWGZ

      ENDDO !ng

C *** Const volume constraint
      nb=NBH(NH_LOC(NH_LOC(0,nx),nx)) !basis fn for hyd. pressure
      DO ns=1,NST(nb)+NAT(nb)
        RE(ns,NH_LOC(NH_LOC(0,nx),nx))=VOLDEF-VOLUND
      ENDDO !ns

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(/'' Residuals for element'',I4,'':'')') ne
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nhx=1,NH_LOC(0,nx)
          nh=NH_LOC(nhx,nx)
          WRITE(OP_STRING,'('' RE(ns,'',I1,''): '',8D11.3,'
     '      //'/(11X,8D11.3))')
     '      nh,(RE(ns,nh),ns=1,NST(NBH(nh))+NAT(NBH(nh)))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

C MPN 28-Aug-95: old code for pressure/vol relationship
C      RE(1,1)=VOLDEF
C      RE(2,1)=VOLUND
C      DVOL=VOLDEF-VOLUND
C      IF(DOP) THEN
CC$      call mp_setlock()
C        WRITE(OP_STRING,'('' Volume change for element '','
C     '    //'I4,'' is '',D12.3)') NE,DVOL
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' Total volumes for element '','
C     '    //'I4,'' is '',D12.3,'' (deformed) and '',D12.3,'
C     '    //''' (undeformed)'')') NE,VOLDEF,VOLUND
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
C      ENDIF

      CALL EXITS('ZERE55')
      RETURN
 9999 CALL ERRORS('ZERE55',ERROR)
      CALL EXITS('ZERE55')
      RETURN 1
      END


