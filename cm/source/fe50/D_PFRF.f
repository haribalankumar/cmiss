      SUBROUTINE D_PFRF(IBT,IDO,INP,IXF,NAN,NBH,NBJ,NBJF,NGAP,NHE,
     '  nhx1,nr,ns1,nx,D_RF,D_ZE,D_ZG,PF,PG,WG,XE,XG,ZE,ZG,
     '  ZG1,ERROR,*)

C#### Subroutine: D_PFRF
C###  Description:
C###    D_PFRF evaluates derivatives wrt geometric variables of
C###    contributions to element residuals RE from the pressure PF(i)
C###    acting on the Xi(3)=IXF face (IXF=iface-1).

C**** Note: GZL & GZU are deformed state metric tensors wrt Xi.
C**** The following code was copied from PFRF and altered and should
C**** be kept in sych with PFRF.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  IXF,NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),NBJF(NJM),
     '  NGAP(NIM,NBM),NHE,nhx1,nr,ns1,nx
      REAL*8 D_RF(32,6),D_ZE(NSM,NHM),D_ZG(NHM,NUM),PF,
     '  PG(NSM,NUM,NGM,NBM),WG(NGM,NBM),XE(NSM,NJM),XG(NJM,NUM),
     '  ZE(NSM,NHM),ZG(NHM,NUM),ZG1(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,NBFF,ng,ng1,ng2,NGI1,NGI2,nh,nhx,NITB,ni,ns,NU1(0:3)
      REAL*8 D(5,5),DELTA_ZE,D_GZ,D_GZU(3,3),D_RGZ,D_RWG,D_SUM,
     '  DXIX(3,3),GZ,GZ1,GZL(3,3),GZL1(3,3),
     '  GZU(3,3),GZU1(3,3),RGZ,RWG,SUM,XI(3)
      DATA NU1/1,2,4,7/
      DATA D/5*0.0d0,-0.288675134594813d0,0.288675134594813d0,3*0.0d0,
     '       -0.387298334620741d0,0.0d0,0.387298334620741d0,2*0.0d0,
     '       -0.430568155797026d0,    -0.169990521792428d0,
     '        0.169990521792428d0,     0.430568155797026d0,0.0d0,
     '       -0.453089922969332d0,    -0.269234655052841d0,0.0d0,
     '        0.269234655052841d0,     0.453089922969332d0/
      DATA DELTA_ZE/1.0d-8/

      CALL ENTERS('D_PFRF',*9999)
      NITB=NIT(NBH(NH_LOC(1,nx)))

      DO i=1,3
        DO j=1,3
          DXIX(i,j)=0.0d0
        ENDDO
      ENDDO

      DO nhx=1,NJ_LOC(NJL_GEOM,0,nr) !Shouldn't this be NH_LOC ? AJP 18/4/97
        nh=NH_LOC(nhx,nx)
        CALL ASSERT(nh.LE.6,'>>Incr dim of D_RF in D_ZERE50, D_PFRF',
     '    ERROR,*9999)
        NBFF=NBJF(nhx)
        DO ns=1,NST(NBFF)+NAT(NBFF)
          D_RF(ns,nh)=0.0d0
        ENDDO !ns
      ENDDO !nhx

      XI(3)=DBLE(IXF)
      NGI1=NGAP(1,NBJF(1))
      NGI2=NGAP(2,NBJF(1))
      ng=0
      DO ng2=1,NGI2
        DO ng1=1,NGI1
          ng=ng+1
          XI(1)=0.5d0+D(ng1,NGI1)
          XI(2)=0.5d0+D(ng2,NGI2)
          CALL XEXW(0,IBT,IDO,INP,NAN,NBJ,nr,XE,XG,XI,ERROR,*9999)
          CALL ZEZW(0,0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIX,ZE,ZG,
     '      XI,ERROR,*9999)
          CALL ZGMG(NBH(NH_LOC(1,nx)),nr,GZ,GZL,GZU,ZG,ERROR,*9999)
          RGZ=DSQRT(GZ)
          RWG=RGZ*WG(ng,NBJF(1))
C         Calc D_ZG, deriv of def coords (ZG), wrt elem coords (ZE)
          D_ZE(ns1,nhx1)=1.0d0 !differentiating wrt ns1,nhx1 elem coord
          CALL ZEZW(0,0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIX,D_ZE,D_ZG,XI,
     '      ERROR,*9999)
          D_ZE(ns1,nhx1)=0.0d0
C         Perturb deformed element coords
          ZE(ns1,nhx1)=ZE(ns1,nhx1)+DELTA_ZE
          CALL ZEZW(0,0,IBT,IDO,INP,NAN,NBH,NHE,nr,nx,DXIX,ZE,ZG1,XI,
     '      ERROR,*9999)
C         Reset deformed element coords
          ZE(ns1,nhx1)=ZE(ns1,nhx1)-DELTA_ZE
          CALL ZGMG(NBH(NH_LOC(1,nx)),nr,GZ1,GZL1,GZU1,ZG1,ERROR,*9999)
C         Finite diff approx to derivatives
          DO I=1,3
            DO J=1,3
              D_GZU(i,j)=(GZU1(i,j)-GZU(i,j))/DELTA_ZE
            ENDDO
          ENDDO
          D_GZ=(GZ1-GZ)/DELTA_ZE
          D_RGZ=D_GZ/(2.d0*DSQRT(GZ))
          D_RWG=D_RGZ*WG(ng,NBJF(1))

C LC 25/2/97 archive section :
C     old MPN 28-Jun-1995: integrals now wrt deformed coords
          DO nhx=1,NJ_LOC(NJL_GEOM,0,nr) !Shouldn't this be NH_LOC ? AJP 18/4/97
            nh=NH_LOC(nhx,nx)
            NBFF=NBJF(nhx)
            SUM=0.0d0
            D_SUM=0.0d0
            DO ni=1,NITB
C MPN 28-Jun-1995: integrals now wrt deformed coords
              SUM=SUM+GZU(ni,3)*ZG(nhx,NU1(ni))
              D_SUM=D_SUM+D_GZU(ni,3)*ZG(nhx,NU1(ni))
     '          +GZU(ni,3)*D_ZG(nhx,NU1(ni))
C old       SUM=SUM+GZU(ni,3)*DZI(nhx,ni)
C old       D_SUM=D_SUM+D_GZU(ni,3)*DZI(nhx,ni)+GZU(ni,3)*D_DZI(nhx,ni)
            ENDDO
            DO ns=1,NST(NBFF)+NAT(NBFF)
              IF(IXF.EQ.1)THEN !positive face
                D_RF(ns,nh)=D_RF(ns,nh)-PF*PG(ns,1,ng,NBFF)
     '            *(D_SUM*RWG+SUM*D_RWG)
              ELSE !negative face
                D_RF(ns,nh)=D_RF(ns,nh)+PF*PG(ns,1,ng,NBFF)
     '            *(D_SUM*RWG+SUM*D_RWG)
              ENDIF
            ENDDO
          ENDDO
        ENDDO !ngi1
      ENDDO !ngi2

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO nhx=1,NJ_LOC(NJL_GEOM,0,nr) !Shouldn't this be NH_LOC ? AJP 18/4/97
          nh=NH_LOC(nhx,nx)
          NBFF=NBJF(nhx)
          WRITE(OP_STRING,
     '      '('' D_RF(ns,'',I2,''): '',5D12.4)')
     '      nh,(D_RF(ns,nh),NS=1,NST(NBFF)+NAT(NBFF))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('D_PFRF')
      RETURN
 9999 CALL ERRORS('D_PFRF',ERROR)
      CALL EXITS('D_PFRF')
      RETURN 1
      END


