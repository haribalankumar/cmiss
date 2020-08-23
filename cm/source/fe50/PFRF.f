      SUBROUTINE PFRF(IBT,IDO,INP,IXF,NAN,NBH,NBJF,NGAP,nr,nx,
     '  PF,PG,RF,WG,ZE,ZG,ERROR,*)

C#### Subroutine: PFRF
C###  Description:
C###    PFRF evaluates contribution RF(ns,nj) to element residuals RE
C###    from the pressure PF(i) acting on the Xi(3)=IXF face
C###    (IXF=iface-1).

C**** Note: GZL & GZU are deformed state metric tensors wrt Xi.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  IXF,NAN(NIM,NAM,NBFM),NBH(NHM),NBJF(NJM),NGAP(NIM,NBM),nx
      REAL*8 PF,PG(NSM,NUM,NGM,NBM),RF(32,6),WG(NGM,NBM),
     '  ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,NBFF,ng,ng1,ng2,NGI1,NGI2,NITB,nh,nhx,ni,
     '  nr,ns,NU1(0:3)
      REAL*8 D(5,5),DXIX(3,3),GZ,GZL(3,3),GZU(3,3),RGZ,RWG,SUM,XI(3)
      DATA NU1/1,2,4,7/
      DATA D/5*0.0d0,-0.288675134594813d0,0.288675134594813d0,3*0.0d0,
     '       -0.387298334620741d0,0.0d0,0.387298334620741d0,2*0.0d0,
     '       -0.430568155797026d0,    -0.169990521792428d0,
     '        0.169990521792428d0,     0.430568155797026d0,  0.0d0,
     '       -0.453089922969332d0,    -0.269234655052841d0,  0.0d0,
     '        0.269234655052841d0,     0.453089922969332d0/

      CALL ENTERS('PFRF',*9999)
      NITB=NIT(NBH(NH_LOC(1,nx)))

      DO i=1,3
        DO j=1,3
          DXIX(i,j)=0.0d0
        ENDDO
      ENDDO

      DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
        nh=NH_LOC(nhx,nx)
        CALL ASSERT(nh.LE.6,'>>Incr dimension of RF in ZERE50, PFRF',
     '    ERROR,*9999)
        NBFF=NBJF(nhx)
        DO ns=1,NST(NBFF)+NAT(NBFF)
          RF(ns,nh)=0.0d0
        ENDDO !ns
      ENDDO !nhx

      XI(3)=DBLE(IXF)
      NGI1=NGAP(1,NBJF(1))
      NGI2=NGAP(2,NBJF(1))
      ng=0
      DO ng2=1,NGI2
        DO ng1=1,NGI1
          ng=ng+1
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,
     '        '('' >>>PFRF diagnostic op at Gauss pt '',I2)') NG
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF !DOP
          XI(1)=0.5d0+D(ng1,NGI1)
          XI(2)=0.5d0+D(ng2,NGI2)
c MPN 2Mar97 not needed
c          CALL XEXW(0,IBT,IDO,INP,NAN,NBJ,nr,XE,XG,XI,ERROR,*9999)
          CALL ZEZW(0,0,IBT,IDO,INP,NAN,NBH,NJ_LOC(NJL_GEOM,0,nr),nr,nx,
     '      DXIX,ZE,ZG,XI,ERROR,*9999)
          CALL ZGMG(NBH(NH_LOC(1,nx)),nr,GZ,GZL,GZU,ZG,ERROR,*9999)
          RGZ=DSQRT(GZ)
          RWG=RGZ*WG(ng,NBJF(1))
c MPN 2Mar97 not needed
c          IF(DOP) THEN
cC$          call mp_setlock()
c            DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
c              nj=NJ_LOC(NJL_GEOM,nhx,nr)
c              WRITE(OP_STRING,'(''  XG('',I1,'',ni): '',4D12.4)')
c     '          nj,(XG(nj,NU1(ni)),ni=0,NITB)
c              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c              WRITE(OP_STRING,'(''  ZG('',I1,'',ni): '',4D12.4)')
c     '          nhx,(ZG(nhx,NU1(ni)),ni=0,NITB)
c              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c            ENDDO !nhx
c            DO mi=1,NITB
c              WRITE(OP_STRING,'('' GZU('',I1,'',ni): '',3D12.4)')
c     '          MI,(GZU(mi,ni),ni=1,NITB)
c              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c            ENDDO !mi
cC$          call mp_unsetlock()
c          ENDIF

C LC 25/2/97 archived section :
C old MPN 28-Jun-1995: integrals now wrt deformed coords
          DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
            nh=NH_LOC(nhx,nx)
            NBFF=NBJF(nhx)
            SUM=0.0d0
            DO ni=1,NITB
C MPN 28-Jun-1995: integrals now wrt deformed coords
              SUM=SUM+GZU(3,ni)*ZG(nhx,NU1(ni))
C old              SUM=SUM+GZU(ni,3)*DZI(nhx,ni)
            ENDDO !ni
            DO ns=1,NST(NBFF)+NAT(NBFF)
              IF(IXF.EQ.1)THEN           !positive face
                RF(ns,nh)=RF(ns,nh)-PF*SUM*PG(ns,1,ng,NBFF)*RWG
              ELSE                       !negative face
                RF(ns,nh)=RF(ns,nh)+PF*SUM*PG(ns,1,ng,NBFF)*RWG
              ENDIF
            ENDDO !ns
          ENDDO !nhx
        ENDDO !ng1
      ENDDO !ng2

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
          nh=NH_LOC(nhx,nx)
          NBFF=NBJF(nhx)
          WRITE(OP_STRING,
     '      '('' RF(ns,'',I2,''): '',5D12.4,/(12X,5D12.4))')
     '      nh,(RF(ns,nh),ns=1,NST(NBFF)+NAT(NBFF))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nhx
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('PFRF')
      RETURN
 9999 CALL ERRORS('PFRF',ERROR)
      CALL EXITS('PFRF')
      RETURN 1
      END


