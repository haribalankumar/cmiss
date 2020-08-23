      SUBROUTINE PLANE(NBH,NBJ,nr,NW,nx,CG,ED,EM,ES,PG,RG,SE,
     '  WG,XE,XG,ERROR,*)

C#### Subroutine: PLANE
C###  Description:
C###    PLANE calculates the element load vector ER and the stiffness
C###    matrix ES for plane stress (NW=11) and plane strain (NW=12).

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b15.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER NBH(NHM),NBJ(NJM),nr,NW,nx
      REAL*8 CG(NMM,NGM),ED(NHM*NSM,NHM*NSM),EM(NHM*NSM,NHM*NSM),
     '  ES(NHM*NSM,NHM*NSM),PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM),
     '  WG(NGM,NBM),XE(NSM,NJM),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ms,nb,ng,NGTB,ns,NSTB
      REAL*8 DENS,CC(3,3),DXIX(3,3),GL(3,3),GU(3,3),PM,PMX,PMY,PN,PNX,
     '  PNY,RWG

      CALL ENTERS('PLANE',*9999)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Plane stress or strain'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      nb=NBH(NH_LOC(1,nx)) !is dependent variable basis
      NGTB=NGT(nb)
      NSTB=NST(nb)

      DO ng=1,NGTB
        CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
        CALL XGMG(0,NIT(NBJ(1)),NBJ(1),nr,DXIX,GL,GU,
     '    RG(ng),XG,ERROR,*9999)
        RWG=RG(ng)*WG(ng,nb)
        CALL CGC11(nr,NW,nx,CC,CG(1,ng),DENS,ERROR,*9999)
        DO ms=1,NSTB
          DO ns=1,NSTB
            PM =PG(ms,1,ng,nb)
            PN =PG(ns,1,ng,nb)
            PMX=PG(ms,2,ng,nb)*DXIX(1,1)+PG(ms,4,ng,nb)*DXIX(2,1)
            PMY=PG(ms,2,ng,nb)*DXIX(1,2)+PG(ms,4,ng,nb)*DXIX(2,2)
            PNX=PG(ns,2,ng,nb)*DXIX(1,1)+PG(ns,4,ng,nb)*DXIX(2,1)
            PNY=PG(ns,2,ng,nb)*DXIX(1,2)+PG(ns,4,ng,nb)*DXIX(2,2)

            EM(ms,ns)          =EM(ms,ns)          +DENS*PM*PN*RWG
            EM(ms+NSTB,ns+NSTB)=EM(ms+NSTB,ns+NSTB)+DENS*PM*PN*RWG

            ES(ms,ns)          =ES(ms,ns)
     '        +(CC(1,1)*PMX*PNX+CC(3,3)*PMY*PNY)*RWG   !(u,u) term
            ES(ms,ns+NSTB)     =ES(ms,ns+NSTB)
     '        +(CC(1,2)*PMX*PNY+CC(3,3)*PMY*PNX)*RWG   !(u,v) term
            ES(ms+NSTB,ns)     =ES(ms+NSTB,ns)
     '        +(CC(2,1)*PMY*PNX+CC(3,3)*PMX*PNY)*RWG   !(v,u) term
            ES(ms+NSTB,ns+NSTB)=ES(ms+NSTB,ns+NSTB)
     '        +(CC(2,2)*PMY*PNY+CC(3,3)*PMX*PNX)*RWG   !(v,v) term

          ENDDO
        ENDDO
      ENDDO

      DO ms=1,NSTB
        DO ns=1,NSTB
          EM(ms,ns)          =EM(ms,ns)          *SE(ms,nb)*SE(ns,nb)
          EM(ms+NSTB,ns+NSTB)=EM(ms+NSTB,ns+NSTB)*SE(ms,nb)*SE(ns,nb)
          ES(ms,ns)          =ES(ms,ns)          *SE(ms,nb)*SE(ns,nb)
          ES(ms,ns+NSTB)     =ES(ms,ns+NSTB)     *SE(ms,nb)*SE(ns,nb)
          ES(ms+NSTB,ns)     =ES(ms+NSTB,ns)     *SE(ms,nb)*SE(ns,nb)
          ES(ms+NSTB,ns+NSTB)=ES(ms+NSTB,ns+NSTB)*SE(ms,nb)*SE(ns,nb)
        ENDDO
      ENDDO

      IF(DAMPING_FACTOR1.GT.0.d0.OR.DAMPING_FACTOR2.GT.0.d0) THEN
        DO ms=1,2*NSTB
          DO ns=1,2*NSTB
            ED(ms,ns)=DAMPING_FACTOR1*DAMPING_FREQUENCY*EM(ms,ns)
     '               +DAMPING_FACTOR2*DAMPING_FREQUENCY*ES(ms,ns)
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('PLANE')
      RETURN
 9999 CALL ERRORS('PLANE',ERROR)
      CALL EXITS('PLANE')
      RETURN 1
      END


