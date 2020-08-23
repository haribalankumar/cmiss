      SUBROUTINE XPES40(NBH,NBJ,NHE,NKJE,NPF,NPNE,nr,
     '  NVJE,nw,nx,CE,CG,CGE,CP,ED,EM,ER,ES,PG,RG,SE,WG,
     '  XA,XE,XG,XP,YG,UPDATE,ERROR,*)

C#### Subroutine: XPES40
C###  Description:
C###    XPES40 calculates element load vector ER and (if UPDATE=.T.)
C###    element stiffness matrix ES for static (ITYP5(nr,nx)=1) or
C###    dynamic (ITYP5(nr,nx)=2) linear elasticity equations.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'dx00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grow00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp60.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER NBH(NHM),NBJ(NJM),NHE,NKJE(NKM,NNM,NJM),NPF(9,NFM),
     '  NPNE(NNM,NBFM),nr,NVJE(NNM,NBFM,NJM),NW,nx
      REAL*8 CE(NMM),CG(NMM,NGM),CGE(NMM,NGM),
     '  CP(NMM,NPM),ED(NHM*NSM,NHM*NSM),
     '  EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM),
     '  WG(NGM,NBM),XA(NAM,NJM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM)
      CHARACTER ERROR*(*)
      LOGICAL UPDATE
!     Local Variables
      INTEGER il,mh,mhs,mhsg,ms,nb,ng,nh,nhx,nhs,nhsg,ns
      REAL*8 GG(2,2,48,25),RR(2,2,48,25),SUM,X3G(4,3,25)

      CALL ENTERS('XPES40',*9999)

      DO mhsg=1,NHM*NSM
        ER(mhsg)=0.0d0
        IF(UPDATE) THEN
          DO nhsg=1,NHM*NSM
            ED(mhsg,nhsg)=0.0d0
            EM(mhsg,nhsg)=0.0d0
            ES(mhsg,nhsg)=0.0d0
          ENDDO
        ENDIF
      ENDDO

      CALL XPXE(NBJ,NKJE,NPF(1,1),NPNE,nr,NVJE,
     '  SE,XA,XE,XP,ERROR,*9999)

      IF(nw.EQ.7) THEN !shell element: find strain & change of curvature
        CALL X3XG(NBJ,D3PG,XE,X3G,ERROR,*9999)
        CALL PGGR(NBH,NBJ,nw,GG,PG,RG,RR,XE,XG,X3G,ERROR,*9999)
      ENDIF

      CALL CPCG(nw,NBJ(1),NPNE,nr,nx,CE,CG,CGE,CP,PG,ERROR,*9999)
      IF(KTYP60.EQ.1) THEN      !Growth law: so get density from YG(5)
        DO ng=1,NGT(NBJ(1))
          CG(5,ng)=YG(5,ng)                  !density
          CG(1,ng)=CE(1)*CG(5,ng)**2         !Young's mod from density
        ENDDO
      ELSE IF(KTYP60.EQ.2) THEN !Growth law: so get density from YG(5)
        DO ng=1,NGT(NBJ(1))
          CG(5,ng)=YG(5,ng)               !density
          CG(1,ng)=CE(1)*(CG(5,ng)/GROW1)**2 !Young's mod from density
        ENDDO
      ENDIF
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO ng=1,NGT(NBJ(1))
          WRITE(OP_STRING,
     '      '('' ng='',I1,'' CG:'',9E11.3)') ng,(CG(il,ng),il=1,9)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

      IF(UPDATE) THEN
        IF(nw.EQ.1) THEN      !truss or cable element
          CALL TRUSS(NBH,NBJ,NHE,nr,NW,nx,CE,CG,
     '      ED,EM,ER,ES,PG,RG,SE,WG,XE,XG,ERROR,*9999)
        ELSE IF(nw.EQ.2) THEN !batten element
          ERROR='>>Not implemented'
          GOTO 9999
        ELSE IF(nw.EQ.3) THEN !beam element
          CALL BEAM(NBH,NBJ,NHE,nw,CE,CG,
     '      EM,ER,ES,PG,RG,SE,WG,XE,XG,ERROR,*9999)
        ELSE IF(nw.EQ.5) THEN !membrane element
          ERROR='>>Not implemented'
          GOTO 9999
        ELSE IF(nw.EQ.6) THEN !plate element
          CALL PLATE(NBH,NBJ,NHE,nr,NW,nx,CE,CG,
     '      ER,ES,PG,RG,SE,WG,XE,XG,ERROR,*9999)
        ELSE IF(nw.EQ.7) THEN !shell element
          CALL SHELL(NBH,NBJ,NHE,nr,NW,nx,CE,CG,
     '      EM,ER,ES,GG,PG,RG,RR,SE,WG,XE,XG,ERROR,*9999)
        ELSE IF(nw.EQ.8) THEN !shell/fluid interface element
          ERROR='>>Not implemented'
          GOTO 9999
        ELSE IF(nw.EQ.9) THEN !3D element
          CALL ELAS3D(NBH,NBJ,nr,nw,nx,CG,
     '      ED,EM,ER,ES,PG,RG,SE,WG,XE,XG,ERROR,*9999)
C CPB 21/3/96 Adding in plane strain
C        ELSE IF(nw.EQ.11) THEN !plane stress element
C          CALL PLANE(NBH,NBJ,nr,nw,nx,CG,ED,EM,ES,PG,RG,SE,WG,XE,
C     '      XG,ERROR,*9999)
C        ELSE IF(nw.EQ.12) THEN !plane strain element
        ELSE IF(nw.EQ.11.OR.nw.EQ.12) THEN
          CALL PLANE(NBH,NBJ,nr,nw,nx,CG,ED,EM,ES,PG,RG,SE,WG,XE,
     '      XG,ERROR,*9999)
        ELSE IF(nw.EQ.13) THEN !fluid surface element
C         CALL AXISYM(NBH,NBJ,NHE,nw,CE,CG,
C    '      ED,EM,ER,ES,PG,RG,SE,WG,XE,ERROR,*9999)
        ENDIF

        IF(LUMP) THEN
          mhs=0
          DO mh=1,NHE
            DO ms=1,NST(NBH(mh))
              mhs=mhs+1
              nhs=0
              SUM=0.0d0
              DO nhx=1,NHE
                nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
                DO ns=1,NST(NBH(nh))
                  nhs=nhs+1
                  SUM=SUM+EM(mhs,nhs)
                  EM(mhs,nhs)=0.0d0
                ENDDO !ns
              ENDDO !nh
              EM(mhs,mhs)=SUM
            ENDDO !ms
          ENDDO !mh
        ENDIF
      ENDIF !lump

      IF(ITYP5(nr,nx).EQ.2) THEN
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          nb=NBH(NH_LOC(1,nx))
          WRITE(OP_STRING,'(/'' ER(mhs):''/,(1X,12E10.3))')
     '      (ER(mhs),mhs=1,NHE*NST(nb))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          IF(UPDATE) THEN
            DO mhs=1,NHE*NST(nb)
              WRITE(OP_STRING,
     '          '('' ES('',I2,'',nhs):''/,(1X,12E10.3))')
     '          mhs,(ES(mhs,nhs),nhs=1,NHE*NST(nb))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF
CC$        call mp_unsetlock()
        ENDIF
      ENDIF

      CALL EXITS('XPES40')
      RETURN
 9999 CALL ERRORS('XPES40',ERROR)
      CALL EXITS('XPES40')
      RETURN 1
      END


