      SUBROUTINE TRUSS(NBH,NBJ,NHE,nr,NW,nx,
     '  CE,CG,ED,EM,ER,ES,PG,RG,SE,WG,XE,XG,ERROR,*)

C#### Subroutine: TRUSS
C###  Description:
C###    TRUSS calculates the element load vector ER, the stiffness
C###    matrix ES and, if ITYP5(nr,nx)>1, the mass matrix EM and the
C###    damping matrix ED for the truss element (NW=1).

C**** Note: SUM1 accumulates the axial or flexural stiffness matrix
C****       SUM2     "           geometric            "        "
C****       SUM3     "           mass matrix.
C**** RT(nh,nhg) rotates variables nh=1,nhe in the element coor system
C****   to equivalent variables nhg in the global system.
C**** Direction cosines of truss axis are calc.d from truss geometry.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'b15.cmn'
      INCLUDE 'b40.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM),NBJ(NJM),NHE,nr,NW,nx
      REAL*8 CE(NMM),CG(NMM,NGM),ED(NHM*NSM,NHM*NSM),
     '  EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),ES(NHM*NSM,NHM*NSM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM),
     '  WG(NGM,NBM),XE(NSM,NJM),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mb,mh,mhg,mhs,mhsg,ms,MSTOT,mvar,
     '  nb,ng,nh,nhg,nhs,nhsg,nhx,nj,NLOAD,nm,ns,NSTOT,nvar
      REAL*8 AREA,DENS,DX(3),DXIX(3,3),GL(3,3),GU(3,3),
     '  RT(6,6),RWG,SUM,SUM1,SUM2,SUM3,XX,XX2,YMOD

      CALL ENTERS('TRUSS',*9999)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,
     '    '('' NVE='',I1,'' NHE='',I1)') NVE(NW),NHE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' CE(nm): '',10E12.3)') (CE(nm),nm=1,10)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

C *** Calculate rotation matrix RT
      SUM=0.0d0
      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
        nb=NBJ(nj)
        DX(nj)=XE(1+(NNT(nb)-1)*NKT(0,nb),nj)-XE(1,nj)
        SUM=SUM+DX(nj)*DX(nj)
      ENDDO
      XX=DSQRT(SUM) !is length of truss
      XX2=SUM
      DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
        RT(1,nj)=DX(nj)/XX
      ENDDO
      IF(NJ_LOC(NJL_GEOM,0,nr).EQ.1) THEN
        RT(1,1)=1.d0
        RT(1,2)=0.d0
        RT(2,1)=0.d0
        RT(2,2)=1.d0
        IF(NHE.EQ.1) THEN
          RT(2,1)=1.d0
        ELSE IF(NHE.EQ.2) THEN
          RT(3,1)=1.d0
          RT(4,2)=1.d0
          RT(3,2)=0.d0
          RT(4,1)=0.d0
        ENDIF
      ELSE IF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN
        RT(2,1)=-RT(1,2)
        RT(2,2)= RT(1,1)
      ELSE IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN
        DO nj=1,3
          RT(2,nj)=CE(ILT(NW,nr,nx)+nj)
          RT(3,nj)=CE(ILT(NW,nr,nx)+nj+3)
        ENDDO
      ENDIF
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Truss axis  direction cosines:'',3F7.4)')
     '    (RT(1,nj),nj=1,NJ_LOC(NJL_GEOM,0,nr))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      DO ng=1,NGT(NBJ(1)) !To determine RG(ng)
        CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
        CALL XGMG(0,NIT(NBJ(1)),NBJ(1),nr,DXIX,GL,GU,
     '    RG(ng),XG,ERROR,*9999)
      ENDDO
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' RG(ng): '',9E10.3)')
     '    (RG(ng),ng=1,NGT(NBJ(1)))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      NLOAD=NMP(NW)+NGP(NW)
      mhs=0
      DO mvar=1,NVE(NW)
        mh=NHV(mvar,NW)
        mb=NBH(mh)
        DO ms=1,NST(mb)
          mhs=mhs+1
          IF(mh.EQ.1) THEN
            SUM=0.d0
          ELSE
            SUM=0.d0
            DO ng=1,NGT(mb)
              RWG=RG(ng)*WG(ng,mb)
              SUM=SUM+CG(NLOAD+mh,ng)*PG(ms,1,ng,mb)*RWG
            ENDDO
          ENDIF
          SUM=SUM*SE(ms,mb)
          MSTOT=0
          DO mhg=1,NHE
            mhsg=ms+MSTOT
            MSTOT=MSTOT+NST(NBH(NH_LOC(1,nx)))
            ER(mhsg)=ER(mhsg)+RT(mh,mhg)*SUM
          ENDDO
          nhs=0
          DO nvar=1,NVE(NW)
            nhx=NHV(nvar,NW)
            nh=NH_LOC(nhx,nx)
            nb=NBH(nh)
            DO ns=1,NST(nb)
              nhs=nhs+1
              SUM1=0.0d0
              SUM2=0.0d0
              SUM3=0.0d0
              DO ng=1,NGT(nb)
                RWG=RG(ng)*WG(ng,nb)
                AREA=CG(NMP(NW)+1,ng)
                YMOD=CG(1,ng)
                IF(ITYP5(nr,nx).EQ.1) THEN !static
                  DENS=0.d0
                ELSE
                  DENS=CG(ILT(NW,nr,nx),ng)*CG(NMP(NW)+1,ng)*1.0d-3
                ENDIF
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,
     '              '('' YMOD='',E12.3,'' AREA='',E12.3,'' DENS='','
     '              //'E12.3)') YMOD,AREA,DENS
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
                IF(mh.EQ.1.AND.nh.EQ.1) THEN !axial component
                  SUM1=SUM1
     '                +YMOD*AREA*PG(ms,2,ng,MB)*PG(ns,2,ng,nb)*RWG/XX2
                  SUM2=0.d0 !Note what sb included here?
                  SUM3=SUM3+DENS*PG(ms,1,ng,MB)*PG(ns,1,ng,nb)*RWG
                ELSE IF(mh.EQ.NH) THEN       !transverse cpts
                ENDIF
              ENDDO
              SUM1=SUM1*SE(ms,mb)*SE(ns,nb)
              SUM2=SUM2*SE(ms,mb)*SE(ns,nb)
              SUM3=SUM3*SE(ms,mb)*SE(ns,nb)
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,
     '            '('' mvar='',I1,'' mh='',I1,'' ms='',I1,'
     '            //''' nvar='',I1,'' nh='',I1,'' ns='',I1)') mvar,
     '            mh,ms,nvar,nh,ns
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' sum1='',E11.3,'' sum2='',E11.3,'
     '            //''' sum3='',E11.3)') SUM1,SUM2,SUM3
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
              MSTOT=0
              DO mhg=1,NHE
                mhsg=ms+MSTOT
                MSTOT=MSTOT+NST(NBH(NH_LOC(1,nx)))
                NSTOT=0
                DO nhg=1,NHE
                  nhsg=ns+NSTOT
                  NSTOT=NSTOT+NST(NBH(NH_LOC(1,nx)))
c                 IF(ITYP5(nr,nx).EQ.1) THEN !static
c                   ES(mhsg,nhsg)=ES(mhsg,nhsg)+RT(mh,mhg)*(SUM1)*RT(nh,nhg)
c                 ELSE IF(ITYP5(nr,nx).EQ.2.OR.ITYP5(nr,nx).EQ.3) THEN
c                   ES(mhsg,nhsg)=ES(mhsg,nhsg)+RT(mh,mhg)*
c    '                                 (SUM1+CG(NLOAD+1,ng)*SUM2)*RT(nh,nhg)
c                   EM(mhsg,nhsg)=EM(mhsg,nhsg)+RT(mh,mhg)*(SUM3)*RT(nh,nhg)
c                 ELSE IF(ITYP5(nr,nx).EQ.4) THEN
c                 ELSE IF(ITYP5(nr,nx).EQ.5) THEN
c                 ENDIF
                  ES(mhsg,nhsg)=ES(mhsg,nhsg)+RT(mh,mhg)*SUM1
     '                                                  *RT(nh,nhg)
                  EM(mhsg,nhsg)=EM(mhsg,nhsg)+RT(mh,mhg)*SUM3
     '                                                  *RT(nh,nhg)
                  IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                    WRITE(OP_STRING,
     '                '('' ES('',I2,'','',I2,'')='',E11.3)')
     '                mhsg,nhsg,ES(mhsg,nhsg)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      IF(DAMPING_FACTOR1.GT.0.d0.OR.DAMPING_FACTOR2.GT.0.d0) THEN
        nb=NBH(NH_LOC(1,nx))
        DO ms=1,2*NST(nb)
          DO ns=1,2*NST(nb)
            ED(ms,ns)=DAMPING_FACTOR1*DAMPING_FREQUENCY*EM(ms,ns)
     '               +DAMPING_FACTOR2*DAMPING_FREQUENCY*ES(ms,ns)
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('TRUSS')
      RETURN
 9999 CALL ERRORS('TRUSS',ERROR)
      CALL EXITS('TRUSS')
      RETURN 1
      END


