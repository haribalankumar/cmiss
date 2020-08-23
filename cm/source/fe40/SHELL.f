      SUBROUTINE SHELL(NBH,NBJ,NHE,nr,NW,nx,
     '  CE,CG,EM,ER,ES,GG,PG,RG,RR,SE,WG,XE,XG,ERROR,*)

C#### Subroutine: SHELL
C###  Description:
C###    SHELL calculates the element load vector ER, the stiffness
C###    matrix ES and, if ITYP5(nr,nx)>1, the mass matrix EM and the
C###    damping matrix ED for the shell element (NW=7).

C**** Note: SUM1 accumulates the STRAIN TENSOR stiffness matrix.
C****       SUM2      "          CHANGE OF CURVATURE    "        "
C****       SUM3      "          mass matrix.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'b40.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM),NBJ(NJM),NHE,nr,NW,nx
      REAL*8 CE(NMM),CG(NMM,NGM),EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),
     '  ES(NHM*NSM,NHM*NSM),GG(2,2,48,*),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  RR(2,2,48,*),SE(NSM,NBFM),WG(NGM,NBM),XE(NSM,NJM),
     '  XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mb,mh,mhs,mk,MMK,MMS,mn,ms,mvar,
     '  nb,ng,nh,nhs,nhx,nk,NLOAD,nm,nn,NNK,NNS,ns,nvar
      REAL*8 DENS,DXIX(3,3),GL(3,3),GU(3,3),P1,P2,P3,P4,POIS,
     '  RWG,SEM,SEN,SUM,SUM1,SUM2,SUM3,
     '  THIC,YMOD

      CALL ENTERS('SHELL',*9999)
      DO ng=1,NGT(NBJ(1)) !To determine RG(ng)
        CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
        CALL XGMG(0,NIT(NBJ(1)),NBJ(1),nr,DXIX,GL,GU,
     '    RG(ng),XG,ERROR,*9999)
      ENDDO

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' RG(ng): '',9E10.3)') (RG(ng),ng=1,9)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' CE(nm): '',14F6.2)') (CE(nm),nm=1,NMM)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' NVE='',I1,'' NHE='',I1)')
     '    NVE(NW),NHE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      SUM=0.d0
C  Initialise the rotation matrix RT using I
C     DO 4 nn=1,NNT(NBJ(1))
C     DO 4 M=1,3
C     DO 4 N=1,3
C       IF(M.EQ.N) THEN
C         RT(M,N,nn)=1.d0
C       ELSE
C         RT(M,N,nn)=0.d0
C       ENDIF
C4    CONTINUE
C
C  Finds the rotation matrix RT from the local coordinate system
C  to the global coordinate system
C<<<
C>>>  CALL XPRT(NBJ(1),NJ_LOC(NJL_GEOM,0,nr),RT,XE)
C===  CALL GUPUN(NBJ(1),NJ_LOC(NJL_GEOM,0,nr),PUE,XE)
C     IF(DOP) THEN
C       DO 101 nn=1,NNT(NBJ(1))
C       DO 101 I=1,NJ_LOC(NJL_GEOM,0,nr)
C         WRITE(IO4,201) (RT(I,J,nn),J=1,NJ_LOC(NJL_GEOM,0,nr))
C101    CONTINUE
C201    FORMAT(' ROTATION= ',9F7.2)
C     ENDIF

C *** Transform into the physical component to assemble

CCC   CALL GUPU(NBH,GU,PU)
      NLOAD=NMP(NW)+NGP(NW)
      mhs=0
      DO mvar=1,NVE(NW)
        mh=NHV(mvar,NW)
        mb=NBH(mh)
       ms=0
        DO mn=1,NNT(mb)
          DO mk=1,NKT(mn,mb)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,
     '          '('' mvar='',I1,'' mn='',I1,'' mk='',I1)')
     '           mvar,mn,mk
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
           ms=ms+1
            mhs=mhs+1
            SUM=0.d0
            DO ng=1,NGT(mb)
              RWG=RG(ng)*WG(ng,mb)
              SUM=SUM+CG(NLOAD+mh,ng)*PG(ms,1,ng,mb)*RWG
C 26185       SUM=SUM+CG(NLOAD+mh,ng)*PG(ms,1,ng,mb)*RWG/PU(mh,ng)
            ENDDO
C           ER(mhs)=SUM*SE(ms,mb)/PUE(mh,mn)
C
C    Note: Here, assuming bicubic Hermite for the geometric variables
C
            IF(mh.LE.2) THEN
              MMK=mh+1
C             MMS=(mn-1)*NKT(0,NBJ(1))+MMK
C             SEM=SE(MMS,NBJ(1))
              MMS=(mn-1)*NKT(0,NBH(NH_LOC(3,nx)))+MMK
              SEM=SE(MMS,NBH(NH_LOC(3,nx)))
            ELSE
              SEM=1.d0
            ENDIF
            ER(mhs)=SUM*SE(ms,mb)*SEM

C           SUM=SUM*SE(ms,mb)

C            MSTOT=0
C>>>        DO mhg=1,NHE
C>>>          mhsg=ms+MSTOT
C>>>          MSTOT=MSTOT+NST(NBH(NH_LOC(1,nx)))
C>>>          ER(mhsg)=ER(mhsg)+RT(mh,mhg,mn)*SUM
C           ENDDO
            nhs=0
            DO nvar=1,NVE(NW)
              nhx=NHV(nvar,NW)
              nh=NH_LOC(nhx,nx)
              nb=NBH(nh)
C             DO ns=1,NST(nb)
                ns=0
                DO nn=1,NNT(nb)
                  DO nk=1,NKT(nn,nb)
                    IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                    call mp_setlock()
                      WRITE(OP_STRING,'('' nvar='',I1,'' nn='',I1,'
     '                  //''' nk='',I1)') nvar,nn,nk
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                    call mp_unsetlock()
                    ENDIF
                    ns=ns+1
                    nhs=nhs+1
                    SUM1=0.d0
                    SUM2=0.d0
                    SUM3=0.d0
                    DO ng=1,NGT(nb)
                      RWG =RG(ng)*WG(ng,nb)
                      THIC=CG(NMP(NW)+1,ng)
                      YMOD=CG(1,ng)
                      POIS=CG(2,ng)
                      P1=YMOD*THIC*(1.d0-POIS)/(1.d0-POIS**2)
                      P2=YMOD*THIC*POIS/(1.d0-POIS**2)
                      P3=YMOD*THIC**3/(12.d0*(1.d0+POIS))
                      P4=YMOD*THIC**3*POIS/(12.d0*(1.d0-POIS**2))
                      IF(ITYP5(nr,nx).EQ.1) THEN
                        DENS=0.d0
                      ELSE
                        DENS=CG(ILT(NW,nr,nx),ng)*CG(NMP(NW)+1,ng)*1.D-3
                      ENDIF
                      SUM1=SUM1+(P1*(GG(1,1,mhs,ng)*GG(1,1,nhs,ng)
     '                              +GG(1,2,mhs,ng)*GG(2,1,nhs,ng)
     '                              +GG(2,1,mhs,ng)*GG(1,2,nhs,ng)
     '                              +GG(2,2,mhs,ng)*GG(2,2,nhs,ng))
     '                          +P2*(GG(1,1,mhs,ng)*GG(1,1,nhs,ng)
     '                              +GG(1,1,mhs,ng)*GG(2,2,nhs,ng)
     '                              +GG(2,2,mhs,ng)*GG(1,1,nhs,ng)
     '                              +GG(2,2,mhs,ng)*GG(2,2,nhs,ng)))*RWG
C261 '                             /PU(mh,ng)/PU(nh,ng)
                      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                      call mp_setlock()
                        WRITE(OP_STRING,'('' sum1='',E12.3)') SUM1
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                      call mp_unsetlock()
                      ENDIF
                      SUM2=SUM2+(P3*(RR(1,1,mhs,ng)*RR(1,1,nhs,ng)
     '                              +RR(1,2,mhs,ng)*RR(2,1,nhs,ng)
     '                              +RR(2,1,mhs,ng)*RR(1,2,nhs,ng)
     '                              +RR(2,2,mhs,ng)*RR(2,2,nhs,ng))
     '                          +P4*(RR(1,1,mhs,ng)*RR(1,1,nhs,ng)
     '                              +RR(1,1,mhs,ng)*RR(2,2,nhs,ng)
     '                              +RR(2,2,mhs,ng)*RR(1,1,nhs,ng)
     '                              +RR(2,2,mhs,ng)*RR(2,2,nhs,ng)))*RWG
C261 '                             /PU(mh,ng)/PU(nh,ng)
                      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                      call mp_setlock()
                        WRITE(OP_STRING,'('' sum2='',E12.3)') SUM2
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                      call mp_unsetlock()
                      ENDIF
                      SUM3=SUM3+DENS*PG(ms,1,ng,MB)*PG(ns,1,ng,nb)*RWG
                      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                      call mp_setlock()
                        WRITE(OP_STRING,'('' sum3='',E12.3)') SUM3
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                      call mp_unsetlock()
                      ENDIF
                    ENDDO
C                   SUM1=SUM1*SE(ms,mb)*SE(ns,nb)/(PUE(mh,mn)
C    '                                            *PUE(nh,nn))
C                   SUM2=SUM2*SE(ms,mb)*SE(ns,nb)/(PUE(mh,mn)
C    '                                            *PUE(nh,nn))
C                   SUM3=SUM3*SE(ms,mb)*SE(ns,nb)/(PUE(mh,mn)
C    '                                            *PUE(nh,nn))
C
C    Note: Here, assuming bicubic Hermite for the geometric variables
C
                    IF(nh.LE.2) THEN
                      NNK=nh+1
C                     NNS=(nn-1)*NKT(0,NBJ(1))+NNK
C                     SEN=SE(NNS,NBJ(1))
                      NNS=(nn-1)*NKT(0,NBH(NH_LOC(3,nx)))+NNK
                      SEN=SE(NNS,NBH(NH_LOC(3,nx)))
                    ELSE
                      SEN=1.d0
                    ENDIF
                    SUM1=SUM1*SE(ms,mb)*SE(ns,nb)*(SEM*SEN)
                    SUM2=SUM2*SE(ms,mb)*SE(ns,nb)*(SEM*SEN)
                    SUM3=SUM3*SE(ms,mb)*SE(ns,nb)*(SEM*SEN)
                    ES(mhs,nhs)=SUM1+SUM2
                    EM(mhs,nhs)=SUM3
C>>>                MSTOT=0
C>>>                DO mhg=1,NHE
C>>>                  mhsg=ms+MSTOT
C>>>                  MSTOT=MSTOT+NST(NBH(NH_LOC(1,nx)))
C>>>                  NSTOT=0
C>>>                  DO nhg=1,NHE
C>>>                    nhsg=ns+NSTOT
C>>>                    NSTOT=NSTOT+NST(NBH(NH_LOX(1,nx)))
C>>>                    IF(ITYP5(nr,nx).EQ.1) THEN
C>>>                      ES(mhsg,nhsg)=ES(mhsg,nhsg)
C>>> '                        +RT(mh,mhg,MN)*(    SUM1 + SUM2     )*RT(nh,nhg,nn)
C                       ELSE IF(ITYP5(nr,nx).EQ.2.OR.ITYP5(nr,nx).EQ.3) THEN
C                         ES(mhsg,nhsg)=ES(mhsg,nhsg)
C    '                        +RT(mh,mhg)*(SUM1+CG(NLOAD+1,ng)*SUM2)*RT(nh,nhg)
C                         EM(mhsg,nhsg)=EM(mhsg,nhsg)
C    '                        +RT(mh,mhg)*(        SUM3            )*RT(nh,nhg)
C                       ELSE IF(ITYP5(nr,nx).EQ.4) THEN
C                         ES(mhsg,nhsg)=ES(mhsg,nhsg)
C    '                        +RT(mh,mhg)*(        SUM1            )*RT(nh,nhg)
C                         EM(mhsg,nhsg)=EM(mhsg,nhsg)
C    '                        +RT(mh,mhg)*(        SUM2            )*RT(nh,nhg)
C                       ENDIF
C                     ENDDO
C                   ENDDO
C                   IF(DOP) THEN
C                     WRITE(IO4,151)         mvar,mh,MS,nvar,nh,ns,SUM1,SUM2,SUM3
C151                  FORMAT(' mvar=',I1,' mh=',I1,' ms=',I2,' nvar=',I1,' nh=',I1,
C    '                   ' ns=',I2,' sum1=',E10.3,' sum2=',E10.3,' sum3=',E10.3)
C                   ENDIF
C***              ES(ms,ns)=SUM1+SUM2
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      CALL EXITS('SHELL')
      RETURN
 9999 CALL ERRORS('SHELL',ERROR)
      CALL EXITS('SHELL')
      RETURN 1
      END


