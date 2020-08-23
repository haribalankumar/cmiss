      SUBROUTINE EVAERO(NBH,NEELEM,NHE,NHP,NKH,NKJ,NPL,NPNODE,
     '  NVHP,NYNE,NYNP,
     '  CE,DL,DLL,WG,XIG,XP,YP,ZA,ZP,STRING,ERROR,*)

C#### Subroutine: EVAERO
C###  Description:
C###    EVAERO evaluates aerofoil lift, drag and circulation and
C###    wake pressure difference.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'aero00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),NKJ(NJM,NPM),
     '  NPL(5,0:3,NLM),NPNODE(0:NP_R_M,0:NRM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM,NXM),DL(3,NLM),DLL(3,NLM),
     '  WG(NGM,NBM),XIG(NIM,NGM,NBM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER STRING*(MXCH),ERROR*(*)
!     Local Variables
      INTEGER i,IBEG,IEND,IFROMC,iXi1,N3CO,nb,nc,ne,ng,ni,nj,nl,NL1,NL2,
     '  no_aero,no_wake,NP1,NP2,nr,nv,nx,OP_LEVEL
      REAL*8 ARC_increm,ARC_length,DRAG,LIFT,X(3),Xi1,Xn(4,3),Z,Zn(4),
     '  PL1,PH3,PHI_WAKE(50),PRESS(50),PRESS_AERO(11,20),
     '  dSdXi(1),dXdXi(3,2),dXndXi(4,3),dZdXi(2),dZndXi(4),X_prev(3),
     '  PHI,PRESSURE,P_trail1,P_trail2,SUM,VELOCITY
      CHARACTER FILE*(MXCH),SURFACE*6
      LOGICAL CBBREV,OPFILE
      EXTERNAL PL1

      CALL ENTERS('EVAERO',*9999)

      nc=1 !temporary cpb 22/11/94
      nx=1 !temporary cpb 22/11/94

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate aerofoil<;FILENAME>
C###  Parameter:        <region #[1]>
C###    Specify the region numbers to evaluate.
C###  Parameter:        <output LEVEL#[0]>
C###    Specify a level of information output required
C###  Description:
C###    Evaluated the aerofoil lift, drag, circulations and
C###      wake presure difference
        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<output LEVEL#[0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVAERO',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opaero','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF

        IF(CBBREV(CO,'OUTPUT',1,noco+1,NTCO,N3CO)) THEN
          OP_LEVEL=IFROMC(CO(N3CO+1))
        ELSE
          OP_LEVEL=0
        ENDIF

        CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),
     '    NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,
     '    ERROR,*9999)

!Update line lengths (to allow for wake iterations and aerofoil shape
!changes - Should update DL(1) and DL(2)?
c       DO nl=1,NLT
c         CALL ARCLEN(IDO,NBJ,NEL(0,nl),nl,NPL(1,0,nl),NPNE,
c    '      DL,XP,ERROR,*9999)
c       ENDDO

c       DO nb=1,NBFT
c         IF(nb.EQ.NBH(NH_LOC(1,nx),1,NEELEM(1,nr))) THEN !reset -ve DLs
c           !Any -ve DL is reset +ve to avoid putting a -ve on SE
c           !(The -ve DL is used to allow a thick aerofoil to have
c           !increasing Xi(1) on upper and lower surfaces while
c           !having a continuous y derivative at the leading edge).
c           DO nl=1,NLT
c             DO i=1,3
c               DLL(i,nl)=DABS(DL(i,nl))
c             ENDDO
c           ENDDO
c           CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,NXI,DLL,
c    '        SE,ERROR,*9999)
c         ELSE
c           CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,NXI,DL,
c    '        SE,ERROR,*9999)
c         ENDIF
c       ENDDO
c Added 29May95 PJH
        DO nl=1,NLT
          DO i=1,3
            DLL(i,nl)=DABS(DL(i,nl))
          ENDDO
        ENDDO

! Aerofoil calculations
        IF(OP_LEVEL.GE.1) THEN
          WRITE(OP_STRING,'(/'' Aerofoil surface lines: '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

!     Calc. circulation from phi difference at trailing edge
        NL1=NL_AERO(NL_AERO(0,1),1) !last line on upper surface
        NL2=NL_AERO(NL_AERO(0,2),2) !last line on lower surface
        NP_aero_TE1=NPL(3,1,NL1) !is trailing edge node # on upper surface
        NP_aero_TE2=NPL(3,1,NL2) !is trailing edge node # on lower surface
        CIRCULATION=ZP(1,2,1,NP_aero_TE1,nc)-ZP(1,1,1,NP_aero_TE1,nc)
        IF(OP_LEVEL.GE.1) THEN
          WRITE(OP_STRING,'(/'' Circulation = '',E12.5)') CIRCULATION
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

!     Calc. trailing edge velocity diff from two versions at TE
        TE_VELOC_DIFF=ZP(2,2,1,NP_aero_TE1,nc)
     '    -ZP(2,1,1,NP_aero_TE1,nc)

!        IF(OP_LEVEL.GE.1) THEN
          WRITE(OP_STRING,'(/''   Circulation = '',E12.5)') CIRCULATION
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' TE veloc diff = '',E12.5)')TE_VELOC_DIFF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
!        ENDIF

        TOT_DRAG   = 0.d0
        TOT_LIFT   = 0.d0
        AERO_PERIM = 0.d0

        DO no_aero=1,NL_AERO(0,1)+NL_AERO(0,2)
          IF(no_aero.LE.NL_AERO(0,1)) THEN !upper surface
            SURFACE='UPPER'
            nl=NL_AERO(no_aero,1)
            ne=NE_AERO(no_aero,1)
          ELSE                             !lower surface
            SURFACE='LOWER'
            nl=NL_AERO(no_aero-NL_AERO(0,1),2)
            ne=NE_AERO(no_aero-NL_AERO(0,1),2)
          ENDIF

          NP1=NPL(2,1,nl) !is node# at Xi=0 end of line
          NP2=NPL(3,1,nl) !is node# at Xi=1 end of line

          REF_VELOC    =CE(1,ne,nx) !is free stream veloc in ne
          FLUID_DENSITY=CE(2,ne,nx) !is fluid density in ne

          DO nj=1,NJT
            Xn(1,nj)=XP(1,1,nj,NP1)
            Xn(2,nj)=XP(1,1,nj,NP2)
          ENDDO

          IF(SURFACE(1:5).EQ.'UPPER'.AND.NP2.EQ.NP_aero_TE1) THEN
            nv=2
          ELSE
            nv=1
          ENDIF

          Zn(1)=ZP(1,nv,1,NP1,nc)
          Zn(2)=ZP(1,nv,1,NP2,nc)

          DO nj=1,NJT
            nb=NPL(1,nj,nl)  !is type of basis for x(j)
            IF(nb.EQ.4) THEN !calc nodal derivs for cubic Hermite
              dXndXi(1,nj)=XP(NPL(4,1,nl),1,nj,NP1)*DL(1,nl)
              dXndXi(2,nj)=XP(NPL(5,1,nl),1,nj,NP2)*DL(2,nl)
            ENDIF
          ENDDO

          IF(NKH(NH_LOC(1,nx),NP1,1,nr).GT.1) THEN
            dZndXi(1)=ZP(2,nv,1,NP1,nc)*DLL(1,nl)
          ENDIF
          IF(NKH(NH_LOC(1,nx),NP2,1,nr).GT.1) THEN
            dZndXi(2)=ZP(2,nv,1,NP2,nc)*DLL(2,nl)
          ENDIF
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(/'' dZndXi(1)='',E12.3,'' dZndXi(2)='','
     '        //'E12.3)') dZndXi(1),dZndXi(2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF

          IF(OP_LEVEL.GE.2) THEN
            IF(no_aero.LE.NL_AERO(0,1)) THEN !upper surface
              WRITE(OP_STRING,'(/'' Line '',I4,'
     '          //''' Upper surface Nodes: '',2I5)') nl,NP1,NP2
            ELSE                             !lower surface
              WRITE(OP_STRING,'(/'' Line '',I4,'
     '          //''' Lower surface Nodes: '',2I5)') nl,NP1,NP2
            ENDIF
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

! Loop over aerofoil surface Xi coords in increments of 0.1
          ni=1 !temporary
          ng=0
          DRAG=0.d0
          LIFT=0.d0
C         DO Xi1=0.d0,1.d0,0.1d0
          DO iXi1=0,10
            Xi1=DBLE(i)/10.0d0

            SUM=0.d0
            DO nj=1,NJT
              nb=NPL(1,nj,nl)  !is type of basis for x(j)
              IF(nb.EQ.1) THEN      !linear Lagrange
                X(nj)       =PL1(1,1,Xi1)*Xn(1,nj)+PL1(2,1,Xi1)*Xn(2,nj)
                dXdXi(nj,ni)=PL1(1,2,Xi1)*Xn(1,nj)+PL1(2,2,Xi1)*Xn(2,nj)
              ELSE IF(nb.EQ.4) THEN !cubic Hermite
                X(nj)       =PH3(1,1,1,Xi1)*Xn(1,nj)
     '                      +PH3(2,1,1,Xi1)*Xn(2,nj)
     '                      +PH3(1,2,1,Xi1)*dXndXi(1,nj)
     '                      +PH3(2,2,1,Xi1)*dXndXi(2,nj)
                dXdXi(nj,ni)=PH3(1,1,2,Xi1)*Xn(1,nj)
     '                      +PH3(2,1,2,Xi1)*Xn(2,nj)
     '                      +PH3(1,2,2,Xi1)*dXndXi(1,nj)
     '                      +PH3(2,2,2,Xi1)*dXndXi(2,nj)
              ENDIF
              SUM=SUM+dXdXi(nj,ni)**2.d0
            ENDDO
            dSdXi(ni)=DSQRT(SUM)

            IF((no_aero.EQ.1.OR.no_aero.EQ.NL_AERO(0,1)+1)
     '        .AND.DABS(Xi1).LE.1.0D-6) THEN
              ARC_length = 0.d0
              DO nj=1,NJT
                X_prev(nj)=X(nj)
              ENDDO
            ELSE
              ARC_increm = 0.d0
              DO nj=1,NJT
                ARC_increm=ARC_increm+(X(nj)-X_prev(nj))**2.d0
                X_prev(nj)=X(nj)
              ENDDO
              ARC_increm=DSQRT(ARC_increm)
              ARC_length = ARC_length + ARC_increm
            ENDIF !no_aero

            IF(NKH(NH_LOC(1,nx),NP1,1,nr).EQ.1) THEN !linear Lagrange
              Z        =PL1(1,1,Xi1)*Zn(1)+PL1(2,1,Xi1)*Zn(2)
              dZdXi(ni)=PL1(1,2,Xi1)*Zn(1)+PL1(2,2,Xi1)*Zn(2)
            ELSE !cubic Hermite
              Z=         PH3(1,1,1,Xi1)*Zn(1)
     '                  +PH3(2,1,1,Xi1)*Zn(2)
     '                  +PH3(1,2,1,Xi1)*dZndXi(1)
     '                  +PH3(2,2,1,Xi1)*dZndXi(2)
              dZdXi(ni)= PH3(1,1,2,Xi1)*Zn(1)
     '                  +PH3(2,1,2,Xi1)*Zn(2)
     '                  +PH3(1,2,2,Xi1)*dZndXi(1)
     '                  +PH3(2,2,2,Xi1)*dZndXi(2)
            ENDIF !nkh

            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' x(1)='',E10.3,'' x(2)='',E10.3,'
     '          //''' dSdXi(1)='',E10.3)') X(1),X(2),dSdXi(1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' z='',E10.3,'' dZdXi(1)='',E10.3)')
     '          Z,dZdXi(1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF !dop

            PHI=Z                      !velocity potential
            VELOCITY=dZdXi(1)/dSdXi(1) !fluid velocity
            PRESSURE=.5d0*FLUID_DENSITY*(REF_VELOC**2.d0-VELOCITY**2.d0)

            IF(OP_LEVEL.GE.2) THEN
              WRITE(OP_STRING,'('' Xi1='',F4.2,'' s='',E10.3,'
     '          //''' Phi='',E10.3,'' Velocity = '',E10.3,'
     '          //''' Pressure = '',E10.3)')
     '          Xi1,ARC_length,PHI,VELOCITY,PRESSURE
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF

c Normal vector not reqd since project p onto dx for lift & dy for drag
c           CALL ORTHOG(NJT,dXdXi(1,1),dXdXi(1,2),NORMAL,ERROR,*9999)

c           IF(DOP) THEN
c             WRITE(OP_STRING,'('' normal vector: '',3E10.3)')
c    '          (NORMAL(nj),nj=1,NJT)
c             CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c           ENDIF

            IF(Xi1.LT.0.05d0.OR.Xi1.GT.0.95d0) THEN
c             LIFT=LIFT+0.5d0*PRESSURE*NORMAL(2)*dXdXi(1,1)*0.1d0
c             DRAG=DRAG+0.5d0*PRESSURE*NORMAL(1)*dXdXi(2,1)*0.1d0
              IF(SURFACE(1:6).EQ.'UPPER') THEN
                IF(NKJ(1,NP1).EQ.1) THEN !linear Lagrange
                  LIFT=LIFT-0.5d0*PRESSURE*dXdXi(1,1)*0.1d0
                  DRAG=DRAG+0.5d0*PRESSURE*dXdXi(2,1)*0.1d0
                ELSE !cubic Hermite (Xi increasing downstream)
c                 write(*,'('' upper surface cubic Hermite'')')
                  LIFT=LIFT-0.5d0*PRESSURE*dXdXi(1,1)*0.1d0
                  DRAG=DRAG+0.5d0*PRESSURE*dXdXi(2,1)*0.1d0
                ENDIF
              ELSE IF(SURFACE(1:6).EQ.'LOWER') THEN
                IF(NKJ(1,NP1).EQ.1) THEN !linear Lagrange
                  LIFT=LIFT+0.5d0*PRESSURE*dXdXi(1,1)*0.1d0
                  DRAG=DRAG-0.5d0*PRESSURE*dXdXi(2,1)*0.1d0
c               ELSE !cubic Hermite (Xi increasing upstream)
c                 write(*,'('' lower surface cubic Hermite'')')
c                 LIFT=LIFT-0.5d0*PRESSURE*dXdXi(1,1)*0.1d0
c                 DRAG=DRAG+0.5d0*PRESSURE*dXdXi(2,1)*0.1d0
                ELSE !cubic Hermite (Xi increasing downstream)
                  LIFT=LIFT+0.5d0*PRESSURE*dXdXi(1,1)*0.1d0
                  DRAG=DRAG-0.5d0*PRESSURE*dXdXi(2,1)*0.1d0
                ENDIF
              ENDIF
            ELSE
c             LIFT=LIFT  +  PRESSURE*NORMAL(2)*dXdXi(1,1)*0.1d0
c             DRAG=DRAG  +  PRESSURE*NORMAL(1)*dXdXi(2,1)*0.1d0
              IF(SURFACE(1:6).EQ.'UPPER') THEN
                IF(NKJ(1,NP1).EQ.1) THEN !linear Lagrange
                  LIFT=LIFT - PRESSURE*dXdXi(1,1)*0.1d0
                  DRAG=DRAG + PRESSURE*dXdXi(2,1)*0.1d0
                ELSE !cubic Hermite (Xi increasing downstream)
c                 write(*,'('' upper surface cubic Hermite'')')
                  LIFT=LIFT - PRESSURE*dXdXi(1,1)*0.1d0
                  DRAG=DRAG + PRESSURE*dXdXi(2,1)*0.1d0
                ENDIF
              ELSE IF(SURFACE(1:6).EQ.'LOWER') THEN
                IF(NKJ(1,NP1).EQ.1) THEN !linear Lagrange
                  LIFT=LIFT + PRESSURE*dXdXi(1,1)*0.1d0
                  DRAG=DRAG - PRESSURE*dXdXi(2,1)*0.1d0
c               ELSE !cubic Hermite (Xi increasing upstream)
c                 write(*,'('' lower surface cubic Hermite'')')
c                 LIFT=LIFT - PRESSURE*dXdXi(1,1)*0.1d0
c                 DRAG=DRAG + PRESSURE*dXdXi(2,1)*0.1d0
                ELSE !cubic Hermite (Xi increasing downstream)
                  LIFT=LIFT + PRESSURE*dXdXi(1,1)*0.1d0
                  DRAG=DRAG - PRESSURE*dXdXi(2,1)*0.1d0
                ENDIF
              ENDIF
            ENDIF
          ENDDO !iXi1

          IF(OP_LEVEL.GE.2) THEN
            WRITE(OP_STRING,'(/'' Drag = '',E12.5)') DRAG
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'( '' Lift = '',E12.5)') LIFT
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'( '' DL(3,nl) = '',E12.5)') DL(3,nl)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

          TOT_DRAG=TOT_DRAG+DRAG
          TOT_LIFT=TOT_LIFT+LIFT
          AERO_PERIM=AERO_PERIM+DL(3,nl)

! Loop over aerofoil Gauss point positions for stress calculations
          ni=1 !temporary
          DRAG=0.d0
          LIFT=0.d0

          IF(CALL_COUP) THEN !coupling to sail defined
            DO ng=1,NGT(NB_AERO_PRESS)
              Xi1=XIG(ni,ng,NB_AERO_PRESS)

              SUM=0.d0
              DO nj=1,NJT
                nb=NPL(1,nj,nl)  !is type of basis for x(j)
                IF(nb.EQ.1) THEN      !linear Lagrange
                  X(nj)       =PL1(1,1,Xi1)*Xn(1,nj)
     '                        +PL1(2,1,Xi1)*Xn(2,nj)
                  dXdXi(nj,ni)=PL1(1,2,Xi1)*Xn(1,nj)
     '                        +PL1(2,2,Xi1)*Xn(2,nj)
                ELSE IF(nb.EQ.4) THEN !cubic Hermite
                  X(nj)       =PH3(1,1,1,Xi1)*Xn(1,nj)
     '                        +PH3(2,1,1,Xi1)*Xn(2,nj)
     '                        +PH3(1,2,1,Xi1)*dXndXi(1,nj)
     '                        +PH3(2,2,1,Xi1)*dXndXi(2,nj)
                  dXdXi(nj,ni)=PH3(1,1,2,Xi1)*Xn(1,nj)
     '                        +PH3(2,1,2,Xi1)*Xn(2,nj)
     '                        +PH3(1,2,2,Xi1)*dXndXi(1,nj)
     '                        +PH3(2,2,2,Xi1)*dXndXi(2,nj)
                ENDIF
                SUM=SUM+dXdXi(nj,ni)**2.d0
              ENDDO
              dSdXi(ni)=DSQRT(SUM)

              IF(NKH(NH_LOC(1,nx),NP1,1,nr).EQ.1) THEN !linear Lagrange
                Z        =PL1(1,1,Xi1)*Zn(1)+PL1(2,1,Xi1)*Zn(2)
                dZdXi(ni)=PL1(1,2,Xi1)*Zn(1)+PL1(2,2,Xi1)*Zn(2)
              ELSE !cubic Hermite
                Z=         PH3(1,1,1,Xi1)*Zn(1)
     '                    +PH3(2,1,1,Xi1)*Zn(2)
     '                    +PH3(1,2,1,Xi1)*dZndXi(1)
     '                    +PH3(2,2,1,Xi1)*dZndXi(2)
                dZdXi(ni)= PH3(1,1,2,Xi1)*Zn(1)
     '                    +PH3(2,1,2,Xi1)*Zn(2)
     '                    +PH3(1,2,2,Xi1)*dZndXi(1)
     '                    +PH3(2,2,2,Xi1)*dZndXi(2)
              ENDIF

              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' x(1)='',E10.3,'' x(2)='',E10.3,'
     '            //''' dSdXi(1)='',E10.3)') X(1),X(2),dSdXi(1)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' z='',E10.3,'
     '            //''' dZdXi(1)='',E10.3)') Z,dZdXi(1)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF

              PHI=Z                      !velocity potential
              VELOCITY=dZdXi(1)/dSdXi(1) !fluid velocity
              PRESSURE=0.5d0*FLUID_DENSITY*(REF_VELOC**2.d0-
     '          VELOCITY**2.d0)

              IF(OP_LEVEL.GE.2) THEN
                WRITE(OP_STRING,'('' Xi1='',F4.2,'' Phi='',E10.3,'
     '            //''' Velocity = '',E10.3,'' Pressure = '',E10.3)')
     '            Xi1,PHI,VELOCITY,PRESSURE
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF

              LIFT=LIFT + PRESSURE*dXdXi(1,1)*WG(ng,NB_AERO_PRESS)
              DRAG=DRAG + PRESSURE*dXdXi(2,1)*WG(ng,NB_AERO_PRESS)

              !Record press at aerofoil Gauss points for stress analysis
              PRESS_AERO(ng,no_aero)=PRESSURE
            ENDDO !ng

            IF(OP_LEVEL.GE.2) THEN
              WRITE(OP_STRING,'(/'' Drag from Gauss points = '','
     '          //'E12.5)') DRAG
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'( '' Lift from Gauss points = '','
     '          //'E12.5)') LIFT
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF

          ENDIF !call_coup

        ENDDO !no_aero

        IF(DABS(TOT_LIFT).GT.1.0D-6) THEN
          LIFT_COEFF=0.5d0*FLUID_DENSITY*REF_VELOC**2.d0*AERO_PERIM
     '                /DABS(TOT_LIFT)
        ELSE
          LIFT_COEFF=0.d0
        ENDIF

        IF(OP_LEVEL.GE.1) THEN
          WRITE(OP_STRING,'(/'' Aerofoil drag = '',E12.5)') TOT_DRAG
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' Circumference = '',E12.5)') AERO_PERIM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' Lift coeff.   = '',E12.5)') LIFT_COEFF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        WRITE(OP_STRING,'( '' Aerofoil lift = '',E12.5)') TOT_LIFT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(CALL_COUP) THEN !coupling to sail defined
! Record pressure difference at Gauss points for stress analysis
          DO no_aero=1,NL_AERO(0,1)
            DO ng=1,NGT(NB_AERO_PRESS)
              PRESS_DIFF_AERO(ng,no_aero)=
     '                          PRESS_AERO(ng,no_aero+NL_AERO(0,1))
     '                         -PRESS_AERO(ng,no_aero)
            ENDDO
            IF(OP_LEVEL.GE.1) THEN
              WRITE(OP_STRING,'('' Press diff.s at Gauss pts '','
     '          //'''for aero line '',I5,'': '',7E12.5)')
     '          NL_AERO(no_aero,1),
     '          (PRESS_DIFF_AERO(ng,no_aero),ng=1,NGT(NB_AERO_PRESS))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
        ENDIF !call_coup

! Wake calculations
        IF(OP_LEVEL.GE.2) THEN
          WRITE(OP_STRING,'(/'' Wake surface lines: '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

        DO no_wake=1,NL_WAKE(0,1)
          nl=NL_WAKE(no_wake,1)

          NP1=NPL(2,1,nl) !is node # at Xi=0 end of line
          NP2=NPL(3,1,nl) !is node # at Xi=1 end of line

          DO nj=1,NJT
            Xn(1,nj)=XP(1,1,nj,NP1)
            Xn(2,nj)=XP(1,1,nj,NP2)
          ENDDO
          Zn(1)=ZP(1,1,1,NP1,nc)
          Zn(2)=ZP(1,1,1,NP2,nc)

          PHI_WAKE(no_wake)=Zn(2) !is veloc pot at Xi=1 end of line

          DO nj=1,NJT
            nb=NPL(1,nj,nl)  !is type of basis for x(j)
            IF(nb.EQ.4) THEN !calc nodal derivs for cubic Hermite
              dXndXi(1,nj)=XP(NPL(4,1,nl),1,nj,NP1)*DL(1,nl)
              dXndXi(2,nj)=XP(NPL(5,1,nl),1,nj,NP2)*DL(2,nl)
            ENDIF
          ENDDO

          IF(NKH(NH_LOC(1,nx),NP1,1,nr).GT.1) THEN
            dZndXi(1)=ZP(2,1,1,NP1,nc)*DLL(1,nl)
          ENDIF
          IF(NKH(NH_LOC(1,nx),NP2,1,nr).GT.1) THEN
            dZndXi(2)=ZP(2,1,1,NP2,nc)*DLL(2,nl)
          ENDIF

          IF(OP_LEVEL.GE.2) THEN
            WRITE(OP_STRING,'(/'' Line '',I4,'
     '          //''' Wake nodes: '',2I5)') nl,NP1,NP2
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

! Loop over surface Xi coords in increments of 0.1
          ni=1 !temporary
C         DO Xi1=0.d0,1.d0,0.1d0
          DO iXi1=0,10
            Xi1=DBLE(iXi1)/10.0d0

            SUM=0.d0
            DO nj=1,NJT
              nb=NPL(1,nj,nl)  !is type of basis for x(j)
              IF(nb.EQ.1) THEN      !linear Lagrange
                X(nj)       =PL1(1,1,Xi1)*Xn(1,nj)+PL1(2,1,Xi1)*Xn(2,nj)
                dXdXi(nj,ni)=PL1(1,2,Xi1)*Xn(1,nj)+PL1(2,2,Xi1)*Xn(2,nj)
              ELSE IF(nb.EQ.4) THEN !cubic Hermite
                X(nj)       =PH3(1,1,1,Xi1)*Xn(1,nj)
     '                      +PH3(2,1,1,Xi1)*Xn(2,nj)
     '                      +PH3(1,2,1,Xi1)*dXndXi(1,nj)
     '                      +PH3(2,2,1,Xi1)*dXndXi(2,nj)
                dXdXi(nj,ni)=PH3(1,1,2,Xi1)*Xn(1,nj)
     '                      +PH3(2,1,2,Xi1)*Xn(2,nj)
     '                      +PH3(1,2,2,Xi1)*dXndXi(1,nj)
     '                      +PH3(2,2,2,Xi1)*dXndXi(2,nj)
              ENDIF
              SUM=SUM+dXdXi(nj,ni)**2.d0
            ENDDO
            dSdXi(ni)=DSQRT(SUM)

            IF((no_wake.EQ.1).AND.DABS(Xi1).LE.1.0D-6) THEN
              ARC_length = 0.d0
              DO nj=1,NJT
                X_prev(nj)=X(nj)
              ENDDO
            ELSE
              ARC_increm = 0.d0
              DO nj=1,NJT
                ARC_increm=ARC_increm+(X(nj)-X_prev(nj))**2.d0
                X_prev(nj)=X(nj)
              ENDDO
              ARC_increm=DSQRT(ARC_increm)
              ARC_length = ARC_length + ARC_increm
            ENDIF

            IF(NKH(NH_LOC(1,nx),NP1,1,nr).EQ.1) THEN !linear Lagrange
              Z        =PL1(1,1,Xi1)*Zn(1)+PL1(2,1,Xi1)*Zn(2)
              dZdXi(ni)=PL1(1,2,Xi1)*Zn(1)+PL1(2,2,Xi1)*Zn(2)
            ELSE !cubic Hermite
              Z=         PH3(1,1,1,Xi1)*Zn(1)
     '                  +PH3(2,1,1,Xi1)*Zn(2)
     '                  +PH3(1,2,1,Xi1)*dZndXi(1)
     '                  +PH3(2,2,1,Xi1)*dZndXi(2)
              dZdXi(ni)= PH3(1,1,2,Xi1)*Zn(1)
     '                  +PH3(2,1,2,Xi1)*Zn(2)
     '                  +PH3(1,2,2,Xi1)*dZndXi(1)
     '                  +PH3(2,2,2,Xi1)*dZndXi(2)
            ENDIF

            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' x(1)='',E10.3,'' x(2)='',E10.3,'
     '          //''' dSdXi(1)='',E10.3)') X(1),X(2),dSdXi(1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' z='',E10.3,'' dZdXi(1)='',E10.3)')
     '          Z,dZdXi(1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF

            PHI=Z                      !velocity potential
            VELOCITY=dZdXi(1)/dSdXi(1) !fluid velocity
            PRESSURE=.5d0*FLUID_DENSITY*(REF_VELOC**2.d0-VELOCITY**2.d0)

            IF(OP_LEVEL.GE.2) THEN
              WRITE(OP_STRING,'('' Xi1='',F4.2,'' s='',E10.3,'
     '          //''' Phi='',E10.3,'' Velocity = '',E10.3,'
     '          //''' Pressure = '',E10.3)')
     '          Xi1,ARC_length,PHI,VELOCITY,PRESSURE
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF

            IF(XI1.GT.0.95d0) THEN !store pressure at Xi1=1
              PRESS(no_wake)=PRESSURE
            ENDIF
            IF(no_wake.EQ.1.AND.XI1.LT.0.05d0) THEN !trailing edge press
              P_trail1=PRESSURE
            ENDIF
          ENDDO !iXi1
        ENDDO !no_wake

        DO no_wake=1,NL_WAKE(0,1)
          !Calculate PHI difference across wake (lower-upper)
          dPHI_resid(no_wake)     =PHI_WAKE(no_wake+NL_WAKE(0,1))
     '                            -PHI_WAKE(no_wake)
          !Calculate PHI_diff - PHI_diff_Trailing_Edge
          dPHI_resid(no_wake)     =dPHI_resid(no_wake)-CIRCULATION
          WRITE(OP_STRING,'('' dPHI_residual '','
     '      //'''for wake line '',I5,'' is '',E12.5)')
     '      NL_WAKE(no_wake,1),dPHI_resid(no_wake)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          !Calculate pressure difference across wake
          PRESS_DIFF_WAKE(no_wake)=PRESS(no_wake+NL_WAKE(0,1))
     '                            -PRESS(no_wake)
          IF(OP_LEVEL.GE.1) THEN
            WRITE(OP_STRING,'('' Pressure difference at Xi1=1 '','
     '        //'''for wake line '',I5,'' is '',E12.5)')
     '        NL_WAKE(no_wake,1),PRESS_DIFF_WAKE(no_wake)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO

! Calc pressure difference at trailing edge
        p_trail2=0.0d0 !AJP 19/9/95 what should p_trail2 be (not set) ????
        PRESS_DIFF_WAKE(1)=P_trail1-P_trail2
        WRITE(OP_STRING,'('' TE press diff = '',E11.4)')
     '    PRESS_DIFF_WAKE(NL_WAKE(0,1)+1)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('EVAERO')
      RETURN
 9999 CALL ERRORS('EVAERO',ERROR)
      CALL EXITS('EVAERO')
      RETURN 1
      END


