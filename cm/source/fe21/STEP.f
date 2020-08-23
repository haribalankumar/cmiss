      SUBROUTINE STEP(IBT,ISCLOC,ISEG,ISLINE,ISLINO,NBH,NEELEM,
     '  NHE,NHP,NKH,NKJ,NLLIST,NPL,NPNODE,NVHP,NXLIST,NYNE,NYNP,
     '  DL,XP,YP,ZA,ZP,CSEG,STRING,ERROR,*)

C#### Subroutine: STEP
C###  Description:
C###    STEP steps through time varying nodal values.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'four00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'gks000.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'moti00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),ISCLOC(NWM),ISEG(*),ISLINE(NWM,2*NGRSEGM),
     '  ISLINO(NWM),NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKJ(NJM,NPM),
     '  NLLIST(0:NLM),NPL(5,0:3,NLM),
     '  NPNODE(0:NP_R_M,0:NRM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 DL(3,NLM),XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,INDEX,INDEX_POLYLINE,iw,IWK(6),N3CO,
     '  nb,nc,nh,nhx,NKEF,
     '  nks,no_coeff,no_coeffs,nocycle,nonode,nostep,np,nr,NT_COEFFS,
     '  NTCYCLE,NTIW,NTSTEP,nv,nx,nxc
      REAL*8 CLOCK,CYCLE_TIME,DISPLACEMENT,FLOW,RFROMC,PERIOD,PF1,TIME,
     '  TIME_INCREMENT
      LOGICAL CBBREV,OPFILE
      CHARACTER FILE*(MXCH)

      CALL ENTERS('STEP',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM step<;FILENAME>
C###  Parameter:    <cycle #CYCLES[1]>
C###    The number of cycles to complete
C###  Parameter:    <step #STEPS[10]>
C###    The number of steps in each cycle
C###  Parameter:    <update WSS[none]>
C###    Updates the GX graphics
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use
C###  Description:
C###    Steps through time-varying nodal soultions.

        OP_STRING(1)=STRING(1:IEND)//';<FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<cycle #CYCLES[1]>'
        OP_STRING(3)=BLANK(1:15)//'<step #STEPS[10]>'
        OP_STRING(4)=BLANK(1:15)//'<update WSS[none]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM step<;FILENAME>
C###  Parameter:    <from INITIAL_TIME#[T0]>
C###    The initial time
C###  Parameter:    <to FINAL_TIME#[T1]>
C###    The final time
C###  Parameter:    <update WSS#[none]>
C###    Updates the GX graphics
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use
C###  Description:
C###    Steps through time-varying nodal soultions.

        OP_STRING(1)=STRING(1:IEND)//';<FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<from INITIAL_TIME#[T0]>'
        OP_STRING(3)=BLANK(1:15)//'<to FINAL_TIME#[T1]>'
        OP_STRING(4)=BLANK(1:15)//'<update WSS[none]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','STEP',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IO4=IOFILE1
          CALL OPENF(IO4,'DISK',FILE(IBEG:IEND)//'.step','NEW',
     '     'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          IO4=IOOP
          OPFILE=.FALSE.
        ENDIF

        CALL ASSERT(CALL_MOTI,'>>Motion not defined',ERROR,*9999)
        CALL ASSERT(OMEGA.GT.0.d0,'>>Omega is zero',ERROR,*9999)

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        IF(CBBREV(CO,'CYCLE',1,noco+1,NTCO,N3CO)) THEN
          NTCYCLE=IFROMC(CO(N3CO+1))
        ELSE
          NTCYCLE=1
        ENDIF
        IF(CBBREV(CO,'STEP',1,noco+1,NTCO,N3CO)) THEN
          NTSTEP=IFROMC(CO(N3CO+1))
        ELSE
          NTSTEP=10
        ENDIF
        IF(CBBREV(CO,'FROM',1,noco+1,NTCO,N3CO)) THEN
          T0=RFROMC(CO(N3CO+1))
        ENDIF
        IF(CBBREV(CO,'TO',1,noco+1,NTCO,N3CO)) THEN
          T1=RFROMC(CO(N3CO+1))
        ENDIF
        IF(CBBREV(CO,'UPDATE',1,noco+1,NTCO,N3CO)) THEN
          UPVU=.TRUE.
          CALL PARSIL(CO(N3CO+1),4,NTIW,IWK,ERROR,*9999)
        ELSE
          UPVU=.FALSE.
        ENDIF

C CPB 22/11/94 Adding NX_LOC
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        PERIOD=2.0d0*PI/OMEGA
        TIME_INCREMENT=PERIOD/NTSTEP
        TIME=0.0d0
        nr=1 !temporarily

        IF(ITYP2(nr,nx).EQ.5.AND.ITYP3(nr,nx).EQ.2.
     '    AND.KTYP58(nr).EQ.3) THEN
          !Lung gas transport with Fourier basis for flow
          NT_COEFFS=IBT(2,NIT(NB_MOTION),NB_MOTION)
          DO nocycle=1,NTCYCLE !loops over cycles
            CYCLE_TIME=0.0d0
            DO nostep=1,NTSTEP !loops over steps within cycle
              CYCLE_TIME=CYCLE_TIME+TIME_INCREMENT
              TIME=TIME+TIME_INCREMENT
              FLOW=0.d0
              DO no_coeffs=1,NT_COEFFS
                FLOW=FLOW+PF1(no_coeffs,1,TIME)*FLOW_COEFFS(no_coeffs)
              ENDDO
              WRITE(OP_STRING,'('' Flow at time '',E12.3,'' is '','
     '          //'E12.3)') TIME,FLOW
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C              CALL MESH_FLOW(NEELEM,nr,NXI,CE(1,1,nx),VOL_CH,
C     '          ERROR,*9999)
              IF(UPVU) THEN
                iw=1
                CALL ACWK(iw,1,ERROR,*9999)
                IF(ISCLOC(iw).GT.0) THEN        !Update clock
                  CLOCK=CYCLE_TIME/PERIOD*2.d0*PI
                  CALL SGCLOC(0,ISCLOC(iw),ISEG,iw,CSEG,CLOCK,
     '              ERROR,*9999)
                ENDIF
                IF(ISLINE(iw,NTLINE).GT.0) THEN !Update line
                  INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','RED')
                  CALL SGLINE(index,ISEG,ISLINE(iw,NTLINE),ISLINO(iw),
     '              iw,NLLIST,NTLINE,NPL,nr,nx,CSEG,'UNDEFORMED',DL,
     '              'SOLID',.TRUE.,XP,ZP,ERROR,*9999)
                ENDIF
                CALL DAWK(iw,1,ERROR,*9999)
              ENDIF
            ENDDO
          ENDDO

        ELSE
          CALL YPZP(MOTION_IY,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '      NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '      YP(1,1,nx),ZA,ZP,ERROR,*9999)
          DO nocycle=1,NTCYCLE !loops over cycles
            CYCLE_TIME=0.0d0
            nc=1 !Temporary AJP 18-12-91
            DO nostep=1,NTSTEP !loops over steps within cycle
              CYCLE_TIME=CYCLE_TIME+TIME_INCREMENT
              TIME=TIME+TIME_INCREMENT
              DO nr=1,NRT
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  DO nhx=1,NHP(np,nr,nx)
                    nh=NH_LOC(nhx,nx)
                    nb=NBH(nh,nc,NEELEM(1,nr))     !Fourier basis number
                    DO nv=1,NVHP(nh,np,nc,nr)
                      DO nks=1,NKJ(nh,np)      !spatial dof number
C                       interpolate node np in time
                        DISPLACEMENT=0.d0
                        DO no_coeff=1,IBT(2,NIT(nb),nb)
                          NKEF=(no_coeff-1)*NKJ(nh,np)+nks !Fourier basis dof num.
                          DISPLACEMENT=DISPLACEMENT+
     '                     PF1(no_coeff,1,CYCLE_TIME)
     '                     *ZP(NKEF,nv,nh,np,nc)
                        ENDDO
                        IF(KTYP58(nr).EQ.2)THEN
                          ZP(nks,nv,nh,np,nc)=XP(nks,nv,nh,np)
     '                      +DISPLACEMENT
                        ELSE
                          ZP(nks,nv,nh,np,nc)=DISPLACEMENT
                        ENDIF
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

            nr=1 !temporarily
            IF(UPVU) THEN
              iw=1
              CALL ACWK(iw,1,ERROR,*9999)
              IF(ISCLOC(iw).GT.0) THEN        !Update clock
                CLOCK=CYCLE_TIME/PERIOD*2.d0*PI
                CALL SGCLOC(0,ISCLOC(iw),ISEG,iw,CSEG,CLOCK,ERROR,*9999)
              ENDIF
              IF(ISLINE(iw,NTLINE).GT.0) THEN !Update line
                CALL SGLINE(0,ISEG,ISLINE(iw,NTLINE),ISLINO(iw),iw,
     '            NLLIST,NTLINE,NPL,nr,nx,CSEG,'DEFORMED',
     '            DL,'SOLID',.TRUE.,XP,ZP,
     '            ERROR,*9999)
              ENDIF
              CALL DAWK(iw,1,ERROR,*9999)
            ENDIF
          ENDDO
        ENDIF
        IF(OPFILE) THEN
          CALL CLOSEF(IO4,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('STEP')
      RETURN
 9999 CALL ERRORS('STEP',ERROR)
      IF(OPFILE) CLOSE(UNIT=IO4)
      CALL EXITS('STEP')
      RETURN 1
      END


