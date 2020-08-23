      SUBROUTINE MECH_FLOW(IBT,IDO,INP,NAN,NBH,NBJ,
     '  NEELEM,NENQ,NEP,NHE,NHP,NKH,NKHE,NKJE,
     '  NPF,NPNE,NPNODE,NQNE,NQS,NQXI,nr_mech,nr_flow,NVHE,NVHP,NVJE,
     '  NW,nx_flow,nx_mech,NYNE,NYNP,NYNQ,
     '  CE,CG,CGE,CP,CQ,CURVCORRECT,FEXT,PG,RG,SE,XA,
     '  XIP,XP,XQ,YG,YP,YQ,ZA,ZP,ERROR,*)

C#### Subroutine: MECH_FLOW
C###  Description:
C###    MECH_FLOW calculates the pressure exerted on vessels embeded
C###    in a host  medium. This is done by transforming the 2nd
C###    Piola-Kirkhoff stress tensor into a local vessel coordinate
C###    system via the host element xi coordinate system. The local
C###    vessel coordinate system has one axis aligned with the
C###    vessel direction and the other two normal to the vessel wall.
C###    Thus the pressure becomes the average of these second two
C###    stress components (see Coronary Flow Mechanics-PhD thesis N. P.
C###    Smith for further details). The xi positions of each grid point
C###     and  the host mesh element number are a calculated where the
C###    option to calculate the xi positions is chosen in the .ipcoup
C###    file. This routine is parralised using simple scheduling over
C###    the grid point loop and thus for large problems with
C###    correspondingly large numbers of grid points is scalable up to
C###     32 processors.


      IMPLICIT NONE

      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'time02.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),NEP(NPM),NHE(NEM),NHP(NPM),
     '  NKH(NHM,NPM,NCM),NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NQS(NEQM),NQXI(0:NIM,NQSCM),NQNE(NEQM,NQEM),
     '  nr_mech,nr_flow,NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),
     '  NW(NEM,3),nx_flow,nx_mech,NYNQ(NHM,NQM,0:NRCM,NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),
     '  CP(NMM,NPM),CQ(NMM,NQM),
     '  CURVCORRECT(2,2,NNM,NEM),FEXT(NIFEXTM,NGM,NEM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XIP(NIM,NPM),
     '  XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),
     '  YG(NIYGM,NGM,NEM),YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER  ne_flow,NITOT, !SMAR009 22/12/98 ne,ne_curent,SCHEME
     '  noelem,nq,ny_art,ny_lambda_art,
     '  ny_lambda_vien,ny_vien
      REAL*8  REF(3)
c     SMAR009 22/12/98 dydNu(3,3),dydx(3,3),MAT_VECTOR(3,3),
c     '  DIR_DEFORM(3),DXIZN(3,3),DZNXI(3,3),
c     '  NODE1(3),NODE2(3),
c     '  PST(3),,RG2D,RGZ,RGZ2D,RM(3,3)
      INTEGER*4 XE_PTR,XG_PTR,ZE_PTR,ZG_PTR
      REAL ELAPSED_TIME,TIME_START1(1),TIME_STOP(1)
      LOGICAL ERROR_FLAG !SMAR009 22/12/98 FOUND


      CALL ENTERS('MECH_FLOW',*9999)


      CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '  nr_mech,NVHP,nx_mech,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
      NITOT=NJ_LOC(NJL_GEOM,0,nr_mech)
      REF(1)=1.0d0
      REF(2)=0.0d0
      REF(3)=0.0d0
C reference direction used to find orthogonal directions to element
      ERROR_FLAG=.FALSE.

      CALL ASSERT(NIQM.GE.10,'>>Increase NIQM',ERROR,*9999)
      CALL ASSERT(NMM.GE.8,'>>Increase NMM',ERROR,*9999)


      CALL CPU_TIMER(CPU_USER,TIME_START1)

      DO nq=NQR(1,nr_flow),NQR(2,nr_flow)
        ny_art=NYNQ(1,nq,0,nx_flow)
        ny_vien=NYNQ(4,nq,0,nx_flow)
        ny_lambda_art=NYNQ(3,nq,0,nx_flow)
        ny_lambda_vien=NYNQ(6,nq,0,nx_flow)
        YQ(ny_art,9,1,nx_flow)=YQ(ny_art,2,1,nx_flow)
        YQ(ny_vien,9,1,nx_flow)=YQ(ny_vien,2,1,nx_flow)
        YQ(ny_lambda_art,9,1,nx_flow)=YQ(ny_lambda_art,2,1,nx_flow)
        YQ(ny_lambda_vien,9,1,nx_flow)=YQ(ny_lambda_vien,2,1,nx_flow)
        YQ(ny_art,2,1,nx_flow)=0.0d0
        YQ(ny_vien,2,1,nx_flow)=0.0d0
        YQ(ny_lambda_art,2,1,nx_flow)=0.0d0
        YQ(ny_lambda_vien,2,1,nx_flow)=0.0d0
      ENDDO

C stores the values of the trace from the previous time step in the
C YQ array

cC$OMP PARALLEL DO
cC$&     PRIVATE(noelem,ne_flow,XE_PTR,XG_PTR,ZE_PTR,ZG_PTR)
cC$&     SCHEDULE(GUIDED)


      DO noelem=1,NEELEM(0,nr_flow)
C parralises by looping over the  elements in the flow region
        IF(.NOT.ERROR_FLAG) THEN
          ne_flow=NEELEM(noelem,nr_flow)

          XE_PTR=0
          XG_PTR=0
          ZE_PTR=0
          ZG_PTR=0

          CALL ALLOCATE_MEMORY(NSM*NJM,1,DPTYPE,XE_PTR,
     '      MEM_INIT,ERROR,*100)
          CALL ALLOCATE_MEMORY(NJM*NUM,1,DPTYPE,XG_PTR,
     '      MEM_INIT,ERROR,*100)
          CALL ALLOCATE_MEMORY(NSM*NHM,1,DPTYPE,ZE_PTR,
     '      MEM_INIT,ERROR,*100)
          CALL ALLOCATE_MEMORY(NHM*NUM,1,DPTYPE,ZG_PTR,
     '      MEM_INIT,ERROR,*100)

C allocates the memory for the arrays passed into MECH_FLOW_DYNAM

          CALL MECH_FLOW_DYNAM(IBT,IDO,INP,NAN,NBH,NBJ,
     '      NEELEM,NENQ,NEP,ne_flow,NHE,NKHE,NKJE,NITOT,
     '      NPF,NPNE,NQNE,NQS,NQXI,nr_mech,
     '      NVHE,NVJE,
     '      NW,nx_flow,nx_mech,NYNQ,
     '      CE,CG,CGE,CP,CQ,CURVCORRECT,FEXT,PG,
     '      REF,RG,SE,XA,%VAL(XE_PTR),
     '      %VAL(XG_PTR),XIP,XP,XQ,YG,YQ,ZA,%VAL(ZE_PTR),%VAL(ZG_PTR),
     '      ZP,ERROR,*100)

          CALL FREE_MEMORY(XE_PTR,ERROR,*100)
          CALL FREE_MEMORY(XG_PTR,ERROR,*100)
          CALL FREE_MEMORY(ZE_PTR,ERROR,*100)
          CALL FREE_MEMORY(ZG_PTR,ERROR,*100)

C frees the memory

          GO TO 102
C             This statement is designed to be skipped if no error
C             occurs. However if a error occurs within a subroutine
C             the alternate return points to line 100 to set the flag
 100      CONTINUE
cC$OMP CRITICAL(MECH_FLOW)
         ERROR_FLAG=.TRUE.
          WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
          CALL WRITES(IOER,OP_STRING,ERROR,*101)
          WRITE(OP_STRING,'(/'' >>An error occurred - '
     '      //'results may be unreliable!'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101      CONTINUE
cC$OMP END CRITICAL(MECH_FLOW)
 102      CONTINUE

        ENDIF !error_flag
      ENDDO !noelem
cC$OMP END PARALLEL DO

      CALL ASSERT(.NOT.ERROR_FLAG,'>>An error occurred during '
     '  //'trace calculations',ERROR,*9999)

      CALL CPU_TIMER(CPU_USER,TIME_STOP)

      ELAPSED_TIME=TIME_STOP(1)-TIME_START1(1)

      WRITE(OP_STRING,'(/''total CPU time for trace '
     '  //'calculation: '',D15.7,'' s'')') ELAPSED_TIME
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)

      CALL EXITS('MECH_FLOW')
      RETURN
 9999 CALL ERRORS('MECH_FLOW',ERROR)
      CALL EXITS('MECH_FLOW')
      RETURN 1
      END



