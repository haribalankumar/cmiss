      SUBROUTINE EVFLOW(NBJ,NDP,NEELEM,NELIST,NENP,NORD,NPLIST,
     &  NPNE,NPNODE,
     &  NRLIST,NVJE,NVJP,NXI,BBM,CE,XAB,XP,ZD,STRING,ERROR,*)

C#### Subroutine: EVFLOW
C###  Description:
C###    EVFLOW solves pulmonary 1D flow distributions. This subroutine
C###    does the problem set up.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn' !MEM_INIT
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn' 
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'mach00.inc' !DP_TYPE
      INCLUDE 'valu00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NDP(NDM),NEELEM(0:NE_R_M,0:NRM),
     &  NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),NORD(5,NE_R_M),
     &  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     &  NRLIST(0:NRM),
     &  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 BBM(2,NEM),CE(NMM,NEM,NXM),XAB(NORM,NEM),
     &  XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local variables
      REAL*8 constrict_by,COV,CW,dPmax,dPmin,dt,ERR_USER,
     &  FlowIN,FlowOUT,I_TO_E_RATIO,Pmax,Pmin,
     &  MeanCompliance,Ppl_step,dt_init,dt_max,PressureIN,refvol,
     &  RMaxMean,RMinMean,
     &  R1,R2,radius,T_interval,undef,volume_target
      INTEGER BELOW,C_ORDER,ERR,Gdirn,IBEG,IEND,IBEG1,IBEG2,IEND1,IEND2,
     &  IFROMC,IPFILE,ITER_USER,MECHANICS_FILETYPE,nb,NBREATHS,ne,
     &  ne0,NJJ_FLOW,NJJ_PRESSURE,NJJ_SV,NJJ_VRATIO,nn,noelem,np,nr,
     &  nsteps,N3CO,nv,nv1,nv2,N1LIST,order
      LOGICAL ALL_REGIONS,CALCU,COMPLIANCE_BC,DIAG_OP,NORMALISE,
     &  UNIFORM,CBBREV,FILIO,FIRST_ORDER,GENER,INITIAL,
     &  INLIST,MOUSE,PATHLENGTHS,PEAK,PRINT,READ_VOLUMES,SETUP
      CHARACTER EXTEND_FRC*1,FILE*(MXCH),filename*255,P_TYPE*50,STATUS*3
!     Functions
      REAL*8 RFROMC

      CALL ENTERS('EVFLOW',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM evaluate flow 
C###  Description:
C###    Solves pulmonary 1D flow distributions in an airway geometry.
C###  Parameter:      <UNIFORM>
C###    The flow distribution is proportional to the subtended
C###    volumes.
C###  Parameter:      <TARGET (#)>
C###    Set target tidal volume in mm3, rather than specifying 
C###    change in Ppl.
C###  Parameter:      <PPL_STEP (#) [-980.0]>
C###    Set change in Ppl over inspiration.
C###  Parameter:      <INTERVAL (#) [2.0]>
C###    Specify breath interval (inspiration + expiration).
C###  Parameter:      <TIME_STEP (#)>
C###    Specify time step.
C###  Parameter:      <NBREATHS (#) [5]>
C###    Specify number of breaths. 
C###  Parameter:      <SET_VOLUMES> 
C###                         <COV> <RMAX (#)[1]> <RMIN (#) [1]> 
C###                         <REFERENCE (#) [0.5]> <DIRECTION (#) [3]>
C###    Set initial distribution of acinar volume ratios (ratio of 
C###    deformed to undeformed volume) to be linear in the direction
C###    of gravity with random isogravitational heterogeneity.   
C###    COV = coefficient of variation. RMAX = maximum ratio.
C###    RMIN = minimum ratio. REFERENCE = ratio of reference volume to 
C###    start volume. DIRECTION=gravitational direction (x=1,y=2,z=3). 
C###  Parameter:      <READ_VOLUMES>
C###                         <UNDEFORMED (#) [1.0]>
C###                         <DATAPOINT/FIELD (#) [5]>
C###    Initial acinar volume ratio distribution is read from 
C###    ipdata file. UNDEFORMED = undeformed volume of an acinus 
C###    (generally = reference lung volume divided by #acini).
C###    Choose to read in either mechanics data from ipdata or ipfiel 
C###    file and set nj index (i.e. index of ZD or XP).
C###  Parameter:      <CHESTWALL (#) [2040]>
C###    Set chest wall compliance (L/cmH2O).
C###  Parameter:      <FILENAME (name)>
C###    Set output filename (for *.extime, *.shear, *.resist files).
C###  Parameter:      <INPRESSURE (#)>
C###    Specify pressure boundary condition at trachea.
C###  Parameter:      <CONSTRICT (#) [1.0]>
C###                         <RANDOM>
C###                         <BELOW (#)>
C###                         <ORDER (#)>
C###     Specify degree of airway constriction and distribution.
C###     RANDOM = random distribution. BELOW = constrict all airways 
C###     below the given Horsfield order. ORDER = constrict all airways 
C###     in the given Horsfield order.
        OP_STRING(1)=BLANK(1:15)//'<UNIFORM>'
        OP_STRING(2)=BLANK(1:15)//'<TARGET (#)>'
        OP_STRING(3)=BLANK(1:15)//'<PPL_STEP (#) [-980.0]>'
        OP_STRING(4)=BLANK(1:15)//'<INTERVAL (#) [2.0]>'
        OP_STRING(5)=BLANK(1:15)//'<TIME_STEP (#)>'
        OP_STRING(6)=BLANK(1:15)//'<NBREATHS (#) [5]>'
        OP_STRING(7)=BLANK(1:15)//'<SET_VOLUMES>'
        OP_STRING(8)=BLANK(1:20)//'<COV>'
        OP_STRING(9)=BLANK(1:20)//'<RMAX (#)[1]>'
        OP_STRING(10)=BLANK(1:20)//'<RMIN (#) [1]>'
        OP_STRING(11)=BLANK(1:20)//'<REFERENCE (#) [0.5]>'
        OP_STRING(12)=BLANK(1:20)//'<DIRECTION (#) [3]>'
        OP_STRING(13)=BLANK(1:15)//'<READ_VOLUMES>'
        OP_STRING(14)=BLANK(1:20)//'<UNDEFORMED (#) [1.0]>'
        OP_STRING(15)=BLANK(1:20)//'<DATAPOINT/FIELD (#)>'
        OP_STRING(16)=BLANK(1:15)//'<CHESTWALL (#) [2040]>'
        OP_STRING(17)=BLANK(1:15)//'<FILENAME (name)>'
        OP_STRING(17)=BLANK(1:15)//'<INPRESSURE (#)>'
        OP_STRING(18)=BLANK(1:15)//'<CONSTRICT (#) [1.0]>'
        OP_STRING(20)=BLANK(1:20)//'<RANDOM>'
        OP_STRING(21)=BLANK(1:20)//'<BELOW (#)>'
        OP_STRING(22)=BLANK(1:20)//'<ORDER (#)>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe31','doc','EVFLOW',ERROR,*9999)
      ELSE

      UNIFORM=.FALSE.
      SETUP=.FALSE.
      
      IF(CBBREV(CO,'UNIFORM',3,noco+1,NTCO,N3CO)) UNIFORM=.TRUE.
      IF(CBBREV(CO,'SETUP',3,noco+1,NTCO,N3CO)) SETUP=.TRUE.


      nr=1 !WARNING!!!!!!!!!!!! hardcoded to region 1

C read in a file, and/or get parameters from the command line
      CALL PARSE_QUALIFIERS('CDLPRW',noco,1,CO,COQU,
     '  CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
      CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
      IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
      CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'flow',
     '  STATUS,ERR,ERROR,*9999)

      CALL IPFLOW(Gdirn,ITER_USER,MECHANICS_FILETYPE,NBREATHS,NEELEM,
     &  NORD,NPLIST,NPNE,NVJE,nr,COV,CW,dt,dt_init,dt_max,ERR_USER,
     &  FlowIN,FlowOUT,I_TO_E_RATIO,MeanCompliance,Pmax,Pmin,
     &  Ppl_step,PressureIN,
     &  T_interval,
     &  refvol,RMaxMean,RMinMean,volume_target,undef,XP,DIAG_OP,
     &  NORMALISE,PATHLENGTHS,PRINT,READ_VOLUMES,UNIFORM,EXTEND_FRC,
     &  filename,
     &  P_TYPE,ERROR,*9999)

      CALL CLOSEF(IFILE,ERROR,*9999)

      PressureIN=PressureIN+PEEP

      IF(CBBREV(CO,'ZERO_FLOW',3,noco+1,NTCO,N3CO)) THEN
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)
      ELSE
        NELIST(0)=0
        NPLIST(0) = 0
      ENDIF

      IF(CBBREV(CO,'DEPOSITION',3,noco+1,NTCO,N3CO)) THEN
        UNIFORM=.TRUE. ! Super Hack, export particles
        T_interval=-1.0d0
      ENDIF

      IF(CBBREV(CO,'PEAK',3,noco+1,NTCO,N3CO)) PEAK=.TRUE.
      IF(CBBREV(CO,'CONSTRICT',3,noco+1,NTCO,N3CO)) THEN
        constrict_by=RFROMC(CO(N3CO+1))
        IF(CBBREV(CO,'BELOW',3,noco+1,NTCO,N3CO)) THEN
          BELOW=IFROMC(CO(N3CO+1))
        ENDIF
        IF(CBBREV(CO,'ORDER',3,noco+1,NTCO,N3CO)) THEN
          C_ORDER=IFROMC(CO(N3CO+1))
        ENDIF

        DO noelem=1,NEELEM(0,1)
          ne=NEELEM(noelem,1)
          order=NORD(2,ne) !Horsfield order
          IF(BELOW.GT.0)THEN
            IF(order.LE.BELOW)THEN !constrict
              DO nn=1,2
                np=NPNE(nn,1,ne) !which node at start and end of element
                nv=NVJE(nn,1,nj_radius,ne) !which version of node in the element
                XP(1,nv,nj_radius,np)=XP(1,nv,nj_radius,np)*constrict_by
              ENDDO !nn
            ENDIF !order
          ELSE IF(C_ORDER.GT.0)THEN
            IF(order.EQ.C_ORDER)THEN !constrict
              DO nn=1,2
                np=NPNE(nn,1,ne) !which node at start and end of element
                nv=NVJE(nn,1,nj_radius,ne) !which version of node in the element
                XP(1,nv,nj_radius,np)=XP(1,nv,nj_radius,np)*constrict_by
              ENDDO !nn
            ENDIF !order
          ENDIF
        ENDDO !noelem
      ENDIF !constrict


      CALL EVFLOW_DYNAM(Gdirn,ITER_USER,MECHANICS_FILETYPE,NBJ,
     &  NBREATHS,NDP,
     &  NEELEM(0,1),NENP(1,0,1),NORD,NPLIST,NPNE,NPNODE(0,1),NVJE,
     &  NVJP,NXI,COV,CW,dt,dt_init,dt_max,ERR_USER,FlowIN,FlowOUT,
     &  I_TO_E_RATIO,MeanCompliance,Pmax,Pmin,Ppl_step,PressureIN,
     &  T_interval,refvol,RMaxMean,RMinMean,
     &  volume_target,BBM,CE(1,1,1),undef,XAB,XP,ZD,COMPLIANCE_BC,
     &  DIAG_OP,FIRST_ORDER,INITIAL,
     &  NORMALISE,PATHLENGTHS,PEAK,PRINT,READ_VOLUMES,SETUP,
     &  UNIFORM,EXTEND_FRC,filename,P_TYPE,ERROR,*9999)
      

      ENDIF !CO='?'
      
      CALL EXITS('EVFLOW')
      RETURN
 9999 CALL ERRORS('EVFLOW',ERROR)
      CALL EXITS('EVFLOW')
      RETURN 1
      END
      
